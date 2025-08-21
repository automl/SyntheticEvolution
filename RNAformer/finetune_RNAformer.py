import os

os.environ['TORCH_DISTRIBUTED_DEBUG'] = 'DETAIL'
import sys
import argparse
import pathlib
import torch
import logging
import torch.cuda
import yaml
import random
import math
import numpy as np
import pandas as pd
import lightning as L
from lightning.fabric.loggers import CSVLogger, TensorBoardLogger
from tqdm import tqdm

from RNAformer.model.RNAformer import RNAFormer
from RNAformer.utils.configuration import Config, read_unknown_args
from RNAformer.utils.folder_manager import get_experiment_folder

logger = logging.getLogger(__name__)

from evaluate_RNAformer import prepare_RNA_sample, evaluate_RNAformer
from RNAformer.utils.data.rna import CollatorRNA

from RNAformer.utils.optim.adam_cpr_fast import AdamCPR, group_parameters_for_cpr_optimizer


def trainable_params(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def total_params(model):
    return sum(p.numel() for p in model.parameters())


def crop_mat(batch, random_crop_mat):
    for idx in range(batch['mask'].shape[0]):

        mat_size = batch['mask'].shape[1]

        if mat_size - random_crop_mat <= 1:
            continue

        crop_x = np.random.randint(0, mat_size - random_crop_mat)
        crop_y = np.random.randint(0, mat_size - random_crop_mat)

        crop_mask = torch.zeros_like(batch['mask'][idx])
        crop_mask[crop_x:crop_x + random_crop_mat, crop_y:crop_y + random_crop_mat] = 1

        batch['mask'][idx] = batch['mask'][idx] & crop_mask

    return batch


def ignore_partial_mat(batch, random_ignore_mat):
    trg_mat = batch['trg_mat']

    with torch.no_grad():
        filter = torch.nn.Conv2d(1, 1, 3, 1, 1, bias=False).to(trg_mat.device)
        torch.nn.init.constant_(filter.weight, 1.0)
        filter_trg = filter(trg_mat.unsqueeze(1).float()).squeeze(1)
        filter_trg = filter_trg.bool()

    rand_mat = torch.rand_like(trg_mat.float()) < random_ignore_mat
    rand_mat = torch.logical_or(rand_mat, filter_trg)
    batch['mask'] = torch.logical_and(batch['mask'], rand_mat)
    return batch


if __name__ == '__main__':

    print("CUDA AVAILABLE", torch.cuda.is_available())
    print("CUDA DEVICES", torch.cuda.device_count())

    parser = argparse.ArgumentParser(description='Finetune RNAformer')
    parser.add_argument('-c', '--config', type=str, help='config file name')

    args, unknown_args = parser.parse_known_args()
    config_name = args.config

    if not config_name.endswith('.yml'):
        config_name += '.yml'

    config_file = os.path.join("configs", config_name)
    with open(config_file, 'r') as f:
        config_dict = yaml.load(f, Loader=yaml.Loader)
    config_dict = read_unknown_args(unknown_args, config_dict)

    if isinstance(config_dict["seed"], str):
        seed_list = [int(i) for i in config_dict["seed"].split("-")]
    else:
        seed_list = [config_dict["seed"]]

    ft_config = Config(config_dict=config_dict)
    expt_dir = get_experiment_folder(**ft_config.experiment, new_folder=True, count_folder=False)
    ft_config.save_config(expt_dir)

    seq_vocab_size = 5
    model_dir = pathlib.Path(ft_config.model_dir)

    if 'config.yml' in os.listdir(model_dir):
        config = Config(config_file=model_dir / 'config.yml')
    else:
        config_name = ft_config.model_ckpt.split(".")[0] + ".yml"
        config_name = config_name.replace("state_dict", "config")
        config = Config(config_file=model_dir / config_name)

    checkpoint = torch.load(model_dir / ft_config.model_ckpt)
    if 'state_dict' in checkpoint:
        state_dict = checkpoint['state_dict']
    else:
        state_dict = checkpoint
    for key in list(state_dict.keys()):
        if key.startswith('model.'):
            state_dict[key[6:]] = state_dict.pop(key)

    config.RNAformer.max_len = ft_config.maxLength
    config.RNAformer.seq_vocab_size = seq_vocab_size
    flash_attn = config.RNAformer.flash_attn

    model = RNAFormer(config.RNAformer)
    model.load_state_dict(state_dict, strict=True)

    trainable_p = trainable_params(model)
    total_p = total_params(model)
    print(
        f"Total params: {total_p:,d} || Trainable params: {trainable_p:,d} || Fraction %: {100 * trainable_p / total_p:.2f}")

    if ft_config.cycling:
        model.initialize_cycling(ft_config.cycling, recycle_norm=False)

    print("Load Datasat")

    pad_index = 0
    ignore_index = -100
    collator = CollatorRNA(pad_index, ignore_index)

    file_train = os.path.join(ft_config.data_dir, ft_config.data_file)
    df = pd.read_pickle(file_train)
    df_train = df[df['set'] == 'train']
    df_valid = df[df['set'] == 'valid']
    raw_samples_train = df_train.to_dict(orient="records")
    raw_samples_valid = df_valid.to_dict(orient="records")

    train_samples = []
    for sample in tqdm(raw_samples_train):
        if len(sample["sequence"]) > ft_config.maxLength:
            continue
        else:
            if any(["T" == s for s in sample["sequence"]]):
                sample["sequence"] = ["U" if s == "T" else s for s in sample["sequence"]]
            sample = prepare_RNA_sample(sample)
            train_samples.append(sample)

    train_loader = torch.utils.data.DataLoader(train_samples, batch_size=ft_config.effective_batch_size,
                                               shuffle=True, num_workers=0, collate_fn=collator, drop_last=True)

    valid_samples = []
    for sample in tqdm(raw_samples_valid):
        if len(sample["sequence"]) > ft_config.maxLength:
            continue
        else:
            if any(["T" == s for s in sample["sequence"]]):
                sample["sequence"] = ["U" if s == "T" else s for s in sample["sequence"]]
            sample = prepare_RNA_sample(sample)
            valid_samples.append(sample)

    valid_loader = torch.utils.data.DataLoader(valid_samples, batch_size=ft_config.effective_batch_size,
                                               shuffle=False, num_workers=0, collate_fn=collator, drop_last=False)

    print("Train Samples", len(train_samples))

    L.seed_everything(ft_config.seed)

    strategy = L.fabric.strategies.DDPStrategy(
        find_unused_parameters=False,
        static_graph=False,
        process_group_backend="nccl"
    )

    tb_logger = TensorBoardLogger(expt_dir, name="ft")

    fabric = L.Fabric(devices=ft_config.devices, precision=ft_config.precision, strategy=strategy,
                      loggers=tb_logger)

    fabric.launch()

    if ft_config.use_cpr:
        parameters = group_parameters_for_cpr_optimizer(model)
        optimizer = AdamCPR(parameters, lr=ft_config.learning_rate, kappa_init_param=ft_config.cpr_param,
                            kappa_init_method=ft_config.cpr_init)
    else:
        optimizer = torch.optim.AdamW(model.parameters(), lr=ft_config.learning_rate,
                                      weight_decay=ft_config.weight_decay)

    model = fabric.setup_module(model)
    optimizer = fabric.setup_optimizers(optimizer)

    model = model.train()

    test_sets = pd.read_pickle(ft_config.test_set_file)
    print("train_samples", len(train_samples))

    train_loader = fabric.setup_dataloaders(train_loader)
    valid_loader = fabric.setup_dataloaders(valid_loader)


    def get_cosine_annealing_with_warmup(learning_rate, num_training_steps, num_warmup_steps, decay_factor):
        def lr_lambda(current_step: int):
            training_steps = num_training_steps - num_warmup_steps
            if current_step <= num_warmup_steps:
                return learning_rate * float(current_step) / float(max(1, num_warmup_steps))
            else:
                # return learning_rate
                cosine_decay = max(0.0, (1 + math.cos(
                    math.pi * (current_step - num_warmup_steps) / float(max(1, training_steps)))) / 2)
                return learning_rate * (decay_factor + (1 - decay_factor) * cosine_decay)

        return lambda current_step: lr_lambda(current_step)


    def get_constant_with_warmup(learning_rate, num_training_steps, num_warmup_steps):
        def lr_lambda(current_step: int):
            if current_step <= num_warmup_steps:
                return learning_rate * float(current_step) / float(max(1, num_warmup_steps))
            else:
                return learning_rate

        return lambda current_step: lr_lambda(current_step)


    def get_up_with_warmup(learning_rate, num_training_steps, num_warmup_steps, decay_factor):
        def lr_lambda(current_step: int):
            training_steps = num_training_steps - num_warmup_steps
            if current_step <= num_warmup_steps:
                return learning_rate * float(current_step) / float(max(1, num_warmup_steps))
            else:
                increase = float(current_step - num_warmup_steps) / float(training_steps - num_warmup_steps)
                return learning_rate * (1 + decay_factor * increase)

        return lambda current_step: lr_lambda(current_step)


    gradient_accumulation_iters = ft_config.batch_size // ft_config.effective_batch_size // ft_config.devices

    num_optim_steps = ft_config.max_train_steps
    num_warmup_steps = ft_config.warmup_steps
    decay_factor = 0.01

    if ft_config.lr_scheduler == "cosine":
        lr_scheduler = get_cosine_annealing_with_warmup(ft_config.learning_rate, num_optim_steps, num_warmup_steps,
                                                        decay_factor)
    elif ft_config.lr_scheduler == "constant":
        lr_scheduler = get_constant_with_warmup(ft_config.learning_rate, num_optim_steps, num_warmup_steps)
    elif ft_config.lr_scheduler == "up":
        lr_scheduler = get_up_with_warmup(ft_config.learning_rate, num_optim_steps, num_warmup_steps, decay_factor)

    optim_step_count = 0
    step_count = 0
    epoch = 0

    loss_train = torch.nn.BCEWithLogitsLoss(reduction='none')

    model.train()

    for epoch in range(ft_config.max_epochs):

        valid_loss = 0
        for batch_idx, batch in tqdm(enumerate(valid_loader), total=len(valid_loader), disable=fabric.world_size > 1):
            pdb_sample = torch.ones_like(batch['pdb_sample'])

            logits_mat, mask = model(batch['src_seq'], batch['length'], pdb_sample, n_cycles=0)

            pred = logits_mat[batch['mask']][:, 0]
            target = batch['trg_mat'][batch['mask']].float()
            loss = loss_train(pred, target)
            loss = loss[~torch.isnan(loss)]
            loss = torch.mean(loss)
            valid_loss += loss.item()

            fabric.log("valid/loss", loss.item(), optim_step_count)
        if optim_step_count % 10 == 0:
            print(f"validation epoch: {epoch} - loss {valid_loss / batch_idx:.4f}")

        random.shuffle(train_samples)
        train_loader = torch.utils.data.DataLoader(train_samples, batch_size=ft_config.effective_batch_size,
                                                   shuffle=True, num_workers=0, collate_fn=collator, drop_last=True)
        train_loader = fabric.setup_dataloaders(train_loader)
        for batch_idx, batch in tqdm(enumerate(train_loader), total=len(train_loader), disable=fabric.world_size > 1):

            if ft_config.random_crop_mat:
                batch = crop_mat(batch, ft_config.random_crop_mat)

            if ft_config.random_ignore_mat:
                batch = ignore_partial_mat(batch, ft_config.random_ignore_mat)

            is_accumulating = batch_idx % gradient_accumulation_iters != 0

            with fabric.no_backward_sync(model, enabled=is_accumulating):
                if model.cycling and ft_config.cycling:
                    if ft_config.cycle_warm_up > 0:
                        max_cycle = min(max(1, epoch // ft_config.cycle_warm_up), model.cycle_steps)
                    else:
                        max_cycle = model.cycle_steps

                    n_cycles = torch.randint(0, max_cycle + 1, [1])
                    n_cycles = fabric.all_gather(n_cycles)
                    n_cycles = n_cycles[0].item()
                    fabric.log(f"train/n_cycles", n_cycles, optim_step_count)
                else:
                    n_cycles = 0

                pdb_sample = torch.ones_like(batch['pdb_sample'])

                logits_mat, mask = model(batch['src_seq'], batch['length'], pdb_sample, n_cycles=n_cycles)

                fabric.barrier()

                pred = logits_mat[batch['mask']][:, 0]
                target = batch['trg_mat'][batch['mask']].float()
                loss = loss_train(pred, target)
                loss = loss[~torch.isnan(loss)]
                loss = torch.mean(loss)

                fabric.backward(loss)

            if not is_accumulating and step_count != 0:
                for param_group in optimizer.param_groups:
                    param_group["lr"] = lr_scheduler(optim_step_count)
                fabric.clip_gradients(model, optimizer, clip_val=ft_config.grad_clip)
                optimizer.step()
                optimizer.zero_grad()

                fabric.log("epoch/step_count", epoch, step_count)
                fabric.log("epoch/optim_step_count", epoch, optim_step_count)
                fabric.log("train/loss", loss.item(), step_count)
                for param_group in optimizer.param_groups:
                    fabric.log("optim/step_count", optim_step_count, step_count)
                    fabric.log("optim/epoch", epoch, optim_step_count)
                    fabric.log("optim/learning_rate", param_group["lr"], optim_step_count)

                optim_step_count += 1

            if optim_step_count % 10 == 0:
                print(
                    f"epoch: {epoch} - iteration: {batch_idx} - optim step: {optim_step_count} - loss {loss.item():.4f}")
            step_count += 1

            if optim_step_count >= ft_config.max_train_steps:
                break
        if optim_step_count >= ft_config.max_train_steps:
            break

        epoch += 1
        fabric.barrier()

    print(f"performance after {epoch + 1} epoch of fine tuning")
    test_results, all_samples = evaluate_RNAformer(model, test_sets, eval_synthetic=False,
                                                   eval_bprna=False)

    for _, row in test_results.iterrows():
        fabric.log(f"eval/{row['test_set']}", row['f1_score'], optim_step_count)
    print(test_results.to_markdown())

    if ft_config.save_samples and fabric.local_rank == 0:
        torch.save(all_samples, os.path.join(model_dir, f"finetune_test_samples.pt"))
        test_results.to_csv(os.path.join(expt_dir, f"finetune_results.csv"))

        state = {"model": model}
        fabric.save(os.path.join(expt_dir, f"finetune_model.ckpt"), state)

    sys.exit()
