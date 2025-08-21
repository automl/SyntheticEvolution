from typing import List
import os, sys, socket
import argparse, collections, yaml

import logging
import torch.cuda
import pytorch_lightning as pl
import numpy as np

from RNAformer.pl_module.datamodule_rna import DataModuleRNA
from RNAformer.pl_module.rna_folding_trainer import RNAFoldingTrainer

from RNAformer.utils.configuration import Config
from RNAformer.utils.instantiate import instantiate
from RNAformer.utils.folder_manager import get_experiment_folder


def bold(msg):
    return f"\033[1m{msg}\033[0m"


def main(cfg):
    """
    Launch pretraining
    """

    torch.set_float32_matmul_precision('high')

    if os.environ.get("GLOBAL_RANK") is None or os.environ.get("GLOBAL_RANK") == 0:
        is_rank_zero = True
        rank = 0
    else:
        is_rank_zero = False
        rank = os.environ.get("GLOBAL_RANK")


    if cfg.resume_training:
        exp_folder = get_experiment_folder(**cfg.experiment, new_folder=False, count_folder=False)
        if (exp_folder / "last.ckpt").exists():
            do_resume_training = True
        else:
            do_resume_training = False
    else:
        do_resume_training = False
        exp_folder = get_experiment_folder(**cfg.experiment, new_folder=is_rank_zero)

    logger = logging.getLogger(__name__)

    if is_rank_zero:
        cfg.save_config(exp_folder)

        logging.basicConfig(
            format="[%(asctime)s][%(levelname)s][%(name)s] - %(message)s",
            datefmt="%d/%m/%Y %H:%M:%S",
            level=logging.INFO,
            handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(exp_folder / "logfile.txt")],
        )

        logger.info(bold("######################################################"))
        logger.info(bold("########          START   TRAINING          ##########"))
        logger.info(bold("######################################################"))

        logger.info(f"########  Project:    {cfg.experiment.project_name}")
        logger.info(f"########  Session:    {cfg.experiment.session_name}")
        logger.info(f"########  Experiment: {cfg.experiment.experiment_name}")
        logger.info(f"save logs and checkpoints in: {exp_folder.as_posix()}")

        logger.info(bold("############### CONFIGURATION"))
        logger.info("RNA Task args")
        logger.info(cfg.rna_data)
        logger.info("Trainer args")
        logger.info(cfg.trainer)
        logger.info("Train args")
        logger.info(cfg.train)
        logger.info("Deepspeed args")
        logger.info(cfg.deepspeed)
        logger.info("Optimizer args")
        logger.info(cfg.train.optimizer)
        logger.info("RNAformer args")
        logger.info(cfg.RNAformer)

    # Set seed before initializing model
    np.random.seed(cfg.train.seed)
    torch.manual_seed(cfg.train.seed)
    torch.cuda.manual_seed_all(cfg.train.seed)

    logger.info(bold(f"############### LOAD DATA on rank {rank}"))

    rna_config = {**cfg.rna_data}
    data_module = DataModuleRNA(**rna_config, logger=logger)

    cfg.RNAformer.max_len = cfg.rna_data.max_len
    cfg.RNAformer.seq_vocab_size = data_module.seq_vocab_size

    model_module = RNAFoldingTrainer(
        cfg_train=cfg.train,
        cfg_model=cfg.RNAformer,
        py_logger=logger,
        data_module=data_module,
    )

    if is_rank_zero:
        def count_parameters(parameters):
            return sum(p.numel() for p in parameters if p.requires_grad)
        logger.info(f"#### trainable_parameters {count_parameters(model_module.parameters())}")


    if cfg.resume_training:
        logger.info(bold(f"############### RESUME TRAINING on rank {rank}"))

    logger.info(f'#### Load logger on rank {rank}')
    training_logger = pl.loggers.tensorboard.TensorBoardLogger(
        save_dir=exp_folder,
        name="",
        version="tb",
        prefix="",
    )

    logger.info(f"#### Load callbacks on rank {rank}")
    callbacks: List[pl.Callback] = []
    if "callbacks" in cfg:
        for cb_name, cb_conf in cfg.callbacks.items():
            if cb_conf is not None and "_target_" in cb_conf:
                logger.info(f"Instantiating callback <{cb_name}>")
                if "dirpath" in cb_conf:
                    cb_conf["dirpath"] = exp_folder
                callbacks.append(instantiate(cb_conf))

    if cfg.train.swa:
        callbacks.append(pl.callbacks.StochasticWeightAveraging(swa_epoch_start=120,
                                                                annealing_epochs=100,
                                                                swa_lrs=cfg.train.swa))

    logger.info(f'#### Load strategy on rank {rank}')
    strategy = pl.strategies.DDPStrategy(
        find_unused_parameters=False,   # TODO: check if this is needed
        static_graph=False,
        process_group_backend="nccl"
    )

    logger.info(bold(f"############### TRAINER on rank {rank}"))

    trainer = instantiate(cfg.trainer, instance=pl.Trainer,
                          callbacks=callbacks,
                          strategy=strategy,
                          logger=training_logger,
                          )


    logger.info(f"Starting training on rank {rank}")
    trainer.fit(
        model=model_module, datamodule=data_module,  ckpt_path=exp_folder / 'last.ckpt' if do_resume_training else None
    )

    if is_rank_zero:
        logger.info(f"Saving model to {exp_folder} on rank {rank}")
        trainer.save_checkpoint(exp_folder / "final_weights.ckpt", weights_only=True)
        logger.info(f"Finished saving model weights on rank {rank}")
    logger.info(f"Wait on barrier: rank {rank}")
    torch.distributed.barrier()

    logger.info("End training!")


if __name__ == "__main__":

    from functools import reduce  # forward compatibility for Python 3
    import operator

    def update(d, u):
        for k, v in u.items():
            if isinstance(v, collections.abc.Mapping):
                d[k] = update(d.get(k, {}), v)
            else:
                d[k] = v
        return d


    def getFromDict(dataDict, mapList):
        return reduce(operator.getitem, mapList, dataDict)

    def setInDict(dataDict, mapList, value):
        getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value

    def convert_string_value(value):
        if value in ('false', 'False'):
            value = False
        elif value in ('true', 'True'):
            value = True
        else:
            try:
                value = int(value)
            except:
                try:
                    value = float(value)
                except:
                    pass
        return value


    if socket.gethostname() == "tower":
        default_config_name = "default_config.yml"
    else:
        default_config_name = "juwels_config.yml"

    parser = argparse.ArgumentParser(description='Train RNAformer')
    parser.add_argument('-c', '--config', type=str, default=default_config_name, help='config file name')

    args, unknown_args = parser.parse_known_args()

    config_name = args.config
    if not config_name.endswith('.yml'):
        config_name += '.yml'

    config_file = os.path.join("configs", config_name)
    with open(config_file, 'r') as f:
        config_dict = yaml.load(f, Loader=yaml.Loader)

    for arg in unknown_args:
        if '=' in arg:
            keys = arg.split('=')[0].split('.')
            value = convert_string_value(arg.split('=')[1])
            print(keys, value)
            setInDict(config_dict, keys, value)
        else:
            raise UserWarning(f"argument unknown: {arg}")

    config = Config(config_dict=config_dict)

    main(cfg=config)
