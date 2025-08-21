import argparse
import logging
import os
import pathlib
from collections import defaultdict

import torch
import torch.cuda
import pandas as pd

from RNAformer.model.RNAformer import RNAFormer
from RNAformer.utils.configuration import Config

logger = logging.getLogger(__name__)


def sequence2index_vector(sequence, mapping):
    int_sequence = list(map(mapping.get, sequence))
    if None in int_sequence:
        raise ValueError("Invalid sequence", sequence, int_sequence)
    return torch.LongTensor(int_sequence)


def prepare_RNA_sample(input_sample, inference=False):
    seq_vocab = ['A', 'C', 'G', 'U', 'N']
    seq_stoi = dict(zip(seq_vocab, range(len(seq_vocab))))

    seq_stoi['Y'] = seq_stoi['C']
    seq_stoi['R'] = seq_stoi['A']
    seq_stoi['W'] = seq_stoi['A']
    seq_stoi['S'] = seq_stoi['C']
    seq_stoi['K'] = seq_stoi['G']
    seq_stoi['M'] = seq_stoi['A']
    seq_stoi['B'] = seq_stoi['C']
    seq_stoi['D'] = seq_stoi['A']
    seq_stoi['H'] = seq_stoi['A']
    seq_stoi['V'] = seq_stoi['A']
    seq_stoi['T'] = seq_stoi['U']

    sequence = input_sample["sequence"]
    length = len(sequence)

    src_seq = sequence2index_vector(sequence, seq_stoi)
    torch_sample = {}
    torch_sample['src_seq'] = src_seq.clone()
    torch_sample['length'] = torch.LongTensor([length])[0]

    if not inference:
        pos1id = input_sample["pos1id"]
        pos2id = input_sample["pos2id"]
        pdb_sample = int(input_sample['is_pdb'])

        torch_sample['pos1id'] = torch.LongTensor(pos1id)
        torch_sample['pos2id'] = torch.LongTensor(pos2id)
        torch_sample['pdb_sample'] = torch.LongTensor([pdb_sample])[0]
        
        trg_mat = torch.LongTensor(length, length).fill_(0)
        trg_mat[pos1id, pos2id] = 1
        trg_mat[pos2id, pos1id] = 1
        torch_sample['trg_mat'] = trg_mat
        torch_sample['trg_seq'] = src_seq.clone()

    if inference:
        idx = input_sample["Id"]
        # breakpoint()
        if type(idx) is not int:
            torch_sample['Id'] = torch.tensor([ord(c) for c in idx], dtype=torch.long)
        else:
            torch_sample['Id'] = torch.LongTensor([idx])[0]

    return torch_sample


def evaluate_RNAformer(model, test_sets, eval_synthetic=False, eval_bprna=False, n_cycles=None):
    model = model.cuda()

    # check GPU can do bf16
    if torch.cuda.is_bf16_supported():
        model = model.bfloat16()
    else:
        model = model.half()

    model.eval()

    test_results = []
    all_samples = {}

    with torch.no_grad():
        for test_set, test_df in test_sets.items():

            if not eval_synthetic and test_set == "synthetic_test":
                continue
            if not eval_bprna and test_set == "bprna_ts0":
                continue

            print("Evaluating on", test_set, test_df.shape)
            metrics = defaultdict(list)

            inference_samples = []

            for id, sample_raw in test_df.iterrows():

                sample = prepare_RNA_sample(sample_raw)

                sequence = sample['src_seq'].unsqueeze(0).cuda()
                src_len = torch.LongTensor([sequence.shape[-1]]).cuda()
                if torch.cuda.is_bf16_supported():
                    pdb_sample = torch.FloatTensor([[1]]).bfloat16().cuda()
                else:
                    pdb_sample = torch.FloatTensor([[1]]).half().cuda()

                if n_cycles is not None:
                    n_cycles = n_cycles

                elif model.cycling:
                    n_cycles = model.cycle_steps
                else:
                    n_cycles = 0

                if n_cycles > 0:
                    model.cycling = True
                logits, pair_mask = model(sequence, src_len, pdb_sample, n_cycles=n_cycles)
                
                pred_mat = torch.sigmoid(logits[0, :, :, -1]) > 0.5
                true_mat = sample['trg_mat'].float().cuda()

                save_sample = {}
                for nkey in ['sequence', 'pk', 'has_multiplet', 'has_pk', 'set', 'has_nc']:
                    save_sample[nkey] = sample_raw[nkey]
                save_sample['true_mat'] = sample['trg_mat']
                save_sample['pred_mat'] = pred_mat.detach().cpu()
                save_sample['logits_mat'] = torch.sigmoid(logits[0, :, :, 0]).detach().cpu()
                inference_samples.append(save_sample)

                solved = torch.equal(true_mat, pred_mat).__int__()

                metrics['solved'].append(torch.tensor([solved], dtype=true_mat.dtype, device=true_mat.device))

                tp = torch.logical_and(pred_mat, true_mat).sum()
                tn = torch.logical_and(torch.logical_not(pred_mat), torch.logical_not(true_mat)).sum()
                fp = torch.logical_and(pred_mat, torch.logical_not(true_mat)).sum()
                fn = torch.logical_and(torch.logical_not(pred_mat), true_mat).sum()
                assert pred_mat.size().numel() == tp + tn + fp + fn
                accuracy = tp / pred_mat.size().numel()
                precision = tp / (1e-4 + tp + fp)
                recall = tp / (1e-4 + tp + fn)
                f1_score = 2 * tp / (1e-4 + (2 * tp + fp + fn))
                mcc = (tp * tn - fp * fn) / (1e-4 + torch.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
                metrics['accuracy'].append(accuracy)
                metrics['precision'].append(precision)
                metrics['recall'].append(recall)
                metrics['f1_score'].append(f1_score)
                metrics['mcc'].append(mcc)

            num_samples = len(metrics['mcc'])
            for key, value_list in metrics.items():
                metrics[key] = torch.stack(value_list).mean().item()
            metrics['num_samples'] = num_samples

            metrics['test_set'] = test_set
            test_results.append(metrics)
            all_samples[test_set] = inference_samples

    test_results = pd.DataFrame(test_results)

    return test_results, all_samples


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Evaluate RNAformer')
    parser.add_argument('-n', '--model_name', type=str, default='state_dict.ckpt')
    parser.add_argument('-m', '--model_dir', type=str, default='.')
    parser.add_argument('-s', '--save_predictions', action='store_true')
    parser.add_argument('-y', '--eval_synthetic', action='store_true')
    parser.add_argument('-b', '--eval_eval_bprna', action='store_true')
    parser.add_argument('-c', '--cycling', type=int, default=0)

    args, unknown_args = parser.parse_known_args()
    model_dir = args.model_dir

    seq_vocab_size = 5
    model_dir = pathlib.Path(model_dir)

    config = Config(config_file=model_dir / 'config.yml')

    if not hasattr(config.RNAformer, "flash_attn"):
        config_dict = config.get_dict()
        print("Updating config", config_dict)

        config_dict['RNAformer']['flash_attn'] = False
        config_dict["RNAformer"]["pos_embedding"] = True
        config_dict["RNAformer"]["learn_ln"] = True
        config_dict["RNAformer"]["rotary_emb"] = False
        config_dict["RNAformer"]["pdb_flag"] = False
        config_dict["RNAformer"]["binary_output"] = False
        config_dict["RNAformer"]["binary_output"] = False

        config = Config(config_dict=config_dict)

    if not hasattr(config.RNAformer, "pos_embedding"):
        config_dict = config.get_dict()
        config_dict["RNAformer"]["pos_embedding"] = True
        config_dict["RNAformer"]["learn_ln"] = True
        config_dict["RNAformer"]["rotary_emb"] = False
        config_dict["RNAformer"]["pdb_flag"] = False
        config_dict["RNAformer"]["binary_output"] = False
        config_dict["RNAformer"]["binary_output"] = False
        config = Config(config_dict=config_dict)

    config.RNAformer.seq_vocab_size = seq_vocab_size
    config.RNAformer.max_len = 600
    flash_attn = config.RNAformer.flash_attn

    checkpoint = torch.load(model_dir / args.model_name)
    if 'state_dict' in checkpoint:
        state_dict = checkpoint['state_dict']
    elif 'model' in checkpoint:
        state_dict = checkpoint['model']
    else:
        state_dict = checkpoint

    for key in list(state_dict.keys()):
        if key.startswith('model.'):
            state_dict[key[6:]] = state_dict.pop(key)

    if 'seq2mat_embed.src_embed_1.embed_pair_pos.weight' in state_dict:
        config.RNAformer.max_len = state_dict['seq2mat_embed.src_embed_1.embed_pair_pos.weight'].shape[1]
    else:
        config.RNAformer.max_len = state_dict['seq2mat_embed.src_embed_1.weight'].shape[1]

    if args.cycling:
        config.RNAformer.cycling = args.cycling
        n_cycles = args.cycling

    if "recycle_pair_norm.weight" in state_dict:
        recycle_norm = True
    else:
        recycle_norm = False

    model = RNAFormer(config.RNAformer, recycle_norm=recycle_norm)
    model.load_state_dict(state_dict, strict=True)

    if args.cycling and args.cycling > 0:
        model.cycle_steps = args.cycling


    def count_parameters(parameters):
        return sum(p.numel() for p in parameters)


    print(f"Model size: {count_parameters(model.parameters()):,d}")

    test_sets = pd.read_pickle("datasets/test_sets.pkl")
    test_results, all_samples = evaluate_RNAformer(model, test_sets, eval_synthetic=args.eval_synthetic,
                                                   eval_bprna=args.eval_eval_bprna, n_cycles=n_cycles)

    if args.save_predictions:
        torch.save(all_samples, os.path.join(model_dir, f"{args.model_name}_test_samples.pt"))

    print(test_results.to_markdown())
