import argparse
import logging
import pathlib
import time

import torch
import torch.cuda
import pandas as pd
import numpy as np

from pathlib import Path


from RNAformer.model.RNAformer import RNAFormer
from RNAformer.utils.configuration import Config
from evaluate_RNAformer import prepare_RNA_sample

logger = logging.getLogger(__name__)

def filter_sequences(sequences_df):
    valid_bases = {'A', 'U', 'C', 'G'}
    erronous_elements = set()

    # Function to identify erronous elements in each sequence
    def check_sequence(sequence):
        erronous = set([item for item in sequence if item not in valid_bases])
        return all(item in valid_bases for item in sequence), erronous

    # Apply check on each sequence and filter out valid sequences
    results = sequences_df['sequence'].apply(check_sequence)
    mask = results.apply(lambda x: x[0])
    sequences_df.loc[:, 'erronous_elements'] = results.apply(lambda x: x[1])

    # Aggregating all erronous elements found across all sequences
    for err_set in sequences_df.loc[~mask, 'erronous_elements']:
        erronous_elements.update(err_set)

    valid_df = sequences_df[mask]
    erronous_df = sequences_df[~mask].drop(columns=['erronous_elements'])

    # Printing the erronous elements and their counts
    if erronous_elements:
        print("Erronous elements found in sequences:", erronous_elements)
        element_count = {element: 0 for element in erronous_elements}
        for element in sequences_df.loc[~mask, 'sequence'].explode():
            if element in erronous_elements:
                element_count[element] += 1
        print("Counts of erronous elements:", element_count)
    return valid_df, erronous_df

def infer_RNAformer(model, path, n_cycles=None):
    model = model.cuda()

    # model = model.cpu()

    # check GPU can do bf16
    if torch.cuda.is_bf16_supported():
        model = model.bfloat16()
    else:
        model = model.half()

    model.eval()
    # breakpoint()

    sequences_df = pd.read_pickle(path)
    # sequences_df = sequences_df["pdb_ts"]
    # breakpoint()
    # valid_df, erronous_df = filter_sequences(sequences_df)
    base_path = path.split(".")[:-1]
    pred_path = ".".join(base_path) + "_predictions.pkl"
    # breakpoint()
    # if not erronous_df.empty:
    #     err_path = ".".join(base_path) + "_erronous.pkl"
    #     erronous_df.to_pickle(err_path) 
    #     print(f"Erronous samples saved to {err_path}")

    start_time = time.time()
    with torch.no_grad():       
        print(f"Infering {sequences_df.shape[0]} RNA sequences")
        inference_samples = []

        for _, sample_raw in sequences_df.iterrows():

            sample = prepare_RNA_sample(sample_raw, inference=True)
            sequence = sample['src_seq'].unsqueeze(0).cuda()

            save_sample = {}
            save_sample['sequence'] = sequence.detach().cpu().numpy() 
            save_sample['Id'] = sample['Id']

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

            logits, _ = model(sequence, src_len, pdb_sample, n_cycles=n_cycles)
            
            pred_mat = torch.sigmoid(logits[0, :, :, -1]) > 0.5
            save_sample['pred_mat'] = pred_mat.detach().cpu().numpy() 
            inference_samples.append(save_sample)

    inference_df = pd.DataFrame(inference_samples)
    inference_df.to_pickle(pred_path)
    elapsed_time = time.time() - start_time
    print(f"Predictions saved to {pred_path}")
    print(f"Inference completed in {elapsed_time:.2f} seconds")
    # Auto-convert to compact format (Id, sequence, pairs)
    try:
        convert_predictions_file(Path(pred_path))
    except Exception as e:
        print(f"[post] Conversion failed: {e}")



def _decode_id(id_field):
    """
    Robustly turn various ID encodings into a human-readable string.

    Accepts:
      - str -> str
      - scalar (torch/numpy/python int) -> "int"
      - 1D array/list of ints/bytes -> ascii if possible, else join numbers
      - bytes/bytearray -> ascii (errors ignored)
    """
    # 1) Already a string?
    if isinstance(id_field, str):
        return id_field

    # 2) Bytes-ish?
    if isinstance(id_field, (bytes, bytearray, memoryview)):
        try:
            return bytes(id_field).decode("utf-8", errors="ignore")
        except Exception:
            try:
                return bytes(id_field).decode("ascii", errors="ignore")
            except Exception:
                return repr(id_field)

    # 3) Torch tensor?
    try:
        import torch  # type: ignore
        if isinstance(id_field, torch.Tensor):
            if id_field.ndim == 0:
                # scalar -> "int"
                return str(int(id_field.item()))
            # convert to numpy for unified handling
            id_field = id_field.detach().cpu().numpy()
    except Exception:
        pass

    # 4) Numpy scalar?
    try:
        import numpy as np  # type: ignore
        if isinstance(id_field, np.generic) and id_field.ndim == 0:  # scalar
            return str(int(id_field.item()))
    except Exception:
        pass

    # 5) Sequence of ints? Try ASCII decode first.
    try:
        # normalize to flat list of ints
        if hasattr(id_field, "tolist"):
            arr = id_field.tolist()
        else:
            arr = list(id_field)

        # flatten single nesting like [[...]] or (..,)
        while isinstance(arr, (list, tuple)) and len(arr) == 1 and isinstance(arr[0], (list, tuple)):
            arr = arr[0]

        # coerce elements to ints (handle tiny tensors, numpy scalars, etc.)
        nums = []
        for v in arr:
            try:
                vv = int(getattr(v, "item", lambda: v)())
            except Exception:
                continue
            nums.append(vv)

        if nums:
            # If they look like bytes, decode to text
            if all(0 <= n <= 255 for n in nums):
                try:
                    return bytes(nums).decode("utf-8", errors="ignore").strip("\x00")
                except Exception:
                    return bytes(nums).decode("ascii", errors="ignore").strip("\x00")
            # If they look like printable ASCII code points, decode via chr
            if all(32 <= n <= 126 for n in nums):
                try:
                    return "".join(chr(n) for n in nums)
                except Exception:
                    pass
            # Fallback: join the numbers
            return "_".join(str(n) for n in nums)
    except Exception:
        pass

    # 6) Last resort: string-repr
    return str(id_field)

def _decode_seq(seq_array):
    """
    Accepts saved numpy array like shape (1, L) or (L,), returns string of A/C/G/U/N.
    Unknown indices map to 'N'.
    """
    seq_vocab = ['A', 'C', 'G', 'U', 'N']
    seq_itos = dict(enumerate(seq_vocab))
    # squeeze first dimension if present
    seq_array = np.asarray(seq_array)
    if seq_array.ndim == 2 and seq_array.shape[0] == 1:
        seq_array = seq_array[0]
    out = []
    for idx in seq_array:
        try:
            ii = int(idx)
        except Exception:
            try:
                ii = int(idx.item())
            except Exception:
                ii = None
        out.append(seq_itos.get(ii, 'N'))
    return "".join(out)

def _mat2pairs(matrix, symmetric=True):
    """
    Convert predicted contact/paired matrix (bool/0-1) to list of index pairs (i, j).
    If symmetric=True, each undirected pair appears once (i <= j).
    """
    mat = np.asarray(matrix)
    nnz = np.argwhere(mat == 1) if mat.dtype != bool else np.argwhere(mat)
    if symmetric:
        return list({tuple(sorted(map(int, pair))) for pair in nnz})
    else:
        return [tuple(map(int, pair)) for pair in nnz]

def convert_predictions_file(predictions_path: Path) -> Path:
    """
    Reads <*_predictions.pkl>, writes <*_processed.pkl> with columns: Id, sequence, pairs.
    Returns output path.
    """
    predictions_path = Path(predictions_path)
    df = pd.read_pickle(predictions_path)

    # Robust decoding
    if 'Id' in df.columns:
        df.loc[:, 'Id'] = df['Id'].apply(_decode_id)
    if 'sequence' in df.columns:
        df.loc[:, 'sequence'] = df['sequence'].apply(_decode_seq)
    if 'pred_mat' in df.columns:
        df.loc[:, 'pairs'] = df['pred_mat'].apply(lambda x: _mat2pairs(x, symmetric=True))
    else:
        df.loc[:, 'pairs'] = [[] for _ in range(len(df))]

    outpath = predictions_path.with_name(predictions_path.stem.replace('_predictions', '') + '_processed.pkl')
    df[['Id', 'sequence', 'pairs']].to_pickle(outpath)

    # Pretty echo for quick CLI feedback
    print(f"[post] Processed dataframe written to: {outpath.resolve()}")
    try:
        print(df[['Id', 'sequence', 'pairs']].head(3))
    except Exception:
        pass
    


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Inference RNAformer')
    parser.add_argument('-n', '--model_name', type=str, default='state_dict.pth')
    parser.add_argument('-m', '--model_dir', type=str, default='.')
    parser.add_argument('-p', '--path', type=str)
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
    config.RNAformer.max_len = 500 # for data
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

    infer_RNAformer(model, args.path, n_cycles=n_cycles)