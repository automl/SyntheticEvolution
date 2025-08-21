import torch
import argparse

import pandas as pd

def ckpt_to_pth(ckpt_path):
    checkpoint = torch.load(ckpt_path, map_location=torch.device('cpu'))
    
    if 'state_dict' in checkpoint:
        state_dict = checkpoint['state_dict']
    elif 'model_state' in checkpoint:
        state_dict = checkpoint['model_state']
    else:
        state_dict = checkpoint

    # breakpoint()
    pth_path = ckpt_path.split(".")[0]+".pth"
    torch.save(state_dict, pth_path)
    print(f"Model's state dictionary has been saved to {pth_path}")

def pth_to_txt(pth_path):
    """
    Load a .pth file and write its contents in a human-readable format to a text file.

    Args:
    - pth_path (str): Path to the .pth file.
    """
    # Load the state dictionary from the .pth file
    state_dict = torch.load(pth_path, map_location=torch.device('cpu'))

    txt_path = pth_path.split(".")[0]+"_256.txt"
    with open(txt_path, 'w') as file:
        # Write each parameter in the state dictionary to the file
        for param_name, tensor in state_dict.items():
            if "model" in param_name:
                temp_name=param_name.split(".")[1:]
                param_name=(".").join(temp_name)
            file.write(f"{param_name}:\n")
            
            if tensor.numel() < 257:  # Check if the number of elements in the tensor is manageable
                file.write(f"{tensor.tolist()}\n")
            else:
                file.write(f"Shape: {tensor.shape}\n")  # Otherwise, just write the shape of the tensor
            
            file.write("\n")  # Add a newline for better separation between parameters

def process_data(df, filename):
    # Initialize new columns
    df['pred_pairs'] = None
    df['true_pairs'] = None
    
    # Iterate over rows in DataFrame
    for index, row in df.iterrows():
        # Extracting pred_pairs
        true_positions = list(zip(*row['pred_mat'].nonzero()))
        df.at[index, 'pred_pairs'] = true_positions
        
        # Extracting true_pairs from pos1id and pos2id
        pairs = [(int(pos1.item()), int(pos2.item())) for pos1, pos2 in zip(row['pos1id'], row['pos2id'])]
        df.at[index, 'true_pairs'] = pairs
    
    save_pairs_to_file(df, filename)
    return df

def save_pairs_to_file(df, filename):
    with open(filename, 'w') as file:
        for index, row in df.iterrows():
            # Get maximum length of pairs to align the outputs
            max_length = max(len(row['pred_pairs']), len(row['true_pairs']))
            
            file.write(f"Sample {index + 1} Pairs:\n")
            file.write(f"{'Pred Pairs':30s} | True Pairs\n")
            file.write("-" * 60 + "\n")
            
            for i in range(max_length):
                # Safely access the index of each list or provide an empty output
                pred_pair = row['pred_pairs'][i] if i < len(row['pred_pairs']) else ""
                true_pair = row['true_pairs'][i] if i < len(row['true_pairs']) else ""
                file.write(f"{str(pred_pair):30s} | {str(true_pair)}\n")
            
            file.write("\n")

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--mode", type=str)
    parser.add_argument("--path", type=str)
    parser.add_argument("--df_path", type=str)
    args=parser.parse_args()

    if args.mode == "pth":
        ckpt_to_pth(args.path)
    elif args.mode == "txt":
        pth_to_txt(args.path)
    elif args.mode == "check":
        df = pd.read_pickle(args.df_path)
        process_data(df, args.path)
    else:
        raise ValueError("Available modes are 'pth' and 'txt'.")