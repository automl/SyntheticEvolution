import os
import statistics
import argparse

def read_list_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            # This assumes one number per line or comma-separated numbers
            return [float(num.strip()) for num in file.read().replace('[', '').replace(']', '').split(',')]
    except Exception as e:
        print(f"Failed to read {file_path}: {str(e)}")
        exit(1)

def calc_mean_std_and_write(list1, list2, filename, directory):
    if len(list1) != len(list2):
        return "Error: The lists do not match in size."

    if not os.path.exists(directory):
        os.makedirs(directory)
    filepath = os.path.join(directory, filename)
    results = []

    for i in range(len(list1)):
        mean = statistics.mean([list1[i], list2[i]])
        std = statistics.stdev([list1[i], list2[i]]) if len({list1[i], list2[i]}) > 1 else 0
        results.append(f"{mean:.7f} +/- {std:.7f}")

    with open(filepath, 'w') as file:
        for result in results:
            file.write(result + "\n")
    return f"Results successfully written to {filepath}"

def main():
    parser = argparse.ArgumentParser(description="Calculate means and stds of list pairs and write to file.")
    parser.add_argument("--filename", type=str, default="bio_RNAformer_layers_0_attn_pair_col_input_norm_weight.txt", help="The filename where the results will be saved.")
    parser.add_argument("--directory", type=str, default="stats", help="The directory where the file will be saved.")
    args = parser.parse_args()

    list1_path = os.path.join(args.directory, 'list1.txt')
    list2_path = os.path.join(args.directory, 'list2.txt')

    list1 = read_list_from_file(list1_path)
    list2 = read_list_from_file(list2_path)
    message = calc_mean_std_and_write(list1, list2, args.filename, args.directory)
    print(message)

if __name__ == "__main__":
    main()
