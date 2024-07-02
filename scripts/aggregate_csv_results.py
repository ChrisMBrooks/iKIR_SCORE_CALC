import argparse
import os
import pandas as pd

def parse_arguments() -> dict:
    parser = argparse.ArgumentParser(
        description="Script to retrieve call calculate functional iKIR score from HLA and KIR genotypes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-id",
        "--input_dir",
        help="Directory of the path containing the files to be combined as str.",
        required=True,
        type=str,
    )

    optional.add_argument(
        "-dd",
        "--drop_duplicates",
        help="Flag to indicate whether to drop duplicates.",
        required=False,
        type=bool,
        action="store_true"
    )

    required.add_argument(
        "-od",
        "--output_dir",
        help="Output directory, full path filename as str.",
        required=True,
        type=str,
    )

    required.add_argument(
        "-of",
        "--output_filename",
        help="Output filename as .csv.",
        required=True,
        type=str,
        default="functional_ikir_scoring.csv"
    )
    
    arguments = vars(parser.parse_args())
    return arguments

def main():
    args = parse_arguments()
    os.makedirs(args["output_dir"], exist_ok=True)

    frames = []
    for filename in os.listdir(args["input_dir"]):
        full_path_filename = os.path.join(args["input_dir"], filename)
        if os.path.isfile(full_path_filename) and full_path_filename.endswith(".csv"):
            frame = pd.read_csv(full_path_filename, index_col=0)
            frames.append(frame)

    output_frame = pd.concat(frames, ignore_index=True)
    if args["drop_duplicates"]:
        output_frame = output_frame.drop_duplicates()

    output_filename = os.path.join(args["output_dir"], args["output_filename"])
    output_frame.to_csv(output_filename)
    print(output_frame)
            
print("Starting...")
main()
print("Complete.")