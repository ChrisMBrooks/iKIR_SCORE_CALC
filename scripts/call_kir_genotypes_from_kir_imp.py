import argparse
import os
import numpy as np
import pandas as pd

def parse_arguments() -> dict:
    parser = argparse.ArgumentParser(
        description="Script to call KIR genotypes from KIR*IMP results (e.g. imputations.csv).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-if",
        "--input_file",
        help="Filename for the KIR*IMP results file (imputations.csv). Full path filename as .csv.",
        required=True,
        type=str,
    )

    required.add_argument(
        "-od",
        "--output_dir",
        help="Output directory, full path filename as str.",
        required=True,
        type=str,
    )

    optional.add_argument(
        "-of",
        "--output_filename",
        help="Output filename as .csv.",
        required=False,
        type=str,
        default="kir_genotype_calls.csv"
    )

    optional.add_argument(
        "-qt",
        "--quality_threshold",
        help="Quality threshold used to call a variant, as float.",
        required=False,
        type=float,
        default=0.5
    )
    
    arguments = vars(parser.parse_args())
    return arguments

def import_tabular_file(filename:str, index_col:int = None) -> pd.DataFrame:
    file_extension = filename.split(".")[-1]
    if file_extension == "parquet":
        frame = pd.read_parquet(filename)
        return frame
    elif file_extension == "csv":
        frame = pd.read_csv(filename, index_col=index_col)
        return frame
    else:
        raise Exception("File type not recognised. Please us parquet or csv.")

def call_kir3dl1(kir3dl1ex4:bool, kir3dl1ex9:bool) -> bool:
    if kir3dl1ex4 or kir3dl1ex9:
        return True
    elif kir3dl1ex4 is None or kir3dl1ex9 is None:
        return None
    else:
        return False

def call_haplotype_gene(haplotype:pd.Series, quality_threshold:float) -> bool:
    if float(haplotype["posterior_probability"]) > quality_threshold:
        if int(haplotype["imputed_type"]) > 0:
            return True
        else:
            return False
    return None

def call_kir_genotype(raw_kir_imp_results:pd.DataFrame, required_kir_genes:list, quality_threshold:float) -> pd.DataFrame:
    records = []
    for subject_id in np.unique(raw_kir_imp_results["subject_id"]):
        kir_genotype_calls = {}
        subject_subset = raw_kir_imp_results[raw_kir_imp_results["subject_id"] == subject_id]
        for kir_gene in required_kir_genes:
            new_kir_key = kir_gene.lower()
            hap_key_1 = "{0}_hap_{1}".format(new_kir_key, 1)
            hap_key_2 = "{0}_hap_{1}".format(new_kir_key, 2)

            kir_subset = subject_subset[subject_subset["locus"] == kir_gene]

            haplotype_1_status = call_haplotype_gene(kir_subset.iloc[0], quality_threshold)
            haplotype_2_status = call_haplotype_gene(kir_subset.iloc[1], quality_threshold)

            kir_genotype_calls[hap_key_1] = haplotype_1_status
            kir_genotype_calls[hap_key_2] = haplotype_2_status

            if haplotype_1_status or haplotype_2_status:
                kir_genotype_calls[new_kir_key] = True
            elif haplotype_1_status is None or haplotype_2_status is None:
                kir_genotype_calls[new_kir_key] = None
            else:
                kir_genotype_calls[new_kir_key] = False

        # Final tasks
        kir_genotype_calls['kir3dl1'] = call_kir3dl1(
            kir_genotype_calls['kir3dl1ex4'], 
            kir_genotype_calls['kir3dl1ex9']
        )
        kir_genotype_calls["subject_id"] = subject_id
        kir_genotype_calls["applied_threshold"] = quality_threshold
        records.append(kir_genotype_calls)

    kir_genotype_calls = pd.DataFrame(records)

    return kir_genotype_calls

def compute_frequency_stats(kir_genotype_calls:pd.DataFrame) -> pd.DataFrame:
    desired_ikirs = ["kir2dl1", "kir2dl2", "kir2dl3", "kir3dl1"]
    frequencies = {key:0.0 for key in desired_ikirs}
    kir_genotype_calls = kir_genotype_calls.dropna()

    for ikir in desired_ikirs:
        frequency = (kir_genotype_calls[ikir] == True).sum()/kir_genotype_calls.shape[0]
        frequencies[ikir] = frequency

    frequencies = pd.DataFrame([frequencies])
    return frequencies

def format_and_export_results(kir_genotype_calls:pd.DataFrame, kir_frequencies:pd.DataFrame, output_filename:str, output_dir:str):
    # Format and Export Primary Results
    output_filename = os.path.join(output_dir, output_filename)
    columns = list(kir_genotype_calls.columns)
    columns.remove("subject_id")
    columns.remove("applied_threshold")
    columns = sorted(columns)
    columns = ["subject_id", "applied_threshold"] + sorted(columns)
    kir_genotype_calls = kir_genotype_calls[columns]
    kir_genotype_calls.to_csv(output_filename)
    print(kir_genotype_calls)

    # Export Frequencies
    output_filename = "kir_genotype_frequencies.csv"
    output_filename = os.path.join(output_dir, output_filename)
    kir_frequencies.to_csv(output_filename)
    print("KIR Population Frequencies:")
    print(kir_frequencies)

def main():
    args = parse_arguments()
    os.makedirs(args["output_dir"], exist_ok=True)

    # Declare Constants 
    required_kir_genes = [
        "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR3DL1ex4", "KIR3DL1ex9", 
        "KIR2DS2", "KIR2DP1", "KIR3DP1", "KIR2DL4", "KIR3DS1", 
        "KIR2DL5", "KIR2DS3", "KIR2DS5", "KIR2DS1"
    ]

    # Import Raw KIR*IMP Data
    raw_kir_imp_results = import_tabular_file(args["input_file"])
    columns = ['subject_id', 'subject_id2', 'haplotype_id', 'locus', 'imputed_type', 'posterior_probability']
    raw_kir_imp_results = pd.DataFrame(raw_kir_imp_results.values, columns=columns)

    kir_genotype_calls = call_kir_genotype(
        raw_kir_imp_results=raw_kir_imp_results, 
        required_kir_genes = required_kir_genes,
        quality_threshold=args["quality_threshold"]
    )

    # Compute QC  Stats
    kir_frequencies = compute_frequency_stats(kir_genotype_calls)

    # Format and Export Results
    format_and_export_results(
        kir_genotype_calls = kir_genotype_calls, 
        kir_frequencies = kir_frequencies, 
        output_filename = args["output_filename"], 
        output_dir = args["output_dir"]
    )

print("Starting...")
main()
print("Complete.")