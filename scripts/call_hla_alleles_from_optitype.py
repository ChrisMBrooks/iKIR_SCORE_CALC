import argparse
import os
import pandas as pd
import numpy as np

def parse_arguments() -> dict:
    parser = argparse.ArgumentParser(
        description="Script to retrieve call hla allele variants from technical replicates with quality score.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-if",
        "--input_file",
        help="Filename of the Optitype HLA allele typing results. Full path filename as .parquet.",
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
        default="hla_allele_calls.csv"
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

def find_best_allele(quality_scores:dict, quality_threshold:float) -> tuple:
    best_allele = None
    best_quality = 0.0
    for allele_name in quality_scores:
        current_quality_score = quality_scores[allele_name]
        if  current_quality_score > best_quality:
            best_quality = current_quality_score
            best_allele = allele_name
    if best_quality > quality_threshold:
        return best_allele, best_quality
    else:
        return None, None

def get_hla_allele_calls_from_tech_reps(hla_allele_subset:pd.DataFrame, quality_threshold:float) -> dict:
    hla_allele_subset = hla_allele_subset.copy()
    hla_allele_subset["hla_gene0"] = hla_allele_subset["hla_gene"].str.slice(0,1)
    hla_allele_subset = hla_allele_subset.sort_values(by="hla_allele")
    hla_allele_subset = hla_allele_subset.sort_values(by="hla_gene")

    hla_allele_quality = {}
    hla_allele_calls = {}

    for hla_gene in ["A1", "A2", "B1", "B2", "C1", "C2"]:
        subset = hla_allele_subset[hla_allele_subset["hla_gene"] == hla_gene].copy()
        uniques_alleles = np.unique(subset["hla_allele"].dropna())
        hla_allele_quality[hla_gene] = {key:0 for key in uniques_alleles}
        subset["quality_score"] =  subset["quality_score"]/ subset["quality_score"].sum()

        for index, row in subset.iterrows():
            hla_gene = row["hla_gene"]
            hla_allele = row["hla_allele"]
            if hla_allele:
                hla_allele_quality[hla_gene][hla_allele] += row["quality_score"]

        hla_allele_calls[hla_gene] = find_best_allele(hla_allele_quality[hla_gene], quality_threshold)[0]

    
    return hla_allele_calls

def call_hla_alleles(hla_alleles_by_sample:pd.DataFrame, quality_threshold:float) -> pd.DataFrame:

    subject_ids = np.unique(hla_alleles_by_sample["subject_id"])
    records = []
    for subject_id in subject_ids:
        technical_replicates = hla_alleles_by_sample[hla_alleles_by_sample["subject_id"] == subject_id]
        record = get_hla_allele_calls_from_tech_reps(technical_replicates, quality_threshold)
        record["subject_id"] = subject_id
        records.append(record)

    motif_posession_by_subject = pd.DataFrame(records)
    motif_posession_by_subject["quality_threshold"] = quality_threshold
    return motif_posession_by_subject

def format_and_export(hla_allele_calls:pd.DataFrame, output_filename:str, output_dir:str):
    hla_genes = ["A1", "A2", "B1", "B2", "C1", "C2"]

    for column_key in hla_genes:
        new_key = "{0}_{1}".format(column_key[0], column_key[1])
        hla_allele_calls[new_key] = hla_allele_calls[column_key]
    
    columns = list(hla_allele_calls.columns)
    for column_key in hla_genes:
        columns.remove(column_key)
    columns.remove("subject_id")
    columns.remove("quality_threshold")
    print(columns)

    columns = ["subject_id", "quality_threshold"] + sorted(columns)
    hla_allele_calls = hla_allele_calls[columns]
    columns = [key.lower() for key in columns]
    hla_allele_calls = pd.DataFrame(hla_allele_calls.values, columns=columns)

    output_filename = os.path.join(output_dir, output_filename)
    hla_allele_calls.to_csv(output_filename)

    print(hla_allele_calls)

def main():
    args = parse_arguments()
    os.makedirs(args["output_dir"], exist_ok=True)

    # Import Consolidated Results
    optitype_hla_results = pd.read_parquet(args["input_file"])

    # Compute HLA Variant Calls
    hla_allele_calls = call_hla_alleles(optitype_hla_results, args["quality_threshold"])

    # Format And Export
    format_and_export(hla_allele_calls, args["output_filename"], args["output_dir"])

print("Starting...")
main()
print("Complete.")