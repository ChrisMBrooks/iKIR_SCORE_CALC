import argparse
import os
import math
import pandas as pd
import numpy as np

def parse_arguments() -> dict:
    parser = argparse.ArgumentParser(
        description="Script to call KIR genotypes from T1K Results (consolidated) results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-if",
        "--input_file",
        help="Filename for the (consolidated) T1K results file. Full path filname as .parquet.",
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
        "-trqt",
        "--tech_rep_quality_threshold",
        help="Quality threshold used to call a variant, as float.",
        required=False,
        type=float,
        default=0.5
    )

    optional.add_argument(
        "-pvqt",
        "--p_val_quality_threshold",
        help="Quality threshold used to call a variant, as float.",
        required=False,
        type=float,
        default=0.1
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

def get_genotype_call_from_tech_reps(
    kir_gene_subset:pd.DataFrame, 
    tech_rep_quality_threshold:float
) -> bool:
        
        kir_gene_subset = kir_gene_subset.copy()
        replicate_count = kir_gene_subset.shape[0]

        # NB - Nans are ignored from both counts. 
        gene_present = (kir_gene_subset["allele_count"] > 0).sum() > (replicate_count*tech_rep_quality_threshold)
        gene_absent = (kir_gene_subset["allele_count"] == 0).sum() > (replicate_count*tech_rep_quality_threshold)

        if gene_present and gene_absent:
            return None
        elif gene_present:
            return True
        elif gene_absent:
            return False
        else:
            return None
        
def call_kir_genotypes_from_gex(
    raw_kir_from_gex:pd.DataFrame, 
    required_kir_genes:list, 
    pvalue_quality_threshold:float, 
    tech_rep_quality_threshold:float
) -> pd.DataFrame:

    # Log Transform the P Value to confirm with reported value
    pvalue_quality_threshold = -1 * math.log10(pvalue_quality_threshold)

    # Replace -1 with infinity, to facilitate filtering.
    mod_kir_from_gex = raw_kir_from_gex.copy()
    mod_kir_from_gex["quality_score"] = mod_kir_from_gex["quality_score"].astype(object)
    mod_kir_from_gex.loc[mod_kir_from_gex["allele_count"].astype(int) == 0, "quality_score"] = 60
    mod_kir_from_gex.loc[mod_kir_from_gex["quality_score"].astype(float) <= pvalue_quality_threshold, "allele_count"] = None

    records = []
    for subject_id in np.unique(raw_kir_from_gex["subject_id"]):
        record = {key:np.nan for key in required_kir_genes}
        subject_subset = mod_kir_from_gex[mod_kir_from_gex["subject_id"] == subject_id]
        for kir_gene in required_kir_genes:
            kir_gene_subset = subject_subset[subject_subset["kir_gene"].str.lower() == kir_gene].copy()
            record[kir_gene] = get_genotype_call_from_tech_reps(kir_gene_subset, tech_rep_quality_threshold)
        
        record["subject_id"] = subject_id
        record["quality_threshold"] = tech_rep_quality_threshold
        record["pval_threshold"] = pvalue_quality_threshold
        records.append(record)
    
    kir_genotypes = pd.DataFrame(records)

    return kir_genotypes

def compute_frequency_stats(kir_genotype_calls:pd.DataFrame, required_kir_genes:list) -> pd.DataFrame:
    frequencies = {key:0.0 for key in required_kir_genes}
    kir_genotype_calls = kir_genotype_calls.dropna()

    for ikir in required_kir_genes:
        frequency = (kir_genotype_calls[ikir] == True).sum()/kir_genotype_calls.shape[0]
        frequencies[ikir] = frequency

    frequencies = pd.DataFrame([frequencies])
    return frequencies

def format_and_export(
        kir_genotype_calls:pd.DataFrame, 
        kir_frequencies:pd.DataFrame, 
        required_kir_genes:list,
        output_filename:str, output_dir:str
):
    
    kir_genotype_calls = kir_genotype_calls.dropna()
    kir_genotype_calls = kir_genotype_calls[["subject_id", "quality_threshold", "pval_threshold"] + sorted(required_kir_genes)]

    output_filename = os.path.join(output_dir, output_filename)
    kir_genotype_calls.to_csv(output_filename)
    print(kir_genotype_calls)

    output_filename = "kir_genotype_frequencies.csv"
    output_filename = os.path.join(output_dir, output_filename)
    kir_frequencies.to_csv(output_filename)
    print(kir_frequencies)

def main():
    args = parse_arguments()
    os.makedirs(args["output_dir"], exist_ok=True)

    # Declare Constants
    required_kir_genes = [
        "kir2dl1", "kir2dl2", "kir2dl3", "kir3dl1", 
        "kir3ds1",  "kir2dl4", "kir3dl2"
    ]

    # Import Data
    raw_t1k_kir_results = import_tabular_file(args["input_file"])

    # Perform KIR typing across patients for all technical replicates
    kir_genotype_calls = call_kir_genotypes_from_gex(
        raw_kir_from_gex = raw_t1k_kir_results, 
        required_kir_genes = required_kir_genes,
        pvalue_quality_threshold = args["p_val_quality_threshold"],
        tech_rep_quality_threshold = args["tech_rep_quality_threshold"]
    )

    # Calculate Summary Stats
    kir_frequencies = compute_frequency_stats(
        kir_genotype_calls = kir_genotype_calls, 
        required_kir_genes = required_kir_genes
    )

    # Export Results
    format_and_export(
        kir_genotype_calls = kir_genotype_calls, 
        kir_frequencies = kir_frequencies, 
        required_kir_genes = required_kir_genes,
        output_filename = args["output_filename"], 
        output_dir = args["output_dir"]
    )

print("Starting...")
main()
print("Complete.")