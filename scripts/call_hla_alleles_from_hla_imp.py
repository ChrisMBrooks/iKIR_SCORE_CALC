import argparse
import os
import numpy as np
import pandas as pd

def parse_arguments() -> dict:
    parser = argparse.ArgumentParser(
        description="Script to call hla allele variants from HLA*IMP3 results (e.g. imputations.csv).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-if",
        "--input_file",
        help="Filename for the HLA*IMP3 results file (imputations.csv). Full path filename as .csv.",
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

def format_hla_shortcode(shortcode:str, gene_letter:str) -> str:
    if len(shortcode) < 4:
        shortcode = "0{0}:{1}".format(shortcode[0], shortcode[1:])
    else:
        shortcode = "{0}:{1}".format(shortcode[0:2], shortcode[2:])

    hla_allele = "{0}*{1}".format(gene_letter.upper(), shortcode)
    return hla_allele

def call_haplotype_allele(haplotype:pd.Series, quality_threshold:float, gene_letter:str) -> bool:
    if haplotype["imputed_type"] is not None:
        if float(haplotype["posterior_probability"]) > quality_threshold:
            hla_allele = format_hla_shortcode(str(haplotype["imputed_type"]), gene_letter)
            return hla_allele
    return None

def call_hla_genotype(raw_hla_imp_results:pd.DataFrame, required_hla_genes:list, quality_threshold:float) -> pd.DataFrame:
    records = []
    for subject_id in np.unique(raw_hla_imp_results["subject_id"]):
        hla_allele_calls = {}
        subject_subset = raw_hla_imp_results[raw_hla_imp_results["subject_id"] == subject_id]
        for hla_gene in required_hla_genes:
            gene_letter = hla_gene.replace("HLA","").lower()
            hap_key_1 = "{0}_{1}".format(gene_letter, 1)
            hap_key_2 = "{0}_{1}".format(gene_letter, 2)

            hla_subset = subject_subset[subject_subset["locus"] == hla_gene]
            hla_allele_calls[hap_key_1] = call_haplotype_allele(hla_subset.iloc[0], quality_threshold, gene_letter)
            hla_allele_calls[hap_key_2] = call_haplotype_allele(hla_subset.iloc[1], quality_threshold, gene_letter)
        
        hla_allele_calls["subject_id"] = subject_id
        hla_allele_calls["applied_threshold"] = quality_threshold
        records.append(hla_allele_calls)
    hla_allele_calls = pd.DataFrame(records)
    return hla_allele_calls

def compute_frequency_stats(hla_allele_calls:pd.DataFrame) -> pd.DataFrame:
    
    all_hla_alleles = ["a_1", "a_2", "b_1", "b_2", "c_1", "c_2"]
    all_hla_alleles = [hla_allele_calls[key].dropna().values for key in all_hla_alleles]
    all_hla_alleles = np.hstack(all_hla_alleles)
    desired_hla_alleles = np.unique(all_hla_alleles)

    frequencies = {key:0.0 for key in desired_hla_alleles}

    for allele_name in all_hla_alleles:
        frequencies[allele_name] +=1/(hla_allele_calls.shape[0]*2)

    for allele_name in frequencies:
         frequencies[allele_name] = 1 - (1.0-frequencies[allele_name])**2
    
    frequencies = pd.DataFrame([frequencies]).iloc[0].sort_values(ascending=False)
    return frequencies

def format_and_export_results(hla_allele_calls:pd.DataFrame, hla_frequencies:pd.DataFrame, output_filename:str, output_dir:str):
    # Format and Export Primary Results
    output_filename = "hla_allele_calls.csv"
    output_filename = os.path.join(output_dir, output_filename)
    columns = list(hla_allele_calls.columns)
    columns.remove("subject_id")
    columns.remove("applied_threshold")
    columns = sorted(columns)
    columns = ["subject_id", "applied_threshold"] + sorted(columns)
    hla_allele_calls = hla_allele_calls[columns]
    hla_allele_calls.to_csv(output_filename)
    print(hla_allele_calls)

    # Export Frequencies
    output_filename = "hla_phenotype_frequencies.csv"
    output_filename = os.path.join(output_dir, output_filename)
    hla_frequencies.to_csv(output_filename)
    print("hla Population Frequencies:")
    print(hla_frequencies)

def main():
    # Parse Inputs
    args = parse_arguments()
    os.makedirs(args["output_dir"], exist_ok=True)

    # Declare Constants
    required_hla_genes = [
        'HLAA', 'HLAB', 'HLAC', 
        'HLADPA1', 'HLADPB1', 'HLADQA1', 'HLADQB1', 
        'HLADRB1','HLADRB3', 'HLADRB4', 'HLADRB5'
    ]

    # Load Data
    raw_hla_imp_results = pd.read_csv(args["input_file"])
    columns = ['subject_id', 'subject_id2', 'haplotype_id', 'locus', 'imputed_type', 'posterior_probability']
    raw_hla_imp_results = pd.DataFrame(raw_hla_imp_results.values, columns=columns)

    # Call Alleles
    hla_allele_calls = call_hla_genotype(
        raw_hla_imp_results = raw_hla_imp_results, 
        required_hla_genes = required_hla_genes, 
        quality_threshold = args["quality_threshold"]
    )

    # Compute Stats
    hla_allele_frequencies = compute_frequency_stats(hla_allele_calls)

    # Export Rsults
    format_and_export_results(
        hla_allele_calls = hla_allele_calls, 
        hla_frequencies = hla_allele_frequencies, 
        output_filename = args["output_filename"],
        output_dir = args["output_dir"]
    )

print("Starting...")
main()
print("Complete.")
