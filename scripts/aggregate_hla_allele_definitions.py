import argparse
import json
import os
import re
import requests
import time
import numpy as np
import pandas as pd

def parse_arguments() -> dict:
    parser = argparse.ArgumentParser(
        description="Script to retrieve the hla allele details required to compute the iKIR score.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-i",
        "--input_file",
        help="Input configuration filename. Full path filename as .parquet.",
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
        default="ipd_hla_allele_definitions.csv"
    )
    
    arguments = vars(parser.parse_args())
    return arguments

def load_json(filename:str) -> dict:
    payload = {}
    with open(filename, 'r') as file_o: 
        raw_text = file_o.read()
        payload = json.loads(raw_text)
    return payload

def format_hla_allele_name(allele_tuple:list) -> str:
    short_code_prefix = str(int(allele_tuple[1]))

    if len(short_code_prefix) < 2:
        allele_tuple[1] = "0{0}".format(short_code_prefix)

    allele_name = "{0}*{1}:{2}".format(
        allele_tuple[0], 
        allele_tuple[1], 
        allele_tuple[2]
    )

    return allele_name

def get_hla_allele_options(allele_name:str) -> list:
    url_template = "https://www.ebi.ac.uk/cgi-bin/ipd/api/allele?limit=10&project=HLA&query=startsWith(name,\"{}\")"
    url = url_template.format(allele_name)

    request = requests.get(url)
    
    records = []
    
    if request.status_code == 200:
        payload = request.json()
        if len(payload["data"]) > 0:
            for list_object in payload["data"]:
                if "accession" in list_object and "name" in list_object:
                    record = {
                        "ebi_id": list_object['accession'],
                        "ebi_name": list_object['name']
                    }
                    records.append(record)
    return records

def get_hla_prot_sequence(ebi_id:str) -> str:

    url_template = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=imgthlapro;id={}&format=fasta&style=raw"
    url = url_template.format(ebi_id)

    ebi_hla_prot_request = requests.get(url)

    raw_strand = ""
    if ebi_hla_prot_request.status_code == 200:
        raw_fasta_prot_txt = ebi_hla_prot_request.text

        regex_defintion = r"(>)(\w+:)(.+\n)([ACDEFGHIKLMNPQRSTVWY\n]+)"
        uniprot_parser = re.compile(regex_defintion)
        matches = uniprot_parser.findall(raw_fasta_prot_txt)

        if len(matches) > 0:
            raw_strand=matches[0][-1].replace('\n','').strip()
    
    return raw_strand

def get_signal_peptide_length(ebi_id:str) -> int:

    url_template = 'https://www.ebi.ac.uk/cgi-bin/ipd/api/allele/{}'

    url = url_template.format(ebi_id)
    request = requests.get(url)

    if request.status_code == 200:
        payload = request.json()
        if "feature" in payload:
            peptides = payload['feature']['protein']
            signal_peptide_lengths = [peptide['length'] for peptide in peptides if peptide['type'] == 'signal']
            if len(signal_peptide_lengths) > 0:
                return int(signal_peptide_lengths[0])

    print("No signal peptide length data available for {0}. Assuming length of 24aa.".format(ebi_id))
    return 24

def get_unique_hla_alleles(hla_alleles:pd.DataFrame) -> list:

    hla_alleles = np.unique(hla_alleles["hla_allele"].dropna())

    unique_allele_tuples = []
    for hla_allele in hla_alleles:
        star_pieces = hla_allele.split("*")
        colon_pieces = star_pieces[1].split(":")
        unique_allele_tuples.append([star_pieces[0].upper(), colon_pieces[0], colon_pieces[1]])

    return unique_allele_tuples

def is_conventional_allele(ebi_name:str) -> bool:
    alternative_expressions = ['N', 'L', 'S', 'C', 'A', 'Q']
    for letter in alternative_expressions:
        if letter in ebi_name[1:]:
            return False
    return True

def find_record_from_options(ebi_hla_options, allele_tuple:tuple) -> list:
    for ebi_hla_option in ebi_hla_options:
        if "ebi_id" in ebi_hla_option and is_conventional_allele(ebi_hla_option["ebi_name"]):
            ebi_hla_option["prot_seq"] = get_hla_prot_sequence(ebi_hla_option["ebi_id"])
            ebi_hla_option["signal_peptide_len"] = get_signal_peptide_length(ebi_hla_option["ebi_id"])

            if ebi_hla_option["prot_seq"]:
                    
                ebi_hla_option["hla_gene"] = allele_tuple[0]
                ebi_hla_option["short_code"] = "{0}{1}".format(allele_tuple[1],allele_tuple[2])

                record = [
                    ebi_hla_option["ebi_id"], 
                    ebi_hla_option["ebi_name"], 
                    ebi_hla_option["hla_gene"], 
                    ebi_hla_option["short_code"], 
                    ebi_hla_option["prot_seq"],
                    ebi_hla_option["signal_peptide_len"]
                ]

                return record
    return None

def compute_motif_posession(
    protein_sequence:str, 
    motif_definitions:dict, 
    offset = 24
) -> dict:

    motif_status = {}
    for key in motif_definitions:
        criterea = motif_definitions[key]
        motif_status[key] = True
        for critereon in criterea:
            amino_acid = protein_sequence[critereon[0]-1+offset]
            condition = amino_acid in critereon[1]
            if not condition:
                motif_status[key] = False
                break
            else: 
                pass
    return motif_status

def retrieve_hla_meta_data_from_ipd(alleles, motif_definitions:dict) -> pd.DataFrame:
    records = []
    failed_records = []
    for allele_tuple in alleles:
        allele_name = format_hla_allele_name(allele_tuple = allele_tuple)
        ebi_hla_options = get_hla_allele_options(allele_name = allele_name)
        ebi_record = find_record_from_options(
            ebi_hla_options = ebi_hla_options, 
            allele_tuple = allele_tuple,
        )

        # Compute Motif Posession
        protein_sequence = ebi_record[4]
        signal_peptide_len = ebi_record[5]
        motif_posession_statuses = compute_motif_posession(
            protein_sequence, 
            motif_definitions = motif_definitions, 
            offset = signal_peptide_len
        )

        ebi_record = ebi_record + [motif_posession_statuses[key] for key in motif_definitions]

        if ebi_record:
            records.append(ebi_record)
        else:
            print("Failed to find record for {}".format(allele_tuple))
            ebi_record = [allele_tuple[0], "{0}{1}".format(allele_tuple[1],allele_tuple[2]), allele_name]
            failed_records.append(ebi_record)
        time.sleep(.1)
    
    column_names = ["ebi_id", "ebi_name", "hla_gene", "short_code", "protein_sequence", "signal_peptide_len"]
    column_names = column_names + list(motif_definitions.keys())
    hla_allele_records = pd.DataFrame(records, columns=column_names)
    hla_allele_records["short_code"] = hla_allele_records["short_code"].astype('int')

    failed_records_df = pd.DataFrame(
        failed_records, 
        columns=["hla_gene", "short_code", "allele_name"]
    )

    return hla_allele_records, failed_records_df

def format_and_export(hla_allele_records:pd.DataFrame, failed_records:pd.DataFrame, output_filename:str, output_dir:str):
    output_filename = os.path.join(output_dir, output_filename)
    hla_allele_records.to_csv(output_filename)
    print(hla_allele_records)

    if failed_records.shape[0] > 0:
        print("Some HLA alleles could not be reconciled against the ipd. Details are reported in the failed_records.csv file.")
        output_filename = os.path.join(output_dir, "failed_records.csv")
        failed_records.to_csv(output_filename)

def main():
    args = parse_arguments()
    os.makedirs(args["output_dir"], exist_ok=True)

    # Import Data
    motif_definitions = load_json("ref_data/hla_ligand_motif_definitions.json")
    hla_allele_typing = pd.read_parquet(args["input_file"])

    # Format Input Data
    alleles = get_unique_hla_alleles(hla_allele_typing)

    # Retrieve Meta Data from IPD
    hla_allele_records, failed_records = retrieve_hla_meta_data_from_ipd(
        alleles = alleles, 
        motif_definitions = motif_definitions
    )

    # Export Results
    format_and_export(
        hla_allele_records = hla_allele_records, 
        failed_records = failed_records, 
        output_filename = args["output_filename"], 
        output_dir = args["output_dir"]
    )

print('Starting...')
main()
print('Complete.')

