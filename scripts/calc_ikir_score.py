import argparse
import json
import os
import numpy as np
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
        "-kf",
        "--kir_file",
        help="Filename of the KIR genotyping results. Full path filename as .csv.",
        required=True,
        type=str,
    )

    required.add_argument(
        "-hf",
        "--hla_file",
        help="Filename of the HLA allele typing results. Full path filename as .csv.",
        required=True,
        type=str,
    )

    required.add_argument(
        "-rf",
        "--ref_file",
        help="Filename of hla allele definition file. Full path filename as .csv.",
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
        default="functional_ikir_scoring.csv"
    )
    
    arguments = vars(parser.parse_args())
    return arguments

def load_json(filename:str) -> dict:
    payload = {}
    with open(filename, 'r') as file_o: 
        raw_text = file_o.read()
        payload = json.loads(raw_text)
    return payload

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

def get_pertinent_motif_status(allele_name:str, hla_allele_definitions:pd.DataFrame) -> dict:
    desired_motifs = ["c1", "c2", "bw4", "bw6"]
    if allele_name and not pd.isna(allele_name):
        motif_statuses = hla_allele_definitions[hla_allele_definitions["ebi_name"].str.startswith(allele_name)].iloc[0][desired_motifs].to_dict()    
    else:
        motif_statuses = {motif_name:None for motif_name in desired_motifs}
    return motif_statuses

def evaluate_ligand_match_criteria(relevant_alleles:dict, ligand_match_defintion:dict, hla_allele_definitions:pd.DataFrame) -> bool:
    
    hla_gene_letter = ligand_match_defintion["gene"][0].lower()
    keys = ["{0}_{1}".format(hla_gene_letter, i) for i in [1,2]]

    # If no allele information, return None
    if relevant_alleles[keys[0]] is None and relevant_alleles[keys[1]] is None:
        return None
    
    # Loop through pertinent alleles, if match return true, else false
    match_status = []
    for key in keys:
        allele_name = relevant_alleles[key]
        if allele_name:
            allele_definition = get_pertinent_motif_status(
                allele_name=allele_name, 
                hla_allele_definitions=hla_allele_definitions
            )
            required_allele_prefix = ligand_match_defintion["gene"]
            required_motif = ligand_match_defintion["motif"]

            if allele_name and not pd.isna(allele_name) and allele_name.startswith(required_allele_prefix) and allele_definition[required_motif]:
                match_status.append(True)
            else:
                match_status.append(False)
        else:
            match_status.append(None)
    
    if match_status[0] or match_status[1]:
        return True
    elif match_status[0] is None or match_status[1] is None:
        return None
    else:
        return False

def add_ligand_possession_status(
    hla_allele_calls:pd.DataFrame, 
    hla_allele_definitions:pd.DataFrame, 
    ligand_matching_criteria:dict
) -> pd.DataFrame:
    
    hla_genes = ["a_1", "a_2", "b_1", "b_2", "c_1", "c_2"]
    ligand_possession_statuses = []
    for index, row in hla_allele_calls.iterrows():
        ligand_possession_status = {}
        
        relevant_alleles = row[hla_genes].to_dict()

        for ligand_name in ligand_matching_criteria:
            ligand_match_defintion = ligand_matching_criteria[ligand_name]
            ligand_possession_status[ligand_name] = evaluate_ligand_match_criteria(
                relevant_alleles = relevant_alleles, 
                ligand_match_defintion = ligand_match_defintion,
                hla_allele_definitions = hla_allele_definitions,
            )
        ligand_possession_status["subject_id"] = row["subject_id"]
        ligand_possession_statuses.append(ligand_possession_status)

    ligand_possession_statuses = pd.DataFrame(ligand_possession_statuses)
    hla_allele_calls = hla_allele_calls.merge(ligand_possession_statuses, on="subject_id", how="left")

    return hla_allele_calls

def get_functional_ligand_affinity(row, matching_criteria:list) -> float:
    none_flag = False
    for matching_criterion in matching_criteria:
        required_ligand = matching_criterion["ligand"]
        ligand_status = row[required_ligand]
        ligation_affinity = matching_criterion["affinity"]
        if ligand_status is None:
             none_flag = True
        elif ligand_status:
            return ligation_affinity
    if none_flag:
        return None
    else:
        return 0.0

def has_nones(input_dict:dict) -> bool:
    for key,value in input_dict.items():
        if value is None:
            return True
    return False

def is_functional(row, kir_gene:str, matching_criteria:list) -> bool:
    none_flag = False
    for matching_criterion in matching_criteria:
        kir_status = row[kir_gene]
        if kir_status is None:
            return None
        elif kir_status:
            required_ligand = matching_criterion["ligand"]
            ligand_status = row[required_ligand]
            if ligand_status is None:
                none_flag = True
            elif ligand_status:
                return True
            else:
                continue
        else:
            return False
        
    # Finally, if any misisng information, can't be sure of true negative, so return None.
    if none_flag:
        return None
    else:
        return False

def compute_functional_kir_status(genotype_information:pd.DataFrame, kir_ligand_matching_criteria:dict) -> pd.DataFrame:
    all_functional_kir_status = []
    for index, row in genotype_information.iterrows():

        functional_kir_statuses = {}
        for kir_gene in kir_ligand_matching_criteria:
            relevant_matching_criteria = kir_ligand_matching_criteria[kir_gene]
            func_kir_key = "func_{0}".format(kir_gene)
            functional_kir_statuses[func_kir_key] = is_functional(row, kir_gene, relevant_matching_criteria)
        functional_kir_statuses["subject_id"] = row["subject_id"]
        all_functional_kir_status.append(functional_kir_statuses)

    # Combine Results     
    all_functional_kir_status = pd.DataFrame(all_functional_kir_status)
    genotype_information = genotype_information.merge(all_functional_kir_status, how="left", on="subject_id").copy()

    return genotype_information

def compute_functional_ikir_stats(genotype_information:pd.DataFrame, kir_ligand_matching_criteria:dict) -> pd.DataFrame:
    ikir_genes = kir_ligand_matching_criteria.keys()
    functional_ikirs = ["func_{0}".format(kir_gene) for kir_gene in ikir_genes]
    records = []
    for index, row in genotype_information.iterrows():
        ikir_stats = {}
        kir_statuses = row[ikir_genes].to_dict()
        if has_nones(kir_statuses):
            ikir_stats["ikir_count"] = None
        else:
            ikir_stats["ikir_count"] = row[ikir_genes].values.sum()

        functional_ikir_statuses = row[functional_ikirs].to_dict()
        if has_nones(functional_ikir_statuses):
            ikir_stats["func_ikir_score"] = None
            ikir_stats["func_ikir_count"] = None
        else:
            ikir_stats["func_ikir_count"] = row[functional_ikirs].values.sum()
            func_ikir_score = 0.0

            for func_ikir_gene in functional_ikir_statuses:
                if functional_ikir_statuses[func_ikir_gene]:
                    func_ikir_score += get_functional_ligand_affinity(row, kir_ligand_matching_criteria[func_ikir_gene[5:]])
            ikir_stats["func_ikir_score"] = func_ikir_score
        ikir_stats["subject_id"] = row["subject_id"]
        records.append(ikir_stats)

    ikir_stats = pd.DataFrame(records)

    genotype_information = genotype_information.merge(ikir_stats, how="left", on="subject_id")
    return genotype_information

def format_and_export(genotype_information:pd.DataFrame, output_filename:str, output_dir:str):
    
    item_to_move_to_front  = ["subject_id", "func_ikir_score", "func_ikir_count", "ikir_count"] 
    columns = list(genotype_information.columns)
    for key in item_to_move_to_front:
        columns.remove(key)

    columns = item_to_move_to_front + sorted(columns)
    genotype_information = genotype_information[columns]
    columns = [column_name.lower() for column_name in columns]
    genotype_information = pd.DataFrame(genotype_information.values, columns=columns)
    
    output_filename = os.path.join(output_dir, output_filename)
    genotype_information.to_csv(output_filename)
    print(genotype_information)

def main():
    args = parse_arguments()
    os.makedirs(args["output_dir"], exist_ok=True)

    # Load Data
    kir_ligand_matching_criteria = load_json("ref_data/kir_ligand_matching_criteria.json")
    ligand_matching_criteria = load_json("ref_data/hla_ligand_definitions.json")
    kir_genotype_calls = import_tabular_file(args["kir_file"], index_col=0)
    hla_allele_calls = import_tabular_file(args["hla_file"], index_col=0)
    hla_allele_definitions = import_tabular_file(args["ref_file"], index_col=0)

    # Agument Data with Reference Data (Compute Motif Posession)
    hla_allele_calls = add_ligand_possession_status(
        hla_allele_calls=hla_allele_calls, 
        hla_allele_definitions=hla_allele_definitions, 
        ligand_matching_criteria=ligand_matching_criteria
    )

    # Merge Data
    genotype_information = kir_genotype_calls.merge(
        right=hla_allele_calls, how="inner", on="subject_id"
    )

    # Compute iKIR Score
    genotype_information = compute_functional_kir_status(
        genotype_information = genotype_information, 
        kir_ligand_matching_criteria = kir_ligand_matching_criteria
    )

    # Compute iKIR Count and iKIR Count
    genotype_information = compute_functional_ikir_stats(
        genotype_information = genotype_information, 
        kir_ligand_matching_criteria = kir_ligand_matching_criteria
    )

    # Export Results
    format_and_export(
        genotype_information = genotype_information,
        output_filename = args["output_filename"],
        output_dir = args["output_dir"] 
    )

print("Starting...")
main()
print("Complete.")