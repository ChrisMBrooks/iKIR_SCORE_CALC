# snakemake --cores 8 --use-conda --conda-frontend conda --printshellcmds
# snakemake --cores 8 --use-conda --conda-frontend conda --rerun-incomplete
# snakemake --cores 8 --use-conda --conda-frontend conda --dry-run --printshellcmds
# snakemake --cores 8 --use-conda --conda-frontend conda --keep-incomplete
# snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf

import json, os, sys
def load_pipeline_config(config_filename:str) -> dict:
    try:
        full_path = os.path.join(os.getcwd(), config_filename)
        f = open(full_path)
        config = json.load(f)
        return config
    except Exception as e:
        print('Failed to load pipeline config. Exiting...')
        sys.exit(1)

# Import Pipeline Configuration File
CONFIG_FILENAME = "pipeline.config.json"
CONFIG = load_pipeline_config(CONFIG_FILENAME)

# Declare Constants
PROJECT = CONFIG["project"]
RAW_HLA_INPUT_FILENAME = CONFIG["hla_allele_calling"]["raw_input_filenames"][0]
RAW_KIR_INPUT_FILENAME = CONFIG["kir_genotype_calling"]["raw_input_filenames"][0]

# Master Rule

rule ikir_score_ready:
    input:
        config_file = CONFIG_FILENAME, 
        ikir_score_file = expand(
            "outputs/{project}/functional_ikir_scoring.csv",
            project = PROJECT
        )

# Static Imports

include: "rules/aggregate_hla_allele_definitions.smk"
include: "rules/calc_ikir_score.smk"

# Conditional Imports

if CONFIG["hla_allele_calling"]["source"].lower() == CONFIG["rule_config"]["hla_imp_calling"]["name"].lower():
    include: "rules/call_hla_alleles_from_hla_imp.smk"
elif CONFIG["hla_allele_calling"]["source"].lower() == CONFIG["rule_config"]["optitype_calling"]["name"].lower():
    include: "rules/call_hla_alleles_from_optitype.smk"

if CONFIG["kir_genotype_calling"]["source"].lower() == CONFIG["rule_config"]["kir_imp_calling"]["name"].lower():
    include: "rules/call_kir_genotypes_from_kir_imp.smk"
elif CONFIG["kir_genotype_calling"]["source"].lower() == CONFIG["rule_config"]["t1k_calling"]["name"].lower():
    include: "rules/call_kir_genotype_from_t1k.smk"


