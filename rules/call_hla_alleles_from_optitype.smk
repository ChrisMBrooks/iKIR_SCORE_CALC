rule call_hla_alleles_from_optitype:
    input:
        script = "scripts/call_hla_alleles_from_optitype.py",
        hla_file = "{0}".format(RAW_HLA_INPUT_FILENAME)
    params:
        output_dir = "outputs/{project}/hla",
        output_filename = "hla_allele_calls.csv",
        quality_threshold = CONFIG["rule_config"]["optitype_calling"]["quality_threshold"]
    output:
        output_filename = "outputs/{project}/hla/hla_allele_calls.csv"
    conda: "../envs/pandas_env.yaml"
    shell:
        """
            python {input.script} \
            --input_file {input.hla_file} \
            --output_filename {params.output_filename} \
            --output_dir {params.output_dir} \
            --quality_threshold {params.quality_threshold}
        """
