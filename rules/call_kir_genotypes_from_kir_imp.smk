rule call_kir_genotypes_from_kir_imp:
    input:
        script = "scripts/call_kir_genotypes_from_kir_imp.py",
        kir_file = "{0}".format(RAW_KIR_INPUT_FILENAME)
    params:
        output_dir = "outputs/{project}/kir",
        output_filename = "kir_genotype_calls.csv",
        quality_threshold = CONFIG["rule_config"]["kir_imp_calling"]["quality_threshold"]
    output:
        output_filename = "outputs/{project}/kir/kir_genotype_calls.csv"
    conda: "../envs/pandas_env.yaml"
    shell:
        """
            python {input.script} \
            --input_file {input.kir_file} \
            --output_filename {params.output_filename} \
            --output_dir {params.output_dir} \
            --quality_threshold {params.quality_threshold}
        """
