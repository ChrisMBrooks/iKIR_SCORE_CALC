rule call_kir_genotype_from_t1k:
    input:
        script = "scripts/call_kir_genotype_from_t1k.py",
        kir_file = "{0}".format(RAW_KIR_INPUT_FILENAME)
    params:
        output_dir = "outputs/{project}/kir",
        output_filename = "kir_genotype_calls.csv",
        p_val_quality_threshold = CONFIG["rule_config"]["t1k_calling"]["p_val_quality_threshold"],
        tech_rep_quality_threshold = CONFIG["rule_config"]["t1k_calling"]["tech_rep_quality_threshold"]
    output:
        output_filename = "outputs/{project}/kir/kir_genotype_calls.csv"
    conda: "../envs/pandas_env.yaml"
    shell:
        """
            python {input.script} \
            --input_file {input.kir_file} \
            --output_filename {params.output_filename} \
            --output_dir {params.output_dir} \
            --p_val_quality_threshold {params.p_val_quality_threshold} \
            --tech_rep_quality_threshold {params.tech_rep_quality_threshold}
        """
