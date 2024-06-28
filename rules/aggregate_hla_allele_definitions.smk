rule aggregate_hla_allele_definitions:
    input:
        script = "scripts/aggregate_hla_allele_definitions.py",
        hla_file = "{0}".format(RAW_HLA_INPUT_FILENAME)
    params:
        output_dir = "outputs/{project}/hla",
        output_filename = "ipd_hla_allele_definitions.csv"
    output:
        output_filename = "outputs/{project}/hla/ipd_hla_allele_definitions.csv"
    conda: "../envs/pandas_env.yaml"
    shell:
        """
            python {input.script} \
            --input_file {input.hla_file} \
            --output_filename {params.output_filename} \
            --output_dir {params.output_dir} 
        """
