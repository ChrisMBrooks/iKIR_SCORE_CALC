rule calc_ikir_score:
    input:
        script = "scripts/calc_ikir_score.py",
        kir_file = "outputs/{project}/kir/kir_genotype_calls.csv",
        hla_file = "outputs/{project}/hla/hla_allele_calls.csv",
        ref_file = "outputs/{project}/hla/ipd_hla_allele_definitions.csv"
    params:
        output_dir = "outputs/{project}",
        output_filename = "functional_ikir_scoring.csv"
    output:
        output_filename = "outputs/{project}/functional_ikir_scoring.csv"
    conda: "../envs/pandas_env.yaml"
    shell:
        """
            python {input.script} \
            --kir_file {input.kir_file} \
            --hla_file {input.hla_file} \
            --ref_file {input.ref_file} \
            --output_filename {params.output_filename} \
            --output_dir {params.output_dir} 
        """