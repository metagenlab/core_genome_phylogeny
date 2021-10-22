

checkpoint download_assemblies:
    conda:
        "../envs/af.yml"
    output:
        "genomes-assemblies-summary.tsv"
    params:
        input_table = config["input_table_af"],
        conda_prefix = "/home/tpillone/miniconda3",
        ncbi_key = config["NCBI_key"],
        ncbi_email = config["NCBI_email"],
        complete_assemblies = config["complete_assemblies"],
        filter_rank = config["Rank_to_filter_by"],
        n_by_rank = config["n_by_rank"]
    shell:
        """
        af run --nolock --input-table {params.input_table} --output_prefix genomes --conda-prefix {params.conda_prefix} --ncbi_key {params.ncbi_key} --ncbi_email {params.ncbi_email} --complete_assemblies {params.complete_assemblies} --filter_rank {params.filter_rank} --n_by_rank {params.n_by_rank}
        """
