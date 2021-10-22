
rule plot_phylogeny: 
    input:
        assembly_table = "genomes-assemblies-summary.tsv",
        tree = "checkm_phylogeny/tree.nwk",
        local = expand("assembly_gz/genomes/local_{sample}.fna", sample=sample_table["label"].to_list()),
        downloaded = aggregate_genomes, 
    output:
        "plots/phylo.svg"
    script: "scripts/plot_phylo.py"
