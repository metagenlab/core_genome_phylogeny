

import pandas
import gzip
import shutil
import yaml 

pipeline_path = workflow.basedir





singularity_envs = yaml.safe_load(open(os.path.join(workflow.basedir,  "envs/sing_envs.yml"), 'r'))

sample_table = pandas.read_csv(config["input_table_local"], sep = '\t', header=0)


rule decompress_assemblies:
    input:
        '{any}.fna.gz'
    output:
        '{any}.fna'
    shell:
        "zcat {input[0]} > {output[0]}"

rule copy_local_genomes:
    output:
        # assembly_gz/genomes/local_1015432923.fna.gz
        expand("assembly_gz/genomes/local_{sample}.fna.gz", sample=sample_table["label"].to_list())
    run:
        for n, row in sample_table.iterrows():
            f_in = open(row["fna_path"], "rb")
            with gzip.open( f'assembly_gz/genomes/local_{row["label"]}.fna.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

# input function for the rule aggregate
def aggregate_genomes(wildcards):
    # decision based on content of output file
    df = pandas.read_csv(checkpoints.download_assemblies.get().output[0],
                         sep="\t", 
                         header=0,
                         index_col=0)
    assembly_list = [f"assembly_gz/genomes/{assembly}_genomic.fna" for assembly in df.AssemblyNames]
    # "InSilico/PCR/{sample}/Primers.fasta.mux.{sample}.fna.amplicons.nr"
    return assembly_list

rule merge_fasta:
    input:
        aggregate_genomes, expand("genomes/local_{sample}.fna", sample=sample_table["label"].to_list())
    output:
        "merged.fna"
    shell:
        """
        cat {input} > merged.fna
        """


include:
    "rules/assembly_finder.smk"
include:
    "rules/checkm.smk"

include:
    "rules/plots.smk" 