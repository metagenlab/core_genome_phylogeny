
rule checkm_analyse:
    conda:
        "../envs/checkm.yml"
    threads:
        8
    resources:
        load_checkm=100
    params:
        markers = os.path.join(pipeline_path, 'data/bacteria.ms')
    input:
        local = expand("assembly_gz/genomes/local_{sample}.fna", sample=sample_table["label"].to_list()),
        downloaded = aggregate_genomes, 
    output:
        log = "checkm_analyse/checkm_log.txt",
    shell:
        """
        checkm data setRoot /data/databases/checkm
        checkm  analyze -x fna {params.markers} assembly_gz/genomes/ checkm_analyse -t 8 --nt > {output[0]}_save
        echo OK > {output[0]}
        # checkm analyze --genes -x faa_gz bacteria.ms faa output -t 8 --nt
        """

checkpoint checkm_gene_fasta_list:
    input:
        "checkm_analyse/checkm_log.txt",
    output:
        "checkm_analyse/bins/bin_list.tab",
    shell:
        """
        ls checkm_analyse/bins/*/genes.faa >> {output[0]}
        """

rule checkm_qa:
    conda:
        "../envs/checkm.yml"
    threads:
        8
    params:
        markers = os.path.join(pipeline_path, 'data/bacteria.ms')
    input:
        checkm_log = "checkm_analyse/checkm_log.txt"
    output:
        qa_results = "checkm_analyse/qa_results.txt"
    shell:
        """
        checkm data setRoot /data/databases/checkm
        checkm qa {params[0]} checkm_analyse -o 2 --tab_table > {output[0]}
        """

rule checkm_marker_set:
    conda:
        "../envs/checkm.yml"
    threads:
        8
    params:
        markers = os.path.join(pipeline_path, 'data/bacteria.ms')
    input:
        qa_results = "checkm_analyse/qa_results.txt"
    output:
        checkm_table = "checkm_analyse/markers.tab"
    shell:
        """
        checkm data setRoot /data/databases/checkm
        checkm qa {params[0]} checkm_analyse -o 5 --tab_table > {output[0]}
        """

def aggregate_checkm_bins(wildcards):
    
    lst = []
    with open(checkpoints.checkm_gene_fasta_list.get().output[0]) as f:
        for row in f:
            lst.append(row.strip())
    return lst
    '''
    print("ok")

    checkpoint_output = checkpoints.checkm_analyse.get(**wildcards).output[1]

    print("checkpoint_output", checkpoint_output)

    return expand("checkm_analyse/bins/{i}/genes.faa",
           taxid=wildcards.taxid,
           i=glob_wildcards(os.path.join(checkpoint_output, "{i}/genes.faa")).i)
    '''


checkpoint checkm_markers_fastas:
    conda:
        "../envs/pandas-openpyxl-pronto-xlrd.yml"
    threads:
        1
    input:
        checkm_fastas = aggregate_checkm_bins,
        checkm_table = "checkm_analyse/markers.tab"
    params:
        markers = os.path.join(pipeline_path, 'data/bacteria.ms')
    output:
        marker_fastas = directory("checkm_marker_fastas")
    script: "scripts/get_checkm_markers_fastas.py"


rule align_marker_fasta_with_mafft:
    conda:
        "../envs/mafft.yml"
    input:
        fasta = "checkm_marker_fastas/{marker}.faa"
    output:
        alignment = "checkm_marker_alignments/{marker}_mafft.faa"
    shell:
        """
        mafft --quiet --anysymbol --amino --auto --maxiterate 1000 {input[0]} > {output[0]}
        """

def aggregate_checkm_markers_alignments(wildcards):
    checkpoint_output = checkpoints.checkm_markers_fastas.get(**wildcards).output[0]
    print("checkpoint_output", checkpoint_output)
    return expand("checkm_marker_alignments/{i}_mafft.faa",i=glob_wildcards(os.path.join(checkpoint_output, "{i}.faa")).i)

rule concatenate_checkm_fastas:
    conda:
        "../envs/pandas-openpyxl-pronto-xlrd.yml"
    input:
        marker_fastas = aggregate_checkm_markers_alignments
    output:
        alignment = "checkm_phylogeny/concatenated_alignment.faa"
    script: "scripts/concat_align.py"


rule build_phylogeny_with_fasttree:
    conda:
        "../envs/fasttree.yml"
    input:
        alignment = "checkm_phylogeny/concatenated_alignment.faa"
    output:
        tree = "checkm_phylogeny/tree.nwk"
    shell:
        """
        fasttree {input[0]} > {output[0]}
        """
