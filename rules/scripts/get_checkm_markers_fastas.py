
import pandas
import re
from Bio import SeqIO
import os

checkm_table = snakemake.input["checkm_table"]
checkm_fastas = snakemake.input["checkm_fastas"]
marker_set = snakemake.params["markers"]

print("checkm_table", checkm_table)
print("checkm_fastas", checkm_fastas)
print("marker_set", marker_set)

os.mkdir("checkm_marker_fastas")

def parse_checkm_marker_table(checkm_table):
    marker_table = pandas.read_csv(checkm_table, delimiter='\t', header=0, index_col=0)

    marker2sample2records = {}

    for n, row in marker_table.iterrows():
        sample_id = row.name
        marker = row[0]
        gene_id = row[1]
        if marker not in marker2sample2records:
            marker2sample2records[marker] = {}
        if sample_id not in marker2sample2records[marker]:
            marker2sample2records[marker][sample_id] = [gene_id]
        else:
            marker2sample2records[marker][sample_id].append(gene_id)
    # remove hits with multiple copies
    for marker in marker2sample2records:
        for sample in list(marker2sample2records[marker].keys()):
            if len(marker2sample2records[marker][sample]) > 1:
                marker2sample2records[marker].pop(sample)
    return marker2sample2records




def parse_fastas(fasta_files):

    sample2genes = {}
    for fasta_file in fasta_files:
        print("one fasta", fasta_file)
        sample = fasta_file.split("/")[2]
        sample_dico = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        sample2genes[sample] = sample_dico
    return sample2genes

marker_genes = parse_checkm_marker_table(checkm_table)

print("marker_genes", marker_genes)

fasta_dico = parse_fastas(checkm_fastas)


marker2records = {}
for marker in marker_genes:
    marker_records = []
    for sample in marker_genes[marker]:
        record = fasta_dico[str(sample)][marker_genes[marker][sample][0]]
        record.name = str(sample)
        record.id = str(sample)
        record.description = ""
        record.seq = record.seq[0:-1]
        marker_records.append(record)
        # expand("databases/virulence/complete_genomes/{{taxid}}/checkm_marker_fastas/{marker}.faa", markler = mylist)
    with open("checkm_marker_fastas/%s.faa" % (marker), "w") as f:
        SeqIO.write(marker_records, f, "fasta")
