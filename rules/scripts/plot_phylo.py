import pandas
import os 
from Bio import SeqIO 
from Bio.SeqUtils import GC

genomes = snakemake.input["downloaded"] + snakemake.input["local"]

leaf2GC = {}
leaf2genome_size = {}
for genome in genomes:
    print("local genome"), genome
    concat_seq = ''
    records = SeqIO.parse(genome, "fasta")
    for record in records:
        concat_seq+=record.seq 
    gc = GC(concat_seq)
    leaf2GC[genome.split(".fna")[0].split("/")[-1].split("_genomic")[0]] = gc
    leaf2genome_size[genome.split(".fna")[0].split("/")[-1].split("_genomic")[0]] = round(len(concat_seq)/1000000,2)
print("leaf2GC", leaf2GC)

df = pandas.read_csv(snakemake.input[0], sep="\t", header=0).set_index(["AssemblyNames"])
from metagenlab_libs import ete_phylo 

ete_tree = ete_phylo.EteTool(snakemake.input[1])

for leave in  ete_tree.tree.iter_leaves():
    leave.name = leave.name.split("_genomic")[0]

leaves = [i.name for i in ete_tree.tree.iter_leaves()]

leaf2assembly = {i:i for i in leaves}

leaf2family = df["family"].to_dict()
leaf2genus = df["genus"].to_dict()
leaf2species = df["species"].to_dict()

print(leaf2assembly)


ete_tree.add_text_face(leaf2family, 
                        header_name="family",
                        color_scale=True)


ete_tree.add_text_face(leaf2genus, 
                        header_name="genus",
                        color_scale=True)

ete_tree.add_text_face(leaf2species, 
                        header_name="species",
                        color_scale=False)

ete_tree.add_text_face(leaf2assembly, 
                        header_name="Assembly",
                        color_scale=False)

ete_tree.add_simple_barplot(leaf2GC, 
                   "GC",
                   color=False,
                   show_values=True,
                   substract_min=True,
                   max_value=False)

ete_tree.add_simple_barplot(leaf2genome_size, 
                   "Genome size",
                   color=False,
                   show_values=True,
                   substract_min=True,
                   max_value=False)

print("ok")
ete_tree.remove_dots()
os.environ['QT_QPA_PLATFORM']='offscreen'
print("writing")
ete_tree.tree.render(snakemake.output[0],tree_style=ete_tree.tss, w=183, units="mm")