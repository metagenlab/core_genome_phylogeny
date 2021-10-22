#!/usr/bin/env python 

from ete3 import Tree
import pandas

p = "/scratch/hdd/airflow/execution/diag_pipelines/2021_08_30-1638_VMH6M2/taxonomy/GTDB/gtdbtk.bac120.classify.tree"

t = Tree(p, quoted_node_names=True, format=1)
#print(t)
s = 'f__Chlamydiaceae'
s = '1015432923_5248'
node = t.search_nodes(name=s)[0]

assembly_table = pandas.read_csv("/data/databases/GTDB/release202/taxonomy/gtdb_taxonomy.tsv", sep="\t", header=None)
assembly_table.columns = ["assembly_accession", "taxonomy"]
taxonomy = assembly_table["taxonomy"].str.split(";", n = -1, expand = True)
taxonomy.columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]

print(taxonomy.query('family == "f__Enterobacteriaceae"').shape)

while node:
    print(node.name)
    if 'f__' in node.name:
        print (node.name)
        leaves = [leaf.name for leaf in node.iter_leaves()]
        print(len(leaves))
        break
    node = node.up