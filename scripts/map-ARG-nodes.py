#! /usr/bin/env python

import pandas as pd
df1 = pd.read_csv("merged.abunds.node.tpm.tsv",sep="\t")
df1 = df1.rename(columns={'ID': 'node_id'})
df2 = pd.read_csv("fam.tsv",sep="\t")
df2 = df2.drop(columns=["gene_symbol","parent_node_id","Dbxref","hmm_id","family_name","hmm_tc1","hmm_tc2","blastrule_complete_ident","blastrule_complete_wp_coverage","blastrule_complete_br_coverage","blastrule_partial_ident","blastrule_partial_wp_coverage","blastrule_partial_br_coverage","reportable"])
df2 = df2.drop_duplicates(["node_id"])
m1 = pd.merge(df1,df2,on="node_id",how="inner")
m1.to_csv("merged.abunds.tpm.nodes.csv", index=False, header=True)
