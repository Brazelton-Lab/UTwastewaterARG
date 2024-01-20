#! /usr/bin/env python

import pandas as pd
import glob

# merge the VFDB result files from multiple samples
dfs = []
for f in glob.glob("*.VFDB.tsv"):
	colnames=["VFDB_Result","Query Title","Pident","Bitscore","Subject Sequence Length","e-value","Query Sequence Length","Start of Alignment in Subject","End of Alignment in Query","Start of Alignment in Query","End of Alignment in Query.1"]
	data = pd.read_csv(f, sep="\t", names=colnames)

	base = f.replace(".VFDB.tsv","")
	
	# extract contig name from "Query Title" column of VFDB
	s = data["Query Title"].str.split("_", expand=True)
	data["Contig"] = base + "_" + s[0] + "_" + s[1]
	data['Sample'] = base
	
	data = data.drop(["Pident","e-value","Start of Alignment in Query","Subject Sequence Length","Query Sequence Length","Start of Alignment in Subject","Start of Alignment in Subject","End of Alignment in Query","End of Alignment in Query.1"], axis=1)
	
	# VFDB has multiple hits per predicted protein
	# use pandas to take the row with the max bit score
	dg = data.sort_values("Bitscore").groupby("Query Title", as_index=False).first()

	dfs.append(dg)	

final = pd.concat(dfs, ignore_index=True)
final.to_csv("VFDB.final.csv", index=False, header=True)
