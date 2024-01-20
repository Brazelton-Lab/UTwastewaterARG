#! /usr/bin/env python

import pandas as pd
import glob


# merge AMR results from multiple samples
# need to use the *.amrfinder.tsv.b6 files that retain the contig names and are ordered by node ids (original .amrfinder.tsv results files are ordered by gene symbol)

dfs = []
for f in glob.glob("amrfinder_results/*.amrfinder.tsv.b6"):
	colnames = ["Protein identifier","node_id","pid","bitscore","nan1","nan2","nan3","nan4","nan5","nan6","nan7","nan8"]
	AMR = pd.read_csv(f, sep="\t", names=colnames)
	AMR = AMR.drop(list(AMR.filter(regex = 'nan')), axis=1)
	AMR = AMR.drop(["pid","bitscore"], axis=1)

	base = f.replace(".amrfinder.tsv.b6","")
	base = base.replace("amrfinder_results/","")
	
	# extract contig name from "Protein identifier" column of AMR
	s = AMR["Protein identifier"].str.split("_", expand=True)
	AMR["Contig"] = base + "_" + s[0] + "_" + s[1]
	AMR['Sample'] = base
	dfs.append(AMR)
final = pd.concat(dfs, ignore_index=True)
final.to_csv("AMR.final.csv", index=False, header=True)
