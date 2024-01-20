#! /usr/bin/env python

import pandas as pd
import glob

# merge the mobileOG result files from multiple samples
dfs = []
for f in glob.glob("*.mobileOG.Alignment.Out.csv"):
	data = pd.read_csv(f)
	data = data.rename(columns={'Unnamed: 0': 'seqID'})
	
	base = f.replace(".mobileOG.Alignment.Out.csv","")
	
	# extract contig name from "Query Title" column of mobileOG
	s = data["Query Title"].str.split("_", expand=True)
	data["Contig"] = base + "_" + s[0] + "_" + s[1]
	data['Sample'] = base

	data = data.drop(["seqID","Pident","Unique_ORF","e-value","ORF_Start_Stop_Strands","Partial Tag","Start Codon","RBS Motif","RBS Spacer","GC Content","Start of Alignment in Query","Start of Alignment in Subject","End of Alignment in Query","End of Alignment in Query.1","Subject Sequence Length","Query Sequence Length","Sense or Antisense Strand","Prodigal ID","Prodigal Designated Contigs","ORF_Start","ORF_End"], axis=1)

	# mobileOG has multiple hits per predicted protein
	# use pandas to take the row with the max bit score
	dg = data.sort_values("Bitscore").groupby("Query Title", as_index=False).first()
	dg = dg.rename(columns={"Sequence Title": "mobileOG_Result"})
	
	dfs.append(dg)

final = pd.concat(dfs, ignore_index=True)
final.to_csv("mobileOG.final.csv", index=False, header=True)
