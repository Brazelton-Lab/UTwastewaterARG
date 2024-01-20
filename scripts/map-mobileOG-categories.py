#! /usr/bin/env python

import pandas as pd
import glob

# first, gather all of the mobileOG result files from individual samples
# I can't use the mobileOG.final.csv file that was generated earlier in pipeline because that one takes only the best hit per predicted protein, which can differ from the best hit chosen by count_features when there are multiple hits with the same score
dfs = []
for f in glob.glob("mobileOG_results/*.mobileOG.Alignment.Out.csv"):
	data = pd.read_csv(f)
	data = data.rename(columns={'Unnamed: 0': 'seqID'})
	data = data.rename(columns={"Sequence Title": "mobileOG_Result"})
	df = data.filter(["mobileOG_Result","mobileOG ID","Gene Name","Best Hit Accession ID","Major mobileOG Category","Minor mobileOG Category","Source Database","Evidence Type"])
	dfs.append(df)

final = pd.concat(dfs, ignore_index=True)
final = final.drop_duplicates(["mobileOG_Result"])

# merged table from count_features and compare_features
df1 = pd.read_csv("merged.mobileOG.abunds.tpm.tsv",sep="\t")
df1 = df1.rename(columns={'ID': 'mobileOG_Result'})

# extract mobileOG ID from full string that contains a bunch of weird characters that will mess up the final merger
s = df1["mobileOG_Result"].str.split("|", expand=True)
df1["mobileOG ID"] = s[0]

m1 = pd.merge(df1,final,on="mobileOG ID",how="inner")
m1.to_csv("merged.abunds.mobileOG.cats.tpm.csv", index=False, header=True)
