#! /usr/bin/env python3

import sys
DB = sys.argv[1]
FAM = DB + "/fam.tab"

import pandas as pd
df = pd.read_csv(FAM,sep="\t")
df = df.rename(columns={'#node_id': 'node_id'})
df['Dbxref'] = df['hmm_id']
df = df.fillna('NA')
df.to_csv("fam.tsv", index=False, header=True, sep="\t", )

