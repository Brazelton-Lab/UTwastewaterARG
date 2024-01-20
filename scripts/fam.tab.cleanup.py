#! /usr/bin/env python3

import pandas as pd
df = pd.read_csv("fam.tab",sep="\t")
df = df.rename(columns={'#node_id': 'node_id'})
df['Dbxref'] = df['hmm_id']
df = df.fillna('NA')
df.to_csv("fam.tsv", index=False, header=True, sep="\t", )

