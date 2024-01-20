#! /usr/bin/env python

# merge AMRfinder, mobileOG, and VFDB results into one table
# the goal is to merge on contigs, in order to see which contigs encode multiple results

import pandas as pd

AMR = "AMR.final.cats.csv"
mobileOG = "mobileOG.final.csv"
VFDB = "VFDB.final.csv"

df1 = pd.read_csv(AMR)
df2 = pd.read_csv(mobileOG)
df3 = pd.read_csv(VFDB)

m1 = pd.merge(df1,df2,on="Contig",how="outer")
m2 = pd.merge(m1,df3,on="Contig",how="outer")

m2['Sample'].update(m2.pop('Sample_x'))
m2['Sample'].update(m2.pop('Sample_y'))
#m2['Bitscore'].update(m2.pop('Bitscore_x'))
#m2['Bitscore'].update(m2.pop('Bitscore_y'))
#m2['Query_Title'].update(m2.pop('Query_Title_x'))
#m2['Query_Title'].update(m2.pop('Query_Title_y'))

m2.to_csv("AMR.mobileOG.VFDB.merged.csv", index=False, header=True)
