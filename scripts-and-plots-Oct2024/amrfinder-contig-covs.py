#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import sys
import argparse

# user provides names of amrfinder table and contig cov table
parser = argparse.ArgumentParser(description='merge a table of AMRfinder results with a table of contig coverages')

parser.add_argument('-a','--amrfinder', help='table of AMRfinder results with a column named "Protein identifier", tab-delimited', required=True)
parser.add_argument('-c','--contigs', help='table of contig coverages with a column named "Contig", tab-delimited', required=True)
args = parser.parse_args()

# amrfinder results
amr = pd.read_csv(args.amrfinder, sep="\t")

amr[['temp_contig','contig_length','ID']] = amr['Protein identifier'].str.split('_',expand=True)
amr['Contig'] = amr['temp_contig'] + '_' + amr['contig_length']
amr = amr.drop(columns=['temp_contig', 'contig_length','ID'])

# contig coverages
cc = pd.read_csv(args.contigs, sep="\t")

# merged
m = pd.merge(amr, cc, on="Contig", how="inner")

outfilename = args.amrfinder.rsplit(".", 1 )[ 0 ] + '.contig_covs.csv'
m.to_csv(outfilename, index=False)
