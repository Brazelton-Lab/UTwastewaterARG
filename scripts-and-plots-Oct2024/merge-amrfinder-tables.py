#!/usr/bin/env python3

import pandas as pd
import os

extension = ".amrfinder.contig_covs.csv"

dfs = [pd.read_csv(f, index_col=0, dtype={1:'string'})
       for f in os.listdir(os.getcwd()) if f.endswith(extension)]

merged = pd.concat(dfs, sort=False)

merged.to_csv("amrfinder-contig-cov-merged.csv")
