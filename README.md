# UTwastewaterARG
Identification of ARGs and MGEs in wastewater metagenomes

This repository contains a complete pipeline for identification of antibiotic resistance genes (ARGs) and mobile genetic elements (MGEs) in metagenomes. It was carried out on a dedicated server in the Brazelton Lab at the University of Utah. 

The batch script "utwaste-pipeline-master.sh" uses an assembly-based approach, beginning with quality control of the original sequencing reads and ending with the calculation of sequencing coverage of each predicted ARG and MGE in each metagenome. 

## Dependencies
### Various python tools written by the Brazelton lab, documented in the following repositories:
https://github.com/Brazelton-Lab/seq-qc

https://github.com/Brazelton-Lab/seq-annot

https://github.com/Brazelton-Lab/lab_scripts

### Main software dependencies
bbduk https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/

Megahit https://github.com/voutcn/megahit

Prodigal https://github.com/hyattpd/Prodigal

AMRFinderPlus https://github.com/ncbi/amr

Bowtie2 https://github.com/BenLangmead/bowtie2

samtools https://github.com/samtools/samtools

diamond https://github.com/bbuchfink/diamond

mobileOG https://github.com/clb21565/mobileOG-db/tree/main

VFDB http://www.mgc.ac.cn/VFs/download.htm


## Citations
The overall assembly-based approach, sequence coverage calculations, and the seq-qc and seq-annot packages were described in:
Thornton, CT, Tanner, W, VanDerslice, J, Brazelton, WJ (2020) Localized effect of treated wastewater effluent on the resistome of an urban watershed. GigaScience. 9 (11): giaa125, https://doi.org/10.1093/gigascience/giaa125

In addition, each of the tools listed above should be cited if you use them in your own work.
