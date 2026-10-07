# user must specify input filenames and path to amrfinder database
# provide raw reads with --skip_none, or else provide qc'd reads for the rest of the pipeline

# example for each input option:

# raw reads as input, don't skip anything:
# python utwaste-pipeline-v2-Sep2025.py --skip_none --forward 2105670_S4_L001_R1_001.fastq.gz --reverse 2105670_S4_L001_R2_001.fastq.gz --db ~/miniconda3/envs/amrfinder/share/amrfinderplus/data/latest

# qc'd reads as input, skip fastp and run assembly and everything downstream:
# python utwaste-pipeline-v2-Sep2025.py --skip_fastp --forward 2105670_S4.forward.atrim.nophix.fastp.fq.gz --reverse 2105670_S4.reverse.atrim.nophix.fastp.fq.gz --db ~/miniconda3/envs/amrfinder/share/amrfinderplus/data/latest

# assembly and qc'd reads as input, skip fastp, skip assembly, and run prodigal and everything downstream:
# python utwaste-pipeline-v2-Sep2025.py --skip_assembly --assembly 2105670_S4.contigs.fa --forward 2105670_S4.forward.atrim.nophix.fastp.fq.gz --reverse 2105670_S4.reverse.atrim.nophix.fastp.fq.gz --db ~/miniconda3/envs/amrfinder/share/amrfinderplus/data/latest

# already-annotated assembly, its predicted proteins, and qc'd reads as input, skip fastp, skip assembly, skip prodigal, and run amrfinder and everything downstream:
# python utwaste-pipeline-v2-Sep2025.py --skip_prodigal --assembly 2105670_S4.renamed.fa --faa 2105670_S4.faa --forward 2105670_S4.forward.atrim.nophix.fastp.fq.gz --reverse 2105670_S4.reverse.atrim.nophix.fastp.fq.gz --db ~/miniconda3/envs/amrfinder/share/amrfinderplus/data/latest

# already-annotated assembly, its amrfinder results, and qc'd reads as input, skip everything except calculation of coverages:
# python utwaste-pipeline-v2-Sep2025.py --skip_amrfinder --assembly 2105670_S4.renamed.fa --amr 2105670_S4.amrfinder.tsv --forward 2105670_S4.forward.atrim.nophix.fastp.fq.gz --reverse 2105670_S4.reverse.atrim.nophix.fastp.fq.gz --db ~/miniconda3/envs/amrfinder/share/amrfinderplus/data/latest

# scripts required for the workflow
SCRIPTS_PATH="scripts"

import os
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Pipeline from raw reads to assembly and predicted ARGs. Options for skipping any of the steps.')

# flags for skip options
parser.add_argument('-s0','--skip_none', action='store_true', help="don't skip anything, run the whole pipeline starting from raw reads", required=False)
parser.add_argument('-s1','--skip_fastp', action='store_true', help="skip fastp step because I am providing reads that are already qc'd", required=False)
parser.add_argument('-s2','--skip_assembly', action='store_true', help="skip fastp and assembly steps because I am providing an assembly as input", required=False)
parser.add_argument('-s3','--skip_prodigal', action='store_true', help="skip fastp, assembly, and prodigal steps because I am providing an already-annotated assembly including .gff and .faa files at the specified path", required=False)
parser.add_argument('-s4','--skip_amrfinder', action='store_true', help="skip fastp, assembly, prodigal, and amrfinder steps because I am providing amrfinder results in the $AMR_PATH", required=False)

# arguments for input options
parser.add_argument('-d','--db', help="path to folder containing amrfinder database downloaded with amrfinder -u", required='--skip_amrfinder' not in sys.argv)
parser.add_argument('-f','--forward', help='input fastq file - forward reads', required=True)
parser.add_argument('-r','--reverse', help='input fastq file - reverse reads', required=True)
parser.add_argument('-a','--assembly', help="input assembly FASTA instead of reads", required='--skip_assembly' in sys.argv or "--skip_prodigal" in sys.argv or "--skip_amrfinder" in sys.argv)
parser.add_argument('-g','--faa', help="path to folder containing .gff and .faa files as input", required='--skip_prodigal' in sys.argv)
parser.add_argument('-m','--amr', help="path to folder containing amrfinder results as input", required='--skip_amrfinder' in sys.argv)

args = parser.parse_args()

# check that at least one skip flag was specified
if not (args.skip_none or args.skip_fastp or args.skip_assembly or args.skip_prodigal or args.skip_amrfinder):
    parser.error('No skip option specified; see --help')

# run fastp if requested by user
FASTP = SCRIPTS_PATH + "/fastp.sh"
if args.skip_none: 
	print("running fastp with", args.forward, "and", args.reverse)
	subprocess.check_call([FASTP, args.forward, args.reverse])
	FORWARD = args.forward.replace("_L001_R1_001.fastq.gz", ".forward.atrim.nophix.fastp.fq.gz")
	REVERSE = args.reverse.replace("_L001_R2_001.fastq.gz", ".reverse.atrim.nophix.fastp.fq.gz")

# run megahit if requested by user
if args.skip_fastp: FORWARD = args.forward
if args.skip_fastp: REVERSE: REVERSE = args.reverse
MEGAHIT = SCRIPTS_PATH + "/megahit.sh"
if args.skip_none or args.skip_fastp:
	print("running megahit with", FORWARD, "and", REVERSE)
	subprocess.check_call([MEGAHIT, FORWARD, REVERSE])
	ASS = os.path.basename(FORWARD).replace(".forward.atrim.nophix.fastp.fq.gz", "")
	ASS = ASS + "_megahit/" + ASS + ".contigs.fa"

# run prodigal if requested by user
if args.skip_assembly: ASS = args.assembly
PRODIGAL = SCRIPTS_PATH + "/prodigal.sh"
if args.skip_none or args.skip_fastp or args.skip_assembly:
	print("running prodigal with", ASS)
	subprocess.check_call([PRODIGAL, ASS])
	FAA = os.path.basename(ASS).replace(".fa", ".faa")

# run amrfinder if requested by user
if args.skip_prodigal: FAA = args.faa
AMRFINDER = SCRIPTS_PATH + "/amrfinder.sh"
if args.skip_none or args.skip_fastp or args.skip_assembly or args.skip_prodigal:
	print("running amrfinder with", FAA)
	DB = args.db
	subprocess.check_call([AMRFINDER, FAA, DB])
	AMR = FAA.replace(".faa", ".amrfinder.tsv")

# run coverm and cleanup files
if args.skip_amrfinder: AMR = args.amr
try: ASS = os.path.basename(ASS).replace(".fa", ".renamed.fa")
except: ASS = args.assembly
FORWARD = args.forward
REVERSE = args.reverse
print("running coverm with", AMR, ",", ASS, ",", FORWARD, ", and", REVERSE)
COVERM = SCRIPTS_PATH + "/coverm.sh"
subprocess.check_call([COVERM, AMR, ASS, FORWARD, REVERSE])

print("pipeline completed!")
