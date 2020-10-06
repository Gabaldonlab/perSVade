#!/usr/bin/env python


# This script tests how to generate random simulations of SVs using only RSVsim, with no inHouse code. it can be either random or inserting a bedpe file

##### DEFINE ENVIRONMENT #######

# module imports
import os
import sys
import pandas as pd
from Bio import SeqIO


# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False    
    threads = 4
else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    threads = 48



# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions
import sv_functions as fun

# define dirs
CurDir = "%s/scripts/perSVade/perSVade_repository/testing/Cglabrata_testingSimulations"%ParentDir
outdir_genomes_and_annotations = "%s/scripts/perSVade/perSVade_repository/testing/genomes_and_annotations"%ParentDir
reference_genome = "%s/Candida_glabrata.fasta"%outdir_genomes_and_annotations

############################

bedpe_breakpoints = "%s/M12_WT_Cglabrata_SVs.bedpe"%CurDir
fun.simulate_SVs_in_genome(reference_genome, "mito_C_glabrata_CBS138", CurDir, nvars=10, bedpe_breakpoints=bedpe_breakpoints)

"""
argp = add_argument(argp, "--regions_bed", help="Path to the bed file where the simulations should be generated")

argp = add_argument(argp, "--number_Ins", default=10, help="The number of insertions to generate")
argp = add_argument(argp, "--number_Inv", default=10, help="The number of inversions to generate")
argp = add_argument(argp, "--number_Del", default=10, help="The number of deletions to generate")
argp = add_argument(argp, "--number_Tra", default=10, help="The number of translocations to generate")
argp = add_argument(argp, "--number_Dup", default=10, help="The number of duplications to generate")

argp = add_argument(argp, "--len_shortest_chr", default=100000, help="The number of duplications to generate")
argp = add_argument(argp, "--percCopiedIns", default=0.5, help="The fraction of INS that are copy-and-paste")
argp = add_argument(argp, "--percBalancedTrans", default=1, help="The fraction of TRA that are balanced")
argp = add_argument(argp, "--max_time_rearrangement", default=120, help="The maximum number of seconds which a rearrangement will take. This is important because sometimes the simulateSV function gets stuck when simulating mtDNA variation")
"""