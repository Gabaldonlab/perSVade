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
#bedpe_breakpoints = None
final_svtype_to_svfile, final_rearranged_genome = fun.simulate_SVs_in_genome(reference_genome, "mito_C_glabrata_CBS138", CurDir, nvars=100, bedpe_breakpoints=bedpe_breakpoints, replace=False)

