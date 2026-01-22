#!/usr/bin/env python

# Calculates genmap idx for genome, necessary for mappability calculations. Should be run on perSVade_env

# imports
import os, sys, re, time
import pandas as pd
import copy as cp
import multiprocessing as multiproc
import shutil
import numpy as np
import itertools

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# import functions
import sv_functions as fun

# parse args
ref_genome, outdir = sys.argv[1:]

# get as full path ome paths
outdir = fun.get_fullpath(outdir)

# log
print("getting mappability into %s."%outdir)
fun.make_folder(outdir)

# define the final file
final_file = outdir+"/finished_file.txt"
if not fun.file_is_empty(final_file):
    print("WARNING: already ran, exiting...")
    sys.exit(0)

# link
ref_genome_linked = "%s/genome.fasta"%outdir
fun.soft_link_files(ref_genome, ref_genome_linked)

# first create the index
genome_dir = "/".join(ref_genome_linked.split("/")[0:-1])
genome_name = ref_genome_linked.split("/")[-1]
idx_folder = "%s/%s_genmap_idx"%(genome_dir, genome_name.split(".")[0])
genmap_std_file = "%s.genmap.std"%ref_genome_linked

# get the genome index files, for the whole genome
print("Getting the index files...")
expected_idx_files = ["index.info.concat", "index.lf.drv.sbl", "index.sa.val", "index.txt.limits", "index.lf.drv", "index.info.limits", "index.rev.lf.drp"]
if any([fun.file_is_empty("%s/%s"%(idx_folder, x)) for x in expected_idx_files]):

    # remove previously generated indices
    if os.path.isdir(idx_folder): shutil.rmtree(idx_folder)
    fun.run_cmd("%s index -F %s -I %s -v > %s 2>&1"%(fun.genmap, ref_genome_linked, idx_folder, genmap_std_file))

# clean and exit
for f in [genmap_std_file, ref_genome_linked]: fun.remove_file(f)
open(final_file, "w").write("finished...")
print("finished")


