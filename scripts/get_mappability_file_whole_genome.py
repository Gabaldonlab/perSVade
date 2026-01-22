#!/usr/bin/env python

# calculates the mappability file per chrom and per whole genome. Should be run on perSVade_env

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
ref_genome, threads, final_file, outdir = sys.argv[1:]

# get as full path ome paths
outdir = fun.get_fullpath(outdir)

# log
print("getting mappability into %s."%outdir)
fun.make_folder(outdir)

# define the final file
if not fun.file_is_empty(final_file):
    print("WARNING: already ran, exiting...")
    sys.exit(0)

# process inputs
threads = int(threads)

# make tmpdir
tmpdir = "%s/tmp"%outdir
fun.make_folder(tmpdir)

# link
ref_genome_linked = "%s/genome.fasta"%tmpdir
fun.soft_link_files(ref_genome, ref_genome_linked)

# running mappability
fun.printing_verbose_mode = True
fun.generate_genome_mappability_file(ref_genome_linked, replace=False, threads=threads)

# keep some files
for f in ["genome.fasta.mappability_per_position_per_chromosome", "genome.fasta.mappability_per_position.bed"]:
	os.rename("%s/%s"%(tmpdir, f), "%s/%s"%(outdir, f))

# clean and exit
fun.delete_folder(tmpdir)
open(final_file, "w").write("finished...")
print("finished")


