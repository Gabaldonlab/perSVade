#!/usr/bin/env python

# This script takes an srr and runs prefetch, fastqdump and trimmomatic on them

##### DEFINE ENVIRONMENT #######

# module imports
import argparse, os
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter
import copy as cp
import pickle
import string
import shutil 
from Bio import SeqIO
import random
import sys
from shutil import copyfile

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import sv_functions as fun

#######

##### CMD ARGS #####

description = """
Downloads an srr 
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# general args
parser.add_argument("--srr", dest="srr", required=True, type=str, help="The SRR")
parser.add_argument("--outdir", dest="outdir", required=True, type=str, help="The output dir")
parser.add_argument("--stop_after_prefetch", dest="stop_after_prefetch", action="store_true", default=False, help="Stop after running the prefetch")
parser.add_argument("--threads", dest="threads", default=16, type=int, help="Number of threads, Default: 16")
parser.add_argument("--replace", dest="replace", action="store_true", default=False, help="Replace existing files")

opt = parser.parse_args()

####################

# define the final reads
final_trimmed_reads1 = "%s/%s_trimmed_reads_1.fastq.gz"%(opt.outdir, opt.srr)
final_trimmed_reads2 = "%s/%s_trimmed_reads_2.fastq.gz"%(opt.outdir, opt.srr)

if any([fun.file_is_empty(f) for f in [final_trimmed_reads1, final_trimmed_reads2]]) or opt.replace is True:

    # make outdir
    fun.make_folder(opt.outdir)

    # run prefetch
    print("running prefetch")
    SRRfile = "%s/%s.srr"%(opt.outdir, opt.srr)
    SRRfile = fun.download_srr_with_prefetch(opt.srr, SRRfile, replace=opt.replace)

    # stop after prefetch
    if opt.stop_after_prefetch is True: 
        print("Exiting after prefetch obtention")
        sys.exit(0)

    # get parallel fastqdump
    print("running parallel fastqdump")
    raw_reads1, raw_reads2 = fun.run_parallelFastqDump_on_prefetched_SRRfile(SRRfile, replace=opt.replace, threads=opt.threads)

    # run trimmomatic
    print("running trimmomatic")
    reads1, reads2 = fun.run_trimmomatic(raw_reads1, raw_reads2, replace=opt.replace, threads=opt.threads)

    # define the SRRfile
    SRRfile = "%s/%s.srr"%(opt.outdir, opt.srr)

    # delete all files that are not the trimmed reads
    for f in os.listdir(opt.outdir): 
        file = "%s/%s"%(opt.outdir, f)

        if file not in {reads1, reads2, SRRfile}: 
            print("removing %s"%file)
            fun.remove_file(file)
            fun.delete_folder(file)

    # move the files
    os.rename(reads1, final_trimmed_reads1)
    os.rename(reads2, final_trimmed_reads2)




