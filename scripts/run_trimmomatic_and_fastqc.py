#!/usr/bin/env python

# This script takes some reads and runs fastqc and trimmomatic for them


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

description = """
Runs fastqc and trimmomatic on an a set of paired WGS samples.
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)


# general args
parser.add_argument("-f1", "--fastq1", dest="fastq1", default=None, help="fastq_1 file. Option required to obtain bam files. It can be 'auto', in which case a set of 10M reads will be generated.")
parser.add_argument("-f2", "--fastq2", dest="fastq2", default=None, help="fastq_2 file. Option required to obtain bam files. It can be 'auto', in which case a set of 10M reads will be generated.")
parser.add_argument("-thr", "--threads", dest="threads", default=16, type=int, help="Number of threads, Default: 16")
parser.add_argument("--replace", dest="replace", action="store_true", default=False, help="Replace existing files")

opt = parser.parse_args()

fun.run_trimmomatic(opt.fastq1, opt.fastq2, replace=opt.replace, threads=opt.threas)

