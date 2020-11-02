#!/usr/bin/env python

# This file is for running trimmomatic. It will not repeat anything unless stated

import argparse
import os

parser = argparse.ArgumentParser(description="Runs trimmomatic for paired end sequencing given the full path to the raw reads and the trimmed reads. Woon't repeat unless stated")

# MANDATORY ARGS
parser.add_argument("-rr1", "--raw_reads1", dest="raw_reads1", required=True, type=str, help="The path to the raw reads 1")
parser.add_argument("-rr2", "--raw_reads2", dest="raw_reads2", required=True, type=str, help="The path to the raw reads 2")
parser.add_argument("-tr1", "--trimmed_reads1", dest="trimmed_reads1", required=True, type=str, help="The path to the trimmed reads1")
parser.add_argument("-tr2", "--trimmed_reads2", dest="trimmed_reads2", required=True, type=str, help="The path to the trimmed reads2")
parser.add_argument("-ad", "--adapters_filename", dest="adapters_filename", required=True, type=str, help="The path to the fasta with adapters")

# OPTIONAL
parser.add_argument("--repeat", dest="repeat", action="store_true", default=False, help="Replace existing files")
parser.add_argument("--number_threads", dest="number_threads", default=2,  type=int, help="number threads")

# clipping options
parser.add_argument("--seed_mismatches", dest="seed_mismatches", default=2,  type=int, help="ILLUMINACLIP seed mismatches")
parser.add_argument("--pal_clip_tshd", dest="pal_clip_tshd", default=30,  type=int, help="ILLUMINACLIP palindrome clip threshold")
parser.add_argument("--simple_clip_tshd", dest="simple_clip_tshd", default=10,  type=int, help="ILLUMINACLIP simple clip threshold")
parser.add_argument("--leading", dest="leading", default=15,  type=int, help="LEADING value")
parser.add_argument("--trailing", dest="trailing", default=15,  type=int, help="TRAILING value")
parser.add_argument("--window_size", dest="window_size", default=4,  type=int, help="window size")
parser.add_argument("--window_req_quality", dest="window_req_quality", default=15,  type=int, help="window required quality")
parser.add_argument("--min_len", dest="min_len", default=36,  type=int, help="min length")

opt = parser.parse_args()

##### FUNCTIONS #######

def run_cmd(cmd):
    out_stat = os.system(cmd); 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd, out_stat))

def file_is_empty(path): 
    
    """ask if a file is empty or does not exist """
    
    if not os.path.isfile(path):
        return_val = True
    elif os.stat(path).st_size==0:
        return_val = True
    else:
        return_val = False
            
    return return_val

#######################

import sys

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# import functions
import sv_functions as fun

# paths 
JAVA = "%s/bin/java"%EnvDir
TRIMMOMATIC = "%s/share/trimmomatic/trimmomatic.jar"%EnvDir 

# parse CMD
raw_reads1 = opt.raw_reads1; raw_reads2 = opt.raw_reads2; trimmed_reads1 = opt.trimmed_reads1; trimmed_reads2 = opt.trimmed_reads2; 

# define the unpaired reads
unpaired_trimmed_reads1 = "%s_unpaired.fq.gz"%(trimmed_reads1.rstrip(".fq.gz"))
unpaired_trimmed_reads2 = "%s_unpaired.fq.gz"%(trimmed_reads2.rstrip(".fq.gz"))

# get the tmp files of each outputfile
trimmed_reads1_tmp = "%s.tmp.fq.gz"%(trimmed_reads1.rstrip(".fq.gz")); trimmed_reads2_tmp = "%s.tmp.fq.gz"%(trimmed_reads2.rstrip(".fq.gz"))
unpaired_trimmed_reads1_tmp = "%s.tmp.fq.gz"%(unpaired_trimmed_reads1.rstrip(".fq.gz")) ; unpaired_trimmed_reads2_tmp = "%s.tmp.fq.gz"%(unpaired_trimmed_reads2.rstrip(".fq.gz"));

# remove any previously existing temporary files
for rdir in {trimmed_reads1_tmp, unpaired_trimmed_reads1_tmp, trimmed_reads2_tmp, unpaired_trimmed_reads2_tmp}:
    if os.path.isfile(rdir): os.unlink(rdir)

# run trimmomatic and generate temporary files
trim_cmd = "%s -jar %s PE -phred33 -threads %i %s %s %s %s %s %s ILLUMINACLIP:%s:%i:%i:%i LEADING:%i TRAILING:%i SLIDINGWINDOW:%i:%i MINLEN:%i TOPHRED33"%(JAVA, TRIMMOMATIC, opt.number_threads, raw_reads1, raw_reads2, trimmed_reads1_tmp, unpaired_trimmed_reads1_tmp, trimmed_reads2_tmp, unpaired_trimmed_reads2_tmp, opt.adapters_filename, opt.seed_mismatches, opt.pal_clip_tshd, opt.simple_clip_tshd, opt.leading, opt.trailing,  opt.window_size, opt.window_req_quality, opt.min_len); run_cmd(trim_cmd)

# check that the reads are correct
fun.check_that_paired_reads_are_correct(trimmed_reads1_tmp, trimmed_reads2_tmp)

# rename the temporary files so that you don't repeat.
os.rename(trimmed_reads1_tmp, trimmed_reads1)
os.rename(trimmed_reads2_tmp, trimmed_reads2)
os.rename(unpaired_trimmed_reads1_tmp, unpaired_trimmed_reads1)
os.rename(unpaired_trimmed_reads2_tmp, unpaired_trimmed_reads2)

