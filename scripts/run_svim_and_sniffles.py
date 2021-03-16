#!/usr/bin/env python

# This is a script that runs SVIM and SNIFFLES in the perSVade_env for a set of reads


# module imports
import sys
import argparse, os
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)
import sv_functions as fun

### CMD LINE ###

description = """
Takes a file with ONT reads and runs SVIM and SINFFLES on them
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
parser.add_argument("--ref", dest="ref", required=True, help="Reference genome. Has to end with .fasta.")
parser.add_argument("--input_reads", dest="input_reads", required=True, help="The input reads")
parser.add_argument("--outdir", dest="outdir", required=True, help="The outdir")
parser.add_argument("--aligner", dest="aligner", default="ngmlr", help="The aligner for SVIM, it can be ngmlr or minimap2")
parser.add_argument("--threads", dest="threads", default=4, type=int, help="Number of threads, Default: 4")
parser.add_argument("--replace", dest="replace", action="store_true", help="Replace existing files")
parser.add_argument("--verbose", dest="verbose", action="store_true", help="verbosity mode")
opt = parser.parse_args()

################

# make outdir
if opt.replace is True: fun.delete_folder(opt.outdir)
fun.make_folder(opt.outdir)

# define the print mode
fun.printing_verbose_mode = opt.verbose

# put the input reads under outdir
reads_dir = "%s/reads"%opt.outdir; fun.make_folder(reads_dir)
dest_input_reads = "%s/%s"%(reads_dir, fun.get_file(opt.input_reads))
fun.soft_link_files(opt.input_reads, dest_input_reads)
opt.input_reads = dest_input_reads

# run fastqdump if the provided reads are an SRR file
if opt.input_reads.endswith(".srr"):
	print("Assuming that the input_reads %s are an SRR file obtained with prefetch"%opt.input_reads)
	opt.input_reads = fun.run_parallelFastqDump_on_prefetched_SRRfile_nanopore(opt.input_reads, replace=opt.replace, threads=opt.threads)

# run porechop
fun.print_if_verbose("Running porechop")
#trimmed_reads = fun.run_porechop(opt.input_reads,  replace=opt.replace, threads=opt.threads)
trimmed_reads = fun.run_porechop_parallel(opt.input_reads, replace=opt.replace, threads=opt.threads, reads_per_chunk=8000)

# run svim (and also get bam)
fun.print_if_verbose("Running SVIM")
outdir_svim = "%s/svim_output"%opt.outdir
sorted_bam = fun.run_svim(trimmed_reads, opt.ref, outdir_svim,  threads=opt.threads, replace=opt.replace, aligner=opt.aligner, is_nanopore=True)

# run SNIFFLES
fun.print_if_verbose("Running SNIFFLES")
outdir_sniffles = "%s/sniffles_output"%opt.outdir
fun.run_sniffles(sorted_bam, outdir_sniffles, opt.replace, opt.threads)

# write the final file
final_file = "%s/ONT_SV_calling_finished.txt"%opt.outdir
open(final_file, "w").write("SVIM and SNIFFLES run finished")
print("SVIM and SNIFFLES finished")
