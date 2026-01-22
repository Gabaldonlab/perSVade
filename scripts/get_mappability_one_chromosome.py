#!/usr/bin/env python

# Gets mappability one chrom

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
ref_genome, outdir, chromosome, threads, idx_folder = sys.argv[1:]
outdir = fun.get_fullpath(outdir)
ref_genome = fun.get_fullpath(ref_genome)
idx_folder = fun.get_fullpath(idx_folder)
threads = int(threads)

# log
print("getting mappability into %s for chrom %s"%(outdir, chromosome))
fun.make_folder(outdir)

# make tmp
tmpdir = outdir+"/tmp"
fun.make_folder(tmpdir)

# define the final file
final_file = outdir+"/finished_file.txt"
if not fun.file_is_empty(final_file):
    print("WARNING: already ran, exiting...")
    sys.exit(0)

# link genome
ref_genome_linked = "%s/genome.fasta"%tmpdir
fun.soft_link_files(ref_genome, ref_genome_linked)

# redefine perSVade tmp dir
os.environ["PERSVADE_TMPDIR"] = tmpdir+"/running_cmds"
fun.make_folder(os.environ["PERSVADE_TMPDIR"])

# get a df with windows of the genome of 0.1Mb, to run per window, keeping only windows of the target chrom
fun.index_genome(ref_genome_linked)
df_windows = fun.get_df_windows(ref_genome_linked, 1000000)
df_windows["chromosome"] = df_windows.chromosome.apply(str)
df_windows = df_windows.sort_values(by=["chromosome", "start", "end"])
df_windows_chrom = df_windows[df_windows.chromosome==chromosome].copy()

# checks
if len(df_windows_chrom)==0: raise ValueError("invalid")

# get the mappability per windows file for whole chrom
bed_region = "%s/bed_region_%s.bed"%(tmpdir, chromosome)
pd.DataFrame({"chromosome":[chromosome], "start":[df_windows_chrom.start.iloc[0]], "end":[df_windows_chrom.end.iloc[-1]]}).to_csv(bed_region, sep="\t", header=False, index=False)

# get the mappability per windows file
map_outfile_chrom = "%s/map_per_window_%s.bed"%(tmpdir, chromosome)
genmap_std_file = "%s.map_running.std"%map_outfile_chrom

if fun.file_is_empty(map_outfile_chrom):
    map_outfile_chrom_tmp = "%s.tmp.bed"%map_outfile_chrom
    fun.remove_file(map_outfile_chrom_tmp)
    fun.run_cmd("%s map -E 2 -K 30 -I %s -O %s -b --threads %i -v --selection %s > %s 2>&1"%(fun.genmap, idx_folder, ".bed".join(map_outfile_chrom_tmp.split(".bed")[0:-1]), threads, bed_region, genmap_std_file))
    os.rename(map_outfile_chrom_tmp, map_outfile_chrom)

# add regions with 0.0 mappability
print("getting bed with 0 mappability")
map_outfile_chrom_0_map = "%s.0map.bed"%map_outfile_chrom
if fun.file_is_empty(map_outfile_chrom_0_map):
    map_outfile_chrom_0_map_tmp = "%s.tmp"%map_outfile_chrom_0_map
    fun.run_cmd("""%s subtract -a %s -b %s | awk -v OFS='\\t' '{print $1, $2, $3, "-", "0.0"}' > %s"""%(fun.bedtools, bed_region, map_outfile_chrom, map_outfile_chrom_0_map_tmp))
    os.rename(map_outfile_chrom_0_map_tmp, map_outfile_chrom_0_map)

# generate the merged file of both
print("get merged bed")
map_outfile_chrom_merged = "%s.merged.bed"%map_outfile_chrom
if fun.file_is_empty(map_outfile_chrom_merged):
    map_outfile_chrom_merged_tmp = "%s.tmp"%map_outfile_chrom_merged
    fun.run_cmd("cat %s %s | sort -k1,1 -k2,2n -k3,3n > %s"%(map_outfile_chrom_0_map, map_outfile_chrom, map_outfile_chrom_merged_tmp))
    os.rename(map_outfile_chrom_merged_tmp, map_outfile_chrom_merged)

# split file in chunks
print("getting chunks...")
chunk_size = 10000
tmpdir_chrom = "%s/chrom_%s_csize_%i"%(tmpdir, chromosome, chunk_size)
fun.make_folder(tmpdir_chrom)

chunks_dir = "%s/chunks_windows"%tmpdir_chrom
if not os.path.isdir(chunks_dir):

    # log chunks
    nlines_merged_bed = fun.get_nlines_file_no_writing(map_outfile_chrom_merged)
    nchunks = len(list(fun.chunks(list(range(int(nlines_merged_bed))), chunk_size)))

    # split reads
    chunks_dir_tmp = "%s_tmp"%chunks_dir
    fun.delete_folder(chunks_dir_tmp); fun.make_folder(chunks_dir_tmp)
    for suffix_len in range(3, 100):
        if (24**suffix_len) > nchunks: break
    fun.run_cmd("split -a %i -l %i %s %s/chunk."%(suffix_len, chunk_size, map_outfile_chrom_merged, chunks_dir_tmp))

    # keep
    for f in os.listdir(chunks_dir_tmp): fun.check_correct_filename_dots(f, 2)
    os.rename(chunks_dir_tmp, chunks_dir)

# for each chunk, get the mappability per position in parallel
print("getting map per pos in parallel...")
chunks_dir_map = "%s/mappability_per_pos"%tmpdir_chrom; fun.make_folder(chunks_dir_map)
nchunks = len(os.listdir(chunks_dir))
inputs_fn = []
fields_mappability = ["chromosome", "position", "unique_map_score", "is_uniquely_mappable"]
for Ic, chunk_name in enumerate(sorted(os.listdir(chunks_dir))):
    mappability_file_chunk = "%s/%s_map_per_pos.bed"%(chunks_dir_map, chunk_name)
    inputs_fn.append((Ic, nchunks, fields_mappability, mappability_file_chunk, "%s/%s"%(chunks_dir, chunk_name), chromosome))

with multiproc.Pool(threads) as pool:
    pool.starmap(fun.get_mappability_per_pos_one_chunk, inputs_fn, chunksize=1) 
    pool.close()
    pool.terminate()

# create map_outfile_long_chrom
print("appending files...")
map_outfile_long_chrom = "%s/mappability_chrom_%s.bed"%(outdir, chromosome)
map_outfile_long_chrom_tmp = "%s.tmp"%map_outfile_long_chrom
fun.remove_file(map_outfile_long_chrom_tmp)
open(map_outfile_long_chrom_tmp, "w").write("\t".join(fields_mappability)+"\n")
fun.run_cmd("cat %s/* | sort -k2,2n >> %s"%(chunks_dir_map, map_outfile_long_chrom_tmp))

# checks
print("checking file...")
expected_chrom_len = df_windows_chrom.end.iloc[-1] - df_windows_chrom.start.iloc[0]
nlines_chrom = fun.get_int_output_bash_command("wc -l %s"%map_outfile_long_chrom_tmp)
if expected_chrom_len!=(nlines_chrom-1): raise ValueError("# lines are not correct")
if fun.get_int_output_bash_command("tail -n +2 %s | head -1 | cut -f2"%map_outfile_long_chrom_tmp) != 0: raise ValueError("pos1 should be 0")
if fun.get_int_output_bash_command("tail -n1 %s | cut -f2"%map_outfile_long_chrom_tmp) != expected_chrom_len-1: raise ValueError("end should be chrom len")
if fun.get_int_output_bash_command("cut -f2 %s | sort | uniq | wc -l"%map_outfile_long_chrom_tmp) != (expected_chrom_len+1): raise ValueError("there should be one pos per line")

# clean
os.rename(map_outfile_long_chrom_tmp, map_outfile_long_chrom)
fun.delete_folder(tmpdir)
open(final_file, "w").write("finished...")
print("finished")

