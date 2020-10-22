#!/usr/bin/env python

# This script is to plot a genome browser from several perSVade results

######################################################
###############  DEFINE ENVIRONMENT ##################
######################################################

import sys

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import sv_functions as fun
import graphics_functions as gfun

# This is a script to run vep
import argparse, os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import Seq 
import Bio.Seq as Seq
import Bio.SeqRecord as SeqRecord
import copy as cp
import itertools
from argparse import RawTextHelpFormatter

######################################################
######################################################
######################################################

### ARGUMENTS ###

description = """
Takes a table with the paths to several genome variation analyses. It generates an .html file with the genome browser. These are the accepted fields in --input_data table, which determine which type of data will be drawn:

    - sorted_bam: the path to the sorted bam. This is useful to draw coverage.
    - sampleID (mandatory)
    - SV_CNV_vcf: a vcf containing SV and CNV info from perSVade.
    - SV_CNV_var_annotation. The annotation with VEP of SV_CNV_vcf
    - smallVars_vcf: a vcf containing the small variants.
    - smallVars_var_annotation: The annotation with VEP of smallVars_vcf
    - bgcolor: The background color of each sample (as a color or HEX). It can be a set of ',' sepparated colors, in which case several colors will be drawn as columns.
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# mandatory args
parser.add_argument("--input_data", dest="input_data", required=True, type=str, help="The input .tab file, which should contain one line per sample with the required information")
parser.add_argument("--outdir", dest="outdir", required=True, type=str, help="The outfile")
parser.add_argument("--reference_genome", dest="reference_genome", required=True, type=str, help="The reference genome")
parser.add_argument("--gff", dest="gff", required=True, type=str, help="The gff file")

# target regions
parser.add_argument("--target_regions", dest="target_regions", default=None, type=str, help="A bed file with target regions. Genes overlapping these target regions will be drawn. It should have 'chromosome', 'start' and 'end' . If none of --target_regions or --target_genes are specified, all variants will be drawn.")
parser.add_argument("--target_genes", dest="target_genes", default=None, type=str, help="A string of comma sepparated genes ro draw. Only these genes and variants related to them will be shown. If not indicated, all genes will be shown. If none of --target_regions or --target_genes are specified, all variants will be drawn.")
parser.add_argument("--sample_group_labels", dest="sample_group_labels", default="sample", type=str, help="A comma-sepparated string of the labels related to the different cathegories of 'bgcolor' from inout_data.")

parser.add_argument("-thr", "--threads", dest="threads", default=16, type=int, help="Number of threads, Default: 16")
parser.add_argument("--replace", dest="replace", action="store_true", help="Replace existing files")
parser.add_argument("--only_affected_genes", dest="only_affected_genes", action="store_true", help="add only the affected genes in the browser")
parser.add_argument("--vcf_fields_onHover", dest="vcf_fields_onHover", default="all", type=str, help="A comma-sepparated string of the interesting fields of the vcf to show. If you want fields from the 'INFO', set them as 'INFO_<field>'.")



opt = parser.parse_args()

##################

print("getting visualization for %s into %s"%(opt.input_data, opt.outdir))

# get the visualization dataframe
df = pd.read_csv(opt.input_data, sep="\t")

# sebug non-uniqueIDs
if len(set(df.sampleID))!=len(df): raise ValueError("sampleID should be unique in --input_data")

######### GET ALL THE DATA UNDER OUTDIR #########

# define the cahce files dir
fun.make_folder(opt.outdir)
data_dir = "%s/cache_files"%opt.outdir
if opt.replace is True: fun.delete_folder(data_dir)
fun.make_folder(data_dir)

# gff
new_gff = "%s/annotations.gff"%data_dir
fun.soft_link_files(opt.gff, new_gff)
opt.gff = new_gff

# ref
new_reference_genome = "%s/reference_genome.fasta"%data_dir
fun.soft_link_files(opt.reference_genome, new_reference_genome)
opt.reference_genome = new_reference_genome

# get the sorted bams
if "sorted_bam" in df.keys():

    sample_ID_to_newSortedBam = {}
    for I, r in df.iterrows():
        new_sorted_bam = "%s/%s_aligned_reads_sorted.bam"%(data_dir, r["sampleID"])
        new_sorted_bai = "%s/%s_aligned_reads_sorted.bam.bai"%(data_dir, r["sampleID"])

        fun.soft_link_files(r["sorted_bam"], new_sorted_bam)
        fun.soft_link_files("%s.bai"%r["sorted_bam"], new_sorted_bai)

        # keep
        sample_ID_to_newSortedBam[r["sampleID"]] = new_sorted_bam

    df["sorted_bam"] = df.sampleID.apply(lambda x: sample_ID_to_newSortedBam[x])

#################################################


# add the bgcolor as a random palette if not there
if "bgcolor" not in df.keys(): 

    sample_to_color, palette_sample = gfun.get_value_to_color(set(df.sampleID), palette="hls", type_color="hex")
    df["bgcolor"] = df.sampleID.apply(lambda x: sample_to_color[x])

# get the bgcolor as a list
df["bgcolor"] = df.bgcolor.apply(lambda x: x.split(","))

# load the gff df
df_gff = gfun.load_gff3_intoDF(opt.gff)

# check that the fields of the gff are as expected
expected_features = {'pseudogene', 'rRNA', 'tRNA', 'repeat_region', 'chromosome', 'CDS', 'mRNA', 'gene', 'exon', 'centromere', 'ncRNA', 'long_terminal_repeat'} # long_terminal_repeat
missing_features = set(df_gff.feature).difference(expected_features)
if len(missing_features)>0: raise ValueError("Features %s are not expected in the gff"%missing_features)


# get the gene info into df

###### DEFINE THE TARGET REGIONS AND GENES ######

# define all genes
all_genes = set(df_gff[df_gff.feature.isin({"pseudogene", "gene"})].upmost_parent)

# define the target genes. Define the target regions
if opt.target_genes is not None: target_genes = set(opt.target_genes.split(","))
else: target_genes = set()

# define the target regions
if opt.target_regions is not None: target_regions = pd.read_csv(opt.target_regions, sep="\t")
else: df_regions = pd.DataFrame(columns=["chromosome", "start", "end"])

# if both are empty, just get them all
if opt.target_regions is None and opt.target_genes is None:

    target_genes = all_genes
    target_regions = pd.DataFrame({I : {"chromosome":chrom, "start":0, "end":length} for I, (chrom, length) in enumerate(fun.get_chr_to_len(opt.reference_genome).items())}).transpose()

# check that the fields are correct
target_regions = target_regions[["chromosome", "start", "end"]]

# check that all the genes are in the gff
missing_genes = target_genes.difference(all_genes)
if len(missing_genes)>0: raise ValueError("%s are missing genes"%missing_genes)

######################################################

######## if there is only one sample, add another as blank ########
"""
if len(df)==1:

    data_dict = {k : df[k].iloc[0] for k in df.keys()}
    data_dict["sampleID"] = "blank_to_skip"
    data_dict["bgcolor"] = ["red"]*len(data_dict["bgcolor"])
    df = df.append(pd.DataFrame({0:data_dict}).transpose())

"""

###################################################################


# define the samples_colors_df
sample_group_labels = opt.sample_group_labels.split(",")
samples_colors_df = pd.DataFrame({sample_label : {sampleID : bgcolor[I] for sampleID, bgcolor in dict(df.set_index("sampleID")["bgcolor"]).items()} for I, sample_label in enumerate(sample_group_labels)})[sample_group_labels]

# define the vcf_fields_onHover
if opt.vcf_fields_onHover!="all": opt.vcf_fields_onHover = set(opt.vcf_fields_onHover.split(","))


# get the browser
filename = "%s/genome_variation_browser.html"%opt.outdir
gfun.get_genome_variation_browser(df, samples_colors_df, target_regions, target_genes, df_gff, filename, data_dir, opt.reference_genome, threads=opt.threads, sample_group_labels=opt.sample_group_labels.split(","), only_affected_genes=opt.only_affected_genes, vcf_fields_onHover=opt.vcf_fields_onHover)
print("genome variation browser was written into %s"%filename)











