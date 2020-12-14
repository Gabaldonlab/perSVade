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

    - sorted_bam: the path to the sorted bam. This is useful to draw coverage. This script will create a <sorted_bam>.coverage_per_region.tab to skip repeating
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
parser.add_argument("--fraction_y_domain_by_gene_browser", dest="fraction_y_domain_by_gene_browser", default=0.3, type=float, help="The fraction of the yaxis taken by the gene browser")
parser.add_argument("--replace", dest="replace", action="store_true", help="Replace existing files")
parser.add_argument("--only_affected_genes", dest="only_affected_genes", action="store_true", help="add only the affected genes in the browser")
parser.add_argument("--vcf_fields_onHover", dest="vcf_fields_onHover", default="all", type=str, help="A comma-sepparated string of the interesting fields of the vcf to show. If you want fields from the 'INFO', set them as 'INFO_<field>'.")
parser.add_argument("-mchr", "--mitochondrial_chromosome", dest="mitochondrial_chromosome", default="mito_C_glabrata_CBS138", type=str, help="The name of the mitochondrial chromosome. This is important if you have mitochondrial proteins for which to annotate the impact of nonsynonymous variants, as the mitochondrial genetic code is different. This should be the same as in the gff. If there is no mitochondria just put 'no_mitochondria'. If there is more than one mitochindrial scaffold, provide them as comma-sepparated IDs.")
parser.add_argument("--gff_annotation_fields", dest="gff_annotation_fields", default="upmost_parent,ANNOTATION_product", type=str, help="A comma-sepparated string of the interesting fields of the gff (it can include fields in the annotation by starting with 'ANNOTATION_') to add to the browser. By default, it will draw the upmost_parent of each feature, which is usually the ID of the corresponding gene.")
parser.add_argument("--interesting_features", dest="interesting_features", default="all", type=str, help="A comma-sepparated string of the interesting features of the gff to draw.")

parser.add_argument("--coverage_range", dest="coverage_range", default="0,2", type=str, help="A comma-sepparated string of the minimum and max coverage.")



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

#################################################


# add the bgcolor as a random palette if not there
if "bgcolor" not in df.keys(): 

    sample_to_color, palette_sample = fun.get_value_to_color(set(df.sampleID), palette="Paired", type_color="hex")
    df["bgcolor"] = df.sampleID.apply(lambda x: sample_to_color[x])

# get the bgcolor as a list
df["bgcolor"] = df.bgcolor.apply(lambda x: x.split(","))

# load the gff df
df_gff = fun.load_gff3_intoDF(opt.gff)

# check that the fields of the gff are as expected
expected_features = {'pseudogene', 'rRNA', 'tRNA', 'repeat_region', 'chromosome', 'CDS', 'mRNA', 'gene', 'exon', 'centromere', 'ncRNA', 'long_terminal_repeat', 'polyA_site', 'region'} # long_terminal_repeat
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


# define the samples_colors_df
sample_group_labels = opt.sample_group_labels.split(",")
samples_colors_df = pd.DataFrame({sample_label : {sampleID : bgcolor[I] for sampleID, bgcolor in dict(df.set_index("sampleID")["bgcolor"]).items()} for I, sample_label in enumerate(sample_group_labels)})[sample_group_labels]

# define the vcf_fields_onHover (by default it will only display certain fields)
if opt.vcf_fields_onHover!="all": opt.vcf_fields_onHover = set(opt.vcf_fields_onHover.split(","))
else: opt.vcf_fields_onHover = {"#CHROM", "POS", "INFO_BREAKEND_overlaps_repeats", "INFO_BREAKEND_real_AF", "INFO_BREAKENDIDs", "INFO_BREAKEND_coordinates", "INFO_BREAKEND_FILTER"}

# define the gff_annotation_fields. These are extra fields of the gff_df that are written
gff_annotation_fields = set(opt.gff_annotation_fields.split(","))

# define the interesting features
if opt.interesting_features=="all": interesting_features = set(df_gff.feature)
else: interesting_features = set(opt.interesting_features.split(","))

# define  coverage_range
min_cov, max_cov = [float(x) for x in opt.coverage_range.split(",")]

# get the browser
filename = "%s/genome_variation_browser.html"%opt.outdir
gfun.get_genome_variation_browser(df, samples_colors_df, target_regions, target_genes, df_gff, filename, data_dir, opt.reference_genome, threads=opt.threads, sample_group_labels=opt.sample_group_labels.split(","), only_affected_genes=opt.only_affected_genes, vcf_fields_onHover=opt.vcf_fields_onHover, replace=opt.replace, mitochondrial_chromosome=opt.mitochondrial_chromosome, gff_annotation_fields=gff_annotation_fields, fraction_y_domain_by_gene_browser=opt.fraction_y_domain_by_gene_browser, interesting_features=interesting_features, min_cov=min_cov, max_cov=max_cov)
print("genome variation browser was written into %s"%filename)





"""
these are the sv_cnv .vcf all INFO tags

{'INFO_QUAL', 'INFO_len_inserted_sequence_min', 'INFO_any_overlaps_repeats', 'INFO_bpIDs', 'INFO_best_FILTER', '#CHROM', 'INFO_length_inexactHomology_mean', 'POS', 'INFO_END', 'INFO_length_event_mean', 'INFO_BREAKEND_overlaps_repeats', 'INFO_SVTYPE', 'INFO_length_microHomology_max', 'INFO_BREAKEND_allele_frequency_SmallEvent', 'INFO_FILTER', 'INFO_allele_frequency_SmallEvent_mean', 'INFO_RELCOVERAGE_NEIGHBOR', 'INFO_length_event_max', 'INFO_allele_frequency', 'INFO_allele_frequency_max', 'INFO_BREAKPOINTIDs', 'INFO_worse_FILTER', 'INFO_real_AF_min', 'INFO_length_inexactHomology_min', 'INFO_QUAL_min', 'INFO_real_AF_max', 'INFO_len_inserted_sequence_max', 'INFO_QUAL_mean', 'INFO_allele_frequency_mean', 'INFO_length_event_min', 'INFO_length_inexactHomology', 'INFO_BREAKEND_QUAL', 'INFO_length_microHomology_mean', 'INFO_RELCOVERAGE', 'INFO_has_poly16GC', 'INFO_BREAKEND_FILTER', 'INFO_REGION_SPEARMANP', 'INFO_length_microHomology', 'ALT', 'INFO_real_AF_mean', 'INFO_length_inexactHomology_max', 'INFO_BREAKEND_length_microHomology', 'INFO_all_FILTERs', 'INFO_allele_frequency_min', 'INFO_BREAKEND_length_inexactHomology', 'INFO_allele_frequency_SmallEvent_max', 'INFO_REGION_PEARSONP', 'ID', 'INFO_BREAKEND_has_poly16GC', 'INFO_REGION_ABS_PEARSONR', 'INFO_BREAKEND_coordinates', 'INFO_REGION_ABS_SPEARMANR', 'INFO_length_microHomology_min', 'INFO_BREAKEND_real_AF', 'INFO_overlaps_repeats', 'INFO_allele_frequency_SmallEvent', 'INFO_allele_frequency_SmallEvent_min', 'INFO_BREAKEND_len_inserted_sequence', 'INFO_BREAKEND_allele_frequency', 'INFO_BREAKEND_length_event', 'INFO_real_AF', 'INFO_len_inserted_sequence_mean', 'INFO_BREAKENDIDs', 'INFO_BPS_TYPE', 'INFO_variantID', 'INFO_BREAKPOINTID', 'INFO_QUAL_max'}

"""





