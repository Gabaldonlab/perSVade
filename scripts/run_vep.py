#!/usr/bin/env python

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

# define binaries that are in the EnvDir
bgzip = "%s/bin/bgzip"%EnvDir
tabix = "%s/bin/tabix"%EnvDir
bedtools = "%s/bin/bedtools"%EnvDir
vep = "%s/bin/vep"%EnvDir

######################################################
######################################################
######################################################


### ARGUMENTS ###

parser = argparse.ArgumentParser(description="Takes a vcf and writes the output of VEP. It also replaces the annotation of the mitochondrial genes, as indicated by those genes that are in mitochondrial_chromosome and also the ones in the nuclear genome, as indicated by gDNA_code")

# general args
parser.add_argument("--input_vcf", dest="input_vcf", required=True, type=str, help="The input vcf file.")
parser.add_argument("--outfile", dest="outfile", required=True, type=str, help="The outfile")
parser.add_argument("--ref", dest="ref", required=True, type=str, help="The reference genome")
parser.add_argument("--gff", dest="gff", required=True, type=str, help="The gff file")
parser.add_argument("--mitochondrial_chromosome", dest="mitochondrial_chromosome", required=True, type=str, help="The name of mitochondrial chromosomes,  which can be comma sepparated if there are more than 1.")
parser.add_argument("--mito_code", dest="mito_code", required=False, default=3, type=int, help="The code of translation of ncbi of mitochondrial proteins. Fungal mitochondrial by default.")
parser.add_argument("--gDNA_code", dest="gDNA_code", required=False, default=1, type=int, help="The code of translation of ncbi of nuclear genes. Standard by default. C. albicans has 12")
parser.add_argument("--log_file_all_cmds", dest="log_file_all_cmds", default=None, help="An existing log_file_all_cmds to store the cmds")

opt = parser.parse_args()

# print the command line to run this

# add the cmds file
fun.log_file_all_cmds = opt.log_file_all_cmds
if fun.log_file_all_cmds is not None and fun.file_is_empty(fun.log_file_all_cmds): raise ValueError("The provided --log_file_all_cmds %s should exist"%fun.log_file_all_cmds)

#############################################
################ RUNNING VEP ################
#############################################

# get the gff tabixed and sorted
gff_clean = "%s_clean.gff"%opt.gff
gff_clean_compressed = "%s_clean.gz"%opt.gff
gff_clean_compressed_tbi = "%s.tbi"%gff_clean_compressed

if fun.file_is_empty(gff_clean_compressed_tbi):

    # remove previous files
    fun.remove_file(gff_clean)
    fun.remove_file(gff_clean_compressed)
    fun.remove_file(gff_clean_compressed_tbi)

    print("compressing gff before running vep")

    # eliminate strange lines,chromosomes and compress
    fun.run_cmd("%s sort -i %s | egrep -v '^#' | egrep -v $'\tchromosome\t' > %s"%(bedtools, opt.gff, gff_clean))
    fun.run_cmd("%s -c %s > %s"%(bgzip, gff_clean, gff_clean_compressed))

    # index with tabix
    fun.run_cmd("%s %s"%(tabix, gff_clean_compressed))


###### test that the gff is correct ######

# load df
gff_fields = ["chromosome", "source", "type_feature", "start", "end", "score", "strand", "phase", "attributes"]
df_gff = pd.read_csv(gff_clean, header=None, sep="\t", names=gff_fields)

# check that all type_feature are correct
print("checking that the specifications are correct according to https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gfftypes")
allowed_fields = {"aberrant_processed_transcript", "CDS", "C_gene_segment", "D_gene_segment", "exon", "gene", "J_gene_segment", "lincRNA", "lincRNA_gene", "miRNA", "miRNA_gene", "mRNA", "mt_gene", "ncRNA", "NMD_transcript_variant", "primary_transcript", "processed_pseudogene", "processed_transcript", "pseudogene", "pseudogenic_transcript", "RNA", "rRNA", "rRNA_gene", "snoRNA", "snoRNA_gene", "snRNA", "snRNA_gene", "supercontig", "transcript", "tRNA", "VD_gene_segment", "V_gene_segment"}

not_allowed_features = set(df_gff[~df_gff.type_feature.isin(allowed_fields)].type_feature)
if len(not_allowed_features)>0:  print("WARNING: %s are features not parsed by VEP in the 3d col of the gff"%not_allowed_features)

# check that all mRNAs have exons
print("If something goes wrong you may want to validate that your GFF is ok (gt gff3validator gff)")
n_exons = sum(df_gff.type_feature=="exon")
n_mRNAs = sum(df_gff.type_feature=="mRNA")
if n_exons<n_mRNAs: raise ValueError("Each mRNA must have at least one exon in the provided gff")

##########################################

# define the outfile of vep raw
outfile_vep_raw = "%s.raw.tbl"%opt.outfile
outfile_vep_raw_tmp = "%s.tmp"%outfile_vep_raw

if fun.file_is_empty(outfile_vep_raw):
    fun.remove_file(outfile_vep_raw_tmp)

    print("running vep. If this fails, you may double check that the provided gff is consistent with https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gfftypes")

    # test that there are variants in the input
    nlines_vep_input = len([l for l in open(opt.input_vcf, "r").readlines() if not l.startswith("#")])
    if nlines_vep_input==0: raise ValueError("The input of vep is empty")

    # define the backbone_cmd
    vep_std = "%s.std"%outfile_vep_raw_tmp

    cmd = '%s --input_file %s --format "vcf" --output_file %s --fasta %s -v --force_overwrite --tab --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,Extra" --gff %s > %s 2>&1'%(vep, opt.input_vcf, outfile_vep_raw_tmp, opt.ref, gff_clean_compressed, vep_std) # --plugin Blosum62
    fun.run_cmd(cmd)

    # check that there are no errors in the output
    if any([any({e in l.upper() for e in {"EXCEPTION", "ERROR"}}) and (not l.upper().startswith("WARNING")) for l in open(vep_std, "r").readlines()]): raise ValueError("There was an error running vep. Check %s to see what happened. You may want to validate that your GFF is ok (in http://genometools.org/cgi-bin/gff3validator.cgi)"%vep_std) 

    # print all the warnings found in vep_std
    for l in open(vep_std, "r").readlines(): 
        if "WARNING" in l.upper(): print("VEP "+l)

    # check that <10% of the variants were not annotated
    nlines_vep_output = len([l for l in open(outfile_vep_raw_tmp, "r").readlines() if not l.startswith("#")])

    if (nlines_vep_output/nlines_vep_input) < 0.9: raise ValueError("There is less than 90 perecent of annotated vars")

    fun.remove_file(vep_std)

    # rename
    os.rename(outfile_vep_raw_tmp, outfile_vep_raw)


#############################################
#############################################
#############################################

#############################################
######## CORRECTING THE GENETIC CODE ########
#############################################

print("Correcting proteins for the genetic_code")

# get into df
vep_df = pd.read_csv(outfile_vep_raw, sep="\t", header=len([x for x in open(outfile_vep_raw, "r") if x.startswith("##")]), na_values=fun.vcf_strings_as_NaNs, keep_default_na=False)

# print a warning if all consequences are intergenic
if all(vep_df.Consequence=="intergenic_variant"): fun.print_with_runtime("!!!! WARNING !!!!! All the variants annotated by VEP are intergenic for a chunk of variants (those from '%s'). This suggests that the input files may not be properly formatted."%opt.input_vcf)

# modify the consequences for BND rows
vep_df["Consequence"] = vep_df.apply(fun.get_corrected_Consequence_for_vep_r, axis=1)

# define the expected vars
all_expected_consequences = {'stop_gained', 'intron_variant', 'upstream_gene_variant', '5_prime_UTR_variant', 'inframe_insertion', 'synonymous_variant', 'non_coding_transcript_exon_variant', 'intergenic_variant', 'protein_altering_variant', 'coding_sequence_variant', 'downstream_gene_variant', '3_prime_UTR_variant', 'missense_variant', 'splice_region_variant', 'splice_acceptor_variant', 'inframe_deletion', 'stop_lost', 'non_coding_transcript_variant', 'start_retained_variant', 'frameshift_variant', 'stop_retained_variant', 'start_lost', 'incomplete_terminal_codon_variant', 'splice_donor_variant', 'transcript_amplification', 'transcript_ablation', 'feature_elongation', 'feature_truncation'}

all_expected_consequences = all_expected_consequences.union({"%s_BND"%cons for cons in all_expected_consequences})

# check that all the found consequences are in the expected ones
all_found_consequences = set.union(*vep_df.Consequence.apply(lambda x: set(x.split(","))))
not_expected_consequences = all_found_consequences.difference(all_expected_consequences)

if len(not_expected_consequences)>0: raise ValueError("%s are not expected consequences. Maybe you are using a different VEP version"%not_expected_consequences)

#print("\n\n", vep_df[vep_df.Consequence.apply(lambda x: "incomplete_terminal_codon_variant" in x)], "\n\n")

# define the mito chromosomes as a set
mitochondrial_chromosomes_set = set(opt.mitochondrial_chromosome.split(","))

# define the types of variants that are affected by the genetic code
genCode_affected_vars = {'stop_retained_variant', 'inframe_deletion', 'inframe_insertion', 'frameshift_variant', 'synonymous_variant', 'missense_variant', 'stop_gained', 'stop_lost', 'protein_altering_variant'}

# define the idxs of each type of genes
typeGenes_to_idx = {"mito": vep_df.apply(lambda row: row["Location"].split(":")[0] in mitochondrial_chromosomes_set and len(set(row["Consequence"].split(",")).intersection(genCode_affected_vars))>0, axis=1),

                    "nuclear": vep_df.apply(lambda row: row["Location"].split(":")[0] not in mitochondrial_chromosomes_set and len(set(row["Consequence"].split(",")).intersection(genCode_affected_vars))>0, axis=1)}


# define the code of each type of genes
typeGenes_to_code = {"mito":opt.mito_code, "nuclear":opt.gDNA_code}

# define a dataframe of non-affected rows in the VEP df
all_df = vep_df[~(typeGenes_to_idx["mito"]) & ~(typeGenes_to_idx["nuclear"])]

# go through each type of genes
for typeGenes, idx_affected_rows in typeGenes_to_idx.items():

    # define the affected df
    affected_df = vep_df[idx_affected_rows]    

    # define the genCode
    genCode = typeGenes_to_code[typeGenes]

    # define the stop codons
    stop_codons = set([codon for codon in ["".join(s) for s in itertools.product(["A", "C", "T", "G"], repeat=3)] if str(Seq.Seq(codon).translate(table = genCode))=="*"])

    # modify the rows if there's something to modify
    if len(affected_df)>0: 
        affected_df[["Amino_acids","Consequence"]] = affected_df.apply(lambda r: fun.modify_DF_cols(r, genCode, stop_codons, genCode_affected_vars), axis=1)

    # keep
    all_df = all_df.append(affected_df)


# add the BLOSUM62 score for missense variants
all_df["BLOSUM62_score"] = all_df.apply(fun.get_BLOSUM62_score_missense_variants, axis=1)

#### CHECK THAT ALL THE VARIANTS HAVE BEEN ANNOTATED ####

# check that the relationship between the VEP Uploaded_var and merged_vcf_all is 1:1
uploaded_variation = set(all_df["#Uploaded_variation"])
all_variants = set(fun.get_df_and_header_from_vcf(opt.input_vcf)[0]["ID"])

if len(uploaded_variation.difference(all_variants))>0: 
    raise ValueError("There are some uploaded variations that can't be found in all_variants. This could be because there are multiallelic variants, which are not handled properly by this module. These are the unexpected vars: %s"%(uploaded_variation.difference(all_variants)))

# deinfe the unnanotated vars as those that are not in the VEP output and are also not missing 
missing_vars = all_variants.difference(uploaded_variation)
unnanotated_vars = {v for v in missing_vars if v.split("/")[-1]!="*"}

if len(unnanotated_vars)>0: 
    print("WARNING: There are some variants that have not been annotated with VEP:\n%s\n (%i/%i in total)"%("\n".join(unnanotated_vars), len(unnanotated_vars), len(all_variants)))

# raise ValueError if there are more than 10% of unannotated vars
n_all_vars = len(all_variants)
n_unnanotated_vars = len(unnanotated_vars)
#if n_unnanotated_vars/n_all_vars>0.1: raise ValueError("There are more than 10 perecent of unannotated variants by VEP")
if (n_unnanotated_vars/n_all_vars)>0.1: print("WARNING: There are more than 10 pct of unannotated variants by VEP in this chunk (%s). This could be due to errors or because there is no annotations in this chromosome."%opt.input_vcf)

########################################################

print("%i/%i vars were annotated, and all of them have IDs matching the input vcf"%(n_all_vars-n_unnanotated_vars, n_all_vars))

# write to the same as outfile
all_df.to_csv(opt.outfile, sep="\t", index=False, header=True)



#############################################
#############################################
#############################################

