#!/usr/bin/env Rscript

# This script takes the vcf output of gridss, and generates some files for further downstream analysis, only the simple version, without IDs or anything

# To install packages I had to use:
# conda install bioconductor-structuralvariantannotation 
# conda install -c bioconda bioconductor-variantannotation
# conda install -c r r-stringr

# load libraries
library(rtracklayer)
library(StructuralVariantAnnotation)

# parse cmd arguments
input_vcf = commandArgs(trailingOnly = TRUE)[1]
#input_vcf = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/output.vcf" # this is tot test
output_bedpe = paste(input_vcf, "bedpe", sep=".")

# read vcf
vcf = readVcf(input_vcf)

# convert vcf to "ranges" objects
bpgr = breakpointRanges(vcf, info_columns=NULL) # contains all breakends that have a partner

###### BEDPE GENERATION ######

# generate bedpe dataframe
bedpe_df = breakpointgr2bedpe(bpgr)

# change the scientific format for locations (i.e.: avoid 1e+2.5)
bedpe_df$start1 <- gsub(" ", "", format(bedpe_df$start1, scientific = FALSE)) 
bedpe_df$start2 <-gsub(" ", "", format(bedpe_df$start2, scientific = FALSE)) 
bedpe_df$end1 <-gsub(" ", "", format(bedpe_df$end1, scientific = FALSE)) 
bedpe_df$end2 <-gsub(" ", "", format(bedpe_df$end2, scientific = FALSE)) 

# write 
write.table(bedpe_df, file=output_bedpe, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

############################
