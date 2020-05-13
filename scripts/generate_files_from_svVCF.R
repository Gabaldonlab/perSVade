#!/usr/bin/env Rscript

# This script takes the vcf output of gridss, and generates some files for further downstream analysis

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

###########  FUNCTIONS ##########

get_filters = function(bnd_row, bpgr){
  
  # takes a breakpond ID and a ranges object of breakpoints. It returns the merged FILTER of all the partners of bnd_id, in order of appearance
  
  # get the id
  bnd_id = bnd_row["name"]
  
  # get own filter
  own_filter = bpgr[bnd_id]$FILTER
  
  # get filter of the partner
  partner_filter = bpgr[bpgr[bnd_id]$partner]$FILTER
  
  # join and return
  return(paste(own_filter, partner_filter, sep="||"))
  
}

get_ALTSs = function(bnd_row, bpgr){
  
  # takes a breakpond row and a ranges object of breakpoints. It returns the merged ALT of all the partners of bnd_id, in order of appearance
  
  # get the id
  bnd_id = bnd_row["name"]
  
  # get own ALT
  own_alt = bpgr[bnd_id]$ALT
  
  # get alt of the partner
  partner_ALT = bpgr[bpgr[bnd_id]$partner]$ALT
  
  # join and return
  return(paste(own_alt, partner_ALT, sep="||"))
  
}

get_IDs = function(bnd_row, bpgr){
  
  # takes a breakpond row and a ranges object of breakpoints. It returns the merged breakpoint IDs of all the partners of bnd_id, in order of appearance
  
  # get the id
  bnd_id = bnd_row["name"]

  # get ID of the partner
  partner_id = bpgr[bnd_id]$partner
  
  # join and return
  return(paste(bnd_id, partner_id, sep="||"))

}
################################

# read vcf
vcf = readVcf(input_vcf)

# convert vcf to "ranges" objects
bpgr = breakpointRanges(vcf, info_columns=NULL) # contains all breakends that have a partner

###### BEDPE GENERATION ######

# generate bedpe dataframe
bedpe_df = breakpointgr2bedpe(bpgr)

# add the filter tags of all the breakends
bedpe_df$FILTERS = apply(bedpe_df, 1, function(x) get_filters(x, bpgr))

# add the ALT allele of the breakends
bedpe_df$ALTs = apply(bedpe_df, 1, function(x) get_ALTSs(x, bpgr))

# add the ID of both breakends
bedpe_df$IDs = apply(bedpe_df, 1, function(x) get_IDs(x, bpgr))

# change the scientific format for locations (i.e.: avoid 1e+2.5)
bedpe_df$start1 <- gsub(" ", "", format(bedpe_df$start1, scientific = FALSE)) 
bedpe_df$start2 <-gsub(" ", "", format(bedpe_df$start2, scientific = FALSE)) 
bedpe_df$end1 <-gsub(" ", "", format(bedpe_df$end1, scientific = FALSE)) 
bedpe_df$end2 <-gsub(" ", "", format(bedpe_df$end2, scientific = FALSE)) 

# write 
write.table(bedpe_df, file=output_bedpe, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

############################
