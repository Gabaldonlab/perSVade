#!/usr/bin/env Rscript

# This script takes as input the vcf output of gridss and it generates a file that has the simple vcf

# load libraries
library(VariantAnnotation, quiet=TRUE)
library(StructuralVariantAnnotation, quiet=TRUE)
library(stringr, quiet=TRUE)

# parse cmd arguments
input_vcf = commandArgs(trailingOnly = TRUE)[1]
#input_vcf = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/output.vcf" # this is tot test
output_vcf = paste(input_vcf, "withSimpleEventType.vcf", sep=".")

#' Simple SV type classifier
simpleEventType <- function(gr) {
  return(ifelse(seqnames(gr) != seqnames(partner(gr)), "CTX", # inter-chromosomosal
                ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
                       ifelse(strand(gr) == strand(partner(gr)), "INV",
                              ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
                                     "DUP")))))
}

# load the vcf
vcf <- readVcf(input_vcf, "version vcf")

# adds some info 
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
  row.names=c("SIMPLE_TYPE"),
  Number=c("1"),
  Type=c("String"),
  Description=c("Simple event type annotation based purely on breakend position and orientation."))), "DataFrame"))

gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen
writeVcf(vcf, output_vcf) 
