#!/usr/bin/env Rscript

# This script runs aneufinder on an input table were there is the coverage per bins of the genome

# define environment
library(argparser, quietly=TRUE)
library(AneuFinder, quietly=TRUE)
library(GenomicRanges, quietly=TRUE)

########## ARGS ##########

options(warn=1)

argp = arg_parser("Perfroms CNV calling with ANEUPLOIDY for a set of chromosomes.")

# mandatory
argp = add_argument(argp, "--coverage_table", help="A table with the input files. This should have 'seqnames', 'start', 'end', 'counts' ")
argp = add_argument(argp, "--outfile", help="The outfile where to write the outfiles.")
argp = add_argument(argp, "--threads", help="The number of threads.")

# optional
argp = add_argument(argp, "--R", default=10, help="The maximum number of random permutations to use in eachiteration of the permutation test (see e.divisive). Increase this value to increase accuracy on the cost of speed")
argp = add_argument(argp, "--sig_lvl", default=0.1, help="The level at which to sequentially test if a proposed change point is statistically significant (see e.divisive). Increase this value to find more breakpoints.")

opt = parse_args(argp)

##########################

# load the coverage df as a GRs
df_coverage = read.table(opt$coverage_table, sep="\t", header=TRUE)
gr_coverage = GRanges(df_coverage)

# define parameters
states = c("zero-inflation", paste0(0:4, "-somy")) # define CN between 0 and 4
threads = as.integer(opt$threads)

# define optional parms
R = as.integer(opt$R)
sig_lvl = as.numeric(opt$sig_lvl)

# print
print("running findCNVs with the following parms")
print(R)
print(sig_lvl)
print(threads)

# run
cnv_calls = findCNVs(gr_coverage, most.frequent.state="2-somy", states=states, ID="AneuFinderID", strand="*", num.threads=threads, verbosity=4, method="edivisive", R=R, sig.lvl=sig_lvl) 

# save
cnv_calls = data.frame(cnv_calls$segments)
write.table(cnv_calls, opt$outfile, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)

# remove all R objects
rm(list = ls())

