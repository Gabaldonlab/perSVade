#!/usr/bin/env Rscript

# This script runs hmmcopy on an input table were there is the coverage per bins of the genome

# define environment
library(argparser, quietly=TRUE)
library(HMMcopy)

# print the traceback on exit
#options(error=function()traceback(2))
options(warn=1)

# parse cmd line args
argp = arg_parser("Perfroms CNV calling with HMMcopy for one chromosome")

argp = add_argument(argp, "--coverage_table", help="A table with the input files")
argp = add_argument(argp, "--outfile", help="The outfile where to write the outdirs")

# parameters
argp = add_argument(argp, "--e", default=0.9999999, help="The e parameters of the HMMcopy")
argp = add_argument(argp, "--mu", default="-0.458558247,-0.215877601,-0.002665686,0.191051578,0.347816046,1.664333241", help="The mu parameters of the HMMcopy")
argp = add_argument(argp, "--lambda", default=20, help="The lambda parameters of the HMMcopy")
argp = add_argument(argp, "--nu", default=2.1, help="The nu parameters of the HMMcopy")
argp = add_argument(argp, "--kappa", default="50,50,700,100,50,50", help="The kappa parameters of the HMMcopy")
argp = add_argument(argp, "--m", default="-0.458558247,-0.215877601,-0.002665686,0.191051578,0.347816046,1.664333241", help="The m parameters of the HMMcopy")
argp = add_argument(argp, "--eta", default=50000.0, help="The eta parameters of the HMMcopy")
argp = add_argument(argp, "--gamma", default=3, help="The gamma parameters of the HMMcopy")
argp = add_argument(argp, "--S", default=0.02930164, help="The S parameters of the HMMcopy")
argp = add_argument(argp, "--strength", default="1e7", help="The strength parameters of the HMMcopy")

opt = parse_args(argp)

# load df
df_coverage = read.table(opt$coverage_table, sep="\t", header=TRUE)

# add columns that have corrected fields
df_coverage = correctReadcount(df_coverage)

# get the default parameters to run HMMsegment
params_HMMsegment = HMMsegment(df_coverage, getparam = TRUE)

# define the parameters according to the inputs
params_HMMsegment$strength = as.numeric(opt$strength)
params_HMMsegment$e = opt$e
params_HMMsegment$mu = as.numeric(strsplit(opt$mu, ",")[[1]])
params_HMMsegment$lambda = opt$lambda
params_HMMsegment$nu = opt$nu
params_HMMsegment$kappa = as.numeric(strsplit(opt$kappa, ",")[[1]])
params_HMMsegment$m = as.numeric(strsplit(opt$m, ",")[[1]])
params_HMMsegment$eta = opt$eta
params_HMMsegment$gamma = opt$gamma
params_HMMsegment$S = opt$S

# classify into copy number variation even
CN_segments = HMMsegment(df_coverage, params_HMMsegment)

# add the segments
df_coverage$state_CNsegment = CN_segments$state # this goes from 1 to 6

# write
write.table(df_coverage, opt$outfile, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)

print("HMM copy worked well")

