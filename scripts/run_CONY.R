#!/usr/bin/env Rscript

# This script runs CONY on a sorted bam, together with several files that are necessary for a proper running of CONY. This is parallelised

# define environment
library(argparser, quietly=TRUE)

# print the traceback on exit
#options(error=function()traceback(2))
options(warn=1)

# parse cmd line args
argp = arg_parser("Perfroms CNV calling with CONY for one chromosome")

argp = add_argument(argp, "--chromosome", help="The name of the chromosome")
argp = add_argument(argp, "--sample_name", help="The name of the sample in the sorted_bam")
argp = add_argument(argp, "--libraries_CONY", help="The path to the CONY.R libraries")
argp = add_argument(argp, "--outdir", help="Set the current working dir")
argp = add_argument(argp, "--ploidy", help="The ploidy")
argp = add_argument(argp, "--fragment_len", help="The fragment_len related to CONY paralelisation")

opt = parse_args(argp)


######## PROCESS INPUTS ########

# activate the libraries
source(opt$libraries_CONY)

# set the outdir
setwd(opt$outdir)

################################

######## CONY PIPELINE ########

# get UsedRD, which is useful to remove the windows with  CN=0 or CN=1
print("running UsedRD")
UsedRD(CRDMethod="SumUp", AnaMethod="Single", TargetChr=opt$chromosome, SampleName=opt$sample_name)

# exit if the usedRD did not yield any regions to use
usedRD_output = paste0("./CONY.3-TempRegion.", opt$chromosome, ".", opt$sample_name, ".SumUp.Single.UsedRD.txt", sep="")
usedRD_df = read.table(usedRD_output, sep=" ", header=TRUE)
if (length(rownames(usedRD_df))==0) {

  print("There are no used regions. exiting")
  quit(status=0)

}

# estimate the parameters used to define each of the copies. NCN=5 means that there will be ranges from 1-to-5 CN
print("running EstPar")

ploidy = as.integer(opt$ploidy)
if (ploidy==1){ncopies = 4}
if (ploidy==2){ncopies = 4}

EstPar(CRDMethod="SumUp", AnaMethod="Single", TargetChr=opt$chromosome, SampleName=opt$sample_name, NCN=ncopies)

# run the pipeline
print("running CONY")

# define the fragment length related to the chromosome length (inferred from )
fragment_len = as.integer(opt$fragment_len)
RunCONY(CRDMethod="SumUp", AnaMethod="Single", TargetChr=opt$chromosome, SampleName=opt$sample_name, RunTime = 300000, BurnN = 5000, RTN = 1000, BCPoint = 20, FragLength=fragment_len)



# get the result
print("running the integration of results")
ComResult(CRDMethod="SumUp", AnaMethod="Single", TargetChr=opt$chromosome, SampleName=opt$sample_name)

###############################

print("CONY worked properly")
