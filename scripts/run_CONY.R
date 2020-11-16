#!/usr/bin/env Rscript

# This script runs CONY on a sorted bam, together with several files that are necessary for a proper running of CONY.

# define environment
library(argparser, quietly=TRUE)

# print the traceback on exit
#options(error=function()traceback(2))
options(warn=1)

# parse cmd line args
argp = arg_parser("Perfroms CNV calling with CONY for one chromosome")

argp = add_argument(argp, "--reference_genome", help="Path to the reference_genome")
argp = add_argument(argp, "--sorted_bam", help="The sorted bam with a .bai and removed duplicates")
argp = add_argument(argp, "--regions_file", help="A 1-based regions_file to analyze, which is a .tab file with seqname, start, end.")
argp = add_argument(argp, "--mpileup_file", help="The mpileup_file, which contains position, 1-based, and read depth")
argp = add_argument(argp, "--chromosome", help="The name of the chromosome")
argp = add_argument(argp, "--sample_name", help="The name of the sample in the sorted_bam")
argp = add_argument(argp, "--libraries_CONY", help="The path to the CONY.R libraries")
argp = add_argument(argp, "--outdir", help="Set the current working dir")
argp = add_argument(argp, "--window_size", default=100, help="The window size on which to run CONY")

opt = parse_args(argp)


######## PROCESS INPUTS ########

# activate the libraries
source(opt$libraries_CONY)

# get the target df
regions_df = read.table(opt$regions_file, sep="\t", header=TRUE)

# set the outdir
setwd(opt$outdir)

################################

######## CONY PIPELINE ########

# run WindowInfo, which calculates the information of each window
print("running WindowInfo")
CONY.TempRegion = as.data.frame(WindowInfo(target.df=regions_df, RefFaFileName=opt$reference_genome, WindowSize=opt$window_size))
colnames(CONY.TempRegion) = c("seqname", "start", "end", "width", "nonAmb", "GC")

# change the types of variables 
CONY.TempRegion$start = as.integer(as.character(CONY.TempRegion$start))
CONY.TempRegion$end = as.integer(as.character(CONY.TempRegion$end))
CONY.TempRegion$width = as.integer(as.character(CONY.TempRegion$width))
CONY.TempRegion$nonAmb = as.numeric(as.character(CONY.TempRegion$nonAmb))
CONY.TempRegion$GC = as.numeric(as.character(CONY.TempRegion$GC))

#outputs a window information file with 6 columns. It includes the name of chromosome (seqname), start position (start), end position (end), width length (width), percentage of indefinable base (nonAmb), and GC percentage (GC) for

# CalRD: calculate the coverage per windows
print("running CalRD")
carRD_df = CalRD(TempRegion=CONY.TempRegion, CRDMethod="SumUp", SampleBamFileName=opt$sorted_bam, MPileCountFileName=opt$mpileup_file,SampleName=opt$sample_name, TargetChr=opt$chromosome, WindowSize=opt$window_size)

# adjusting the coverage per GC content and nonAmb
print("running AdjRD")
AdjRD(CRDMethod= "SumUp", TargetChr=opt$chromosome, SampleName=opt$sample_name)

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
EstPar(CRDMethod="SumUp", AnaMethod="Single", TargetChr=opt$chromosome, SampleName=opt$sample_name, NCN=4)

# run the pipeline
print("running CONY")

# define the fragment length related to the chromosome length (inferred from )
fragment_len = min(c(as.integer(regions_df$end * 0.05) + 1, 500000))

RunCONY(CRDMethod="SumUp", AnaMethod="Single", TargetChr=opt$chromosome, SampleName=opt$sample_name, RunTime = 300000, BurnN = 5000, RTN = 1000, BCPoint = 20, FragLength=fragment_len)



# get the result
print("running the integration of results")
ComResult(CRDMethod="SumUp", AnaMethod="Single", TargetChr=opt$chromosome, SampleName=opt$sample_name)

###############################

print("CONY worked properly")
