#!/usr/bin/env Rscript

# This script runs CONY on a sorted bam, together with several files that are necessary for a proper running of CONY.

# define environment
library(argparser, quietly=TRUE)

# print the traceback on exit
#options(error=function()traceback(2))

# parse cmd line args
argp = arg_parser("Perfroms CNV calling with CONY")

argp = add_argument(argp, "--reference_genome", help="Path to the reference_genome")
argp = add_argument(argp, "--sorted_bam", help="The sorted bam with a .bai and removed duplicates")
argp = add_argument(argp, "--regions_file", help="A 1-based regions_file to analyze, which is a .tab file with seqname, start, end.")
argp = add_argument(argp, "--libraries_CONY", help="The path to the CONY.R libraries")
argp = add_argument(argp, "--window_size", default=100, help="The window size on which to run CONY")

opt = parse_args(argp)


######## PROCESS INPUTS ########

# activate the libraries
source(opt$libraries_CONY)

# get the target df
regions_df = read.table(opt$regions_file, sep="\t", header=TRUE)

################################

######## CONY PIPELINE ########

print("running WindowInfo")
CONY.TempRegion = WindowInfo(target.df = regions_df, RefFaFileName=opt$reference_genome, WindowSize=opt$window_size)

print(CONY.TempRegion)

adkjhgkadaj

###############################


print(regions_df)

lajdhjdkha


#argp = add_argument(argp, "--number_Ins", default=0, help="The number of insertions to generate")
#argp = add_argument(argp, "--number_Inv", default=0, help="The number of inversions to generate")
#argp = add_argument(argp, "--number_Del", default=0, help="The number of deletions to generate")
#argp = add_argument(argp, "--number_Tra", default=0, help="The number of translocations to generate")
#argp = add_argument(argp, "--number_Dup", default=0, help="The number of duplications to generate")

#argp = add_argument(argp, "--len_shortest_chr", default=1000, help="Len of the shortest chrom")
#argp = add_argument(argp, "--percCopiedIns", default=0.5, help="The fraction of INS that are copy-and-paste")
#argp = add_argument(argp, "--percBalancedTrans", default=1, help="The fraction of TRA that are balanced")
#argp = add_argument(argp, "--replace", flag=TRUE, default=FALSE, help="Replace genomes that a")

  
# these are the real from cmd args:
input_genome = argv$input_genome
outdir = argv$outdir
regions_tra_bed = argv$regions_tra_bed
len_shortest_chr = argv$len_shortest_chr

# get the genome object, only taking the ID (first space)
genome_obj = readDNAStringSet(input_genome)
names(genome_obj) = lapply(names(genome_obj), function(x) (strsplit(x, "(\t)|( )")[[1]][1]))
chromosomes = names(genome_obj) 

# if there is only one chromosome there cannot be translocations
if (length(chromosomes) == 1) { argv$number_Tra = 0 }

# make the outdir
if (!dir.exists(outdir)){dir.create(outdir)}
  
# initialize the rearranged genome
rearranged_genome_generated = FALSE # a tag that defines that the rearranged genome has been created

# define the fraction of longest chromosome that will be considered
all_fraction_shortest_chr_to_consider = c(0.1, 0.07, 0.05, 0.03, 0.02, 0.01, 0.005)
#all_fraction_shortest_chr_to_consider = c(0.005)

# define the fraction of nevents that will be considered
all_fraction_n_events = c(1, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01, 0.05)
#all_fraction_n_events = c(1)

# define the max_time_rearrangement
all_max_time_rearrangement = c(200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 10000)

# first iterate through the number of events. It is more important to get as much events as possible
for (fraction_n_events in all_fraction_n_events) {

  # define the number of each event
  number_Ins = as.integer(argv$number_Ins*fraction_n_events)
  number_Inv = as.integer(argv$number_Inv*fraction_n_events)
  number_Del = as.integer(argv$number_Del*fraction_n_events)
  number_Tra = as.integer(argv$number_Tra*fraction_n_events)
  number_Dup = as.integer(argv$number_Dup*fraction_n_events)

  # now iterate on the maxiumum length of the events
  for (fraction_shortest_chr_to_consider in all_fraction_shortest_chr_to_consider) {
    print("fraction_n_events:"); print(fraction_n_events)
    print("fraction_shortest_chr_to_consider:"); print(fraction_shortest_chr_to_consider)

    # iterate through several time rearrangements
    for (max_time_rearrangement in all_max_time_rearrangement) {
      print("max_time_rearrangement:"); print(max_time_rearrangement)

      # define the sizes of these evebts
      print("getting the size of the SVs")
      size = as.integer(len_shortest_chr*fraction_shortest_chr_to_consider)

      # define length of each sv
      size_Del = estimateSVSizes(n=number_Del, minSize=50, maxSize=size, default="deletions", hist=FALSE)
      size_Ins = estimateSVSizes(n=number_Ins, minSize=50, maxSize=size, default="insertions", hist=FALSE)
      size_Inv = estimateSVSizes(n=number_Inv, minSize=50, maxSize=size, default="inversions", hist=FALSE)
      size_Dup = estimateSVSizes(n=number_Dup, minSize=50, maxSize=size, default="inversions", hist=FALSE)

      # adjust the number of tra according to the number of chromosomes in the 
      max_tra = length(chromosomes)-1
      if (number_Tra>max_tra) {number_Tra = max_tra}

      # new
      rearranged_genome = withTimeout(try(simulateSV(output=NA, genome=genome_obj, chrs=chromosomes, dels=number_Del, ins=number_Ins, invs=number_Inv, trans=number_Tra, dups=number_Dup, sizeDels=size_Del, sizeIns=size_Ins, sizeInvs=size_Inv,  sizeDups=size_Dup, percCopiedIns=argv$percCopiedIns, percBalancedTrans=argv$percBalancedTrans, bpFlankSize=0, percSNPs=0, indelProb=0, maxIndelSize=0, repeatBias=FALSE, bpSeqSize=100, random=TRUE, verbose=TRUE, maxDups=1)), timeout=max_time_rearrangement)

      # if it passes the simulation then break
      if (class(rearranged_genome)[1]=="DNAStringSet"){ rearranged_genome_generated=TRUE; break; }
    }

    # if it passes the simulation then break
    if (class(rearranged_genome)[1]=="DNAStringSet"){ break }

  }

  # if it passes the simulation then break
  if (class(rearranged_genome)[1]=="DNAStringSet"){ break }

}
  
# at the end check that the genome has been created 
if (!rearranged_genome_generated){stop("It has not been possible to generate the rearranged genome. You may try different combinations of fraction_shortest_chr_to_consider")}
  
# save the report
report_genome_rearranging = "#This is a report of how the genome was rearranged:\tnumber_Ins\tnumber_Inv\tnumber_Del\tnumber_Tra\tlongest_event\tfraction_shortest_chromosome\n"

report = sprintf("%i\t%i\t%i\t%i\t%i\t%.3f\n", number_Ins, number_Inv, number_Del, number_Tra, size, fraction_shortest_chr_to_consider)
report_genome_rearranging = paste(report_genome_rearranging, report, sep="")

# write output
if (number_Ins>0){
  write.table(metadata(rearranged_genome)$insertions, file=paste(outdir, "insertions.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}
if (number_Del>0){
  write.table(metadata(rearranged_genome)$deletions, file=paste(outdir, "deletions.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}
if (number_Inv>0){
  write.table(metadata(rearranged_genome)$inversions, file=paste(outdir, "inversions.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}
if (number_Dup>0){
  write.table(metadata(rearranged_genome)$tandemDuplications, file=paste(outdir, "tandemDuplications.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

# only write translocations if any generated
if (number_Tra>0){
  write.table(metadata(rearranged_genome)$translocations, file=paste(outdir, "translocations.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

# join genomes amd write
#writeXStringSet(rearranged_genome, filepath = paste(outdir, "rearranged_genome.fasta", sep="/"), format = "fasta")

# write the report
writeLines(report_genome_rearranging, con=paste(outdir, "report_generating_genome.tab", sep="/"))



