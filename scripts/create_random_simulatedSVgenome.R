#!/usr/bin/env Rscript

# This script takes a genome and outputs a simulated genome and breakpoint info under ---outdir. I was using version   - bioconductor-rsvsim=1.24.0 of bioconductor at 15/05/2020. To update it I had to tweak conda, so that it may be worth to come back to this version

# define environment
library(argparser, quietly=TRUE)
library(RSVSim, quietly=TRUE) # doc in https://rdrr.io/bioc/RSVSim/man/simulateSV.html
library(emdbook, quietly=TRUE)
library(rtracklayer, quietly=TRUE)
library(GenomicRanges, quietly=TRUE)
library(R.utils, quietly=TRUE)

# print the traceback on exit
#options(error=function()traceback(2))

# parse cmd line args
argp = arg_parser("Takes a genome and generates a simulated genome with rearrangements with the known rearrangements in outdir. It will generate these rearrangements and the rearranged genome under outdir, only in the regions that are provided (--regions_bed)")

argp = add_argument(argp, "--input_genome", help="Path to the genome where to generate the SV")
argp = add_argument(argp, "--outdir", help="Path to the directory where to write all the files")
argp = add_argument(argp, "--regions_bed", help="Path to the bed file where the simulations should be generated")

argp = add_argument(argp, "--number_Ins", default=0, help="The number of insertions to generate")
argp = add_argument(argp, "--number_Inv", default=0, help="The number of inversions to generate")
argp = add_argument(argp, "--number_Del", default=0, help="The number of deletions to generate")
argp = add_argument(argp, "--number_Tra", default=0, help="The number of translocations to generate")
argp = add_argument(argp, "--number_Dup", default=0, help="The number of duplications to generate")

argp = add_argument(argp, "--len_shortest_chr", default=1000, help="Len of the shortest chrom")
argp = add_argument(argp, "--percCopiedIns", default=0.5, help="The fraction of INS that are copy-and-paste")
argp = add_argument(argp, "--percBalancedTrans", default=1, help="The fraction of TRA that are balanced")
argp = add_argument(argp, "--max_time_rearrangement", default=120, help="The maximum number of seconds which a rearrangement will take. This is important because sometimes the simulateSV function gets stuck when simulating mtDNA variation")
#argp = add_argument(argp, "--replace", flag=TRUE, default=FALSE, help="Replace genomes that a")

argv = parse_args(argp)

# these are to simulate on one example
#input_genome = "/home/mschikora/samba/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_chromosomes.fasta"
#outdir = "/home/mschikora/samba/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_chromosomes_100dupinsinvdel_50tra_0.5percCopiedIns_0.75percBalancedTrans"
#regions_bed = "/home/mschikora/samba/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_chromosomes_repeats_regions_+-1kb.bed"
#regions_bed = "all_regions"

# define the number of translocations
#n_translocations = argv$number_Tra
#max_trans =  length(chromosomes)-1
#if (n_translocations > max_trans) { n_translocations = max_trans }
  
# these are the real from cmd args:
input_genome = argv$input_genome
outdir = argv$outdir
regions_bed = argv$regions_bed
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
    
    # define the sizes of these evebts
    print("getting the size of the SVs")
    size = as.integer(len_shortest_chr*fraction_shortest_chr_to_consider)

    # define length of each sv
    size_Del = estimateSVSizes(n=number_Del, minSize=50, maxSize=size, default="deletions", hist=FALSE)
    size_Ins = estimateSVSizes(n=number_Ins, minSize=50, maxSize=size, default="insertions", hist=FALSE)
    size_Inv = estimateSVSizes(n=number_Inv, minSize=50, maxSize=size, default="inversions", hist=FALSE)
    size_Dup = estimateSVSizes(n=number_Dup, minSize=50, maxSize=size, default="inversions", hist=FALSE)

    # try to run simulations, depending on the provided bed file
    gr_regions = import(regions_bed)

    # adjust the number of tra according to the number of chromosomes in the 
    max_tra = length(unique(seqnames(gr_regions)))-1
    if (number_Tra>max_tra) {number_Tra = max_tra}

    rearranged_genome = withTimeout(try(simulateSV(output=NA, genome=genome_obj, chrs=chromosomes, dels=number_Del, ins=number_Ins, invs=number_Inv, trans=number_Tra, dups=number_Dup, sizeDels=size_Del, sizeIns=size_Ins, sizeInvs=size_Inv,  sizeDups=size_Dup, percCopiedIns=argv$percCopiedIns, percBalancedTrans=argv$percBalancedTrans, bpFlankSize=0, percSNPs=0, indelProb=0, maxIndelSize=0, repeatBias=FALSE, bpSeqSize=100, random=TRUE, verbose=TRUE, regionsDels=gr_regions, regionsIns=gr_regions, regionsInvs=gr_regions, regionsTrans=gr_regions, regionsDups=gr_regions, maxDups=3)), timeout=argv$max_time_rearrangement)

    #rearranged_genome = simulateSV(output=NA, genome=genome_obj, chrs=chromosomes, dels=number_Del, ins=number_Ins, invs=number_Inv, dups=number_Dup, trans=number_Tra, sizeDels=size_Del, sizeIns=size_Ins, sizeInvs=size_Inv, sizeDups=size_Dup, maxDups=2, percCopiedIns=argv$percCopiedIns, percBalancedTrans=argv$percBalancedTrans, bpFlankSize=0, percSNPs=0, indelProb=0, maxIndelSize=0, repeatBias=FALSE, bpSeqSize=100, random=TRUE, verbose=TRUE, regionsDels=gr_regions, regionsIns=gr_regions, regionsInvs=gr_regions, regionsDups=gr_regions, regionsTrans=gr_regions)
      
    # if it passes the simulation then break
    if (class(rearranged_genome)[1]=="DNAStringSet"){ rearranged_genome_generated=TRUE; break; }
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



