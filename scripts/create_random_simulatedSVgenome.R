#!/usr/bin/env Rscript

# This script takes a genome and outputs a simulated genome and breakpoint info under ---outdir

# define environment
library(argparser, quietly=TRUE)
library(RSVSim, quietly=TRUE) # doc in https://rdrr.io/bioc/RSVSim/man/simulateSV.html
library(emdbook, quietly=TRUE)
library(rtracklayer, quietly=TRUE)
library(GenomicRanges, quietly=TRUE)
library(R.utils, quietly=TRUE)

# parse cmd line args
argp = arg_parser("Takes a genome and generates a simulated genome with rearrangements with the known rearrangements in outdir. It will generate these rearrangements and the rearranged genome under outdir, only in the regions that are provided (--regions_bed)")

argp = add_argument(argp, "--input_genome", help="Path to the genome where to generate the SV")
argp = add_argument(argp, "--outdir", help="Path to the directory where to write all the files")
argp = add_argument(argp, "--regions_bed", help="Path to the bed file where the simulations should be generated")

argp = add_argument(argp, "--number_Dup", default=100, help="The number of duplications to generate")
argp = add_argument(argp, "--number_Ins", default=100, help="The number of insertions to generate")
argp = add_argument(argp, "--number_Inv", default=100, help="The number of inversions to generate")
argp = add_argument(argp, "--number_Del", default=100, help="The number of deletions to generate")
argp = add_argument(argp, "--number_Tra", default=100, help="The number of translocations to generate")

argp = add_argument(argp, "--percCopiedIns", default=0.5, help="The fraction of INS that are copy-and-paste")
argp = add_argument(argp, "--percBalancedTrans", default=1, help="The fraction of TRA that are balanced")
argp = add_argument(argp, "--max_time_rearrangement", default=30, help="The maximum number of seconds which a rearrangement will take. This is important because sometimes the simulateSV function gets stuck when simulating mtDNA variation")
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

# define length of each sv
number_Dup = argv$number_Dup
number_Ins = argv$number_Ins
number_Inv = argv$number_Inv
number_Del = argv$number_Del
number_Tra = argv$number_Tra

# re-define the number of translocations if it is too high
#n_translocations = argv$number_Tra
#max_trans =  length(chromosomes)-1
#if (n_translocations > max_trans) { n_translocations = max_trans }



# get the genome object, only taking the ID (first space)
genome_obj = readDNAStringSet(input_genome)
names(genome_obj) = lapply(names(genome_obj), function(x) (strsplit(x, "(\t)|( )")[[1]][1]))
chromosomes = names(genome_obj) 

# if there is only one chromosome there cannot be translocations
if (length(chromosomes) == 1) { number_Tra = 0 }


# make the outdir
if (!dir.exists(outdir)){dir.create(outdir)}

# initialize a log of how the genome was generated
report_genome_rearranging = "#This is a report of how the genome was rearranged:\ntype_genome\tn_InsInvDelDup\tn_Tra\tlongest_event\tfraction_shortest_chromosome\n"

# define the length of the shortest chromosome
len_shortest_chr = min(width(genome_obj[chromosomes]))
  
# initialize the rearranged genome
rearranged_genome_generated = FALSE # a tag that defines that the rearranged genome has been created

# define the fraction of longest chromosome that will be considered
all_fraction_shortest_chr_to_consider = c(0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01, 0.005)

# now iterate on the maxiumum length of the events
for (fraction_shortest_chr_to_consider in all_fraction_shortest_chr_to_consider) {
  print("fraction_shortest_chr_to_consider:"); print(fraction_shortest_chr_to_consider)
  
  # define the sizes of these evebts
  print("getting the size of the SVs")
  size = as.integer(len_shortest_chr*fraction_shortest_chr_to_consider)
  size_Del = estimateSVSizes(n=number_Del, minSize=50, maxSize=size, default="deletions", hist=FALSE)
  size_Ins = estimateSVSizes(n=number_Ins, minSize=50, maxSize=size, default="insertions", hist=FALSE)
  size_Dup = estimateSVSizes(n=number_Dup, minSize=50, maxSize=size, default="tandemDuplications", hist=FALSE)
  size_Inv = estimateSVSizes(n=number_Inv, minSize=50, maxSize=size, default="inversions", hist=FALSE)
  
  # try to run simulations, depending on the provided bed file
  print("getting bed"); print(regions_bed)
  gr_regions = import(regions_bed)

  print("running simulation")
  rearranged_genome = withTimeout(try(simulateSV(output=NA, genome=genome_obj, chrs=chromosomes, dels=number_Del, ins=number_Ins, invs=number_Inv, dups=number_Dup, trans=number_Tra, sizeDels=size_Del, sizeIns=size_Ins, sizeInvs=size_Inv, sizeDups=size_Dup, maxDups=3, percCopiedIns=argv$percCopiedIns, percBalancedTrans=argv$percBalancedTrans, bpFlankSize=0, percSNPs=0, indelProb=0, maxIndelSize=0, repeatBias=FALSE, bpSeqSize=100, random=TRUE, verbose=TRUE, regionsDels=gr_regions, regionsIns=gr_regions, regionsInvs=gr_regions, regionsDups=gr_regions, regionsTrans=gr_regions)), timeout=argv$max_time_rearrangement)
    
  # if it passes the simulation then break
  if (class(rearranged_genome)[1]=="DNAStringSet"){ rearranged_genome_generated=TRUE; break; }
}
  
# at the end check that the genome has been created 
if (!rearranged_genome_generated){stop("It has not been possible to generate the rearranged genome. You may try different combinations of fraction_shortest_chr_to_consider")}
  
# save the report
report = sprintf("%s\t%i\t%i\t%i\t%.3f\n", genome_type, n_events, n_translocations, size, fraction_shortest_chr_to_consider)
report_genome_rearranging = paste(report_genome_rearranging, report, sep="")

# join all the dfs 
write.table(metadata(rearranged_genome)$insertions, file=paste(outdir, "insertions.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(metadata(rearranged_genome)$deletions, file=paste(outdir, "deletions.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(metadata(rearranged_genome)$inversions, file=paste(outdir, "inversions.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(metadata(rearranged_genome)$tandemDuplications, file=paste(outdir, "tandemDuplications.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(metadata(rearranged_genome)$translocations, file=paste(outdir, "translocations.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# join genomes amd write
writeXStringSet(rearranged_genome, filepath = paste(outdir, "rearranged_genome.fasta", sep="/"), format = "fasta")

# write the report
writeLines(report_genome_rearranging, con=paste(outdir, "report_generating_genome.tab", sep="/"))
