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
argp = arg_parser("Takes a genome and outputs a simulated genome and breakpoint info under ---outdir. This script generates SV that are up to 10% of the shortest chromosome (for gDNA) or 5% of the chromosome (for mtDNA). For mtDNA, only 20% of the --number_DupInsInvDel will be generated. The distribution of lengths of the SV will be simulated from a beta distribution, which has been observed in the human genome (there are less large variants than short). This script writes under outdir the rearranged genome and the tables with the different rearrangements. ")

argp = add_argument(argp, "--input_genome", help="Path to the genome where to generate the SV")
argp = add_argument(argp, "--outdir", help="Path to the directory where to write all the files")
argp = add_argument(argp, "--regions_bed", default="all_regions", help="Target regions in bed format. If not provided, the program will be run in all regions ")
argp = add_argument(argp, "--number_DupInsInvDel", default=50, help="The number of non-imbalanced events to generate. Ins are regions that are cut and pasted (and/or) repeated in the genome. The cmd percCopiedIns specifies a float that indicates which fraction of these are copy-and-paste.")
argp = add_argument(argp, "--percCopiedIns", default=0.5, help="The fraction of INS that are copy-and-paste")
argp = add_argument(argp, "--number_Tra", default=50, help="The number of translocations to generate. There are fewer by default")
argp = add_argument(argp, "--percBalancedTrans", default=0.75, help="The fraction of TRA that are balanced")
argp = add_argument(argp, "--mitochondrial_chromosome", default="mito_C_glabrata_CBS138", help="The mito chromosome. This script will obviate it when deciding the length of the largest event (which will be 0.5 of the size of the smallest nuclear chromosome)")
argp = add_argument(argp, "--max_time_rearrangement", default=60, help="The maximum number of seconds which a rearrangement will take. This is important because sometimes the simulateSV function gets stuck when simulating mtDNA variation")
#argp = add_argument(argp, "--replace", flag=TRUE, default=FALSE, help="Replace genomes that a")

argv = parse_args(argp)

# these are to simulate on one example
#input_genome = "/home/mschikora/samba/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_chromosomes.fasta"
#outdir = "/home/mschikora/samba/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_chromosomes_100dupinsinvdel_50tra_0.5percCopiedIns_0.75percBalancedTrans"
#regions_bed = "/home/mschikora/samba/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_chromosomes_repeats_regions_+-1kb.bed"
#regions_bed = "all_regions"

# these are the real from cmd args:
input_genome = argv$input_genome
outdir = argv$outdir
regions_bed = argv$regions_bed

# get the genome object, only taking the ID (first space)
genome_obj = readDNAStringSet(input_genome)
names(genome_obj) = lapply(names(genome_obj), function(x) (strsplit(x, "(\t)|( )")[[1]][1]))

# get the chromosomes of each type
gDNA_chromosomes = grep(argv$mitochondrial_chromosome, names(genome_obj), value=TRUE, invert=TRUE)
mtDNA_chromosomes =  grep(argv$mitochondrial_chromosome, names(genome_obj), value=TRUE, invert=FALSE)

# make the outdir
if (!dir.exists(outdir)){dir.create(outdir)}

# initialize dataframes with the metadata
df_ins = data.frame()
df_del = data.frame()
df_inv = data.frame()
df_dup = data.frame()
df_tra = data.frame()

# initialize a log of how the genome was generated
report_genome_rearranging = "#This is a report of how the genome was rearranged:\ntype_genome\tn_InsInvDelDup\tn_Tra\tlongest_event\tfraction_shortest_chromosome\n"

# go through each type of genome and generate rearrangements
for (genome_type in c("gDNA", "mtDNA")) {
  
  # get a list with the chromosomes of each genome
  if (genome_type=="gDNA"){
    chromosomes = gDNA_chromosomes
  }else {
    chromosomes = mtDNA_chromosomes
  }
  
  # define the length of the shortest chromosome
  len_shortest_chr = min(width(genome_obj[chromosomes]))
  
  # define the number of translocations
  n_translocations = argv$number_Tra
  max_trans =  length(chromosomes)-1
  if (n_translocations > max_trans) { n_translocations = max_trans }
  
  # initialize the rearranged genome
  rearranged_genome_generated = FALSE # a tag that defines that the rearranged genome has been created
  
  # try several combinations of the number of events (not translocations) and the size,until you get one that works. You will start with the biggest possible and keep decreasing
  if (genome_type=="gDNA"){
    all_fraction_n_events = c(1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01, 0.05)
  }else {
    all_fraction_n_events = c(0.1, 0.05, 0.01) # for mtDNA simulate less rearrangements
  }
  
  all_fraction_shortest_chr_to_consider = c(0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01, 0.005)

  # first iterate through the number of events. It is more important to get as much events as possible
  for (fraction_n_events in all_fraction_n_events) {
    
    # now iterate on the maxiumum length of the events
    for (fraction_shortest_chr_to_consider in all_fraction_shortest_chr_to_consider) {
      print(genome_type); print("fraction_n_events:"); print(fraction_n_events); print("fraction_shortest_chr_to_consider:"); print(fraction_shortest_chr_to_consider)
      
      # define the number of events
      n_events = as.integer(argv$number_DupInsInvDel*fraction_n_events)
      
      # define the sizes of these evebts
      size = as.integer(len_shortest_chr*fraction_shortest_chr_to_consider)
      size_Del = estimateSVSizes(n=n_events, minSize=50, maxSize=size, default="deletions", hist=FALSE)
      size_Ins = estimateSVSizes(n=n_events, minSize=50, maxSize=size, default="insertions", hist=FALSE)
      size_Dup = estimateSVSizes(n=n_events, minSize=50, maxSize=size, default="tandemDuplications", hist=FALSE)
      size_Inv = estimateSVSizes(n=n_events, minSize=50, maxSize=size, default="inversions", hist=FALSE)
      
      # try to run simulations, depending on the provided bed file
      if (regions_bed=="all_regions"){
        
        rearranged_genome = withTimeout(try(simulateSV(output=NA, genome=genome_obj, chrs=chromosomes, dels=n_events, ins=n_events, invs=n_events, dups=n_events, trans=n_translocations, sizeDels=size_Del, sizeIns=size_Ins, sizeInvs=size_Inv, sizeDups=size_Dup, maxDups=5, percCopiedIns=argv$percCopiedIns, percBalancedTrans=argv$percBalancedTrans, bpFlankSize=0, percSNPs=0, indelProb=0, maxIndelSize=0, repeatBias=FALSE, bpSeqSize=100, random=TRUE, verbose=TRUE)), timeout=argv$max_time_rearrangement)
        
      }else{
        
        gr_regions = import(regions_bed)
        rearranged_genome = withTimeout(try(simulateSV(output=NA, genome=genome_obj, chrs=chromosomes, dels=n_events, ins=n_events, invs=n_events, dups=n_events, trans=n_translocations, sizeDels=size_Del, sizeIns=size_Ins, sizeInvs=size_Inv, sizeDups=size_Dup, maxDups=5, percCopiedIns=argv$percCopiedIns, percBalancedTrans=argv$percBalancedTrans, bpFlankSize=0, percSNPs=0, indelProb=0, maxIndelSize=0, repeatBias=FALSE, bpSeqSize=100, random=TRUE, verbose=TRUE, regionsDels=gr_regions, regionsIns=gr_regions, regionsInvs=gr_regions, regionsDups=gr_regions, regionsTrans=gr_regions)), timeout=argv$max_time_rearrangement)
        
      }
      
      # if it passes the simulation then break
      if (class(rearranged_genome)[1]=="DNAStringSet"){ rearranged_genome_generated=TRUE; break; }
    }
    
    # if it passes the simulation then break
    if (class(rearranged_genome)[1]=="DNAStringSet"){ break }
  }
  
  # at the end check that the genome has been created 
  if (!rearranged_genome_generated){stop("It has not been possible to generate the rearranged genome. You may try different combinations of fraction_n_events and fraction_shortest_chr_to_consider")}
  
  # save the report
  report = sprintf("%s\t%i\t%i\t%i\t%.3f\n", genome_type, n_events, n_translocations, size, fraction_shortest_chr_to_consider)
  report_genome_rearranging = paste(report_genome_rearranging, report, sep="")

  # save genomes
  if (genome_type=="gDNA"){
    gDNA_rearranged_genome = rearranged_genome
  }else {
    mtDNA_rearranged_genome = rearranged_genome
  }
  
  # keep the events
  df_ins = rbind(df_ins, metadata(rearranged_genome)$insertions)
  df_del = rbind(df_del, metadata(rearranged_genome)$deletions)
  df_inv = rbind(df_inv, metadata(rearranged_genome)$inversions)
  df_dup = rbind(df_dup, metadata(rearranged_genome)$tandemDuplications)
  df_tra = rbind(df_tra, metadata(rearranged_genome)$translocations)
  
}

# join all the dfs 
write.table(df_ins, file=paste(outdir, "insertions.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(df_del, file=paste(outdir, "deletions.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(df_inv, file=paste(outdir, "inversions.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(df_dup, file=paste(outdir, "tandemDuplications.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(df_tra, file=paste(outdir, "translocations.tab", sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# join genomes amd write
final_rearranged_genome = gDNA_rearranged_genome
final_rearranged_genome[mtDNA_chromosomes[1]] = mtDNA_rearranged_genome
writeXStringSet(final_rearranged_genome, filepath = paste(outdir, "rearranged_genome.fasta", sep="/"), format = "fasta")

# write the report
writeLines(report_genome_rearranging, con=paste(outdir, "report_generating_genome.tab", sep="/"))

