#!/usr/bin/env Rscript

# this script generates an SV genome under outdir given the precise location of each SV

# define environment
library(argparser, quietly=TRUE)
library(RSVSim, quietly=TRUE) # doc in https://rdrr.io/bioc/RSVSim/man/simulateSV.html
library(emdbook, quietly=TRUE)
library(rtracklayer, quietly=TRUE)
library(GenomicRanges, quietly=TRUE)
library(R.utils, quietly=TRUE)

# parse cmd line args
argp = arg_parser("Takes a genome and generates a simulated genome with rearrangements in outdir. It will generate these rearrangements. It requires all the desired rearrangements files")

argp = add_argument(argp, "--input_genome", help="Path to the genome where to generate the SV")
argp = add_argument(argp, "--output_genome", help="Path to the rearranged genome")

argp = add_argument(argp, "--insertions_file", help="The path to the insertions")
argp = add_argument(argp, "--translocations_file", help="The path to the translocations")
argp = add_argument(argp, "--deletions_file", help="The path to the deletions")
argp = add_argument(argp, "--tandemDuplications_file", help="The path to the tandemDuplications")
argp = add_argument(argp, "--inversions_file", help="The path to the inversions")

argv = parse_args(argp)

# load examples for debugging
#argv$insertions_file = "/home/mschikora/samba/scripts/perSVade/perSVade_repository/testing/Cglabrata/RUN1_CST34_2G_FLZ_subsampled/outdir_perSVade/SVdetection_output/test_FindSVinAssembly/simulation_0/final_simulated_SVs/insertions.tab"
#argv$deletions_file = "/home/mschikora/samba/scripts/perSVade/perSVade_repository/testing/Cglabrata/RUN1_CST34_2G_FLZ_subsampled/outdir_perSVade/SVdetection_output/test_FindSVinAssembly/simulation_0/final_simulated_SVs/deletions.tab"
#argv$inversions_file = "/home/mschikora/samba/scripts/perSVade/perSVade_repository/testing/Cglabrata/RUN1_CST34_2G_FLZ_subsampled/outdir_perSVade/SVdetection_output/test_FindSVinAssembly/simulation_0/final_simulated_SVs/inversions.tab"
#argv$tandemDuplications_file = "/home/mschikora/samba/scripts/perSVade/perSVade_repository/testing/Cglabrata/RUN1_CST34_2G_FLZ_subsampled/outdir_perSVade/SVdetection_output/test_FindSVinAssembly/simulation_0/final_simulated_SVs/tandemDuplications.tab"
#argv$translocations_file = "/home/mschikora/samba/scripts/perSVade/perSVade_repository/testing/Cglabrata/RUN1_CST34_2G_FLZ_subsampled/outdir_perSVade/SVdetection_output/test_FindSVinAssembly/simulation_0/final_simulated_SVs/translocations.tab"
#argv$input_genome = "/home/mschikora/samba/scripts/perSVade/perSVade_repository/testing/Cglabrata/RUN1_CST34_2G_FLZ_subsampled/outdir_perSVade/reference_genome_dir/reference_genome.fasta"


# load the regions into gr_objects

# insertions
insertions_df = read.table(argv$insertions_file, header=TRUE)[,c("ChrA","StartA","EndA","ChrB","StartB","EndB","Copied")]
colnames(insertions_df) = c("chr","start","end","chrB","startB","endB","Copied")
regionsIns = GRanges(insertions_df)

# change trans to be 0-based
tra_df = read.table(argv$translocations_file, header=TRUE)[,c("ChrA","StartA","EndA","ChrB","StartB","EndB", "Balanced")]
colnames(tra_df) = c("chr","start","end","chrB","startB","endB","Balanced")
regionsTrans = GRanges(tra_df)

# load the single-region vars
regionsDels = GRanges(read.table(argv$deletions_file, header=TRUE)[,c("Chr","Start","End")])
regionsInvs = GRanges(read.table(argv$inversions_file, header=TRUE)[,c("Chr","Start","End")])
regionsDups = GRanges(read.table(argv$tandemDuplications_file, header=TRUE)[,c("Chr","Start","End","Duplications")])

# get the genome as an object
genome_obj = readDNAStringSet(argv$input_genome)
names(genome_obj) = lapply(names(genome_obj), function(x) (strsplit(x, "(\t)|( )")[[1]][1]))

# get the rearranged genome
rearranged_genome = simulateSV(output=NA, genome=genome_obj, random=FALSE, verbose=TRUE, regionsDels=regionsDels, regionsInvs=regionsInvs, regionsIns=regionsIns, regionsTrans=regionsTrans, regionsDups=regionsDups)

# write the rearranged genome
writeXStringSet(rearranged_genome, filepath = argv$output_genome, format = "fasta")
