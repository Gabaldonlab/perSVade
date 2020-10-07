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

# get the genome as an object
genome_obj = readDNAStringSet(argv$input_genome)
names(genome_obj) = lapply(names(genome_obj), function(x) (strsplit(x, "(\t)|( )")[[1]][1]))
all_chromosomes = names(genome_obj)

# insertions
print("reading insertions")
if (!is.na(argv$insertions_file)){

	insertions_df = read.table(argv$insertions_file, header=TRUE)[,c("ChrA","StartA","EndA","ChrB","StartB","EndB","Copied")]
	colnames(insertions_df) = c("chr","start","end","chrB","startB","endB","Copied")
	insertions_df$Copied = as.logical(insertions_df$Copied)
	regionsIns = GRanges(insertions_df)

} else { 

	chr = c(all_chromosomes[1]); start = c(1); end = c(1); chrB = c(all_chromosomes[1]); startB = c(1); endB = c(1); Copied = c(TRUE)
	regionsIns = GRanges(data.frame(chr, start, end, chrB, startB, endB, Copied))

}


print("reading translocations")
if (!is.na(argv$translocations_file)) {

	tra_df = read.table(argv$translocations_file, header=TRUE)[,c("ChrA","StartA","EndA","ChrB","StartB","EndB", "Balanced")]
	colnames(tra_df) = c("chr","start","end","chrB","startB","endB","Balanced")
	tra_df$Balanced = as.logical(tra_df$Balanced)
	regionsTrans = GRanges(tra_df)

} else { 

	chr = c(all_chromosomes[1]); start = c(1); end = c(1); chrB = c(all_chromosomes[1]); startB = c(1); endB = c(1); Balanced = c(TRUE)
	regionsTrans = GRanges(data.frame(chr, start, end, chrB, startB, endB, Balanced))
}

# load the single-region vars
print("reading small vars")

if (!is.na(argv$deletions_file)) {
	regionsDels = GRanges(read.table(argv$deletions_file, header=TRUE)[,c("Chr","Start","End")])
} else { 
	chr = c(all_chromosomes[1]); start = c(1); end = c(1);
	regionsDels = GRanges(data.frame(chr, start, end))
}

if (!is.na(argv$inversions_file)) {
	regionsInvs = GRanges(read.table(argv$inversions_file, header=TRUE)[,c("Chr","Start","End")])
} else { 
	chr = c(all_chromosomes[1]); start = c(1); end = c(1);
	regionsInvs = GRanges(data.frame(chr, start, end))
 }

if (!is.na(argv$tandemDuplications_file)) {
	regionsDups = GRanges(read.table(argv$tandemDuplications_file, header=TRUE)[,c("Chr","Start","End")])
} else { 
	chr = c(all_chromosomes[1]); start = c(1); end = c(1);
	regionsDups = GRanges(data.frame(chr, start, end))
}

# get the rearranged genome
rearranged_genome = simulateSV(output=NA, genome=genome_obj, random=FALSE, verbose=TRUE, regionsDels=regionsDels, regionsInvs=regionsInvs, regionsIns=regionsIns, regionsDups=regionsDups, regionsTrans=regionsTrans, bpSeqSize=50, maxDups=4, percBalancedTrans=1)

# test that all the insetions are cut-and-paste (this is so because after this they should be changed to copy-paste )
if (!is.na(argv$insertions_file)){
	print("checking that insertions are ok")

	copy_column_insertions = unique(metadata(rearranged_genome)$insertions$Copied)
	if (copy_column_insertions!=c(FALSE)){ stop("ERROR: All the insertions should be cut-and-paste in the simulation. The copy-and-paste will be generated afterwards with a custom script") }
}

# write the rearranged genome
writeXStringSet(rearranged_genome, filepath = argv$output_genome, format = "fasta")
