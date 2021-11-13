#!/usr/bin/env python

# This script iterates through all species and calculates the insert sizes for the bam files used in the simulations.

##### DEFINE ENVIRONMENT #######

# module imports
import os
import sys
import pandas as pd

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    running_in_cluster = False    
    threads = 4
else:
    running_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
        
# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions
print("importing functions")
import sv_functions as fun
fun.printing_verbose_mode = False

# import testing functions
sys.path.insert(0, "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir)
import testing_functions as test_fun

# get the cluster name
if running_in_cluster is True:

    cluster_name = fun.get_current_clusterName_mareNostrum()
    if cluster_name=="MN4": threads = 48
    elif cluster_name=="Nord3": threads = 16
    else: raise ValueError("cluster could not be identified")

# define paths
perSVade_py = "%s/perSVade.py"%perSVade_dir

# define dirs
outdir_testing_human = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_humanGoldenSet"%ParentDir
outdir_testing_otherSpecies = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies"%ParentDir
CurDir = "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir
outdir_genomes_and_annotations = "%s/scripts/perSVade/perSVade_repository/testing/genomes_and_annotations"%ParentDir

################################

# print the insert size len for several species
for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in test_fun.species_Info:
	print(spName)

	# define the outdir of several samples
	if spName=="Candida_glabrata": outdir_severalSamples = "%s/%s_%s/findingRealSVs_providedCloseReads/all_realVars"%(outdir_testing_otherSpecies, taxID, spName)
	else: outdir_severalSamples = "%s/%s_%s/findingRealSVs_automaticFindingOfCloseReads/all_realVars"%(outdir_testing_otherSpecies, taxID, spName)

	# go through each sample
	for s in os.listdir(outdir_severalSamples):

		# debug
		if not s.startswith("shortReads"): continue

		# get the bam file
		sorted_bam = "%s/%s/aligned_reads.bam.sorted"%(outdir_severalSamples, s)

		# calculate the insert size
		median_insert_size, median_insert_size_sd  = fun.get_insert_size_distribution(sorted_bam, replace=False, threads=threads)
		fun.print_with_runtime("%s, %s. median insert size is %i, absolute deviation of %i"%(spName, s, median_insert_size, median_insert_size_sd))

# print for human
outdir_severalSamples_human = "%s/running_on_hg38/findingRealSVs_providedCloseReads/all_realVars"%(outdir_testing_human)

# go through each sample
for s in os.listdir(outdir_severalSamples_human):

	# debug
	if not s.startswith("shortReads"): continue

	# get the bam file
	sorted_bam = "%s/%s/aligned_reads.bam.sorted"%(outdir_severalSamples_human, s)

	# calculate the insert size
	median_insert_size, median_insert_size_sd  = fun.get_insert_size_distribution(sorted_bam, replace=False, threads=threads)
	fun.print_with_runtime("Homo_sapiens, %s. median insert size is %i, absolute deviation of %i"%(s, median_insert_size, median_insert_size_sd))

# all the insert sizes make sense