#!/usr/bin/env python

# This is a script that runs the plots of the testing of perSVade on several species

#%% DEFINE ENVIRONMENT

# module imports
import os
import sys
import pandas as pd
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False    
    threads=4

else:
    run_in_cluster = True    
    threads=48
    ParentDir = "/gpfs/projects/bsc40/mschikora"
        
# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions
print("importing functions")
import sv_functions as fun

# import testing functions
sys.path.insert(0, "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir)
import testing_functions as test_fun

# define paths
perSVade_py = "%s/perSVade.py"%perSVade_dir

# define dirs
outdir_testing = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies"%ParentDir; fun.make_folder(outdir_testing)
outdir_testing_GoldenSet = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies_goldenSet"%ParentDir; fun.make_folder(outdir_testing)
outdir_testing_human = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_humanGoldenSet/running_on_hg38"%ParentDir

CurDir = "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir; fun.make_folder(outdir_testing)
genomes_and_annotations_dir = "%s/genomes_and_annotations"%CurDir

# define plots dir
PlotsDir = "%s/plots"%CurDir; fun.make_folder(PlotsDir)
ProcessedDataDir = "%s/processed_data"%PlotsDir; fun.make_folder(ProcessedDataDir)

#%% GET PROCESSED DFs

# get a cross accuracy df from the real SVs, both based on the golden set and the human real set of SVs
df_cross_accuracy_benchmark_realSVs = test_fun.get_cross_accuracy_df_realSVs(CurDir, ProcessedDataDir, threads=threads, replace=False)

# get cross-accuracy measurements testing on simulations (it already includes the human hg38 as testing)
df_cross_accuracy_benchmark = test_fun.get_cross_accuracy_df_several_perSVadeSimulations(outdir_testing, outdir_testing_human, genomes_and_annotations_dir, replace=False)

# load golden set df that comes from the testing (it does not include the human testing)
df_goldenSetAccuracy = test_fun.get_accuracy_df_goldenSet(outdir_testing_GoldenSet)

# load used parameters (this already includes the human hg38)
df_parameters_used = test_fun.get_used_parameters_testing_several_species(outdir_testing, outdir_testing_human)

#%% PLOT USED RESOURCES

# it plots the used memory and time in the testing on simulations
threads = 48
test_fun.plot_used_resources_testing_on_simulations(CurDir, ProcessedDataDir, PlotsDir, threads)

#%% PLOT USED PARAMETERS

filename = "%s/used_parameters_testing_several_species.pdf"%PlotsDir
df_plot = test_fun.get_heatmaps_used_parameters(df_parameters_used, filename)

#%% ACCURACY ON SIMULATIONS VS DEFAULT

# define the accuracy_f
accuracy_f = "precision"; width_multiplier = 1.6; legend_deviation = 2.1 
accuracy_f = "recall"; width_multiplier = 1.6; legend_deviation = 2.1 
#accuracy_f = "Fvalue"; width_multiplier = 2.8; legend_deviation = 1.9

# all data
fileprefix = "%s/accuracy_simulations_vs_default"%PlotsDir
test_fun.get_crossbenchmarking_distributions_default_and_best(df_cross_accuracy_benchmark, fileprefix, accuracy_f=accuracy_f, width_multiplier=width_multiplier, legend_deviation=legend_deviation)

#%% CROSS BENCHMARKING PLOT SIMULATIONS

# define the accuracy_f
accuracy_f = "precision"; # it could be Fvalue, precision or recall

# all data
fileprefix = "%s/all_cross_accuracy"%PlotsDir
test_fun.generate_heatmap_accuracy_of_parameters_on_test_samples(df_cross_accuracy_benchmark, fileprefix, replace=False, threads=4, accuracy_f=accuracy_f, svtype="integrated", col_cluster = False, row_cluster = False, show_only_species_and_simType=True)

#%% CROSS BENCHMARKING PLOT SIMULATIONS ONLY UNIFORM SVs

df_plot = df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.parms_typeSimulations.isin({"fast", "uniform"})) & (df_cross_accuracy_benchmark.test_typeSimulations=="uniform")]

# all data
fileprefix = "%s/all_cross_accuracy_onlyUniform"%PlotsDir
test_fun.generate_heatmap_accuracy_of_parameters_on_test_samples(df_plot, fileprefix, replace=False, threads=4, accuracy_f="Fvalue", svtype="integrated", col_cluster = False, row_cluster = False, show_only_species_and_simType=True, multiplier_width_colorbars=1.5, show_only_species=True)

#%% CROSSACCURACY DISTRIBUTION PER CATHEGORIES

# define the accuracy_f
accuracy_f = "Fvalue"; # it could be Fvalue, precision or recall

# plots the cross accuracy in a distributions-like manner
fileprefix = "%s/cross_accuracy_distribution"%PlotsDir
test_fun.get_crossbenchmarking_distributions_differentSetsOfParameters(df_cross_accuracy_benchmark, fileprefix, accuracy_f=accuracy_f, svtype="integrated")

#%% CROSSACCURACY REAL SVs HEATMAP

# define the accuracy_f
accuracy_f = "Fvalue"; # it could be Fvalue, precision or recall

# all data
fileprefix = "%s/all_cross_accuracy_realSVs"%PlotsDir
test_fun.generate_heatmap_accuracy_of_parameters_on_test_samples_realSVs(df_cross_accuracy_benchmark_realSVs, fileprefix, replace=False, threads=4, accuracy_f=accuracy_f, svtype="integrated", col_cluster = False, row_cluster = False)

#%% CROSSACCURACY DISTRIBUTION REAL SVs PER CATHEGORIES

# define the accuracy_f
accuracy_f = "Fvalue"; # it could be Fvalue, precision or recall

# plots the cross accuracy in a distributions-like manner
fileprefix = "%s/cross_accuracy_distribution_realSVs"%PlotsDir
ax = test_fun.get_crossbenchmarking_distributions_differentSetsOfParameters_realSVs(df_cross_accuracy_benchmark_realSVs, fileprefix, accuracy_f=accuracy_f, svtype="deletions")

# plot scatterplot as compared to optimised parameters
#fileprefix = "%s/cross_accuracy_scatter_realSVs"%PlotsDir
#ax = test_fun.get_crossbenchmarking_distributions_differentSetsOfParameters_realSVs_scatter(df_cross_accuracy_benchmark_realSVs, fileprefix, accuracy_f=accuracy_f, svtype="integrated")

#%% GOLDEN ACCURACY BAR PLOTS (ONLY DEFAULT AND OPTIMISED PARAMETERS)

fileprefix = "%s/goldenSetAccuracy"%PlotsDir
test_fun.plot_goldenSet_accuracy_barplots(df_goldenSetAccuracy, fileprefix, accuracy_f="Fvalue", svtype="integrated")

#%% GOLDEN ACCURACY LINE PLOT (ONLY DEFAULT AND OPTIMISED PARAMETERS)

fileprefix = "%s/goldenSetAccuracy_lineplot"%PlotsDir
test_fun.plot_goldenSet_accuracy_lineplots(df_goldenSetAccuracy, fileprefix, accuracy_f="recall", svtype="integrated")

#%% PRINT CONTENT OF THE FINAL SV_CATHEGORY

typeSim = "realSVs"

# prints the content of the SV dfs
for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in test_fun.species_Info:

    outdir_species = "%s/%s_%s/testing_Accuracy/%s"%(outdir_testing, taxID, spName, typeSim)
    for sampleID in os.listdir(outdir_species):
        
        SV_CNV_file = "%s/%s/SVcalling_output/SV_and_CNV_variant_calling.vcf"%(outdir_species, sampleID)
        SV_CNV = fun.get_vcf_df_with_INFO_as_single_fields(fun.get_df_and_header_from_vcf(SV_CNV_file)[0])
        
        adkhdjg

#%% PRINT OUT
print("testing several species finished")
sys.exit(0)
