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

else:
    run_in_cluster = True    
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

CurDir = "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir; fun.make_folder(outdir_testing)
genomes_and_annotations_dir = "%s/genomes_and_annotations"%CurDir

# define plots dir
PlotsDir = "%s/plots"%CurDir; fun.make_folder(PlotsDir)
ProcessedDataDir = "%s/processed_data"%PlotsDir; fun.make_folder(ProcessedDataDir)

#%% GET PROCESSED DFs

# get cross-accuracy measurements
df_cross_accuracy_benchmark = test_fun.get_cross_accuracy_df_several_perSVadeSimulations(outdir_testing, genomes_and_annotations_dir, replace=False)

# load golden set df
df_goldenSetAccuracy = test_fun.get_accuracy_df_goldenSet(outdir_testing_GoldenSet)

# load used parameters 
df_parameters_used = test_fun.get_used_parameters_testing_several_species(outdir_testing)

#%% PLOT USED PARAMETERS

filename = "%s/used_parameters_testing_several_species.pdf"%PlotsDir
df_plot = test_fun.get_heatmaps_used_parameters(df_parameters_used, filename)

#%% CROSS BENCHMARKING PLOT

# all data
fileprefix = "%s/all_cross_accuracy"%PlotsDir
test_fun.generate_heatmap_accuracy_of_parameters_on_test_samples(df_cross_accuracy_benchmark, fileprefix, replace=False, threads=4, accuracy_f="Fvalue", svtype="integrated", col_cluster = False, row_cluster = False, show_only_species_and_simType=True)

#%% CROSS BENCHMARKING PLOT ONLY RANDOM SVs

df_plot = df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.parms_typeSimulations.isin({"fast", "uniform"})) & (df_cross_accuracy_benchmark.test_typeSimulations=="uniform")]

# all data
fileprefix = "%s/all_cross_accuracy_onlyUniform"%PlotsDir
test_fun.generate_heatmap_accuracy_of_parameters_on_test_samples(df_plot, fileprefix, replace=False, threads=4, accuracy_f="Fvalue", svtype="integrated", col_cluster = False, row_cluster = False, show_only_species_and_simType=True, multiplier_width_colorbars=1.5, show_only_species=True)


#%% CROSS-ACCURACY JITTER PLOTS PER CATHEGORY
    
fileprefix = "%s/cross_accuracy_distribution"%PlotsDir
ax  = test_fun.get_crossaccuracy_distributions(df_cross_accuracy_benchmark, fileprefix, accuracy_f="precision",  svtype="integrated")
    


#%% CROSS-ACCURACY JITTER PLOTS PER CATHEGORY ONE SINGLE TYPE
    
fileprefix = "%s/cross_accuracy_distribution_onlyUniform_SameRun"%PlotsDir
test_fun.plot_accuracy_distributions_sameRun_bySpecies(df_cross_accuracy_benchmark, fileprefix, accuracy_f="Fvalue", all_types_simulations={"fast", "uniform"})
    
#%% GOLDEN ACCURACY BAR PLOTS

fileprefix = "%s/goldenSetAccuracy"%PlotsDir
test_fun.plot_goldenSet_accuracy_barplots(df_goldenSetAccuracy, fileprefix, accuracy_f="Fvalue", svtype="integrated")

#%% GOLDEN ACCURACY LINE PLOT

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

#%% 
#%% 
#%% 
#%% 
#%% 
#%% 
#%% 
#%% 
#%%
#%% CODE GRAVEYARD

##%%  CROSS BENCHMARKING SEPPARATE

# each parms species
plots_dir_parms_species = "%s/parms_only_one_species"%PlotsDir; fun.make_folder(plots_dir_parms_species)
for parms_species in ['Arabidopsis_thaliana', 'Candida_albicans', 'Candida_glabrata', 'Cryptococcus_neoformans', 'Drosophila_melanogaster']:
    df_sp = df_cross_accuracy_benchmark[df_cross_accuracy_benchmark.parms_species.isin({parms_species, "none"})]
    
    fileprefix = "%s/%s"%(plots_dir_parms_species, parms_species)
    test_fun.generate_heatmap_accuracy_of_parameters_on_test_samples(df_sp, fileprefix, replace=False, threads=4, accuracy_f="Fvalue", svtype="integrated", col_cluster = False, row_cluster = False)

# each test species
plots_dir_test_species = "%s/test_only_one_species"%PlotsDir; fun.make_folder(plots_dir_test_species)
for test_species in ['Arabidopsis_thaliana', 'Candida_albicans', 'Candida_glabrata', 'Cryptococcus_neoformans', 'Drosophila_melanogaster']:
    df_sp = df_cross_accuracy_benchmark[df_cross_accuracy_benchmark.test_species.isin({test_species, "none"})]
    
    fileprefix = "%s/%s"%(plots_dir_test_species, test_species)
    test_fun.generate_heatmap_accuracy_of_parameters_on_test_samples(df_sp, fileprefix, replace=False, threads=4, accuracy_f="Fvalue", svtype="integrated", col_cluster = False, row_cluster = False)
    
# only one species against itself
plots_dir_species = "%s/only_one_species"%PlotsDir; fun.make_folder(plots_dir_species)
for species in ['Arabidopsis_thaliana', 'Candida_albicans', 'Candida_glabrata', 'Cryptococcus_neoformans', 'Drosophila_melanogaster']:
    df_sp = df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.test_species.isin({species, "none"})) & (df_cross_accuracy_benchmark.parms_species.isin({species, "none"}))]
    
    fileprefix = "%s/%s"%(plots_dir_species, species)
    test_fun.generate_heatmap_accuracy_of_parameters_on_test_samples(df_sp, fileprefix, replace=False, threads=4, accuracy_f="precision", svtype="integrated", col_cluster = False, row_cluster = False)
    
########### COLLECT DATA ON EACH SAMPLE ########### 

# we will collect data of the simulations accuracy and the number of each type of dat

# define files
all_df_simulations_accuracy_file = "%s/all_df_simulations_accuracy.tab"%ProcessedDataDir
all_df_nvariants_file = "%s/all_df_nvariants.tab"%ProcessedDataDir


if any([fun.file_is_empty(f) for f in [all_df_simulations_accuracy_file, all_df_nvariants_file]]):
    print("getting dfs")

    species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138"),
                    ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314"),
                    ("5207", "Cryptococcus_neoformans", 1, "CP003834.1"),
                    ("746128", "Aspergillus_fumigatus", 1, "CM016889.1"),
                    ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1"),
                    ("7227", "Drosophila_melanogaster", 2, "KJ947872.2")]

    # int dfs
    all_df_simulations_accuracy = pd.DataFrame()
    all_df_nvariants = pd.DataFrame()

    # define the samples to skip 
    samples_to_skip = {"sample2608267_ERR3514862", "sample38785_SRR7904950"}

    # go through each species
    for taxID, spName, ploidy, mitochondrial_chromosome in species_Info:
        print(taxID, spName)

        # define teh name
        short_spName = "%s. %s"%(spName.split("_")[0][0], spName.split("_")[1])


        # go through each type of optimisation
        for type_optimisation in ["fast", "realSVs", "uniform"]:

            # define the outdir
            outdir_typeOptimisation = "%s/%s_%s/testing_Accuracy/%s"%(outdir_testing, taxID, spName, type_optimisation)

            # go through each sample
            for sampleID in os.listdir(outdir_typeOptimisation):

                # skip unfinished samples
                if sampleID in samples_to_skip: continue

                # get the outdir
                outdir_sampleID = "%s/%s"%(outdir_typeOptimisation, sampleID)

                ######### SIMULATIONS ACCURACY #########

                if type_optimisation!="fast":

                    # define the df that has the accuracy
                    accuray_file = "%s/SVdetection_output/parameter_optimisation/benchmarking_all_filters_for_all_genomes_and_ploidies/df_cross_benchmark_best.tab"%outdir_sampleID
            
                    # get the accuracy df
                    df_accuracy = pd.read_csv(accuray_file, sep="\t")

                    # add the metadata of this sample
                    df_accuracy["spName"] = short_spName
                    df_accuracy["type_optimisation"] = type_optimisation
                    df_accuracy["sampleID"] = sampleID

                    all_df_simulations_accuracy = all_df_simulations_accuracy.append(df_accuracy)

                #########################################

                ######## NUMBER VARIANTS ########

                # get the svtype to svfile
                svtype_to_svfile = fun.get_svtype_to_svfile_from_perSVade_outdir(outdir_sampleID)

                # init dict
                data_dict = {}
                for Isv, svtype in enumerate({"insertions", "deletions", "tandemDuplications", "translocations", "inversions", "remaining"}):

                    # add the number of each svtype
                    if svtype not in svtype_to_svfile: nsvs = 0
                    else: nsvs = int(subprocess.check_output("wc -l %s"%(svtype_to_svfile[svtype]), shell=True).split()[0]) - 1

                    data_dict[Isv] = {"svtype":svtype, "nvars": nsvs}

                # get df
                df_nvars = pd.DataFrame(data_dict).transpose()
                df_nvars["spName"] = short_spName
                df_nvars["type_optimisation"] = type_optimisation
                df_nvars["sampleID"] = sampleID

                all_df_nvariants = all_df_nvariants.append(df_nvars)

                #################################

    # save
    fun.save_df_as_tab(all_df_nvariants, all_df_nvariants_file)
    fun.save_df_as_tab(all_df_simulations_accuracy, all_df_simulations_accuracy_file)


# load dfs
all_df_nvariants = fun.get_tab_as_df_or_empty_df(all_df_nvariants_file)
all_df_simulations_accuracy = fun.get_tab_as_df_or_empty_df(all_df_simulations_accuracy_file)

###################################################

############## PLOT THE NUMBER OF VARIANTS ##############

print("plotting number of vars")

# make a folder for the nvars
nvars_dir = "%s/nvars"%PlotsDir; fun.make_folder(nvars_dir)

# a swarmplot where the x is the species, the y is the accuracy and the hue is the svtype
for type_optimisation in ["uniform", "realSVs", "fast"]:

    # get the df
    df_plot = all_df_nvariants[all_df_nvariants.type_optimisation==type_optimisation]

    # add a pseudocount
    df_plot.nvars += 1

    fig = plt.figure(figsize=(8,3))
    ax = sns.boxplot(x="spName", y="nvars", data=df_plot, hue="svtype", boxprops=dict(alpha=.45))
    ax = sns.swarmplot(x="spName", y="nvars", hue="svtype", data=df_plot, dodge=True)

    ax.set_xlabel("species")
    ax.set_ylabel("number of variants")
    ax.set_yscale("log")
    ax.set_title("variant number for %s simulations"%(type_optimisation))
    ax.legend(bbox_to_anchor=(1, 1))
    for item in ax.get_xticklabels(): item.set_rotation(60)


    filename = "%s/%s.pdf"%(nvars_dir, type_optimisation)
    fig.savefig(filename, bbox_inches='tight')
    plt.close(fig)

#########################################################

############## PLOT THE SIMULATIONS ACCURACY ##############

print("plotting simulations accuracy")

# make a folder for the simulations accuracy
simAccuracy_dir = "%s/simulations_accuracy"%PlotsDir; fun.make_folder(simAccuracy_dir)

# a swarmplot where the x is the species, the y is the accuracy and the hue is the svtype
for type_optimisation in ["uniform", "realSVs"]:

    # get the df
    df_plot = all_df_simulations_accuracy[all_df_simulations_accuracy.type_optimisation==type_optimisation]

    for accuracy_measure in ["precision", "recall", "Fvalue"]:

        fig = plt.figure(figsize=(8,3))
        ax = sns.boxplot(x="spName", y=accuracy_measure, data=df_plot, hue="test_svtype", boxprops=dict(alpha=.45))
        ax = sns.swarmplot(x="spName", y=accuracy_measure, hue="test_svtype", data=df_plot, dodge=True)

        ax.set_xlabel("species")
        ax.set_title("%s for %s simulations"%(accuracy_measure, type_optimisation))
        ax.legend(bbox_to_anchor=(1, 1))
        for item in ax.get_xticklabels(): item.set_rotation(60)


        filename = "%s/%s_%s.pdf"%(simAccuracy_dir, type_optimisation, accuracy_measure)
        fig.savefig(filename, bbox_inches='tight')
        plt.close(fig)

###########################################################


print("plots correctly performed")
