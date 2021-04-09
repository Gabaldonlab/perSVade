#!/usr/bin/env python

# This is a script that runs the testing of perSVade on several species. THere are functions and more info in testing_functions.py

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

# import testing functions
sys.path.insert(0, "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir)
import testing_functions as test_fun

# get the cluster name
if running_in_cluster is True:

    cluster_name = fun.get_current_clusterName_mareNostrum()
    if cluster_name=="MN4": threads = 24
    elif cluster_name=="Nord3": threads = 16
    else: raise ValueError("cluster could not be identified")

# define paths
perSVade_py = "%s/perSVade.py"%perSVade_dir

# define dirs
outdir_testing = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies"%ParentDir; fun.make_folder(outdir_testing)
CurDir = "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir; fun.make_folder(outdir_testing)
outdir_genomes_and_annotations = "%s/scripts/perSVade/perSVade_repository/testing/genomes_and_annotations"%ParentDir

################################

# define the table for C. glabrata
close_shortReads_table_Cglabrata = "%s/scripts/perSVade/perSVade_repository/testing/Cglabrata_table_short_reads.tab"%ParentDir

# define the type of run
run_in_cluster = True
if running_in_cluster is False: run_in_cluster = False

# define a dir with the STDs of the normal run's testing accuracy
all_STDs_dir = "%s/all_STDs_testingAccuracySeveralSpecies"%outdir_testing; fun.make_folder(all_STDs_dir)

# go through each species
for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in test_fun.species_Info:
    print(taxID, spName)

    #if spName=="Candida_glabrata": continue # debug

    # define  the genome and annotations
    genome = "%s/%s.fasta"%(outdir_genomes_and_annotations, spName)
    gff = "%s/%s.gff"%(outdir_genomes_and_annotations, spName)

    # create an outdir
    outdir_perSVade = "%s/%s_%s"%(outdir_testing, taxID, spName); fun.make_folder(outdir_perSVade)

    # get the repeats for this genome
    previous_repeats_table = fun.get_repeat_maskerDF(genome, threads=threads, replace=False)[1]

    # this is testing the whole perSVade pipeline on 3 runs of 3 close taxIDs to the reference genome. It will run only SV calling.

    # record the used resources in this run (this should be only implemented when there are no running jobs)
    test_fun.keep_STDfiles_nord3Runs_testingAccuracy(all_STDs_dir, outdir_perSVade, spName) 

    # define the table with short reads
    if spName=="Candida_glabrata": close_shortReads_table = close_shortReads_table_Cglabrata
    else: close_shortReads_table = "auto"

    # define the ploidy
    if ploidy==1: simulation_ploidies = "haploid"
    elif ploidy==2: simulation_ploidies = "diploid_hetero"
    else: raise ValueError("ploidy %i is not valid"%ploidy)

    # get the blastn of the genome against itself
    simulate_SVs_arround_HomologousRegions_maxEvalue = 0.00001
    simulate_SVs_arround_HomologousRegions_queryWindowSize = 500
    simulate_SVs_arround_HomologousRegions_previousBlastnFile = fun.get_blastn_regions_genome_against_itself(genome, simulate_SVs_arround_HomologousRegions_maxEvalue, simulate_SVs_arround_HomologousRegions_queryWindowSize, False, threads)

    # get the reads from SRA. 3 samples, 3 runs per sample. Process with the. --verbose
    cmd = "%s --ref %s --threads %i -o %s --close_shortReads_table %s --target_taxID %s --n_close_samples 3 --nruns_per_sample 1 -f1 skip -f2 skip --mitochondrial_chromosome %s --testAccuracy --verbose --max_coverage_sra_reads %i --gff %s --nsimulations 2 --skip_CNV_calling --simulation_ploidies %s --previous_repeats_table %s --StopAfter_testAccuracy --simulate_SVs_arround_HomologousRegions_maxEvalue %.10f --simulate_SVs_arround_HomologousRegions_queryWindowSize %i --simulate_SVs_arround_HomologousRegions_previousBlastnFile %s --StopAfter_testAccuracy_perSVadeRunning"%(perSVade_py, genome, threads, outdir_perSVade, close_shortReads_table, taxID, mitochondrial_chromosome, max_coverage_sra_reads, gff, simulation_ploidies, previous_repeats_table, simulate_SVs_arround_HomologousRegions_maxEvalue, simulate_SVs_arround_HomologousRegions_queryWindowSize, simulate_SVs_arround_HomologousRegions_previousBlastnFile)

    """ 
    Relevant args
    --skip_SVcalling
    --skip_CNV_calling: This is always used because this testing is about SV calling
    --StopAfter_readObtentionFromSRA
    --StopAfter_obtentionOFcloseSVs

    # --StopAfter_testAccuracy_perSVadeRunning --slurm_constraint, --StopAfter_obtentionOFcloseSVs --gff %s. Need to add the ploidy (-p ploidy) min_CNV_size # replace_SV_CNVcalling_and_optimisation --skip_cleaning_simulations_files_and_parameters --skip_repeat_analysis, --StopAfterPrefecth_of_reads, --StopAfter_sampleIndexingFromSRA

    """
    
    # add options depending on the machine
    if run_in_cluster is True: cmd += " --job_array_mode job_array"
    else: cmd += " --job_array_mode local"

    cmd_output = "%s/cmd_testing.std"%outdir_perSVade
    print("running std into %s"%cmd_output)
    fun.run_cmd("%s > %s 2>&1"%(cmd, cmd_output)) # run with stdout
    #fun.run_cmd(cmd); continue # run locally

    ###### RUN JOB ARRAYS ######

    # get the jobs file to run
    all_lines_jobfile = [l for l in open(cmd_output, "r").readlines() if l.startswith("You need to successfully run all jobs in")]

    if len(all_lines_jobfile)==1 and run_in_cluster is True:

        jobs_filename = [x for x in all_lines_jobfile[-1].split() if x.startswith("/gpfs/projects/bsc40/mschikora")][0]

        # define parameters
        name = "%s_jobs"%spName
     
        # run jobs
        if cluster_name=="MN4": 

            queue = "bsc_ls"
            time = "48:00:00"
            nodes = 3

            fun.run_jobarray_file_MN4_greasy(jobs_filename, name, time=time, queue=queue, threads_per_job=threads, nodes=nodes)

        elif cluster_name=="Nord3": 

            # define resources
            if spName in {"Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans"}:
                time = "48:00:00"
                RAM_per_thread = 2000

            elif spName in {"Arabidopsis_thaliana", "Drosophila_melanogaster"}:
                time = "48:00:00"
                RAM_per_thread = 4000 # first 3600

            else: raise ValueError("%s is not a valid spName"%spName)

            # run
            fun.run_jobarray_file_Nord3(jobs_filename, name, time=time, queue="bsc_ls", threads_per_job=threads, RAM_per_thread=RAM_per_thread, max_njobs_to_run=10000)

    elif len(all_lines_jobfile)!=0: raise ValueError("something went wrong")

    ############################

    #if taxID=="5476": adkjhdakg # stop after C. albicans

print("the testing of several species finsihed susccesffully")




