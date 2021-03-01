#!/usr/bin/env python

# this is a script to test the running of perSVade on a golden set (ONT reads). It should be run when testing_several_species.py is finished

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
    if cluster_name=="MN4": threads = 48
    elif cluster_name=="Nord3": threads = 16
    else: raise ValueError("cluster could not be identified")

# define paths
perSVade_py = "%s/perSVade.py"%perSVade_dir

# define dirs
outdir_testing = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies_goldenSet"%ParentDir; fun.make_folder(outdir_testing)
CurDir = "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir; fun.make_folder(outdir_testing)
outdir_genomes_and_annotations = "%s/scripts/perSVade/perSVade_repository/testing/genomes_and_annotations"%ParentDir

################################

# define the taxIDs that have no
taxIDs_with_noON_overalpping = {"746128"}

# define the golden set dir of C. glabrata
goldenSet_dir_Cglabrata = "%s/scripts/perSVade/perSVade_repository/testing/Cglabrata_goldenSetReads_BG2"%ParentDir

# define the run in cluster (and debug)
run_in_cluster = True
if running_in_cluster is False: run_in_cluster = False

# go through each species
for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in test_fun.species_Info:
    print(taxID, spName)

    #if spName=="Candida_glabrata": continue # debug

    # define  the genome and annotations
    genome = "%s/%s.fasta"%(outdir_genomes_and_annotations, spName)
    gff = "%s/%s.gff"%(outdir_genomes_and_annotations, spName)

    # create an outdir
    outdir_perSVade = "%s/%s_%s"%(outdir_testing, taxID, spName); fun.make_folder(outdir_perSVade)

    ###### TEST SVIM AND SNIFFLES NORMALISATION ######
    """
    svtype_to_file_svim = fun.get_svim_output_as_perSVade("%s/testing_goldenSetAccuracy/ONT_SV_calling/svim_output"%outdir_perSVade)
    print(svtype_to_file_svim)
    jagjdadjgadj
    svtype_to_file_sniffles = fun.get_sniffles_output_as_perSVade("%s/testing_goldenSetAccuracy/ONT_SV_calling/sniffles_output"%outdir_perSVade)
    adkhadhjgdahda
    """

    ##################################################

    # get the repeats for this genome
    previous_repeats_table = fun.get_repeat_maskerDF(genome, threads=threads, replace=False)[1]

    # define the goldenSet_dir
    if spName=="Candida_glabrata": goldenSet_dir = goldenSet_dir_Cglabrata
    else: goldenSet_dir = "auto"

    # define the real bedepe breakpoints
    real_bedpe_breakpoints = "%s/"%outdir_perSVade

    # get the golden set running 
    if taxID in taxIDs_with_noON_overalpping: continue
    cmd = "%s --ref %s --threads %i -o %s --target_taxID %s --real_bedpe_breakpoints %s -f1 skip -f2 skip --mitochondrial_chromosome %s --gff %s --goldenSet_dir %s --skip_SVcalling --verbose --nsimulations 2 --simulation_ploidies haploid,diploid_hetero --previous_repeats_table %s --QC_and_trimming_reads --StopAfter_goldenSetAnalysis "%(perSVade_py, genome, threads, outdir_perSVade, taxID, real_bedpe_breakpoints, mitochondrial_chromosome, gff, goldenSet_dir, previous_repeats_table)

    """
    StopAfter_goldenSetAnalysis
    StopAfterPrefecth_of_reads_goldenSet
	StopAfter_goldenSetAnalysis_readObtention
	StopAfter_goldenSetAnalysis_readTrimming
    StopAfter_sampleIndexingFromSRA

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
            nodes = 2

            fun.run_jobarray_file_MN4_greasy(jobs_filename, name, time=time, queue=queue, threads_per_job=threads, nodes=nodes)

        elif cluster_name=="Nord3": 

            queue = "bsc_ls"; 
            RAM_per_thread = 5000; # 1800 or 5000 
            time = "48:00:00" # per job

            fun.run_jobarray_file_Nord3(jobs_filename, name, time=time, queue=queue, threads_per_job=threads, RAM_per_thread=RAM_per_thread, max_njobs_to_run=10000)

    elif len(all_lines_jobfile)!=0: raise ValueError("something went wrong")

    ############################

print("Golden set analysis worked")