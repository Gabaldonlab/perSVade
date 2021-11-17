#!/usr/bin/env python

# This script contains functions useful for the testing of perSVade on several samples

######### ENV ########

import os
import sys
import pandas as pd
import numpy as np
import multiprocessing as multiproc
import seaborn as sns
from matplotlib import gridspec
import subprocess
import itertools
from matplotlib.lines import Line2D
import copy as cp

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if not os.path.exists(ParentDir): 
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    run_in_cluster = True
else:
    run_in_cluster = False

# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions from perSVade
print("importing functions")
import sv_functions as fun

import matplotlib
import matplotlib.pyplot as plt

"""
if run_in_cluster is False: 
    try: matplotlib.use('TkAgg')
    except: print("setting TkAgg does not work") 
"""

######################

##### GENERAL INFO #####

"""
This is how the genomes were obtained:

C. glabrata: reference genome from CGD: the latest version by 12/03/2019, which is v_s02-m07-r35 

C. albicans: 

    ref genome CGD: http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_version_A22-s07-m01-r110_chromosomes.fasta.gz

    gff from CGD: http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_version_A22-s07-m01-r110_features.gff

    From here I keep 'haplotype A' for both files

C. neoformans: ref genome from GenBank GCA_000149245.3

A. fumigatus: 

    gDNA from reference NCBI:

    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.gff.gz
    
    mtDNA from https://www.ncbi.nlm.nih.gov/nuccore/CM016889.1

A. thaliana: ref genome from GenBank GCA_000001735.2

D. melanogaster: ref genome from GenBank GCA_000001215.4

D. rerio: ref genome from GenBank removing the alternate haplotypes. (this is GCA_000002035.4)

H. sapiens: ref genome from GenBank removing the alternate haplotypes. (this is GCA_000001405.28)

For C. glabrata I got the nanopore reads from ~/../mmarcet/nanopore/GABALDON02/assembly_files/BG2/nanopore.reads.pass.fastq.gz and the short reads from Ewa's experiment in RUN4_BG2_SRA_WT

"""

###############

####### DEFINE VARIABLES ##########

#("7955", "Danio_rerio", 2, "NC_002333.2")]
#("9606", "Homo_sapiens", 2, "NC_012920.1")]


species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", 10000000000000000),
                ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314", 10000000000000000),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]

species_Info_WithHumanHg38 = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", 10000000000000000),
                              ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314", 10000000000000000),
                              ("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000),
                              ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                              ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30),
                              ("9606", "Homo_sapiens", 2, "chrM", 10000000000000)]


species_to_taxID = {x[1]:int(x[0]) for x in species_Info}

###################################

def update_df_resources_nord3Runs_testingAccuracy(df_resources_file, outdir_perSVade, spName, all_STDs_dir):

    """Takes the outdir of perSVade for one species and add to df_resources the data of that runID"""

    #print("getting resources")

    # get the df_resources
    if not fun.file_is_empty(df_resources_file): df_resources = fun.get_tab_as_df_or_empty_df(df_resources_file)
    else: df_resources = pd.DataFrame()

    # define the roundID
    if len(df_resources)==0 or sum(df_resources.spName==spName)==0: roundID = 0
    else: roundID = max(df_resources[df_resources.spName==spName].roundID)+1

    # init the df resources of this round
    dict_resources = {}; 
    if len(df_resources)==0: uniqueID = 0
    else: uniqueID = max(df_resources.uniqueID)+1

    # keep the resources of the previous run (this is assumed to be ran in Nord3)
    for typeSim in ["fast", "uniform", "realSVs"]:
        outdir_testAccuracy = "%s/testing_Accuracy/%s"%(outdir_perSVade, typeSim)
        if not os.path.isdir(outdir_testAccuracy): continue

        # go through each sample (note that each sample is only useful when )
        for sampleID in os.listdir(outdir_testAccuracy): 

            # define files
            outdir_sample = "%s/%s"%(outdir_testAccuracy, sampleID)
            final_file = "%s/perSVade_finished_file.txt"%(outdir_sample)

            # define whether some files where generated
            perSVade_was_ran = any([f not in {"aligned_reads.bam.sorted", "aligned_reads.bam.sorted.bai"} for f in os.listdir(outdir_sample)])

            # only keep info about perSVade runs
            if perSVade_was_ran is True:

                # get the time
                time_fields = {"time_GeneralProcessing", "time_alignment", "time_all", "time_obtentionCloseSVs", "time_SVcalling", "time_SVandCNVcalling", "time_smallVarsCNV"}
                if fun.file_is_empty(final_file):
                    perSVade_finished = False
                    time_dict = {t : -1 for t in time_fields}

                else:
                    perSVade_finished = True
                    time_dict = {l.split(":")[0] : float(l.strip().split(":")[1]) for l in open(final_file, "r").readlines() if l.startswith("time_")}

                # if this run is already in the dataset, skip
                if perSVade_finished is True and len(df_resources)>0 and any((df_resources.spName==spName) & (df_resources.typeSim==typeSim) & (df_resources.sampleID==sampleID) & (df_resources.perSVade_finished==True)): continue

                # identify the command ID for this STD. 
                stddir = "%s/testing_Accuracy/STDfiles"%outdir_perSVade
                commandIDs = [int(f.split("command.")[1]) for f in os.listdir(stddir) if f.startswith("command.") and "/".join(open("%s/%s"%(stddir,f), "r").readlines()[0].strip().split("--outdir ")[1].split()[0].split("/")[-2:])=="%s/%s"%(typeSim, sampleID)]

                if len(commandIDs)>1: raise ValueError("There are these possible command IDs: %s. This is the stddir: %s, This is the expected ID: %s/%s"%(commandIDs, stddir, typeSim, sampleID))

                # if there are no cmds it means that you lost the calculus of time, but you can still keep the final time. This is important because it could be useful to get the time of the missed command
                if len(commandIDs)==0 and perSVade_finished is True:

                    #print("WARNING: THere are no commands for %s-%s, although perSVade perSVade_finished is True."%(sampleID, typeSim))

                    dict_resources[uniqueID] = {"typeSim":typeSim, "sampleID":sampleID, "jobID":-1, "commandID":-1, "std_file":"no_file", "requested_mem_Mb":-1, "max_mem_Mb":-1, "cpu_time":-1, "run_time":-1, "perSVade_finished":perSVade_finished, "spName":spName, "roundID":roundID, "uniqueID":uniqueID}
                    for f in time_fields: dict_resources[uniqueID][f] = time_dict[f]
                    uniqueID+=1

                elif len(commandIDs)==0: raise ValueError("WARNING: THere are no commands for %s-%s. perSVade_finished:%s. stddir:%s"%(sampleID, typeSim, perSVade_finished, stddir))

                # if there are commands, save
                else:
                
                    commandID = commandIDs[0]

                    # get the jobID
                    main_stdout_lines = [l.strip() for l in open("%s/%s_jobs_stdout.txt"%(stddir, spName), "r").readlines()]
                    all_jobIDs = {int(l.split("Subject: Job ")[1].split("[")[0]) for l in main_stdout_lines if l.startswith("Subject: Job")}
                    if len(all_jobIDs)!=1: raise ValueError("There are more than 1 job IDS: %s"%all_jobIDs)
                    jobID = next(iter(all_jobIDs))

                    # move the STD to all_STDs_dir
                    origin_std = "%s/accuracyRealSVs.%i.out"%(stddir, commandID)
                    dest_std = "%s/%s_job%i_command%i_%s_%s_std.txt"%(all_STDs_dir, spName, jobID, commandID, typeSim, sampleID)
                    fun.rsync_file(origin_std, dest_std)

                    # get the resource consumption
                    jobID_line = [Iline for Iline,l in enumerate(main_stdout_lines) if l.startswith("Subject: Job %i[%i]:"%(jobID, commandID))][0]
                    requested_mem_Mb = [float(l.split("Total Requested Memory :")[1].split()[0]) for l in main_stdout_lines[jobID_line:] if "Total Requested Memory :" in l and "MB" in l][0]
                    max_mem_Mb = [float(l.split("Max Memory :")[1].split()[0]) for l in main_stdout_lines[jobID_line:] if "Max Memory :" in l and "MB" in l][0]
                    cpu_time = [float(l.split("CPU time :")[1].split()[0]) for l in main_stdout_lines[jobID_line:] if "CPU time :" in l and "sec" in l][0]

                    # get the elapsed time in seconds
                    start_date =  [l.strip().split()[3:] for l in main_stdout_lines[jobID_line:] if l.startswith("Started at")][0]
                    end_date =  [l.strip().split()[4:] for l in main_stdout_lines[jobID_line:] if l.startswith("Results reported on")][0]

                    s_month = start_date[0]
                    s_day = int(start_date[1])
                    s_h, s_min, s_sec = [float(x) for x in start_date[2].split(":")]
                    s_year = int(start_date[3])

                    e_month = end_date[0]
                    e_day = int(end_date[1])
                    e_h, e_min, e_sec = [float(x) for x in end_date[2].split(":")]
                    e_year = int(end_date[3])  

                    # define the number of days that distanced the start and the end. Each month is particular
                    if s_year==e_year and s_year==2021 and s_month=="Feb" and e_month=="Mar": transcurred_days = e_day - (s_day-28)

                    # most cases
                    else: 
                        transcurred_days = e_day-s_day
                        if s_month!=e_month: raise ValueError("s_month %s is different to e_month %s"%(s_month, e_month))  
                        if s_year!=e_year: raise ValueError("s_year %s is different to e_year %s"%(s_year, e_year))  

                    # get the total time
                    run_time =  transcurred_days*(24*3600) + (e_h-s_h)*3600 + (e_min-s_min)*60 + (e_sec-s_sec)

                    # checks
                    for x in [cpu_time, max_mem_Mb, requested_mem_Mb, run_time]:
                        if pd.isna(x) or x<=0.0: raise ValueError("there is an error with the parsing of the output: %s"%x)

                
                    # keep as a dict if this job has not already been saved
                    if len(df_resources)==0 or not any((df_resources.spName==spName) & (df_resources.typeSim==typeSim) & (df_resources.sampleID==sampleID) & (df_resources.jobID==jobID) & (df_resources.commandID==commandID) & (df_resources.max_mem_Mb==max_mem_Mb) & (df_resources.perSVade_finished==perSVade_finished) & (df_resources.cpu_time==cpu_time)):

                        dict_resources[uniqueID] = {"typeSim":typeSim, "sampleID":sampleID, "jobID":jobID, "commandID":commandID, "std_file":dest_std, "requested_mem_Mb":requested_mem_Mb, "max_mem_Mb":max_mem_Mb, "cpu_time":cpu_time, "run_time":run_time, "perSVade_finished":perSVade_finished, "spName":spName, "roundID":roundID, "uniqueID":uniqueID}
                        for f in time_fields: dict_resources[uniqueID][f] = time_dict[f]

                        uniqueID+=1

            else: print("perSVade was never ran for sample %s typeSim %s"%(sampleID, typeSim))

    # save the updated df_resources
    df_resources = df_resources.append(pd.DataFrame(dict_resources).transpose())
    if len(set(df_resources.uniqueID))!=len(df_resources): raise ValueError("uniqueID is not unique")

    # make the std_file as only the file (as this can be ran from different systems)
    def get_empty_std_file_to_no_file(x):
        if pd.isna(x) or x=="": return "no_file"
        else: return x
    df_resources["std_file"] = df_resources.std_file.apply(get_empty_std_file_to_no_file).apply(fun.get_file)

    # add the pair of sample and typeSim
    df_resources["species_sampleID_sim"] = df_resources.spName + "_" + df_resources.sampleID + "_" + df_resources.typeSim

    # check that there is as much one row with perSVade finished
    species_sampleID_sim_to_nFinished = df_resources.groupby("species_sampleID_sim").apply(lambda df_sim: sum(df_sim.perSVade_finished))
    if not all(species_sampleID_sim_to_nFinished.isin({0, 1})): raise ValueError("There are some sims with more than 1 finished file")

    # print the samples that are missing and the latest std
    missing_species_sampleID_sim = list(species_sampleID_sim_to_nFinished[species_sampleID_sim_to_nFinished==0].index)
    #for x in missing_species_sampleID_sim: print(df_resources[df_resources.species_sampleID_sim==x].iloc[-1]["std_file"])

    # save
    fun.save_df_as_tab(df_resources, df_resources_file)

    return df_resources, roundID


def keep_STDfiles_nord3Runs_testingAccuracy(all_STDs_dir, outdir_perSVade, spName, type_testing="several_species"):

    """This function records all the new files from outdir_perSVade testing accuracy"""

    fun.make_folder(all_STDs_dir)

    # define the dirs
    testing_Accuracy_dir = "%s/testing_Accuracy"%outdir_perSVade
    stddir = "%s/STDfiles"%testing_Accuracy_dir

    # get the current STD file metadata. This indicates the run
    if type_testing=="several_species":

        stdout_report = "%s/%s_jobs_stdout.txt"%(stddir, spName)
        stderr_report = "%s/%s_jobs_stderr.txt"%(stddir, spName)

    elif type_testing=="human":

        stdout_report = "%s/%s_stdout.txt"%(stddir, spName)
        stderr_report = "%s/%s_stderr.txt"%(stddir, spName)


    if fun.file_is_empty(stdout_report): return

    jobIDs = {int(l.split("Subject: Job ")[1].split("[")[0]) for l in open(stdout_report, "r").readlines() if l.startswith("Subject:")}
    if len(jobIDs)!=1: raise ValueError("There has to be some jobID")
    jobID = next(iter(jobIDs))

    taskIDs =  {int(l.split("Subject: Job ")[1].split("[")[1].split("]")[0]) for l in open(stdout_report, "r").readlines() if l.startswith("Subject:")}
    if taskIDs!=set(range(1, max(taskIDs)+1)): raise ValueError("The task IDs should be from 1 to the max")

    nlines = len(open(stdout_report, "r").readlines())
        
    # define some destination files
    prefix_stds = "%s/%s_job=%i_tasks=[1-%i]_%ilines"%(all_STDs_dir, spName, jobID, max(taskIDs), nlines)
    dest_stdout_report = "%s_stdout_report.txt"%prefix_stds
    dest_stderr_report = "%s_stderr_report.txt"%prefix_stds
    dest_jobs_file = "%s_jobsfile.txt"%prefix_stds
    dest_jobEndingStatus_file = "%s_jobEnding.txt"%prefix_stds

    # if the dest report has not been written, continue
    if fun.file_is_empty(dest_stdout_report):
        print("getting STDs")

        # define the jobs file
        jobs_file = "%s/jobs.testingRealDataAccuracy"%(testing_Accuracy_dir)
        if fun.file_is_empty(jobs_file): raise ValueError("The jobs file should exist")

        # define a file for each job that indicates whether it finished or it was ran
        dict_data = {}
        for taskID, task_str in enumerate(open(jobs_file, "r").readlines()):

            # get the outdir
            outdir_task = task_str.split("--outdir ")[1].split()[0].replace("/gpfs/projects/bsc40/mschikora", ParentDir)

            # define the exit stats
            perSVade_was_ran = os.path.isdir(outdir_task) and any([f not in {"aligned_reads.bam.sorted", "aligned_reads.bam.sorted.bai"} for f in os.listdir(outdir_task)])
            perSVade_finished = (not fun.file_is_empty("%s/perSVade_finished_file.txt"%(outdir_task)))

            dict_data[taskID] = {"taskID":taskID+1, "outdir_task":outdir_task, "perSVade_was_ran":perSVade_was_ran, "perSVade_finished":perSVade_finished, "dest_stdout":fun.get_file(dest_stdout_report), "jobs_file":fun.get_file(dest_jobs_file)}

        df_data = pd.DataFrame(dict_data).transpose()
        fun.save_df_as_tab(df_data, dest_jobEndingStatus_file)

        # at the end keep the report with all files
        fun.rsync_file(jobs_file, dest_jobs_file)
        fun.rsync_file(stderr_report, dest_stderr_report)
        fun.rsync_file(stdout_report, dest_stdout_report)
        
def get_goldenSetTable_Cglabrata(CurDir):

    """Generates a table with the golden set reads for C. glabrata"""

    # define file
    table_file = "%s/goldenSet_table_Cglabrata_ONTreads.tab"%CurDir

    print("getting C. glabrata golden set table")

    # init
    data_dict = {}

    # define the parent dir
    ParentDir = "%s/samba_bsc40"%(os.getenv("HOME")); # local
    if not os.path.exists(ParentDir): ParentDir = "/gpfs/projects/bsc40"

    # move the reads of Marina's assemblies
    marina_dir = "%s/current/mmarcet/nanopore/GABALDON02/assembly_files"%ParentDir
    dest_reads_dir = "%s/Cglabrata_reads"%CurDir; fun.make_folder(dest_reads_dir)
    for sampleID in ["B1012M", "BG2", "CST110", "EB0911Sto", "M12", "P35-2"]:

        # define origin files
        origin_nanoporeFile = "%s/%s/nanopore.pass.fastq.gz"%(marina_dir, sampleID)
        if fun.file_is_empty(origin_nanoporeFile): origin_nanoporeFile = "%s/%s/nanopore.reads.pass.fastq.gz"%(marina_dir, sampleID)
        if fun.file_is_empty(origin_nanoporeFile): origin_nanoporeFile = "%s/%s/nanopore.fastq.gz"%(marina_dir, sampleID)

        illumina_reads_dir = "%s/%s/illumina"%(marina_dir, sampleID)

        # define the origin illumina based on the reads
        if sampleID!="BG2":

            origin_reads1 = ["%s/%s"%(illumina_reads_dir, f) for f in os.listdir(illumina_reads_dir) if "_1." in f][0]
            origin_reads2 = ["%s/%s"%(illumina_reads_dir, f) for f in os.listdir(illumina_reads_dir) if "_2." in f][0]

        else:

            origin_reads1 = "%s/mschikora/Cglabrata_antifungals/data/raw_reads/RUN4_BG2_SRA_WT_R1.fq.gz"%ParentDir
            origin_reads2 = "%s/mschikora/Cglabrata_antifungals/data/raw_reads/RUN4_BG2_SRA_WT_R2.fq.gz"%ParentDir




        # dest files
        dest_nanopore_file = "%s/%s_nanopore.fastq.gz"%(dest_reads_dir, sampleID)
        dest_reads_1 = "%s/%s_illumina_1.fastq.gz"%(dest_reads_dir, sampleID)
        dest_reads_2 = "%s/%s_illumina_2.fastq.gz"%(dest_reads_dir, sampleID)

        # rsync reads only if missing
        if any([fun.file_is_empty(f) for f in [dest_nanopore_file, dest_reads_1, dest_reads_2]]): 

            if any([fun.file_is_empty(f) for f in [origin_nanoporeFile, origin_reads1, origin_reads2]]): 

                print(origin_nanoporeFile)
                print(origin_reads1)
                print(origin_reads2)

                raise ValueError("There origin files should be defined in %s"%sampleID)
            # rsync
            fun.rsync_file(origin_nanoporeFile, dest_nanopore_file)
            fun.rsync_file(origin_reads1, dest_reads_1)
            fun.rsync_file(origin_reads2, dest_reads_2)

        data_dict[sampleID] = {"sampleID":sampleID, "short_reads_1":dest_reads_1, "short_reads_2":dest_reads_2, "long_reads":dest_nanopore_file}

    df = pd.DataFrame(data_dict).transpose()
    fun.save_df_as_tab(df, table_file)

    return table_file

def get_df_benchmark_from_file(Idf, df_benchmarking_file, parms_row, test_row, parameters_df_metadata, test_df_metadata, all_keys):

    """Takes a df benchmarking file and returns the df"""

    print("df %i"%Idf)

    # load the df
    df_benchmark = fun.get_tab_as_df_or_empty_df(df_benchmarking_file) 

    # add metadata fields
    for x in parameters_df_metadata: df_benchmark["parms_%s"%x] = parms_row[x]
    for x in test_df_metadata: df_benchmark["test_%s"%x] = test_row[x]

    # keep
    return df_benchmark[all_keys]

def get_df_accuracy_of_parameters_on_test_samples(parameters_df, test_df, outdir, replace=False, threads=4, threads_integration=48, remove_SVs_overlapping_simple_repeats=False):

    """This function tests the accuracy of each of the parameters df on the simulations of test_df. It runs each comparison in a sepparate perSVade job in the cluster"""

    if replace is True: fun.delete_folder(outdir)
    fun.make_folder(outdir)

    # define the metadata of each df
    parameters_df_metadata = [k for k in parameters_df.keys() if k not in {"parameters_json"}]
    test_df_metadata = [k for k in test_df.keys() if k not in {"sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix", "interesting_svtypes"}]

    # define the outdir
    outdir_cross_benchmark_files = "%s/tmp_files"%outdir; fun.make_folder(outdir_cross_benchmark_files)

    # define the benchmarking file
    df_benchmark_all_file = "%s/cross_benchmarking_parameters.tab"%outdir

    if fun.file_is_empty(df_benchmark_all_file):
        print("running get_df_accuracy_of_parameters_on_test_samples")

        # initialize the df of the benchmarking
        benchmarking_fields = ['FN', 'FP', 'Fvalue', 'TP', 'nevents', 'precision', 'recall', 'svtype']
        all_keys = ["parms_%s"%x for x in parameters_df_metadata] + ["test_%s"%x for x in test_df_metadata] + benchmarking_fields

        list_dfs_benchmark_files = []

        # init cmds
        all_cmds = []

        # map each parameterID to the equivalent parameters
        parameters_df["parameters_json_dict"] = parameters_df.parameters_json.apply(fun.get_parameters_from_json)
        parmID_to_equal_parmIDs = {parmID : {other_parmID for other_parmID, r_other in parameters_df.iterrows() if r_other["parameters_json_dict"]==r["parameters_json_dict"]} for parmID, r in parameters_df.iterrows()}

        # define the parameter to run for each parameter
        parmID_to_parmIDtoRun = {pID : sorted(equal_pIDs)[0] for pID, equal_pIDs in parmID_to_equal_parmIDs.items()}

        # iterate through each set of parameters and testing
        Idf = 0
        for numeric_parameter_index, (Irow, parms_row) in enumerate(parameters_df.iterrows()):
            print(numeric_parameter_index, Irow, parmID_to_parmIDtoRun[Irow])

            # define the running parm from parmID_to_parmIDtoRun (avoid duplications)
            running_parmID = parmID_to_parmIDtoRun[Irow]
            unique_parms_row = parameters_df.loc[running_parmID]

            print(unique_parms_row)
            print(parameters_df_metadata)

            outdir_parms = "%s/parms_%s"%(outdir_cross_benchmark_files, "_".join(unique_parms_row[parameters_df_metadata]))
            fun.make_folder(outdir_parms)

            for numeric_test_index, (Itest, test_row) in enumerate(test_df.iterrows()):
            
                # define the final df
                df_benchmarking_file = "%s/testON_%s_benchmarking_df.tab"%(outdir_parms, "_".join(test_row[test_df_metadata]))

                # define cmd
                testing_script = "%s/scripts/perSVade/perSVade_repository/testing/get_accuracy_parameters_on_sorted_bam.py"%ParentDir

                if fun.file_is_empty(df_benchmarking_file): 

                    # define cmd
                    cmd_test = "%s --reference_genome %s --df_benchmarking_file %s --sorted_bam %s --gridss_vcf %s --svtables_prefix %s --threads %i --parameters_json %s --verbose --mitochondrial_chromosome %s"%(testing_script, test_row.reference_genome, df_benchmarking_file, test_row.sorted_bam, test_row.gridss_vcf, test_row.svtables_prefix, threads, unique_parms_row.parameters_json, test_row.mitochondrial_chromosome)

                    if "interesting_svtypes" in set(test_row.keys()): cmd_test += " --interesting_svtypes %s"%(test_row.interesting_svtypes)
                    if remove_SVs_overlapping_simple_repeats is True: cmd_test += " --remove_SVs_overlapping_simple_repeats"

                    # keep cmd and continue
                    all_cmds.append(cmd_test)
                    continue

                # append
                list_dfs_benchmark_files.append((Idf, df_benchmarking_file, parms_row, test_row, parameters_df_metadata, test_df_metadata, all_keys))
                Idf += 1

        # get unique cmds
        all_cmds_unique = sorted(set(all_cmds))
        print("There are %i jobs to run. These are %i unique jobs"%(len(all_cmds), len(all_cmds_unique)))

        # run cmds
        if run_in_cluster is False: # 'or True is debug' 
            for cmd in all_cmds_unique: 
                fun.run_cmd(cmd)
                raiseErrorAfterFirstTry

        # run in the cluster
        elif len(all_cmds_unique)>0:

            # get jobs file
            print("submitting %i cmds to the cluster"%len(all_cmds_unique))
            jobs_filename = "%s/jobs.getting_crossAccuracy"%outdir
            open(jobs_filename, "w").write("\n".join(all_cmds_unique))
            fun.generate_jobarray_file(jobs_filename, "gettingCloseShortReads")

            # submit to the cluster
            fun.run_jobarray_file_MN4_greasy(jobs_filename, "getting_crossAccuracy", time="02:00:00", queue="debug", threads_per_job=threads, nodes=8) # max 8 nodes (often ran on 6 nodes)
            #fun.run_jobarray_file_MN4_greasy(jobs_filename, "getting_crossAccuracy", time="10:00:00", queue="bsc_ls", threads_per_job=threads, nodes=3)

            # exit before it starts
            print("You need to run all the jobs from %s"%jobs_filename)
            sys.exit(0)

        # get the df for a function
        with multiproc.Pool(threads_integration) as pool:

            # run in parallel porechop runs for each chunk
            list_dfs_benchmark = pool.starmap(get_df_benchmark_from_file, list_dfs_benchmark_files)

            # close the pool
            pool.close()
            pool.terminate()

        # merge
        df_benchmark_all = pd.concat(list_dfs_benchmark)

        # save
        fun.save_df_as_tab(df_benchmark_all, df_benchmark_all_file)

    # clean
    fun.delete_folder(outdir_cross_benchmark_files)
    fun.delete_folder("%s/STDfiles"%outdir)

    # load
    df_benchmark_all = fun.get_tab_as_df_or_empty_df(df_benchmark_all_file)

    return df_benchmark_all

def get_accuracy_df_goldenSet(outdir_testing_GoldenSet):

    """This function loads into a single df the accuracy results of the golden set bencjmarking"""

    df_accuracy_all = pd.DataFrame()
    for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in species_Info:
        print(spName)

        # define outir
        outdir_species = "%s/%s_%s/testing_goldenSetAccuracy/perSVade_vs_longReads"%(outdir_testing_GoldenSet, taxID, spName)

        for file in os.listdir(outdir_species):
            

            sampleID = file.split(".tab")[0].split("_")[-1]
            df_accuracy = pd.read_csv("%s/%s"%(outdir_species, file), sep="\t")

            # keep
            df_accuracy["species"] = spName
            df_accuracy["sampleID"] = sampleID

            df_accuracy_all = df_accuracy_all.append(df_accuracy)

    return df_accuracy_all


def plot_goldenSet_accuracy_lineplots(df, fileprefix, accuracy_f="Fvalue", svtype="integrated"):

    """
    Plots the accuracy of the golden sets as lineplots

    """


    #youshouldremaketheplotssothattheyconsiderTheFactThatPrecision_and_recall_are_calculated_differently

    # filter get a df plot 
    df_all = df[(df.svtype==svtype) & (df.type_SVs_longReads=="all_SVs") & (df.remaining_treatment=="drop")]

    df_recall = df_all[(df_all.threshold_fractionParms==0.50) & (df_all.tol_bp==50) & (df_all.pct_overlap==0.75)][["recall", "species", "comparisonID", "sampleID"]]

    df_precision = df_all[(df_all.threshold_fractionParms==min(df_all.threshold_fractionParms)) & (df_all.tol_bp==50) & (df_all.pct_overlap==0.75)][["precision", "species", "comparisonID", "sampleID"]]

    #df_precision = df_all[(df_all.threshold_fractionParms==min(df_all.threshold_fractionParms)) & (df_all.tol_bp==500) & (df_all.pct_overlap==0.50)][["precision", "species", "comparisonID", "sampleID"]]


    df_plot = df_precision.merge(df_recall, on=["species", "comparisonID", "sampleID"], how="left", validate="one_to_one").set_index("species", drop=False)

    df_plot["Fvalue"] = (2*df_plot.precision*df_plot.recall) / (df_plot.precision+df_plot.recall)

    # define vars
    typeRun_to_color = {"perSVade-uniform":"blue", "perSVade-fast":"gray", "perSVade-realSVs":"red", 'perSVade-arroundHomRegions':"black"}
    sorted_typeRuns  = ["perSVade-arroundHomRegions", "perSVade-uniform", "perSVade-realSVs"]
    sorted_species =  ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana", "Drosophila_melanogaster"]

    typeRun_to_typeRunText = {"perSVade-realSVs":"known SVs", "perSVade-uniform":"random", "perSVade-arroundHomRegions":"homologous SVs"}

    # plot init
    nrows = len(sorted_species)
    ncols = len(sorted_typeRuns)
    fig = plt.figure(figsize=(ncols*1.3, nrows*1.3)); I = 1


    for Is, species in enumerate(sorted_species):
        print(species)

        # define the df
        df_species = df_plot[(df_plot.species==species)]

        # define the species
        ylim = [min(df_species[accuracy_f])-0.05, max(df_species[accuracy_f])+0.05]
        for Ir, typeRun in enumerate(sorted_typeRuns):

            ax = plt.subplot(nrows, ncols, I); I+=1

            # get the df
            df = df_species[(df_species.comparisonID.isin({typeRun, "perSVade-fast"}))]
            if len(df)!=len(set(df.sampleID))*2: raise ValueError("df is not correct")

            # add the xval
            compID_to_number = {"perSVade-fast":0, typeRun:1}
            df["x_value"] = df.comparisonID.map(compID_to_number) 
            df = df.sort_values(by="x_value")

            # plot 
            ax = sns.lineplot(data=df, x="x_value", y=accuracy_f, style="sampleID", color=typeRun_to_color[typeRun], markers=True, dashes=False)
            
            # adjust
            ax.get_legend().remove()    
            if Is==0: ax.set_title("%s"%(typeRun.split("-")[1].replace("arroundHomRegions", "HomRegions")))
            ax.set_ylim(ylim)
            ax.set_xlim([-0.45, 1.45])

            if Is==(nrows-1): 
                ax.set_xticks([0,1])
                ax.set_xticklabels(["default parms.", typeRun_to_typeRunText[typeRun]])
                for label in ax.get_xticklabels(): label.set_rotation(90)

            else: 
                ax.set_xticks([])

            ax.set_xlabel("")
            if Ir!=0:
                ax.set_yticklabels([])
                ax.set_ylabel("")

            else: ax.set_ylabel("%s. %s\n%s"%(species.split("_")[0][0], species.split("_")[1], accuracy_f))

    plt.subplots_adjust(wspace=0.15, hspace=0.15) 
    filename = "%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype)
    print("saving %s"%filename)
    plt.subplots_adjust(wspace=0.05, hspace=0.15) 

    fig.savefig(filename, bbox_inches='tight')


def plot_goldenSet_accuracy_barplots(df, fileprefix, accuracy_f="Fvalue", svtype="integrated"):

    """
    Plots the accuracy of the golden sets as barplots by defining true variants those that have tshd_parms.

    """

    # filter the df
    df_plot = df[(df.svtype==svtype) & (df.comparisonID!="perSVade-arroundRepeats") & (df.threshold_fractionParms==0.5310526315789474) & (df.type_SVs_longReads=="all_SVs") & (df.remaining_treatment=="drop") & (df.tol_bp==50) & (df.pct_overlap==0.75)].set_index("species", drop=False)

    # define vars
    typeRun_to_color = {"perSVade-uniform":"blue", "perSVade-fast":"gray", "perSVade-realSVs":"red", 'perSVade-arroundHomRegions':"black"}

    ylim = [0, max(df_plot[accuracy_f])+0.05]

    # plot init
    fig = plt.figure(figsize=(6*0.5, 4*2)) 


    # C. glabrata
    ax = plt.subplot2grid(shape=(6, 6), loc=(0, 0), colspan=6) # loc is (row, col) of the first, colspa
    sns.barplot(x="sampleID", y=accuracy_f, data=df_plot.loc["Candida_glabrata"], hue="comparisonID", palette=typeRun_to_color)
    ax.legend(bbox_to_anchor=(1, 1))
    ax.set_xticklabels([])
    ax.set_xlabel("")
    ax.set_ylim(ylim)
    ax.set_ylabel("")
    ax.set_title("%s for %s SVs"%(accuracy_f, svtype))

    # C. albicans
    ax = plt.subplot2grid(shape=(6, 6), loc=(1, 0), colspan=1) # loc is (row, col) of the first, colspa
    sns.barplot(x="sampleID", y=accuracy_f, data=df_plot.loc["Candida_albicans"], hue="comparisonID", palette=typeRun_to_color)
    ax.set_xlabel("")
    ax.set_xticklabels([])
    ax.get_legend().remove()
    ax.set_ylim(ylim)
    ax.set_ylabel("")

    # C. neoformans
    ax = plt.subplot2grid(shape=(6, 6), loc=(1, 1), colspan=5) # loc is (row, col) of the first, colspa
    sns.barplot(x="sampleID", y=accuracy_f, data=df_plot.loc["Cryptococcus_neoformans"], hue="comparisonID", palette=typeRun_to_color)
    ax.set_xlabel("")
    ax.set_xticklabels([])
    ax.get_legend().remove()
    ax.set_ylim(ylim)
    ax.set_ylabel("")
    ax.set_yticklabels([])


    # Drosophila
    ax = plt.subplot2grid(shape=(6, 6), loc=(2, 0), colspan=6) # loc is (row, col) of the first, colspa
    sns.barplot(x="sampleID", y=accuracy_f, data=df_plot.loc["Drosophila_melanogaster"], hue="comparisonID", palette=typeRun_to_color)
    ax.set_xlabel("")
    ax.set_xticklabels([])
    ax.get_legend().remove()    
    ax.set_ylim(ylim)
    ax.set_ylabel("")

    # A. thaliana
    ax = plt.subplot2grid(shape=(6, 6), loc=(3, 0), colspan=4) # loc is (row, col) of the first, colspa
    sns.barplot(x="sampleID", y=accuracy_f, data=df_plot.loc["Arabidopsis_thaliana"], hue="comparisonID", palette=typeRun_to_color)
    ax.set_xlabel("")
    ax.set_xticklabels([])
    ax.get_legend().remove()    
    ax.set_ylim(ylim)
    ax.set_ylabel("")


    #if I!=0: ax.get_legend().remove()

    filename = "%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype)
    print("saving %s"%filename)
    plt.subplots_adjust(wspace=0.2, hspace=0.15) 

    fig.savefig(filename, bbox_inches='tight')

def get_used_parameters_testing_several_species(outdir_testing, outdir_testing_human):

    """Takes the outdir testing and outputs a df with the used parameters."""

    parameters_dict = {}; I=0

    # map each species to a type of simulations and to outdir_species_simulations
    spName_to_typeSimulations_to_outdir_species_simulations = {spName : {typeSimulations : "%s/%s_%s/testing_Accuracy/%s"%(outdir_testing, taxID, spName, typeSimulations) for typeSimulations in ["arroundHomRegions", "uniform", "realSVs", "fast"]} for (taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads) in species_Info}

    spName_to_typeSimulations_to_outdir_species_simulations["Homo_sapiens"] = {typeSimulations : "%s/testing_Accuracy/%s"%(outdir_testing_human, typeSimulations) for typeSimulations in ["uniform", "realSVs", "fast"]}

    # create the used parameters df
    for spName, typeSimulations_to_outdir_species_simulations in spName_to_typeSimulations_to_outdir_species_simulations.items(): 
        for typeSimulations, outdir_species_simulations in typeSimulations_to_outdir_species_simulations.items():
                
                # define samples
                all_sampleIDs = [x for x in os.listdir(outdir_species_simulations)]

                # go through each sampleID
                for sampleID in all_sampleIDs:
                    print(spName, typeSimulations, sampleID)
                 
                    # define the parameters
                    parameters_json = "%s/%s/simulations_files_and_parameters/final_parameters.json"%(outdir_species_simulations, sampleID)

                    # load them
                    gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup = fun.get_parameters_from_json(parameters_json)

                    # modify the filters dict
                    gridss_filters_dict = {"GRIDSS %s"%filt_name : filt_val for filt_name, filt_val in gridss_filters_dict.items()}

                    # add them to dict
                    parameters_dict[I] = {**{"species":spName, "typeSimulations":typeSimulations, "sampleID":sampleID, "CLOVE max_rel_coverage_del":max_rel_coverage_to_consider_del, "CLOVE min_rel_coverage_dup":min_rel_coverage_to_consider_dup}, **gridss_filters_dict}; I+=1

    df_parameters = pd.DataFrame(parameters_dict).transpose()

    return df_parameters

def get_heatmaps_used_parameters(df_parameters, filename):

    """takes a df from get_used_parameters_testing_several_species and plots the used parameters as a heatmap. The rows will be the different types of parameters and the cols the different species, samples. The colormap will indicate how many samples in a given species/type simulation yields a given parameter """

    # sort df
    df_parameters = df_parameters.sort_values(by=["species", "typeSimulations"])
    df_parameters["species_simulations"] = df_parameters.species + "||||" + df_parameters.typeSimulations

    # filterout fast
    df_parameters = df_parameters[df_parameters.typeSimulations!="fast"]

    # define vars
    all_typesSimulations = ["arroundHomRegions", "uniform", "realSVs"] # fast
    all_sorted_species = ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana", "Drosophila_melanogaster", "Homo_sapiens"]

    species_to_color = {'none': 'gray', 'Drosophila_melanogaster': 'darkorange', 'Arabidopsis_thaliana': 'olive', 'Cryptococcus_neoformans': 'lightcoral', 'Candida_albicans': 'magenta', 'Candida_glabrata': 'lightseagreen', "Homo_sapiens":"brown"}

    typeSimulations_to_color = {"uniform":"blue", "realSVs":"red", "arroundHomRegions":"black", "fast":"gray"}

    cathegory_to_colors_dict = {"species" : species_to_color,
                            "typeSimulations" : typeSimulations_to_color}

    all_filter_names = [k for k in df_parameters.keys() if k.split()[0] in {"GRIDSS", "CLOVE"}]

    # init a list of weights
    yticklabels_weights = []

    # init df
    df_square_all = pd.DataFrame()

    # one plot for each filter
    for filter_name in all_filter_names:

        # get the values
        sorted_filter_vals = sorted(set(df_parameters[filter_name]))
        if len(sorted_filter_vals)==1: continue

        # get df square where the rows are the values of this filter_name, the columns are the
        def get_series_filter_val_to_n_samples(df_s): return pd.Series({filt_val : sum(df_s[filter_name]==filt_val) for filt_val in sorted_filter_vals})
        df_square = df_parameters.groupby("species_simulations").apply(get_series_filter_val_to_n_samples).transpose()

        # rename the index
        df_square.index = [str(x) for x in df_square.index]

        # add a blank
        df_square_all = df_square_all.append(pd.DataFrame({c : {"" : -1} for c in df_square.columns}))

        # add the label
        df_square_all = df_square_all.append(pd.DataFrame({c : {"%s"%filter_name : -1} for c in df_square.columns}))

        # append df_square
        df_square_all = df_square_all.append(df_square)

        yticklabels_weights += (["normal", "bold"] + ["normal"]*len(df_square))
     
    # sort by species
    species_to_order =  {'none': 0, "Candida_glabrata":1, "Candida_albicans":2, "Cryptococcus_neoformans":3, "Arabidopsis_thaliana":4, "Drosophila_melanogaster":5, "Homo_sapiens":6}
    typeSimulations_to_order = {"fast":0, "uniform":1, "realSVs":2, "arroundHomRegions":3}
    col_to_order = {c : (species_to_order[c.split("||||")[0]], typeSimulations_to_order[c.split("||||")[1]])  for c in df_square_all.columns}
    sorted_cols = sorted(df_square_all.columns, key=(lambda x: col_to_order[x]))
    df_square_all = df_square_all[sorted_cols]

    # generate the cols colors df
    def get_colors_series(idx):
        # type_keys can be parms or test

        # get the color dicts
        keys = ["species", "typeSimulations"]

        # get the content
        idx_content = idx.split("||||")

        # define the series
        field_to_color = {keys[I] : cathegory_to_colors_dict[keys[I]][c] for I,c in enumerate(idx_content)}

        return pd.Series(field_to_color)

    col_colors_df = pd.Series(df_square_all.columns, index=df_square_all.columns).apply(get_colors_series)

    # get colormap
    cm = sns.clustermap(df_square_all, col_cluster=False, row_cluster=False, row_colors=None, col_colors=col_colors_df, cbar_kws={'label': "# samples", "ticks":[0, 1, 2, 3]}, xticklabels=False, square=False, cmap=["white", "whitesmoke", "silver", "dimgray", "black"],  yticklabels=True, linewidth=0.1) 

    # set the weight
    for I, fw in enumerate(yticklabels_weights): cm.ax_heatmap.get_yticklabels()[I].set_weight(fw) 

    # set the ticks
    cm.ax_heatmap.tick_params(axis='both', which='both', length=0)


    # adjust the heatmap to fit the row colors
    cmap_pos = cm.ax_col_colors.get_position()
    size_box = cmap_pos.height/2
    new_hm_height = size_box*len(df_square_all)

    hm_pos = cm.ax_heatmap.get_position()
    cm.ax_heatmap.set_position([hm_pos.x0, cmap_pos.y0-new_hm_height-(size_box*0.0), hm_pos.width, new_hm_height])

    # do not rotate yaxis
    plt.setp(cm.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)

    # add the legend of the heatmap
    legend_elements = [("species", "white")] + [(s, species_to_color[s]) for s in all_sorted_species] + [("", "white"), ("typeSimulations", "white")] + [(s, typeSimulations_to_color[s]) for s in all_typesSimulations] 

    legend_elements = [Line2D([0], [0], color=color, lw=7, label=label, alpha=1.0) for label, color in legend_elements]
    cm.ax_heatmap.legend(handles=legend_elements, loc='right', bbox_to_anchor=(1.9, 0.8))

    # ste title
    cm.ax_col_colors.set_title("parameters chosen from simulations")

    print("saving %s"%filename)
    cm.savefig(filename, bbox_inches='tight')


def get_df_cross_accuracy_benchmark_withExtraFields(df_cross_accuracy_benchmark):

    """Takes a df_cross_accuracy_benchmark and adds some extra fields that will be useful for plotting the accuracies"""

    print("adding extra fields to df_cross_accuracy_benchmark")

    # add sample and run IDs (human samples have the same sample and run)
    human_samples = {"HG002run1", "CHMrun1", "NA12878run1"}
    def get_run(sampleID):
        
        if sampleID in human_samples: return sampleID
        else:
            if sampleID=="fast": return "fast"
            else: return sampleID.split("_")[1]

    def get_sample(sampleID): 

        if sampleID in human_samples: return sampleID
        else: return sampleID.split("_")[0]

    df_cross_accuracy_benchmark["parms_sample"] = df_cross_accuracy_benchmark.parms_sampleID.apply(get_sample)
    df_cross_accuracy_benchmark["parms_run"] = df_cross_accuracy_benchmark.parms_sampleID.apply(get_run)
    df_cross_accuracy_benchmark["test_sample"] = df_cross_accuracy_benchmark.test_sampleID.apply(get_sample)
    df_cross_accuracy_benchmark["test_run"] = df_cross_accuracy_benchmark.test_sampleID.apply(get_run)

    # add the type of comparison
    print("running get_type_comparison")
    def get_type_comparison(r):

        # trained on fast
        if r["parms_species"]=="none": return "fast"

        # trained on the same paramteres as tested
        elif r["parms_species"]==r["test_species"] and r["parms_sampleID"]==r["test_sampleID"] and r["parms_typeSimulations"]==r["test_typeSimulations"]: return "same_run_and_simulation"

        # trained on the same specues and simulation type
        elif r["parms_species"]==r["test_species"] and r["parms_typeSimulations"]==r["test_typeSimulations"]: return "same_species_and_simulation"

        # trained on the same species
        elif r["parms_species"]==r["test_species"]: return "same_species"

        # trained on other speceis
        elif r["parms_species"]!=r["test_species"]: return "different_species"

        else: raise ValueError("r is not valid")

    df_cross_accuracy_benchmark["type_comparison"] = df_cross_accuracy_benchmark.apply(get_type_comparison, axis=1)
    
    # add the intra-species sample and run ID
    species_to_sample_to_numericSample = {}
    species_to_sample_to_run_to_numericRun = {}

    df = df_cross_accuracy_benchmark[["parms_species", "parms_sample", "parms_run"]].drop_duplicates()
    for species in sorted(set(df.parms_species)):
        df_species = df[df.parms_species==species]

        for Is, sample in enumerate(sorted(set(df_species.parms_sample))):
            df_sample = df[df.parms_sample==sample]
            species_to_sample_to_numericSample.setdefault(species, {}).setdefault(sample, Is)

            for Ir, run in enumerate(sorted(set(df_sample.parms_run.values))): species_to_sample_to_run_to_numericRun.setdefault(species, {}).setdefault(sample, {}).setdefault(run, Ir)


    print("get_parms_numeric_sample")
    def get_parms_numeric_sample(r): return species_to_sample_to_numericSample[r.parms_species][r.parms_sample]
    df_cross_accuracy_benchmark["parms_numeric_sample"] = df_cross_accuracy_benchmark.apply(get_parms_numeric_sample, axis=1)

    print("get_test_numeric_sample")
    def get_test_numeric_sample(r): return species_to_sample_to_numericSample[r.test_species][r.test_sample]
    df_cross_accuracy_benchmark["test_numeric_sample"] = df_cross_accuracy_benchmark.apply(get_test_numeric_sample, axis=1)

    print("get_parms_numeric_run")
    def get_parms_numeric_run(r): return species_to_sample_to_run_to_numericRun[r.parms_species][r.parms_sample][r.parms_run]
    df_cross_accuracy_benchmark["parms_numeric_run"] = df_cross_accuracy_benchmark.apply(get_parms_numeric_run, axis=1)

    print("get_test_numeric_run")
    def get_test_numeric_run(r): return species_to_sample_to_run_to_numericRun[r.test_species][r.test_sample][r.test_run]
    df_cross_accuracy_benchmark["test_numeric_run"] = df_cross_accuracy_benchmark.apply(get_test_numeric_run, axis=1)

    return df_cross_accuracy_benchmark



def get_cross_accuracy_df_several_perSVadeSimulations(outdir_testing, outdir_testing_human, genomes_and_annotations_dir, replace=False):

    """This function tests how each of the perSVade configurations works on the others. It runs one job for each type of simulations, and it iterates through them inside of the job. It takes the runs on human hg38 as the human data set."""


    # define the final outdir
    outdir_cross_accuracy = "%s/cross_accuracy_calculations"%outdir_testing
    df_benchmark_all_file = "%s/cross_benchmarking_parameters.tab"%outdir_cross_accuracy

    if fun.file_is_empty(df_benchmark_all_file) or replace is True:

        ###### GET PARAMETERS DF AND TESTING DF #######

        # the parameters_df. The first cols are metadata (like sampleID, runID and optimisation type) and the others are things necessary for runnning gridss: and the path to the parameters_json
        parameters_df_dict = {}

        # test_df: This is info on which to test the running of gridss+clove. It contains metadata cols (sampleID, runID, optimisation type (real, uniform), simName, ploidy, svtype) and data to run the optimisation on (sorted_bam, gridss_vcf, reference_genome, mitochondrial_chromosome)
        test_df_dict = {}

        # map each species to a type of simulations and to outdir_species_simulations
        spName_to_typeSimulations_to_outdir_species_simulations = {spName : {typeSimulations : "%s/%s_%s/testing_Accuracy/%s"%(outdir_testing, taxID, spName, typeSimulations) for typeSimulations in ["arroundHomRegions", "uniform", "realSVs"]} for (taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads) in species_Info}

        spName_to_typeSimulations_to_outdir_species_simulations["Homo_sapiens"] = {typeSimulations : "%s/testing_Accuracy/%s"%(outdir_testing_human, typeSimulations) for typeSimulations in ["uniform", "realSVs"]}

        # create the used parameters df
        for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in species_Info_WithHumanHg38:
            for typeSimulations in ["arroundHomRegions", "uniform", "realSVs"]:

                # skip the human arroundHomRegions
                if spName=="Homo_sapiens" and typeSimulations=="arroundHomRegions": continue

                # define outir
                outdir_species_simulations = spName_to_typeSimulations_to_outdir_species_simulations[spName][typeSimulations]

                # define samples and runs
                all_sampleIDs = [x for x in os.listdir(outdir_species_simulations)]

                # go through each sampleID
                for sampleID in os.listdir(outdir_species_simulations):
                    print(spName, typeSimulations, sampleID)
                 
                    # define the outdir of the run
                    sampleID_simulations_files = "%s/%s/simulations_files_and_parameters"%(outdir_species_simulations, sampleID)
                    if not os.path.isdir(sampleID_simulations_files): raise ValueError("%s does not exist"%sampleID_simulations_files)

                    # keep parameters
                    parameters_json = "%s/final_parameters.json"%sampleID_simulations_files
                    parameters_df_dict[(spName, sampleID, typeSimulations)] = {"species":spName, "typeSimulations":typeSimulations, "sampleID":sampleID, "parameters_json":parameters_json}

                    # define simulation ploidies
                    if ploidy==1: simulation_ploidy = "haploid"
                    elif ploidy==2: simulation_ploidy = "diploid_hetero"

                    # go through each 
                    for simName in ["sim1", "sim2"]:

                        # define things
                        sorted_bam = "%s/reads_%s_%s.bam"%(sampleID_simulations_files, simName, simulation_ploidy)
                        gridss_vcf = "%s/gridss_vcf_%s_%s.vcf"%(sampleID_simulations_files, simName, simulation_ploidy)
                        svtables_prefix =  "%s/SVs_%s"%(sampleID_simulations_files, simName)

                        # define the reference genome (note that it has to be the hg38 for human)
                        if spName=="Homo_sapiens": reference_genome = "%s/../data/hg38.fa.corrected.fasta"%outdir_testing_human
                        else: reference_genome = "%s/%s.fasta"%(genomes_and_annotations_dir, spName)

                        # check that the genome is there
                        if fun.file_is_empty(reference_genome): raise ValueError("%s does not exist"%reference_genome)
                        

                        # add to dict
                        test_df_dict[(spName, sampleID, typeSimulations, simName, simulation_ploidy)] = {"species":spName, "sampleID":sampleID, "typeSimulations":typeSimulations, "simName":simName, "simulation_ploidy":simulation_ploidy, "sorted_bam":sorted_bam, "gridss_vcf":gridss_vcf, "reference_genome":reference_genome, "mitochondrial_chromosome":mitochondrial_chromosome, "svtables_prefix":svtables_prefix}


        # add the fast parameters
        parameters_json_fast = "%s/5478_Candida_glabrata/testing_Accuracy/fast/BG2_ANI/simulations_files_and_parameters/final_parameters.json"%outdir_testing
        parameters_df_dict[("none", "fast", "fast")] = {"species":"none", "typeSimulations":"fast", "sampleID":"fast", "parameters_json":parameters_json_fast}

        # get the dfs
        parameters_df = pd.DataFrame(parameters_df_dict).transpose()[["species", "sampleID", "typeSimulations", "parameters_json"]]
        test_df = pd.DataFrame(test_df_dict).transpose()[["species", "sampleID", "typeSimulations", "simName", "simulation_ploidy", "sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix"]]


        # plot the cross-accuracy
        print("getting cross-accuracy df")
        df_cross_accuracy_benchmark = get_df_accuracy_of_parameters_on_test_samples(parameters_df, test_df, outdir_cross_accuracy, replace=replace)

    df_cross_accuracy_benchmark = fun.get_tab_as_df_or_empty_df(df_benchmark_all_file)

    ##################################

    # add extra fields for plotting accuracies
    df_cross_accuracy_benchmark = get_df_cross_accuracy_benchmark_withExtraFields(df_cross_accuracy_benchmark)

    return df_cross_accuracy_benchmark

def get_parameters_df_cross_accuracy_df_realSVs(CurDir, humanSample_to_refGenomeID):

    """Gets a df with the parameters to be considered for the cross-accuracy analysis from the real df"""

    # init dict
    parameters_df_dict = {}

    for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in species_Info_WithHumanHg38:
        for typeSimulations in ["arroundHomRegions", "uniform", "realSVs"]:

            # skip the human arroundHomRegions
            if spName=="Homo_sapiens" and typeSimulations=="arroundHomRegions": continue

            # for non human species where there was a ON-based definition of 'golden set SVs'
            if spName!="Homo_sapiens":

                # define a dir where all the files are
                outdir_perSVadeCalling = "%s/outdirs_testing_severalSpecies_goldenSet/%s_%s/testing_goldenSetAccuracy/perSVade_SV_calling"%(CurDir, taxID, spName)

                # iterate through the contents of the dir that are related to typeSimulations
                for f in os.listdir(outdir_perSVadeCalling):

                    # continue if it is not this typeSimulations
                    if not f.startswith("perSVade_calling_%s_"%typeSimulations): continue

                    # define the parameters and sampleID
                    sampleID = "_".join(f.split("_")[3:])
                    parameters_json = "%s/%s/simulations_files_and_parameters/final_parameters.json"%(outdir_perSVadeCalling, f)
                    if fun.file_is_empty(parameters_json): raise ValueError("%s should exist"%parameters_json)

                    # add to the dict
                    parameters_df_dict[(spName, sampleID, typeSimulations)] = {"species":spName, "typeSimulations":typeSimulations, "sampleID":sampleID, "parameters_json":parameters_json}

            # for human samples, where the directory structure is a bit different
            else:

                for sampleID, refGenomeID in humanSample_to_refGenomeID.items():

                    # define the parameters
                    parameters_json = "%s/outdirs_testing_humanGoldenSet/running_on_%s/testing_Accuracy/%s/%s/simulations_files_and_parameters/final_parameters.json"%(CurDir, refGenomeID, typeSimulations, sampleID)
                    if fun.file_is_empty(parameters_json): raise ValueError("%s should exist"%parameters_json)

                    # add to the dict
                    parameters_df_dict[(spName, sampleID, typeSimulations)] = {"species":spName, "typeSimulations":typeSimulations, "sampleID":sampleID, "parameters_json":parameters_json}

    # add the fast parameters
    parameters_json_fast = "%s/outdirs_testing_severalSpecies/5478_Candida_glabrata/testing_Accuracy/fast/BG2_ANI/simulations_files_and_parameters/final_parameters.json"%CurDir
    if fun.file_is_empty(parameters_json): raise ValueError("%s should exist"%parameters_json)

    parameters_df_dict[("none", "fast", "fast")] = {"species":"none", "typeSimulations":"fast", "sampleID":"fast", "parameters_json":parameters_json_fast}

    parameters_df = pd.DataFrame(parameters_df_dict).transpose()[["species", "sampleID", "typeSimulations", "parameters_json"]]

    return parameters_df

def generate_highConfidence_and_all_SVs_files_ONcalling_severalParameter_combinations(dict_paths, outdir, reference_genome, replace, threads, run_in_parallel=True):

    """This function generates the SV files from different combinations of filtering the ON SV calls. It takes parts of fun.get_df_accuracy_perSVade_vs_longReads_one_sample. It will generate a folder with high confidence and one with low confidence thresholds"""

    # define final files
    final_file = "%s/ONbased_SVfiles_generated.txt"%outdir
    fun.make_folder(outdir)

    if fun.file_is_empty(final_file) or replace is True:
        print("generating SV files from ON reads")

        # define parms
        min_svim_QUAL = 2

        # get the svim and sniffles dfs
        svim_df = fun.get_svim_as_df(dict_paths["svim_outdir"], reference_genome, min_svim_QUAL)
        sniffles_df = fun.get_sniffles_as_df(dict_paths["sniffles_outdir"], reference_genome)

        from collections import Counter
        print("svim GT:", Counter(svim_df.GT))
        print("sniffles GT:", Counter(sniffles_df.GT))

        # define SVIM and SNIFFLES filters (same number each)
        all_min_QUAL_svim = [3, 5, 7, 10, 12, 15, 20, 30, 40, 50] # max of 100
        all_min_RE_sniffles = np.linspace(min(sniffles_df.INFO_RE), np.percentile(sniffles_df.INFO_RE, 90), len(all_min_QUAL_svim))

        # define necessary things for clove runs
        sorted_bam_longReads =  "%s/aligned_reads.sorted.bam"%dict_paths["svim_outdir"]

        # define general parameters
        filter_IMPRECISE_sniffles = True # remove the SVs with a tag of 'IMPRECISE'
        remaining_treatment = "drop" # do not consider the 'remaining' SVs
        tol_bp = 50 # standard to define the overlaps
        pct_overlap = 0.75 # standard to define the overlaps

        # define the type_SVs_longReads based on the species
        species = "_".join(reference_genome.split("/")[-2].split("_")[0:2])
        if species in {"Candida_glabrata", "Cryptococcus_neoformans"}: type_SVs_longReads = "haploid_SVs"
        elif species in {"Candida_albicans", "Arabidopsis_thaliana", "Drosophila_melanogaster"}: type_SVs_longReads = "all_SVs"
        else: raise ValueError("%s is not valid"%species)

        # define an outdir for this type of overlaps
        outdir_SVIMandSNIFFLEScalling = "%s/svim_and_snifflesCalling_%s_%s_%s_%s"%(outdir, type_SVs_longReads, filter_IMPRECISE_sniffles, remaining_treatment, tol_bp); fun.make_folder(outdir_SVIMandSNIFFLEScalling)

        # define the inputs of 
        svim_inputs = [(min_QUAL_svim, tol_bp, "svim", "QUAL", svim_df, outdir_SVIMandSNIFFLEScalling, sorted_bam_longReads, reference_genome, type_SVs_longReads, filter_IMPRECISE_sniffles, remaining_treatment) for min_QUAL_svim in all_min_QUAL_svim]

        sniffles_inputs = [(min_RE_sniffles, tol_bp, "sniffles", "INFO_RE", sniffles_df, outdir_SVIMandSNIFFLEScalling, sorted_bam_longReads, reference_genome, type_SVs_longReads, filter_IMPRECISE_sniffles, remaining_treatment) for min_RE_sniffles in all_min_RE_sniffles]

        inputs_fn = svim_inputs + sniffles_inputs
        inputs_fn = [tuple(list(x)+[I]) for I,x in enumerate(inputs_fn)]

        # run fn
        print("getting SVIM and SNIFFLES sv_to_svdf for %i filters and %i threads"%(len(inputs_fn), threads))
        
        if run_in_parallel is False: list_tuples_svtype_to_svDF_ON = list(map(lambda x: fun.get_svtype_to_svDF_withFiltering(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]), inputs_fn))

        else:
            with multiproc.Pool(threads) as pool:
                list_tuples_svtype_to_svDF_ON = pool.starmap(fun.get_svtype_to_svDF_withFiltering, inputs_fn)
                
                pool.close()
                pool.terminate()

        # get the commonID between svim and sniffles
        ON_svtype_to_svDF = fun.process_list_tuples_svtype_to_svDF_ON_to_add_commonID(cp.deepcopy(list_tuples_svtype_to_svDF_ON), outdir_SVIMandSNIFFLEScalling,  tol_bp, pct_overlap)

        # generate a folder with the high-confidence and all SVs 
        print("writing SVs")
        for threshold_fractionParms, typeSVs in [(0.5, "high_confidence_ONvars"), (0.0, "all_ONvars")]:

            # get the filtered dict
            ON_svtype_to_svDF_filt = {svtype : svDF[svDF.fraction_filters_withSV>=threshold_fractionParms] for svtype, svDF in ON_svtype_to_svDF.items()}

            # write each of the SVs
            outdir_SVs = "%s/%s"%(outdir, typeSVs); fun.make_folder(outdir_SVs)
            for svtype, svDF in ON_svtype_to_svDF_filt.items(): fun.save_df_as_tab(svDF, "%s/SVs_%s.tab"%(outdir_SVs, svtype))

        # at the end generate the file
        open(final_file, "w").write("ONreads generated")


def get_df_benchmark_cross_accuracy_ONbasedGoldenSet(CurDir, parameters_df, outdir_cross_accuracy_ONbased, threads, replace):

    """This function runs each of the parameters_df on two subsets of SVs (high-confidence for recall and all SVs for precision) for all other samples. It writes the results under outdir_cross_accuracy_ONbased. 

    Note that the test_df should end with [["species", "sampleID", "sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix", "interesting_svtypes"]]"""

    # define outdirs
    fun.make_folder(outdir_cross_accuracy_ONbased)
    df_benchmark_ONbased_file = "%s/df_benchmark_ONbased.tab"%outdir_cross_accuracy_ONbased

    if fun.file_is_empty(df_benchmark_ONbased_file) or replace is True:
        print("getting df-cross-accuracy ONbased on %i threads"%threads)

        ######## GENERATE ALL THE ON-BASED 'REAL' SV CALLS AS DONE IN THE GOLDEN-SET TESTING IMPLEMENTED IN PERSVADE ########

        # init the testing dfs dict
        typeONvars_to_test_df_dict = {}

        # go through each species and sample
        for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in species_Info:

            # go through all the samples of the golden set testing
            for sampleID in sorted(set(parameters_df[parameters_df.species==spName].sampleID)):
                print(spName, sampleID)

                # define the reference genome
                outdir_sampleID = "%s/%s_%s"%(outdir_cross_accuracy_ONbased, spName, sampleID); fun.make_folder(outdir_sampleID)
                origin_reference_genome = "%s/genomes_and_annotations/%s.fasta"%(CurDir, spName)
                reference_genome = "%s/reference_genome.fasta"%outdir_sampleID
                fun.soft_link_files(origin_reference_genome, reference_genome)

                # generate all the ON-based SVs files as done in fun.get_df_accuracy_perSVade_vs_longReads_one_sample
                goldenSetDir = "%s/outdirs_testing_severalSpecies_goldenSet/%s_%s/testing_goldenSetAccuracy"%(CurDir, taxID, spName)
                dict_paths = {"svim_outdir":"%s/ONT_SV_calling/%s/svim_output"%(goldenSetDir, sampleID),
                              "sniffles_outdir":"%s/ONT_SV_calling/%s/sniffles_output"%(goldenSetDir, sampleID)}

                generate_highConfidence_and_all_SVs_files_ONcalling_severalParameter_combinations(dict_paths, outdir_sampleID, reference_genome, replace, threads)

                # define files for the test_df. Note that we take one single run from perSVade
                perSVade_outdir = "%s/perSVade_SV_calling/perSVade_calling_uniform_%s"%(goldenSetDir, sampleID)
                sorted_bam = "%s/aligned_reads.bam.sorted"%perSVade_outdir
                gridss_vcf = "%s/SVdetection_output/final_gridss_running/gridss_output.raw.vcf"%perSVade_outdir # this already has the withSimpleEventType

                # obtain the CollectInsertSizeMetrics file
                fun.get_insert_size_distribution(sorted_bam, replace=replace, threads=threads) 

                # obtain the coverage per window
                destination_dir = "%s.calculating_windowcoverage"%sorted_bam
                min_chromosome_len = 100000
                fun.window_l = fun.get_perSVade_window_l(reference_genome, mitochondrial_chromosome, min_chromosome_len)
                coverage_file = fun.generate_coverage_per_window_file_parallel(origin_reference_genome, destination_dir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=True, delete_bams=True, threads=threads) # get the per window coverage

                # fill the test_df_dict:
                for typeSVs in ["high_confidence_ONvars", "all_ONvars"]:

                    # define the dictionary with all the fields for the test df
                    #perSVade_outdir
                    dict_data = {"species":spName, "sampleID":sampleID, "sorted_bam":sorted_bam, "gridss_vcf":gridss_vcf, "reference_genome":origin_reference_genome, "mitochondrial_chromosome":mitochondrial_chromosome, "svtables_prefix": "%s/%s/SVs"%(outdir_sampleID, typeSVs), "interesting_svtypes":"insertions,deletions,translocations,inversions,tandemDuplications"}

                    # define the files
                    typeONvars_to_test_df_dict.setdefault(typeSVs, {}).setdefault((spName, sampleID), dict_data)

        #####################################################################################################################


        ######## GET THE ACCURACY DF FOR EACH TYPE OF ON SVs ###########

        # init a dict that maps each type of ON vars to the cross-accuracy df
        typeONvars_to_crossAccuracy_df = {}

        for typeONvars, test_df_dict in typeONvars_to_test_df_dict.items():
            print("running get_df_accuracy_of_parameters_on_test_samples on %s"%typeONvars)

            # define the test_df
            test_df = pd.DataFrame(test_df_dict).transpose()[["species", "sampleID", "sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix", "interesting_svtypes"]]

            # define the cross-accuracy benchmark
            outdir_cross_accuracy_typeONvars = "%s/cross_accuracy_%s"%(outdir_cross_accuracy_ONbased, typeONvars)
            typeONvars_to_crossAccuracy_df[typeONvars] = get_df_accuracy_of_parameters_on_test_samples(parameters_df, test_df, outdir_cross_accuracy_typeONvars, replace=replace, remove_SVs_overlapping_simple_repeats=True)

        # merge the two
        print("integrating and saving")
        merge_fields = ["parms_species", "parms_sampleID", "parms_typeSimulations", "test_species", "test_sampleID", "svtype"]
        df_benchmark_ONbased = typeONvars_to_crossAccuracy_df["high_confidence_ONvars"].merge(typeONvars_to_crossAccuracy_df["all_ONvars"], on=merge_fields, validate="one_to_one", suffixes=("_highConf", "_allVars"))

        # add fields, derived from different sets of SVs
        df_benchmark_ONbased["recall"] = df_benchmark_ONbased.recall_highConf # the recall relates to the high confidence
        df_benchmark_ONbased["precision"] = df_benchmark_ONbased.precision_allVars # precision is calculated on all vars

        def calculate_Fvalue_from_r(r):

            if r.precision<=0.0 or r.recall<=0.0: Fvalue = 0.0
            else: Fvalue = (2*r.precision*r.recall)/(r.precision+r.recall)

            return Fvalue

        df_benchmark_ONbased["Fvalue"] = df_benchmark_ONbased.apply(calculate_Fvalue_from_r, axis=1)

        if any(pd.isna(df_benchmark_ONbased.Fvalue)): raise ValueError("there are nans in the Fvalue")

        ################################################################

        # save
        fun.save_df_as_tab(df_benchmark_ONbased, df_benchmark_ONbased_file)


    df_benchmark_ONbased = fun.get_tab_as_df_or_empty_df(df_benchmark_ONbased_file)

    return df_benchmark_ONbased

def get_svtables_prefix_human_goldenSets(human_goldenSet_dir, knownSVs_dir, sampleID, reference_genome):

    """This function generates under knownSVs_dir the files with known SVs (correctly formated), returning the prefix and the interesting_svtypes"""

    print("getting get_svtables_prefix_human_goldenSets for %s"%sampleID)

    ############ GENERATING SV FILES ################

    # define all the chromosomes
    all_chromosomes = set(fun.get_chr_to_len(reference_genome))

    # define the prefix
    svtables_prefix = "%s/SVs"%knownSVs_dir

    # NA12878 has only deletions. In the NatComm review they ignore deletions close to ENCODE DAC list (not done here). The downloaded file needs some refactoring to adjust the chromosomes
    if sampleID=="NA12878run1":

        # get the deletions df and modify fields
        deletions_df = fun.get_tab_as_df_or_empty_df("%s/data/NA12878_deletions.bed"%human_goldenSet_dir)
        deletions_df["Chr"] = ["chr%s"%x for x in deletions_df.Chr]
        deletions_df["ID"] = ["deletion%i"%I for I in range(len(deletions_df))]

        # keep only deletions >50bp
        deletions_df["svlen"] = deletions_df.End - deletions_df.Start
        deletions_df = deletions_df[deletions_df.svlen>=50]

        # write as deletions
        fun.save_df_as_tab(deletions_df[fun.svtype_to_fieldsDict["deletions"]["all_fields"]], "%s_deletions.tab"%svtables_prefix)

    # HG002run1 comes from a GiAB project. All kinds of SVs (They exclude variants that fall in the high-confidenceTier1regionsdefined region). We will only use variants in th  high-confidence set. The perSVade output can only be adapated to the deletions, as the HG002run1 DUP mainly show the sequence, but it is not clear if it is a tandem duplication

    elif sampleID=="HG002run1":

        # get the vcf loaded
        vcf_df = fun.get_vcf_df_with_INFO_as_single_fields(fun.get_df_and_header_from_vcf("%s/data/HG002_SVs_Tier1_v0.6.vcf"%human_goldenSet_dir)[0])

        # get only the deletions
        deletions_df = vcf_df[(vcf_df.INFO_SVTYPE=="DEL") & (vcf_df.INFO_REPTYPE.isin({"SIMPLEDEL", "SUBSDEL"}))]

        # calculate the len
        deletions_df["REF_len"] = deletions_df.REF.apply(len)
        deletions_df["ALT_len"] = deletions_df.ALT.apply(len)
        if any(deletions_df.ALT_len>=deletions_df.REF_len): raise ValueError("The ALT seq should be shorter")

        deletions_df["svlen"] = deletions_df.REF_len

        # keep only deletions >50bp
        deletions_df = deletions_df[deletions_df.svlen>=50]

        # add fields
        deletions_df["Start"] = deletions_df.POS-1
        deletions_df["End"] = deletions_df.Start + deletions_df.svlen
        deletions_df["Chr"] = ["chr%s"%x for x in deletions_df["#CHROM"]]

        # write as deletions
        fun.save_df_as_tab(deletions_df[fun.svtype_to_fieldsDict["deletions"]["all_fields"]], "%s_deletions.tab"%svtables_prefix)

    # CHMrun1 correspond to two syntetically merged haploid cell lines. The union of all variants forund in CHM1 and CHM13 is taken as the through set. I will only keep deletions and inversions

    elif sampleID=="CHMrun1":

        svtypes_CHM = ["inversions", "deletions"]
        positions_fields = ["Chr", "Start", "End"]
        svtype_to_svDF = {svtype : pd.DataFrame(columns=positions_fields) for svtype in svtypes_CHM}

        # go through each of the diploids
        for chmID in ["CHM1", "CHM13"]:

            # load the vcf 
            vcf_df = fun.get_vcf_df_with_INFO_as_single_fields(fun.get_df_and_header_from_vcf("%s/data/%s_SVs.annotated.vcf"%(human_goldenSet_dir, chmID))[0])

            # remove variants involving strange chroms
            vcf_df = vcf_df[~(vcf_df["#CHROM"].apply(lambda c: c.endswith("_random"))) & ~(vcf_df["#CHROM"].apply(lambda c: c.startswith("chrUn_")))]

            # remove variants whose length is below 50
            vcf_df = vcf_df[vcf_df.INFO_SVLEN>=50]

            # go through each svtype and get the svDF
            for svtype in svtypes_CHM:

                svDF = vcf_df[vcf_df.INFO_SVTYPE==(svtype[0:-1])]
                svDF["Chr"] = svDF["#CHROM"]
                svDF["Start"] = svDF.POS-1
                svDF["End"] = svDF.Start + svDF.INFO_SVLEN

                svtype_to_svDF[svtype] = svtype_to_svDF[svtype].append(svDF[positions_fields])

        # write the concatenated SVs
        for svtype, svDF in svtype_to_svDF.items():

            # get the filtered svDF with ID
            svDF_filt = svDF.sort_values(by=positions_fields).drop_duplicates()
            svDF_filt["ID"] = ["%s%i"%(svtype, I) for I in range(len(svDF_filt))]

            # write
            fun.save_df_as_tab(svDF_filt[fun.svtype_to_fieldsDict[svtype]["all_fields"]], "%s_%s.tab"%(svtables_prefix, svtype))

    else: raise ValueError("%s was not taken into consideration"%sampleID)

    #################################################

    ################ DEBUGGING AND SAVING ##############

    # define the interesting svtypes and check integrity of files
    interesting_svtypes_list = []

    for svtype, fieldsDict in fun.svtype_to_fieldsDict.items():

        # get the svfile (if existing)
        svfile = "%s_%s.tab"%(svtables_prefix, svtype)
        if fun.file_is_empty(svfile): continue

        # load df
        svDF = fun.get_tab_as_df_or_empty_df(svfile)

        # keep
        interesting_svtypes_list.append(svtype)

        # check that the chromosomes are correct
        for chr_f in fieldsDict["chromosome_fields"]:

            strange_chromosomes = set(svDF[chr_f]).difference(all_chromosomes)
            if len(strange_chromosomes)>0: raise ValueError("There are some strange chromosomes in the deletions_df: %s"%strange_chromosomes)

        # check that the fileds are correct
        if list(svDF.keys())!=fieldsDict["all_fields"]: raise ValueError("The fields %s are not correct for %s"%(svDF.keys(), svtype))

    # debug
    if len(interesting_svtypes_list)==0: raise("There should be some SVs")

    # join
    interesting_svtypes = ",".join(interesting_svtypes_list)

    ####################################################

    # return
    return svtables_prefix, interesting_svtypes

def get_df_benchmark_cross_accuracy_humanGoldenSet(CurDir, parameters_df, outdir_cross_accuracy_human, humanSample_to_refGenomeID, threads, replace):

    """This function runs each of the parameters_df on the human dataset of cross-accuracy. This is analogous to  get_df_benchmark_cross_accuracy_ONbasedGoldenSet.

    Note that the test_df should end with [["species", "sampleID", "sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix", "interesting_svtypes"]]"""

    # define outdirs
    fun.make_folder(outdir_cross_accuracy_human)
    df_benchmark_human_file = "%s/df_benchmark_human.tab"%outdir_cross_accuracy_human

    if fun.file_is_empty(df_benchmark_human_file) or replace is True:
        print("getting df-cross-accuracy human on %i threads"%threads)

        # define general things
        human_goldenSet_dir = "%s/outdirs_testing_humanGoldenSet"%CurDir
        refGenomeID_to_mtChromosome = {"hg38":"chrM", "hg19":"chrMT"}

        # init dict with the test data
        test_df_dict = {}

        for sampleID, refGenomeID in humanSample_to_refGenomeID.items():

            # define the general things
            perSVade_outdir = "%s/running_on_%s/testing_Accuracy/uniform/%s"%(human_goldenSet_dir, refGenomeID, sampleID)
            sorted_bam = "%s/aligned_reads.bam.sorted"%perSVade_outdir
            gridss_vcf = "%s/SVdetection_output/final_gridss_running/gridss_output.raw.vcf"%perSVade_outdir
            reference_genome = "%s/data/%s.fa.corrected.fasta"%(human_goldenSet_dir, refGenomeID)
            mitochondrial_chromosome = refGenomeID_to_mtChromosome[refGenomeID]

            # obtain the CollectInsertSizeMetrics file
            fun.get_insert_size_distribution(sorted_bam, replace=replace, threads=threads) 

            # obtain the coverage per window
            destination_dir = "%s.calculating_windowcoverage"%sorted_bam
            min_chromosome_len = 100000
            fun.window_l = fun.get_perSVade_window_l(reference_genome, mitochondrial_chromosome, min_chromosome_len)
            coverage_file = fun.generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=True, delete_bams=True, threads=threads) # get the per window coverage

            # note that the reference genomes have chr1 ... chrX ... chrM/chrMT, which needs to be adapted in the sv files

            # make a knownSVs dir
            knownSVs_dir = "%s/knownSVs_%s"%(outdir_cross_accuracy_human, sampleID); fun.make_folder(knownSVs_dir)

            # get the known SV files
            svtables_prefix, interesting_svtypes = get_svtables_prefix_human_goldenSets(human_goldenSet_dir, knownSVs_dir, sampleID, reference_genome)      

            # add to the test_df_dict
            test_df_dict[sampleID] = {"species":"Homo_sapiens", "sampleID":sampleID, "sorted_bam":sorted_bam, "gridss_vcf":gridss_vcf, "reference_genome":reference_genome, "mitochondrial_chromosome":mitochondrial_chromosome, "svtables_prefix":svtables_prefix, "interesting_svtypes":interesting_svtypes}

        # run the cross accuracy
        test_df = pd.DataFrame(test_df_dict).transpose()[["species", "sampleID", "sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix", "interesting_svtypes"]]

        outdir_cross_accuracy_human_generatingFiles = "%s/cross_accuracy_files"%outdir_cross_accuracy_human
        df_benchmark_human = get_df_accuracy_of_parameters_on_test_samples(parameters_df, test_df, outdir_cross_accuracy_human_generatingFiles, replace=replace, threads=8, remove_SVs_overlapping_simple_repeats=True)

        # save
        fun.save_df_as_tab(df_benchmark_human, df_benchmark_human_file)

    # return
    df_benchmark_human = fun.get_tab_as_df_or_empty_df(df_benchmark_human_file)
    return df_benchmark_human

def get_cross_accuracy_df_realSVs(CurDir, ProcessedDataDir, threads=4, replace=False):

    """This function is similar to  get_cross_accuracy_df_several_perSVadeSimulations, but working with the samples used for the real SV testing. Note that for some human datasets we only have some SVs. In addition, for the ON 'golden set testing' I have the precision and the recall calculated from a different set of ON-based SVs, which should be also considered."""

    # define outdirs
    outdir_cross_accuracy = "%s/cross_accuracy_calculations_realSVs"%CurDir; fun.make_folder(outdir_cross_accuracy)
    df_benchmark_all_file = "%s/cross_benchmarking_parameters.tab"%outdir_cross_accuracy

    if fun.file_is_empty(df_benchmark_all_file) or replace is True:

        # map each human dataset to the reference genome
        humanSample_to_refGenomeID = {"NA12878run1":"hg19", "HG002run1":"hg19", "CHMrun1":"hg38"}

        # get the training parameters df (only one per each sample, and taking only the corresponding reference genome from the human datasets)
        parameters_df = get_parameters_df_cross_accuracy_df_realSVs(CurDir, humanSample_to_refGenomeID)

        # get the df_bechmark for the human golden set
        outdir_cross_accuracy_human = "%s/cross_accuracy_human"%outdir_cross_accuracy
        df_benchmark_human = get_df_benchmark_cross_accuracy_humanGoldenSet(CurDir, parameters_df, outdir_cross_accuracy_human, humanSample_to_refGenomeID, threads, replace)

        # get the df_benchmark for the ON-based golden set SVs (non-human species)
        outdir_cross_accuracy_ONbased = "%s/cross_accuracy_ONbased"%outdir_cross_accuracy
        df_benchmark_ONbased = get_df_benchmark_cross_accuracy_ONbasedGoldenSet(CurDir, parameters_df, outdir_cross_accuracy_ONbased, threads, replace)

        # merge and save
        df_benchmark_all = df_benchmark_human.append(df_benchmark_ONbased)
        fun.save_df_as_tab(df_benchmark_all, df_benchmark_all_file)

    # load
    df_benchmark_all = fun.get_tab_as_df_or_empty_df(df_benchmark_all_file)

    ########### ADD METADATA #############

    # add the type of comparison
    print("running get_type_comparison")
    def get_type_comparison(r):

        # trained on fast
        if r["parms_species"]=="none": return "fast"

        # trained on the same sample as tested
        elif r["parms_species"]==r["test_species"] and r["parms_sampleID"]==r["test_sampleID"]: return "same_sample"

        # trained on the same species
        elif r["parms_species"]==r["test_species"]: return "same_species"

        # trained on other speceis
        elif r["parms_species"]!=r["test_species"]: return "different_species"

        else: raise ValueError("r is not valid")

    df_benchmark_all["type_comparison"] = df_benchmark_all.apply(get_type_comparison, axis=1)

    # keep df so that there are


    ######################################

    return df_benchmark_all


def get_cross_accuracy_df_realSVs_onlyHuman(CurDir, ProcessedDataDir, threads=4, replace=False):

    """This function is similar to  get_cross_accuracy_df_realSVs, but working only on the human SVs and optimizing on the simulations' parameters."""

    # define outdirs
    outdir_testing = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies"%ParentDir
    outdir_cross_accuracy = "%s/cross_accuracy_calculations_realSVs_onlyHuman_trainOnSimulationsDatasets"%CurDir; fun.make_folder(outdir_cross_accuracy)
    df_benchmark_all_file = "%s/cross_benchmarking_parameters.tab"%outdir_cross_accuracy

    if fun.file_is_empty(df_benchmark_all_file) or replace is True:

        ###### GET PARAMETERS DF OF NON-HUMAN #######

        # the parameters_df. The first cols are metadata (like sampleID, runID and optimisation type) and the others are things necessary for runnning gridss: and the path to the parameters_json
        parameters_df_dict = {}

        # map each species to a type of simulations and to outdir_species_simulations
        spName_to_typeSimulations_to_outdir_species_simulations = {spName : {typeSimulations : "%s/%s_%s/testing_Accuracy/%s"%(outdir_testing, taxID, spName, typeSimulations) for typeSimulations in ["arroundHomRegions", "uniform", "realSVs"]} for (taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads) in species_Info}

        # create the used parameters df
        for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in species_Info:
            for typeSimulations in ["arroundHomRegions", "uniform", "realSVs"]:

                # define outir
                outdir_species_simulations = spName_to_typeSimulations_to_outdir_species_simulations[spName][typeSimulations]

                # go through each sampleID
                for sampleID in os.listdir(outdir_species_simulations):
                    print(spName, typeSimulations, sampleID)
                    if sampleID.startswith("."): continue
                 
                    # define the outdir of the run
                    sampleID_simulations_files = "%s/%s/simulations_files_and_parameters"%(outdir_species_simulations, sampleID)
                    if not os.path.isdir(sampleID_simulations_files): raise ValueError("%s does not exist"%sampleID_simulations_files)

                    # keep parameters
                    parameters_json = "%s/final_parameters.json"%sampleID_simulations_files
                    parameters_df_dict[(spName, sampleID, typeSimulations)] = {"species":spName, "typeSimulations":typeSimulations, "sampleID":sampleID, "parameters_json":parameters_json}

        # add the fast parameters
        parameters_json_fast = "%s/5478_Candida_glabrata/testing_Accuracy/fast/BG2_ANI/simulations_files_and_parameters/final_parameters.json"%outdir_testing
        parameters_df_dict[("none", "fast", "fast")] = {"species":"none", "typeSimulations":"fast", "sampleID":"fast", "parameters_json":parameters_json_fast}

        # get the dfs
        parameters_df_notHuman = pd.DataFrame(parameters_df_dict).transpose()[["species", "sampleID", "typeSimulations", "parameters_json"]]

        #############################################

        ####### ADD THE HUMAN-OPTIMIZED PARAMETERS #######

        # map each human dataset to the reference genome to take the parameters that fit each df
        humanSample_to_refGenomeID = {"NA12878run1":"hg19", "HG002run1":"hg19", "CHMrun1":"hg38"}
        spName = "Homo_sapiens"

        # init dict
        parameters_df_dict = {}

        # go through each type of simulations
        for typeSimulations in ["uniform", "realSVs"]:

            # for human samples, where the directory structure is a bit different
            for sampleID, refGenomeID in humanSample_to_refGenomeID.items():

                # define the parameters
                parameters_json = "%s/outdirs_testing_humanGoldenSet/running_on_%s/testing_Accuracy/%s/%s/simulations_files_and_parameters/final_parameters.json"%(CurDir, refGenomeID, typeSimulations, sampleID)
                if fun.file_is_empty(parameters_json): raise ValueError("%s should exist"%parameters_json)

                # add to the dict
                parameters_df_dict[(spName, sampleID, typeSimulations)] = {"species":spName, "typeSimulations":typeSimulations, "sampleID":sampleID, "parameters_json":parameters_json}

        parameters_df_Human = pd.DataFrame(parameters_df_dict).transpose()[["species", "sampleID", "typeSimulations", "parameters_json"]]

        #################################################

        # define the merged parameters
        parameters_df = parameters_df_notHuman.append(parameters_df_Human)

        # get the df_bechmark for the human golden set
        outdir_cross_accuracy_human = "%s/cross_accuracy_human"%outdir_cross_accuracy
        df_benchmark_all = get_df_benchmark_cross_accuracy_humanGoldenSet(CurDir, parameters_df, outdir_cross_accuracy_human, humanSample_to_refGenomeID, threads, replace)

        # merge and save
        fun.save_df_as_tab(df_benchmark_all, df_benchmark_all_file)

    # load
    df_benchmark_all = fun.get_tab_as_df_or_empty_df(df_benchmark_all_file)

    ########### ADD METADATA #############

    # add the type of comparison
    print("running get_type_comparison")
    def get_type_comparison(r):

        # trained on fast
        if r["parms_species"]=="none": return "fast"

        # trained on the same sample as tested
        elif r["parms_species"]==r["test_species"] and r["parms_sampleID"]==r["test_sampleID"]: return "same_sample"

        # trained on the same species
        elif r["parms_species"]==r["test_species"]: return "same_species"

        # trained on other speceis
        elif r["parms_species"]!=r["test_species"]: return "different_species"

        else: raise ValueError("r is not valid")

    df_benchmark_all["type_comparison"] = df_benchmark_all.apply(get_type_comparison, axis=1)

    ######################################

    return df_benchmark_all


def generate_heatmap_accuracy_of_parameters_on_test_samples(df_benchmark, fileprefix, replace=False, threads=4, accuracy_f="Fvalue", svtype="integrated", col_cluster = False, row_cluster = False, show_only_species_and_simType=False, multiplier_width_colorbars=3, show_only_species=False):

    """
    This function takes a df where each row is one set of training parameters and test data svtype, together with the accuracy records. It generates a heatmap were the rows are each of the training parameters and the cols are the test samples.
    """

    print("plotting cross-accuracy")

    # define graphics
    #species_to_color = {'none': 'gray', 'Drosophila_melanogaster': 'black', 'Arabidopsis_thaliana': 'gray', 'Cryptococcus_neoformans': 'lightcoral', 'Candida_albicans': 'blue', 'Candida_glabrata': 'cyan'}

    #species_to_color = {'none': 'gray', 'Drosophila_melanogaster': 'darkorange', 'Arabidopsis_thaliana': 'olive', 'Cryptococcus_neoformans': 'lightcoral', 'Candida_albicans': 'magenta', 'Candida_glabrata': 'lightseagreen'}

    species_to_color = {'none': 'gray', 'Drosophila_melanogaster': 'darkorange', 'Arabidopsis_thaliana': 'olive', 'Cryptococcus_neoformans': 'lightcoral', 'Candida_albicans': 'magenta', 'Candida_glabrata': 'lightseagreen', "Homo_sapiens":"brown"}


    #typeSimulations_to_color = {"uniform":"blue", "realSVs":"red", "arroundRepeats":"black", "arroundHomRegions":"olive", "fast":"gray"}
    typeSimulations_to_color = {"uniform":"blue", "realSVs":"red", "arroundHomRegions":"black", "fast":"gray"}
    numericSample_to_color = {"0":"white", "1":"gray", "2":"black"}
    numericRun_to_color = {"0":"lightcoral", "1":"brown", "2":"r"}
    simName_to_color = {"sim1":"white", "sim2":"black"}


    cathegory_to_colors_dict = {"parms_species" : species_to_color,
                                "parms_typeSimulations" : typeSimulations_to_color,
                                "parms_numeric_sample":numericSample_to_color,
                                "parms_numeric_run":numericRun_to_color,
                                "test_species" : species_to_color,
                                "test_typeSimulations" : typeSimulations_to_color,
                                "test_numeric_sample":numericSample_to_color,
                                "test_numeric_run":numericRun_to_color,
                                "test_simName":simName_to_color}


    # keep only the df with the svtype
    df_plot = df_benchmark[(df_benchmark.svtype==svtype)] 


    # define square df
    parms_keys = ["parms_species", "parms_typeSimulations", "parms_numeric_sample", "parms_numeric_run"]
    test_keys = ["test_species", "test_typeSimulations", "test_numeric_sample", "test_numeric_run", "test_simName"]
    df_plot["parms_idx"] = df_plot.apply(lambda r: "||||".join([str(r[k]) for k in parms_keys]), axis=1)
    df_plot["test_idx"] = df_plot.apply(lambda r: "||||".join([str(r[k]) for k in test_keys]), axis=1)
    df_square = df_plot[["parms_idx", "test_idx", accuracy_f]].pivot(index='parms_idx', columns='test_idx', values=accuracy_f)

    # sort by species
    print("sorting")
    species_to_order =  {'none': 0, "Candida_glabrata":1, "Candida_albicans":2, "Cryptococcus_neoformans":3, "Arabidopsis_thaliana":4, "Drosophila_melanogaster":5, "Homo_sapiens":6}
    index_to_order = {c : species_to_order[c.split("||||")[0]]  for c in df_square.index}
    col_to_order = {c : species_to_order[c.split("||||")[0]]  for c in df_square.columns}

    sorted_index = sorted(df_square.index, key=(lambda x: index_to_order[x]))
    sorted_cols = sorted(df_square.columns, key=(lambda x: col_to_order[x]))

    df_square = df_square.loc[sorted_index, sorted_cols]

    # add the label
    type_comparison_to_label = {"same_run_and_simulation":"*", "fast":"", "same_species_and_simulation":"", "same_species":"", "different_species":""}
    df_plot["label"] = df_plot.type_comparison.apply(lambda x: type_comparison_to_label[x])
    df_annotations = df_plot[["parms_idx", "test_idx", "label"]].pivot(index='parms_idx', columns='test_idx', values="label").loc[sorted_index, sorted_cols]

    # define dicts mapping objects
    type_keys_to_keys = {"parms":parms_keys, "test":test_keys}

    # generate the cols colors df
    def get_colors_series(idx, type_keys="parms"):
        # type_keys can be parms or test

        # get the color dicts
        keys = type_keys_to_keys[type_keys]

        # get the content
        idx_content = idx.split("||||")

        # define the series
        field_to_color = {keys[I] : cathegory_to_colors_dict[keys[I]][c] for I,c in enumerate(idx_content)}

        return pd.Series(field_to_color)
    
    row_colors_df = pd.Series(df_square.index, index=df_square.index).apply(lambda x: get_colors_series(x, type_keys="parms"))
    col_colors_df = pd.Series(df_square.columns, index=df_square.columns).apply(lambda x: get_colors_series(x, type_keys="test"))

    # keep only simulation type and species twice
    if show_only_species_and_simType is True:

        row_colors_df = row_colors_df[["parms_species", "parms_typeSimulations"]]
        col_colors_df = col_colors_df[["test_species", "test_typeSimulations"]]

    # keep only species colors
    if show_only_species is True:

        row_colors_df = row_colors_df[["parms_species"]]
        col_colors_df = col_colors_df[["test_species"]]

    # get the plot
    filename = "%s_cross_accuracy_%s_%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype, col_cluster, row_cluster)
    print("getting %s"%filename)

    # define the figure size
    figsize = (int(len(df_square.columns)*0.03), int(len(df_square)*0.03))

    fun.plot_clustermap_with_annotation(df_square, row_colors_df, col_colors_df, filename, title="cross accuracy", col_cluster=col_cluster, row_cluster=row_cluster, colorbar_label=accuracy_f, adjust_position=True, legend=True, idxs_separator_pattern="||||", texts_to_strip={"L001"}, default_label_legend="control", df_annotations=df_annotations, cmap=sns.color_palette("RdBu_r", 50), ylabels_graphics_df=None, grid_lines=False, figsize=figsize, multiplier_width_colorbars=multiplier_width_colorbars, vmax=1.0, vmin=0.0, size_annot=12)


def get_value_to_color(values, palette="mako", n=100, type_color="rgb", center=None):

    """TAkes an array and returns the color that each array has. Checj http://seaborn.pydata.org/tutorial/color_palettes.html"""

    # get the colors
    colors = sns.color_palette(palette, n)

    # change the colors
    if type_color=="rgb": colors = colors
    elif type_color=="hex": colors = [rgb_to_hex(c) for c in colors]
    else: raise ValueError("%s is not valid"%palette)

    # if they are strings
    if type(list(values)[0])==str:

        palette_dict = dict(zip(values, colors))
        value_to_color = palette_dict

    # if they are numbers
    else:

        # map eaqually distant numbers to colors
        if center==None:
            min_palette = min(values)
            max_palette = max(values)
        else: 
            max_deviation = max([abs(fn(values)-center) for fn in [min, max]])
            min_palette = center - max_deviation
            max_palette = center + max_deviation

        all_values_palette = list(np.linspace(min_palette, max_palette, n))
        palette_dict = dict(zip(all_values_palette, colors))

        # get value to color
        value_to_color = {v : palette_dict[find_nearest(all_values_palette, v)] for v in values}

    return value_to_color, palette_dict

def get_df_benchmark_one_svtype_onlySamples_with_min_SVs(df, min_n_SVs):

    """Gets a benchmarking df of one SV type and filters out the samples that have less than min_n_SVs"""

    # filter out df
    df = df[(df.nevents>=min_n_SVs) | (df.nevents_highConf>=min_n_SVs)]

    # keep the rows with the samples as in test_sampleID
    valid_sampleIDs = set(df.test_sampleID).union({"fast"})
    df = df[df.parms_sampleID.isin(valid_sampleIDs)]

    return df

def generate_heatmap_accuracy_of_parameters_on_test_samples_realSVs(df_benchmark, fileprefix, replace=False, threads=4, accuracy_f="Fvalue", svtype="integrated", col_cluster = False, row_cluster = False, multiplier_width_colorbars=3, show_only_species=False, min_n_SVs=10):

    """
    This is similar to generate_heatmap_accuracy_of_parameters_on_test_samples, but adapted to the realSVs
    """

    print("plotting cross-accuracy")

    # define graphics
    species_to_color = {'none': 'gray', 'Drosophila_melanogaster': 'darkorange', 'Arabidopsis_thaliana': 'olive', 'Cryptococcus_neoformans': 'lightcoral', 'Candida_albicans': 'magenta', 'Candida_glabrata': 'lightseagreen', "Homo_sapiens":"brown"}

    typeSimulations_to_color = {"uniform":"blue", "realSVs":"red", "arroundHomRegions":"black", "fast":"gray"}

    all_sampleIDs = sorted(set(df_benchmark.parms_sampleID).union(set(df_benchmark.test_sampleID)))
    sampleID_to_color = get_value_to_color(all_sampleIDs)[0]

    cathegory_to_colors_dict = {"parms_species" : species_to_color,
                                "parms_typeSimulations" : typeSimulations_to_color,
                                "test_species" : species_to_color,
                                "parms_sampleID":sampleID_to_color,
                                "test_sampleID":sampleID_to_color}

    # keep only the df with the svtype
    df_plot = df_benchmark[(df_benchmark.svtype==svtype)]

    # keep only the samples where there are >=min_n_SVs
    df_plot = get_df_benchmark_one_svtype_onlySamples_with_min_SVs(df_plot, min_n_SVs)

    # define square df
    parms_keys = ["parms_species", "parms_typeSimulations", "parms_sampleID"]
    test_keys = ["test_species", "test_sampleID"]
    df_plot = df_plot.sort_values(by=(parms_keys + test_keys))

    df_plot["parms_idx"] = df_plot.apply(lambda r: "||||".join([str(r[k]) for k in parms_keys]), axis=1)
    df_plot["test_idx"] = df_plot.apply(lambda r: "||||".join([str(r[k]) for k in test_keys]), axis=1)
    df_square = df_plot[["parms_idx", "test_idx", accuracy_f]].pivot(index='parms_idx', columns='test_idx', values=accuracy_f)

    # sort by species
    print("sorting")
    species_to_order =  {'none': 0, "Candida_glabrata":1, "Candida_albicans":2, "Cryptococcus_neoformans":3, "Arabidopsis_thaliana":4, "Drosophila_melanogaster":5, "Homo_sapiens":6}
    index_to_order = {c : species_to_order[c.split("||||")[0]]  for c in df_square.index}
    col_to_order = {c : species_to_order[c.split("||||")[0]]  for c in df_square.columns}

    sorted_index = sorted(df_square.index, key=(lambda x: index_to_order[x]))
    sorted_cols = sorted(df_square.columns, key=(lambda x: col_to_order[x]))

    df_square = df_square.loc[sorted_index, sorted_cols]

    # add the label
    type_comparison_to_label = {"same_sample":"*", "fast":"", "same_species":"", "different_species":""}
    df_plot["label"] = df_plot.type_comparison.apply(lambda x: type_comparison_to_label[x])
    df_annotations = df_plot[["parms_idx", "test_idx", "label"]].pivot(index='parms_idx', columns='test_idx', values="label").loc[sorted_index, sorted_cols]

    # define dicts mapping objects
    type_keys_to_keys = {"parms":parms_keys, "test":test_keys}

    # generate the cols colors df
    def get_colors_series(idx, type_keys="parms"):
        # type_keys can be parms or test

        # get the color dicts
        keys = type_keys_to_keys[type_keys]

        # get the content
        idx_content = idx.split("||||")

        # define the series
        field_to_color = {keys[I] : cathegory_to_colors_dict[keys[I]][c] for I,c in enumerate(idx_content)}

        return pd.Series(field_to_color)
    
    row_colors_df = pd.Series(df_square.index, index=df_square.index).apply(lambda x: get_colors_series(x, type_keys="parms"))
    col_colors_df = pd.Series(df_square.columns, index=df_square.columns).apply(lambda x: get_colors_series(x, type_keys="test"))

    # keep only simulation type and species
    row_colors_df = row_colors_df[["parms_species", "parms_typeSimulations"]]
    col_colors_df = col_colors_df[["test_species"]]

    # keep only species colors
    if show_only_species is True:

        row_colors_df = row_colors_df[["parms_species"]]
        col_colors_df = col_colors_df[["test_species"]]

    # get the plot
    filename = "%s_cross_accuracy_%s_%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype, col_cluster, row_cluster)
    print("getting %s"%filename)

    # define the figure size
    figsize = (int(len(df_square.columns)*0.03), int(len(df_square)*0.03))

    fun.plot_clustermap_with_annotation(df_square, row_colors_df, col_colors_df, filename, title="cross accuracy", col_cluster=col_cluster, row_cluster=row_cluster, colorbar_label=accuracy_f, adjust_position=True, legend=True, idxs_separator_pattern="||||", texts_to_strip={"L001"}, default_label_legend="control", df_annotations=df_annotations, cmap=sns.color_palette("RdBu_r", 50), ylabels_graphics_df=None, grid_lines=False, figsize=figsize, multiplier_width_colorbars=multiplier_width_colorbars, vmax=1.0, vmin=0.0, size_annot=12)



def get_sorted_bam_by_readName(bam, threads=4, replace=False):

    """Sorts a bam file into sorted_bam"""

    sorted_bam = "%s.sorted"%bam
    sorted_bam_tmp = "%s.tmp"%sorted_bam

    # define outdir
    outdir = fun.get_dir(bam)

    if fun.file_is_empty(sorted_bam) or replace is True:

        # remove all temporary files generated previously in samtools sort (they'd make a new sort to be an error)
        for outdir_file in os.listdir(outdir): 
            fullfilepath = "%s/%s"%(outdir, outdir_file)
            if outdir_file.startswith(fun.get_file(bam)) and ".tmp." in outdir_file: fun.remove_file(fullfilepath)

        print("sorting bam")
        fun.run_cmd("%s sort -n --threads %i -o %s %s"%(fun.samtools, threads, sorted_bam_tmp, bam))

        os.rename(sorted_bam_tmp, sorted_bam)

    return sorted_bam

def get_fastqgz_from_bam(bamfile, threads=4, replace=False, already_sorted_by_readName=False):

    """This function takes a bamfile and writes the reads """

    # get the sorted bam and indexed bam
    if already_sorted_by_readName is False: sorted_bam = get_sorted_bam_by_readName(bamfile, threads=threads, replace=replace)
    else: sorted_bam = bamfile
    #index_bam = get_index_bam(sorted_bam, threads=threads, replace=replace)

    # define the reads
    reads1_gz = "%s.reads1.fastq.gz"%sorted_bam
    reads2_gz = "%s.reads2.fastq.gz"%sorted_bam

    if fun.file_is_empty(reads1_gz) or fun.file_is_empty(reads2_gz) or replace is True:

        # define the tmp reads
        reads1_tmp = "%s.reads1.tmp.fastq"%sorted_bam
        reads2_tmp = "%s.reads2.tmp.fastq"%sorted_bam
        reads1_tmp_gz = "%s.gz"%reads1_tmp
        reads2_tmp_gz = "%s.gz"%reads2_tmp

        # remove tmp files
        for f in [reads1_tmp, reads2_tmp, reads1_tmp_gz, reads2_tmp_gz]: fun.remove_file(f)

        # get reads
        print("getting reads from %s"%sorted_bam)
        bedtools_std = "%s.std"%reads1_tmp
        fun.run_cmd("%s bamtofastq -i %s -fq %s -fq2 %s > %s 2>&1"%(fun.bedtools, sorted_bam, reads1_tmp, reads2_tmp, bedtools_std))

        # gzip
        for f in [reads1_tmp, reads2_tmp]: fun.run_cmd("pigz --fast %s"%f)

        # rename
        os.rename(reads1_tmp_gz, reads1_gz)
        os.rename(reads2_tmp_gz, reads2_gz)

    return reads1_gz, reads2_gz


def get_crossbenchmarking_distributions_default_and_best(df_cross_accuracy_benchmark, fileprefix, accuracy_f="Fvalue", width_multiplier=2.8, legend_deviation=1.9):

    """Plots the svtype vs accuracy. The rows are different species. The cols are different types of simulations  """

    # keep only default and same sample
    df_cross_accuracy_benchmark = df_cross_accuracy_benchmark[df_cross_accuracy_benchmark.type_comparison.isin({"fast", "same_run_and_simulation"})]

    # define things
    sorted_species = ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana", "Drosophila_melanogaster", "Homo_sapiens"]
    sorted_simTypes = ["uniform", "realSVs", "arroundHomRegions"]
    simType_to_CorrectSimType = {"uniform":"random", "realSVs":"real", "arroundHomRegions":"homologous"}

    # sort by svtype
    sorted_svtypes = ["deletions", "tandemDuplications", "inversions", "insertions", "translocations", "integrated"]
    svtype_to_orderI = dict(zip(sorted_svtypes, range(len(sorted_svtypes))))
    df_cross_accuracy_benchmark["svtype_I"] = df_cross_accuracy_benchmark.svtype.apply(lambda x: svtype_to_orderI[x])

    # add the parameters
    typeComparison_to_newTypeComparison = {"fast":"default", "same_run_and_simulation":"optimized"}
    df_cross_accuracy_benchmark["parameters"] = df_cross_accuracy_benchmark.type_comparison.apply(lambda x: typeComparison_to_newTypeComparison[x])


    # define graphics
    parms_to_color = {"default":"gray", "optimized":"red"}

    # create small svtype
    svtype_to_shortSVtype = {"deletions":"del", "tandemDuplications":"tan", "insertions":"ins", "translocations":"tra", "inversions":"inv", "integrated":"all"}

    df_cross_accuracy_benchmark["svtype"] = df_cross_accuracy_benchmark.svtype.apply(lambda x: svtype_to_shortSVtype[x])


    # init fig
    nrows = len(sorted_species)
    ncols = len(sorted_simTypes)

    I = 1
    fig = plt.figure(figsize=(ncols*width_multiplier, nrows*1.3))

    for Ir, test_species in enumerate(sorted_species):
        for Ic, test_simType in enumerate(sorted_simTypes):
            print(test_species, test_simType)

            # get df
            df_plot =  df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.test_species==test_species) & (df_cross_accuracy_benchmark.test_typeSimulations==test_simType)].sort_values(by=["svtype_I"])

            #do not add empty
            if len(df_plot)==0: continue

            # get the subplot
            ax = plt.subplot(nrows, ncols, I); I+=1

            # get the stripplot
            ax = sns.stripplot(data=df_plot, x="svtype", y=accuracy_f, hue="parameters", palette=parms_to_color,  dodge=True,  linewidth=.3, edgecolor="black", size=4, alpha=.9)

            # set the ylims
            ax.set_ylim([-0.05, 1.1])

            # add title
            if Ir==0: ax.set_title(simType_to_CorrectSimType[test_simType])

            # define ylabels
            if Ir==2: ax.set_ylabel("%s\n%s. %s"%(accuracy_f, test_species.split("_")[0][0], test_species.split("_")[1]))
            else: ax.set_ylabel("%s. %s"%(test_species.split("_")[0][0], test_species.split("_")[1]))

            # remove the xticklabels except in the last row
            if Ir==(nrows-1) and Ic==1: ax.set_xlabel("SV type")
            else: ax.set_xlabel("")
            
            if test_species!=sorted_species[-1] and not (test_species=="Drosophila_melanogaster" and test_simType=="arroundHomRegions"): ax.set_xticklabels([])
            else: 
                plt.setp(ax.get_xticklabels(), rotation=90, fontsize=10)

            # add the 0.75 lines
            plt.axhline(0.75, linestyle="--", linewidth=.7, color="k")

            # remove the yticks unless it is the first column
            if Ic!=0: 
                ax.set_yticklabels([])
                ax.set_yticks([])
                ax.set_ylabel("")

            # get the legen only in the last box
            #if Ir==0 and Ic==(ncols-1): ax.legend(bbox_to_anchor=(1, 1), title="parameters") 
            if Ir==(nrows-1) and Ic==(ncols-2): ax.legend(bbox_to_anchor=(legend_deviation, 0.5), title="parameters") 

            elif len(df_plot)>0: ax.get_legend().remove()

    # spaces
    plt.subplots_adjust(wspace=0.05, hspace=0.1)
    filename = "%s_%s.pdf"%(fileprefix, accuracy_f)
    fig.savefig(filename, bbox_inches='tight')



def get_crossbenchmarking_distributions_differentSetsOfParameters(df_cross_accuracy_benchmark, fileprefix, accuracy_f="Fvalue", svtype="integrated"):

    """Takes the cross benchmarking dataset and plots on the X different simulation types. Each row is one species, and the color is the type of parameter"""


    # keep one df
    df_cross_accuracy_benchmark = df_cross_accuracy_benchmark[df_cross_accuracy_benchmark.svtype==svtype]

    # define parms
    sorted_species = ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana", "Drosophila_melanogaster", "Homo_sapiens"]
    sorted_simTypes = ["uniform", "realSVs", "arroundHomRegions"]

    # add the order of the type of teh comparison
    sorted_typeComparisons = ["fast", "different_species", "same_species", "same_species_and_simulation", "same_run_and_simulation"]
    typeComparison_to_orderI = dict(zip(sorted_typeComparisons, range(len(sorted_typeComparisons))))
    df_cross_accuracy_benchmark["type_comparison_I"] = df_cross_accuracy_benchmark.type_comparison.apply(lambda x: typeComparison_to_orderI[x])

    # add a shorter to type comparison
    typeComparison_to_newTypeComparison = {"fast":"default", "different_species": "different spp", "same_species":"same spp", "same_species_and_simulation":"same simulation", "same_run_and_simulation":"same sample"}
    df_cross_accuracy_benchmark["training parameters"] = df_cross_accuracy_benchmark.type_comparison.apply(lambda x: typeComparison_to_newTypeComparison[x])

    # define graphics
    trainingParmsToColor = {"default":"gray", "different spp": "greenyellow", "same spp":"lime", "same simulation":"green", "same sample":"black"}


    # init fig
    nrows = len(sorted_species)
    ncols = len(sorted_simTypes)
    I = 1
    fig = plt.figure(figsize=(ncols*0.8, nrows*1.3))

    for Ir, test_species in enumerate(sorted_species):
        for Ic, test_simType in enumerate(sorted_simTypes):
            print(test_species, test_simType)

            # get df
            df_plot =  df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.test_species==test_species) & (df_cross_accuracy_benchmark.test_typeSimulations==test_simType)].sort_values(by=["type_comparison_I"])

            # get the subplot
            ax = plt.subplot(nrows, ncols, I); I+=1

            # get the stripplot
            if len(df_plot)>0: ax = sns.stripplot(data=df_plot, x="test_typeSimulations", y=accuracy_f, hue="training parameters", palette=trainingParmsToColor,  dodge=True,  linewidth=.01, edgecolor="gray", size=3, alpha=.6)

            # set the ylims
            ax.set_ylim([-0.05, 1.1])

            # remove title
            ax.set_title("")
            if Ir==2: ax.set_ylabel("%s\n%s. %s"%(accuracy_f, test_species.split("_")[0][0], test_species.split("_")[1]))
            else: ax.set_ylabel("%s. %s"%(test_species.split("_")[0][0], test_species.split("_")[1]))

            # remove the xticklabels except in the last row
            ax.set_xlabel("")
            ax.set_xticks([])
            if test_species!=sorted_species[-1]: ax.set_xticklabels([])
            else: 
                ax.set_xlabel(test_simType, rotation=75)

            # remove the yticks unless it is the first column
            if Ic!=0: 
                ax.set_yticklabels([])
                ax.set_yticks([])
                ax.set_ylabel("")

            # get the legen only in the first box
            if Ir==0 and Ic==(ncols-1): ax.legend(bbox_to_anchor=(1, 1)) 
            elif len(df_plot)>0: ax.get_legend().remove()

    # spaces
    plt.subplots_adjust(wspace=0.00, hspace=0.1)
    filename = "%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype)
    fig.savefig(filename, bbox_inches='tight')


def get_crossbenchmarking_distributions_differentSetsOfParameters_realSVs_scatter(df_cross_accuracy_benchmark, fileprefix, accuracy_f="Fvalue", svtype="integrated", min_n_SVs=10):

    """Plots a scatterplot of maximum accuracy (y) vs accuray of different parameters (x). The color would be the type of parameters. The color would be the type of parameters and the symbol the sampleID. The rows will be different specues. The cols are different types of training parameters."""

    # keep one df
    df_cross_accuracy_benchmark = cp.deepcopy(df_cross_accuracy_benchmark[df_cross_accuracy_benchmark.svtype==svtype])

    # get only minSVs
    df_cross_accuracy_benchmark = get_df_benchmark_one_svtype_onlySamples_with_min_SVs(df_cross_accuracy_benchmark, min_n_SVs)


    # define parms
    sorted_species = ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana", "Drosophila_melanogaster", "Homo_sapiens"]

    # add the order of the type of teh comparison
    sorted_typeComparisons = ["fast", "different_species", "same_species", "same_sample"]
    typeComparison_to_orderI = dict(zip(sorted_typeComparisons, range(len(sorted_typeComparisons))))
    df_cross_accuracy_benchmark["type_comparison_I"] = df_cross_accuracy_benchmark.type_comparison.apply(lambda x: typeComparison_to_orderI[x])

    # define the type of training simulations
    sorted_training_typeSim = ["uniform", "realSVs", "arroundHomRegions"]

    # add a shorter to type comparison
    typeComparison_to_newTypeComparison = {"fast":"default", "different_species": "different spp", "same_species":"same spp", "same_sample":"same sample"}
    df_cross_accuracy_benchmark["training parameters"] = df_cross_accuracy_benchmark.type_comparison.apply(lambda x: typeComparison_to_newTypeComparison[x])

    # define graphics
    trainingParmsToColor = {"default":"gray", "different spp": "m", "same spp":"black"}
    #trainingParmsToLetter = {"default":"d", "different spp": "\neq", "same spp":"~", "same sample":"="}

    # add fields
    df_cross_accuracy_benchmark["test_typeSVs"] = "real SVs"

    # init fig
    nrows = len(sorted_species)
    ncols = 3
    I = 1
    fig = plt.figure(figsize=(ncols*2.8, nrows*2.8))

    for Ir, test_species in enumerate(sorted_species):

        # define the ylim
        #df_species =  df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.test_species==test_species)].sort_values(by=["type_comparison_I"])
        #lims = [min(df_species[accuracy_f])-0.05, max(df_species[accuracy_f])+0.05]
        #ylim = [-0.05, 1.05]
        lims = [-0.05, max(df_cross_accuracy_benchmark[accuracy_f])+0.05]

        for Ic, parms_typeSim in enumerate(sorted_training_typeSim):


            # get df
            df_plot =  df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.test_species==test_species) & (df_cross_accuracy_benchmark.parms_typeSimulations.isin({"fast", parms_typeSim}))].sort_values(by=["type_comparison_I"])

            # get the subplot
            ax = plt.subplot(nrows, ncols, I); I+=1

            # skip if empty
            if sum(df_plot["training parameters"]=="same sample")==0: 
                ax.set_xlabel("")
                ax.set_xticks([])
                ax.set_ylabel("")
                ax.set_yticks([])

                continue

            # get the df_plot of accuracy of the same sample
            df_plot_same_sample = df_plot[df_plot["training parameters"]=="same sample"]
            accuracy_same_sample_f = "%s optimized same sample"%accuracy_f 
            df_plot_same_sample[accuracy_same_sample_f] = df_plot_same_sample[accuracy_f]

            # add to df_plot the accuracy of the same sample
            df_plot = df_plot.merge(df_plot_same_sample[[accuracy_same_sample_f, "test_sampleID"]], on="test_sampleID", validate="many_to_one")
            df_plot = df_plot[df_plot["training parameters"]!="same sample"]
            if any(pd.isna(df_plot[accuracy_same_sample_f])): raise ValueError("There can't be NaNs in accuracy_same_sample_f")

            # create a df that has the 'training parameters' 
            ax = sns.scatterplot(data=df_plot, x=accuracy_f, y=accuracy_same_sample_f, style="test_sampleID", edgecolor=[trainingParmsToColor[x] for x in df_plot["training parameters"]], alpha=1, facecolor="none", linewidth=1.5)

            # set lims
            #lims = [-0.05, 1.05]
            #all_vals = df_plot[accuracy_f].append(df_plot[accuracy_same_sample_f])
            #lims = [min(all_vals)-.05, max(all_vals)+.05]

            ax.set_ylim(lims)
            ax.set_xlim(lims)

            # add horizontal line
            plt.plot([0, 1], [0, 1], linewidth=.7, color="gray", linestyle="--")


            # legend
            ax.get_legend().remove()
            #if Ir==0 and Ic==2: ax.legend(bbox_to_anchor=(1, 1)) 
            #elif len(df_plot)>0: ax.get_legend().remove()


            #if len(df_plot)>0: ax = sns.stripplot(data=df_plot, x="test_typeSVs", y=accuracy_f, hue="training parameters", palette=trainingParmsToColor,  dodge=True,  linewidth=.01, edgecolor="gray", size=3, alpha=.75, edgecolor="none")

            #if len(df_plot)>0: ax = sns.stripplot(data=df_plot, x="test_typeSVs", y=accuracy_f, hue="training parameters", palette=trainingParmsToColor,  dodge=True,  linewidth=.01, edgecolor="gray", size=3, alpha=.75)



            
def get_crossbenchmarking_distributions_differentSetsOfParameters_realSVs(df_cross_accuracy_benchmark, fileprefix, accuracy_f="Fvalue", svtype="integrated", min_n_SVs=10):

    """This is like get_crossbenchmarking_distributions_differentSetsOfParameters but for real SVs"""

    # imports
    from matplotlib.ticker import FormatStrFormatter

    # keep one df
    df_cross_accuracy_benchmark = cp.deepcopy(df_cross_accuracy_benchmark[df_cross_accuracy_benchmark.svtype==svtype])

    # keep only the samples with min_n_SVs
    df_cross_accuracy_benchmark = get_df_benchmark_one_svtype_onlySamples_with_min_SVs(df_cross_accuracy_benchmark, min_n_SVs)


    # define parms
    sorted_species = ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana", "Drosophila_melanogaster", "Homo_sapiens"]

    # add the order of the type of teh comparison
    sorted_typeComparisons = ["fast", "different_species", "same_species", "same_sample"]
    typeComparison_to_orderI = dict(zip(sorted_typeComparisons, range(len(sorted_typeComparisons))))
    df_cross_accuracy_benchmark["type_comparison_I"] = df_cross_accuracy_benchmark.type_comparison.apply(lambda x: typeComparison_to_orderI[x])

    # define the type of training simulations
    sorted_training_typeSim = ["uniform", "realSVs", "arroundHomRegions"]

    # add a shorter to type comparison
    typeComparison_to_newTypeComparison = {"fast":"default", "different_species": "different spp", "same_species":"same spp", "same_sample":"same sample"}
    df_cross_accuracy_benchmark["training parameters"] = df_cross_accuracy_benchmark.type_comparison.apply(lambda x: typeComparison_to_newTypeComparison[x])

    # define graphics
    trainingParmsToColor = {"default":"gray", "different spp": "greenyellow", "same spp":"lime", "same sample":"black"}
    trainingParmsToLetter = {"default":"d", "different spp": r'$\neq$', "same spp":"~", "same sample":"="}
    traininTypeSim_to_correctName = {"uniform":"random", "realSVs":"known", "arroundHomRegions":"homologous"}

    # add fields
    df_cross_accuracy_benchmark["test_typeSVs"] = "real SVs"

    # init fig
    nrows = len(sorted_species)
    ncols = 3
    I = 1
    fig = plt.figure(figsize=(ncols*0.8, nrows*1.8))


    for Ir, test_species in enumerate(sorted_species):

        # define the ylim
        df_species =  df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.test_species==test_species)].sort_values(by=["type_comparison_I"])
        #ylim = [min(df_species[accuracy_f])-0.05, max(df_species[accuracy_f])+0.05]
        ylim = [-0.05, 1.05]

        # print the number of events of each sample
        if test_species=="Homo_sapiens": field_nevents = "nevents"
        else: field_nevents = "nevents_highConf"
        df_species[field_nevents] = df_species[field_nevents].apply(int)

        print(df_species[["test_species", "test_sampleID", field_nevents]].drop_duplicates())

        for Ic, parms_typeSim in enumerate(sorted_training_typeSim):

            # get df
            df_plot =  df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.test_species==test_species) & (df_cross_accuracy_benchmark.parms_typeSimulations.isin({"fast", parms_typeSim}))].sort_values(by=["type_comparison_I"])

            # get the subplot
            ax = plt.subplot(nrows, ncols, I); I+=1

            # skip if empty
            if sum(df_plot["training parameters"]=="same sample")==0: 
                ax.set_xlabel("")
                ax.set_xticks([])
                ax.set_ylabel("")
                ax.set_yticks([])

                continue



            # get the stripplot
            #if len(df_plot)>0: ax = sns.stripplot(data=df_plot, x="test_typeSVs", y=accuracy_f, hue="training parameters", palette=trainingParmsToColor,  dodge=True,  linewidth=.01, edgecolor="gray", size=3, alpha=.75)

            # add lines pairing the same test_sampleID
            sorted_training_parameters = [x for x in ["default", "different spp", "same spp", "same sample"] if x in set(df_plot["training parameters"])]

            sorted_samples = sorted(set(df_plot.test_sampleID))
            sampleID_to_color = fun.get_value_to_color(sorted_samples, palette="tab10", type_color="rgb", center=None)[0]


            for Ifrom in range(len(sorted_training_parameters)-1):
                Ito = Ifrom+1

                # define the parameters
                parm_from = sorted_training_parameters[Ifrom]
                parm_to = sorted_training_parameters[Ito]

                # go through all samples
                for sampleID in set(df_plot.test_sampleID):
                    
                    # get the df of the from and to
                    df_from = df_plot[(df_plot.test_sampleID==sampleID) & (df_plot["training parameters"]==parm_from)]
                    df_to = df_plot[(df_plot.test_sampleID==sampleID) & (df_plot["training parameters"]==parm_to)]

                    # go through each
                    for from_y in df_from[accuracy_f].values:
                        for to_y in df_to[accuracy_f].values:

                            # get the xfrom to 
                            #from_x = np.mean(np.array(ax.collections[Ifrom].get_offsets()).T[0])
                            #to_x = np.mean(np.array(ax.collections[Ito].get_offsets()).T[0])

                            plt.plot([Ifrom, Ito], [from_y, to_y], linestyle="-", linewidth=.7, color=sampleID_to_color[sampleID], alpha=.8)

            # set the ylims
            ax.set_ylim(ylim)

            # add the lines
            for y in [0.5, 0.75]: plt.axhline(y, color="gray", linestyle="--", linewidth=.7)

            # title for the first row
            if Ir==0: ax.set_title(traininTypeSim_to_correctName[parms_typeSim],  rotation=75)

            # ylabel (only first col)
            if Ic==0:

                if Ir==3: ax.set_ylabel("%s\n%s. %s"%(accuracy_f, test_species.split("_")[0][0], test_species.split("_")[1]))
                else: ax.set_ylabel("%s. %s"%(test_species.split("_")[0][0], test_species.split("_")[1]))

                ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

            else:
                ax.set_yticklabels([])
                ax.set_ylabel("")
                ax.set_yticks([])

            # remove the xticklabels except in the last row
            ax.set_xlabel("")
            #ax.set_xticks([])
            ax.set_xticks(list(range(len(sorted_training_parameters))))
            ax.set_xticklabels([trainingParmsToLetter[p] for p in sorted_training_parameters])

            #if test_species!=sorted_species[-1]: ax.set_xticklabels([])
            #elif Ic==1: 
            if Ic==1 and test_species==sorted_species[-1]: ax.set_xlabel("real SVs", rotation=0)

            # get the legen only in the first box
            #if Ir==0 and Ic==2: ax.legend(bbox_to_anchor=(1, 1)) 
            #elif len(df_plot)>0: ax.get_legend().remove()

    # spaces
    print("saving")
    plt.subplots_adjust(wspace=0.1, hspace=0.3)
    filename = "%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype)
    fig.savefig(filename, bbox_inches='tight')

def get_crossaccuracy_distributions(df_cross_accuracy_benchmark, fileprefix, accuracy_f="Fvalue", svtype="integrated"):

    """Prints a boxplot for each species and the accuracy on different types of simulations depending on the training parameters. Each row is a tested species and the X is the type of parameters"""

    # keep one df
    df_cross_accuracy_benchmark = df_cross_accuracy_benchmark[df_cross_accuracy_benchmark.svtype==svtype]

    # define parms
    sorted_species = ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana", "Drosophila_melanogaster", "Homo_sapiens"]
    typeSimulations_to_color = {"uniform":"blue", "realSVs":"red", "arroundHomRegions":"black", "fast":"gray"}
    type_comparison_to_xvalue = {"fast":0, "same_run_and_simulation":1, "same_species_and_simulation":2, "same_species":3, "different_species":4}
    xvalue_to_typeCompLabel = {0:"default parms.", 1:"same sample", 2:"same simulation", 3:"same species", 4:"different species"}

    df_cross_accuracy_benchmark["training parameters"] = df_cross_accuracy_benchmark.type_comparison.apply(lambda x: type_comparison_to_xvalue[x])
    all_xvalues = sorted(set(df_cross_accuracy_benchmark["training parameters"]))

    # init fig
    nrows = len(sorted_species)
    fig = plt.figure(figsize=(5*0.6, nrows*1.3))

    for I, test_species in enumerate(sorted_species):

        # get df
        df_plot =  df_cross_accuracy_benchmark[df_cross_accuracy_benchmark.test_species==test_species]

        # get the subplot
        ax = plt.subplot(nrows, 1, I+1)

        # get the violin
        #ax = sns.violinplot(data=df_plot, x="training parameters", y=accuracy_f, hue="test_typeSimulations", palette=typeSimulations_to_color)
        #ax = sns.boxplot(data=df_plot, x="training parameters", y=accuracy_f, hue="test_typeSimulations", palette=typeSimulations_to_color)

        ax = sns.stripplot(data=df_plot, x="training parameters", y=accuracy_f, hue="test_typeSimulations", palette=typeSimulations_to_color,  dodge=True,  linewidth=.01, edgecolor="gray", size=3, alpha=.3)
        #ax = sns.stripplot(x="training parameters", y=accuracy_f, hue="test_typeSimulations", data=df_plot, dodge=True, linewidth=.01, edgecolor="gray", size=3, palette=typeSimulations_to_color)
        #ax = sns.swarmplot(x="training parameters", y=accuracy_f, hue="test_typeSimulations", data=df_plot, linewidth=.01, edgecolor="gray", size=3, palette=typeSimulations_to_color)

        # change the edge colors
        """
        for i,box in enumerate(ax.artists):
            #box.set_edgecolor(box.get_facecolor())
            box.set_edgecolor("gray")
            box.set_linewidth(0.02)
        """

        # remove title
        ax.set_title("")
        ax.set_ylabel("%s. %s\n%s"%(test_species.split("_")[0][0], test_species.split("_")[1], accuracy_f))

        # change the xticklabels
        if I!=(nrows-1):
            ax.set_xticklabels("")
            ax.set_xlabel("")
            ax.set_xticks(all_xvalues)
        else:
            ax.set_xticklabels([xvalue_to_typeCompLabel[x] for x in all_xvalues])
            for label in ax.get_xticklabels(): label.set_rotation(90)

        ax.set_ylim([-0.05, 1.1])

        # add axvlines
        for x in [0.5, 1.5, 2.5, 3.5]: plt.axvline(x, linewidth=.9, color="gray", linestyle="--")


        # get the legen only in the first box
        ax.legend(bbox_to_anchor=(1, 1))
        if I!=0: ax.get_legend().remove()

    # spaces
    plt.subplots_adjust(wspace=0.01, hspace=0.01)
    filename = "%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype)
    fig.savefig(filename, bbox_inches='tight')
    #plt.close(fig)

def plot_accuracy_distributions_sameRun_bySpecies(df_cross_accuracy_benchmark, fileprefix, accuracy_f="Fvalue", all_types_simulations={"fast", "uniform"}, all_types_comparisons={'fast', 'same_run_and_simulation'}, svtype="integrated", ylim=[0,1.1]):

    """This plots on the X the test species and the hue is the type_comparison"""


    # get the df
    df_plot = df_cross_accuracy_benchmark[(df_cross_accuracy_benchmark.parms_typeSimulations.isin(all_types_simulations)) & (df_cross_accuracy_benchmark.test_typeSimulations.isin(all_types_simulations)) & (df_cross_accuracy_benchmark.type_comparison.isin(all_types_comparisons)) & (df_cross_accuracy_benchmark.svtype==svtype)]

    # define vars
    palette = {"fast":"silver", "same_run_and_simulation":"black"}

    # plot
    fig = plt.figure(figsize=(5, 3))


    # add violin
    ax = sns.violinplot(x="test_species", y=accuracy_f, hue="type_comparison", data=df_plot, palette=palette, linewidth=.5, color=palette, dodge=True, split=False)

    for art in ax.get_children():
        if type(art)==matplotlib.collections.PolyCollection: art.set_alpha(0.4)


    #ax = sns.stripplot(data=df_plot, x="test_species", y=accuracy_f, hue="type_comparison", palette=palette,  dodge=True,  linewidth=.01, edgecolor="gray", size=3, alpha=.4)
    ax = sns.swarmplot(data=df_plot, x="test_species", y=accuracy_f, hue="type_comparison", palette=palette,  dodge=True,  linewidth=.01, edgecolor="dimgray", size=6, alpha=.9)

    


    for label in ax.get_xticklabels(): label.set_rotation(90)
    ax.legend(bbox_to_anchor=(1, 1))
    ax.set_ylim(ylim)

    filename = "%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype)
    fig.savefig(filename, bbox_inches='tight')


def rsync_file_shh(origin, target):

    """Copy a file with tmp and rsync"""

    fun.run_cmd("rsync -v %s %s"%(origin, target))


def run_parallelFastqDump_fromSRR_pairedReads_localComputer(srr, outdir, replace=False, threads=4):

    """Runs parallel fastqdump from an srr into outdir"""

    # make dirs
    fun.make_folder(outdir)

    # define the reads
    final_reads1 = "%s/%s_1.fastq.gz"%(outdir, srr)
    final_reads2 = "%s/%s_2.fastq.gz"%(outdir, srr)

    if fun.file_is_empty(final_reads1) or fun.file_is_empty(final_reads2) or replace is True:

        # check that you are locally
        if not str(subprocess.check_output("uname -a", shell=True)).startswith("b'Linux bscls063 4.12.14-lp150.12.48-default"): raise ValueError("You are not running in the BSC local machine") 

        # define the local outdir to download
        local_outdir = "/data/projects/downloading_%s"%srr; 
        fun.make_folder(local_outdir)
        os.chdir(local_outdir)

        # download prefetched SRA file
        prefetch_outdir = "%s/downloading_SRAfile"%local_outdir
        SRAfile = "%s/%s/%s.sra"%(prefetch_outdir, srr, srr)

        prefetch_finished_file = "%s/prefetch_finished.txt"%local_outdir
        if fun.file_is_empty(prefetch_finished_file) or replace is True:

            fun.delete_folder(prefetch_outdir)
            prefetch_std = "%s/prefetch_std.txt"%(local_outdir)
            print("running prefetch. STD in %s"%prefetch_std)
            fun.run_cmd("%s --output-directory %s --max-size 500G --progress %s > %s 2>&1"%(fun.prefetch, prefetch_outdir, srr, prefetch_std))
            open(prefetch_finished_file, "w").write("prefetch finished")    

        # clean
        for f in os.listdir(local_outdir): 
            if f not in {"downloading_SRAfile", srr, "prefetch_finished.txt", "prefetch_std.txt"}: fun.delete_file_or_folder("%s/%s"%(local_outdir, f))

        # run fastqdump parallel into the tmpdir 
        stdfile = "%s/std_fastqdump.txt"%local_outdir
        print("running fastq dump in parallel. The std is in %s"%stdfile)
        fun.run_cmd("%s -s %s -t %i -O %s --tmpdir %s --split-3 --gzip -vvv > %s 2>&1 "%(fun.parallel_fastq_dump, SRAfile, threads, local_outdir, local_outdir, stdfile))

        # check that the fastqdump is correct
        std_lines = open(stdfile, "r").readlines()
        any_error = any(["ERROR" in l.upper() for l in std_lines])
        if any_error:
            raise ValueError("Something went wrong with the fastqdump. Check the log in %s"%stdfile)

        # move reads to outdir
        local_reads1 = "%s/%s_1.fastq.gz"%(local_outdir, srr)
        local_reads2 = "%s/%s_2.fastq.gz"%(local_outdir, srr)

        # define options for rsync
        gpfs_parentDir = "/gpfs/projects/bsc40/mschikora"
        cluster_name = "bsc40395@dt01.bsc.es"

        print("rsyncing reads1")
        clusterDir_reads1 = "%s:%s/%s"%(cluster_name, fun.get_dir(final_reads1).replace(ParentDir, gpfs_parentDir), fun.get_file(final_reads1))
        rsync_file_shh(local_reads1, clusterDir_reads1)

        print("rsyncing reads2")
        clusterDir_reads2 = "%s:%s/%s"%(cluster_name, fun.get_dir(final_reads2).replace(ParentDir, gpfs_parentDir), fun.get_file(final_reads2))
        rsync_file_shh(local_reads2, clusterDir_reads2)

        # clean
        fun.delete_folder(local_outdir)

    return final_reads1, final_reads2


def merge_reads_into_one_file(readsA, readsB, merged_reads, replace=False):

    """Takes 2 reads and merges them into a single file"""

    if fun.file_is_empty(merged_reads) or replace is True:

        merged_reads_tmp = "%s.tmp.fastq.gz"%merged_reads
        print("merging reads into %s"%merged_reads)

        fun.run_cmd("cat %s %s > %s"%(readsA, readsB, merged_reads_tmp))
        os.rename(merged_reads_tmp, merged_reads)



def get_correct_human_genome(raw_genome, type_genome="hg19"):

    """This function takes a human genome and rewrites it so that it is correct"""

    # define the corrected_genome
    corrected_genome = "%s.corrected.fasta"%raw_genome

    if fun.file_is_empty(corrected_genome):
        print("getting corrected genome for %s"%type_genome)

        # load genome and keep only the important scaffolds
        if type_genome=="hg19": chr_to_seq = {seq.id : seq for seq in fun.SeqIO.parse(raw_genome, "fasta") if "_" not in seq.id and seq.id!="chrM"} # remove all the chromosomes with a '_' in it. I also remove chrM, which is the old version of the mtDNA

        elif type_genome=="hg38": chr_to_seq = {seq.id : seq for seq in fun.SeqIO.parse(raw_genome, "fasta") if "_" not in seq.id} # same as hg19, but chromosome is called M
        
        else: raise ValueError("%s is not valid"%type_genome)

        # get the length
        genome_size = sum(map(len, chr_to_seq.values()))
        print("The genome %s has %.3f Gbp"%(type_genome, genome_size/1e9))


        for c in sorted(chr_to_seq.keys()): print(c)

        # write
        corrected_genome_tmp = "%s.tmp.fasta"%corrected_genome
        fun.SeqIO.write(list(chr_to_seq.values()), corrected_genome_tmp, "fasta")
        os.rename(corrected_genome_tmp, corrected_genome)

    return corrected_genome


def get_blastn_regions_genome_against_itself_equal_regions(windows_multifasta, replace=False, threads=4, window_size=500):

    """Takes a multifasta and returns a blastn of regions that are the same"""

    blastn_outfile = "%s.equal_sequences_like_blastn.tab"%windows_multifasta
    if fun.file_is_empty(blastn_outfile) or replace is True:

        print("getting %s"%blastn_outfile)

        # load multifasta in as a seq_df
        def get_id_from_seq(seq): return seq.id
        def get_str_from_seq(seq): return str(seq.seq).upper()

        indices = list(map(get_id_from_seq, fun.SeqIO.parse(windows_multifasta, "fasta")))
        sequences = list(map(get_str_from_seq, fun.SeqIO.parse(windows_multifasta, "fasta")))

        seq_df = pd.DataFrame({"sequence":sequences, "seq_name":indices}).sort_values(by="sequence")

        # get a series that maps each sequence to a list of pairs of different indices
        def get_pairs_homRegions(df_s): return list(itertools.combinations(df_s.seq_name, 2))
        seq_to_pairsEqualRegions = seq_df.groupby("sequence").apply(get_pairs_homRegions).reset_index(drop=True).apply(set)

        # get all the pairs of identical regions
        all_pairs_identical = pd.Series(list(set.union(*seq_to_pairsEqualRegions)))
        #all_pairs_identical = pd.Series(list({("chrX||59229000||59229500", "chrX||59229000||59229500"), ("chrX||59229000||59229500", "chr3||131069500||131070000")})) # example to debug

        if len(all_pairs_identical)==0: raise ValueError("There are no identical regions")

        # get a df with the pairs of identical regions
        def get_r_blastn_from_pairIDs(pair):

            query_chromosome, qstart, qend = pair[0].split("||")
            subject_chromosome, sstart, send = pair[1].split("||")

            return pd.Series({"qseqid":pair[0], "sseqid":pair[1], "query_chromosome":query_chromosome, "subject_chromosome":subject_chromosome, "qstart":qstart, "sstart":sstart, "qend":qend, "send":send})

        blast_df = all_pairs_identical.apply(get_r_blastn_from_pairIDs)

        # remove regions that are unnecessary
        blast_df = blast_df[(blast_df.qseqid!=blast_df.sseqid)]

        # change to ints
        blast_df["qstart"] = blast_df.qstart.apply(int)
        blast_df["qend"] = blast_df.qend.apply(int)
        blast_df["sstart"] = blast_df.sstart.apply(int)
        blast_df["send"] = blast_df.send.apply(int)

        # add uniform fields
        blast_df["qlen"] = window_size
        blast_df["slen"] = window_size
        blast_df["evalue"] = 0.0 
        blast_df["bitscore"] = 100.0
        blast_df["score"] = 100.0 
        blast_df["length"] =  window_size
        blast_df["pident"] = 100.0 
        blast_df["nident"] = window_size 
        blast_df["qcovs"] = 100.0
        blast_df["sstart_0based"] =  blast_df.sstart
        blast_df["subject_like_qseqid"] = blast_df.sseqid
        blast_df["query_start"] = blast_df.qstart
        blast_df["query_end"] = blast_df.qend

        # add the centers
        blast_df["query_center"] = blast_df.query_start + (int(window_size/2))
        blast_df["subject_center"] = blast_df.sstart + (int(window_size/2))

        # save df
        fun.save_df_as_tab(blast_df, blastn_outfile)

    return blastn_outfile

def print_df_keys(df):

    """print df keys"""

    for k in df.keys(): 
        if any(~pd.isna(df[k])): print( "\n", k, "\n", df[k],  "\n")


def get_flagstat_file_from_bamfile(bamfile, threads=4, replace=False):

    """runs flagstat on bamfile"""

    flagstat_file = "%s.flagstat"%bamfile; flagstat_file_tmp = "%s.tmp"%flagstat_file

    # get the total n reads
    if fun.file_is_empty(flagstat_file) or replace is True:

        read_count_stderr = "%s.generating.stderr"%flagstat_file
        print("calculating n reads. The stderr is in %s"%read_count_stderr)
        fun.run_cmd("%s flagstat --threads %i %s > %s 2>%s"%(fun.samtools, threads, bamfile, flagstat_file_tmp, read_count_stderr))
        fun.remove_file(read_count_stderr)

        os.rename(flagstat_file_tmp, flagstat_file)

    return flagstat_file


def get_df_resources_simulations(CurDir, threads):

    """Generates a df with the resources used by each simulation in perSVade. It takes the 5 testing species and the hg38 human datasets. This was run on 16 threads and 6000 Mb RAM/thread (1800 Mb RAM/thread) in Nord3."""


    ########### DEFINE A DF WITH ALL THE CONCATENATED PATHS TO STDOUTs ################

    # add the human ones
    df_stdouts_human = pd.DataFrame()

    human_dir = "%s/outdirs_testing_humanGoldenSet/all_STDs_testingAccuracySeveralSpecies"%CurDir
    for f in [x for  x in os.listdir(human_dir) if x.startswith("hg38_") and x.endswith("_jobEnding.txt")]: 
        df_stdouts_human = df_stdouts_human.append(fun.get_tab_as_df_or_empty_df("%s/%s"%(human_dir, f)))

    df_stdouts_human["species"] = "Homo_sapiens"
    df_stdouts_human["dest_stdout"] = ["%s/%s"%(human_dir, f) for f in df_stdouts_human.dest_stdout]
    df_stdouts_human["jobs_file"] = ["%s/%s"%(human_dir, f) for f in df_stdouts_human.jobs_file]


    # add all the testing species
    df_stdouts_testing = pd.DataFrame()

    testingSpecies_dir = "%s/outdirs_testing_severalSpecies/all_STDs_testingAccuracySeveralSpecies"%CurDir
    for f in [x for  x in os.listdir(testingSpecies_dir) if x.endswith("_jobEnding.txt")]: 
        df_stdouts_testing = df_stdouts_testing.append(fun.get_tab_as_df_or_empty_df("%s/%s"%(testingSpecies_dir, f)))


    df_stdouts_testing["species"] = df_stdouts_testing.jobs_file.apply(lambda f: "_".join(f.split("_")[0:2]))
    df_stdouts_testing["dest_stdout"] = ["%s/%s"%(testingSpecies_dir, f) for f in df_stdouts_testing.dest_stdout]
    df_stdouts_testing["jobs_file"] = ["%s/%s"%(testingSpecies_dir, f) for f in df_stdouts_testing.jobs_file]

    # merge and add fields
    df_stdouts = df_stdouts_human.append(df_stdouts_testing)
    df_stdouts["outdir_task"] = df_stdouts.outdir_task.apply(lambda x: x.replace("/gpfs/projects/bsc40/mschikora", ParentDir))
    df_stdouts["type_simulation"] = df_stdouts.outdir_task.apply(lambda x: x.split("/")[-2])
    df_stdouts["sampleID"] = df_stdouts.outdir_task.apply(lambda x: x.split("/")[-1])

    ################################################################################### 

    ######### DEFINE THE RESOURCES DF ##########

    # create a df were each line is a combination of species, sampleID and type_simulation. Note that taskID indicates the number in the stdout related to this command

    def add_time_and_memory_to_df_stdouts_r(r):

        # define the lines in the stdout
        main_stdout_lines = [l.strip() for l in open(r.dest_stdout, "r").readlines()]

        # get the job ID
        all_jobIDs = {int(l.split("Subject: Job ")[1].split("[")[0]) for l in main_stdout_lines if l.startswith("Subject: Job")}
        if len(all_jobIDs)!=1: raise ValueError("There are more than 1 job IDS: %s"%all_jobIDs)
        jobID = next(iter(all_jobIDs))

        # get the line in stdout with related to this task
        jobID_line = [Iline for Iline,l in enumerate(main_stdout_lines) if l.startswith("Subject: Job %i[%i]:"%(jobID, r.taskID))][0]
        requested_mem_Mb = [float(l.split("Total Requested Memory :")[1].split()[0]) for l in main_stdout_lines[jobID_line:] if "Total Requested Memory :" in l and "MB" in l][0]
        max_mem_Mb = [float(l.split("Max Memory :")[1].split()[0]) for l in main_stdout_lines[jobID_line:] if "Max Memory :" in l and "MB" in l][0]
        cpu_time = [float(l.split("CPU time :")[1].split()[0]) for l in main_stdout_lines[jobID_line:] if "CPU time :" in l and "sec" in l][0]

        # get the elapsed time in seconds
        start_date =  [l.strip().split()[3:] for l in main_stdout_lines[jobID_line:] if l.startswith("Started at")][0]
        end_date =  [l.strip().split()[4:] for l in main_stdout_lines[jobID_line:] if l.startswith("Results reported on")][0]

        s_month = start_date[0]
        s_day = int(start_date[1])
        s_h, s_min, s_sec = [float(x) for x in start_date[2].split(":")]
        s_year = int(start_date[3])

        e_month = end_date[0]
        e_day = int(end_date[1])
        e_h, e_min, e_sec = [float(x) for x in end_date[2].split(":")]
        e_year = int(end_date[3])  

        # define the number of days that distanced the start and the end. Each month is particular
        if s_year==e_year and s_year==2021 and s_month=="Feb" and e_month=="Mar": transcurred_days = e_day - (s_day-28)

        # most cases
        else: 
            transcurred_days = e_day-s_day
            if s_month!=e_month: raise ValueError("s_month %s is different to e_month %s"%(s_month, e_month))  
            if s_year!=e_year: raise ValueError("s_year %s is different to e_year %s"%(s_year, e_year))  

        # get the total time
        run_time =  transcurred_days*(24*3600) + (e_h-s_h)*3600 + (e_min-s_min)*60 + (e_sec-s_sec)

        # checks
        for x in [cpu_time, max_mem_Mb, requested_mem_Mb, run_time]:
            if pd.isna(x) or x<=0.0: raise ValueError("there is an error with the parsing of the output: %s"%x)

        # add things to r
        r["run_time"] = run_time
        r["max_mem_Mb"] = max_mem_Mb
        r["cpu_time"] = cpu_time
        r["requested_mem_Mb"] = requested_mem_Mb

        return r

    def get_resources_dict(df):

        # keep indexes
        species, sampleID, type_simulation = df.name

        # add the time and the row
        df = df.apply(add_time_and_memory_to_df_stdouts_r, axis=1)

        # check that there is only one row with 'perSVade_finished'
        if sum(df.perSVade_finished)!=1: raise ValueError("There should be only one job where perSVade was finished")

        # check that the requested memory was always the same
        if len(set(df.requested_mem_Mb))!=1: raise ValueError("the requested memory was different across runs")

        # check that the outdir task is the same
        if len(set(df.outdir_task))!=1: raise ValueError("the outdir task was different across runs")

        # define the dict with all the data
        data_dict = {"species":species, "sampleID":sampleID, "type_simulation":type_simulation, "run time (h)":sum(df.run_time)/3600, "max memory (Gb)":max(df.max_mem_Mb)/1000, "cpu_time":sum(df.cpu_time), "requested_mem_Mb":df.requested_mem_Mb.iloc[0],  "outdir_task":df.outdir_task.iloc[0]}

        return data_dict

    df_resources = df_stdouts[df_stdouts.perSVade_was_ran].groupby(["species", "sampleID", "type_simulation"]).apply(get_resources_dict).reset_index(drop=True).apply(pd.Series)

    ############################################

    ####### ADD FIELDS FOR PLOTTING #########
    print("adding extra fields for plotting")

    # add the genome size in bp
    species_to_genome_size = {species :  sum(fun.get_chr_to_len("%s/genomes_and_annotations/%s.fasta"%(CurDir, species)).values())  for species in set(df_resources.species)}
    df_resources["genome size (Mb)"] = df_resources.species.apply(lambda x: species_to_genome_size[x])/1e6

    # add the number of reads 
    def get_number_pairs_from_outdir_task(o):

        flagstat_lines = open(get_flagstat_file_from_bamfile("%s/aligned_reads.bam.sorted"%o, threads=threads), "r").readlines()

        return [int(l.split()[0]) for l in flagstat_lines if "mapped" in l][0]

    df_resources["mapped pairs"] = df_resources.outdir_task.apply(get_number_pairs_from_outdir_task)

    #########################################

    return df_resources

def plot_used_resources_testing_on_simulations(CurDir, ProcessedDataDir, PlotsDir, threads):

    """This function plots the resources used in simulations, sepparated by simulation type, taking the hg38 as example."""

    ########## GET THE DATASET WITH THE USED RESOURCES ########

    df_resources_file = "%s/used_resources_simulations.py"%ProcessedDataDir
    if fun.file_is_empty(df_resources_file):

        print("getting resources df")
        df_resources = get_df_resources_simulations(CurDir, threads)

        fun.save_object(df_resources, df_resources_file)

    df_resources = fun.load_object(df_resources_file)


    ###########################################################

    ######## PLOT THE SCATTER WITH TIME AND MAX MEMORY CONSUMED ######

    print("plotting")

    # define graphics
    typeSimulations_to_color = {"uniform":"blue", "realSVs":"red", "arroundHomRegions":"black", "fast":"gray"}
    species_to_marker = {"Candida_glabrata":"o", "Candida_albicans":"<", "Cryptococcus_neoformans":"s", "Arabidopsis_thaliana":"^", "Drosophila_melanogaster":"P", "Homo_sapiens":"*"}

    # sort for plotting 
    df_resources = df_resources.sort_values(by=["species", "type_simulation"])

    print(df_resources[["species", "genome size (Mb)"]].drop_duplicates())
    print(df_resources[["species", "mapped pairs"]].drop_duplicates().sort_values(by=["mapped pairs"]))


    # init figure
    fig = plt.figure(figsize=(6,6))

    # go through different subplots
    for I, (xfield, yfield) in enumerate([("genome size (Mb)", "run time (h)"), ("mapped pairs", "run time (h)"), ("genome size (Mb)", "max memory (Gb)"), ("mapped pairs", "max memory (Gb)")]):

        Ip = I+1

        # init the subplot
        ax = plt.subplot(2,2,Ip)
        ax = sns.scatterplot(data=df_resources, x=xfield, y=yfield, style="species", alpha=.95, markers=species_to_marker, edgecolor="white", linewidth=0.02, hue="type_simulation", palette=typeSimulations_to_color)

        #ax = sns.scatterplot(data=df_resources, x=xfield, y=yfield, style="species", alpha=.95, markers=species_to_marker, edgecolor="gray", linewidth=0.9, hue="type_simulation", palette=typeSimulations_to_color, edgecolor=[typeSimulations_to_color[ts] for ts in df_resources.type_simulation], facecolor="white")

        # set scales
        ax.set_xscale("log")
        if yfield=="run time (h)": ax.set_yscale("log")


        # remove axes
        if Ip in {1,2}: 
            #ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_xlabel("")

        if Ip in {2,4}: 
            ax.set_ylabel("")
            ax.set_yticklabels([])

        # add title 
        if Ip==1: ax.set_title("perSVade resource consumption")

        # add legend
        if Ip==2: ax.legend(bbox_to_anchor=(2.1, 0), loc="right", borderaxespad=0.)
        else: ax.get_legend().remove()

    plt.subplots_adjust(wspace=0.01, hspace=0.01)
    fig.savefig("%s/resource_consumption_Nord3_simulations.pdf"%PlotsDir, bbox_inches='tight')

    ##################################################################

svtype_to_positionField_to_chromosome = {"inversions": {"Start":"Chr", "End":"Chr"},
                                         "tandemDuplications": {"Start":"Chr", "End":"Chr"},
                                         "deletions": {"Start":"Chr", "End":"Chr"},
                                         "insertions": {"StartA":"ChrA", "EndA":"ChrA", "StartB":"ChrB"},
                                         "translocations_noINV": {"EndA":"ChrA", "EndB":"ChrB"},
                                         "translocations_INV": {"EndA":"ChrA", "StartB":"ChrB"},
                                         "remaining": {"POS":"#CHROM", "START":"CHR2", "END":"CHR2"}}

def get_sv_dict_without_simple_repeats(sv_dict, repeats_file, outdir, reference_genome):

    """Gets an svdict without variants that overlap repeats by <50 bp in any of the relevant positions"""

    # make the outdir
    fun.delete_folder(outdir)
    fun.make_folder(outdir)

    # get the chromosome to len
    chr_to_len = fun.get_chr_to_len(reference_genome)

    # get the repeats df of simple SVs
    repeats_df = fun.get_tab_as_df_or_empty_df(repeats_file)
    repeats_df = repeats_df[repeats_df["type"].isin({"Low_complexity", "Simple_repeat"})]
    simple_repeats_table = "%s/simple_repeats.tab"%outdir
    fun.save_df_as_tab(repeats_df, simple_repeats_table)

    # init modified sv_dict
    filt_sv_dict = {}

    # go through each svfile
    for svtype, svfile in sv_dict.items():

        # define the name of the
        svfile_filt = "%s/%s_noSimpleRepeats.tab"%(outdir, svtype)

        # load the svDF and add an ID
        svDF = fun.get_tab_as_df_or_empty_df(svfile).reset_index(drop=True)
        initial_fields = cp.deepcopy(list(svDF.keys()))
        svDF["varID"] = svtype + "_" + pd.Series(list(range(len(svDF))), index=svDF.index).apply(str)

        # if there is something
        if len(svDF)>0:

            # create an svdict that maps different parts of the chroms
            if svtype=="translocations": 
                svDF["is_inverted"] = svDF.StartB>0
                print("There are %i/%i inverted balanced translocations"%(sum(svDF.is_inverted), len(svDF)))                
                inner_svtype_to_svDF = {"translocations_INV":svDF[svDF.is_inverted], "translocations_noINV":svDF[~(svDF.is_inverted)]}

            else: inner_svtype_to_svDF = {svtype : svDF}

            # create a vcf_df with all the regions under SV
            vcf_df = pd.DataFrame()
            for inner_svtype, inner_svDF in inner_svtype_to_svDF.items():
                for pos_field, chrom_field in svtype_to_positionField_to_chromosome[inner_svtype].items():
                    vcf_df = vcf_df.append(inner_svDF[[chrom_field, pos_field, "varID"]].rename(columns={chrom_field:"#CHROM", pos_field:"POS"}))

            # add the ID (to be 1-based)
            vcf_df["POS"] = vcf_df.POS.apply(int)+1
            vcf_df = vcf_df.sort_values(by=["#CHROM", "POS"]).reset_index(drop=True)
            vcf_df["ID"] = "region_" + pd.Series(list(range(len(vcf_df))), index=vcf_df.index).apply(str)
            vcf_df = vcf_df[["#CHROM", "POS", "ID", "varID"]].set_index("ID", drop=False)

            # get a vcf_df with the CHROM and POS of svtype_to_positionField_to_chromosome
            vcf_df["overlaps_repeats"] = fun.get_series_variant_in_repeats(vcf_df, simple_repeats_table, replace=False)
            print("There are %i/%i variant regions overlapping repeats"%(sum(vcf_df.overlaps_repeats), len(vcf_df)))

            # keep only svs that don't overlap repets
            varIDs_overlap_repeats = set(vcf_df[vcf_df.overlaps_repeats].varID)
            svDF_filt = svDF[~(svDF.varID.isin(varIDs_overlap_repeats))]
            print("There are %i/%i SVs that don't overlap repeats"%(len(svDF_filt), len(svDF)))

            fun.save_df_as_tab(svDF_filt[initial_fields], svfile_filt)

        # else
        else: fun.soft_link_files(svfile, svfile_filt)

        # keep the filtered SVs
        filt_sv_dict[svtype] = svfile_filt

    return filt_sv_dict