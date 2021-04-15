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
#if run_in_cluster is False: matplotlib.use('TkAgg')

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

"""
species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", 10000000000000000),
                ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314", 10000000000000000),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000),
                ("746128", "Aspergillus_fumigatus", 1, "CM016889.1", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]
"""

"""
species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", 10000000000000000),
                ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314", 10000000000000000),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]
"""

#species_Info = [("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]



species_Info = [("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]

"""
species_Info = [("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30)]
"""

#species_Info = [("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000)]


"""
species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", 10000000000000000),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30)]

"""
"""
species_Info = [("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]

"""

"""
species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", 10000000000000000),
                ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314", 10000000000000000),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000)]

"""

"""
species_Info = [("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]
"""

"""
species_Info = [("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000),
                ("746128", "Aspergillus_fumigatus", 1, "CM016889.1", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]
"""
"""
species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", 10000000000000000),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000)]
"""
"""
species_Info = [
                ("746128", "Aspergillus_fumigatus", 1, "CM016889.1", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30)
                ]
"""

#species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", 10000000000000000)]
#species_Info = [("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000)]
#species_Info = [("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314", 10000000000000000)]

#species_Info = [("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30)]

#species_Info = [("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]

"""
species_Info = [("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]
"""

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


def keep_STDfiles_nord3Runs_testingAccuracy(all_STDs_dir, outdir_perSVade, spName):

    """This function records all the new files from outdir_perSVade testing accuracy"""

    fun.make_folder(all_STDs_dir)

    # define the dirs
    testing_Accuracy_dir = "%s/testing_Accuracy"%outdir_perSVade
    stddir = "%s/STDfiles"%testing_Accuracy_dir

    # get the current STD file metadata. This indicates the run
    stdout_report = "%s/%s_jobs_stdout.txt"%(stddir, spName)
    stderr_report = "%s/%s_jobs_stderr.txt"%(stddir, spName)

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

def get_df_accuracy_of_parameters_on_test_samples(parameters_df, test_df, outdir, replace=False, threads=4, threads_integration=48):

    """This function tests the accuracy of each of the parameters df on the simulations of test_df. It runs each comparison in a sepparate perSVade job in the cluster"""

    if replace is True: fun.delete_folder(outdir)
    fun.make_folder(outdir)

    # define the metadata of each df
    parameters_df_metadata = [k for k in parameters_df.keys() if k not in {"parameters_json"}]
    test_df_metadata = [k for k in test_df.keys() if k not in {"sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix"}]

    # define the outdir
    outdir_cross_benchmark_files = "%s/tmp_files"%outdir; fun.make_folder(outdir_cross_benchmark_files)

    # define the benchmarking file
    df_benchmark_all_file = "%s/cross_benchmarking_parameters.tab"%outdir

    if fun.file_is_empty(df_benchmark_all_file):

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
            outdir_parms = "%s/parms_%s"%(outdir_cross_benchmark_files, "_".join(unique_parms_row[parameters_df_metadata]))
            fun.make_folder(outdir_parms)

            for numeric_test_index, (Itest, test_row) in enumerate(test_df.iterrows()):
            
                # define the final df
                df_benchmarking_file = "%s/testON_%s_benchmarking_df.tab"%(outdir_parms, "_".join(test_row[test_df_metadata]))

                # define cmd
                testing_script = "%s/scripts/perSVade/perSVade_repository/testing/get_accuracy_parameters_on_sorted_bam.py"%ParentDir

                if fun.file_is_empty(df_benchmarking_file): 

                    all_cmds.append("%s --reference_genome %s --df_benchmarking_file %s --sorted_bam %s --gridss_vcf %s --svtables_prefix %s --threads %i --parameters_json %s --verbose --mitochondrial_chromosome %s"%(testing_script, test_row.reference_genome, df_benchmarking_file, test_row.sorted_bam, test_row.gridss_vcf, test_row.svtables_prefix, threads, unique_parms_row.parameters_json, test_row.mitochondrial_chromosome))

                    continue

              
                # append
                list_dfs_benchmark_files.append((Idf, df_benchmarking_file, parms_row, test_row, parameters_df_metadata, test_df_metadata, all_keys))
                Idf += 1

        # get unique cmds
        all_cmds_unique = sorted(set(all_cmds))
        print("There are %i jobs to run. These are %i unique jobs"%(len(all_cmds), len(all_cmds_unique)))

        # run cmds
        if run_in_cluster is False: 
            for cmd in all_cmds_unique: fun.run_cmd(cmd)

        # run in the cluster
        elif len(all_cmds_unique)>0:

            # get jobs file
            print("submitting %i cmds to the cluster"%len(all_cmds_unique))
            jobs_filename = "%s/jobs.getting_crossAccuracy"%outdir
            open(jobs_filename, "w").write("\n".join(all_cmds_unique))
            fun.generate_jobarray_file(jobs_filename, "gettingCloseShortReads")

            # submit to the cluster
            #fun.run_jobarray_file_MN4_greasy(jobs_filename, "getting_crossAccuracy", time="02:00:00", queue="debug", threads_per_job=threads, nodes=8)
            fun.run_jobarray_file_MN4_greasy(jobs_filename, "getting_crossAccuracy", time="10:00:00", queue="bsc_ls", threads_per_job=threads, nodes=3)

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


    youshouldremaketheplotssothattheyconsiderTheFactThatPrecision_and_recall_are_calculated_differently

    # filter the df
    df_plot = df[(df.svtype==svtype) & (df.comparisonID!="perSVade-arroundRepeats") & (df.type_SVs_longReads=="all_SVs") & (df.remaining_treatment=="drop")].set_index("species", drop=False)


    #& (df.threshold_fractionParms==0.4789473684210527) &(df.tol_bp==50) & (df.pct_overlap==0.75)

    # define vars
    typeRun_to_color = {"perSVade-uniform":"blue", "perSVade-fast":"gray", "perSVade-realSVs":"red", 'perSVade-arroundHomRegions':"black"}
    sorted_typeRuns  = ["perSVade-arroundHomRegions", "perSVade-uniform", "perSVade-realSVs"]
    #sorted_species =  ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana", "Drosophila_melanogaster"]
    sorted_species =  ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana"]

    typeRun_to_typeRunText = {"perSVade-realSVs":"known SVs", "perSVade-uniform":"random", "perSVade-arroundHomRegions":"homologous SVs"}

    # plot init
    nrows = len(sorted_species)
    ncols = len(sorted_typeRuns)
    fig = plt.figure(figsize=(ncols*1, nrows*1)); I = 1


    for Is, species in enumerate(sorted_species):

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
            if Is==0: ax.set_title("%s"%(typeRun.split("-")[1]))
            ax.set_ylim(ylim)
            ax.set_xlim([-0.45, 1.45])

            if Is==(nrows-1): 
                ax.set_xticks([0,1])
                ax.set_xticklabels(["default parms.", typeRun_to_typeRunText[typeRun]])
                for label in ax.get_xticklabels(): label.set_rotation(90)

            ax.set_xlabel("")
            if Ir!=0:
                ax.set_yticklabels([])
                ax.set_ylabel("")

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


def get_cross_accuracy_df_several_perSVadeSimulations(outdir_testing, genomes_and_annotations_dir, replace=False):

    """This function tests how each of the perSVade configurations works on the others. It runs one job for each type of simulations, and it iterates through them inside of the job."""


    # define the final outdir
    outdir_cross_accuracy = "%s/cross_accuracy_calculations"%outdir_testing
    df_benchmark_all_file = "%s/cross_benchmarking_parameters.tab"%outdir_cross_accuracy

    if fun.file_is_empty(df_benchmark_all_file) or replace is True:

        ###### GET PARAMETERS DF AND TESTING DF #######

        # the parameters_df. The first cols are metadata (like sampleID, runID and optimisation type) and the others are things necessary for runnning gridss: and the path to the parameters_json
        parameters_df_dict = {}

        # test_df: This is info on which to test the running of gridss+clove. It contains metadata cols (sampleID, runID, optimisation type (real, uniform), simName, ploidy, svtype) and data to run the optimisation on (sorted_bam, gridss_vcf, reference_genome, mitochondrial_chromosome)
        test_df_dict = {}

        for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in species_Info:
            for typeSimulations in ["arroundHomRegions", "arroundRepeats", "uniform", "realSVs"]:

                # define outir
                outdir_species_simulations = "%s/%s_%s/testing_Accuracy/%s"%(outdir_testing, taxID, spName, typeSimulations)

                # define samples and runs
                all_sampleIDs = [x for x in os.listdir(outdir_species_simulations)]

                # go through each sampleID
                for sampleID in os.listdir(outdir_species_simulations):
                 
                    # define the outdir of the run
                    sampleID_simulations_files = "%s/%s/simulations_files_and_parameters"%(outdir_species_simulations, sampleID)

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
                        reference_genome = "%s/%s.fasta"%(genomes_and_annotations_dir, spName)
                        svtables_prefix =  "%s/SVs_%s"%(sampleID_simulations_files, simName)

                        # add to dict
                        test_df_dict[(spName, sampleID, typeSimulations, simName, simulation_ploidy)] = {"species":spName, "sampleID":sampleID, "typeSimulations":typeSimulations, "simName":simName, "simulation_ploidy":simulation_ploidy, "sorted_bam":sorted_bam, "gridss_vcf":gridss_vcf, "reference_genome":reference_genome, "mitochondrial_chromosome":mitochondrial_chromosome, "svtables_prefix":svtables_prefix}


        # add the fast parameters
        parameters_json_fast = "%s/5478_Candida_glabrata/testing_Accuracy/fast/BG2_YPD/simulations_files_and_parameters/final_parameters.json"%outdir_testing
        parameters_df_dict[("none", "fast", "fast")] = {"species":"none", "typeSimulations":"fast", "sampleID":"fast", "parameters_json":parameters_json_fast}

        # get the dfs
        parameters_df = pd.DataFrame(parameters_df_dict).transpose()[["species", "sampleID", "typeSimulations", "parameters_json"]]
        test_df = pd.DataFrame(test_df_dict).transpose()[["species", "sampleID", "typeSimulations", "simName", "simulation_ploidy", "sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix"]]

        # plot the cross-accuracy
        print("getting cross-accuracy df")
        df_cross_accuracy_benchmark = get_df_accuracy_of_parameters_on_test_samples(parameters_df, test_df, outdir_cross_accuracy, replace=replace)

    df_cross_accuracy_benchmark = fun.get_tab_as_df_or_empty_df(df_benchmark_all_file)

    ##################################

    ######### ADD FIELDS TO df_cross_accuracy_benchmark #########
    print("adding extra fields")


    # add sample and run IDs
    def get_run(sampleID):
        
        if sampleID=="fast": return "fast"
        else: return sampleID.split("_")[1]

    def get_sample(sampleID): return sampleID.split("_")[0]

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




    #############################################################

    return df_cross_accuracy_benchmark



def generate_heatmap_accuracy_of_parameters_on_test_samples(df_benchmark, fileprefix, replace=False, threads=4, accuracy_f="Fvalue", svtype="integrated", col_cluster = False, row_cluster = False, show_only_species_and_simType=False):

    """
    This function takes a df where each row is one set of training parameters and test data svtype, together with the accuracy records. It generates a heatmap were the rows are each of the training parameters and the cols are the test samples.
    """

    print("plotting cross-accuracy")

    # define graphics
    #species_to_color = {'none': 'gray', 'Drosophila_melanogaster': 'black', 'Arabidopsis_thaliana': 'gray', 'Cryptococcus_neoformans': 'lightcoral', 'Candida_albicans': 'blue', 'Candida_glabrata': 'cyan'}

    species_to_color = {'none': 'gray', 'Drosophila_melanogaster': 'darkorange', 'Arabidopsis_thaliana': 'olive', 'Cryptococcus_neoformans': 'lightcoral', 'Candida_albicans': 'magenta', 'Candida_glabrata': 'lightseagreen'}


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
    df_plot = df_benchmark[(df_benchmark.svtype==svtype) & (df_benchmark.parms_species!="Drosophila_melanogaster") & (df_benchmark.test_species!="Drosophila_melanogaster")] # hide drosophila
    #df_plot = df_benchmark[(df_benchmark.svtype==svtype)]


    # define square df
    parms_keys = ["parms_species", "parms_typeSimulations", "parms_numeric_sample", "parms_numeric_run"]
    test_keys = ["test_species", "test_typeSimulations", "test_numeric_sample", "test_numeric_run", "test_simName"]
    df_plot["parms_idx"] = df_plot.apply(lambda r: "||||".join([str(r[k]) for k in parms_keys]), axis=1)
    df_plot["test_idx"] = df_plot.apply(lambda r: "||||".join([str(r[k]) for k in test_keys]), axis=1)
    df_square = df_plot[["parms_idx", "test_idx", accuracy_f]].pivot(index='parms_idx', columns='test_idx', values=accuracy_f)

    # sort by species
    print("sorting")
    species_to_order =  {'none': 0, "Candida_glabrata":1, "Candida_albicans":2, "Cryptococcus_neoformans":3, "Arabidopsis_thaliana":4, "Drosophila_melanogaster":5}
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


    # get the plot
    filename = "%s_cross_accuracy_%s_%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype, col_cluster, row_cluster)
    print("getting %s"%filename)

    # define the figure size
    figsize = (int(len(df_square.columns)*0.03), int(len(df_square)*0.03))

    fun.plot_clustermap_with_annotation(df_square, row_colors_df, col_colors_df, filename, title="cross accuracy", col_cluster=col_cluster, row_cluster=row_cluster, colorbar_label=accuracy_f, adjust_position=True, legend=True, idxs_separator_pattern="||||", texts_to_strip={"L001"}, default_label_legend="control", df_annotations=df_annotations, cmap=sns.color_palette("RdBu_r", 50), ylabels_graphics_df=None, grid_lines=False, figsize=figsize, multiplier_width_colorbars=5)




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



def get_crossaccuracy_distributions(df_cross_accuracy_benchmark, fileprefix, accuracy_f="Fvalue", svtype="integrated"):

    """Prints a boxplot for each species and the accuracy on different types of simulations depending on the training parameters. Each row is a tested species and the X is the type of parameters"""

    # keep one df
    df_cross_accuracy_benchmark = df_cross_accuracy_benchmark[df_cross_accuracy_benchmark.svtype==svtype]

    # define parms
    #sorted_species = ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana", "Drosophila_melanogaster"]
    sorted_species = ["Candida_glabrata", "Candida_albicans", "Cryptococcus_neoformans", "Arabidopsis_thaliana"]
    typeSimulations_to_color = {"uniform":"blue", "realSVs":"red", "arroundHomRegions":"black", "fast":"gray"}
    type_comparison_to_xvalue = {"fast":0, "same_run_and_simulation":1, "same_species_and_simulation":1, "same_species":1, "different_species":2}
    xvalue_to_typeCompLabel = {0:"default parms.", 1:"same species", 2:"different species"}
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
        ax.set_ylabel("%s\n%s"%(test_species, accuracy_f))

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
        for x in [0.5, 1.5]: plt.axvline(x, linewidth=.9, color="gray", linestyle="--")


        # get the legen only in the first box
        ax.legend(bbox_to_anchor=(1, 1))
        if I!=0: ax.get_legend().remove()

    # spaces
    plt.subplots_adjust(wspace=0.01, hspace=0.01)
    filename = "%s_%s_%s.pdf"%(fileprefix, accuracy_f, svtype)
    fig.savefig(filename, bbox_inches='tight')
    #plt.close(fig)


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


