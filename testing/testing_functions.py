#!/usr/bin/env python

# This script contains functions useful for the testing of perSVade on several samples

######### ENV ########

import os
import sys
import pandas as pd
import numpy as np

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if not os.path.exists(ParentDir): ParentDir = "/gpfs/projects/bsc40/mschikora"
        
# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions from perSVade
print("importing functions")
import sv_functions as fun

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
species_Info = [("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000),
                ("746128", "Aspergillus_fumigatus", 1, "CM016889.1", 10000000000000000),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]
"""

species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", 10000000000000000),
                ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314", 10000000000000000),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1", 10000000000000000)]
                
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
