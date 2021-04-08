#!/usr/bin/env python

# this script is to evaluate how perSVade works on several human datasets as in cameron 2019

##### DEFINE ENVIRONMENT #######

# module imports
import os
import sys
import pandas as pd

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    running_in_cluster = False    
    threads = 8
else:
    running_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
        
# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions
print("importing functions")
import sv_functions as fun
fun.printing_verbose_mode = True

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
outdir_testing = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_humanGoldenSet"%ParentDir; fun.make_folder(outdir_testing)
CurDir = "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir
DataDir = "%s/data"%outdir_testing
outdir_genomes_and_annotations = "%s/scripts/perSVade/perSVade_repository/testing/genomes_and_annotations"%ParentDir

################################

######### GENERATING INPUT DATA ###############


"""
This is how I obtained the datasets:

NA12878:

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed; mv Personalis_1000_Genomes_deduplicated_deletions.bed NA12878_deletions.bed # SVs

(reference genome h19)

HG002:

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam # bam

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz; gunzip HG002_SVs_Tier1_v0.6.vcf.gz # SVs

(reference genome h19)

CHM1, CHM13:

wget http://eichlerlab.gs.washington.edu/publications/Huddleston2016/structural_variants/CHM13_SVs.annotated.vcf.gz; gunzip CHM13_SVs.annotated.vcf.gz --> 20470
wget http://eichlerlab.gs.washington.edu/publications/Huddleston2016/structural_variants/CHM1_SVs.annotated.vcf.gz; gunzip CHM1_SVs.annotated.vcf.gz --> 20602

(reference genome h38)
(these should be merged into a single callset)

I download the reference genomes from the latest UCSC (at 06/04/2021):

ref. genome h19: 

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz; gunzip hg19.fa.gz

ref. genome h38:

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz; gunzip hg38.fa.gz

"""

# get the reads for CHM1 and CHM13 and merge them
CHM1_13_raw_reads_dir =  "%s/CHM1_13_raw_reads"%DataDir; fun.make_folder(CHM1_13_raw_reads_dir)
for srr in ["ERR1341794", "ERR1341795"]: fun.run_parallelFastqDump_fromSRR_pairedReads_localComputer(srr,  "%s/%s"%(CHM1_13_raw_reads_dir, srr), replace=False, threads=threads)


akhagjdagdag

"""
CHM_r1
CHM_r2
adjhgadjhgda
"""

# get the reads from NA12878
NA12878_raw_reads_dir = "%s/NA12878_raw_reads"%DataDir
srr = "ERR194147"
fun.run_parallelFastqDump_fromSRR_pairedReads(srr, NA12878_raw_reads_dir, replace=False, threads=threads)


kadhgahdgjadg
"""
NA12878_r1
NA12878_r2
adjhgadjgjad
"""


# get the reads from HG002
HG002_bam = "%s/HG002.hs37d5.60x.1.bam"%DataDir
HG002_r1, HG002_r2 = test_fun.get_fastqgz_from_bam(HG002_bam, threads=threads, replace=False)


adadgadjgjhadg

# get repeatmasker for each genome, without repeatmodeller

###############################################

############## TRIM READS ##################




#############################################

















# define the taxIDs that have no
taxIDs_with_noON_overalpping = {"746128"}

# define the run in cluster (and debug)
run_in_cluster = True
if running_in_cluster is False: run_in_cluster = False

# go through each species
for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in test_fun.species_Info:
    print(taxID, spName)

    # define  the genome and annotations
    genome = "%s/%s.fasta"%(outdir_genomes_and_annotations, spName)
    gff = "%s/%s.gff"%(outdir_genomes_and_annotations, spName)

    # create an outdir
    outdir_perSVade = "%s/%s_%s"%(outdir_testing, taxID, spName); fun.make_folder(outdir_perSVade)

    # get the repeats for this genome
    previous_repeats_table = fun.get_repeat_maskerDF(genome, threads=threads, replace=False)[1]

    # get the blastn of the genome against itself
    simulate_SVs_arround_HomologousRegions_maxEvalue = 0.00001
    simulate_SVs_arround_HomologousRegions_queryWindowSize = 500
    simulate_SVs_arround_HomologousRegions_previousBlastnFile = fun.get_blastn_regions_genome_against_itself(genome, simulate_SVs_arround_HomologousRegions_maxEvalue, simulate_SVs_arround_HomologousRegions_queryWindowSize, False, threads)

    # define the goldenSet_dir and the real bedpe breakpoints
    if spName=="Candida_glabrata": 
        goldenSet_table = test_fun.get_goldenSetTable_Cglabrata(CurDir)
        real_bedpe_breakpoints = "%s/outdirs_testing_severalSpecies/%s_%s/findingRealSVs_providedCloseReads/integrated_breakpoints.bedpe"%(CurDir, taxID, spName)

    else: 
        goldenSet_table = "auto"
        real_bedpe_breakpoints = "%s/outdirs_testing_severalSpecies/%s_%s/findingRealSVs_automaticFindingOfCloseReads/integrated_breakpoints.bedpe"%(CurDir, taxID, spName)

    # define the ploidy
    if ploidy==1: simulation_ploidies = "haploid"
    elif ploidy==2: simulation_ploidies = "diploid_hetero"
    else: raise ValueError("ploidy %i is not valid"%ploidy)

    # get the golden set running 
    if taxID in taxIDs_with_noON_overalpping: continue
    cmd = "%s --ref %s --threads %i -o %s --target_taxID %s --real_bedpe_breakpoints %s -f1 skip -f2 skip --mitochondrial_chromosome %s --gff %s --goldenSet_table %s --skip_SVcalling --verbose --nsimulations 2 --simulation_ploidies %s --previous_repeats_table %s --QC_and_trimming_reads --StopAfter_goldenSetAnalysis --simulate_SVs_arround_HomologousRegions_previousBlastnFile %s --simulate_SVs_arround_HomologousRegions_maxEvalue %.10f --simulate_SVs_arround_HomologousRegions_queryWindowSize %i"%(perSVade_py, genome, threads, outdir_perSVade, taxID, real_bedpe_breakpoints, mitochondrial_chromosome, gff, goldenSet_table, simulation_ploidies, previous_repeats_table, simulate_SVs_arround_HomologousRegions_previousBlastnFile, simulate_SVs_arround_HomologousRegions_maxEvalue, simulate_SVs_arround_HomologousRegions_queryWindowSize)

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
    #fun.run_cmd("%s > %s 2>&1"%(cmd, cmd_output)) # run with stdout
    fun.run_cmd(cmd); continue # run locally 
 	
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
            nodes = 4

            fun.run_jobarray_file_MN4_greasy(jobs_filename, name, time=time, queue=queue, threads_per_job=threads, nodes=nodes)

        elif cluster_name=="Nord3": 

            queue = "bsc_ls"; 
            RAM_per_thread = 4000; # 1800 or 5000 
            time = "48:00:00" # per job

            #fun.run_jobarray_file_Nord3(jobs_filename, name, time=time, queue=queue, threads_per_job=threads, RAM_per_thread=RAM_per_thread, max_njobs_to_run=10000)
            fun.run_jobarray_file_Nord3_greasy(jobs_filename, name, time=time, queue=queue, threads_per_job=threads, RAM_per_thread=RAM_per_thread, nodes=4)

    elif len(all_lines_jobfile)!=0: raise ValueError("something went wrong")

    ############################

print("Golden set analysis worked")
