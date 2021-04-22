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
fun.printing_verbose_mode = True

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


# I downloaded the raw reads from the webpage https://www.ebi.ac.uk/ena/browser/view/ERR194147, because the fastqdump did not work

(reference genome h19)

HG002:

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam # bam

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam.bai

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

# define the reference genomes
hg19_genome = test_fun.get_correct_human_genome("%s/hg19.fa"%DataDir, type_genome="hg19")
hg38_genome = test_fun.get_correct_human_genome("%s/hg38.fa"%DataDir, type_genome="hg38")

# get the reads for CHM1 and CHM13 and merge them
CHM1_13_raw_reads_dir =  "%s/CHM1_13_raw_reads"%DataDir; fun.make_folder(CHM1_13_raw_reads_dir)
srr_to_reads = {srr : test_fun.run_parallelFastqDump_fromSRR_pairedReads_localComputer(srr,  "%s/%s"%(CHM1_13_raw_reads_dir, srr), replace=False, threads=threads) for srr in ["ERR1341794", "ERR1341795"]}

CHM_r1 = "%s/merged_CHM_reads_1.fastq.gz"%CHM1_13_raw_reads_dir
CHM_r2 = "%s/merged_CHM_reads_2.fastq.gz"%CHM1_13_raw_reads_dir
test_fun.merge_reads_into_one_file(srr_to_reads["ERR1341794"][0], srr_to_reads["ERR1341795"][0], CHM_r1, replace=False)
test_fun.merge_reads_into_one_file(srr_to_reads["ERR1341794"][1], srr_to_reads["ERR1341795"][1], CHM_r2, replace=False)

# get the reads from NA12878
NA12878_raw_reads_dir = "%s/NA12878_raw_reads"%DataDir # downloaded manually from ERR194147
NA12878_r1 = "%s/ERR194147_1.fastq.gz"%NA12878_raw_reads_dir
NA12878_r2 = "%s/ERR194147_2.fastq.gz"%NA12878_raw_reads_dir

# get the reads from HG002
HG002_bam = "%s/HG002.hs37d5.60x.1.bam"%DataDir
HG002_r1, HG002_r2 = test_fun.get_fastqgz_from_bam(HG002_bam, threads=threads, replace=False, already_sorted_by_readName=False)

###############################################

######### run SV calling  ###########

# define the type of run
run_in_cluster = True
if running_in_cluster is False: run_in_cluster = False

# we will run SV calling on each referenge genome for the three datasets (using as "real SVs" those found in all of them)
for genome_name, genome, mitochondrial_chromosome in [("hg38", hg38_genome, "chrM"), ("hg19", hg19_genome, "chrMT")]:
    print(genome_name)

    # define the outdir
    outdir_perSVade = "%s/running_on_%s"%(outdir_testing, genome_name); fun.make_folder(outdir_perSVade)

    # get the repeats for this genome
    previous_repeats_table = fun.get_repeat_maskerDF(genome, threads=threads, replace=False, use_repeat_modeller=False)[1]

    # get homologous regions that are exactly the same (this is to generate simulate_SVs_arround_HomologousRegions_previousBlastnFile)
    simulate_SVs_arround_HomologousRegions_maxEvalue = 0.00001
    simulate_SVs_arround_HomologousRegions_queryWindowSize = 500

    windows_multifasta = fun.get_multifasta_genome_split_into_windows(genome, simulate_SVs_arround_HomologousRegions_queryWindowSize, False, 10000000)

    simulate_SVs_arround_HomologousRegions_previousBlastnFile = test_fun.get_blastn_regions_genome_against_itself_equal_regions(windows_multifasta, replace=False, threads=threads, window_size=simulate_SVs_arround_HomologousRegions_queryWindowSize)

    # define the table with short reads
    close_shortReads_table_df = pd.DataFrame({Is : {"sampleID":sampleID, "runID":sampleID+"run1", "short_reads1":reads[0], "short_reads2":reads[1]} for Is, (sampleID, reads) in enumerate({"CHM":(CHM_r1, CHM_r2), "NA12878":(NA12878_r1, NA12878_r2),  "HG002":(HG002_r1, HG002_r2)}.items())}).transpose()
    close_shortReads_table = "%s/table_reads_genome_%s.tab"%(DataDir, genome_name)
    fun.save_df_as_tab(close_shortReads_table_df, close_shortReads_table)


    # define the simulation ploidies as diploid
    simulation_ploidies = "diploid_hetero"

    # define the cmd
    cmd = "%s --ref %s --threads %i -o %s --close_shortReads_table %s --n_close_samples 3 --nruns_per_sample 1 -f1 skip -f2 skip --mitochondrial_chromosome %s --testAccuracy --verbose --nsimulations 2 --skip_CNV_calling --simulation_ploidies %s --previous_repeats_table %s --StopAfter_testAccuracy --simulate_SVs_arround_HomologousRegions_maxEvalue %.10f --simulate_SVs_arround_HomologousRegions_queryWindowSize %i --simulate_SVs_arround_HomologousRegions_previousBlastnFile %s --StopAfter_testAccuracy_perSVadeRunning --skip_SV_CNV_calling --min_gb_RAM_required 3"%(perSVade_py, genome, threads, outdir_perSVade, close_shortReads_table, mitochondrial_chromosome, simulation_ploidies, previous_repeats_table, simulate_SVs_arround_HomologousRegions_maxEvalue, simulate_SVs_arround_HomologousRegions_queryWindowSize, simulate_SVs_arround_HomologousRegions_previousBlastnFile)

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
        name = genome_name
     
        # run jobs
        if cluster_name=="MN4": 

            queue = "bsc_ls"
            time = "48:00:00"
            nodes = 3

            fun.run_jobarray_file_MN4_greasy(jobs_filename, name, time=time, queue=queue, threads_per_job=threads, nodes=nodes)

        elif cluster_name=="Nord3": 

            # define resources
            time = "48:00:00"
            RAM_per_thread = 6000 # first 3600

            # run
            fun.run_jobarray_file_Nord3(jobs_filename, name, time=time, queue="bsc_ls", threads_per_job=threads, RAM_per_thread=RAM_per_thread, max_njobs_to_run=10000)


    elif len(all_lines_jobfile)!=0: raise ValueError("something went wrong")

    ############################

    #if taxID=="5476": adkjhdakg # stop after C. albicans

print("the testing of several species finsihed susccesffully")

#####################################



















