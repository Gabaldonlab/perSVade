#!/usr/bin/env python

# This is the perSVade pipeline main script, which shoul dbe run on the perSVade conda environment


##### DEFINE ENVIRONMENT #######

# module imports
import argparse, os
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter
import copy as cp
import pickle
import string
import shutil 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import sys
from shutil import copyfile

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import sv_functions as fun

# packages installed into the conda environment 
samtools = "%s/bin/samtools"%EnvDir
java = "%s/bin/java"%EnvDir

# scripts that are installed under this software
varcall_cnv_pipeline = "%s/varcall_cnv_pipeline.py"%CWD
perSVade_genome_browser = "%s/perSVade_genome_browser.py"%CWD

#######

description = """
Runs perSVade pipeline on an input set of paired end short ends. It is expected to be run on a coda environment and have several dependencies (see https://github.com/Gabaldonlab/perSVade). Some of these dependencies are included in the respository "installation/external_software". These are gridss (tested on version 2.8.1), clove (tested on version 0.17) and NINJA (we installed it from https://github.com/TravisWheelerLab/NINJA/releases/tag/0.95-cluster_only). If you have any trouble with these you can replace them from the source code.
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# general args
parser.add_argument("-r", "--ref", dest="ref", required=True, help="Reference genome. Has to end with .fasta.")
parser.add_argument("-thr", "--threads", dest="threads", default=16, type=int, help="Number of threads, Default: 16")
parser.add_argument("-o", "--outdir", dest="outdir", action="store", required=True, help="Directory where the data will be stored")
parser.add_argument("--replace", dest="replace", action="store_true", help="Replace existing files")
parser.add_argument("-p", "--ploidy", dest="ploidy", default=1, type=int, help="Ploidy, can be 1 or 2")


# replace CNV_calling
parser.add_argument("--replace_SV_CNVcalling_and_optimisation", dest="replace_SV_CNVcalling_and_optimisation", action="store_true", help="Replace everything related to the SV and CNV calling.")

parser.add_argument("--replace_only_SV_CNVcalling", dest="replace_only_SV_CNVcalling", action="store_true", help="Replace the running of SV_CNVcalling, but not the coverage thresholds")

parser.add_argument("--replace_FromGridssRun_final_perSVade_run", dest="replace_FromGridssRun_final_perSVade_run", action="store_true", help="Replace from the clove running in the final gridss+clove running")


# different modules to be executed
parser.add_argument("--testSVgen_from_DefaulReads", dest="testSVgen_from_DefaulReads", default=False, action="store_true", help="This indicates whether to generate a report of how the generation of SV works with the default parameters on random simulations")

parser.add_argument("--close_shortReads_table", dest="close_shortReads_table", type=str, default=None, help="This is the path to a table that has 4 fields: sampleID,runID,short_reads1,short_reads2. These should be WGS runs of samples that are close to the reference genome and some expected SV. Whenever this argument is provided, the pipeline will find SVs in these samples and generate a folder <outdir>/findingRealSVs<>/SVs_compatible_to_insert that will contain one file for each SV, so that they are compatible and ready to insert in a simulated genome. This table will be used if --testAccuracy is specified, which will require at least 3 runs for each sample. It can be 'auto', in which case it will be inferred from the taxID provided by --target_taxID.")

parser.add_argument("--target_taxID", dest="target_taxID", type=int, default=None, help="This is the taxID (according to NCBI taxonomy) to which your reference genome belongs. If provided it is used to download genomes and reads.")

parser.add_argument("--n_close_samples", dest="n_close_samples", default=5, type=int, help="Number of close samples to search in case --target_taxID is provided")

parser.add_argument("--nruns_per_sample", dest="nruns_per_sample", default=3, type=int, help="Number of runs to download for each sample in the case that --target_taxID is specified. ")

parser.add_argument("--real_bedpe_breakpoints", dest="real_bedpe_breakpoints", type=str, default=None, help="A file with the list of 'real' breakpoints arround which to insert the SVs in simulations. It may be created with --close_shortReads_table. If both --real_bedpe_breakpoints and --close_shortReads_table are provided, --real_bedpe_breakpoints will be used, and --close_shortReads_table will have no effect. If none of them are provided, this pipeline will base the parameter optimization on randomly inserted SVs (the default behavior). The coordinates have to be 1-based, as they are ready to insert into RSVsim.")

parser.add_argument("--parameters_json_file", dest="parameters_json_file", type=str, default=None, help="A file with the json parameters to use. This only has effect if --fast_SVcalling is specified")


parser.add_argument("--fast_SVcalling", dest="fast_SVcalling", action="store_true", default=False, help="Run SV calling with a default set of parameters. There will not be any optimisation nor reporting of accuracy. This is expected to work almost as fast as gridss and clove together. If --parameters_json_file, the parameters are substituted by the json parameters.")


# set the minimum SVsize
parser.add_argument("--min_CNVsize_betweenBPs", dest="min_CNVsize_betweenBPs", default=500, type=int, help="The minimum size of a CNV inferred brtween breakpoints in the perSVade running.")

# pipeline skipping options 
parser.add_argument("--skip_SVcalling", dest="skip_SVcalling", action="store_true", default=False, help="Do not run SV calling.")

parser.add_argument("--skip_SV_CNV_calling", dest="skip_SV_CNV_calling", action="store_true", default=False, help="Do not run the SV and CNV calling.")


# options of the long reads-based benchmarking
parser.add_argument("--goldenSet_dir", dest="goldenSet_dir", type=str, default=None, help="This is the path to a directory that has some oxford nanopore reads (should end with'long_reads.fasta') and some  short paired end reads (ending with '1.fastq.gz' and '2.fastq.gz'). These are assumed to be from the exact same sample. If provided, perSVade will call SVs from it using the --nanopore configuration of svim and validate each of the 'real' (if provided), 'uniform' and 'fast' versions from the short reads on it. If you state 'auto', it will look for samples of your --target_taxID in the SRA that are suited. We already provide some automated finding of reads in the SRA for several taxIDs: 3702 (Arabidopsis_thaliana). All the jobs will be run in an squeue as specified by --job_array_mode.")

# pipeline stopping options
parser.add_argument("--StopAfter_readObtentionFromSRA", dest="StopAfter_readObtentionFromSRA", action="store_true", default=False, help="Stop after obtaining reads from SRA.")
parser.add_argument("--StopAfter_sampleIndexingFromSRA", dest="StopAfter_sampleIndexingFromSRA", action="store_true", default=False, help="It will stop after indexing the samples of SRA. You can use this if, for example, your local machine has internet connection and your slurm cluster does not. You can first obtain the SRA indexes in the local machine. And then run again this pipeline without this option in the slurm cluster.")
parser.add_argument("--StopAfter_genomeObtention", dest="StopAfter_genomeObtention", action="store_true", default=False, help="Stop after genome obtention.")
parser.add_argument("--StopAfter_bamFileObtention", dest="StopAfter_bamFileObtention", action="store_true", default=False, help="Stop after obtaining the BAM file of aligned reads.")
parser.add_argument("--StopAfterPrefecth_of_reads", dest="StopAfterPrefecth_of_reads", action="store_true", default=False, help="Stop after obtaining the prefetched .srr file in case close_shortReads_table is 'auto'")
parser.add_argument("--StopAfterPrefecth_of_reads_goldenSet", dest="StopAfterPrefecth_of_reads_goldenSet", action="store_true", default=False, help="Stop after obtaining the prefetched .srr file in case --goldenSet_dir is specified.")
parser.add_argument("--StopAfter_obtentionOFcloseSVs", dest="StopAfter_obtentionOFcloseSVs", action="store_true", default=False, help="Stop after obtaining the real_bedpe_breakpoints ")
parser.add_argument("--StopAfter_repeatsObtention", dest="StopAfter_repeatsObtention", action="store_true", default=False, help="Stop after obtaining  the repeats table")

parser.add_argument("--StopAfter_testAccuracy_perSVadeRunning", dest="StopAfter_testAccuracy_perSVadeRunning", action="store_true", default=False, help="When --testAccuracy is specified, the pipeline will stop after the running of perSVade on all the inputs of --close_shortReads_table with the different configurations.")

# testing options
parser.add_argument("--testAccuracy", dest="testAccuracy", action="store_true", default=False, help="Reports the accuracy  of your calling on the real data, simulations and fastSVcalling for all the WGS runs specified in --close_shortReads_table. ")


# simulation parameter args
parser.add_argument("--nvars", dest="nvars", default=50, type=int, help="Number of variants to simulate for each SVtype.")
parser.add_argument("--nsimulations", dest="nsimulations", default=3, type=int, help="The number of 'replicate' simulations that will be produced.")
parser.add_argument("--simulation_ploidies", dest="simulation_ploidies", type=str, default="auto", help='A comma-sepparated string of the ploidies to simulate for parameter optimisation. It can have any of "haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1" . By default it will be inferred from the ploidy. For example, if you are running on ploidy=2, it will optimise for diploid_hetero and diploid_homo.')
parser.add_argument("--range_filtering_benchmark", dest="range_filtering_benchmark", type=str, default="theoretically_meaningful", help='The range of parameters that should be tested in the SV optimisation pipeline. It can be any of large, medium, small, theoretically_meaningful or single.')

# alignment args
parser.add_argument("-f1", "--fastq1", dest="fastq1", default=None, help="fastq_1 file. Option required to obtain bam files. It can be 'skip'")
parser.add_argument("-f2", "--fastq2", dest="fastq2", default=None, help="fastq_2 file. Option required to obtain bam files. It can be 'skip'")
parser.add_argument("-sbam", "--sortedbam", dest="sortedbam", default=None, help="The path to the sorted bam file, which should have a bam.bai file in the same dir. For example, if your bam file is called 'aligned_reads.bam', there should be an 'aligned_reads.bam.bai' as well. This is mutually exclusive with providing reads. By default, it is assumed that this bam has marked duplicates.")

# machine options
parser.add_argument("--job_array_mode", dest="job_array_mode", type=str, default="local", help="It specifies in how to run the job arrays for,  --testAccuracy, the downloading of reads if  --close_shortReads_table is auto, and the SV calling for the table in --close_shortReads_table. It can be 'local' (runs one job after the other or 'job_array'. If 'job_array' is specified, this pipeline will generate a file with all the jobs to run, and stop. You will have to run these jobs before stepping into downstream analyses.")


# other args
parser.add_argument("-mchr", "--mitochondrial_chromosome", dest="mitochondrial_chromosome", default="mito_C_glabrata_CBS138", type=str, help="The name of the mitochondrial chromosome. This is important if you have mitochondrial proteins for which to annotate the impact of nonsynonymous variants, as the mitochondrial genetic code is different. This should be the same as in the gff. If there is no mitochondria just put 'no_mitochondria'. If there is more than one mitochindrial scaffold, provide them as comma-sepparated IDs.")

# do not clean the outdir
parser.add_argument("--skip_cleaning_outdir", dest="skip_cleaning_outdir", action="store_true", default=False, help="Will NOT remove all the unnecessary files of the perSVade outdir")
parser.add_argument("--skip_cleaning_simulations_files_and_parameters", dest="skip_cleaning_simulations_files_and_parameters", action="store_true", default=False, help="Will NOT remove all the files related to simulations and parameters from the <outdir>/simulations_files_and_parameters")

# arg to run the trimming of the reads
parser.add_argument("--QC_and_trimming_reads", dest="QC_and_trimming_reads", action="store_true", default=False, help="Will run fastq and trimmomatic of reads, and use the trimmed reads for downstream analysis. This option will generate files under the same dir as f1 and f2, so be aware of it.")

# maximum coverage args
parser.add_argument("--max_coverage_sra_reads", dest="max_coverage_sra_reads", default=10000000000000000, type=int, help="This is the maximum coverage allowed in the reads downloaded by SRA. If the datasets have a coverage above this perSVade will randmomly subset reads to match this.")

# min chromosome name
parser.add_argument("--min_chromosome_len", dest="min_chromosome_len", default=100000, type=int, help="The minimum length to consider chromosomes from the provided fasta for calculating the window length. Any chromosomes that shorter than the window length will not be considered in the random SV simulations.")

# the fraction of available memory
parser.add_argument("--fraction_available_mem", dest="fraction_available_mem", default=None, help="The fraction of RAM that is being allocated to this perSVade run. In several steps, this pipeline needs to calculate the available memory (using psutil.virtual_memory()). This returns all the available memory in the computer. If you are running on a fraction of the computers' resources, this calculation is overestimating the available RAM. In such case you can provide the fraction available through this argument. By default, it will calculate the available ram by filling the memory, which may give errors.")

# coverage CNV
parser.add_argument("--min_coverage_duplication", dest="min_coverage_duplication", default="auto", help="The minimum rel. coverage to call a duplicated region. If set to auto, it will be calculated from optimisation")
parser.add_argument("--max_coverage_deletion", dest="max_coverage_deletion", default="auto", help="The maximum rel. coverage to call a deleted region. If set to auto, it will be calculated from optimisation")
parser.add_argument("--min_r_pearson_noFlatRegions", dest="min_r_pearson_noFlatRegions", default="auto", help="The minimum pearson correlation that there should be between position and coverage in a region for being considered as a duplication.")
parser.add_argument("--min_r_spearman_noFlatRegions", dest="min_r_spearman_noFlatRegions", default="auto", help="The minimum spearman correlation that there should be between position and coverage in a region for being considered as a duplication.")

# small VarCall and CNV args
parser.add_argument("--run_smallVarsCNV", dest="run_smallVarsCNV", action="store_true", default=False, help="Will call small variants and CNV.")
parser.add_argument("-gff", "--gff-file", dest="gff", default=None, help="path to the GFF3 annotation of the reference genome. Make sure that the IDs are completely unique for each 'gene' tag. This is necessary for both the CNV analysis (it will look at genes there) and the annotation of the variants.")
parser.add_argument("-caller", "--caller", dest="caller", required=False, default="all", help="SNP caller option to obtain vcf file. options: no/all/HaplotypeCaller/bcftools/freebayes. It can be a comma-sepparated string, like 'HaplotypeCaller,freebayes'")
parser.add_argument("-c", "--coverage", dest="coverage", default=20, type=int, help="minimum Coverage (int)")
parser.add_argument("--minAF_smallVars", dest="minAF_smallVars", default="infer", help="The minimum fraction of reads covering a variant to be called. The default is 'infer', which will set a threshold based on the ploidy. This is only relevant for the final vcfs, where only PASS vars are considered. It can be a number between 0 and 1.")

parser.add_argument("-mcode", "--mitochondrial_code", dest="mitochondrial_code", default=3, type=int, help="The code of the NCBI mitochondrial genetic code. For yeasts it is 3. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")
parser.add_argument("-gcode", "--gDNA_code", dest="gDNA_code", default=1, type=int, help="The code of the NCBI gDNA genetic code. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . For C. albicans it is 12. ")
parser.add_argument("--remove_smallVarsCNV_nonEssentialFiles", dest="remove_smallVarsCNV_nonEssentialFiles", action="store_true", default=False, help="Will remove all the varCall files except the integrated final file and the bam file.")
parser.add_argument("--replace_var_integration", dest="replace_var_integration", action="store_true", help="Replace all the variant integration steps for smallVariantCalling.")
parser.add_argument("--generate_alternative_genome", dest="generate_alternative_genome", default=False, action="store_true", help="Generate an alternative genome in smallVariantCalling.")
parser.add_argument("--skip_cnv_analysis", dest="skip_cnv_analysis", default=False, action="store_true", help="Don't perform the cnv analysis")


# add the CNV calling args
parser.add_argument("--window_size_CNVcalling", dest="window_size_CNVcalling", default=3000, type=int, help="The window size in which the genome will be fragmented for CNV calling.")


# visualization
parser.add_argument("--visualization_results", dest="visualization_results", default=False, action="store_true", help="Visualize the results")


parser.add_argument("--pooled_sequencing", dest="pooled_sequencing", action="store_true", default=False, help="It is a pooled sequencing run, which means that the small variant calling is not done based on ploidy. If you are also running SV calling, you will run parameters optimisation on a sample that has 1-10 pooling.")

# verbosity
parser.add_argument("--verbose", dest="verbose", action="store_true", default=False, help="Whether to print a verbose output. This is important if there are any errors in the run.")

# repeat obtention
parser.add_argument("--consider_repeats_smallVarCall", dest="consider_repeats_smallVarCall", action="store_true", default=False, help="If --run_smallVarsCNV, this option will imply that each small  variant will have an annotation of whether it overlaps a repeat region.")
parser.add_argument("--previous_repeats_table", dest="previous_repeats_table", default=None, help="This may be the path to a file that contains the processed output of RepeatMasker (such as the one output by the function get_repeat_maskerDF). This should be a table with the following header: 'SW_score, perc_div, perc_del, perc_ins, chromosome, begin_repeat, end_repeat, left_repeat, strand, repeat, type, position_inRepeat_begin, position_inRepeat_end, left_positionINrepeat, IDrepeat'. It is created by parsing the tabular output of RepeatMasker and putting into a real .tab format.")

parser.add_argument("--skip_repeat_analysis", dest="skip_repeat_analysis", default=False, action="store_true", help="Skip the inference of repeats. If --previous_repeats_table is provided, this argument will override it.")

# small varCall stop options
parser.add_argument("--StopAfter_smallVarCallSimpleRunning", dest="StopAfter_smallVarCallSimpleRunning", action="store_true", default=False, help="Stop after obtaining the filtered vcf outputs of each program.")


opt = parser.parse_args()

########################################
##### GENERAL PROCESSING OF INPUTS #####
########################################

# if replace is set remove the outdir, and then make it
if opt.replace is True: fun.delete_folder(opt.outdir)
fun.make_folder(opt.outdir)

# define the final file. and exit if it exists
final_file = "%s/perSVade_finished_file.txt"%opt.outdir
if not fun.file_is_empty(final_file): 
    
    print("WARNING: %s exists, suggesting that perSVade was already  run in this folder. Remove this file if you want this command to work. Exiting..."%final_file)
    sys.exit(0)

# define the name as the sample as the first 10 characters of the outdir
name_sample = fun.get_file(opt.outdir)[0:10]
print("Running perSVade into %s"%opt.outdir)

#### REPLACE THE REF GENOME ####

# define where the reference genome will be stored
reference_genome_dir = "%s/reference_genome_dir"%(opt.outdir); fun.make_folder(reference_genome_dir)
new_reference_genome_file = "%s/reference_genome.fasta"%reference_genome_dir

# rewrite the reference genome so that all the chars ar upper
all_chroms_seqRecords = [SeqRecord(Seq(str(seq.seq).upper()), id=seq.id, description="", name="") for seq in SeqIO.parse(opt.ref, "fasta")]
SeqIO.write(all_chroms_seqRecords, new_reference_genome_file, "fasta")
opt.ref = new_reference_genome_file

# check that the mitoChromosomes are in the ref
all_chroms = {s.id for s in SeqIO.parse(opt.ref, "fasta")}
if any([x not in all_chroms for x in opt.mitochondrial_chromosome.split(",")]) and opt.mitochondrial_chromosome!="no_mitochondria":
    raise ValueError("The provided mitochondrial_chromosomes are not in the reference genome.")

# get the genome len
genome_length = sum(fun.get_chr_to_len(opt.ref).values())
print("The genome has %.2f Mb"%(genome_length/1000000 ))

##################################

#### REPLACE THE GFF ####
target_gff = "%s/reference_genome_features.gff"%reference_genome_dir

# copy the gff
if opt.gff is None: print("WARNING: gff was not provided. This will be a problem if you want to annotate small variant calls")
else:

    if fun.file_is_empty(target_gff): fun.soft_link_files(opt.gff, target_gff)

    # change the path
    opt.gff = target_gff

#########################

#### REPLACE THE REPEATS TABLE IF PROVIDED ####

# define the 
repeats_table_file = "%s.repeats.tab"%opt.ref

if opt.previous_repeats_table is not None:
    print("using privided repeats %s"%opt.previous_repeats_table)

    if fun.file_is_empty(opt.previous_repeats_table): raise ValueError("The provided repeats table does not exist")
    # softlink
    fun.soft_link_files(opt.previous_repeats_table, repeats_table_file)

###############################################

if opt.StopAfter_genomeObtention is True: 
    print("Stopping pipeline after the genome obtention.")
    sys.exit(0)

#### define misc args ####

# get the simulation ploidies
if opt.simulation_ploidies!="auto": simulation_ploidies = opt.simulation_ploidies.split(",")

else: 

    # for pooled seq it takes simulated on 1 in 10
    if opt.pooled_sequencing is True: simulation_ploidies = ["diploid_homo", "ref:9_var:1"]
    elif opt.ploidy==1: simulation_ploidies = ["haploid"]
    else: simulation_ploidies = ["diploid_homo", "ref:%i_var:1"%(opt.ploidy-1)]

# set a specific calling for pooled sequencing
if opt.pooled_sequencing is True: print("WARNING: Running on pooled sequencing.  These are the simulated ploidies for the SV calling parameter optimisation:", simulation_ploidies)

# the window length for all operations
valid_chrom_lens = [len_seq for chrom, len_seq  in fun.get_chr_to_len(opt.ref).items() if chrom not in opt.mitochondrial_chromosome.split(",") and len_seq>=opt.min_chromosome_len]
if len(valid_chrom_lens)==0: raise ValueError("There are no chromosomes to calculate the window_l. Decrease --min_chromosome_len.")
print("There are %i chromosomes to calculate window length"%len(valid_chrom_lens))
fun.window_l = int(np.median(valid_chrom_lens)*0.05) + 1
if pd.isna(fun.window_l): fun.window_l = 1000
print("using a window length of %i"%fun.window_l)

# define the verbosity. If opt.verbose is False, none of the 'print' statements of sv_functions will have an effect
fun.printing_verbose_mode = opt.verbose

# defin the fraction of available mem
fun.fraction_available_mem = opt.fraction_available_mem

# define the min_CNVsize_betweenBPs
fun.min_CNVsize_betweenBPs = opt.min_CNVsize_betweenBPs

# redefine the real threads
real_available_threads = fun.get_available_threads(opt.outdir)
if opt.threads>real_available_threads:  print("WARNING: There are %i available threads, and you required %i."%(real_available_threads, opt.threads))

print("Running with %i Gb of RAM and %i cores"%(int(fun.get_availableGbRAM(opt.outdir)), opt.threads))

# change the default parameters if specified
if opt.parameters_json_file is not None:

    gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup = fun.get_parameters_from_json(opt.parameters_json_file)

    fun.default_filtersDict_gridss = gridss_filters_dict
    fun.default_gridss_blacklisted_regions = gridss_blacklisted_regions
    fun.default_gridss_maxcoverage = gridss_maxcoverage
    fun.default_max_rel_coverage_to_consider_del = max_rel_coverage_to_consider_del
    fun.default_min_rel_coverage_to_consider_dup = min_rel_coverage_to_consider_dup


# get the gff info
if opt.gff is not None: correct_gff, gff_with_biotype = fun.get_correct_gff_and_gff_with_biotype(opt.gff, replace=opt.replace)

# get the repeats table
if opt.skip_repeat_analysis is False:

    print("getting repeats")
    repeats_df, repeats_table_file = fun.get_repeat_maskerDF(opt.ref, threads=opt.threads, replace=opt.replace)

else:

    print("skipping the repeats analysis")
    fun.write_repeats_table_file(repeats_table_file)

if opt.StopAfter_repeatsObtention is True:
    print("Stopping after the obtention of repeats")
    sys.exit(0)

#############################

########################################
########################################
########################################


#####################################
############# BAM FILE ##############
#####################################

if not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]):

    ##### DEFINE THE SORTED BAM #####

    # define files that may be used in many steps of the pipeline
    if opt.sortedbam is None:

        bamfile = "%s/aligned_reads.bam"%opt.outdir
        sorted_bam = "%s.sorted"%bamfile
        index_bam = "%s.bai"%sorted_bam

    else:

        # debug the fact that you prvided reads and bam. You should just provide one
        if any([not x is None for x in {opt.fastq1, opt.fastq2}]): raise ValueError("You have provided reads and a bam, you should only provide one")

        # get the files
        sorted_bam = opt.sortedbam
        index_bam = "%s.bai"%sorted_bam

    ###################################

    # normal alignment of provided reads
    if all([not x is None for x in {opt.fastq1, opt.fastq2}]):

        # if the reads have to be QC and trimmed:
        if opt.QC_and_trimming_reads is True: 
            print("running trimming and QC of the reads")
            opt.fastq1, opt.fastq2 = fun.run_trimmomatic(opt.fastq1, opt.fastq2, replace=opt.replace, threads=opt.threads)

        print("WORKING ON ALIGNMENT")
        fun.run_bwa_mem(opt.fastq1, opt.fastq2, opt.ref, opt.outdir, bamfile, sorted_bam, index_bam, name_sample, threads=opt.threads, replace=opt.replace)

        fun.clean_sorted_bam_coverage_per_window_files(sorted_bam)

    else: print("Warning: No fastq file given, assuming that you provided a bam file")

    # check that all the important files exist
    if any([fun.file_is_empty(x) for x in {sorted_bam, index_bam}]): raise ValueError("You need the sorted and indexed bam files in ")

#####################################
#####################################
#####################################

###########################################
############# NECESSARY FILES #############
###########################################

# First create some files that are important for any program

# Create a reference dictionary
fun.create_sequence_dict(opt.ref, replace=opt.replace)

# Index the reference
fun.index_genome(opt.ref, replace=opt.replace)

#### calculate coverage per windows of window_l ####

if not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]):

    destination_dir = "%s.calculating_windowcoverage"%sorted_bam
    coverage_file = fun.generate_coverage_per_window_file_parallel(opt.ref, destination_dir, sorted_bam, windows_file="none", replace=opt.replace, run_in_parallel=True, delete_bams=True)

####################################################



###########################################
###########################################
###########################################

if opt.StopAfter_bamFileObtention is True: 
    print("Stopping pipeline after the bamfile obtention.")
    sys.exit(0)

#####################################
##### STRUCTURAL VARIATION ##########
#####################################


##### find a dict that maps each svtype to a file with a set of real SVs (real_svtype_to_file) #####
all_svs = {'translocations', 'insertions', 'deletions', 'inversions', 'tandemDuplications'}

if opt.real_bedpe_breakpoints is not None and opt.fast_SVcalling is False: 
    print("using the set of real variants from %s"%opt.real_bedpe_breakpoints)
    real_bedpe_breakpoints = opt.real_bedpe_breakpoints

elif opt.fast_SVcalling is False and opt.close_shortReads_table is not None:
    
    # the table was provided
    if opt.close_shortReads_table!="auto": 

        print("finding the set of compatible SVs from %s"%opt.close_shortReads_table)

        # define the outdir for the real vars
        outdir_finding_realVars = "%s/findingRealSVs_providedCloseReads"%opt.outdir

    # a taxID was provided, which overrides the value of opt.genomes_withSV_and_shortReads_table
    else:

        print("finding close genomes or reads for close taxIDs in the SRA database for taxID %s"%opt.target_taxID)

        # define the outdir for the real vars
        outdir_finding_realVars = "%s/findingRealSVs_automaticFindingOfCloseReads"%opt.outdir; fun.make_folder(outdir_finding_realVars)

        # define the outdir where the close genomes whould be downloaded
        outdir_getting_closeReads = "%s/getting_closeReads"%outdir_finding_realVars; fun.make_folder(outdir_getting_closeReads)

        opt.close_shortReads_table = fun.get_close_shortReads_table_close_to_taxID(opt.target_taxID, opt.ref, outdir_getting_closeReads, opt.ploidy, n_close_samples=opt.n_close_samples, nruns_per_sample=opt.nruns_per_sample, replace=opt.replace, threads=opt.threads, job_array_mode=opt.job_array_mode, StopAfter_sampleIndexingFromSRA=opt.StopAfter_sampleIndexingFromSRA, StopAfterPrefecth_of_reads=opt.StopAfterPrefecth_of_reads, max_coverage_sra_reads=opt.max_coverage_sra_reads)

    # skip the running of the pipeline 
    if opt.StopAfter_readObtentionFromSRA:
        print("Stopping pipeline after the reads obtention from SRA")
        sys.exit(0) 

    # skip pipeline running if you have to stop after prefetch of reads
    if opt.StopAfterPrefecth_of_reads:
        print("Stopping pipeline after the prefetch of reads")
        sys.exit(0) 

    # get the real SVs
    real_bedpe_breakpoints = fun.get_compatible_real_bedpe_breakpoints(opt.close_shortReads_table, opt.ref, outdir_finding_realVars, replace=opt.replace, threads=opt.threads, max_nvars=opt.nvars, mitochondrial_chromosome=opt.mitochondrial_chromosome, job_array_mode=opt.job_array_mode, parameters_json_file=opt.parameters_json_file)

else: 
    print("Avoiding the simulation of real variants. Only inserting randomSV.")

    # define the set of vars as empty. This will trigger the random generation of vars
    real_bedpe_breakpoints = None


if opt.StopAfter_obtentionOFcloseSVs: 
    print("stopping pipeline after obtention of close SVs")
    sys.exit(0)



###################################################################################################

# test accuracy on real data
if opt.testAccuracy is True:  

    print("testing accuracy on simulations and real variants (if provided)")

    # test that you have provided a opt.close_shortReads_table
    if opt.close_shortReads_table is None or opt.fast_SVcalling is True: 
        raise ValueError("You have to specify a --close_shortReads_table and not run in --fast_SVcalling to test the accuracy of the pipeline on several datasets (--testAccuracy)")

    fun.report_accuracy_realSVs(opt.close_shortReads_table, opt.ref, "%s/testing_Accuracy"%opt.outdir, real_bedpe_breakpoints, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, job_array_mode=opt.job_array_mode, StopAfter_testAccuracy_perSVadeRunning=opt.StopAfter_testAccuracy_perSVadeRunning, skip_cleaning_simulations_files_and_parameters=opt.skip_cleaning_simulations_files_and_parameters, skip_cleaning_outdir=opt.skip_cleaning_outdir, parameters_json_file=opt.parameters_json_file, gff=opt.gff, replace_FromGridssRun_final_perSVade_run=opt.replace_FromGridssRun_final_perSVade_run, fraction_available_mem=opt.fraction_available_mem, replace_SV_CNVcalling_and_optimisation=opt.replace_SV_CNVcalling_and_optimisation, replace_only_SV_CNVcalling=opt.replace_only_SV_CNVcalling)


# get the golden set
if opt.goldenSet_dir is not None:

    needstoberefactoredwithnewstructure

    outdir_goldenSet = "%s/testing_goldenSetAccuracy"%opt.outdir
    fun.report_accuracy_golden_set(opt.goldenSet_dir, outdir_goldenSet, opt.ref, real_svtype_to_file, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, job_array_mode=opt.job_array_mode, StopAfter_sampleIndexingFromSRA=opt.StopAfter_sampleIndexingFromSRA, StopAfterPrefecth_of_reads=opt.StopAfterPrefecth_of_reads_goldenSet, target_taxID=opt.target_taxID, parameters_json_file=opt.parameters_json_file, fraction_available_mem=opt.fraction_available_mem)

# run the actual perSVade function optimising parameters
if opt.skip_SVcalling is False and not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]):

    SVdetection_outdir = "%s/SVdetection_output"%opt.outdir
    outdir_gridss_final = fun.run_GridssClove_optimising_parameters(sorted_bam, opt.ref, SVdetection_outdir, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, fast_SVcalling=opt.fast_SVcalling, real_bedpe_breakpoints=real_bedpe_breakpoints, replace_FromGridssRun_final_perSVade_run=opt.replace_FromGridssRun_final_perSVade_run)


print("structural variation analysis with perSVade finished")

#####################################
#####################################
#####################################

#####################################
###### SV and CNV ANNOTATION ########
#####################################

if opt.skip_SVcalling is False and not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]) and opt.skip_SV_CNV_calling is False:

    # define the df_bedpe and df_gridss
    df_gridss = fun.get_svtype_to_svfile_and_df_gridss_from_perSVade_outdir(opt.outdir, opt.ref)[1]

    # run CNVcalling 
    cnv_calling_outdir = "%s/CNV_calling"%opt.outdir
    fun.run_CNV_calling(sorted_bam, opt.ref, cnv_calling_outdir, opt.threads, opt.replace, opt.mitochondrial_chromosome, df_gridss, opt.window_size_CNVcalling, opt.ploidy)












    finished_CNV_calling

    # delete key files if replace_SV_CNVcalling
    if opt.replace_SV_CNVcalling_and_optimisation is True: fun.remove_files_SV_CNVcalling(opt.outdir)

    #### get CNV parameters ####

    # get the default parameters
    if opt.min_coverage_duplication=="auto" or opt.max_coverage_deletion=="auto" or opt.min_r_pearson_noFlatRegions=="auto" or opt.min_r_spearman_noFlatRegions=="auto":

        auto_min_coverage_duplication, auto_max_coverage_deletion, auto_min_r_pearson_noFlatRegions, auto_min_r_spearman_noFlatRegions = fun.get_automatic_coverage_thresholds(opt.outdir, opt.fast_SVcalling, opt.ploidy, opt.nsimulations, simulation_ploidies, opt.mitochondrial_chromosome, opt.ref, threads=opt.threads, replace=opt.replace)


    # set real vars
    if opt.min_coverage_duplication=="auto": min_coverage_duplication = auto_min_coverage_duplication
    else: min_coverage_duplication = float(opt.min_coverage_duplication)
        
    if opt.max_coverage_deletion=="auto": max_coverage_deletion = auto_max_coverage_deletion
    else: max_coverage_deletion = float(opt.max_coverage_deletion)

    if opt.min_r_pearson_noFlatRegions=="auto": min_r_pearson_noFlatRegions = auto_min_r_pearson_noFlatRegions
    else: min_r_pearson_noFlatRegions = float(opt.min_r_pearson_noFlatRegions)
        
    if opt.min_r_spearman_noFlatRegions=="auto": min_r_spearman_noFlatRegions = auto_min_r_spearman_noFlatRegions
    else: min_r_spearman_noFlatRegions = float(opt.min_r_spearman_noFlatRegions)

    ############################

    # get the variant calling 
    outdir_var_calling = "%s/SVcalling_output"%opt.outdir
    print("getting all SVs into one VCF. Regions with a cov>%.3f will be treated as DUP, and regions with cov <%.3f will be treated as DEL. Regions with a correlation between position and coverage >%.3f (pearson) and >%.3f (spearman) will not be considered as CNV, as they may be related to the smiley-face effect."%(min_coverage_duplication, max_coverage_deletion, min_r_pearson_noFlatRegions, min_r_spearman_noFlatRegions))

    # delete only the CNV calling, but not the optimisation
    if opt.replace_only_SV_CNVcalling is True: fun.delete_folder(outdir_var_calling)

    SV_CNV_vcf = fun.get_vcf_all_SVs_and_CNV(opt.outdir, outdir_var_calling, sorted_bam, opt.ref, opt.ploidy, replace=opt.replace, threads=opt.threads, mitochondrial_chromosome=opt.mitochondrial_chromosome, mito_code=opt.mitochondrial_code, gDNA_code=opt.gDNA_code, max_coverage_deletion=max_coverage_deletion, min_coverage_duplication=min_coverage_duplication, min_r_pearson_noFlatRegions=min_r_pearson_noFlatRegions, min_r_spearman_noFlatRegions=min_r_spearman_noFlatRegions)

    print("the SV and CNV calling vcf can be found in %s"%SV_CNV_vcf)

    # get variant annotation
    if opt.gff is not None:

        print("annotating SV, CNV variants with VEP")
        SV_CNV_vcf_annotated = fun.annotate_SVs_inHouse(SV_CNV_vcf, gff_with_biotype, opt.ref, replace=opt.replace, threads=opt.threads, mitochondrial_chromosome=opt.mitochondrial_chromosome, mito_code=opt.mitochondrial_code, gDNA_code=opt.gDNA_code)

        print("annotated SV vcf can be found in %s"%SV_CNV_vcf_annotated)
    
    else: print("WARNING: Skipping SV annotation because -gff was not provided.")

#####################################
#####################################
#####################################

# stop after the generation of SV and CNV calls


#####################################
###### SMALL VARS AND CNV ###########
#####################################

if opt.run_smallVarsCNV:

    # define an outdir
    outdir_varcall = "%s/smallVars_CNV_output"%opt.outdir

    # define the basic cmd
    varcall_cmd = "%s -r %s --threads %i --outdir %s -p %i -sbam %s --caller %s --coverage %i --mitochondrial_chromosome %s --mitochondrial_code %i --gDNA_code %i --minAF_smallVars %s"%(varcall_cnv_pipeline, opt.ref, opt.threads, outdir_varcall, opt.ploidy, sorted_bam, opt.caller, opt.coverage, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code, opt.minAF_smallVars)

    # add options
    if opt.replace is True: varcall_cmd += " --replace"
    if opt.gff is not None: varcall_cmd += " -gff %s"%opt.gff
    if opt.remove_smallVarsCNV_nonEssentialFiles is True: varcall_cmd += " --remove_smallVarsCNV_nonEssentialFiles"
    if opt.StopAfter_smallVarCallSimpleRunning is True: varcall_cmd += " --StopAfter_smallVarCallSimpleRunning"
    if opt.replace_var_integration is True: varcall_cmd += " --replace_var_integration"
    if opt.pooled_sequencing is True: varcall_cmd += " --pooled_sequencing"
    if opt.consider_repeats_smallVarCall is True: varcall_cmd += " --repeats_table %s"%repeats_table_file
    if opt.generate_alternative_genome is True: varcall_cmd += " --generate_alternative_genome"
    if opt.skip_cnv_analysis is True: varcall_cmd += " --skip_cnv_analysis"

    # run
    if __name__ == '__main__': fun.run_cmd(varcall_cmd)

    # define the small variants vcf
    small_vars_vcf = "%s/variants_atLeast1PASS_ploidy%i.vcf"%(outdir_varcall, opt.ploidy)

    # define the variant annotation
    small_vars_var_annotation = "%s/variant_annotation_ploidy%i.tab"%(outdir_varcall, opt.ploidy)


#####################################
#####################################
#####################################


#####################################
##### VARIANTS VISUALIZATION ########
#####################################

if (opt.skip_SVcalling is False or run_smallVarsCNV is True) and not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]) and opt.visualization_results is True and opt.gff is not None:

    # visualize the results with the browser

    # create a table with the data to visualize the browser for this sample. This will include many things
    dict_data = {"sorted_bam":sorted_bam, "sampleID":name_sample}
    if opt.skip_SVcalling is False: 
        dict_data["SV_CNV_vcf"] = SV_CNV_vcf
        dict_data["SV_CNV_var_annotation"] = SV_CNV_vcf_annotated

    if opt.run_smallVarsCNV is True: 
        dict_data["smallVars_vcf"] = small_vars_vcf
        dict_data["smallVars_var_annotation"] = small_vars_var_annotation

    # make df
    df_visualization = pd.DataFrame({0:dict_data}).transpose()

    # save to file
    outdir_visualization = "%s/variant_visualization"%opt.outdir; fun.make_folder(outdir_visualization)
    visualization_file = "%s/visualization_data.tab"%outdir_visualization
    df_visualization.to_csv(visualization_file, sep="\t", index=False, header=True)

    # get the visualization cmd
    cmd_visualization = "%s --input_data %s --outdir %s --reference_genome %s --gff %s --threads %i --only_affected_genes --mitochondrial_chromosome %s"%(perSVade_genome_browser, visualization_file, outdir_visualization, opt.ref, opt.gff, opt.threads, opt.mitochondrial_chromosome)
    if opt.replace is True: cmd_visualization += " --replace"

    fun.run_cmd(cmd_visualization)

#####################################
#####################################
#####################################

# at the end you want to clean the outdir to keep only the essential files
if opt.skip_cleaning_outdir is False: fun.clean_perSVade_outdir(opt.outdir)

# generate a file that indicates whether the gridss run is finished
open(final_file, "w").write("perSVade_finished finished...")

print("perSVade Finished correctly")


