#!/usr/bin/env python

# This is the perSVade pipeline main script, which shoul dbe run on the perSVade conda environment


##### DEFINE ENVIRONMENT #######

# general module imports
import argparse, os
from argparse import RawTextHelpFormatter
import copy as cp
import pickle
import string
import shutil 
import random
import sys
from shutil import copyfile
import time

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import sv_functions as fun

# import persvade-specific modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np

# packages installed into the conda environment 
samtools = "%s/bin/samtools"%EnvDir
java = "%s/bin/java"%EnvDir

# scripts that are installed under this software
varcall_cnv_pipeline = "%s/varcall_cnv_pipeline.py"%CWD
perSVade_genome_browser = "%s/perSVade_genome_browser.py"%CWD

#######

description = """
Runs perSVade pipeline on an input set of paired end short ends. See https://github.com/Gabaldonlab/perSVade for more information.
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

#### general, mandatory args ####

mandatory_args = parser.add_argument_group("MANDATORY ARGUMENTS")

mandatory_args.add_argument("-r", "--ref", dest="ref", required=True, help="Reference genome. Has to end with .fasta.")
mandatory_args.add_argument("-o", "--outdir", dest="outdir", action="store", required=True, help="Directory where the data will be stored")
mandatory_args.add_argument("-mchr", "--mitochondrial_chromosome", dest="mitochondrial_chromosome", required=True, type=str, help="The name of the mitochondrial chromosome. This is important if you have mitochondrial proteins for which to annotate the impact of nonsynonymous variants, as the mitochondrial genetic code is different. This should be the same as in the gff. If there is no mitochondria just put 'no_mitochondria'. If there is more than one mitochindrial scaffold, provide them as comma-sepparated IDs.")

#################################

#### modularity of the pipeline ####
module_inputs = parser.add_argument_group("PIPELINE MODULES")
module_inputs.add_argument("--type_variant_calling", dest="type_variant_calling", default=None, help="Specify which modules to execute. By default, perSVade runs only coverage-based CNV and SV calling (with random-based parameter optimisation). This command can be any combination of 'SV', 'coverageCNV' and 'small_vars' in a comma comma-sepparated manner. For example, if you want to execute SV and small variant calling you can specify --type_variant_calling SV,small_vars. Note that this command overrides the effect of --run_smallVarsCNV, --remove_smallVarsCNV_nonEssentialFiles, --skip_SVcalling, --skip_CNV_calling. If you want to execute this argument you need to specify some sequence input (either reads or bam file).")

####################################

###### inputs #######

seq_inputs = parser.add_argument_group("SEQUENCE INPUTS")

seq_inputs.add_argument("-f1", "--fastq1", dest="fastq1", default=None, help="fastq_1 file. Option required to obtain bam files. It can be 'skip', which will tell the pipeline to not use any fastq or bam input.")
seq_inputs.add_argument("-f2", "--fastq2", dest="fastq2", default=None, help="fastq_2 file. Option required to obtain bam files. It can be 'skip', which will tell the pipeline to not use any fastq or bam input.")

seq_inputs.add_argument("--input_SRRfile", dest="input_SRRfile", default=None, help="An input srr file that can be provided instead of the fastq files. If this is provided the pipeline will run fastqdump on the reads, and also run fastqc and trimmomatic on them. This may be dangerous, since you may supervise the downloading and quality control of the reads.")

seq_inputs.add_argument("-sbam", "--sortedbam", dest="sortedbam", default=None, help="The path to the sorted bam file, which should have a bam.bai file in the same dir. For example, if your bam file is called 'aligned_reads.bam', there should be an 'aligned_reads.bam.bai' as well. This is mutually exclusive with providing reads. By default, it is assumed that this bam has marked duplicates.")

seq_inputs.add_argument("--downsampled_coverage", dest="downsampled_coverage", default=None, type=float, help="A float indicating whether to downsample to a specific coverage for faster running. This can be useful for testing some argument. For example, '--downsampled_coverage 5.0' would downsample randomly your input reads or bam to an average coverage of 5x.")

seq_inputs.add_argument("--other_perSVade_outdirs_sameReadsANDalignment", dest="other_perSVade_outdirs_sameReadsANDalignment",  default=None, help="A comma-sepparated set of full paths to perSVade outdirs that can be used to replace the sorted bam and reads dir.")


#####################

########### resource allocation ##############

resources_args = parser.add_argument_group("RESOURCES")

resources_args.add_argument("-thr", "--threads", dest="threads", default=16, type=int, help="Number of threads, Default: 16")

resources_args.add_argument("--fraction_available_mem", dest="fraction_available_mem", default=None, type=float, help="The fraction of RAM that is being allocated to this perSVade run. In several steps, this pipeline needs to calculate the available memory (using psutil.virtual_memory()). This returns all the available memory in the computer. If you are running on a fraction of the computers' resources, this calculation is overestimating the available RAM. In such case you can provide the fraction available through this argument. By default, it will calculate the available ram by filling the memory, which may give errors. It is highly reccommended that you provide this option. If you want to use all the allocated memory you should specify --fraction_available_mem 1.0")

resources_args.add_argument("--fractionRAM_to_dedicate", dest="fractionRAM_to_dedicate", type=float,  default=0.5, help="This is the fraction of the available memory that will be used by several java programs that require a heap size. By default we set this to 0.5 to not overload the system")

resources_args.add_argument("--tmpdir", dest="tmpdir", default=None, help="A full path to a directory where to write intermediate files. This is useful if you are running on a cluster that has some directories that have higher writing speed than others.")

resources_args.add_argument("--min_gb_RAM_required", dest="min_gb_RAM_required", default=2, type=int, help="The minimum number of RAM required to run perSVade. This can be related to --fraction_available_mem")



##############################################

#### general, optional args #### 

general_optional_args = parser.add_argument_group("GENERAL OPTIONAL ARGUMENTS")

general_optional_args.add_argument("-p", "--ploidy", dest="ploidy", default=1, type=int, help="Ploidy, can be 1 or 2")

general_optional_args.add_argument("-gff", "--gff-file", dest="gff", default=None, help="path to the GFF3 annotation of the reference genome. Make sure that the IDs are completely unique for each 'gene' tag. This is necessary for both the CNV analysis (it will look at genes there) and the annotation of the variants.")

general_optional_args.add_argument("-mcode", "--mitochondrial_code", dest="mitochondrial_code", default=3, type=int, help="The code of the NCBI mitochondrial genetic code. For yeasts it is 3. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. The information of this website may be wrong, so you may want to double check with the literature.")

general_optional_args.add_argument("-gcode", "--gDNA_code", dest="gDNA_code", default=1, type=int, help="The code of the NCBI gDNA genetic code. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . For C. albicans it is 12. The information of this website may be wrong, so you may want to double check with the literature.")


general_optional_args.add_argument("--QC_and_trimming_reads", dest="QC_and_trimming_reads", action="store_true", default=False, help="Will run fastq and trimmomatic of reads, and use the trimmed reads for downstream analysis. This option will generate files under the same dir as f1 and f2, so be aware of it. This is a rather dangerous option, sinc you may want to do the quality control of the reads before running perSVade.")

general_optional_args.add_argument("--min_chromosome_len", dest="min_chromosome_len", default=100000, type=int, help="The minimum length to consider chromosomes from the provided fasta for calculating the window length. Any chromosomes that shorter than the window length will not be considered in the random SV simulations.")

general_optional_args.add_argument("--verbose", dest="verbose", action="store_true", default=False, help="Whether to print a verbose output. This is important if there are any errors in the run.")

general_optional_args.add_argument("--previous_repeats_table", dest="previous_repeats_table", default=None, help="This may be the path to a file that contains the processed output of RepeatMasker (such as the one output by the function get_repeat_maskerDF). This should be a table with the following header: 'SW_score, perc_div, perc_del, perc_ins, chromosome, begin_repeat, end_repeat, left_repeat, strand, repeat, type, position_inRepeat_begin, position_inRepeat_end, left_positionINrepeat, IDrepeat'. It is created by parsing the tabular output of RepeatMasker and putting into a real .tab format.")

################################

######### SV calling and optimisation vars ############

SVcalling_args = parser.add_argument_group("SV CALLING")

SVcalling_args.add_argument("--parameters_json_file", dest="parameters_json_file", type=str, default=None, help="A file with the json parameters to use. This only has effect if --fast_SVcalling is specified")

SVcalling_args.add_argument("--fast_SVcalling", dest="fast_SVcalling", action="store_true", default=False, help="Run SV calling with a default set of parameters. There will not be any optimisation nor reporting of accuracy. This is expected to work almost as fast as gridss and clove together. If --parameters_json_file, the parameters are substituted by the json parameters.")

SVcalling_args.add_argument("--simulation_chromosomes", dest="simulation_chromosomes", type=str, default=None, help="A comma-sepparated set of chromosomes (i.e.: chr1,chr2,chr3) in which to perform simulations. By default it takes all the chromosomes.")

SVcalling_args.add_argument("--nvars", dest="nvars", default=50, type=int, help="Number of variants to simulate for each SVtype.")

SVcalling_args.add_argument("--nsimulations", dest="nsimulations", default=2, type=int, help="The number of 'replicate' simulations that will be produced.")

SVcalling_args.add_argument("--simulation_ploidies", dest="simulation_ploidies", type=str, default="auto", help='A comma-sepparated string of the ploidies to simulate for parameter optimisation. It can have any of "haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1" . By default it will be inferred from the ploidy. For example, if you are running on ploidy=2, it will optimise for diploid_hetero.')

SVcalling_args.add_argument("--range_filtering_benchmark", dest="range_filtering_benchmark", type=str, default="theoretically_meaningful", help='The range of parameters that should be tested in the SV optimisation pipeline. It can be any of large, medium, small, theoretically_meaningful, theoretically_meaningful_NoFilterRepeats or single. ')

SVcalling_args.add_argument("--simulate_SVs_around_repeats", dest="simulate_SVs_around_repeats", action="store_true", default=False, help="Simulate SVs around repeats. This requires that there are some repeats inferred. This option will generate a simulated set of breakpoints around repeats, if possible of the same family, and with random orientations.")

SVcalling_args.add_argument("--simulate_SVs_around_HomologousRegions", dest="simulate_SVs_around_HomologousRegions", action="store_true", default=False, help="Simulate SVs around regions that have high similarity.")

SVcalling_args.add_argument("--simulate_SVs_around_HomologousRegions_queryWindowSize", dest="simulate_SVs_around_HomologousRegions_queryWindowSize",  default=500, type=int, help="The window size used for finding regions with high similarity. This only works if --simulate_SVs_around_HomologousRegions is specified.")

SVcalling_args.add_argument("--simulate_SVs_around_HomologousRegions_maxEvalue", dest="simulate_SVs_around_HomologousRegions_maxEvalue",  default=1e-5, type=float, help="The maximum evalue by which two regions will be said to have high similarity. This only works if --simulate_SVs_around_HomologousRegions is specified.")

SVcalling_args.add_argument("--simulate_SVs_around_HomologousRegions_minPctOverlap", dest="simulate_SVs_around_HomologousRegions_minPctOverlap",  default=50, type=int, help="The minimum percentage of overlap between two homologous regions so that there can be a brèkpoint drawn in between. This only works if --simulate_SVs_around_HomologousRegions is specified.")

SVcalling_args.add_argument("--simulate_SVs_around_HomologousRegions_previousBlastnFile", dest="simulate_SVs_around_HomologousRegions_previousBlastnFile",  default=None, help="A file that contains the output of blastn of regions of the genome against itself. It will be put under refetence genome dir")

########################################################

#### arguments related to the definition of "real SVs". All are optional #### 

realSVs_args = parser.add_argument_group("SV CALLING. FINDING REAL SVs")

realSVs_args.add_argument("--close_shortReads_table", dest="close_shortReads_table", type=str, default=None, help="This is the path to a table that has 4 fields: sampleID,runID,short_reads1,short_reads2. These should be WGS runs of samples that are close to the reference genome and some expected SV. Whenever this argument is provided, the pipeline will find SVs in these samples and generate a folder <outdir>/findingRealSVs<>/SVs_compatible_to_insert that will contain one file for each SV, so that they are compatible and ready to insert in a simulated genome. This table will be used if --testAccuracy is specified, which will require at least 3 runs for each sample. It can be 'auto', in which case it will be inferred from the taxID provided by --target_taxID.")

realSVs_args.add_argument("--real_bedpe_breakpoints", dest="real_bedpe_breakpoints", type=str, default=None, help="A file with the list of 'real' breakpoints around which to insert the SVs in simulations. It may be created with --close_shortReads_table. If both --real_bedpe_breakpoints and --close_shortReads_table are provided, --real_bedpe_breakpoints will be used, and --close_shortReads_table will have no effect. If none of them are provided, this pipeline will base the parameter optimization on randomly inserted SVs (the default behavior). The coordinates have to be 1-based, as they are ready to insert into RSVsim.")

realSVs_args.add_argument("--job_array_mode", dest="job_array_mode", type=str, default="local", help="It specifies in how to run the job arrays for,  --testAccuracy, the downloading of reads if  --close_shortReads_table is auto, and the SV calling for the table in --close_shortReads_table. It can be 'local' (runs one job after the other or 'job_array'. If 'job_array' is specified, this pipeline will generate a file with all the jobs to run, and stop. You will have to run these jobs before stepping into downstream analyses.")

##############################################################################

##### read depth-based CNV callings #####

CNV_args = parser.add_argument_group("READ DEPTH-BASED CNV CALLING")

CNV_args.add_argument("--min_CNVsize_coverageBased", dest="min_CNVsize_coverageBased", default=300, type=int, help="The minimum size of a CNV inferred from coverage.")

CNV_args.add_argument("--window_size_CNVcalling", dest="window_size_CNVcalling", default=100, type=int, help="The window size in which the genome will be fragmented for CNV calling.")

CNV_args.add_argument("--cnv_calling_algs", dest="cnv_calling_algs", default="HMMcopy,AneuFinder", type=str, help="A comma-sepparated string thatindicates which programs should be used for the CNV calling. It can be any of HMMcopy,CONY,AneuFinder. We note that CONY does not work well for small chromosomes or large binned windows.")

CNV_args.add_argument("--bg_sorted_bam_CNV", dest="bg_sorted_bam_CNV", default=None, help="This is a sorted bam (with duplicated marked) that is taken as a 'reference' background in the CNV calling. By default, perSVade corrects the coverage by GC content, mappability and the distance to the telomere. If --bg_sorted_bam_CNV, the coverage will be normalised by the coverage of this sorted bam.")

#########################################

######### small variant calling #######

smallVars_args = parser.add_argument_group("SMALL VARIANT CALLING AND COVERAGE PER GENE CALCULATION")

smallVars_args.add_argument("--run_smallVarsCNV", dest="run_smallVarsCNV", action="store_true", default=False, help="Will call small variants and CNV.")

smallVars_args.add_argument("-caller", "--caller", dest="caller", required=False, default="all", help="SNP caller option to obtain vcf file. options: no/all/HaplotypeCaller/bcftools/freebayes. It can be a comma-sepparated string, like 'HaplotypeCaller,freebayes'")

smallVars_args.add_argument("-c", "--coverage", dest="coverage", default=20, type=int, help="minimum Coverage (int) of a position to be considered for small variant calling. This parameter should be related to the coverage of your library. You may be careful with setting it too low (i.e. <15) as it will yield many false positive calls. It is reasonable to check how other similar studies set this parameter.")

smallVars_args.add_argument("--minAF_smallVars", dest="minAF_smallVars", default="infer", help="The minimum fraction of reads covering a variant to be called. The default is 'infer', which will set a threshold based on the ploidy. This is only relevant for the final vcfs, where only PASS vars are considered. It can be a number between 0 and 1.")

smallVars_args.add_argument("--window_freebayes_bp", dest="window_freebayes_bp", default=10000, type=int, help="freebayes is run in parallel by chunks of the genome. This cmd specifies the window (in bp) in which freebayes regions are split to. If you increase this number the splitting will be in larger chunks of the genome.")

smallVars_args.add_argument("--pooled_sequencing", dest="pooled_sequencing", action="store_true", default=False, help="It is a pooled sequencing run, which means that the small variant calling is not done based on ploidy. If you are also running SV calling, you will run parameters optimisation on a sample that has 1-10 pooling.")

smallVars_args.add_argument("--run_ploidy2_ifHaploid", dest="run_ploidy2_ifHaploid", action="store_true", default=False, help="If ploidy==1, run also in diploid configuration. This is useful because there may be diploid mutations in duplicated regions.")

smallVars_args.add_argument("--consider_repeats_smallVarCall", dest="consider_repeats_smallVarCall", action="store_true", default=False, help="If --run_smallVarsCNV, this option will imply that each small  variant will have an annotation of whether it overlaps a repeat region.")

smallVars_args.add_argument("--remove_smallVarsCNV_nonEssentialFiles", dest="remove_smallVarsCNV_nonEssentialFiles", action="store_true", default=False, help="Will remove all the varCall files except the integrated final file and the bam file.")

#######################################

##### skipping some steps #####

skipping_args = parser.add_argument_group("SKIP SOME STEPS")

skipping_args.add_argument("--skip_repeat_analysis", dest="skip_repeat_analysis", default=False, action="store_true", help="Skip the inference of repeats. If --previous_repeats_table is provided, this argument will override it.")

skipping_args.add_argument("--skip_SVcalling", dest="skip_SVcalling", action="store_true", default=False, help="Do not run SV calling.")

skipping_args.add_argument("--skip_SV_CNV_calling", dest="skip_SV_CNV_calling", action="store_true", default=False, help="Do not run the integration of SV and CNV calling into a final vcf. Note that this integration will only take place if SV calling is performed.")

skipping_args.add_argument("--skip_CNV_calling", dest="skip_CNV_calling", action="store_true", default=False, help="Do not run the calling of CNV based on coverage. This is reccommended if you have samples with smiley face in coverage in a fragmented genome, which does not allow for a proper normalisation of coverage.")

skipping_args.add_argument("--skip_cnv_per_gene_analysis", dest="skip_cnv_per_gene_analysis", default=False, action="store_true", help="Don't perform the cnv analysis where we calculate coverage per gene and per equal windows of the genome. This refers to the per-gene CNV analysis, not the per-window.")

skipping_args.add_argument("--skip_marking_duplicates", dest="skip_marking_duplicates", default=False, action="store_true", help="Don't mark the duplicate reads in the .bam file.")


###############################

##### stopping options ######

stopping_args = parser.add_argument_group("STOP AFTER SOME OPERATIONS")

stopping_args.add_argument("--StopAfter_bamFileObtention", dest="StopAfter_bamFileObtention", action="store_true", default=False, help="Stop after obtaining the BAM file of aligned reads.")

stopping_args.add_argument("--StopAfter_obtentionOFcloseSVs", dest="StopAfter_obtentionOFcloseSVs", action="store_true", default=False, help="Stop after obtaining the real_bedpe_breakpoints ")

stopping_args.add_argument("--StopAfter_repeatsObtention", dest="StopAfter_repeatsObtention", action="store_true", default=False, help="Stop after obtaining  the repeats table")

stopping_args.add_argument("--StopAfter_smallVarCallSimpleRunning", dest="StopAfter_smallVarCallSimpleRunning", action="store_true", default=False, help="Stop after obtaining the filtered vcf outputs of each program.")

stopping_args.add_argument("--StopAfter_smallVarCall", dest="StopAfter_smallVarCall", action="store_true", default=False, help="Stop after running for small variant calls and calculation of per gene coverage.")

#############################

#### replacing options ####

replacing_args = parser.add_argument_group("REPLACE SOME STEPS")

replacing_args.add_argument("--replace", dest="replace", action="store_true", help="Replace existing files")

replacing_args.add_argument("--replace_SV_CNVcalling", dest="replace_SV_CNVcalling", action="store_true", help="Replace everything related to the SV and CNV calling.")

replacing_args.add_argument("--replace_FromGridssRun_final_perSVade_run", dest="replace_FromGridssRun_final_perSVade_run", action="store_true", help="Replace from the clove running in the final gridss+clove running")

replacing_args.add_argument("--replace_var_integration", dest="replace_var_integration", action="store_true", help="Replace all the variant integration steps for smallVariantCalling.")

replacing_args.add_argument("--replace_var_annotation", dest="replace_var_annotation", action="store_true", help="Replace all the variant annotation steps for smallVariantCalling, SV and CNV calling.")

###########################

#### debugging options. For testing only ####

debug_args = parser.add_argument_group("ONLY TO BE USED FOR DEVELOPMENT. USE AT YOUR OWN RISK!")

debug_args.add_argument("--target_taxID", dest="target_taxID", type=int, default=None, help="This is the taxID (according to NCBI taxonomy) to which your reference genome belongs. If provided it is used to download genomes and reads.")

debug_args.add_argument("--goldenSet_max_n_samples", dest="goldenSet_max_n_samples", default=6, type=int, help="Number of runs to run golden set analysis on.")

debug_args.add_argument("--goldenSet_table", dest="goldenSet_table", type=str, default=None, help="This is the path to a table that has ONT reads for several samples, each line one sample. It can be 'auto', which will download as many datasets as possible not exceeding goldenSet_max_n_samples")

debug_args.add_argument("--n_close_samples", dest="n_close_samples", default=5, type=int, help="Number of close samples to search in case --target_taxID is provided")

debug_args.add_argument("--nruns_per_sample", dest="nruns_per_sample", default=3, type=int, help="Number of runs to download for each sample in the case that --target_taxID is specified. ")

debug_args.add_argument("--StopAfter_readObtentionFromSRA", dest="StopAfter_readObtentionFromSRA", action="store_true", default=False, help="Stop after obtaining reads from SRA.")

debug_args.add_argument("--StopAfter_sampleIndexingFromSRA", dest="StopAfter_sampleIndexingFromSRA", action="store_true", default=False, help="It will stop after indexing the samples of SRA. You can use this if, for example, your local machine has internet connection and your slurm cluster does not. You can first obtain the SRA indexes in the local machine. And then run again this pipeline without this option in the slurm cluster.")

debug_args.add_argument("--StopAfterPrefecth_of_reads", dest="StopAfterPrefecth_of_reads", action="store_true", default=False, help="Stop after obtaining the prefetched .srr file in case close_shortReads_table is 'auto'")

debug_args.add_argument("--StopAfterPrefecth_of_reads_goldenSet", dest="StopAfterPrefecth_of_reads_goldenSet", action="store_true", default=False, help="Stop after obtaining the prefetched .srr file in case --goldenSet_table is specified.")

debug_args.add_argument("--StopAfter_testAccuracy_perSVadeRunning", dest="StopAfter_testAccuracy_perSVadeRunning", action="store_true", default=False, help="When --testAccuracy is specified, the pipeline will stop after the running of perSVade on all the inputs of --close_shortReads_table with the different configurations.")

debug_args.add_argument("--StopAfter_testAccuracy", dest="StopAfter_testAccuracy", action="store_true", default=False, help="When --testAccuracy is specified, the pipeline will stop after testing the accuracy.")

debug_args.add_argument("--StopAfter_goldenSetAnalysis", dest="StopAfter_goldenSetAnalysis", action="store_true", default=False, help="When --goldenSet_table is specified, the pipeline will stop after running the golden set analysis.")

debug_args.add_argument("--StopAfter_replace_SV_CNVcalling", dest="StopAfter_replace_SV_CNVcalling", action="store_true", help="Stop after the removal of files for repeating the CNV calling.")

debug_args.add_argument("--testAccuracy", dest="testAccuracy", action="store_true", default=False, help="Reports the accuracy  of your calling on the real data, simulations and fastSVcalling for all the WGS runs specified in --close_shortReads_table. ")

debug_args.add_argument("--testAccuracy_skipHomRegionsSimulation", dest="testAccuracy_skipHomRegionsSimulation", action="store_true", default=False, help="If --testAccuracy is provided, skip the simulation of SVs around hom regions ")

debug_args.add_argument("--correct_close_shortReads_table_by_nRunsAndSamples", dest="correct_close_shortReads_table_by_nRunsAndSamples", action="store_true", default=False, help="Corrects --close_shortReads_table keeping only nruns and nsamples. ")

debug_args.add_argument("--skip_cleaning_outdir", dest="skip_cleaning_outdir", action="store_true", default=False, help="Will NOT remove all the unnecessary files of the perSVade outdir")

debug_args.add_argument("--max_coverage_sra_reads", dest="max_coverage_sra_reads", default=10000000000000000, type=int, help="This is the maximum coverage allowed in the reads downloaded by SRA. If the datasets have a coverage above this perSVade will randmomly subset reads to match this.")

debug_args.add_argument("--replace_addingCNstate_to_smallVars", dest="replace_addingCNstate_to_smallVars", action="store_true", default=False, help="Replace the step of adding the Copy Number of each variant.")

debug_args.add_argument("--generate_alternative_genome", dest="generate_alternative_genome", default=False, action="store_true", help="Generate an alternative genome in smallVariantCalling.")

debug_args.add_argument("--visualization_results", dest="visualization_results", default=False, action="store_true", help="Visualize the results")

debug_args.add_argument("--keep_simulation_files", dest="keep_simulation_files", action="store_true", default=False, help="Keeps the simulation files from perSVade.")

#############################################

opt = parser.parse_args()


########################################
##### GENERAL PROCESSING OF INPUTS #####
########################################

# define the modules
if opt.type_variant_calling is not None:

    # get the modules and debug
    modules_set = set(opt.type_variant_calling.split(","))
    strange_modules = modules_set.difference({"SV", "coverageCNV", "small_vars"})
    if len(strange_modules)>0: raise ValueError("There are some strange modules: %s. Provide a valid --type_variant_calling."%strange_modules)

    # check that there is some sequence input
    if any([x=="skip" for x in {opt.fastq1, opt.fastq2}]): raise ValueError("If you want to execute some modules with --type_variant_calling you have to provide reads (-f1 and -f2) or a bam file (-sbam). If you specified '-f1 skip' and/or '-f2 skip' you are telling perSVade to not use any sequence input, which is incomaptible with --type_variant_calling. ")

    # override arguments depending on the combinations
    if "SV" in modules_set: opt.skip_SVcalling = False
    else: opt.skip_SVcalling = True

    if "coverageCNV" in modules_set: opt.skip_CNV_calling = False
    else: opt.skip_CNV_calling = True

    if "small_vars" in modules_set: opt.run_smallVarsCNV = True
    else: opt.run_smallVarsCNV = False 

    # redefine some args
    if opt.run_smallVarsCNV is True: opt.remove_smallVarsCNV_nonEssentialFiles = True

# start time of processing
start_time_GeneralProcessing =  time.time()
start_time_all =  time.time()

# if replace is set remove the outdir, and then make it
if opt.replace is True: fun.delete_folder(opt.outdir)
fun.make_folder(opt.outdir)

# define the final file. and exit if it exists
final_file = "%s/perSVade_finished_file.txt"%opt.outdir
if opt.replace or opt.replace_SV_CNVcalling or opt.replace_FromGridssRun_final_perSVade_run or opt.replace_var_integration: fun.remove_file(final_file)

if not fun.file_is_empty(final_file): 
    
    print("WARNING: %s exists, suggesting that perSVade was already  run in this folder. Remove this file if you want this command to work. Exiting..."%final_file)
    sys.exit(0)

# define the name as the sample as the first 10 characters of the outdir
name_sample = fun.get_sampleName_from_perSVade_outdir(opt.outdir)
print("Running perSVade into %s. The name_sample is %s"%(opt.outdir, name_sample))

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

# check that the simulation_chromosomes are in all_chroms
if opt.simulation_chromosomes is not None: 
    strange_chroms = set(opt.simulation_chromosomes.split(",")).difference(all_chroms)
    if len(strange_chroms)>0: raise ValueError("The --simulation_chromosomes argument was not properly set. There are some chromosome IDs (%s) not found in the provided reference genome"%strange_chroms)

# get the genome len
genome_length = sum(fun.get_chr_to_len(opt.ref).values())
print("The genome has %.2f Mb"%(genome_length/1000000 ))

# Index the reference
fun.index_genome(opt.ref, replace=opt.replace)

##################################

#### REPLACE THE GFF ####
target_gff = "%s/reference_genome_features.gff"%reference_genome_dir

# if you want to repeat the annotation, remove all the files in the reference genome dir that are related to the annotation
if opt.replace_var_annotation is True:
    for f in os.listdir(reference_genome_dir): 

        file_path = "%s/%s"%(reference_genome_dir, f)
        if file_path.startswith(target_gff): fun.remove_file(file_path)

# copy the gff
if opt.gff is None: print("WARNING: gff was not provided. This will be a problem if you want to annotate small variant calls")
else:

    if fun.file_is_empty(target_gff): fun.soft_link_files(opt.gff, target_gff)

    # change the path
    opt.gff = target_gff

#########################

#### REPLACE THE REPEATS TABLE IF PROVIDED ####

# define the repeats table. This will work for any downstream analysis of this
repeats_table_file = "%s.repeats.tab"%opt.ref

if opt.previous_repeats_table is not None:
    print("using privided repeats %s"%opt.previous_repeats_table)

    if fun.file_is_empty(opt.previous_repeats_table): raise ValueError("The provided repeats table does not exist")
    
    # softlink
    fun.remove_file(repeats_table_file)
    fun.soft_link_files(opt.previous_repeats_table, repeats_table_file)

###############################################


# GENERATE READS FROM OTHER perSVade_outdirs
if opt.other_perSVade_outdirs_sameReadsANDalignment is not None: fun.link_files_from_other_perSVade_outdirs_reads_and_alignment(opt.outdir, opt.other_perSVade_outdirs_sameReadsANDalignment)

#### define misc args ####

# define a log file for perSVade
fun.log_file_all_cmds = "%s/all_cmds_ran.txt"%opt.outdir
if fun.file_is_empty(fun.log_file_all_cmds): open(fun.log_file_all_cmds, "w").write("# These are all the cmds:\n")

# get the simulation ploidies
if opt.simulation_ploidies!="auto": simulation_ploidies = opt.simulation_ploidies.split(",")

else: 

    # for pooled seq it takes simulated on 1 in 10
    if opt.pooled_sequencing is True: simulation_ploidies = ["diploid_homo", "ref:9_var:1"]
    elif opt.ploidy==1: simulation_ploidies = ["haploid"]
    else: simulation_ploidies = ["ref:%i_var:1"%(opt.ploidy-1)]

# define the CNV calling algs
cnv_calling_algs = set(opt.cnv_calling_algs.split(","))
all_expected_cnv_calling_algs = {"HMMcopy", "AneuFinder", "CONY"}
if len(cnv_calling_algs.difference(all_expected_cnv_calling_algs))>0: raise ValueError("the cnv calling algs should be in %s"%all_expected_cnv_calling_algs)

# set a specific calling for pooled sequencing
if opt.pooled_sequencing is True: print("WARNING: Running on pooled sequencing.  These are the simulated ploidies for the SV calling parameter optimisation:", simulation_ploidies)

# the window length for all operations
fun.window_l = fun.get_perSVade_window_l(opt.ref, opt.mitochondrial_chromosome, opt.min_chromosome_len)

# define the verbosity. If opt.verbose is False, none of the 'print' statements of sv_functions will have an effect
fun.printing_verbose_mode = opt.verbose

# defin the fraction of available mem
fun.fraction_available_mem = opt.fraction_available_mem
if opt.fraction_available_mem is None: print("WARNING: You did not specify how much RAM should be used through --fraction_available_mem. perSVade will calculate this by filling the memory, which may be dangerous. If you want to use all the allocated memory you should specify --fraction_available_mem 1.0")

# define the fraction of RAM to dedicate
if opt.fractionRAM_to_dedicate>0.95: raise ValueError("You are using >95 pct of the systems RAM, which is dangerous")
fun.fractionRAM_to_dedicate = opt.fractionRAM_to_dedicate

# define the min_CNVsize_coverageBased
fun.min_CNVsize_coverageBased = opt.min_CNVsize_coverageBased

# redefine the real threads
real_available_threads = fun.get_available_threads(opt.outdir)
if opt.threads>real_available_threads:  print("WARNING: There are %i available threads, and you required %i."%(real_available_threads, opt.threads))

available_Gb_RAM = fun.get_availableGbRAM(opt.outdir)
if available_Gb_RAM<opt.min_gb_RAM_required: raise ValueError("You should be running with at least %i Gb of RAM. There are just %.3f in this system"%(opt.min_gb_RAM_required, available_Gb_RAM))
print("Running with %.3f Gb of RAM and %i cores"%(available_Gb_RAM, opt.threads))

# change the default parameters if specified
if opt.parameters_json_file is not None:

    gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup = fun.get_parameters_from_json(opt.parameters_json_file)

    fun.default_filtersDict_gridss = gridss_filters_dict
    fun.default_gridss_blacklisted_regions = gridss_blacklisted_regions
    fun.default_gridss_maxcoverage = gridss_maxcoverage
    fun.default_max_rel_coverage_to_consider_del = max_rel_coverage_to_consider_del
    fun.default_min_rel_coverage_to_consider_dup = min_rel_coverage_to_consider_dup

# test whether the gff is correct
if opt.gff is not None: 
    correct_gff, gff_with_biotype = fun.get_correct_gff_and_gff_with_biotype(opt.gff, replace=opt.replace)
    fun.check_that_gff_is_correct(opt.gff, opt.ref, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code, opt.threads, opt.replace)

# check that the tmpdir exists
if opt.tmpdir is not None:
    if not os.path.isdir(opt.tmpdir): raise ValueError("The folder that you specified with --tmpdir does not exist")

# get the repeats table
if opt.skip_repeat_analysis is False:

    print("getting repeats")
    repeats_df, repeats_table_file = fun.get_repeat_maskerDF(opt.ref, threads=opt.threads, replace=opt.replace)

else:

    print("skipping the repeats analysis")
    fun.write_repeats_table_file(repeats_table_file)


if opt.StopAfter_repeatsObtention is True:
    print("Stopping after the obtention of repeats. They can be found in %s"%repeats_table_file)
    sys.exit(0)

####### GENERATE A real_bedpe_breakpoints AROUND REPEATS IF NECESSARY ##########

if opt.simulate_SVs_around_repeats is True:
    print("simulating around repeats")
    
    # debug
    if opt.real_bedpe_breakpoints is not None: raise ValueError("You can not specify --simulate_SVs_around_repeats and --real_bedpe_breakpoints. You may definir either of both.")
    if opt.skip_repeat_analysis is True: raise ValueError("You should not skip the repeats analysis (with --skip_repeat_analysis) if you want to simulate SVs around repeats.")
    if opt.simulate_SVs_around_HomologousRegions is True: raise ValueError("You can not specify --simulate_SVs_around_repeats and --simulate_SVs_around_HomologousRegions.")


    # override the real_bedpe_breakpoints with the bedpe comming from repeats
    opt.real_bedpe_breakpoints = fun.get_bedpe_breakpoints_around_repeats(repeats_table_file, replace=opt.replace, max_breakpoints=(opt.nvars*5*10000), max_breakpoints_per_repeat=1, threads=opt.threads, max_repeats=(opt.nvars*5*2500))

################################################################################

####### GENERATE real_bedpe_breakpoints AROUND HOMOLOGOUS REGIONS #########

if opt.simulate_SVs_around_HomologousRegions is True:
    print("simulating around Homologous regions")

    # debug
    if opt.real_bedpe_breakpoints is not None: raise ValueError("You can not specify --simulate_SVs_around_repeats and --real_bedpe_breakpoints. You may definie either of both.") 
    if opt.simulate_SVs_around_repeats is True: raise ValueError("You can not specify --simulate_SVs_around_repeats and --simulate_SVs_around_HomologousRegions.")

    # get a file that contains the blastn of some regions of the genome
    if opt.simulate_SVs_around_HomologousRegions_previousBlastnFile is None: blastn_file = fun.get_blastn_regions_genome_against_itself(opt.ref, opt.simulate_SVs_around_HomologousRegions_maxEvalue, opt.simulate_SVs_around_HomologousRegions_queryWindowSize, opt.replace, opt.threads)

    else: blastn_file = opt.simulate_SVs_around_HomologousRegions_previousBlastnFile

    # define the bedpe breakpoints around the homologous regions
    bedpe_breakpoints = "%s/breakpoints_aroundHomRegions_wsize=%ibp_maxEval=%s_minQcovS=%i.bedpe"%(opt.outdir, opt.simulate_SVs_around_HomologousRegions_queryWindowSize, opt.simulate_SVs_around_HomologousRegions_maxEvalue, opt.simulate_SVs_around_HomologousRegions_minPctOverlap)

    opt.real_bedpe_breakpoints = fun.get_bedpe_breakpoints_around_homologousRegions(blastn_file, bedpe_breakpoints, replace=opt.replace, threads=opt.threads, max_eval=opt.simulate_SVs_around_HomologousRegions_maxEvalue, query_window_size=opt.simulate_SVs_around_HomologousRegions_queryWindowSize, min_qcovs=opt.simulate_SVs_around_HomologousRegions_minPctOverlap, max_n_hits=(opt.nvars*5*2500))

############################################################################

#############################

# end time of processing
end_time_GeneralProcessing =  time.time()

########################################
########################################
########################################


#####################################
############# BAM FILE ##############
#####################################

start_time_alignment =  time.time()

if not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]):

    ##### DEFINE THE SORTED BAM #####

    # define the dest bam files
    bamfile = "%s/aligned_reads.bam"%opt.outdir
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam

    # if you provided a sorted bam it should be placed into sorted_bam
    if opt.sortedbam is not None:

        # debug the fact that you prvided reads and bam. You should just provide one
        if any([not x is None for x in {opt.fastq1, opt.fastq2}]): raise ValueError("You have provided reads and a bam, you should only provide one")

        # downsample the bam for testing if speciefied
        if opt.downsampled_coverage is not None: 
            print("WARNING: You are running perSVade on a downsampled bam (to %.3fx)."%opt.downsampled_coverage)
            opt.sortedbam = fun.get_downsampled_bam_to_specified_coverage(opt.outdir, opt.sortedbam, opt.downsampled_coverage, opt.ref, replace=opt.replace, threads=opt.threads) 

        # get the linked bam files
        fun.soft_link_files(fun.get_fullpath(opt.sortedbam), sorted_bam)
        fun.soft_link_files(fun.get_fullpath(opt.sortedbam)+".bai", sorted_bam+".bai")
        
    ###################################

    # normal alignment of provided reads
    if all([not x is None for x in {opt.fastq1, opt.fastq2}]) or opt.input_SRRfile is not None:

        # if you have reads, you may want to downsample them
        if all([not x is None for x in {opt.fastq1, opt.fastq2}]) and opt.downsampled_coverage is not None:

            # downsample the reads
            print("WARNING: You are running perSVade on downsampled reads (to %.3fx)."%opt.downsampled_coverage)
            opt.fastq1, opt.fastq2 = fun.get_paired_downsampled_reads(opt.fastq1, opt.fastq2, "%s/downsampled_reads"%opt.outdir, opt.downsampled_coverage, opt.ref, replace=opt.replace, threads=opt.threads)

        # if the reads have to be QC and trimmed:
        if opt.QC_and_trimming_reads is True or opt.input_SRRfile is not None: 

            # define the reads dir
            reads_dir = "%s/reads"%opt.outdir; fun.make_folder(reads_dir)

            # define the raw reads under reads dir
            dest_fastq1 = "%s/raw_reads1.fastq.gz"%reads_dir
            dest_fastq2 = "%s/raw_reads2.fastq.gz"%reads_dir

            # define the trimmed reads
            trimmed_reads1 = "%s.trimmed.fastq.gz"%dest_fastq1
            trimmed_reads2 = "%s.trimmed.fastq.gz"%dest_fastq2

            if fun.file_is_empty(trimmed_reads1) or fun.file_is_empty(trimmed_reads2) or opt.replace is True:
                print("running trimming and QC of the reads")

                # get reads from prefetch
                if opt.input_SRRfile is not None and opt.fastq1 is None and opt.fastq2 is None:
                    print("Getting raw reads from SRR file")

                    dest_input_SRRfile = "%s/input_SRRfile.srr"%reads_dir
                    fun.soft_link_files(opt.input_SRRfile, dest_input_SRRfile)
                    opt.fastq1, opt.fastq2 = fun.run_parallelFastqDump_on_prefetched_SRRfile(dest_input_SRRfile, replace=opt.replace, threads=opt.threads)

                # get the reads under outdir
                fun.soft_link_files(opt.fastq1, dest_fastq1)
                fun.soft_link_files(opt.fastq2, dest_fastq2)

            # trim reads
            opt.fastq1, opt.fastq2 = fun.run_trimmomatic(dest_fastq1, dest_fastq2, replace=opt.replace, threads=opt.threads)

            # clean
            for f in os.listdir(reads_dir): 
                if f not in {fun.get_file(opt.fastq1), fun.get_file(opt.fastq2)}: fun.delete_file_or_folder("%s/%s"%(reads_dir, f))


        # deifine if marking duplicates (default yes)
        if opt.skip_marking_duplicates is True: bwa_mem_MarkDuplicates = False
        else: bwa_mem_MarkDuplicates = True

        print("WORKING ON ALIGNMENT")
        fun.run_bwa_mem(opt.fastq1, opt.fastq2, opt.ref, opt.outdir, bamfile, sorted_bam, index_bam, name_sample, threads=opt.threads, replace=opt.replace, tmpdir_writingFiles=opt.tmpdir, MarkDuplicates=bwa_mem_MarkDuplicates)
        fun.clean_sorted_bam_coverage_per_window_files(sorted_bam)
        
    else: print("Warning: No fastq file given, assuming that you provided a bam file")

    # check that all the important files exist
    if any([fun.file_is_empty(x) for x in {sorted_bam, index_bam}]): raise ValueError("You need the sorted and indexed bam files in ")

end_time_alignment =  time.time()

#####################################
#####################################
#####################################

###########################################
############# NECESSARY FILES #############
###########################################

# First create some files that are important for any program

# Create a reference dictionary
fun.create_sequence_dict(opt.ref, replace=opt.replace)

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

start_time_obtentionCloseSVs =  time.time()

##### find a set of 'real_bedpe_breakpoints' that will be used for the simulations #####

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

    # redefine close_shortReads_table to match n_close_samples and nruns_per_sample
    if opt.correct_close_shortReads_table_by_nRunsAndSamples is True: opt.close_shortReads_table = fun.get_redefined_close_shortReads_table_with_meaningful_samples(opt.close_shortReads_table, n_close_samples=opt.n_close_samples, nruns_per_sample=opt.nruns_per_sample)

    # skip the running of the pipeline 
    if opt.StopAfter_readObtentionFromSRA:
        print("Stopping pipeline after the reads obtention from SRA")
        sys.exit(0) 

    # skip pipeline running if you have to stop after prefetch of reads
    if opt.StopAfterPrefecth_of_reads:
        print("Stopping pipeline after the prefetch of reads")
        sys.exit(0) 

    # get the real SVs
    real_bedpe_breakpoints = fun.get_compatible_real_bedpe_breakpoints(opt.close_shortReads_table, opt.ref, outdir_finding_realVars, replace=opt.replace, threads=opt.threads, max_nvars=opt.nvars, mitochondrial_chromosome=opt.mitochondrial_chromosome, job_array_mode=opt.job_array_mode, parameters_json_file=opt.parameters_json_file, tmpdir=opt.tmpdir, skip_marking_duplicates=opt.skip_marking_duplicates)

else: 
    print("Avoiding the simulation of real variants. Only inserting randomSV.")

    # define the set of vars as empty. This will trigger the random generation of vars
    real_bedpe_breakpoints = None


if opt.StopAfter_obtentionOFcloseSVs: 
    print("stopping pipeline after obtention of close SVs")
    sys.exit(0)

end_time_obtentionCloseSVs =  time.time()

# test that the real_bedpe_breakpoints are correct
if real_bedpe_breakpoints is not None and not os.path.isfile(real_bedpe_breakpoints): raise ValueError("The provided real bedpe breakpoints %s are incorrect."%real_bedpe_breakpoints)

###################################################################################################

# test accuracy on real data
if opt.testAccuracy is True:  

    print("testing accuracy on simulations and real variants (if provided)")

    # test that you have provided a opt.close_shortReads_table
    if opt.close_shortReads_table is None or opt.fast_SVcalling is True: 
        raise ValueError("You have to specify a --close_shortReads_table or --real_bedpe_breakpoints and not run in --fast_SVcalling to test the accuracy of the pipeline on several datasets (--testAccuracy)")

    ### RUN PERSVADE ###
    print("testing perSVade on several samples. --tmpdir is %s"%opt.tmpdir)

    dict_perSVade_outdirs = fun.report_accuracy_realSVs_perSVadeRuns(opt.close_shortReads_table, opt.ref, "%s/testing_Accuracy"%opt.outdir, real_bedpe_breakpoints, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, job_array_mode=opt.job_array_mode, parameters_json_file=opt.parameters_json_file, gff=opt.gff, replace_FromGridssRun_final_perSVade_run=opt.replace_FromGridssRun_final_perSVade_run, fraction_available_mem=opt.fraction_available_mem, replace_SV_CNVcalling=opt.replace_SV_CNVcalling, skip_CNV_calling=opt.skip_CNV_calling, outdir_finding_realVars=outdir_finding_realVars, simulate_SVs_around_HomologousRegions_previousBlastnFile=opt.simulate_SVs_around_HomologousRegions_previousBlastnFile, simulate_SVs_around_HomologousRegions_maxEvalue=opt.simulate_SVs_around_HomologousRegions_maxEvalue, simulate_SVs_around_HomologousRegions_queryWindowSize=opt.simulate_SVs_around_HomologousRegions_queryWindowSize, skip_SV_CNV_calling=opt.skip_SV_CNV_calling, tmpdir=opt.tmpdir, simulation_chromosomes=opt.simulation_chromosomes, testAccuracy_skipHomRegionsSimulation=opt.testAccuracy_skipHomRegionsSimulation)

    if opt.StopAfter_testAccuracy_perSVadeRunning is True: 
        print("You already ran all the configurations of perSVade. Stopping after the running of perSVade on testAccuracy")
        sys.exit(0)

    ####################   

    ### REPORT ACCURACY SINGLE SAMPLE ###
    #youhavetoaddcodeof_codeGraveyard_report_accuracy_realSVs


    if opt.StopAfter_testAccuracy is True: 
        print(" Stopping after the running of perSVade on testAccuracy")
        sys.exit(0)

    ##################################### 


# get the golden set
if opt.goldenSet_table is not None:

    # run jobs golden set testing
    outdir_goldenSet = "%s/testing_goldenSetAccuracy"%opt.outdir
    dict_paths_goldenSetAnalysis = fun.report_accuracy_golden_set_runJobs(opt.goldenSet_table, outdir_goldenSet, opt.ref, real_bedpe_breakpoints, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, job_array_mode=opt.job_array_mode, StopAfterPrefecth_of_reads=opt.StopAfterPrefecth_of_reads_goldenSet, target_taxID=opt.target_taxID, parameters_json_file=opt.parameters_json_file, fraction_available_mem=opt.fraction_available_mem, verbose=opt.verbose, simulate_SVs_around_HomologousRegions_previousBlastnFile=opt.simulate_SVs_around_HomologousRegions_previousBlastnFile, simulate_SVs_around_HomologousRegions_maxEvalue=opt.simulate_SVs_around_HomologousRegions_maxEvalue, simulate_SVs_around_HomologousRegions_queryWindowSize=opt.simulate_SVs_around_HomologousRegions_queryWindowSize, max_n_samples=opt.goldenSet_max_n_samples)

    # plot the accurac
    fun.report_accuracy_golden_set_reportAccuracy(dict_paths_goldenSetAnalysis, outdir_goldenSet, opt.ref, threads=opt.threads, replace=opt.replace)

    if opt.StopAfter_goldenSetAnalysis is True: 
        print(" Stopping after the running of golden-set analysis")
        sys.exit(0) 

start_time_SVcalling =  time.time()

# run the actual perSVade function optimising parameters
if opt.skip_SVcalling is False and not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]):

    SVdetection_outdir = "%s/SVdetection_output"%opt.outdir
    outdir_gridss_final = fun.run_GridssClove_optimising_parameters(sorted_bam, opt.ref, SVdetection_outdir, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, fast_SVcalling=opt.fast_SVcalling, real_bedpe_breakpoints=real_bedpe_breakpoints, replace_FromGridssRun_final_perSVade_run=opt.replace_FromGridssRun_final_perSVade_run, tmpdir=opt.tmpdir, simulation_chromosomes=opt.simulation_chromosomes)

end_time_SVcalling =  time.time()

print("structural variation analysis with perSVade finished")

#####################################
#####################################
#####################################

###########################
###### CNV CALLING ########
###########################

# run CNVcalling by getting df_CNV_coverage
minimal_CNV_fields = ["chromosome", "merged_relative_CN", "start", "end", "CNVid", "median_coverage", "median_coverage_corrected", "SVTYPE"] + ["median_relative_CN_%s"%x for x in cnv_calling_algs]

run_CNV_calls = (opt.skip_CNV_calling is False and not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]))
if run_CNV_calls is True: 
    print("RUNNING COVERAGE-BASED CNV CALLING...")

    # make folder
    cnv_calling_outdir = "%s/CNV_calling"%opt.outdir
    fun.make_folder(cnv_calling_outdir)

    # run CNV calling
    df_CNV_coverage = fun.run_CNV_calling(sorted_bam, opt.ref, cnv_calling_outdir, opt.threads, opt.replace, opt.mitochondrial_chromosome, opt.window_size_CNVcalling, opt.ploidy, bg_sorted_bam_CNV=opt.bg_sorted_bam_CNV, cnv_calling_algs=cnv_calling_algs)

else: df_CNV_coverage = pd.DataFrame(columns=minimal_CNV_fields)

###########################
###########################
###########################

#####################################
###### SV and CNV ANNOTATION ########
#####################################

start_time_SVandCNVcalling =  time.time()

run_SV_CNV_calling = (opt.skip_SVcalling is False and not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]) and opt.skip_SV_CNV_calling is False)
if run_SV_CNV_calling is True:

    print("running CNV calling per window and integrating to SV calling")

    # define outdirs
    outdir_var_calling = "%s/SVcalling_output"%opt.outdir

    # remove folders if there is some replacement to be done. Remove
    if opt.replace_SV_CNVcalling is True: fun.delete_folder(outdir_var_calling)

    # stop after the removal
    if opt.StopAfter_replace_SV_CNVcalling is True: 
        print("exitting after the --replace_SV_CNVcalling action")
        sys.exit(0)

    # make folder
    fun.make_folder(outdir_var_calling)
    
    # get the variant calling 
    SV_CNV_vcf = fun.get_vcf_all_SVs_and_CNV(opt.outdir, outdir_var_calling, sorted_bam, opt.ref, opt.ploidy, df_CNV_coverage, opt.window_size_CNVcalling, cnv_calling_algs, replace=opt.replace, threads=opt.threads, mitochondrial_chromosome=opt.mitochondrial_chromosome)

    print("the SV and CNV calling vcf can be found in %s"%SV_CNV_vcf)

    # get variant annotation
    if opt.gff is not None:

        # remove the annotated vcf if needed 
        if opt.replace_var_annotation is True: 
            fun.remove_file("%s_annotated_VEP.tab"%SV_CNV_vcf)
            fun.remove_file("%s_annotated_VEP.tab.raw.tbl.tmp_summary.html"%SV_CNV_vcf)

        print("annotating SV, CNV variants with VEP")
        SV_CNV_vcf_annotated = fun.annotate_SVs_inHouse(SV_CNV_vcf, gff_with_biotype, opt.ref, replace=opt.replace, threads=opt.threads, mitochondrial_chromosome=opt.mitochondrial_chromosome, mito_code=opt.mitochondrial_code, gDNA_code=opt.gDNA_code)

        print("annotated SV vcf can be found in %s"%SV_CNV_vcf_annotated)
    
    else: print("WARNING: Skipping SV annotation because -gff was not provided.")

end_time_SVandCNVcalling =  time.time()

#####################################
#####################################
#####################################

# stop after the generation of SV and CNV calls

#####################################
###### SMALL VARS AND CNV ###########
#####################################

start_time_smallVarsCNV =  time.time()

if opt.run_smallVarsCNV:

    # define an outdir
    outdir_varcall = "%s/smallVars_CNV_output"%opt.outdir

    # define the basic cmd
    varcall_cmd = "%s -r %s --threads %i --outdir %s -sbam %s --caller %s --coverage %i --mitochondrial_chromosome %s --mitochondrial_code %i --gDNA_code %i --minAF_smallVars %s --window_freebayes_bp %i --log_file_all_cmds %s"%(varcall_cnv_pipeline, opt.ref, opt.threads, outdir_varcall, sorted_bam, opt.caller, opt.coverage, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code, opt.minAF_smallVars, opt.window_freebayes_bp, fun.log_file_all_cmds)

    # add options
    if opt.replace is True: varcall_cmd += " --replace"
    if opt.gff is not None: varcall_cmd += " -gff %s"%opt.gff
    if opt.StopAfter_smallVarCallSimpleRunning is True: varcall_cmd += " --StopAfter_smallVarCallSimpleRunning"
    if opt.replace_var_integration is True: varcall_cmd += " --replace_var_integration"
    if opt.pooled_sequencing is True: varcall_cmd += " --pooled_sequencing"
    if opt.consider_repeats_smallVarCall is True: varcall_cmd += " --repeats_table %s"%repeats_table_file
    if opt.generate_alternative_genome is True: varcall_cmd += " --generate_alternative_genome"
    if opt.skip_cnv_per_gene_analysis is True: varcall_cmd += " --skip_cnv_analysis"

    # define which ploidies to run
    if opt.ploidy==1 and opt.run_ploidy2_ifHaploid is False: ploidies_varcall = [1]
    if opt.ploidy==1 and opt.run_ploidy2_ifHaploid is True: ploidies_varcall = [1, 2]
    else: ploidies_varcall = [opt.ploidy]

    # run for each ploidy
    for ploidy_varcall in ploidies_varcall:

        # remove the annotation-derived outputs
        if opt.replace_var_annotation is True: 
            fun.remove_file("%s/variant_annotation_ploidy%i.tab"%(outdir_varcall, ploidy_varcall))
            fun.delete_folder("%s/CNV_results"%outdir_varcall)
            
        # run the variant calling command
        varcall_cmd += " -p %i"%ploidy_varcall
        if __name__ == '__main__': fun.run_cmd(varcall_cmd)

        # regenerate the variant calling file according to run_SV_CNV_calling
        if run_CNV_calls is True: fun.get_small_variant_calling_withCNstate("%s/variant_calling_ploidy%i.tab"%(outdir_varcall, ploidy_varcall), df_CNV_coverage, replace=(opt.replace or opt.replace_addingCNstate_to_smallVars))
  
    # define the small variants vcf
    small_vars_vcf = "%s/variants_atLeast1PASS_ploidy%i.vcf"%(outdir_varcall, opt.ploidy)

    # define the variant annotation
    small_vars_var_annotation = "%s/variant_annotation_ploidy%i.tab"%(outdir_varcall, opt.ploidy)

    # clean the varcall dir if specified
    if opt.remove_smallVarsCNV_nonEssentialFiles is True: fun.remove_smallVarsCNV_nonEssentialFiles_severalPloidies(outdir_varcall, ploidies_varcall)

end_time_smallVarsCNV =  time.time()

if opt.StopAfter_smallVarCall is True:
    print("WARNING: Ending after the running of small variant calling...")
    sys.exit(0)

#####################################
#####################################
#####################################


#####################################
##### VARIANTS VISUALIZATION ########
#####################################

if (opt.skip_SVcalling is False or opt.run_smallVarsCNV is True) and not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]) and opt.visualization_results is True and opt.gff is not None:

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

# keep the simulation files
if opt.keep_simulation_files is True: fun.keep_simulation_files_for_perSVade_outdir(opt.outdir, replace=opt.replace, n_simulated_genomes=opt.nsimulations, simulation_ploidies=simulation_ploidies)

# at the end you want to clean the outdir to keep only the essential files
if opt.skip_cleaning_outdir is False: fun.clean_perSVade_outdir(opt.outdir)

# define times
end_time_all = time.time()

# generate a file that indicates whether the gridss run is finished
fun.generate_final_file_report(final_file, start_time_GeneralProcessing, end_time_GeneralProcessing, start_time_alignment, end_time_alignment, start_time_all, end_time_all, start_time_obtentionCloseSVs, end_time_obtentionCloseSVs, start_time_SVcalling, end_time_SVcalling, start_time_SVandCNVcalling, end_time_SVandCNVcalling, start_time_smallVarsCNV, end_time_smallVarsCNV)

print("perSVade Finished correctly")

