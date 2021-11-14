#!/usr/bin/env python

######### define environment ##########

# this is to test that perSVade was properly installed. It takes ~45 minutes on a system of 4 cores, 32Gb of RAM

# module imports
import sys
import os

# get the cwd were all the scripts are 
test_dir = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, test_dir)
scripts_dir = "%s/../../scripts"%test_dir; sys.path.insert(0, scripts_dir)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import test_functions as test_fun

# define the testing inputs (based on the first 2 genes of two chromosomes of Cglabrata)
testing_inputs_dir = "%s/testing_inputs"%test_dir
test_ref_genome = "%s/reduced_genome.fasta"%testing_inputs_dir
test_mutated_genome = "%s/reduced_genome_mutated.fasta"%testing_inputs_dir
test_gff = "%s/reduced_annotation.gff"%testing_inputs_dir

# load the functions (test if you can import python packages)
import sv_functions as fun
fun.print_with_runtime("loading python packages worked successfully")
fun.printing_verbose_mode = False # this may be changed

# define the testing inuts dir 
testing_outputs_dir = "%s/testing_outputs"%test_dir # this is the normal place
#testing_outputs_dir = "%s/Desktop/testing_perSVade_outputs"%(os.getenv("HOME")) # this is to be faster
test_output_perSVade = "%s/perSVade_output"%testing_outputs_dir
outdir_small_variantCalling = "%s/smallVars_CNV_output"%test_output_perSVade

fun.print_with_runtime("all output files will be written to %s"%testing_outputs_dir)

# delete and cretae outdir
#fun.delete_folder(testing_outputs_dir)
fun.make_folder(testing_outputs_dir)
fun.make_folder(test_output_perSVade)

fun.print_with_runtime("setting all the files")

# redefine the reference genome location
ref_genome = "%s/reduced_genome.fasta"%testing_outputs_dir
if fun.file_is_empty(ref_genome): fun.run_cmd("cp %s %s"%(test_ref_genome, ref_genome))

# define reads that come from simulating reads on ref_genome with simulated SVs
sim_reads1 = "%s/all_reads1.correct.fq.gz"%testing_inputs_dir
sim_reads2 = "%s/all_reads2.correct.fq.gz"%testing_inputs_dir

# define a close_shortReads_table
close_shortReads_table_df = fun.pd.DataFrame({I : {"sampleID":"sample_%i"%I, "short_reads1":sim_reads1, "short_reads2":sim_reads2} for I in range(2)}).transpose()
close_shortReads_table = "%s/close_shortReads_table.tab"%testing_outputs_dir
fun.save_df_as_tab(close_shortReads_table_df, close_shortReads_table)

# define randomly subsampled C. glabrata reads
cglab_reads1 = "%s/Cglabrata_reads/sampled_readsR1_first100k.fq.gz"%testing_inputs_dir
cglab_reads2 = "%s/Cglabrata_reads/sampled_readsR2_first100k.fq.gz"%testing_inputs_dir

# redefine the gff
gff = "%s/reduced_annotations.gff"%testing_outputs_dir
if fun.file_is_empty(gff): fun.run_cmd("cp %s %s"%(test_gff, gff))

# redefine the mutated genome location
mut_genome = "%s/mutated_genome.fasta"%testing_outputs_dir
if fun.file_is_empty(mut_genome): fun.run_cmd("cp %s %s"%(test_mutated_genome, mut_genome))

# define an example calbicans varCall_outout
Calbicans_varCall_outdir = "%s/varcalling_output_Calbicans_SRR2088862"%testing_inputs_dir

# define the Calbicans genome
inputs_Calbicans_genome = "%s/Candida_albicans.fasta"%testing_inputs_dir
Calbicans_genome = "%s/Candida_albicans.fasta"%testing_outputs_dir
fun.soft_link_files(inputs_Calbicans_genome, Calbicans_genome)

# define the Calbicans repeats
inputs_Calbicans_repeats = "%s/Candida_albicans.fasta.repeats.tab"%testing_inputs_dir
Calbicans_repeats = "%s/Candida_albicans.fasta.repeats.tab"%testing_outputs_dir
fun.soft_link_files(inputs_Calbicans_repeats, Calbicans_repeats)

# define the MERV genome
inputs_MERS_genome = "%s/MERS_CoV_genome.fasta"%testing_inputs_dir
MERS_genome = "%s/MERS_CoV_genome.fasta"%testing_outputs_dir
fun.soft_link_files(inputs_MERS_genome, MERS_genome)

# define the C. albicans chromosome
inputs_Calbicans_chr1 = "%s/Candida_albicans_Ca22chr1A_C_albicans_SC5314.fasta"%testing_inputs_dir
Calbicans_chr1 = "%s/Candida_albicans_Ca22chr1A_C_albicans_SC5314.fasta"%testing_outputs_dir
fun.soft_link_files(inputs_Calbicans_chr1, Calbicans_chr1)

inputs_Calbicans_chr1_2_6 = "%s/Candida_albicans_chr1_2_6.fasta"%testing_inputs_dir
Calbicans_chr1_2_6 = "%s/Candida_albicans_chr1_2_6.fasta"%testing_outputs_dir
fun.soft_link_files(inputs_Calbicans_chr1_2_6, Calbicans_chr1_2_6)

# define the threads
threads = fun.multiproc.cpu_count()

########################################

######## TEST DIFFERENT MODULES #########

# align reads with SVs
outdir_align_reads_SVs = "%s/align_reads_sim_SVs"%testing_outputs_dir
fun.run_cmd("%s align_reads --threads %i --fraction_available_mem 1.0 -f1 %s -f2 %s -o %s --ref %s --min_chromosome_len 100"%(fun.perSVade_modules, threads, sim_reads1, sim_reads2, outdir_align_reads_SVs, ref_genome))

# run CNV calling
outdir_CNVcalling = "%s/call_CNVs"%testing_outputs_dir
fun.run_cmd("%s call_CNVs --threads %i --fraction_available_mem 1.0 -o %s --ref %s --min_chromosome_len 10000 -sbam %s/aligned_reads.bam.sorted --mitochondrial_chromosome mito_C_glabrata_CBS138 -p 1 --cnv_calling_algs HMMcopy,CONY,AneuFinder --min_CNVsize 500 --window_size_CNVcalling 500 --verbose"%(fun.perSVade_modules, threads, outdir_CNVcalling, ref_genome, outdir_align_reads_SVs))

# SV calling
outdir_SVcalling = "%s/call_SVs"%testing_outputs_dir
fun.run_cmd("%s call_SVs --threads %i --fraction_available_mem 1.0 -o %s --ref %s --min_chromosome_len 10000 -sbam %s/aligned_reads.bam.sorted --mitochondrial_chromosome mito_C_glabrata_CBS138 --SVcalling_parameters default --repeats_file skip"%(fun.perSVade_modules, threads, outdir_SVcalling, ref_genome, outdir_align_reads_SVs))

# integrate SV and CNV calls
outdir_merged_calling = "%s/integrate_SV_and_CNV"%testing_outputs_dir
fun.run_cmd("%s integrate_SV_CNV_calls --threads %i --fraction_available_mem 1.0 -o %s --ref %s --min_chromosome_len 10000  --mitochondrial_chromosome mito_C_glabrata_CBS138 -sbam %s/aligned_reads.bam.sorted -p 1 --outdir_callSVs %s --outdir_callCNVs %s --verbose --repeats_file skip"%(fun.perSVade_modules, threads, outdir_merged_calling, ref_genome, outdir_align_reads_SVs, outdir_SVcalling, outdir_CNVcalling))

print(outdir_merged_calling)
djhgdajgad

# get regions with known SVs
outdir_known_regions = "%s/known_SVs"%testing_outputs_dir
fun.run_cmd("%s find_knownSVs_regions --threads %i --fraction_available_mem 1.0 -o %s --ref %s --min_chromosome_len 10000 --mitochondrial_chromosome mito_C_glabrata_CBS138 --SVcalling_parameters default --repeats_file skip --close_shortReads_table %s --skip_marking_duplicates"%(fun.perSVade_modules, threads, outdir_known_regions, ref_genome, close_shortReads_table))

# run parameter optimization on pre-defined regions (knownSVs) and 2 chromosomes
outdir_parameter_optimization_known = "%s/parameter_optimization_known"%testing_outputs_dir
fun.run_cmd("%s optimize_parameters --threads %i --fraction_available_mem 1.0 -o %s --ref %s --min_chromosome_len 1000 --sortedbam %s/aligned_reads.bam.sorted --mitochondrial_chromosome mito_C_glabrata_CBS138 --repeats_file skip --regions_SVsimulations %s/knownSVs_breakpoints.bedpe --simulation_ploidies haploid --nsimulations 1 --nvars 5 --simulation_chromosomes ChrA_C_glabrata_CBS138,ChrB_C_glabrata_CBS138 --verbose"%(fun.perSVade_modules, threads, outdir_parameter_optimization_known, ref_genome, outdir_align_reads_SVs, outdir_known_regions))

print(outdir_parameter_optimization_known)

adkjghaj
close_shortReads_table

# get homologous regions in C. albicans
outdir_homRegions = "%s/find_hom_regions"%testing_outputs_dir
fun.run_cmd("%s find_homologous_regions --threads %i --fraction_available_mem 1.0 -o %s --ref %s --min_chromosome_len 1000"%(fun.perSVade_modules, threads, outdir_homRegions, Calbicans_chr1_2_6))

# check that that the database has been created
repeat_masker_db = "%s/Libraries/RepeatMasker.lib.nsq"%(fun.repeatmasker_dir) 
if fun.file_is_empty(repeat_masker_db): raise ValueError("%s is missing. Check that you ran ./installation/setup_environment.sh"%repeat_masker_db)

# get repeats for a large genome
#outdir_repeats = "%s/repeats_infer_Calbicans"%testing_outputs_dir
#fun.run_cmd("%s infer_repeats --threads %i --fraction_available_mem 1.0 -o %s --ref %s --min_chromosome_len 1000 --verbose"%(fun.perSVade_modules, threads, outdir_repeats, Calbicans_chr1_2_6))

# get repeats small genome
outdir_repeats_fast = "%s/repeats_infer_Cglab"%testing_outputs_dir
fun.run_cmd("%s infer_repeats --threads %i --fraction_available_mem 1.0 -o %s --ref %s --min_chromosome_len 100"%(fun.perSVade_modules, threads, outdir_repeats_fast, ref_genome))

# run parameter optimisation based on random regions
outdir_parameter_optimization = "%s/parameter_optimization"%testing_outputs_dir
fun.run_cmd("%s optimize_parameters --threads %i --fraction_available_mem 1.0 -o %s --ref %s --min_chromosome_len 1000 --sortedbam %s/aligned_reads.bam.sorted --mitochondrial_chromosome mito_C_glabrata_CBS138 --repeats_file %s/combined_repeats.tab --regions_SVsimulations random --simulation_ploidies diploid_hetero --nsimulations 1 --nvars 5"%(fun.perSVade_modules, threads, outdir_parameter_optimization, ref_genome, outdir_align_reads_SVs, outdir_repeats_fast))

# trimming of reads
outdir_trimmed_reads = "%s/read_trimmingCglab"%testing_outputs_dir
fun.run_cmd("%s trim_reads_and_QC --threads %i --fraction_available_mem 1.0 -f1 %s -f2 %s -o %s"%(fun.perSVade_modules, threads, cglab_reads1, cglab_reads2, outdir_trimmed_reads))

# align reads
outdir_align_reads = "%s/align_reads_Cglab"%testing_outputs_dir
fun.run_cmd("%s align_reads --threads %i --fraction_available_mem 1.0 -f1 %s -f2 %s -o %s --ref %s --min_chromosome_len 100"%(fun.perSVade_modules, threads, cglab_reads1, cglab_reads2, outdir_align_reads, ref_genome))






print(outdir_homRegions)
jdhgjhda
#test_fun.test_get_repeat_maskerDF(Calbicans_chr1_2_6, replace=False)

#test_fun.test_get_repeat_maskerDF(ref_genome, replace=False) # this will not work for repeat masker


print(outdir_repeats)


#########################################

fun.print_with_runtime("\n\n---\nSUCCESS: perSVade was properly installed\n---\n\n")
adkjahkjadd




JDAVJGAJGDA


######## TEST SPECIFIC PARTS ##########

if len(sys.argv)>1:
	part_test = sys.argv[1]

	# check that the picard environment is fine
	if part_test=="picard_env": test_fun.test_picard_env(testing_inputs_dir, testing_outputs_dir)

	else: raise ValueError("%s is not valid"%part_test)

	# exit
	print("testing %s worked well"%part_test)
	sys.exit(0)

#######################################

###### GENERAL TESTING ########

# calculate the memory that is available
#fun.print_with_runtime("running testings on %s Gb of RAM"%(fun.get_availableGbRAM(testing_outputs_dir)))

# check that that the database has been created
repeat_masker_db = "%s/Libraries/RepeatMasker.lib.nsq"%(fun.repeatmasker_dir) 
if fun.file_is_empty(repeat_masker_db): raise ValueError("%s is missing. Check that you ran ./installation/setup_environment.sh"%repeat_masker_db)

# test that the environment can be recreated. This is only important from the developer's package
#test_fun.test_conda_env_generation(testing_outputs_dir, replace=True)

# test repeat masker obtention
test_fun.test_get_repeat_maskerDF(ref_genome, replace=False) # this will not work for repeat masker

# test read simulation by simulating reads from the mutated genome
r1_mutGenome, r2_mutGenome = test_fun.test_read_simulation_and_get_reads(mut_genome)

# test bwa mem, samtools and picard
sorted_bam_mutGenome = test_fun.test_bwa_mem_and_get_bam(r1_mutGenome, r2_mutGenome, ref_genome)

# test the joining and unjoining of multiallleles. This is a test for bcftools 1.10
outdir_testing_CalbicansVarCall = "%s/testing_CalbicansVarCall"%testing_outputs_dir
test_fun.test_processing_varcalling(Calbicans_varCall_outdir, Calbicans_genome, outdir_testing_CalbicansVarCall, sorted_bam_mutGenome)

# test bcftools, freebayes, gatk4, mosdepth, vep by running the small variant calling pipeline
test_fun.test_smallVarCall_CNV_running(sorted_bam_mutGenome, outdir_small_variantCalling, ref_genome, gff, replace=False)

# generate a genome that is rearranged 
rearranged_genome = test_fun.test_rearranging_genome_random(ref_genome, replace=False)

# generate simulated reads from this rearranged genome, and also a bam file
r1_svGenome, r2_svGenome = test_fun.test_read_simulation_and_get_reads(rearranged_genome, replace=False)
sorted_bam_svGenome = test_fun.test_bwa_mem_and_get_bam(r1_svGenome, r2_svGenome, ref_genome, replace=False)

# test whether you can run the gridss and clove pipeline
outdir_testing_gridss_clove = "%s/testing_gridss_clove_pipeline_default_parms"%(testing_outputs_dir)
test_fun.test_gridss_clove_pipeline(sorted_bam_svGenome, ref_genome, outdir_testing_gridss_clove, replace=False)

# test whether the parameter optimisation pipeline works well to find the rearrangements
outdir_testing_parameterOptimisation = "%s/testing_parameter_optimisation_pipeline"%(testing_outputs_dir)
test_fun.test_parameter_optimisation_perSVade(sorted_bam_svGenome, ref_genome, outdir_testing_parameterOptimisation, replace=False)

# test the integration of variants from previous sequencing data
Cglabrata_genome = "%s/Candida_glabrata.fasta"%testing_inputs_dir
Cglabrata_subsampled_reads_dir = "%s/Cglabrata_reads"%testing_inputs_dir
Cglabrata_subsampled_perSVade_outdir = "%s/Cglabrata_reads_subsampled_perSVade_outdir"%testing_outputs_dir
Cglabrata_repeats = "%s/Candida_glabrata.repeats.tab"%testing_inputs_dir
relaxed_parms = "%s/perSVade_parameters_relaxed.json"%testing_inputs_dir
real_bedpe_breakpoints = test_fun.test_realSVgeneration(Cglabrata_subsampled_reads_dir, Cglabrata_subsampled_perSVade_outdir, Cglabrata_repeats, Cglabrata_genome, relaxed_parms, replace=False)

# test the parameter optimisation based on real_bedpe_breakpoints
#test_fun.test_parameter_optimisation_perSVade_real(Cglabrata_subsampled_reads_dir, Cglabrata_subsampled_perSVade_outdir, Cglabrata_repeats, Cglabrata_genome, real_bedpe_breakpoints, replace=False)

# test CNV calling with each of the three programs (CONY is not tested because it should not be used by default)

### TESTING THINGS THAT ARE DISPENSABLE ###

# test the querying of the SRA database and downloading, and trimming of reads and also gztool
#outdir_SRAdb_query_downloading_and_readTrimming = "%s/testing_SRAdb_query_downloading_and_readTrimming_MERS"%testing_outputs_dir
#test_fun.test_SRAdb_query_downloading_and_readTrimming(outdir_SRAdb_query_downloading_and_readTrimming, MERS_genome, 1335626, replace=False)

# test the repeat masker obtention for a long chromosome 1, 2 and 6
#test_fun.test_get_repeat_maskerDF(Calbicans_chr1_2_6, replace=False)


fun.print_with_runtime("\n\n---\nSUCCESS: perSVade was properly installed\n---\n\n")

###############################

############################################

fun.print_with_runtime("\n\n---\nSUCCESS: perSVade was properly installed\n---\n\n")




