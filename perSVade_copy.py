#!/usr/bin/env python

# This is the perSVade pipeline main script, which shoul dbe run on the perSVade conda environment


##### DEFINE ENVIRONMENT #######


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

mandatory_args.add_argument("-r", "--ref", dest="ref", required=True, help="Reference genome. Has to end with .fasta.")
mandatory_args.add_argument("-mchr", "--mitochondrial_chromosome", dest="mitochondrial_chromosome", required=True, type=str, help="The name of the mitochondrial chromosome. This is important if you have mitochondrial proteins for which to annotate the impact of nonsynonymous variants, as the mitochondrial genetic code is different. This should be the same as in the gff. If there is no mitochondria just put 'no_mitochondria'. If there is more than one mitochindrial scaffold, provide them as comma-sepparated IDs.")

#################################


###### inputs #######

seq_inputs = parser.add_argument_group("SEQUENCE INPUTS")

seq_inputs.add_argument("-f1", "--fastq1", dest="fastq1", default=None, help="fastq_1 file. Option required to obtain bam files. It can be 'skip', which will tell the pipeline to not use any fastq or bam input.")
seq_inputs.add_argument("-f2", "--fastq2", dest="fastq2", default=None, help="fastq_2 file. Option required to obtain bam files. It can be 'skip', which will tell the pipeline to not use any fastq or bam input.")


seq_inputs.add_argument("-sbam", "--sortedbam", dest="sortedbam", default=None, help="The path to the sorted bam file, which should have a bam.bai file in the same dir. For example, if your bam file is called 'aligned_reads.bam', there should be an 'aligned_reads.bam.bai' as well. This is mutually exclusive with providing reads. By default, it is assumed that this bam has marked duplicates.")


#####################

########### resource allocation ##############

resources_args = parser.add_argument_group("RESOURCES")

resources_args.add_argument("--tmpdir", dest="tmpdir", default=None, help="A full path to a directory where to write intermediate files. This is useful if you are running on a cluster that has some directories that have higher writing speed than others.")


general_optional_args = parser.add_argument_group("GENERAL OPTIONAL ARGUMENTS")

general_optional_args.add_argument("-p", "--ploidy", dest="ploidy", default=1, type=int, help="Ploidy, can be 1 or 2")

general_optional_args.add_argument("-gff", "--gff-file", dest="gff", default=None, help="path to the GFF3 annotation of the reference genome. Make sure that the IDs are completely unique for each 'gene' tag. This is necessary for both the CNV analysis (it will look at genes there) and the annotation of the variants.")

general_optional_args.add_argument("-mcode", "--mitochondrial_code", dest="mitochondrial_code", default=3, type=int, help="The code of the NCBI mitochondrial genetic code. For yeasts it is 3. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. The information of this website may be wrong, so you may want to double check with the literature.")

general_optional_args.add_argument("-gcode", "--gDNA_code", dest="gDNA_code", default=1, type=int, help="The code of the NCBI gDNA genetic code. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . For C. albicans it is 12. The information of this website may be wrong, so you may want to double check with the literature.")


general_optional_args.add_argument("--QC_and_trimming_reads", dest="QC_and_trimming_reads", action="store_true", default=False, help="Will run fastq and trimmomatic of reads, and use the trimmed reads for downstream analysis. This option will generate files under the same dir as f1 and f2, so be aware of it. This is a rather dangerous option, sinc you may want to do the quality control of the reads before running perSVade.")

general_optional_args.add_argument("--min_chromosome_len", dest="min_chromosome_len", default=100000, type=int, help="The minimum length to consider chromosomes from the provided fasta for calculating the window length. Any chromosomes that shorter than the window length will not be considered in the random SV simulations.")


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





########################################
##### GENERAL PROCESSING OF INPUTS #####
########################################





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

#### define misc args ####


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




# define the min_CNVsize_coverageBased
fun.min_CNVsize_coverageBased = opt.min_CNVsize_coverageBased

# change the default parameters if specified
if opt.parameters_json_file is not None:

    gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup = fun.get_parameters_from_json(opt.parameters_json_file)

    fun.default_filtersDict_gridss = gridss_filters_dict
    fun.default_gridss_blacklisted_regions = gridss_blacklisted_regions
    fun.default_gridss_maxcoverage = gridss_maxcoverage
    fun.default_max_rel_coverage_to_consider_del = max_rel_coverage_to_consider_del
    fun.default_min_rel_coverage_to_consider_dup = min_rel_coverage_to_consider_dup

# test whether the gff is correct
if opt.gff is not None: fun.check_that_gff_is_correct(opt.gff, opt.ref, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code, opt.threads, opt.replace)

# check that the tmpdir exists
if opt.tmpdir is not None:
    if not os.path.isdir(opt.tmpdir): raise ValueError("The folder that you specified with --tmpdir does not exist")

# get the repeats table
if opt.skip_repeat_analysis is False:

else:

    print("skipping the repeats analysis")
    fun.write_repeats_table_file(repeats_table_file)



####### GENERATE real_bedpe_breakpoints AROUND HOMOLOGOUS REGIONS #########

if opt.simulate_SVs_around_HomologousRegions is True:
    print("simulating around Homologous regions")


    # get a file that contains the blastn of some regions of the genome
    if opt.simulate_SVs_around_HomologousRegions_previousBlastnFile is None: blastn_file = fun.get_blastn_regions_genome_against_itself(opt.ref, opt.simulate_SVs_around_HomologousRegions_maxEvalue, opt.simulate_SVs_around_HomologousRegions_queryWindowSize, opt.replace, opt.threads)

    else: blastn_file = opt.simulate_SVs_around_HomologousRegions_previousBlastnFile

    # define the bedpe breakpoints around the homologous regions
    bedpe_breakpoints = "%s/breakpoints_aroundHomRegions_wsize=%ibp_maxEval=%s_minQcovS=%i.bedpe"%(opt.outdir, opt.simulate_SVs_around_HomologousRegions_queryWindowSize, opt.simulate_SVs_around_HomologousRegions_maxEvalue, opt.simulate_SVs_around_HomologousRegions_minPctOverlap)

    opt.real_bedpe_breakpoints = fun.get_bedpe_breakpoints_around_homologousRegions(blastn_file, bedpe_breakpoints, replace=opt.replace, threads=opt.threads, max_eval=opt.simulate_SVs_around_HomologousRegions_maxEvalue, query_window_size=opt.simulate_SVs_around_HomologousRegions_queryWindowSize, min_qcovs=opt.simulate_SVs_around_HomologousRegions_minPctOverlap, max_n_hits=(opt.nvars*5*2500))

############################################################################



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


        # get the linked bam files
        fun.soft_link_files(fun.get_fullpath(opt.sortedbam), sorted_bam)
        fun.soft_link_files(fun.get_fullpath(opt.sortedbam)+".bai", sorted_bam+".bai")
        
    ###################################



#####################################
#####################################
#####################################



#####################################
##### STRUCTURAL VARIATION ##########
#####################################


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
#####################################

# keep the simulation files
if opt.keep_simulation_files is True: fun.keep_simulation_files_for_perSVade_outdir(opt.outdir, replace=opt.replace, n_simulated_genomes=opt.nsimulations, simulation_ploidies=simulation_ploidies)

# at the end you want to clean the outdir to keep only the essential files
if opt.skip_cleaning_outdir is False: fun.clean_perSVade_outdir(opt.outdir)
