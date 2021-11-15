#!/usr/bin/env python

######### define environment ##########

# module imports
import sys
import os
import traceback
from Bio import SeqIO
from Bio.Seq import Seq
import random

# get the cwd were all the scripts are 
test_dir = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, test_dir)
scripts_dir = "%s/../../scripts"%test_dir; sys.path.insert(0, scripts_dir)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])
EnvName = EnvDir.split("/")[-1]
AnacondaDir = "/".join(sys.executable.split("/")[0:-4])

# define script names
varcall_cnv_pipeline = "%s/varcall_cnv_pipeline.py"%scripts_dir

# import functions
import sv_functions as fun
import multiprocessing as multiproc

# change the verbosity
fun.printing_verbose_mode = True

# define the threads
threads = multiproc.cpu_count()


######################################

def test_get_repeat_maskerDF(test_genome, replace=False):

    """Tests the generation of repeats"""

    try: 

        fun.print_with_runtime("trying to generate repeats for %s"%test_genome)

        # define the ref genome
        df_repeats, repeat_masker_outfile_default = fun.get_repeat_maskerDF(test_genome, threads=threads, replace=replace)

        # test
        expected_fields = {'perc_ins', 'perc_del', 'type', 'end_repeat', 'perc_div', 'position_inRepeat_begin', 'repeat', 'left_repeat', 'IDrepeat', 'strand', 'left_positionINrepeat', 'SW_score', 'chromosome', 'begin_repeat', 'position_inRepeat_end'}

        if set(list(df_repeats.keys()))!=expected_fields: raise ValueError("something went wrong with the repeats generation")

        fun.print_with_runtime("repeats were generated correctly for %s"%test_genome)

    except Exception as err:


        fun.print_with_runtime("\n\n---\nWARNING: The running of RepeatModeller and/or RepeatMasker did not work. This is expected to happen in old linux systems with a GLIBC version <2.12. You may provide the already calculated repeats through the command --previous_repeats_table. This is returned by the function of get_repeat_maskerDF from sv_functions.py . \n---\n\n")

        fun.print_with_runtime("---\nThis is the error:")

        traceback.print_tb(err.__traceback__)
        fun.print_with_runtime(err)
        fun.print_with_runtime("---\n")

def test_conda_env_generation(outdir, replace=False):

    """This function exports the current perSVade_env to a .yml file, and generates a conda file"""

    # define the file that indicates that the enviornment is correct
    correct_env_file = "%s/correct_env.txt"%outdir

    # define a test_env_name
    test_env_name = "%s_test"%EnvName

    if fun.file_is_empty(correct_env_file) or replace is True:

        # remove previous env
        fun.print_with_runtime("removing previous env")
        try: fun.run_cmd("conda remove -y -n %s --all"%test_env_name)
        except: fun.print_with_runtime("%s does not exist"%test_env_name)

        # export file
        fun.print_with_runtime("creating %s yml"%EnvName)
        yml_file = "%s/%s.yml"%(outdir, test_env_name)
        fun.run_cmd("conda env export --no-builds --from-history -n %s --file %s"%(EnvName, yml_file))

        # create environment
        fun.print_with_runtime("re-generating as %s"%test_env_name)
        fun.run_cmd("conda env create --file %s --name %s"%(yml_file, test_env_name))
        
        # test that the activation works
        fun.print_with_runtime("activating %s"%test_env_name)
        cmd_activation = "source %s/etc/profile.d/conda.sh && conda activate %s && python -c 'import sys; sys.path.insert(0, \"%s\"); import sv_functions as fun'"%(AnacondaDir, test_env_name, fun.get_fullpath(scripts_dir))
        fun.run_cmd(cmd_activation)

        # remove file
        fun.print_with_runtime("removing envs")
        fun.remove_file(yml_file)

        # remove env
        fun.run_cmd("conda remove -y -n %s --all"%test_env_name)

        # create file stating that the env is correct
        open(correct_env_file, "w").write("env is correct")

    fun.print_with_runtime("%s can be correctly regenerated"%EnvName)

def test_read_simulation_and_get_reads(genome, window_l=2000, npairs=50000, read_length=150, median_insert_size=250, median_insert_size_sd=50, threads=threads, replace=False):

    """ 
    Takes a genome and simulates reads for it, saving them under <genome>_simulating_reads 
    """

    # define the outdir
    outdir = "%s_simulating_reads"%genome; 
    outdir_reads = "%s/getting_reads"%outdir; 
    
    # remove the outdirs if replace is True
    if replace is True: 
        fun.delete_folder(outdir)
        fun.delete_folder(outdir_reads)

    # make folders 
    fun.make_folder(outdir)
    fun.make_folder(outdir_reads)

    # define the expected reads
    reads1 = "%s/all_reads1.correct.fq.gz"%outdir_reads
    reads2 = "%s/all_reads2.correct.fq.gz"%outdir_reads

    if any([fun.file_is_empty(f) for f in [reads1, reads2]]):

        # run index the genome
        fun.run_cmd("%s faidx %s"%(fun.samtools, genome))

        # get the windows df
        windows_bed = "%s/windows_file.bed"%outdir
        fun.run_cmd("%s makewindows -g %s.fai -w %i > %s"%(fun.bedtools, genome, window_l, windows_bed))
        df_windows = fun.pd.read_csv(windows_bed, sep="\t", header=-1, names=["chromosome", "start", "end"])
        df_windows["predicted_relative_coverage"] = fun.random.sample(list(fun.np.linspace(0.5, 2, 10000)), len(df_windows))

        # simulate reads
        fun.simulate_readPairs_per_window(df_windows, genome, npairs, outdir_reads, read_length, median_insert_size, median_insert_size_sd, replace=False, threads=threads) 

    fun.print_with_runtime("read simulation works well")
    return reads1, reads2

def test_bwa_mem_and_get_bam(r1, r2, ref_genome, replace=False):

    """Runs bwa mem on the reads and returns the sorted bam with marked duplicates"""

    # define the outdir
    outdir = "%s/aligning_reads_against_%s"%(fun.get_dir(r1), fun.get_file(ref_genome))

    # if replace is True, delete the outdir
    if replace is True: fun.delete_folder(outdir)

    # make de outdir
    fun.make_folder(outdir)

    # define the inputs of bam
    bamfile = "%s/aligned_reads.bam"%outdir
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam
    name_sample = "test_sample"

    # run
    fun.print_with_runtime("aligning reads")
    fun.run_bwa_mem(r1, r2, ref_genome, outdir, bamfile, sorted_bam, index_bam, name_sample, threads=threads, replace=False, MarkDuplicates=True)

    return sorted_bam


def test_processing_varcalling(smallVars_input_outdir, reference_genome, outdir, sorted_bam, replace=False, threads=threads):

    """This function takes a varcall file were all the variant calling has been performed and checks that the processing of vcfs works in varcall_cnv_pipeline. sorted_bam is just a fake sorted bam not to repeat the pipeline running"""

    # get full paths
    outdir = fun.get_fullpath(outdir)
    smallVars_input_outdir = fun.get_fullpath(smallVars_input_outdir)
    reference_genome = fun.get_fullpath(reference_genome)

    # cp the files under outdir
    fun.make_folder(outdir)

    target_smallVars_input_outdir = "%s/smallVars_input_outdir"%outdir
    target_smallVars_input_outdir_tmp = "%s.tmp"%target_smallVars_input_outdir
    if not os.path.isdir(target_smallVars_input_outdir) or replace is True:
        fun.run_cmd("cp -r %s %s "%(smallVars_input_outdir, target_smallVars_input_outdir_tmp))
        os.rename(target_smallVars_input_outdir_tmp, target_smallVars_input_outdir)

    # final file 
    final_file = "%s/variants_atLeast3PASS_ploidy2.vcf"%target_smallVars_input_outdir

    if fun.file_is_empty(final_file) or replace is True:

        cmd = "%s -r %s -o %s -p 2 -sbam %s -caller all -c 5 -mchr no_mitochondria -mcode 3 -gcode 1 --repeats_table %s.repeats.tab --remove_smallVarsCNV_nonEssentialFiles -thr %i --skip_cnv_analysis"%(varcall_cnv_pipeline, reference_genome, target_smallVars_input_outdir, sorted_bam, reference_genome, threads) 

        fun.run_cmd(cmd)

    fun.print_with_runtime("you can run successfully the variant processing")


def test_smallVarCall_CNV_running(sorted_bam, outdir, ref_genome, gff, threads=threads, mitochondrial_chromosome="mito_C_glabrata_CBS138", replace=False):

    """Takes a sorted bam (shuld have some mutations) and runs the variant calling pipeline on it"""

    # if replace is True, remove the outdir
    if replace is True: fun.delete_folder(outdir)

    # make the outdir
    fun.make_folder(outdir)

    # get the repeats
    repeats_table = "%s.repeats.tab"%ref_genome
    fun.write_repeats_table_file(repeats_table)

    for pooled_seq in [False]: # this may be also [False, True] to test pooled seq

        outdir_varCall = "%s/varcall_pooledSeq_%s"%(outdir, str(pooled_seq))

        # define the final file
        final_file = "%s/variant_annotation_ploidy2.tab"%outdir_varCall

        if fun.file_is_empty(final_file) or replace is True:
            fun.print_with_runtime("running on pooled_seq=%s. If pooled_seq is True, this may take a bit because a lot of variants will be considered"%pooled_seq)

            # define the cmd
            cmd = "%s -r %s -o %s -p 2 -sbam %s -caller all -c 5 -mchr %s -mcode 3 -gcode 1 --repeats_table %s --remove_smallVarsCNV_nonEssentialFiles -gff %s -thr %i"%(varcall_cnv_pipeline, ref_genome, outdir_varCall, sorted_bam, mitochondrial_chromosome, repeats_table, gff, threads) 

            # add pooled seq
            if pooled_seq is True: cmd += " --pooled_sequencing"

            fun.run_cmd(cmd)

    fun.print_with_runtime("small variant calling and CNV of genes works")

def test_SRAdb_query_downloading_and_readTrimming(outdir, reference_genome, target_taxID, replace=False, threads=threads):

    """This function runs get_close_shortReads_table_close_to_taxID for the MERS coronavirus and taking the lowest coverage reads. This tests that sra tools, entrez tools, trimmomatic and fastqc work well.

    5476 is C. albicans
    1335626 is MERS
    """

    # make the outdir
    fun.make_folder(outdir)


    # set ploidy to 1
    ploidy=1

    try:
        # run with 'get_lowest_coverage_possible=True', which will take the lowest coverage datasets
        close_shortReads_table = fun.get_close_shortReads_table_close_to_taxID(target_taxID, reference_genome, outdir, ploidy, n_close_samples=1, nruns_per_sample=1, replace=replace, threads=threads, min_fraction_reads_mapped=0.0, coverage_subset_reads=0.1, min_coverage=5, job_array_mode="local", StopAfter_sampleIndexingFromSRA=False, StopAfterPrefecth_of_reads=False, get_lowest_coverage_possible=True)

        # check
        df_close_shortReads_table = fun.pd.read_csv(close_shortReads_table, sep="\t")

        if set(df_close_shortReads_table.keys())!={'short_reads2', 'short_reads1', 'runID', 'sampleID'} or len(df_close_shortReads_table)!=1: raise ValueError("The close_shortReads_table %s was not created as expected"%close_shortReads_table)

        fun.print_with_runtime("The system to query the SRA database, dowload and trim reads works")

    except Exception as err:


        fun.print_with_runtime("\n\n---\nWARNING: The connection to SRA did not work. This means that the automated obtention of reads of close species for benchmarking (involving the arguments --target_taxID, --n_close_samples, --nruns_per_sample or --goldenSet_dir) may fail. You can download the reads on your own and provide them with --close_shortReads_table. This can be also due to network problems at this moment. \n---\n\n")

        fun.print_with_runtime("---\nThis is the error:")

        traceback.print_tb(err.__traceback__)
        fun.print_with_runtime(err)
        fun.print_with_runtime("---\n")

def test_rearranging_genome_random(ref_genome, replace=False, threads=threads, mitochondrial_chromosome="mito_C_glabrata_CBS138", nvars=5):


    """This function takes a reference genome and simulates random variation on it, returning the rearranged genome in fasta format"""


    # define the outdir
    outdir = "%s.testing_rearranged_genome_generation"%(ref_genome); fun.make_folder(outdir)

    sim_svtype_to_svfile, rearranged_genome = fun.simulate_SVs_in_genome(ref_genome, mitochondrial_chromosome, outdir, nvars=nvars, bedpe_breakpoints=None, replace=replace)

    fun.print_with_runtime("The generation of a genome with randomly-inserted SVs works")

    return rearranged_genome


def test_gridss_clove_pipeline(sorted_bam, reference_genome, outdir, threads=threads, replace=False):

    """Tests that the default gridss and clove pipeline can be obtained"""

    fun.make_folder(outdir)

    # define the median coverage (to be recalculated)
    median_coverage = -1

    fun.print_with_runtime("running gridss+clove pipeline on %i threads"%threads)
    
    SV_dict, df_gridss = fun.run_gridssClove_given_filters(sorted_bam, reference_genome, outdir, median_coverage, replace=replace, threads=threads, gridss_blacklisted_regions="", gridss_VCFoutput="", gridss_maxcoverage=50000, median_insert_size=250, median_insert_size_sd=25, gridss_filters_dict=fun.default_filtersDict_gridss, run_in_parallel=True, max_rel_coverage_to_consider_del=0.2, min_rel_coverage_to_consider_dup=1.8, replace_FromGridssRun=replace)

    fun.print_with_runtime("you could run the gridss + clove pipeline succesfully")





def test_parameter_optimisation_perSVade(sorted_bam, reference_genome, outdir, threads=threads, replace=False):

    """This pipeline will test the parameter optimisation features of perSVade into outdir. It is expected to work for C.glabrata"""

    if fun.file_is_empty("%s/perSVade_finished_file.txt"%outdir):

        #cmd = "%s -r %s -thr %i -o %s -sbam %s --nvars 5 --simulation_ploidies haploid --range_filtering_benchmark theoretically_meaningful -mchr mito_C_glabrata_CBS138 --min_chromosome_len 100 --nsimulations 1 --skip_repeat_analysis --skip_CNV_calling"%(fun.perSVade_py, reference_genome, threads, outdir, sorted_bam)

        cmd = "%s -r %s -thr %i -o %s -sbam %s --nvars 5 --simulation_ploidies haploid --range_filtering_benchmark theoretically_meaningful -mchr mito_C_glabrata_CBS138 --min_chromosome_len 100 --nsimulations 1 --skip_repeat_analysis --window_size_CNVcalling 200 --cnv_calling_algs HMMcopy,AneuFinder"%(fun.perSVade_py, reference_genome, threads, outdir, sorted_bam)

        if fun.printing_verbose_mode is True: cmd += " --verbose"

        fun.run_cmd(cmd)

    fun.print_with_runtime("parameter optimisation worked successfully")


def test_realSVgeneration(reads_dir, outdir, repeats, reference_genome, relaxed_parms, replace=False, threads=threads):

    """This function takes a reads dir (from C. glabrata), generates a table that is passed to perSVade to generate the real variants. It returns the integrated bedpe"""


    # define the real bedpe breakpoints
    real_bedpe_breakpoints = "%s/findingRealSVs_providedCloseReads/integrated_breakpoints.bedpe"%outdir

    if fun.file_is_empty(real_bedpe_breakpoints) or replace is True:

        # define the close_shortReads_table
        reads_dict = {"first100k": {"sampleID":"first100k", "runID":"first100k", "short_reads1":"sampled_readsR1_first100k.fq.gz", "short_reads2":"sampled_readsR2_first100k.fq.gz"},

                      "last100k": {"sampleID":"last100k", "runID":"last100k", "short_reads1":"sampled_readsR1_last100k.fq.gz", "short_reads2":"sampled_readsR2_last100k.fq.gz"}}

        reads_df = fun.pd.DataFrame(reads_dict).transpose()
        for f in ["short_reads1", "short_reads2"]: reads_df[f] = reads_dir + "/" + reads_df[f]
        fun.make_folder(outdir)
        close_shortReads_table = "%s/subsampled_shortReads.tab"%outdir
        reads_df.to_csv(close_shortReads_table, sep="\t", index=False, header=True)

        # run persvade
        cmd = "%s -r %s -thr %i -o %s -f1 skip -f2 skip --nvars 500 -mchr mito_C_glabrata_CBS138 --min_chromosome_len 10000 --StopAfter_obtentionOFcloseSVs --close_shortReads_table %s --previous_repeats_table %s --parameters_json_file %s --skip_CNV_calling"%(fun.perSVade_py, reference_genome, threads, outdir, close_shortReads_table, repeats, relaxed_parms)

        if fun.printing_verbose_mode is True: cmd += " --verbose"
        fun.run_cmd(cmd)

    fun.print_with_runtime("generation of real SVs worked. This also validates that the parameter json loading did work")

    return real_bedpe_breakpoints


def test_parameter_optimisation_perSVade_real(reads_dir, outdir, repeats, reference_genome, real_bedpe_breakpoints, replace=False, threads=threads):

    """This function runs perSVade optimisation based on real SVs defined by real_bedpe_breakpoints."""

    fun.print_with_runtime("testing parameter optimisation based on real SVs")

    if fun.file_is_empty("%s/perSVade_finished_file.txt"%outdir) or replace is True:

        # define the reads
        reads1 = "%s/sampled_readsR1_first100k.fq.gz"%reads_dir
        reads2 = "%s/sampled_readsR2_first100k.fq.gz"%reads_dir

        # run persvade
        cmd = "%s -r %s -thr %i -o %s -f1 %s -f2 %s --nvars 10 -mchr mito_C_glabrata_CBS138 --min_chromosome_len 10000 --real_bedpe_breakpoints %s --previous_repeats_table %s --simulation_ploidies haploid --range_filtering_benchmark large --nsimulations 1 -p 1 --skip_repeat_analysis"%(fun.perSVade_py, reference_genome, threads, outdir, reads1, reads2, real_bedpe_breakpoints, repeats)


        if fun.printing_verbose_mode is True: cmd += " --verbose"
        fun.run_cmd(cmd)

    fun.print_with_runtime("the parameter optimisation based on real SVs worked")

def test_greasy():

    """This function tests whether greasy can be found in the path"""

    try:

        fun.run_cmd("module load greasy")
        fun.print_with_runtime("greasy module can be loaded")

        fun.run_cmd("module load greasy && which greasy")
        fun.print_with_runtime("greasy is in the path")

        fun.run_cmd("which sbatch")
        fun.print_with_runtime("sbatch is in the path")

        fun.print_with_runtime("greasy can be used for running parallel perSVade jobs")

    except Exception as err:

        fun.print_with_runtime("\n\n---\nWARNING: greasy is not installed properly in your system. This means that setting '--job_array_mode greasy' will fail. You can set '--job_array_mode local' and run jobs sequentially and not in parallel. '--job_array_mode greasy' will only work on machines that use SLURM for managing jobs and greasy installed (and callable with a command like 'module load greasy && greasy <jobs_file>', and '<jobs_file>' is a file where each line corresponds to a command to be executed in a sepparate SLURM job). If you have a cluster environment without greasy, you can run in '--job_array_mode greasy' and manage the jobs on your own. \n---\n\n")


        fun.print_with_runtime("---\nThis is the error:")

        traceback.print_tb(err.__traceback__)
        fun.print_with_runtime(err)
        fun.print_with_runtime("---\n")

def test_picard_env(testing_inputs_dir, testing_outputs_dir):

    """Makes sure that the picard environment was properly generated"""

    # define inputs
    sorted_bam = "%s/readsWithSVs_against_reducedGenome.sorted.bam"%testing_outputs_dir 
    fun.soft_link_files("%s/readsWithSVs_against_reducedGenome.sorted.bam"%testing_inputs_dir , sorted_bam)

    ref_genome = "%s/reduced_genome.fasta"%testing_outputs_dir
    fun.soft_link_files("%s/reduced_genome.fasta"%testing_inputs_dir , ref_genome)

    # get the insert size distribution, which relies on picard
    fun.get_insert_size_distribution(sorted_bam, replace=False, threads=threads)

    # get the sequence dict
    fun.create_sequence_dict(ref_genome, replace=False)

    fun.print_with_runtime("picard_env was properly generated")


def get_mutated_seqrecord(seq, mutation_rate=0.02):

    """Gets a seqrecord and mutates it"""

    # keep seq
    import copy as cp
    import numpy as np
    modified_seq = cp.deepcopy(seq)

    # define a function that takes a base and returns a random one in some cases
    all_bases = ["A", "C", "G", "T"]
    def get_random_base(input_base):
        if random.random()<=mutation_rate: return random.choice(all_bases)
        else: return input_base.upper()

    # get mutated seq
    modified_seq.seq = Seq("".join(list(map(get_random_base, str(seq.seq)))))

    # calculate the number of modified bases
    mutated_bases = sum(np.array(list(modified_seq.seq))!=np.array(list(seq.seq)))

    # log
    print("There are %i/%i mutated bases in %s"%(mutated_bases, len(seq.seq), seq.id))

    return modified_seq

def get_reads_with_small_variants(reference_genome, outdir, threads):

    """Takes a reference genome and generates some reads with mutations"""  

    # get outdir
    #fun.delete_folder(outdir)
    fun.make_folder(outdir)

    # get mutated reference genome
    mut_reference_genome = "%s/mutated_ref_genome.fasta"%outdir
    if fun.file_is_empty(mut_reference_genome):
        print("getting mutated genome")

        # get mutated seqs
        mutated_seqs = [get_mutated_seqrecord(seq) for seq in  SeqIO.parse(reference_genome, "fasta")]

        # write them
        mut_reference_genome_tmp  = "%s.tmp"%mut_reference_genome
        SeqIO.write(mutated_seqs, mut_reference_genome_tmp, "fasta")
        os.rename(mut_reference_genome_tmp, mut_reference_genome)

    # get the mutated reads
    print("getting the mutated reads")
    genome_len = sum(fun.get_chr_to_len(mut_reference_genome).values())
    read_len = 150
    coverage = 30
    npairs = (coverage*genome_len)/read_len

    mut_reads1, mut_reads2 = fun.simulate_testing_reads_on_genome(mut_reference_genome, window_l=4000, npairs=npairs, read_length=read_len, median_insert_size=150, median_insert_size_sd=15, threads=threads, replace=False)

    return mut_reads1, mut_reads2