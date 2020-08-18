#!/usr/bin/env python

######### define environment ##########

# module imports
import sys
import os

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
fun.printing_verbose_mode = False

# define the threads
threads = multiproc.cpu_count()


######################################

def test_get_repeat_maskerDF(test_genome, replace=False):

    """Tests the generation of repeats"""

    # define the ref genome
    df_repeats, repeat_masker_outfile_default = fun.get_repeat_maskerDF(test_genome, threads=threads, replace=replace)

    # test
    expected_fields = {'perc_ins', 'perc_del', 'type', 'end_repeat', 'perc_div', 'position_inRepeat_begin', 'repeat', 'left_repeat', 'IDrepeat', 'strand', 'left_positionINrepeat', 'SW_score', 'chromosome', 'begin_repeat', 'position_inRepeat_end'}

    if set(list(df_repeats.keys()))!=expected_fields: raise ValueError("something went wrong with the repeats generation")

    print("repeats were generated correctly")

def test_conda_env_generation(outdir, replace=False):

    """This function exports the current perSVade_env to a .yml file, and generates a conda file"""

    # define the file that indicates that the enviornment is correct
    correct_env_file = "%s/correct_env.txt"%outdir

    # define a test_env_name
    test_env_name = "%s_test"%EnvName

    if fun.file_is_empty(correct_env_file) or replace is True:

        # remove previous env
        print("removing previous env")
        try: fun.run_cmd("conda remove -y -n %s --all"%test_env_name)
        except: print("%s does not exist"%test_env_name)

        # export file
        print("creating %s yml"%EnvName)
        yml_file = "%s/%s.yml"%(outdir, test_env_name)
        fun.run_cmd("conda env export --no-builds --from-history -n %s --file %s"%(EnvName, yml_file))

        # create environment
        print("re-generating as %s"%test_env_name)
        fun.run_cmd("conda env create --file %s --name %s"%(yml_file, test_env_name))
        
        # test that the activation works
        print("activating %s"%test_env_name)
        cmd_activation = "source %s/etc/profile.d/conda.sh && conda activate %s && python -c 'import sys; sys.path.insert(0, \"%s\"); import sv_functions as fun'"%(AnacondaDir, test_env_name, fun.get_fullpath(scripts_dir))
        fun.run_cmd(cmd_activation)

        # remove file
        print("removing envs")
        fun.remove_file(yml_file)

        # remove env
        fun.run_cmd("conda remove -y -n %s --all"%test_env_name)

        # create file stating that the env is correct
        open(correct_env_file, "w").write("env is correct")

    print("%s can be correctly regenerated"%EnvName)

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

    print("read simulation works well")
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
    print("aligning reads")
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

    print("you can run successfully the variant processing")


def test_smallVarCall_CNV_running(sorted_bam, outdir, ref_genome, gff, threads=threads, mitochondrial_chromosome="mito_C_glabrata_CBS138", replace=False):

    """Takes a sorted bam (shuld have some mutations) and runs the variant calling pipeline on it"""

    # if replace is True, remove the outdir
    if replace is True: fun.delete_folder(outdir)

    # make the outdir
    fun.make_folder(outdir)

    # get the repeats
    repeats_table = fun.get_repeat_maskerDF(ref_genome, threads=threads, replace=False)[1]

    for pooled_seq in [False]: # this may be also [False, True] to test pooled seq

        outdir_varCall = "%s/varcall_pooledSeq_%s"%(outdir, str(pooled_seq))

        # define the final file
        final_file = "%s/variant_annotation_ploidy2.tab"%outdir_varCall

        if fun.file_is_empty(final_file) or replace is True:
            print("running on pooled_seq=%s. This may take a bit because a lot of variants will be considered"%pooled_seq)

            # define the cmd
            cmd = "%s -r %s -o %s -p 2 -sbam %s -caller all -c 5 -mchr %s -mcode 3 -gcode 1 --repeats_table %s --remove_smallVarsCNV_nonEssentialFiles -gff %s -thr %i"%(varcall_cnv_pipeline, ref_genome, outdir_varCall, sorted_bam, mitochondrial_chromosome, repeats_table, gff, threads) 

            # add pooled seq
            if pooled_seq is True: cmd += " --pooled_sequencing"

            fun.run_cmd(cmd)

    print("small variant calling and CNV of genes works")


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
        close_shortReads_table = fun.get_close_shortReads_table_close_to_taxID(target_taxID, reference_genome, outdir, ploidy, n_close_samples=2, nruns_per_sample=1, replace=replace, threads=threads, min_fraction_reads_mapped=0.0, coverage_subset_reads=0.1, min_coverage=5, job_array_mode="local", StopAfter_sampleIndexingFromSRA=False, queue_jobs="debug", max_ncores_queue=768, time_read_obtention="02:00:00", StopAfterPrefecth_of_reads=False, get_lowest_coverage_possible=True)


        # check
        df_close_shortReads_table = fun.pd.read_csv(close_shortReads_table, sep="\t")

        if set(df_close_shortReads_table.keys())!={'short_reads2', 'short_reads1', 'runID', 'sampleID'} or len(df_close_shortReads_table)!=2: raise ValueError("The close_shortReads_table %s was not created as expected"%close_shortReads_table)

        print("The system to query the SRA database, dowload and trim reads works")

    except:

        print("\n\n---\nWARNING: The connection to SRA did not work. This means that the automated obtention of reads of close species for benchmarking (involving the arguments --target_taxID, --n_close_samples, --nruns_per_sample or --goldenSet_dir) may fail. You can download the reads on your own and provide them with --close_shortReads_table. This can be also due to network problems at this moment. \n---\n\n")
        
def test_rearranging_genome_random(ref_genome, replace=False, threads=threads, mitochondrial_chromosome="mito_C_glabrata_CBS138", nvars=5):


    """This function takes a reference genome and simulates random variation on it, returning the rearranged genome in fasta format"""


    # define the outdir
    outdir = "%s.testing_rearranged_genome_generation"%(ref_genome); fun.make_folder(outdir)

    sim_svtype_to_svfile, rearranged_genome = fun.rearrange_genomes_simulateSV(ref_genome, outdir, replace=replace, nvars=nvars, mitochondrial_chromosome=mitochondrial_chromosome)


    print("The generation of a genome with randomly-inserted SVs works")

    return rearranged_genome


def test_gridss_clove_pipeline(sorted_bam, reference_genome, outdir, threads=threads, replace=False):

    """Tests that the default gridss and clove pipeline can be obtained"""

    fun.make_folder(outdir)

    # define the median coverage (to be recalculated)
    median_coverage = -1

    print("running gridss+clove pipeline on %i threads"%threads)
    
    SV_dict, df_gridss = fun.run_gridssClove_given_filters(sorted_bam, reference_genome, outdir, median_coverage, replace=replace, threads=threads, gridss_blacklisted_regions="", gridss_VCFoutput="", gridss_maxcoverage=50000, median_insert_size=250, median_insert_size_sd=25, gridss_filters_dict=fun.default_filtersDict_gridss, run_in_parallel=True, max_rel_coverage_to_consider_del=0.2, min_rel_coverage_to_consider_dup=1.8, replace_FromGridssRun=replace)

    print("you could run the gridss + clove pipeline succesfully")


def test_parameter_optimisation_perSVade(sorted_bam, reference_genome, outdir, threads=threads, replace=False):

    """This pipeline will test the parameter optimisation features of perSVade into outdir. It is expected to work for C.glabrata"""

    cmd = "%s -r %s -thr %i -o %s -sbam %s --nvars 5 --simulation_ploidies haploid,diploid_hetero --range_filtering_benchmark theoretically_meaningful -mchr mito_C_glabrata_CBS138"%(fun.perSVade_py, reference_genome, threads, outdir, sorted_bam)

    fun.run_cmd(cmd)

    print("parameter optimisation worked successfully")

def test_greasy():

    """This function tests whether greasy can be found in the path"""

    try:

        fun.run_cmd("module load greasy")
        print("greasy module can be loaded")

        fun.run_cmd("module load greasy && which greasy")
        print("greasy is in the path")

        fun.run_cmd("which sbatch")
        print("sbatch is in the path")

        print("greasy can be used for running parallel perSVade jobs")

    except:

        print("\n\n---\nWARNING: greasy is not installed properly in your system. This means that setting '--job_array_mode greasy' will fail. You can set '--job_array_mode local' and run jobs sequentially and not in parallel. '--job_array_mode greasy' will only work on machines that use SLURM for managing jobs and greasy installed (and callable with a command like 'module load greasy && greasy <jobs_file>', and '<jobs_file>' is a file where each line corresponds to a command to be executed in a sepparate SLURM job) \n---\n\n")





