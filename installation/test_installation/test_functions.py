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

######################################

def test_get_repeat_maskerDF(test_genome):

    """Tests the generation of repeats"""

    # define the ref genome
    df_repeats, repeat_masker_outfile_default = fun.get_repeat_maskerDF(test_genome, threads=4, replace=False)

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

def test_read_simulation_and_get_reads(genome, window_l=2000, npairs=50000, read_length=150, median_insert_size=250, median_insert_size_sd=50, threads=4):

    """ 
    Takes a genome and simulates reads for it, saving them under <genome>_simulating_reads 
    """

    # define the outdir
    outdir = "%s_simulating_reads"%genome; fun.make_folder(outdir)
    outdir_reads = "%s/getting_reads"%outdir; fun.make_folder(outdir_reads)

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
        fun.simulate_readPairs_per_window(df_windows, genome, npairs, outdir_reads, read_length, median_insert_size, median_insert_size_sd, replace=False, threads=4) 

    print("read simulation works well")
    return reads1, reads2

def test_bwa_mem_and_get_bam(r1, r2, ref_genome):

    """Runs bwa mem on the reads and returns the sorted bam with marked duplicates"""

    # define the outdir
    outdir = "%s/aligning_reads_against_%s"%(fun.get_dir(r1), fun.get_file(ref_genome))
    fun.make_folder(outdir)

    # define the inputs of bam
    bamfile = "%s/aligned_reads.bam"%outdir
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam
    name_sample = "test_sample"

    # run
    print("aligning reads")
    fun.run_bwa_mem(r1, r2, ref_genome, outdir, bamfile, sorted_bam, index_bam, name_sample, threads=4, replace=False, MarkDuplicates=True)

    return sorted_bam


def test_processing_varcalling(smallVars_input_outdir, reference_genome, outdir, sorted_bam, replace=False, threads=4):

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

    target_reference_genome = "%s/reference_genome.fasta"%outdir
    target_reference_genome_tmp = "%s.tmp"%target_reference_genome
    if fun.file_is_empty(target_reference_genome) or replace is True:
        fun.run_cmd("cp %s %s "%(reference_genome, target_reference_genome_tmp))
        os.rename(target_reference_genome_tmp, target_reference_genome)

    target_repeats = "%s.repeats.tab"%target_reference_genome
    target_repeats_tmp = "%s.tmp"%target_repeats
    if fun.file_is_empty(target_repeats) or replace is True:
        fun.run_cmd("cp %s.repeats.tab %s "%(reference_genome, target_repeats_tmp))
        os.rename(target_repeats_tmp, target_repeats)

    # final file 
    final_file = "%s/variants_atLeast3PASS_ploidy2.vcf"%target_smallVars_input_outdir

    if fun.file_is_empty(final_file) or replace is True:

        cmd = "%s -r %s -o %s -p 2 -sbam %s -caller all -c 5 -mchr no_mitochondria -mcode 3 -gcode 1 --repeats_table %s --remove_smallVarsCNV_nonEssentialFiles -thr %i --skip_cnv_analysis"%(varcall_cnv_pipeline, target_reference_genome, target_smallVars_input_outdir, sorted_bam, target_repeats, threads) 

        fun.run_cmd(cmd)

    print("you can run successfully the variant processing")


def test_smallVarCall_CNV_running(sorted_bam, outdir, ref_genome, gff, threads=4, mitochondrial_chromosome="mito_C_glabrata_CBS138", replace=False):

    """Takes a sorted bam (shuld have some mutations) and runs the variant calling pipeline on it"""

    # make the outdir
    fun.make_folder(outdir)

    # get the repeats
    repeats_table = fun.get_repeat_maskerDF(ref_genome, threads=4, replace=False)[1]

    for pooled_seq in [False, True]:

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


def test_SRAdb_query_downloading_and_readTrimming(outdir, reference_genome, target_taxID, replace=False, threads=4):

    """This function runs get_close_shortReads_table_close_to_taxID for the MERS coronavirus and taking the lowest coverage reads. This tests that sra tools, entrez tools, trimmomatic and fastqc work well.

    5476 is C. albicans
    1335626 is MERS
    """

    # make the outdir
    fun.make_folder(outdir)


    # set ploidy to 1
    ploidy=1

    # run with 'get_lowest_coverage_possible=True', which will take the lowest coverage datasets
    close_shortReads_table = fun.get_close_shortReads_table_close_to_taxID(target_taxID, reference_genome, outdir, ploidy, n_close_samples=2, nruns_per_sample=1, replace=replace, threads=threads, min_fraction_reads_mapped=0.0, coverage_subset_reads=0.1, min_coverage=5, job_array_mode="local", StopAfter_sampleIndexingFromSRA=False, queue_jobs="debug", max_ncores_queue=768, time_read_obtention="02:00:00", StopAfterPrefecth_of_reads=False, get_lowest_coverage_possible=True)


    # check
    df_close_shortReads_table = fun.pd.read_csv(close_shortReads_table, sep="\t")

    if set(df_close_shortReads_table.keys())!={'short_reads2', 'short_reads1', 'runID', 'sampleID'} or len(df_close_shortReads_table)!=2: raise ValueError("The close_shortReads_table %s was not created as expected"%close_shortReads_table)

    print("The system to query the SRA database, dowload and trim reads works")
    


