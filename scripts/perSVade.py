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
picard = "%s/share/picard-2.18.26-0/picard.jar"%EnvDir
samtools = "%s/bin/samtools"%EnvDir
java = "%s/bin/java"%EnvDir

# scripts that are installed under this software
varcall_cnv_pipeline = "%s/varcall_cnv_pipeline.py"%CWD

#######

description = """
Runs perSVade pipeline on an input set of paired end short ends. It is expected to be run on a coda environment and have several dependencies (see https://github.com/Gabaldonlab/perSVade). Some of these dependencies are included in the respository "installation/external_software". These are gridss (tested on version 2.8.1), clove (tested on version 0.17) and NINJA (we installed it from https://github.com/TravisWheelerLab/NINJA/releases/tag/0.95-cluster_only). If you have any trouble with these you can replace them from the source code.
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# general args
parser.add_argument("-r", "--ref", dest="ref", required=True, help="Reference genome. Has to end with .fasta. It can also be 'auto', in which case the pipeline will download the reference genome from GenBank, if the taxID is specified with --target_taxID. It can also be a specific GenBank accession to be downloaded, like  GCA_000149245.3.")
parser.add_argument("-thr", "--threads", dest="threads", default=16, type=int, help="Number of threads, Default: 16")
parser.add_argument("-o", "--outdir", dest="outdir", action="store", required=True, help="Directory where the data will be stored")
parser.add_argument("--replace", dest="replace", action="store_true", help="Replace existing files")
parser.add_argument("-p", "--ploidy", dest="ploidy", default=1, type=int, help="Ploidy, can be 1 or 2")

# different modules to be executed
parser.add_argument("--testSVgen_from_DefaulReads", dest="testSVgen_from_DefaulReads", default=False, action="store_true", help="This indicates whether to generate a report of how the generation of SV works with the default parameters on random simulations")

parser.add_argument("--close_shortReads_table", dest="close_shortReads_table", type=str, default=None, help="This is the path to a table that has 4 fields: sampleID,runID,short_reads1,short_reads2. These should be WGS runs of samples that are close to the reference genome and some expected SV. Whenever this argument is provided, the pipeline will find SVs in these samples and generate a folder <outdir>/findingRealSVs<>/SVs_compatible_to_insert that will contain one file for each SV, so that they are compatible and ready to insert in a simulated genome. This table will be used if --testRealDataAccuracy is specified, which will require at least 3 runs for each sample. It can be 'auto', in which case it will be inferred from the taxID provided by --target_taxID.")

parser.add_argument("--target_taxID", dest="target_taxID", type=int, default=None, help="This is the taxID (according to NCBI taxonomy) to which your reference genome belongs. If provided it is used to download genomes and reads.")

parser.add_argument("--n_close_samples", dest="n_close_samples", default=5, type=int, help="Number of close samples to search in case --target_taxID is provided")

parser.add_argument("--nruns_per_sample", dest="nruns_per_sample", default=3, type=int, help="Number of runs to download for each sample in the case that --target_taxID is specified. ")

parser.add_argument("--SVs_compatible_to_insert_dir", dest="SVs_compatible_to_insert_dir", type=str, default=None, help="A directory with one file for each SV that can be inserted into the reference genome in simulations. It may be created with --close_shortReads_table. If both --SVs_compatible_to_insert_dir and --close_shortReads_table are provided, --SVs_compatible_to_insert_dir will be used, and --close_shortReads_table will have no effect. If none of them are provided, this pipeline will base the parameter optimization on randomly inserted SVs (the default behavior). The coordinates have to be 1-based, as they are ready to insert into RSVsim.")

parser.add_argument("--fast_SVcalling", dest="fast_SVcalling", action="store_true", default=False, help="Run SV calling with a default set of parameters. There will not be any optimisation nor reporting of accuracy. This is expected to work almost as fast as gridss and clove together.")

# pipeline skipping options 
parser.add_argument("--skip_SVcalling", dest="skip_SVcalling", action="store_true", default=False, help="Do not run SV calling.")

# pipeline stopping options
parser.add_argument("--StopAfter_readObtentionFromSRA", dest="StopAfter_readObtentionFromSRA", action="store_true", default=False, help="Stop after obtaining reads from SRA.")
parser.add_argument("--StopAfter_genomeObtention", dest="StopAfter_genomeObtention", action="store_true", default=False, help="Stop after genome obtention.")
parser.add_argument("--StopAfter_bamFileObtention", dest="StopAfter_bamFileObtention", action="store_true", default=False, help="Stop after obtaining the BAM file of aligned reads.")


# testing options
parser.add_argument("--testRealDataAccuracy", dest="testRealDataAccuracy", action="store_true", default=False, help="Reports the accuracy  of your calling on the real data for all the WGS runs specified in --close_shortReads_table. ")

parser.add_argument("--testSimulationsAccuracy", dest="testSimulationsAccuracy", action="store_true", default=False, help="Reports the accuracy  of your calling on the simulations that are uniform, based on realSVs, and without optimisation. ")


# simulation parameter args
parser.add_argument("--nvars", dest="nvars", default=50, type=int, help="Number of variants to simulate. Note that the number of balanced translocations inserted in simulations will be always as maximum the number of gDNA chromosome-pairs implicated.")

parser.add_argument("--nsimulations", dest="nsimulations", default=2, type=int, help="The number of 'replicate' simulations that will be produced.")

parser.add_argument("--simulation_ploidies", dest="simulation_ploidies", type=str, default="haploid,diploid_hetero", help='A comma-sepparated string of the ploidies to simulate for parameter optimisation. It can have any of "haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1" ')

parser.add_argument("--range_filtering_benchmark", dest="range_filtering_benchmark", type=str, default="theoretically_meaningful", help='The range of parameters that should be tested in the SV optimisation pipeline. It can be any of large, medium, small, theoretically_meaningful or single.')

# alignment args
parser.add_argument("-f1", "--fastq1", dest="fastq1", default=None, help="fastq_1 file. Option required to obtain bam files. It can be 'auto', in which case a set of 10M reads will be generated.")
parser.add_argument("-f2", "--fastq2", dest="fastq2", default=None, help="fastq_2 file. Option required to obtain bam files. It can be 'auto', in which case a set of 10M reads will be generated.")
parser.add_argument("-sbam", "--sortedbam", dest="sortedbam", default=None, help="The path to the sorted bam file, which should have a bam.bai file in the same dir. This is mutually exclusive with providing reads")
parser.add_argument("--run_qualimap", dest="run_qualimap", action="store_true", default=False, help="Run qualimap for quality assessment of bam files. This may be inefficient sometimes because of the ")

# machine options
parser.add_argument("--run_in_slurm", dest="run_in_slurm", action="store_true", default=False, help="If provided, it will run in a job arrays the works of --testSimulationsAccuracy and --testRealDataAccuracy. This requires this script to be run on a slurm-based cluster.")

# other args
parser.add_argument("-mchr", "--mitochondrial_chromosome", dest="mitochondrial_chromosome", default="mito_C_glabrata_CBS138", type=str, help="The name of the mitochondrial chromosome. This is important if you have mitochondrial proteins for which to annotate the impact of nonsynonymous variants, as the mitochondrial genetic code is different. This should be the same as in the gff. If there is no mitochondria just put 'no_mitochondria'. If there is more than one mitochindrial scaffold, provide them as comma-sepparated IDs. It can be auto, in which case the ref genome will be parsed and automatically identified.")


# small VarCalk and CNV args
parser.add_argument("--run_smallVarsCNV", dest="run_smallVarsCNV", action="store_true", default=False, help="Will call small variants and CNV.")
parser.add_argument("-gff", "--gff-file", dest="gff", default=None, help="path to the GFF3 annotation of the reference genome. Make sure that the IDs are completely unique for each 'gene' tag. This is necessary for both the CNV analysis (it will look at genes there) and the annotation of the variants.")
parser.add_argument("-caller", "--caller", dest="caller", required=False, default="all", help="SNP caller option to obtain vcf file. options: no/all/HaplotypeCaller/bcftools/freebayes.")
parser.add_argument("-c", "--coverage", dest="coverage", default=20, type=int, help="minimum Coverage (int)")
parser.add_argument("-mcode", "--mitochondrial_code", dest="mitochondrial_code", default=3, type=int, help="The code of the NCBI mitochondrial genetic code. For yeasts it is 3. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")
parser.add_argument("-gcode", "--gDNA_code", dest="gDNA_code", default=1, type=int, help="The code of the NCBI gDNA genetic code. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . For C. albicans it is 12. ")

opt = parser.parse_args()

########################################
##### GENERAL PROCESSING OF INPUTS #####
########################################

# if replace is set remove the outdir, and then make it
if opt.replace is True: fun.delete_folder(opt.outdir)
fun.make_folder(opt.outdir)

# define the name that will be used as tag, it is the name of the outdir, without the full path
name_sample = opt.outdir.split("/")[-1]; print("working on %s"%name_sample)

# define where the reference genome will be stored
reference_genome_dir = "%s/reference_genome_dir"%(opt.outdir); fun.make_folder(reference_genome_dir)
new_reference_genome_file = "%s/reference_genome.fasta"%reference_genome_dir

if fun.file_is_empty(new_reference_genome_file) or opt.replace is True:

    # download the reference genome from GenBank given the taxID and also the gff annotation
    if opt.ref=="auto": opt.ref, opt.gff = fun.get_reference_genome_from_GenBank(opt.target_taxID, new_reference_genome_file, replace=opt.replace)

    # get by GenBank annotation
    elif opt.ref.startswith("GCA_"): opt.ref, opt.gff = fun.get_reference_genome_from_GenBank(opt.ref, new_reference_genome_file, replace=opt.replace)

    # just move the ref genome in the outdir
    else:

        # move the reference genome into the outdir, so that every file is written under outdir
        try: run_cmd("rm %s"%new_reference_genome_file)
        except: pass
        fun.run_cmd("ln -s %s %s"%(opt.ref, new_reference_genome_file))

    if opt.gff is None: print("WARNING: gff was not provided. This will be a problem if you want to annotate small variant calls")

    # check that the mitoChromosomes are in the ref
    all_chroms = {s.id for s in SeqIO.parse(opt.ref, "fasta")}
    if any([x not in all_chroms for x in opt.mitochondrial_chromosome.split(",")]) and opt.mitochondrial_chromosome!="no_mitochondria":
        raise ValueError("The provided mitochondrial_chromosomes are not in the reference genome.")

opt.ref = new_reference_genome_file

# if the gff is not none, move it to the same dir of the reference genome
target_gff = "%s/reference_genome_features.gff"%reference_genome_dir
if opt.gff is not None:

    if fun.file_is_empty(target_gff): 
        try: fun.run_cmd("rm %s"%target_gff)
        except: pass

        fun.run_cmd("ln -s %s %s"%(opt.gff, target_gff))
    
    opt.gff = target_gff

# get the genome len
genome_length = sum(fun.get_chr_to_len(opt.ref).values())
print("The genome has %.2f Mb"%(genome_length/1000000 ))

if opt.StopAfter_genomeObtention is True: 
    print("Stopping pipeline after the genome obtention.")
    sys.exit(0)


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

########################################
########################################
########################################


#####################################
############# BAM FILE ##############
#####################################

##### YOU NEED TO RUN THE BAM FILE #####

# if you stated auto in the reads, generate a 30x coverage bam file
if any([x=="auto" for x in {opt.fastq1, opt.fastq2}]):

    # define the number of reads as a function of the coverage
    read_length = 150
    total_nread_pairs = int((genome_length*30)/read_length)
    print("Simulating %.2fM reads"%(total_nread_pairs/1000000))

    sorted_bam, index_bam = fun.get_simulated_bamFile(opt.outdir, opt.ref, replace=opt.replace, threads=opt.threads, total_nread_pairs=total_nread_pairs, read_length=read_length)
    print("using simulated bam file from %s"%sorted_bam)


# normal alignment of provided reads
elif all([not x is None for x in {opt.fastq1, opt.fastq2}]):

    print("WORKING ON ALIGNMENT")
    fun.run_bwa_mem(opt.fastq1, opt.fastq2, opt.ref, opt.outdir, bamfile, sorted_bam, index_bam, name_sample, threads=opt.threads, replace=opt.replace)


else: print("Warning: No fastq file given, assuming that you provided a bam file")

#####################################
#####################################
#####################################


###########################################
############# NECESSARY FILES #############
###########################################

# check that all the important files exist
if any([fun.file_is_empty(x) for x in {sorted_bam, index_bam}]): raise ValueError("You need the sorted and indexed bam files in ")

#### bamqc
if opt.run_qualimap is True:
    
    bamqc_outdir = "%s/bamqc_out"%opt.outdir
    if fun.file_is_empty("%s/qualimapReport.html"%bamqc_outdir) or opt.replace is True:
        print("Running bamqc to analyze the bam alignment")
        qualimap_std = "%s/std.txt"%bamqc_outdir
        try: bamqc_cmd = "%s bamqc -bam %s -outdir %s -nt %i > %s 2>&1"%(qualimap, sorted_bam, bamqc_outdir, opt.threads, qualimap_std); fun.run_cmd(bamqc_cmd)
        except: print("WARNING: qualimap failed likely due to memory errors, check %s"%qualimap_std)

# First create some files that are important for any program

# Create a reference dictionary
rstrip = opt.ref.split(".")[-1]
dictionary = "%sdict"%(opt.ref.rstrip(rstrip)); tmp_dictionary = "%s.tmp"%dictionary; 
if fun.file_is_empty(dictionary) or opt.replace is True:

    # remove any previously created tmp_file
    if not fun.file_is_empty(tmp_dictionary): os.unlink(tmp_dictionary)

    print("Creating picard dictionary")
    cmd_dict = "%s -jar %s CreateSequenceDictionary R=%s O=%s TRUNCATE_NAMES_AT_WHITESPACE=true"%(java, picard, opt.ref, tmp_dictionary); fun.run_cmd(cmd_dict)   
    os.rename(tmp_dictionary , dictionary)

# Index the reference
if fun.file_is_empty("%s.fai"%opt.ref) or opt.replace is True:
    print ("Indexing the reference...")
    cmd_indexRef = "%s faidx %s"%(samtools, opt.ref); fun.run_cmd(cmd_indexRef) # This creates a .bai file of the reference

###########################################
###########################################
###########################################

if opt.StopAfter_bamFileObtention is True: 
    print("Stopping pipeline after the bamfile obtention.")
    sys.exit(0)

#####################################
##### STRUCTURAL VARIATION ##########
#####################################


#### define general args ####

# the simulation ploidies as a list
simulation_ploidies = opt.simulation_ploidies.split(",")

# the window length for all operations
fun.window_l = int(np.median([len_seq for chrom, len_seq  in fun.get_chr_to_len(opt.ref).items() if chrom not in opt.mitochondrial_chromosome.split(",")])*0.05) + 1

print("using a window length of %i"%fun.window_l)

#############################

#### test how well the finding of SVs in an assembly works ####
if opt.testSVgen_from_DefaulReads:

    outdir_test_FindSVinAssembly = "%s/test_FindSVfromDefaultSimulations"%opt.outdir
    if __name__ == '__main__': fun.test_SVgeneration_from_DefaultParms(opt.ref, outdir_test_FindSVinAssembly, sorted_bam, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, nvars=opt.nvars)

###############################################################

##### find a dict that maps each svtype to a file with a set of real SVs (real_svtype_to_file) #####
all_svs = {'translocations', 'insertions', 'deletions', 'inversions', 'tandemDuplications'}

if opt.SVs_compatible_to_insert_dir is not None and opt.fast_SVcalling is False: 
    print("using the set of real variants from %s"%opt.SVs_compatible_to_insert_dir)

    # if it is already predefined
    real_svtype_to_file = {svtype : "%s/%s.tab"%(opt.SVs_compatible_to_insert_dir, svtype) for svtype in all_svs if not fun.file_is_empty("%s/%s.tab"%(opt.SVs_compatible_to_insert_dir, svtype))}

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

        opt.close_shortReads_table = fun.get_close_shortReads_table_close_to_taxID(opt.target_taxID, opt.ref, outdir_getting_closeReads, sorted_bam, n_close_samples=opt.n_close_samples, nruns_per_sample=opt.nruns_per_sample, replace=opt.replace, threads=opt.threads, run_in_slurm=opt.run_in_slurm)

        # skip the running of the pipeline 
        if opt.StopAfter_readObtentionFromSRA:
            print("Stopping pipeline after the reads obtention from SRA")
            sys.exit(0) 

    # get the real SVs
    real_svtype_to_file = fun.get_compatible_real_svtype_to_file(opt.close_shortReads_table, opt.ref, outdir_finding_realVars, replace=opt.replace, threads=opt.threads, max_nvars=opt.nvars, mitochondrial_chromosome=opt.mitochondrial_chromosome, run_in_slurm=opt.run_in_slurm)

else: 
    print("Avoiding the simulation of real variants. Only inserting randomSV.")

    # define the set of vars as empty. This will trigger the random generation of vars
    real_svtype_to_file = {}

###################################################################################################

# test the accuracy on each of the simulations types
if opt.testSimulationsAccuracy is True: fun.report_accuracy_simulations(sorted_bam, opt.ref, "%s/testing_SimulationsAccuracy"%opt.outdir, real_svtype_to_file, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars)

# test accuracy on real data
if opt.testRealDataAccuracy is True:  fun.report_accuracy_realSVs(opt.close_shortReads_table, opt.ref, "%s/testing_RealSVsAccuracy"%opt.outdir, real_svtype_to_file, outdir_finding_realVars, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, run_in_slurm=opt.run_in_slurm)

# run the actual perSVade function optimising parameters
if opt.skip_SVcalling is False:

    SVdetection_outdir = "%s/SVdetection_output"%opt.outdir
    fun.run_GridssClove_optimising_parameters(sorted_bam, opt.ref, SVdetection_outdir, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, fast_SVcalling=opt.fast_SVcalling, real_svtype_to_file=real_svtype_to_file)

print("structural variation analysis with perSVade finished")

#####################################
#####################################
#####################################


#####################################
###### SMALL VARS AND CNV ###########
#####################################

if opt.run_smallVarsCNV:

    # define an outdir
    outdir_varcall = "%s/smallVars_CNV_output"%opt.outdir

    # define the basic cmd
    varcall_cmd = "%s -r %s --threads %i --outdir %s -p %i -sbam %s --caller %s --coverage %i --mitochondrial_chromosome %s --mitochondrial_code %i --gDNA_code %i "%(varcall_cnv_pipeline, opt.ref, opt.threads, outdir_varcall, opt.ploidy, sorted_bam, opt.caller, opt.coverage, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code)

    # add options
    if opt.replace is True: varcall_cmd += " --replace"
    if opt.gff is not None: varcall_cmd += " -gff %s"%opt.gff

    # run
    fun.run_cmd(varcall_cmd)

#####################################
#####################################
#####################################



print("perSVade Finished")


