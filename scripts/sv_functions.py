#!/usr/bin/env python

# this contains all the functions related to the definition of the environment

######################################################
###############  DEFINE ENVIRONMENT ##################
######################################################

# module imports
import os
import re
import string
import random
import pandas as pd
import numpy as np
import sys
from collections import ChainMap
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle
import itertools
import copy as cp
import re
import shutil
from datetime import date
import multiprocessing as multiproc
import scipy.stats
import numpy.polynomial.polynomial as poly
from scipy.optimize import curve_fit
from sklearn import linear_model
from statsmodels.stats import multitest
import time
from sklearn.metrics import r2_score
import time
from collections import Counter, defaultdict
import inspect
import collections
from shutil import copyfile

warnings.simplefilter(action='ignore', category=pd.core.common.SettingWithCopyWarning) # avoid the slicing warning
#pd.options.mode.chained_assignment = 'raise'

# load a specific matplotlib library for cluster envs
try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

except: import matplotlib.pyplot as plt
import seaborn as sns

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import graphics_functions as graph_fun
import smallVarsCNV_functions as fun

# EnvDir executables. These are all the ones installed with conda
JAVA = "%s/bin/java"%EnvDir
bcftools = "%s/bin/bcftools"%EnvDir
bgzip = "%s/bin/bgzip"%EnvDir
tabix = "%s/bin/tabix"%EnvDir
bedtools = "%s/bin/bedtools"%EnvDir
bwa = "%s/bin/bwa"%EnvDir
samtools = "%s/bin/samtools"%EnvDir
gatk = "%s/bin/gatk"%EnvDir # this is gatk4
freebayes = "%s/bin/freebayes"%EnvDir
genmap = "%s/bin/genmap"%EnvDir
repeat_masker = "%s/bin/RepeatMasker"%EnvDir
vcffilter = "%s/bin/vcffilter"%EnvDir
freebayes_parallel = "%s/bin/freebayes-parallel"%EnvDir
fasta_generate_regions_py = "%s/bin/fasta_generate_regions.py"%EnvDir
kmercountexact = "%s/bin/kmercountexact.sh"%EnvDir
wgsim = "%s/bin/wgsim"%EnvDir
picard_exec = "%s/bin/picard"%EnvDir
minimap2 = "%s/bin/minimap2"%EnvDir
svim = "%s/bin/svim"%EnvDir
bbmap_reformat_sh = "%s/bin/reformat.sh"%EnvDir
mosdepth = "%s/bin/mosdepth"%EnvDir
repeat_modeller_BuildDatabase = "%s/bin/BuildDatabase"%EnvDir
repeat_modeller = "%s/bin/RepeatModeler"%EnvDir
perl = "%s/bin/perl"%EnvDir
makeblastdb = "%s/bin/makeblastdb"%EnvDir
repeatmoder_dir = "%s/share/RepeatModeler"%EnvDir
repeatmasker_dir = "%s/share/RepeatMasker"%EnvDir
abblast_dir = "%s/bin"%EnvDir
cdhit_dir = "%s/bin"%EnvDir
genometools_dir = "%s/bin"%EnvDir
ltr_retriever_dir = "%s/bin"%EnvDir
mafft_dir = "%s/bin"%EnvDir
recon_dir = "%s/bin"%EnvDir
rmblast_dir = "%s/bin"%EnvDir
rscout_dir = "%s/bin"%EnvDir
trf_prgm_dir = "%s/bin/trf"%EnvDir

# executables that are provided in the repository
external_software = "%s/../installation/external_software"%CWD
gridss_run = "%s/gridss.sh"%external_software
gridss_jar = "%s/gridss-2.8.1-gridss-jar-with-dependencies.jar"%external_software
clove = "%s/clove-0.17-jar-with-dependencies.jar"%external_software
ninja_dir = "%s/NINJA-0.95-cluster_only/NINJA"%external_software

# scripts that are of this pipeline
create_random_simulatedSVgenome_R = "%s/create_random_simulatedSVgenome.R"%CWD

annotate_simpleEvents_gridssVCF_R = "%s/annotate_simpleEvents_gridssVCF.R"%CWD
analyze_svVCF = "%s/generate_files_from_svVCF.R"%CWD
analyze_svVCF_simple = "%s/generate_files_from_svVCF_simple.R"%CWD


######################################################
######################################################

####################################
######## DEFINE VARIABLES ##########
####################################

# define the strings that have to be considered as NaN in the VCF parsing
vcf_strings_as_NaNs = ['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'n/a', 'nan', 'null']

####################################
####################################
####################################

def clean_reference_genome_windows_files(reference_genome):

    """Cleans all the files under reference_genome that are windows files and bed's """

    print("removing windows files")
    ref_dir = fun.get_dir(reference_genome)
    ref_name = fun.get_file(reference_genome)

    for file in os.listdir(ref_dir):
        if file.startswith(ref_name) and "windows" in file and "bp.bed" in file : fun.remove_file("%s/%s"%(ref_dir, file))


def get_affected_region_bed_for_SVdf(svDF, svtype, interesting_chromosomes, add_interval_bp=2000):

    """This function takes a df with SVs and returns the regions that are affected in bed format, which depends on the svtype. It adds an interval arround this values.

    only consider svs were all the chroms are in interesting_chromosomes"""

    # if the df is empty just return an empty df
    if len(svDF)==0: affected_region_bed_df = pd.DataFrame()
    else:

        pass

    print(svDF)

    youshoulddefinethebedregionsofSVtypes

    affected_region_bed_df

    return affected_region_bed_df

def get_bed_df_not_overlapping_with_SVs(all_regions_bed_df, svtypes, svtype_to_svfile, bed_regions_prefix, distance_between_SVs=2000):

    """This function takes a bed df where all the regions of a genomes. It goes through a set of svs and returns the bed were there are regions where there are no SVs (at least by distance_between_SVs). It writes files under outdir"""

    # write the all regions
    all_regions_bed = "%s_all_regions.bed"%bed_regions_prefix
    all_regions_bed_df.to_csv(all_regions_bed, sep="\t", index=False, header=False)

    # define the interesting chromosomes
    interesting_chromosomes = set(all_regions_bed_df.chromosome)    

    # initialize a df with all the regions with SVs
    regions_with_SV_bed_df = pd.DataFrame(columns=["chromosome", "start", "end"])

    # intialize a dict that has the number of SVs
    svtype_to_nSVs = {svtype : 0 for svtype in svtypes}

    # go through each svtype and add regions to regions_with_SV_bed_df
    for svtype in svtypes:

        # whenever the svtype is in svtype_to_svfile it has to be added to regions_with_SV_bed_df
        if svtype in svtype_to_svfile: 

            # get the df
            svDF = pd.read_csv(svtype_to_svfile[svtype], sep="\t")

            # keep the number of SVs
            svtype_to_nSVs[svtype] += len(svDF)

            # append the bed regions
            regions_with_SV_bed_df = regions_with_SV_bed_df.append(get_affected_region_bed_for_SVdf(svDF, svtype, interesting_chromosomes, add_interval_bp=distance_between_SVs))

    # write the bed with the regions with SV
    regions_with_SV_bed = "%s_sv_regions.bed"%bed_regions_prefix
    regions_with_SV_bed_df.to_csv(regions_with_SV_bed, sep="\t", index=False, header=False)

    # get the regions in all_regions_bed that are not in regions_with_SV_bed
    regions_without_SV_bed = "%s_noSV_regions.bed"%bed_regions_prefix
    fun.run_cmd("%s subtract -a %s -b %s > %s"%(bedtools, all_regions_bed, regions_with_SV_bed, regions_without_SV_bed))

    return regions_without_SV_bed, svtype_to_nSVs

def rearrange_genomes_simulateSV(reference_genome, outdir, replace=True, nvars=100, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulated_svtype_to_svfile={}, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}):

    """Runs a simulation of nvars SVs of each type into the reference genome. It will be sepparated by gDNA and mtDNA. mtDNA will only include 5% of the gDNA vars. Everything is written to outdir. simulated_svtype_to_svfile is a dictionary that maps each svtype to a file with it. This function will insert these SVs plus the remaining ones up to nvars (which will be placed randomly in the remaining spots on the genome). Note that only balanced translocations will be simulated, as unbalanced translocations are hard to bechmark based on coverage. """

    # define the different types of chromosomes. Note that the mtDNA chromosomes will be simulated appart
    all_chromosomes = {s.id for s in SeqIO.parse(reference_genome, "fasta")}
    if mitochondrial_chromosome!="no_mitochondria": mtDNA_chromosomes = set(mitochondrial_chromosome.split(","))
    else: mtDNA_chromosomes = set()
    gDNA_chromosomes = all_chromosomes.difference(mtDNA_chromosomes)

    # map the chromosome to the length
    chrom_to_len = {s.id : len(s.seq) for s in SeqIO.parse(reference_genome, "fasta")}

    # go through each of the mtDNA and gDNA
    for type_genome, chroms in [("mtDNA", mtDNA_chromosomes), ("gDNA", gDNA_chromosomes)]:

        # if there are chroms just continue
        if len(chroms)==0: continue

        # if the genome is mtDNA you shoudl simulate less vars
        if type_genome=="gDNA": vars_to_simulate = nvars
        else: vars_to_simulate = int(nvars*0.05) + 1

        # define the outdir
        genome_outdir = "%s/simulation_%s"%(outdir, type_genome); fun.make_folder(genome_outdir)

        # get the genome 
        genome_file = "%s/genome.fasta"%genome_outdir
        SeqIO.write([c for c in SeqIO.parse(reference_genome, "fasta") if c.id in chroms], genome_file, "fasta")

        # define a bed file with all the data
        all_regions_bed_df = pd.DataFrame({chrom: {"start":0, "end":chrom_to_len[chrom]} for chrom in chroms}).transpose()
        all_regions_bed_df["chromosome"] = all_regions_bed_df.index

        # get the regions without SV where simulations should be placed
        bed_regions_prefix = "%s/bed_regions"%genome_outdir
        regions_without_SV_bed, svtype_to_nSVs = get_bed_df_not_overlapping_with_SVs(all_regions_bed_df[["chromosome", "start", "end"]], svtypes, simulated_svtype_to_svfile, bed_regions_prefix)

        # simulate random SVs into regions without previous SVs 
        random_sim_dir = "%s/random_SVs"%genome_outdir; fun.delete_folder(random_sim_dir); fun.make_folder(random_sim_dir)

        # get the cmd of the simulation
        backbone_cmd = "%s --input_genome %s --outdir %s"%(create_random_simulatedSVgenome_R, genome_file, random_sim_dir)

        # add the number of each SV that should be added
        for svtype, number_alreadyGeneratedSVs in svtype_to_nSVs.items():



        #

        std_rearranging_genome = "%s/simulation_std.txt"%random_sim_dir
        run_cmd("%s --input_genome %s --outdir %s --regions_bed %s --mitochondrial_chromosome %s --max_time_rearrangement %i  --percBalancedTrans %.2f > %s 2>&1"%(simulateSV_R, reference_genome, sim_outdir, sim_type_to_regionsBed[simulation_type], mitochondrial_chromosome, max_time_rearrangement, percBalancedTrans, std_rearranging_genome))


        argp = arg_parser("Takes a genome and generates a simulated genome with rearrangements with the known rearrangements in outdir. It will generate these rearrangements and the rearranged genome under outdir, only in the regions that are provided (--regions_bed)")

argp = add_argument(argp, "--number_Dup", default=100, help="The number of duplications to generate")
argp = add_argument(argp, "--number_Ins", default=100, help="The number of insertions to generate")
argp = add_argument(argp, "--number_Inv", default=100, help="The number of inversions to generate")
argp = add_argument(argp, "--number_Del", default=100, help="The number of deletions to generate")
argp = add_argument(argp, "--number_Tra", default=100, help="The number of translocations to generate")







        adljhdaldah




        print(all_regions_bed)

        slkjfsslfhj






def run_GridssClove_optimising_parameters(sorted_bam, reference_genome, outdir, threads=4, replace=False, window_l=1000, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_types=["uniform", "biased_towards_repeats"], target_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], replace_covModelObtention=False, range_filtering_benchmark="theoretically_meaningful", known_genomes_withSV_and_shortReads_table=True, expected_ploidy=1, nvars=100):

    """
    Takes some aligned reads and runs the GridssPipeline optimising the parameters of GRIDSS filtering. These are the different parameters of the function:

    - sorted_bam: the path to a sorted and indexed bam where we want to find the SV
    - reference_genome: the fasta of the reference genome
    - outdir: a directory where all the sample-specific files will be written. All the reference_genome-specific files will be written in the same dir as the reference_genome is
    - window_l: the length of the windows to generate the coverage model
    - n_simulated_genomes is the number of genomes that will be simulated as simulation replicates
    - mitochondrial_chromosome is the name of the mtDNA chromosome in the reference genome. This is passed to some parts of the pipeline to consider differently the gDNA and the mtDNA (the SV simulation). It can be comma-sepparated. It can be "no_mitochondria" if there is no mitochondrial chromosome
    - simulation_types indicates the types of SV performed
    - target_ploidies indicates which poplulations or ploidies have to be simulated. 2ref_1sv means that we will simulate a genome that has 2 reference genomes and 1 genome under structural variation (this can be a population sequencing)
    - replace_covModelObtention indicates whether the process of predicting coverage from seq features has to be replaced
    - range_filtering_benchmark indicates which type of simulation will be performed, it can be "large", "medium", "small", "single", "theoretically_meaningful". This is passed to benchmark_GridssClove_for_knownSV.
    - known_genomes_withSV_and_shortReads_table is a file with a table that has three fields: ID,assembly,shoort_reads1,short_reads2 . This can be, for example a set of NANOPORE assemblies of Candida glabrata and the corresponding short reads' sequencing in YPD
    - expected_ploidy is a number that states the expected ploidy. 
    - nvars deteremines the number of SVs to simulate in each simulation of each type. The mtDNA will get 5% of these.

    This is the workflow:

    1) Test the pipeline that finds SV in a given genome assembly. Generates nvars of each type for n_simulated_genomes and tests  the accuracy of each simulation type on it.

    """

    # prepare files
    fun.make_folder(outdir)

    ##### test how well the finding of SVs in an assembly works #####
    outdir_test_FindSVinAssembly = "%s/test_FindSVinAssembly"%outdir; fun.make_folder(outdir_test_FindSVinAssembly)

    # go through each simulation
    for simID in range(n_simulated_genomes):

        # define outdir 
        outdir_sim = "%s/simulation_%i"%(outdir_test_FindSVinAssembly, simID); fun.make_folder(outdir_sim)

        # generate genome with simulated SVs
        rearrange_genomes_simulateSV(reference_genome, outdir_sim, replace=replace, nvars=nvars, mitochondrial_chromosome=mitochondrial_chromosome)







    #simType_to_rearrangedGenome = rearrange_genomes_simulateSV(reference_genome, simulation_outdir, repeats_1kbwindow_bed=repeats_regions_bed, replace=replace, simulation_types=simulation_types, mitochondrial_chromosome=mitochondrial_chromosome) # !!!! i have to change the representation of translocations with different 3' and 5' end



    #################################################################



    


    # at the end clean the reference genome
    clean_reference_genome_windows_files(reference_genome)






