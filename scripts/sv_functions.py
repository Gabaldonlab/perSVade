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
create_targeted_simulatedSVgenome_R = "%s/create_targeted_simulatedSVgenome.R"%CWD

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


def get_affected_region_bed_for_SVdf(svDF, svtype, interesting_chromosomes, add_interval_bp=1000):

    """This function takes a df with SVs and returns the regions that are affected in bed format, which depends on the svtype. It adds an interval arround this values. It also returns the number 

    only consider svs were ANY of the chroms are in interesting_chromosomes"""


    # load df
    if type(svDF)==str: svDF = pd.read_csv(svDF, sep="\t")

    # if the df is empty just return an empty df
    if len(svDF)==0: 
        affected_region_bed_df = pd.DataFrame(columns=["chromosome", "start", "end"])
        nSVs_in_interesting_chromosomes = 0
    else:

        # simple svtypes
        if svtype in {"deletions", "inversions", "tandemDuplications"}: 

            # everything goes to the bed
            affected_region_bed_df = svDF.rename(columns={"Chr":"chromosome", "Start":"start", "End":"end"})[["chromosome", "start", "end"]]

            # all SVs count
            nSVs_in_interesting_chromosomes = sum(affected_region_bed_df.chromosome.isin(interesting_chromosomes))


        # insertions (the b region will only be included as start)
        elif svtype=="insertions":

            # count the SVs
            nSVs_in_interesting_chromosomes = sum(svDF.apply(lambda r: r["ChrA"] in interesting_chromosomes or r["ChrB"] in interesting_chromosomes, axis=1))

            # initialize with the origin chromosome
            affected_region_bed_df = svDF.rename(columns={"ChrA":"chromosome", "StartA":"start", "EndA":"end"})[["chromosome", "start", "end"]]

            # add the end only with the start
            chrB_df = svDF.rename(columns={"ChrB":"chromosome", "StartB":"start"}); chrB_df["end"] = chrB_df["start"]
            affected_region_bed_df = affected_region_bed_df.append(chrB_df[["chromosome", "start", "end"]], sort=True)

        elif svtype=="translocations":

            # count the SVs
            nSVs_in_interesting_chromosomes = sum(svDF.apply(lambda r: r["ChrA"] in interesting_chromosomes or r["ChrB"] in interesting_chromosomes, axis=1))

            # initialize with the second position in the As chromosome (where the A bp is)
            affected_region_bed_df = svDF.rename(columns={"ChrA":"chromosome", "EndA":"bp_pos"})[["chromosome", "bp_pos"]]

            # add the position in the A chromosome that is not 1 (where the B bp is)
            chrB_df = svDF.rename(columns={"ChrB":"chromosome", "StartB":"start", "EndB":"end"})
            chrB_df["bp_pos"] = chrB_df.apply(lambda r: min({r["start"], r["end"]}.difference({1})), axis=1)

            # define the start and the end to be the bp pos
            chrB_df["start"] = chrB_df.bp_pos
            chrB_df["end"] = chrB_df.bp_pos

            # merge both regions
            affected_region_bed_df = affected_region_bed_df.append(chrB_df[["chromosome", "start", "end"]], sort=True)


        else: raise ValueError("%s is not valid"%(svtype))


        # add the confidence arround each region
        affected_region_bed_df["start"] = (affected_region_bed_df.start - add_interval_bp).apply(lambda x: max([0,x]))
        affected_region_bed_df["end"] = affected_region_bed_df.start + add_interval_bp


    return affected_region_bed_df[["chromosome", "start", "end"]], nSVs_in_interesting_chromosomes

def get_bed_df_not_overlapping_with_SVs(all_regions_bed_df, svtypes, svtype_to_svfile, bed_regions_prefix, distance_between_SVs=1000):

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

            regions_with_SV_bed_svtype_df, nSVs = get_affected_region_bed_for_SVdf(svtype_to_svfile[svtype], svtype, interesting_chromosomes, add_interval_bp=distance_between_SVs)

            # keep the number of SVs
            svtype_to_nSVs[svtype] += nSVs

            # append the bed regions
            regions_with_SV_bed_df = regions_with_SV_bed_df.append(regions_with_SV_bed_svtype_df, sort=True)

    # write the bed with the regions with SV
    regions_with_SV_bed = "%s_sv_regions.bed"%bed_regions_prefix
    regions_with_SV_bed_df[["chromosome", "start", "end"]].to_csv(regions_with_SV_bed, sep="\t", index=False, header=False)

    # get the regions in all_regions_bed that are not in regions_with_SV_bed
    regions_without_SV_bed = "%s_noSV_regions.bed"%bed_regions_prefix
    fun.run_cmd("%s subtract -a %s -b %s > %s"%(bedtools, all_regions_bed, regions_with_SV_bed, regions_without_SV_bed))

    return regions_without_SV_bed, svtype_to_nSVs


def format_translocation_row_simulateSV(r, chr_to_len, start_pos=1):

    """balanced translocations have chromosomes sorted alphabetically and the A region has 1-BP and B region has BP-end. Start pos indicates where the start is"""

    # initiaize edited row to put the content
    er = {}

    # edit balanced translocations
    if r["Balanced"] is True:

        # define the important fields
        important_fields = {"Chr", "Start", "End", "Size", "BpSeq"}.intersection(set([k.rstrip("A").rstrip("B") for k in r.keys()]))

        # map each END to the coordinates
        coords_df = cp.deepcopy(pd.DataFrame({end : {field :  r["%s%s"%(field, end)] for field in important_fields} for end in {"A", "B"}}).transpose()) # columns are the coords and the index is each end

        # add the breakpoint
        coords_df["breakpoint"] = coords_df.apply(lambda c: sorted({c["End"], c["Start"]}.difference({start_pos, chr_to_len[c["Chr"]]}))[0], axis=1)
        coords_df = coords_df.set_index("Chr", drop=False)

        # add the segment which is copied (3end or 5end)
        start_to_segment = {**{start_pos:"5end"}, **{s:"3end" for s in set(coords_df["Start"]).difference({start_pos})}}
        coords_df["segment"] = coords_df.Start.apply(lambda x: start_to_segment[x])

        # get the sorted chromosome
        chrA_sorted, chrB_sorted = sorted(coords_df["Chr"])
        er["ChrA"] = chrA_sorted
        er["ChrB"] = chrB_sorted

        # if the ends are the same, define the 5' ends to be joined, in a way that is alphabeticaally sorted by chromosome
        if set(coords_df.segment)=={"3end"} or set(coords_df.segment)=={"5end"}:
            er["StartA"] = start_pos
            er["StartB"] = start_pos
            er["EndA"] = coords_df.loc[chrA_sorted, "breakpoint"]
            er["EndB"] = coords_df.loc[chrB_sorted, "breakpoint"]

        # when there are different ends, you first put the 5' of the chrA_sorted, which is replaced with the 3' end of the chrB_sorte
        else:
            er["StartA"] = start_pos
            er["StartB"] = coords_df.loc[chrB_sorted, "breakpoint"]
            er["EndA"] = coords_df.loc[chrA_sorted, "breakpoint"]
            er["EndB"] = chr_to_len[chrB_sorted]

        # calculate the size
        if "SizeA" in r.keys():
            er["SizeA"] = er["EndA"] - er["StartA"] + start_pos
            er["SizeB"] = er["EndB"] - er["StartB"] + start_pos

        # add the bpseq as XXX
        if "BpSeqA" in r.keys():
            er["BpSeqA"] = "N"*len(r["BpSeqA"])
            er["BpSeqB"] = "N"*len(r["BpSeqB"])

    else: er = {field : r[field] for field in ["StartA", "EndA", "SizeA", "ChrA", "StartB", "EndB", "SizeB", "ChrB", "BpSeqA", "BpSeqB"] if field in r.keys()}

    # change the fields
    for k, v in er.items():
        if pd.isna(v): er[k] = ""

    # add extra things
    if "Name" in r.keys(): er["Name"] = r["Name"]
    er["Balanced"] = str(r["Balanced"]).upper()

    # add X for empty BpSeqA 
    if "BpSeqA" in er.keys():
        for letter in {"A", "B"}:
            if er["BpSeq%s"%letter]=="": er["BpSeq%s"%letter] = "X"

    # return series
    return pd.Series(er)

def rewrite_translocations_uniformizedFormat_simulateSV(translocations_file, reference_genome):

    """Takes the translocations file output by simulateSV and rewrites it in a way that the balanced translocations have chromosomes sorted alphabetically and the A region has 1-BP and B region has BP-end"""

    print("rewriting %s"%translocations_file)

    # keep the unmodified version
    fun.run_cmd("cp %s %s.unmodified"%(translocations_file, translocations_file))

    # define chromosome_to_length
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

    # load translocations
    df = pd.read_csv(translocations_file, sep="\t")

    # define a function that takes a translocation and formats it
    df_corrected = df.apply(lambda r: format_translocation_row_simulateSV(r, chr_to_len), axis=1)

    df_corrected = df_corrected[[c for c in ["Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB", "Balanced", "BpSeqA", "BpSeqB"] if c in df_corrected.columns]]

    # write the translocations
    df_corrected.to_csv(translocations_file, sep="\t", header=True, index=False)

def change_EmptyString_to_X(string):

    """Takes an empty string and replaces it with X"""

    if pd.isna(string) or len(string)==0 or string=="": return "X"
    else: return string


def rewrite_insertions_uniformizedFormat_simulateSV(insertions_file):

    """Takes an insertions file and rewrites the fact that the BpSeq is empty"""

    print("rewriting %s"%insertions_file)

    # load df
    df = pd.read_csv(insertions_file, sep="\t")

    # correct BpSeq
    for field in ["BpSeqA", "BpSeqB_5prime", "BpSeqB_3prime"]: 
        if field in df.keys(): df[field] = df[field].apply(change_EmptyString_to_X)
    df["Copied"] = df.Copied.apply(lambda x: str(x).upper())

    # write as corrected
    df.to_csv(insertions_file, sep="\t", header=True, index=False)

def rearrange_genomes_simulateSV(reference_genome, outdir, replace=True, nvars=50, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulated_svtype_to_svfile={}, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}):

    """Runs a simulation of nvars SVs of each type into the reference genome. It will be sepparated by gDNA and mtDNA. mtDNA will only include 5% of the gDNA vars. Everything is written to outdir. simulated_svtype_to_svfile is a dictionary that maps each svtype to a file with it. This function will insert these SVs plus the remaining ones up to nvars (which will be placed randomly in the remaining spots on the genome). Note that only balanced translocations will be simulated, as unbalanced translocations are hard to bechmark based on coverage. """

    # define the final outdirs
    final_simulated_SVs_dir = "%s/final_simulated_SVs"%(outdir); 
    final_rearranged_genome = "%s/rearranged_genome.fasta"%final_simulated_SVs_dir
    final_rearranged_genome_finalFile = "%s.performed"%(final_rearranged_genome)

    if fun.file_is_empty(final_rearranged_genome_finalFile) or replace is True:

        # make the folder again
        fun.delete_folder(outdir)
        fun.make_folder(outdir)
        fun.make_folder(final_simulated_SVs_dir)

        # define the different types of chromosomes. Note that the mtDNA chromosomes will be simulated appart
        all_chromosomes = {s.id for s in SeqIO.parse(reference_genome, "fasta")}
        if mitochondrial_chromosome!="no_mitochondria": mtDNA_chromosomes = set(mitochondrial_chromosome.split(","))
        else: mtDNA_chromosomes = set()
        gDNA_chromosomes = all_chromosomes.difference(mtDNA_chromosomes)

        # map the chromosome to the length
        chrom_to_len = {s.id : len(s.seq) for s in SeqIO.parse(reference_genome, "fasta")}

        # initialize a df where each svtype is mapped against a 
        final_svtype_to_svDF = {}
        for svtype in svtypes:
            if svtype in simulated_svtype_to_svfile: 

                # get the df and change the names
                svDF = pd.read_csv(simulated_svtype_to_svfile[svtype], sep="\t")
                svDF["Name"] = svDF.Name + "_realSV"

                # keep all vals but the non real ones
                final_svtype_to_svDF[svtype] = svDF[[c for c in svDF.keys() if "BpSeq" not in c]]

            else: final_svtype_to_svDF[svtype] = pd.DataFrame()


        ###### GENERATE ALL THE RANDOM SIMULATIONS THAT ARE NECESSARY TO ADD ON simulated_svtype_to_svfile ########
        print("generating random simulations")

        # go through each of the mtDNA and gDNA
        for type_genome, chroms in [("mtDNA", mtDNA_chromosomes), ("gDNA", gDNA_chromosomes)]:
            print(type_genome)

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
            all_regions_bed_df = all_regions_bed_df[["chromosome", "start", "end"]]

            # get the regions without SV where simulations should be placed
            bed_regions_prefix = "%s/bed_regions"%genome_outdir
            regions_without_SV_bed, svtype_to_nSVs = get_bed_df_not_overlapping_with_SVs(all_regions_bed_df, svtypes, simulated_svtype_to_svfile, bed_regions_prefix)

            # simulate random SVs into regions without previous SVs 
            random_sim_dir = "%s/random_SVs"%genome_outdir; 

            #### GET THE RANDOM INS,INV,DEL,TRA ####

            if any([fun.file_is_empty("%s/%s.tab"%(random_sim_dir, svtype)) for svtype in {"insertions", "deletions", "translocations", "inversions"}]) or replace is True:

                # make and delete the folder
                fun.delete_folder(random_sim_dir); fun.make_folder(random_sim_dir)

                # get the cmd of the simulation
                randomSV_cmd = "%s --input_genome %s --outdir %s --regions_bed %s"%(create_random_simulatedSVgenome_R, genome_file, random_sim_dir, regions_without_SV_bed)

                # add the number of each SV that should be added
                svtype_to_arg = {"insertions":"number_Ins", "deletions":"number_Del", "inversions":"number_Inv", "translocations":"number_Tra", "tandemDuplications":"number_Dup"}
                #svtype_to_arg = {"insertions":"number_Ins", "deletions":"number_Del", "inversions":"number_Inv", "translocations":"number_Tra"}
            
                for svtype, number_alreadyGeneratedSVs in svtype_to_nSVs.items(): 
                    if svtype not in svtype_to_arg: continue

                    randomSV_cmd += " --%s %i"%(svtype_to_arg[svtype], max([1, (vars_to_simulate-svtype_to_nSVs[svtype])]))

                # run the random simulation
                std_rearranging_genome = "%s/simulation_std.txt"%random_sim_dir
                fun.run_cmd("%s > %s 2>&1"%(randomSV_cmd, std_rearranging_genome))

                # edit the translocations so that the balanced ones are sorted
                translocations_file = "%s/translocations.tab"%random_sim_dir
                if fun.file_is_empty(translocations_file): open(translocations_file, "w").write("\t".join(["Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB", "Balanced", "BpSeqA", "BpSeqB"])) # this needs to be 

            ########################################

            # add the simulations into simulated_svtype_to_svDF
            for svtype in final_svtype_to_svDF.keys():
                svDF = final_svtype_to_svDF[svtype]

                # get the new sv
                new_svDF = pd.read_csv("%s/%s.tab"%(random_sim_dir, svtype), sep="\t")
                new_svDF = new_svDF[[c for c in new_svDF.keys() if "BpSeq" not in c]]

                # add the name
                new_svDF["Name"] = new_svDF.Name + "_sim_%s"%type_genome

                # append 
                final_svtype_to_svDF[svtype] = svDF.append(new_svDF, sort=True)

        ###############################################################################################################

        ####### generate a rearranged genome with all the simulations in final_svtype_to_svDF #########
        print("inserting these random simulations")

        # initialize a cmd to create the simulated genome
        targetSV_cmd = "%s --input_genome %s --output_genome %s"%(create_targeted_simulatedSVgenome_R, reference_genome, final_rearranged_genome)

        # write the SVs into files and add to the cmd
        for svtype, svDF in final_svtype_to_svDF.items(): 

            # write file
            svfile = "%s/%s.tab"%(final_simulated_SVs_dir, svtype)
            svDF.to_csv(svfile, sep="\t", header=True, index=False)

            # get cmd
            targetSV_cmd += " --%s_file %s"%(svtype, svfile)

        # run the cmd
        std_rearranging_genome = "%s/simulation_std.txt"%final_simulated_SVs_dir
        fun.run_cmd("%s > %s 2>&1"%(targetSV_cmd, std_rearranging_genome))

        # edit the insertions
        insertions_file = "%s/insertions.tab"%final_simulated_SVs_dir
        rewrite_insertions_uniformizedFormat_simulateSV(insertions_file)

        # rewrite the variants so that they are optimal for comparison 
        translocations_file = "%s/translocations.tab"%final_simulated_SVs_dir
        rewrite_translocations_uniformizedFormat_simulateSV(translocations_file, genome_file)

        # write a file that indicates that this has finsihed
        open(final_rearranged_genome_finalFile, "w").write("finsihed")

    ###############################################################################################

    # return the genome and the set of SV dict
    final_svtype_to_svfile = {svtype : "%s/%s.tab"%(final_simulated_SVs_dir, svtype) for svtype in svtypes}

    return final_svtype_to_svfile, final_rearranged_genome



def generate_tables_of_SV_between_genomes_gridssClove(query_genome, reference_genome, replace=False, threads=4, coverage=10, insert_size=250, read_lengths=[kb*1000 for kb in [0.5, 0.7, 0.9, 1, 1.3, 1.5, 2]], error_rate=0.0, expected_ploidy=1):

    """Takes a bam file with aligned reads or genomes and generates calls, returning a dict that maps variation type to variants"""

    # first run svim under the outdir of the aligned reads
    working_dir = "%s/findingSVlongReadsSimulation_ouptut_%s_against_%s"%(get_dir(query_genome), get_file(query_genome), get_file(reference_genome))
    print("generating svtables into %s"%working_dir)
    if replace is True: delete_folder(working_dir); 
    make_folder(working_dir)

    ########## SIMULATE LONG READS AS IF THEY WERE PAIRED. UNIFORM #############

    # define the output
    ID = ("%ix_i%i_rl%s"%(coverage, insert_size, "-".join([get_int_or_float_as_text(x/1000) for x in read_lengths])))
    
    all_reads_1 = "%s/uniformSim_%s_all_reads1.fasta"%(working_dir, ID)
    all_reads_2 = "%s/uniformSim_%s_all_reads2.fasta"%(working_dir, ID)

    all_reads_1_fqgz = "%s/uniformSim_%s_all_reads1.fq.gz"%(working_dir, ID)
    all_reads_2_fqgz = "%s/uniformSim_%s_all_reads2.fq.gz"%(working_dir, ID)

    if any([file_is_empty(x) for x in [all_reads_1_fqgz, all_reads_2_fqgz]]) or replace is True:
        print("simulating long paired reads")

        # first simulate reads in parallel for each chromosome, as if they were paired
        inputs_fn = [(chr_obj, coverage, insert_size, read_lengths) for chr_obj in SeqIO.parse(query_genome, "fasta")]

        with multiproc.Pool(multiproc.cpu_count()) as pool:
            all_reads_list_tuples = pool.starmap(simulate_pairedEndReads_per_chromosome_uniform, inputs_fn)
            pool.close()

        #all_reads_list_tuples = list(map(lambda x: simulate_pairedEndReads_per_chromosome_uniform(x[0], x[1], x[2]), inputs_fn))
        reads_objects_1 = make_flat_listOflists([r[0] for r in all_reads_list_tuples])
        reads_objects_2 = make_flat_listOflists([r[1] for r in all_reads_list_tuples])
        print("There are %i reads"%len(reads_objects_1))

        # get the reads writen
        print("writing reads")
        all_reads_1_tmp = "%s.tmp"%all_reads_1; all_reads_2_tmp = "%s.tmp"%all_reads_2

        SeqIO.write(reads_objects_1, all_reads_1, "fasta")
        SeqIO.write(reads_objects_2, all_reads_2, "fasta")

        # get them as fastq gzipped
        all_reads_1_fqgz = convert_fasta_to_fqgz(all_reads_1, replace=replace)
        all_reads_2_fqgz = convert_fasta_to_fqgz(all_reads_2, replace=replace)

    # run bwa mem
    bamfile = "%s/uniformSim_%s_aligned_reads.bam"%(working_dir, ID)
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam

    run_bwa_mem(all_reads_1_fqgz, all_reads_2_fqgz, reference_genome, working_dir, bamfile, sorted_bam, index_bam, name_sample="uniformSim_%s"%ID, threads=threads, replace=replace)

    # plot the coverage across genome
    plot_coverage_across_genome_pairedEndReads(sorted_bam, reference_genome, replace=replace, window_l=10000)

    ############################################################################
    
    # define the gridss filters according to the freq, which is related to the expected ploidy
    gridss_filters_dict = default_filtersDict_gridss
    gridss_filters_dict["min_af"] = max([0.0, (1/expected_ploidy)*0.25])
    print("Filtering out when any AF is below %.3f"%(gridss_filters_dict["min_af"]))

    # run the gridss and clove pipeline with high-confidence parameters
    gridss_outdir = "%s/%s_gridss_outdir"%(working_dir, ID)
    SV_dict, df_gridss =  run_gridssClove_given_filters(sorted_bam, reference_genome, gridss_outdir, coverage, replace=replace, threads=threads, median_insert_size=insert_size, gridss_filters_dict=gridss_filters_dict, replace_FromGridssRun=False) # DEBUG. The replace_FromGridssRun in True would be debugging is to replace from the GRIDSS run step


    # remove all the chromosomal bam files and coverage measurements
    for file in os.listdir(working_dir):
        filepath = "%s/%s"%(working_dir, file)

        if filepath.startswith("%s."%sorted_bam) and filepath!=index_bam and filepath!="%s.coverage_per_window.tab"%sorted_bam: remove_file(filepath)

    return SV_dict



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
        sim_svtype_to_svfile, rearranged_genome = rearrange_genomes_simulateSV(reference_genome, outdir_sim, replace=replace, nvars=nvars, mitochondrial_chromosome=mitochondrial_chromosome)





        khvkhg



        # get the variants from simulating reads
        predicted_tablesSV = generate_tables_of_SV_between_genomes_gridssClove(rearranged_genome, reference_genome, replace=replace, threads=threads, expected_ploidy=expected_ploidy)








    #simType_to_rearrangedGenome = rearrange_genomes_simulateSV(reference_genome, simulation_outdir, repeats_1kbwindow_bed=repeats_regions_bed, replace=replace, simulation_types=simulation_types, mitochondrial_chromosome=mitochondrial_chromosome) # !!!! i have to change the representation of translocations with different 3' and 5' end



    #################################################################



    


    # at the end clean the reference genome
    clean_reference_genome_windows_files(reference_genome)






