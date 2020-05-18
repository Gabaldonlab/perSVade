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
import igraph

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


# define default parameters for gridss filtering
default_filtersDict_gridss = {"min_Nfragments":0, "min_af":0.25, "wrong_FILTERtags":("",), "filter_polyGC":False, "filter_noSplitReads":False, "filter_noReadPairs":False, "maximum_strand_bias":1.0, "maximum_microhomology":1000000000, "maximum_lenght_inexactHomology":100000000, "range_filt_DEL_breakpoints":(0,1), "min_length_inversions":0, "dif_between_insert_and_del":0, "max_to_be_considered_small_event":1, "wrong_INFOtags":('IMPRECISE',), "min_size":50, "min_af_EitherSmallOrLargeEvent":0.25} # the minimum af is 0.25 to include both heterozygous and homozygous vars as default



####################################
####################################
####################################

def id_generator(size=10, chars=string.ascii_uppercase + string.digits, already_existing_ids=set()):

    """ already_existing_ids is a set that indicates whihc IDs can't be picked """

    ID = ''.join(random.choice(chars) for _ in range(size))
    while ID in already_existing_ids:
        ID = ''.join(random.choice(chars) for _ in range(size))

    return ID

def run_cmd(cmd):

    """Runs a cmd under the VarCall_CNV_env environment, as defined in CONDA_ACTIVATING_CMD"""

    out_stat = os.system(cmd); 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd, out_stat))

def get_dir(filename): return "/".join(filename.split("/")[0:-1])

def get_file(filename): return filename.split("/")[-1]

def save_object(obj, filename):
    
    """ This is for saving python objects """
    
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    
    """ This is for loading python objects """
        
    return pickle.load(open(filename,"rb"))

def file_is_empty(path): 
    
    """ask if a file is empty or does not exist """
    
    if not os.path.isfile(path):
        return_val = True
    elif os.stat(path).st_size==0:
        return_val = True
    else:
        return_val = False
            
    return return_val

def remove_file(f):

    if os.path.isfile(f): os.unlink(f)

def delete_folder(f):

    if os.path.isdir(f): shutil.rmtree(f)

def make_folder(f):

    if not os.path.isdir(f): os.mkdir(f)

def delete_file_or_folder(f):

    """Takes a path and removes it"""

    if os.path.isdir(f): shutil.rmtree(f)
    if os.path.isfile(f): os.unlink(f)

def run_bwa_mem(fastq1, fastq2, ref, outdir, bamfile, sorted_bam, index_bam, name_sample, threads=1, replace=False):

    """Takes a set of files and runs bwa mem getting sorted_bam and index_bam"""

    if file_is_empty(sorted_bam) or file_is_empty(index_bam) or replace is True:

        #index fasta
        index_files = ["%s.%s"%(ref, x) for x in ["amb", "ann", "bwt", "pac", "sa"]]

        if any([file_is_empty(x) for x in index_files]) or replace is True:
            print("Indexing fasta")

            # create a branch reference, which will have a tag that is unique to this run. This is important since sometimes you run this pipeline in parallel, and this may give errors in fasta indexing.
            branch_ref = "%s.%s.fasta"%(ref, id_generator())
            shutil.copy2(ref, branch_ref)

            # run indexing in the copy
            cmd_indexFasta = "%s index %s"%(bwa, branch_ref); run_cmd(cmd_indexFasta) # creates a set of indexes of fasta
            index_files_branch = ["%s.%s"%(branch_ref, x) for x in ["amb", "ann", "bwt", "pac", "sa"]]

            # rename each of the indices so that it matches the ref format
            for branchIDX, realIDX in dict(zip(index_files_branch, index_files)).items(): os.rename(branchIDX, realIDX)

            # rm the branch
            os.unlink(branch_ref)

        #BWA MEM --> get .sam
        samfile = "%s/aligned_reads.sam"%outdir;
        if file_is_empty(samfile) or replace is True:

            # remove previuous generated temporary file
            if os.path.isfile("%s.tmp"%samfile): os.unlink("%s.tmp"%samfile)

            print("Running bwa mem")
            cmd_bwa = '%s mem -R "@RG\\tID:%s\\tSM:%s" -t %i %s %s %s > %s.tmp'%(bwa, name_sample, name_sample, threads, ref, fastq1, fastq2, samfile); run_cmd(cmd_bwa)
            os.rename("%s.tmp"%samfile , samfile)

        # convert to bam 
        if file_is_empty(bamfile) or replace is True:
            print("Converting to bam")
            cmd_toBAM = "%s view -Sbu %s > %s.tmp"%(samtools, samfile, bamfile); run_cmd(cmd_toBAM)
            os.rename("%s.tmp"%bamfile , bamfile)

            # remove the sam
            print("Removing sam"); os.unlink(samfile)

        # sorting bam
        if file_is_empty(sorted_bam) or replace is True:
            print("Sorting bam")

            # remove all temporary files generated previously in samtools sort (they'd make a new sort to be an error)
            for outdir_file in os.listdir(outdir): 
                fullfilepath = "%s/%s"%(outdir, outdir_file)
                if outdir_file.startswith("aligned_reads") and ".tmp." in outdir_file: os.unlink(fullfilepath)

            # sort
            cmd_sort = "%s sort --threads %i -o %s.tmp %s"%(samtools, threads, sorted_bam, bamfile); run_cmd(cmd_sort)
            os.rename("%s.tmp"%sorted_bam , sorted_bam)

            # remove the the raw bam file
            print("Removing unsorted bam"); os.unlink(bamfile)

    # indexing bam
    if file_is_empty(index_bam) or replace is True:
        print("Indexing bam")
        cmd_indexBam = "%s index -@ %i %s"%(samtools, threads, sorted_bam); run_cmd(cmd_indexBam)   # creates a .bai of sorted_bam


    print("ALIGNMENT STEP WAS CORRECTLY PERFORMED")

def make_flat_listOflists(LoL):

    return list(itertools.chain.from_iterable(LoL))

def clean_reference_genome_windows_files(reference_genome):

    """Cleans all the files under reference_genome that are windows files and bed's """

    print("removing windows files")
    ref_dir = get_dir(reference_genome)
    ref_name = get_file(reference_genome)

    for file in os.listdir(ref_dir):
        if file.startswith(ref_name) and "windows" in file and "bp.bed" in file : remove_file("%s/%s"%(ref_dir, file))

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
    run_cmd("%s subtract -a %s -b %s > %s"%(bedtools, all_regions_bed, regions_with_SV_bed, regions_without_SV_bed))

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
    run_cmd("cp %s %s.unmodified"%(translocations_file, translocations_file))

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

    if file_is_empty(final_rearranged_genome_finalFile) or replace is True:

        # make the folder again
        delete_folder(outdir)
        make_folder(outdir)
        make_folder(final_simulated_SVs_dir)

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
            genome_outdir = "%s/simulation_%s"%(outdir, type_genome); make_folder(genome_outdir)

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

            if any([file_is_empty("%s/%s.tab"%(random_sim_dir, svtype)) for svtype in {"insertions", "deletions", "translocations", "inversions"}]) or replace is True:

                # make and delete the folder
                delete_folder(random_sim_dir); make_folder(random_sim_dir)

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
                run_cmd("%s > %s 2>&1"%(randomSV_cmd, std_rearranging_genome))

                # edit the translocations so that the balanced ones are sorted
                translocations_file = "%s/translocations.tab"%random_sim_dir
                if file_is_empty(translocations_file): open(translocations_file, "w").write("\t".join(["Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB", "Balanced", "BpSeqA", "BpSeqB"])) # this needs to be 

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
        run_cmd("%s > %s 2>&1"%(targetSV_cmd, std_rearranging_genome))

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

def get_int_or_float_as_text(number):

    """Formats numbers"""

    if int(number)==float(number): return "%i"%number
    else: 
        float_parts = ("%.1f"%number).split(".")

        if float_parts[0]=="0": return ".%s"%float_parts[1]
        else: return "%.1f"%number

def convert_fasta_to_fqgz(fasta_file, replace=False, remove_fasta=True):

    """Takes a fasta file and converts it to fq"""

    # define fastq name
    prefix = fasta_file.rstrip(".fasta").rstrip(".fa")
    fastqgz = "%s.fq.gz"%prefix
    fastqgz_tmp = "%s.tmp.fq.gz"%prefix

    # convert into fastq and gzip
    if file_is_empty(fastqgz) or replace is True:
        print("generating %s"%fastqgz)

        # convert
        print("running reformat")
        run_cmd("%s in=%s out=%s qfake=50 overwrite=true"%(bbmap_reformat_sh, fasta_file, fastqgz_tmp))

        # remove the fasta
        if remove_fasta is True: os.unlink(fasta_file)

        os.rename(fastqgz_tmp, fastqgz)

    return fastqgz

def get_mosdepth_coverage_per_windows_output_likeBamStats(fileprefix, sorted_bam, windows_bed, replace=False, extra_threads=0, chromosome_id=""):

    """This function uses mosdepth to get the coverage for some regions in bed for a sorted_bam """

    # get outfiles
    fileprefix_tmp = "%s.tmp"%fileprefix
    regions_file = "%s.regions.bed.gz"%fileprefix 
    regions_file_tmp = "%s.regions.bed.gz"%fileprefix_tmp 
    thresholds_file = "%s.thresholds.bed.gz"%fileprefix 
    thresholds_file_tmp = "%s.thresholds.bed.gz"%fileprefix_tmp

    if file_is_empty(regions_file) or file_is_empty(thresholds_file) or replace is True:

        print("running mosdepth into %s"%fileprefix)
        
        # change the end, setting it to -1, and also sorting
        windows_1_based = "%s.1_based.bed"%windows_bed
        run_cmd(""" awk '{print $1 "\t" ($2+1) "\t" ($3)}' %s | sort -k1,1 -k2,2n > %s"""%(windows_bed, windows_1_based))

        # get the cmd
        cmd = "%s --threads %i --by %s --no-per-base --fast-mode --thresholds 1 --use-median %s %s"%(mosdepth, extra_threads, windows_1_based, fileprefix_tmp, sorted_bam) # mosdepth does not look at internal cigar operations or correct mate overlaps (recommended for most use-cases). It is also faster

        # add the chromosome_id if provided
        if chromosome_id!="": cmd = cmd.replace("--use-median", "--use-median --chrom %s"%chromosome_id)

        # run 
        run_cmd(cmd)

        # remove the 1-based file
        remove_file(windows_1_based)
 
        # remove innecessary files
        for sufix in ["mosdepth.global.dist.txt", "mosdepth.region.dist.txt", "mosdepth.summary.txt", "thresholds.bed.gz.csi", "regions.bed.gz.csi"]: remove_file("%s.%s"%(fileprefix_tmp, sufix))

        # keep
        os.rename(regions_file_tmp, regions_file)
        os.rename(thresholds_file_tmp, thresholds_file)

    print("arranging mosdepth output into df")

    # get as dfs
    df_regions = pd.read_csv(regions_file, sep="\t", header=-1, names=["#chrom",  "start", "end", "mediancov_1"]).drop_duplicates(subset=["#chrom",  "start", "end"])
    df_thresholds = pd.read_csv(thresholds_file, sep="\t").drop_duplicates(subset=["#chrom",  "start", "end"])

    # add the number of basepairs in each region that are covered by at least one
    try: df = df_regions.merge(df_thresholds, on=["#chrom",  "start", "end"], validate="one_to_one").rename(columns={"1X":"nbp_more_than_1x"})
    except:
        print("!!!!!!!!regions:%s, \n thresholds:%s, \n prefix:%s"%(regions_file, thresholds_file, fileprefix))
        raise ValueError("There was an error with joining the mosdepth outputs")

    def convert_to_float(x):
        
        """Takes a number and converts to float, and 0 if NaN"""

        if pd.isna(x): return 0.0
        else: return float(x)


    # add some trivial info
    df["length"] = df.end - df.start
    df["nocoveragebp_1"] = df.length - df.nbp_more_than_1x
    df["percentcovered_1"] = 100 - ((df.nocoveragebp_1/df.length) * 100).apply(convert_to_float)

    # get as 0-based
    df["start"] = df.start - 1

    # delete unnecessary files
    for f in [regions_file, thresholds_file]: remove_file(f)

    return df

def get_coverage_per_window_for_chromosomeDF(chromosome_id, destination_dir, windows_bed, sorted_bam, replace, window_l):

    """Takes a chromosome id, a destination dir where to write files, a windows file (provided by generate_coverage_per_window_file_parallel) and a sorted bam and it generates a dataframe with the coverage stats"""

    # define the output coverage file
    print("running coverage calculation for %s"%chromosome_id)
        
    # generate a randomID
    randID = id_generator(25)

    # define a file for the coverage
    windows_bed_chromsome = "%s.%s.%s.bed"%(windows_bed, chromosome_id, randID)
    run_cmd("grep $'%s\t' %s > %s"%(chromosome_id, windows_bed, windows_bed_chromsome))

    # if there is nothing, return an empty df
    bamstats_fields = ["#chrom", "start", "end", "length", "mediancov_1", "nocoveragebp_1", "percentcovered_1"]
    if file_is_empty(windows_bed_chromsome): return pd.DataFrame(columns=bamstats_fields)

    # calculate extra threads
    extra_threads = multiproc.cpu_count() - 1

    # define a file prefix on which to calculate the coverage
    mosdepth_outprefix = "%s.mosdepth_output"%windows_bed_chromsome

    # get a df that has all the
    df_coverage =  get_mosdepth_coverage_per_windows_output_likeBamStats(mosdepth_outprefix, sorted_bam, windows_bed_chromsome, replace=replace, extra_threads=extra_threads, chromosome_id=chromosome_id)

    remove_file(windows_bed_chromsome)

    return df_coverage[bamstats_fields]



def generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file="none", replace=False, window_l=1000, run_in_parallel=True, delete_bams=True):

    """Takes a reference genome and a sorted bam and runs a calculation of coverage per window (with bamstats04_jar) in parallel for sorted_bam, writing results under ddestination_dir. if window_file is provided then it is used. If not, it generates a file with non overlappping windows of length window_l"""

    # in the case that you have provided a window file
    if windows_file=="none":

        make_folder(destination_dir)

        # first generate the windows file
        windows_file = "%s.windows%ibp.bed"%(reference_genome, window_l)
        run_cmd("%s makewindows -g %s.fai -w %i > %s"%(bedtools, reference_genome, window_l, windows_file))

        # define the file
        coverage_file = "%s/coverage_windows_%ibp.tab"%(destination_dir, window_l)

        # define the chromosomes
        all_chromosome_IDs = [seq.id for seq in SeqIO.parse(reference_genome, "fasta")]

    # in the case you have provied a window file
    elif not file_is_empty(windows_file):

        # rename windows_file so that it is sorted
        pd.read_csv(windows_file, sep="\t", header=None, names=["chromosome", "start", "end"]).sort_values(by=["chromosome", "start", "end"]).to_csv(windows_file, header=False, index=False, sep="\t")

        # create the coverage file
        coverage_file = "%s.coverage_provided_windows.tab"%windows_file

        # define the destination dir personalized for the given windows file
        destination_dir = "%s.coverage_measurement_destination"%windows_file; make_folder(destination_dir)

        # define the chromosomes
        all_chromosome_IDs = sorted(set(pd.read_csv(windows_file, sep="\t", header=None, names=["chromosome", "start", "end"])["chromosome"]))

        # debug the fact that there is nothing to analyze
        if len(all_chromosome_IDs)==0: raise ValueError("There are no chromosomes in %s, so that no coverage can be calculated"%windows_file)

    else: raise ValueError("The provided windows_file %s does not exist"%windows_file)

    # generate the coverage file in parallel    
    if file_is_empty(coverage_file) or replace is True:

        # get the chromosomal dfs
        inputs_run = [(ID, destination_dir, windows_file, sorted_bam, replace, window_l) for ID in all_chromosome_IDs]

        if run_in_parallel is True:
        #if False is True: # DEBUG!!! ALWAYS NOT RUN IN PARALLEL # never running in parallel
            
            try:

                # initialize the pool class with the available CPUs --> this is syncronous parallelization
                pool = multiproc.Pool(multiproc.cpu_count())

                # run in parallel the coverage generation, which returns a list of dataframes, each with one chromosome
                chromosomal_dfs = pool.starmap(get_coverage_per_window_for_chromosomeDF, inputs_run)

                # close the pool
                pool.close(); pool.terminate(); pool.join()
                
            except KeyboardInterrupt:
                
                pool.close(); pool.terminate(); pool.join()
                raise ValueError("Keyboard Interrupt")

        else: chromosomal_dfs = list(map(lambda x: get_coverage_per_window_for_chromosomeDF(x[0], x[1], x[2], x[3], x[4], x[5]), inputs_run))

        # merge the dfs
        all_df = pd.DataFrame()
        for df in chromosomal_dfs: all_df = all_df.append(df, sort=True)

        # remove chromosomal files:
        for ID in all_chromosome_IDs: remove_file("%s/%s_coverage_windows%ibp.tab"%(destination_dir, ID, window_l))

        # check that it is not empty
        if len(all_df)==0: raise ValueError("There is no proper coverage calculation for %s on windows %s"%(sorted_bam, windows_file))

        # write
        all_df.to_csv(coverage_file, sep="\t", header=True, index=False)

        # at the end remove all the bam files # at some point I commented the lines below and I don't know why
        """
        if delete_bams is True:
            print("removing chromosomal bamfiles")

            for chrom in all_chromosome_IDs: 
                sorted_bam_chr = "%s.%s.bam"%(sorted_bam, chrom)
                os.unlink(sorted_bam_chr); os.unlink("%s.bai"%sorted_bam_chr)
        """

    return coverage_file

def plot_coverage_across_genome_pairedEndReads(sorted_bam, reference_genome, window_l=10000, replace=False):

    """Takes a sorted_bam and plots the coverage for windows of the genome"""

    print("plotting coverage across genome")

    # get coverage df  
    calculate_coverage_dir = "%s.calculating_windowcoverage"%sorted_bam; make_folder(calculate_coverage_dir)
    coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, calculate_coverage_dir, sorted_bam, replace=replace, window_l=window_l), "\t")
    all_chromosomes = sorted(set(coverage_df["#chrom"]))

    # plot, for each chromosome, the coverage
    fig = plt.figure(figsize=(6*len(all_chromosomes), 8)); I=1
    for yfield in ["mediancov_1", "percentcovered_1"]:
        for chromosome in all_chromosomes:

            ax = plt.subplot(2, len(all_chromosomes), I); I+=1
            df_c = coverage_df[coverage_df["#chrom"]==chromosome]
            sns.lineplot(x="start", y=yfield, data=df_c)

            ax.set_title(chromosome)


    #fig.tight_layout()  # otherwise the right y-label is slightly 
    filename="%s.coverage.pdf"%(sorted_bam)
    fig.savefig(filename, bbox_inches='tight');
    #if is_cluster is False: plt.show()
    plt.close(fig)



def simulate_pairedEndReads_per_chromosome_uniform(chr_obj, coverage, insert_size, read_lengths, max_fraction_chromosome=0.1):

    """This function takes a chromosome object (SeqRecord) and it generates reads that are as long as read_lengths, in a window that increases. It returns a list of objects, each of which has a chromosomeID_readID"""

    # define a function that takes a start, a length and an insert size and returns formated reads
    def get_paired_reads_list(chr_obj, startRead, len_read, insert_size, extraID):

        # define the end of the read
        endRead = startRead + len_read

        # define the coordinates of the pair
        startPair = endRead + insert_size
        endPair = startPair + len_read

        if endRead>len_chromosome or endPair>len_chromosome: raise ValueError("The calculation of read ends was not correct")

        # define ID of the read
        ID = "%s_readStart%i_%s"%(chr_obj.id, startRead, extraID)

        # keep the read and the pair
        read = chr_obj[startRead:endRead]
        read.id = "%s/1"%ID; read.name = ""; read.description = ""

        pair = chr_obj[startPair:endPair].reverse_complement()
        pair.id = "%s/2"%ID; pair.name = ""; pair.description = ""

        return [read, pair]

    # define overall metrics
    len_chromosome = len(chr_obj)
    read_lengths = [int(x) for x in sorted(read_lengths) if x<=int(len_chromosome*max_fraction_chromosome)]
    per_readLength_coverage = int((coverage/len(read_lengths))/2) + 1  # coverage that has to be achieved by each length of reads

    # initialize the reads
    reads_1 = []
    reads_2 = []

    # go through each read_length
    for len_read in read_lengths:

        # define the size of the increasing window that has to be reached to generate the per_readLength_coverage
        len_window_increment = int(len_read/per_readLength_coverage)

        # debug 
        if len_read<len_window_increment: raise ValueError("the read length has to be longer that the window increment")

        # go through each window, considering the pair
        startRead = 0
        while (startRead + 2*len_read + insert_size)<len_chromosome:

            # get the read pairs
            reads = get_paired_reads_list(chr_obj, startRead, len_read, insert_size, extraID="readLen%i_intrachromosomal_read"%len_read)
            reads_1.append(reads[0])
            reads_2.append(reads[1])

            # increment for the next window
            startRead += len_window_increment

        # add some reads at the end and the start
        for I in range(per_readLength_coverage):

            # get the reads for the first read
            reads = get_paired_reads_list(chr_obj, 0, len_read, insert_size, extraID="readLen%i_firstRead%i"%(len_read, I))
            reads_1.append(reads[0])
            reads_2.append(reads[1])

            # define coordinates of the last read
            start_lastRead = len_chromosome - (2*len_read + insert_size)

            # keep 
            reads = get_paired_reads_list(chr_obj, start_lastRead, len_read, insert_size, extraID="readLen%i_lastRead%i"%(len_read, I))
            reads_1.append(reads[0])
            reads_2.append(reads[1])



    return (reads_1, reads_2)


###################################################################################################
################################## RUN GRIDSS AND CLOVE PIPELINE ##################################
###################################################################################################

def run_gridss_and_annotateSimpleType(sorted_bam, reference_genome, outdir, replace=False, threads=4, blacklisted_regions="", maxcoverage=50000, max_threads=16):

    """Runs gridss for the sorted_bam in outdir, returning the output vcf. blacklisted_regions is a bed with regions to blacklist"""

    # define the output
    gridss_VCFoutput = "%s/gridss_output.vcf"%outdir; 
    if gridss_VCFoutput.split(".")[-1]!="vcf": raise ValueError("gridss needs a .vcf file. this is not the case for"%gridss_VCFoutput)

    if file_is_empty(gridss_VCFoutput) or replace is True:

        # define other files
        gridss_assemblyBAM = "%s/gridss_assembly.bam"%outdir
        gridss_tmpdir = "%s/gridss_tmp"%outdir

        # if the blacklisted_regions does not exist, just create an empty file
        if file_is_empty(blacklisted_regions): 

            blacklisted_regions = "%s/empty_regions.bed"%outdir; open(blacklisted_regions, "w").write("")

        print("blacklisting %s\n"%blacklisted_regions)
        
        # change the number of threads if more than max_threads, which is the optimum for gridss (8 is the optimum)
        if threads>max_threads: threads =  max_threads # this is to optimise for the reccommended level of parallelism

        # define the out and error of gridss
        #gridss_std = "%s/gridss_run_std.txt"%outdir
        gridss_std = "stdout"
        
        max_tries = 2
        for Itry in range(max_tries):
            print("running gridss try %i"%(Itry+1))
            try: 
                # delete previous files
                delete_folder(gridss_tmpdir); make_folder(gridss_tmpdir)
                remove_file(gridss_assemblyBAM)

                # define the heap size, which depends on the cloud or not
                #jvmheap = "27.5g" # this is the default
                jvmheap = "20g" # this works in MN. This can be changed to fit the machine

                # run
                print("running gridss on %s jvmheap"%jvmheap)

                gridss_cmd = "%s --jar %s --reference %s -o %s --assembly %s --threads %i --workingdir %s --maxcoverage %i --blacklist %s --jvmheap %s %s"%(gridss_run, gridss_jar, reference_genome, gridss_VCFoutput, gridss_assemblyBAM, threads, gridss_tmpdir, maxcoverage, blacklisted_regions, jvmheap, sorted_bam)
                if gridss_std!="stdout": gridss_cmd += " > %s 2>&1"%gridss_std
                run_cmd(gridss_cmd)

                break

            except: print("WARNING: GRIDSS failed on try %i"%(Itry+1))

        # if it did not finish correctly finish it
        if file_is_empty(gridss_VCFoutput): raise ValueError("gridss did not finish correctly after %i tries. Check the log in %s and the std in %s"%(max_tries, gridss_tmpdir, gridss_std))

        delete_folder(gridss_tmpdir); remove_file(gridss_assemblyBAM)

    # annotated simple events, the ones that can be predicted from each breakpoint, but only if there are some predicted events
    gridss_VCFoutput_with_simple_event = "%s.withSimpleEventType.vcf"%gridss_VCFoutput
    simple_event_std = "%s/simple_event_annotation.std"%outdir
    n_breakends = len([l for l in open(gridss_VCFoutput, "r").readlines() if not l.startswith("#")])
    if (file_is_empty(gridss_VCFoutput_with_simple_event) or replace is True) and n_breakends>0 : run_cmd("%s %s > %s 2>&1"%(annotate_simpleEvents_gridssVCF_R, gridss_VCFoutput, simple_event_std))

    return gridss_VCFoutput_with_simple_event

def load_single_sample_VCF(path):

    """Loads a vcf with a single sample into a df."""

    # load the df
    df = pd.read_csv(path, sep="\t", header = len([l for l in open(path, "r") if l.startswith("##")]))

    # change the name of the last column, which is the actual sample data
    df = df.rename(columns={df.keys()[-1]: "DATA"})

    # get the filter as set and filter
    df["FILTER_set"] = df.FILTER.apply(lambda x: set(x.split(";")))
    #df = df[df.FILTER_set.apply(lambda x: len(x.intersection(interesting_filterTAGs))>0)]

    # check that there

    ### INFO COLUMN ####

    def get_real_value(string):

        try: return ast.literal_eval(string)
        except: return string

    def get_dict_fromRow(row):

        # get as a list, considering to put appart whatever has an "="
        list_row = row.split(";")
        rows_with_equal = [x.split("=") for x in list_row if "=" in x]
        rows_without_equal = [x for x in list_row if "=" not in x]

        # add the values with an "=" to the dictionary
        final_dict = {"INFO_%s"%x[0] : get_real_value(x[1]) for x in rows_with_equal}

        # add the others collapsed
        final_dict["INFO_misc"] = ";".join(rows_without_equal)

        return final_dict

    # add a column that has a dictionary with the info fields
    df["INFO_as_dict"] = df.INFO.apply(get_dict_fromRow)
    all_INFO_fields = sorted(list(set.union(*df.INFO_as_dict.apply(lambda x: set(x)))))

    # add them as sepparated columns
    def get_from_dict_orNaN(value, dictionary):

        if value in dictionary: return dictionary[value]
        else: return np.nan

    for f in all_INFO_fields: df[f] = df.INFO_as_dict.apply(lambda d: get_from_dict_orNaN(f, d))
    df.pop("INFO_as_dict")

    #########################

    ### FORMAT COLUMN ###

    # check that there is only one format
    all_formats = set(df.FORMAT)
    if len(all_formats)!=1: raise ValueError("There should be only one format in the FORMAT file")

    # get as dictionary
    format_list = next(iter(all_formats)).split(":")

    def getINTorFLOATdictionary(data):

        # get dictionary
        dictionary = dict(zip(format_list, data.split(":")))

        # change the type
        final_dict = {}
        for k, v in dictionary.items():

            if v==".": final_dict[k] = "."
            elif "/" in v or "," in v: final_dict[k] = v
            elif "." in v: final_dict[k] = float(v)
            else: final_dict[k] = int(v)

        return final_dict

    df["FORMAT_as_dict"] = df.DATA.apply(getINTorFLOATdictionary)

    # get as independent fields
    for f in sorted(format_list): df["DATA_%s"%f] = df.FORMAT_as_dict.apply(lambda d: d[f])
    df.pop("FORMAT_as_dict")

    #########################

    # calculate allele frequencies (if AD is provided)
    if "DATA_AD" in df.keys():

        df["allele_freqs"] = df.apply(lambda r: [int(reads_allele)/r["DATA_DP"] for reads_allele in r["DATA_AD"].split(",")], axis=1)

        # assign a boolean whether it is heterozygius
        df["is_heterozygous_coverage"] = df.allele_freqs.apply(lambda x: any([freq>=0.25 and freq<=0.75 for freq in x]))
        df["is_heterozygous_GT"] = df.DATA_GT.apply(lambda x: len(set(x.split("/")))>1)

    return df

def getNaN_to_0(x):

    if pd.isna(x): return 0.0
    else: return x

def add_info_to_gridssDF(df, expected_fields={"allele_frequency", "allele_frequency_SmallEvent", "other_coordinates", "other_chromosome", "other_position", "other_orientation",  "inserted_sequence", "len_inserted_sequence", "length_event", "has_poly16GC", "length_inexactHomology", "length_microHomology"}, median_insert_size=500, median_insert_size_sd=50):

    """This function takes a gridss df and returns the same adding expected_fields"""

    #print("adding info to gridss df")
    if len(set(df.keys()).intersection(expected_fields))!=len(expected_fields):

        df = cp.deepcopy(df)

        # add the allele frequencies
        #print("adding allele freq")
        df["allele_frequency"] = df.apply(lambda r: np.divide(r["DATA_VF"] , (r["DATA_VF"] + r["DATA_REF"] + r["DATA_REFPAIR"])), axis=1).apply(getNaN_to_0)
        df["allele_frequency_SmallEvent"] = df.apply(lambda r: np.divide(r["DATA_VF"] , (r["DATA_VF"] + r["DATA_REF"])), axis=1).apply(getNaN_to_0)

        # add data to calculate the other breakpoint
        #print("get position of other bend")
        def get_other_position(r):
            
            if "." in r["ALT"]: return "%s:%i"%(r["#CHROM"], r["POS"])
            else: return re.split("\]|\[", r["ALT"])[1]

        df["other_coordinates"] = df.apply(get_other_position, axis=1)
        df["other_chromosome"] = df.other_coordinates.apply(lambda x: x.split(":")[0])
        df["other_position"] = df.other_coordinates.apply(lambda x: int(x.split(":")[1]))

        # add the orientation of the other
        def get_other_orientation(alt): 

            if "]" in alt: return "]"
            elif "[" in alt: return "["
            else: return "unk"

        df["other_orientation"] = df.ALT.apply(get_other_orientation)

        # get the inserted sequence
        def get_inserted_seq(ALT):
            
            if "." in ALT: return ""
            else: return [x for x in re.split("\]|\[", ALT) if all([base.upper() in {"A", "C", "T", "G", "N"} and len(base)>0 for base in set(x)])][0]
            
        df["inserted_sequence"] = df.ALT.apply(get_inserted_seq)
        df["len_inserted_sequence"] = df.inserted_sequence.apply(len)

        def get_length(r):

            
            # same chromosome
            if r["#CHROM"]==r["other_chromosome"]: return abs(r["other_position"] - r["POS"])
            
            # different chromosomes
            return 100000000000

        df["length_event"] = df[["#CHROM", "other_chromosome", "POS", "other_position"]].apply(get_length, axis=1)
        
    
        # add if the df has a polyG tag
        df["has_poly16GC"] = df.ALT.apply(lambda x: ("G"*16) in x or ("C"*16) in x)

        # add the range of inexact homology
        def add_inexact_homology_length(IHOMPOS):
            
            if pd.isna(IHOMPOS): return  0
            else: return int(IHOMPOS.split(",")[1]) - int(IHOMPOS.split(",")[0])

        df["length_inexactHomology"] = df.INFO_IHOMPOS.apply(add_inexact_homology_length)

        # add the length of homology
        def get_homology_length(HOMLEN):
            
            if pd.isna(HOMLEN): return 0
            else: return int(HOMLEN)

        df["length_microHomology"] = df.INFO_HOMLEN.apply(get_homology_length)

        # add the actual allele_frequency, which depends on the median_insert_size and the median_insert_size_sd. If the breakend is longer than the insert size it is a large event
        maxiumum_insert_size = median_insert_size + median_insert_size_sd
        def get_allele_freq(r):
            
            if r["length_event"]>maxiumum_insert_size: return r["allele_frequency"]
            else: return r["allele_frequency_SmallEvent"]
        
        df["real_AF"] = df.apply(get_allele_freq, axis=1)

        if any(pd.isna(df["real_AF"])): raise ValueError("There are NaNs in the real_AF field")


    return df

def get_gridssDF_filtered(df, min_Nfragments=8, min_af=0.005, wrong_INFOtags={"IMPRECISE"}, wrong_FILTERtags={"NO_ASSEMBLY"}, filter_polyGC=True, filter_noSplitReads=True, filter_noReadPairs=True, maximum_strand_bias=0.95, maximum_microhomology=50, maximum_lenght_inexactHomology=50, range_filt_DEL_breakpoints=[100, 800], min_length_inversions=40, dif_between_insert_and_del=5, max_to_be_considered_small_event=1000, min_size=50, add_columns=True, min_af_EitherSmallOrLargeEvent=0.0 ):

    """Takes a gridss df (each line is a breakend and it has the field INFO_SIMPLE_TYPE) and it returns the same df bt with a personalFILTER field, according to all the passed filters .
    Below are the filters with the meaning:

    max_length_SmallDupDel = 1000 # the maximum length to be considered a small duplication or deletion
    min_Nfragments = 8 # the minimum number of reads that should support the variant
    min_af = 0.005 # the minimum allele frequency for a var to be considered. Heterozygous variants will have much less tan 0.5 frequency
    wrong_INFOtags = {"IMPRECISE"} # the tags in miscINFO to filter out
    wrong_FILTERtags = {"NO_ASSEMBLY"} # the tags in FILTER to avoid
    filter_polyGC = True # eliminate breakpoints that have a long G or C tract
    filter_noSplitReads = True # eliminate breakpoints that have no split reads supporting them
    filter_noReadPairs = True # eliminate brreakpints that do not have any discordant read pairs
    maximum_strand_bias = 0.95 # all vars with a SB above this will be removed
    maximum_microhomology = 50 # all vars with a microhomology above this are removed
    maximum_lenght_inexactHomology = 50 # all vars that are longer than min_len and have an indexact homology above this number will be filtered
    range_filt_DEL_breakpoints = [100, 800] # remove all variants that are DEL with a length within this range if their inexact microhomology is >= 6bp
    min_length_inversions= 40 # remove all INV variants that are below this length and have >= 6bp of microhomology
    dif_between_insert_and_del = 5 # any small deletion that has an insertion sequence longer than length_event-dif_between_insert_and_del will be removed
    max_to_be_considered_small_event is the maxiumum length that small events is considered
    min_size is the minimum length that a breakend should have to be considered
    min_af_EitherSmallOrLargeEvent is the minimum allele frequency regardless of if it is a small or a large event, which depends on the insert size that is not always constant. The idea is that any breakend with any AF (calculated as VF/VF+REF+REFPAIR or VF/VF+REF) above min_af_EitherSmallOrLargeEvent will be kept

    add_columns indicates whether to add columns to the df, which may have been done before

    The default values are the ones used in the gridss-purple-linx pipeline to generate the somatic callset

    """

    ######## ADD COLUMNS TO THE DF FOR FURTHER CALCULATION ##########
    if add_columns is True: df = add_info_to_gridssDF(df)

    # define whether the variant is a small duplication or insertion. These have special filters
    df["is_small_DupDel"] = (df.INFO_SIMPLE_TYPE.isin({"DEL", "DUP"})) & (df.length_event<=max_to_be_considered_small_event)

    # get sets
    wrong_INFOtags = set(wrong_INFOtags)
    wrong_FILTERtags = set(wrong_FILTERtags)

    ############ APPLY THE FILTERS ###########

    idx = ((df.length_event>=min_size) &
       (df.DATA_VF>=min_Nfragments) &  
       (df.real_AF>=min_af) & 
       ( (df.allele_frequency>=min_af_EitherSmallOrLargeEvent) | (df.allele_frequency_SmallEvent>=min_af_EitherSmallOrLargeEvent) ) & 
       ~(df.INFO_misc.isin(wrong_INFOtags)) & 
       (df.FILTER.apply(lambda f: len(set(f.split(";")).intersection(wrong_FILTERtags))==0) ) & 
       ( ((df.INFO_SB.apply(float)<maximum_strand_bias) | pd.isna(df.INFO_SB.apply(float))) | ~(df.is_small_DupDel) ) &
       (df.length_microHomology<maximum_microhomology) &
       ((df.length_inexactHomology<maximum_lenght_inexactHomology) | (df.is_small_DupDel)) &
       ~((df.INFO_SIMPLE_TYPE=="DEL") & (df.length_inexactHomology>=6) & (df.length_event>=range_filt_DEL_breakpoints[0]) & (df.length_event<=range_filt_DEL_breakpoints[1])) & 
       ~((df.INFO_SIMPLE_TYPE=="INV") & (df.length_event<min_length_inversions) & (df.length_microHomology>=6)) &       
       ~((df.INFO_SIMPLE_TYPE=="DEL") & (df.length_event<=max_to_be_considered_small_event) & (df.len_inserted_sequence>=(df.length_event-dif_between_insert_and_del)) )       
       )

    if filter_polyGC: idx = idx & ~(df.has_poly16GC)
    if filter_noSplitReads: idx = idx & (~(df.DATA_SR==0) | ~(df.is_small_DupDel))
    if filter_noReadPairs: idx = idx & (~(df.DATA_RP==0) | (df.is_small_DupDel))

    # return the filtered df

    return df[idx]

def get_gridssDF_filtered_from_filtersDict(df_gridss, filters_dict):

    """Takes a df gridss and returns the filtered one, according to filters_dict"""

    # debug the fact that there is no min_af_EitherSmallOrLargeEvent
    if "min_af_EitherSmallOrLargeEvent" not in filters_dict: filters_dict["min_af_EitherSmallOrLargeEvent"] = 0.0

    # get the filtered df
    df_filt = get_gridssDF_filtered(df_gridss, min_Nfragments=filters_dict["min_Nfragments"], min_af=filters_dict["min_af"], wrong_INFOtags=filters_dict["wrong_INFOtags"], wrong_FILTERtags=filters_dict["wrong_FILTERtags"], filter_polyGC=filters_dict["filter_polyGC"], filter_noSplitReads=filters_dict["filter_noSplitReads"], filter_noReadPairs=filters_dict["filter_noReadPairs"], maximum_strand_bias=filters_dict["maximum_strand_bias"], maximum_microhomology=filters_dict["maximum_microhomology"], maximum_lenght_inexactHomology=filters_dict["maximum_lenght_inexactHomology"], range_filt_DEL_breakpoints=filters_dict["range_filt_DEL_breakpoints"], min_length_inversions=filters_dict["min_length_inversions"], dif_between_insert_and_del=filters_dict["dif_between_insert_and_del"], max_to_be_considered_small_event=filters_dict["max_to_be_considered_small_event"], min_size=filters_dict["min_size"], add_columns=False, min_af_EitherSmallOrLargeEvent=filters_dict["min_af_EitherSmallOrLargeEvent"] )

    return df_filt



def get_bedpe_from_svVCF(svVCF, outdir, replace=False, only_simple_conversion=False):

    """Takes a svVCF and writes a bedpe file. outdir is where to write files"""

    # generate the bedpe file 
    bedpe_file = "%s.bedpe"%svVCF

    if file_is_empty(bedpe_file) or replace is True:

        remove_file(bedpe_file)

        # write one file or another if the file has no records
        len_vcf_records = len([l for l in open(svVCF, "r").readlines() if not l.startswith("#")])

        if len_vcf_records>0:

            if only_simple_conversion is True:

                print("Getting files for svVCF file. Mainly generating a bedpe file for breakpoints with some extra info, but simple.")
                r_stdout = "%s/svVCF_analysis_log.out"%outdir
                run_cmd("%s %s > %s 2>&1"%(analyze_svVCF_simple, svVCF, r_stdout))                

            else:
                print("Getting files for svVCF file. Mainly generating a bedpe file for breakpoints with some extra info.")
                r_stdout = "%s/svVCF_analysis_log.out"%outdir
                run_cmd("%s %s > %s 2>&1"%(analyze_svVCF, svVCF, r_stdout))

        else: open(bedpe_file, "w").write("no_vcf_records\n")

    return bedpe_file


def get_bedpe_df_with_added_feats(df_bedpe, df_gridss):

    """This function takes a bedpe file and a df gridss and returns de bedpe as df with the real_af of each of the breakends"""

    # add the eventID to the gridss df
    df_gridss["eventID"] = df_gridss.INFO_EVENT.apply(lambda x: x+"o")

    # map each eventID to the chromosome and the af
    eventID_to_chrom_to_pos_to_af = df_gridss.groupby("eventID").apply(lambda df_e: df_e.set_index(["#CHROM", "POS"])["real_AF"])

    # map each event to an af
    df_bedpe["af_1"] = df_bedpe.apply(lambda r: [af for pos,af in dict(eventID_to_chrom_to_pos_to_af.loc[(r["name"], r["chrom1"])]).items() if pos>=r["start1"] and pos<=r["end1"]][0], axis=1)
    df_bedpe["af_2"] = df_bedpe.apply(lambda r: [af for pos,af in dict(eventID_to_chrom_to_pos_to_af.loc[(r["name"], r["chrom2"])]).items() if pos>=r["start2"] and pos<=r["end2"]][0], axis=1)

    # get the ploidy
    df_bedpe["bp1_heterozygous"] = df_bedpe.af_1<0.75
    df_bedpe["bp2_heterozygous"] = df_bedpe.af_2<0.75

    ##### add the 5'->3' added regions #####
    def get_tuple_5_to_3_positions(r):

        """Takes a row of the bedpe and returns a tuple where there is the 5' position and the 3' attached. This depends on the orientation. See the clove paper to understand what the orientations mean. Only tandem duplictaions are avoided"""

        # get data
        chrom1 = r["chrom1"]
        pos1 = r["start1"]
        strand1 = r["strand1"]
        chrom2 = r["chrom2"]
        pos2 = r["start2"]
        strand2 = r["strand2"]

        # define for chroms that are the same
        if chrom1==chrom2:

            # first check that the second position is 3'
            if pos2<=pos1: raise ValueError("Your bedpe may not be properly formated as pos2<=pos1")

            # for TANDEM duplications just add the next position
            if strand1=="-" and strand2=="+": tuple_5_to_3_positions = [(chrom2, pos2), (chrom2, pos2+1)]

            # INV1, DEL and INV2 are the same
            else: tuple_5_to_3_positions = [(chrom1, pos1), (chrom2, pos2)]

        # define for chroms that are different
        else: tuple_5_to_3_positions = [(chrom1, pos1), (chrom2, pos2)]

        return tuple_5_to_3_positions

    
    df_bedpe["5_to_3_positions_tuple_withoutTAN"] = df_bedpe.apply(get_tuple_5_to_3_positions, axis=1)

    ########################################

    return df_bedpe


def get_genomeGraph_object(genome, df_bedpe, df_gridss_filt, genomeGraph_outfileprefix, replace=False):

    """This function takes a bedpe and generates an undirected graph where each position is mapped to the 3' position in the genome. Tandem duplications are not considered as they may get weird loops"""

    # define the files
    genomeGraph_outfile = "%s.graph.py"%genomeGraph_outfileprefix
    genomeGraph_positions_df = "%s.df_positions.py"%genomeGraph_outfileprefix

    if any([file_is_empty(x) for x in {genomeGraph_outfile, genomeGraph_positions_df}]) or replace is True:
    #if True: # debug

        # map each chromosome to an offset
        chrom_to_lenSeq = {seq.id : len(seq.seq) for seq in SeqIO.parse(genome, "fasta")}

        # deffine an offset for each chromosome, which is necessary to keep all the positions of the genome as independent numbers
        chrom_to_offset = {}
        current_offset = 0
        for chrom, seqLen in chrom_to_lenSeq.items():
            chrom_to_offset[chrom] = current_offset
            current_offset+=seqLen

        # create the graph
        genome_graph = igraph.Graph(directed=False)

        # add one vertex (node) for each position in the genome
        npositions = sum(chrom_to_lenSeq.values())
        genome_graph.add_vertices(npositions)

        # define the edges that are the end of chromosomes
        chromosome_start_nodes = {offset for chrom, offset in chrom_to_offset.items()}
        chromosome_end_nodes = {(offset + chrom_to_lenSeq[chrom] - 1) for chrom, offset in chrom_to_offset.items()}

        # define the connections (edges, parts of the genomes that are linked) from the breakpoint
        all_edges = [] 

        # if a df bedpe is provided, calculate the edges resulting from the breakpoints
        if df_bedpe is not None and len(df_bedpe)>0:

            # get the 5' to 3' positions and ploidy of the breakpoints into bedpe
            df_bedpe = get_bedpe_df_with_added_feats(df_bedpe, df_gridss_filt)
            print("getting graph-based genome")

            for (bp1_heterozygous, bp2_heterozygous, bp_positions_tuple) in df_bedpe[["bp1_heterozygous", "bp2_heterozygous", "5_to_3_positions_tuple_withoutTAN"]].values:

                # initialize a list of the existing edges in the form of [(chr1, pos1), (chr2, pos2)]
                existing_joined_positions = []

                # add the breakpoint edge
                existing_joined_positions.append(bp_positions_tuple)

                # get the positions of the breakpoint
                chrom1, pos1 = bp_positions_tuple[0]
                chrom2, pos2 = bp_positions_tuple[1]

                # add the map between each position and next one if it is heterozygous
                if bp1_heterozygous: existing_joined_positions.append([(chrom1, pos1), (chrom1, pos1+1)])
                if bp2_heterozygous: existing_joined_positions.append([(chrom2, pos2), (chrom2, pos2+1)])

                # add the existing mapped positions
                all_edges += [(chrom_to_offset[pos1[0]]+pos1[1], chrom_to_offset[pos2[0]]+pos2[1]) for pos1, pos2 in existing_joined_positions]

        # add the non-mapped positions to the 3' ones
        positions_already_added = {edge[0] for edge in all_edges}
        all_positions = set(range(npositions))
        remaining_positions = all_positions.difference(positions_already_added.union(chromosome_end_nodes))
        all_edges += [(pos, pos+1) for pos in remaining_positions]

        # add the edges to the graph
        genome_graph.add_edges(all_edges)
        print("genome graph got")

        # get the real ends of the chromosomes
        sorted_positions = sorted(all_positions)
        pos_to_nNeighbors = pd.Series(dict(zip(sorted_positions, map(lambda x: len(genome_graph.neighbors(x, mode="ALL")), sorted_positions))))

        # debug
        if any(pos_to_nNeighbors<1): raise ValueError("there are some unnconected nodes in the graph genome")
        if any(pos_to_nNeighbors>100000): raise ValueError("there are some very highly connected regions in the graph genome")
        if any(pd.isna(pos_to_nNeighbors)): raise ValueError("there are some NaNs in the graph genome")

        real_chromosome_end_nodes = set(pos_to_nNeighbors[pos_to_nNeighbors==1].index)
        print("There are %i telomeric nodes in the graph genome"%len(real_chromosome_end_nodes))

        # generate a df that maps each position to the real position
        positions_real = []
        chromosomes_real = []
        for chrom, lenChrom in chrom_to_lenSeq.items():
            positions_real += list(range(lenChrom))
            chromosomes_real += [chrom]*lenChrom

        df_positions = pd.DataFrame()
        df_positions["chromosome"] =  chromosomes_real
        df_positions["real_position"] =  positions_real
        df_positions["offset"] = df_positions.chromosome.apply(lambda x: chrom_to_offset[x])
        df_positions["graph_position"] = df_positions.real_position + df_positions.offset
        df_positions["is_end_of_chr"] = df_positions.graph_position.isin(real_chromosome_end_nodes)

        if set(df_positions.graph_position)!=all_positions: raise ValueError("There is a bad graph calculation of the positions")

        # save
        save_object(genome_graph, genomeGraph_outfile)
        save_object(df_positions, genomeGraph_positions_df)

    else:
        print("loading graph genome")
        genome_graph = load_object(genomeGraph_outfile)
        df_positions = load_object(genomeGraph_positions_df)

    return genome_graph, df_positions


def run_clove_filtered_bedpe(bedpe_file, outfile, sorted_bam, replace=False, median_coverage=10, median_coverage_dev=1, check_coverage=True):

    """ Takes a bedpe file and a sorted bam and generates the clove output. df_cov is a df of the first 4 columns of the mpileup output. If sorted_bam and median_coverage are provided, there will be a check for TAN and DEL of coverage 

    every TAN or DEL outside median_coverage +- median_coverage_dev will be filtered out"""

    if file_is_empty(outfile) or replace is True:

        print("running clove")
        outfile_tmp = "%s.tmp"%outfile
        clove_std = "%s.std"%outfile
        cmd = "%s -jar %s -i %s BEDPE -b %s -o %s -c %i %i"%(JAVA, clove, bedpe_file, sorted_bam, outfile_tmp, median_coverage, median_coverage_dev)

        if check_coverage is False:

            print("avoid checking coverage")
            cmd += " -r "

        print("running clove with cmd:", cmd)


        # add the std
        cmd = "%s > %s 2>&1"%(cmd, clove_std)

        run_cmd(cmd)
        os.rename(outfile_tmp, outfile)

        # remove at the end
        remove_file(clove_std)

    #print("clove finsihed correctly")

def get_distance_to_telomere_series(df_chromosome_position, genome_graph, df_positions_graph):

    """This function takes a df with chromosome and position and returns the distance to the telomere according to the genome graph"""

    print("getting distance to the telomere")

    # rename the df to have chrom and pos
    df_chromosome_position = df_chromosome_position.rename(columns=dict(zip(df_chromosome_position.columns, ["chromosome", "position"])))

    # add the graph positions
    df_chromosome_position = df_chromosome_position.merge(df_positions_graph, left_on=["chromosome", "position"], right_on=["chromosome", "real_position"], how="left", validate="one_to_one")

    # define all the positions
    all_positions = sorted(set(df_chromosome_position.graph_position))

    # define the positions that are ends of chromosomes
    chrom_end_positions = sorted(set(df_positions_graph[df_positions_graph.is_end_of_chr].graph_position))

    # calculate the distance from each position to the 
    print("calculating shortest paths in genome graph from %i end positions to %i positions. This may take a lot"%(len(chrom_end_positions), len(all_positions)))
    shortestPath_lengths_df = pd.DataFrame(genome_graph.shortest_paths(source=chrom_end_positions, target=all_positions, mode="IN"), columns=all_positions, index=chrom_end_positions)
    distance_telomere_series = shortestPath_lengths_df.apply(min, axis=0).apply(int)
    print("distance calculated")

    #distance_telomere_series = pd.Series([1]*len(all_positions), index=all_positions)  # debug

    # reorder so that it fits the input df
    distance_telomere_series = pd.Series(distance_telomere_series.loc[df_chromosome_position.graph_position])
    distance_telomere_series.index = df_chromosome_position.index

    # return ordered as in df_chromosome_position
    return distance_telomere_series


def generate_nt_content_file(genome, target_nts="GC", replace=False):

    """Takes a genome and outputs a file with chromosome, position and 1 or 0 regarding if any of the target_nts is the same in the genome. This is 0-based"""

    target_nt_content_file = "%s.%scontent.tab"%(genome, target_nts)

    if file_is_empty(target_nt_content_file) or replace is True:
        print("getting GC content df per position")

        print("calculating %s content"%target_nts)

        # initialize a dict to create a genome df
        genome_dict = {"chromosome":[], "position":[], "base":[]}
        for seq in SeqIO.parse(genome, "fasta"): 

            # load with the sequence content at each position
            sequence_str = str(seq.seq).upper()
            genome_dict["chromosome"] += [seq.id]*len(sequence_str)
            genome_dict["position"] += list(range(0, len(sequence_str)))
            genome_dict["base"] += list(sequence_str)

        # get into df
        genome_df = pd.DataFrame(genome_dict)

        # add if the base is in target_nts
        target_nts_set = set(target_nts.upper())
        genome_df["is_in_%s"%target_nts] = genome_df["base"].isin(target_nts_set).apply(int)

        # write
        target_nt_content_file_tmp = "%s.tmp"%target_nt_content_file
        genome_df.to_csv(target_nt_content_file_tmp, sep="\t", header=True, index=False)
        os.rename(target_nt_content_file_tmp, target_nt_content_file)

    return target_nt_content_file

def get_df_with_GCcontent(df_windows, genome, gcontent_outfile, replace=False):

    """This function takes a df with windows of the genome and adds the gc content for each window, writing a file under gcontent_outfile"""

    print("Getting GC content")

    if file_is_empty(gcontent_outfile) or replace is True:

        # resort
        df_windows = df_windows.sort_values(by=["chromosome", "start", "end"])

        # get the GC content file for each position
        gc_content_outfile_perPosition = generate_nt_content_file(genome, replace=replace, target_nts="GC")
        gc_df = pd.read_csv(gc_content_outfile_perPosition, sep="\t")[["chromosome", "position", "is_in_GC"]].sort_values(by=["chromosome", "position"])

        # define a df where each position is one row and it has the start_window as an add
        df_windows["length"] = df_windows.end - df_windows.start
        positions = make_flat_listOflists(list(df_windows.apply(lambda r: list(range(r["start"], r["end"])), axis=1)))
        start_windows = make_flat_listOflists(list(df_windows.apply(lambda r: [r["start"]]*r["length"], axis=1)))
        chromosomes = make_flat_listOflists(list(df_windows.apply(lambda r: [r["chromosome"]]*r["length"], axis=1)))
        df_positions = pd.DataFrame({"position":positions, "chromosome":chromosomes, "start_window":start_windows})

        # add the positions to the gc df
        gc_df = gc_df.merge(df_positions, on=["chromosome", "position"], how="right")        

        # calculate the GC content and add to df
        startWindow_to_gc = gc_df[["start_window", "is_in_GC"]].groupby("start_window").mean()["is_in_GC"]
        df_windows["GCcontent"] = list(startWindow_to_gc.loc[df_windows.start])

        # at the end save
        save_object(df_windows, gcontent_outfile)

    else: df_windows = load_object(gcontent_outfile)


    return df_windows

def get_distanceToTelomere_chromosome_GCcontent_to_coverage_fn(df_coverage_train, genome, genome_graph, df_positions_graph, outdir, mitochondrial_chromosome="mito_C_glabrata_CBS138", replace=False):

    """This function takes a training df_coverage (with windows of a genome) and returns a lambda function that takes GC content, chromosome and  distance to the telomere and returns coverage according to the model.

    Your genome graph can also be a linear genome, you just have to create it without considering breakpoints"""
    print("getting coverage-predictor function")

    # rename the training df
    df = df_coverage_train.rename(columns={"#chrom":"chromosome", "mediancov_1":"coverage"})

    # add the distance to the telomere
    df_with_distance_to_telomere_file = "%s/df_with_distance_to_telomere_file.py"%outdir
    if file_is_empty(df_with_distance_to_telomere_file) or replace is True:

        df["middle_position"] = (df.start + (df.end - df.start)/2).apply(int)
        df["distance_to_telomere"] = get_distance_to_telomere_series(df[["chromosome", "middle_position"]], genome_graph, df_positions_graph)

        save_object(df, df_with_distance_to_telomere_file)

    else: df = load_object(df_with_distance_to_telomere_file)

    # add the gc content
    gcontent_outfile = "%s/GCcontent.py"%outdir
    df = get_df_with_GCcontent(df, genome, gcontent_outfile, replace=replace)

    # define the set of each type of chromosomes
    all_chromosomes = set(df.chromosome)
    if mitochondrial_chromosome!="no_mitochondria": mtDNA_chromosomes = set(mitochondrial_chromosome.split(","))
    else: mtDNA_chromosomes = set()
    gDNA_chromosomes = all_chromosomes.difference(mtDNA_chromosomes)

    ######## find the coeficients for each chromosome #########

    # map each chromosome to the coefs of the quadratic fit that explains coverage form the distance to the telomere and also the coefs of the GC content explaining the resiudal of this fit
    chrom_to_coefType_to_coefs = {}

    # go through each type of genome
    for type_genome, chroms in [("mtDNA", mtDNA_chromosomes), ("gDNA", gDNA_chromosomes)]:
        print("investigating %s"%type_genome)

        # define the training df
        df_g = df[df.chromosome.isin(chroms)]

        # define the relative coverage of each window of this genome
        df_g["relative_coverage"] = df_g.coverage / np.median(df_g.coverage)

        # go through each chrom and identify that it is duplicated if the quadratic fit from the prediction of the distance to the telomere suggests a minimum of >=1.6
        duplicated_chroms = set()
        for chrom in chroms:

            # get df of this chrom
            df_c = df_g[df_g.chromosome==chrom]

            # fit a quadratic fit
            coefs = poly.polyfit(df_c.distance_to_telomere, df_c.relative_coverage, 2)
            predicted_coverage = poly.polyval(df_c.distance_to_telomere, coefs)

            if min(predicted_coverage)>1.6: duplicated_chroms.add(chrom)

        # define the training set for the modelling
        df_correct = df_g[(df_g.relative_coverage<=4) & (df_g.relative_coverage>0.1) & ~(df_g.chromosome.isin(duplicated_chroms))]

        # now fit the model that predicts coverage from the distance to the telomere
        coefs_dist_to_telomere = poly.polyfit(df_correct.distance_to_telomere, df_correct.coverage, 2)

        # get the residual variation in coverage
        df_correct["coverage_from_dist_to_telomere"] = poly.polyval(df_correct.distance_to_telomere, coefs_dist_to_telomere)
        df_correct["residualCoverage_from_dist_to_telomere"] = df_correct.coverage - df_correct.coverage_from_dist_to_telomere

        # get a quadratic fit that predicts coverage from GC content
        coefs_GCcontent = poly.polyfit(df_correct.GCcontent, df_correct.residualCoverage_from_dist_to_telomere, 2)
        df_correct["residualCoverage_from_dist_to_telomere_from_GC_content"] = poly.polyval(df_correct.GCcontent, coefs_GCcontent)

        df_correct["coverage_from_dist_to_telomere_and_GC_content"] = df_correct["coverage_from_dist_to_telomere"] + df_correct["residualCoverage_from_dist_to_telomere_from_GC_content"]

        # save the coefficients
        for chrom in chroms: chrom_to_coefType_to_coefs[chrom] = {"dist_telomere":coefs_dist_to_telomere, "GCcontent":coefs_GCcontent}

        # plot the coverage for each of the chromosomes
        fig = plt.figure(figsize=(7, len(chroms)*5))
        for I, chrom in enumerate(chroms):

            # initialize a subplot, where each row is one chromosome
            ax = plt.subplot(len(chroms), 1, I+1)

            # get df of this chrom
            df_c = df_correct[df_correct.chromosome==chrom]

            # make a line plot for the real coverage
            plt.scatter(df_c.start, df_c.coverage, marker="o", color="gray", label="data")

            # make a line for the prediction from the distance to the telomere
            plt.plot(df_c.start, df_c.coverage_from_dist_to_telomere, linestyle="-", color="blue", label="pred_dist_telomere")

            # make a line for the prediction for both
            plt.plot(df_c.start, df_c.coverage_from_dist_to_telomere_and_GC_content, linestyle="-", color="red", label="pred_dist_and_gc_content")

            # add a line with the distance to the telomere
            #plt.plot(df_c.start, df_c.distance_to_telomere, linestyle="-", color="green", label="dist_telomere")

            ax.legend()
            ax.set_ylabel("coverage")
            ax.set_xlabel("position (bp)")
            ax.set_title(chrom)

        # save
        outfile = "%s/coverage_modelling_%s.pdf"%(outdir, type_genome)
        #print("saving %s"%outfile)
        fig.savefig(outfile, bbox_inches="tight")

        # get the rsquare of the model
        r2 = r2_score(df_correct.coverage, df_correct.coverage_from_dist_to_telomere_and_GC_content)
        print("The rsquare for %s is %.3f"%(type_genome, r2))

    ###############################################################

    # define the function that takes a tuple of (distToTelomere, chromosome and GCcontent) and returns the predicted relative coverage
    final_function = (lambda dist_telomere, chrom, GCcontent:  # this is suposed to be the tuple

                        (poly.polyval([dist_telomere], chrom_to_coefType_to_coefs[chrom]["dist_telomere"]) + # from the dist to tel
                        poly.polyval([GCcontent], chrom_to_coefType_to_coefs[chrom]["GCcontent"]))[0] # residual predicted from GC

                     )

    # check that it works
    df_correct["cov_predicted_from_final_lambda"] = df_correct.apply(lambda r: final_function(r["distance_to_telomere"], r["chromosome"], r["GCcontent"]), axis=1)

    if any(((df_correct["coverage_from_dist_to_telomere_and_GC_content"]-df_correct["cov_predicted_from_final_lambda"]).apply(abs))>0.01): raise ValueError("error in lambda function generation for coverage")

      
    return final_function


def get_clove_output(output_vcf_clove):

    """Gets the raw output of clove and returns a df with it"""

    # load df
    df = pd.read_csv(output_vcf_clove, skiprows=list(range(len([line for line in open(output_vcf_clove, "r", encoding='utf-8', errors='ignore') if line.startswith("##")]))), sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)

    # get FORMAT into several cells
    INFOfields_data = pd.DataFrame(dict(df.INFO.apply(lambda x: {content.split("=")[0] : content.split("=")[1]  for content in  make_flat_listOflists([y.split(";") for y in x.split("; ")])}))).transpose()
    df = df.merge(INFOfields_data, left_index=True, right_index=True, validate="one_to_one")

    # debug that some of the CHR2 are NaN or without N
    if any([pd.isna(x) for x in set(df.CHR2)]): raise ValueError("The CLOVE output parsing was incorrect for %s"%output_vcf_clove)
    if any([pd.isna(x) for x in set(df.END)]): raise ValueError("The CLOVE output parsing was incorrect for %s"%output_vcf_clove)
    if any([pd.isna(x) for x in set(df.SVTYPE)]): raise ValueError("The CLOVE output parsing was incorrect for %s"%output_vcf_clove)


    # change the END by -1 if it is NaN
    def getNaN_to_minus1(x):

        if pd.isna(x): return -1
        else: return int(x)

    if "START" in df.keys(): df["START"] = df.START.apply(getNaN_to_minus1)
    else: df["START"] = [-1]*len(df)

    df["END"] = df.END.apply(getNaN_to_minus1)

    return df

def get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, windows_file, replace=False, run_in_parallel=True, delete_bams=False):

    """This function takes a windows file and a bam, and it runs generate_coverage_per_window_file_parallel but only for regions that are not previously calculated"""
    #print("calculating coverage per windows on %s"%sorted_bam)

    # define bams
    outdir_bam = "/".join(sorted_bam.split("/")[0:-1])
    sorted_bam_name = sorted_bam.split("/")[-1]

    # define the query_windows
    query_windows_df = pd.read_csv(windows_file, sep="\t").set_index(["chromosome", "start", "end"], drop=False)
    if len(query_windows_df)==0: return pd.DataFrame()

    # define the file were the coverage will be calculated
    calculated_coverage_file_prefix = "%s.coverage_per_window.tab"%sorted_bam_name

    # define all previous files, all those that do not end with tmp 
    previously_calculated_windows_coverage_files = ["%s/%s"%(outdir_bam, x) for x in os.listdir(outdir_bam) if x.startswith(calculated_coverage_file_prefix) and "temporary_file" not in x]

    ######### define previously measured regions ##########

    if len(previously_calculated_windows_coverage_files)>0:

        # get into a df all the previously calculated coverages
        df_previosuly_calculated_coverages = pd.concat([pd.read_csv(file, sep="\t").set_index(["chromosome", "start", "end"], drop=False) for file in previously_calculated_windows_coverage_files], sort=True)

        # define the windows_to_measure_df as those that are not in df_previosuly_calculated_coverages
        windows_still_to_measure = set(query_windows_df.index).difference(set(df_previosuly_calculated_coverages.index))

        if len(windows_still_to_measure)==0: windows_to_measure_df = pd.DataFrame()
        else: windows_to_measure_df = query_windows_df.loc[list(windows_still_to_measure)]

    else: 
        windows_to_measure_df = query_windows_df
        df_previosuly_calculated_coverages = pd.DataFrame()

    ###################################

    ######## measure regions, keeping in tmp ######

    if len(windows_to_measure_df)>0:

        # write a bed with the regions to measure
        bed_windows_to_measure = "%s/%s.%s.temporary_file"%(outdir_bam, calculated_coverage_file_prefix, id_generator(20))
        windows_to_measure_df[["chromosome", "start", "end"]].to_csv(bed_windows_to_measure, sep="\t", header=False, index=False)

        # get the coverage df
        destination_dir = "%s.calculating_windowcoverage"%sorted_bam
        coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file=bed_windows_to_measure, replace=replace, window_l=1000, run_in_parallel=run_in_parallel, delete_bams=delete_bams), sep="\t").rename(columns={"#chrom":"chromosome"}).set_index(["chromosome", "start", "end"], drop=False)

        # save file
        coverage_df_file = "%s/%s.%s"%(outdir_bam, calculated_coverage_file_prefix, id_generator(20)); coverage_df_file_tmp = "%s.temporary_file"%coverage_df_file
        coverage_df.to_csv(coverage_df_file_tmp, sep="\t", header=True, index=False)

        # remove the bed
        remove_file(bed_windows_to_measure)
        remove_file("%s.coverage_provided_windows.tab"%bed_windows_to_measure)

        # rename to keep
        os.rename(coverage_df_file_tmp, coverage_df_file)

    else: coverage_df = pd.DataFrame()

    #############################################

    # append the dataframes
    df_coverage_all = df_previosuly_calculated_coverages.append(coverage_df, sort=True).loc[list(query_windows_df.index)]
    df_coverage_all.index = list(range(len(df_coverage_all)))

    # drop duplicates
    df_coverage_all = df_coverage_all.drop_duplicates(subset=["chromosome", "start", "end", "mediancov_1", "percentcovered_1"], keep='first') 

    return df_coverage_all

def get_int(x):

    # def get nans to -1

    try: return int(x)
    except: return -1

def get_coverage_list_relative_to_predictedFromTelomereAndGCcontent(df_cov, genome, distToTel_chrom_GC_to_coverage_fn, genome_graph, df_positions_graph, outdir, real_coverage_field="mediancov_1", replace=False):

    """This function takes a df with coverage for some windows of the genome and coverage. The returned list is in the same order as the df_cov"""


    print("getting coverage relative to the one predicted from seq features")

    # make the outdir
    make_folder(outdir)

    # rename the training df and copy
    df = cp.deepcopy(df_cov.rename(columns={"#chrom":"chromosome", real_coverage_field:"coverage"}))
    df.index = list(range(len(df)))
    initial_index = list(df.index)

    # add the distance to the telomere
    df_with_distance_to_telomere_file = "%s/df_with_distance_to_telomere_file.py"%outdir
    if file_is_empty(df_with_distance_to_telomere_file) or replace is True:

        df["middle_position"] = (df.start + (df.end - df.start)/2).apply(int)
        df["distance_to_telomere"] = get_distance_to_telomere_series(df[["chromosome", "middle_position"]], genome_graph, df_positions_graph)

        save_object(df, df_with_distance_to_telomere_file)

    else: df = load_object(df_with_distance_to_telomere_file)

    # add the gc content
    gcontent_outfile = "%s/GCcontent.py"%outdir
    df = get_df_with_GCcontent(df, genome, gcontent_outfile, replace=replace)

    # predict genome from the sequence features 
    df["cov_predicted_from_features"] = df.apply(lambda r: distToTel_chrom_GC_to_coverage_fn(r["distance_to_telomere"], r["chromosome"], r["GCcontent"]), axis=1)

    # get the relative to the 
    df["cov_rel_to_predFromFeats"] = df.coverage/df.cov_predicted_from_features
    if any(pd.isna(df["cov_rel_to_predFromFeats"])): raise ValueError("There was a problem with the prediction from features")

    # get the final list
    return list(df.loc[initial_index, "cov_rel_to_predFromFeats"])

def get_clove_output_with_coverage_forTANDEL(outfile_clove, reference_genome, sorted_bam, distToTel_chrom_GC_to_coverage_fn, genome_graph, df_positions_graph, replace=False, run_in_parallel=False, delete_bams=False):

    """Takes the output of clove and adds the coverage of the TAN and DEL, or -1"""

    # first load clove into a df
    df_clove = get_clove_output(outfile_clove)

    if len(df_clove)>0:

        # now write a bed with the TANDEL regions
        bed_TANDEL_regions = "%s.TANDEL.bed"%outfile_clove
        df_TANDEL = df_clove[df_clove.SVTYPE.isin({"TAN", "DEL"})][["#CHROM", "POS", "END"]].rename(columns={"#CHROM":"chromosome", "POS":"start", "END":"end"})
        
        if len(df_TANDEL)>0:

            df_TANDEL.to_csv(bed_TANDEL_regions, sep="\t", header=True, index=False)

            # get a file that has the coverage of these windows
            coverage_df = get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, bed_TANDEL_regions, replace=replace, run_in_parallel=run_in_parallel, delete_bams=delete_bams)

            # get the coverage relative to prediction from features
            outdir_rel_cov_calculation = "%s_calculatig_rel_coverage"%outfile_clove
            coverage_df["coverage_rel_to_predFromFeats"] = get_coverage_list_relative_to_predictedFromTelomereAndGCcontent(coverage_df, reference_genome, distToTel_chrom_GC_to_coverage_fn, genome_graph, df_positions_graph, outdir_rel_cov_calculation, real_coverage_field="mediancov_1", replace=replace)

        else: coverage_df = pd.DataFrame(columns=["chromosome", "end", "length", "mediancov_1", "nocoveragebp_1", "percentcovered_1", "start"])

        # merge
        merged_df = df_clove.merge(coverage_df, how="left", left_on=["#CHROM", "POS", "END"], right_on=["chromosome", "start", "end"], validate="many_to_one")

        # change types of fields
        merged_df["POS"] = merged_df.POS.apply(get_int)
        merged_df["END"] = merged_df.END.apply(get_int)
        merged_df["START"] = merged_df.START.apply(get_int)

        return merged_df 

    else: return pd.DataFrame()


def get_covfilter_cloveDF_row_according_to_SVTYPE(r, maxDELcoverage=0.1, minDUPcoverage=1.9, coverage_field="coverage_rel_to_predFromFeats"):

    # define a function that takes a row of the dataframe and does the filterinf can be any that is in coverage

    if r["SVTYPE"]=="DEL":

        if r[coverage_field]<=maxDELcoverage: return "PASS"
        else: return "FAIL" 

    elif r["SVTYPE"]=="TAN":

        if r[coverage_field]>=minDUPcoverage: return "PASS"
        else: return "FAIL" 

    else: return "PASS"

def run_gridssClove_given_filters(sorted_bam, reference_genome, working_dir, median_coverage, replace=True, threads=4, gridss_blacklisted_regions="", gridss_VCFoutput="", gridss_maxcoverage=50000, median_insert_size=500, median_insert_size_sd=0, gridss_filters_dict=default_filtersDict_gridss, tol_bp=50, threshold_p_unbalTRA=0.7, run_in_parallel=True, max_rel_coverage_to_consider_del=0.1, min_rel_coverage_to_consider_dup=1.5, replace_FromGridssRun=False, include_breakpoints_in_genomeGraph=True, mitochondrial_chromosome="mito_C_glabrata_CBS138", type_coverage_to_filterTANDEL="coverage_rel_to_predFromFeats"):

    """This function runs gridss and clove with provided filtering and parameters. This can be run at the end of a parameter optimisation process. It returns a dict mapping each SV to a table, and a df with the gridss.

    coverage_field is the field where clove is filtered to detect CNV. It can be relative_coverage or relative_coverage_dist_to_telomere.

    include_breakpoints_in_genomeGraph indicates whether breakpoints should be included in the genome graph

    type_coverage_to_filterTANDEL is a field that determines what type of coverage is used to filter deletions and tandem duplications in the CLOVE output:

    - mediancov_1 indicates that it should be the abslute and raw coverage
    - coverage_rel_to_predFromFeats means that it should be  based on the coverage relative to the predicted from the features


    """

    print("running gridss and clove with given parameter")
    make_folder(working_dir)

    # edit the replace, regarding if filtering from the run of GRIDSS
    if replace is True and replace_FromGridssRun is False: replace_FromGridssRun = True

    # first obtain the gridss output if it is not provided
    if file_is_empty(gridss_VCFoutput) or replace is True: gridss_VCFoutput = run_gridss_and_annotateSimpleType(sorted_bam, reference_genome, working_dir, replace=replace, threads=threads, blacklisted_regions=gridss_blacklisted_regions, maxcoverage=gridss_maxcoverage)

    #################################################
    ##### GET A LIST OF FILTERED BREAKPOINTS ########
    #################################################
    
    # get the output of gridss into a df
    print("getting gridss")
    df_gridss = add_info_to_gridssDF(load_single_sample_VCF(gridss_VCFoutput), median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd) # this is a dataframe with some info

    # filter according to gridss_filters_dict
    print("filtering gridss")
    df_gridss_filt = get_gridssDF_filtered_from_filtersDict(df_gridss, gridss_filters_dict)
    gridss_VCFoutput_filt = "%s.filtered_default.vcf"%(gridss_VCFoutput)

    # write the filtered gridss vcf
    df_text = df_gridss_filt[list(df_gridss_filt.keys())[0:10]].to_csv(sep="\t", header=True, index=False)
    header_lines = "".join([l for l in open(gridss_VCFoutput, "r").readlines() if l.startswith("##")])
    open(gridss_VCFoutput_filt, "w").write("%s%s"%(header_lines, df_text))

    # get bedpe into a file, only with the minimal fields
    bedpe_with_adds = get_bedpe_from_svVCF(gridss_VCFoutput_filt, working_dir, replace=replace_FromGridssRun, only_simple_conversion=True) # REPLACE debug

    # check that it is not empty
    if open(bedpe_with_adds, "r").readlines()[0].startswith("no_vcf_records"): return {}, df_gridss

    # write the raw fields
    raw_bedpe_file = "%s.raw.bedpe"%bedpe_with_adds
    bedpe_fields = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]
    df_bedpe = pd.read_csv(bedpe_with_adds, sep="\t")[bedpe_fields]
    df_bedpe.to_csv(raw_bedpe_file, sep="\t", header=False, index=False)
    print("there are %i breakpoints"%len(df_bedpe))


    ##### GRAPH GENOME OPERATIONS #####

    # get a graph of the genome
    genomeGraph_outfileprefix = "%s.genomeGraph_incluingBPs%s"%(raw_bedpe_file, include_breakpoints_in_genomeGraph)
    if include_breakpoints_in_genomeGraph is True:
        df_bedpe_arg = df_bedpe
        df_gridss_filt_arg = df_gridss_filt
    else:
        df_bedpe_arg = None
        df_gridss_filt_arg = None

    genome_graph, df_positions_graph = get_genomeGraph_object(reference_genome, df_bedpe_arg, df_gridss_filt_arg, genomeGraph_outfileprefix, replace=replace)

    # get a function that takes the GC content, chromosome and distance to the telomere and returns coverage. This is actually a lambda function
    outdir_coverage_calculation = "%s/coverage_per_regions2kb_incluingBPs%s"%(working_dir, include_breakpoints_in_genomeGraph); make_folder(outdir_coverage_calculation)
    df_coverage_train = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_coverage_calculation, sorted_bam, windows_file="none", replace=replace, window_l=2000), sep="\t")

    distToTel_chrom_GC_to_coverage_fn = get_distanceToTelomere_chromosome_GCcontent_to_coverage_fn(df_coverage_train, reference_genome, genome_graph, df_positions_graph, outdir_coverage_calculation, mitochondrial_chromosome=mitochondrial_chromosome, replace=replace)

    ###################################

    #################################################
    #################################################
    #################################################

    # run clove without checking filtering
    outfile_clove = "%s.clove.vcf"%(raw_bedpe_file)
    run_clove_filtered_bedpe(raw_bedpe_file, outfile_clove, sorted_bam, replace=replace_FromGridssRun, median_coverage=median_coverage, median_coverage_dev=1, check_coverage=False) #  REPLACE debug

    # add the filter of coverage to the clove output
    df_clove = get_clove_output_with_coverage_forTANDEL(outfile_clove, reference_genome, sorted_bam, distToTel_chrom_GC_to_coverage_fn,  genome_graph, df_positions_graph, replace=replace_FromGridssRun, run_in_parallel=run_in_parallel, delete_bams=run_in_parallel)

    # define the coverage filtering based on the type_coverage_to_filterTANDEL
    if type_coverage_to_filterTANDEL=="mediancov_1":
        maxDELcoverage = int(max_rel_coverage_to_consider_del*median_coverage)
        minDUPcoverage = int(min_rel_coverage_to_consider_dup*median_coverage) 

    elif type_coverage_to_filterTANDEL=="coverage_rel_to_predFromFeats":
        maxDELcoverage = max_rel_coverage_to_consider_del
        minDUPcoverage = min_rel_coverage_to_consider_dup

    df_clove["coverage_FILTER"] = df_clove.apply(lambda r: get_covfilter_cloveDF_row_according_to_SVTYPE(r, maxDELcoverage=maxDELcoverage, minDUPcoverage=minDUPcoverage, coverage_field="coverage_rel_to_predFromFeats"), axis=1)

    print(df_clove)




    khadkjghadkjdakghj

    # annotated clove 
    fileprefix = "%s.structural_variants"%outfile_clove
    remaining_df_clove, svtype_to_SVtable = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, fileprefix, reference_genome, sorted_bam, tol_bp=tol_bp, replace=replace_FromGridssRun, median_coverage=median_coverage, svtypes_to_consider={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}, threshold_p_unbalTRA=threshold_p_unbalTRA, run_in_parallel=run_in_parallel)

    # merge the coverage files in one
    merge_coverage_per_window_files_in_one(sorted_bam)

    return svtype_to_SVtable, df_gridss


###################################################################################################
###################################################################################################
###################################################################################################

def generate_tables_of_SV_between_genomes_gridssClove(query_genome, reference_genome, replace=False, threads=4, coverage=10, insert_size=250, read_lengths=[kb*1000 for kb in [0.5, 0.7, 0.9, 1]], error_rate=0.0, gridss_min_af=0.25, mitochondrial_chromosome="mito_C_glabrata_CBS138"):

    """Takes a bam file with aligned reads or genomes and generates calls, returning a dict that maps variation type to variants
    - aligner can be minimap2 or ngmlr"""

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
    gridss_filters_dict["min_af"] = gridss_min_af
    print("Filtering out when any AF is below %.3f"%(gridss_filters_dict["min_af"]))

    # define the median coverage per region
    outdir_coverage_calculation = "%s/coverage_per_regions"%working_dir; make_folder(outdir_coverage_calculation)
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_coverage_calculation, sorted_bam, windows_file="none", replace=replace, window_l=10000), sep="\t")
    median_coverage = np.median(coverage_df.mediancov_1)
    print("The median coverage is %i"%median_coverage)

    # run the gridss and clove pipeline with high-confidence parameters
    gridss_outdir = "%s/%s_gridss_outdir"%(working_dir, ID)
    SV_dict, df_gridss =  run_gridssClove_given_filters(sorted_bam, reference_genome, gridss_outdir, median_coverage, replace=replace, threads=threads, median_insert_size=insert_size, gridss_filters_dict=gridss_filters_dict, replace_FromGridssRun=False) # DEBUG. The replace_FromGridssRun in True would be debugging is to replace from the GRIDSS run step

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
    make_folder(outdir)

    ##### test how well the finding of SVs in an assembly works #####
    outdir_test_FindSVinAssembly = "%s/test_FindSVinAssembly"%outdir; make_folder(outdir_test_FindSVinAssembly)

    print("WORKING ON THE VALIDATION THAT WE CAN FIND READS IN AN ASSEMBLY")

    # go through each simulation
    for simID in range(n_simulated_genomes):

        # define outdir 
        outdir_sim = "%s/simulation_%i"%(outdir_test_FindSVinAssembly, simID); make_folder(outdir_sim)

        # generate genome with simulated SVs
        sim_svtype_to_svfile, rearranged_genome = rearrange_genomes_simulateSV(reference_genome, outdir_sim, replace=replace, nvars=nvars, mitochondrial_chromosome=mitochondrial_chromosome)

        # get the variants from simulating reads. Always ploidy 1 to get homozygous SVs
        predicted_svtype_to_svfile = generate_tables_of_SV_between_genomes_gridssClove(rearranged_genome, reference_genome, replace=replace, threads=threads, mitochondrial_chromosome=mitochondrial_chromosome)

        finsihedrunningGridssClove








    #simType_to_rearrangedGenome = rearrange_genomes_simulateSV(reference_genome, simulation_outdir, repeats_1kbwindow_bed=repeats_regions_bed, replace=replace, simulation_types=simulation_types, mitochondrial_chromosome=mitochondrial_chromosome) # !!!! i have to change the representation of translocations with different 3' and 5' end



    #################################################################



    


    # at the end clean the reference genome
    clean_reference_genome_windows_files(reference_genome)






