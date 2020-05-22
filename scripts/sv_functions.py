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
from ete3 import Tree

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
JolyTree_sh = "%s/bin/JolyTree.sh"%EnvDir

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


# define default parameters for gridss filtering. This has changed from v0
default_filtersDict_gridss = {"min_Nfragments":5, "min_af":0.25, "wrong_FILTERtags":("NO_ASSEMBLY",), "filter_polyGC":True, "filter_noSplitReads":False, "filter_noReadPairs":False, "maximum_strand_bias":0.99, "maximum_microhomology":50, "maximum_lenght_inexactHomology":50, "range_filt_DEL_breakpoints":(0, 1), "min_length_inversions":40, "dif_between_insert_and_del":5, "max_to_be_considered_small_event":1000, "wrong_INFOtags":('IMPRECISE',), "min_size":50, "min_af_EitherSmallOrLargeEvent":0.25} # the minimum af is 0.25 to include both heterozygous and homozygous vars as default

# define lists of filters ordered from less conservative to most conservative
g_all_FILTER_tags = ("ASSEMBLY_ONLY", "NO_ASSEMBLY", "ASSEMBLY_TOO_FEW_READ", "ASSEMBLY_TOO_SHORT", "INSUFFICIENT_SUPPORT", "LOW_QUAL", "REF", "SINGLE_ASSEMBLY")
g_meaningful_FILTER_tags = ("NO_ASSEMBLY", "INSUFFICIENT_SUPPORT", "LOW_QUAL")

g_min_Nfragments_l = [5, 8, 10, 15, 20, 30]
g_min_af_l = [0.0, 0.05, 0.1, 0.2, 0.5, 0.9, 0.99]
g_min_af_EitherSmallOrLargeEvent_l = [0.0, 0.05, 0.1, 0.2, 0.5, 0.9, 0.99]
g_wrong_FILTERtags_l = [("",), ("NO_ASSEMBLY",), ("NO_ASSEMBLY", "INSUFFICIENT_SUPPORT"), ("NO_ASSEMBLY", "LOW_QUAL"), ("LOW_QUAL", "INSUFFICIENT_SUPPORT"), g_meaningful_FILTER_tags, g_all_FILTER_tags] 
g_filter_polyGC_l = [False, True]
g_filter_noSplitReads_l = [False, True]
g_filter_noReadPairs_l = [False, True]
g_maximum_strand_bias_l = [0.99, 0.95, 0.9]    
g_maximum_microhomology_l = [200, 100, 50, 10]
g_maximum_lenght_inexactHomology_l = [200, 100, 50, 10]
g_range_filt_DEL_breakpoints_l = [(0,1), (200, 700), (100, 800), (50, 900)]
g_min_length_inversions_l = [40, 50, 60, 70]
g_dif_between_insert_and_del_l = [0, 5, 10, 20]
g_max_to_be_considered_small_event_l = [100, 200, 500, 1000, 1500]
g_wrong_INFOtags_l = [('NONE',), ("IMPRECISE",)]
g_min_size_l = [50]

# map each filter to the ordered list
g_filterName_to_filtersList = {"min_Nfragments":g_min_Nfragments_l, "min_af":g_min_af_l, "wrong_FILTERtags":g_wrong_FILTERtags_l, "filter_polyGC":g_filter_polyGC_l, "filter_noSplitReads":g_filter_noSplitReads_l, "filter_noReadPairs":g_filter_noReadPairs_l, "maximum_strand_bias":g_maximum_strand_bias_l, "maximum_microhomology":g_maximum_microhomology_l, "maximum_lenght_inexactHomology":g_maximum_lenght_inexactHomology_l, "range_filt_DEL_breakpoints":g_range_filt_DEL_breakpoints_l, "min_length_inversions":g_min_length_inversions_l, "dif_between_insert_and_del":g_dif_between_insert_and_del_l, "max_to_be_considered_small_event":g_max_to_be_considered_small_event_l, "wrong_INFOtags":g_wrong_INFOtags_l, "min_size":g_min_size_l, "min_af_EitherSmallOrLargeEvent":g_min_af_EitherSmallOrLargeEvent_l}

# map each value of each filter list to a value depending on the position
g_filterName_to_filterValue_to_Number = {filterName : dict(zip(filtersList, range(len(filtersList)))) for filterName, filtersList in g_filterName_to_filtersList.items()}

# define a dict that maps each svtype to the fields that are important to define the overlaps
svtype_to_fieldsDict = {"inversions": {"equal_fields": ["Chr"], 
                                        "approximate_fields": ["Start", "End"],
                                        "chromField_to_posFields": {"Chr":{"start": "Start", "end": "End"}},
                                        "all_fields": ["Chr", "Start", "End", "ID"],
                                        "position_fields": ["Start", "End"],
                                        "positionField_to_chromosome": {"Start":"Chr", "End":"Chr"}
                                        }, 

                        "tandemDuplications": {"equal_fields": ["Chr"], 
                                               "approximate_fields": ["Start", "End"],
                                                "chromField_to_posFields": {"Chr":{"start": "Start", "end": "End"}},
                                                "all_fields": ["Chr", "Start", "End", "ID"],
                                                "position_fields": ["Start", "End"],
                                                "positionField_to_chromosome": {"Start":"Chr", "End":"Chr"}
                                                }, 

                        "deletions": {"equal_fields": ["Chr"], 
                                      "approximate_fields": ["Start", "End"],
                                      "chromField_to_posFields": {"Chr":{"start": "Start", "end": "End"}},
                                      "all_fields": ["Chr", "Start", "End", "ID"],
                                      "position_fields": ["Start", "End"],
                                      "positionField_to_chromosome": {"Start":"Chr", "End":"Chr"}
                                     }, 

                        "insertions": {"equal_fields": ["ChrA", "ChrB", "Copied"], 
                                        "approximate_fields": ["StartA", "EndA", "StartB", "EndB"],
                                        "chromField_to_posFields": {"ChrA":{"start": "StartA", "end": "EndA"}},
                                        "all_fields": ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Copied", "ID"],
                                        "position_fields": ["StartA", "EndA", "StartB", "EndB"],
                                        "positionField_to_chromosome": {"StartA":"ChrA", "EndA":"ChrA", "StartB":"ChrB", "EndB":"ChrB"}
                                        }, 

                        "translocations": {"equal_fields": ["ChrA", "ChrB", "Balanced"], 
                                           "approximate_fields": ["StartA", "EndA", "StartB", "EndB"],
                                           "chromField_to_posFields": {},
                                           "all_fields": ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Balanced", "ID"],
                                           "position_fields": ["StartA", "EndA", "StartB", "EndB"],
                                           "positionField_to_chromosome": {"StartA":"ChrA", "EndA":"ChrA", "StartB":"ChrB", "EndB":"ChrB"}
                                          }, 

                        "remaining": {"equal_fields": ["#CHROM", "CHR2", "SVTYPE"], 
                                      "approximate_fields": ["POS", "START", "END"],
                                      "chromField_to_posFields": {"CHR2":{"start":"START", "end":"END"}},
                                      "all_fields": ["#CHROM", "POS", "CHR2", "START", "END", "SVTYPE", "ID"],
                                      "position_fields": ["POS", "START", "END"],
                                      "positionField_to_chromosome": {"POS":"#CHROM", "START":"CHR2", "END":"CHR2"}

                                      }}

####################################
####################################
####################################

def id_generator(size=10, chars=string.ascii_uppercase + string.digits, already_existing_ids=set()):

    """ already_existing_ids is a set that indicates whihc IDs can't be picked """

    ID = ''.join(random.choice(chars) for _ in range(size))
    while ID in already_existing_ids:
        ID = ''.join(random.choice(chars) for _ in range(size))

    return ID


def find_nearest(a, a0):

    """Element in nd array `a` closest to the scalar value `a0`"""
    
    # Debug elements that are inf
    if a0 not in [np.inf, -np.inf]:
        a = np.array(a)
        idx = np.abs(a - a0).argmin()
        closest_in_a = a.flat[idx]
        
    elif a0==np.inf:
        closest_in_a = max(a)
        
    elif a0==-np.inf:
        closest_in_a = min(a)        

    return closest_in_a

def chunks(l, n):
    
    """Yield successive n-sized chunks from a list l"""
    
    for i in range(0, len(l), n):
        yield l[i:i + n]

def get_dict_as_tuple(dictionary):

    """Takes a dict and converts it to a sorted tuple"""

    return tuple([(k, dictionary[k]) for k in sorted(dictionary.keys())])

def union_empty_sets(set_iterable):

    """It returns the union of an iterable of empty sets"""

    return set.union(*list(set_iterable) + [set()])

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

def get_chr_to_len(genome):

    # define chromosome_to_length for a genome
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(genome, "fasta")}

    return chr_to_len

def get_uniqueVals_df(df): return set.union(*[set(df[col]) for col in df.columns])

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

def get_affected_region_bed_for_SVdf(svDF, svtype, interesting_chromosomes, add_interval_bp=1000, first_position_idx=0, translocations_type="breakpoint_pos", chr_to_len={}):

    """This function takes a df with SVs and returns the regions that are affected in bed format, which depends on the svtype. It adds an interval arround this values. It also returns the number 

    only consider svs were ANY of the chroms are in interesting_chromosomes.

    The bed for the translocations includes the left and right regions of all positions.

    translocations_type can be 'start_and_end_pos' or 'whole_arms' or 'whole_chromosomes', which will imply that only the regions arround breakpoints or the whole arms will be returned. chr_to_len is only needed when you have whole_arms_2orientations"""


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

            if len(chr_to_len)==0: raise ValueError("chr_to_len has to be full for translocations")

            # count the SVs
            nSVs_in_interesting_chromosomes = sum(svDF.apply(lambda r: r["ChrA"] in interesting_chromosomes or r["ChrB"] in interesting_chromosomes, axis=1))

            # get only the breakpoint pos
            if translocations_type=="breakpoint_pos":

                # go through each position
                affected_region_bed_df = pd.DataFrame()
                for chrom_f, pos_f in [("ChrA", "StartA"), ("ChrA", "EndA"), ("ChrB", "StartB"), ("ChrB", "EndB")]:

                    # go through each region
                    for chrom, pos in svDF[[chrom_f, pos_f]].values:

                        # do not consider the first positions and last ones
                        if pos>first_position_idx and pos<chr_to_len[chrom]:

                            # else add as a target region
                            affected_region_bed_df = affected_region_bed_df.append(pd.DataFrame({0: {"chromosome":chrom, "start":pos, "end":pos}}).transpose())

            # get the translocations bed as the breakpont, start and end regions of chromosomes
            elif translocations_type=="start_and_end_pos":

                # go through each position
                affected_region_bed_df = pd.DataFrame()
                for chrom_f, pos_f in [("ChrA", "StartA"), ("ChrA", "EndA"), ("ChrB", "StartB"), ("ChrB", "EndB")]:

                    # go through each region
                    for chrom, pos in svDF[[chrom_f, pos_f]].values:

                        # else add as a target region
                        affected_region_bed_df = affected_region_bed_df.append(pd.DataFrame({0: {"chromosome":chrom, "start":pos, "end":pos}}).transpose())

            # get the whole arms of the chromosome as bed
            elif translocations_type=="whole_arms":

                # go through each position
                affected_region_bed_df = pd.DataFrame()
                for chrom_f, start_f, end_f in [("ChrA", "StartA", "EndA"), ("ChrB", "StartB", "EndB")]:

                    # go through each region
                    for chrom, start, end in svDF[[chrom_f, start_f, end_f]].values:

                        # else add as a target region
                        affected_region_bed_df = affected_region_bed_df.append(pd.DataFrame({0: {"chromosome":chrom, "start":start, "end":end}}).transpose())

            # get all the chromosome as affected region
            elif translocations_type=="whole_chromosomes": 

                affected_region_bed_df = pd.DataFrame({I: {"chromosome":chromosome, "start":1, "end":chr_to_len[chromosome]} for I, chromosome in enumerate(get_uniqueVals_df(svDF[["ChrA", "ChrB"]])) }).transpose().drop_duplicates()

        else: raise ValueError("%s is not valid"%(svtype))


        # add the confidence arround each region
        affected_region_bed_df["start"] = (affected_region_bed_df.start - add_interval_bp).apply(lambda x: max([first_position_idx,x]))
        affected_region_bed_df["end"] = affected_region_bed_df.end + add_interval_bp

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

    # keep an uncorrected version
    run_cmd("cp %s %s.unmodified"%(insertions_file, insertions_file))

    # write as corrected
    df.to_csv(insertions_file, sep="\t", header=True, index=False)


def get_breakpoint_positions_df_in_svDF(svDF):

    """Takes and svDF and returns a df with the chromosomes and positions of the breakpoimts"""

    # small vars
    if "Chr" in svDF.keys(): 
        df_start = svDF[["Chr", "Start"]].rename(columns={"Chr":"Chr", "Start":"Pos"})
        df_end = svDF[["Chr", "End"]].rename(columns={"Chr":"Chr", "End":"Pos"})

        df_pos = df_start.append(df_end)

    # bedpe vars
    else:
        df_startA = svDF[["ChrA", "StartA"]].rename(columns={"ChrA":"Chr", "StartA":"Pos"})
        df_startB = svDF[["ChrB", "StartB"]].rename(columns={"ChrB":"Chr", "StartB":"Pos"})
        df_endA = svDF[["ChrA", "EndA"]].rename(columns={"ChrA":"Chr", "EndA":"Pos"})
        df_endB = svDF[["ChrB", "EndB"]].rename(columns={"ChrB":"Chr", "EndB":"Pos"})

        df_pos = pd.concat([d for d in [df_startA, df_startB, df_endA, df_endB]])

    return df_pos

def transform_cut_and_paste_to_copy_and_paste_insertions(reference_genome, rearranged_genome, insertions_file, svtype_to_svDF):

    """ This function takes a rearranged genome and reinserts the copy-and-paste insertions where they should be """

    print("reinserting-copy-and-paste insertions into %s"%insertions_file)

    # load df and keep the copy-and-paste insertions
    df = pd.read_csv(insertions_file, sep="\t")
    df = df[df.Copied]

    if len(df)>0:

        # define an unmodified genome
        rearranged_genome_unmodified = "%s.unmodified.fasta"%rearranged_genome
        rearranged_genome_unmodified_tmp = "%s.tmp"%rearranged_genome_unmodified

        if file_is_empty(rearranged_genome_unmodified):

            # if the unmodified tmps is writen, replace the rearranged_genome with it
            if not file_is_empty(rearranged_genome_unmodified_tmp): os.rename(rearranged_genome_unmodified_tmp, rearranged_genome)

            # get the rearranged genome seq
            chr_to_rearrangedSeq = {seq.id: str(seq.seq) for seq in SeqIO.parse(rearranged_genome, "fasta")}
            all_rearranged_chromosomes_together = "".join(chr_to_rearrangedSeq.values())

            # get the seq
            chr_to_refSeq = {seq.id: str(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

            # define the length of each chrom
            chr_to_lenSeq = {chrom : len(seq) for chrom, seq in chr_to_refSeq.items()}

            # define all the positions with breakpoints
            df_positions = pd.concat([get_breakpoint_positions_df_in_svDF(svDF) for svtype, svDF in svtype_to_svDF.items()])
            chr_to_bpPositions = dict(df_positions.groupby("Chr").apply(lambda df_c: set(df_c["Pos"])))

            # add the ends of the chromosome, and convert to np array
            for chrom, lenSeq in chr_to_lenSeq.items(): 

                chr_to_bpPositions[chrom].update({1, lenSeq})
                chr_to_bpPositions[chrom] = np.array(sorted(chr_to_bpPositions[chrom]))

            # add the closest breakpoint position of ChrA in the reference
            df["closest_5'breakpoint_position"] = df.apply(lambda r: find_nearest(chr_to_bpPositions[r["ChrA"]][chr_to_bpPositions[r["ChrA"]]<(r["StartA"])], r["StartA"]), axis=1)

            df["closest_3'breakpoint_position"] = df.apply(lambda r: find_nearest(chr_to_bpPositions[r["ChrA"]][chr_to_bpPositions[r["ChrA"]]>(r["EndA"])], r["EndA"]), axis=1)

            # get the 5' sequence (from one position after the closest breakpoint to the position before the breakpoint)
            df["5'sequence"] = df.apply(lambda r: chr_to_refSeq[r["ChrA"]][r["closest_5'breakpoint_position"]:r["StartA"]-1], axis=1)

            # get the 3' sequence (from the position after End to the position before the closest breakpoint)
            df["3'sequence"] = df.apply(lambda r: chr_to_refSeq[r["ChrA"]][r["EndA"]:r["closest_3'breakpoint_position"]-1], axis=1)

            # get the deleted sequence (from the start to the end)
            df["deleted_sequence"] = df.apply(lambda r: chr_to_refSeq[r["ChrA"]][r["StartA"]-1:r["EndA"]], axis=1)

            # change the chromosome seq in the sequence 
            for I, (chrA, seq5, seq3, del_seq) in enumerate(df[["ChrA", "5'sequence", "3'sequence", "deleted_sequence"]].values):
                print("copy-paste-insertion %i.."%I)

                # all seq
                ref_seq = seq5+del_seq+seq3

                # conformation in the rearranged chromosome
                rearranged_seq = seq5+seq3

                # check that the rearranged seq appears once in the genome and the ref seq in the ref genome. And they do not cross.
                chrA_refSeq = chr_to_refSeq[chrA]
                if not(chrA_refSeq.count(ref_seq)==1 and chrA_refSeq.count(rearranged_seq)==0 and all_rearranged_chromosomes_together.count(rearranged_seq)==1 and all_rearranged_chromosomes_together.count(ref_seq)==0): raise ValueError("The sequence is not unique")

                # go through each chrom of the rearranged seqs
                for chrom in chr_to_rearrangedSeq.keys():

                    # get the rearranged sequence
                    seq = cp.deepcopy(chr_to_rearrangedSeq[chrom])

                    # if the rearrangement sequence is in this chromosome, change it
                    if rearranged_seq in seq: 

                        # update the chr_to_rearrangedSeq so that it contains the reference sequence (copied)
                        chr_to_rearrangedSeq[chrom] = seq.replace(rearranged_seq, ref_seq)
                        break

            # get the rearranged genome into the file
            seq_records_list = [SeqRecord(Seq(seq), id=chrom, name=chrom, description=chrom) for chrom, seq in chr_to_rearrangedSeq.items()]

            # write the unmodified one
            print("writing")
            run_cmd("cp %s %s.tmp"%(rearranged_genome, rearranged_genome_unmodified_tmp))
            os.rename("%s.tmp"%rearranged_genome_unmodified_tmp, rearranged_genome_unmodified_tmp)

            # write the modified genome
            SeqIO.write(seq_records_list, rearranged_genome, "fasta")

            # write the modified genome
            os.rename(rearranged_genome_unmodified_tmp, rearranged_genome_unmodified)

        else: print("the insertions have already been modified")


def get_bed_df_not_overlapping_with_translocations_allChromARM(target_regions_bed, translocations_file, outdir, chr_to_len):

    """This function takes a bed file with target regions and subtracts the regions of the translocations file that are affected. It returns this bed file where translocations can be placed"""

    # get the interesting regions
    target_regions_df = pd.read_csv(target_regions_bed, sep="\t", header=None, names=["chromosome", "start", "end"])
    interesting_chromosomes = set(target_regions_df.chromosome)

    # get the bed df from the translocations
    regions_with_tra_bed_df, nSVs = get_affected_region_bed_for_SVdf(translocations_file, "translocations", interesting_chromosomes, first_position_idx=1, translocations_type="whole_chromosomes", chr_to_len=chr_to_len)

    # write the bed with the regions with SV
    regions_with_tra_bed = "%s/regions_with_translocations_whole_arms.bed"%outdir
    regions_with_tra_bed_df[["chromosome", "start", "end"]].to_csv(regions_with_tra_bed, sep="\t", index=False, header=False)

    # get the regions in all_regions_bed that are not in regions_with_SV_bed
    regions_without_tra_bed = "%s/noTRA_regions.bed"%outdir

    run_cmd("%s subtract -a %s -b %s > %s"%(bedtools, target_regions_bed, regions_with_tra_bed, regions_without_tra_bed))

    return regions_without_tra_bed

def check_consistency_of_svtype_to_svDF(svtype_to_svDF, all_chromosomes, chr_to_len):

    """Checks whether any of the breakpoints overlap with the others and reports those vars that do"""

 
    print("checking consistency of svtype to svDF")
    df_bed_allRegions = pd.DataFrame()

    # go through each df
    for svtype, svDF in svtype_to_svDF.items():

        if type(svDF)==str: svDF = pd.read_csv(svDF, sep="\t")

        svDF = svDF.set_index("ID", drop=False)

        for varID, sv_series in svDF.iterrows():

            # define series as df
            sv_series_df = pd.DataFrame({0 : sv_series}).transpose()

            # get the bed of this var
            sv_bed, nSVs = get_affected_region_bed_for_SVdf(sv_series_df, svtype, all_chromosomes, chr_to_len=chr_to_len, translocations_type="breakpoint_pos", first_position_idx=1)
            sv_bed["varID"] = varID
            sv_bed["svtype"] = svtype

            # get if there is any overlap with df_bed_allRegions
            regions_overlapping = df_bed_allRegions.apply(lambda rprevious: any(sv_bed.apply(lambda rnew: get_is_overlapping_query_vs_target_region(rprevious, rnew), axis=1)), axis=1)

            # if there is any region matching with the previous, continue, if not, keep
            if any(regions_overlapping): 
                # raise error if they do not come both from simulation
                if not (all(svDF.loc[{varID}, "ID"].apply(lambda x: "_sim_" in x)) and "_sim_" in varID):

                    print("%s has these overlapping regions:\n"%varID, df_bed_allRegions[regions_overlapping])
                    print("The actual var is \n", svDF.loc[varID, svtype_to_fieldsDict[svtype]["all_fields"]],"\n")
                    raise ValueError("there are overlapping regions")

            # add the bed to the regions matching
            df_bed_allRegions = df_bed_allRegions.append(sv_bed)


def generate_rearranged_genome_from_svtype_to_svDF(reference_genome, svtype_to_svDF, outdir, replace=False):

    """This function generates a rearranged genome for the provided svDFs, writing their fileds under outdir"""

    make_folder(outdir)

    # define the rearranged genome, and generated if not already done
    rearranged_genome = "%s/rearranged_genome.fasta"%outdir
    rearranged_genome_finalFile = "%s.performed"%(rearranged_genome)

    if file_is_empty(rearranged_genome_finalFile) or replace is True:

        # initialize a cmd to create the simulated genome
        targetSV_cmd = "%s --input_genome %s --output_genome %s"%(create_targeted_simulatedSVgenome_R, reference_genome, rearranged_genome)

        for svtype, svDF in svtype_to_svDF.items():

            # shift the insertions by 15 bp so that they are not at the beginning of the chrom
            if svtype=="insertions": svDF["StartA"] = svDF["StartA"] + 10

            # keep only the interesting svs
            svDF = svDF[svtype_to_fieldsDict[svtype]["all_fields"]]

            # write file
            svfile = "%s/%s.tab"%(outdir, svtype)
            svDF.to_csv(svfile, sep="\t", header=True, index=False)

            # add the generation of SVs into targetSV_cmd
            targetSV_cmd += " --%s_file %s"%(svtype, svfile)

        # run the cmd
        #std_rearranging_genome = "%s/simulation_std.txt"%outdir
        std_rearranging_genome = "stdout"

        if std_rearranging_genome!="stdout": run_cmd("%s > %s 2>&1"%(targetSV_cmd, std_rearranging_genome))
        else: run_cmd(targetSV_cmd)

        # transform the cut-and-paste insertions to copy-and-paste, whenever necessary
        insertions_file = "%s/insertions.tab"%outdir
        transform_cut_and_paste_to_copy_and_paste_insertions(reference_genome, rearranged_genome, insertions_file, svtype_to_svDF)

        # edit the insertions
        insertions_file = "%s/insertions.tab"%outdir
        rewrite_insertions_uniformizedFormat_simulateSV(insertions_file)

        # rewrite the variants so that they are optimal for comparison 
        translocations_file = "%s/translocations.tab"%outdir
        if not file_is_empty(translocations_file): rewrite_translocations_uniformizedFormat_simulateSV(translocations_file, reference_genome)

        # write a file that indicates that this has finsihed
        open(rearranged_genome_finalFile, "w").write("finsihed")

def get_random_svtype_to_svDF(reference_genome, mitochondrial_chromosome, outdir, nvars=200, replace=False, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}, check_random_genome_generation=False):

    """This function generates nvars into the reference genome splitting by gDNA and mtDNA with files written under outdir. It returns the randomly drawn variants and no-genome"""


    print("generating random simulations")

    # initialize a df that will contain the randomly-simulated vars
    random_svtype_to_svDF = {svtype : pd.DataFrame() for svtype in svtypes}

    # define the different types of chromosomes. Note that the mtDNA chromosomes will be simulated appart
    all_chromosomes = {s.id for s in SeqIO.parse(reference_genome, "fasta")}
    if mitochondrial_chromosome!="no_mitochondria": mtDNA_chromosomes = set(mitochondrial_chromosome.split(","))
    else: mtDNA_chromosomes = set()
    gDNA_chromosomes = all_chromosomes.difference(mtDNA_chromosomes)

    # map the chromosome to the length
    chrom_to_len = {s.id : len(s.seq) for s in SeqIO.parse(reference_genome, "fasta")}

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

        # define a bed file with all the regions
        all_regions_bed_df = pd.DataFrame({chrom: {"start":1, "end":chrom_to_len[chrom]} for chrom in chroms}).transpose()
        all_regions_bed_df["chromosome"] = all_regions_bed_df.index
        all_regions_bed_df = all_regions_bed_df[["chromosome", "start", "end"]]
        all_regions_bed = "%s/all_regions_index1.bed"%genome_outdir
        all_regions_bed_df.to_csv(all_regions_bed, sep="\t", header=False, index=False)

        # simulate random SVs into regions without previous SVs 
        random_sim_dir = "%s/random_SVs"%genome_outdir

        #### GET THE RANDOM INS,INV,DEL,TRA ####

        if any([file_is_empty("%s/%s.tab"%(random_sim_dir, svtype)) for svtype in {"insertions", "deletions", "translocations", "inversions", "tandemDuplications"}]) or replace is True:

            print("generating random SVs")

            # make and delete the folder
            delete_folder(random_sim_dir); make_folder(random_sim_dir)

            # get the cmd of the simulation
            randomSV_cmd = "%s --input_genome %s --outdir %s --regions_bed %s"%(create_random_simulatedSVgenome_R, genome_file, random_sim_dir, all_regions_bed)

            # add the number of each SV that should be added
            svtype_to_arg = {"insertions":"number_Ins", "deletions":"number_Del", "inversions":"number_Inv", "translocations":"number_Tra", "tandemDuplications":"number_Dup"}
        
            for svtype, arg in svtype_to_arg.items(): 
                if svtype not in svtype_to_arg: continue

                randomSV_cmd += " --%s %i"%(arg, vars_to_simulate)

            # run the random simulation
            #std_rearranging_genome = "%s/simulation_std.txt"%random_sim_dir
            std_rearranging_genome = "stdout"
            if std_rearranging_genome!="stdout": run_cmd("%s > %s 2>&1"%(randomSV_cmd, std_rearranging_genome))
            else: run_cmd(randomSV_cmd)

            # edit the translocations so that the balanced ones are sorted
            translocations_file = "%s/translocations.tab"%random_sim_dir
            if file_is_empty(translocations_file): open(translocations_file, "w").write("\t".join(["Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB", "Balanced", "BpSeqA", "BpSeqB"])) # this needs to be 

            # edit the insertions 
            insertions_file = "%s/insertions.tab"%random_sim_dir
            rewrite_insertions_uniformizedFormat_simulateSV(insertions_file)

        ########################################

        # add the simulations into random_svtype_to_svDF
        for svtype in random_svtype_to_svDF.keys():
            svDF = random_svtype_to_svDF[svtype]

            # get the new sv
            new_svDF = pd.read_csv("%s/%s.tab"%(random_sim_dir, svtype), sep="\t")
            new_svDF = new_svDF[[c for c in new_svDF.keys() if "BpSeq" not in c]]

            # add the name
            new_svDF["ID"] = new_svDF.Name + "_sim_%s"%type_genome

            # append 
            random_svtype_to_svDF[svtype] = svDF.append(new_svDF, sort=True)


    ####### test that you can insert the randomly simulated variants into the genome #######

    if check_random_genome_generation is True:

        # get the outdir
        outdir_randomVars_rearranging_genome = "%s/randomVars_rearranging_genome"%outdir; make_folder(outdir_randomVars_rearranging_genome)

        # get the rearranged genome
        generate_rearranged_genome_from_svtype_to_svDF(reference_genome, random_svtype_to_svDF, outdir_randomVars_rearranging_genome, replace=replace)

    ########################################################################################

    return random_svtype_to_svDF


def target_svDFseries_overlaps_bed_df(svDFseries, df_bed, svtype, first_position_idx, translocations_type, chr_to_len):

    """Returns whether a series of an svDF overlaps any of the regions in df_bed"""

    # get the interesting chromosomes
    interesting_chromosomes = set(chr_to_len)

    # get svDFseries as df
    svDF = pd.DataFrame({0: svDFseries}).transpose()

    # get the bed
    sv_bed = get_affected_region_bed_for_SVdf(svDF, svtype, interesting_chromosomes, first_position_idx=first_position_idx, translocations_type=translocations_type, chr_to_len=chr_to_len)[0]

    # get if there is any overlap
    regions_overlapping = df_bed.apply(lambda rprevious: any(sv_bed.apply(lambda rnew: get_is_overlapping_query_vs_target_region(rprevious, rnew), axis=1)), axis=1)

    return any(regions_overlapping)

def rearrange_genomes_simulateSV(reference_genome, outdir, replace=False, nvars=50, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulated_svtype_to_svfile={}, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}):

    """Runs a simulation of nvars SVs of each type into the reference genome. It will be sepparated by gDNA and mtDNA. mtDNA will only include 5% of the gDNA vars. Everything is written to outdir. simulated_svtype_to_svfile is a dictionary that maps each svtype to a file with it. This function will insert these SVs plus the remaining ones up to nvars (which will be placed randomly in the remaining spots on the genome). Note that only balanced translocations will be simulated, as unbalanced translocations are hard to bechmark based on coverage. 

    Keep in mind that all RSVSim are 1-based coordinates"""


    # change the simulated simulated_svtype_to_svfile to not include translocations
    #simulated_svtype_to_svfile = {svtype : svfile for svtype, svfile in simulated_svtype_to_svfile.items() if svtype!="translocations"} # debug    

    # only random var
    #simulated_svtype_to_svfile = {}

    # map each chrom to a len
    chr_to_len = get_chr_to_len(reference_genome)

    # check that the vars provided are consistent
    real_svtype_to_svDF = {svtype : pd.read_csv(file, sep="\t") for svtype, file in simulated_svtype_to_svfile.items()}
    check_consistency_of_svtype_to_svDF(real_svtype_to_svDF, set(chr_to_len), chr_to_len)

    # define the final outdirs
    final_simulated_SVs_dir = "%s/final_simulated_SVs"%(outdir); 
    final_rearranged_genome = "%s/rearranged_genome.fasta"%final_simulated_SVs_dir
    final_rearranged_genome_finalFile = "%s.performed"%(final_rearranged_genome)


    #if file_is_empty(final_rearranged_genome_finalFile) or replace is True:
    if True:

        # get random simulations for number of variants (twice the nvars if it is for future merge with simulated_svtype_to_svfile)
        if len(simulated_svtype_to_svfile)==0: random_nvars = nvars
        else: random_nvars = 2*nvars
        random_svtype_to_svDF = get_random_svtype_to_svDF(reference_genome, mitochondrial_chromosome, outdir, nvars=random_nvars, replace=replace, svtypes=svtypes)

        # get a df with all the bed regions of the real vars
        all_real_SVs_bed_df = pd.concat([get_affected_region_bed_for_SVdf(svDF, svtype, set(chr_to_len), first_position_idx=1, translocations_type="breakpoint_pos", chr_to_len=chr_to_len)[0] for svtype, svDF in real_svtype_to_svDF.items()])

        # initialize the final_svtype_to_svDF, which will contain vars for both real and simulated regions 
        final_svtype_to_svDF = {}

        # add random vars to final_svtype_to_svDF if they do not overlap with the bed regions
        svtypes = svtypes
        #svtypes = svtypes.difference({"translocations"})
        for svtype in svtypes:
            print("generating %s"%svtype)

            # define the real svtype
            if svtype in real_svtype_to_svDF: real_svDF = real_svtype_to_svDF[svtype]
            else: real_svDF = pd.DataFrame(columns=svtype_to_fieldsDict[svtype]["all_fields"])

            # define the random svtypes
            random_svDF = random_svtype_to_svDF[svtype]

            # add to the random SVtype if it overlaps any real svtypes
            random_svDF["overlaps_realSV"] = random_svDF.apply(lambda r: target_svDFseries_overlaps_bed_df(r, all_real_SVs_bed_df, svtype, 1, "breakpoint_pos", chr_to_len), axis=1)


            # translocations deserve a more complex processing
            if svtype=="translocations":

                ######### CHECK THAT THE REAL SVs are ready to INSERT #########

                print("WARING: RSVSim cannot rearrange-multiple times a chromosome, so that maybe the 'realVars' are not ready to insert. If RSVSim fails under this warning it may be because of this. And it would be necessary to implement a filtering of real vars taking this into account")

                ###############################################################

                #

            random_svDF = random_svDF[~random_svDF.overlaps_realSV]

            # get the concatenated df
            svDF = real_svDF.append(random_svDF, sort=True).iloc[0:nvars]

            # add to the final set
            final_svtype_to_svDF[svtype] = svDF

        # check the consistency of the genome
        check_consistency_of_svtype_to_svDF(final_svtype_to_svDF, set(chr_to_len), chr_to_len)

        # get the rearranged genome and simulations
        print("rearranging genome with real + random SVs")
        make_folder(final_simulated_SVs_dir)
        generate_rearranged_genome_from_svtype_to_svDF(reference_genome, final_svtype_to_svDF, final_simulated_SVs_dir, replace=True)

        bzvdmbmvznbnmvvdbz



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
    print("calculating shortest paths in genome graph from %i end positions to %i positions. This may take a lot if there are a lot of end positions"%(len(chrom_end_positions), len(all_positions)))
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

    """This function takes a df with windows of the genome and adds the gc content for each window, writing a file under gcontent_outfile. It will only do those that have been already measured"""

    print("Getting GC content")

    if file_is_empty(gcontent_outfile) or replace is True:

        # define the initial index
        initial_index = list(df_windows.index)

        # resort
        df_windows = df_windows.sort_values(by=["chromosome", "start", "end"]).set_index(["chromosome", "start", "end"], drop=False)

        print("getting GC content for %i new windows"%len(df_windows))

        # get the GC content file for each position
        gc_content_outfile_perPosition = generate_nt_content_file(genome, replace=replace, target_nts="GC")
        gc_df = pd.read_csv(gc_content_outfile_perPosition, sep="\t")[["chromosome", "position", "is_in_GC"]].sort_values(by=["chromosome", "position"])


        # define a df where each position is one row and it has the start_window as an add
        df_windows["length"] = df_windows.end - df_windows.start
        positions = make_flat_listOflists(list(df_windows.apply(lambda r: list(range(r["start"], r["end"])), axis=1)))
        start_windows = make_flat_listOflists(list(df_windows.apply(lambda r: [r["start"]]*r["length"], axis=1)))
        end_windows = make_flat_listOflists(list(df_windows.apply(lambda r: [r["end"]]*r["length"], axis=1)))
        chromosomes = make_flat_listOflists(list(df_windows.apply(lambda r: [r["chromosome"]]*r["length"], axis=1)))
        df_positions = pd.DataFrame({"position":positions, "chromosome":chromosomes, "start_window":start_windows, "end_window":end_windows})

        # add the positions to the gc df
        gc_df = gc_df.merge(df_positions, on=["chromosome", "position"], how="right")        

        # calculate the GC content and add to df
        window_to_gc = gc_df[["chromosome", "start_window", "end_window", "is_in_GC"]].groupby(["chromosome", "start_window", "end_window"]).mean()["is_in_GC"]
     
        # get into df_windows
        df_windows["GCcontent"] = list(window_to_gc.loc[df_windows.index])

        # at the end save the df windows
        df_windows.index = initial_index
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

        outfile = "%s/coverage_modelling_%s.pdf"%(outdir, type_genome)

        if file_is_empty(outfile) or replace is True:

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
    final_list = list(df.loc[initial_index, "cov_rel_to_predFromFeats"])

    return final_list

def get_target_region_row(r, region, breakpoint_positions, maxPos, max_relative_len_neighbors=2):

    """This function takes a row of a df with windows and returns a row of the 'region', which can be 5 or 3. It will span to the next breakpoint or the max_relative_len_neighbors"""

    # sort the breakpoint positions
    breakpoint_positions = np.array(sorted(breakpoint_positions))

    # get the maximum region length
    max_region_length = (r["end"]-r["start"])*max_relative_len_neighbors

    # define the minimum region length
    min_region_len = 100

    # if the start is 1 , set the 5' region to the 3' region
    if region=="5" and r["start"]==1: region="3"

    # get the 5' region
    if region=="5":

        # the end is the start of the region
        end = r["start"]-1

        ##### define the start #####

        # get the breakpoint positions that are before the end, if any
        bps_before_end = breakpoint_positions[breakpoint_positions<(end-min_region_len)]

        if len(bps_before_end)>0: previous_bp = sorted(bps_before_end)[-1]
        else: previous_bp = 0

        # only keep the previous_bp as start if the length is not above max_region_length
        length = end - previous_bp
        if length<=max_region_length: start = previous_bp
        else: start = end - max_region_length

        ############################

    # get the 3' region
    elif region=="3":

        # the start is the end of the region
        start = r["end"]+1

        ##### define the end #####

        # get the breakpoint positions that are after the start, if any
        bps_after_start = breakpoint_positions[breakpoint_positions>(start+min_region_len)]

        if len(bps_after_start)>0: next_bp = sorted(bps_after_start)[0]
        else: next_bp = maxPos

         # only keep the next_bp as start if the length is not above max_region_length
        length = next_bp - start
        if length<=max_region_length: end = next_bp
        else: end = start + max_region_length

        ##########################

    # define general names
    region_name = "%s_region"%region

    # return a series of all important fields
    return pd.Series({"chromosome":r["chromosome"], "start":start, "end":end, "region_name":region_name})

def get_df_with_coverage_per_windows_relative_to_neighbor_regions(df_windows, bed_windows_prefix, reference_genome, sorted_bam, df_clove, median_coverage, replace=True, run_in_parallel=True, delete_bams=True):

    """Takes a df with windows of the genome and returns it with the coverage and the relative to the genome. It returns a df with several relative coverage measures."""

    print("getting coverage relative to neighbors")

    # get the coverage to len
    chrom_to_maxPos= {seq.id : len(seq.seq)-1 for seq in SeqIO.parse(reference_genome, "fasta")}

    # map the chromosome to the positions with breakpoints
    chrom_to_bpPositions_CHROM = dict(df_clove.groupby("#CHROM").apply(lambda df_c: get_uniqueVals_df(df_c[["POS"]])))
    chrom_to_bpPositions_CHR2 = dict(df_clove.groupby("CHR2").apply(lambda df_c: get_uniqueVals_df(df_c[["START", "END"]]).difference({-1})))
    chrom_to_bpPositions = {chrom:set() for chrom in chrom_to_maxPos}
    for chromDict in [chrom_to_bpPositions_CHROM, chrom_to_bpPositions_CHR2]:
        for chrom, bpPositions in chromDict.items():
            chrom_to_bpPositions[chrom].update(bpPositions)

    # initialize a df windows for the 5' and 3' regions
    df_windows = df_windows.sort_values(by=["chromosome", "start", "end"]).drop_duplicates(subset=["chromosome", "start", "end"])
    all_sorted_chromosomes = list(df_windows.chromosome)
    all_df_windows = cp.deepcopy(df_windows)

    # go through each region and get a coverage df
    for region in ["target", "5", "3"]: 

        if region=="target":
            df_region = df_windows
            df_region["region_name"] = "target_region"

        else:

            # get a df with the regions
            df_region = df_windows.apply(lambda r: get_target_region_row(r, region, chrom_to_bpPositions[r["chromosome"]], chrom_to_maxPos[r["chromosome"]]), axis=1)

        # get the coverage df
        bed_file = "%s.%s.bed"%(bed_windows_prefix, region)
        df_region.to_csv(bed_file, sep="\t", header=True, index=False)

        coverage_df = get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, bed_file, replace=replace, run_in_parallel=run_in_parallel, delete_bams=delete_bams)

        # add to df region
        df_region = df_region.merge(coverage_df[["chromosome", "start", "end", "mediancov_1"]], on=["chromosome", "start", "end"], validate="many_to_one", how="left").sort_values(by=["chromosome", "start", "end"])

        # check that the regions are the same
        if all_sorted_chromosomes!=list(df_region.chromosome): raise ValueError("something went wrong with the merging of the line above")

        # add to all df
        median_coverage_list = list(df_region["mediancov_1"])
        all_df_windows["%s_coverage"%region] = median_coverage_list

        # add relative parms
        all_df_windows["relative_coverage_%s"%region] = all_df_windows["%s_coverage"%region]/median_coverage

    # get the coverage relative to the regions
    for region in ["5", "3"]: all_df_windows["coverage_rel_to_%s"%region] = all_df_windows["target_coverage"]/all_df_windows["%s_coverage"%region]

    # get estimate of both relative coverages
    all_df_windows["mean_rel_coverage_to_neighbor"] = (all_df_windows["coverage_rel_to_5"]+all_df_windows["coverage_rel_to_3"])/2
    all_df_windows["closestTo1_rel_coverage_to_neighbor"] = all_df_windows.apply(lambda r: find_nearest([r["coverage_rel_to_5"], r["coverage_rel_to_3"]], 1), axis=1)

    return all_df_windows

def get_clove_output_with_coverage(outfile_clove, reference_genome, sorted_bam, median_coverage, replace=False, run_in_parallel=True, delete_bams=True):

    """Takes the output of clove and adds the coverage of the TAN, DEL and INS-like features or -1. The added coverage is the closest-to-1 relative coverage of the the 5' and 3' regions to each region. Each 5' or 3' regions spans from the end of the target region to the next breakpoint, given that the region is maximum 2x the target region.  """

    # first load clove into a df
    df_clove = get_clove_output(outfile_clove)

    # define the SVtypes
    all_svtypes = set(df_clove.SVTYPE)
    tanDEL_svtypes = {"TAN", "DEL"}
    ins_svtypes = {"CID", "CIT", "DUP", "TRA"}
    remaining_svtypes = all_svtypes.difference(tanDEL_svtypes.union(ins_svtypes))

    if len(df_clove)>0:

        # get the regions of TANDEL and INS that should be renamed
        df_TANDEL = df_clove[df_clove.SVTYPE.isin(tanDEL_svtypes)].rename(columns={"#CHROM":"chromosome", "POS":"start", "END":"end"})[["chromosome", "start", "end"]]

        df_INS = df_clove[df_clove.SVTYPE.isin(ins_svtypes)].rename(columns={"CHR2":"chromosome", "START":"start", "END":"end"})[["chromosome", "start", "end"]]

        df_TANDELINS = df_TANDEL.append(df_INS)
        
        if len(df_TANDELINS)>0:

            # define a bed with these regions
            bed_TANDELINS_regions = "%s.TANDELINS.bed"%outfile_clove

            # define a coverage df
            coverage_df = get_df_with_coverage_per_windows_relative_to_neighbor_regions(df_TANDELINS, bed_TANDELINS_regions, reference_genome, sorted_bam, df_clove, median_coverage, replace=replace, run_in_parallel=run_in_parallel, delete_bams=delete_bams)
         
        else: coverage_df = pd.DataFrame(columns=['chromosome', 'start', 'end', 'target_coverage', 'relative_coverage_target', '5_coverage', 'relative_coverage_5', '3_coverage', 'relative_coverage_3', 'coverage_rel_to_5', 'coverage_rel_to_3', 'mean_rel_coverage_to_neighbor', 'closestTo1_rel_coverage_to_neighbor'])


        # initialize merge
        merged_df = pd.DataFrame()

        svtypeName_to_coordFields = {"tanDel":["#CHROM", "POS", "END"], "remaining":["#CHROM", "POS", "END"], "ins":["CHR2", "START", "END"], }
        # merge by SVtypes
        for svtypeName, svtypes in [("tanDel", tanDEL_svtypes), ("ins", ins_svtypes), ("remaining", remaining_svtypes)]:

            # get the coord fields
            coord_fields = svtypeName_to_coordFields[svtypeName]

            # interesting df clove
            svtype_df_clove = df_clove[df_clove.SVTYPE.isin(svtypes)]

            # get the merged df
            if svtypeName in {"tanDel", "ins"}: 

                # the df with the actual coverage
                merged_df_svtype = svtype_df_clove.merge(coverage_df, left_on=coord_fields, right_on=["chromosome", "start", "end"], validate="one_to_one", how="left")

            else: 

                # initialize as clove
                merged_df_svtype = svtype_df_clove

                # add nan fields
                for field in set(merged_df.columns).difference(set(merged_df_svtype.columns)): merged_df_svtype[field] = np.nan

            # add
            merged_df = merged_df.append(merged_df_svtype, sort=True)

        # change types of fields
        merged_df["POS"] = merged_df.POS.apply(get_int)
        merged_df["END"] = merged_df.END.apply(get_int)
        merged_df["START"] = merged_df.START.apply(get_int)

        # set index
        merged_df.index = list(range(len(merged_df)))

        return merged_df 

    else: return pd.DataFrame()

def get_covfilter_cloveDF_row_according_to_SVTYPE(r, max_rel_coverage_to_consider_del=0.1, min_rel_coverage_to_consider_dup=1.9, coverage_field="mean_rel_coverage_to_neighbor"):

    # define a function that takes a row of the dataframe and does the filterinf can be any that is in coverage

    potential_duplication_fields = {"CID", "CIT", "DUP", "TRA", "TAN"} # note that this includes insertions and TAN 

    if r["SVTYPE"]=="DEL":

        if r[coverage_field]<=max_rel_coverage_to_consider_del: return "PASS"
        else: return "FAIL" 

    elif r["SVTYPE"] in potential_duplication_fields:

        if r[coverage_field]>=min_rel_coverage_to_consider_dup: return "PASS"
        else: return "FAIL" 

    else: return "PASS"

def get_bedpeDF_for_clovebalTRA_5with5_or_3with3(df_clove, tol_bp=50):

    """Takes a clove df and returns a dataframe with events that could be balanced translocations, in the sense that they are ITX1 and ITX2 events very close (according to tol_bp) """

    # get the ITX dataframes
    df_ITX1 = df_clove[df_clove.SVTYPE=="ITX1"]
    df_ITX2 = df_clove[df_clove.SVTYPE=="ITX2"]

    # initialize a dict
    data_dict = {}

    # go through each combination
    for I1 in df_ITX1.index:
        I1s = df_ITX1.loc[I1]

        for I2 in df_ITX2.index:
            I2s = df_ITX2.loc[I2]

            # ask whether this combination can be close
            if I1s["#CHROM"]==I2s["#CHROM"] and I1s["CHR2"]==I2s["CHR2"] and abs(I1s["POS"]-I2s["POS"])<=tol_bp and abs(I1s["END"]-I2s["END"])<=tol_bp:

                # keep in a way that it is in the same order as in the clove VCF
                data_dict[(I1, I2)] = {"ChrA":I1s["#CHROM"], "StartA":0, "EndA":I1s["POS"], "ChrB":I1s["CHR2"], "StartB":0, "EndB":I1s["END"], "Balanced":True} # all the starts are 0s

    df = pd.DataFrame(data_dict).transpose()

    return df

def get_bedpeDF_for_clovebalTRA_5with3_or_3with5_INVTXbreakpoints(df_clove, chr_to_len, tol_bp=50):

    """Takes a clove df and returns a dataframe with events that could be balanced translocations, in the sense that they are INVTX1 and INVTX2 events very close (according to tol_bp) """

    # get the INVTX2 dataframes
    df_INVTX1 = df_clove[df_clove.SVTYPE=="INVTX1"]
    df_INVTX2 = df_clove[df_clove.SVTYPE=="INVTX2"]

    # initialize a dict
    data_dict = {}

    # go through each combination
    for I1 in df_INVTX1.index:
        I1s = df_INVTX1.loc[I1]

        for I2 in df_INVTX2.index:
            I2s = df_INVTX2.loc[I2]

            # ask whether this combination can be close
            if I1s["#CHROM"]==I2s["#CHROM"] and I1s["CHR2"]==I2s["CHR2"] and abs(I1s["POS"]-I2s["POS"])<=tol_bp and abs(I1s["END"]-I2s["END"])<=tol_bp:

                # keep in a way that it is in the same order as in the clove VCF
                data_dict[(I1, I2)] = {"ChrA":I1s["#CHROM"], "StartA":0, "EndA":I1s["POS"], "ChrB":I1s["CHR2"], "StartB":I1s["END"], "EndB":chr_to_len[I1s["CHR2"]], "Balanced":True} # all the starts are 0s

    df = pd.DataFrame(data_dict).transpose()

    return df

def get_bedpe_for_clovebalTRA_5with3(r, chr_to_len):

    """Takes a row of the df_balTRA_5with3 df and returns a bedpe row, sorted"""

    return pd.Series({"ChrA":r["#CHROM"], "StartA":0, "EndA":r["POS"], "ChrB":r["CHR2"], "StartB":r["END"], "EndB":chr_to_len[r["CHR2"]]-1, "Balanced":True})

def write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, fileprefix, reference_genome, sorted_bam, tol_bp=50, replace=False, svtypes_to_consider={"insertions", "deletions", "inversions", "translocations", "tandemDuplications", "remaining"}, run_in_parallel=False, define_insertions_based_on_coverage=False):

    """Takes a clove dataframe and writes the different SV into several files, all starting with fileprefix. it returns a dict mapping each SVtype to the file with the bed or bedpe containing it. tol_bp indicates the basepairs that are considered as tolerated to be regarded as 'the same event' 

    consider_TANDEL indicates whether to write, it requires the coverage_FILTER to PASS.

    only balanced translocations are considered.

    Thesse are the SVTYPE fields:

    DEL "Deletion"
    TAN "Tandem Duplication"
    INV "Inversion"
    INS "Insertion" 
    DUP "Complex Duplication"
    TRA "Complex Translocation"
    ##ALT=<ID=CIV,Description="Complex Inversion">
    ##ALT=<ID=CVT,Description="Complex Inverted Translocation">
    ##ALT=<ID=CVD,Description="Complex Inverted Duplication">
    ##ALT=<ID=CIT,Description="Complex Interchromosomal Translocation">
    ##ALT=<ID=CID,Description="Complex Interchromosomal Duplication">
    ##ALT=<ID=IVT,Description="Complex Inverted Interchromosomal Translocation">
    ##ALT=<ID=IVD,Description="Complex Inverted Interchromosomal Duplication">

    define_insertions_based_on_coverage indicates whether to filter the insertions based on coverage. If so, it will set as copied insertions those that have a coverage above 

    """
    print("getting SVs from clove")

    # initialize as a copy
    df_clove = cp.deepcopy(df_clove)

    # initialize the final dict
    svtype_to_svfile = {}

    # initialize the considered idxs
    df_clove.index = list(range(len(df_clove)))
    considered_idxs = []

    # map each index to the ID
    cloveIDX_to_ID = dict(df_clove.ID)

    # GENERAL THINGS #############

    # add length chromosome
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

    # define the svtypes that are straightforwardly classified
    cloveSVtypes_easy_classification = {"CID", "CIT", "DUP", "TRA", "CIV", "IVD"} # note that DEL and TAN are left out because they may not directly be assignable to DEL and TAN

    ###############################


    ###### TRANSLOCATIONS:  A segment from the 5’ or 3’ end of one chromosome A is exchanged with the 5’ or 3’ end of another chromosome B. ONLY BALANCED TRANSLOCATIONS ARE CONSIDERED#######
    if any([x in set(df_clove.SVTYPE) for x in {"ITX1", "ITX2", "IVD", "INVTX1", "INVTX2"}]) and "translocations" in svtypes_to_consider:

        # balanced translocations 5with5
        df_balTRA_5with5_or_3with3 = get_bedpeDF_for_clovebalTRA_5with5_or_3with3(df_clove, tol_bp=tol_bp) # here the index is not balanced
        considered_idxs += make_flat_listOflists(df_balTRA_5with5_or_3with3.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

        # balanced translocations 5with3 (these are the ones with an IVD field, assigned by clove)
        df_balTRA_5with3_IVD = df_clove[(df_clove.SVTYPE=="IVD") & ((df_clove.START - df_clove.END)<=tol_bp)].apply(lambda r: get_bedpe_for_clovebalTRA_5with3(r, chr_to_len), axis=1)
        considered_idxs += list(df_balTRA_5with3_IVD.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

        # balanced translocations 5with3 where there are two close INVTX breakpoints
        df_balTRA_5with3_INVTXbps = get_bedpeDF_for_clovebalTRA_5with3_or_3with5_INVTXbreakpoints(df_clove, chr_to_len, tol_bp=tol_bp).apply(lambda r: get_bedpe_for_clovebalTRA_5with3(r, chr_to_len), axis=1)
        considered_idxs += list(df_balTRA_5with3_INVTXbps.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

        # merge both
        df_balTRA_5with3 = df_balTRA_5with3_IVD.append(df_balTRA_5with3_INVTXbps)

        # merge together and add some fields
        important_fields = ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Balanced"]
        translocations_dfs =  [df_balTRA_5with5_or_3with3, df_balTRA_5with3]
        if any([len(d)>0 for d in translocations_dfs]): 

            # non empty df
            df_tra = pd.concat([d[important_fields] for d in translocations_dfs if all([f in d.keys() for f in important_fields])], sort=True)
            
            # forrmat the translocation table in the same way as simulateSV
            df_tra = df_tra.apply(lambda r: format_translocation_row_simulateSV(r, chr_to_len, start_pos=0), axis=1)[important_fields]
            df_tra["Balanced"] = df_tra.Balanced.apply(lambda x: str(x).upper())

        else: df_tra = pd.DataFrame(columns=important_fields)


        # add the ID to the df
        def get_ID(idx):
            if type(idx)==int: return cloveIDX_to_ID[idx]
            else: return "+".join([cloveIDX_to_ID[x] for x in idx])

        df_tra["ID"] = [get_ID(idx) for idx in df_tra.index]
        important_fields += ["ID"]

        # write, rewrite, and keep
        bedpe_translocations = "%s.translocations.bedpe.withBlancedINFO"%fileprefix
        df_tra[important_fields].to_csv(bedpe_translocations, sep="\t", header=True, index=False)

        # keep
        svtype_to_svfile["translocations"] = bedpe_translocations

        print("There are %i translocations"%len(df_tra))

    #############################

    ####### INVERSIONS ##########
    if "CIV" in set(df_clove.SVTYPE) and "inversions" in svtypes_to_consider:

        df_inversions = df_clove[df_clove.SVTYPE=="CIV"][["#CHROM", "POS", "END", "ID"]].rename(columns={"#CHROM":"Chr", "POS":"Start", "END":"End"})
        inversions_bed = "%s.inversions.bed"%fileprefix
        df_inversions.to_csv(inversions_bed, sep="\t", header=True, index=False)

        considered_idxs += list(df_inversions.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]
        svtype_to_svfile["inversions"] = inversions_bed

        print("There are %i inversions"%len(df_inversions))


    #############################

    ####### INSERTIONS ########
    if any([x in set(df_clove.SVTYPE) for x in {"CID", "CIT", "DUP", "TRA"}]) and "insertions" in svtypes_to_consider:

        # get the df
        df_ins = df_clove[df_clove.SVTYPE.isin({"CID", "CIT", "DUP", "TRA"})][["#CHROM", "POS", "CHR2", "START", "END", "ID", "SVTYPE", "coverage_FILTER"]]

        # rename
        df_ins = df_ins.rename(columns={"#CHROM":"ChrB", "POS":"StartB", "CHR2":"ChrA", "START":"StartA", "END":"EndA"})
        
        # add the end. It is formated in a way that the insertion length is equivalent to the inserted fragment
        df_ins["EndB"] = df_ins.StartB + (df_ins.EndA - df_ins.StartA)

        # add whether it is copied
        if define_insertions_based_on_coverage is False:
            
            svtype_to_isCopied = {"CID":"TRUE", "CIT":"FALSE", "DUP":"TRUE", "TRA":"FALSE"}
            df_ins["Copied"] = df_ins.SVTYPE.apply(lambda x: svtype_to_isCopied[x])

        # based on the coverage field
        else: df_ins["Copied"] = (df_ins.coverage_FILTER=="PASS").apply(lambda x: str(x).upper())

        # write as bedpe
        important_fields = ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Copied", "ID"]
        bedpe_insertions = "%s.insertions.bedpe.withCopiedINFO"%fileprefix
        df_ins[important_fields].to_csv(bedpe_insertions, sep="\t", header=True, index=False)

        # keep
        svtype_to_svfile["insertions"] = bedpe_insertions
        considered_idxs += list(df_ins.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

        print("There are %i insertions, %i of which are copy-and-paste"%(len(df_ins), sum(df_ins.Copied=="TRUE")))

    ############################

    ###### DEL and TAN #######

    typeSV_to_tag = {"deletions":"DEL", "tandemDuplications":"TAN"}
    for typeSV, tag  in typeSV_to_tag.items():

        # check that it exists
        if tag in set(df_clove.SVTYPE) and typeSV in svtypes_to_consider:

            # write the file
            bed_filename = "%s.%s.bed"%(fileprefix, typeSV)
            df_svtype = df_clove[(df_clove.SVTYPE==tag) & (df_clove.coverage_FILTER=="PASS")].rename(columns={"#CHROM":"Chr", "POS":"Start", "END":"End"})
            df_svtype[["Chr", "Start", "End", "ID"]].to_csv(bed_filename, sep="\t", header=True, index=False)

            # keep
            svtype_to_svfile[typeSV] = bed_filename
            considered_idxs += list(df_svtype.index)

            print("There are %i %s"%(len(df_svtype), typeSV))

    # keep only df_clove that has not been already used, which should be done after each step
    df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

    ###########################


    # write the remaining events which are not easily assignable
    df_notAssigned = df_clove[~(df_clove.SVTYPE.isin(cloveSVtypes_easy_classification))]
    df_notAssigned_file = "%s.remaining.tab"%(fileprefix)
    df_notAssigned[["ID", "#CHROM", "POS", "CHR2", "START", "END", "SVTYPE"]].to_csv(df_notAssigned_file, sep="\t", header=True, index=False)
    svtype_to_svfile["remaining"] = df_notAssigned_file

    print("There are %i remaining SVs"%len(df_notAssigned))

    # at the end make sure that the considered idxs are unique
    if len(considered_idxs)!=len(set(considered_idxs)): 
        print(fileprefix, considered_idxs)
        raise ValueError("ERROR: Some clove events are assigned to more than one cathegory. Check the insertions and translocations calling")
        #print("WARNING: Some clove events are assigned to more than one cathegory. Check the insertions and translocations calling")

    # return the df_clove and the remaining SVs
    return df_clove, svtype_to_svfile

def merge_coverage_per_window_files_in_one(bamfile, bam_sufix=".coverage_per_window.tab"):

    """This function takes all files that start with bamfile and end with coverage_per_window, """

    print("merging coverage tables")

    # define prefixes
    bam_dir = get_dir(bamfile)
    fileprefix = get_file(bamfile) + bam_sufix

    # remove dirs
    dirs_to_remove = ["%s/%s"%(bam_dir, f) for f in os.listdir(bam_dir) if os.path.isdir("%s/%s"%(bam_dir, f)) and f.startswith(fileprefix) and len(os.listdir("%s/%s"%(bam_dir, f)))==0] 
    for f in dirs_to_remove: delete_folder(f)

    # unite files
    files_prefix = ["%s/%s"%(bam_dir, f) for f in os.listdir(bam_dir) if not file_is_empty("%s/%s"%(bam_dir, f)) and "temporary_file" not in f and f.startswith(fileprefix)]
    df_all = pd.concat([pd.read_csv(f, sep="\t") for f in files_prefix])

    # write into one
    integrated_file = bamfile+bam_sufix
    df_all.to_csv(integrated_file, sep="\t", header=True, index=False)

    # remove other files
    for f in files_prefix: 
        if f!=integrated_file: remove_file(f)


def run_gridssClove_given_filters(sorted_bam, reference_genome, working_dir, median_coverage, replace=True, threads=4, gridss_blacklisted_regions="", gridss_VCFoutput="", gridss_maxcoverage=50000, median_insert_size=250, median_insert_size_sd=0, gridss_filters_dict=default_filtersDict_gridss, tol_bp=50, run_in_parallel=True, max_rel_coverage_to_consider_del=0.2, min_rel_coverage_to_consider_dup=1.8, replace_FromGridssRun=False, define_insertions_based_on_coverage=False):

    """This function runs gridss and clove with provided filtering and parameters. This can be run at the end of a parameter optimisation process. It returns a dict mapping each SV to a table, and a df with the gridss.

    coverage_field is the field where clove is filtered to detect CNV. It can be relative_coverage or relative_coverage_dist_to_telomere.

    include_breakpoints_in_genomeGraph indicates whether breakpoints should be included in the genome graph

    type_coverage_to_filterTANDEL is a field that determines what type of coverage is used to filter deletions and tandem duplications in the CLOVE output:

    - mediancov_1 indicates that it should be the abslute and raw coverage
    - coverage_rel_to_predFromFeats means that it should be  based on the coverage relative to the predicted from the features
    - define_insertions_based_on_coverage indicates whether insertions should be defined as "Copied" if they have a coverage >min_rel_coverage_to_consider_dup relative to the nighbor regions



    """

    print("running gridss and clove with given parameter with %.2f min_rel_coverage_to_consider_dup"%min_rel_coverage_to_consider_dup)
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

    ###################################

    #################################################
    #################################################
    #################################################

    # run clove without checking filtering
    outfile_clove = "%s.clove.vcf"%(raw_bedpe_file)
    run_clove_filtered_bedpe(raw_bedpe_file, outfile_clove, sorted_bam, replace=replace_FromGridssRun, median_coverage=median_coverage, median_coverage_dev=1, check_coverage=False) #  REPLACE debug

    #replace_FromGridssRun = True # debug

    # add the filter of coverage to the clove output
    df_clove = get_clove_output_with_coverage(outfile_clove, reference_genome, sorted_bam, median_coverage, replace=replace_FromGridssRun, run_in_parallel=run_in_parallel, delete_bams=run_in_parallel)

    # define the coverage filtering based on the type_coverage_to_filterTANDEL
    df_clove["coverage_FILTER"] = df_clove.apply(lambda r: get_covfilter_cloveDF_row_according_to_SVTYPE(r, max_rel_coverage_to_consider_del=max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup=min_rel_coverage_to_consider_dup, coverage_field="mean_rel_coverage_to_neighbor"), axis=1)

    # annotated clove 
    fileprefix = "%s.structural_variants"%outfile_clove

    remaining_df_clove, svtype_to_SVtable = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, fileprefix, reference_genome, sorted_bam, tol_bp=tol_bp, replace=replace_FromGridssRun, svtypes_to_consider={"insertions", "deletions", "inversions", "translocations", "tandemDuplications", "remaining"}, run_in_parallel=run_in_parallel, define_insertions_based_on_coverage=define_insertions_based_on_coverage)

    # merge the coverage files in one
    merge_coverage_per_window_files_in_one(sorted_bam)

    return svtype_to_SVtable, df_gridss


###################################################################################################
###################################################################################################
###################################################################################################

def generate_tables_of_SV_between_genomes_gridssClove(query_genome, reference_genome, replace=False, threads=4, coverage=40, insert_size=250, read_lengths=[kb*1000 for kb in [0.5, 1, 1.5, 2]], error_rate=0.0, gridss_min_af=0.25):

    """Takes a bam file with aligned reads or genomes and generates calls, returning a dict that maps variation type to variants
    - aligner can be minimap2 or ngmlr.

    [0.5, 0.7, 0.9, 1]

    """


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

    return SV_dict, df_gridss

def get_is_matching_predicted_and_known_rows(rk, rp, equal_fields, approximate_fields, chromField_to_posFields, tol_bp=50, pct_overlap=0.75):

    """Takes a row of a knownID (rk) and a predictedID (rp) and returns a boolean indicating if they match. These rk and rp can be any dict-like structures that have the expected equal_fields and so."""

    # ask if the equal fields match
    equal_fields_match = all([rp[f]==rk[f] for f in equal_fields])

    # ask if the approximate fields match
    approximate_fields_match = all([abs(rp[f]-rk[f])<=tol_bp for f in approximate_fields])

    # ask if all the chromField_to_posFields overlap by more than pct_overlap
    overlapping_regions_list = [True] # stores tested regions overlap. It is initalized with a True so that if there are no chromField_to_posFields it is True
    for chromField, posFields in chromField_to_posFields.items(): 

        # define the regions
        chrom_k = rk[chromField]
        start_k = rk[posFields["start"]]
        end_k = rk[posFields["end"]]
        len_k = end_k-start_k

        chrom_p = rp[chromField]
        start_p = rp[posFields["start"]]
        end_p = rp[posFields["end"]]
        len_p = end_p-start_p

        # make sure that the end is after the start
        if len_k<=0 or len_p<=0: raise ValueError("start is after end")

        # get the length of the longest region
        longest_region_len = max([len_k, len_p])

        # define the boundaries of the overlap
        start_overlap_pos = max([start_p, start_k])
        end_overlap_pos = min([end_p, end_k])
        len_overlap = end_overlap_pos-start_overlap_pos

        # add that the region overalps by more than pct_overlap
        overlapping_regions_list.append( chrom_k==chrom_p and (len_overlap/longest_region_len)>=pct_overlap)

    regions_overlap = all(overlapping_regions_list)

    return  equal_fields_match and approximate_fields_match and regions_overlap

def get_SVbenchmark_dict(df_predicted, df_known, equal_fields=["Chr"], approximate_fields=["Start", "End"], chromField_to_posFields={}):

    """Takes dfs for known and predicted SVs and returns a df with the benchmark. approximate_fields are fields that have to overlap at least by tolerance_bp. It returns a dict that maps each of the benchmark fields to the value. pct_overlap is the percentage of overlap between each of the features in approximate_fields.

    chromField_to_posFields is a dict that maps each chromosome field to the start and end fields that need to be interrogated by pct_overlap"""

    # define the ID fields 
    if "Name" in df_known.keys(): known_IDfield = "Name"
    else: known_IDfield = "ID"

    predicted_IDfield = "ID"


    # get the predictedIDs as those that have the same equal_fields and overlap in all approximate_fields
    if len(df_predicted)>0: 

        df_known["predictedSV_IDs"] = df_known.apply(lambda rk: set(df_predicted[df_predicted.apply(lambda rp: get_is_matching_predicted_and_known_rows(rk, rp, equal_fields, approximate_fields, chromField_to_posFields), axis=1)][predicted_IDfield]), axis=1)

    else: df_known["predictedSV_IDs"] = [set()]*len(df_known)

    # calculate the length
    df_known["predictedSV_IDs_len"] = df_known["predictedSV_IDs"].apply(len)

    # check that there is only 0 or 1 IDs matching
    if any([x not in {0,1} for x in set(df_known.predictedSV_IDs_len)]): 
        print("WARNING: There are some predictedIDs that match more than one variant")
        print(df_known)

    # define a df with the predicted events
    df_known_matching = df_known[df_known.predictedSV_IDs_len>0]

    # define sets of IDspredictedSV_IDs_len
    all_known_IDs = set(df_known[known_IDfield])
    if len(df_predicted)>0: all_predicted_IDs = set(df_predicted[predicted_IDfield])
    else: all_predicted_IDs = set()

    true_positives_knownIDs = set(df_known_matching[known_IDfield])
    false_negatives_knownIDs = all_known_IDs.difference(true_positives_knownIDs)
    
    if len(df_known_matching)==0:  true_positives_predictedIDs = set()
    else: true_positives_predictedIDs = set.union(*df_known_matching.predictedSV_IDs)
    false_positives_predictedIDs = all_predicted_IDs.difference(true_positives_predictedIDs)

    # calculate stats
    TP = len(true_positives_knownIDs)
    TP_predictedIDs = len(true_positives_predictedIDs)
    FP = len(false_positives_predictedIDs)
    FN = len(false_negatives_knownIDs)
    nevents = len(all_known_IDs)
    if nevents==0: precision=1.0; recall=1.0
    else:
        if TP==0 and FP==0: precision =  0.0
        else: precision = TP/(TP + FP)
        recall = TP/(TP + FN)
        
    if precision<=0.0 or recall<=0.0: Fvalue = 0.0
    else: Fvalue = (2*precision*recall)/(precision+recall)

    # convert set to str
    def set_to_str(set_obj): return "||".join(set_obj)

    # get dict
    return {"TP":TP, "FP":FP, "FN":FN, "Fvalue":Fvalue, "nevents":nevents, "precision":precision, "recall":recall, "TP_predictedIDs":TP_predictedIDs, "true_positives_knownIDs":set_to_str(true_positives_knownIDs), "false_negatives_knownIDs":set_to_str(false_negatives_knownIDs), "true_positives_predictedIDs":set_to_str(true_positives_predictedIDs), "false_positives_predictedIDs":set_to_str(false_positives_predictedIDs)}

def get_represenative_filtersDict_for_filtersDict_list(filtersDict_list, type_filters="most_conservative"):

    """Takes a lis, each position with a list of filters like passed to get_tupleBreakpoints_for_filters_GRIDSS and returns a representative dict, according to less_conservative"""

    # map a score for each dict
    score_to_dict = {sum([g_filterName_to_filterValue_to_Number[fname][find_nearest(g_filterName_to_filtersList[fname], fvalue)] for fname, fvalue in filtersDict.items()]) : filtersDict for filtersDict in filtersDict_list}

    # get the dict with the min or max score, depedning on the approach
    type_filters_to_getPositionFun = {"less_conservative":min, "most_conservative":max}

    return score_to_dict[type_filters_to_getPositionFun[type_filters](score_to_dict.keys())]

def get_changing_fields_in_df_benchmark(df):

    """This function takes a df such as the input of get_best_most_conservative_row_df_benchmark and returns a set with the keys that are different across rows"""

    changing_fields = set()

    # go through each field
    for f in df.columns:

        # filters dict are treated specially
        if f=="filters_dict": 
            if len(set(df[f].apply(get_dict_as_tuple)))>1: changing_fields.add(f)

        elif len(set(df[f]))>1: changing_fields.add(f)

    return changing_fields

def get_best_most_conservative_row_df_benchmark(df_benchmark):

    """Takes a df_benchmark, and returns the row with the most conservative row, given that it has the highest Fvalue. The least conservative is made in a step wise way filtering several things one after the other"""

    # get the maximum df
    df_best = df_benchmark[df_benchmark.Fvalue==max(df_benchmark.Fvalue)]
    if len(df_best)==1: return df_best.iloc[0]
    
    # get the df with the highest precision
    df_best = df_best[df_best.precision==max(df_best.precision)]
    if len(df_best)==1: return df_best.iloc[0]

    # get the one with the highest recall
    df_best = df_best[df_best.recall==max(df_best.recall)]
    if len(df_best)==1: return df_best.iloc[0]

    # get the most conservative set of filters for gridss
    most_conservative_filtersDict_tuple = get_dict_as_tuple(get_represenative_filtersDict_for_filtersDict_list(df_best.filters_dict, type_filters="most_conservative"))
    df_best = df_best[df_best.filters_dict.apply(get_dict_as_tuple)==most_conservative_filtersDict_tuple]
    if len(df_best)==1: return df_best.iloc[0]

    # get the minimum clove_max_rel_coverage_to_consider_del
    df_best = df_best[df_best.clove_max_rel_coverage_to_consider_del==min(df_best.clove_max_rel_coverage_to_consider_del)]
    if len(df_best)==1: return df_best.iloc[0]

    # get the max clove_min_rel_coverage_to_consider_dup
    df_best = df_best[df_best.clove_min_rel_coverage_to_consider_dup==max(df_best.clove_min_rel_coverage_to_consider_dup)]
    if len(df_best)==1: return df_best.iloc[0]

    # get filters with min tolerated cov
    if "gridss_maxcoverage" in df_best.keys():
        df_best = df_best[df_best.gridss_maxcoverage==min(df_best.gridss_maxcoverage)]
        if len(df_best)==1: return df_best.iloc[0]

    # if any, take the ones filtering regions in gridss
    if "gridss_regionsToIgnoreBed" in df_best.keys():
        
        if any(df_best.gridss_regionsToIgnoreBed!=""):    
            df_best = df_best[df_best.gridss_regionsToIgnoreBed!=""]
            if len(df_best)==1: return df_best.iloc[0]


    # at the end just return the best one
    print("Warning: there is no single type of filtering that can fullfill all the requirements") 


    # if you didn't find a single best, raise error
    print("\nthis is the best df:\n", df_best, "printing the non equal fields across all rows:\n")
    changing_fields = get_changing_fields_in_df_benchmark(df_best)
    for f in changing_fields:
        print("\t", f)
        for Irow in range(len(df_best)): print("\t\t", df_best[f].iloc[Irow])


    raise ValueError("There is not a single best filtering")


all_svs = {'translocations', 'insertions', 'deletions', 'inversions', 'tandemDuplications', 'remaining'}
def get_integrated_benchmarking_fields_series_for_setFilters_df(df):

    """This function takes a grouped per-filterSet df and returns a row with all the integrated accuracy measurements. The filters of gridss that are best for each SV may vary. If so we will take the most conservative filters of all of the fiters that are best for each SV."""

    # get a df where each row is one df
    df_best_filters = df.groupby("svtype").apply(get_best_most_conservative_row_df_benchmark)

    # debug when there are filters_dict
    if "filters_dict" in set(df_best_filters.keys()):

        if len(set(df_best_filters["filters_dict"].apply(get_dict_as_tuple)))!=1: 
            pass
            #raise ValueError("There are more than 1 filtersDict")

    # initialize a dict that will contain all the integrated filters
    integrated_benchmarking_results_dict = {}

    # get the numeric vals
    for f in ["FN", "FP", "TP", "nevents"]: integrated_benchmarking_results_dict[f] = sum(df_best_filters[f])

    # get through the event IDs 
    for f in ['TP_predictedIDs', 'false_negatives_knownIDs', 'false_positives_predictedIDs', 'true_positives_knownIDs', 'true_positives_predictedIDs']:  integrated_benchmarking_results_dict[f] = "||".join(df_best_filters[f].apply(str))

    # add the calculation of accuracy statistics
    TP = integrated_benchmarking_results_dict["TP"]
    FP = integrated_benchmarking_results_dict["FP"]
    FN = integrated_benchmarking_results_dict["FN"]
    nevents = integrated_benchmarking_results_dict["nevents"]

    if nevents==0: precision=1.0; recall=1.0
    else:
        if TP==0 and FP==0: precision =  0.0
        else: precision = TP/(TP + FP)
        recall = TP/(TP + FN)
        
    if precision<=0.0 or recall<=0.0: Fvalue = 0.0
    else: Fvalue = (2*precision*recall)/(precision+recall)

    integrated_benchmarking_results_dict["precision"] = precision
    integrated_benchmarking_results_dict["recall"] = recall
    integrated_benchmarking_results_dict["Fvalue"] = Fvalue

    # add other fields
    integrated_benchmarking_results_dict["svtype"] = "integrated"

    # add the fileds corresponding to when there are filters dicts
    if "filters_dict" in set(df_best_filters.keys()): 

        integrated_benchmarking_results_dict["filters_dict"] = get_represenative_filtersDict_for_filtersDict_list(list(df_best_filters["filters_dict"]), type_filters="most_conservative")
        integrated_benchmarking_results_dict["clove_max_rel_coverage_to_consider_del"] = df_best_filters.loc["deletions", "clove_max_rel_coverage_to_consider_del"]
        integrated_benchmarking_results_dict["clove_min_rel_coverage_to_consider_dup"] = df_best_filters.loc["tandemDuplications", "clove_min_rel_coverage_to_consider_dup"]
        integrated_benchmarking_results_dict["threshold_p_unbalTRA"] = df_best_filters.loc["translocations", "threshold_p_unbalTRA"]

        integrated_benchmarking_results_dict["median_insert_size"] = df_best_filters.loc["deletions", "median_insert_size"]
        integrated_benchmarking_results_dict["median_insert_size_sd"] = df_best_filters.loc["deletions", "median_insert_size_sd"]
        integrated_benchmarking_results_dict["sorted_bam"] = df_best_filters.loc["deletions", "sorted_bam"]
        integrated_benchmarking_results_dict["median_coverage"] = df_best_filters.loc["deletions", "median_coverage"]

    return pd.Series(integrated_benchmarking_results_dict)

def benchmark_processedSVs_against_knownSVs_inHouse(svtype_to_predsvfile, know_SV_dict, fileprefix, replace=False, add_integrated_benchmarking=True):

    """Takes two dictionaries that map some SVfiles. It runs, for all the types in svtype_to_predsvfile, a benchmarking against the known ones, writing a file under fileprefix. It returns a df of this benchmark, created with functions written here. It returns as matching events those that have an overlap of at least 50 bp.

    know_SV_dict marks the expected events. If they do not exist you have 0 accuracy.

    add_integrated_benchmarking indicates whether to perform a global benchmarking (not only per svtype).

    The 'analysis_benchmarking' feature was removed from this version """

    # map cmplex events to a boolean field
    complexEvt_to_boolField = {"translocations":"Balanced", "insertions":"Copied"}


    # initialize benchmark dict
    benchmark_dict = {}

    # go through each type of event
    for svtype in know_SV_dict:
        print("benchmarking %s"%svtype)

        # load dataframes
        if svtype in svtype_to_predsvfile.keys(): 
            df_predicted = pd.read_csv(svtype_to_predsvfile[svtype], sep="\t")
        else: df_predicted = pd.DataFrame()

        df_known = pd.read_csv(know_SV_dict[svtype], sep="\t")

        # define the fields to find overlaps
        equal_fields = svtype_to_fieldsDict[svtype]["equal_fields"]
        approximate_fields = svtype_to_fieldsDict[svtype]["approximate_fields"]
        chromField_to_posFields = chromField_to_posFields[svtype]["approximate_fields"]

        # get the dict of the benchmark
        dict_benchmark_svtype = get_SVbenchmark_dict(df_predicted, df_known, equal_fields=equal_fields, approximate_fields=approximate_fields, chromField_to_posFields=chromField_to_posFields)
        dict_benchmark_svtype["svtype"] = svtype

        # keep
        benchmark_dict[svtype] = dict_benchmark_svtype

    # get the benchmarking
    df_benchmark = pd.DataFrame(benchmark_dict).transpose()

    ###### perform integrated benchmarking #####
    if add_integrated_benchmarking is True:

        # get a per-filt row
        integrated_df_benchmark = pd.DataFrame({"integrated": get_integrated_benchmarking_fields_series_for_setFilters_df(df_benchmark)}).transpose()

        # keep
        df_benchmark = df_benchmark.append(integrated_df_benchmark[list(df_benchmark.keys())])
    
    ###########################################

    return df_benchmark

"""
- known_genomes_withSV_and_shortReads_table is a file with a table that has three fields: ID,assembly,shoort_reads1,short_reads2 . This can be, for example a set of NANOPORE assemblies of Candida glabrata and the corresponding short reads' sequencing in YPD
"""

def plot_bars_single_df_benchmark(df_benchmark, filename):

    """Takes a df_benchmark with different SV types and plots the results into benchmark"""

    palette_dict = {"precision":"black", "recall":"magenta", "Fvalue":"cyan", "FN":"magenta", "FP": "black", "TP":"cyan", "nevents":"gray"}

    # make a long df where you have a "value" field
    df_long = pd.melt(df_benchmark, id_vars=["svtype"], value_vars=["precision", "recall", "Fvalue", "FN", "FP", "TP", "nevents"])

    # make two subplots, one for each type of data
    fig = plt.figure(figsize=(5, 7))

    for I, variable_fields in enumerate([["precision", "recall", "Fvalue"], ["FP","FN","TP","nevents"]]):

        # initialize subpplot
        ax =  plt.subplot(2, 1, I+1)

        # get df and get barplot
        df = df_long[df_long.variable.isin(variable_fields)]
        sns.barplot(x="svtype", y="value", data=df, hue="variable", palette=palette_dict)

        # add a line
        if I==0: 
            for Y in [0.9, 0.95, 1.0]: plt.axhline(Y, color="gray", linewidth=0.9, linestyle="--")

        # change rotation
        for label in ax.get_xticklabels(): label.set_rotation(45)

    #plt.show()
    fig.tight_layout()  # otherwise the right y-label is slightly 
    fig.savefig(filename, bbox_inches='tight');
    plt.close(fig)

def test_SVgeneration_from_assembly(reference_genome, outdir, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", nvars=100):

    """This function reports how well the finding of SV from a genome works from random simulations. Writing under outdir"""

    # define the output
    precision_and_recall_filename = "%s/precision_and_recall_SVgeneration_from_assembly.pdf"%outdir
    if file_is_empty(precision_and_recall_filename) or replace is True:

        # initialize the start time
        pipeline_start_time = time.time()

        # prepare files
        make_folder(outdir)

        print("WORKING ON THE VALIDATION THAT WE CAN FIND READS IN AN ASSEMBLY")

        # initialize a df that will contain the benchmarking
        all_df_benchmark_longReads = pd.DataFrame()

        # go through each simulation
        for simID in range(n_simulated_genomes):
            print("working on simulation %i"%simID)

            # define outdir 
            outdir_sim = "%s/simulation_%i"%(outdir, simID); make_folder(outdir_sim)

            # generate genome with simulated SVs
            sim_svtype_to_svfile, rearranged_genome = rearrange_genomes_simulateSV(reference_genome, outdir_sim, replace=replace, nvars=nvars, mitochondrial_chromosome=mitochondrial_chromosome)

            # get the variants from simulating reads. Always ploidy 1 to get homozygous SVs
            predicted_svtype_to_svfile, df_gridss = generate_tables_of_SV_between_genomes_gridssClove(rearranged_genome, reference_genome, replace=replace, threads=threads, mitochondrial_chromosome=mitochondrial_chromosome)

            # get a df of benchmarking
            fileprefix = "%s/rearranged_genome_benchmarking_SV"%outdir_sim
            df_benchmark_longReads = benchmark_processedSVs_against_knownSVs_inHouse(predicted_svtype_to_svfile, sim_svtype_to_svfile, fileprefix, replace=replace)

            # keep
            df_benchmark_longReads["simID"] = [simID]*len(df_benchmark_longReads)
            all_df_benchmark_longReads = all_df_benchmark_longReads.append(df_benchmark_longReads)

        # plot the benchmarking
        plot_bars_single_df_benchmark(all_df_benchmark_longReads, precision_and_recall_filename)

        # at the end clean the generation
        clean_reference_genome_windows_files(reference_genome)

        print("--- the testing of SV generation from an assembly took %s seconds in %i cores ---"%(time.time() - pipeline_start_time, threads))


def get_speciesTree_multipleGenomes_JolyTree(input_dir_withGenomes, outdir, threads=4, replace=False):

    """This function generates a species tree under outdir with all the genomes (files ending with fasta) in input_dir_withGenomes. It returns the newick file with the tree"""

    # make the outdir
    make_folder(outdir)


    # define the outprefix and the expected species tree file
    outprefix = "%s/outputJolyTree"%outdir
    species_treefile = "%s.nwk"%outprefix

    if file_is_empty(species_treefile) or replace is True:

        # move all the fasta files in input_dir_withGenomes into input_dir
        input_dir = "%s/input_genomes"%outdir; make_folder(input_dir)
        for file in os.listdir(input_dir_withGenomes):
            origin_file = "%s/%s"%(input_dir_withGenomes, file)

            # if it is a fasta file, softlink to input_dir
            if file.split(".")[-1] in {"fasta", "fa"} and os.path.isfile(origin_file):
                dest_file = "%s/%s"%(input_dir, file)
                if file_is_empty(dest_file): run_cmd("ln -s %s %s"%(origin_file, dest_file))


        # run JolyTree
        print("running JolyTree to get species tree")
        run_cmd("%s -i %s -b %s -t %i"%(JolyTree_sh, input_dir, outprefix, threads))

    return species_treefile

def ask_if_overlapping_breakends_in_parents_withEqualChromosomes(r, parents_gridss_df, tol_bp):

    """Returns a boolean that indicates whether ther is any breakend in parents that overlaps with the breakedn in r by less than tol_bp bp, where the orientation of the brekends is expected to be the same"""

    # get the df in parents_gridss_df that fullfils the overlap conditions
    df_parents_overlap = parents_gridss_df[(parents_gridss_df.other_orientation==r["other_orientation"]) 
                                         & ((parents_gridss_df.POS-r["POS"]).apply(abs)<=tol_bp) 
                                         & ((parents_gridss_df.other_position-r["other_position"]).apply(abs)<=tol_bp) ]

    # return True if there is some overlap
    if len(df_parents_overlap)>0: return True
    else: return False

def get_eventIDs_already_in_parents(sampleID, parentIDs, sampleID_to_dfGRIDSS, tol_bp=50):

    """This function returns a set with the eventIDs that are in sample but also in any of the parents. The idea is that an event is considered to be in a parent if any of the breakends overlaps with a breakend in the parent, where the orientation is the same and the position is less or equal far appart than tol_bp """

    # define the gridss df of the sample
    sample_gridss_df = sampleID_to_dfGRIDSS[sampleID]
    if len(parentIDs)==0: return set()

    # stack all the parents' gridss 
    parents_gridss_df = pd.concat([sampleID_to_dfGRIDSS[p] for p in parentIDs])[["#CHROM", "POS", "other_chromosome", "other_position", "other_orientation"]]

    # add a boolean that indicates whether the breakend can be found in the parents, for set of pairs of chromosomes
    all_chromosomal_combinations = sample_gridss_df[["#CHROM", "other_chromosome"]].sort_values(by=["#CHROM", "other_chromosome"]).drop_duplicates().values
    final_sample_gridss_df = pd.DataFrame() # this will contain the boolean

    # go through all chromosomal combinations appearing in the sample
    print("getting overlaps in breakpoints between samples")
    for I, (chrom, other_chromosome) in enumerate(all_chromosomal_combinations):

        # get the dfs with these chromosomal combinations
        chr_sample_gridss_df = sample_gridss_df[(sample_gridss_df["#CHROM"]==chrom) & (sample_gridss_df["other_chromosome"]==other_chromosome)]
        chr_parents_gridss_df = parents_gridss_df[(parents_gridss_df["#CHROM"]==chrom) & (parents_gridss_df["other_chromosome"]==other_chromosome)][["POS", "other_position", "other_orientation"]]

        # add the boolean that indicates whether the breakend is in the parents
        chr_sample_gridss_df["is_in_parents"] = chr_sample_gridss_df[["POS", "other_position", "other_orientation"]].apply(lambda r: ask_if_overlapping_breakends_in_parents_withEqualChromosomes(r, chr_parents_gridss_df, tol_bp), axis=1)

        # keep 
        final_sample_gridss_df = final_sample_gridss_df.append(chr_sample_gridss_df)

    # define the eventIDs
    final_sample_gridss_df["eventID"] = final_sample_gridss_df.INFO_EVENT.apply(lambda x: x+"o")

    # define the events that are in the parents
    all_eventIDs = set(final_sample_gridss_df.eventID)
    eventIDs_already_in_parents = set(final_sample_gridss_df[final_sample_gridss_df.is_in_parents].eventID)
    print("There are %i of %i breakpoints already in the parents"%(len(eventIDs_already_in_parents), len(all_eventIDs)))

    return eventIDs_already_in_parents

def format_svDF(svDF, svtype="unknown", interesting_chromosomes="all", sampleName="sampleX"):

    """Takes an svDF and formats the ID if suspicious. estimate_fn_af is the function that states which breakend af has to be taken"""


    # return the same if empty
    if len(svDF)==0: return svDF

    if set(svDF.ID)=={""}: svDF["ID"] = svDF.Name

    # add a uniqueID
    svDF["uniqueID"] = ["%s_%s_%i"%(sampleName, svtype, I+1) for I in range(len(svDF))]

    # filter and get only chromosomes that are in interesting_chromosomes
    if interesting_chromosomes!="all":

        # define the chromosomal fields
        chromosome_fields = {"Chr", "ChrA", "ChrB", "#CHROM", "CHR2"}.intersection(set(svDF.keys()))

        # get the df where the interesting_chromosomes are involved
        svDF = svDF[svDF.apply(lambda r: all([r[c] in interesting_chromosomes for c in chromosome_fields]), axis=1)]

    def get_estimate_AF_for_breakends(bends_metadata_dict, AF_field="real_AF", estimate_fn=min):

        """Takes a dict that maps each breakpoints to a list of metadata of each breakend, and returns the estimate_fn AF observed in any breakpoint"""

        estimate_AF =  estimate_fn([estimate_fn([bend_info[AF_field] for bend_info in list_breakend_info]) for list_breakend_info in bends_metadata_dict.values()])

        return estimate_AF

    # get the estimated minimum and maxium allele freq
    for estimate_fn_name, estimate_fn in [("min", min), ("max", max), ("mean", np.mean)]:
        svDF["estimate_AF_%s"%estimate_fn_name] = svDF.bends_metadata_dict.apply(get_estimate_AF_for_breakends, estimate_fn=estimate_fn)

    return svDF

def get_sampleID_to_svtype_to_svDF_filtered(sampleID_to_svtype_to_file, sampleID_to_dfGRIDSS, sampleID_to_parentIDs={}, breakend_info_to_keep=['#CHROM', 'POS', 'other_coordinates', 'allele_frequency', 'allele_frequency_SmallEvent', 'real_AF', 'FILTER', 'inserted_sequence', 'has_poly16GC', 'length_inexactHomology', 'length_microHomology']):

    """This function takes a dictionary that maps sampleIDs to svtpes and the corresponding files (the ones returned in the SV-calling pipeline) and the corresponding gridss dataframes, returning a sampleID_to_svtype_to_svDF, after removal of the variants that are in sampleID_to_parentIDs. The idea is that if any breakend is shared between the sample and any of the parents it is removed and so does any sv that relies on this breakend. Other arguments:

    All across this function, eventID reflects the INFO_EVENT with an extra fields

    """

    # debug empty 
    if len(sampleID_to_parentIDs)==0: sampleID_to_parentIDs = {s:set() for s in sampleID_to_svtype_to_file}

    # map each sampleID to the eventIDs (with a final o, as formated in the sampleID_to_svtype_to_file) that are already in the parents
    all_samples = sorted(sampleID_to_parentIDs.keys())
    get_eventIDs_already_in_parents_inputs = [(s, sampleID_to_parentIDs[s], sampleID_to_dfGRIDSS) for s in all_samples]

    # get the overlapping events with a map
    map_eventIDs_already_in_parents = list(map(lambda x: get_eventIDs_already_in_parents(x[0], x[1], x[2]), get_eventIDs_already_in_parents_inputs))    

    sampleID_to_eventIDs_alreadyInParents = dict(zip(all_samples, map_eventIDs_already_in_parents))


    # load the dfs of the sampleID_to_svtype_to_file and add the info from sampleID_to_dfGRIDSS, keeping the info of the breakpoints that are not in the parent
    sampleID_to_svtype_to_svDF = {}
    for sample, svtype_to_file in sampleID_to_svtype_to_file.items():
        print(sample)

        # get the gridss df
        df_gridss = sampleID_to_dfGRIDSS[sample]

        # add the eventID, as formated in the SVdict
        df_gridss["eventID"] = df_gridss.INFO_EVENT.apply(lambda x: x+"o")

        # define the events that are already in the parents
        all_breakpoints = set(df_gridss["eventID"]).union({""})
        breakpoints_already_in_parents = sampleID_to_eventIDs_alreadyInParents[sample]
        print("There are %i of %i breakpoints already in the parents"%(len(breakpoints_already_in_parents), len(all_breakpoints)))

        # go through each svtype
        for svtype, file in svtype_to_file.items():

            if file!="":

                # get as df
                df_sv = pd.read_csv(file, sep="\t")

                # get for empty df
                if "ID" not in df_sv.keys(): df_sv["ID"] = [""]*len(df_sv)

                # get only the svs that are not in the parents
                df_sv["IDs_set"] = df_sv.ID.apply(lambda x: set(re.split("\+|\-", x)))

                if len(df_sv)>0:

                    # check that all of them are in the breakpoints called by gridss
                    all_breakpoints_in_sv = set.union(*df_sv.IDs_set)
                    breakpoints_not_in_df_gridss = all_breakpoints_in_sv.difference(all_breakpoints)
                    if len(breakpoints_not_in_df_gridss)>0: 
                        print("These are breakpoints not in the df_gridss: %s"%(breakpoints_not_in_df_gridss))
                        raise ValueError("There are breakpoints not in the df_gridss, suggesting some errors, such as the fact that you are loading ")

                    # keep only the df with IDs that are not in the parents
                    df_sv = df_sv[df_sv.IDs_set.apply(lambda ids: len(ids.intersection(breakpoints_already_in_parents))==0)]
                    if len(df_sv)==0: 
                        print("There are no new SVs in %s, %s"%(sample, svtype))
                        continue

                    # add the metadata of the breakends in breakend_info_to_keep from df_gridss
                    if len(df_gridss)>0:
                        df_gridss_sv = df_gridss[df_gridss.eventID.isin(set.union(*df_sv.IDs_set))].set_index("eventID", drop=False)
                        df_gridss_sv["metadata_dict"] = df_gridss_sv.apply(lambda r: {col: r[col] for col in breakend_info_to_keep}, axis=1)

                        df_sv["bends_metadata_dict"] = df_sv.IDs_set.apply(lambda ids: {I:[metadata_dict for metadata_dict in df_gridss_sv.loc[I, "metadata_dict"]] for I in ids})

                    else: df_sv["bends_metadata_dict"] = [{}]*len(df_sv)

            # empty, there are no such events
            else: df_sv = pd.DataFrame()

            # keep
            sampleID_to_svtype_to_svDF.setdefault(sample ,{}).setdefault(svtype, df_sv)

    # add the metadata that is necessary and also unique IDs and allele freqs
    sampleID_to_svtype_to_svDF = {sID : {sv : format_svDF(df, svtype=sv, sampleName=sID) for sv, df in SV_dict.items()} for sID, SV_dict in sampleID_to_svtype_to_svDF.items()}

    return sampleID_to_svtype_to_svDF

def add_svID_to_IDtoSVTYPEtoDF(ID_to_svtype_to_svDF):

    """Takes an ID_to_svtype_to_svDF and adds the svID if it overlaps with other vars"""

    # define all svtypes
    all_svtypes = all_svtypes = set.union(*[set(x.keys()) for x in ID_to_svtype_to_svDF.values()])

    # go through each svtype
    for svtype in all_svtypes:
        print("getting svID for %s"%svtype)

        # define the fields to filter on
        equal_fields = svtype_to_fieldsDict[svtype]["equal_fields"]
        approximate_fields = svtype_to_fieldsDict[svtype]["approximate_fields"]
        chromField_to_posFields = svtype_to_fieldsDict[svtype]["chromField_to_posFields"]

        # initialize a dict that maps each svID to a row with its' coordinates
        svID_to_svIDrow = {0:{}}

        # go through each sample looking for 
        for ID in ID_to_svtype_to_svDF:

            # if the svtype exists
            if svtype in ID_to_svtype_to_svDF[ID] and len(ID_to_svtype_to_svDF[ID][svtype])>0:

                # get the df
                df = ID_to_svtype_to_svDF[ID][svtype].set_index("uniqueID")

                # initialize a dict that maps each uniqueID to a svID
                uniqueID_to_svID = {}

                # go through each SV of this df
                for uniqueID, query_row in df.iterrows():

                    # initialize the found svID
                    svID_found = False

                    # go through each row of the svIDs looking if any is matching
                    for previous_svID, svID_row in svID_to_svIDrow.items():
                        if len(svID_row)==0: continue

                        # if there is an overlap with the svID_row, keep the previous_svID 
                        if get_is_matching_predicted_and_known_rows(query_row, svID_row, equal_fields, approximate_fields, chromField_to_posFields): 

                            uniqueID_to_svID[uniqueID] = previous_svID
                            svID_found = True
                            break

                    # if there was no overlapping sv, initialize it
                    if svID_found is False:

                        svID = max(svID_to_svIDrow)+1
                        uniqueID_to_svID[uniqueID] = svID
                        svID_to_svIDrow[svID] = query_row

                # add the svID to the df
                ID_to_svtype_to_svDF[ID][svtype]["svID"] = ID_to_svtype_to_svDF[ID][svtype].uniqueID.apply(lambda x: "%s_%i"%(svtype, uniqueID_to_svID[x]))

def get_svIDs_inMoreThan1species(ID_to_svtype_to_svDF):

    """This function takes an ID_to_svtype_to_svDF were there is an svID and adds a boolean to each df indicating whether the svID is in more than 1 species. """

    # get all svtypes
    all_svtypes = set.union(*[set(x.keys()) for x in ID_to_svtype_to_svDF.values()])

    # initialize set
    svIDs_inMoreThan1species = set()

    # iterate through SVs
    for svtype in all_svtypes:
        print("assessing monophily for %s"%svtype)

        # get the IDs 
        ID_to_svIDs = {ID : set(ID_to_svtype_to_svDF[ID][svtype].svID)  for ID in ID_to_svtype_to_svDF if svtype in ID_to_svtype_to_svDF[ID] and len(ID_to_svtype_to_svDF[ID][svtype])>0}

        # fill empty IDs
        for ID in set(ID_to_svtype_to_svDF).difference(set(ID_to_svIDs)): ID_to_svIDs[ID] = set()

        # get all the svIDs
        all_svIDs = union_empty_sets(ID_to_svIDs.values())

        # go through each ID
        for svID in all_svIDs:

            # calculate the number of IDs with this svID. If it is above 1, get
            nIDs_with_svID = len({ID for ID, svIDs in ID_to_svIDs.items() if svID in svIDs})
            if nIDs_with_svID>=2: svIDs_inMoreThan1species.add(svID)


    print("There are %i vars in more than 1 species"%len(svIDs_inMoreThan1species))

    return svIDs_inMoreThan1species

    
def prune_IDtoSVTYPEtoDF_keeping_HighConfidenceVars(ID_to_svtype_to_svDF, species_tree, min_af=0.4):

    """This function takes a df that maps and ID to an svtypes to df and only keeps those that are High Confidence. These are vars that are either:

    - found in >1 IDs
    - The minimum allele frequency of all breakends is above min_af and the filter is PASS"""

    # first add the svID to each var. This is an ID that represents this variant across all sampleIDs, so that if there is a var overlapping across 2 samples it will have the same svID

    add_svID_to_IDtoSVTYPEtoDF(ID_to_svtype_to_svDF)

    # get the svIDs that follow the species phylogeny
    svIDs_inMoreThan1species = get_svIDs_inMoreThan1species(ID_to_svtype_to_svDF)

    # prune the df to keep only high-confidence vars
    for ID, svtype_to_svDF in ID_to_svtype_to_svDF.items():
        for svtype, svDF in svtype_to_svDF.items():

            if len(svDF)==0:  ID_to_svtype_to_svDF[ID][svtype] = svDF
            else:   

                # add whether the var is in more than 1 species
                svDF["sv_in_more_than1_spp"] = svDF.svID.isin(svIDs_inMoreThan1species)

                # add whether it is high confidence
                svDF["all_bends_PASS"] = svDF.apply(lambda r: set.union(*[set([bend_info["FILTER"] for bend_info in list_breakend_info]) for list_breakend_info in r["bends_metadata_dict"].values()])=={"PASS"}, axis=1)

                svDF["PASS_and_highAF"] = (svDF.estimate_AF_min>=min_af) & (svDF.all_bends_PASS)

                svDF["high_confidence"] = (svDF.PASS_and_highAF) | (svDF.sv_in_more_than1_spp)

                # keep the df that has high confidence
                ID_to_svtype_to_svDF[ID][svtype] = svDF[svDF.high_confidence]

def get_is_overlapping_query_vs_target_region(q, r):

    """This function takes two 'bed'-like regions and returns whether they are overlapping by some extent """

    return (q["chromosome"]==r["chromosome"]) and ((r["start"]<=q["start"]<=r["end"]) or (r["start"]<=q["end"]<=r["end"]) or (q["start"]<=r["start"]<=q["end"]) or (q["start"]<=r["end"]<=q["end"]))

def get_ID_to_svtype_to_svDF_for_setOfGenomes_highConfidence(genomes_withSV_and_shortReads_table, reference_genome, outdir, species_treefile=None, replace=False, threads=4):

    """Generates a dict that maps each sample in genomes_withSV_and_shortReads_table to an svtype and a DF with all the info about several vars. It only gets the high-confidence vars"""

    # load the df
    df_genomes = pd.read_csv(genomes_withSV_and_shortReads_table, sep="\t")

    # define an outdir that will store all the real_vars
    make_folder(outdir)
    all_realVars_dir = "%s/all_realVars"%(outdir)
    if replace is True: delete_folder(all_realVars_dir)
    make_folder(all_realVars_dir)

    # define the final object that contains them all
    ID_to_svtype_to_svDF_file = "%s/ID_to_svtype_to_svDF.py"%all_realVars_dir

    if file_is_empty(ID_to_svtype_to_svDF_file) or replace is True:

        # initialize dicts that keep them all together
        all_sampleID_to_svtype_to_file = {}
        all_sampleID_to_dfGRIDSS = {}
       
        # generate all real vars
        for ID, assembly in df_genomes[["ID", "assembly"]].values:
            print("getting vars for %s"%ID)

            # softlink the genome
            dest_genomeFile = "%s/genome_%s.fasta"%(all_realVars_dir, ID)
            if file_is_empty(dest_genomeFile): run_cmd("ln -s %s %s"%(assembly, dest_genomeFile))

            # find the real vars
            svtype_to_svfile, df_gridss = generate_tables_of_SV_between_genomes_gridssClove(dest_genomeFile, reference_genome, replace=replace, threads=threads)

            # keep 
            all_sampleID_to_svtype_to_file[ID] =  svtype_to_svfile
            all_sampleID_to_dfGRIDSS[ID] = df_gridss

        # get a df that has all the info for each SV, and then the df with allele freq, metadata and 
        ID_to_svtype_to_svDF = get_sampleID_to_svtype_to_svDF_filtered(all_sampleID_to_svtype_to_file, all_sampleID_to_dfGRIDSS)

        # build the species tree
        if species_treefile is None: 
            outdir_species_tree = "%s/output_speciesTree_JolyTree"%outdir
            species_treefile = get_speciesTree_multipleGenomes_JolyTree(all_realVars_dir, outdir_species_tree, threads=threads, replace=replace)

        # get the names so that they are IDs
        species_tree = Tree(species_treefile)
        for l in species_tree.get_leaves(): l.name = l.name.split("_")[1]
        print("This is the considered species tree:\n", species_tree)

        # preseve only the high confidence vars from ID_to_svtype_to_svDF
        prune_IDtoSVTYPEtoDF_keeping_HighConfidenceVars(ID_to_svtype_to_svDF, species_tree)

        print(ID_to_svtype_to_svDF)

        # save
        save_object(ID_to_svtype_to_svDF, ID_to_svtype_to_svDF_file)

    else: ID_to_svtype_to_svDF = load_object(ID_to_svtype_to_svDF_file)

    return ID_to_svtype_to_svDF

def add1_unless_it_is_minus1(x):

    """Takes an int and adds 1 unless it is -1"""

    if x==-1: return x
    else: return x+1 

def set_position_to_max(pos, maxPos):

    """Sets a position to a maximum"""

    if pos>maxPos: return maxPos
    else: return pos

def get_compatible_real_svtype_to_file(genomes_withSV_and_shortReads_table, reference_genome, outdir, species_treefile=None, replace=False, threads=4, max_nvars=100):

    """This function generates a dict of svtype to the file for SVs that are compatible and ready to insert into the reference_genome. All the files are written into outdir. Only a set of 'high-confidence' SVs are reported, which are those that have ether a minimum allele frequency above 0.75 and all breakpoints with 'PASS' or are found in >2 genomes and following the phylogeny of the species tree. This is inferred with JolyTree if not provided. At the end, this pipeline reports a set of compatible SVs, that are ready to insert into RSVsim (so that the coordinates are 1-based). 

    Rememeber that this will not include 'remaining' events, as these can't be inserted 
    
    For transloctaions, there will only be one translocation allowed for each arm. This is how RSVSim allows the simulation.

    max_nvars is the maximum number of variants of each type that will be generated
    """

    # initialize the start time
    pipeline_start_time = time.time()

    # get all the high-confidence real variants
    ID_to_svtype_to_svDF = get_ID_to_svtype_to_svDF_for_setOfGenomes_highConfidence(genomes_withSV_and_shortReads_table, reference_genome, outdir, species_treefile=species_treefile, replace=replace, threads=threads)

    # define the df with the realVars info
    all_realVars_dir = "%s/all_realVars"%(outdir)

    ########## GET THE FILES OF COMPATIBLE SVs ##########

    # create a folder where to write the high-confidence vars
    highConfidenceVars_perGenome_dir = "%s/highConfidenceVars_perGenome"%all_realVars_dir
    make_folder(highConfidenceVars_perGenome_dir)

    # create a folder where you have the set of compatible vars
    SVs_compatible_to_insert_dir = "%s/SVs_compatible_to_insert"%outdir; make_folder(SVs_compatible_to_insert_dir)

    # initialize the dict that will contain the final files
    compatible_real_svtype_to_file = {}

    # initialize a bed_df with all the regions
    df_bed_allRegions = pd.DataFrame()

    # define all chromosomes as interesting
    all_chromosomes = {seq.id for seq in SeqIO.parse(reference_genome, "fasta")}

    # go through each svtype
    all_svtypes = set.union(*[set(x.keys()) for x in ID_to_svtype_to_svDF.values()]).difference({"remaining"})
    all_svtypes = [s for s in ["insertions", "translocations", "inversions", "tandemDuplications", "deletions"] if s in all_svtypes]

    # map each chromosome to the len
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

    for svtype in all_svtypes:
        outfile_compatible_SVs = "%s/%s.tab"%(SVs_compatible_to_insert_dir, svtype)

        if file_is_empty(outfile_compatible_SVs) or replace is True:
        #if True:

            print("writing consensus %s"%svtype)

            # initalize a df with all the compatible svDFs
            compatible_svDF = pd.DataFrame()

            # go through each ID
            for ID in ID_to_svtype_to_svDF.keys():

                # debug empty dfs
                if svtype not in ID_to_svtype_to_svDF[ID] or len(ID_to_svtype_to_svDF[ID][svtype])==0: continue

                svDF = ID_to_svtype_to_svDF[ID][svtype].set_index("uniqueID", drop=False)

                # write the svDF to the high-confidence dir
                sampleDir = "%s/%s"%(highConfidenceVars_perGenome_dir, ID); make_folder(sampleDir)
                svDF.to_csv("%s/%s.tab"%(sampleDir, svtype), sep="\t", header=True, index=False)

                # go though each variant
                for varID, sv_series in svDF.iterrows():

                    # define series as df
                    sv_series_df = pd.DataFrame({0 : sv_series}).transpose()

                    # get the bed of this var (for translocations only consider the positions)
                    sv_bed, nSVs = get_affected_region_bed_for_SVdf(sv_series_df, svtype, all_chromosomes, translocations_type="breakpoint_pos", chr_to_len=chr_to_len)

                    # get if there is any overlap with df_bed_allRegions
                    regions_overlapping = df_bed_allRegions.apply(lambda rprevious: any(sv_bed.apply(lambda rnew: get_is_overlapping_query_vs_target_region(rprevious, rnew), axis=1)), axis=1)

                    # if there is any region matching with the previous, continue, if not, keep
                    if any(regions_overlapping): continue
                    else:

                        # add the bed to the regions matching
                        df_bed_allRegions = df_bed_allRegions.append(sv_bed)

                        # add to the compatible SVs
                        compatible_svDF = compatible_svDF.append(sv_series_df)

            # get only the important fields
            compatible_svDF = compatible_svDF[svtype_to_fieldsDict[svtype]["all_fields"]]
            compatible_svDF.index = list(range(len(compatible_svDF)))

            # define the maping between the position field and the chromosome field
            posF_to_chrF = svtype_to_fieldsDict[svtype]["positionField_to_chromosome"]



            # add +1 to all the positions (this is because RSVSim requires 1-based positions), also that the last position of the chromosome is not exceeded
            for f in svtype_to_fieldsDict[svtype]["position_fields"]: compatible_svDF[f] = compatible_svDF.apply(lambda r: set_position_to_max(add1_unless_it_is_minus1(r[f]), chr_to_len[r[posF_to_chrF[f]]]), axis=1)

            # for translocations, define the set of maximum nvars based on the number of chromosomes
            if svtype=="translocations": real_max_nvars = min([max_nvars, len(all_chromosomes)-1])
            else: real_max_nvars = max_nvars

            # if this number exceeds the number of variants it will chop the df
            if len(compatible_svDF)>real_max_nvars: compatible_svDF = compatible_svDF.iloc[0:real_max_nvars]

            # write the compatible svDF into the final set of vars
            compatible_svDF.to_csv(outfile_compatible_SVs, sep="\t", header=True, index=False)

        # keep 
        compatible_real_svtype_to_file[svtype] = outfile_compatible_SVs

    #####################################################

    print("--- the generation of real SVs took %s seconds in %i cores ---"%(time.time() - pipeline_start_time, threads))

    return compatible_real_svtype_to_file

def get_insert_size_distribution(sorted_bam, replace=False, threads=4):

    """Takes a bam file of aligned paired end reads and retuns the mean and SD insert size of the library."""

    # define outfiles
    hist_file = "%s.histogram_insertsizes.pdf"%sorted_bam
    outfile = "%s.CollectInsertSizeMetrics.out"%sorted_bam; outfile_tmp = "%s.tmp"%outfile

    # run
    if file_is_empty(outfile) or replace is True: 
        remove_file(outfile_tmp)
        run_cmd("%s CollectInsertSizeMetrics HISTOGRAM_FILE=%s INPUT=%s OUTPUT=%s"%(picard_exec, hist_file, sorted_bam, outfile_tmp))
        os.rename(outfile_tmp, outfile)

    # get stats
    wrong_foot_lines = [l for l in open(outfile, "r").readlines() if len(l.split("\t"))==2 and not l.startswith("## METRICS CLASS")]
    df = pd.read_csv(outfile, sep="\t", skip_blank_lines=False, header=6, skipfooter=len(wrong_foot_lines)+2, engine='python').iloc[0]

    return (df["MEDIAN_INSERT_SIZE"], df["MEDIAN_ABSOLUTE_DEVIATION"])


def get_windows_infoDF_with_predictedFromFeatures_coverage(genome, distToTel_chrom_GC_to_coverage_fn, expected_coverage_per_bp, replace=False, window_l=1000, threads=4, make_plots=True):

    """This function gets a genome and returns a df for windows of the genome and the relative coverage predicted from distToTel_chrom_GC_to_coverage_fn"""

    windows_infoDF_file = "%s_windows_with_predictedRelCov_from_features.py"%genome

    if file_is_empty(windows_infoDF_file) or replace is True:
        print("getting relCov predicted from feats")

        ##### get the windows df ####

        # get the file
        windows_file = "%s.windows%ibp.bed"%(genome, window_l)
        run_cmd("%s makewindows -g %s.fai -w %i > %s"%(bedtools, genome, window_l, windows_file))

        # get into df 
        df = pd.read_csv(windows_file, sep="\t", header=None, names=["chromosome", "start", "end"])

        #############################

        # get the genome graph with no bps
        genomeGraph_outfileprefix = "%s_genomeGraph_withoutBPs.py"%genome
        genome_graph, df_positions_graph = get_genomeGraph_object(genome, None, None, genomeGraph_outfileprefix, replace=replace)

        # get the distance to the telomere in the df
        df["middle_position"] = (df.start + (df.end - df.start)/2).apply(int)
        df["distance_to_telomere"] = get_distance_to_telomere_series(df[["chromosome", "middle_position"]], genome_graph, df_positions_graph)

        # add the gc content
        gcontent_outfile = "%s_GCcontent.py"%windows_file
        df = get_df_with_GCcontent(df, genome, gcontent_outfile, replace=replace)

        # predict genome from the sequence features 
        df["cov_predicted_from_features"] = df.apply(lambda r: distToTel_chrom_GC_to_coverage_fn(r["distance_to_telomere"], r["chromosome"], r["GCcontent"]), axis=1)

        # get to non-0 values
        def get_to_non0_cov(cov, tsh=0.1):

            if cov<tsh: return tsh
            else: return cov

        # add the relative coverage
        df["predicted_relative_coverage"] = (df["cov_predicted_from_features"] / expected_coverage_per_bp).apply(get_to_non0_cov)

        # plot per chrom
        if make_plots is True:

            all_chromosomes = df.chromosome.unique()
            print("making plots")

            fig = plt.figure(figsize=(10, len(all_chromosomes)*4.5))
            for I, chrom in enumerate(all_chromosomes):

                # get df
                df_c = df[df.chromosome==chrom]

                # print each of the predictions
                ax = plt.subplot(len(all_chromosomes), 1, I+1)
                sns.lineplot(x="start", y="predicted_relative_coverage", data=df_c, linewidth=2, color="blue")

                ax.set_title(chrom)

            fig.tight_layout()  # otherwise the right y-label is slightly 
            filename="%s_predicted_relative_coverage.pdf"%(windows_file)
            fig.savefig(filename, bbox_inches='tight');
            plt.close(fig)

        # save
        save_object(df, windows_infoDF_file)

    else: df = load_object(windows_infoDF_file)

    return df

def get_read_length(bamfile, threads=4, nreads=5000, replace=False):

    """Calculates the readlength for a bamfile"""

    readlen_dist_file = "%s.read_length_dist_first%ireads.txt"%(bamfile, nreads); readlen_dist_file_tmp = "%s.tmp"%readlen_dist_file
    if file_is_empty(readlen_dist_file) or replace is True:

        print("The following command will throw a warning stating that 'samtools view: writing to standard output failed: Broken pipe'. This is because the output of samtools view is piped, which is expected.")
        cmd = "%s view --threads %i %s | head -n %i | cut -f10 | perl -ne 'chomp;print length($_) . \"\n\"' | sort > %s"%(samtools, threads, bamfile, nreads, readlen_dist_file_tmp)
        run_cmd(cmd)

        os.rename(readlen_dist_file_tmp, readlen_dist_file)

    return int(np.median([int(l.strip()) for l in open(readlen_dist_file, "r").readlines()]))


def count_number_read_pairs(bamfile, replace=False, threads=4):

    """counts the total number of reads of a bamfile"""

    read_count_file = "%s.flagstat"%bamfile; read_count_file_tmp = "%s.tmp"%read_count_file

    # get the total n reads
    if file_is_empty(read_count_file) or replace is True:

        print("calculating n reads")
        run_cmd("%s flagstat --threads %i %s > %s"%(samtools, threads, bamfile, read_count_file_tmp))
        os.rename(read_count_file_tmp, read_count_file)

    return [int(l.split()[0]) for l in open(read_count_file, "r").readlines() if " read1" in l][0]



def run_wgsim_pairedEnd_for_window(genome, chromosome, start, end, readPairs, read_length, outdir,  median_insert_size, median_insert_size_sd, replace, error_rate=0.02):

    """Takes a region of a genome and generates a fasta file for it under outdir. Then it runs wgsim to simulate total_readPairs. It returns a tupple of the (read1, read2) in .gz . It assumes that the error rate is 0.02, which was benchmarked in https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-184.

    This is run on haplotype mode (which means that all introduced vars are homozygous (although I am not introducing any vars))"""

    # get the time

    # define the ID
    ID = "%s_%i_%i_readlen%i_nreads%i_insertsize%i+-%i_error%.3f"%(chromosome, start, end, read_length, readPairs, median_insert_size, median_insert_size_sd, error_rate)

    # define region
    len_region = end - start

    # if the region is smaller than the read length, just make it a little bit longer than the readlen + insert size + sd
    min_length_region = int((read_length + median_insert_size + median_insert_size_sd)*1.5)
    if len_region<min_length_region: 
        print("WARNING: The length of the region is %i and the length of the reads is %i. Setting the length of the region to %i. We will try to set the region to the 5' of the chromosome, and to the 3' if not possible."%(len_region, read_length, min_length_region))

        # get the length of the chromosome
        len_chr = [len(chrom) for chrom in SeqIO.parse(genome, "fasta") if chrom.id==chromosome][0]

        # define the region to extend
        length_to_extend = min_length_region - len_region

        # if the start can be extended
        if (start - length_to_extend)>=0: start = start - length_to_extend

        # if the end can be extended
        elif (end + length_to_extend)<=len_chr: end = end + length_to_extend

        else: raise ValueError("Region %i cannot be corrected according to read length, maybe because the chromosome is to short. Check that %s is ok"%(ID, chromosome))

    # first get the region
    region_fasta = "%s/%s.fasta"%(outdir, ID)
    if file_is_empty(region_fasta) or replace is True: SeqIO.write([chrom[start:end] for chrom in SeqIO.parse(genome, "fasta") if chrom.id==chromosome], region_fasta, "fasta")

    # define the distance between ends of the reads, which depends on the length of the region
    if len_region<=(median_insert_size/2):

        # for small regions, set the median insert size to be half of the region, and the SD to be the 10% of the median
        median_insert_size = int(len_region/2)+1
        median_insert_size_sd = int(median_insert_size/10)+1

    # now run wgsim
    fastq_1 = "%s/%s_read1.fq"%(outdir, ID); fastq_2 = "%s/%s_read2.fq"%(outdir, ID); 
    fastq_1_tmp = "%s.tmp"%fastq_1; fastq_2_tmp = "%s.tmp"%fastq_2; 

    fastqgz_1 = "%s.gz"%fastq_1; fastqgz_2 = "%s.gz"%fastq_2;
    fastqgz_1_tmp = "%s.gz"%fastq_1_tmp; fastqgz_2_tmp = "%s.gz"%fastq_2_tmp;

 
    if any([file_is_empty(x) for x in {fastqgz_1, fastqgz_2}]) or replace is True:

        #print("simulating %s"%ID)

        # define the stderr of the reads
        std = "%s/%s_std.txt"%(outdir, ID)

        run_cmd("%s -e %.2f -N %i -1 %i -2 %i -r 0.0 -R 0.0 -X 0.0 -h -d %i -s %i %s %s %s > %s 2>&1"%(wgsim, error_rate, readPairs, read_length, read_length, median_insert_size, median_insert_size_sd, region_fasta, fastq_1_tmp, fastq_2_tmp, std))

        # check that the generated files are not empty, this may make the parallelization to fail
        for f in [fastq_1_tmp, fastq_2_tmp]:
            if os.stat(f).st_size==0: print("!!WARNING: No reads could be generated for region %s into %s. It is likely too short. This may cause problems with multiprocessing."%(ID, f))

        # delete previously generated gzips
        for f in [fastqgz_1_tmp, fastqgz_2_tmp]: remove_file(f)
        
        # compress        
        #print("%s performing gzip"%ID)
        for f in [fastq_1_tmp, fastq_2_tmp]: run_cmd("gzip %s"%f)
        os.rename(fastqgz_1_tmp, fastqgz_1); os.rename(fastqgz_2_tmp, fastqgz_2); 

        # check that everything is fine in the std and remove it
        if all([l.startswith("[wgsim") for l in open(std, "r").readlines()]): remove_file(std)
        else: raise ValueError("There may have been an error in generating some reads. Check %s"%std)


    # remove the fasta file
    remove_file(region_fasta)
    #print("%s returning"%ID)

    return (fastqgz_1, fastqgz_2)

def run_wgsim_pairedEnd_per_windows_in_parallel(df_windows, genome, outdir, read_length,  median_insert_size, median_insert_size_sd, replace=False, max_n_windows_at_once=500, error_rate=0.02):

    """Takes a dataframe with windows of ["chromosome", "start", "end", "readPairs"] and writes, under outdir, two fastq.gz files of simulated paired end  reads. The parallel runs are written in subfolder under outir which is finally removed. It returns a tuple of the two fastq.gz files.

    max_n_windows_at_once indicates the number of windows that can be processed at once. """

    # make the outdir if it does not exist
    make_folder(outdir)
    print("running simulation in a parallelized manner. Each window will be run in a different core. Running on chunks of %i windows"%max_n_windows_at_once)

    # define the final outdirs
    allChunks_fastqgz_1 = "%s/allChunks_reads1.fq.gz"%outdir; allChunks_fastqgz_2 = "%s/allChunks_reads2.fq.gz"%outdir

    if any([file_is_empty(f) for f in [allChunks_fastqgz_1, allChunks_fastqgz_2]]) or replace is True:

        # define the important files
        df_windows = df_windows[["chromosome", "start", "end", "readPairs"]]
        df_windows.index = list(range(len(df_windows)))

        # define chunks of indices
        chunks_indices = list(chunks(list(df_windows.index), max_n_windows_at_once))

        # initialize chunks files
        chunks_fastq_files = []

        # go through slices of max_n_windows_at_once of the df_windows
        for I, idxs in enumerate(chunks_indices):
            print("Working on chunk %i of %i"%(I+1, len(chunks_indices)))

            # define the dataframe of this run
            df_chunk = df_windows.loc[idxs]

            # define the final output
            chunk_all_fastqgz_1 = "%s/chunk%i_ofmax%i_all_reads1.fq.gz"%(outdir, I, max_n_windows_at_once); chunk_all_fastqgz_2 = "%s/chunk%i_ofmax%i_all_reads2.fq.gz"%(outdir, I, max_n_windows_at_once); 

            # define the location of the intermediate files
            parallel_files_outdir = "%s/parallel_files_generated_chunk%i_ofmax%i"%(outdir, I, max_n_windows_at_once)

            if any([file_is_empty(f) for f in [chunk_all_fastqgz_1, chunk_all_fastqgz_2]]) or replace is True:

                # define a folder where to store the parallel runs of this function
                delete_folder(parallel_files_outdir)
                make_folder(parallel_files_outdir)

                print("opening multiprocessing")

                # initialize the pool
                start_time =  time.time()
                with  multiproc.Pool(multiproc.cpu_count()) as pool:

                    # define the list of args
                    args = [(genome, chromosome, start, end, readPairs, read_length, parallel_files_outdir, median_insert_size, median_insert_size_sd, replace, error_rate) for chromosome, start, end, readPairs in df_chunk.values]

                    read1_read2_tuples = pool.starmap(run_wgsim_pairedEnd_for_window, args)
                    pool.close()
                    pool.terminate()
                    pool.join()
         
                # concatenate them all
                print("Integrating all in one for chunk %i..."%(I+1))
                chunk_all_fastqgz_1_tmp = "%s.tmp"%chunk_all_fastqgz_1; chunk_all_fastqgz_2_tmp = "%s.tmp"%chunk_all_fastqgz_2; 

                run_cmd("cat %s/*_read1.fq.gz > %s"%(parallel_files_outdir, chunk_all_fastqgz_1_tmp)) 
                run_cmd("cat %s/*_read2.fq.gz > %s"%(parallel_files_outdir, chunk_all_fastqgz_2_tmp))

                # rename to keep
                os.rename(chunk_all_fastqgz_1_tmp, chunk_all_fastqgz_1); os.rename(chunk_all_fastqgz_2_tmp, chunk_all_fastqgz_2)

                # print the time it took to generate these files
                print("--- creating reads for this chunk took %s seconds ---"%(time.time() - start_time))

            # remoove intermediate files of the fastq generation
            print("removing files")
            delete_folder(parallel_files_outdir)

            # keep the fastq files
            chunks_fastq_files += [chunk_all_fastqgz_1, chunk_all_fastqgz_2]

        # integrate all the chunks in one
        allChunks_fastqgz_1_tmp = "%s.tmp"%allChunks_fastqgz_1; allChunks_fastqgz_2_tmp = "%s.tmp"%allChunks_fastqgz_2
        chunks_read1_files = "%s/chunk*_ofmax%i_all_reads1.fq.gz"%(outdir, max_n_windows_at_once); chunks_read2_files = "%s/chunk*_ofmax%i_all_reads2.fq.gz"%(outdir, max_n_windows_at_once)

        run_cmd("cat %s > %s"%(chunks_read1_files, allChunks_fastqgz_1_tmp)) 
        run_cmd("cat %s > %s"%(chunks_read2_files, allChunks_fastqgz_2_tmp))

        # remove chunks' fastq's
        for f in chunks_fastq_files: remove_file(f)

        # at the end rename so that everything is clean
        os.rename(allChunks_fastqgz_1_tmp, allChunks_fastqgz_1)
        os.rename(allChunks_fastqgz_2_tmp, allChunks_fastqgz_2)

    print("All the reads simulated")

    return (allChunks_fastqgz_1, allChunks_fastqgz_2) 

def simulate_readPairs_per_window(df_windows, genome, npairs, outdir, read_length,  median_insert_size, median_insert_size_sd, replace=False, threads=4):

    """Simulates npairs of short reads in a way that is consistent with the predicted_relative_coverage of the df_windows"""

    # make the outdir if not there
    make_folder(outdir)

    # define the outputdirs
    all_fastqgz_1 = "%s/all_reads1.fq.gz"%outdir; all_fastqgz_2 = "%s/all_reads2.fq.gz"%outdir
    outdir_whole_chromosomes = "%s/whole_chromosomes_simulation_with_min_reads"%outdir
    outdir_per_windows = "%s/per_window_simulation_with_extra_reads"%outdir

    if any([file_is_empty(f) for f in [all_fastqgz_1, all_fastqgz_2]]) or replace is True:

        # add the chromosome lengths
        chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(genome, "fasta")}
        df_windows["len_chromosome"] = df_windows.chromosome.apply(lambda x: chr_to_len[x])

        # assign the number of readpairs to each window
        uniform_n_readPairs_per_base = npairs / sum(set(df_windows.len_chromosome))
        df_windows["window_length"] = df_windows.end - df_windows.start
        df_windows["total_readPairs"] = (df_windows.window_length*uniform_n_readPairs_per_base) * df_windows.predicted_relative_coverage

        # generate a df_windows that have sliding windows, so that each window appears twice and covers the intersection of df_windows. each window will get half of the reads that are expected
        df_windows_sliding = pd.DataFrame()

        for chromosome in df_windows.chromosome.unique():
            print("geting sliding window for %s"%chromosome)

            # get the df with a unique numeric index
            df_c = df_windows[df_windows.chromosome==chromosome][["chromosome", "start", "end", "total_readPairs", "window_length"]].sort_values(by=["chromosome", "start"])
            df_c.index = list(range(len(df_c)))

            ##### define df_c_halfs, which has non_overlapping windows ####

            # get a df that has halfs of the windows, only from the first to the forelast ones --> this lacks the half of the first window and half of the last window
            df_c_halfs = cp.deepcopy(df_c).iloc[0:-1][["chromosome", "start", "end", "window_length"]]

            # define the starts
            starts = list(df_c_halfs.apply(lambda r: r["start"] + int(r["window_length"]*0.5), axis=1)) # start at +half of the window
            starts[0] = 0 # the first start is 0

            # define the ends (half of the next window)
            ends = [df_c_halfs.loc[I, "end"] + int(df_c.loc[I+1, "window_length"]*0.5) for I in df_c_halfs.index]
            ends[-1] = df_c.iloc[-1].end

            # add to df
            df_c_halfs["start"] = starts; df_c_halfs["end"] = ends
            df_c_halfs["window_length"] = df_c_halfs.end - df_c_halfs.start

            #################################################################

            # add to df_halfs the expected coverage as a function of the overlapping windows 
            def get_reads_from_overlap_with_other_windows(r):

                """Takes a row and of the windows df and returns the reads as a function of the overlap with df_c"""

                # get the index of this window, which is -500 than df_c
                I = r.name

                # for the first window, take the first window and half of the second
                if I==0: return (df_c.loc[0].total_readPairs + df_c.loc[1].total_readPairs/2)

                # for the last window, get half of the forelast window and all the last one 
                if I==last_w: return (df_c.iloc[-2].total_readPairs/2 + df_c.iloc[-1].total_readPairs)

                # for the others, take half and half
                else: return (df_c.loc[I, "total_readPairs"]/2 + df_c.loc[I+1, "total_readPairs"]/2)

            last_w = df_c_halfs.index[-1]
            df_c_halfs["total_readPairs"] = df_c_halfs[["start", "end"]].apply(get_reads_from_overlap_with_other_windows, axis=1)
            df_c_halfs = df_c_halfs[list(df_c.keys())]

            # get the concatenated df
            df_c_both = df_c_halfs.append(df_c, sort=True).sort_values(by=["chromosome", "start"])

            # get half of the total_readPairs, as you have overlapping region
            df_c_both["total_readPairs"] = (df_c_both.total_readPairs/2).apply(int) + 1 # I also add a pseudocount

            # record the minimum 
            min_read_pairs = min(df_c_both.total_readPairs)
            df_c_both["min_readPairs"] = min_read_pairs
            df_c_both["extra_readPairs"] = (df_c_both.total_readPairs - min_read_pairs) + 1 # also add a pseudocount to have non-empty regions

            # append both dfs in a sorted manner
            df_windows_sliding = df_windows_sliding.append(df_c_both, sort=True)

        # get chr as index
        df_windows_sliding = df_windows_sliding.set_index("chromosome", drop=False)

        # first simulate reads for the whole chromosome, so that we simulate the minimum number of reads for each chromosome. This is useful to have paired end reads that span the whole chromosome, not only regions
        chrom_to_len = {s.id: len(s.seq) for s in SeqIO.parse(genome, "fasta")}
        df_windows_chromosomes_minCov = pd.DataFrame({I : {"start": 0, "end": len_chr, "chromosome": chrom, "readPairs": sum(df_windows_sliding.loc[chrom, "min_readPairs"])} for I, (chrom, len_chr) in enumerate(chrom_to_len.items())}).transpose()
        wholeChr_fastqgz1, wholeChr_fastqgz2 = run_wgsim_pairedEnd_per_windows_in_parallel(df_windows_chromosomes_minCov, genome, outdir_whole_chromosomes, read_length, median_insert_size, median_insert_size_sd, replace=replace)

        # now simulate the reads for the regions, only considering the extra_readPairs
        df_windows_sliding["readPairs"] = df_windows_sliding.extra_readPairs        
        perWindow_fastqgz1, perWindow_fastqgz2 = run_wgsim_pairedEnd_per_windows_in_parallel(df_windows_sliding, genome, outdir_per_windows, read_length, median_insert_size, median_insert_size_sd, replace=replace)

        # now concatenate all the simulated reads
        print("Integrating all in one...")
        all_fastqgz_1_tmp = "%s.tmp"%all_fastqgz_1; all_fastqgz_2_tmp = "%s.tmp"%all_fastqgz_2; 

        run_cmd("cat %s %s > %s"%(wholeChr_fastqgz1, perWindow_fastqgz1, all_fastqgz_1_tmp))
        run_cmd("cat %s %s > %s"%(wholeChr_fastqgz2, perWindow_fastqgz2, all_fastqgz_2_tmp))

        # rename to keep
        os.rename(all_fastqgz_1_tmp, all_fastqgz_1); os.rename(all_fastqgz_2_tmp, all_fastqgz_2)

    # remove previously generated folders
    delete_folder(outdir_whole_chromosomes); delete_folder(outdir_per_windows)

    # return the simulated reads
    return (all_fastqgz_1, all_fastqgz_2) 

def simulate_and_align_PairedReads_perWindow(df_windows, genome_interest, reference_genome, npairs, read_length, outdir, median_insert_size, median_insert_size_sd, replace=False, threads=4):

    """Takes a dataframe with windows of the genome, which also has a predicted_relative_coverage (which indicates by how much should the coverage be multiplied in this window). This function generates a fastq (and deletes it at the end), aligns it and returns the bam. All files are written under outdir. It returns the aligned bamfile. All the chromosomes are simulated as linear."""

    print("Simulating reads and aliging them for %s"%genome_interest)
    start_time = time.time()

    # make folder if it does not exist
    make_folder(outdir)

    # define final outputs
    sim_bamfile = "%s/aligned_reads%ipairs_readlen%i_insertsize%i+-%i.bam"%(outdir, npairs, read_length, median_insert_size, median_insert_size_sd)
    sim_sorted_bam = "%s.sorted"%sim_bamfile
    sim_index_bam = "%s.bai"%sim_sorted_bam

    if file_is_empty(sim_sorted_bam) or replace is True:

        # sort by chromosome and start
        df_windows = df_windows.sort_values(by=["chromosome", "start"])
      
        ##### generate the reads for the genome_interest #####
        outdir_reads = "%s/simulated_reads_%ipairs_readlen%i_insertsize%i+-%i"%(outdir, npairs, read_length, median_insert_size, median_insert_size_sd) 
        read1_fastqgz, read2_fastqgz = simulate_readPairs_per_window(df_windows, genome_interest, npairs, outdir_reads, read_length, median_insert_size, median_insert_size_sd, replace=replace, threads=threads) 

        ######################################################

        ##### align the reads and retun the bam ######
        run_bwa_mem(read1_fastqgz, read2_fastqgz, reference_genome, outdir, sim_bamfile, sim_sorted_bam, sim_index_bam, name_sample="simulations_reference_genome", threads=threads, replace=replace)

        # remove the fastq files
        print("deleting reads")
        delete_folder(outdir_reads)

    # record the time consumed
    print("--- generating %i reads from %s and aligning them  took %s seconds ---"%(npairs, genome_interest, time.time() - start_time))

    return sim_sorted_bam

def get_best_parameters_for_GridssClove_run(sorted_bam, reference_genome, outdir, threads=4, replace=False, window_l=5000, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", expected_ploidy=1, nvars=100, real_svtype_to_file={}, median_insert_size=250, median_insert_size_sd=0):

    """This finds the optimum parameters for running GRIDSS clove and returns them. The parameters are equivalent to the run_GridssClove_optimising_parameters function"""


    # define plots dir
    PlotsDir = "%s/plots"%outdir; make_folder(PlotsDir)

    ###### MODELLING COVERAGE ######
    print("modelling coverage of the sample")

    # get a graph of the genome without any breakpoints (thus the None, None)
    genomeGraph_outfileprefix = "%s/genomeGraph_withoutBPs"%(outdir)
    genome_graph, df_positions_graph = get_genomeGraph_object(reference_genome, None, None, genomeGraph_outfileprefix, replace=replace)

    # get a function that takes the GC content, chromosome and distance to the telomere and returns coverage. This is actually a lambda function
    outdir_coverage_calculation = "%s/coverage_per_regions%ibb"%(outdir, window_l); make_folder(outdir_coverage_calculation)
    df_coverage_train = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_coverage_calculation, sorted_bam, windows_file="none", replace=replace, window_l=window_l), sep="\t")

    distToTel_chrom_GC_to_coverage_fn = get_distanceToTelomere_chromosome_GCcontent_to_coverage_fn(df_coverage_train, reference_genome, genome_graph, df_positions_graph, outdir_coverage_calculation, mitochondrial_chromosome=mitochondrial_chromosome, replace=replace)

    print("coverage model obtained")

    ################################

    ############ GENERAL OPERATIONS THAT WILL BE NEEDED FOR ALL THE STEPS #####

    # the dir and genome names
    genome_dir = "/".join(reference_genome.split("/")[0:-1])
    genome_name = reference_genome.split("/")[-1].split(".")[0]

    # map each chromosome to length
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

    # count the length od the reads
    read_length = get_read_length(sorted_bam, threads=threads, replace=replace)
    print("The median read length is %i"%read_length)

    # count total number of reads
    total_nread_pairs = count_number_read_pairs(sorted_bam, replace=False, threads=threads)
    #total_nread_pairs  = 1000000 # this is to debug the simulation pipeline
    expected_coverage_per_bp = int((total_nread_pairs*read_length) / sum(chr_to_len.values())) +  1 # the expected coverage per position with pseudocount
    print("There are %i read pairs in your library. The expected coverage is %ix."%(total_nread_pairs, expected_coverage_per_bp))

    # get the info of the reference genome with predictions of coverage per window
    df_REFgenome_info = get_windows_infoDF_with_predictedFromFeatures_coverage(reference_genome, distToTel_chrom_GC_to_coverage_fn, expected_coverage_per_bp, replace=replace, window_l=window_l, threads=threads)

    # simulate reads for the reference if you are not only simulating haploid
    if set(simulation_ploidies)!={"haploid"}: 
        outdir_ref = "%s/simulation_reference_genome_%ibp_windows"%(outdir, window_l)
        simulated_reference_bam_file = simulate_and_align_PairedReads_perWindow(df_REFgenome_info, reference_genome, reference_genome, total_nread_pairs, read_length, outdir_ref, median_insert_size, median_insert_size_sd, replace=replace, threads=threads)

    ##############################################################################

    ################ SIMULATION PIPELINE ################ 

    # initialize a df with all the benchmarking data
    df_benchmark_all = pd.DataFrame()
    genomeID_to_knownSVdict = {}

    df_benchmark_all_file = "%s/df_benchmark_all.py"%outdir
    genomeID_to_knownSVdict_file= "%s/genomeID_to_knownSVdict.py"%outdir

    if file_is_empty(df_benchmark_all_file) or file_is_empty(genomeID_to_knownSVdict_file) or replace is True:

        # go throigh each simulation (these are technical replicates of the pipeline)
        for simulation_ID in range(1, n_simulated_genomes+1):
            print("working on simulation %i"%simulation_ID)

            # get an outdir where all the simulations of this ID will be stored
            simulation_outdir = "%s/simulation_%i"%(outdir, simulation_ID); make_folder(simulation_outdir)

            # get the simulated SVs, which are an integration of 
            sim_svtype_to_svfile, rearranged_genome = rearrange_genomes_simulateSV(reference_genome, simulation_outdir, replace=replace, nvars=nvars, mitochondrial_chromosome=mitochondrial_chromosome, simulated_svtype_to_svfile=real_svtype_to_file, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"})

            kjhckgadkjgkadg



    ####################################################








    lndalshjldjasjjskdllkjskjlsakjlasdkjladsjklasdjklasdjkljklsjklsad

    return gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup


def run_GridssClove_optimising_parameters(sorted_bam, reference_genome, outdir, threads=4, replace=False, window_l=5000, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", expected_ploidy=1, nvars=100, fast_SVcalling=False, real_svtype_to_file={}):

    """
    Takes some aligned reads and runs the GridssPipeline optimising the parameters of GRIDSS filtering. These are the different parameters of the function:

    - sorted_bam: the path to a sorted and indexed bam where we want to find the SV
    - reference_genome: the fasta of the reference genome
    - outdir: a directory where all the sample-specific files will be written. All the reference_genome-specific files will be written in the same dir as the reference_genome is
    - window_l: the length of the windows to generate the coverage model
    - n_simulated_genomes is the number of genomes that will be simulated as simulation replicates
    - mitochondrial_chromosome is the name of the mtDNA chromosome in the reference genome. This is passed to some parts of the pipeline to consider differently the gDNA and the mtDNA (the SV simulation). It can be comma-sepparated. It can be "no_mitochondria" if there is no mitochondrial chromosome
    - simulation_ploidies indicates which poplulations or ploidies have to be simulated. 2ref_1sv means that we will simulate a genome that has 2 reference genomes and 1 genome under structural variation (this can be a population sequencing)
    - range_filtering_benchmark indicates which type of simulation will be performed, it can be "large", "medium", "small", "single", "theoretically_meaningful". This is passed to benchmark_GridssClove_for_knownSV.
    - expected_ploidy is a number that states the expected ploidy. 
    - nvars deteremines the number of SVs to simulate in each simulation of each type. The mtDNA will get 5% of these. There will be as maximum len(gDNA chromosomes)-1 balanced translocations, so that each chromosomal arm is only implicated once.
    - fast_SVcalling runs the SV on a predefined set of parameters, without optimisation
    - real_svtype_to_file is a dict that maps each type of 

    """

    # make the outdir
    make_folder(outdir)

    # initialize the start time
    pipeline_start_time = time.time()

    ###### CALCULATE GENERAL PARMS ######

    # clean the reference genome windows files
    clean_reference_genome_windows_files(reference_genome)

    # calculate the insert size
    median_insert_size, median_insert_size_sd  = get_insert_size_distribution(sorted_bam, replace=replace, threads=threads)
    print("The median insert size is %i, with an absolute deviation of %i"%(median_insert_size, median_insert_size_sd))

    #####################################

    ###### GET PARAMETERS ######

    # get the default running and filtering parameters
    if fast_SVcalling is False:
        print("getting optimised parameters")

        parameter_optimisation_dir = "%s/parameter_optimisation"%outdir; make_folder(parameter_optimisation_dir)

        gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup = get_best_parameters_for_GridssClove_run(sorted_bam, reference_genome, parameter_optimisation_dir, threads=threads, replace=replace, window_l=window_l, n_simulated_genomes=n_simulated_genomes, mitochondrial_chromosome=mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=range_filtering_benchmark, expected_ploidy=expected_ploidy, nvars=nvars, real_svtype_to_file=real_svtype_to_file, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd)

    # get the parameters from an optimisation
    else: 
        print("running with default un-optimised parameters")

        # get the default parameters 
        gridss_blacklisted_regions = ""
        gridss_maxcoverage = 50000
        gridss_filters_dict = default_filtersDict_gridss
        max_rel_coverage_to_consider_del = 0.1
        min_rel_coverage_to_consider_dup = 1.8

    ###########################

    ###### FINAL GRIDSS-CLOVE RUNNING ######

    # define the final outdir 
    outdir_gridss_final = "%s/final_gridss_running"%outdir; make_folder(outdir_gridss_final)
    final_gridss_vcf = "%s/output_gridss.vcf"%outdir_gridss_final

    # define the median coverage across window_l windows of the genome
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_gridss_final, sorted_bam, windows_file="none", replace=replace, window_l=window_l), sep="\t")
    median_coverage = np.median(coverage_df[~coverage_df["#chrom"].isin(mitochondrial_chromosome.split(","))].mediancov_1)
    print("The median coverage is %i"%median_coverage)

    # run the pipeline
    print("running final gridss with parameters...")
    run_gridssClove_given_filters(sorted_bam, reference_genome, outdir_gridss_final, median_coverage, replace=replace, threads=threads, gridss_blacklisted_regions=gridss_blacklisted_regions, gridss_VCFoutput=final_gridss_vcf, gridss_maxcoverage=gridss_maxcoverage, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, gridss_filters_dict=gridss_filters_dict, run_in_parallel=True, max_rel_coverage_to_consider_del=max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup=min_rel_coverage_to_consider_dup, replace_FromGridssRun=replace, define_insertions_based_on_coverage=False)

    ########################################

    ##### PIPELINE ENDING OPERATIONS ##### 

    # at the end, remove all the mosdepth and windows files under the reference
    clean_reference_genome_windows_files(reference_genome)
    
    print("GRIDSS pipeline finished correctly")

    print("--- the gridss pipeline optimising parameters took %s seconds in %i cores ---"%(time.time() - pipeline_start_time, threads))

    # at the end remove all the files that are unnecessary

    # simulations
    #for simulation_ID in range(1, n_simulated_genomes+1): delete_folder("%s/simulation_%i"%(outdir, simulation_ID))

    # remove the cross-benchmarking files
    #delete_folder("%s/cross_benchmarking_files"%outdir_benchmarking)

    # generate a file that indicates whether the gridss run is finished
    #final_file = "%s/gridss_finished_file_final_gridss_final_with_window_size.txt"%outdir
    #open(final_file, "w").write("gridss finished...")

    ######################################




