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
import time
from sklearn.metrics import r2_score
from collections import Counter, defaultdict
import inspect
import collections
from shutil import copyfile
import igraph
from ete3 import Tree, NCBITaxa
import urllib
from subprocess import STDOUT, check_output
import subprocess
import subprocess, datetime, signal
import json
import sklearn
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import scipy.stats as stats

#### UNIVERSAL FUNCTIONS ####

def get_fullpath(x):

    """Takes a path and substitutes it bu the full path"""

    if x.startswith("/"): return x
    elif x.startswith("."): path = "/".join(x.split("/")[1:])
    else: path = x

    return os.getcwd() + "/" + path


def get_date_and_time():

    """Gets the date of today"""

    current_day = date.today().strftime("%d%m%Y")
    current_time = time.strftime("%H%M%S", time.localtime())

    return "%s_%s"%(current_day, current_time)

########################

# packages to remove (potentially)
#import psutil
#from statsmodels.stats import multitest

# chnage the errirs ti report
warnings.simplefilter(action='ignore', category=pd.core.common.SettingWithCopyWarning) # avoid the slicing warning
#pd.options.mode.chained_assignment = 'raise'
warnings.simplefilter('ignore', sklearn.metrics.regression.UndefinedMetricWarning)
warnings.simplefilter('ignore', np.RankWarning)

np.seterr(divide='ignore', invalid='ignore')



# load a specific matplotlib library for cluster envs
try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

except: import matplotlib.pyplot as plt

import seaborn as sns

# get the cwd were all the scripts are 
CWD = get_fullpath("/".join(__file__.split("/")[0:-1])); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# define the conda base dir and the env name
CondaDir =  "/".join(sys.executable.split("/")[0:-4])
EnvName = EnvDir.split("/")[-1]


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
repeat_masker = "%s/bin/RepeatMasker"%EnvDir
vcffilter = "%s/bin/vcffilter"%EnvDir
freebayes_parallel = "%s/bin/freebayes-parallel"%EnvDir
fasta_generate_regions_py = "%s/bin/fasta_generate_regions.py"%EnvDir
wgsim = "%s/bin/wgsim"%EnvDir
picard_exec = "%s/bin/picard"%EnvDir
svim = "%s/bin/svim"%EnvDir
bbmap_reformat_sh = "%s/bin/reformat.sh"%EnvDir
mosdepth = "%s/bin/mosdepth"%EnvDir
repeat_modeller_BuildDatabase = "%s/bin/BuildDatabase"%EnvDir
perl = "%s/bin/perl"%EnvDir
makeblastdb = "%s/bin/makeblastdb"%EnvDir
abblast_dir = "%s/bin"%EnvDir
cdhit_dir = "%s/bin"%EnvDir
genometools_dir = "%s/bin"%EnvDir
ltr_retriever_dir = "%s/bin"%EnvDir
mafft_dir = "%s/bin"%EnvDir
recon_dir = "%s/bin"%EnvDir
rmblast_dir = "%s/bin"%EnvDir
rscout_dir = "%s/bin"%EnvDir
trf_prgm_dir = "%s/bin/trf"%EnvDir
ninja_dir = "%s/bin"%EnvDir
esearch = "%s/bin/esearch"%EnvDir
efetch = "%s/bin/efetch"%EnvDir
prefetch = "%s/bin/prefetch"%EnvDir
fastqdump = "%s/bin/fastq-dump"%EnvDir
parallel_fastq_dump = "%s/bin/parallel-fastq-dump"%EnvDir
FASTQC = "%s/bin/fastqc"%EnvDir
porechop = "%s/bin/porechop"%EnvDir
seqtk = "%s/bin/seqtk"%EnvDir
repeatmoder_dir = "%s/share/RepeatModeler"%EnvDir
repeat_modeller = "%s/bin/RepeatModeler"%EnvDir
repeatmasker_dir = "%s/share/RepeatMasker"%EnvDir
fasta_generate_regions_py = "%s/bin/fasta_generate_regions.py"%EnvDir
pigz = "%s/bin/pigz"%EnvDir
unpigz = "%s/bin/unpigz"%EnvDir

# executables that are provided in the repository
external_software = "%s/../installation/external_software"%CWD
gridss_run = "%s/gridss.sh"%external_software
gridss_jar = "%s/gridss-2.9.2-gridss-jar-with-dependencies.jar"%external_software
clove = "%s/clove-0.17-jar-with-dependencies.jar"%external_software
gztool = "%s/gztool"%external_software

# define the bcftools=1.10 by activating the conda env
#SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
# CONDA_ACTIVATING_CMD = "conda activate %s;"%EnvName
#bcftools_latest = "%s && conda activate %s_bcftools_1.10.2_env && bcftools"%(SOURCE_CONDA_CMD, EnvName)
#bcftools_latest_cmd = "%s && conda activate %s_bcftools_1.10.2_env && bcftools"%(SOURCE_CONDA_CMD, EnvName)


# scripts that are of this pipeline
create_random_simulatedSVgenome_R = "%s/create_random_simulatedSVgenome.R"%CWD
create_targeted_simulatedSVgenome_R = "%s/create_targeted_simulatedSVgenome.R"%CWD
annotate_simpleEvents_gridssVCF_R = "%s/annotate_simpleEvents_gridssVCF.R"%CWD
analyze_svVCF = "%s/generate_files_from_svVCF.R"%CWD
analyze_svVCF_simple = "%s/generate_files_from_svVCF_simple.R"%CWD
TRIMMOMATIC = "%s/run_trimmomatic.py"%CWD 
perSVade_py = "%s/perSVade.py"%CWD 
run_trimmomatic_and_fastqc_py = "%s/run_trimmomatic_and_fastqc.py"%CWD
get_trimmed_reads_for_srr_py = "%s/get_trimmed_reads_for_srr.py"%CWD
run_vep = "%s/run_vep.py"%CWD

######################################################
######################################################

####################################
######## DEFINE VARIABLES ##########
####################################

# define the strings that have to be considered as NaN in the VCF parsing
vcf_strings_as_NaNs = ['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'n/a', 'nan', 'null']

# define default parameters for gridss filtering. This has changed from v0
default_filtersDict_gridss = {"min_Nfragments":5, "min_af":0.25, "wrong_FILTERtags":("NO_ASSEMBLY",), "filter_polyGC":True, "filter_noSplitReads":False, "filter_noReadPairs":False, "maximum_strand_bias":0.99, "maximum_microhomology":50, "maximum_lenght_inexactHomology":50, "range_filt_DEL_breakpoints":(0, 1), "min_length_inversions":40, "dif_between_insert_and_del":5, "max_to_be_considered_small_event":1000, "wrong_INFOtags":('IMPRECISE',), "min_size":50, "min_af_EitherSmallOrLargeEvent":0.25, "min_QUAL":0, "filter_overlappingRepeats":False} # the minimum af is 0.25 to include both heterozygous and homozygous vars as default

# define other default parameters
default_gridss_blacklisted_regions = ""
default_gridss_maxcoverage = 50000
default_max_rel_coverage_to_consider_del = 0.1
default_min_rel_coverage_to_consider_dup = 1.8

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
g_min_QUAL_l = [0, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
g_filter_overlappingRepeats_l = [False, True]

# map each filter to the ordered list
g_filterName_to_filtersList = {"min_Nfragments":g_min_Nfragments_l, "min_af":g_min_af_l, "wrong_FILTERtags":g_wrong_FILTERtags_l, "filter_polyGC":g_filter_polyGC_l, "filter_noSplitReads":g_filter_noSplitReads_l, "filter_noReadPairs":g_filter_noReadPairs_l, "maximum_strand_bias":g_maximum_strand_bias_l, "maximum_microhomology":g_maximum_microhomology_l, "maximum_lenght_inexactHomology":g_maximum_lenght_inexactHomology_l, "range_filt_DEL_breakpoints":g_range_filt_DEL_breakpoints_l, "min_length_inversions":g_min_length_inversions_l, "dif_between_insert_and_del":g_dif_between_insert_and_del_l, "max_to_be_considered_small_event":g_max_to_be_considered_small_event_l, "wrong_INFOtags":g_wrong_INFOtags_l, "min_size":g_min_size_l, "min_af_EitherSmallOrLargeEvent":g_min_af_EitherSmallOrLargeEvent_l, "min_QUAL":g_min_QUAL_l, "filter_overlappingRepeats":g_filter_overlappingRepeats_l}

# map each value of each filter list to a value depending on the position
g_filterName_to_filterValue_to_Number = {filterName : dict(zip(filtersList, range(len(filtersList)))) for filterName, filtersList in g_filterName_to_filtersList.items()}

# define the default window_l
window_l = 10000

# define the verbosity. If False, none of the print statements will have an effect.
printing_verbose_mode = True

# define the wrong SRRs (could not be downloaded)
wrong_SRRs = {"ERR2163728", "SRR8373447"}

# define the fraction of RAM to dedicate to the pipeline
fractionRAM_to_dedicate = 0.8

svtype_to_color={"tandemDuplications": "gray", "deletions": "black", "inversions": "blue", "translocations": "olive", "insertions": "red", "remaining":"magenta", "integrated":"c"}

# define a dict that maps each svtype to the fields that are important to define the overlaps
svtype_to_fieldsDict = {"inversions": {"equal_fields": ["Chr"], 
                                        "approximate_fields": ["Start", "End"],
                                        "chromField_to_posFields": {"Chr":{"start": "Start", "end": "End"}},
                                        "all_fields": ["Chr", "Start", "End", "ID"],
                                        "position_fields": ["Start", "End"],
                                        "positionField_to_chromosome": {"Start":"Chr", "End":"Chr"},
                                        "chromosome_fields": ["Chr"]
                                        }, 

                        "tandemDuplications": {"equal_fields": ["Chr"], 
                                               "approximate_fields": ["Start", "End"],
                                                "chromField_to_posFields": {"Chr":{"start": "Start", "end": "End"}},
                                                "all_fields": ["Chr", "Start", "End", "ID"],
                                                "position_fields": ["Start", "End"],
                                                "positionField_to_chromosome": {"Start":"Chr", "End":"Chr"},
                                                "chromosome_fields": ["Chr"]
                                                }, 

                        "deletions": {"equal_fields": ["Chr"], 
                                      "approximate_fields": ["Start", "End"],
                                      "chromField_to_posFields": {"Chr":{"start": "Start", "end": "End"}},
                                      "all_fields": ["Chr", "Start", "End", "ID"],
                                      "position_fields": ["Start", "End"],
                                      "positionField_to_chromosome": {"Start":"Chr", "End":"Chr"},
                                      "chromosome_fields": ["Chr"]
                                     }, 

                        "insertions": {"equal_fields": ["ChrA", "ChrB", "Copied"], 
                                        "approximate_fields": ["StartA", "EndA", "StartB", "EndB"],
                                        "chromField_to_posFields": {"ChrA":{"start": "StartA", "end": "EndA"}},
                                        "all_fields": ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Copied", "ID"],
                                        "position_fields": ["StartA", "EndA", "StartB", "EndB"],
                                        "positionField_to_chromosome": {"StartA":"ChrA", "EndA":"ChrA", "StartB":"ChrB", "EndB":"ChrB"},
                                        "chromosome_fields": ["ChrA", "ChrB"]
                                        }, 

                        "translocations": {"equal_fields": ["ChrA", "ChrB", "Balanced"], 
                                           "approximate_fields": ["StartA", "EndA", "StartB", "EndB"],
                                           "chromField_to_posFields": {},
                                           "all_fields": ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Balanced", "ID"],
                                           "position_fields": ["StartA", "EndA", "StartB", "EndB"],
                                           "positionField_to_chromosome": {"StartA":"ChrA", "EndA":"ChrA", "StartB":"ChrB", "EndB":"ChrB"},
                                           "chromosome_fields": ["ChrA", "ChrB"]
                                          }, 

                        "remaining": {"equal_fields": ["#CHROM", "CHR2", "SVTYPE"], 
                                      "approximate_fields": ["POS", "START", "END"],
                                      "chromField_to_posFields": {"CHR2":{"start":"START", "end":"END"}},
                                      "all_fields": ["#CHROM", "POS", "CHR2", "START", "END", "SVTYPE", "ID"],
                                      "position_fields": ["POS", "START", "END"],
                                      "positionField_to_chromosome": {"POS":"#CHROM", "START":"CHR2", "END":"CHR2"},
                                      "chromosome_fields": ["#CHROM", "CHR2"]

                                      }}

# define the golden set reads
taxID_to_srrs_goldenSet = {3702: {"short_reads":"ERR2173372", "long_reads":"ERR2173373"}} # arabidopsis thaliana sample

####################################
####################################
####################################


def print_if_verbose(*x):

    """This function overrides the print statement"""

    if printing_verbose_mode is True: print(*x)

def get_aa(codons, genetic_code):

    """Takes a string of codons and returns the AA, according to the genetic_code"""

    # if there are no codons
    if codons=="-": return "-"

    # if it is a protein
    elif len(codons)%3==0: return str(Seq(codons).translate(table = genetic_code))

    # if not
    else: 
        if len(codons)<3: return "X"
        else: return (str(Seq("".join(list(chunks(codons, 3))[0:-1])).translate(table = genetic_code)) + "X")

def modify_DF_cols(row, genetic_code, stop_codons, genCode_affected_vars):

    """Takes a row of a VEP df according to the genetic_code and modifies the necessary rows"""

    if row["Codons"]=="-": codons_ref = "-"; codons_alt = "-"

    else:
        # get the codons
        codons_ref, codons_alt = row["Codons"].split("/")
        codons_ref = codons_ref.upper(); codons_alt = codons_alt.upper(); 

    # get refs
    aa_ref = get_aa(codons_ref, genetic_code); aa_alt = get_aa(codons_alt, genetic_code);

    # define the consequences
    consequences = set()

    # retainig stops
    if codons_ref in stop_codons and codons_alt in stop_codons: consequences.add('stop_retained_variant')

    # indels
    if (len(aa_ref.strip("-"))>len(aa_alt.strip("-")) and len(codons_ref)%3==0 and len(codons_alt)%3==0) or (aa_ref not in {"*", "X"} and aa_alt=="-"): consequences.add('inframe_deletion')
    if (len(aa_ref.strip("-"))<len(aa_alt.strip("-")) and len(codons_ref)%3==0 and len(codons_alt)%3==0) or (aa_alt not in {"*", "X"} and aa_ref=="-"): consequences.add('inframe_insertion')

    # frameshift
    if "X" in aa_alt: consequences.add('frameshift_variant')

    # SNPs
    if len(codons_ref)==3 and len(codons_alt)==3 and aa_ref!="*" and aa_alt!="*":
        if aa_ref==aa_alt: consequences.add('synonymous_variant')
        else: consequences.add('missense_variant')

    # stop codon info
    nStops_ref = aa_ref.count("*"); nStops_alt = aa_alt.count("*"); 
    if nStops_ref<nStops_alt: consequences.add('stop_gained')
    if nStops_ref>nStops_alt: consequences.add('stop_lost')

    # protein altering. This is equivalent to the thing of the row already. 
    if "protein_altering_variant" in row["Consequence"]: consequences.add('protein_altering_variant')

    # add any consequence that is not in any of the genCode_affected_vars
    consequences.update(set(row["Consequence"].split(",")).difference(genCode_affected_vars))

    # return the aa and the consequences
    return pd.Series({"Amino_acids": "%s/%s"%(aa_ref, aa_alt), "Consequence": ",".join(list(consequences))})

def id_generator(size=10, chars=string.ascii_uppercase + string.digits, already_existing_ids=set()):

    """ already_existing_ids is a set that indicates whihc IDs can't be picked """

    ID = ''.join(random.choice(chars) for _ in range(size))
    while ID in already_existing_ids:
        ID = ''.join(random.choice(chars) for _ in range(size))

    return ID

def find_nearest(a, a0):
    """Element in nd array `a` closest to the scalar value `a0`"""
    
    if type(a0) in {float, int}:

        # Debug elements that are inf
        if a0 not in [np.inf, -np.inf]:
            a = np.array(a)
            idx = np.abs(a - a0).argmin()
            closest_in_a = a.flat[idx]
            
        elif a0==np.inf:
            closest_in_a = max(a)
            
        elif a0==-np.inf:
            closest_in_a = min(a)        

    elif type(a0) in {bool, str, tuple}:

        # return the value in a that has the exact same value
        if a0 in a: closest_in_a = a0
        else: raise ValueError("%s can't be found in %s"%(a0, a))

    else: raise ValueError("%s is not a valid type for find_nearest function"%(a0))

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

    # define the way to export the path to be first the EnvDir
    #path_definition_cmd = "export PATH=%s/bin:$PATH"%EnvDir
    #cmd_to_run = "%s && %s"%(path_definition_cmd, cmd) # the cmd with the path definition
    cmd_to_run = cmd # the raw cmd
    out_stat = os.system(cmd_to_run) 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd, out_stat))


def run_bcftools_latest(kwargs):

    """This funcion runs bcftools 1.10.2 with the provided kwargs"""

    # define the cmds
    SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
    bcftools_cmd = "%s && conda activate %s_bcftools_1.10.2_env && bcftools"%(SOURCE_CONDA_CMD, EnvName)

    # define the running
    run_cmd("bash -c '%s %s'"%(bcftools_cmd, kwargs))

def get_dir(filename): return "/".join(filename.split("/")[0:-1])

def get_file(filename): return filename.split("/")[-1]

def save_object(obj, filename):
    
    """ This is for saving python objects """

    filename_tmp = "%s.tmp"%filename
    remove_file(filename_tmp)
    
    with open(filename_tmp, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

    os.rename(filename_tmp, filename)

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

    if os.path.isfile(f): 

        try: run_cmd("rm %s > /dev/null 2>&1"%f)
        except: pass

def delete_folder(f):

    if os.path.isdir(f): shutil.rmtree(f)


def make_folder(f):

    if not os.path.isdir(f): os.mkdir(f)

def delete_file_or_folder(f):

    """Takes a path and removes it"""

    if os.path.isdir(f): shutil.rmtree(f)
    if os.path.isfile(f): os.unlink(f)

def get_chr_to_len(genome):

    chr_to_len_file = "%s.chr_to_len.py"%genome
    chr_to_len_file_tmp = "%s.tmp"%chr_to_len_file

    if file_is_empty(chr_to_len_file):

        remove_file(chr_to_len_file_tmp)

        # define chromosome_to_length for a genome
        chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(genome, "fasta")}

        # save
        save_object(chr_to_len, chr_to_len_file_tmp)
        os.rename(chr_to_len_file_tmp, chr_to_len_file)

    else: chr_to_len = load_object(chr_to_len_file)

    return chr_to_len

def get_uniqueVals_df(df): return set.union(*[set(df[col]) for col in df.columns])




def get_date():

    """Gets the date of today"""

    today = date.today()

    return today.strftime("%d/%m/%Y")

def extract_BEDofGENES_of_gff3(gff, bed, replace=False, reference=""):

    """Takes the full path to a gff annotation and writes a bed file to bed, which only has genes' info. 
    The bed will be 1-indexed, which is not usually the case. Bed files are 0 inexed. 

    It also writes a file called 'bed.regions_lengthwindow_l' that has the regions of - and + window_l of the start and end of the gene, respectively.
    The regions bed is returned by this function.

    """

    # load the gff
    df_gff3 = pd.read_csv(gff, skiprows=list(range(len([line for line in open(gff, "r") if line.startswith("#")]))), sep="\t", names=["chromosome", "source", "type_feature", "start", "end", "score", "strand", "phase", "attributes"])

    # define the regions you are interested in 
    interesting_features = {"gene"}
    df_gff3 = df_gff3[df_gff3.type_feature.isin(interesting_features)]

    # define a function that takes attribues and returns ID
    def get_ID_gff(attributes):

        # get the ID
        IDlist = [x.lstrip("ID=") for x in attributes.split(";") if x.startswith("ID=")]

        # check that the ID is correct
        if len(IDlist)!=1: 
            print_if_verbose(IDlist, attributes)
            raise ValueError("IDlist has to be of length 1")

        # get the ID
        ID = IDlist[0]

        # add the part if it is there
        if ";part=" in attributes: ID += "_part%s"%([x.lstrip("part=").split("/")[0] for x in attributes.split(";") if x.startswith("part=")][0])

        return ID

    df_gff3["ID"] = df_gff3.attributes.apply(get_ID_gff)

    # get the bed of the interesting features
    df_bed = df_gff3[["chromosome", "start", "end", "ID"]]
    if file_is_empty(bed) or replace is True: df_bed.to_csv(path_or_buf=bed, sep="\t", index=False, header=False)

    # check that the ID assignation is correct
    if len(set(df_bed.index))!=len(df_bed): raise ValueError("The gene IDs unique IDs in the provided gffs does not match the unique numer of gene IDs. This may be because of unexpected formatting of the 'attributes' field of the provided gff file")

    #### get the bed of the regions sorrounding the genes ####

    # map each chromosome to it's length
    chromosome_to_lenght = {chrRec.id : len(str(chrRec.seq)) for chrRec in SeqIO.parse(reference, "fasta")}

    # define two functions that take the start and the end and return the left and right bounds of the gene +- window
    def get_left_bound_of_gene(start):

        """Takes only the start"""

        left_bound = start - window_l
        if left_bound<1: return 1
        else: return left_bound

    def get_right_bound_of_gene(chr_end_tup):

        """Takes a tuple of chromosome and end"""

        chromosome, end = chr_end_tup
        chromosome_len = chromosome_to_lenght[chromosome]

        right_bound = end + window_l
        if right_bound>chromosome_len: return chromosome_len
        else: return right_bound

    # get the regions
    df_bed = df_bed.rename(index=str, columns={"start":"gene_start", "end":"gene_end", "ID":"gene_ID"}) # rename the column names
    df_bed["region_start"] = df_bed.gene_start.apply(get_left_bound_of_gene)

    df_bed["chromosome_end_tuple"] = [(chrom, end) for chrom, end in df_bed[["chromosome", "gene_end"]].values]
    df_bed["region_end"] = df_bed.chromosome_end_tuple.apply(get_right_bound_of_gene)

    # write always
    regions_filename = "%s.regions_lenght%i"%(bed, window_l)
    df_bed[["chromosome", "region_start", "region_end", "gene_ID"]].to_csv(path_or_buf=regions_filename, sep="\t", index=False, header=False)

    # return for further usage
    return regions_filename


def get_availableGbRAM():

    """This function returns a float with the available memory in your system"""

    # define the memory available on several ways

    # check if meminfo exists
    try: 
        if len(open("/proc/meminfo", "r").readlines())>0:
            meminfo_exists = True
        else: meminfo_exists = False
    except: meminfo_exists = False

    # slurm cluster environment
    if "SLURM_MEM_PER_CPU" in os.environ: 

        # the available memory is the number of CPUs x the number of mem per CPU
        mem_per_cpu = int(os.environ["SLURM_MEM_PER_CPU"])/1000 # the variable is in Mb
        ncpus = int(os.environ["SLURM_CPUS_PER_TASK"])
        available_mem = mem_per_cpu*ncpus

    # the /proc/meminfo file
    elif meminfo_exists:

        lines_availableMem = [float(l.split()[1])/1000000 for l in open("/proc/meminfo", "r").readlines() if l.startswith("MemAvailable:") and l.strip().endswith("kB")]

        if len(lines_availableMem)!=1: raise ValueError("there are more than one correct lines when calculating the memory consumption.")

        available_mem = lines_availableMem[0]

    else: raise ValueError("The available memory cannot be determined")

    return available_mem


def write_coverage_per_gene_mosdepth_and_parallel(sorted_bam, reference_genome, cnv_outdir, bed, gene_to_coverage_file, replace=False):

    """Takes a bam, a bed (1-indexed) and a file were the gene-to-coverage should be written, and writes a file with the gene_to_coverage info """

    if file_is_empty(gene_to_coverage_file) or replace is True:
        print_if_verbose("generating %s"%gene_to_coverage_file)

        # get the first three cols, important for the coverage calc
        cnv_bed = "%s/%s"%(cnv_outdir, get_file(bed))
        cnv_bed_stderr = "%s.generating.stderr"%cnv_bed
        print_if_verbose("getting the first 3 cols of a bed. The stderr can be found in %s"%cnv_bed_stderr)
        run_cmd("cut -f1-3 %s > %s 2>%s"%(bed, cnv_bed, cnv_bed_stderr))
        remove_file(cnv_bed_stderr)

        # get the coverage file for the first three
        coverage_file = generate_coverage_per_window_file_parallel(reference_genome, cnv_outdir, sorted_bam, windows_file=cnv_bed, replace=replace, run_in_parallel=True, delete_bams=True)

        # add the ID
        df_bed = pd.read_csv(bed, sep="\t", header=-1, names=["chromosome", "start", "end", "ID"])
        df_coverage = pd.read_csv(coverage_file, sep="\t")
        df_coverage_with_ID = df_coverage.merge(df_bed, how="right", left_on=["#chrom", "start", "end"], right_on=["chromosome", "start", "end"], right_index=False, left_index=False, validate="one_to_many")[['chromosome', 'start', 'end', 'length', 'mediancov_1', 'nocoveragebp_1', 'percentcovered_1', 'ID']].rename(columns={"mediancov_1":"median_reads_per_gene"})

        # add fields
        df_coverage_with_ID["fraction_covered_by_MoreThan1read"] = df_coverage_with_ID.percentcovered_1 / 100
        df_coverage_with_ID["relative_coverage"] = df_coverage_with_ID.median_reads_per_gene / np.median(df_coverage_with_ID.median_reads_per_gene)

        # replace 
        coverage_file_with_ID = "%s.withID"%coverage_file
        df_coverage_with_ID.to_csv(coverage_file_with_ID, sep="\t", index=False, header=True)

        print_if_verbose("writing %s"%gene_to_coverage_file)
        os.rename(coverage_file_with_ID, gene_to_coverage_file)


def run_gatk_HaplotypeCaller(outdir_gatk, ref, sorted_bam, ploidy, threads, coverage, replace=False):

    """Runs haplotype caller under outdir and returns the filename of the filtered results"""

    # make the outdir if not there
    if not os.path.isdir(outdir_gatk): os.mkdir(outdir_gatk)

    # run GATK
    gatk_out = "%s/output.raw.vcf"%outdir_gatk; gatk_out_tmp = "%s.tmp"%gatk_out
    if file_is_empty(gatk_out) or replace is True:
        print_if_verbose("Running GATK HaplotypeCaller in ploidy mode...")
        
        gatk_std = "%s.running.std"%gatk_out
        print_if_verbose("the log of gatk can be found in %s"%gatk_std)
        gatk_cmd = "%s HaplotypeCaller -R %s -I %s -O %s -ploidy %i --genotyping-mode DISCOVERY --emit-ref-confidence NONE --stand-call-conf 30 --native-pair-hmm-threads %i > %s 2>&1"%(gatk, ref, sorted_bam, gatk_out_tmp, ploidy, threads, gatk_std)
        run_cmd(gatk_cmd)
        os.rename(gatk_out_tmp, gatk_out)
        remove_file(gatk_std)

        # rename the index as well
        os.rename("%s.tmp.idx"%gatk_out, "%s.idx"%gatk_out)

    # variant filtration. There's a field called filter that has the FILTER argument
    gatk_out_filtered = "%s/output.filt.vcf"%outdir_gatk; gatk_out_filtered_tmp = "%s.tmp"%gatk_out_filtered
    if file_is_empty(gatk_out_filtered) or replace is True:
        print_if_verbose("Running GATK HaplotypeCaller Variant filtration...")

        # this depends on the ploidy. If ploidy is 2 you don't want to filter out heterozygous positions
        if ploidy==1: filterHeterozygous = '-G-filter-name "heterozygous" -G-filter "isHet == 1"'
        else: filterHeterozygous = ''

        gatk_filtered_std = "%s.filtering.std"%gatk_out
        print_if_verbose("the log of GATK filtering is in %s"%gatk_filtered_std)
        gatk_filt_cmd = '%s VariantFiltration -V %s -O %s -cluster 5 -window 20 %s --filter-name "BadDepthofQualityFilter" -filter "DP <= %i || QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" > %s 2>&1'%(gatk, gatk_out, gatk_out_filtered_tmp, filterHeterozygous , coverage, gatk_filtered_std)

        run_cmd(gatk_filt_cmd)
        os.rename(gatk_out_filtered_tmp, gatk_out_filtered)
        remove_file(gatk_filtered_std)

        # rename the index as well
        os.rename("%s.tmp.idx"%gatk_out_filtered, "%s.idx"%gatk_out_filtered)

   
    # return the filtered file
    return gatk_out_filtered

def run_bwa_mem(fastq1, fastq2, ref, outdir, bamfile, sorted_bam, index_bam, name_sample, threads=1, replace=False, MarkDuplicates=True):

    """Takes a set of files and runs bwa mem getting sorted_bam and index_bam. skip_MarkingDuplicates will not mark duplicates"""

    if file_is_empty(sorted_bam) or file_is_empty(index_bam) or replace is True:

        #index fasta
        index_files = ["%s.%s"%(ref, x) for x in ["amb", "ann", "bwt", "pac", "sa"]]

        if any([file_is_empty(x) for x in index_files]) or replace is True:

            # create a branch reference, which will have a tag that is unique to this run. This is important since sometimes you run this pipeline in parallel, and this may give errors in fasta indexing.
            branch_ref = "%s.%s.fasta"%(ref, id_generator())
            shutil.copy2(ref, branch_ref)

            # run indexing in the copy
            indexFasta_std = "%s.std.txt"%branch_ref
            print_if_verbose("indexing fasta. The std can be found in %s"%indexFasta_std)
            cmd_indexFasta = "%s index %s > %s 2>&1"%(bwa, branch_ref, indexFasta_std); run_cmd(cmd_indexFasta) # creates a set of indexes of fasta
            index_files_branch = ["%s.%s"%(branch_ref, x) for x in ["amb", "ann", "bwt", "pac", "sa"]]

            # rename each of the indices so that it matches the ref format
            for branchIDX, realIDX in dict(zip(index_files_branch, index_files)).items(): os.rename(branchIDX, realIDX)

            # rm the branch
            os.unlink(branch_ref)
            os.unlink(indexFasta_std)

        #BWA MEM --> get .sam
        samfile = "%s/aligned_reads.sam"%outdir;
        if (file_is_empty(samfile) and file_is_empty(bamfile)) or replace is True:

            # remove previuous generated temporary file
            if os.path.isfile("%s.tmp"%samfile): os.unlink("%s.tmp"%samfile)

            bwa_mem_stderr = "%s.tmp.stderr"%samfile
            print_if_verbose("running bwa mem. The std is in %s"%bwa_mem_stderr)
            cmd_bwa = '%s mem -R "@RG\\tID:%s\\tSM:%s" -t %i %s %s %s > %s.tmp 2>%s'%(bwa, name_sample, name_sample, threads, ref, fastq1, fastq2, samfile, bwa_mem_stderr); run_cmd(cmd_bwa)
            os.rename("%s.tmp"%samfile , samfile)
            remove_file(bwa_mem_stderr)

        # convert to bam 
        if file_is_empty(bamfile) or replace is True:
            bamconversion_stderr = "%s.tmp.stderr"%bamfile
            print_if_verbose("Converting to bam. The std is in %s"%bamconversion_stderr)
            cmd_toBAM = "%s view -Sbu %s > %s.tmp 2>%s"%(samtools, samfile, bamfile, bamconversion_stderr); run_cmd(cmd_toBAM)

            # remove the sam
            os.unlink(samfile)
            remove_file(bamconversion_stderr)
            os.rename("%s.tmp"%bamfile , bamfile)


        # define the sorted_bam with Dups
        sorted_bam_noMarkDups = "%s.noMarkDups"%sorted_bam
        sorted_bam_noMarkDups_tmp = "%s.tmp"%sorted_bam_noMarkDups

        # sorting bam
        if file_is_empty(sorted_bam_noMarkDups) or replace is True:
            print_if_verbose("Sorting bam")

            # remove all temporary files generated previously in samtools sort (they'd make a new sort to be an error)
            for outdir_file in os.listdir(outdir): 
                fullfilepath = "%s/%s"%(outdir, outdir_file)
                if outdir_file.startswith("aligned_reads") and ".tmp." in outdir_file: os.unlink(fullfilepath)

            # sort
            bam_sort_std = "%s.tmp.sortingBam_std.txt"%sorted_bam
            print_if_verbose("the sorting bam std is in %s"%bam_sort_std)
            cmd_sort = "%s sort --threads %i -o %s %s > %s 2>&1"%(samtools, threads, sorted_bam_noMarkDups_tmp, bamfile, bam_sort_std); run_cmd(cmd_sort)

            # rename
            remove_file(bam_sort_std)
            os.rename(sorted_bam_noMarkDups_tmp, sorted_bam_noMarkDups)

        # mark duplicates or not, depending on MarkDuplicates
        if file_is_empty(sorted_bam) or replace is True:

            if MarkDuplicates is True:

                print_if_verbose("marking duplicates")
                # mark duplicates
                sorted_bam_MarkedDuplicates = get_sortedBam_with_duplicatesMarked(sorted_bam_noMarkDups, threads=threads, replace=replace)

                # remove the the raw bam file
                remove_file("%s.bai"%sorted_bam_MarkedDuplicates)

                # replace
                os.rename(sorted_bam_MarkedDuplicates, sorted_bam)

            else: os.rename(sorted_bam_noMarkDups, sorted_bam)

        # remove unnecessary files
        remove_file(bamfile)
        remove_file(sorted_bam_noMarkDups)

    # indexing bam
    if file_is_empty(index_bam) or replace is True:
        bam_index_std = "%s.indexingBam_std.txt"%sorted_bam
        print_if_verbose("indexing bam. The std is in %s"%bam_index_std)
        cmd_indexBam = "%s index -@ %i %s > %s 2>&1"%(samtools, threads, sorted_bam, bam_index_std); run_cmd(cmd_indexBam)   # creates a .bai of sorted_bam
        remove_file(bam_index_std)


def make_flat_listOflists(LoL):

    return list(itertools.chain.from_iterable(LoL))

def clean_reference_genome_windows_files(reference_genome):

    """Cleans all the files under reference_genome that are windows files and bed's """

    print_if_verbose("removing windows files")
    ref_dir = get_dir(reference_genome)
    ref_name = get_file(reference_genome)

    for file in os.listdir(ref_dir):
        if file.startswith(ref_name) and "windows" in file and "bp.bed" in file : remove_file("%s/%s"%(ref_dir, file))

def get_affected_region_bed_for_SVdf(svDF, svtype, interesting_chromosomes, add_interval_bp=1000, first_position_idx=0, translocations_type="breakpoint_pos", chr_to_len={}, insertions_type="only_one_chrB_breakpoint"):

    """This function takes a df with SVs and returns the regions that are affected in bed format, which depends on the svtype. It adds an interval arround this values. It also returns the number 

    only consider svs were ANY of the chroms are in interesting_chromosomes.

    The bed for the translocations includes the left and right regions of all positions.

    translocations_type can be 'start_and_end_pos' or 'whole_arms' or 'whole_chromosomes', which will imply that only the regions arround breakpoints or the whole arms will be returned. chr_to_len is only needed when you have whole_arms_2orientations

    
    insertions_type defines how to define the bed for insertions
    """


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

            if insertions_type=="only_one_chrB_breakpoint":

                # initialize with the origin chromosome
                affected_region_bed_df = svDF.rename(columns={"ChrA":"chromosome", "StartA":"start", "EndA":"end"})[["chromosome", "start", "end"]]

                # add the end only with the start
                chrB_df = svDF.rename(columns={"ChrB":"chromosome", "StartB":"start"}); chrB_df["end"] = chrB_df["start"]
                affected_region_bed_df = affected_region_bed_df.append(chrB_df[["chromosome", "start", "end"]], sort=True)

            elif insertions_type=="start_and_end_chrB":

                # initialize with the origin chromosome
                affected_region_bed_df = svDF.rename(columns={"ChrA":"chromosome", "StartA":"start", "EndA":"end"})[["chromosome", "start", "end"]]

                # add the end only with the start
                chrB_df = svDF.rename(columns={"ChrB":"chromosome", "StartB":"start", "EndB":"end"})
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
    regions_without_SV_bed_stderr = "%s.generating.stderr"%regions_without_SV_bed
    print_if_verbose("generating regions_without_SV_bed. The stderr can be found in %s"%regions_without_SV_bed_stderr)
    run_cmd("%s subtract -a %s -b %s > %s 2>%s"%(bedtools, all_regions_bed, regions_with_SV_bed, regions_without_SV_bed, regions_without_SV_bed_stderr))
    remove_file(regions_without_SV_bed_stderr)

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

        # add the ID
        if "ID" in r.keys(): er["ID"] = r["ID"]

    else: er = {field : r[field] for field in ["StartA", "EndA", "SizeA", "ChrA", "StartB", "EndB", "SizeB", "ChrB", "BpSeqA", "BpSeqB", "ID"] if field in r.keys()}

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

    print_if_verbose("rewriting %s"%translocations_file)

    # keep the unmodified version
    copying_translocations_std = "%s.copying.std"%translocations_file
    print_if_verbose("Copying translocations. The std is in %s"%copying_translocations_std)
    run_cmd("cp %s %s.unmodified > %s 2>&1"%(translocations_file, translocations_file, copying_translocations_std))
    remove_file(copying_translocations_std)

    # define chromosome_to_length
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

    # load translocations
    df = pd.read_csv(translocations_file, sep="\t")

    # define a function that takes a translocation and formats it
    df_corrected = df.apply(lambda r: format_translocation_row_simulateSV(r, chr_to_len), axis=1)

    df_corrected = df_corrected[[c for c in ["Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB", "Balanced", "BpSeqA", "BpSeqB", "ID"] if c in df_corrected.columns]]

    # write the translocations
    df_corrected.to_csv(translocations_file, sep="\t", header=True, index=False)

def change_EmptyString_to_X(string):

    """Takes an empty string and replaces it with X"""

    if pd.isna(string) or len(string)==0 or string=="": return "X"
    else: return string

def rewrite_insertions_uniformizedFormat_simulateSV(insertions_file):

    """Takes an insertions file and rewrites the fact that the BpSeq is empty"""

    print_if_verbose("rewriting %s"%insertions_file)

    # load df
    df = pd.read_csv(insertions_file, sep="\t")

    # correct BpSeq
    for field in ["BpSeqA", "BpSeqB_5prime", "BpSeqB_3prime"]: 
        if field in df.keys(): df[field] = df[field].apply(change_EmptyString_to_X)
    df["Copied"] = df.Copied.apply(lambda x: str(x).upper())

    # keep an uncorrected version
    copying_insertions_std = "%s.copying.std"%insertions_file
    print_if_verbose("copying insertions. The std is in %s"%copying_insertions_std)
    run_cmd("cp %s %s.unmodified > %s 2>&1"%(insertions_file, insertions_file, copying_insertions_std))
    remove_file(copying_insertions_std)

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

    """ This function takes a rearranged genome and reinserts the copy-and-paste insertions where they should be. The copy-and-paste insertions will be set to cut-and-paste if they can't be processed by this pipeline. This means that the insertions_file will be overwritten """

    print_if_verbose("reinserting-copy-and-paste insertions into %s"%insertions_file)

    # load df and keep the copy-and-paste insertions
    df_all = pd.read_csv(insertions_file, sep="\t")

    # define all the initial ids
    initial_IDs = cp.deepcopy(set(df_all.ID))

    # check that the keys are the ones of the simulation by RSVsim
    insertions_fields = ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Copied", "ID"]
    expected_keys = set(insertions_fields)
    if set(df_all.keys())!=expected_keys: raise ValueError("insertons_file %s does not have the expected headers"%insertions_file)

    # check that the IDs are unique
    if len(set(df_all.ID))!=len(df_all): raise ValueError("the IDs should be unique")

    # define the df of copied insertions (which have to be changed)
    df = df_all[df_all.Copied]

    if len(df)>0:

        # define an unmodified genome
        rearranged_genome_unmodified = "%s.unmodified.fasta"%rearranged_genome
        rearranged_genome_unmodified_tmp = "%s.tmp"%rearranged_genome_unmodified

        if file_is_empty(rearranged_genome_unmodified):

            # initialize a df that will contain the final insertions (with the cut-and-paste, that are correct)
            df_final_insertions = df_all[~df_all.Copied] ### NEWLY ADDED LINE

            # if the unmodified tmps is writen, replace the rearranged_genome with it
            if not file_is_empty(rearranged_genome_unmodified_tmp): os.rename(rearranged_genome_unmodified_tmp, rearranged_genome)

            # get the rearranged genome seq
            chr_to_rearrangedSeq = {seq.id: str(seq.seq).upper() for seq in SeqIO.parse(rearranged_genome, "fasta")}
            all_rearranged_chromosomes_together = "".join(chr_to_rearrangedSeq.values())

            # get the seq
            chr_to_refSeq = {seq.id: str(seq.seq).upper()  for seq in SeqIO.parse(reference_genome, "fasta")}

            # define the length of each chrom
            chr_to_lenSeq = {chrom : len(seq) for chrom, seq in chr_to_refSeq.items()}

            # define all the positions with breakpoints
            df_positions = pd.concat([get_breakpoint_positions_df_in_svDF(svDF) for svtype, svDF in svtype_to_svDF.items()])
            chr_to_bpPositions = dict(df_positions.groupby("Chr").apply(lambda df_c: set(df_c["Pos"])))

            # add the ends of the chromosome, and convert to np array
            for chrom, lenSeq in chr_to_lenSeq.items(): 

                if chrom in chr_to_bpPositions.keys(): chr_to_bpPositions[chrom].update({1, lenSeq})
                else: chr_to_bpPositions[chrom] = {1, lenSeq}
                
                chr_to_bpPositions[chrom] = np.array(sorted(chr_to_bpPositions[chrom]))

            # add the closest breakpoint position of ChrA in the reference
            df["closest_5'breakpoint_position"] = df.apply(lambda r: find_nearest(chr_to_bpPositions[r["ChrA"]][chr_to_bpPositions[r["ChrA"]]<(r["StartA"])], r["StartA"]), axis=1).apply(int)

            ######### debugging code #########
            """
            for Irow, r in df.iterrows():
                print(Irow, r["ChrA"], "len", chr_to_lenSeq[r["ChrA"]])
                print(chr_to_bpPositions[r["ChrA"]])

                print(chr_to_bpPositions[r["ChrA"]]>=(r["EndA"]))
                print(chr_to_bpPositions[r["ChrA"]]>=(r["EndA"]))

                print(chr_to_bpPositions[r["ChrA"]][chr_to_bpPositions[r["ChrA"]]>=(r["EndA"])])

                print(find_nearest(chr_to_bpPositions[r["ChrA"]][chr_to_bpPositions[r["ChrA"]]>=(r["EndA"])], r["EndA"]))
            """
            ##################################

            df["closest_3'breakpoint_position"] = df.apply(lambda r: find_nearest(chr_to_bpPositions[r["ChrA"]][chr_to_bpPositions[r["ChrA"]]>=(r["EndA"])], r["EndA"]), axis=1).apply(int)

            ######### debugging code #########
            """
            for Irow, r in df.iterrows():
                
                print(Irow, r["ChrA"], "len", chr_to_lenSeq[r["ChrA"]])
                print(r["StartA"])
                print(r["closest_5'breakpoint_position"])
                print(chr_to_refSeq[r["ChrA"]][r["closest_5'breakpoint_position"]:r["StartA"]-1])
            """
            ##################################

            # get the 5' sequence (from one position after the closest breakpoint to the position before the breakpoint)
            df["5'sequence"] = df.apply(lambda r: chr_to_refSeq[r["ChrA"]][r["closest_5'breakpoint_position"]:r["StartA"]-1], axis=1)

            # get the 3' sequence (from the position after End to the position before the closest breakpoint)
            df["3'sequence"] = df.apply(lambda r: chr_to_refSeq[r["ChrA"]][r["EndA"]:r["closest_3'breakpoint_position"]-1], axis=1)

            # get the deleted sequence (from the start to the end)
            df["deleted_sequence"] = df.apply(lambda r: chr_to_refSeq[r["ChrA"]][r["StartA"]-1:r["EndA"]], axis=1)

            # change the chromosome seq in the sequence 
            for I, row in df.iterrows():
                print_if_verbose("copy-paste-insertion %i/%i.."%(I+1, len(df_all)))

                # define critical vars
                chrA = row["ChrA"]
                seq5 = row["5'sequence"]
                seq3 = row["3'sequence"]
                del_seq = row["deleted_sequence"]

                # all seq
                ref_seq = seq5+del_seq+seq3

                # conformation in the rearranged chromosome
                rearranged_seq = seq5+seq3

                # check that the rearranged seq appears once in the genome and the ref seq in the ref genome. And they do not cross. If it does, just keep the genome as it was, and set this insertion to false
                chrA_refSeq = chr_to_refSeq[chrA]
                if not(chrA_refSeq.count(ref_seq)==1 and all_rearranged_chromosomes_together.count(rearranged_seq)==1 and all_rearranged_chromosomes_together.count(ref_seq)==0): # chrA_refSeq.count(rearranged_seq)==0

                    print_if_verbose("WARNING: insertion %i is not unique. setting as cut-and-paste instead of copy-and-paste"%(I+1))
                    row["Copied"] = False

                else:
                    # change the location
                    # go through each chrom of the rearranged seqs
                    for chrom in chr_to_rearrangedSeq.keys():

                        # get the rearranged sequence
                        seq = cp.deepcopy(chr_to_rearrangedSeq[chrom])

                        # if the rearrangement sequence is in this chromosome, change it
                        if rearranged_seq in seq: 

                            # update the chr_to_rearrangedSeq so that it contains the reference sequence (copied)
                            chr_to_rearrangedSeq[chrom] = seq.replace(rearranged_seq, ref_seq)
                            break

                # keep this row
                row_df = pd.DataFrame({I : row}).transpose()[insertions_fields]
                df_final_insertions = df_final_insertions.append(row_df)

            # get the rearranged genome into the file
            seq_records_list = [SeqRecord(Seq(seq), id=chrom, name=chrom, description=chrom) for chrom, seq in chr_to_rearrangedSeq.items()]

            # write the unmodified one
            copying_rearranged_genome_std = "%s.copying.std"%rearranged_genome
            print_if_verbose("writing rearranged_genome. The std is in %s"%(copying_rearranged_genome_std))
            run_cmd("cp %s %s.tmp > %s 2>&1"%(rearranged_genome, rearranged_genome_unmodified_tmp, copying_rearranged_genome_std))
            os.rename("%s.tmp"%rearranged_genome_unmodified_tmp, rearranged_genome_unmodified_tmp)
            remove_file(copying_rearranged_genome_std)

            # write the modified genome
            SeqIO.write(seq_records_list, rearranged_genome, "fasta")

            # test that the IDs are the same as at the beggining
            if set(df_final_insertions.ID)!=initial_IDs: raise ValueError("The IDs have changed along the rearrangement of the genome")

            # write the insertions file
            insertions_file_tmp = "%s.tmp"%insertions_file
            df_final_insertions.to_csv(insertions_file_tmp, sep="\t", header=True, index=False)
            os.rename(insertions_file_tmp, insertions_file)

            # write the modified genome
            os.rename(rearranged_genome_unmodified_tmp, rearranged_genome_unmodified)

        else: print_if_verbose("the insertions have already been modified")


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

    bed_non_overlapping_stderr = "%s.generating.stderr"%regions_without_tra_bed
    print_if_verbose("Getting regions without translocations. The stderr is in %s"%bed_non_overlapping_stderr)
    run_cmd("%s subtract -a %s -b %s > %s 2>%s"%(bedtools, target_regions_bed, regions_with_tra_bed, regions_without_tra_bed, bed_non_overlapping_stderr))
    remove_file(bed_non_overlapping_stderr)

    return regions_without_tra_bed

def check_consistency_of_svtype_to_svDF(svtype_to_svDF, all_chromosomes, chr_to_len):

    """Checks whether any of the breakpoints overlap with the others and reports those vars that do"""
 
    print_if_verbose("checking consistency of svtype to svDF")
    df_bed_allRegions = pd.DataFrame()

    # go through each df
    for svtype, svDF in svtype_to_svDF.items():

        if type(svDF)==str: svDF = pd.read_csv(svDF, sep="\t")

        svDF = svDF.set_index("ID", drop=False)

        for varID, sv_series in svDF.iterrows():

            # define series as df
            sv_series_df = pd.DataFrame({0 : sv_series}).transpose()

            # get the bed of this var
            sv_bed, nSVs = get_affected_region_bed_for_SVdf(sv_series_df, svtype, all_chromosomes, chr_to_len=chr_to_len, translocations_type="breakpoint_pos", first_position_idx=1, add_interval_bp=0) # the interval should be 0 because we want to keep only true overlaps 
            sv_bed["varID"] = varID
            sv_bed["svtype"] = svtype

            # get if there is any overlap with df_bed_allRegions
            regions_overlapping = df_bed_allRegions.apply(lambda rprevious: any(sv_bed.apply(lambda rnew: get_is_overlapping_query_vs_target_region(rprevious, rnew), axis=1)), axis=1)

            # if there is any region matching with the previous, continue, if not, keep
            if any(regions_overlapping): 

                # raise error if they do not come both from simulation or if any of the vars comes from a translocation
                if not (all(svDF.loc[{varID}, "ID"].apply(lambda x: "_sim_" in x)) and "_sim_" in varID) or "translocation" in varID or any(svDF.loc[{varID}, "ID"].apply(lambda x: "translocation" in x)):

                    print_if_verbose(regions_overlapping)

                    print_if_verbose("%s has these overlapping regions:\n"%varID, df_bed_allRegions[regions_overlapping])
                    print_if_verbose("The actual var is \n", svDF.loc[varID, svtype_to_fieldsDict[svtype]["all_fields"]],"\n")
                    
                    raise ValueError("there are overlapping regions where they should not")

            # add the bed to the regions matching
            df_bed_allRegions = df_bed_allRegions.append(sv_bed)


def get_ChrB_bp_pos_translocations(r, chr_to_len):

    """Takes a row of a translocations df and returns the breakpoint position"""

    if r["StartB"]==1: return r["EndB"]
    elif r["EndB"]==chr_to_len[r["ChrB"]]: return r["StartB"]
    else: raise ValueError("r is not properly formatted")

def get_orientation_translocation(r, chr_to_len):

    """Gets the orientation of the translocation"""

    if r["StartB"]==1: return "5_to_5"
    elif r["EndB"]==chr_to_len[r["ChrB"]]: return "5_to_3"
    else: raise ValueError("r is not properly formatted")

def get_svDF_in_coords_of_rearranged_genome(svDF, reference_genome, rearranged_genome, svtype, svtype_to_svDF):

    """Takes an svDF and returns it with the coordinates matching those of the rearranged genome (by unique sequence). Those events that can't be mapped will be discarded from the returned svDF. These should be 1-based and the chromosomes in one place should match the chromosomes in the other"""

    # if it is empty, just return it as it is
    if len(svDF)==0: return svDF

    # get the rearranged genome seq
    chr_to_rearrangedSeq = {seq.id: str(seq.seq).upper() for seq in SeqIO.parse(rearranged_genome, "fasta")}
    #all_rearranged_chromosomes_together = "".join(chr_to_rearrangedSeq.values())

    # get the seq
    chr_to_refSeq = {seq.id: str(seq.seq).upper() for seq in SeqIO.parse(reference_genome, "fasta")}

    # define the length of each chrom
    chr_to_ref_lenSeq = {chrom : len(seq) for chrom, seq in chr_to_refSeq.items()}
    chr_to_rearranged_lenSeq = {chrom : len(seq) for chrom, seq in chr_to_rearrangedSeq.items()}

    # define all the positions with breakpoints (this also includes the breakpoints of this svDF). These are positions of the refGenome
    df_positions = pd.concat([get_breakpoint_positions_df_in_svDF(df) for svtype, df in svtype_to_svDF.items()])
    chr_to_bpPositions = dict(df_positions.groupby("Chr").apply(lambda df_c: set(df_c["Pos"])))

    # add the ends of the chromosome, and convert to np array
    for chrom, lenSeq in chr_to_ref_lenSeq.items(): 

        if chrom in chr_to_bpPositions.keys(): chr_to_bpPositions[chrom].update({1, lenSeq})
        else: chr_to_bpPositions[chrom] = {1, lenSeq}

        chr_to_bpPositions[chrom] = np.array(sorted(chr_to_bpPositions[chrom]))


    # check that the ID is unique 
    if len(set(svDF.ID))!=len(svDF): raise ValueError("IDs are not unique")

    ##### PIPELINE DIFFERENTIAL FOR EACH SVTYPE #####

    if svtype=="translocations":

        # make sure that the format is correct
        if set(svDF["StartA"])!={1}: raise ValueError("This svDF is not properly formatted")

        # set the index to be the ID
        svDF = svDF.set_index("ID", drop=False)

        # add if it is 5_to_5 or 5_t_3
        svDF["orientation"] = svDF.apply(lambda r: get_orientation_translocation(r, chr_to_ref_lenSeq), axis=1) 

        # add the position of each breakpoint
        svDF["ChrA_bp_pos"] = svDF["EndA"]
        svDF["ChrB_bp_pos"] = svDF.apply(lambda r: get_ChrB_bp_pos_translocations(r, chr_to_ref_lenSeq), axis=1)

        # at the positions of the closest breakpoints and the corresponding sequences arround the breakpoints
        for chrom in ["ChrA", "ChrB"]:

            # define the breakpoint position field and the sequences
            bp_pos_fiel = "%s_bp_pos"%chrom
            seq_5_field = "%s_5seq"%chrom # this will include until the position before the breakpoint
            seq_3_field = "%s_3seq"%chrom # this will start on the position of the breakpoint
            seq_field = "%s_seq"%chrom # the whole sequence

            # add the closest breakpoint position of chrom in the reference
            svDF["%s_closest_5'bp_pos"%chrom] = svDF.apply(lambda r: find_nearest(chr_to_bpPositions[r[chrom]][chr_to_bpPositions[r[chrom]]<(r[bp_pos_fiel])], r[bp_pos_fiel]), axis=1).apply(int)

            svDF["%s_closest_3'bp_pos"%chrom] = svDF.apply(lambda r: find_nearest(chr_to_bpPositions[r[chrom]][chr_to_bpPositions[r[chrom]]>(r[bp_pos_fiel])], r[bp_pos_fiel]), axis=1).apply(int)

            # add the sequences 

            # 5' seq starts at the position after the breakpoint and ends including the breakpoint position (this may give errors)
            svDF[seq_5_field] = svDF.apply(lambda r: chr_to_refSeq[r[chrom]][r["%s_closest_5'bp_pos"%chrom] : int(r[bp_pos_fiel]-1)], axis=1)

            # 3' seq starts right after the breakpoint and spans until the position before the nex breakpoint
            svDF[seq_3_field] = svDF.apply(lambda r: chr_to_refSeq[r[chrom]][int(r[bp_pos_fiel]-1) : (r["%s_closest_3'bp_pos"%chrom]-1)], axis=1)

            # the merged seqs
            svDF[seq_field] = svDF[seq_5_field] + svDF[seq_3_field]


        # initialize the df svDF_rearrangedCoords
        svDF_rearrangedCoords = pd.DataFrame(columns=svtype_to_fieldsDict[svtype]["all_fields"])

        # go through each SVdf and add to svDF_rearrangedCoords if the sequences are unique
        for ID, sv_row  in svDF.iterrows():

            # check that the 3' seq is unique
            if not ( chr_to_rearrangedSeq[sv_row["ChrA"]].count(sv_row["ChrA_3seq"])==1 and chr_to_rearrangedSeq[sv_row["ChrB"]].count(sv_row["ChrB_3seq"])==1 ): 

                print_if_verbose("The 3' seq appears %i times in chrA "%chr_to_rearrangedSeq[sv_row["ChrA"]].count(sv_row["ChrA_3seq"]))
                print_if_verbose("The 3' seq appears %i times in chrB "%chr_to_rearrangedSeq[sv_row["ChrB"]].count(sv_row["ChrB_3seq"]))
                #print_if_verbose("The 3' seq of chrA is %s"%sv_row["ChrA_3seq"])
                #print_if_verbose("The 3' seq of chrB is %s"%sv_row["ChrB_3seq"])
                print_if_verbose("WARNING: The sequences for %s are not unique enough to find the position of the bp in the rearranged genome"%ID)
                continue

            # define general parameters of the rearranged genome
            ChrA = sv_row["ChrA"]
            ChrB = sv_row["ChrB"]
            Balanced = sv_row["Balanced"]
            StartA = 1

            # define the breakpoint positions in 1-based coordinates (the find() returns 0 if not found)
            ChrA_bp_pos = chr_to_rearrangedSeq[ChrA].find(sv_row["ChrA_3seq"]) + 1
            ChrB_bp_pos = chr_to_rearrangedSeq[ChrB].find(sv_row["ChrB_3seq"]) + 1

            if any([x==0 for x in {ChrA_bp_pos, ChrB_bp_pos}]): raise ValueError("The breakpoints can't be at 0")

            # define the other coordinates
            EndA = ChrA_bp_pos

            if sv_row["orientation"]=="5_to_5": 
                StartB = 1
                EndB = ChrB_bp_pos

            elif sv_row["orientation"]=="5_to_3":
                StartB = ChrB_bp_pos
                EndB = chr_to_rearranged_lenSeq[ChrB]

            else: raise ValueError('sv_row["orientation"] is incorrect') 

            # add to df
            dict_var = {"ChrA":ChrA, "StartA":StartA, "EndA":EndA, "ChrB":ChrB, "StartB":StartB, "EndB":EndB, "Balanced":Balanced, "ID":ID}
            svDF_rearrangedCoords = svDF_rearrangedCoords.append(pd.DataFrame({ID: dict_var}).transpose()[svtype_to_fieldsDict[svtype]["all_fields"]])

        print_if_verbose("You have been able to remap the positions for %i/%i translocations"%(len(svDF_rearrangedCoords), len(svDF)))

        if len(svDF)>5 and len(svDF_rearrangedCoords)==0: raise ValueError("The remapping of translocations did not work")

    else: raise ValueError("This has not been developed for %s"%svtype)

    return svDF_rearrangedCoords

def get_genomeGraph_object_5to3_noBreakpoints(genome, genomeGraph_outfileprefix, replace=False, check_genome=False):

    """This function takes a genome and generates a directed graph where each node is a position in the genome and the edges are 5->3 relationships. It is saved under genomeGraph_outfileprefix"""

    # define the files
    genomeGraph_outfile = "%s.graph.py"%genomeGraph_outfileprefix
    genomeGraph_positions_df = "%s.df_positions.py"%genomeGraph_outfileprefix

    if any([file_is_empty(x) for x in {genomeGraph_outfile, genomeGraph_positions_df}]) or replace is True:
        print_if_verbose("getting genome graph")

        # map each chromosome to an offset
        chrom_to_lenSeq = get_chr_to_len(genome)

        # deffine an offset for each chromosome, which is necessary to keep all the positions of the genome as independent numbers
        chrom_to_offset = {}
        current_offset = 0
        for chrom, seqLen in chrom_to_lenSeq.items():
            chrom_to_offset[chrom] = current_offset
            current_offset+=seqLen

        # create the graph
        genome_graph = igraph.Graph(directed=True)

        # add one vertex (node) for each position in the genome
        npositions = sum(chrom_to_lenSeq.values())
        genome_graph.add_vertices(npositions)

        # define the edges that are the end of chromosomes
        chromosome_start_nodes = {offset for chrom, offset in chrom_to_offset.items()}
        chromosome_end_nodes = {(offset + chrom_to_lenSeq[chrom] - 1) for chrom, offset in chrom_to_offset.items()}

        # define the edges mappping each position to the next one, but not the chromosome_end ones
        all_positions = set(range(npositions))
        non_end_positons = all_positions.difference(chromosome_end_nodes)
        all_edges = [(pos, pos+1) for pos in non_end_positons]

        # add the edges to the graph
        genome_graph.add_edges(all_edges)
        print_if_verbose("genome graph got")

        # get the real ends of the chromosomes (regardless of the connected regions)
        if check_genome is True:

            sorted_positions = sorted(all_positions)
            pos_to_nNeighbors = pd.Series(dict(zip(sorted_positions, map(lambda x: len(genome_graph.neighbors(x, mode="ALL")), sorted_positions))))

            real_chromosome_end_nodes = set(pos_to_nNeighbors[pos_to_nNeighbors==1].index)
            print_if_verbose("There are %i telomeric nodes in the graph genome"%len(real_chromosome_end_nodes))

        print_if_verbose("getting positions df")

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
        df_positions["is_start_of_chr"] = df_positions.graph_position.isin(chromosome_start_nodes)
        df_positions["is_end_of_chr"] = df_positions.graph_position.isin(chromosome_end_nodes)

        print_if_verbose("getting nucleotide")

        # map each chromosome to a sequence
        chr_to_seq = {seq.id: str(seq.seq).upper() for seq in SeqIO.parse(genome, "fasta")}

        # map each chrom to each position to a seq
        print_if_verbose("getting positions dict")
        chrom_to_pos_to_seq = {chrom : dict(zip(range(len(seq)) , seq)) for chrom, seq in chr_to_seq.items()}

        chrom_AND_pos_to_seq = {}
        for chrom, pos_to_seq in chrom_to_pos_to_seq.items():
            for pos, seq in pos_to_seq.items(): chrom_AND_pos_to_seq["%s_%i"%(chrom, pos)] = seq

        print_if_verbose("adding to df")
        df_positions["chrom_AND_pos"] = df_positions.chromosome + "_" + df_positions.real_position.apply(str)
        df_positions["nucleotide"] = df_positions.chrom_AND_pos.map(chrom_AND_pos_to_seq)
        nucleotide_to_compNucleotide = {"A":"T", "C":"G", "T":"A", "G":"C", "N":"N", "W":"W", "S":"S", "M":"K", "K":"M", "R":"Y", "Y":"R", "B":"V", "D":"H", "H":"D", "V":"B", "Z":"Z"}
        df_positions["complement_nucleotide"] = df_positions.nucleotide.apply(lambda x: x.upper()).map(nucleotide_to_compNucleotide)        

        if set(df_positions.graph_position)!=all_positions: raise ValueError("There is a bad graph calculation of the positions")

        if any(pd.isna(df_positions.nucleotide)): raise ValueError("There should be no NaNs in the sequence")
        if any(pd.isna(df_positions.complement_nucleotide)): 
            print_if_verbose("These are the nts in the sequence:", set(df_positions.nucleotide))
            print_if_verbose("These are the nts not in the sequence:", set(df_positions.nucleotide).difference(set(nucleotide_to_compNucleotide)))
            raise ValueError("There should be no NaNs in the complementary seq. This suggests that there are strange nulceotides in the genome which have not been considered. perSVade currently considers only %s"%set(nucleotide_to_compNucleotide))

        # save
        save_object(genome_graph, genomeGraph_outfile)
        save_object(df_positions, genomeGraph_positions_df)

    else:
        print_if_verbose("loading graph genome")
        genome_graph = load_object(genomeGraph_outfile)
        df_positions = load_object(genomeGraph_positions_df)

    return genome_graph, df_positions

def get_graph_subcomponents(graph):

    """Takes a graph and returns a list of unique subcomponents (ordered tuples)"""

    components_set = set()

    for node in graph.vs.indices:

        # get the component
        all_graph_positions = tuple(sorted(set(graph.subcomponent(node, mode="ALL"))))

        # add
        components_set.add(all_graph_positions)

    return components_set

def get_clusters_translocated_chromosomes(translocations_df):

    """Takes a translocations df and returns a list of sets, each of which contains a set of chromosomes that are translocated with each other"""

    # get all chromosomes
    all_chromosomes = sorted(set(translocations_df["ChrA"]).union(translocations_df["ChrB"]))
    cID_to_chromosome = dict(zip(range(len(all_chromosomes)), all_chromosomes))
    chrom_to_cID = {chromosome : cID for cID, chromosome in cID_to_chromosome.items()}

    # initialize a graph
    chrom_graph = igraph.Graph(directed=False)

    # add one vertex for each chromosome
    chrom_graph.add_vertices(len(cID_to_chromosome))

    # add connections for each pair of rearranged chromosomes
    edges = [(chrom_to_cID[chrA], chrom_to_cID[chrB]) for chrA, chrB in translocations_df[["ChrA", "ChrB"]].values]
    chrom_graph.add_edges(edges)

    # get all the components of this graph
    chromosomes_set_of_tuples = get_graph_subcomponents(chrom_graph)

    # get as list of sets
    clusters_translocated_chromosomes = [set([cID_to_chromosome[ID] for ID in set_tuples]) for set_tuples in chromosomes_set_of_tuples]

    return clusters_translocated_chromosomes


def map_number_inversions_to_isComplement(ninversions):

    """takes the number of inversions that a position in the genome has and returns whether the complementary nucleotide should be considered"""

    if ninversions%2==0: return False
    else: return True


def write_genome_graph_to_fasta(genome_graph, df_positions, outfile_fasta, translocations_df, df_inverted_positions, replace=False):

    """This function writes a  genome graoh and the associated df_positions into the outfile fasta. It makes some prints to validate that it is correct"""

    print_if_verbose("writing genome graph to %s "%outfile_fasta)

    if file_is_empty(outfile_fasta) or replace is True:

        # set the index of the inverted positions to be the one of the graph_position
        df_inverted_positions = df_inverted_positions.set_index("graph_position", drop=True)

        # get the graph position as the index of the chromosome
        df_positions = df_positions.set_index("graph_position", drop=False)

        # add the nucleotide has to be complementary
        print_if_verbose("adding whether the complementary nt should be obtained")
        df_positions["number_inversions"] = df_inverted_positions.sum(axis=1).loc[df_positions.index]
        df_positions["is_complement_nt"] = df_positions.number_inversions.apply(map_number_inversions_to_isComplement)

        print_if_verbose("These are the number of inversions:", set(df_positions.number_inversions))

        # get the final nucleotide
        print_if_verbose("getting the final nucleotide")
        series_nucleotide = df_positions[~df_positions.is_complement_nt]["nucleotide"]
        series_complement_nucleotide = df_positions[df_positions.is_complement_nt]["complement_nucleotide"]
        df_positions["final_nucleotide"] = series_nucleotide.append(series_complement_nucleotide).loc[df_positions.index]

        print_if_verbose("There are %i positions that should be complemented "%(sum(df_positions["is_complement_nt"])))

        # get the clusters of translocated chromosomes
        print_if_verbose("getting clusters of translocated chroms")
        clusters_translocated_chromosomes = get_clusters_translocated_chromosomes(translocations_df)

        # add content to each of the df_positions
        df_positions["n_succesors"] = df_positions.graph_position.apply(lambda x: genome_graph.successors(x)).apply(len)
        df_positions["n_predecessors"] = df_positions.graph_position.apply(lambda x: genome_graph.predecessors(x)).apply(len)
        df_positions["is_start"] = df_positions.n_predecessors==0
        df_positions["is_end"] = df_positions.n_succesors==0

        # debug
        if any(df_positions["n_succesors"]>1) or any(df_positions["n_predecessors"]>1) or any((df_positions["is_start"]) & df_positions["is_end"]): raise ValueError("Your graph is incorrect")

        # define the end and start positions of the chromosome
        start_positions = set(df_positions[df_positions.is_start]["graph_position"])
        end_positions = set(df_positions[df_positions.is_end]["graph_position"])
        all_positions = set(df_positions["graph_position"])

        # define all the chromosome names, and check that their length is equal to the number of start_pos
        all_chromosomes = set(df_positions.chromosome)
        if len(all_chromosomes)!=len(start_positions) or len(all_chromosomes)!=len(end_positions): raise ValueError("The number of chromosomes and the one of end or start positions does not match")

        # map each end position to a chromosome
        startPos_to_chrom = dict(df_positions[df_positions.is_start]["chromosome"])
        endPos_to_chrom = dict(df_positions[df_positions.is_end]["chromosome"])

        # initialize all the found positions (to check that the genome is a linear graph)
        all_matched_positions = set()
        all_matched_end_positions = set()

        #initialize an iteratior of seq records
        all_chromosomes_SeqRecords = []

        print_if_verbose("There are %i start and %i end positions"%(len(start_positions), len(end_positions)))

        # go through each start position and find the end
        for start_pos in start_positions:
            print_if_verbose("working on start pos %i"%start_pos)

            # get all the positions interconnected
            all_chr_positions = set(genome_graph.subcomponent(start_pos, mode="OUT"))

            # define the end pos
            end_pos = end_positions.intersection(all_chr_positions)
            if len(end_pos)!=1: raise ValueError("There can't be more or less than 1 end position")
            end_pos = next(iter(end_pos))

            # get the paths from start to end
            all_paths = genome_graph.get_shortest_paths(start_pos, end_pos, mode="OUT")
            if len(all_paths)!=1: raise ValueError("There has to be only one path from start to end")
            chr_path = all_paths[0]

            # get the sequence
            chrom_seq = "".join(df_positions.loc[chr_path, "final_nucleotide"])

            # get the ID randomly of the start and end IDs
            start_chrom = startPos_to_chrom[start_pos]
            end_chrom = endPos_to_chrom[end_pos]
            if start_chrom==end_chrom: print_if_verbose("The start and end chromosomes are the same for %s"%start_chrom)

            # take the chromosome from the available chromosomes, if there are none, from the chromosomes that are in there
            chromosomes_this_chr = set(df_positions.loc[chr_path, "chromosome"])
            available_chromosomes = chromosomes_this_chr.intersection(all_chromosomes)
            
            # if there is no such intersection pick any chromosome from the same cluster
            if len(available_chromosomes)==0: 
                print_if_verbose("looking for chromosomes in the similar cluster")

                available_chromosomes = set.union(*[cluster for cluster in clusters_translocated_chromosomes if len(chromosomes_this_chr.intersection(cluster))>0]).intersection(all_chromosomes)

            chrom_ID = next(iter(available_chromosomes))

            # add to the SeqRecord
            all_chromosomes_SeqRecords.append(SeqRecord(Seq(chrom_seq), id=chrom_ID, name=chrom_ID, description=chrom_ID))

            # delete this from all chromosomes
            all_chromosomes = all_chromosomes.difference({chrom_ID})

            # record which positions have been added to a chromosome
            all_matched_positions.update(all_chr_positions)
            all_matched_end_positions.add(end_pos)

        # check that each end and each position is correctly placed
        if all_matched_positions!=all_positions: raise ValueError("Not all positions have been included in the genome")
        if all_matched_end_positions!=end_positions: raise ValueError("Not all end positions have been included")
        if len(all_chromosomes)!=0: raise ValueError("Not all chromosomes were picked")

        print_if_verbose("writing %s"%outfile_fasta)
        SeqIO.write(all_chromosomes_SeqRecords, outfile_fasta, "fasta")

def insert_translocations_into_rearranged_genome(reference_genome, input_rearranged_genome, output_rearranged_genome, svDF, translocations_file, svtype_to_svDF, replace=False):

    """This function takes a rearranged genome and insert the translocations generating output_rearranged_genome. This substitutes the translocations_file in case that some translocations cannot be inserted. svtype_to_svDF should contain translocations. The svDF should have 1s on it. svDF has to be in 1-based coordinates"""

    print_if_verbose("inserting translocations inHouse")

    # test that all are balanced translocations
    if not all(svDF.Balanced): raise ValueError("This has not been developed for unbalanced translocations")

    # get the svDF in coordinates of the rearranged genome
    svDF_rearrangedCoords = get_svDF_in_coords_of_rearranged_genome(svDF, reference_genome, input_rearranged_genome, "translocations", svtype_to_svDF)

    if len(svDF_rearrangedCoords)>0:

        # add fields to the rearrangedCoords df
        chr_to_rearranged_len = get_chr_to_len(input_rearranged_genome)
        svDF_rearrangedCoords["orientation"] = svDF_rearrangedCoords.apply(lambda r: get_orientation_translocation(r, chr_to_rearranged_len), axis=1) 
        svDF_rearrangedCoords["ChrA_bp_pos"] = svDF_rearrangedCoords["EndA"]
        svDF_rearrangedCoords["ChrB_bp_pos"] = svDF_rearrangedCoords.apply(lambda r: get_ChrB_bp_pos_translocations(r, chr_to_rearranged_len), axis=1)

        # rewrite positions so that they are 0-based (so that each position is the real one)
        for pos_f in ["StartA", "EndA", "StartB", "EndB", "ChrA_bp_pos", "ChrB_bp_pos"]: svDF_rearrangedCoords[pos_f] = svDF_rearrangedCoords[pos_f] - 1

    # get the rearranged genome as a directed liner graph
    genomeGraph_outfileprefix = "%s.linearDirectedGraph"%input_rearranged_genome
    genome_graph, df_positions = get_genomeGraph_object_5to3_noBreakpoints(input_rearranged_genome, genomeGraph_outfileprefix, replace=replace) # everything here is 0-based locations

    # define the index of the positons df as the chromosome and real_position
    df_positions = df_positions.set_index(["chromosome", "real_position"], drop=False)

    # define the file that will store the rearranged graph
    genomeGraph_final_file = "%s.linearDirectedGraph.graph.py"%output_rearranged_genome
    feasible_translocation_IDs_file = "%s.feasible_translocation_IDs.py"%translocations_file
    df_inverted_positions_file = "%s.df_inverted_positions.py"%output_rearranged_genome

    if file_is_empty(genomeGraph_final_file) or file_is_empty(feasible_translocation_IDs_file) or file_is_empty(df_inverted_positions_file) or replace is True:

        # initialize the set of feasible IDs
        feasible_translocation_IDs = set()

        # initialize a df that has as index all the positions in the genome
        df_inverted_positions = pd.DataFrame({"graph_position":list(df_positions.graph_position)})

        # inialize the positions that have been already inverted
        already_inverted_positions = set()

        # go through each translocation and modify the corresponding edges in the genome graph
        for ID, r in svDF_rearrangedCoords.iterrows():
            print_if_verbose("generating %s in the graph"%ID)

            # define coords in the location of the graph
            chrA_bp = df_positions.loc[(r["ChrA"], r["ChrA_bp_pos"]), "graph_position"]
            chrB_bp = df_positions.loc[(r["ChrB"], r["ChrB_bp_pos"]), "graph_position"]

            # define the next positions and before ones
            chrA_bp_3pos = genome_graph.successors(chrA_bp)[0]
            chrA_bp_5pos = genome_graph.predecessors(chrA_bp)[0]

            chrB_bp_3pos = genome_graph.successors(chrB_bp)[0]
            chrB_bp_5pos = genome_graph.predecessors(chrB_bp)[0]

            # check that the chrA and chrB are different, and skip if it is not the case. Sometimes successive balanced translocations make others being in the same chromosome
            all_chrA_positions = cp.deepcopy(set(genome_graph.subcomponent(chrA_bp, mode="ALL")))
            all_chrB_positions = cp.deepcopy(set(genome_graph.subcomponent(chrB_bp, mode="ALL")))
            all_chrAandB_positions = all_chrA_positions.union(all_chrB_positions)

            if len(all_chrA_positions.intersection(all_chrB_positions))>0: 
                print_if_verbose("%s is not between different chromosomes (maybe because the multiple rearrangements performed), skipping..."%ID)
                continue

            # if any of the chrB or chrB positions was already inverted, skip this translocation
            if len(already_inverted_positions.intersection(all_chrAandB_positions))>0: 
                print_if_verbose("The chromosome A or B were already inverted once, skipping...")
                continue

            ##### 5_to_5 positions are straightforward because they don't alter the order of anything #####
            if r["orientation"]=="5_to_5":

                # delete the WT unions in chromosome A and chromosome B
                genome_graph.delete_edges([(chrA_bp, chrA_bp_3pos), (chrB_bp, chrB_bp_3pos)])

                # add the union between chrA and the position after chromosome B
                genome_graph.add_edges([(chrA_bp, chrB_bp_3pos)])

                # add the union between chromosome B and the next position of chrA
                genome_graph.add_edges([(chrB_bp, chrA_bp_3pos)])

            ###############################################################################################

            ####### 5_to_3 orientations require changing the orientation of many nodes #######
            elif r["orientation"]=="5_to_3":
                print_if_verbose("inverted translocation")


                # chrA 5' is united with chromB 5' and chrB 3' is united with chrB 3'. We will change the orientation of the connections in chromosome B so that they follow this

                # get the connection bewteen each chrB node and its positions
                pos_to_pos5 = cp.deepcopy({pos : genome_graph.predecessors(pos)[0] for pos in all_chrB_positions if len(genome_graph.predecessors(pos))>0})
                pos_to_pos3 = cp.deepcopy({pos5 : pos for pos, pos5 in pos_to_pos5.items()})


                # keep the positions that will be inverted (those of chromosome)              
                already_inverted_positions.update(all_chrB_positions)

                # add the chrB positions to the inverted genome
                df_inverted_positions[ID] = df_inverted_positions.graph_position.isin(all_chrB_positions)

                # check that the length
                if not len(all_chrB_positions)==(len(pos_to_pos5)+1): raise ValueError("something went wrong with the graph")

                #### delete edges ####

                # delete all the connections in chromosome B
                genome_graph.delete_edges([(pos, pos3) for pos, pos3 in pos_to_pos3.items()])

                # delete the edge in chrA
                genome_graph.delete_edges([(chrA_bp, chrA_bp_3pos)]) # this does not work sometimes

                #######################

                #### add simple edges ####

                # add the edge between the ChrA pos and the ChrB pos
                genome_graph.add_edges([(chrA_bp, chrB_bp)])

                # add the edge between ChrB+1 and ChrA+1
                genome_graph.add_edges([(chrB_bp_3pos, chrA_bp_3pos)])

                ###########################

                # add the 3'->5' connections in chrB, unless it is chrB_bp_3pos
                genome_graph.add_edges([(pos, pos5) for pos, pos5 in pos_to_pos5.items() if pos!=chrB_bp_3pos])

            ##################################################################################

            else: raise ValueError("orientation is not correct")


            # keep the translocation for further printing
            feasible_translocation_IDs.add(ID)

        # debug the fact that there are no translocations
        if len(svDF)>0 and len(feasible_translocation_IDs)==0: raise ValueError("There are no translocations generated, there must be an error.")

        # save
        save_object(genome_graph, genomeGraph_final_file)
        save_object(feasible_translocation_IDs, feasible_translocation_IDs_file)
        save_object(df_inverted_positions, df_inverted_positions_file)

    else:  

        print_if_verbose("loading graph and feasible translocation IDs and inverted positions df")
        genome_graph = load_object(genomeGraph_final_file)
        feasible_translocation_IDs = load_object(feasible_translocation_IDs_file)
        df_inverted_positions = load_object(df_inverted_positions_file)

    # write fasta for the rearranged genome
    write_genome_graph_to_fasta(genome_graph, df_positions, output_rearranged_genome, svDF, df_inverted_positions, replace=replace)

    #### replace translocations file ####

    # delete the translocations file because it has to be rewritten as the end point of this function
    remove_file(translocations_file)
 
    # write into translocations file (this is important so that in case any of the translocations could not be mapped from the sequence)
    svDF_final = svDF[svDF.ID.isin(feasible_translocation_IDs)][svtype_to_fieldsDict["translocations"]["all_fields"]]
    svDF_final.to_csv(translocations_file, sep="\t", header=True, index=False)

    #####################################

    print_if_verbose("There are %i/%i translocations that are feasible"%(len(svDF_final), len(svDF)))

   

def generate_rearranged_genome_from_svtype_to_svDF(reference_genome, svtype_to_svDF, outdir, replace=False):

    """This function generates a rearranged genome for the provided svDFs, writing their fileds under outdir"""    

    # define the rearranged genome, and generated if not already done
    rearranged_genome = "%s/rearranged_genome.fasta"%outdir
    rearranged_genome_InsInvDelTan = "%s/rearranged_genome_InsInvDelTan.fasta"%outdir
    rearranged_genome_InsInvDelTan_tmp = "%s.tmp.fasta"%rearranged_genome_InsInvDelTan
    rearranged_genome_finalFile = "%s.performed"%(rearranged_genome)

    if file_is_empty(rearranged_genome_finalFile) or replace is True:

        # generate the rearranged genome with INS,DEL,INV,TAN
        if file_is_empty(rearranged_genome_InsInvDelTan) or replace is True:

            # remove the outdir
            delete_folder(outdir); make_folder(outdir)

            # initialize a cmd to create the simulated genome
            targetSV_cmd = "%s --input_genome %s --output_genome %s"%(create_targeted_simulatedSVgenome_R, reference_genome, rearranged_genome_InsInvDelTan_tmp)

            for svtype, svDF in svtype_to_svDF.items():

                # shift the insertions by 15 bp so that they are not at the beginning of the chrom
                if svtype=="insertions": svDF["StartA"] = svDF["StartA"] + 10

                # keep only the interesting svs
                svDF = svDF[svtype_to_fieldsDict[svtype]["all_fields"]]

                # write file
                svfile = "%s/%s.tab"%(outdir, svtype)
                svDF.to_csv(svfile, sep="\t", header=True, index=False)

                # skip translocations from simulation
                if svtype=="translocations": continue

                # add the generation of SVs into targetSV_cmd
                targetSV_cmd += " --%s_file %s"%(svtype, svfile)

            # run the cmd
            std_rearranging_genome = "%s/simulation_std.txt"%outdir
            #std_rearranging_genome = "stdout" # debug
            print_if_verbose("rearranging genome. The std is in %s"%std_rearranging_genome)

            if std_rearranging_genome!="stdout": run_cmd("%s > %s 2>&1"%(targetSV_cmd, std_rearranging_genome))
            else: run_cmd(targetSV_cmd)

            # transform the cut-and-paste insertions to copy-and-paste, whenever necessary
            insertions_file = "%s/insertions.tab"%outdir
            if not file_is_empty(insertions_file):

                # transform the insertions
                transform_cut_and_paste_to_copy_and_paste_insertions(reference_genome, rearranged_genome_InsInvDelTan_tmp, insertions_file, svtype_to_svDF)

                # edit the insertions so that they are in the correct format
                rewrite_insertions_uniformizedFormat_simulateSV(insertions_file)


            remove_file(std_rearranging_genome)

            # rename the genome 
            os.rename(rearranged_genome_InsInvDelTan_tmp, rearranged_genome_InsInvDelTan)

        
        # generate translocations
        translocations_file = "%s/translocations.tab"%outdir
        if "translocations" in svtype_to_svDF: insert_translocations_into_rearranged_genome(reference_genome, rearranged_genome_InsInvDelTan, rearranged_genome, svtype_to_svDF["translocations"], translocations_file, svtype_to_svDF)

        # rewrite the variants so that they are optimal for comparison. This is important to re-sort the chromosomes if necessary
        print_if_verbose("rewriting %s"%translocations_file)
        if not file_is_empty(translocations_file): rewrite_translocations_uniformizedFormat_simulateSV(translocations_file, reference_genome)

        # write a file that indicates that this has finsihed
        open(rearranged_genome_finalFile, "w").write("finsihed")



def get_translocations_randomly_placed_in_target_regions(target_regions_bed, translocations_file, chr_to_len, nvars=100, only_5_to_3=False):

    """Writes nvars randomly placed translocations in target_regions_bed, and writes them to translocations_file. It will draw as maximum number of translocations as possbile. Half of them will be inverted and half in the same orientation. All of them are balanced."""

    print_if_verbose("getting randomly inserted translocations")

    # get the bed into a df
    target_regions_df = pd.read_csv(target_regions_bed, sep="\t", names=["chromosome", "start", "end"], header=None).drop_duplicates()
    target_regions_df.index = list(range(len(target_regions_df)))

    # add the breakpoint region in the middle
    target_regions_df["bp_pos"] = (target_regions_df.start + (target_regions_df.end-target_regions_df.start)/2).apply(int)

    # initialize a dict that will be used for the tra df
    varID_to_colName_to_value = {}

    # keep simulating translocations until you have nvars
    nvars_simulated = 0

    while nvars_simulated<nvars:

        # if the target df has less than 2 vars with a different chromosome, drop
        if len(set(target_regions_df.chromosome))<2: break

        # pick a regionA
        all_regions = list(target_regions_df.index)
        regionA = target_regions_df.loc[random.choice(all_regions)]

        # get the regions from different chromosomes
        target_regions_df_B = target_regions_df[target_regions_df.chromosome!=regionA["chromosome"]]
        all_regions_B = list(target_regions_df_B.index)
        regionB = target_regions_df_B.loc[random.choice(all_regions_B)]

        # get the dict for this tra
        tra_dict = {"ChrA":regionA["chromosome"], "StartA":1, "EndA":regionA["bp_pos"], "ChrB":regionB["chromosome"], "Balanced":True}

        # define if inverted or not, which defines the orientation of chrB
        if only_5_to_3 is True: is_inverted = False
        else: 
            # is_inverted = bool(random.randrange(0, 2)) 50% each
            is_inverted = random.uniform(0, 1)>=0.8 # 80% inverted translocations

        if is_inverted is True: 
            tra_dict["StartB"] = regionB["bp_pos"]
            tra_dict["EndB"] = chr_to_len[regionB["chromosome"]]

        else:
            tra_dict["StartB"] = 1
            tra_dict["EndB"] = regionB["bp_pos"]

        # add the 
        varID_to_colName_to_value[nvars_simulated] = tra_dict

        # delete both regions from the possibilities
        target_regions_df = target_regions_df.drop(regionA.name)
        target_regions_df = target_regions_df.drop(regionB.name)

        # update the number of simulated
        nvars_simulated+=1

    # get df
    tra_df = pd.DataFrame(varID_to_colName_to_value).transpose()[["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Balanced"]]
    tra_df["Name"] = ["translocation_%i"%(I+1) for I in range(len(tra_df))]
    tra_df["ID"] = tra_df["Name"]

    print_if_verbose("writing %s"%translocations_file)
    tra_df.to_csv(translocations_file, sep="\t")

def get_random_svtype_to_svDF(reference_genome, mitochondrial_chromosome, outdir, nvars=200, replace=False, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}, check_random_genome_generation=False, only_5_to_3_translocations=False):

    """This function generates nvars into the reference genome splitting by gDNA and mtDNA with files written under outdir. It returns the randomly drawn variants and no-genome"""


    print_if_verbose("generating random simulations")

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
        print_if_verbose(type_genome)

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

        #### GET THE RANDOM INS,INV,DEL ####

        if any([file_is_empty("%s/%s.tab"%(random_sim_dir, svtype)) for svtype in {"insertions", "deletions", "inversions", "tandemDuplications"}.intersection(svtypes)]) or replace is True:

            print_if_verbose("generating random SVs")

            # make and delete the folder
            delete_folder(random_sim_dir); make_folder(random_sim_dir)

            # get the cmd of the simulation
            randomSV_cmd = "%s --input_genome %s --outdir %s --regions_bed %s"%(create_random_simulatedSVgenome_R, genome_file, random_sim_dir, all_regions_bed)

            # add the number of each SV that should be added
            svtype_to_arg = {"insertions":"number_Ins", "deletions":"number_Del", "inversions":"number_Inv", "tandemDuplications":"number_Dup"}
        
            for svtype, arg in svtype_to_arg.items(): 
                if svtype not in svtype_to_arg or svtype not in svtypes: continue

                randomSV_cmd += " --%s %i"%(arg, vars_to_simulate)

            # run the random simulation
            std_rearranging_genome = "%s/simulation_std.txt"%random_sim_dir
            #std_rearranging_genome = "stdout"
            print_if_verbose("getting random SVs. The std is in %s"%std_rearranging_genome)
            if std_rearranging_genome!="stdout": run_cmd("%s > %s 2>&1"%(randomSV_cmd, std_rearranging_genome))
            else: run_cmd(randomSV_cmd)
            remove_file(std_rearranging_genome)

            # edit the insertions 
            insertions_file = "%s/insertions.tab"%random_sim_dir
            rewrite_insertions_uniformizedFormat_simulateSV(insertions_file)

        ########################################

        # define the translocations file
        translocations_file = "%s/translocations.tab"%random_sim_dir

        if file_is_empty(translocations_file) or replace is True:

            ##### CREATE RANDOMLY PLACED TRANSLOCATIONS #####
            if len(chroms)>1 and "translocations" in svtypes: 
                print_if_verbose("generating non-overlapping translocations")

                # get a bed with the previous variants' locations
                InsInvDelTan_bed_df = pd.concat([get_affected_region_bed_for_SVdf("%s/%s.tab"%(random_sim_dir, svtype), svtype, chroms, first_position_idx=1, translocations_type="breakpoint_pos")[0] for svtype in {"insertions", "deletions", "inversions", "tandemDuplications"}.intersection(svtypes)]).sort_values(by=["chromosome", "start", "end"])

                InsInvDelTan_bed = "%s/InsInvDelTan_regions.bed"%genome_outdir
                InsInvDelTan_bed_df.to_csv(InsInvDelTan_bed, sep="\t", header=False, index=False)

                # get a bed file where to place the randomly chosen 
                noInsInvDelTan_bed = "%s/noInsInvDelTan_regions.bed"%genome_outdir
                noInsInvDelTan_bed_stderr = "%s.generating.stderr"%noInsInvDelTan_bed
                print_if_verbose("the stderr is in %s "%noInsInvDelTan_bed_stderr)
                run_cmd("%s subtract -a %s -b %s > %s 2>%s"%(bedtools, all_regions_bed, InsInvDelTan_bed, noInsInvDelTan_bed, noInsInvDelTan_bed_stderr))
                remove_file(noInsInvDelTan_bed_stderr)

                # get the translocations randomly placed
                get_translocations_randomly_placed_in_target_regions(noInsInvDelTan_bed, translocations_file, chrom_to_len, nvars=nvars*2, only_5_to_3=only_5_to_3_translocations)
            #################################################

            # write an empty translocations file
            else: open(translocations_file, "w").write("\t".join(["Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB", "Balanced", "BpSeqA", "BpSeqB"])) # this needs to be

        # add the simulations into random_svtype_to_svDF
        for svtype in random_svtype_to_svDF.keys():
            svDF = random_svtype_to_svDF[svtype]

            # get the new sv
            new_svDF = pd.read_csv("%s/%s.tab"%(random_sim_dir, svtype), sep="\t")
            new_svDF = new_svDF[[c for c in new_svDF.keys() if "BpSeq" not in c]]

            # add the name
            new_svDF["ID"] = new_svDF.Name + "_sim_%s"%type_genome

            # define the merged df


            # append 
            random_svtype_to_svDF[svtype] = svDF.append(new_svDF, sort=True)


    ####### test that you can insert the randomly simulated variants into the genome #######

    if check_random_genome_generation is True:
        print_if_verbose("checking random genome generation")

        # check the consistency of the generated vars. This takes a lot of time
        #check_consistency_of_svtype_to_svDF(random_svtype_to_svDF, set(chrom_to_len), chrom_to_len)

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

def rearrange_genomes_simulateSV(reference_genome, outdir, replace=False, nvars=50, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulated_svtype_to_svfile={}, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}, check_consistency=False):

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
    if check_consistency is True: check_consistency_of_svtype_to_svDF(real_svtype_to_svDF, set(chr_to_len), chr_to_len)

    # define the final outdirs
    final_simulated_SVs_dir = "%s/final_simulated_SVs"%(outdir); 
    final_rearranged_genome = "%s/rearranged_genome.fasta"%final_simulated_SVs_dir
    final_rearranged_genome_finalFile = "%s.performed"%(final_rearranged_genome)


    if file_is_empty(final_rearranged_genome_finalFile) or replace is True:
    #if True:

        #delete_folder(final_simulated_SVs_dir); make_folder(final_simulated_SVs_dir)

        # get random simulations for number of variants (twice the nvars if it is for future merge with simulated_svtype_to_svfile)
        if len(simulated_svtype_to_svfile)==0: random_nvars = nvars
        else: random_nvars = 2*nvars
        random_svtype_to_svDF = get_random_svtype_to_svDF(reference_genome, mitochondrial_chromosome, outdir, nvars=random_nvars, replace=replace, svtypes=svtypes)

        # get the real vars so that no var covers the whole chromosome if it is a tandem duplication, as this is circularization and you we don't want to simulate it
        for svtype, svDF in real_svtype_to_svDF.items():

            if svtype=="tandemDuplications": svDF = svDF[~svDF.apply(lambda r: r["Start"]<(chr_to_len[r["Chr"]]*0.02) and r["End"]>(chr_to_len[r["Chr"]]*0.98), axis=1)]

            real_svtype_to_svDF[svtype] = svDF

        # get a df with all the bed regions of the real vars
        list_affected_bed_regions = [get_affected_region_bed_for_SVdf(svDF, svtype, set(chr_to_len), first_position_idx=1, translocations_type="breakpoint_pos", chr_to_len=chr_to_len)[0] for svtype, svDF in real_svtype_to_svDF.items()]
        if len(list_affected_bed_regions)>0: all_real_SVs_bed_df = pd.concat(list_affected_bed_regions)
        else: all_real_SVs_bed_df = pd.DataFrame()

        # initialize the final_svtype_to_svDF, which will contain vars for both real and simulated regions 
        final_svtype_to_svDF = {}

        # add random vars to final_svtype_to_svDF if they do not overlap with the bed regions
        for svtype in svtypes:
            print_if_verbose("generating %s"%svtype)

            # define the real svtype
            if svtype in real_svtype_to_svDF: real_svDF = real_svtype_to_svDF[svtype]
            else: real_svDF = pd.DataFrame(columns=svtype_to_fieldsDict[svtype]["all_fields"])

            # define the random svtypes
            random_svDF = random_svtype_to_svDF[svtype]

            if len(random_svDF)>0:

                # add to the random SVtype if it overlaps any real svtypes
                random_svDF["overlaps_realSV"] = random_svDF.apply(lambda r: target_svDFseries_overlaps_bed_df(r, all_real_SVs_bed_df, svtype, 1, "breakpoint_pos", chr_to_len), axis=1)

                random_svDF = random_svDF[~random_svDF.overlaps_realSV]

            # for translocations, get the nvars to be twice. This is because some of them may not be feasible
            if svtype=="translocations": real_nvars = nvars*2
            else: real_nvars = nvars

            # get the concatenated df
            svDF = real_svDF.append(random_svDF, sort=True).iloc[0:real_nvars]
   
            # add all the mitochondrial variants (they may be excluded with real_nvars)
            mitochondrial_chromosomes_set = set(mitochondrial_chromosome.split(","))
            random_svDF_mito = random_svDF[random_svDF.apply(lambda r: any([r[f] in mitochondrial_chromosomes_set for f in svtype_to_fieldsDict[svtype]["chromosome_fields"]]), axis=1)].iloc[0:real_nvars]

            svDF = svDF.append(random_svDF_mito, sort=True).drop_duplicates(subset=svtype_to_fieldsDict[svtype]["all_fields"])

            # add to the final set
            final_svtype_to_svDF[svtype] = svDF

        # check the consistency of the genome
        if check_consistency is True: check_consistency_of_svtype_to_svDF(final_svtype_to_svDF, set(chr_to_len), chr_to_len)

        # get the rearranged genome and simulations
        print_if_verbose("rearranging genome with real + random SVs")
        make_folder(final_simulated_SVs_dir)
        generate_rearranged_genome_from_svtype_to_svDF(reference_genome, final_svtype_to_svDF, final_simulated_SVs_dir, replace=replace)

    # define the map between each svtype and the file that defines it
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
        print_if_verbose("generating %s"%fastqgz)

        # convert
        reformatting_std = "%s.generating.std"%fastqgz_tmp
        print_if_verbose("running reformat. The std is in %s"%reformatting_std)
        run_cmd("%s in=%s out=%s qfake=50 overwrite=true > %s 2>&1"%(bbmap_reformat_sh, fasta_file, fastqgz_tmp, reformatting_std))

        # remove the fasta
        if remove_fasta is True: os.unlink(fasta_file)

        remove_file(reformatting_std)
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
        
        # change the end, setting it to -1, and also sorting
        windows_1_based = "%s.1_based.bed"%windows_bed
        windows_1_based_stderr = "%s.generating.stderr"%windows_1_based
        print_if_verbose("getting 1-based bed file. The stderr is in %s"%windows_1_based_stderr)
        run_cmd(""" awk '{print $1 "\t" ($2+1) "\t" ($3)}' %s | sort -k1,1 -k2,2n > %s 2>%s"""%(windows_bed, windows_1_based, windows_1_based_stderr))
        remove_file(windows_1_based_stderr)

        # get the cmd
        mosdepth_std = "%s.generating.std"%fileprefix_tmp
        print_if_verbose("running mosdepth. The std is in %s "%mosdepth_std)
        cmd = "%s --threads %i --by %s --no-per-base --fast-mode --thresholds 1 --use-median %s %s > %s 2>&1"%(mosdepth, extra_threads, windows_1_based, fileprefix_tmp, sorted_bam, mosdepth_std) # mosdepth does not look at internal cigar operations or correct mate overlaps (recommended for most use-cases). It is also faster

        # add the chromosome_id if provided
        if chromosome_id!="": cmd = cmd.replace("--use-median", "--use-median --chrom %s"%chromosome_id)

        # run 
        run_cmd(cmd)
        remove_file(mosdepth_std)

        # remove the 1-based file
        remove_file(windows_1_based)
 
        # remove innecessary files
        for sufix in ["mosdepth.global.dist.txt", "mosdepth.region.dist.txt", "mosdepth.summary.txt", "thresholds.bed.gz.csi", "regions.bed.gz.csi"]: remove_file("%s.%s"%(fileprefix_tmp, sufix))

        # keep
        os.rename(regions_file_tmp, regions_file)
        os.rename(thresholds_file_tmp, thresholds_file)

    # get as dfs
    df_regions = pd.read_csv(regions_file, sep="\t", header=-1, names=["#chrom",  "start", "end", "mediancov_1"]).drop_duplicates(subset=["#chrom",  "start", "end"])
    df_thresholds = pd.read_csv(thresholds_file, sep="\t").drop_duplicates(subset=["#chrom",  "start", "end"])

    # add the number of basepairs in each region that are covered by at least one
    try: df = df_regions.merge(df_thresholds, on=["#chrom",  "start", "end"], validate="one_to_one").rename(columns={"1X":"nbp_more_than_1x"})
    except:
        print_if_verbose("!!!!!!!!regions:%s, \n thresholds:%s, \n prefix:%s"%(regions_file, thresholds_file, fileprefix))
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
        
    # generate a randomID
    randID = id_generator(25)

    # define a file for the coverage
    windows_bed_chromsome = "%s.%s.%s.bed"%(windows_bed, chromosome_id, randID)
    egrepping_windows_stderr = "%s.generaing.stderr"%windows_bed_chromsome
    print_if_verbose("Running egrep on chromosome. The stderr is in %s"%egrepping_windows_stderr)
    run_cmd("egrep '%s\t' %s > %s 2>%s"%(chromosome_id, windows_bed, windows_bed_chromsome, egrepping_windows_stderr))
    remove_file(egrepping_windows_stderr)

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


def generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file="none", replace=False, run_in_parallel=True, delete_bams=True):

    """Takes a reference genome and a sorted bam and runs a calculation of coverage per window (with bamstats04_jar) in parallel for sorted_bam, writing results under ddestination_dir. if window_file is provided then it is used. If not, it generates a file with non overlappping windows of length window_l"""

    # in the case that you have provided a window file
    if windows_file=="none":

        make_folder(destination_dir)

        # first generate the windows file
        windows_file = "%s.windows%ibp.bed"%(reference_genome, window_l)
        windows_file_stderr = "%s.generating.stderr"%windows_file
        print_if_verbose("genearting windows_file. The stderr is in %s"%windows_file_stderr)
        run_cmd("%s makewindows -g %s.fai -w %i > %s 2>%s"%(bedtools, reference_genome, window_l, windows_file, windows_file_stderr))
        remove_file(windows_file_stderr)

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
        coverage_file_tmp = "%s.tmp"%coverage_file
        all_df.to_csv(coverage_file_tmp, sep="\t", header=True, index=False)

        # at the end remove all the bam files # at some point I commented the lines below and I don't know why
        if delete_bams is True and run_in_parallel is True:
            print_if_verbose("removing chromosomal bamfiles")

            for chrom in all_chromosome_IDs: 
                sorted_bam_chr = "%s.%s.bam"%(sorted_bam, chrom)
                remove_file(sorted_bam_chr); remove_file("%s.bai"%sorted_bam_chr)

        # rename
        os.rename(coverage_file_tmp, coverage_file)

    return coverage_file

def plot_coverage_across_genome_pairedEndReads(sorted_bam, reference_genome, replace=False):

    """Takes a sorted_bam and plots the coverage for windows of the genome"""

    print_if_verbose("plotting coverage across genome")

    # get coverage df  
    calculate_coverage_dir = "%s.calculating_windowcoverage"%sorted_bam; make_folder(calculate_coverage_dir)
    coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, calculate_coverage_dir, sorted_bam, replace=replace), "\t")
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

def soft_link_files(origin, target):

    """This function takes an origin file and makes it accessible through a link (target)"""

    if file_is_empty(target):

        # rename as full paths
        origin = get_fullpath(origin)
        target = get_fullpath(target)

        # remove previous lisqnk
        try: run_cmd("rm %s > /dev/null 2>&1"%target)
        except: pass

        soft_linking_std = "%s.softlinking.std"%(target)
        print_if_verbose("softlinking. The std is in %s"%soft_linking_std)
        run_cmd("ln -s %s %s > %s 2>&1"%(origin, target, soft_linking_std))
        remove_file(soft_linking_std)

def run_gridss_and_annotateSimpleType(sorted_bam, reference_genome, outdir, replace=False, threads=4, blacklisted_regions="", maxcoverage=50000):

    """Runs gridss for the sorted_bam in outdir, returning the output vcf. blacklisted_regions is a bed with regions to blacklist"""

    # define the output
    gridss_VCFoutput = "%s/gridss_output.vcf"%outdir
    if gridss_VCFoutput.split(".")[-1]!="vcf": raise ValueError("gridss needs a .vcf file. this is not the case for"%gridss_VCFoutput)

    # there may be a previous gridss vcf, called 'gridss_output.raw.vcf' that is equally valid. Rename it here
    previous_gridss_VCFoutput = "%s/gridss_output.raw.vcf"%outdir
    if not file_is_empty(previous_gridss_VCFoutput) and file_is_empty(gridss_VCFoutput) and replace is False: 
        os.rename(previous_gridss_VCFoutput, gridss_VCFoutput)

    # softlink the sorted_bam under outdir so that the naming is not too long
    sorted_bam_renamed = "%s/aligned_reads.sorted.bam"%outdir
    index_bam_renamed = "%s/aligned_reads.sorted.bam.bai"%outdir
    soft_link_files(sorted_bam, sorted_bam_renamed)
    soft_link_files("%s.bai"%sorted_bam, index_bam_renamed)

    if file_is_empty(gridss_VCFoutput) or replace is True:

        # define other files
        gridss_assemblyBAM = "%s/gridss_assembly.bam"%outdir
        gridss_tmpdir = "%s/gridss_tmp"%outdir

        # if the blacklisted_regions does not exist, just create an empty file
        if file_is_empty(blacklisted_regions): 

            blacklisted_regions = "%s/empty_regions.bed"%outdir; open(blacklisted_regions, "w").write("")

        print_if_verbose("blacklisting %s\n"%blacklisted_regions)
        
        # define the out and error of gridss
        gridss_std = "%s/gridss_run_std.txt"%outdir
        
        max_tries = 2
        for Itry in range(max_tries):
            print_if_verbose("running gridss try %i. The std can be found in %s"%(Itry+1, gridss_std))
            try: 
                # delete previous files
                delete_folder(gridss_tmpdir); make_folder(gridss_tmpdir)
                remove_file(gridss_assemblyBAM)

                # define the ram available
                allocated_ram = get_availableGbRAM()*fractionRAM_to_dedicate
                print_if_verbose("running gridss with %iGb of RAM"%allocated_ram)

                # define the heap size, which depends on the cloud or not
                #jvmheap = "27.5g" # this is the default
                #jvmheap = "20g" # this works in MN. This can be changed to fit the machine
                jvmheap = "%ig"%min([31, int(allocated_ram)]) # this is automatically adjusted for the given machine. Note that due to Java's use of Compressed Oops, specifying a max heap size of between 32-48GB effectively reduces the memory available to GRIDSS so is strongly discouraged.

                # define the maxiumum number of threads so that each thread has 8Gb of ram (jvmheap)
                max_threads = max([1, int(allocated_ram/8 - 1)]) 
                if threads>max_threads: threads =  max_threads # this is to optimise for the reccommended level of parallelism

                # run
                print_if_verbose("running gridss on %s jvmheap and %i threads"%(jvmheap, threads))

                gridss_cmd = "%s --jar %s --reference %s -o %s --assembly %s --threads %i --workingdir %s --maxcoverage %i --blacklist %s --jvmheap %s %s > %s 2>&1"%(gridss_run, gridss_jar, reference_genome, gridss_VCFoutput, gridss_assemblyBAM, threads, gridss_tmpdir, maxcoverage, blacklisted_regions, jvmheap, sorted_bam_renamed, gridss_std)
                run_cmd(gridss_cmd)

                break

            except: print_if_verbose("WARNING: GRIDSS failed on try %i"%(Itry+1))

        # if it did not finish correctly finish it
        if file_is_empty(gridss_VCFoutput): raise ValueError("gridss did not finish correctly after %i tries. Check the log in %s and the std in %s"%(max_tries, gridss_tmpdir, gridss_std))

        # keep the full log
        lines_with_full_log = [l for l in open(gridss_std, "r").readlines() if "Full log file is" in l]
        if len(lines_with_full_log)!=1: raise ValueError("there should only be one log file")
        full_gridss_log = lines_with_full_log[0].split("Full log file is:")[1].split()[0]
        dest_full_gridss_log = "%s/gridss_full_log.txt"%outdir
        os.rename(full_gridss_log, dest_full_gridss_log)

        # delete unnecessary files
        delete_folder(gridss_tmpdir); remove_file(gridss_assemblyBAM)

    # annotated simple events, the ones that can be predicted from each breakpoint, but only if there are some predicted events
    gridss_VCFoutput_with_simple_event = "%s.withSimpleEventType.vcf"%gridss_VCFoutput
    simple_event_std = "%s/simple_event_annotation.std"%outdir
    print_if_verbose("annotating simple events in the grids output. The std is in %s"%simple_event_std)
    n_breakends = len([l for l in open(gridss_VCFoutput, "r").readlines() if not l.startswith("#")])
    if (file_is_empty(gridss_VCFoutput_with_simple_event) or replace is True) and n_breakends>0 : run_cmd("%s %s > %s 2>&1"%(annotate_simpleEvents_gridssVCF_R, gridss_VCFoutput, simple_event_std))

    remove_file(simple_event_std)

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

def add_info_to_gridssDF(df, reference_genome, expected_fields={"allele_frequency", "allele_frequency_SmallEvent", "other_coordinates", "other_chromosome", "other_position", "other_orientation",  "inserted_sequence", "len_inserted_sequence", "length_event", "has_poly16GC", "length_inexactHomology", "length_microHomology", "overlaps_repeats"}, median_insert_size=500, median_insert_size_sd=50):

    """This function takes a gridss df and returns the same adding expected_fields. This adds if not already done."""

    if len(set(df.keys()).intersection(expected_fields))!=len(expected_fields):

        df = cp.deepcopy(df)
        print_if_verbose("Adding info to df_gridss")


        ####### ADD WHETHER THE VARIANTS OVERLAP REPEATS #######

        # get the repeats table
        repeats_table = get_repeat_maskerDF(reference_genome, threads=1, replace=False)[1]
        df["overlaps_repeats"] = get_series_variant_in_repeats(df, repeats_table, replace=False)

        ########################################################

        # add the allele frequencies
        df["allele_frequency"] = df.apply(lambda r: np.divide(r["DATA_VF"] , (r["DATA_VF"] + r["DATA_REF"] + r["DATA_REFPAIR"])), axis=1).apply(getNaN_to_0)
        df["allele_frequency_SmallEvent"] = df.apply(lambda r: np.divide(r["DATA_VF"] , (r["DATA_VF"] + r["DATA_REF"])), axis=1).apply(getNaN_to_0)

        # add data to calculate the other breakpoint

        # check that the chromosomes do not start or end with a ".", which would give an error in the calculation of the length of the event
        if any([c.startswith(".") or c.endswith(".")  for c in set(df["#CHROM"])]): raise ValueError("If chromosome names start or end with a '.' the parsing of the gridss output is incorrect.")

        def get_other_position(r):
            
            # in case that this is just a bp, without any other chromsomes
            if r["ALT"].startswith(".") or r["ALT"].endswith("."): return "%s:%i"%(r["#CHROM"], r["POS"])

            # in case it is a breakpoint
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

def get_gridssDF_filtered(df, reference_genome, min_Nfragments=8, min_af=0.005, wrong_INFOtags={"IMPRECISE"}, wrong_FILTERtags={"NO_ASSEMBLY"}, filter_polyGC=True, filter_noSplitReads=True, filter_noReadPairs=True, maximum_strand_bias=0.95, maximum_microhomology=50, maximum_lenght_inexactHomology=50, range_filt_DEL_breakpoints=[100, 800], min_length_inversions=40, dif_between_insert_and_del=5, max_to_be_considered_small_event=1000, min_size=50, add_columns=True, min_af_EitherSmallOrLargeEvent=0.0, min_QUAL=0, filter_overlappingRepeats=False ):

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
    min_QUAL is the minimum value of the 'QUAL' field
    filter_overlappingRepeats indicates whether to filter out the breakends that are overlapping any repeat

    add_columns indicates whether to add columns to the df, which may have been done before

    The default values are the ones used in the gridss-purple-linx pipeline to generate the somatic callset

    """

    ######## ADD COLUMNS TO THE DF FOR FURTHER CALCULATION ##########
    if add_columns is True: df = add_info_to_gridssDF(df, reference_genome)

    # define whether the variant is a small duplication or insertion. These have special filters
    df["is_small_DupDel"] = (df.INFO_SIMPLE_TYPE.isin({"DEL", "DUP"})) & (df.length_event<=max_to_be_considered_small_event)

    # get sets
    wrong_INFOtags = set(wrong_INFOtags)
    wrong_FILTERtags = set(wrong_FILTERtags)

    ############ APPLY THE FILTERS ###########

    idx = ((df.length_event>=min_size) &
       (df.QUAL>=min_QUAL) &
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
    if filter_overlappingRepeats: idx = idx & ~(df.overlaps_repeats)

    # return the filtered df
    return df[idx]

def get_gridssDF_filtered_from_filtersDict(df_gridss, filters_dict, reference_genome):

    """Takes a df gridss and returns the filtered one, according to filters_dict"""

    # debug the fact that there is no min_af_EitherSmallOrLargeEvent
    if "min_af_EitherSmallOrLargeEvent" not in filters_dict: filters_dict["min_af_EitherSmallOrLargeEvent"] = 0.0

    # get the filtered df
    df_filt = get_gridssDF_filtered(df_gridss, reference_genome, min_Nfragments=filters_dict["min_Nfragments"], min_af=filters_dict["min_af"], wrong_INFOtags=filters_dict["wrong_INFOtags"], wrong_FILTERtags=filters_dict["wrong_FILTERtags"], filter_polyGC=filters_dict["filter_polyGC"], filter_noSplitReads=filters_dict["filter_noSplitReads"], filter_noReadPairs=filters_dict["filter_noReadPairs"], maximum_strand_bias=filters_dict["maximum_strand_bias"], maximum_microhomology=filters_dict["maximum_microhomology"], maximum_lenght_inexactHomology=filters_dict["maximum_lenght_inexactHomology"], range_filt_DEL_breakpoints=filters_dict["range_filt_DEL_breakpoints"], min_length_inversions=filters_dict["min_length_inversions"], dif_between_insert_and_del=filters_dict["dif_between_insert_and_del"], max_to_be_considered_small_event=filters_dict["max_to_be_considered_small_event"], min_size=filters_dict["min_size"], add_columns=False, min_af_EitherSmallOrLargeEvent=filters_dict["min_af_EitherSmallOrLargeEvent"], min_QUAL=filters_dict["min_QUAL"], filter_overlappingRepeats=filters_dict["filter_overlappingRepeats"] )

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

            # define the std
            r_stdout = "%s/svVCF_analysis_log.out"%outdir

            if only_simple_conversion is True:

                print_if_verbose("Getting files for svVCF file. Mainly generating a bedpe file for breakpoints with some extra info, but simple. The std can be found in %s"%r_stdout)
                run_cmd("%s %s > %s 2>&1"%(analyze_svVCF_simple, svVCF, r_stdout))                

            else:
                
                print_if_verbose("Getting files for svVCF file. Mainly generating a bedpe file for breakpoints with some extra info. The std can be found in %s"%r_stdout)
                run_cmd("%s %s > %s 2>&1"%(analyze_svVCF, svVCF, r_stdout))

            remove_file(r_stdout)

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
        print_if_verbose("getting genome graph object")

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
            print_if_verbose("getting graph-based genome")

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
        print_if_verbose("genome graph got")

        # get the real ends of the chromosomes
        if df_bedpe is not None:

            print_if_verbose("get the ends of the chromosome")
            sorted_positions = sorted(all_positions)
            pos_to_nNeighbors = pd.Series(dict(zip(sorted_positions, map(lambda x: len(genome_graph.neighbors(x, mode="ALL")), sorted_positions))))

            # debug
            if any(pos_to_nNeighbors<1): raise ValueError("there are some unnconected nodes in the graph genome")
            if any(pos_to_nNeighbors>100000): raise ValueError("there are some very highly connected regions in the graph genome")
            if any(pd.isna(pos_to_nNeighbors)): raise ValueError("there are some NaNs in the graph genome")

            real_chromosome_end_nodes = set(pos_to_nNeighbors[pos_to_nNeighbors==1].index)
            print_if_verbose("There are %i telomeric nodes in the graph genome"%len(real_chromosome_end_nodes))

            # clean
            del pos_to_nNeighbors

        # generate a df that maps each position to the real position
        print_if_verbose("defining df_positions")
        positions_real = []
        chromosomes_real = []
        is_end_of_chr = []
        for chrom, lenChrom in chrom_to_lenSeq.items():
            positions_real += list(range(lenChrom))
            chromosomes_real += [chrom]*lenChrom
            is_end_of_chr += ([True] + [False]*(lenChrom-2) + [True])

        print_if_verbose("adding lists to df_positions")
        df_positions = pd.DataFrame()
        df_positions["chromosome"] =  chromosomes_real
        df_positions["real_position"] =  positions_real
        df_positions["offset"] = df_positions.chromosome.apply(lambda x: chrom_to_offset[x])
        df_positions["graph_position"] = df_positions.real_position + df_positions.offset

        # add the is_end_of_chr depending on the 
        if df_bedpe is not None: df_positions["is_end_of_chr"] = df_positions.graph_position.isin(real_chromosome_end_nodes)
        else: df_positions["is_end_of_chr"] = is_end_of_chr
        print_if_verbose("There are %i telomeric nodes in the graph genome"%sum(df_positions["is_end_of_chr"]))

        if set(df_positions.graph_position)!=all_positions: raise ValueError("There is a bad graph calculation of the positions")

        # save
        print_if_verbose("saving graph object that has %.10f Mb"%(sys.getsizeof(genome_graph)/1000000))
        genomeGraph_outfile_tmp = "%s.tmp"%genomeGraph_outfile
        save_object(genome_graph, genomeGraph_outfile_tmp)
        os.rename(genomeGraph_outfile_tmp, genomeGraph_outfile)

        print_if_verbose("saving df positions objects")
        genomeGraph_positions_df_tmp = "%s.tmp"%genomeGraph_positions_df
        save_object(df_positions, genomeGraph_positions_df_tmp)
        os.rename(genomeGraph_positions_df_tmp, genomeGraph_positions_df)

    else:
        print_if_verbose("loading graph genome")
        genome_graph = load_object(genomeGraph_outfile)
        df_positions = load_object(genomeGraph_positions_df)

    return genome_graph, df_positions


def run_clove_filtered_bedpe(bedpe_file, outfile, sorted_bam, replace=False, median_coverage=10, median_coverage_dev=1, check_coverage=True):

    """ Takes a bedpe file and a sorted bam and generates the clove output. df_cov is a df of the first 4 columns of the mpileup output. If sorted_bam and median_coverage are provided, there will be a check for TAN and DEL of coverage 

    every TAN or DEL outside median_coverage +- median_coverage_dev will be filtered out"""

    if file_is_empty(outfile) or replace is True:

        print_if_verbose("running clove")
        outfile_tmp = "%s.tmp"%outfile
        clove_std = "%s.std"%outfile
        cmd = "%s -jar %s -i %s BEDPE -b %s -o %s -c %i %i"%(JAVA, clove, bedpe_file, sorted_bam, outfile_tmp, median_coverage, median_coverage_dev)

        if check_coverage is False:

            print_if_verbose("avoid checking coverage")
            cmd += " -r "

        print_if_verbose("running clove. The std can ve found in %s"%clove_std)

        # add the std
        cmd = "%s > %s 2>&1"%(cmd, clove_std)

        run_cmd(cmd)
        remove_file(clove_std)
        os.rename(outfile_tmp, outfile)
        
def get_distance_to_telomere_series(df_chromosome_position, genome_graph, df_positions_graph):

    """This function takes a df with chromosome and position and returns the distance to the telomere according to the genome graph"""

    print_if_verbose("getting distance to the telomere")

    # rename the df to have chrom and pos
    df_chromosome_position = df_chromosome_position.rename(columns=dict(zip(df_chromosome_position.columns, ["chromosome", "position"])))

    # add the graph positions
    df_chromosome_position = df_chromosome_position.merge(df_positions_graph, left_on=["chromosome", "position"], right_on=["chromosome", "real_position"], how="left", validate="one_to_one")

    # define all the positions
    all_positions = sorted(set(df_chromosome_position.graph_position))

    # define the positions that are ends of chromosomes
    chrom_end_positions = sorted(set(df_positions_graph[df_positions_graph.is_end_of_chr].graph_position))

    # calculate the distance from each position to the 
    print_if_verbose("calculating shortest paths in genome graph from %i end positions to %i positions. This may take a lot if there are a lot of end positions"%(len(chrom_end_positions), len(all_positions)))
    shortestPath_lengths_df = pd.DataFrame(genome_graph.shortest_paths(source=chrom_end_positions, target=all_positions, mode="IN"), columns=all_positions, index=chrom_end_positions)
    distance_telomere_series = shortestPath_lengths_df.apply(min, axis=0).apply(int)
    print_if_verbose("distance calculated")

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
        print_if_verbose("getting GC content df per position")

        print_if_verbose("calculating %s content"%target_nts)

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

    print_if_verbose("Getting GC content")

    if file_is_empty(gcontent_outfile) or replace is True:

        # define the initial index
        initial_index = list(df_windows.index)

        # resort
        df_windows = df_windows.sort_values(by=["chromosome", "start", "end"]).set_index(["chromosome", "start", "end"], drop=False)

        print_if_verbose("getting GC content for %i new windows"%len(df_windows))

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
    print_if_verbose("getting coverage-predictor function")

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
        print_if_verbose("investigating %s"%type_genome)

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

        # if the filtering is useless, use all the df. This is a way to skip the modelling.
        if len(df_correct)==0: df_correct = df_g

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
            fig.savefig(outfile, bbox_inches="tight")

        # get the rsquare of the model
        r2 = r2_score(df_correct.coverage, df_correct.coverage_from_dist_to_telomere_and_GC_content)
        print_if_verbose("The rsquare for %s is %.3f"%(type_genome, r2))

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

    if len(df)==0: return df

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
        coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file=bed_windows_to_measure, replace=replace, run_in_parallel=run_in_parallel, delete_bams=delete_bams), sep="\t").rename(columns={"#chrom":"chromosome"}).set_index(["chromosome", "start", "end"], drop=False)

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


    print_if_verbose("getting coverage relative to the one predicted from seq features")

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

    # define general names
    region_name = "%s_region"%region

    # define the minimum region length
    min_region_len = 100

    #### readjust for regions that are close to the telomere ####

    # if the start is 1 , set the 5' region to the 3' region
    if region=="5" and r["start"]<min_region_len: region="3"

    # if the region in the 3' is not larger than  enough, just define the 5' region
    elif region=="3" and (maxPos-r["end"])<min_region_len: region="5"

    #############################################################

    # if the region spans the whole chromosome, just set whole region as the 'region'
    if r["start"]<=min_region_len and (maxPos-r["end"])<=min_region_len: 
        start = r["start"]
        end = r["end"]

    # get the 5' region
    elif region=="5":

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

    # if the start is after the end, exit
    if start>=end: 
        print_if_verbose(r, region_name, maxPos, start, end)
        raise ValueError("start after end")

    # return a series of all important fields
    return pd.Series({"chromosome":r["chromosome"], "start":start, "end":end, "region_name":region_name})

def get_df_with_coverage_per_windows_relative_to_neighbor_regions(df_windows, bed_windows_prefix, reference_genome, sorted_bam, df_clove, median_coverage, replace=True, run_in_parallel=True, delete_bams=True):

    """Takes a df with windows of the genome and returns it with the coverage and the relative to the genome. It returns a df with several relative coverage measures."""

    print_if_verbose("getting coverage relative to neighbors")

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
    if len(df_clove)>0: all_svtypes = set(df_clove.SVTYPE)
    else: all_svtypes = set()

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
                merged_df_svtype = svtype_df_clove.merge(coverage_df, left_on=coord_fields, right_on=["chromosome", "start", "end"], how="left")

                # check that it is correct
                if len(merged_df_svtype)!=len(svtype_df_clove) or any(pd.isna(merged_df_svtype.target_coverage)): raise ValueError("There is an error with the merge")
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

    # initialize the already used breakpoints
    already_used_bps = set()

    # go through each combination
    for I1 in df_ITX1.index:
        I1s = df_ITX1.loc[I1]

        for I2 in df_ITX2.index:
            I2s = df_ITX2.loc[I2]

            # ask whether this combination can be close
            if I1s["#CHROM"]==I2s["#CHROM"] and I1s["CHR2"]==I2s["CHR2"] and abs(I1s["POS"]-I2s["POS"])<=tol_bp and abs(I1s["END"]-I2s["END"])<=tol_bp:

                # if any of the bp has been used, skip
                if I1 in already_used_bps or I2 in already_used_bps: continue

                # keep the already used bps
                already_used_bps.update({I1, I2})

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

    # initialize the already used breakpoints
    already_used_bps = set()

    # go through each combination
    for I1 in df_INVTX1.index:
        I1s = df_INVTX1.loc[I1]

        for I2 in df_INVTX2.index:
            I2s = df_INVTX2.loc[I2]

            # ask whether this combination can be close
            if I1s["#CHROM"]==I2s["#CHROM"] and I1s["CHR2"]==I2s["CHR2"] and abs(I1s["POS"]-I2s["POS"])<=tol_bp and abs(I1s["END"]-I2s["END"])<=tol_bp:

                # if any of the bp has been used, skip
                if I1 in already_used_bps or I2 in already_used_bps: continue

                # keep the already used bps
                already_used_bps.update({I1, I2})

                # keep in a way that it is in the same order as in the clove VCF
                data_dict[(I1, I2)] = {"ChrA":I1s["#CHROM"], "StartA":0, "EndA":I1s["POS"], "ChrB":I1s["CHR2"], "StartB":I1s["END"], "EndB":chr_to_len[I1s["CHR2"]], "Balanced":True} # all the starts are 0s

    df = pd.DataFrame(data_dict).transpose()

    return df

def get_bedpe_for_clovebalTRA_5with3(r, chr_to_len):

    """Takes a row of the df_balTRA_5with3 df and returns a bedpe row, sorted"""

    # define the fields
    if "#CHROM" in r.keys():
        chrA_field = "#CHROM"
        endA_field = "POS"
        chrB_field = "CHR2"
        startB_field = "END"

    else:
        chrA_field = "ChrA"
        endA_field = "EndA"
        chrB_field = "ChrB"
        startB_field = "StartB"

    return pd.Series({"ChrA":r[chrA_field], "StartA":0, "EndA":r[endA_field], "ChrB":r[chrB_field], "StartB":r[startB_field], "EndB":chr_to_len[r[chrB_field]]-1, "Balanced":True})

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
    print_if_verbose("getting SVs from clove")

    # initialize as a copy
    df_clove_initial = cp.deepcopy(df_clove)
    df_clove = cp.deepcopy(df_clove)

    # initialize the final dict
    svtype_to_svfile = {}

    # initialize the considered idxs
    df_clove.index = list(range(len(df_clove)))
    considered_idxs = []

    # initialize a df with the considered IDX types
    series_considered_idxs = pd.Series()

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
        series_considered_idxs["balTRA_5with5_or_3with3"] = make_flat_listOflists(df_balTRA_5with5_or_3with3.index)

        # balanced translocations 5with3 (these are the ones with an IVD field, assigned by clove)
        df_balTRA_5with3_IVD = df_clove[(df_clove.SVTYPE=="IVD") & ((df_clove.START - df_clove.END)<=tol_bp)].apply(lambda r: get_bedpe_for_clovebalTRA_5with3(r, chr_to_len), axis=1)
        considered_idxs += list(df_balTRA_5with3_IVD.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]
        series_considered_idxs["balTRA_5with3_IVD"] = list(df_balTRA_5with3_IVD.index)

        # balanced translocations 5with3 where there are two close INVTX breakpoints
        df_balTRA_5with3_INVTXbps = get_bedpeDF_for_clovebalTRA_5with3_or_3with5_INVTXbreakpoints(df_clove, chr_to_len, tol_bp=tol_bp).apply(lambda r: get_bedpe_for_clovebalTRA_5with3(r, chr_to_len), axis=1)
        considered_idxs += make_flat_listOflists(df_balTRA_5with3_INVTXbps.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]
        series_considered_idxs["balTRA_5with3_INVTXbps"] = make_flat_listOflists(df_balTRA_5with3_INVTXbps.index)

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

        print_if_verbose("There are %i translocations"%len(df_tra))

    #############################

    ####### INVERSIONS ##########
    if "CIV" in set(df_clove.SVTYPE) and "inversions" in svtypes_to_consider:

        df_inversions = df_clove[df_clove.SVTYPE=="CIV"][["#CHROM", "POS", "END", "ID"]].rename(columns={"#CHROM":"Chr", "POS":"Start", "END":"End"})
        inversions_bed = "%s.inversions.bed"%fileprefix
        df_inversions.to_csv(inversions_bed, sep="\t", header=True, index=False)

        considered_idxs += list(df_inversions.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]
        svtype_to_svfile["inversions"] = inversions_bed
        series_considered_idxs["inversions"] = list(df_inversions.index)


        print_if_verbose("There are %i inversions"%len(df_inversions))


    #############################

    ####### INSERTIONS ########
    if any([x in set(df_clove.SVTYPE) for x in {"CID", "CIT", "DUP", "TRA"}]) and "insertions" in svtypes_to_consider:

        # get the df
        df_ins = df_clove[df_clove.SVTYPE.isin({"CID", "CIT", "DUP", "TRA"})][[f for f in ["#CHROM", "POS", "CHR2", "START", "END", "ID", "SVTYPE", "coverage_FILTER"] if f in df_clove.keys()]]

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
        series_considered_idxs["insertions"] = list(df_ins.index)

        print_if_verbose("There are %i insertions, %i of which are copy-and-paste"%(len(df_ins), sum(df_ins.Copied=="TRUE")))

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
            series_considered_idxs[typeSV] = list(df_svtype.index)

            print_if_verbose("There are %i %s"%(len(df_svtype), typeSV))

    # keep only df_clove that has not been already used, which should be done after each step
    df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

    ###########################


    # write the remaining events which are not easily assignable
    df_notAssigned = df_clove[~(df_clove.SVTYPE.isin(cloveSVtypes_easy_classification))]
    df_notAssigned_file = "%s.remaining.tab"%(fileprefix)
    df_notAssigned[["ID", "#CHROM", "POS", "CHR2", "START", "END", "SVTYPE"]].to_csv(df_notAssigned_file, sep="\t", header=True, index=False)
    svtype_to_svfile["remaining"] = df_notAssigned_file

    print_if_verbose("There are %i remaining SVs"%len(df_notAssigned))

    # at the end make sure that the considered idxs are unique
    if len(considered_idxs)!=len(set(considered_idxs)): 

        print_if_verbose(series_considered_idxs)

        print_if_verbose("These are the IDXs considered more than once, with the number:", [(k,v) for k,v in Counter(considered_idxs).items() if v!=1])

        # save the df_clove
        df_clove_filename = "%s_df_clove_initial.py"%fileprefix
        save_object(df_clove_initial, df_clove_filename)

        cmds_write = ["import sys; sys.path.insert(0, '%s')"%CWD,
                      "import sv_functions as fun",
                      "df_clove = fun.load_object('%s')"%df_clove_filename,
                      "fun.write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, '%s', '%s', '%s', tol_bp=%i, define_insertions_based_on_coverage='%s')"%(fileprefix, reference_genome, sorted_bam, tol_bp, define_insertions_based_on_coverage)]

        print_if_verbose("ERROR log: \n---saving the untouched df_clove into %s---\n\n. You can  with the following cmds to reproduce the ERROR:---\n%s\n---"%(df_clove_filename, "\n".join(cmds_write)))

        raise ValueError("ERROR: Some clove events are assigned to more than one cathegory. Check the insertions and translocations calling")

    # return the df_clove and the remaining SVs
    return df_clove, svtype_to_svfile

def merge_coverage_per_window_files_in_one(bamfile, bam_sufix=".coverage_per_window.tab"):

    """This function takes all files that start with bamfile and end with coverage_per_window, """

    print_if_verbose("merging coverage tables")

    # define prefixes
    bam_dir = get_dir(bamfile)
    fileprefix = get_file(bamfile) + bam_sufix

    # remove dirs
    dirs_to_remove = ["%s/%s"%(bam_dir, f) for f in os.listdir(bam_dir) if os.path.isdir("%s/%s"%(bam_dir, f)) and f.startswith(fileprefix) and len(os.listdir("%s/%s"%(bam_dir, f)))==0] 
    for f in dirs_to_remove: delete_folder(f)

    # unite files
    files_prefix = ["%s/%s"%(bam_dir, f) for f in os.listdir(bam_dir) if not file_is_empty("%s/%s"%(bam_dir, f)) and "temporary_file" not in f and f.startswith(fileprefix)]

    # if there are no files, just skip the writing of the 'coverage_per_window.tab' file
    if len(files_prefix)==0: return

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
    print_if_verbose("running gridss and clove")

    # if the median_coverage is -1, calculate it from the provided sorted_bam
    if median_coverage==-1: 

        print_if_verbose("recalculating median coverage")
        destination_dir = "%s.calculating_windowcoverage"%sorted_bam
        coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=run_in_parallel, delete_bams=True), sep="\t")

        median_coverage = np.median(coverage_df.mediancov_1)

    print_if_verbose("running gridss and clove with given parameter with %.2f min_rel_coverage_to_consider_dup"%min_rel_coverage_to_consider_dup)
    make_folder(working_dir)

    # edit the replace, regarding if filtering from the run of GRIDSS
    if replace is True and replace_FromGridssRun is False: replace_FromGridssRun = True

    # if there is a 

    # first obtain the gridss output if it is not provided
    if file_is_empty(gridss_VCFoutput) or replace is True: gridss_VCFoutput = run_gridss_and_annotateSimpleType(sorted_bam, reference_genome, working_dir, replace=replace, threads=threads, blacklisted_regions=gridss_blacklisted_regions, maxcoverage=gridss_maxcoverage)

    #################################################
    ##### GET A LIST OF FILTERED BREAKPOINTS ########
    #################################################

    
    # get the output of gridss into a df
    print_if_verbose("getting gridss")
    df_gridss = add_info_to_gridssDF(load_single_sample_VCF(gridss_VCFoutput), reference_genome, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd) # this is a dataframe with some info

    # filter according to gridss_filters_dict
    print_if_verbose("filtering gridss")
    df_gridss_filt = get_gridssDF_filtered_from_filtersDict(df_gridss, gridss_filters_dict, reference_genome)
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
    print_if_verbose("there are %i breakpoints"%len(df_bedpe))

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

    if len(df_clove)==0: return {}, df_gridss

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


def get_taxid2name(taxIDs):

    """Takes an iterable of taxIDs and returns a dict mapping each of them to the scientific name"""

    taxIDs = list(taxIDs)

    ncbi = NCBITaxa()
    taxid2name = ncbi.get_taxid_translator(taxIDs)

    return taxid2name

def get_allOxfordNanopore_runInfo_fromSRA_forDivision(fileprefix, taxID_division, reference_genome, taxIDs_to_exclude=set(), replace=False, min_coverage=30):

    """This function tages all the SRRs that have WGS for the required taxIDs"""

    SRA_runInfo_df_file = "%s.SRA_runInfo_df.py"%fileprefix

    if file_is_empty(SRA_runInfo_df_file) or replace is True:
        print_if_verbose("Getting WGS for %i"%taxID_division)

        # define the WGS fastq filters
        ONT_filters = '("biomol dna"[Properties] AND "strategy wgs"[Properties] AND "platform oxford nanopore"[Properties] AND ("strategy wgs"[Properties] OR "strategy wga"[Properties] OR "strategy wcs"[Properties] OR "strategy clone"[Properties] OR "strategy finishing"[Properties] OR "strategy validation"[Properties]))'
    
        # define the esearch query
        esearch_query = "(txid%i[Organism:exp]) AND %s"%(taxID_division, ONT_filters) 

        # add the taxIDs_to_exclude 
        if len(taxIDs_to_exclude)>0: esearch_query += " NOT (%s)"%(" OR ".join(["(txid%i[Organism:exp])"%ID for ID in taxIDs_to_exclude]))

        print_if_verbose("This is the esearch query (you can check it against the SRA):\n\n %s \n\n"%esearch_query)

        # get esearch
        esearch_outfile = "%s.esearch_output.txt"%fileprefix
        esearch_stderr = "%s.stderr"%esearch_outfile
        efetch_outfile = "%s.efetch_output.txt"%fileprefix
        efetch_stderr = "%s.stderr"%efetch_outfile

        columns_efetch = "Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash".split(",")

        # get the esearch
        print_if_verbose("running esearch. The stderr is in %s"%esearch_stderr)
        run_cmd("%s -db sra -query '%s' > %s 2>%s"%(esearch, esearch_query, esearch_outfile, esearch_stderr))
        remove_file(esearch_stderr)

        # get the number of rows
        nresults = [int(l.split("<Count>")[1].split("<")[0]) for l in open(esearch_outfile, "r").readlines() if "<Count>" in l and "</Count>" in l][0]
        print_if_verbose("There are %i results in esearch"%nresults)

        if nresults==0: SRA_runInfo_df = pd.DataFrame(columns=columns_efetch)
        else:

            # run efetch
            print_if_verbose("running efetch. The stderr is in %s"%efetch_stderr)
            run_cmd("cat %s | %s -db sra --format runinfo | egrep -v '^Run' | egrep 'https' > %s 2>%s"%(esearch_outfile, efetch, efetch_outfile, efetch_stderr))
            remove_file(efetch_stderr)


            # get into df
            SRA_runInfo_df = pd.read_csv(efetch_outfile, sep=",", header=None, names=columns_efetch)

        save_object(SRA_runInfo_df, SRA_runInfo_df_file)

    else: SRA_runInfo_df = load_object(SRA_runInfo_df_file)

    # return if empty 
    if len(SRA_runInfo_df)==0: return SRA_runInfo_df

    # filter out the undesired taxIDs
    SRA_runInfo_df = SRA_runInfo_df[~SRA_runInfo_df.TaxID.isin(taxIDs_to_exclude)]


    # define the expected coverage
    length_genome = sum(get_chr_to_len(reference_genome).values())
    SRA_runInfo_df["expected_coverage"] = SRA_runInfo_df.apply(lambda r: (r["spots"]*r["avgLength"])/length_genome,axis=1)

    # plot the number of spots
    filename = "%s.distribution_parameters.pdf"%fileprefix
    if file_is_empty(filename):

        fig = plt.figure(figsize=(5,13))
        for I, field in enumerate(["spots", "spots_with_mates", "avgLength", "InsertSize", "size_MB", "expected_coverage"]):
            ax = plt.subplot(6, 1, I+1)
            sns.distplot(SRA_runInfo_df[field], kde=False, rug=True)
            ax.set_xlabel(field)

        fig.tight_layout()  # otherwise the right y-label is slightly 
        fig.savefig(filename, bbox_inches='tight');
        plt.close(fig)

    # keep only those that have at least a coverage of min_coverage    
    SRA_runInfo_df = SRA_runInfo_df[SRA_runInfo_df["expected_coverage"]>=min_coverage]

    print_if_verbose("There are %i SRRs ready to use with at least %ix coverage"%(len(SRA_runInfo_df), min_coverage))

    return SRA_runInfo_df

def get_allWGS_runInfo_fromSRA_forDivision(fileprefix, taxID_division, reference_genome, taxIDs_to_exclude=set(), replace=False, min_coverage=30):

    """This function tages all the SRRs that have WGS for the required taxIDs"""

    SRA_runInfo_df_file = "%s.SRA_runInfo_df.py"%fileprefix

    if file_is_empty(SRA_runInfo_df_file) or replace is True:
        print_if_verbose("Getting WGS for %i"%taxID_division)

        # define the WGS fastq filters
        WGS_filters = '("biomol dna"[Properties] AND "strategy wgs"[Properties] AND "library layout paired"[Properties] AND "platform illumina"[Properties] AND "strategy wgs"[Properties] OR "strategy wga"[Properties] OR "strategy wcs"[Properties] OR "strategy clone"[Properties] OR "strategy finishing"[Properties] OR "strategy validation"[Properties])'

        # define the esearch query
        esearch_query = "(txid%i[Organism:exp]) AND %s"%(taxID_division, WGS_filters) 

        # add the taxIDs_to_exclude 
        if len(taxIDs_to_exclude)>0: esearch_query += " NOT (%s)"%(" OR ".join(["(txid%i[Organism:exp])"%ID for ID in taxIDs_to_exclude]))

        print_if_verbose("This is the esearch query (you can check it against the SRA):\n\n %s \n\n"%esearch_query)

        # get esearch
        esearch_outfile = "%s.esearch_output.txt"%fileprefix
        esearch_stderr = "%s.stderr"%esearch_outfile
        efetch_outfile = "%s.efetch_output.txt"%fileprefix
        efetch_stderr = "%s.stderr"%efetch_outfile

        columns_efetch = "Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash".split(",")

        # get the esearch
        print_if_verbose("running esearch. The stderr is in %s"%esearch_stderr)
        run_cmd("%s -db sra -query '%s' > %s 2>%s"%(esearch, esearch_query, esearch_outfile, esearch_stderr))
        remove_file(esearch_stderr)

        # get the number of rows
        nresults = [int(l.split("<Count>")[1].split("<")[0]) for l in open(esearch_outfile, "r").readlines() if "<Count>" in l and "</Count>" in l][0]
        print_if_verbose("There are %i results in esearch"%nresults)

        if nresults==0: SRA_runInfo_df = pd.DataFrame(columns=columns_efetch)
        else:

            # run efetch
            print_if_verbose("running efetch. The stderr is in %s"%efetch_stderr)
            run_cmd("cat %s | %s -db sra --format runinfo | egrep -v '^Run' | egrep 'https' > %s 2>%s"%(esearch_outfile, efetch, efetch_outfile, efetch_stderr))
            remove_file(efetch_stderr)

            # get into df
            SRA_runInfo_df = pd.read_csv(efetch_outfile, sep=",", header=None, names=columns_efetch)

        save_object(SRA_runInfo_df, SRA_runInfo_df_file)

    else: SRA_runInfo_df = load_object(SRA_runInfo_df_file)

    # return if empty 
    if len(SRA_runInfo_df)==0: return SRA_runInfo_df

    # filter out the undesired taxIDs
    SRA_runInfo_df = SRA_runInfo_df[~SRA_runInfo_df.TaxID.isin(taxIDs_to_exclude)]

    # return if empty 
    if len(SRA_runInfo_df)==0: return SRA_runInfo_df

    # define the expected coverage
    length_genome = sum(get_chr_to_len(reference_genome).values())
    SRA_runInfo_df["expected_coverage"] = SRA_runInfo_df.apply(lambda r: (r["spots_with_mates"]*r["avgLength"])/length_genome,axis=1)

    # plot the number of spots
    filename = "%s.distribution_parameters.pdf"%fileprefix
    if file_is_empty(filename):

        fig = plt.figure(figsize=(5,13))
        for I, field in enumerate(["spots", "spots_with_mates", "avgLength", "InsertSize", "size_MB", "expected_coverage"]):
            ax = plt.subplot(6, 1, I+1)
            sns.distplot(SRA_runInfo_df[field], kde=False, rug=True)
            ax.set_xlabel(field)

        fig.tight_layout()  # otherwise the right y-label is slightly 
        fig.savefig(filename, bbox_inches='tight');
        plt.close(fig)

    # keep only those that have at least a coverage of min_coverage    
    SRA_runInfo_df = SRA_runInfo_df[SRA_runInfo_df["expected_coverage"]>=min_coverage]

    print_if_verbose("There are %i SRRs ready to use with at least %ix coverage"%(len(SRA_runInfo_df), min_coverage))

    return SRA_runInfo_df

def get_allWGS_runInfo_fromSRA_forTaxIDs(fileprefix, taxIDs, reference_genome, replace=False, min_coverage=30):

    """This function tages all the SRRs that have WGS for the required taxIDs"""

    SRA_runInfo_df_file = "%s.SRA_runInfo_df.py"%fileprefix

    if file_is_empty(SRA_runInfo_df_file) or replace is True:

        # define the WGS fastq filters
        WGS_filters = '("biomol dna"[Properties] AND "strategy wgs"[Properties] AND "library layout paired"[Properties] AND "platform illumina"[Properties] AND "strategy wgs"[Properties] OR "strategy wga"[Properties] OR "strategy wcs"[Properties] OR "strategy clone"[Properties] OR "strategy finishing"[Properties] OR "strategy validation"[Properties])'


        # get the name of these taxIDs
        taxid2name = get_taxid2name(taxIDs)

        # define the esearch query
        organism_filters = " OR ".join(['"%s"[orgn:__txid%i]'%(name.split()[0], taxID) for taxID, name in taxid2name.items()])
        esearch_query = "(%s) AND %s"%(organism_filters, WGS_filters) 

        # get esearch
        efetch_outfile = "%s.efetch_output.txt"%fileprefix

        columns_efetch = "Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash".split(",")

        # if there are no runs, it will run an error
        esearch_efetch_stderr = "%s.generating.stderr"%efetch_outfile
        print_if_verbose("Querying the SRA database. This will throw an error if there are no results. The stderr is in %s"%esearch_efetch_stderr)
        run_cmd("%s -db sra -query '%s' | %s -db sra --format runinfo | egrep -v '^Run' | egrep 'https' > %s 2>%s"%(esearch, esearch_query, efetch, efetch_outfile, esearch_efetch_stderr))
        remove_file(esearch_efetch_stderr)

        SRA_runInfo_df = pd.read_csv(efetch_outfile, sep=",", header=None, names=columns_efetch)

        save_object(SRA_runInfo_df, SRA_runInfo_df_file)

    else: SRA_runInfo_df = load_object(SRA_runInfo_df_file)


    # define the expected coverage
    length_genome = sum(get_chr_to_len(reference_genome).values())
    SRA_runInfo_df["expected_coverage"] = SRA_runInfo_df.apply(lambda r: (r["spots_with_mates"]*r["avgLength"])/length_genome,axis=1)

    # plot the number of spots
    filename = "%s.distribution_parameters.pdf"%fileprefix
    fig = plt.figure(figsize=(5,13))
    for I, field in enumerate(["spots", "spots_with_mates", "avgLength", "InsertSize", "size_MB", "expected_coverage"]):
        ax = plt.subplot(6, 1, I+1)
        sns.distplot(SRA_runInfo_df[field], kde=False, rug=True)
        ax.set_xlabel(field)

    fig.tight_layout()  # otherwise the right y-label is slightly 
    fig.savefig(filename, bbox_inches='tight');
    plt.close(fig)

    # keep only those that have at least a coverage of min_coverage    
    SRA_runInfo_df = SRA_runInfo_df[SRA_runInfo_df["expected_coverage"]>=min_coverage]

    for field in ["AssemblyName", "SampleType", "TaxID"]: pass

    print_if_verbose("There are %i SRRs ready to use with at least %ix coverage"%(len(SRA_runInfo_df), min_coverage))

    return SRA_runInfo_df

def download_srr_subsetReads_onlyFastqDump(srr, download_dir, subset_n_reads=1000000):

    """This function downloads a subset of reads with fastqdump"""

    # define the previous folder
    folder_upstream_download_dir = get_dir(download_dir); make_folder(folder_upstream_download_dir)

    # define a tmp_folder
    download_dir_tmp = "%s_tmp"%download_dir
    make_folder(download_dir_tmp)

    # define the reads
    reads1 = "%s/%s_1.fastq.gz"%(download_dir, srr)
    reads2 = "%s/%s_2.fastq.gz"%(download_dir, srr)

    # define the std
    fastqdump_std = "%s/std.txt"%download_dir_tmp

    for Itry in range(3):

        if file_is_empty(reads1) or file_is_empty(reads2):

            # define previous runs
            delete_folder(download_dir_tmp)
            make_folder(download_dir_tmp)

            # run dump
            print_if_verbose("running fastqdump. The std is in %s"%fastqdump_std)
            fastqdump_cmd = "%s --split-3 --gzip --maxSpotId %i --outdir %s %s > %s 2>&1"%(fastqdump, subset_n_reads, download_dir_tmp, srr, fastqdump_std)
            
            try:
                run_cmd(fastqdump_cmd)

                # move the tmp to the final
                os.rename(download_dir_tmp, download_dir)

            except: print_if_verbose("fastqdump did not work this time. Retrying up to 3 times")

    if any([file_is_empty(x) for x in [reads1, reads2]]): raise ValueError("You could not download properly %s. Try running again"%srr)

    remove_file(fastqdump_std)

    return reads1, reads2

def run_freebayes_for_region(region, outvcf_folder, ref, sorted_bam, ploidy, coverage, replace, pooled_sequencing):

    """Takes a region (chromosome:start-end) and the fasta file and an outvcf and runs freebayes on it"""

    # define the output vcf file
    outvcf = "%s/%s_freebayes.vcf"%(outvcf_folder, region); outvcf_tmp = "%s.tmp.vcf"%outvcf
    print_if_verbose("running freebayes for %s"%region)

    # remove previously existing files
    if file_is_empty(outvcf) or replace is True:
        remove_file(outvcf_tmp)

        # run freebayes
        freebayes_std = "%s.std"%outvcf_tmp

        if pooled_sequencing is True:
            print_if_verbose("running for pooled data")
            run_cmd("%s -f %s --haplotype-length -1 --use-best-n-alleles 20 --min-alternate-count %i --min-alternate-fraction 0 --pooled-continuous -b %s -v %s --region %s > %s 2>&1"%(freebayes, ref, coverage, sorted_bam, outvcf_tmp, region, freebayes_std))
        else: 
            print_if_verbose("running unpooled sequencing")
            run_cmd("%s -f %s -p %i --min-coverage %i -b %s --haplotype-length -1 -v %s --region %s > %s 2>&1"%(freebayes, ref, ploidy, coverage, sorted_bam, outvcf_tmp, region, freebayes_std))

        # remove the intermediate files
        remove_file(freebayes_std)

        # rename
        os.rename(outvcf_tmp, outvcf)

    # return the vcfs
    return outvcf


def leftTrimVariant(pos, ref, alt, onlyOneBp=True):

    """Takes a variant with a position, ref and alt alleles and returns the trimmed to the left. onlyOneBp means that it will only trim one basepair. Fo some reason this is how the vcflib and gatk normalisations have the trims. """

    # get the lengths
    len_ref = len(ref); len_alt = len(alt)

    # those that start with the same base can be trimmed
    if ref[0]==alt[0]:

        len_overlap = 0
        # go through each position and keep how long is the overlap on the left
        if onlyOneBp is True: rangelength = 1 # only the first bp
        else: rangelength = min([len_ref, len_alt]) # all

        for i in range(rangelength):
            if ref[i]==alt[i]: len_overlap+=1
            else: break

        # define the modified position, which increments depending on the length of the overlap
        mod_pos = pos + len_overlap

        # define the ref
        if len_overlap==len_ref: mod_ref = "-"
        else: mod_ref = ref[len_overlap:]

        # define the alt the same way
        if len_overlap==len_alt: mod_alt = "-"
        else: mod_alt = alt[len_overlap:]

    else: mod_pos, mod_ref, mod_alt = pos, ref, alt # the same

    return mod_pos, mod_ref, mod_alt



def load_vep_table_intoDF(vep_filename):

    """Takes the tabular output of ensembl vep and returns a table where the index is (chr, pos, ref, alt)"""

    # load normal df
    vep_df = pd.read_csv(vep_filename, sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)

    # keep the orignial uploaded variation
    vep_df["#Uploaded_variation_original"] = vep_df["#Uploaded_variation"]

    # change the uploaded variation so that it includes only the allele position
    vep_df["#Uploaded_variation"] = vep_df["#Uploaded_variation"].apply(lambda x: x.split("/")[0]) + "/" + vep_df.Allele

    # add the index of the alternative allele which corresponds to this var. For example, if the gt is 1|1|1, and this var has a 1, it will mean that it is the gt of the middle.
    print_if_verbose("geting GT index")
    vep_df["GT_index"] = vep_df[["#Uploaded_variation_original", "Allele"]].apply(lambda r: [Iallele+1 for Iallele, alt in enumerate(r["#Uploaded_variation_original"].split("/")[1:]) if r["Allele"]==alt][0], axis=1)

    #/home/mschikora/samba/Cglabrata_antifungals/VarCall/VarCallOutdirs/RUN1_EF1620_7F_ANIFLZ_VarCallresults/freebayes_ploidy1_out/output.filt.norm_vcflib.vcf_annotated.tab

    def get_idx_from_row(row):

        """Thakes a row of vep_df and returns the index"""

        chrom = row["Location"].split(":")[0]
        pos_str = row["#Uploaded_variation"].split("_")[-2]
        ref = row["#Uploaded_variation"].split("_")[-1].split("/")[0]
        pos = int(pos_str)
        alt = row["Allele"]

        return (chrom, pos, ref, alt)

    vep_df.index = vep_df.apply(get_idx_from_row, axis=1)
    vep_df["chromosome"] = [x[0] for x in vep_df.index] 
    vep_df["position"] = [x[1] for x in vep_df.index]
    vep_df["ref"] = [x[2] for x in vep_df.index]
    vep_df["alt"] = [x[3] for x in vep_df.index]

    # check that the index and uploaded var are the same
    if len(set(vep_df.index))!=len(set(vep_df["#Uploaded_variation"])): raise ValueError("each chrom/pos/ref/alt combination does not match the uploaded variation")

    return vep_df



def load_vcf_intoDF_GettingFreq_AndFilter(vcf_file):
    
    """Takes a vcf and loads it into a pandas dataframe. It assumes that it is a vcf with only one sample, at the end. It returns the allele frequency of all the representations of the variants. At least those that have been tested for the VEP output. These are important filters:
    
    # for all
    DP (from DATA, one for all alleles)
    QUAL (one for all)

    # specific from freebayes
    SAF (from INFO, one for each alleles)
    SAR (from INFO, one for each alleles)
    RPR (from INFO, one for each allele)
    RPL (from INFO, one for each allele)

    # specific from HC:
    QD (from INFO)
    MQ (from INFO)
    FS (from INFO)
    MQRankSum (from INFO)
    ReadPosRankSum (from INFO)

    """

    # get the df (avoid NA as a default NaN)

    df = pd.read_csv(vcf_file, skiprows=list(range(len([line for line in open(vcf_file, "r", encoding='utf-8', errors='ignore') if line.startswith("##")]))), sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)

    # set the index to be a tuple of (chromosome, location, ref, alt)
    df["CHROM_POS_REF_ALT"] = [tuple(x) for x in df[["#CHROM", "POS", "REF", "ALT"]].values]; df = df.set_index("CHROM_POS_REF_ALT")

    # add a colum that will result from the merging of FORMAT and the last column (which are the values of FORMAT)
    data_colname = list(df.keys())[-1]
    df["METADATA"] = [dict(zip(x[0].split(":"), x[1].split(":"))) for x in df[["FORMAT", data_colname]].values]
    features = df.iloc[0].METADATA.keys()

    # add as columns all the fetaures
    for feature in features: 

        # go through each data record
        data = []
        for rec in df.METADATA:

            if feature in rec: data.append(rec[feature])
            else: data.append("")
        df[feature] = data

    ##### ADD INFO COLS #####

    interesting_INFO_filters = {"SAF", "SAR", "RPR", "RPL", "QD", "MQ", "FS", "MQRankSum", "ReadPosRankSum"}
    df["INFO_dict"] = df.INFO.apply(lambda x: {item.split("=")[0] : item.split("=")[1].split(",") for item in x.split(";") if item.split("=")[0] in interesting_INFO_filters})

    def get_INFOfilt(r, INFO_filt):
        
        if INFO_filt in set(r["INFO_dict"].keys()): return r["INFO_dict"][INFO_filt]
        else: return [np.nan]*len(r["ALT"].split(","))

    for INFO_filt in interesting_INFO_filters: df["INFO_%s_list"%INFO_filt] = df.apply(lambda r: get_INFOfilt(r, INFO_filt), axis=1)

    #########################


    # a function that avoids nans
    def returnNonNAN(float_val):
        if pd.isna(float_val): return 0.0000001
        else: return float_val 

    # calculate the real allelle frequency
    def calculate_altAllele_freq(x):

        """Takes AD and returns a tuple with the alternagtive alleles frequencies"""
        alleles_reads = [int(nr) for nr in x.split(",")]
        sum_alleles = sum(alleles_reads)

        return [returnNonNAN(np.divide(altAlleleReads, sum_alleles)) for altAlleleReads in alleles_reads[1:]]

    df["alternative_allelle_frequencies"] = df.AD.apply(calculate_altAllele_freq)

    # get the frequencies into a dictionary
    print_if_verbose("calculating features for multialleles")

    # define the in

    def get_dicts_metadata(row):

        """Takes a row of df and returns a dictionary that contains var_to_frequency, var_to_filter and var_to_GT for the variant of the row, split by ",", and the left-trimmed one. It also returns a dict mapping each var to the filters that are important. 
        """

        # initialize dicts
        var_to_frequency = {} 
        var_to_filter = {}
        var_to_GT = {}
        var_to_filters = {} # these are some potentially interesting filters, where each allele gets one

        # get the tuple var
        var = row.name
        chrom, pos, ref, alt = var

        # get all alternative alleles
        alt_alleles = alt.split(",")

        # get all the alternative frequencies
        altAllele_to_freq = dict(zip(alt_alleles, row["alternative_allelle_frequencies"]))
        INFOfilt_to_altAllele_to_val = {INFOfilt : dict(zip(alt_alleles , row["INFO_%s_list"%INFOfilt])) for INFOfilt in interesting_INFO_filters}

        # map them to freqs and tags
        for real_alt, freq in altAllele_to_freq.items(): 

            # add the interesting filters from the INFO field, QUAL and DP
            filters_str = ";".join(["%s=%s"%(INFOfilt, INFOfilt_to_altAllele_to_val[INFOfilt][real_alt]) for INFOfilt in interesting_INFO_filters] + ["QUAL=%s"%row["QUAL"], "DP=%i"%int(row["DP"])])
         
            # first the direct representation of the variant
            var_to_frequency[(chrom, pos, ref, real_alt)] = freq
            var_to_filter[(chrom, pos, ref, real_alt)] = row["FILTER"]
            var_to_GT[(chrom, pos, ref, real_alt)] = str(row["GT"])
            var_to_filters[(chrom, pos, ref, real_alt)] = filters_str

            # now the left-most trimmed
            mod_pos, mod_ref, mod_alt = leftTrimVariant(pos, ref, real_alt)

            # keep if changed
            if (mod_pos, mod_ref, mod_alt)!=(pos, ref, real_alt): 
                untrimmed_var = (chrom, pos, ref, real_alt)

                var_to_frequency[(chrom, mod_pos, mod_ref, mod_alt)] = var_to_frequency[untrimmed_var]
                var_to_filter[(chrom, mod_pos, mod_ref, mod_alt)] = var_to_filter[untrimmed_var]
                var_to_GT[(chrom, mod_pos, mod_ref, mod_alt)] = var_to_GT[untrimmed_var]
                var_to_filters[(chrom, mod_pos, mod_ref, mod_alt)] = var_to_filters[untrimmed_var]

        return {"var_to_frequency":var_to_frequency, "var_to_filter":var_to_filter, "var_to_GT":var_to_GT, "var_to_filters":var_to_filters}

    series_dicts = df.apply(get_dicts_metadata, axis=1)

    # get the union of all dicts
    print_if_verbose("Geting final dicts")
    var_to_frequency = dict(j for i in series_dicts.apply(lambda d: d["var_to_frequency"]) for j in i.items())
    var_to_filter = dict(j for i in series_dicts.apply(lambda d: d["var_to_filter"]) for j in i.items())
    var_to_GT = dict(j for i in series_dicts.apply(lambda d: d["var_to_GT"]) for j in i.items())
    var_to_filters = dict(j for i in series_dicts.apply(lambda d: d["var_to_filters"]) for j in i.items())

    # get only the columns that you want to keep the real vcf
    return df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",  data_colname]], var_to_frequency, var_to_filter, var_to_GT, var_to_filters


def run_freebayes_parallel_regions(outdir_freebayes, ref, sorted_bam, ploidy, coverage, threads=4, replace=False, pooled_sequencing=False):

    """It parallelizes over the provided threads of the system"""

    # make the dir if not already done
    make_folder(outdir_freebayes)

    #run freebayes
    freebayes_output ="%s/output.raw.vcf"%outdir_freebayes; freebayes_output_tmp = "%s.tmp"%freebayes_output
    if file_is_empty(freebayes_output) or replace is True:

        print_if_verbose("running freebayes in parallel with %i threads"%(threads))

        # define the regions file
        window_fb = 10000
        regions_file = "%s/regions_genome_%ibp.tab"%(outdir_freebayes, window_fb)
        regions_file_stderr = "%s.generating.stderr"%regions_file
        print_if_verbose("getting regions for freebayes run in parallel. The stderr is in %s"%regions_file_stderr)
        run_cmd("%s %s.fai %i > %s 2>%s"%(fasta_generate_regions_py, ref, window_fb, regions_file, regions_file_stderr))
        remove_file(regions_file_stderr)
        regions = [l.strip() for l in open(regions_file, "r").readlines()]

        # remove the previous tmp file
        remove_file(freebayes_output_tmp)

        # make a dir to store the vcfs
        regions_vcfs_dir = "%s/regions_vcfs"%outdir_freebayes; make_folder(regions_vcfs_dir)

        # define the inputs of the function
        inputs_fn = [(region, regions_vcfs_dir, ref, sorted_bam, ploidy, coverage, replace, pooled_sequencing) for region in regions]

        # initialize the pool class with the available CPUs --> this is asyncronous parallelization
        with multiproc.Pool(threads) as pool:

            # run in parallel the freebayes generation for all the 
            regions_vcfs = pool.starmap(run_freebayes_for_region, inputs_fn)

            # close the pool
            pool.close()
            pool.terminate()

        # go through each of the chromosomal vcfs and append to a whole df
        print_if_verbose("appending all vcfs of the individual regions together")
        all_df = pd.DataFrame()
        all_header_lines = []
        for Iregion, vcf in enumerate(regions_vcfs):
            print_if_verbose("Working on region %i/%i"%(Iregion+1, len(regions_vcfs)))

            # load the df keeping the header lines
            header_lines = [l for l in open(vcf, "r") if l.startswith("##")]
            df = pd.read_csv(vcf, sep="\t", header = len(header_lines))

            # define the vcf header
            vcf_header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            sample_header = [c for c in df.columns if c not in vcf_header][0]
            vcf_header.append(sample_header)

            # keep header that is unique
            all_header_lines.append("".join([line for line in header_lines if line.split("=")[0] not in {"##reference", "##commandline", "##fileDate"}]))
            
            # append to the previous df
            all_df = all_df.append(df[vcf_header], sort=True)

        # sort and remove duplicate entries
        all_df = all_df.sort_values(by=["#CHROM", "POS"]).drop_duplicates(subset=["#CHROM", "POS", "REF", "ALT"])

        # check that all headers are the same
        if len(set(all_header_lines))!=1: 
            print_if_verbose("These are the header lines: ", set(all_header_lines))
            print_if_verbose("There are %i unique headers"%len(set(all_header_lines)))
            raise ValueError("Not all headers are the same in the individual chromosomal vcfs. This may indicate a problem with parallelization of freebayes")

        # write the file
        open(freebayes_output_tmp, "w").write(all_header_lines[0] + all_df[vcf_header].to_csv(sep="\t", index=False, header=True))

        # remove tmp vcfs
        delete_folder(regions_vcfs_dir)

        # rename
        os.rename(freebayes_output_tmp, freebayes_output)

    # filter the freebayes by quality
    freebayes_filtered = "%s/output.filt.vcf"%outdir_freebayes; freebayes_filtered_tmp = "%s.tmp"%freebayes_filtered
    if file_is_empty(freebayes_filtered) or replace is True:

        if pooled_sequencing is True: soft_link_files(freebayes_output, freebayes_filtered)
        else:

            filtering_stderr = "%s.generating.stderr"%freebayes_filtered_tmp
            print_if_verbose("filtering freebayes. The stderr is in %s"%filtering_stderr)
            cmd_filter_fb = '%s -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" --tag-pass PASS %s > %s 2>%s'%(vcffilter, freebayes_output, freebayes_filtered_tmp, filtering_stderr); run_cmd(cmd_filter_fb)
            remove_file(filtering_stderr)
            os.rename(freebayes_filtered_tmp, freebayes_filtered)

    return freebayes_filtered




def get_set_adapter_fastqc_report(fastqc_report):

    """Takes a fastqc report and returns a set with the adapters"""
    with open(fastqc_report, "r") as fd:

        # look for the line with the content of the table
        for line in fd:
            if line.startswith("</style>"):

                if "No overrepresented sequences" in line: overrepresented_seqs = set()
                else:
                    # split by overrepresented seqs
                    after_OR_seqs = line.split("Overrepresented sequences")[2]
                    overrepresented_seqs = set([x for x in re.split("</?td>", after_OR_seqs) if all([char in {"A", "C", "T", "G"} for char in x]) and len(x)>0])
                break

    return overrepresented_seqs

def run_trimmomatic(reads1, reads2, replace=False, threads=1):

    """Trims the reads and returns the paths to the trimmed reads"""

    # define the trimmmed reads dir
    trimmed_reads1 = "%s.trimmed.fastq.gz"%reads1
    trimmed_reads2 = "%s.trimmed.fastq.gz"%reads2

    if file_is_empty(trimmed_reads1) or file_is_empty(trimmed_reads2) or replace is True:

        # initialize all the html files
        all_html_files = []

        # run fastqc to get the adapters
        for reads in [reads1, reads2]:

            fastqc_dir = "%s_fastqc_dir"%(reads); make_folder(fastqc_dir)
            html_files = ["%s/%s"%(fastqc_dir, x) for x in os.listdir(fastqc_dir) if x.endswith(".html")]

            if len(html_files)==0 or replace is True:

                std_fastqc = "%s/std.txt"%fastqc_dir
                print_if_verbose("running fastqc. The std is in %s"%std_fastqc)
                run_cmd("%s -o %s --threads %i --extract --java %s %s > %s 2>&1"%(FASTQC, fastqc_dir, threads, JAVA, reads, std_fastqc))
                remove_file(std_fastqc)

            # get again the html files
            html_files = ["%s/%s"%(fastqc_dir, x) for x in os.listdir(fastqc_dir) if x.endswith(".html")]

            all_html_files.append(html_files[0])


        # get the adapters from the fastqc report
        adapters = set.union(*[get_set_adapter_fastqc_report(x) for x in all_html_files])

        # write adapters to fasta
        all_seqs = []
        existing_ids = set()
        for adapter in adapters:
            ID = id_generator(already_existing_ids=existing_ids); existing_ids.add(ID)
            all_seqs.append(SeqRecord(Seq(adapter), id=ID, name="", description=""))

        adapters_filename = "%s/adapters.fasta"%get_dir(reads1)
        SeqIO.write(all_seqs, adapters_filename, "fasta")

        # run trimmomatic
        if file_is_empty(trimmed_reads1) or file_is_empty(trimmed_reads2) or replace is True:

            std_trimmomatic = "%s.trimmomatic_std.txt"%trimmed_reads1
            print_if_verbose("running trimmomatic. The std is in %s"%std_trimmomatic)
            trim_cmd = "%s --number_threads %i -rr1 %s -rr2 %s -tr1 %s -tr2 %s -ad %s > %s 2>&1"%(TRIMMOMATIC, threads, reads1, reads2, trimmed_reads1, trimmed_reads2, adapters_filename, std_trimmomatic)

            run_cmd(trim_cmd)
            remove_file(std_trimmomatic)

        # check that the reads are correct
        check_that_paired_reads_are_correct(trimmed_reads1, trimmed_reads2)

    return trimmed_reads1, trimmed_reads2

def run_porechop(raw_reads, replace=False, threads=4):

    """This function takes some raw nanopore reads and runs porechop on them to get the trimmed reads, returning them. All files are written under raw_reads """

    # define the trimmed reads
    trimmed_reads = "%s.trimmed.fastq.gz"%raw_reads
    trimmed_reads_tmp = "%s.trimmed.tmp.fastq.gz"%raw_reads

    if file_is_empty(trimmed_reads) or replace is True:

        # delete previous files
        remove_file(trimmed_reads_tmp)

        # run cmd
        porechop_std = "%s.std.txt"%trimmed_reads
        print_if_verbose("running porechop. The std is in %s"%porechop_std)
        run_cmd("%s -i %s -o %s --threads %i > %s 2>&1"%(porechop, raw_reads, trimmed_reads_tmp, threads, porechop_std))
        remove_file(porechop_std)

        # rename
        os.rename(trimmed_reads_tmp, trimmed_reads)


    return trimmed_reads

def get_coverage_for_reads(reads1, reads2, reference_genome, outdir, threads=4, replace=False):

    """This function runs the coverage calculation for sets of reads """

    # get the trimmed reads
    trimmed_reads1, trimmed_reads2 = run_trimmomatic(reads1, reads2, replace=replace, threads=threads)

    # get the aligned reads
    bamfile = "%s/aligned_reads.bam"%outdir
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam

    run_bwa_mem(trimmed_reads1, trimmed_reads2, reference_genome, outdir, bamfile, sorted_bam, index_bam, "nameSample", threads=threads, replace=replace)
      
    # define the parallel running of mosdepth 
    if threads==1: run_in_parallel=False
    else: run_in_parallel = True

    # get the coverage df
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=True), sep="\t")

    # define stats
    mean_coverage = np.mean(coverage_df.mediancov_1)
    fraction_genome_covered = np.mean(coverage_df.percentcovered_1)/100

    return mean_coverage, fraction_genome_covered


def get_sorted_bam_for_untrimmed_reads(reads1, reads2, reference_genome, outdir, threads=4, replace=False):

    """Takes some untrimmed reads and returns the sorted_bam"""

    # get the trimmed reads
    print_if_verbose("running trimmomatic")
    trimmed_reads1, trimmed_reads2 = run_trimmomatic(reads1, reads2, replace=replace, threads=threads)

    # get the aligned reads
    bamfile = "%s/aligned_reads.bam"%outdir
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam
    print_if_verbose("running bwa mem into %s"%sorted_bam)

    run_bwa_mem(trimmed_reads1, trimmed_reads2, reference_genome, outdir, bamfile, sorted_bam, index_bam, "nameSample", threads=threads, replace=replace)

    return sorted_bam



def get_fraction_readPairsMapped(bamfile, replace=False, threads=4):

    """Returns the fraction of reads mappend for a bam file"""

    # get the npairs, which already generates the flagstat
    npairs = count_number_read_pairs(bamfile, replace=replace, threads=threads)

    # get the fraction mapped
    flagstat_file = "%s.flagstat"%bamfile

    if npairs==0: fraction_mapped = 0.0
    else: fraction_mapped = [float(l.split("mapped (")[1].split("%")[0])/100 for l in open(flagstat_file, "r").readlines() if "mapped (" in l][0]

    return fraction_mapped


def downsample_bamfile_keeping_pairs(bamfile, fraction_reads=0.1, replace=True, threads=4, name="sampleX", sampled_bamfile=None):

    """Takes a sorted and indexed bam and samples a fraction_reads randomly. This is a fast process so that it can be repeated many times"""

    # define the outfile
    seed = random.choice(list(range(300)))

    # define sampled_bamfile if not provided
    if sampled_bamfile is None:
        sampled_bamfile = "%s.%ipct_reads_seed%i_%s.bam"%(bamfile, int(fraction_reads*100), seed, name); remove_file(sampled_bamfile)


    # get fraction as str
    fraction_str = str(fraction_reads).split(".")[1]

    # run the sampling
    sampled_bamfile_unedited = "%s.unedited"%sampled_bamfile
    sampled_bamfile_unedited_std = "%s.generating.std"%sampled_bamfile_unedited
    print_if_verbose("running samtools view. The std is in %s"%sampled_bamfile_unedited_std)
    run_cmd("%s view -o %s -s %i.%s --threads %i -b %s > %s 2>&1"%(samtools, sampled_bamfile_unedited, seed, fraction_str, threads, bamfile, sampled_bamfile_unedited_std))
    remove_file(sampled_bamfile_unedited_std)

    # now rename the reads, so that the read name gets a prefix called name
    sampled_bamfile_generating_stderr = "%s.generating.stderr"%sampled_bamfile
    print_if_verbose("getting the sampled bam. The stderr is in %s"%sampled_bamfile_generating_stderr)
    run_cmd("%s view -h --threads %i %s | sed -e 's/^/%s_/' | sed 's/^%s_@/@/' | %s view --threads %i -bSh > %s 2>%s"%(samtools, threads, sampled_bamfile_unedited, name, name, samtools, threads, sampled_bamfile, sampled_bamfile_generating_stderr))
    remove_file(sampled_bamfile_generating_stderr)

    # at the end remove the unedited
    remove_file(sampled_bamfile_unedited)

    return sampled_bamfile


def download_srr_with_prefetch(srr, SRRfile, replace=False):

    """This function downloads an srr file for an srr if not already done"""

    # define the downloading dir
    downloading_dir = get_dir(SRRfile)

    # make the downloading dir
    make_folder(downloading_dir)

    # define the std files
    prefetch_std = "%s.std.txt"%SRRfile
    prefetch_std_copy = "%s.copy"%prefetch_std

    # try a couple of times
    for Itry in range(2):

        if file_is_empty(SRRfile) or replace is True:
            print_if_verbose("running prefetch for %s"%srr)

            # remove the locks of previous runs
            for file in ["%s/%s"%(downloading_dir, f) for f in os.listdir(downloading_dir) if ".lock" in f or ".tmp." in f]: remove_file(file)

            # remove the actual srr
            remove_file(SRRfile)

            # run prefetch
            print_if_verbose("running prefetch. The std can be found in %s"%prefetch_std)
            try: run_cmd("%s -o %s --max-size 500G --progress 1 %s > %s 2>&1"%(prefetch, SRRfile, srr, prefetch_std))
            except: print_if_verbose("prefetch did not work for %s"%srr)

            # test that the std of prefetch states that there are no unresolved dependencies
            std_lines = open(prefetch_std, "r").readlines()
            successful_download = any(["was downloaded successfully" in l for l in std_lines])
            no_dependencies_left = any(["has 0 unresolved dependencies" in l for l in std_lines])
            has_dependencies_line = any(["unresolved dependencies" in l for l in std_lines])

            if not successful_download or (has_dependencies_line and not no_dependencies_left): 
                run_cmd("cp %s %s"%(prefetch_std, prefetch_std_copy))
                print_if_verbose("prefetch did not work for %s. Test the log in %s"%(srr, prefetch_std_copy))
                remove_file(SRRfile)

    # check that the prefetch works 
    if file_is_empty(SRRfile): 
        raise ValueError("prefetch did not work for %s"%srr)

    remove_file(prefetch_std)
    remove_file(prefetch_std_copy)

    return SRRfile


def get_n_pairs_in_fastqgz(file):

    """Takes a fastqgz file and returns the number of reads"""

    # get the number of lines into file
    file_wc = "%s.wc"%file
    file_wc_tmp = "%s.tmp"%file_wc
    if file_is_empty(file_wc):

        unpigz_stderr = "%s.generating.stderr"%file_wc_tmp
        print_if_verbose("calculating # reads for %s. The stderr is in %s"%(file, unpigz_stderr))

        run_cmd("%s -c %s | wc -l > %s 2>%s"%(unpigz, file, file_wc_tmp, unpigz_stderr))
        remove_file(unpigz_stderr)

        os.rename(file_wc_tmp, file_wc)

    # get the number
    nlines = int(open(file_wc, "r").readlines()[0].strip())

    # check that it is multiple of 4
    if nlines%4!=0: raise ValueError("nlines %i is not valid"%nlines)

    # get the number of lines between 4
    return nlines/4

def readIDs_are_correct(readIDs):

    """Takes an iterable of fastq IDs and returns whether they look good"""

    return all([r.startswith("@") for r in readIDs])

def get_last_reads_fastqgz_file(file, nreads=100):

    """Takes a fastq.gz file and returns the read names of the last nreads"""

    # get the last reads
    file_last_reads = "%s.lastReads%iReads"%(file, nreads)
    file_last_reads_tmp = "%s.tmp"%file_last_reads

    # try twice
    for Itry in range(2):

        try:

            if file_is_empty(file_last_reads):
                print_if_verbose("getting last reads for %s"%file)

                # remove the gz index file
                remove_file("%si"%file)
                gztool_stderr = "%s.generating.stderr"%file_last_reads_tmp
                print_if_verbose("running gztool. The stderr is in %s"%gztool_stderr)
                cmd_ztail = "%s -t %s -v 0 | tail -%i | sed -n '1~4p' | cut -f1 -d ' ' > %s 2>%s"%(gztool, file, int(nreads*4), file_last_reads_tmp, gztool_stderr)
                run_cmd(cmd_ztail)
                remove_file(gztool_stderr)

                os.rename(file_last_reads_tmp, file_last_reads)

            # get as list
            last_reads = [l.strip() for l in open(file_last_reads, "r").readlines()]

            if not readIDs_are_correct(last_reads): remove_file(file_last_reads) # remove the file and try again

        except: print_if_verbose("gztool did not work. This may be due to stochastic errors or impossible installation")

    # if the reads are incorrect, get them by another, slower way
    if file_is_empty(file_last_reads) or not readIDs_are_correct(last_reads): 

        # getting last reads in a more traditional way
        tail_zcat_stderr = "%s.generating.stderr"%file_last_reads_tmp
        print_if_verbose("getting last reads for %s in a slower way (zcat + tail). The stderr is in %s"%(file, tail_zcat_stderr))
        run_cmd("zcat %s | tail -%i | sed -n '1~4p' | cut -f1 -d ' ' > %s 2>%s"%(file, int(nreads*4), file_last_reads_tmp, tail_zcat_stderr))

        os.rename(file_last_reads_tmp, file_last_reads)

        last_reads = [l.strip() for l in open(file_last_reads, "r").readlines()]

        if not readIDs_are_correct(last_reads): 

            print_if_verbose("These are the wrong reads:", last_reads)
            raise ValueError("There are some incorrect reads in %s"%file)

    return last_reads

def check_that_paired_reads_are_correct(reads1, reads2):

    """This function takes some paired end reads and throws an error if the last 500 reads do not have exactly the same reads. And also that the read counts are the same"""

    # get the number of reads
    """
    nreads1 = get_n_pairs_in_fastqgz(reads1)
    nreads2 = get_n_pairs_in_fastqgz(reads2)

    if nreads1!=nreads2: raise ValueError("%s and %s have a different number of reads"%(reads1, reads2))
    """

    # get the last reads
    last_reads1 = get_last_reads_fastqgz_file(reads1)
    last_reads2 = get_last_reads_fastqgz_file(reads2)

    if last_reads1!=last_reads2: raise ValueError("%s and %s do not have a proper read pairing"%(reads1, reads2))

def run_parallelFastqDump_on_prefetched_SRRfile(SRRfile, replace=False, threads=4):

    """This function runs parallel fastqdump on a tmp dir in SRR file"""


    # define the outdir and tmp dir
    outdir = get_dir(SRRfile)
    tmpdir = "%s_downloading_tmp"%SRRfile

    # define the final reads
    reads1 = "%s_1.fastq.gz"%SRRfile
    reads2 = "%s_2.fastq.gz"%SRRfile

    print_if_verbose("These are the expected final reads:\n%s\n%s. You can create these files and there will not be any downloading"%(reads1, reads2))

    if any([file_is_empty(f) for f in [reads1, reads2]]) or replace is True:

        # delete previous files
        print_if_verbose(tmpdir)
        for f in [reads1, reads2]: remove_file(f)
        delete_folder(tmpdir); make_folder(tmpdir)

        # test if fastqdump works. Run for 15 seconds into test_dir and check if any files have been generated
        print_if_verbose("testing fastqdump")
        test_dir = "%s/testing"%tmpdir
        delete_folder(test_dir); make_folder(test_dir)

        # run cmd and kill after 30 seconds
        testing_fastqdump_parallel_std = "%s.testing_fastqdumpParallel.std"%SRRfile
        print_if_verbose("running the testing of the parallel fastqdump. The std is in %s"%testing_fastqdump_parallel_std)
        test_cmd = "timeout 30 %s -s %s -t %i -O %s --tmpdir %s --split-3 --gzip > %s 2>&1"%(parallel_fastq_dump, SRRfile, threads, test_dir, test_dir, testing_fastqdump_parallel_std)
        try: run_cmd(test_cmd)
        except: print_if_verbose("fastqdump did not work before timeout. This is fine.")

        any_fastqgz_generated = any([any([file.endswith("fastq.gz") for file in f[2]]) for f in os.walk(test_dir)])

        # define the std of fastdump
        stdfile = "%s/std_fastqdump.txt"%tmpdir

        if any_fastqgz_generated is True:

            # run fastqdump parallel into the tmpdir 
            print_if_verbose("running fastq dump in parallel. The std is in %s"%stdfile)
            run_cmd("%s -s %s -t %i -O %s --tmpdir %s --split-3 --gzip > %s 2>&1 "%(parallel_fastq_dump, SRRfile, threads, tmpdir, tmpdir, stdfile))

            # check that the fastqdump is correct
            std_lines = open(stdfile, "r").readlines()
            written_lines_as_threads = len([l for l in std_lines if l.startswith("Written")])==threads
            any_error = any(["ERROR" in l.upper() for l in std_lines])
            if not written_lines_as_threads or any_error:
                raise ValueError("Something went wrong with the fastqdump. Check the log in %s"%stdfile)

            # define the tmp reads
            tmp_reads1 = "%s/%s"%(tmpdir, get_file(reads1))
            tmp_reads2 = "%s/%s"%(tmpdir, get_file(reads2))

        else: 
            print_if_verbose("running normal fastqdump. This is slower. The std is in %s"%stdfile)

            # run fastqdump
            srr = get_file(SRRfile).split(".")[0]
            run_cmd("%s --split-3 --gzip --outdir %s %s > %s 2>&1"%(fastqdump, tmpdir, srr, stdfile))

            # define the tmp reads
            tmp_reads1 = "%s/%s_1.fastq.gz"%(tmpdir, srr)
            tmp_reads2 = "%s/%s_2.fastq.gz"%(tmpdir, srr)

        # check that the read pairs are correct
        check_that_paired_reads_are_correct(tmp_reads1, tmp_reads2)

        # remove unneccessary logs
        remove_file(testing_fastqdump_parallel_std)
        remove_file(stdfile)


        # rename
        os.rename(tmp_reads1, reads1)
        os.rename(tmp_reads2, reads2)

    # delete the tmpdir
    delete_folder(tmpdir)

    check_that_paired_reads_are_correct(reads1, reads2)

    return reads1, reads2


def run_parallelFastqDump_on_prefetched_SRRfile_nanopore(SRRfile, replace=False, threads=4):

    """This function is equivalent to run_parallelFastqDump_on_prefetched_SRRfile but returning only one fastq. """

    # define the outdir and tmp dir
    outdir = get_dir(SRRfile)
    tmpdir = "%s_downloading_tmp"%SRRfile

    # define the final reads
    reads = "%s.fastq.gz"%SRRfile

    if file_is_empty(reads) or replace is True:

        # delete previous files
        remove_file(reads)
        delete_folder(tmpdir); make_folder(tmpdir)

        # run fastqdump parallel into the tmpdir 
        stdfile = "%s/std_fastqdump.txt"%tmpdir
        print_if_verbose("running fastqdump in parallel for ONT prefetched .srr file. The std is in %s"%stdfile)
        run_cmd("%s -s %s -t %i -O %s --tmpdir %s --gzip > %s 2>&1"%(parallel_fastq_dump, SRRfile, threads, tmpdir, tmpdir, stdfile))

        # check that the fastqdump is correct
        std_lines = open(stdfile, "r").readlines()
        written_lines_as_threads = len([l for l in std_lines if l.startswith("Written")])==threads
        any_error = any(["ERROR" in l.upper() for l in std_lines])
        if not written_lines_as_threads or any_error:
            raise ValueError("Something went wrong with the fastqdump. Check the log in %s"%stdfile)

        remove_file(stdfile)

        # rename the reads
        tmp_reads = "%s/%s"%(tmpdir, get_file(reads))
        os.rename(tmp_reads, reads)

    # delete the tmpdir
    delete_folder(tmpdir)

    return reads

def run_svim(reads, reference_genome, outdir,  threads=4, replace=False, min_sv_size=50, max_sv_size=1000000000000000000000, aligner="ngmlr", is_nanopore=True, minimum_depth=5):

    """Takes some reads and a reference genome and runs svim. The reads should be in fastq.gz"""

    # get the name of the sorted bam
    sorted_bam_long = "%s/%s.%s.coordsorted.bam"%(outdir, get_file(reads).rstrip(".gz"), aligner)
    sorted_bam_long_idx = "%s.bai"%sorted_bam_long

    # change the name so that it is shorter, this is good for making further folders
    sorted_bam_short = "%s/aligned_reads.sorted.bam"%outdir
    sorted_bam_short_idx = "%s.bai"%sorted_bam_short

    #if any([not os.path.isfile(x) for x in svType_to_file.values()]) or file_is_empty(sorted_bam_short) or file_is_empty(sorted_bam_short_idx) or replace is True:

    # run svim
    if file_is_empty(sorted_bam_short) or file_is_empty(sorted_bam_short_idx) or replace is True:
     
        # make the folder
        delete_folder(outdir); make_folder(outdir)

        # define the std
        svim_std = "%s/running_svim.std"%outdir
        print_if_verbose("running svim. The std is in %s"%svim_std)

        # run svim with few filters
        svim_cmd = "%s reads %s %s %s --min_sv_size %i --max_sv_size %i --cores %i --aligner %s --minimum_depth %s --min_mapq 0 --skip_genotyping > %s 2>&1"%(svim, outdir, reads, reference_genome, min_sv_size, max_sv_size, threads, aligner, minimum_depth, svim_std)
        if is_nanopore is True: svim_cmd += " --nanopore"
        run_cmd(svim_cmd)
        remove_file(svim_std)
        
        os.rename(sorted_bam_long, sorted_bam_short)
        os.rename(sorted_bam_long_idx, sorted_bam_short_idx)


    # calculate the coverage
    destination_dir = "%s.calculating_windowcoverage"%sorted_bam_short
    coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam_short, windows_file="none", replace=replace, run_in_parallel=True, delete_bams=True), sep="\t")

    median_coverage = np.median(coverage_df.mediancov_1)
    fraction_genome_covered = np.mean(coverage_df.percentcovered_1)/100

    print_if_verbose("The median coverage is %.2f and the fraction of the genome covered is %.3f"%(median_coverage, fraction_genome_covered))

    # define the outfiles
    svType_to_file = {svtype : "%s/candidates/candidates_%s.corrected.bed"%(outdir, svtype) for svtype in {"breakends", "deletions", "int_duplications_dest", "int_duplications_source", "inversions", "novel_insertions", "tan_duplications_dest", "tan_duplications_source"}}

    if any([not os.path.isfile(x) for x in svType_to_file.values()]) or replace is True:

        #### ADD HEADER TO TABLES ####

        # define the column names
        col3_Bnd_IntDup_TanDup = "svtype;partner_dest;std_pos_across_cluster;std_span_across_cluster"
        col3_Del_Inv_Ins = "svtype;std_pos_across_cluster;std_span_across_cluster"

        svtype_to_col3_name = {"breakends":col3_Bnd_IntDup_TanDup, "deletions":col3_Del_Inv_Ins, "int_duplications_dest":col3_Bnd_IntDup_TanDup, "int_duplications_source":col3_Bnd_IntDup_TanDup, "inversions":col3_Del_Inv_Ins, "novel_insertions":col3_Del_Inv_Ins, "tan_duplications_dest":col3_Bnd_IntDup_TanDup, "tan_duplications_source":col3_Bnd_IntDup_TanDup}

        colnamesDict_InsDelTanInv = {0:"Chr", 1:"Start", 2:"End", 4:"score", 5:"evidence_deleted_origin", 6:"signatures_making_this_candidate"}
        colnamesDict_Bnd = {0:"Chr", 1:"Start", 2:"End", 4:"score", 5:"signatures_making_this_candidate"}

        # rewrite the candidates adding header
        candidates_dir = "%s/candidates"%outdir
        for file in os.listdir(candidates_dir):
            svtype = file.split(".")[0].split("candidates_")[1]

            # define the colnames for this svtype
            if svtype=="breakends": colnames_dict = colnamesDict_Bnd
            else: colnames_dict = colnamesDict_InsDelTanInv
            colnames_dict[3] = svtype_to_col3_name[svtype]

            # get the df
            filename = "%s/%s"%(candidates_dir, file)
            df = pd.read_csv(filename, sep="\t", header=-1).rename(columns=colnames_dict)

            # write
            df.to_csv(svType_to_file[svtype], sep="\t", index=False, header=True)

        ############################

    return svType_to_file, sorted_bam_short, median_coverage

def download_srr_parallelFastqDump(srr, destination_dir, is_paired=True, threads=4, replace=False):

    """Takes an SRR accession and downloads the fastq files into destination_dir with prefetch"""

    isPaired_to_fastqIDXs = {True: [1,2], False: [1]}

    make_folder(destination_dir)

    # define a tmp dir
    downloading_dir= "%s/downloading"%destination_dir # this is temporary and will be empty
    tmp_dir = "%s/tmp"%downloading_dir

    # define the expected files
    expected_files = ["%s/%s_%i.fastq.gz"%(destination_dir, srr, N) for N in isPaired_to_fastqIDXs[is_paired]]

    # download into downloading_dir
    if any([file_is_empty(x) for x in expected_files]) or replace is True: 
        print_if_verbose("downloading %s"%srr)

        # remove the tmp dir if it exists
        delete_folder(tmp_dir)

        # make the dirs newly
        for f in [downloading_dir, tmp_dir]: make_folder(f)

        # first run prefetch if not already done
        SRRfile = "%s/%s.srr"%(downloading_dir, srr)
        download_srr_with_prefetch(srr, SRRfile)

        # download into fastq split files
        std_fqdump_parallel = "%s.running_fqDumpParallel.std"%SRRfile
        print_if_verbose("running fastqdump in parallel. The std is in %s"%std_fqdump_parallel)
        run_cmd("%s -s %s -t %i -O %s --tmpdir %s --split-3 --gzip > %s 2>&1"%(parallel_fastq_dump, SRRfile, threads, downloading_dir, tmp_dir, std_fqdump_parallel))
        remove_file(std_fqdump_parallel)

        # move the fastq files into destination
        for N in isPaired_to_fastqIDXs[is_paired]: os.rename("%s/%s.srr_%i.fastq.gz"%(downloading_dir, srr, N), "%s/%s_%i.fastq.gz"%(destination_dir, srr, N))

    # at the end remove the downloading dir
    delete_folder(downloading_dir)

def get_SRA_runInfo_df_with_mapping_data(SRA_runInfo_df, reference_genome, outdir, replace=False, threads=4, coverage_subset_reads=10):

    """This function takes a df with info about SRA runs and returns it with fraction_reads_mapped, fraction_genome_covered and mean_coverage"""

    make_folder(outdir)

    # define the outdir that will store the seq data
    seq_data_dir = "%s/seq_data"%outdir; make_folder(seq_data_dir)

    # change index
    SRA_runInfo_df = SRA_runInfo_df.set_index("Run", drop=False)

    # calculate the length of the genome
    length_genome = sum(get_chr_to_len(reference_genome).values())

    # add the subset_n_reads depending on the coverage and the read length
    SRA_runInfo_df["subset_n_reads"] = (length_genome*coverage_subset_reads / SRA_runInfo_df.avgLength).apply(int) + 1

    ###### DOWNLOAD THE SUBSET OF READS WITH FASTQDUMP ####

    inputs_download_srr = [(srr, "%s/%s/reads_dir"%(seq_data_dir, srr), SRA_runInfo_df.loc[srr, "subset_n_reads"]) for srr in SRA_runInfo_df.Run]

    print_if_verbose("Downloading %i sra datasets or %i threads"%(len(inputs_download_srr), threads))
    with multiproc.Pool(threads) as pool:
        list_reads_tuples = pool.starmap(download_srr_subsetReads_onlyFastqDump, inputs_download_srr)
        pool.close()
        pool.terminate()

  
    ###################################################################

    ##### GET FRAC READS MAPPED ##### 

    # initialize lists
    fraction_reads_mapped_list = []

    # get the files
    for I, (reads1, reads2) in enumerate(list_reads_tuples):

        # define the outdir
        srr = SRA_runInfo_df.iloc[I]["Run"]
        outdir_srr = "%s/%s"%(seq_data_dir, srr)

        # get the sorted_bam
        sorted_bam = get_sorted_bam_for_untrimmed_reads(reads1, reads2, reference_genome, outdir_srr, threads=threads, replace=replace)

        # get fraction reads mapped
        fraction_mapped = get_fraction_readPairsMapped(sorted_bam)

        #except: fraction_mapped = 0.0

        # get the fraction of reads mapped
        fraction_reads_mapped_list.append(fraction_mapped)

    SRA_runInfo_df["fraction_reads_mapped"] = fraction_reads_mapped_list

    ############################################

    return SRA_runInfo_df

def get_taxID_or_BioSample(r, target_taxID):

    """Takes a row of the all_SRA_runInfo_df and returns the TaxID. If the taxID is the target_taxID it returns the biosample"""

    if r["TaxID"]==target_taxID: return r["BioSample"]
    else: return r["TaxID"]

def get_SRA_runInfo_df(target_taxID, n_close_samples, nruns_per_sample, outdir, reference_genome, min_coverage, replace, threads, coverage_subset_reads, min_fraction_reads_mapped, get_lowest_coverage_possible=False):

    """This function mines the SRA to find n_close_samples and nruns_per_sample, returning the necessary df """

    ######## UPDATE NCBI TAXONOMY ########

    print_if_verbose("Getting genomes for taxID into %s"%(outdir))

    # load the NCBI taxonomy database and upgrade it if not already done
    print_if_verbose("getting NCBI taxonomy database")

    ncbi = NCBITaxa()
    ncbiTaxa_updated_file = "%s/ncbiTaxa_updated.txt"%outdir
    if file_is_empty(ncbiTaxa_updated_file) or replace is True: 

        # update
        ncbi.update_taxonomy_database()

        # write file
        open(ncbiTaxa_updated_file, "w").write("NCBItaxa updated\n")

    #######################################

    # get the outdir were to store the IDs 
    outdir_gettingID = "%s/getting_sample_IDs"%outdir; make_folder(outdir_gettingID)

    total_nruns = n_close_samples*nruns_per_sample
    print_if_verbose("Looking for %i runs"%total_nruns)

    # initialize a set that defines the runs of the previous node
    runs_previous_nodes = set()

    # initialize the final df
    final_SRA_runInfo_df = pd.DataFrame()

    # define all potentially interesting taxIDs close to the target_taxIDs
    for nancestorNodes in range(1, 100): # one would mean to consider only IDs that are under the current species
        print_if_verbose("Considering %i ancestor nodes"%nancestorNodes)

        # initialize a df that stores the tested SRRs for this division
        df_division_tested = pd.DataFrame()

        # create a folder for this number of ancestors
        outdir_ancestors = "%s/all_runsWithWGS_arround_target_taxID_%i_considering%iAncestors"%(outdir, target_taxID, nancestorNodes); make_folder(outdir_ancestors)

        # get the ancestor
        ancestor_taxID = int(ncbi.get_lineage(target_taxID)[-nancestorNodes])

        # get the runs for this division excluding the target taxID
        print_if_verbose("Getting WGS info")
        fileprefix = "%s/output"%(outdir_ancestors)
        all_SRA_runInfo_df = get_allWGS_runInfo_fromSRA_forDivision(fileprefix, ancestor_taxID, reference_genome, taxIDs_to_exclude={target_taxID}, replace=False, min_coverage=min_coverage).set_index("Run", drop=False)

        # exclude the wrong SRRs
        all_SRA_runInfo_df = all_SRA_runInfo_df[~all_SRA_runInfo_df.Run.isin(wrong_SRRs)]

        # if it is empty, continue
        if len(all_SRA_runInfo_df)==0: continue

        # define the runs with target taxID
        runs_target_taxID = set(all_SRA_runInfo_df[all_SRA_runInfo_df.TaxID==target_taxID].Run)

        # get the tree
        tree = ncbi.get_descendant_taxa(ancestor_taxID, collapse_subspecies=False, return_tree=True, intermediate_nodes=True)
        #print("You can visualize the tree in http://etetoolkit.org/treeview/: \n\n%s\n\n"%(tree.write()))

        if type(tree)==list:
            interesting_taxIDs = set(tree)
            taxID_to_distanceToTarget = {tree[0] : 0}

        else:

            # define interesting taxIDs (the leafs and the species names that may be intermediate)
            interesting_taxIDs = {int(x) for x in set(tree.get_leaf_names()).union({n.name for n in tree.traverse() if n.rank=="species"})}.difference({target_taxID})

            # map the distance between each leave and the target
            taxID_to_distanceToTarget = {taxID : tree.get_distance(str(target_taxID), str(taxID)) for taxID in interesting_taxIDs}

        # map names
        taxID_to_sciName = ncbi.get_taxid_translator(interesting_taxIDs)

        # add the sciName
        all_SRA_runInfo_df["sci_name"] = all_SRA_runInfo_df.TaxID.map(taxID_to_sciName)

        # get the taxIDs sorted by the distance (so that the closest )
        interesting_taxIDs_sorted = sorted(interesting_taxIDs, key=(lambda x: taxID_to_distanceToTarget[x]))

        # add the number of runs that each taxID has 
        taxID_to_nRuns = Counter(all_SRA_runInfo_df.TaxID)
        all_SRA_runInfo_df["nRuns_with_taxID"] = all_SRA_runInfo_df.TaxID.apply(lambda x: taxID_to_nRuns[x])
        all_SRA_runInfo_df = all_SRA_runInfo_df[all_SRA_runInfo_df.nRuns_with_taxID>=nruns_per_sample]

        # if you did not find anything get to farther ancestors
        if len(set(all_SRA_runInfo_df.TaxID).intersection(interesting_taxIDs))<n_close_samples: continue

        # iterate through each taxID from closest to farthest
        for Itax, taxID in enumerate(interesting_taxIDs_sorted):

            # get the df for this taxID
            df_taxID = all_SRA_runInfo_df[all_SRA_runInfo_df.TaxID==taxID].sort_values(by="expected_coverage", ascending=get_lowest_coverage_possible)
            if len(df_taxID)==0: continue

            # define the length of the final_SRA_runInfo_df taxIDs
            if len(final_SRA_runInfo_df)==0: nFinal_taxIDs = 0
            else: nFinal_taxIDs = len(set(final_SRA_runInfo_df.TaxID))

            print_if_verbose("working on taxon %s. It is at %i nodes of the target. It is taxID %i/%i. final_SRA_runInfo_df has %i/%i taxIDs"%(taxID_to_sciName[taxID], taxID_to_distanceToTarget[taxID], Itax+1, len(interesting_taxIDs_sorted), nFinal_taxIDs, n_close_samples))

            # go through several slices of the taxID df
            for end_indx in range(nruns_per_sample, len(df_taxID)+1):

                # get the df of a slice
                df_slice = df_taxID.iloc[0:end_indx]

                # add the mapping data 
                df_slice = get_SRA_runInfo_df_with_mapping_data(df_slice, reference_genome, outdir_gettingID, replace=replace, threads=threads, coverage_subset_reads=coverage_subset_reads)

                # keep
                df_division_tested = df_division_tested.append(df_slice)

                # filter            
                idx = ((df_slice.fraction_reads_mapped>=min_fraction_reads_mapped))

                print_if_verbose(df_slice[["fraction_reads_mapped"]])
                print_if_verbose("There are %i/%i files that passed the filters"%(sum(idx), len(idx)), min_fraction_reads_mapped)

                df_slice = df_slice[idx]

                # if there is something left keep it
                if len(df_slice)>=nruns_per_sample: 

                    # keep to the final df
                    final_SRA_runInfo_df = final_SRA_runInfo_df.append(df_slice.iloc[0:nruns_per_sample])

                    # break the iteration through this taxID
                    break

                # if none of the df_slice are good and you have tried at least 3 times, break
                if len(df_slice)==0 and len(idx)>=3: break

                # if there are more than 5 idxs and more than half are bad, break
                if len(idx)>=6 and sum(idx)/len(idx)<0.5: break

            # once you already got enough data, break
            if len(final_SRA_runInfo_df)==total_nruns: break

        # if there are no new nodes, break
        runs_in_this_node = set(final_SRA_runInfo_df.Run)
        if len(runs_in_this_node.difference(runs_previous_nodes))==0: break
        runs_previous_nodes.update(runs_in_this_node)

        # if more than 90% of the samples tested in this division are below 
        fraction_runs_tested_correct = sum(df_division_tested.fraction_reads_mapped>=min_fraction_reads_mapped)/len(df_division_tested)
        print_if_verbose("There are %.3f of runs analyzed that map correctly"%fraction_runs_tested_correct)
        if fraction_runs_tested_correct<0.1: break


        # if you already found the IDs, break
        if len(final_SRA_runInfo_df)==total_nruns: break

    # debug
    if len(final_SRA_runInfo_df)!=total_nruns: raise ValueError("You could not find any datasets in SRA that would be useful")

    # change the names
    final_SRA_runInfo_df["sampleID"] = final_SRA_runInfo_df.TaxID

    return final_SRA_runInfo_df

def connected_to_network(test_host='http://google.com'):

    """Returns whether you are connected to the network"""

    try: 

        urllib.request.urlopen(test_host)
        return True

    except: return False

def close_shortReads_table_is_correct(close_shortReads_table):

    """Takes a file that is the close_shortReads_table and returns whether all the reads exist"""

    # if it is empty return false
    if file_is_empty(close_shortReads_table): return False

    # load as df
    close_shortReads_table_df = pd.read_csv(close_shortReads_table, sep="\t")

    # check that all the reads exist
    reads_files = set(close_shortReads_table_df["short_reads1"]).union(close_shortReads_table_df["short_reads2"])

    if any([file_is_empty(f) for f in reads_files]): return False
    else: return True

def get_close_shortReads_table_close_to_taxID(target_taxID, reference_genome, outdir, ploidy, n_close_samples=3, nruns_per_sample=3, replace=False, threads=4, min_fraction_reads_mapped=0.1, coverage_subset_reads=0.1, min_coverage=30, job_array_mode="local", StopAfter_sampleIndexingFromSRA=False, queue_jobs="debug", max_ncores_queue=768, time_read_obtention="02:00:00", StopAfterPrefecth_of_reads=False, get_lowest_coverage_possible=False):

    """
    This function takes a taxID and returns the close_shortReads_table that is required to do optimisation of parameters
    """

    ##### GET SRA_runInfo_df##### 

    # define the total number of runs that you need
    SRA_runInfo_df_file = "%s/final_SRA_runInfo_df.py"%outdir

    if file_is_empty(SRA_runInfo_df_file) or replace is True:
        print_if_verbose("getting SRRs")

        # get the df
        SRA_runInfo_df = get_SRA_runInfo_df(target_taxID, n_close_samples, nruns_per_sample, outdir, reference_genome, min_coverage, replace, threads, coverage_subset_reads, min_fraction_reads_mapped, get_lowest_coverage_possible=get_lowest_coverage_possible)

        # load df
        save_object(SRA_runInfo_df, SRA_runInfo_df_file)

    else: SRA_runInfo_df = load_object(SRA_runInfo_df_file)

    #print_if_verbose("these are the samples chosen:\n:", SRA_runInfo_df[["Run", "sampleID", "SampleName", "sci_name", "fraction_reads_mapped"]].sort_values("sampleID"))

    if StopAfter_sampleIndexingFromSRA is True: 
        print_if_verbose("stopping after generation of SRA_runInfo_df into %s"%SRA_runInfo_df_file)
        exit(0)

    ##############################

    ###### GETTING THE FINAL DATASETS ######

    # define the path to the final table
    close_shortReads_table = "%s/close_shortReads_table.tbl"%outdir

    # if it exists 

    # define the dir were you will store the reads
    reads_dir = "%s/reads"%outdir; make_folder(reads_dir)

    if not close_shortReads_table_is_correct(close_shortReads_table) or replace is True:
        print_if_verbose("getting short reads table")

        # change things
        SRA_runInfo_df["sampleID"] = "sample" + SRA_runInfo_df.sampleID.apply(str)
        SRA_runInfo_df["runID"] = SRA_runInfo_df.sampleID + "_" + SRA_runInfo_df.Run

        # initialize the final dict
        srr_to_readsDict = {}

        # initialize the cmds
        all_cmds = []

        for srr in SRA_runInfo_df.Run:
            print_if_verbose("trimming %s reads"%srr)

            # make the outdir for this
            outdir_srr = "%s/%s"%(reads_dir, srr); make_folder(outdir_srr)

            # define files
            trimmed_reads1 = "%s/%s_trimmed_reads_1.fastq.gz"%(outdir_srr, srr)
            trimmed_reads2 = "%s/%s_trimmed_reads_2.fastq.gz"%(outdir_srr, srr)
            SRRfile = "%s/%s.srr"%(outdir_srr, srr)

            # define the finalisation file
            if StopAfterPrefecth_of_reads is True: final_file = SRRfile
            else: final_file = trimmed_reads2 # the last one generated

            # get the cmd if necessary
            if file_is_empty(final_file) and file_is_empty(trimmed_reads2) or replace is True:

                # define the cmd and add it
                cmd = "%s --srr %s --outdir %s --threads %i"%(get_trimmed_reads_for_srr_py, srr, outdir_srr, threads)
                if StopAfterPrefecth_of_reads is True: cmd += " --stop_after_prefetch"

                all_cmds.append(cmd)

            # keep the files
            srr_to_readsDict[srr] = {"short_reads1":trimmed_reads1, "short_reads2":trimmed_reads2}

        # if there are all_cmds, run them in a job array
        if len(all_cmds)>0:
            print_if_verbose("Getting reads into %s"%(reads_dir))

            # if you are in local, run them in parallel. Each job in one threads
            if job_array_mode=="local" and StopAfterPrefecth_of_reads is True: 

                print_if_verbose("running prefetch in parallel")
                with multiproc.Pool(multiproc.cpu_count()) as pool:
                    pool.starmap(run_cmd, [(x,) for x in all_cmds])
                    pool.close()

            # local
            elif job_array_mode=="local" and StopAfterPrefecth_of_reads is False: 

                print_if_verbose("running all files one after the other")
                for cmd in all_cmds: run_cmd(cmd)

            # if you are in a greasy environment
            elif job_array_mode=="greasy": 

                read_downloading_dir = "%s/read_downloading_files"%outdir; make_folder(read_downloading_dir)
                print_if_verbose("Submitting %i downloading reads jobs to greasy. All files will be stored in %s"%(len(all_cmds), read_downloading_dir))

                jobs_filename = "%s/jobs.getting_SRAdatasets"%read_downloading_dir
                open(jobs_filename, "w").write("\n".join(all_cmds))

                generate_jobarray_file_greasy(jobs_filename, walltime=time_read_obtention,  name="getting_SRAdatasets", queue=queue_jobs, sbatch=True, ncores_per_task=threads, number_tasks_to_run_at_once="all", max_ncores_queue=max_ncores_queue)

                # exit before it starts
                print_if_verbose("You need to wait until the read obtention is is finsihed")
                sys.exit(0)

            else: raise ValueError("%s is not valid"%job_array_mode)




        if StopAfterPrefecth_of_reads is True: 
            print_if_verbose("Stopping after prefetch of the reads. You still need to re-run this pipeline to get the actual reads.")
            sys.exit(0)

        # add to the df
        for f in ["short_reads1", "short_reads2"]: SRA_runInfo_df[f] = SRA_runInfo_df.Run.apply(lambda srr: srr_to_readsDict[srr][f])

        # get the final df
        final_df = SRA_runInfo_df[["sampleID", "runID", "short_reads1", "short_reads2"]]

        # write
        final_df.to_csv(close_shortReads_table, sep="\t", header=True, index=False)

    close_shortReads_table_df = pd.read_csv(close_shortReads_table, sep="\t")

    if StopAfterPrefecth_of_reads is True: 
        print_if_verbose("Stopping after prefetch of the reads. You still need to re-run this pipeline to get the actual reads.")
        sys.exit(0)


    # debug
    if not close_shortReads_table_is_correct(close_shortReads_table): raise ValueError("%s has empty reads files"%close_shortReads_table)
    
    # define the important files
    important_files = set(close_shortReads_table_df["short_reads1"]).union(close_shortReads_table_df["short_reads2"])

    # remove the files that are not the reads
    for srr in os.listdir(reads_dir):
        srr_dir = "%s/%s"%(reads_dir, srr)
        for file in os.listdir(srr_dir):
            filepath = "%s/%s"%(srr_dir, file)

            if filepath not in important_files: remove_file(filepath)

    # remove all the other unimportant files and folders
    for f in os.listdir(outdir):

        if f not in {"close_shortReads_table.tbl", "final_SRA_runInfo_df.py", "reads"}:

            path = "%s/%s"%(outdir, f)
            delete_folder(path)
            remove_file(path)

    # last check
    if not close_shortReads_table_is_correct(close_shortReads_table): raise ValueError("%s has empty reads files"%close_shortReads_table)
 
    #########################################

    return close_shortReads_table

def generate_tables_of_SV_between_genomes_gridssClove(query_genome, reference_genome, replace=False, threads=4, coverage=30, insert_size=500, read_lengths=[kb*1000 for kb in [0.3, 0.5, 1, 1.5, 2, 2.5]], error_rate=0.0, gridss_min_af=0.25):

    """Takes a bam file with aligned reads or genomes and generates calls, returning a dict that maps variation type to variants
    - aligner can be minimap2 or ngmlr.

    [0.5, 0.7, 0.9, 1]

    """


    # first run svim under the outdir of the aligned reads
    working_dir = "%s/findingSVlongReadsSimulation_ouptut_%s_against_%s"%(get_dir(query_genome), get_file(query_genome), get_file(reference_genome))
    print_if_verbose("generating svtables into %s"%working_dir)
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
        print_if_verbose("simulating long paired reads")

        # first simulate reads in parallel for each chromosome, as if they were paired
        inputs_fn = [(chr_obj, coverage, insert_size, read_lengths) for chr_obj in SeqIO.parse(query_genome, "fasta")]

        with multiproc.Pool(multiproc.cpu_count()) as pool:
            all_reads_list_tuples = pool.starmap(simulate_pairedEndReads_per_chromosome_uniform, inputs_fn)
            pool.close()

        #all_reads_list_tuples = list(map(lambda x: simulate_pairedEndReads_per_chromosome_uniform(x[0], x[1], x[2]), inputs_fn))
        reads_objects_1 = make_flat_listOflists([r[0] for r in all_reads_list_tuples])
        reads_objects_2 = make_flat_listOflists([r[1] for r in all_reads_list_tuples])
        print_if_verbose("There are %i reads"%len(reads_objects_1))

        # get the reads writen
        print_if_verbose("writing reads")
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
    plot_coverage_across_genome_pairedEndReads(sorted_bam, reference_genome, replace=replace)

    ############################################################################
    
    # define the gridss filters according to the freq, which is related to the expected ploidy
    gridss_filters_dict = default_filtersDict_gridss
    gridss_filters_dict["min_af"] = gridss_min_af
    print_if_verbose("Filtering out when any AF is below %.3f"%(gridss_filters_dict["min_af"]))

    # define the median coverage per region
    outdir_coverage_calculation = "%s/coverage_per_regions"%working_dir; make_folder(outdir_coverage_calculation)
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_coverage_calculation, sorted_bam, windows_file="none", replace=replace), sep="\t")
    median_coverage = np.median(coverage_df.mediancov_1)
    print_if_verbose("The median coverage is %i"%median_coverage)

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
    if equal_fields_match is False: return False

    # ask if the approximate fields match
    approximate_fields_match = all([abs(rp[f]-rk[f])<=tol_bp for f in approximate_fields])
    if approximate_fields_match is False: return False

    # if the rknown is all 0s, return False
    if "POS" in rk.keys():
        if all([rk[f]==0 for f in ["POS", "START", "END"]]): return False

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
        print_if_verbose("WARNING: There are some predictedIDs that match more than one variant")
        #print_if_verbose(df_known)

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

def get_represenative_filtersDict_for_filtersDict_list(filtersDict_list, type_filters="less_conservative"):

    """Takes a lis, each position with a list of filters like passed to get_tupleBreakpoints_for_filters_GRIDSS and returns a representative dict, according to less_conservative"""

    # map a score for each dict
    score_to_dict = {sum([g_filterName_to_filterValue_to_Number[fname][find_nearest(g_filterName_to_filtersList[fname], fvalue)] for fname, fvalue in filtersDict.items()]) : filtersDict for filtersDict in filtersDict_list}

    # get the dict with the min or max score, depedning on the approach
    type_filters_to_getPositionFun = {"less_conservative":min, "most_conservative":max}

    return score_to_dict[type_filters_to_getPositionFun[type_filters](score_to_dict.keys())]

def get_changing_fields_in_df_benchmark(df):

    """This function takes a df such as the input of get_best_less_conservative_row_df_benchmark and returns a set with the keys that are different across rows"""

    changing_fields = set()

    # go through each field
    for f in df.columns:

        # filters dict are treated specially
        if f=="filters_dict": 
            if len(set(df[f].apply(get_dict_as_tuple)))>1: changing_fields.add(f)

        elif len(set(df[f]))>1: changing_fields.add(f)

    return changing_fields

def get_best_less_conservative_row_df_benchmark(df_benchmark):

    """Takes a df_benchmark, and returns the row with the less conservative row, given that it has the highest Fvalue. The least conservative is made in a step wise way filtering several things one after the other"""

    # get the maximum df
    df_best = df_benchmark[df_benchmark.Fvalue==max(df_benchmark.Fvalue)]
    if len(df_best)==1: return df_best.iloc[0]
    
    # get the df with the highest precision
    df_best = df_best[df_best.precision==max(df_best.precision)]
    if len(df_best)==1: return df_best.iloc[0]

    # get the one with the highest recall
    df_best = df_best[df_best.recall==max(df_best.recall)]
    if len(df_best)==1: return df_best.iloc[0]

    # get the least conservative set of filters for gridss
    less_conservative_filtersDict_tuple = get_dict_as_tuple(get_represenative_filtersDict_for_filtersDict_list(df_best.filters_dict, type_filters="less_conservative"))
    df_best = df_best[df_best.filters_dict.apply(get_dict_as_tuple)==less_conservative_filtersDict_tuple]
    if len(df_best)==1: return df_best.iloc[0]

    # get the maximum clove_max_rel_coverage_to_consider_del
    df_best = df_best[df_best.clove_max_rel_coverage_to_consider_del==max(df_best.clove_max_rel_coverage_to_consider_del)]
    if len(df_best)==1: return df_best.iloc[0]

    # get the minimum clove_min_rel_coverage_to_consider_dup
    df_best = df_best[df_best.clove_min_rel_coverage_to_consider_dup==min(df_best.clove_min_rel_coverage_to_consider_dup)]
    if len(df_best)==1: return df_best.iloc[0]

    # get filters with maximum tolerated cov
    if "gridss_maxcoverage" in df_best.keys():
        df_best = df_best[df_best.gridss_maxcoverage==max(df_best.gridss_maxcoverage)]
        if len(df_best)==1: return df_best.iloc[0]

    # if any, take the ones without filtering any regions in gridss
    if "gridss_regionsToIgnoreBed" in df_best.keys():
        
        if any(df_best.gridss_regionsToIgnoreBed==""):    
            df_best = df_best[df_best.gridss_regionsToIgnoreBed==""]
            if len(df_best)==1: return df_best.iloc[0]

    # get the one with the maxiumum threshold defining unbalanced translocations
    df_best = df_best[df_best.threshold_p_unbalTRA==min(df_best.threshold_p_unbalTRA)]
    if len(df_best)==1: return df_best.iloc[0]

    # at the end just return the best one
    print_if_verbose("Warning: there is no single type of filtering that can fullfill all the requirements") 


    # if you didn't find a single best, raise error
    print_if_verbose("\nthis is the best df:\n", df_best, "printing the non equal fields across all rows:\n")
    changing_fields = get_changing_fields_in_df_benchmark(df_best)
    for f in changing_fields:
        print_if_verbose("\t", f)
        for Irow in range(len(df_best)): print_if_verbose("\t\t", df_best[f].iloc[Irow])


    raise ValueError("There is not a single best filtering")

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
    print_if_verbose("Warning: there is no single type of filtering that can fullfill all the requirements") 

    # if you didn't find a single best, raise error
    print_if_verbose("\nthis is the best df:\n", df_best, "printing the non equal fields across all rows:\n")
    changing_fields = get_changing_fields_in_df_benchmark(df_best)
    for f in changing_fields:
        print_if_verbose("\t", f)
        for Irow in range(len(df_best)): print_if_verbose("\t\t", df_best[f].iloc[Irow])


    raise ValueError("There is not a single best filtering")

all_svs = {'translocations', 'insertions', 'deletions', 'inversions', 'tandemDuplications', 'remaining'}
def get_integrated_benchmarking_fields_series_for_setFilters_df(df):

    """This function takes a grouped per-filterSet df and returns a row with all the integrated accuracy measurements. The filters of gridss that are best for each SV may vary. If so we will take the less conservative filters of all of the fiters that are best for each SV."""

    # get a df where each row is one df
    df_best_filters = df.groupby("svtype").apply(get_best_less_conservative_row_df_benchmark)

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

        integrated_benchmarking_results_dict["filters_dict"] = get_represenative_filtersDict_for_filtersDict_list(list(df_best_filters["filters_dict"]), type_filters="less_conservative")
        integrated_benchmarking_results_dict["clove_max_rel_coverage_to_consider_del"] = df_best_filters.loc["deletions", "clove_max_rel_coverage_to_consider_del"]
        integrated_benchmarking_results_dict["clove_min_rel_coverage_to_consider_dup"] = df_best_filters.loc["tandemDuplications", "clove_min_rel_coverage_to_consider_dup"]

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
        print_if_verbose("benchmarking %s"%svtype)

        # load dataframes
        if svtype in svtype_to_predsvfile.keys(): 
            df_predicted = pd.read_csv(svtype_to_predsvfile[svtype], sep="\t")
        else: df_predicted = pd.DataFrame()

        df_known = pd.read_csv(know_SV_dict[svtype], sep="\t")

        # define the fields to find overlaps
        equal_fields = svtype_to_fieldsDict[svtype]["equal_fields"]
        approximate_fields = svtype_to_fieldsDict[svtype]["approximate_fields"]
        chromField_to_posFields = svtype_to_fieldsDict[svtype]["chromField_to_posFields"]

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


def plot_boxplots_allele_freqs(sampleID_to_svtype_to_svDF, filename):

    """This function takes a sampleID_to_svtype_to_svDF and draws a boxplot were the x is the svtype and the y is the allele frequency, for different allele frequencies as hue. Different plots are several IDs """

    print_if_verbose("getting boxplot allele frequencies for sampleID_to_svtype_to_svDF")

    sampleID_to_svtype_to_svDF = cp.deepcopy(sampleID_to_svtype_to_svDF)

    # make two subplots, one for each type of data
    fig = plt.figure(figsize=(10, 4*len(sampleID_to_svtype_to_svDF)))

    # go through each ID
    for I, (ID, svtype_to_svDF) in enumerate(sampleID_to_svtype_to_svDF.items()):

        # sepparate insertions into copy-paste and cut-paste
        svDF_insertions = svtype_to_svDF["insertions"]
        svtype_to_svDF["insertions_copy"] = svDF_insertions[svDF_insertions.Copied]
        svtype_to_svDF["insertions_cut"] = svDF_insertions[~svDF_insertions.Copied]
        del svtype_to_svDF["insertions"]

        # initialize subpplot
        ax =  plt.subplot(len(sampleID_to_svtype_to_svDF), 1, I+1)

        # initialize a df that will contain the allele freqs and the svtype
        df = pd.DataFrame()

        for svtype, svDF in svtype_to_svDF.items():
            for afEstimate in ["estimate_AF_min", "estimate_AF_max", "estimate_AF_mean"]:

                df_af = svDF[[afEstimate]].rename(columns={afEstimate:"af"})

                df_af["svtype"] = svtype
                df_af["af_estimate"] = afEstimate

                df = df.append(df_af, sort=False)

        # get boxplot
        #bp = sns.boxplot(x="svtype", y="af", data=df, hue="af_estimate", notch=True, boxprops=dict(alpha=.99), linewidth=0.5)
        jit = sns.swarmplot(x="svtype", y="af", data=df, hue="af_estimate", dodge=True, size=5, edgecolor="black", linewidth=0.5)

        ax.set_title(ID)

        # add hlines
        for y in [0.25, 0.5, 0.75, 0.9, 1]: plt.axhline(y, color="k", linewidth=.2, linestyle="--")

    # get figure
    fig.tight_layout() 
    fig.savefig(filename, bbox_inches='tight');
    plt.close(fig)

def test_SVgeneration_from_DefaultParms(reference_genome, outdir, sample_sorted_bam, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", nvars=100):

    """This function reports how well the finding of SV from reads works from random simulations. Writing under outdir"""

    # define the output
    precision_and_recall_filename = "%s/precision_and_recall_SVgeneration_from_reads.pdf"%(outdir)
    allele_frequency_boxplots_filename = "%s/allele_frequency_boxplots_SVgeneration_from_reads.pdf"%(outdir)
    if file_is_empty(precision_and_recall_filename) or file_is_empty(allele_frequency_boxplots_filename) or replace is True:

        # initialize the start time
        pipeline_start_time = time.time()

        # prepare files
        make_folder(outdir)

        print_if_verbose("WORKING ON THE VALIDATION THAT WE CAN FIND READS IN AN ASSEMBLY")

        # initialize a df that will contain the benchmarking
        all_df_benchmark_longReads = pd.DataFrame()

        # iniialize dicts
        sampleID_to_svtype_to_file = {}
        sampleID_to_dfGRIDSS = {}

        # go through each simulation
        for simID in range(n_simulated_genomes):
            print_if_verbose("working on simulation %i"%simID)

            # define outdir 
            outdir_sim = "%s/simulation_%i"%(outdir, simID); make_folder(outdir_sim)

            # generate genome with simulated SVs
            sim_svtype_to_svfile, rearranged_genome = rearrange_genomes_simulateSV(reference_genome, outdir_sim, replace=replace, nvars=nvars, mitochondrial_chromosome=mitochondrial_chromosome)

            # get the variants by simulating short reads from the genome
            print_if_verbose("getting SVs from reads")

            # define properties of the run
            chr_to_len = get_chr_to_len(reference_genome)
            median_insert_size, median_insert_size_sd  = get_insert_size_distribution(sample_sorted_bam, replace=replace, threads=threads)
            read_length = get_read_length(sample_sorted_bam, threads=threads, replace=replace)
            total_nread_pairs = count_number_read_pairs(sample_sorted_bam, replace=replace, threads=threads)
            expected_coverage_per_bp = int((total_nread_pairs*read_length) / sum(chr_to_len.values())) +  1 

            # define the function that gets coverage from seq properties
            distToTel_chrom_GC_to_coverage_fn = (lambda x,y,z: expected_coverage_per_bp)

            # get the info of the reference genome with predictions of coverage per window
            df_genome_info = get_windows_infoDF_with_predictedFromFeatures_coverage(rearranged_genome, distToTel_chrom_GC_to_coverage_fn, expected_coverage_per_bp, replace=replace, threads=threads)

            # simulate reads and align them to the reference genome
            outdir_simulation_short_reads = "%s/simulation_shortReads"%(outdir_sim); make_folder(outdir_simulation_short_reads)
            simulated_bam_file = simulate_and_align_PairedReads_perWindow(df_genome_info, rearranged_genome, reference_genome, total_nread_pairs, read_length, outdir_simulation_short_reads, median_insert_size, median_insert_size_sd, replace=replace, threads=threads)

            # call GRIDSS and CLOVE for the simulated reads
            final_run_dir = "%s/final_run_dir"%(outdir_simulation_short_reads); make_folder(final_run_dir)

            predicted_svtype_to_svfile, df_gridss = run_GridssClove_optimising_parameters(simulated_bam_file, reference_genome, final_run_dir, threads=threads, replace=replace, mitochondrial_chromosome=mitochondrial_chromosome, fast_SVcalling=True)

            # get a df of benchmarking
            fileprefix = "%s/rearranged_genome_benchmarking_SV"%outdir_sim
            df_benchmark_longReads = benchmark_processedSVs_against_knownSVs_inHouse(predicted_svtype_to_svfile, sim_svtype_to_svfile, fileprefix, replace=replace)

            # keep
            df_benchmark_longReads["simID"] = [simID]*len(df_benchmark_longReads)
            all_df_benchmark_longReads = all_df_benchmark_longReads.append(df_benchmark_longReads)

            # keep the gridss df and files
            sampleID_to_svtype_to_file[simID] = predicted_svtype_to_svfile
            sampleID_to_dfGRIDSS[simID] = df_gridss
 
        # plot the benchmarking
        plot_bars_single_df_benchmark(all_df_benchmark_longReads, precision_and_recall_filename)

        # get the sampleID_to_svtype_to_svDF
        sampleID_to_svtype_to_svDF = get_sampleID_to_svtype_to_svDF_filtered(sampleID_to_svtype_to_file, sampleID_to_dfGRIDSS)

        # get the boxplots of the allele frequencies
        plot_boxplots_allele_freqs(sampleID_to_svtype_to_svDF, allele_frequency_boxplots_filename)

        # at the end clean the generation
        clean_reference_genome_windows_files(reference_genome)

        print_if_verbose("--- the testing of SV generation from an assembly took %s seconds in %i cores ---"%(time.time() - pipeline_start_time, threads))


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
    print_if_verbose("getting overlaps in breakpoints between samples")
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
    print_if_verbose("There are %i of %i breakpoints already in the parents"%(len(eventIDs_already_in_parents), len(all_eventIDs)))

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

def get_sampleID_to_svtype_to_svDF_filtered(sampleID_to_svtype_to_file, sampleID_to_dfGRIDSS, sampleID_to_parentIDs={}, breakend_info_to_keep=['#CHROM', 'POS', 'other_coordinates', 'allele_frequency', 'allele_frequency_SmallEvent', 'real_AF', 'FILTER', 'inserted_sequence', 'has_poly16GC', 'length_inexactHomology', 'length_microHomology', 'QUAL', 'overlaps_repeats']):

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
        print_if_verbose(sample)

        # get the gridss df
        df_gridss = sampleID_to_dfGRIDSS[sample]

        # add the eventID, as formated in the SVdict
        df_gridss["eventID"] = df_gridss.INFO_EVENT.apply(lambda x: x+"o")

        # define the events that are already in the parents
        all_breakpoints = set(df_gridss["eventID"]).union({""})
        breakpoints_already_in_parents = sampleID_to_eventIDs_alreadyInParents[sample]
        print_if_verbose("There are %i of %i breakpoints already in the parents"%(len(breakpoints_already_in_parents), len(all_breakpoints)))

        # go through each svtype
        for svtype, file in svtype_to_file.items():

            if file!="":

                # get as df
                df_sv = pd.read_csv(file, sep="\t")

                # get for empty df
                if "ID" not in df_sv.keys(): df_sv["ID"] = [""]*len(df_sv)

                # get only the svs that are not in the parents
                df_sv["IDs_set"] = df_sv.ID.apply(lambda x: set([ID[0:-1]+"o" for ID in re.split("\+|\-", x)]))

                if len(df_sv)>0:

                    # check that all of them are in the breakpoints called by gridss
                    all_breakpoints_in_sv = set.union(*df_sv.IDs_set)
                    breakpoints_not_in_df_gridss = all_breakpoints_in_sv.difference(all_breakpoints)
                    if len(breakpoints_not_in_df_gridss)>0: 
                        print_if_verbose("These are breakpoints not in the df_gridss: %s"%(breakpoints_not_in_df_gridss))
                        raise ValueError("There are breakpoints not in the df_gridss, suggesting some errors, such as the fact that you are loading an incorrect df_gridss")

                    # keep only the df with IDs that are not in the parents
                    df_sv = df_sv[df_sv.IDs_set.apply(lambda ids: len(ids.intersection(breakpoints_already_in_parents))==0)]
                    if len(df_sv)==0: 
                        print_if_verbose("There are no new SVs in %s, %s"%(sample, svtype))
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
        print_if_verbose("getting svID for %s"%svtype)

        # define the fields to filter on
        equal_fields = svtype_to_fieldsDict[svtype]["equal_fields"]
        approximate_fields = svtype_to_fieldsDict[svtype]["approximate_fields"]
        chromField_to_posFields = svtype_to_fieldsDict[svtype]["chromField_to_posFields"]

        # initialize a dict that maps each svID to a row with its' coordinates
        svID_to_svIDrow = {0:{}}

        # go through each sample looking for 
        for ID in ID_to_svtype_to_svDF:
            print_if_verbose(ID)

            # if the svtype exists
            if svtype in ID_to_svtype_to_svDF[ID] and len(ID_to_svtype_to_svDF[ID][svtype])>0:
                print_if_verbose(svtype)

                # get the df
                df = ID_to_svtype_to_svDF[ID][svtype].set_index("uniqueID")
                print_if_verbose("There are %i variants"%len(df))

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
        print_if_verbose("assessing monophily for %s"%svtype)

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


    print_if_verbose("There are %i vars in more than 1 species"%len(svIDs_inMoreThan1species))

    return svIDs_inMoreThan1species

    
def prune_IDtoSVTYPEtoDF_keeping_HighConfidenceVars(ID_to_svtype_to_svDF, df_samples, min_af_TraInvDel=0.5, min_af_Tan=0.4, min_af_Ins=0.25):

    """This function takes a df that maps and ID to an svtypes to df and only keeps those that are High Confidence. These are vars that are either:

    - found in any of the IDs of the same sampleID in df_samples
    - The minimum allele frequency across all the breakends is above the min_af_TraInvDel, min_af_TanInscut or min_af_Inscopy (depending on the svtype) and all breakpends are PASS

    and also that they are not in >=75% of the IDs
    """

    # first add the svID to each var. This is an ID that represents this variant across all sampleIDs, so that if there is a var overlapping across 2 samples it will have the same svID


    add_svID_to_IDtoSVTYPEtoDF(ID_to_svtype_to_svDF)

    # map each ID to the IDs of the same sample (not the sameID)
    ID_to_replicateIDs = {ID : set(df_samples[df_samples.sampleID==df_samples.loc[ID, "sampleID"]].index).difference({ID}) for ID in ID_to_svtype_to_svDF.keys()}

    # map each ID to the svTypes
    ID_to_svIDs = {ID : set.union(*[set(svDF.svID) for svDF in svtype_to_svDF.values() if len(svDF)>0]) for ID, svtype_to_svDF in ID_to_svtype_to_svDF.items()}

    # map each svID to the fraction of IDs that have it
    svID_to_fractionIDsPresent = {svID : sum([svID in svIDs for svIDs in ID_to_svIDs.values()])/len(ID_to_svIDs) for svID in set.union(*ID_to_svIDs.values())}

    # prune the df to keep only high-confidence vars
    for ID, svtype_to_svDF in ID_to_svtype_to_svDF.items():
        for svtype, svDF in svtype_to_svDF.items():

            if len(svDF)==0:  ID_to_svtype_to_svDF[ID][svtype] = svDF
            else:   

                # add wether it is in >75% of all IDs
                svDF["sv_in_75pct_samples"] = svDF.svID.apply(lambda svID: svID_to_fractionIDsPresent[svID]>0.75)

                # add whether the var is in any of the replicate IDs 
                svDF["sv_in_any_replicateIDs"] = svDF.svID.apply(lambda svID: any([svID in ID_to_svIDs[repID] for repID in ID_to_replicateIDs[ID]]))

                # add whether the breakends are all PASS
                svDF["all_bends_PASS"] = svDF.apply(lambda r: set.union(*[set([bend_info["FILTER"] for bend_info in list_breakend_info]) for list_breakend_info in r["bends_metadata_dict"].values()])=={"PASS"}, axis=1)

                # add whether all breakends have a high quality (QUAL>1000)
                svDF["all_bends_highQUAL"] = svDF.apply(lambda r: all(set.union(*[set([bend_info["QUAL"]>1000 for bend_info in list_breakend_info]) for list_breakend_info in r["bends_metadata_dict"].values()])), axis=1)


                # add wthether the minimum AF is above the min_af for this breakend
                if svtype in {"translocations", "deletions", "inversions", "remaining"}: min_af = min_af_TraInvDel
                elif svtype=="tandemDuplications": min_af = min_af_Tan
                elif svtype=="insertions": min_af = min_af_Ins
                else: raise ValueError("%s is not a valid svtype"%svtype)

                svDF["highAF"] = (svDF.estimate_AF_min>=min_af)

                # define those that have high confidence as those that have the correct af and are in all replicateIDs
                svDF["high_confidence"] = ((svDF.sv_in_any_replicateIDs) | ( (svDF.highAF) & ((svDF.all_bends_highQUAL) | (svDF.all_bends_PASS)) )) & ~(svDF.sv_in_75pct_samples)

                # keep the df that has high confidence
                ID_to_svtype_to_svDF[ID][svtype] = svDF[svDF.high_confidence]

def get_is_overlapping_query_vs_target_region(q, r):

    """This function takes two 'bed'-like regions and returns whether they are overlapping by some extent """

    return (q["chromosome"]==r["chromosome"]) and ((r["start"]<=q["start"]<=r["end"]) or (r["start"]<=q["end"]<=r["end"]) or (q["start"]<=r["start"]<=q["end"]) or (q["start"]<=r["end"]<=q["end"]))


def get_svtype_to_svfile_and_df_gridss_from_perSVade_outdir(perSVade_outdir, reference_genome):

    """This function takes from the perSVade outdir the svdict and the df_gridss"""

    # define the paths
    outdir = "%s/SVdetection_output/final_gridss_running"%perSVade_outdir
    gridss_vcf = "%s/gridss_output.raw.vcf"%outdir

    # this means that it is a cleaned dir
    if not file_is_empty(gridss_vcf):

         # get the svtype_to_svfile
        svtype_to_svfile = {svtype : "%s/%s.tab"%(outdir, svtype)  for svtype in {"insertions", "deletions", "tandemDuplications", "translocations", "inversions"}}
        svtype_to_svfile["remaining"] = "%s/unclassified_SVs.tab"%outdir 

    else: 

        # assume that it is a not cleaned dir
        gridss_vcf = "%s/gridss_output.vcf.withSimpleEventType.vcf"%outdir
        svtype_to_svfile = {file.split(".structural_variants.")[1].split(".")[0] : "%s/%s"%(outdir, file) for file in os.listdir(outdir) if ".structural_variants." in file}

    # get the df gridss
    df_gridss = add_info_to_gridssDF(load_single_sample_VCF(gridss_vcf), reference_genome)

    # keep only the ones that exist
    svtype_to_svfile = {svtype : file for svtype, file in svtype_to_svfile.items() if not file_is_empty(file)}
    print_if_verbose("There are %i svfiles"%len(svtype_to_svfile))

    return svtype_to_svfile, df_gridss

def copy_file(origin, target):

    """Copy a file with tmp"""

    target_tmp = "%s.tmp"%target

    if file_is_empty(target):

        shutil.copy2(origin, target_tmp)
        os.rename(target_tmp, target)

def clean_perSVade_outdir(outdir):

    """This function takes an outdir of perSVade and cleans it only keeping the most important files """

    # intialize the filenames
    files_to_remove = []
    file_to_dest_file = {}

    # add the most immediate files
    files_to_remove += [

       # immediate files
       "aligned_reads.bam.sorted.CollectInsertSizeMetrics.out",
       "aligned_reads.bam.sorted.coverage_per_window.tab",
       "aligned_reads.bam.sorted.histogram_insertsizes.pdf",
       "aligned_reads.bam.sorted.tmp.MarkDups.bam.bai",
       "aligned_reads.bam.sorted.tmp.MarkDups.metrics",
       "aligned_reads.bam.sorted.tmp.sortingBam_std.txt",
       "aligned_reads.bam.sorted.noMarkDups.MarkDups.metrics",
       "aligned_reads.bam.sorted.read_length_dist_first5000reads.txt",

       # files under SVdetection
       "SVdetection_output/gridss_finished.txt",
       "SVdetection_output/gridssClove_finished.txt"

    ]

    # add all the temporary files
    files_to_remove += [f for f in os.listdir(outdir) if "temporary_file" in f or f.endswith(".tmp")] 

    ########## FILES IN final_gridss_running  ######### 

    # add the files in the final_gridss_running
    final_gridss_running = "SVdetection_output/final_gridss_running"

    # add files to remove
    files_to_remove_final_gridss_running = [
        "aligned_reads.sorted.bam",
        "aligned_reads.sorted.bam.bai",
        "coverage_windows_%ibp.tab"%(window_l),
        "empty_regions.bed",                       
        "gridss_output.vcf",
        "gridss_output.vcf.idx",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf.TANDELINS.bed.3.bed", 
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf.TANDELINS.bed.5.bed",            
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf.TANDELINS.bed.target.bed",
        "simple_event_annotation.std",
        "svVCF_analysis_log.out"
    ]

    # add the names to change
    file_to_dest_file_final_gridss_running = {
        "gridss_output.vcf.withSimpleEventType.vcf":"gridss_output.raw.vcf",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf":"gridss_output.filt.vcf",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe":"gridss_output.filt.bedpe",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf":"clove_output.vcf",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf.structural_variants.deletions.bed":"deletions.tab",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf.structural_variants.inversions.bed":"inversions.tab",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf.structural_variants.remaining.tab":"unclassified_SVs.tab",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf.structural_variants.tandemDuplications.bed":"tandemDuplications.tab",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf.structural_variants.translocations.bedpe.withBlancedINFO":"translocations.tab",
        "gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf.structural_variants.insertions.bedpe.withCopiedINFO":"insertions.tab"
    }

    # keep
    files_to_remove += ["%s/%s"%(final_gridss_running, f) for f in files_to_remove_final_gridss_running]
    file_to_dest_file = {**file_to_dest_file, **{"%s/%s"%(final_gridss_running, origin) : "%s/%s"%(final_gridss_running, dest) for origin, dest in file_to_dest_file_final_gridss_running.items()}}

    ##################################################

    #### files in reference genome dir #### 
    files_to_remove_reference_genome_dir = ["reference_genome.fasta.repeat_modeler_outdir",
                                            "reference_genome_repeat_masker_outdir",
                                            "reference_genome.fasta.amb",
                                            "reference_genome.fasta.ann",
                                            "reference_genome.fasta.bwt",
                                            "reference_genome.fasta.chr_to_len.py",
                                            "reference_genome.fasta.fai",
                                            "reference_genome.fasta.pac",
                                            "reference_genome.fasta.img",
                                            "reference_genome.fasta.sa",
                                            "reference_genome_features.gff",
                                            "reference_genome.dict",
                                            "reference_genome.fasta",
                                            "reference_genome.fasta.gridsscache",
                                            "reference_genome.fasta_genomeGraph_withoutBPs.py.df_positions.py",
                                            "reference_genome.fasta_genomeGraph_withoutBPs.py.graph.py",
                                            "reference_genome.fasta.GCcontent.tab"
                                            ]

    files_to_remove += ["reference_genome_dir/%s"%f for f in files_to_remove_reference_genome_dir]
    #######################################

    ################## files in parameter_optimisation ##################

    # define the dirs
    parameter_optimisation = "SVdetection_output/parameter_optimisation"

    files_to_remove_parameter_optimisation  = ["genomeGraph_withoutBPs.df_positions.py",
                                               "genomeGraph_withoutBPs.graph.py",
                                               "genomeID_to_knownSVdict.py",
                                               "coverage_per_regions%ibb"%window_l,
                                               "simulation_reference_genome_%ibp_windows"%window_l,
                                               "benchmarking_all_filters_for_all_genomes_and_ploidies/plots"

                                            ]

    file_to_dest_file_parameter_optimisation  = {
    "coverage_per_regions%ibb/coverage_modelling_mtDNA.pdf"%window_l: "plots/coverage_modelling_mtDNA.pdf",
    "coverage_per_regions%ibb/coverage_modelling_gDNA.pdf"%window_l: "plots/coverage_modelling_gDNA.pdf",
    "benchmarking_all_filters_for_all_genomes_and_ploidies/plots/cross_accuracy_heatmaps": "plots/cross_accuracy_heatmaps",
    }  
    
    # make the simulations' SVfiles directiory
    parameter_optimisation_dir = "%s/%s"%(outdir, parameter_optimisation)
    SVfiles_dir = "%s/SVfiles"%parameter_optimisation_dir
    if os.path.isdir(parameter_optimisation_dir): 
        make_folder(SVfiles_dir)

        # go through each simulation
        for simDir in [f for f in os.listdir(parameter_optimisation_dir) if f.startswith("simulation_")]:

            # remove the dir
            files_to_remove_parameter_optimisation.append(simDir)

            # rename the SVfiles
            file_to_dest_file_parameter_optimisation = {**file_to_dest_file_parameter_optimisation,
            **{
              "%s/final_simulated_SVs/deletions.tab"%simDir: "SVfiles/%s_deletions.tab"%simDir,
              "%s/final_simulated_SVs/insertions.tab"%simDir: "SVfiles/%s_insertions.tab"%simDir,
              "%s/final_simulated_SVs/inversions.tab"%simDir: "SVfiles/%s_inversions.tab"%simDir,
              "%s/final_simulated_SVs/tandemDuplications.tab"%simDir: "SVfiles/%s_tandemDuplications.tab"%simDir,
              "%s/final_simulated_SVs/translocations.tab"%simDir: "SVfiles/%s_translocations.tab"%simDir
              }
            }

            # go through each ploidy
            for ploidyDir in [f for f in os.listdir("%s/%s"%(parameter_optimisation_dir, simDir)) if f.startswith("benchmark_GridssClove_")]:

                # define the ploidy
                ploidy = "_".join(ploidyDir.split("_")[2:])

                file_to_dest_file_parameter_optimisation = {**file_to_dest_file_parameter_optimisation,
                **{
                  "%s/%s/plots_benchmark"%(simDir, ploidyDir): "plots/plots_benchmark_%s"%ploidy,
                  }             
                }

    # keep
    files_to_remove += ["%s/%s"%(parameter_optimisation, f) for f in files_to_remove_parameter_optimisation]
    file_to_dest_file = {**file_to_dest_file, **{"%s/%s"%(parameter_optimisation, origin) : "%s/%s"%(parameter_optimisation, dest) for origin, dest in file_to_dest_file_parameter_optimisation.items()}}

    #####################################################################

    ####### REMOVE AND CHANGE FILENAMES #######

    # change name
    for o, d in file_to_dest_file.items():
        origin = "%s/%s"%(outdir, o)
        dest = "%s/%s"%(outdir, d)

        # try for files
        if not file_is_empty(origin): os.rename(origin, dest)   

        # try for directories
        if os.path.isdir(origin) and not os.path.isdir(dest): os.rename(origin, dest)   

    # remove 
    for f in files_to_remove:
        file = "%s/%s"%(outdir, f)
        remove_file(file)
        delete_folder(file)

    ###########################################

def get_ID_to_svtype_to_svDF_for_setOfGenomes_highConfidence(close_shortReads_table, reference_genome, outdir, replace=False, threads=4, mitochondrial_chromosome="mito_C_glabrata_CBS138", job_array_mode="local", max_ncores_queue=768, time_perSVade_running="48:00:00", queue_jobs="bsc_ls", max_nvars=100, name_job_array="getRealSVs"):

    """Generates a dict that maps each sample in genomes_withSV_and_shortReads_table to an svtype and a DF with all the info about several vars. It only gets the high-confidence vars.

    realSV_calling_on can be reads or assembly"""

    # load the df
    df_genomes = pd.read_csv(close_shortReads_table, sep="\t").set_index("runID")

    # define an outdir that will store all the real_vars
    make_folder(outdir)
    all_realVars_dir = "%s/all_realVars"%(outdir)
    if replace is True: delete_folder(all_realVars_dir)
    make_folder(all_realVars_dir)

    # define the final object that contains them all
    ID_to_svtype_to_svDF_file = "%s/ID_to_svtype_to_svDF.py"%all_realVars_dir

    if file_is_empty(ID_to_svtype_to_svDF_file) or replace is True:

        # initialize a list of cmds to run
        all_cmds = []

        # initialize dicts that keep them all together
        all_sampleID_to_svtype_to_file = {}
        all_sampleID_to_dfGRIDSS = {}

        # define the name of the final file name
        final_file_name = "perSVade_finished_file.txt"

        # check if there are some jobs to run SV calling on
        njobs_to_run_SVcalling_on = sum([file_is_empty("%s/shortReads_realVarsDiscovery_%s/%s"%(all_realVars_dir,ID, final_file_name)) for ID, row in df_genomes.iterrows()])

        print_if_verbose("There are %i jobs still to run"%njobs_to_run_SVcalling_on)

        # generate all real vars
        for ID, row in df_genomes.iterrows():
            print_if_verbose(ID)

            # define a file with the df_gridss and svtype_to_svfile
            df_gridss_svtype_to_svfile_objectFile = "%s/df_gridss_svtype_to_svfile_tuple_%s.py"%(all_realVars_dir,ID)

            if file_is_empty(df_gridss_svtype_to_svfile_objectFile) or replace is True:

                # run in the gridss and clove with the fast parameters
                outdir_gridssClove = "%s/shortReads_realVarsDiscovery_%s"%(all_realVars_dir,ID); make_folder(outdir_gridssClove)

                # define the previous important files
                final_file = "%s/%s"%(outdir_gridssClove, final_file_name)

                # define the previous repeats file 
                previous_repeats_table = "%s.repeats.tab"%reference_genome
                if file_is_empty(previous_repeats_table): raise ValueError("%s should exist"%previous_repeats_table)

                # only contine if the final file is not defined
                if file_is_empty(final_file) or replace is True:
                    print_if_verbose("getting vars for %s"%ID)

                    # define the cmd. This is a normal perSvade.py run with the vars of the previous dir  
                    cmd = "python %s -r %s --threads %i --outdir %s  --mitochondrial_chromosome %s --fast_SVcalling --previous_repeats_table %s"%(perSVade_py, reference_genome, threads, outdir_gridssClove, mitochondrial_chromosome, previous_repeats_table)

                    # add arguments depending on the pipeline
                    if replace is True: cmd += " --replace"

                    # add the input
                    all_keys_df = set(df_genomes.keys())

                    # reads
                    if "short_reads1" in all_keys_df and "short_reads2" in all_keys_df: cmd += " -f1 %s -f2 %s"%(row["short_reads1"], row["short_reads2"])

                    # bams
                    elif "sorted_bam" in all_keys_df: cmd += " -sbam %s"%(row["sorted_bam"])

                    else: raise ValueError("The provided close_shortReads_table is not valid") 

                    # if the running in slurm is false, just run the cmd
                    if job_array_mode=="local": run_cmd(cmd)
                    elif job_array_mode=="greasy": 
                        all_cmds.append(cmd)
                        continue

                    else: raise ValueError("%s is not valid"%job_array_mode)


                # if the mode is greasy and there are jobs to run, continue
                if njobs_to_run_SVcalling_on>0 and job_array_mode=="greasy": continue

                # define the svdict
                print_if_verbose("running get_svtype_to_svfile_and_df_gridss_from_perSVade_outdir")
                svtype_to_svfile,  df_gridss = get_svtype_to_svfile_and_df_gridss_from_perSVade_outdir(outdir_gridssClove, reference_genome)

                # save
                print_if_verbose("saving %s"%df_gridss_svtype_to_svfile_objectFile)
                df_gridss_svtype_to_svfile_objectFile_tmp = "%s.tmp"%df_gridss_svtype_to_svfile_objectFile
                save_object((svtype_to_svfile , df_gridss), df_gridss_svtype_to_svfile_objectFile_tmp)
                os.rename(df_gridss_svtype_to_svfile_objectFile_tmp, df_gridss_svtype_to_svfile_objectFile)

            else: svtype_to_svfile , df_gridss = load_object(df_gridss_svtype_to_svfile_objectFile)

            # keep 
            all_sampleID_to_svtype_to_file[ID] =  svtype_to_svfile
            all_sampleID_to_dfGRIDSS[ID] = df_gridss

        # if yoy are running on slurm, get it in a job array
        if job_array_mode=="greasy": 

            if len(all_cmds)>0: 
                print_if_verbose("submitting %i jobs to the cluster for the real data. The files can be monitored from %s"%(len(all_cmds), all_realVars_dir))
                jobs_filename = "%s/jobs.getting_realSVs"%all_realVars_dir
                open(jobs_filename, "w").write("\n".join(all_cmds))

                generate_jobarray_file_greasy(jobs_filename, walltime=time_perSVade_running,  name=name_job_array, queue=queue_jobs, sbatch=True, ncores_per_task=threads, constraint="", number_tasks_to_run_at_once="all", max_ncores_queue=max_ncores_queue )

                print_if_verbose("Exiting... You have to wait until all the jobs in testRealSVs are done. Wait until the jobs are done and rerun this pipeline to continue")
                sys.exit(0)

        # get a df that has all the info for each SV, and then the df with allele freq, metadata and 
        print_if_verbose("running get_sampleID_to_svtype_to_svDF_filtered")
        ID_to_svtype_to_svDF = get_sampleID_to_svtype_to_svDF_filtered(all_sampleID_to_svtype_to_file, all_sampleID_to_dfGRIDSS)

        # keep 5x the maximum number of vars
        ID_to_svtype_to_svDF = {ID : {svtype : svDF.iloc[0:max_nvars*5] for svtype,svDF in svtype_to_svDF.items()} for ID,svtype_to_svDF in ID_to_svtype_to_svDF.items()}

        # preseve only the high confidence vars from ID_to_svtype_to_svDF
        print_if_verbose("running prune_IDtoSVTYPEtoDF_keeping_HighConfidenceVars")
        prune_IDtoSVTYPEtoDF_keeping_HighConfidenceVars(ID_to_svtype_to_svDF, df_genomes)

        # save
        print_if_verbose("saving")
        save_object(ID_to_svtype_to_svDF, ID_to_svtype_to_svDF_file)

    else: ID_to_svtype_to_svDF = load_object(ID_to_svtype_to_svDF_file)

    # remove the bam files and reference files
    print_if_verbose("cleaning files")
    for ID, row in df_genomes.iterrows():

        # run in the gridss and clove with the fast parameters
        outdir_gridssClove = "%s/shortReads_realVarsDiscovery_%s"%(all_realVars_dir,ID)

        for f in ["aligned_reads.bam.sorted", "aligned_reads.bam.sorted.bai", "reference_genome_dir"]:
            file = "%s/%s"%(outdir_gridssClove, f)
            delete_folder(file)
            remove_file(file)

        # clean again
        clean_perSVade_outdir(outdir_gridssClove)

        # remove the file that has the sample
        remove_file("%s/df_gridss_svtype_to_svfile_tuple_%s.py"%(all_realVars_dir,ID))

    return ID_to_svtype_to_svDF

def add1_unless_it_is_minus1(x):

    """Takes an int and adds 1 unless it is -1"""

    if x==-1: return x
    else: return x+1 

def set_position_to_max(pos, maxPos):

    """Sets a position to a maximum"""

    if pos>maxPos: return maxPos
    else: return pos

def get_compatible_real_svtype_to_file(close_shortReads_table, reference_genome, outdir, replace=False, threads=4, max_nvars=100, mitochondrial_chromosome="mito_C_glabrata_CBS138", job_array_mode="local", max_ncores_queue=768, time_perSVade_running="48:00:00", queue_jobs="bsc_ls", name_job_array="getRealSVs"):

    """This function generates a dict of svtype to the file for SVs that are compatible and ready to insert into the reference_genome. All the files are written into outdir. Only a set of 'high-confidence' SVs are reported, which are those that, for each sampleID inclose_shortReads_table, have a reasonable minimum allele frequency and all breakpoints with 'PASS' and are found in all the genomes of the same sampleID.

    At the end, this pipeline reports a set of compatible SVs, that are ready to insert into RSVsim (so that the coordinates are 1-based). 

    Rememeber that this will not include 'remaining' events, as these can't be inserted 
    
    max_nvars is the maximum number of variants of each type that will be generated
    """

    # initialize the start time
    pipeline_start_time = time.time()

    # get all the high-confidence real variants
    ID_to_svtype_to_svDF = get_ID_to_svtype_to_svDF_for_setOfGenomes_highConfidence(close_shortReads_table, reference_genome, outdir, replace=replace, threads=threads, mitochondrial_chromosome=mitochondrial_chromosome, job_array_mode=job_array_mode, max_ncores_queue=max_ncores_queue, time_perSVade_running=time_perSVade_running, queue_jobs=queue_jobs, max_nvars=max_nvars, name_job_array=name_job_array)

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

    # inialize a file that will contain the number of each real SV
    nSVs_statistics_filecontent = "svtype\tnSVs\n" 

    for svtype in all_svtypes:
        outfile_compatible_SVs = "%s/%s.tab"%(SVs_compatible_to_insert_dir, svtype)

        if file_is_empty(outfile_compatible_SVs) or replace is True:
        #if True:

            # initalize a df with all the compatible svDFs
            compatible_svDF = pd.DataFrame(columns=svtype_to_fieldsDict[svtype]["all_fields"])

            # go through each ID
            print_if_verbose("getting %s"%svtype)
            for ID in ID_to_svtype_to_svDF.keys():
                print_if_verbose(ID)

                # debug empty dfs
                if svtype not in ID_to_svtype_to_svDF[ID] or len(ID_to_svtype_to_svDF[ID][svtype])==0: continue

                svDF = ID_to_svtype_to_svDF[ID][svtype].set_index("uniqueID", drop=False).iloc[0:max_nvars]

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
                        df_bed_allRegions = df_bed_allRegions.append(sv_bed, sort=True)

                        # add to the compatible SVs
                        compatible_svDF = compatible_svDF.append(sv_series_df, sort=True)


            if len(compatible_svDF)>0:

                # get only the important fields
                compatible_svDF = compatible_svDF[svtype_to_fieldsDict[svtype]["all_fields"]]
                compatible_svDF.index = list(range(len(compatible_svDF)))

                # define the maping between the position field and the chromosome field
                posF_to_chrF = svtype_to_fieldsDict[svtype]["positionField_to_chromosome"]

                # add +1 to all the positions (this is because RSVSim requires 1-based positions), also that the last position of the chromosome is not exceeded
                for f in svtype_to_fieldsDict[svtype]["position_fields"]: compatible_svDF[f] = compatible_svDF.apply(lambda r: set_position_to_max(add1_unless_it_is_minus1(r[f]), chr_to_len[r[posF_to_chrF[f]]]), axis=1)


                # if this number exceeds the number of variants it will chop the df
                if len(compatible_svDF)>max_nvars: compatible_svDF = compatible_svDF.iloc[0:max_nvars]

            # write the compatible svDF into the final set of vars
            compatible_svDF.to_csv(outfile_compatible_SVs, sep="\t", header=True, index=False)

        # keep 
        compatible_real_svtype_to_file[svtype] = outfile_compatible_SVs

        # write the number
        nVars = len(pd.read_csv(outfile_compatible_SVs, sep="\t"))
        print_if_verbose("Defining %i compatible %s"%(nVars, svtype))
        nSVs_statistics_filecontent += "%s\t%i\n"%(svtype, nVars)

    # write the statistics file
    open("%s/number_SVs.tab"%SVs_compatible_to_insert_dir, "w").write(nSVs_statistics_filecontent)


    #####################################################

    print_if_verbose("--- the generation of real SVs took %s seconds in %i cores ---"%(time.time() - pipeline_start_time, threads))

    return compatible_real_svtype_to_file

def get_insert_size_distribution(sorted_bam, replace=False, threads=4):

    """Takes a bam file of aligned paired end reads and retuns the mean and SD insert size of the library."""

    # define outfiles
    hist_file = "%s.histogram_insertsizes.pdf"%sorted_bam
    outfile = "%s.CollectInsertSizeMetrics.out"%sorted_bam; outfile_tmp = "%s.tmp"%outfile

    # run
    if file_is_empty(outfile) or replace is True: 
        remove_file(outfile_tmp)

        # downsample bam to 1%
        print_if_verbose("getting 1pct of the reads to calculate insert size")
        sampled_bam = downsample_bamfile_keeping_pairs(sorted_bam, fraction_reads=0.01, replace=replace, threads=threads)

        # run the calculation of insert sizes
        picard_insertSize_std = "%s.generating.std"%outfile_tmp
        print_if_verbose("calculating insert size distribution. The std is in %s"%picard_insertSize_std)
        run_cmd("%s CollectInsertSizeMetrics HISTOGRAM_FILE=%s INPUT=%s OUTPUT=%s > %s 2>&1"%(picard_exec, hist_file, sampled_bam, outfile_tmp, picard_insertSize_std))
        remove_file(picard_insertSize_std)
        remove_file(sampled_bam)

        os.rename(outfile_tmp, outfile)

    # get stats
    wrong_foot_lines = [l for l in open(outfile, "r").readlines() if len(l.split("\t"))==2 and not l.startswith("## METRICS CLASS")]
    df = pd.read_csv(outfile, sep="\t", skip_blank_lines=False, header=6, skipfooter=len(wrong_foot_lines)+2, engine='python').iloc[0]

    return (int(float(df["MEDIAN_INSERT_SIZE"])), int(float(df["MEDIAN_ABSOLUTE_DEVIATION"])))


def get_windows_infoDF_with_predictedFromFeatures_coverage(genome, distToTel_chrom_GC_to_coverage_fn, expected_coverage_per_bp, replace=False, threads=4, make_plots=True):

    """This function gets a genome and returns a df for windows of the genome and the relative coverage predicted from distToTel_chrom_GC_to_coverage_fn"""

    windows_infoDF_file = "%s_windows_with_predictedRelCov_from_features.py"%genome

    if file_is_empty(windows_infoDF_file) or replace is True:
        print_if_verbose("getting relCov predicted from feats")

        # index the genome of interest if not already done
        if file_is_empty("%s.fai"%genome) or replace is True: 

            faidx_std = "%s.indexing.std"%genome
            print_if_verbose("running faidx. The std is in %s"%faidx_std)
            run_cmd("%s faidx %s > %s 2>&1"%(samtools, genome, faidx_std))
            remove_file(faidx_std)

        ##### get the windows df ####

        # get the file
        windows_file = "%s.windows%ibp.bed"%(genome, window_l)
        windows_file_stderr = "%s.generating.stderr"%windows_file
        print_if_verbose("running makewindows. The stderr is in %s"%windows_file_stderr)
        run_cmd("%s makewindows -g %s.fai -w %i > %s 2>%s"%(bedtools, genome, window_l, windows_file, windows_file_stderr))
        remove_file(windows_file_stderr)

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
            print_if_verbose("making plots")

            # define the chromosomes that are interesting
            interesting_chromosomes = [c for c in all_chromosomes if len(df[df.chromosome==c])>=5]
            print_if_verbose("There are %i/%i chromosomes with >=5 windows"%(len(interesting_chromosomes), len(all_chromosomes)))

            fig = plt.figure(figsize=(10, len(interesting_chromosomes)*4.5))
            for I, chrom in enumerate(interesting_chromosomes):

                # get df
                df_c = df[df.chromosome==chrom]

                # print each of the predictions
                ax = plt.subplot(len(interesting_chromosomes), 1, I+1)
                sns.lineplot(x="start", y="predicted_relative_coverage", data=df_c, linewidth=2, color="blue")

                ax.set_title(chrom)

            try:

                fig.tight_layout()  # otherwise the right y-label is slightly 
                filename="%s_predicted_relative_coverage.pdf"%(windows_file)
                fig.savefig(filename, bbox_inches='tight');
                plt.close(fig)

            except: print_if_verbose("The rendering of the coverage as image did not work. This is likely because there are too many plots ")

        # save
        save_object(df, windows_infoDF_file)

    else: df = load_object(windows_infoDF_file)

    return df

def get_read_length(bamfile, threads=4, nreads=5000, replace=False):

    """Calculates the readlength for a bamfile"""

    readlen_dist_file = "%s.read_length_dist_first%ireads.txt"%(bamfile, nreads); readlen_dist_file_tmp = "%s.tmp"%readlen_dist_file
    if file_is_empty(readlen_dist_file) or replace is True:

        samtools_read_len_stderr = "%s.generating.stderr"%readlen_dist_file_tmp

        print_if_verbose("Running samtools view. The following command will throw a warning stating that 'samtools view: writing to standard output failed: Broken pipe'. This is because the output of samtools view is piped, which is expected. The stderr is in %s"%samtools_read_len_stderr)
        cmd = "%s view --threads %i %s | head -n %i | cut -f10 | perl -ne 'chomp;print length($_) . \"\n\"' | sort > %s 2>%s"%(samtools, threads, bamfile, nreads, readlen_dist_file_tmp, samtools_read_len_stderr)
        run_cmd(cmd)
        remove_file(samtools_read_len_stderr)

        os.rename(readlen_dist_file_tmp, readlen_dist_file)

    return int(np.median([int(l.strip()) for l in open(readlen_dist_file, "r").readlines()]))


def count_number_read_pairs(bamfile, replace=False, threads=4):

    """counts the total number of reads of a bamfile"""

    read_count_file = "%s.flagstat"%bamfile; read_count_file_tmp = "%s.tmp"%read_count_file

    # get the total n reads
    if file_is_empty(read_count_file) or replace is True:

        read_count_stderr = "%s.generating.stderr"%read_count_file
        print_if_verbose("calculating n reads. The stderr is in %s"%read_count_stderr)
        run_cmd("%s flagstat --threads %i %s > %s 2>%s"%(samtools, threads, bamfile, read_count_file_tmp, read_count_stderr))
        remove_file(read_count_stderr)

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
        print_if_verbose("WARNING: The length of the region is %i and the length of the reads is %i. Setting the length of the region to %i. We will try to set the region to the 5' of the chromosome, and to the 3' if not possible."%(len_region, read_length, min_length_region))

        # get the length of the chromosome
        len_chr = [len(chrom) for chrom in SeqIO.parse(genome, "fasta") if chrom.id==chromosome][0]

        # define the region to extend
        length_to_extend = min_length_region - len_region

        # if the start can be extended
        if (start - length_to_extend)>=0: start = start - length_to_extend

        # if the end can be extended
        elif (end + length_to_extend)<=len_chr: end = end + length_to_extend

        # if the region is longer than the chrom, just set the whole chromosome
        elif len_chr<min_length_region: 
            start = 0
            end = len_chr

        else: raise ValueError("Region %s cannot be corrected according to read length, maybe because the chromosome is to short. Check that %s is ok"%(ID, chromosome))

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
 
    if any([file_is_empty(x) for x in {fastq_1, fastq_2}]) or replace is True:


        # define the stderr of the reads
        std = "%s/%s_std.txt"%(outdir, ID)
        print_if_verbose("simulating %s. The std is in %s"%(ID, std))

        run_cmd("%s -e %.2f -N %i -1 %i -2 %i -r 0.0 -R 0.0 -X 0.0 -h -d %i -s %i %s %s %s > %s 2>&1"%(wgsim, error_rate, readPairs, read_length, read_length, median_insert_size, median_insert_size_sd, region_fasta, fastq_1_tmp, fastq_2_tmp, std)) # bear in mind '-A FLOAT discard if the fraction of ambiguous bases higher than FLOAT [0.05]'

        # check that the generated files are not empty, this may make the parallelization to fail
        for f in [fastq_1_tmp, fastq_2_tmp]:
            if os.stat(f).st_size==0: print_if_verbose("!!WARNING: No reads could be generated for region %s into %s. It is likely too short or it has too many Ns. This may cause problems with multiprocessing."%(ID, f))

        # check that everything is fine in the std and remove it
        if all([l.startswith("[wgsim") or l.startswith("Failed to produce a sequence with insufficient Ns") for l in open(std, "r").readlines()]): remove_file(std)
        else: raise ValueError("There may have been an error in generating some reads. Check %s"%std)

        os.rename(fastq_1_tmp, fastq_1); os.rename(fastq_2_tmp, fastq_2); 

    # remove the fasta file
    remove_file(region_fasta)

    return (fastq_1, fastq_2)

def run_wgsim_pairedEnd_per_windows_in_parallel(df_windows, genome, outdir, read_length,  median_insert_size, median_insert_size_sd, replace=False, max_n_windows_at_once=2000, error_rate=0.02):

    """Takes a dataframe with windows of ["chromosome", "start", "end", "readPairs"] and writes, under outdir, two fastq.gz files of simulated paired end  reads. The parallel runs are written in subfolder under outir which is finally removed. It returns a tuple of the two fastq.gz files.

    max_n_windows_at_once indicates the number of windows that can be processed at once. """

    # make the outdir if it does not exist
    make_folder(outdir)
    print_if_verbose("running simulation in a parallelized manner. Each window will be run in a different core. Running on chunks of %i windows"%max_n_windows_at_once)

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
            print_if_verbose("Working on chunk %i of %i"%(I+1, len(chunks_indices)))

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

                print_if_verbose("opening multiprocessing on %i cores on %i windows"%(multiproc.cpu_count(), len(df_chunk)))

                # initialize the pool
                start_time =  time.time()
                with  multiproc.Pool(multiproc.cpu_count()) as pool:

                    # define the list of args
                    args = [(genome, chromosome, start, end, readPairs, read_length, parallel_files_outdir, median_insert_size, median_insert_size_sd, replace, error_rate) for chromosome, start, end, readPairs in df_chunk.values]

                    read1_read2_tuples = pool.starmap(run_wgsim_pairedEnd_for_window, args)
                    pool.close()
                    pool.terminate()
                    pool.join()

                # run gzip in parallel
                pigz_std = "%s.pigz_compression.std"%chunk_all_fastqgz_1
                print_if_verbose("running pgzip in parallel. The std is in %s"%pigz_std)
                run_cmd("%s --fast %s/* > %s 2>&1"%(pigz, parallel_files_outdir, pigz_std))

                # concatenate them all
                chunk_all_fastqgz_1_tmp = "%s.tmp"%chunk_all_fastqgz_1
                chunk_all_fastqgz_2_tmp = "%s.tmp"%chunk_all_fastqgz_2
                chunk_all_fastqgz_1_stderr = "%s.generating.stderr"%chunk_all_fastqgz_1_tmp
                chunk_all_fastqgz_2_stderr = "%s.generating.stderr"%chunk_all_fastqgz_2_tmp

                print_if_verbose("Integrating all in one for chunk %i. The stderrs are in %s and %s"%(I+1, chunk_all_fastqgz_1_stderr, chunk_all_fastqgz_2_stderr))


                run_cmd("cat %s/*_read1.fq.gz > %s 2>%s"%(parallel_files_outdir, chunk_all_fastqgz_1_tmp, chunk_all_fastqgz_1_stderr)) 
                run_cmd("cat %s/*_read2.fq.gz > %s 2>%s"%(parallel_files_outdir, chunk_all_fastqgz_2_tmp, chunk_all_fastqgz_2_stderr))

                remove_file(chunk_all_fastqgz_1_stderr)
                remove_file(chunk_all_fastqgz_2_stderr)

                # rename to keep
                os.rename(chunk_all_fastqgz_1_tmp, chunk_all_fastqgz_1); os.rename(chunk_all_fastqgz_2_tmp, chunk_all_fastqgz_2)

                # print the time it took to generate these files
                print_if_verbose("--- creating reads for this chunk took %s seconds ---"%(time.time() - start_time))

            # remoove intermediate files of the fastq generation
            print_if_verbose("removing files")
            delete_folder(parallel_files_outdir)

            # keep the fastq files
            chunks_fastq_files += [chunk_all_fastqgz_1, chunk_all_fastqgz_2]

        # integrate all the chunks in one
        allChunks_fastqgz_1_tmp = "%s.tmp"%allChunks_fastqgz_1
        allChunks_fastqgz_2_tmp = "%s.tmp"%allChunks_fastqgz_2
        chunks_read1_files = "%s/chunk*_ofmax%i_all_reads1.fq.gz"%(outdir, max_n_windows_at_once)
        chunks_read2_files = "%s/chunk*_ofmax%i_all_reads2.fq.gz"%(outdir, max_n_windows_at_once)

        allChunks_fastqgz_1_tmp_stderr = "%s.generating.stderr"%allChunks_fastqgz_1_tmp
        allChunks_fastqgz_2_tmp_stderr = "%s.generating.stderr"%allChunks_fastqgz_2_tmp
        print_if_verbose("integrating all the gzipped reads into one.The stderrs are in %s and %s"%(allChunks_fastqgz_1_tmp_stderr, allChunks_fastqgz_2_tmp_stderr))

        run_cmd("cat %s > %s 2>%s"%(chunks_read1_files, allChunks_fastqgz_1_tmp, allChunks_fastqgz_1_tmp_stderr)) 
        run_cmd("cat %s > %s 2>%s"%(chunks_read2_files, allChunks_fastqgz_2_tmp, allChunks_fastqgz_2_tmp_stderr))

        remove_file(allChunks_fastqgz_1_tmp_stderr)
        remove_file(allChunks_fastqgz_2_tmp_stderr)

        # remove chunks' fastq's
        for f in chunks_fastq_files: remove_file(f)

        # at the end rename so that everything is clean
        os.rename(allChunks_fastqgz_1_tmp, allChunks_fastqgz_1)
        os.rename(allChunks_fastqgz_2_tmp, allChunks_fastqgz_2)

    print_if_verbose("All the reads simulated")

    return (allChunks_fastqgz_1, allChunks_fastqgz_2) 

def run_seqtk_rename(origin_fastqgz, dest_fastqgz):

    """This function takes an origin and a dest fastq gz and runs seqtk to rename the reads"""

    if file_is_empty(seqtk): raise ValueError("seqtk is expected to be in %s"%seqtk)

    seqtk_stderr = "%s.generating.stderr"%dest_fastqgz
    print_if_verbose("Running seqtk. The stderr is in %s"%seqtk_stderr)
    run_cmd("%s rename %s read_ | %s -c > %s 2>%s"%(seqtk, origin_fastqgz, pigz, dest_fastqgz, seqtk_stderr))
    remove_file(seqtk_stderr)


def simulate_readPairs_per_window(df_windows, genome, npairs, outdir, read_length,  median_insert_size, median_insert_size_sd, replace=False, threads=4):

    """Simulates npairs of short reads in a way that is consistent with the predicted_relative_coverage of the df_windows. At the end, this will correct the reads by """

    # make the outdir if not there
    make_folder(outdir)

    # define the outputdirs
    all_fastqgz_1 = "%s/all_reads1.fq.gz"%outdir
    all_fastqgz_2 = "%s/all_reads2.fq.gz"%outdir
    all_fastqgz_1_correct = "%s/all_reads1.correct.fq.gz"%outdir
    all_fastqgz_2_correct = "%s/all_reads2.correct.fq.gz"%outdir
    outdir_whole_chromosomes = "%s/whole_chromosomes_simulation_with_min_reads"%outdir
    outdir_per_windows = "%s/per_window_simulation_with_extra_reads"%outdir

    if any([file_is_empty(f) for f in [all_fastqgz_1_correct, all_fastqgz_2_correct]]) or replace is True:

        # define the tmp files
        all_fastqgz_1_correct_tmp = "%s/all_reads1.correct.tmp.fq.gz"%outdir
        all_fastqgz_2_correct_tmp = "%s/all_reads2.correct.tmp.fq.gz"%outdir

        # get the raw reads from wgsim
        if any([file_is_empty(f) for f in [all_fastqgz_1, all_fastqgz_2]]) or replace is True:

            # delete previous folders and files
            #delete_folder(outdir_whole_chromosomes); delete_folder(outdir_per_windows); remove_file(all_fastqgz_1); remove_file(all_fastqgz_2)

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
                print_if_verbose("geting sliding window for %s"%chromosome)

                # get the df with a unique numeric index
                df_c = df_windows[df_windows.chromosome==chromosome][["chromosome", "start", "end", "total_readPairs", "window_length"]].sort_values(by=["chromosome", "start"])
                df_c.index = list(range(len(df_c)))

                if len(df_c)>1:

                    # this only applies when there is more than 1 window

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

                else: 

                    # else all the reads go to this window
                    df_c_both = df_c.sort_values(by=["chromosome", "start"])
                    df_c_both["total_readPairs"] = df_c_both.total_readPairs.apply(int) + 1

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
            df_windows_chromosomes_minCov = pd.DataFrame({I : {"start": 0, "end": len_chr, "chromosome": chrom, "readPairs": sum(df_windows_sliding[df_windows_sliding.chromosome==chrom]["min_readPairs"])} for I, (chrom, len_chr) in enumerate(chrom_to_len.items())}).transpose()
            wholeChr_fastqgz1, wholeChr_fastqgz2 = run_wgsim_pairedEnd_per_windows_in_parallel(df_windows_chromosomes_minCov, genome, outdir_whole_chromosomes, read_length, median_insert_size, median_insert_size_sd, replace=replace)

            # now simulate the reads for the regions, only considering the extra_readPairs
            df_windows_sliding["readPairs"] = df_windows_sliding.extra_readPairs        
            perWindow_fastqgz1, perWindow_fastqgz2 = run_wgsim_pairedEnd_per_windows_in_parallel(df_windows_sliding, genome, outdir_per_windows, read_length, median_insert_size, median_insert_size_sd, replace=replace)

            # now concatenate all the simulated reads
            all_fastqgz_1_tmp = "%s.tmp"%all_fastqgz_1
            all_fastqgz_2_tmp = "%s.tmp"%all_fastqgz_2

            all_fastqgz_1_tmp_stderr = "%s.generating.stderr"%all_fastqgz_1_tmp
            all_fastqgz_2_tmp_stderr = "%s.generating.stderr"%all_fastqgz_2_tmp

            print_if_verbose("Integrating all in one. The stderrs are in %s and %s"%(all_fastqgz_1_tmp_stderr, all_fastqgz_2_tmp_stderr))

            run_cmd("cat %s %s > %s 2>%s"%(wholeChr_fastqgz1, perWindow_fastqgz1, all_fastqgz_1_tmp, all_fastqgz_1_tmp_stderr))
            run_cmd("cat %s %s > %s 2>%s"%(wholeChr_fastqgz2, perWindow_fastqgz2, all_fastqgz_2_tmp, all_fastqgz_2_tmp_stderr))

            remove_file(all_fastqgz_1_tmp_stderr)
            remove_file(all_fastqgz_2_tmp_stderr)

            # rename to keep
            os.rename(all_fastqgz_1_tmp, all_fastqgz_1); os.rename(all_fastqgz_2_tmp, all_fastqgz_2)


        # correct the names with seqtk rename
        print_if_verbose("running seqtk to correct read names from wgsim")
        run_seqtk_rename(all_fastqgz_1, all_fastqgz_1_correct_tmp)
        run_seqtk_rename(all_fastqgz_2, all_fastqgz_2_correct_tmp)

        for f in [all_fastqgz_1_correct_tmp, all_fastqgz_2_correct_tmp]:
            if file_is_empty(f): raise ValueError("%s is empty. Something went wrong with wgsim"%f)

        # keep
        os.rename(all_fastqgz_1_correct_tmp, all_fastqgz_1_correct)
        os.rename(all_fastqgz_2_correct_tmp, all_fastqgz_2_correct)

    # remove previously generated folders
    delete_folder(outdir_whole_chromosomes); delete_folder(outdir_per_windows)

    # delete all reads
    remove_file(all_fastqgz_1)
    remove_file(all_fastqgz_2)

    # return the simulated reads
    return (all_fastqgz_1_correct, all_fastqgz_2_correct) 

def simulate_and_align_PairedReads_perWindow(df_windows, genome_interest, reference_genome, npairs, read_length, outdir, median_insert_size, median_insert_size_sd, replace=False, threads=4):

    """Takes a dataframe with windows of the genome, which also has a predicted_relative_coverage (which indicates by how much should the coverage be multiplied in this window). This function generates a fastq (and deletes it at the end), aligns it and returns the bam. All files are written under outdir. It returns the aligned bamfile. All the chromosomes are simulated as linear."""

    print_if_verbose("Simulating reads and aliging them for %s"%genome_interest)
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
        print_if_verbose("deleting reads")
        delete_folder(outdir_reads)

    # record the time consumed
    print_if_verbose("--- generating %i reads from %s and aligning them  took %s seconds ---"%(npairs, genome_interest, time.time() - start_time))

    return sim_sorted_bam


def get_fractions_reads_for_ploidy(ploidy):

    """Takes a ploidy and returns the fraction_var and fraction_ref."""

    # define the fraction of reads that come from each genome
    if ploidy in {"haploid", "diploid_homo"}:

        fraction_var = 1.0
        fraction_ref = 0.0

    elif ploidy=="diploid_hetero":

        fraction_var = 0.5
        fraction_ref = 0.5 

    else:

        # These are specific combinations of genomes (polyploids or populations)
        try:
            typeGenome_to_ncopies = {x.split(":")[0] : int(x.split(":")[1]) for x in ploidy.split("_")}
            typeGenome_to_fraction_reads = {type_genome: ncopies / sum(typeGenome_to_ncopies.values()) for type_genome, ncopies in typeGenome_to_ncopies.items()}

        except: raise ValueError("ploidy %s is incorrect. The format is ref:2_var:1 , which would be one copy of variant genome for every two copies of reference genome"%ploidy)
        
        fraction_var = typeGenome_to_fraction_reads["var"]
        fraction_ref = typeGenome_to_fraction_reads["ref"]

    return fraction_var, fraction_ref



def merge_2bams(bamA, bamB, merged, threads=4):

    """Merges two bams"""
    remove_file(merged)

    merged_std = "%s.generating.std"%merged
    print_if_verbose("merging %s and %s into %s. The std is in %s"%(bamA, bamB, merged, merged_std))
    run_cmd("%s merge --threads %i %s %s %s > %s 2>&1"%(samtools, threads, merged, bamA, bamB, merged_std))
    remove_file(merged_std)

def sort_bam(bam, sorted_bam, threads=4):

    """Sorts a bam file into sorted_bam"""

    sorting_std = "%s.generating.std"%sorted_bam
    print_if_verbose("sorting bam. The std is in %s"%sorting_std)
    run_cmd("%s sort --threads %i -o %s %s > %s 2>&1"%(samtools, threads, sorted_bam, bam, sorting_std))
    remove_file(sorting_std)

def index_bam(bam, threads=4):

    """indexes bam and creates a .bai file"""

    indexing_std = "%s.generating.std"%bam
    print_if_verbose("indexing bam. The std is in %s"%indexing_std)
    run_cmd("%s index -@ %i %s > %s 2>&1"%(samtools, threads, bam, indexing_std))
    remove_file(indexing_std)


def merge_2bams_sort_and_index(bamfileA, bamfileB, merged_bam, merged_sorted_bam, merged_sorted_bam_index, threads=4):

    """Takes two bams and merges them. writes files under merged_bam_prefix"""

    # first merge
    merge_2bams(bamfileA, bamfileB, merged_bam, threads=threads)

    # sort 
    sort_bam(merged_bam, merged_sorted_bam, threads=threads)

    # index
    index_bam(merged_sorted_bam, threads=threads)

    # remove the unsorted bam
    remove_file(merged_bam)


def get_merged_bamfile_for_ploidy(variant_bamfile, reference_bamfile, ploidy, replace=False, threads=4):

    """This function takes two bam files, one that includes aligned reads for a reference genome and another that includes aligned reads for a variat genome. It subsamples each of the bamfiles so that there is a proportion of reads comming from each reference and variant, and the proportion is indicated by ploidy, which can be any of the strings in the default values of target_ploidies of run_GridssClove_optimising_parameters."""

    print_if_verbose("merging %s and %s in ploidy %s"%(variant_bamfile, reference_bamfile, ploidy))

    # initialize counter
    start_time = time.time()

    # haploid means only variant
    if ploidy=="haploid": merged_sorted_bam = variant_bamfile

    # in all the other cases, do some merging
    else:

        # define the expected files
        merged_bam = "%s.%s.bam"%(variant_bamfile, ploidy)
        merged_sorted_bam = "%s.sorted"%merged_bam
        merged_sorted_bam_index = "%s.bai"%merged_sorted_bam

        if any([file_is_empty(f) for f in [merged_sorted_bam, merged_sorted_bam_index]]) or replace is True:

            # FIRST DEFINE TWO BAMFILES

            # diploid homo is merging a subset of 50% of reads for each variant_bamfile
            if ploidy=="diploid_homo":

                # simulate 50% of reads for each the variant
                bamfileA = downsample_bamfile_keeping_pairs(variant_bamfile, fraction_reads=0.5, replace=replace, threads=threads, name="varGenomeA")
                bamfileB = downsample_bamfile_keeping_pairs(variant_bamfile, fraction_reads=0.5, replace=replace, threads=threads, name="varGenomeB")

            # all the situations where you have a combination of reference and variant genomes
            else:

                # get the fractions of reads
                fraction_var, fraction_ref = get_fractions_reads_for_ploidy(ploidy)

                # simulate the reads for the giiven fraction
                bamfileA = downsample_bamfile_keeping_pairs(variant_bamfile, fraction_reads=fraction_var, replace=replace, threads=threads, name="varGenome")
                bamfileB = downsample_bamfile_keeping_pairs(reference_bamfile, fraction_reads=fraction_ref, replace=replace, threads=threads, name="refGenome")

            # NOW MERGE, SORT and INDEX
            merge_2bams_sort_and_index(bamfileA, bamfileB, merged_bam, merged_sorted_bam, merged_sorted_bam_index, threads=threads)

            # FINALLY REMOVE THE SAMPLED BAMFILES
            for f in [bamfileA, bamfileB]: remove_file(f)

    # get the time the function took to run
    print_if_verbose("--- merging took for ploidy %s reads took %s seconds ---"%(ploidy, time.time() - start_time))

    return merged_sorted_bam

def get_tupleBreakpoints_for_filters_GRIDSS(df_gridss, filters_dict, reference_genome, return_timing=False):

    """ Takes a df_gridss (the output vcf) and a dictionary with filters, returning a tuple of the breakpoints where both breakends have passed the filters."""

    # initialize time
    start_time = time.time()

    # debug the fact that there is no min_af_EitherSmallOrLargeEvent
    if "min_af_EitherSmallOrLargeEvent" not in filters_dict: filters_dict["min_af_EitherSmallOrLargeEvent"] = 0.0

    # get the filtered df
    df_filt = get_gridssDF_filtered(df_gridss, reference_genome, min_Nfragments=filters_dict["min_Nfragments"], min_af=filters_dict["min_af"], wrong_INFOtags=filters_dict["wrong_INFOtags"], wrong_FILTERtags=filters_dict["wrong_FILTERtags"], filter_polyGC=filters_dict["filter_polyGC"], filter_noSplitReads=filters_dict["filter_noSplitReads"], filter_noReadPairs=filters_dict["filter_noReadPairs"], maximum_strand_bias=filters_dict["maximum_strand_bias"], maximum_microhomology=filters_dict["maximum_microhomology"], maximum_lenght_inexactHomology=filters_dict["maximum_lenght_inexactHomology"], range_filt_DEL_breakpoints=filters_dict["range_filt_DEL_breakpoints"], min_length_inversions=filters_dict["min_length_inversions"], dif_between_insert_and_del=filters_dict["dif_between_insert_and_del"], max_to_be_considered_small_event=filters_dict["max_to_be_considered_small_event"], min_size=filters_dict["min_size"], add_columns=False, min_af_EitherSmallOrLargeEvent=filters_dict["min_af_EitherSmallOrLargeEvent"], min_QUAL=filters_dict["min_QUAL"], filter_overlappingRepeats=filters_dict["filter_overlappingRepeats"] )

    # get the breakpoints that have both breakends
    correct_breakpoints = tuple(sorted([bp for bp, N in Counter(df_filt.breakpointID).items() if N==2]))

    if return_timing: return  (time.time() - start_time)
    else: return correct_breakpoints

def keep_relevant_filters_lists_inparallel(filterName_to_filtersList, df_gridss, reference_genome, type_filtering="keeping_all_filters_that_change",  wrong_INFOtags=("IMPRECISE",), min_size=50):

    """Takes a dictionary that maps the filterName to the list of possible filters. It modifies each of the lists in filterName_to_filtersList in a way that only those values that yield a unique set of breakpoints when being applied in the context of a set breakpoints. The final set of filters taken are the less conservative of each. """

    # define a set of filters that are very unconservative (they take all the breakpoints)
    unconservative_filterName_to_filter = {"min_Nfragments":-1, "min_af":-1, "wrong_FILTERtags":("",), "filter_polyGC":False, "filter_noSplitReads":False, "filter_noReadPairs":False, "maximum_strand_bias":1.1, "maximum_microhomology":1000000000000, "maximum_lenght_inexactHomology":1000000000000, "range_filt_DEL_breakpoints":(0,1), "min_length_inversions":-1, "dif_between_insert_and_del":0, "max_to_be_considered_small_event":1, "wrong_INFOtags":wrong_INFOtags, "min_size":min_size, "min_af_EitherSmallOrLargeEvent":-1, "min_QUAL":0, "filter_overlappingRepeats":False}

    # define an unconservative set of breakpoints
    unconservative_breakpoints = tuple(sorted([bp for bp, N in Counter(df_gridss.breakpointID).items() if N==2]))
    print_if_verbose("There are %i bp in total"%len(unconservative_breakpoints))

    # define a list of filters dict, only changing one of the values in filterName_to_filtersList, recording at the same time an ID
    filters_dict_list = [] # a list of filterDicts
    filter_changing_list = [] # a list of (filterName, value) tuples
    for filterName, filtersList in filterName_to_filtersList.items():
        for filterVal in filtersList:

            # define the changed dict in a copy of the unconservative one
            filters_dict = cp.deepcopy(unconservative_filterName_to_filter)
            filters_dict[filterName] = filterVal

            # keep  
            filter_changing_list.append((filterName, filterVal))
            filters_dict_list.append(filters_dict)

    # run in a map or a pool the obtention of tuples of breakpoints for each parameter combination
    inputs_fn = [(df_gridss, fd, reference_genome) for fd in filters_dict_list]

    # pool
    with  multiproc.Pool(multiproc.cpu_count()) as pool:
        bp_tuples_list = pool.starmap(get_tupleBreakpoints_for_filters_GRIDSS, inputs_fn)
        pool.close()

    # map
    #bp_tuples_list = list(map(lambda x: get_tupleBreakpoints_for_filters_GRIDSS(x[0], x[1], x[2]), inputs_fn))

    # map the filter changing to the dict and the breakpoints
    filterChanging_to_filtersDict = dict(zip(filter_changing_list, filters_dict_list))
    filterChanging_to_breakpoints = dict(zip(filter_changing_list, bp_tuples_list))

    # get those filters that change the df
    if type_filtering=="keeping_all_filters_that_change":

        # go through each parameter combination
        for (filterName, filterValue), breakpoints in filterChanging_to_breakpoints.items():

            # if it is the same, remove from 
            if (breakpoints==unconservative_breakpoints or len(breakpoints)==0) and len(filterName_to_filtersList[filterName])>1: filterName_to_filtersList[filterName].remove(filterValue)

    # get those filters that yield a unique set of BPs
    elif type_filtering=="keeping_filters_that_yield_uniqueBPs":

        # go through each filter type, and edit filtersList to keep the ones that yield a unique set of breakpoints
        for filterName, filtersList in filterName_to_filtersList.items():

            # define the filtersChanginIDs for this filter
            filterChanging_IDs = set([ID for ID in filterChanging_to_breakpoints if ID[0]==filterName])

            # get the lists of breakpoints and filtersDicts
            filtersDict_list = [filterChanging_to_filtersDict[ID] for ID in filterChanging_IDs]
            breakpoints_list = [filterChanging_to_breakpoints[ID] for ID in filterChanging_IDs]

            # map each unique breakpoints to a list of the dicts that change it
            bpointTuple_to_filterDicts = {}
            for bpointTuple, filterDict in zip(breakpoints_list, filtersDict_list): bpointTuple_to_filterDicts.setdefault(bpointTuple, []).append(filterDict)

            # for each set, get the less conservative dict, and so extract the value of filterName
            all_filter_values_unique = set([get_represenative_filtersDict_for_filtersDict_list(list_unique_filter_dicts, type_filters="less_conservative")[filterName] for list_unique_filter_dicts in bpointTuple_to_filterDicts.values()])

            # edit filterslist to keep only the unique filter values
            all_non_unique_vals = cp.deepcopy(set(filtersList).difference(all_filter_values_unique))
            for useless_filter_value in all_non_unique_vals: filtersList.remove(useless_filter_value)

    elif "none": pass

    else: raise ValueError("%s is not valid"%type_filtering)

def write_bedpeANDfilterdicts_for_breakpoints(df_bedpe, breakpoints, filterDicts, outdir):

    """Takes a df_bedpe that is already filtered (only has the fields to write) and it writes the breakpoints into outdir. It also writes a series with the less conservative filter set that gace with the filters that have rise to this bedpe"""

    # define the folder
    make_folder(outdir)

    # write bedpe
    bedpe_file = "%s/filtered_breakpoints.bedpe"%outdir
    df_bedpe[df_bedpe.name.isin(breakpoints)].to_csv(bedpe_file, sep="\t", header=False, index=False)

    # get the name of the folder as a file
    outdir_name = outdir.split("/")[-1]

    # get the less conservative filterdict
    less_conservative_filtersDict = get_represenative_filtersDict_for_filtersDict_list(filterDicts, type_filters="less_conservative")
    save_object(less_conservative_filtersDict, "%s/less_conservative_filtersDict.py"%outdir)

def write_breakpoints_for_parameter_combinations_and_get_filterIDtoBpoints_gridss(df_gridss, df_bedpe, outdir, reference_genome, range_filtering="theoretically_meaningful", expected_AF=1.0, replace=False):

    """Gets, for a range of filters defined byrange_filtering, a dictionary that maps a string defining all these filters to a df that has the filtered breakpoints (bedpe) . If range_filtering is large, we expect ~13 h to run on 48 cores for the Candida glabrata genome"""

    # define files that will be written at the end of this function
    filtersID_to_breakpoints_file  = "%s/filtersID_to_breakpoints_file.py"%outdir

    print_if_verbose("getting lists of bedpe breakpoints")

    if any([file_is_empty(f) for f in [filtersID_to_breakpoints_file]]) or replace is True:

        ##### WORK WITH DF_BEDPE ########

        # map each brekendID to the breakpointID
        all_bend_IDs = set.union(*df_bedpe.IDs_set)
        bendID_to_bpointID = {bendID : df_bedpe[df_bedpe.IDs_set.apply(lambda IDs_set: bendID in IDs_set)].iloc[0]["name"] for bendID in df_gridss.ID if bendID in all_bend_IDs}

        # get only breakends that are in bendID_to_bpointID, meaning that there are 2 breakends
        df_gridss_twoBreakEnds = df_gridss[df_gridss.ID.isin(bendID_to_bpointID)]
        df_gridss_twoBreakEnds["breakpointID"] = df_gridss_twoBreakEnds.ID.apply(lambda x: bendID_to_bpointID[x])
        df_gridss_twoBreakEnds = df_gridss_twoBreakEnds.set_index("ID")[["INFO_SIMPLE_TYPE", "length_event", "allele_frequency", "allele_frequency_SmallEvent", "DATA_VF", "INFO_misc", "FILTER", "INFO_SB", "length_microHomology", "length_inexactHomology", "len_inserted_sequence", "has_poly16GC", "DATA_SR", "DATA_RP", "breakpointID", "real_AF", "QUAL", "overlaps_repeats"]]

        # check that all breakpoints have two breakends in the df
        if set(Counter(df_gridss_twoBreakEnds["breakpointID"]).values())!={2}: raise ValueError("Not all breakpoints have 2 brekends")
        
        #################################

        # initialize all filter tags
        all_FILTER_tags = ("ASSEMBLY_ONLY", "NO_ASSEMBLY", "ASSEMBLY_TOO_FEW_READ", "ASSEMBLY_TOO_SHORT", "INSUFFICIENT_SUPPORT", "LOW_QUAL", "REF", "SINGLE_ASSEMBLY")

        # I discard the following filters with these reasons:

        """
        - ASSEMBLY_ONLY is fine, as it is likely that this is enough evidence
        - ASSEMBLY_TOO_FEW_READ. The read correction is already done with the filter min_Nfragments
        - ASSEMBLY_TOO_SHORT. This is already considered in min_size
        - REF: This is already considered in min_af
        - SINGLE_ASSEMBLY: this is already considered given that we filterout breakpoints where one breakend is discarded
        - SINGLE_SUPPORT: already considered from min_Nfragments_l

        In the gridss-purple-linx paper they filter out if there is NO_ASSEMBLY, in the GenomeBiology benchmark they filter if there is LOW_QUAL;NO_ASSEMBLY. Thus, I always filter by NO_ASSEMBLY
        """
        meaningful_FILTER_tags = ("NO_ASSEMBLY", "INSUFFICIENT_SUPPORT", "LOW_QUAL")

        # define arrays of parameters, depending on range_filtering
        if range_filtering=="large":

            min_Nfragments_l = [0, 1, 2, 3, 4, 5, 8, 10, 15, 20, 30]
            min_af_l = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5, expected_AF*0.9]
            min_af_EitherSmallOrLargeEvent_l = min_af_l
            wrong_FILTERtags_l = [("",), ("NO_ASSEMBLY",), ("NO_ASSEMBLY", "INSUFFICIENT_SUPPORT"), ("NO_ASSEMBLY", "LOW_QUAL"), ("LOW_QUAL", "INSUFFICIENT_SUPPORT"), all_FILTER_tags, meaningful_FILTER_tags] 
            filter_polyGC_l = [True, False]
            filter_noSplitReads_l = [True, False]
            filter_noReadPairs_l = [True, False]
            maximum_strand_bias_l = [0.9, 0.95, 0.99, 1.0]
            maximum_microhomology_l = [10, 50, 100, 200, 1000, 100000000]
            maximum_lenght_inexactHomology_l = [10, 50, 100, 200, 1000, 10000000]
            range_filt_DEL_breakpoints_l = [(100, 800), (50, 900), (200, 700), (0,1)]
            min_length_inversions_l = [0, 40, 50, 60, 70, 10000]
            dif_between_insert_and_del_l = [0, 1, 5, 10, 20, 10000000]
            max_to_be_considered_small_event_l = [1, 100, 200, 500, 1000, 1500, 1000000000]
            min_QUAL_l = [0, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1000000000]
            filter_overlappingRepeats_l = [False, True]

        elif range_filtering=="medium":

            min_Nfragments_l = [5, 8, 10]
            min_af_l = [0.05, 0.1, expected_AF*0.9]
            min_af_EitherSmallOrLargeEvent_l = min_af_l
            wrong_FILTERtags_l = [("",), ("NO_ASSEMBLY",), ("NO_ASSEMBLY", "LOW_QUAL"), all_FILTER_tags] 
            filter_polyGC_l = [True]
            filter_noSplitReads_l = [True]
            filter_noReadPairs_l = [True]
            maximum_strand_bias_l = [0.95]
            maximum_microhomology_l = [10, 50, 100]
            maximum_lenght_inexactHomology_l = [10, 50, 100]
            range_filt_DEL_breakpoints_l = [(100, 800), (50, 900), (200, 700)]
            min_length_inversions_l = [40, 50, 60]
            dif_between_insert_and_del_l = [5, 10, 20]
            max_to_be_considered_small_event_l = [200, 500, 1000]
            min_QUAL_l = [0, 50, 100, 500, 900, 1000, 1000000000]
            filter_overlappingRepeats_l = [False, True]

        elif range_filtering=="small":

            min_Nfragments_l = [5, 10]
            min_af_l = [0.05, expected_AF*0.9]
            min_af_EitherSmallOrLargeEvent_l = min_af_l
            wrong_FILTERtags_l = [("",), all_FILTER_tags] 
            filter_polyGC_l = [True]
            filter_noSplitReads_l = [True]
            filter_noReadPairs_l = [True]
            maximum_strand_bias_l = [0.95]
            maximum_microhomology_l = [50]
            maximum_lenght_inexactHomology_l = [50]
            range_filt_DEL_breakpoints_l = [(100, 800)]
            min_length_inversions_l = [40, 50]
            dif_between_insert_and_del_l = [5, 10]
            max_to_be_considered_small_event_l = [1000]
            min_QUAL_l = [0, 500, 1000, 1000000000]
            filter_overlappingRepeats_l = [False, True]


        elif range_filtering=="theoretically_meaningful":

            min_Nfragments_l = [5, 8, 10, 15, 20, 30]
            min_af_l = [0.05, 0.1, 0.2, 0.5, expected_AF*0.9]
            min_af_EitherSmallOrLargeEvent_l = min_af_l
            wrong_FILTERtags_l = [("",), ("NO_ASSEMBLY",), ("NO_ASSEMBLY", "INSUFFICIENT_SUPPORT"), ("NO_ASSEMBLY", "LOW_QUAL"), ("LOW_QUAL", "INSUFFICIENT_SUPPORT"), all_FILTER_tags, meaningful_FILTER_tags] 
            filter_polyGC_l = [True]
            filter_noSplitReads_l = [True, False]
            filter_noReadPairs_l = [True, False]
            maximum_strand_bias_l = [0.9, 0.95]
            maximum_microhomology_l = [100, 200, 1000, 100000000]
            maximum_lenght_inexactHomology_l = [100, 200, 1000, 10000000]
            range_filt_DEL_breakpoints_l = [(100, 800), (50, 900), (200, 700), (0,1)]
            min_length_inversions_l = [0, 50, 1000000000]
            dif_between_insert_and_del_l = [0, 5, 10, 1000000000]
            max_to_be_considered_small_event_l = [100, 200, 500, 1000, 1500, 1000000000]
            min_QUAL_l = [0, 100, 300, 500, 800, 1000, 1000000000]
            filter_overlappingRepeats_l = [False, True]

        elif range_filtering=="single":

            min_Nfragments_l = [8]
            min_af_l = [0.05]
            min_af_EitherSmallOrLargeEvent_l = min_af_l
            wrong_FILTERtags_l = [("NO_ASSEMBLY",)] 
            filter_polyGC_l = [True]
            filter_noSplitReads_l = [True]
            filter_noReadPairs_l = [True]
            maximum_strand_bias_l = [0.95]
            maximum_microhomology_l = [50]
            maximum_lenght_inexactHomology_l = [50]
            range_filt_DEL_breakpoints_l = [(100, 800)]
            min_length_inversions_l = [40]
            dif_between_insert_and_del_l = [5]
            max_to_be_considered_small_event_l = [1000]
            min_QUAL_l = [0]
            filter_overlappingRepeats_l = [False]

        else: raise ValueError("%s is not a valid range_filtering parameter, it has to be 'large', 'medium', 'small' or 'single' "%range_filtering)

        # define filters that are always the same
        wrong_INFOtags = ("IMPRECISE",)
        min_size = 50

        # map the filters through a dict
        filterName_to_filtersList = {"min_Nfragments":min_Nfragments_l, "min_af":min_af_l, "wrong_FILTERtags":wrong_FILTERtags_l, "filter_polyGC":filter_polyGC_l, "filter_noSplitReads":filter_noSplitReads_l, "filter_noReadPairs":filter_noReadPairs_l, "maximum_strand_bias":maximum_strand_bias_l, "maximum_microhomology":maximum_microhomology_l, "maximum_lenght_inexactHomology":maximum_lenght_inexactHomology_l, "range_filt_DEL_breakpoints":range_filt_DEL_breakpoints_l, "min_length_inversions":min_length_inversions_l, "dif_between_insert_and_del":dif_between_insert_and_del_l, "max_to_be_considered_small_event":max_to_be_considered_small_event_l, "min_af_EitherSmallOrLargeEvent":min_af_EitherSmallOrLargeEvent_l, "min_QUAL":min_QUAL_l, "filter_overlappingRepeats":filter_overlappingRepeats_l}

        # edit the filter list, to keep only those that, when applied, change the called breakpoints
        keep_relevant_filters_lists_inparallel(filterName_to_filtersList, df_gridss_twoBreakEnds, reference_genome, type_filtering="keeping_filters_that_yield_uniqueBPs", wrong_INFOtags=wrong_INFOtags, min_size=min_size) # it can also be keeping_all_filters_that_change or keeping_filters_that_yield_uniqueBPs or none

        # initialize objects to store the filtering
        I = 1
        filters_dict_list = []

        # go through each range of filters
        print_if_verbose("generating dictionaries of filters")
        for min_Nfragments in min_Nfragments_l:
          for min_af in min_af_l:
            for min_af_EitherSmallOrLargeEvent in min_af_EitherSmallOrLargeEvent_l:
                for wrong_FILTERtags in wrong_FILTERtags_l:
                  for filter_polyGC in filter_polyGC_l:
                    for filter_noSplitReads in filter_noSplitReads_l:
                      for filter_noReadPairs in filter_noReadPairs_l:
                        for maximum_strand_bias in maximum_strand_bias_l:
                          for maximum_microhomology in maximum_microhomology_l:
                            for maximum_lenght_inexactHomology in maximum_lenght_inexactHomology_l:
                              for range_filt_DEL_breakpoints in range_filt_DEL_breakpoints_l:
                                for min_length_inversions in min_length_inversions_l:
                                  for dif_between_insert_and_del in dif_between_insert_and_del_l:
                                    for max_to_be_considered_small_event in max_to_be_considered_small_event_l:
                                        for min_QUAL in min_QUAL_l:
                                            for filter_overlappingRepeats in filter_overlappingRepeats_l:

                                                I+=1

                                                # get the parameters_dict
                                                filters_dict = dict(min_Nfragments=min_Nfragments, min_af=min_af, wrong_INFOtags=wrong_INFOtags, wrong_FILTERtags=wrong_FILTERtags, filter_polyGC=filter_polyGC, filter_noSplitReads=filter_noSplitReads, filter_noReadPairs=filter_noReadPairs, maximum_strand_bias=maximum_strand_bias, maximum_microhomology=maximum_microhomology, maximum_lenght_inexactHomology=maximum_lenght_inexactHomology, range_filt_DEL_breakpoints=range_filt_DEL_breakpoints, min_length_inversions=min_length_inversions, dif_between_insert_and_del=dif_between_insert_and_del, max_to_be_considered_small_event=max_to_be_considered_small_event, min_size=min_size, min_af_EitherSmallOrLargeEvent=min_af_EitherSmallOrLargeEvent, min_QUAL=min_QUAL, filter_overlappingRepeats=filter_overlappingRepeats)

                                                # keep
                                                filters_dict_list.append(filters_dict)
        
        print_if_verbose("There are %i combinations of parameters"%I)

        # first try for some combinations, which will give you the timing 
        times = [get_tupleBreakpoints_for_filters_GRIDSS(df_gridss_twoBreakEnds, filters_dict, reference_genome, return_timing=True) for filters_dict in random.sample(filters_dict_list, min(10, len(filters_dict_list)))]
        ncores = multiproc.cpu_count()
        print_if_verbose("Obtaining the list of tuples of breakpoints will take arround %.2f minutes on %i cores"%(((np.mean(times)*I)/ncores)/60, ncores))


        # obtain the list of tuples for each parameter combintaion
        with  multiproc.Pool(multiproc.cpu_count()) as pool:
            bp_tuples_list = pool.starmap(get_tupleBreakpoints_for_filters_GRIDSS, [(df_gridss_twoBreakEnds, fd, reference_genome) for fd in filters_dict_list])
            pool.close()

        # map each tuple o bp to the dicts of parameters that gave it
        bpointTuple_to_filterDicts = {}
        for bpointTuple, filterDict in zip(bp_tuples_list, filters_dict_list): 
            if len(bpointTuple)>0: bpointTuple_to_filterDicts.setdefault(bpointTuple, []).append(filterDict)
        print_if_verbose("There are %i sets of breakpoints that can be created with %i combinations of parameters"%(len(bpointTuple_to_filterDicts), I))

        # map each tuple pf breakpoints to an ID that will be saved
        bpoints_to_ID = dict(zip(bpointTuple_to_filterDicts.keys(), map(lambda I: "filters_%i"%I, range(len(bpointTuple_to_filterDicts)))))

        # generate under otdir all the breakpoints from df_bedpe
        print_if_verbose("writing bedpefiles")
        bedpe_fields = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]
        df_bedpe = df_bedpe[bedpe_fields]
        
        # run generation
        inputs_function = [(df_bedpe, bpoints, filterDicts, "%s/%s"%(outdir, bpoints_to_ID[bpoints])) for bpoints, filterDicts in bpointTuple_to_filterDicts.items()]

        with multiproc.Pool(multiproc.cpu_count()) as pool:
            pool.starmap(write_bedpeANDfilterdicts_for_breakpoints, inputs_function)
            pool.close()

        # save the map between each filter 
        print_if_verbose("writing files")
        filtersID_to_breakpoints = dict(zip(bpoints_to_ID.values(), bpoints_to_ID.keys()))
        save_object(filtersID_to_breakpoints, filtersID_to_breakpoints_file)

    else: filtersID_to_breakpoints = load_object(filtersID_to_breakpoints_file)

    # return the dataframe with all the parameter combinations and the filter
    return filtersID_to_breakpoints

def benchmark_bedpe_with_knownSVs(bedpe, know_SV_dict, reference_genome, sorted_bam, median_coverage, replace=False, ID_benchmark="defaultID", delete_intermediate_files=True):

    """Takes the full path to a bedpe file and generates files, under the same directory, that indicate the benchmarking."""

    # write files under the bedpe outdir
    outdir = "/".join(bedpe.split("/")[0:-1])
    bedpe_filename = bedpe.split("/")[-1]

    # get the benchmark file
    benchmark_df_filename = "%s/df_benchmarking_allParms.py"%outdir
    #remove_file(benchmark_df_filename) # debug

    if file_is_empty(benchmark_df_filename) or replace is True:
        print_if_verbose("benchmarking")

        # first run clove without checking for coverage deviations
        outfile_clove = "%s.clove.vcf"%(bedpe)
        run_clove_filtered_bedpe(bedpe, outfile_clove, sorted_bam, replace=replace, median_coverage=10, median_coverage_dev=1, check_coverage=False)

        # now convert it to a df that has also the coverage for TANDEL REGIONS
        df_clove = get_clove_output_with_coverage(outfile_clove, reference_genome, sorted_bam, median_coverage, replace=replace, run_in_parallel=True, delete_bams=False)

        # add an arti

        # initialize benchmark_df
        df_benchmark_all = pd.DataFrame()

        ##### BENCHMARK INSERTIONS, INVERSIONS and TRANSLOCATIONS ##########
        
        # get files in a way that it is similar to RSVSim, only for complex variants
        fileprefix = "%s/insertionsANDinversionsANDtranslocations"%(outdir)
        remaining_df_clove, svtype_to_predsvfile = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, fileprefix, reference_genome, sorted_bam, replace=replace, svtypes_to_consider={"insertions", "inversions", "translocations"})
        svtype_to_predsvfile_InsInvTra = {svtype : svtype_to_predsvfile[svtype] for svtype in {"insertions", "inversions", "translocations"} if svtype in svtype_to_predsvfile}
        know_SV_dict_InsInvTra = {svtype : know_SV_dict[svtype] for svtype in {"insertions", "inversions", "translocations"} if svtype in know_SV_dict}

        # benchmark (and write missing events with overlaps)
        df_benchmark_InsInvTra = benchmark_processedSVs_against_knownSVs_inHouse(svtype_to_predsvfile_InsInvTra, know_SV_dict_InsInvTra, fileprefix, replace=replace, add_integrated_benchmarking=False) # this time you don't want the 'integrated thing'
        
        # keep
        df_benchmark_all = df_benchmark_all.append(df_benchmark_InsInvTra, sort=True)

        # add fields that are necessary to compare with the TANDEL df
        df_benchmark_all["clove_max_rel_coverage_to_consider_del"] = [-1]*len(df_benchmark_all)
        df_benchmark_all["clove_min_rel_coverage_to_consider_dup"] = [-1]*len(df_benchmark_all)


        #########################################################################

        # initialize a benchmark df for tandel
        df_benchmark_TANDEL = pd.DataFrame()

        ##### deletions #################

        # get df
        df_DEL = df_clove[df_clove.SVTYPE=="DEL"]

        # go through different deletion ranges that define true deletion
        max_rel_coverage_to_consider_del_l = [0.0, 0.0001, 0.001, 0.01, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0]
        for max_rel_coverage_to_consider_del in max_rel_coverage_to_consider_del_l:

            # count the maximum coverage for the deletion
            maxDELcoverage = int(max_rel_coverage_to_consider_del*median_coverage)

            # define ID
            coveragefiltID = "maxDELcoverage%i"%(maxDELcoverage)

            # add the filter 
            if len(df_DEL)>0: df_DEL["coverage_FILTER"] = df_DEL.apply(lambda r: get_covfilter_cloveDF_row_according_to_SVTYPE(r, max_rel_coverage_to_consider_del=max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup=0), axis=1)
            else: df_DEL["coverage_FILTER"] = []

            # get a dict svtype_to_svfile
            fileprefix = "%s/%s"%(outdir, coveragefiltID)
            remaining_df_clove, svtype_to_predsvfile = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_DEL, fileprefix, reference_genome, sorted_bam, replace=replace, svtypes_to_consider={"deletions"})

            # benchmark
            know_SV_dict_DEL = {svtype : know_SV_dict[svtype] for svtype in {"deletions"} if svtype in know_SV_dict}

            df_benchmark = benchmark_processedSVs_against_knownSVs_inHouse(svtype_to_predsvfile, know_SV_dict_DEL, fileprefix, replace=replace, add_integrated_benchmarking=False)

            df_benchmark["clove_max_rel_coverage_to_consider_del"] = [max_rel_coverage_to_consider_del]*len(df_benchmark)
            df_benchmark["clove_min_rel_coverage_to_consider_dup"] = [-1]*len(df_benchmark)

            # keep
            df_benchmark_TANDEL = df_benchmark_TANDEL.append(df_benchmark, sort=True)

        ################################

        ##### tandem duplications ######

        # get df
        df_TAN = df_clove[df_clove.SVTYPE=="TAN"]

        # go through different TAN ranges that define true TAN
        min_rel_coverage_to_consider_dup_l = [0.0, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.5]
        for min_rel_coverage_to_consider_dup in min_rel_coverage_to_consider_dup_l:

            # count the maximum coverage for the tan
            minDUPcoverage = int(min_rel_coverage_to_consider_dup*median_coverage)

            # define ID
            coveragefiltID = "minDUPcoverage%i"%(minDUPcoverage)

            # add the filter 
            if len(df_TAN)>0: df_TAN["coverage_FILTER"] = df_TAN.apply(lambda r: get_covfilter_cloveDF_row_according_to_SVTYPE(r, max_rel_coverage_to_consider_del=1000000, min_rel_coverage_to_consider_dup=min_rel_coverage_to_consider_dup), axis=1)
            else: df_TAN["coverage_FILTER"] = []

            # get a dict svtype_to_svfile
            fileprefix = "%s/%s"%(outdir, coveragefiltID)
            remaining_df_clove, svtype_to_predsvfile = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_TAN, fileprefix, reference_genome, sorted_bam, replace=replace, svtypes_to_consider={"tandemDuplications"})

            # benchmark
            know_SV_dict_TAN = {svtype : know_SV_dict[svtype] for svtype in {"tandemDuplications"} if svtype in know_SV_dict}

            df_benchmark = benchmark_processedSVs_against_knownSVs_inHouse(svtype_to_predsvfile, know_SV_dict_TAN, fileprefix, replace=replace, add_integrated_benchmarking=False)
            df_benchmark["clove_max_rel_coverage_to_consider_del"] = [-1]*len(df_benchmark)
            df_benchmark["clove_min_rel_coverage_to_consider_dup"] = [min_rel_coverage_to_consider_dup]*len(df_benchmark)

            # keep
            df_benchmark_TANDEL = df_benchmark_TANDEL.append(df_benchmark, sort=True)

        ################################

        # keep
        df_benchmark_all = df_benchmark_all.append(df_benchmark_TANDEL, sort=True)

        # append both and return
        df_benchmark_all["benchmarkID"] = [ID_benchmark]*len(df_benchmark_all)
        df_benchmark_all["bedpe"] = [bedpe]*len(df_benchmark_all) # I changed this line at some point because it was raising integers

        # add the parameters that yielded this dict
        filters_dict = load_object("%s/less_conservative_filtersDict.py"%outdir)
        df_benchmark_all["filters_dict"] = [filters_dict]*len(df_benchmark_all)

        # save in disk
        print_if_verbose("saving into %s"%benchmark_df_filename)
        save_object(df_benchmark_all, benchmark_df_filename)

    else: df_benchmark_all = load_object(benchmark_df_filename)

    # delete intermediate fields
    if delete_intermediate_files is True:

        filenames_to_keep = {bedpe_filename, "%s.clove.vcf.TANDEL.bed.coverage_provided_windows.tab"%bedpe_filename, "unbalanced_translocations_5with5_or_3with3.bed.coverage_provided_windows.tab", "%s.clove.vcf"%bedpe_filename, "uniform_filters_series.py", "variable_filters_df.py", "df_benchmarking_allParms.py", "less_conservative_filtersDict.py"}

        #for file in os.listdir(outdir):
        #    if file not in filenames_to_keep and "benchmark_analysis_" not in file: remove_file("%s/%s"%(outdir, file))


    return df_benchmark_all



def merge_tables_into_file(list_table_files, outfile):

    """Takes a list of table files and merges them into outfile"""
    
    if len(list_table_files)>0:

        df = pd.concat([pd.read_csv(x, sep="\t") for x in list_table_files if os.path.isfile(x)], sort=True)

        # check that the outilfe is not in the list_table_files
        if outfile in list_table_files: outfile += ".%s"%(id_generator(25))
        df.to_csv(outfile, sep="\t", header=True, index=False)

        # remove previously generated files
        for f in list_table_files: remove_file(f)


def makePlots_gridsss_benchmarking_oneGenome(df_benchmark, PlotsDir, plots={"histogram", "scatter_PRvsRC", "scatter_PRvsRCa_eachSVtype", "Fscore_correlation_scatter", "Fscore_correlation_mat"}):

    """Takes a dataframe such as the output of benchmark_GridssClove_for_knownSV and writes plots under PlotsDir, as specified in plots. These are several """

    print_if_verbose("performing plots into %s"%PlotsDir)

    make_folder(PlotsDir)

    # map each svtype to a marker
    svtype_to_color = {"tandemDuplications": "gray", "deletions": "black", "inversions": "blue", "translocations": "olive", "insertions": "red", "remaining":"magenta"}
    svtype_to_marker = {"tandemDuplications": "P", "deletions": "s", "inversions": "^", "translocations": "D", "insertions": "o", "remaining":"v"}

    # define things
    all_events = {'inversions', 'translocations', 'deletions', 'tandemDuplications', 'insertions'}

    # only consider the benchmarking types if there are events
    all_events = [e for e in {'inversions', 'translocations', 'deletions', 'tandemDuplications', 'insertions', 'remaining'} if len(df_benchmark[df_benchmark.svtype==e])>0]

    # change the vals to floats
    for field in ["precision", "recall", "Fvalue"]: df_benchmark[field] = df_benchmark[field].apply(float)

    # a histogram with the distribution of precision, recall and Fvalue
    if "histogram" in plots:

        fig = plt.figure(figsize=(12, 4))

        # go through each field to plot
        for I, field in enumerate(["precision", "recall", "Fvalue"]):

            # initialize subplot
            ax = plt.subplot(1, 3, I+1)

            # add hists for each svtype
            for svtype in set(df_benchmark.svtype): sns.distplot(list(df_benchmark[df_benchmark.svtype==svtype][field]), hist=True, kde=False, rug=True, color=svtype_to_color[svtype], label=svtype, kde_kws=dict(linewidth=3))

            ax.set_xlabel(field)
            ax.set_ylabel("n parameter combinations")


        fig.tight_layout()  # otherwise the right y-label is slightly 
        filename="%s/histogram_PR_RC_Fvalue.pdf"%(PlotsDir)
        fig.savefig(filename, bbox_inches='tight');
        #if is_cluster is False: plt.show()
        plt.close(fig)

    # a scatterplot correlating precision an recall
    if "scatter_PRvsRC" in plots:

        fig = plt.figure(figsize=(12, 4))

        # first subplot, raw values
        ax = plt.subplot(1, 3, 1)
        sns.scatterplot(x="recall", y="precision", data=df_benchmark, hue="svtype", palette=svtype_to_color)

        # raw with alpha
        ax = plt.subplot(1, 3, 2)
        sns.scatterplot(x="recall", y="precision", data=df_benchmark, hue="svtype", palette=svtype_to_color, alpha=0.15, edgecolors=None)

        # color according to Fvalue
        ax = plt.subplot(1, 3, 3)
        cmap = sns.cubehelix_palette(dark=.3, light=.8, as_cmap=True)

        sns.scatterplot(x="recall", y="precision", data=df_benchmark, hue="Fvalue", palette=cmap, edgecolors=None, style="svtype", markers=svtype_to_marker)
        #sns.scatterplot(x="recall", y="precision", data=df_benchmark, hue="Fvalue", edgecolors=None, style="svtype", markers=svtype_to_marker)


        fig.tight_layout()  # otherwise the right y-label is slightly 
        filename="%s/scatter_PRvsRCvsFvalue.pdf"%(PlotsDir)
        fig.savefig(filename, bbox_inches='tight');
        #if is_cluster is False: plt.show()
        plt.close(fig)

    # a PRvsRC plot for each svtype
    if "scatter_PRvsRCa_eachSVtype" in plots:

        all_svtypes = sorted(set(df_benchmark.svtype))

        fig = plt.figure(figsize=(4*len(all_svtypes), 4))

        for I, svtype in enumerate(all_svtypes):

            # color according to Fvalue
            ax = plt.subplot(1, len(all_svtypes), I+1)
            cmap = sns.cubehelix_palette(dark=.3, light=.6, as_cmap=True)

            # this sometimes fails with the Fvalue as hue
            try: sns.scatterplot(x="recall", y="precision", data=df_benchmark[df_benchmark.svtype==svtype], hue="Fvalue", palette=cmap, edgecolors=None, style="svtype", markers=svtype_to_marker)
            except: sns.scatterplot(x="recall", y="precision", data=df_benchmark[df_benchmark.svtype==svtype], edgecolors=None, style="svtype", markers=svtype_to_marker)

            plt.axvline(1, color="black", linestyle="--", linewidth=1.0)
            plt.axhline(1, color="black", linestyle="--", linewidth=1.0)
            plt.plot([0, 1], color="black", linestyle="--", linewidth=1.0)
            

            ax.set_title(svtype)

        fig.tight_layout()  # otherwise the right y-label is slightly 
        filename="%s/scatter_PRvsRCa_eachSVtype.pdf"%(PlotsDir)
        fig.savefig(filename, bbox_inches='tight');
        #if is_cluster is False: plt.show()
        plt.close(fig)


    # a scatterplot correlating, the Fscores of each pairwise comparison of SVTYPE, for this you need to define a dataframe that has a unique identifier for each unique [svtype, gridss_maxcoverage, gridss_regionsToIgnoreBed], and taking the maxiumum Fscore

    filter_set_ID = ["gridss_maxcoverage", "gridss_regionsToIgnoreBed", "benchmarkID"]
    unique_set_ID = filter_set_ID + ["svtype"]
    df_Fvalue = df_benchmark[unique_set_ID + ["Fvalue"]].groupby(unique_set_ID).max() # this has index the unique_set_ID
    df_Fvalue["svtype"] = [x[3] for x in df_Fvalue.index]
    df_Fvalue["filterID"] = [tuple(x[0:3]) for x in df_Fvalue.index]

    # define, for each svtype combination, the Fscores correlated
    svX_to_svY_to_vals = {}
    svX_to_svY_to_SameMaxAccuracyFilters = {}

    all_svtypes = list(set(df_benchmark.svtype))
    for Ix, svtype_X in enumerate(all_svtypes):
        for svtype_Y in all_svtypes[Ix:]:

            # same idx
            if svtype_Y==svtype_X: 

                svX_to_svY_to_SameMaxAccuracyFilters.setdefault(svtype_X, {}).setdefault(svtype_Y, "")
                svX_to_svY_to_SameMaxAccuracyFilters.setdefault(svtype_Y, {}).setdefault(svtype_X, "")
                continue

            # get the Fscores X and Y, corresponding to the same filters
            svXfilt_to_values = dict(zip(df_Fvalue[df_Fvalue.svtype==svtype_X].filterID, df_Fvalue[df_Fvalue.svtype==svtype_X].Fvalue))
            svYfilt_to_values = dict(zip(df_Fvalue[df_Fvalue.svtype==svtype_Y].filterID, df_Fvalue[df_Fvalue.svtype==svtype_Y].Fvalue))

            svX_to_svY_to_vals.setdefault(svtype_X, {}).setdefault(svtype_Y, pd.DataFrame([[svXfilt_to_values[filtID], svYfilt_to_values[filtID]] for filtID in svYfilt_to_values], columns=["X", "Y"]))

            # record whether the filter ID of the X and the Y are the same
            maxXfilt = max(svXfilt_to_values.items(), key=(lambda x:x[1]))[0]
            maxYfilt = max(svYfilt_to_values.items(), key=(lambda x:x[1]))[0]

            sameMaxFilt_to_str = {True:"*", False:""}
            svX_to_svY_to_SameMaxAccuracyFilters.setdefault(svtype_X, {}).setdefault(svtype_Y, sameMaxFilt_to_str[maxXfilt==maxYfilt])
            svX_to_svY_to_SameMaxAccuracyFilters.setdefault(svtype_Y, {}).setdefault(svtype_X, sameMaxFilt_to_str[maxXfilt==maxYfilt])

    # scatterplot
    if "Fscore_correlation_scatter" in plots:

        fig = plt.figure(figsize=(7, 7))

        # go through each combination and plot
        for svX, svY_to_vals in svX_to_svY_to_vals.items():
            for svY, vals in svY_to_vals.items():

                corr, p = scipy.stats.spearmanr(vals, axis=0)
                if pd.isna(corr): corr = 0.0; p = 1.0
                ax = sns.scatterplot(x="X", y="Y", data=vals, color=svtype_to_color[svX], marker=svtype_to_marker[svY], label="X:%s, Y:%s, r=%.2f, p=%.3f"%(svX, svY, corr, p))

        ax.set_xlabel("Fvalue X")
        ax.set_ylabel("Fvalue Y")


        fig.tight_layout()  # otherwise the right y-label is slightly 
        filename="%s/Fscore_correlation_scatter.pdf"%(PlotsDir)
        fig.savefig(filename, bbox_inches='tight');
        #if is_cluster is False: plt.show()
        plt.close(fig)

    # s
    if "Fscore_correlation_mat" in plots:

        fig = plt.figure(figsize=(4, 4))

        svX_to_svY_to_corr = {}
        svX_to_svY_to_p = {}

        # go through each combination and get the correlation data
        for svX, svY_to_vals in svX_to_svY_to_vals.items():
            for svY, vals in svY_to_vals.items():

                corr, p = scipy.stats.spearmanr(vals, axis=0)
                if pd.isna(corr): corr = 0.0; p = 1.0

                # correlation
                svX_to_svY_to_corr.setdefault(svX, {}).setdefault(svY, corr)
                svX_to_svY_to_corr.setdefault(svY, {}).setdefault(svX, corr)

        # get as dfs
        df_corr = pd.DataFrame(svX_to_svY_to_corr); df_corr = df_corr[list(df_corr.index)]
        df_labels = pd.DataFrame(svX_to_svY_to_SameMaxAccuracyFilters)[list(df_corr.index)].loc[list(df_corr.index)].applymap(str)

        # Generate a mask for the upper triangle
        mask = np.zeros_like(df_corr, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True

        # Generate a custom diverging colormap
        cmap = sns.cubehelix_palette(dark=.3, light=.6, as_cmap=True)

        # Draw the heatmap with the mask and correct aspect ratio
        sns.heatmap(df_corr, cmap=cmap, square=True, linewidths=.5, cbar_kws={"shrink": .5, "label":"correlation, *: same best filtering"}, mask=mask, annot=df_labels, fmt = '')

        fig.tight_layout()  # otherwise the right y-label is slightly 
        filename="%s/Fscore_correlation_mat.pdf"%(PlotsDir)
        fig.savefig(filename, bbox_inches='tight');
        #if is_cluster is False: plt.show()
        plt.close(fig)

def benchmark_GridssClove_for_knownSV(sample_bam, reference_genome, know_SV_dict, outdir, range_filtering="theoretically_meaningful", expected_AF=1.0, replace=False, threads=4, median_insert_size=500, median_insert_size_sd=50, mitochondrial_chromosome="mito_C_glabrata_CBS138", run_in_parallel=False):

    """Runs a benchmarking for several combinations of filters of a GridsssClove pipeline of a given bam file (sample_bam), writing files under outdir. The known SV are provided as a dictionary that maps each type of SV to a path where a table with the SVs are known .

    range_filtering indicates which type of simulation will be performed, it can be "large", "medium", "small", "single" and correlates with the range of parameters to use.
    expected_AF is the expected allele frequency, for a haploid it should be 1.0 

    median_insert_size is used to define small breakends to calculate their allele frequency

    window_l is used to define the coverage for winows of regions"""


    ###### DEFINE GENERAL THINGS

    start_time = time.time()

    # define the median coverage of regions
    print_if_verbose("getting coverage")
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir, sample_bam, windows_file="none", replace=replace), sep="\t")
    median_coverage = np.median(coverage_df[~coverage_df["#chrom"].isin(mitochondrial_chromosome.split(","))].mediancov_1); print_if_verbose("The median coverage is %i"%median_coverage)

    ##################
    
    # define the combinations of gridss parameters
    maxcoverage_list = [50000]
    ignore_regions_list = [False]

    # initialize a df with all the benchmarking
    all_benchmarking_df = pd.DataFrame()

    # run gridss in several different ways:
    I = 1
    for maxcoverage in maxcoverage_list:
        for ignore_regions in ignore_regions_list:

            # make an outdir for this gridss run
            id_gridss_run = "max%ix_ignoreRegions%s"%(maxcoverage, ignore_regions)
            gridss_outdir = "%s/benchmark_%s"%(outdir, id_gridss_run); make_folder(gridss_outdir)

            # define the black listed regins
            if ignore_regions is True: raise ValueError("This has not been configured yet to blacklist regions from gridss running")
            else: blacklisted_regions = ""

            # get the gridss outputs
            gridss_VCFoutput = run_gridss_and_annotateSimpleType(sample_bam, reference_genome, gridss_outdir, replace=replace, threads=threads, blacklisted_regions=blacklisted_regions, maxcoverage=maxcoverage)
                  
            bedpe_with_adds = get_bedpe_from_svVCF(gridss_VCFoutput, gridss_outdir, replace=replace)

            # get into dfs with generally interesting info
            df_gridss = add_info_to_gridssDF(load_single_sample_VCF(gridss_VCFoutput), reference_genome, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd) # this is a dataframe with some extra info
            df_bedpe = pd.read_csv(bedpe_with_adds, sep="\t")
            df_bedpe["IDs_set"] = df_bedpe.IDs.apply(lambda x: set(x.split("||")))

            # write the breakpoints. The point of this is that with many parameter combinations we may yield the same breakpoints, so that it's more efficient to create them first
            outdir_parameter_combinations = "%s/several_parameter_combinations_filter_%s_af%.2f"%(gridss_outdir, range_filtering, expected_AF)
            #delete_folder(outdir_parameter_combinations) # DEBUG
            make_folder(outdir_parameter_combinations)
            filtersID_to_breakpoints = write_breakpoints_for_parameter_combinations_and_get_filterIDtoBpoints_gridss(df_gridss, df_bedpe, outdir_parameter_combinations, reference_genome, range_filtering=range_filtering, expected_AF=expected_AF, replace=replace) # this is a dataframe with all the filter combinations and the map between filterID and the actual filtering

            # define the paths to the breakpoints
            paths_to_bedpe_breakpoints = ["%s/%s/filtered_breakpoints.bedpe"%(outdir_parameter_combinations, filterID) for filterID in filtersID_to_breakpoints]

            # define inputs of the benchmarking pipeline
            inputs_benchmarking_pipeline = [(bedpe, know_SV_dict, reference_genome, sample_bam, median_coverage, replace, bedpe.split("/")[-2], True) for bedpe in paths_to_bedpe_breakpoints]

            if run_in_parallel is True:

                # initialize the list of benchmarking dfs
                all_benchmarking_dfs = []

                # go through each chunk of ncpus
                for Ichunk, chunk_inputs_benchmarking_pipeline in enumerate(chunks(inputs_benchmarking_pipeline, threads)):
                    print_if_verbose("working on chunk %i"%Ichunk)

                    # get the parallelized obtention of data
                    print_if_verbose("getting benchmarking for each set of filters in parallel")
                    with multiproc.Pool(threads) as pool:
                        all_benchmarking_dfs += pool.starmap(benchmark_bedpe_with_knownSVs, chunk_inputs_benchmarking_pipeline)
                        pool.close()
                        pool.terminate()
            else:

                all_benchmarking_dfs = list(map(lambda x: benchmark_bedpe_with_knownSVs(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]), inputs_benchmarking_pipeline))

            #### delete the bamfiles for the chromosomes, which are created in benchmark_bedpe_with_knownSVs ######
            all_chromosome_IDs = set([s.id for s in SeqIO.parse(reference_genome, "fasta")])
            outdir_bam = get_dir(sample_bam)
            name_bam = get_file(sample_bam)
            for chrom in all_chromosome_IDs: 
                files_to_remove = ["%s/%s"%(outdir_bam, x) for x in os.listdir(outdir_bam) if x.startswith(name_bam) and chrom in x]
                for f in files_to_remove: remove_file(f)
            #######################################################################################################

            ###### remove temporary files and folders #####
            for f in ["%s/%s"%(outdir_bam, x) for x in os.listdir(outdir_bam) if "temporary_file" in x]: remove_file(f); delete_folder(f)

            ##### merge the tables of coverage regions and remove intermediate files ######
            tables_coverage_windows = ["%s/%s"%(outdir_bam, x) for x in os.listdir(outdir_bam) if "temporary_file" not in x and "%s.coverage_per_window.tab"%name_bam in x]
            outfile = "%s.coverage_per_window.tab.%s"%(sample_bam, get_date().replace("/", "_"))
            merge_tables_into_file(tables_coverage_windows, outfile)

            # merge all the coverage files generated in one
            merge_coverage_per_window_files_in_one(sample_bam)

            #################################################

            # define the fields
            fields_benchmarking_df = list(all_benchmarking_dfs[0].keys())

            # concatenate and add fields
            print_if_verbose("concatenating dfs")
            benchmarking_df = pd.concat([d[fields_benchmarking_df] for d in all_benchmarking_dfs], sort=True)
            benchmarking_df["gridss_maxcoverage"] = [maxcoverage]*len(benchmarking_df)
            benchmarking_df["gridss_regionsToIgnoreBed"] = [blacklisted_regions]*len(benchmarking_df)
            benchmarking_df["gridss_VCFoutput"] = [gridss_VCFoutput]*len(benchmarking_df)
            benchmarking_df["median_coverage"] =  [median_coverage]*len(benchmarking_df)

            # keep
            if I==1: all_benchmarking_df = all_benchmarking_df.append(benchmarking_df, sort=True)
            else: all_benchmarking_df = all_benchmarking_df[list(benchmarking_df.keys())].append(benchmarking_df, sort=True)

            I+=1

    # debug strange calculations
    if any(pd.isna(all_benchmarking_df.precision)) or any(pd.isna(all_benchmarking_df.recall)): raise ValueError("There are NaNs in the precision or recall measures") 

    # make plots of the benchmarking
    PlotsDir = "%s/plots_benchmark"%outdir
    print_if_verbose("making plots into %s"%PlotsDir)
    makePlots_gridsss_benchmarking_oneGenome(all_benchmarking_df, PlotsDir)

    # get the time
    print_if_verbose("----It took %s seconds to run the whole benchmarking of one set of SV----"%(time.time() - start_time))

    return all_benchmarking_df

def get_df_accuracy_for_train_filer(r, outdir, test_gridss_info_dict, sorted_bam, reference_genome, median_coverage, replace, median_insert_size, median_insert_size_sd, test_SVdict):

    """define a function that takes a row of df_filters_train and returns a series with the accuracy values of each filter"""

    # define outdir
    working_dir = "%s/train_on_%s_%s_%s"%(outdir, r["genomeID"], r["ploidy"], r["svtype"]); make_folder(working_dir)

    # define the file
    df_benchmark_filename = "%s/df_benchmark.tab"%working_dir

    if file_is_empty(df_benchmark_filename) or replace is True:

        # define the gridss_VCFoutput based on test_gridss_info_dict_under_outdir
        gridss_VCFoutput = test_gridss_info_dict[r["gridss_regionsToIgnoreBed"]][r["gridss_maxcoverage"]]

        # make a link under working_dir
        gridss_VCFoutput_underWorkDir = "%s/gridss_output.vcf"%(working_dir)
        print_if_verbose("testing...", gridss_VCFoutput_underWorkDir)
        if file_is_empty(gridss_VCFoutput_underWorkDir) or replace is True: soft_link_files(gridss_VCFoutput, gridss_VCFoutput_underWorkDir)

        # get the svs
        predicted_svtype_to_SVtable, df_gridss = run_gridssClove_given_filters(sorted_bam, reference_genome, working_dir, median_coverage, replace=replace, threads=multiproc.cpu_count(), gridss_blacklisted_regions=r["gridss_regionsToIgnoreBed"], gridss_VCFoutput=gridss_VCFoutput_underWorkDir, gridss_maxcoverage=r["gridss_maxcoverage"], median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, gridss_filters_dict=r["filters_dict"], tol_bp=50, run_in_parallel=True, max_rel_coverage_to_consider_del=r["clove_max_rel_coverage_to_consider_del"], min_rel_coverage_to_consider_dup=r["clove_min_rel_coverage_to_consider_dup"], replace_FromGridssRun=False)

        # get the benchmarking df
        fileprefix = "%s.benchmarking"%working_dir
        df_benchmark_filtN = benchmark_processedSVs_against_knownSVs_inHouse(predicted_svtype_to_SVtable, test_SVdict, fileprefix, replace=replace, add_integrated_benchmarking=True)

        # add the metdadata
        df_benchmark_filtN["train_genomeID"] = r["genomeID"]
        df_benchmark_filtN["train_ploidy"] = r["ploidy"]
        df_benchmark_filtN["train_svtype"] = r["svtype"]

        # save
        df_benchmark_filtN.to_csv(df_benchmark_filename, sep="\t", header=True, index=False)

    else: df_benchmark_filtN = pd.read_csv(df_benchmark_filename, sep="\t")

    return df_benchmark_filtN

def get_benchmarking_df_for_testSVs_from_trainSV_filterSets(test_SVdict, outdir, df_filters_train, test_gridss_info_dict, genomeID, ploidy, sorted_bam, reference_genome, median_coverage, median_insert_size, median_insert_size_sd, replace):

    """This function takes a  set of test SVdict and it tries all the filters in df_filters_train on the gridss vcfs in gridss_info_dict, writing files under outdir. It returns a df with the accuracy for each svtype and filters from df_filters_train tested on SVdict. This will only train on 'integrated' and test on all SVtypes  """

    start_time = time.time()

    # check that the df_filters_train contains unique vals for each genomeID, ploidy and svtype
    if len(df_filters_train)!=len(df_filters_train[["genomeID", "ploidy", "svtype"]].drop_duplicates()): raise ValueError('df_filters_train does not contain unique vals for "genomeID", "ploidy", "svtype"')

    # define the df_benchmark
    df_benchmark_all_filename = "%s/df_benchmark_all.tab"%outdir
    print_if_verbose("working on %s"%df_benchmark_all_filename)

    if file_is_empty(df_benchmark_all_filename) or replace is True:

        # keep only the integrated train set. If this is commented it does not happen
        #df_filters_train = df_filters_train[df_filters_train.svtype=="integrated"]

        df_benchmark = pd.concat(list(df_filters_train.apply(lambda r: get_df_accuracy_for_train_filer(r, outdir, test_gridss_info_dict, sorted_bam, reference_genome, median_coverage, replace, median_insert_size, median_insert_size_sd, test_SVdict), axis=1)))

        # add metadata
        df_benchmark["test_genomeID"] = genomeID
        df_benchmark["test_ploidy"] = ploidy
        df_benchmark["test_svtype"] = df_benchmark.svtype

        # save
        print_if_verbose("saving %s"%df_benchmark_all_filename)
        df_benchmark.to_csv(df_benchmark_all_filename, sep="\t", header=True, index=False)

    else: df_benchmark = pd.read_csv(df_benchmark_all_filename, sep="\t")

    print_if_verbose("----It took %s seconds to run the whole benchmarking of one set of test filters----"%(time.time() - start_time))

    return df_benchmark


######################################################
################# GRAPHICS FUNCTIONS #################
######################################################

def plot_clustermap_with_annotation(df, row_colors_df, col_colors_df, filename, title="clustermap", col_cluster=False, row_cluster=False, colorbar_label="default label", adjust_position=True, legend=True, idxs_separator_pattern="_", texts_to_strip={"L001"}, default_label_legend="control", df_annotations=None, cmap=sns.color_palette("RdBu_r", 50), ylabels_graphics_df=None, grid_lines=True, add_to_legend_x=1, figsize=None):

    """Takes a df were the index is the annotation and the cols are samples. It will be saved under filename. ylabels_graphics_df can be a df containing fontweight and color for each index value in df"""

    # define the yticklabels
    if ylabels_graphics_df is not None: yticklabels = True
    else: yticklabels = False

    # change the values of the df to floats
    df = df.applymap(float)

    # get the ordered xlabels and ylabels
    cm_labels = sns.clustermap(df, col_cluster=col_cluster, row_cluster=row_cluster, yticklabels=True, xticklabels=True)
    # plt.close()
    xlabels = [x.get_text() for x in cm_labels.ax_heatmap.get_xticklabels()]
    ylabels = [y.get_text() for y in cm_labels.ax_heatmap.get_yticklabels()]

    # decide whether to add annotations
    if df_annotations is not None:  annot = df_annotations.loc[ylabels][xlabels]
    else: annot = False

    # define the line color from the last item in the df
    if grid_lines is True: 
        linecolor = list(col_colors_df.loc[xlabels][list(col_colors_df.keys())[-1]])
        linewidths = 1.5

    else: 
        linecolor = "gray"
        linewidths = 0

    if adjust_position is True:

        # define figsize 
        figsize = (len(df.columns)*0.3, len(df)*0.35)

        # get the clustermap
        cm = sns.clustermap(df, col_cluster=col_cluster, row_cluster=row_cluster, row_colors=row_colors_df, col_colors=col_colors_df, cbar_kws={'label': colorbar_label}, xticklabels=False, square=False, figsize=figsize, cmap=cmap, annot=annot, fmt="", annot_kws={"size": 6.5}, linecolor=linecolor, linewidths=linewidths, yticklabels=yticklabels) # figsize=figsize, linecolor=linecolor, linewidths=linewidths, yticklabels=yticklabels

        # move the heatmap to the right
        hm_pos = cm.ax_heatmap.get_position()
        cm.ax_heatmap.set_position([hm_pos.x0, hm_pos.y0, hm_pos.width, hm_pos.height]); hm_pos = cm.ax_heatmap.get_position()

        # adjust the row colorbar, proportional to the row colorbar
        width_row_colorbar = (hm_pos.width/len(col_colors_df)) * len(row_colors_df.columns)
        rc_pos = cm.ax_row_colors.get_position()
        cm.ax_row_colors.set_position([hm_pos.x0 - width_row_colorbar, rc_pos.y0, width_row_colorbar, rc_pos.height]); rc_pos = cm.ax_row_colors.get_position()

        # adjust the col colorbar proporitonal to the col colors
        height_col_colorbar = (hm_pos.height/len(df)) * len(col_colors_df.columns)
        cc_pos = cm.ax_col_colors.get_position()
        cm.ax_col_colors.set_position([hm_pos.x0, hm_pos.y0 + hm_pos.height, hm_pos.width, height_col_colorbar]); cc_pos = cm.ax_col_colors.get_position()

        # adjust the row dendrogram
        if row_cluster is True: width_row_dendrogram = (hm_pos.width/len(df.columns))*4
        else: width_row_dendrogram = 0
        cm.ax_row_dendrogram.set_position([rc_pos.x0-width_row_dendrogram-0.001, rc_pos.y0, width_row_dendrogram, rc_pos.height])
        row_dendrogram_pos = cm.ax_row_dendrogram.get_position()

        
        # adjust the position of the colorbar to the left of the col dednrogram
        
        rdendro_pos = cm.ax_row_dendrogram.get_position()
        height_colorbar = hm_pos.height/2
        width_colorbar = rc_pos.width/2
        cm.cax.set_position([rdendro_pos.x0 - width_row_dendrogram - hm_pos.width*0.5, rdendro_pos.y0 + (rdendro_pos.height-rdendro_pos.y0)/2 - height_colorbar/2, width_colorbar, height_colorbar])
        
        # adjust position of the col dendrogram
        cdendro_pos = cm.ax_col_dendrogram.get_position()
        if col_cluster is True: height_col_dendrogram = (hm_pos.height/len(df))*4
        else: height_col_dendrogram = 0
        cm.ax_col_dendrogram.set_position([cc_pos.x0, cc_pos.y0+cc_pos.height, cc_pos.width, height_col_dendrogram])
    
    else:

        # define figsize 
        #figsize = (len(df.columns)*0.3, len(df)*0.35)

        # get the clustermap
        cm = sns.clustermap(df, col_cluster=col_cluster, row_cluster=row_cluster, row_colors=row_colors_df, col_colors=col_colors_df, cbar_kws={'label': colorbar_label}, xticklabels=False, square=True, cmap=cmap, annot=annot, fmt="", annot_kws={"size": 6}, linecolor=linecolor, linewidths=linewidths, yticklabels=yticklabels, figsize=figsize) # figsize=figsize, linecolor=linecolor, 


        ########### adjsut square position to the hm ###########

        # define the size of each square
        hm_pos = cm.ax_heatmap.get_position()
        size_square = hm_pos.width / len(df.columns)

        # adjust the col colorbar
        height_col_colorbar = size_square * len(col_colors_df.columns)
        cm.ax_col_colors.set_position([hm_pos.x0, hm_pos.y0 + hm_pos.height + size_square*0.2, hm_pos.width, height_col_colorbar]); cc_pos = cm.ax_col_colors.get_position()

        # adjust the row colorbar
        width_row_colorbar = size_square * len(row_colors_df.columns)
        cm.ax_row_colors.set_position([hm_pos.x0 - width_row_colorbar - size_square*0.2, hm_pos.y0, width_row_colorbar, hm_pos.height]); rc_pos = cm.ax_row_colors.get_position()

        # adjust the row dendrogram
        if row_cluster is True: width_row_dendrogram = width_row_colorbar
        else: width_row_dendrogram = 0
        cm.ax_row_dendrogram.set_position([rc_pos.x0-width_row_dendrogram, rc_pos.y0, width_row_dendrogram, rc_pos.height])
        row_dendrogram_pos = cm.ax_row_dendrogram.get_position()


        # adjust position of the col dendrogram
        if col_cluster is True: height_col_dendrogram = height_col_colorbar
        else: height_col_dendrogram = 0
        cm.ax_col_dendrogram.set_position([cc_pos.x0, cc_pos.y0+cc_pos.height, cc_pos.width, height_col_dendrogram])

        ########################################################

    # get the tile
    cm.ax_col_dendrogram.set_title(title)

    # remove the colors from the colorbars
    cm.ax_col_colors.collections[0].set_linewidth(0.1)
    cm.ax_col_colors.collections[0].set_edgecolor("black")
    cm.ax_row_colors.collections[0].set_linewidth(0.1)
    cm.ax_row_colors.collections[0].set_edgecolor("black")

    # add the position of the legend. one for rows, one for cols
    if  legend is True:

        # get the maximum number of boxes in legend
        max_n_boxes_rows = 1 + len(row_colors_df.columns) + len(set.union(*[set(idx.split(idxs_separator_pattern)) for idx in row_colors_df.index]))
        max_n_boxes_cols = 1 + len(col_colors_df.columns) + len(set.union(*[set(idx.split(idxs_separator_pattern)) for idx in col_colors_df.index]))
        max_nboxes = max([max_n_boxes_cols, max_n_boxes_rows])

        # go through rows and cols
        for Ilegend, (legendlabel, colors_df) in enumerate([("rows", row_colors_df), ("cols", col_colors_df)]):

            # add the label for the type of column
            cm.ax_heatmap.bar(1, 0, color="white", label="--%s--"%legendlabel, linewidth=0); I=1

            # go through each of the labels of the cols
            for If, field in enumerate(colors_df.columns):

                # get each of the label and color combinations
                label_and_color_combinations = set(zip(map(lambda x: x.split(idxs_separator_pattern)[If], colors_df.index), colors_df[field]))

                # add bar for the field
                cm.ax_heatmap.bar(1, 0, color="white", label=field, linewidth=0); I+=1

                for label, color in sorted(label_and_color_combinations):

                    # add bar
                    cm.ax_heatmap.bar(1, 0, color=color, label=label, linewidth=0); I+=1

            # add missing bars to the legend
            for extraBarI in range(max_nboxes-I): cm.ax_heatmap.bar(1, 0, color="white", label="-", linewidth=0)

        # add the legend, depending on the location
        #cm.ax_heatmap.legend(bbox_to_anchor=(hm_pos.x0+hm_pos.width+add_to_legend_x, hm_pos.y0, hm_pos.width, hm_pos.height), ncol=2) # loc is the lower left corner bbox_to_anchor
        cm.ax_heatmap.legend(bbox_to_anchor=(2, 2), ncol=2) # loc is the lower left corner bbox_to_anchor

        ###################################



    # add graphics to y labels if provided
    if ylabels_graphics_df is not None:

        for I, (fw, c) in enumerate(ylabels_graphics_df[["fontweight", "color"]].values):

            cm.ax_heatmap.get_yticklabels()[I].set_weight(fw) 
            cm.ax_heatmap.get_yticklabels()[I].set_color(c) 


    # SAVE
    print_if_verbose("saving %s"%filename)
    cm.savefig(filename)

    return cm


def getPlots_filtering_accuracy_across_genomes_and_ploidies(df_cross_benchmark, PlotsDir, simName_to_color={"simulation_1":"black", "simulation_2":"red", "simulation_3":"green"}, simType_to_color={'biased_towards_repeats':"red", 'realData':"black", 'uniform':"blue", "simulated":"blue"}, ploidy_to_color={'consensus_ref': 'gray', 'haploid': 'black', 'diploid_hetero': 'maroon', 'ref:3_var:1': 'red', 'ref:9_var:1': 'lightsalmon', 'ref:99_var:1': 'white'}, svtype_to_color={"tandemDuplications": "gray", "deletions": "black", "inversions": "blue", "translocations": "olive", "insertions": "red", "remaining":"magenta", "integrated":"c"}):

    """This function takes a df that has several training and testing genomes with accuracy in each, and makes several plots to represent the data under outdir """

    #### cross-accuracy heatmap between combinations of parameters ###

    # define training and testing sets
    all_train_ploidies = set(df_cross_benchmark.train_ploidy)
    all_train_ploidies_str = "_".join(sorted(all_train_ploidies))

    all_train_genomeIDs = set(df_cross_benchmark.train_genomeID)
    all_train_genomeIDs_str = "_".join(sorted(all_train_genomeIDs))

    all_test_ploidies = set(df_cross_benchmark.test_ploidy)
    all_test_ploidies_str = "_".join(sorted(all_test_ploidies))

    all_test_genomeIDs = set(df_cross_benchmark.test_genomeID)
    all_test_genomeIDs_str = "_".join(sorted(all_test_genomeIDs))

    all_svtypes = {'integrated', 'translocations', 'tandemDuplications', 'remaining', 'inversions', 'deletions', 'insertions'}
    all_svtypes_str = "-".join(sorted(all_svtypes))

    # define the lists of ploidies

    # training
    #interesting_train_plodies_list = [all_train_ploidies, {"haploid"}, {"haploid", "diploid_hetero"}]
    #interesting_train_plodies_list = [{"haploid"}]
    interesting_train_plodies_list = [all_train_ploidies]

    #interesting_train_genomeIDs_list = [all_train_genomeIDs, {x for x in all_train_genomeIDs if "biased_towards_repeats" in x}, {x for x in all_train_genomeIDs if "uniform" in x}]
    interesting_train_genomeIDs_list = [all_train_genomeIDs]

    #interesting_train_svtypes_list = [all_svtypes] + [{x} for x in all_svtypes]
    interesting_train_svtypes_list = [all_svtypes]
    #interesting_train_svtypes_list = [{x} for x in all_svtypes]

    # testing

    #interesting_test_plodies_list = [all_test_ploidies, {"haploid"}, {"haploid", "diploid_hetero"}]
    interesting_test_plodies_list = [all_test_ploidies]

    #interesting_test_genomeIDs_list = [all_test_genomeIDs, {x for x in all_test_genomeIDs if "biased_towards_repeats" in x}, {x for x in all_test_genomeIDs if "uniform" in x}]
    interesting_test_genomeIDs_list = [all_test_genomeIDs]

    #interesting_test_svtypes_list = [all_svtypes] + [{x} for x in all_svtypes]
    interesting_test_svtypes_list = [all_svtypes]
    #interesting_test_svtypes_list = [{x} for x in all_svtypes]

    # accuracies and testings
    interesting_accuracies = ["Fvalue", "recall", "precision"]
    #interesting_accuracies = ["precision"]


    #row_cluster_list = [False, True]
    row_cluster_list = [True]
    #col_cluster_list = [False, True]
    col_cluster_list = [True]

    # map each accuracy 

    # go through each combination
    for interesting_train_plodies in interesting_train_plodies_list:
        for interesting_train_genomeIDs in interesting_train_genomeIDs_list:
            for interesting_train_svtypes in interesting_train_svtypes_list:
                for interesting_test_plodies in interesting_test_plodies_list:
                    for interesting_test_genomeIDs in interesting_test_genomeIDs_list:
                        for interesting_test_svtypes in interesting_test_svtypes_list:        

                            # get the filtered df
                            df = df_cross_benchmark[(df_cross_benchmark.train_ploidy.isin(interesting_train_plodies)) & (df_cross_benchmark.train_genomeID.isin(interesting_train_genomeIDs)) & (df_cross_benchmark.train_svtype.isin(interesting_train_svtypes)) & (df_cross_benchmark.test_ploidy.isin(interesting_test_plodies)) & (df_cross_benchmark.test_genomeID.isin(interesting_test_genomeIDs)) & (df_cross_benchmark.test_svtype.isin(interesting_test_svtypes)) & (df_cross_benchmark.nevents>=10)]

                            # add a simulation tag

                            # add the train and test indices
                            df["train_idx"] = (df.train_simName + "||||" + df.train_simType + "||||" + df.train_ploidy + "||||" + df.train_svtype)
                            df["test_idx"] = (df.test_simName + "||||" + df.test_simType + "||||" + df.test_ploidy + "||||" + df.test_svtype)

                            # find the train test that has the best Fvalue (the one that is maximum in the)
                            df_square_best = df[["train_idx", "test_idx", "Fvalue"]].pivot(index='train_idx', columns='test_idx', values="Fvalue")
                            train_idx_to_accuracy_df = df_square_best.apply(lambda r: pd.Series({"min_Fscore":min(r), "mean_Fscore":np.mean(r), "max_Fscore":max(r), "inverse_std_Fscore": 1/np.std(r), "median_Fscore":np.median(r)}) , axis=1).sort_values(by=["min_Fscore", "inverse_std_Fscore", "mean_Fscore", "median_Fscore", "max_Fscore"])
                            best_train_idx = train_idx_to_accuracy_df.iloc[-1].name

                            # go through each accuracy measurement
                            for accuracy in interesting_accuracies:

                                # prepare the dataframe so that the rows reflect the trainset and the cols the test set
                                df_square = df[["train_idx", "test_idx", accuracy]].pivot(index='train_idx', columns='test_idx', values=accuracy)
                                df_square_nevents = df[["train_idx", "test_idx", "nevents"]].pivot(index='train_idx', columns='test_idx', values="nevents")

                                # generate the cols colors df
                                def get_colors_series(idx):

                                    # get the different keywords
                                    simName, simType, ploidy, svtype = idx.split("||||")

                                    return pd.Series({"simName":simName_to_color[simName], "simType":simType_to_color[simType], "ploidy":ploidy_to_color[ploidy], "svtype":svtype_to_color[svtype]})
                                
                                row_colors_df = pd.Series(df_square.index, index=df_square.index).apply(get_colors_series)
                                col_colors_df = pd.Series(df_square.columns, index=df_square.columns).apply(get_colors_series)

                                #### get df annotations if the svtype is the same #######

                                # and a + when there are very few cases in the test
                                colIDX_to_rowIDX_to_label = {}
                                for col in df_square.columns: # iterate through testing
                                    svtype_col = col.split("||||")[-1]

                                    for row in df_square.index: # iterate through training
                                        svtype_row = row.split("||||")[-1]

                                        # define the number of events in test
                                        nevents_test = df_square_nevents.loc[row, col]

                                        # if they are the same
                                        if nevents_test<10: label = str(nevents_test)
                                        elif col==row: label = "=" 
                                        elif svtype_col==svtype_row: label = "*"
                                        else: label = ""

                                        # when the row is the same as the best_train_idx, add a best to the label
                                        if row==best_train_idx: 
                                            #label = "%s(B)"%label (a best)

                                            # add the accuracy
                                            if df_square.loc[row, col]==1: label = "1"
                                            else: label = "." + ("%.2f"%df_square.loc[row, col]).split(".")[1]

                                        # keep
                                        colIDX_to_rowIDX_to_label.setdefault(col, {}).setdefault(row, label)

                                df_annotations = pd.DataFrame(colIDX_to_rowIDX_to_label)

                                #########################################################

                                # go through each clustering
                                for row_cluster in row_cluster_list:
                                    for col_cluster in col_cluster_list:

                                        # define the filename

                                        #train vals
                                        train_genomeIDs = "_".join(sorted(interesting_train_genomeIDs)); 
                                        if train_genomeIDs==all_train_genomeIDs_str: train_genomeIDs = "all"
                                        train_ploidies = "_".join(sorted(interesting_train_plodies)); 
                                        if train_ploidies==all_train_ploidies_str: train_ploidies = "all"
                                        train_svtypes = "_".join(sorted(interesting_train_svtypes)); 
                                        if train_svtypes==all_svtypes_str: train_svtypes = "all"

                                        # test vals
                                        test_genomeIDs = "_".join(sorted(interesting_test_genomeIDs)); 
                                        if test_genomeIDs==all_test_genomeIDs_str: test_genomeIDs = "all"
                                        test_ploidies = "_".join(sorted(interesting_test_plodies)); 
                                        if test_ploidies==all_test_ploidies_str: test_ploidies = "all"
                                        test_svtypes = "_".join(sorted(interesting_test_svtypes)); 
                                        if test_svtypes==all_svtypes_str: test_svtypes = "all"

                                        plots_dir_crossaccuracyHeatmaps = "%s/cross_accuracy_heatmaps"%PlotsDir; make_folder(plots_dir_crossaccuracyHeatmaps)
                                        string_title_train = "train_genomes:%s_ploidies:%s_svtypes:%s"%(train_genomeIDs, train_ploidies, train_svtypes)
                                        string_title_test = "test_genomes:%s_ploidies:%s_svtypes:%s"%(test_genomeIDs, test_ploidies, test_svtypes)

                            
                                        filename = "%s/accuracyHeatmap_%s.pdf"%(plots_dir_crossaccuracyHeatmaps, accuracy)
                                        

                                        filename = filename.replace(":", "_")

                                        #filename = "%s/accuracyHeatmap_shortName.pdf"%(plots_dir_crossaccuracyHeatmaps)


                                        #title = "%s\n%s"%(string_title_train, string_title_test)
                                        title ="checking overfitting in SV calling"


                                        # get the dendrogram, either adjusting or not
                                        plot_clustermap_with_annotation(df_square, row_colors_df, col_colors_df, filename, title=title, col_cluster=col_cluster, row_cluster=row_cluster, colorbar_label=accuracy, adjust_position=False, legend=True, idxs_separator_pattern="||||", texts_to_strip={"L001"}, default_label_legend="control", df_annotations=df_annotations, cmap=sns.color_palette("RdBu_r", 50), ylabels_graphics_df=None, grid_lines=True)
    

    # return the best_train_idx
    return dict(zip(["simName", "simType", "ploidy", "svtype"], best_train_idx.split("||||")))

    ##################################################################

######################################################
######################################################
######################################################


def get_and_report_filtering_accuracy_across_genomes_and_ploidies(df_benchmark, genomeID_to_knownSVdict, outdir, PlotsDir, reference_genome, replace=False, consider_integrated_filtering=True, threads=4, run_in_parallel=False):

    """This function takes a df that has the benchmarking info (each line is a set of filtering parameters) and a "genomeID" and "ploidy" fields, which indicate the unique genomes. The idea is to pick, for each genomeID and ploidy combination, the best filtering set (highest Fscore) for the others simulations. If the genomeID is not in df_benchmark.genomeID it will run the whole gridss pipeline for the desired filtering sets, considering known_genomes_withSV_and_shortReads_table """

    # general vars
    
    ##### add an overal accuracy measurement for df_benchmark  ########
    if consider_integrated_filtering is True:

        # define the expected df field
        df_benchmark_with_integratedInfo_file = "%s/df_benchmark_all_with_integratedInfo.py"%outdir


        if file_is_empty(df_benchmark_with_integratedInfo_file) or replace is True:

            # define the fileds related to a filtering ID
            fields_filtering_ID = ['bedpe', 'benchmarkID','genomeID', 'gridss_VCFoutput', 'gridss_maxcoverage', 'gridss_regionsToIgnoreBed', 'ploidy']

            # get as filtered df
            print_if_verbose("getting all the variants integrated for each set of filters")
            df_benchmark_allSVtypes = df_benchmark.groupby(fields_filtering_ID, as_index=True).apply(get_integrated_benchmarking_fields_series_for_setFilters_df)

            # add the indices as fields, to match those of df_benchmark
            for I, field in enumerate(fields_filtering_ID): df_benchmark_allSVtypes[field] =  df_benchmark_allSVtypes.index.get_level_values(I)
            df_benchmark_allSVtypes = df_benchmark_allSVtypes.set_index("svtype", drop=False)

            # append to df benchmark
            df_benchmark = df_benchmark.append(df_benchmark_allSVtypes[list(df_benchmark.keys())])

            # keep
            save_object(df_benchmark, df_benchmark_with_integratedInfo_file)

        else: df_benchmark = load_object(df_benchmark_with_integratedInfo_file)

    #####################################################################

    # get, for each combination of genomeID, ploidy, svtype the best and less conservative filterset into a list
    print_if_verbose("Getting list of best filters for each genome, ploidy and svtype")
    IDs_sepparate_measurements = ["genomeID", "ploidy", "svtype"]
    df_best_filters = df_benchmark.groupby(IDs_sepparate_measurements).apply(get_best_less_conservative_row_df_benchmark)
    for I, field in enumerate(IDs_sepparate_measurements): df_best_filters[field] =  df_best_filters.index.get_level_values(I)

    # define the combinations of regions_to_ignore and max_coverage
    df_regionsIgnore_maxCov = df_best_filters[["gridss_regionsToIgnoreBed", "gridss_maxcoverage"]].drop_duplicates()


    ####### GENERATE A DF WITH THE INFO OF EACH SV SET TO BE TESTED THROUGH df_best_filters ############

    # initialize dicts to run cross-benchmarking
    genomeIDandPlody_to_info = {}

    # initialize a folder were files of the benchmarking will be stored
    cross_benchmarking_files_dir = "%s/cross_benchmarking_files"%outdir; make_folder(cross_benchmarking_files_dir)

    ## add the info for all genomes and ploidies found in df_best_filters
    for genomeID, ploidy in df_best_filters[["genomeID", "ploidy"]].drop_duplicates().values:

        # initialize test dict
        test_gridss_info_dict = {}

        # go through the gridss filterings
        for gridss_regionsToIgnoreBed, gridss_maxcoverage in df_regionsIgnore_maxCov.values:

            # find in df_benchmark the corresponding value
            df_ben_int = df_benchmark[(df_benchmark.genomeID==genomeID) & (df_benchmark.ploidy==ploidy) & (df_benchmark.gridss_maxcoverage==gridss_maxcoverage) & (df_benchmark.gridss_regionsToIgnoreBed==gridss_regionsToIgnoreBed)][["gridss_VCFoutput", "sorted_bam", "median_coverage", "median_insert_size", "median_insert_size_sd"]].drop_duplicates()

            # debug
            if len(df_ben_int)!=1: raise ValueError("There are not only one gridss vcfs with the given genomeID and ploidy")
            gridss_VCFoutput = df_ben_int.gridss_VCFoutput.iloc[0]

            # get the known vars
            knownSVdict = genomeID_to_knownSVdict[genomeID]

            # keep into test_gridss_info_dict
            test_gridss_info_dict.setdefault(gridss_regionsToIgnoreBed, {}).setdefault(gridss_maxcoverage, gridss_VCFoutput)

        # keep in the dict
        knownSVdict = genomeID_to_knownSVdict[genomeID]
        processing_dir = "%s/%s_%s"%(cross_benchmarking_files_dir, genomeID, ploidy); make_folder(processing_dir)
        genomeIDandPlody_to_info[(genomeID, ploidy)] = {"test_SVdict":knownSVdict, "outdir":processing_dir, "df_filters_train":df_best_filters, "test_gridss_info_dict":test_gridss_info_dict, "sorted_bam":df_ben_int.sorted_bam.iloc[0], "median_coverage":df_ben_int.median_coverage.iloc[0], "median_insert_size":df_ben_int.median_insert_size.iloc[0], "median_insert_size_sd":df_ben_int.median_insert_size_sd.iloc[0]}

    ###########################################################################################

    # run a function in parallel that will take a genome and ploidy combination and evaluate the accuracy of all the filters in df_best_filters. first prepare input as list of tuples
    list_inputs = [(d["test_SVdict"], d["outdir"], d["df_filters_train"], d["test_gridss_info_dict"], genomeID, ploidy, d["sorted_bam"], reference_genome, d["median_coverage"], d["median_insert_size"], d["median_insert_size_sd"],  replace) for (genomeID, ploidy), d in genomeIDandPlody_to_info.items()] # 0 is for debug

    # get the cross benchmarking df
    df_cross_benchmark_file = "%s/df_cross_benchmark.py"%outdir

    if file_is_empty(df_cross_benchmark_file) or replace is True:

        if run_in_parallel is True:

            # run in parallel
            with multiproc.Pool(threads) as pool:
                list_cross_benchmarking_dfs = pool.starmap(get_benchmarking_df_for_testSVs_from_trainSV_filterSets, list_inputs) # needs if __name__=="__main__" 
                
                pool.close()
                pool.terminate()

        else:

            list_cross_benchmarking_dfs = list(map(lambda x: get_benchmarking_df_for_testSVs_from_trainSV_filterSets(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]), list_inputs))
        
        # concatenate all the dfs
        df_cross_benchmark = pd.concat(list_cross_benchmarking_dfs)

        # save
        save_object(df_cross_benchmark, df_cross_benchmark_file)

    else: df_cross_benchmark = load_object(df_cross_benchmark_file)


    # add the simulation type for train and test
    def add_simulation_name_and_type(genomeID, tag):

        if "biased_towards_repeats" in genomeID: 
            simName = genomeID.split("_simType_biased_towards_repeats")[0]
            simType = "biased_towards_repeats"

        elif "uniform" in genomeID: 
            simName = genomeID.split("_simType_uniform")[0]
            simType = "uniform"

        elif "realData" in genomeID: 
            simName = genomeID
            simType = "realData"

        else:
            simName = genomeID
            simType = "simulated"

        return pd.Series({"%s_simName"%tag : simName, "%s_simType"%tag : simType})

    df_cross_benchmark[["train_simName", "train_simType"]] = df_cross_benchmark.train_genomeID.apply(lambda x: add_simulation_name_and_type(x, "train"))
    df_cross_benchmark[["test_simName", "test_simType"]] = df_cross_benchmark.test_genomeID.apply(lambda x: add_simulation_name_and_type(x, "test"))

    ######### get the plots and the best filters ###########

    best_filters_dict = getPlots_filtering_accuracy_across_genomes_and_ploidies(df_cross_benchmark, PlotsDir)

    #########################################################

    # get the best filters
    best_filters_series = df_best_filters.loc[(best_filters_dict["simName"], best_filters_dict["ploidy"], best_filters_dict["svtype"])]

    # get the df_cross_benchmark that matches the best series
    df_cross_benchmark_best = df_cross_benchmark[(df_cross_benchmark.train_genomeID==best_filters_dict["simName"]) & (df_cross_benchmark.train_ploidy==best_filters_dict["ploidy"]) & (df_cross_benchmark.train_svtype==best_filters_dict["svtype"]) & (df_cross_benchmark.nevents>5)]

    print_if_verbose("These are the accuracy measurements on all data")
    print_if_verbose(df_cross_benchmark_best[["test_genomeID", "test_ploidy", "Fvalue", "precision", "recall"]].sort_values(by=["Fvalue", "precision"]))

    return df_cross_benchmark_best, best_filters_series

def get_best_parameters_for_GridssClove_run(sorted_bam, reference_genome, outdir, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", nvars=100, real_svtype_to_file={}, median_insert_size=250, median_insert_size_sd=0):

    """This finds the optimum parameters for running GRIDSS clove and returns them. The parameters are equivalent to the run_GridssClove_optimising_parameters function"""


    # define plots dir
    PlotsDir = "%s/plots"%outdir; make_folder(PlotsDir)

    ###### MODELLING COVERAGE ######
    print_if_verbose("modelling coverage of the sample")

    # get a graph of the genome without any breakpoints (thus the None, None)
    genomeGraph_outfileprefix = "%s/genomeGraph_withoutBPs"%(outdir)
    genome_graph, df_positions_graph = get_genomeGraph_object(reference_genome, None, None, genomeGraph_outfileprefix, replace=replace)

    # get a function that takes the GC content, chromosome and distance to the telomere and returns coverage. This is actually a lambda function
    outdir_coverage_calculation = "%s/coverage_per_regions%ibb"%(outdir, window_l); make_folder(outdir_coverage_calculation)
    df_coverage_train = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_coverage_calculation, sorted_bam, windows_file="none", replace=replace), sep="\t")

    distToTel_chrom_GC_to_coverage_fn = get_distanceToTelomere_chromosome_GCcontent_to_coverage_fn(df_coverage_train, reference_genome, genome_graph, df_positions_graph, outdir_coverage_calculation, mitochondrial_chromosome=mitochondrial_chromosome, replace=replace)

    print_if_verbose("coverage model obtained")

    ################################

    ############ GENERAL OPERATIONS THAT WILL BE NEEDED FOR ALL THE STEPS #####

    # the dir and genome names
    genome_dir = "/".join(reference_genome.split("/")[0:-1])
    genome_name = reference_genome.split("/")[-1].split(".")[0]

    # map each chromosome to length
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

    # count the length od the reads
    read_length = get_read_length(sorted_bam, threads=threads, replace=replace)
    print_if_verbose("The median read length is %i"%read_length)

    # count total number of reads
    total_nread_pairs = count_number_read_pairs(sorted_bam, replace=False, threads=threads)
    #total_nread_pairs  = 100000 # this is to debug the simulation pipeline
    expected_coverage_per_bp = int((total_nread_pairs*read_length) / sum(chr_to_len.values())) +  1 # the expected coverage per position with pseudocount
    print_if_verbose("There are %i read pairs in your library. The expected coverage is %ix."%(total_nread_pairs, expected_coverage_per_bp))

    # get the info of the reference genome with predictions of coverage per window
    df_REFgenome_info = get_windows_infoDF_with_predictedFromFeatures_coverage(reference_genome, distToTel_chrom_GC_to_coverage_fn, expected_coverage_per_bp, replace=replace, threads=threads)

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
            print_if_verbose("working on simulation %i"%simulation_ID)

            # get an outdir where all the simulations of this ID will be stored
            simulation_outdir = "%s/simulation_%i"%(outdir, simulation_ID); make_folder(simulation_outdir)

            # get the simulated SVs, which are an integration of 
            sim_svtype_to_svfile, rearranged_genome = rearrange_genomes_simulateSV(reference_genome, simulation_outdir, replace=replace, nvars=nvars, mitochondrial_chromosome=mitochondrial_chromosome, simulated_svtype_to_svfile=real_svtype_to_file, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"})

            # define the genome ID
            genomeID = "simulation_%i"%(simulation_ID)

            # get a df that has genome info and predicted coverage from seq fetaures
            df_genome_info = get_windows_infoDF_with_predictedFromFeatures_coverage(rearranged_genome, distToTel_chrom_GC_to_coverage_fn, expected_coverage_per_bp, replace=replace, threads=threads)

            # get the aligned reads to the reference
            simulation_bam_file = simulate_and_align_PairedReads_perWindow(df_genome_info, rearranged_genome, reference_genome, total_nread_pairs, read_length, simulation_outdir, median_insert_size, median_insert_size_sd, replace=replace, threads=threads)

            # define a path to the known SVs (know_SV_dict should be changed to sim_svtype_to_svfile)

            # add the "remaining" cathegory, as an empty field
            remaining_file = "%s/remaining_sv.tab"%(get_dir(sim_svtype_to_svfile["translocations"]))
            open(remaining_file, "w").write("ID\t#CHROM\tPOS\tCHR2\tSTART\tEND\tSVTYPE\n"+"iii\tccc\t0\tyyy\t0\t0\tzzz\n")
            sim_svtype_to_svfile["remaining"] = remaining_file

            # map each genomeID to the known variants
            genomeID_to_knownSVdict[genomeID] = sim_svtype_to_svfile

            # go through each of the target ploidies and generate the resulting bam files:
            for ploidy in simulation_ploidies:
                print_if_verbose("working on %s"%ploidy)

                # define the final sorted bam depending on the ploidy (which also includes populations)
                ploidy_merged_bam = get_merged_bamfile_for_ploidy(variant_bamfile=simulation_bam_file, reference_bamfile=simulated_reference_bam_file, ploidy=ploidy, replace=replace, threads=threads)

                # calculate the expected fraction of reads comming from each genome
                fraction_var, fraction_ref = get_fractions_reads_for_ploidy(ploidy)

                # write a table and some files with the benchmarking of several filtering strategies of the data
                ploidy_dir = "%s/benchmark_GridssClove_%s"%(simulation_outdir, ploidy); make_folder(ploidy_dir)

                # get a df with a benchmark of many different parameters. This will also report some plots with the 
                benchmarking_df = benchmark_GridssClove_for_knownSV(ploidy_merged_bam, reference_genome, sim_svtype_to_svfile, ploidy_dir, range_filtering=range_filtering_benchmark, expected_AF=fraction_var, replace=replace, threads=threads, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, mitochondrial_chromosome=mitochondrial_chromosome)


                # add some parms and keep
                benchmarking_df["genomeID"] = [genomeID]*len(benchmarking_df)
                benchmarking_df["ploidy"] = [ploidy]*len(benchmarking_df)
                benchmarking_df["sorted_bam"] = [ploidy_merged_bam]*len(benchmarking_df)
                benchmarking_df["median_insert_size"] = [median_insert_size]*len(benchmarking_df)
                benchmarking_df["median_insert_size_sd"] = [median_insert_size_sd]*len(benchmarking_df)
                df_benchmark_all = df_benchmark_all.append(benchmarking_df, sort=True)

        print_if_verbose("GRIDSS simulation finished correctly")

        # save important files
        print_if_verbose("saving important files...")
        save_object(df_benchmark_all, df_benchmark_all_file)
        save_object(genomeID_to_knownSVdict, genomeID_to_knownSVdict_file)


    else:
        print_if_verbose("GRIDSS simulation finished correctly. Loading previous files ...")
        df_benchmark_all = load_object(df_benchmark_all_file)
        genomeID_to_knownSVdict = load_object(genomeID_to_knownSVdict_file)

    ####################################################

    ################### REPORT ACCURACY BETWEEN PARAMETERS OF DIFFERENT OPTIMISATIONS ####################

    # we will take the parameters that work best for all the simulations that we input

    # define the outputs
    outdir_benchmarking = "%s/benchmarking_all_filters_for_all_genomes_and_ploidies"%outdir; make_folder(outdir_benchmarking)
    PlotsDir_benchmarking = "%s/plots"%outdir_benchmarking; make_folder(PlotsDir_benchmarking)


    print_if_verbose("getting report of the accuracies between simulations")

    df_cross_benchmark_best, best_f = get_and_report_filtering_accuracy_across_genomes_and_ploidies(df_benchmark_all, genomeID_to_knownSVdict, outdir_benchmarking, PlotsDir_benchmarking, reference_genome, replace=replace, consider_integrated_filtering=True, threads=threads)

    # write this df
    df_cross_benchmark_best.to_csv("%s/df_cross_benchmark_best.tab"%(outdir_benchmarking), sep="\t", header=True, index=False)

    # define the filters
    gridss_blacklisted_regions = best_f["gridss_regionsToIgnoreBed"]
    gridss_maxcoverage = best_f["gridss_maxcoverage"]
    gridss_filters_dict = best_f["filters_dict"]
    max_rel_coverage_to_consider_del = best_f["clove_max_rel_coverage_to_consider_del"]
    min_rel_coverage_to_consider_dup = best_f["clove_min_rel_coverage_to_consider_dup"]

    ######################################################################################################


    #### remove unnecessary files ####

    # simulations
    #for simulation_ID in range(1, n_simulated_genomes+1): delete_folder("%s/simulation_%i"%(outdir, simulation_ID))

    # remove the cross-benchmarking files
    delete_folder("%s/cross_benchmarking_files"%outdir_benchmarking)


    ####################################

    return gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup, df_cross_benchmark_best



class NpEncoder(json.JSONEncoder):

    """A class to encode the json for numpy objects"""

    def default(self, obj):

        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)

def write_gridss_parameters_as_json(gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup, json_file, replace=False):

    """
    This function takes all the gridss and clove parameters and writes a json under json_file
    """

    if file_is_empty(json_file) or replace is True:
        print_if_verbose("saving parameters into %s"%json_file)

        # get the json
        parameters_dict = {"gridss_blacklisted_regions":gridss_blacklisted_regions,
                           "gridss_maxcoverage":gridss_maxcoverage,
                           "gridss_filters_dict":gridss_filters_dict,
                           "max_rel_coverage_to_consider_del":max_rel_coverage_to_consider_del,
                           "min_rel_coverage_to_consider_dup":min_rel_coverage_to_consider_dup}

        parameters_json = json.dumps(parameters_dict, cls=NpEncoder, indent=4, sort_keys=True)

        # write
        json_file_tmp = "%s.tmp"%json_file
        open(json_file_tmp, "w").write(parameters_json)
        os.rename(json_file_tmp, json_file)

def get_parameters_from_json(json_file):

    """Takes a json file with the parameters and returns each important thing"""

    with open(json_file) as f: parameters_dict = json.load(f)

    # define things
    gridss_blacklisted_regions = parameters_dict["gridss_blacklisted_regions"]
    gridss_maxcoverage = parameters_dict["gridss_maxcoverage"]
    gridss_filters_dict_initial = parameters_dict["gridss_filters_dict"]
    max_rel_coverage_to_consider_del = parameters_dict["max_rel_coverage_to_consider_del"]
    min_rel_coverage_to_consider_dup = parameters_dict["min_rel_coverage_to_consider_dup"]

    # modify the filters
    def convert_tuple_to_list(x):
        if type(x)==list: return tuple(x)
        else: return x

    gridss_filters_dict = {k : convert_tuple_to_list(v) for k,v in gridss_filters_dict_initial.items()}

    return gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup

def run_GridssClove_optimising_parameters(sorted_bam, reference_genome, outdir, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", nvars=100, fast_SVcalling=False, real_svtype_to_file={}, gridss_VCFoutput=""):

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
    - real_svtype_to_file is a dict that maps each type of sv to the file that contains it. It is useful to introduce a previously found set of SVs.
    - gridss_VCFoutput is passed to run_gridssClove_given_filters. It can be useful if you don't want to rerun gridss, which is very expensive.

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
    print_if_verbose("The median insert size is %i, with an absolute deviation of %i"%(median_insert_size, median_insert_size_sd))

    #####################################

    ###### GET PARAMETERS ######

    # get the default running and filtering parameters
    if fast_SVcalling is False:
        print_if_verbose("getting optimised parameters")

        parameter_optimisation_dir = "%s/parameter_optimisation"%outdir; make_folder(parameter_optimisation_dir)

        gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup, df_cross_benchmark_best = get_best_parameters_for_GridssClove_run(sorted_bam, reference_genome, parameter_optimisation_dir, threads=threads, replace=replace, n_simulated_genomes=n_simulated_genomes, mitochondrial_chromosome=mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=range_filtering_benchmark, nvars=nvars, real_svtype_to_file=real_svtype_to_file, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd)

    # get the parameters from an optimisation
    else: 
        print_if_verbose("running with default un-optimised parameters")

        # get the default parameters 
        gridss_blacklisted_regions = default_gridss_blacklisted_regions
        gridss_maxcoverage = default_gridss_maxcoverage
        gridss_filters_dict = default_filtersDict_gridss
        max_rel_coverage_to_consider_del = default_max_rel_coverage_to_consider_del
        min_rel_coverage_to_consider_dup = default_min_rel_coverage_to_consider_dup

    ###########################

    ###### FINAL GRIDSS-CLOVE RUNNING ######

    # define the final outdir 
    outdir_gridss_final = "%s/final_gridss_running"%outdir; make_folder(outdir_gridss_final)

    # write the parameters of the running
    json_file = "%s/perSVade_parameters.json"%outdir_gridss_final
    write_gridss_parameters_as_json(gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup, json_file, replace=replace)

    # define the final vcf dir
    final_gridss_vcf = "%s/output_gridss.vcf"%outdir_gridss_final

    # if there is a provided gridss_VCFoutput, softlink it to the outdir_gridss_final
    if not file_is_empty(gridss_VCFoutput) and file_is_empty(final_gridss_vcf): soft_link_files(gridss_VCFoutput, final_gridss_vcf)

    # define the median coverage across window_l windows of the genome
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_gridss_final, sorted_bam, windows_file="none", replace=replace), sep="\t")
    median_coverage = np.median(coverage_df[~coverage_df["#chrom"].isin(mitochondrial_chromosome.split(","))].mediancov_1)
    print_if_verbose("The median coverage is %i"%median_coverage)

    # run the pipeline
    print_if_verbose("running final gridss with parameters...")
    final_sv_dict, df_gridss = run_gridssClove_given_filters(sorted_bam, reference_genome, outdir_gridss_final, median_coverage, replace=replace, threads=threads, gridss_blacklisted_regions=gridss_blacklisted_regions, gridss_VCFoutput=final_gridss_vcf, gridss_maxcoverage=gridss_maxcoverage, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, gridss_filters_dict=gridss_filters_dict, run_in_parallel=True, max_rel_coverage_to_consider_del=max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup=min_rel_coverage_to_consider_dup, replace_FromGridssRun=replace)

    ########################################

    ##### PIPELINE ENDING OPERATIONS ##### 

    # at the end, remove all the mosdepth and windows files under the reference
    clean_reference_genome_windows_files(reference_genome)
    
    print_if_verbose("GRIDSS pipeline finished correctly")

    print_if_verbose("--- the gridss pipeline optimising parameters took %s seconds in %i cores ---"%(time.time() - pipeline_start_time, threads))

    # generate a file that indicates whether the gridss run is finished
    final_file = "%s/gridssClove_finished.txt"%outdir
    open(final_file, "w").write("gridssClove_finished finished...")

    ######################################

    return final_sv_dict, df_gridss


def plot_report_accuracy_simulations(df_benchmarking, filename):

    """Plots the accuracy of each type of SVcalling on simulations. There will be one subplot for each precision/recall/Fvalue and ploidy combination. The rows will be for ploidies and the cols for accuracy measurements."""

    df_benchmarking_long = pd.DataFrame()

    accuracy_fields = ["precision", "recall", "Fvalue"]
    nonAccuracy_fields = [c for c in df_benchmarking.columns if c not in accuracy_fields]

    # convert to a long format
    for f in accuracy_fields:

        df = df_benchmarking[nonAccuracy_fields + [f]].rename(columns={f : "accuracy"})
        df["type_ac"] = f
        df_benchmarking_long = df_benchmarking_long.append(df)

    # define the palette
    palette = {"uniform":"navy", "realSVs":"red", "fastSV_on_uniform":"cyan", "fastSV_on_realSVs":"magenta"}


    # change the svtype name
    svtype_to_shortSVtype = {"deletions":"del", "tandemDuplications":"tan", "insertions":"ins", "translocations":"tra", "inversions":"inv", "integrated":"all"}
    df_benchmarking_long["svtype"] = df_benchmarking_long.svtype.apply(lambda x: svtype_to_shortSVtype[x])

    g = sns.catplot(x="svtype", y="accuracy", hue="typeParameterOptimisation", col="type_ac", row="ploidy", data=df_benchmarking_long,kind="bar", height=2, aspect=2, ci="sd", palette=palette)

    g.set(ylim=(0.5, 1))

    g.savefig(filename, bbox_inches='tight')

def report_accuracy_simulations(sorted_bam, reference_genome, outdir, real_svtype_to_file, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", nvars=100):

    """This function runs tests the accuracy on uniform simulations and those indicated by real_svtype_to_file. It reports the accuracy of the parameter optimisation on each of these simulations, as well as how the fast method works for each of the different simulation types. """

    # test that there is real data
    if len(real_svtype_to_file)==0: raise ValueError("You need real data if you want to test accuracy")

    # make the outdir
    make_folder(outdir)

    print_if_verbose("Testing accuracy of simulations")

    # calculate the insert size statistics
    median_insert_size, median_insert_size_sd  = get_insert_size_distribution(sorted_bam, replace=replace, threads=threads)

    # define a file that will contain the benchmarking
    df_benchmarking_file = "%s/df_benchmarking.tab"%outdir

    if file_is_empty(df_benchmarking_file) or replace is True:

        # initialize a df that will contain the accuracy of each simulation type. The fields will be genomeID, ploidy, svtype, typeParameterOptimisation (this can be uniform, realSVs, fastSV_on_uniform or fastSV_on_realSVs), Fvalue, precision and recall
        df_benchmarking_fields = ["genomeID", "ploidy", "svtype", "typeParameterOptimisation", "Fvalue", "precision", "recall"]
        df_benchmarking = pd.DataFrame(columns=df_benchmarking_fields)

        # go through each simulation type
        for typeSimulations, svtype_to_svfile in [("uniform", {}), ("realSVs", real_svtype_to_file)]:
            print_if_verbose(typeSimulations)

            # define the parameter optimisation dir
            parameter_optimisation_dir = "%s/parameter_optimisation_%s"%(outdir, typeSimulations); make_folder(parameter_optimisation_dir)

            # get the accuracy of these types of simulations, as well as the best parameters
            gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup, df_cross_benchmark_best = get_best_parameters_for_GridssClove_run(sorted_bam, reference_genome, parameter_optimisation_dir, threads=threads, replace=replace, n_simulated_genomes=n_simulated_genomes, mitochondrial_chromosome=mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=range_filtering_benchmark, nvars=nvars, real_svtype_to_file=svtype_to_svfile, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd)

            # add to the df the simulations paramteres
            df_cross_benchmark_best["typeParameterOptimisation"] = typeSimulations
            df_cross_benchmark_best = df_cross_benchmark_best.rename(columns={"test_genomeID": "genomeID", "test_ploidy":"ploidy"})
            df_benchmarking = df_benchmarking.append(df_cross_benchmark_best[df_benchmarking_fields])

            # go through each simulation and ploidy and run the fastSV calling on it. This will be fast by putting the vcf file under the 
            for simulation_ID in range(1, n_simulated_genomes+1):
                print_if_verbose(simulation_ID)

                # define the genome ID
                genomeID = "simulation_%i"%(simulation_ID)

                # define the known SVs
                knownSV_dict = {svtype : "%s/simulation_%i/final_simulated_SVs/%s.tab"%(parameter_optimisation_dir, simulation_ID, svtype) for svtype in {"insertions", "translocations", "deletions", "tandemDuplications", "inversions"}}

                # go through each of the target ploidies and generate the resulting bam files:
                for ploidy in simulation_ploidies:

                    # define the vcf outdir
                    benchmarking_dir = "%s/simulation_%i/benchmark_GridssClove_%s/benchmark_max%ix_ignoreRegions%s"%(parameter_optimisation_dir, simulation_ID, ploidy, gridss_maxcoverage, gridss_blacklisted_regions!="")
                    gridss_vcf = "%s/gridss_output.vcf.withSimpleEventType.vcf"%benchmarking_dir
                    print_if_verbose("working on %s"%ploidy)

                    # define the outdir
                    outdir_fastCalling = "%s/fastSVcalling_%s_simulation_%i_ploidy_%s"%(outdir, typeSimulations, simulation_ID, ploidy); make_folder(outdir_fastCalling)

                    # define the sorted bam
                    sorted_bam_fastSV = "%s/aligned_reads.sorted.bam"%benchmarking_dir

                    # debug
                    if any([file_is_empty(x) for x in [sorted_bam_fastSV, gridss_vcf]]): raise ValueError("Some files do not exist")

                    # get the fast calling (without repeating)
                    sv_dict, df_gridss = run_GridssClove_optimising_parameters(sorted_bam_fastSV, reference_genome, outdir_fastCalling, threads=threads, replace=replace, mitochondrial_chromosome=mitochondrial_chromosome, fast_SVcalling=True, gridss_VCFoutput=gridss_vcf)

                    # define the benchmarking accuracy
                    fileprefix = "%s/benchmarking"%outdir_fastCalling
                    df_benchmark_fastSVcalling = benchmark_processedSVs_against_knownSVs_inHouse(sv_dict, knownSV_dict, fileprefix, replace=replace, add_integrated_benchmarking=True)

                    # add fields
                    df_benchmark_fastSVcalling["genomeID"] = genomeID
                    df_benchmark_fastSVcalling["ploidy"] = ploidy
                    df_benchmark_fastSVcalling["typeParameterOptimisation"] = "fastSV_on_%s"%typeSimulations

                    # keep
                    df_benchmarking = df_benchmarking.append(df_benchmark_fastSVcalling[df_benchmarking_fields])


        # save
        df_benchmarking.to_csv(df_benchmarking_file, sep="\t", header=True, index=False)

    else: df_benchmarking = pd.read_csv(df_benchmarking_file, sep="\t")


    # make plots to report the accuracy of each simulation type. 
    filename = "%s/accuracy_on_simulations.pdf"%outdir
    plot_report_accuracy_simulations(df_benchmarking, filename)


def generate_jobarray_file_greasy(jobs_filename, walltime="48:00:00",  name="JobArray", queue="bsc_ls", sbatch=False, ncores_per_task=1, constraint="", email="miquel.schikora@bsc.es", number_tasks_to_run_at_once="all", max_ncores_queue=768):
    
    """ This function takes:
        jobs_filename: a path to a file in which each line is a command that has to be executed as an independent job
        name: the name of the jobs array
        walltime is the time in "dd-hh:mm:ss". It is the sum of job times
        ncores_per_task is the number of cores that each job gets
        
        name is the name prefix
        queue can be "debug" or "bsc_ls", use bsc_queues to understand which works best
        rmstd indicates if the previous std has to be removed
        constraint is a list of constraints to pass to sbatch. For example “highmem” is useful for requesting more memory. You cannot submit a job requesting memory parameters, memory is automatically set for each asked cpu (2G/core by default, 8G/core for highmem)

        number_tasks_to_run_at_once are the number of tasks in a job array to run at once
        
        It returns a jobs_filename.run file, which can be sbatch to the cluster directly if sbatch is True
        This is run in the VarCall_CNV_env
    """

    # define the stddir
    outdir = get_dir(jobs_filename)
    stddir = "%s/STDfiles"%outdir; make_folder(stddir)

    # remove all previous files from stddir that start with the same name
    for file in os.listdir(stddir): 
        if file.startswith(name): remove_file("%s/%s"%(stddir, file))
    
    # Get the number of jobs
    n_jobs = len(open(jobs_filename, "r").readlines())

    # if default, number_tasks_to_run_at_once is all, which means that it will try to run all tasks at once
    if number_tasks_to_run_at_once=="all": number_tasks_to_run_at_once = min([int(max_ncores_queue/ncores_per_task), n_jobs])

    # define the constraint only if it is necessary
    if constraint!="": constraint_line = "#SBATCH --constraint=%s"%constraint
    else: constraint_line = ""

    # remove previous rst files
    name_jobs_filename = get_file(jobs_filename)
    for file in os.listdir(get_dir(jobs_filename)): 

        if file.startswith("%s-"%name_jobs_filename) and file.endswith(".rst"): 
            remove_file("%s/%s"%(get_dir(jobs_filename), file))

    # rewrite the jobs_filename so that each std goes to a different file
    std_perJob_prefix = "%s/%s"%(stddir, name)
    jobs_filename_lines = ["%s > %s.%i.out 2>&1"%(l.strip(), std_perJob_prefix, I+1) for I, l in enumerate(open(jobs_filename, "r").readlines())]
    open(jobs_filename, "w").write("\n".join(jobs_filename_lines))

    # define the environment activating parms
    #"echo 'sourcing conda to run pipeline...';",
    #SOURCE_CONDA_CMD,
    #CONDA_ACTIVATING_CMD,

    # "#SBATCH --mail-type=all",
    # "#SBATCH --mail-user=%s"%email,

    # define the std files
    greasy_logfile = "%s/%s_greasy.log"%(stddir, name)
    stderr_file = "%s/%s_stderr.txt"%(stddir, name)
    stdout_file = "%s/%s_stdout.txt"%(stddir, name)

    # define arguments
    arguments = ["#!/bin/sh", # the interpreter
                 "#SBATCH --time=%s"%walltime, # several SBATCH misc commands
                 "#SBATCH --qos=%s"%queue,
                 "#SBATCH --job-name=%s"%name,
                 "#SBATCH --cpus-per-task=%i"%ncores_per_task,
                 "#SBATCH --error=%s"%(stderr_file), # the standard error
                 "#SBATCH --output=%s"%(stdout_file), # the standard output
                 "#SBATCH --get-user-env", # this is to maintain the environment
                 "#SBATCH --workdir=%s"%outdir,
                 "#SBATCH --ntasks=%i"%number_tasks_to_run_at_once,
                 constraint_line,
                 "",
                 "module load greasy",
                 "export GREASY_LOGFILE=%s;"%(greasy_logfile),
                 "echo 'running pipeline';",
                 "greasy %s"%jobs_filename]


    # define and write the run filename
    jobs_filename_run = "%s.run"%jobs_filename
    with open(jobs_filename_run, "w") as fd: fd.write("\n".join(arguments))
    
    # run in cluster if specified
    if sbatch is True: run_cmd("sbatch %s"%jobs_filename_run)

    # get info about the exit status: sacct -j <jobid> --format=JobID,JobName,MaxRSS,Elapsed



def generate_jobarray_file_slurm(jobs_filename, stderr="./STDERR", stdout="./STDOUT", walltime="02:00:00",  name="JobArray", queue="bsc_ls", sbatch=False, ncores_per_task=1, rmstd=True, constraint="", number_tasks_to_run_at_once="all", email="mikischikora@gmail.com"):
    
    """ This function takes:
        jobs_filename: a path to a file in which each line is a command that has to be executed in a sepparate cluster node
        name: the name of the jobs array
        stderr and stdout: paths to directories that will contains the STDERR and STDOUT files
        walltime is the time in "dd-hh:mm:ss"
        memory is the RAM: i.e.: 4G, 2M, 200K
        ncores_per_task is the number of cores that each job gets
        
        name is the name prefix
        queue can be "debug" or "bsc_ls", use bsc_queues to understand which works best
        rmstd indicates if the previous std has to be removed
        constraint is a list of constraints to pass to sbatch. For example “highmem” is useful for requesting more memory. You cannot submit a job requesting memory parameters, memory is automatically set for each asked cpu (2G/core by default, 8G/core for highmem)

        number_tasks_to_run_at_once are the number of tasks in a job array to run at once
        
        It returns a jobs_filename.run file, which can be sbatch to the cluster directly if sbatch is True
        This is run in the VarCall_CNV_env
    """

    def removeSTDfiles(stddir):
        """ Will remove all files in stddir with name"""
        for file in os.listdir(stddir):
            if file.startswith(name): os.unlink("%s/%s"%(stddir, file))

    # prepare the stderr and stdout
    if not os.path.isdir(stderr): os.mkdir(stderr)
    elif rmstd is True: removeSTDfiles(stderr)

    if not os.path.isdir(stdout): os.mkdir(stdout)
    elif rmstd is True: removeSTDfiles(stdout)
    
    # Get the number of jobs
    n_jobs = len(open(jobs_filename, "r").readlines())

    # if default, number_tasks_to_run_at_once is 0, which means that it will try to run all tasks at once
    if number_tasks_to_run_at_once=="all": number_tasks_to_run_at_once = n_jobs

    # define the number of nodes, consider that each node has 48 cpus, you need to request the number of nodes accordingly
    nnodes = int((ncores_per_task/48)+1) # get the non decimal part of a float

    # define the constraint only if it is necessary
    if constraint!="": constraint_line = "#SBATCH --constraint=%s"%constraint
    else: constraint_line = ""

    # define arguments
    arguments = ["#!/bin/sh", # the interpreter
                 "#SBATCH --time=%s"%walltime, # several SBATCH misc commands
                 "#SBATCH --qos=%s"%queue,
                 "#SBATCH --job-name=%s"%name,
                 "#SBATCH --cpus-per-task=%i"%ncores_per_task,
                 "#SBATCH --error=%s/%s_%sA_%sa.err"%(stderr, name, "%", "%"), # the standard error
                 "#SBATCH --output=%s/%s_%sA_%sa.out"%(stdout, name, "%", "%"), # the standard error
                 "#SBATCH --get-user-env", # this is to maintain the environment
                 "#SBATCH --workdir=.",
                 "#SBATCH --array=1-%i%s%i"%(n_jobs, "%", number_tasks_to_run_at_once),
                 #"#SBATCH --array=1-%i"%n_jobs,
                 "#SBATCH --mail-type=all",
                 "#SBATCH --mail-user=%s"%email,
                 constraint_line,
                 "",
                 "echo 'sourcing conda to run pipeline...';",
                 "echo 'running pipeline';",
                 'srun $(head -n $SLURM_ARRAY_TASK_ID %s | tail -n 1)'%jobs_filename]


    # define and write the run filename
    jobs_filename_run = "%s.run"%jobs_filename
    with open(jobs_filename_run, "w") as fd: fd.write("\n".join(arguments))
    
    # run in cluster if specified
    if sbatch is True: out_state = os.system("sbatch %s"%jobs_filename_run); print_if_verbose("%s sbatch out state: %i"%(name, out_state))

    # get info about the exit status: sacct -j <jobid> --format=JobID,JobName,MaxRSS,Elapsed


def plot_accuracy_simulations_from_all_sampleID_to_dfBestAccuracy(all_sampleID_to_dfBestAccuracy, filename):

    """This function takes a dict that maps each combination of typeSimulation||||runID to the df with the best accuracies. There will be one subplot for each precision/recall/Fvalue and ploidy combination. The rows will be for ploidies and the cols for accuracy measurements. """

    print_if_verbose("running plot_accuracy_simulations_from_all_sampleID_to_dfBestAccuracy")

    # generate df
    df_benchmarking = pd.DataFrame()

    for ID, df in all_sampleID_to_dfBestAccuracy.items():

        simID, sampleID = ID.split("||||")
        df["typeParameterOptimisation"] = simID
        df["sampleID"] = sampleID

        # keep
        df_benchmarking = df_benchmarking.append(df)


    # get a df benchamrking long
    df_benchmarking_long = pd.DataFrame()

    accuracy_fields = ["precision", "recall", "Fvalue"]
    nonAccuracy_fields = [c for c in df_benchmarking.columns if c not in accuracy_fields]

    # convert to a long format
    for f in accuracy_fields:

        df = df_benchmarking[nonAccuracy_fields + [f]].rename(columns={f : "accuracy"})
        df["ac"] = f
        df_benchmarking_long = df_benchmarking_long.append(df)

    # change the length of the df
    svtype_to_shortSVtype = {"deletions":"del", "tandemDuplications":"tan", "insertions":"ins", "translocations":"tra", "inversions":"inv", "integrated":"all"}
    df_benchmarking_long["svtype"] = df_benchmarking_long.svtype.apply(lambda x: svtype_to_shortSVtype[x])

    # rename the ploidy
    df_benchmarking_long = df_benchmarking_long.rename(columns={"test_ploidy":"ploidy"})

    # define the minimum accuracy
    min_y = max([min(df_benchmarking_long["accuracy"])-0.1, 0])

    # get a violin catplot
    print_if_verbose("getting %s"%filename)
    palette = {"uniform":"navy", "realSVs":"red"}

    g = sns.catplot(x="svtype", y="accuracy", hue="typeParameterOptimisation", col="ac", row="ploidy", data=df_benchmarking_long,kind="swarm", height=2, aspect=2, palette=palette, dodge=True)
    g.set(ylim=(min_y, 1))

    # add horizontal lines
    #for y in [0.5, 0.7, 0.8, 0.9]: g.map(plt.axhline, y=y, ls='--', c='k')

    g.savefig(filename, bbox_inches='tight')
    
def plot_fraction_overlapping_realSVs(df_benchmarking, filename):

    """This function takes a df that has the fraction_overlapping vars of each svtype for each runID. It plots it for each svtype. """

    print_if_verbose("getting %s"%filename)
    df_benchmarking = cp.deepcopy(df_benchmarking)

    palette = {"uniform":"navy", "realSVs":"red", "fast":"magenta"}

    # define the minimum accuracy
    svtype_to_shortSVtype = {"deletions":"del", "tandemDuplications":"tan", "insertions":"ins", "translocations":"tra", "inversions":"inv", "integrated":"all", "remaining":"rem"}
    df_benchmarking["svtype"] = df_benchmarking.svtype.apply(lambda x: svtype_to_shortSVtype[x])

    fig = plt.figure(figsize=(len(set(df_benchmarking.svtype)), 7))

    label_to_ylabel = {"fraction overlapping SVs": "fraction overlap. SVs ~ precision" , "n SVs":"n SVs ~ recall"}

    for I, y in enumerate(["fraction overlapping SVs", "n SVs"]): # 

        ax = plt.subplot(2, 1, I+1)

        df_benchmarking[y] = df_benchmarking[y].astype(float)

        # get a violin plot
        ax = sns.boxplot(x="svtype", y=y, data=df_benchmarking, hue="simulationID", palette=palette, boxprops=dict(alpha=.45))

        ax = sns.swarmplot(x="svtype", y=y, hue="simulationID", data=df_benchmarking, palette=palette, dodge=True, linewidth=.5, edgecolor="k")

        ax.legend(bbox_to_anchor=(1, 1))
        ax.set_xlabel("")
        ax.set_ylabel(label_to_ylabel[y])

        if I in [1]: ax.get_legend().remove()


    fig.savefig(filename, bbox_inches='tight')
    plt.close(fig)


def plot_fraction_overlapping_realSVs_precision_vs_recall(df_benchmarking, filename):

    """This function takes a df that has the fraction_overlapping vars of each svtype for each runID. It plots it for each svtype as a precision vs recall plot """

    print_if_verbose("getting %s"%filename)
    df_benchmarking = cp.deepcopy(df_benchmarking)

    palette = {"uniform":"navy", "realSVs":"red", "fast":"magenta"}
    simID_to_marker = {"uniform":"o", "realSVs":"s", "fast":"X"}

    # define the minimum accuracy
    nsvtypes = len(set(df_benchmarking.svtype))
    fig = plt.figure(figsize=(nsvtypes*5, 4))

    for I, svtype in enumerate(sorted(set(df_benchmarking.svtype))):
        ax = plt.subplot(1, nsvtypes, I+1)

        df_sv = df_benchmarking[df_benchmarking.svtype==svtype]

        # get a scatter plot
        ax = sns.scatterplot(x="n SVs", y="fraction overlapping SVs", data=df_sv, hue="simulationID", palette=palette, style="simulationID", markers=simID_to_marker)

        # add a shaded error bar with the median and mad
        for simID, color in palette.items():

            # define the data
            df_sim = df_sv[df_sv.simulationID==simID]
            x = df_sim["n SVs"].apply(float)
            y = df_sim["fraction overlapping SVs"].apply(float)

            # median
            median_x = np.median(x); mad_x = stats.median_absolute_deviation(x)
            median_y = np.median(y); mad_y = stats.median_absolute_deviation(y)

            # mean
            median_x = np.mean(x); mad_x = np.std(x)
            median_y = np.mean(y); mad_y = np.std(y)

            # add the patch
            fancybox = mpatches.FancyBboxPatch([median_x-mad_x, median_y-mad_y], 2*mad_x, 2*mad_y, boxstyle=mpatches.BoxStyle("Round", pad=0))
            ax.add_collection(PatchCollection([fancybox], color=color, alpha=0.3))

            # add the median
            plt.scatter(median_x, median_y, s=42, c=color, marker=simID_to_marker[simID], linewidths=.65, edgecolor="black")
        
        # graphic parms
        ax.set_title(svtype)
        if I==0: ax.set_ylabel("fraction overlap. SVs ~ precision")
        ax.set_xlabel("n SVs ~ recall")
        #ax.legend(bbox_to_anchor=(1, 1))
        #if I in [1]: ax.get_legend().remove()

    fig.savefig(filename, bbox_inches='tight')
    plt.close(fig)

def get_sorted_colors():
    
    """Returns a list with all CSS colors sorted"""

    colors = mcolors.CSS4_COLORS

    list_colors = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(color))), name) for name, color in colors.items())
    sorted_colors = [x[1] for x in list_colors]

    return sorted_colors

def get_sampleID_and_runID_to_color(df):

    """Returns a dict with unique sampleIDs and runIDs"""

    # find the sample_and_run_ID_combinations 
    sample_and_run_ID_combinations = set.union(*[{(s,r) for s,r in df[["%s_sampleID"%prefix, "%s_runID"%prefix]].values} for prefix in ["parms", "test"]])
    sampleID_to_runIDs = {x[0]:[] for x in sample_and_run_ID_combinations}
    for s, r in sample_and_run_ID_combinations: sampleID_to_runIDs[s].append(r)

    # get the sorted colors
    sorted_colors = get_sorted_colors()

    # define the location of the sample colors
    sample_colors_idxs = [int(x) for x in np.linspace(0, len(sorted_colors)-1, len(sampleID_to_runIDs))]
    
    # initialize the dicts
    sampleID_to_color = {}
    runID_to_color = {}
    for Isample, (sampleID, runIDs) in enumerate(sampleID_to_runIDs.items()):

        # define the color I
        colorI = sample_colors_idxs[Isample]

        # get the sample color
        sampleID_to_color[sampleID]  = sorted_colors[colorI]

        # start_runID_colorID
        min_start_runID_colorI = 0
        max_start_runID_colorI = max(sample_colors_idxs) - len(runIDs)
        start_runID_colorI = colorI - int(len(runIDs)/2)

        # debug the fact that it is inconsistent
        if start_runID_colorI<min_start_runID_colorI: start_runID_colorI = min_start_runID_colorI
        elif start_runID_colorI>max_start_runID_colorI: start_runID_colorI = max_start_runID_colorI

        # go through each run, so that you get colors arround the sample color
        for Irun, runID in enumerate(runIDs):

            # define the colorI of the run
            runID_colorI = start_runID_colorI +  Irun
            runID_to_color[runID] = sorted_colors[runID_colorI]

    return sampleID_to_color, runID_to_color
 
def generate_heatmap_accuracy_of_parameters_on_test_samples(df_benchmark, plots_dir, replace=False, threads=4):

    """
    This function takes a df where each row is one set of training parameters and test data svtype, together with the accuracy records. It generates a heatmap were the rows are each of the training parameters and the cols are the test samples.
    """

    print_if_verbose("plotting cross-accuracy")

    # define  graphics
    simName_to_color = {"simulation_1":"black", "simulation_2":"gray", "simulation_3":"green"}
    ploidy_to_color = {'consensus_ref': 'gray', 'haploid': 'black', 'diploid_hetero': 'maroon', 'ref:3_var:1': 'red', 'ref:9_var:1': 'lightsalmon', 'ref:99_var:1': 'white'}
    svtype_to_color = {"tandemDuplications": "gray", "deletions": "black", "inversions": "blue", "translocations": "olive", "insertions": "red", "remaining":"magenta", "integrated":"c"}
    typeSimulations_to_color = {"uniform":"blue", "realSVs":"red", "fast":"magenta"}
    sampleID_to_color, runID_to_color = get_sampleID_and_runID_to_color(df_benchmark) # automatic definition of graphics of sampleID and runID

    # map each cathegory to the colors
    cathegory_to_colors_dict = {"parms_sampleID" : sampleID_to_color,
                                "parms_runID" : runID_to_color,
                                "parms_typeSimulations": typeSimulations_to_color,
                                "test_sampleID" : sampleID_to_color,
                                "test_runID" : runID_to_color,
                                "test_typeSimulations": typeSimulations_to_color,
                                "test_simName" : simName_to_color,
                                "test_ploidy" : ploidy_to_color,
                                "svtype": svtype_to_color
                                }

    # define the lists of things
    #interesting_ploidies_list = [set(df_benchmark.test_ploidy), {"haploid"}, {"diploid_hetero"}]
    interesting_ploidies_list = [set(df_benchmark.test_ploidy)]

    #interesting_svtypes_list = [set(df_benchmark.svtype), {"integrated"}]
    interesting_svtypes_list = [{"integrated"}]

    
    #interesting_typeSimulations_list = [set(df_benchmark.parms_typeSimulations), {"fast"}, {"uniform", "realSVs"}]
    interesting_typeSimulations_list = [set(df_benchmark.parms_typeSimulations), {"uniform", "realSVs"}]
    
    interesting_accuracies = ["Fvalue", "precision", "recall"]

    # go through each accuracy measurement
    for accuracy in interesting_accuracies:
        for interesting_ploidies in interesting_ploidies_list:
            for interesting_svtypes in interesting_svtypes_list:
                for interesting_typeSimulations in interesting_typeSimulations_list:

                    # define the tags
                    if len(interesting_ploidies)==1: ploidy_tag = next(iter(interesting_ploidies))
                    else: ploidy_tag = "allPloidies"

                    if len(interesting_svtypes)==1: svtype_tag = next(iter(interesting_svtypes))
                    else: svtype_tag = "allSVtypes"

                    if len(interesting_typeSimulations)==1: typeSimulations_tag = next(iter(interesting_typeSimulations))
                    elif len(interesting_typeSimulations)==2: typeSimulations_tag = "realANDuniform"
                    else: typeSimulations_tag = "allSimulations"

                    # get the filtered df
                    df = df_benchmark[(df_benchmark.test_ploidy.isin(interesting_ploidies)) & (df_benchmark.svtype.isin(interesting_svtypes)) & (df_benchmark.parms_typeSimulations.isin(interesting_typeSimulations))]

                    # add the indices
                    parms_keys = [k for k in df.keys() if k.startswith("parms_")]
                    test_keys = [k for k in df.keys() if k.startswith("test_")] + ["svtype"]
                    df["parms_idx"] = df.apply(lambda r: "||||".join([r[k] for k in parms_keys]), axis=1)
                    df["test_idx"] = df.apply(lambda r: "||||".join([r[k] for k in test_keys]), axis=1)

                    # add the label
                    def get_label(r):

                        if r["parms_sampleID"]==r["test_sampleID"] and r["parms_runID"]==r["test_runID"] and r["parms_typeSimulations"]==r["test_typeSimulations"]: label = "="
                        else: label = ""

                        return label

                    df["label"] = df.apply(get_label, axis=1)

                    # get the square df
                    df_square = df[["parms_idx", "test_idx", accuracy]].pivot(index='parms_idx', columns='test_idx', values=accuracy)

                    # define dicts mapping objects
                    type_keys_to_keys = {"parms":parms_keys, "test":test_keys}


                    # generate the cols colors df
                    def get_colors_series(idx, type_keys="parms"):
                        # type_keys can be parms or test

                        # get the color dicts
                        keys = type_keys_to_keys[type_keys]

                        # get the content
                        idx_content = idx.split("||||")

                        # define the series
                        field_to_color = {keys[I] : cathegory_to_colors_dict[keys[I]][c] for I,c in enumerate(idx_content)}

                        return pd.Series(field_to_color)
                    
                    row_colors_df = pd.Series(df_square.index, index=df_square.index).apply(lambda x: get_colors_series(x, type_keys="parms"))
                    col_colors_df = pd.Series(df_square.columns, index=df_square.columns).apply(lambda x: get_colors_series(x, type_keys="test"))


                    # define the col clustering
                    col_cluster = True
                    row_cluster = True

                    # define the annotations
                    df_annotations = df[["parms_idx", "test_idx", "label"]].pivot(index='parms_idx', columns='test_idx', values="label")

                    # define the filename
                    filename = "%s/cross_accuracy_%s_%s_%s_%s.pdf"%(plots_dir, accuracy, ploidy_tag, svtype_tag, typeSimulations_tag)
                    print_if_verbose("getting %s"%filename)

                    # define the title
                    title = "%s when running the best filters according for each sample/condition (rows) tested on each simulation (columns)"%accuracy

                    # define the figure size
                    figsize = (int(len(df_square.columns)*0.03), int(len(df_square)*0.03))
                    #figsize = None

                    plot_clustermap_with_annotation(df_square, row_colors_df, col_colors_df, filename, title=title, col_cluster=col_cluster, row_cluster=row_cluster, colorbar_label=accuracy, adjust_position=True, legend=True, idxs_separator_pattern="||||", texts_to_strip={"L001"}, default_label_legend="control", df_annotations=df_annotations, cmap=sns.color_palette("RdBu_r", 50), ylabels_graphics_df=None, grid_lines=False, figsize=figsize)

def generate_boxplot_comparing_cross_accuracy(df_benchmark, plots_dir):

    """This function takes a cross-benchmarking df from plot_accuracy_of_parameters_on_test_samples and generates a set of boxplots, each of them with one accuracy measurement. The x will be the svtype, and the hue the type of comparison (fastSVcalling, different_sample, different_run, same_run) """

    print_if_verbose("plotting cross-accuracy boxplot")

    df_benchmark = cp.deepcopy(df_benchmark)

    # define a shortened version of svtype
    svtype_to_shortSVtype = {"deletions":"del", "tandemDuplications":"tan", "insertions":"ins", "translocations":"tra", "inversions":"inv", "integrated":"all", "remaining":"rem"}
    df_benchmark["svtype"] = df_benchmark.svtype.apply(lambda x: svtype_to_shortSVtype[x])

    #df_benchmark = df_benchmark.iloc[0:1000].append(df_benchmark[df_benchmark.parms_typeSimulations=="fast"].iloc[0:100]) # debug

    # add the type of comparison
    def get_type_comparison(r):

        if r["parms_typeSimulations"]=="fast": return "default parameters"
        elif r["parms_sampleID"]!=r["test_sampleID"]: return "optimised on different taxID"
        elif r["parms_sampleID"]==r["test_sampleID"] and r["parms_runID"]!=r["test_runID"]: return "optimised on different run"
        elif r["parms_sampleID"]==r["test_sampleID"] and r["parms_runID"]==r["test_runID"]: return "optimised on same run"
        else: raise ValueError("The row is not valid")

    df_benchmark["type comparison"] = df_benchmark.apply(get_type_comparison, axis=1)

    print_if_verbose("plotting")

    #label_to_ylabel = {"fraction overlapping SVs": "fraction overlap. SVs ~ precision" , "n SVs":"n SVs ~ recall"}

    for parms_tag, interesting_type_comps in [["allComparisons", set(df_benchmark["type comparison"])], ["noDefault", set(df_benchmark["type comparison"]).difference({"default parameters"})]]:

        # filter
        df_plot = df_benchmark[df_benchmark["type comparison"].isin(interesting_type_comps)]

        fig = plt.figure(figsize=(len(set(df_benchmark.svtype)), 8))

        for I, y in enumerate(["precision", "recall", "Fvalue"]): # 
            print_if_verbose(y)

            ax = plt.subplot(3, 1, I+1)

            # get a violin plot
            #ax = sns.boxplot(x="svtype", y=y, data=df_plot, hue="type comparison", boxprops=dict(alpha=.45))
            ax = sns.violinplot(x="svtype", y=y, data=df_plot, hue="type comparison", boxprops=dict(alpha=.9))

            #ax = sns.swarmplot(x="svtype", y=y, hue="type comparison", data=df_plot, dodge=True, linewidth=.5, edgecolor="k")
            ax = sns.stripplot(x="svtype", y=y, hue="type comparison", data=df_plot, dodge=True, linewidth=.1, edgecolor="k", size=2)

            ax.legend(bbox_to_anchor=(1, 1))
            ax.set_xlabel("")
            #ax.set_ylabel(label_to_ylabel[y])

            if I in [1,2]: ax.get_legend().remove()


        # save
        print_if_verbose("saving")
        filename = "%s/cross_benchmarking_boxplots_%s.pdf"%(plots_dir, parms_tag)
        fig.savefig(filename, bbox_inches='tight')
        plt.close(fig)

def plot_accuracy_of_parameters_on_test_samples(parameters_df, test_df, outdir, plots_dir, replace=False, threads=4):

    """
    This function runs the gridss+clove pipeline from several parameters (specified in parameters_df) on other datasets (specified in test_df). All the files will be written under outdir. At the end, a heatmap with the accuracy of each of the parameters on each of the test datasets will be generated. The index of each df should be a unique tuple indicating the cells in the final heatmap.
    """

    if replace is True: delete_folder(outdir)
    make_folder(outdir)

    # define the metadata of each df
    parameters_df_metadata = [k for k in parameters_df.keys() if k not in {"parameters_json"}]
    test_df_metadata = [k for k in test_df.keys() if k not in {"sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix"}]

    # keep only one "fast" if typeSimulations is there:
    if "typeSimulations" in parameters_df.keys():

        parameters_df_noFast = parameters_df[parameters_df.typeSimulations!="fast"]
        parameters_df_Fast = parameters_df[parameters_df.typeSimulations=="fast"].iloc[0]
        parameters_df = parameters_df_noFast.append(parameters_df_Fast)

    # add the parameters as a dict
    parameters_df["parameters_json_dict"] = parameters_df.parameters_json.apply(get_parameters_from_json)

    # map each parameterID to the equivalent parameters
    parmID_to_equal_parmIDs = {parmID : {other_parmID for other_parmID, r_other in parameters_df.iterrows() if r_other["parameters_json_dict"]==r["parameters_json_dict"] and other_parmID!=parmID} for parmID, r in parameters_df.iterrows()}
    for parmID, equal_parmIDs in parmID_to_equal_parmIDs.items(): 
        if len(equal_parmIDs)>0: print_if_verbose("%s and %s have equal parms"%(parmID, equal_parmIDs))

    ######## CREATE A df_benchmarking CONTAINING ALL THE RESULTS ##########

    # define the outdir
    outdir_cross_benchmark_files = "%s/tmp_files"%outdir; make_folder(outdir_cross_benchmark_files)

    # define the benchmarking file
    df_benchmark_file = "%s/benchmarking_parameters.tab"%outdir

    if file_is_empty(df_benchmark_file):

        # initialize the df of the benchmarking
        benchmarking_fields = ['FN', 'FP', 'Fvalue', 'TP', 'nevents', 'precision', 'recall', 'svtype']
        df_benchmark = pd.DataFrame(columns=["parms_%s"%x for x in parameters_df_metadata] + ["test_%s"%x for x in test_df_metadata] + benchmarking_fields)

        for numeric_parameter_index, (Irow, parms_row) in enumerate(parameters_df.iterrows()):
            Irow_str = "_".join(Irow)
            
            # get the parameters
            gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup = get_parameters_from_json(parms_row["parameters_json"])

            for numeric_test_index, (Itest, test_row) in enumerate(test_df.iterrows()):
                Itest_str = "_".join(Itest)
                print_if_verbose("\n\n---------\nusing best parameters by %s (%i/%i)"%(Irow_str, numeric_parameter_index+1, len(parameters_df)))
                print_if_verbose("testing on %s (%i/%i)"%(Itest_str, numeric_test_index+1, len(test_df)))

                # define an outdir and put the gridss vcf there with a softlink
                outdir_cross_benchmark = "%s/%s_parameters_tested_on_%s"%(outdir_cross_benchmark_files, Irow_str, Itest_str); make_folder(outdir_cross_benchmark)

                # define the test data file
                df_benchmark_test_file = "%s/df_benchmark_test.py"%outdir_cross_benchmark

                # generate this file through softlinking of equivalent parameterIDs
                equal_parmIDs = parmID_to_equal_parmIDs[Irow]
                if len(equal_parmIDs)>0:

                    for equal_parmID in equal_parmIDs:

                        # define the equal parmID test dir
                        equal_parmID_outdir_cross_benchmark = "%s/%s_parameters_tested_on_%s"%(outdir_cross_benchmark_files, "_".join(equal_parmID), Itest_str); make_folder(outdir_cross_benchmark)

                        # define the equal parm ID df_benchmark_test 
                        equal_parmID_df_benchmark_test_file = "%s/df_benchmark_test.py"%equal_parmID_outdir_cross_benchmark

                        if not file_is_empty(equal_parmID_df_benchmark_test_file): 
                            print_if_verbose("Taking df benchmark from %s for %s"%("_".join(equal_parmID), Irow_str))
                            soft_link_files(equal_parmID_df_benchmark_test_file, df_benchmark_test_file)

                # get this file if not done
                if file_is_empty(df_benchmark_test_file) or replace is True:

                    # define the outdir
                    gridss_vcf = "%s/gridss_vcf.vcf"%outdir_cross_benchmark
                    soft_link_files(test_row["gridss_vcf"], gridss_vcf)

                    # calculate the median insert sizes
                    median_insert_size, median_insert_size_sd  = get_insert_size_distribution(test_row["sorted_bam"], replace=replace, threads=threads)

                    # get the gridss-clove run
                    sv_dict, df_gridss = run_gridssClove_given_filters(test_row["sorted_bam"], test_row["reference_genome"], outdir_cross_benchmark, -1, replace=replace, threads=threads, gridss_blacklisted_regions=gridss_blacklisted_regions, gridss_VCFoutput=gridss_vcf, gridss_maxcoverage=gridss_maxcoverage, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, gridss_filters_dict=gridss_filters_dict, run_in_parallel=True, max_rel_coverage_to_consider_del=max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup=min_rel_coverage_to_consider_dup, replace_FromGridssRun=replace)

                    # get the known SV dict
                    known_sv_dict = {svtype : "%s_%s.tab"%(test_row["svtables_prefix"], svtype) for svtype in {"insertions", "deletions", "translocations", "inversions", "tandemDuplications"}}
                    if any([file_is_empty(x) for x in known_sv_dict.values()]): raise ValueError("There are some un existing files")

                    # get the benchmarking
                    fileprefix = "%s/benchmarking"%outdir_cross_benchmark
                    df_benchmark_test = benchmark_processedSVs_against_knownSVs_inHouse(sv_dict, known_sv_dict, fileprefix, replace=replace, add_integrated_benchmarking=True)

                    # save
                    print_if_verbose("saving")
                    save_object(df_benchmark_test, df_benchmark_test_file)

                else: df_benchmark_test = load_object(df_benchmark_test_file)

                # add metadata fields
                for x in parameters_df_metadata: df_benchmark_test["parms_%s"%x] = parms_row[x]
                for x in test_df_metadata: df_benchmark_test["test_%s"%x] = test_row[x]

                # keep
                df_benchmark = df_benchmark.append(df_benchmark_test[list(df_benchmark.keys())])

        # save
        print_if_verbose("saving")
        df_benchmark_file_tmp = "%s.tmp"%df_benchmark_file
        df_benchmark.to_csv(df_benchmark_file_tmp, sep="\t", index=False, header=True)
        os.rename(df_benchmark_file_tmp, df_benchmark_file)

    else: df_benchmark = pd.read_csv(df_benchmark_file, sep="\t")

    # delete the folder with all the intermediate files
    delete_folder(outdir_cross_benchmark_files)

    #######################################################################

    # plot boxplot
    generate_boxplot_comparing_cross_accuracy(df_benchmark, plots_dir)

    # plot heatmap of cross accuracy
    generate_heatmap_accuracy_of_parameters_on_test_samples(df_benchmark, plots_dir, replace=replace, threads=threads)

  

def report_accuracy_realSVs(close_shortReads_table, reference_genome, outdir, real_svtype_to_file, SVs_compatible_to_insert_dir, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", nvars=100, job_array_mode="local", max_ncores_queue=48, time_perSVade_running="02:00:00", queue_jobs="debug", StopAfter_testAccuracy_perSVadeRunning=False, skip_cleaning_simulations_files_and_parameters=False, skip_cleaning_outdir=False, slurm_constraint=""):


    """This function runs the SV pipeline for all the datasets in close_shortReads_table with the fastSV, optimisation based on uniform parameters and optimisation based on realSVs (specified in real_svtype_to_file). The latter is skipped if real_svtype_to_file is empty.

    First, it runs perSVade on all parameters without cleaning. At the end it cleans."""

    # this pipeline requires real data and close_shortReads_table that is not none
    if len(real_svtype_to_file)==0: raise ValueError("You need real data if you want to test accuracy")
    if file_is_empty(close_shortReads_table): raise ValueError("You need real data reads if you want to test accuracy")

    # make the outdir
    make_folder(outdir)

    # make a plots dir
    plots_dir = "%s/plots"%outdir; make_folder(plots_dir)

    print_if_verbose("testing the accuracy of perSVade")

    # load the real data table
    df_reads = pd.read_csv(close_shortReads_table, sep="\t").set_index("runID", drop=False)

    # define the outfiles
    all_sampleID_to_dfBestAccuracy_file = "%s/all_sampleID_to_dfBestAccuracy.py"%outdir
    ID_to_svtype_to_svDF_file = "%s/ID_to_svtype_to_svDF.py"%outdir

    if file_is_empty(ID_to_svtype_to_svDF_file) or file_is_empty(all_sampleID_to_dfBestAccuracy_file) or replace is True:

        # initialize a dict that will contain all the data
        all_sampleID_to_svtype_to_file = {}
        all_sampleID_to_dfGRIDSS = {}
        all_sampleID_to_dfBestAccuracy = {}

        # initialize the cmds to run 
        all_cmds = []


        # predefine if some jobs need to be ran
        n_remaining_jobs = sum([sum([file_is_empty("%s/%s/%s/perSVade_finished_file.txt"%(outdir, typeSimulations, runID)) for runID in set(df_reads.runID)]) for typeSimulations in ["uniform", "fast", "realSVs"]])
        print_if_verbose("There are %i remaining jobs"%n_remaining_jobs)

        # go through each run and configuration
        for typeSimulations, svtype_to_svfile, fast_SVcalling in [("uniform", {}, False), ("realSVs", real_svtype_to_file, False), ("fast", {}, True)]:

            # define an outdir for this type of simulations
            outdir_typeSimulations = "%s/%s"%(outdir, typeSimulations); make_folder(outdir_typeSimulations)

            # go though each runID
            for runID in set(df_reads.runID):
                print_if_verbose(typeSimulations, runID)

                # define an outdir for this runID
                outdir_runID = "%s/%s"%(outdir_typeSimulations, runID); make_folder(outdir_runID)

                # define the reads
                r1 = df_reads.loc[runID, "short_reads1"]
                r2 = df_reads.loc[runID, "short_reads2"]

                # define the final file 
                final_file = "%s/perSVade_finished_file.txt"%outdir_runID
                parameters_file = "%s/SVdetection_output/final_gridss_running/perSVade_parameters.json"%outdir_runID

                # define the previous repeats file 
                previous_repeats_table = "%s.repeats.tab"%reference_genome
                if file_is_empty(previous_repeats_table): raise ValueError("%s should exist"%previous_repeats_table)
                
                # only contine if the final file is not defined
                if file_is_empty(final_file) or replace is True:# or file_is_empty(parameters_file):

                    # define the cmd. This is a normal perSvade.py run with the vars of the previous dir  
                    cmd = "python %s -r %s --threads %i --outdir %s --nvars %i --nsimulations %i --simulation_ploidies %s --range_filtering_benchmark %s --mitochondrial_chromosome %s -f1 %s -f2 %s --previous_repeats_table %s --skip_cleaning_outdir"%(perSVade_py, reference_genome, threads, outdir_runID, nvars, n_simulated_genomes, ",".join(simulation_ploidies), range_filtering_benchmark, mitochondrial_chromosome, r1, r2, previous_repeats_table)

                    # add arguments depending on the pipeline
                    if replace is True: cmd += " --replace"
                    if fast_SVcalling is True: cmd += " --fast_SVcalling"
                    if len(svtype_to_svfile)>0: cmd += " --SVs_compatible_to_insert_dir %s"%SVs_compatible_to_insert_dir
                    if printing_verbose_mode is True: cmd += " --verbose"

                    # if the running in slurm is false, just run the cmd
                    if job_array_mode=="local": run_cmd(cmd)
                    elif job_array_mode=="greasy": 
                        all_cmds.append(cmd)
                        continue

                    else: raise ValueError("%s is not valid"%job_array_mode)

                if StopAfter_testAccuracy_perSVadeRunning is True: continue

                # define the svdict and the df_gridss 
                svtype_to_svfile, df_gridss = get_svtype_to_svfile_and_df_gridss_from_perSVade_outdir(outdir_runID, reference_genome)

                # add to dict
                ID = "%s||||%s"%(typeSimulations, runID)
                all_sampleID_to_svtype_to_file[ID] = svtype_to_svfile
                all_sampleID_to_dfGRIDSS[ID] = df_gridss

                # get the df best accuracy
                if typeSimulations!="fast": all_sampleID_to_dfBestAccuracy[ID] = pd.read_csv("%s/SVdetection_output/parameter_optimisation/benchmarking_all_filters_for_all_genomes_and_ploidies/df_cross_benchmark_best.tab"%outdir_runID, sep="\t")

        # if you are not running on slurm, just execute one cmd after the other
        if job_array_mode=="greasy":

            if len(all_cmds)>0: 
                print_if_verbose("submitting %i jobs to the cluster for testing accuracy of perSVade on several combinations of parameters. The files of the submission are in %s"%(len(all_cmds), outdir))
                jobs_filename = "%s/jobs.testingRealDataAccuracy"%outdir
                open(jobs_filename, "w").write("\n".join(all_cmds))

                generate_jobarray_file_greasy(jobs_filename, walltime=time_perSVade_running,  name="testAccuracy", queue=queue_jobs, sbatch=True, ncores_per_task=threads, constraint=slurm_constraint, number_tasks_to_run_at_once="all", max_ncores_queue=max_ncores_queue )

                print_if_verbose("You have to wait under all the jobs in testRealSVs are done")
                sys.exit(0)


        if StopAfter_testAccuracy_perSVadeRunning is True: 
            print_if_verbose("You already ran all the configurations of perSVade. Stopping after the running of perSVade on testAccuracy")
            sys.exit(0)

        print_if_verbose("getting ID_to_svtype_to_svDF")
        ID_to_svtype_to_svDF = get_sampleID_to_svtype_to_svDF_filtered(all_sampleID_to_svtype_to_file, all_sampleID_to_dfGRIDSS)

        # add the 'svID', which is useful to calculate overlaps
        print_if_verbose("adding svID")
        add_svID_to_IDtoSVTYPEtoDF(ID_to_svtype_to_svDF)

        # save
        save_object(ID_to_svtype_to_svDF, ID_to_svtype_to_svDF_file)
        save_object(all_sampleID_to_dfBestAccuracy, all_sampleID_to_dfBestAccuracy_file)

    else: 
        print_if_verbose("loading objects")
        ID_to_svtype_to_svDF = load_object(ID_to_svtype_to_svDF_file)
        all_sampleID_to_dfBestAccuracy = load_object(all_sampleID_to_dfBestAccuracy_file) 

    if StopAfter_testAccuracy_perSVadeRunning is True: 
        print_if_verbose("You already ran all the configurations of perSVade. Stopping after the running of perSVade on testAccuracy")
        sys.exit(0)

    print_if_verbose("getting the plots about accuracy")
    # map each runID to the IDs of the same sample 
    runID_to_replicateIDs = {runID : set(df_reads[df_reads.sampleID==df_reads.loc[runID, "sampleID"]].index).difference({runID}) for runID in df_reads.runID}

    # map each ID to the svIDs 
    ID_to_svIDs = {ID : set.union(*[set(svDF.svID) for svDF in svtype_to_svDF.values() if len(svDF)>0]) for ID, svtype_to_svDF in ID_to_svtype_to_svDF.items()}

    # define the fraction of samples with this ID
    svID_to_fractionIDsPresent = {svID : sum([svID in svIDs for svIDs in ID_to_svIDs.values()])/len(ID_to_svIDs) for svID in set.union(*ID_to_svIDs.values())}

    # define the wrong svIDs as those that are present in >75% of the samples. These are likely errors in the assembly
    wrong_svIDs = {svID for svID, fraction in svID_to_fractionIDsPresent.items() if fraction>0.75}

    # map each runID to the svtype to the svIDs
    ID_to_svtype_to_svIDs = {ID : {svtype : set(svDF.svID).difference(wrong_svIDs) for svtype, svDF in svtype_to_svDF.items() if len(svDF)>0} for ID, svtype_to_svDF in ID_to_svtype_to_svDF.items()}

    # initialize a benchmarking dict
    df_benchmarking_realSVs_dict = {}

    for ID, svtype_to_svIDs in ID_to_svtype_to_svIDs.items():

        # initialize the total numbers
        nSVs_list = []
        fraction_overlapping_list = []

        for svtype, svIDs in svtype_to_svIDs.items():

            # define the IDs
            simulationID, runID = ID.split("||||")

            # only keep data if there are more than 10 svIDs
            if len(svIDs)>=0: 

                # define the svIDs in other runs of the same sample
                list_other_svIDs = [ID_to_svtype_to_svIDs["%s||||%s"%(simulationID, otherRunID)][svtype] for otherRunID in runID_to_replicateIDs[runID] if svtype in ID_to_svtype_to_svIDs["%s||||%s"%(simulationID, otherRunID)]] + [set()]

                # define all the others
                all_other_svIDs = set.union(*list_other_svIDs)

                # define the true even
                intersection_other_svIDs = set.intersection(*list_other_svIDs)

                # define the overlap as the ratio 
                n_SVs = len(svIDs)
                n_overlapping = len(all_other_svIDs.intersection(svIDs))
                if n_SVs>0: overlap_SVs = n_overlapping/len(svIDs)
                else: overlap_SVs = 0

                # define the 'recall' of real vars
                TPs = intersection_other_svIDs.intersection(svIDs)
                #recall = len(TPs)/len(intersection_other_svIDs)

                # keep
                IDdict = "%s||||%s"%(ID, svtype)

                #df_benchmarking_realSVs_dict[IDdict] = {"simulationID":simulationID, "runID":runID, "sampleID":df_reads.loc[runID, "sampleID"], "svtype":svtype, "fraction overlapping SVs":overlap_SVs, "n_overlapping":n_overlapping, "n SVs":n_SVs, "n_HighConfidence_SVs":len(intersection_other_svIDs)}   

                #keep
                nSVs_list.append(n_SVs)
                fraction_overlapping_list.append(overlap_SVs)

                # only keep important things
                df_benchmarking_realSVs_dict[IDdict] = {"simulationID":simulationID, "runID":runID, "sampleID":df_reads.loc[runID, "sampleID"], "svtype":svtype, "fraction overlapping SVs":overlap_SVs, "n SVs":n_SVs}  


        # add the integrated vars
        integrated_n_SVs = sum(nSVs_list)
        integeated_overlap_SVs = np.mean(make_flat_listOflists([[fraction_overlapping_list[Isvtype]]*nSVs for Isvtype, nSVs in enumerate(nSVs_list)]))

        df_benchmarking_realSVs_dict["%s||||integrated"%ID] = {"simulationID":simulationID, "runID":runID, "sampleID":df_reads.loc[runID, "sampleID"], "svtype":"integrated", "fraction overlapping SVs":integeated_overlap_SVs, "n SVs":integrated_n_SVs}  

    df_benchmarking_realSVs = pd.DataFrame(df_benchmarking_realSVs_dict).transpose()

    # plot the fraction of overlapping real SVs as a boxplot
    filename = "%s/accuracy_realSVs.pdf"%plots_dir
    plot_fraction_overlapping_realSVs(df_benchmarking_realSVs, filename)

    # plot the fraction of overlapping real SVs as a precision-vs-recall plot
    filename = "%s/accuracy_realSVs_precisionVSrecall.pdf"%plots_dir
    plot_fraction_overlapping_realSVs_precision_vs_recall(df_benchmarking_realSVs, filename)
 
    # plot the accuracy on simulations
    filename = "%s/accuracy_simulations.pdf"%plots_dir
    plot_accuracy_simulations_from_all_sampleID_to_dfBestAccuracy(all_sampleID_to_dfBestAccuracy, filename)

    ########## plot the cross-accuracy between different parameters ##########

    # initialize a dir that will contain the sorted bams and gridss vcf outputs
    simulations_files_and_parameters_dir = "%s/simulations_files_and_parameters"%outdir
    if replace is True: delete_folder(simulations_files_and_parameters_dir)
    make_folder(simulations_files_and_parameters_dir)

    # This requires the definition of two things:

    # the parameters_df. The first cols are metadata (like sampleID, runID and optimisation type) and the others are things necessary for runnning gridss: and the path to the parameters_json
    parameters_df_dict = {}

    # test_df: This is info on which to test the running of gridss+clove. It contains metadata cols (sampleID, runID, optimisation type (real, uniform), simName, ploidy, svtype) and data to run the optimisation on (sorted_bam, gridss_vcf, reference_genome, mitochondrial_chromosome)
    test_df_dict = {}

    for typeSimulations in ["uniform", "realSVs", "fast"]:

        # go through each sampleID
        for sampleID in sorted(set(df_reads.sampleID)):

            # go though each runID
            for runID in sorted(set(df_reads[df_reads.sampleID==sampleID].runID)):

                # define the outdir of the run
                runID_outdir = "%s/%s/%s/SVdetection_output"%(outdir, typeSimulations, runID)

                # get the parameters
                parameters_json_origin = "%s/final_gridss_running/perSVade_parameters.json"%runID_outdir
                parameters_json_dest = "%s/parameters_%s_%s_%s.json"%(simulations_files_and_parameters_dir, typeSimulations, sampleID, runID)
                if file_is_empty(parameters_json_dest): 


                    print(parameters_json_origin, parameters_json_dest)
                    parameters_json_dest_tmp = "%s.tmp"%parameters_json_dest
                    run_cmd("cp %s %s"%(parameters_json_origin, parameters_json_dest_tmp))
                    os.rename(parameters_json_dest_tmp, parameters_json_dest)

                parameters_df_dict[(sampleID, runID, typeSimulations)] = {"sampleID":sampleID, "runID":runID, "typeSimulations":typeSimulations, "parameters_json":parameters_json_dest}

                # go through additional things
                if typeSimulations!="fast":
                    for Isim in range(n_simulated_genomes):

                        # define the name
                        simName = "simulation_%i"%(Isim+1)
                        sim_outdir = "%s/parameter_optimisation/%s"%(runID_outdir, simName)

                        # go through each ploidy
                        for ploidy in simulation_ploidies:

                            # define the destintaion bam 
                            destination_bam = "%s/reads_%s_%s_%s_sim%i_%s.bam"%(simulations_files_and_parameters_dir, typeSimulations, sampleID, runID, Isim+1, ploidy)

                            # change the place
                            if file_is_empty(destination_bam): 

                                # define the origin bam
                                if ploidy=="haploid": suffix_ploidy = "bam.sorted"
                                else: suffix_ploidy = "bam.sorted.%s.bam.sorted"%ploidy
                                bam_files = ["%s/%s"%(sim_outdir, file) for file in os.listdir(sim_outdir) if file.startswith("aligned_reads") and ".".join(file.split(".")[1:])==suffix_ploidy]
                                if len(bam_files)!=1: raise ValueError("There should be only one bam")
                                origin_bam = bam_files[0]
   
                                # copy the metrics
                                origin_metrics = "%s/%s/%s/aligned_reads.bam.sorted.CollectInsertSizeMetrics.out"%(outdir, typeSimulations, runID)
                                destination_metrics = "%s.CollectInsertSizeMetrics.out"%destination_bam
                                destination_metrics_tmp = "%s.tmp"%destination_metrics
                                if file_is_empty(destination_metrics):
                                    run_cmd("cp %s %s"%(origin_metrics, destination_metrics_tmp))
                                    os.rename(destination_metrics_tmp, destination_metrics)

                                # change the bai
                                if file_is_empty("%s.bai"%destination_bam): os.rename("%s.bai"%origin_bam, "%s.bai"%destination_bam)

                                # change the coverage per window (this is any destination)
                                if file_is_empty("%s.coverage_per_window.tab"%destination_bam): os.rename("%s.coverage_per_window.tab"%origin_bam, "%s.coverage_per_window.tab"%destination_bam)

                                # change the coverage per constant windows
                                calculating_windowcoverage_dir = "%s.calculating_windowcoverage"%destination_bam; make_folder(calculating_windowcoverage_dir)
                                destination_windowcoverage_file = "%s/coverage_windows_%ibp.tab"%(calculating_windowcoverage_dir, window_l)
                                origin_windowcoverage_file = "%s/benchmark_GridssClove_%s/coverage_windows_%ibp.tab"%(sim_outdir, ploidy, window_l)
                                if file_is_empty(destination_windowcoverage_file): os.rename(origin_windowcoverage_file, destination_windowcoverage_file)

                                # rename the bam
                                os.rename(origin_bam, destination_bam)

                            # change the gridss vcf
                            origin_gridss_vcf = "%s/benchmark_GridssClove_%s/benchmark_max50000x_ignoreRegionsFalse/gridss_output.vcf.withSimpleEventType.vcf"%(sim_outdir, ploidy)
                            dest_gridss_vcf = "%s/gridss_vcf_%s_%s_%s_sim%i_%s.vcf"%(simulations_files_and_parameters_dir, typeSimulations, sampleID, runID, Isim+1, ploidy)
                            if file_is_empty(dest_gridss_vcf): os.rename(origin_gridss_vcf, dest_gridss_vcf)

                            # change the location of the simulated SVs
                            svtables_prefix =  "%s/SVs_%s_%s_%s_sim%i"%(simulations_files_and_parameters_dir, typeSimulations, sampleID, runID, Isim+1)
                            for svtype in {"insertions", "deletions", "translocations", "inversions", "tandemDuplications"}:

                                # define the files
                                origin_file = "%s/final_simulated_SVs/%s.tab"%(sim_outdir, svtype)
                                dest_file = "%s_%s.tab"%(svtables_prefix, svtype)
                                dest_file_tmp = "%s.tmp"%dest_file

                                # move
                                if file_is_empty(dest_file):
                                    run_cmd("cp %s %s"%(origin_file, dest_file_tmp))
                                    os.rename(dest_file_tmp, dest_file)

                            # get the name
                            test_df_dict[(sampleID, runID, typeSimulations, simName, ploidy)] = {"sampleID":sampleID, "runID":runID, "typeSimulations":typeSimulations, "simName":simName, "ploidy":ploidy, "sorted_bam":destination_bam, "gridss_vcf":dest_gridss_vcf, "reference_genome":reference_genome, "mitochondrial_chromosome":mitochondrial_chromosome, "svtables_prefix":svtables_prefix}

    # get the dfs
    parameters_df = pd.DataFrame(parameters_df_dict).transpose()[["sampleID", "runID", "typeSimulations", "parameters_json"]]
    test_df = pd.DataFrame(test_df_dict).transpose()[["sampleID", "runID", "typeSimulations", "simName", "ploidy", "sorted_bam", "gridss_vcf", "reference_genome", "mitochondrial_chromosome", "svtables_prefix"]]

    # plot the cross-accuracy
    print_if_verbose("plotting cross-accuracy")
    outdir_cross_accuracy = "%s/cross_accuracy_calculations"%outdir
    plot_accuracy_of_parameters_on_test_samples(parameters_df, test_df, outdir_cross_accuracy, plots_dir, replace=replace, threads=threads)

    ##########################################################################

    ##### VERY IMPORTANT: CLEAN THE OUTPUT #####

    # clean the simulations_files_and_parameters if not stated otherwise
    if skip_cleaning_simulations_files_and_parameters is False: 
        print_if_verbose("cleaning simulations_files_and_parameters_dir")
        delete_folder(simulations_files_and_parameters_dir)

    if skip_cleaning_outdir is False:
        print_if_verbose("cleaning outdir of each run")

        # clean the outdir of each sample
        for typeSimulations in ["uniform", "realSVs", "fast"]:

            # go though each runID
            for runID in set(df_reads.runID):
                print(typeSimulations, runID)

                # define an outdir for this runID and clean
                outdir_runID = "%s/%s/%s"%(outdir, typeSimulations, runID)
                clean_perSVade_outdir(outdir_runID)

                # remove the bam files
                for f in ["aligned_reads.bam.sorted", "aligned_reads.bam.sorted.bai"]: remove_file("%s/%s"%(outdir_runID, f))

    #############################################

def get_simulated_bamFile(outdir, reference_genome, replace=False, threads=4, total_nread_pairs=10000000, read_length=150, median_insert_size=500, median_insert_size_sd=50):

    """This function simulates reads and aligns them for a reference genomes in an uniform way to coverage 50 """

    # simulate a 'windows' df were each chrom is a window
    df_genome_info = pd.DataFrame({chrom : {"chromosome":chrom, "start":0, "end":lenChrom, "predicted_relative_coverage":1.0} for chrom, lenChrom in get_chr_to_len(reference_genome).items()}).transpose()

    # define an outdir
    outdir_sim = "%s/simulatingReadsFromReference"%outdir; make_folder(outdir_sim)

    # get the simulated bam
    simulated_sorted_bam = simulate_and_align_PairedReads_perWindow(df_genome_info, reference_genome, reference_genome, total_nread_pairs, read_length, outdir_sim, median_insert_size, median_insert_size_sd, replace=replace, threads=threads)

    # define the index of this bam
    index_bam = "%s.bai"%simulated_sorted_bam

    return simulated_sorted_bam, index_bam


def get_short_and_long_reads_sameBioSample(outdir, taxID, reference_genome, replace=False, min_coverage=30):

    """This function takes a taxID and looks in the SRA database if there are short paired Illumina reads and MinION ONT reads for the same BioSample, returning the SRRs."""

    make_folder(outdir)

    # get the WGS dataset for all the datasets under this taxID
    fileprefix = "%s/wgs_data"%outdir
    wgs_runInfo_df = get_allWGS_runInfo_fromSRA_forDivision(fileprefix, taxID, reference_genome, replace=False, min_coverage=min_coverage).set_index("Run", drop=False)

    # get the oxford nanopore data
    fileprefix = "%s/ONT_data"%outdir
    ONT_runInfo_df = get_allOxfordNanopore_runInfo_fromSRA_forDivision(fileprefix, taxID, reference_genome, replace=replace, min_coverage=min_coverage)

    # get intersecting biosamples
    intersecting_biosamples = set(wgs_runInfo_df.BioSample).intersection(set(ONT_runInfo_df.BioSample))

    # debug the fact that there are empty datasets
    if len(intersecting_biosamples)==0: raise ValueError("There are no intersecting BioSamples")

    # map each biosample to the ilumina coverage
    interseting_BioSamples = wgs_runInfo_df[wgs_runInfo_df.BioSample.isin(intersecting_biosamples)].sort_values("expected_coverage", ascending=False)["BioSample"]

    # get the best SRR
    for bioSample in interseting_BioSamples:
        print_if_verbose("\n\n", bioSample)

        # get the df for each
        df_wgs = wgs_runInfo_df[wgs_runInfo_df.BioSample==bioSample].sort_values("expected_coverage", ascending=False)
        df_ONT = ONT_runInfo_df[ONT_runInfo_df.BioSample==bioSample].sort_values("expected_coverage", ascending=False)

        wgs_run = df_wgs.iloc[0]["Run"]
        ONT_run = df_ONT.iloc[0]["Run"]

        print_if_verbose("These are the Illumina reads:\n:", df_wgs[["Run", "SampleName"]].iloc[0])
        print_if_verbose("These are the ONT reads:\n:", df_ONT[["Run", "SampleName"]].iloc[0])

        break

    return wgs_run, ONT_run


def report_accuracy_golden_set(goldenSet_dir, outdir, reference_genome, real_svtype_to_file, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", nvars=100, job_array_mode="local", time_read_obtention="02:00:00", time_perSVade_running="48:00:00", queue_jobs="debug", StopAfterPrefecth_of_reads=False, StopAfter_sampleIndexingFromSRA=False, max_ncores_queue=768, target_taxID=None, min_coverage=30):

    """This function takes a directory that has the golden set vars and generates plots reporting the accuracy. If auto, it will find them in the SRA and write them under outdir."""

    print_if_verbose("calculating accuracy for golden set SVcalls")
    make_folder(outdir)

    ### automatic obtention of golden set reads ###
    if goldenSet_dir=="auto":

        # create this dir un
        goldenSet_dir = "%s/automatic_obtention_goldenSetReads"%outdir; make_folder(goldenSet_dir)

        #### define the SRRs to download ####

        # if the tax ID is in taxID_to_srrs_goldenSet, get it
        if target_taxID in taxID_to_srrs_goldenSet: 
            short_reads_srr = taxID_to_srrs_goldenSet[target_taxID]["short_reads"]
            long_reads_srr = taxID_to_srrs_goldenSet[target_taxID]["long_reads"]

        else: short_reads_srr, long_reads_srr = get_short_and_long_reads_sameBioSample("%s/finding_sameBioSample_srrs"%goldenSet_dir, target_taxID, reference_genome, min_coverage=min_coverage, replace=replace)

        #####################################

        # download each of the reads (raw). Stop after fastqdump
        for type_data, srr in [("illumina_paired", short_reads_srr), ("nanopore", long_reads_srr)]:
            print_if_verbose("Getting raw reads for %s"%type_data)

            # define the outdir
            outdir_srr = "%s/%s"%(goldenSet_dir, srr)

            # define the cmd downloading after the fastq-dump
            cmd = "%s --srr %s --outdir %s --threads %i --stop_after_fastqdump --type_data %s"%(get_trimmed_reads_for_srr_py, srr, outdir_srr, threads, type_data)
            if StopAfterPrefecth_of_reads is True: cmd += " --stop_after_prefetch"

            run_cmd(cmd)

        # define the reads
        longReads = "%s/%s/%s.srr.fastq.gz"%(goldenSet_dir, long_reads_srr, long_reads_srr)
        short_reads1 = "%s/%s/%s.srr_1.fastq.gz"%(goldenSet_dir, short_reads_srr, short_reads_srr)
        short_reads2 = "%s/%s/%s.srr_2.fastq.gz"%(goldenSet_dir, short_reads_srr, short_reads_srr)

    #####################################
    else:

        # define the reads, they are suposed to be called like this
        longReads = "%s/long_reads.fastq.gz"%goldenSet_dir
        short_reads1 = "%s/short_reads_1.fastq.gz"%goldenSet_dir
        short_reads2 = "%s/short_reads_2.fastq.gz"%goldenSet_dir

    if any([file_is_empty(f) for f in [longReads, short_reads1, short_reads2]]): raise ValueError("Your golden dir %s should contain long_reads.fasta, short_reads_1.fastq.gz and short_reads_2.fastq.gz"%goldenSet_dir)

    # trim the long reads
    print_if_verbose("running porechop")
    trimmed_long_reads = run_porechop(longReads,  replace=replace, threads=threads)

    # trim the reads
    print_if_verbose("Running trimmomatic")
    trimmed_reads1, trimmed_reads2 = run_trimmomatic(short_reads1, short_reads2, replace=replace, threads=threads)

    # run svim
    print_if_verbose("running svim")
    outdir_svim = "%s/svim_output"%outdir
    svType_to_file_longReads, sorted_bam_short_longReads, median_coverage_longReads = run_svim(trimmed_long_reads, reference_genome, outdir_svim,  threads=threads, replace=replace, aligner="ngmlr", is_nanopore=True)

    # run perSVade with all the types of optimisat


def remove_smallVarsCNV_nonEssentialFiles(outdir, ploidy):

    """Removes all the files in outdir that are not essential. The outdir has to be the one of the VarCall outdir. The files to r"""

    # initialize the files to remove
    files_to_remove = ["%s/CNV_results/gene_to_coverage_genes.tab"%outdir, # the genes coverage
                       "%s/CNV_results/gene_to_coverage_regions.tab"%outdir # the regions coverage
                       ]

    # add the bcftools
    bcftools_dir = "%s/bcftools_ploidy%i_out"%(outdir, ploidy)
    HC_dir = "%s/HaplotypeCaller_ploidy%i_out"%(outdir, ploidy)
    fb_dir = "%s/freebayes_ploidy%i_out"%(outdir, ploidy)

    # go through each dir
    for vcfDir in [bcftools_dir, HC_dir, fb_dir]:

        if os.path.isdir(vcfDir):
            for file in os.listdir(vcfDir):

                if file not in {"output.raw.vcf", "output.filt.vcf"}: files_to_remove.append("%s/%s"%(vcfDir, file))

    # go through the files in the outdir and just keep the essential ones
    for f in os.listdir(outdir):
        file = "%s/%s"%(outdir, f)

        if os.path.isfile(file):

            files_to_keep = {"merged_vcfs_allVars_ploidy%i.vcf"%ploidy,
                             "variant_annotation_ploidy%i.tab"%ploidy,
                             "variant_calling_ploidy%i.tab"%ploidy,

                             "variants_atLeast1PASS_ploidy%i.vcf"%ploidy,
                             "variants_atLeast2PASS_ploidy%i.vcf"%ploidy,
                             "variants_atLeast3PASS_ploidy%i.vcf"%ploidy,

                             "variants_atLeast1PASS_ploidy%i.withMultiAlt.vcf"%ploidy,
                             "variants_atLeast2PASS_ploidy%i.withMultiAlt.vcf"%ploidy,
                             "variants_atLeast3PASS_ploidy%i.withMultiAlt.vcf"%ploidy,

                             "variant_calling_stats_ploidy%i_called.tab"%ploidy,
                             "variant_calling_stats_ploidy%i_PASS.tab"%ploidy
                             }


            if f not in files_to_keep: files_to_remove.append(file)

    for f in files_to_remove: remove_file(f)

def get_sortedBam_with_duplicatesMarked(sorted_bam, threads=4, replace=False):

    """This function takes a sorted bam and returns the equivalent with the duplicates marked with picard MarkDuplicates. It also indexes this bam"""

    # define dirs
    sorted_bam_dupMarked = "%s.MarkDups.bam"%sorted_bam
    sorted_bam_dupMarked_tmp = "%s.MarkDups.tmp.bam"%sorted_bam
    sorted_bam_dupMarked_metrics = "%s.MarkDups.metrics"%sorted_bam

    if file_is_empty(sorted_bam_dupMarked) or replace is True:
        print_if_verbose("marking duplicate reads")

        # define the java memory
        #javaRamGb = int(get_availableGbRAM()*fractionRAM_to_dedicate) # At Broad, we run MarkDuplicates with 2GB Java heap (java -Xmx2g) and 10GB hard memory limit
        #javaRamGb = int(get_availableGbRAM()*0.5) # At Broad, we run MarkDuplicates with 2GB Java heap (java -Xmx2g) and 10GB hard memory limit
        #javaRamGb = int(get_availableGbRAM() - 2) # rule of thumb from GATK
        javaRamGb = 4 # this is from a post from 2011, reccommended for a 170Gb RAM

        # define the MAX_RECORDS_IN_RAM
        MAX_RECORDS_IN_RAM = int(250000*javaRamGb*0.8) # 250,000 reads for each Gb given (for SortSam, I don't know if this will work for Picard tools)

        # define the number of MAX_FILE_HANDLES_FOR_READ_ENDS_MAP
        max_nfilehandles = int(subprocess.check_output("ulimit -n", shell=True))
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = int(max_nfilehandles*0.5) # a little lower than ulimit -n

        # SORTING_COLLECTION_SIZE_RATIO is 0.25 by default. If I have memory issues I can reduce this number.

        markduplicates_std = "%s.markingDuplicates.std"%sorted_bam
        print_if_verbose("running MarkDuplicates with %iGb of RAM and %i MAX_FILE_HANDLES_FOR_READ_ENDS_MAP. The std is in %s"%(javaRamGb, MAX_FILE_HANDLES_FOR_READ_ENDS_MAP, markduplicates_std))

        run_cmd("%s -Xmx%ig MarkDuplicates I=%s O=%s M=%s ASSUME_SORT_ORDER=coordinate MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=%i MAX_RECORDS_IN_RAM=%i > %s 2>&1"%(picard_exec, javaRamGb, sorted_bam, sorted_bam_dupMarked_tmp, sorted_bam_dupMarked_metrics, MAX_FILE_HANDLES_FOR_READ_ENDS_MAP, MAX_RECORDS_IN_RAM, markduplicates_std))
        #REMOVE_DUPLICATES=Boolean

        remove_file(markduplicates_std)

        # keep
        os.rename(sorted_bam_dupMarked_tmp, sorted_bam_dupMarked)

    # index the bam with the duplicate reads
    index_sorted_bam_dupMarked = "%s.bai"%sorted_bam_dupMarked
    if file_is_empty(index_sorted_bam_dupMarked) or replace is True:
        index_bam(sorted_bam_dupMarked, threads=threads)

    return sorted_bam_dupMarked

def write_integrated_smallVariantsTable_as_vcf_old(df, filename, ploidy):

    """This function takes a df table with the ouptut of VEP and writes a vcf only with the variants (no annotation)"""

    print_if_verbose("getting vcf intersection")
    # get a df that has unique vars
    df = cp.deepcopy(df)
    df = df.drop_duplicates(subset="#Uploaded_variation")

    # get the vcf df with info
    df["#CHROM"] = df.chromosome
    df["POS"] = df.position
    df["REF"] = df.ref
    df["ALT"] = df.alt

    # get an empty ID and quality
    df["ID"] = df["#Uploaded_variation"]
    df["QUAL"] = "."

    # get the filter as the number of programs that pass the calling
    df["FILTER"] = df.number_PASS_programs.apply(str) + "xPASS"

    # the sample will contain the genotype and the allele frequency
    df["FORMAT"] = "GT:AF"  

    # define the programs
    programs = ["freebayes", "HaplotypeCaller", "bcftools"]

    # add the PASS programs
    print_if_verbose("getting PASS programs")
    df["PASS_programs"] = df.apply(lambda r: [p for p in programs if r["%s_PASS"%p]], axis=1)
    df["PASS_programs_str"] = df.PASS_programs.apply(lambda x: ",".join(x))

    # add the ploidies
    if ploidy==1: df["GT"] = "."
    elif ploidy==2:

        print_if_verbose("adding ploidies")

        # add the ploidies
        df["all_ploidies"] =  df.apply(lambda r: {"/".join(re.split("/|\|", r["%s_GT"%p])) for p in r["PASS_programs"]}, axis=1)
        df["all_ploidies_len"] = df.all_ploidies.apply(len)

        # add the final ploidy depending on the len
        def get_ploidy_diploid(all_ploidies):

            if len(all_ploidies)==1: return next(iter(all_ploidies))
            else: return "."

        df["GT"] = df.all_ploidies.apply(get_ploidy_diploid)

    else: raise ValueError("There is no proper testing on ploidy %i about the representation"%ploidy)

    # get the AF as the mean of the pass programs
    print_if_verbose("getting allele frequency")
    df["AF"] = df.apply(lambda r: "%.4f"%np.mean([r["%s_fractionReadsCoveringThisVariant"%p] for p in r["PASS_programs"]]), axis=1)

    # add the INFO
    print_if_verbose("getting INFO")
    df["GT_eachProgram"] = df.apply(lambda r: ";".join(["%s_GT=%s"%(p, r["%s_GT"%p]) for p in programs]), axis=1)
    df["INFO"] = "PASSALGS=" + df.PASS_programs_str + ";" + df.GT_eachProgram

    # add to df
    df["SAMPLE"] = df.GT + ":" + df.AF

    # write the final vcf
    vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    df_vcf = df[vcf_fields].sort_values(by=["#CHROM", "POS", "REF"])

    if len(df_vcf)!=len(set(df["#Uploaded_variation"])): raise ValueError("There are some duplicated fields")

    # get the vcf content
    vcf_lines = df_vcf.to_csv(sep="\t", header=True, index=False)

    # get the header
    header_lines = ["##fileformat=VCFv4.2",
                    "##perSVade small variant calling pipeline. This is the merged output of freebayes, GATK Haplotype Caller and bcftools for variants that PASS the filters in at least %i algorithms."%min(df.number_PASS_programs),
                    "##FILTER indicates the number of algorithms were this variant was called and PASSed the filters",
                    "##FORMAT includes the GT (genotype) and AF (allele frequency)",
                    "##GT includes the genotype in the case that all the PASS algorithms called the same GT, and '.' otherwise",
                    "##AF includes the mean fraction of reads calling this variant across PASS alorithms",
                    "##INFO includes the name of the algorithms that called this variant (PASSALGS) and the GT of each of these"
                    ]
    
    filename_tmp = "%s.tmp"%filename
    open(filename_tmp, "w").write("\n".join(header_lines) + "\n" + vcf_lines)
    os.rename(filename_tmp, filename)


def convert_NaN_to0(x):
    
    if pd.isna(x): return 0.0
    else: return x
    
def get_nChroms_with_var(r, p):
    
    GT = r["%s_GT"%p]
    index = r["%s_GT_index"%p]
    
    if GT=="" or index=="": return 0
    else: return GT.count(index)

def get_numberChromosomes_withVar(r):
    
    all_PASS_nChroms = {r["%s_numberChromosomes_withVar"%p] for p in r["PASS_programs"]}
    
    # if there are none, just get a 0
    if len(all_PASS_nChroms)==0: return 0
    
    # if there are is one, take it
    elif len(all_PASS_nChroms)==1: return next(iter(all_PASS_nChroms))
    
    # if there are more than one option get the minimum (most conservative one)
    elif len(all_PASS_nChroms)>1: return min(all_PASS_nChroms)
    
    else: raise ValueError("There are incorrect numbers of all_PASS_nChroms")
    
def get_corrected_INFO(x):
    
    return x.replace("freebayes", "fb").replace("HaplotypeCaller", "HC").replace("bcftools", "bt")
    
def get_row_vcf(df_pos, ploidy):
    
    # takes a df with the position and returns a series that has all the vcf fields
                      
    # initialize the row with unique vals
    r = pd.Series()
    r["#CHROM"] = df_pos["#CHROM"].iloc[0]
    r["POS"] = df_pos["POS"].iloc[0]
    r["REF"] = df_pos["REF"].iloc[0]
    r["QUAL"] = "."
    r["FORMAT"] = "GT:AF" 
    
    # initialize with shared vals
    r["ALT"] = ",".join(df_pos["ALT"])
    r["ID"] = ",".join(df_pos["#Uploaded_variation"])
    r["FILTER"] = ",".join(df_pos.FILTER)
    r["PASS_programs_str"] =  ",".join(df_pos["PASS_programs_str"])
    
    # define the index of each var (so that they start on 1)
    df_pos["var_index"] = list(range(1, len(df_pos)+1))
    df_pos["var_index"] = df_pos["var_index"].apply(str)
    
    # define the allele frequency of each program
    programs = ["freebayes", "HaplotypeCaller", "bcftools"]
    r["AF_programs_str"] = ";".join(["%s_AF="%p + ",".join(df_pos["%s_fractionReadsCoveringThisVariant"%p].apply(lambda x: "%.4f"%x)) for p in programs if "%s_called"%p in df_pos.keys()])
    
    # define the info
    r["INFO"] = get_corrected_INFO("PASSALGS=" + r["PASS_programs_str"] + ";" + r.AF_programs_str)

    # initialize a list that will contain the var_index
    p_GTlist = []
    
    # define a df with the GT info
    df_GT = df_pos[df_pos.numberChromosomes_withVar>0]
    if 0<sum(df_GT.numberChromosomes_withVar)<=ploidy:
        
        p_GTlist += make_flat_listOflists(df_GT.apply(lambda x: [x["var_index"]]*x["numberChromosomes_withVar"], axis=1))
        
    # add missing chromosomes
    p_GTlist = ["0"]*(ploidy-len(p_GTlist)) + p_GTlist
    if len(p_GTlist)!=ploidy: raise ValueError("The genotype was not properly calculated")
        
    r["GT"] = "/".join(p_GTlist)
    
    # add the sample
    r["SAMPLE"] = r.GT + ":" + ",".join(df_pos.AF)
    
    vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    return r[vcf_fields]

def get_int_or0(x):

    """Gets an int from a string, or 0"""

    if x in {"", ".", " "}: return 0
    else: return int(x)
        
def write_integrated_smallVariantsTable_as_vcf(df, filename, ploidy):

    """This function takes a df table with the ouptut of VEP and writes a vcf only with the variants (no annotation)"""

    print_if_verbose("getting vcf intersection")
    
    # get a df that has unique vars
    df = cp.deepcopy(df)
    df = df.drop_duplicates(subset="#Uploaded_variation")
                                
    # get the vcf df with info
    df["#CHROM"] = df.chromosome
    df["POS"] = df.position
    df["REF"] = df.ref
    df["ALT"] = df.alt

    # get the filter as the number of programs that pass the calling
    df["FILTER"] = df.number_PASS_programs.apply(str) + "xPASS"

    # define the programs
    programs = ["freebayes", "HaplotypeCaller", "bcftools"]

    # add the PASS programs
    print_if_verbose("getting PASS programs")
    df["PASS_programs"] = df.apply(lambda r: [p for p in programs if r["%s_PASS"%p]], axis=1)
    df["PASS_programs_str"] = df.PASS_programs.apply(lambda x: "|".join(x))
    
    # get the AF as the mean of the pass programs
    print_if_verbose("getting allele frequency")

    df["AF"] = df.apply(lambda r: "%.4f"%convert_NaN_to0(np.mean([r["%s_fractionReadsCoveringThisVariant"%p] for p in r["PASS_programs"]])), axis=1)
    
    # get the AF for each program
    df["AF_programs"] = df.apply(lambda r: ["%s_AF=%.4f"%(p, r["%s_fractionReadsCoveringThisVariant"%p]) for p in programs], axis=1)
            
    # define the vcffields
    vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    
    # if it is haploid, avoid any issues with the genotyping
    if ploidy==1:
        
        # define the vcf fields
        df["ID"] = df["#Uploaded_variation"]
        df["QUAL"] = "."
        df["FORMAT"] = "GT:AF" 
        
        # add genotyping fields
        df["GT"] = "."
        
        # add to the info the agorithms that PASS the str
        df["AF_programs_str"] = df.AF_programs.apply(lambda x: ";".join(x))
        df["INFO"] = ("PASSALGS=" + df.PASS_programs_str + ";" + df.AF_programs_str).apply(get_corrected_INFO)
        
        # add the sample
        df["SAMPLE"] = df.GT + ":" + df.AF
        
        # the df_vcf is equivalent to the df
        df_vcf = df[vcf_fields].sort_values(by=["#CHROM", "POS", "REF"])
        
    else:
        
        # add the number of chromosomes with this variant according to each genotype
        for p in programs: 
            
            # change vars
            df["%s_GT"%p] = df["%s_GT"%p].apply(lambda x: "/".join(re.split("/|\|", x)))
            df["%s_GTset"%p] = df["%s_GT"%p].apply(lambda x: set(x.split("/")))
            df["%s_GT_index"%p] =  df["%s_GT_index"%p].apply(str)
                        
            # add n chromosomes
            df["%s_numberChromosomes_withVar"%p] = df.apply(lambda r: get_nChroms_with_var(r, p), axis=1)
                    
            # test        
            if sum(df["%s_numberChromosomes_withVar"%p]>0)!=sum(df["%s_called"%p]): 
            
                # check that there are no rows with GT and not called
                df_notCalled_withGT = df[(df["%s_numberChromosomes_withVar"%p]>0) & ~(df["%s_called"%p])] 
                if len(df_notCalled_withGT)>0: raise ValueError("There are some uncalled vars with GT")
                
        # add the numberChromosomes_withVar considering only PASS vars
        df["numberChromosomes_withVar"] = df.apply(get_numberChromosomes_withVar, axis=1) 
       
        # report if there is any POS 
        nvars_consideringREF =  len(df.drop_duplicates(subset=["#CHROM", "POS", "REF"]))
        nvars_not_consideringREF =  len(df.drop_duplicates(subset=["#CHROM", "POS"]))
        if nvars_consideringREF!=nvars_not_consideringREF: print_if_verbose("Warning there are some positions with >1 REF")                                                        
                
        # get the grouped df per chromosome and position  
        print_if_verbose("getting vcf lines")
        df_vcf = df.groupby(["#CHROM", "POS", "REF"], as_index=False).apply(lambda df_pos: get_row_vcf(df_pos, ploidy))
        df_vcf.index = list(range(len(df_vcf)))
                
    # get the vcf content
    vcf_lines = df_vcf[vcf_fields].sort_values(by=["#CHROM", "POS", "REF"]).to_csv(sep="\t", header=True, index=False)

    # get the header
    header_lines = ["##fileformat=VCFv4.2",
                    "##perSVade small variant calling pipeline. This is the merged output of freebayes (fb), GATK Haplotype Caller (HC) and bcftools (bt) for variants that PASS the filters in at least %i algorithms."%min(df.number_PASS_programs),
                    "##FILTER indicates the number of algorithms were this variant was called and PASSed the filters",
                    "##FORMAT includes the GT (genotype) and AF (allele frequency).",
                    "##GT includes the genotype in the case that all the PASS algorithms called the same GT, and the one that implies least varying positions otherwise.",
                    "##AF includes the mean fraction of reads calling this variant across PASS alorithms",
                    "##INFO includes the name of the algorithms that called this variant (PASSALGS) and the AF of each of the programs. Note that for multiallelic positions the ',' indicates each of the alleles in the order of 'ALT'"
                    ]
    
    print_if_verbose("writing %s"%(filename))
    filename_tmp = "%s.tmp"%filename
    open(filename_tmp, "w").write("\n".join(header_lines) + "\n" + vcf_lines)
    os.rename(filename_tmp, filename)

def get_PASS_vcf(vcf, replace=False):

    """This function takes a vcf and writes and returns the PASS one"""

    # define
    pass_vcf = "%s.PASS.vcf"%vcf
    pass_vcf_tmp = "%s.tmp"%pass_vcf

    if file_is_empty(pass_vcf) or replace is True:

        # get the vcf lines with PASS
        linesPASS = "%s.PASSvcflines"%vcf
        linesPASS_stderr = "%s.generating.stderr"%linesPASS
        print_if_verbose("getting the PASS lines. The stderr is in %s"%linesPASS_stderr)
        run_cmd("egrep '\tPASS\t' %s | egrep -v '^#' > %s 2>%s"%(vcf, linesPASS, linesPASS_stderr))

        # get the header
        header = "%s.header"%vcf
        header_stderr = "%s.generating.stderr"%header
        print_if_verbose("getting header. The stderr is in %s"%header_stderr)
        run_cmd("egrep '^#' %s > %s 2>%s"%(vcf, header, header_stderr))
        
        # merge into the tmp
        merging_stderr = "%s.generating.stderr"%pass_vcf_tmp
        print_if_verbose("merging. The stderr is in %s"%merging_stderr)
        run_cmd("cat %s %s > %s 2>%s"%(header, linesPASS, pass_vcf_tmp, merging_stderr))

        # remove intermediate files
        for f in [linesPASS, header, linesPASS_stderr, header_stderr, merging_stderr]: remove_file(f)

        # keep
        os.rename(pass_vcf_tmp, pass_vcf)

    return pass_vcf


def get_gzipped_file(file, replace=False):

    """Takes a file and returns the gzipped version"""

    gz_file = "%s.gz"%file
    gz_file_tmp = "%s.gz_tmp"%file

    if file_is_empty(gz_file) or replace is True:

        # run pigz with the 'gz_vcf' suffix
        pigz_std = "%s.generating.std"%gz_file_tmp
        print_if_verbose("running pigz. The std is in %s"%pigz_std)
        run_cmd("%s --keep --suffix .gz_tmp %s > %s 2>&1"%(pigz, file, pigz_std))
        remove_file(pigz_std)

        # rename 
        os.rename(gz_file_tmp, gz_file)

    return gz_file




def get_vcf_as_df_simple_oneSample(vcf_file):

    """Takes a vcf file and returns as df"""

    # get the df (avoid NA as a default NaN)
    df = pd.read_csv(vcf_file, skiprows=list(range(len([line for line in open(vcf_file, "r", encoding='utf-8', errors='ignore') if line.startswith("##")]))), sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)

    # set the index to be a tuple of (chromosome, location, ref, alt)
    df["CHROM_POS_REF_ALT"] = [tuple(x) for x in df[["#CHROM", "POS", "REF", "ALT"]].values]; df = df.set_index("CHROM_POS_REF_ALT")

    # return an empty df
    if len(df)==0: return pd.DataFrame()

    # add a colum that will result from the merging of FORMAT and the last column (which are the values of FORMAT)
    data_colname = list(df.keys())[-1]
    df["METADATA"] = [dict(zip(x[0].split(":"), x[1].split(":"))) for x in df[["FORMAT", data_colname]].values]
    features = df.iloc[0].METADATA.keys()

    # add as columns all the fetaures
    for feature in features: 

        # go through each data record
        data = []
        for rec in df.METADATA:

            if feature in rec: data.append(rec[feature])
            else: data.append("")
        df[feature] = data

    return df

def get_df_and_header_from_vcf(vcf_file):

    """Takes a vcf file and returns the df and the header lines (as a list)"""

    # define the header
    header_lines = [line.strip() for line in open(vcf_file, "r", encoding='utf-8', errors='ignore') if line.startswith("##")]

    # get the vcf
    df = pd.read_csv(vcf_file, skiprows=len(header_lines), sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)

    return df, header_lines

def get_program_that_called_vcf(vcf):

    """Takes a vcf filename and returns the program that called it"""

    all_programs = []
    for p in ["freebayes", "bcftools", "HaplotypeCaller"]:

        if p in vcf: all_programs.append(p)

    if len(all_programs)!=1: raise ValueError("The sample could not be identified")
    else: return all_programs[0]

def get_GTto0(x):

    """This function returns 0 if GT is provided, and 1 otherwise"""

    if x=="GT": return 0
    else: return 1


def get_consensus_GT_row_multialleles_diploid(r, var_to_GTaf):

    """This function takes a row of a multiallelic vcf gt and a dictionary where each variant (in the ID field) is mappend to a predicted AF from the GT of the monoallelic var"""

    # map each GT to a var
    GTaf_to_GT = {1.0:'1/1', 0.5:'0/1', 0.0:'0/0'}

    # get the split vars
    all_alts = r["ID"].split(";")
    if len(all_alts)==0: raise ValueError("The all_alts can't be empty")

    # if it is monoallelic, return the GT corresponding to the monoallelic GT
    elif len(all_alts)==1: return GTaf_to_GT[var_to_GTaf[all_alts[0]]]

    # if there 2 alleles, find a consensus
    elif len(all_alts)==2: 

        # define the vars
        alt1 = all_alts[0]
        alt2 = all_alts[1]

        # define the GTs
        gt_AF1 = var_to_GTaf[alt1]
        gt_AF2 = var_to_GTaf[alt2]

        # if both are theerozygous
        if gt_AF1==0.5 and gt_AF2==0.5: return "1/2"

        # if alt1 is homozygous
        elif gt_AF1==1.0 and gt_AF2==0.0: return "1/1"

        # if alt2 is homozygous:
        elif gt_AF1==0.0 and gt_AF2==1.0: return "2/2"

        # if alt1 is heterozygous
        elif gt_AF1==0.5 and gt_AF2==0.0: return "0/1"

        # if alt2 is homozygous:
        elif gt_AF1==0.0 and gt_AF2==0.5: return "0/2"

        # if none are called
        elif gt_AF1==0.0 and gt_AF2==0.0: return "0/0"

        # if both are called as homozygous
        elif (gt_AF1==1.0 and gt_AF2==1.0) or (gt_AF1==1.0 and gt_AF2==0.5) or (gt_AF1==0.5 and gt_AF2==1.0): return "."

        else: raise ValueError("the provided gt are not valid") 

    # if there are more than 2 alts, the GT can't be resolved
    else: return "."





def get_vcf_with_joined_multialleles_diploid(input_vcf, output_vcf, reference_genome, replace=False, threads=4):

    """Takes a vcf and joins the multialleles"""

    # define the tmp vcf
    output_vcf_tmp = "%s.tmp"%output_vcf

    if file_is_empty(output_vcf) or replace is True:
        print_if_verbose("joining vcf records in %s"%input_vcf)

        # get the processed input
        input_vcf_only_knownGT = "%s.only_knownGT.vcf"%input_vcf
        generating_knownGT_stderr = "%s.generating.stderr"%input_vcf_only_knownGT
        print_if_verbose("getting only known GT. The stderr is in %s"%generating_knownGT_stderr)
        run_cmd("egrep  '(\tknown_GT\t)|(^#)' %s > %s 2>%s"%(input_vcf, input_vcf_only_knownGT, generating_knownGT_stderr))

        # run the joining
        joining_std = "%s.joining.std"%input_vcf_only_knownGT

        print_if_verbose("running bcftools latest. The std is in %s"%joining_std)
        run_bcftools_latest("norm --check-ref ws --fasta-ref %s --multiallelics +any -o %s --output-type v --threads %i %s > %s 2>&1"%(reference_genome, output_vcf_tmp, threads, input_vcf_only_knownGT, joining_std))
        
        # check that none were changed
        if any([set(l.split()[2].split("/")[1:])!={"0"} for l in open(joining_std, "r").readlines() if "total/split/realigned/skipped" in l]): raise ValueError("some variants changed the format in bcftools norm, which is likely a bug")

        ####### EDITING GT #######
        print_if_verbose("finding consensus genotype")

        # load into vcf
        noMultiAlt_vcf = get_vcf_as_df_simple_oneSample(input_vcf_only_knownGT).set_index("ID", drop=False)
        vcf_df = get_vcf_as_df_simple_oneSample(output_vcf_tmp)

        if len(noMultiAlt_vcf)!=len(set(noMultiAlt_vcf.index)): raise ValueError("There are duplicate records in the split alleles")

        # map each var to the fraction of vars that correspond to GT
        GT_to_fractionReads = {'1/1':1.0, '0/1':0.5, '0/0':0.0}
        var_to_GTaf = dict(noMultiAlt_vcf.GT.apply(lambda x: GT_to_fractionReads[x]))

        # get new GT
        vcf_df["GT"] = vcf_df.apply(lambda r: get_consensus_GT_row_multialleles_diploid(r, var_to_GTaf), axis=1)

        # debug
        df_noGT = vcf_df[vcf_df["GT"]=="."]
        if len(df_noGT)>0: print_if_verbose("WARNING: There are %i/%i loci that can't have an assigned GT. These will be skipped"%(len(df_noGT), len(vcf_df)))

        # get only vcf with known GT
        vcf_df = vcf_df[vcf_df["GT"]!="."]

        #########################

        #### DEBUGGING DATASETS ####

        # check that the inputs are good
        if any([x not in GT_to_fractionReads for x in set(noMultiAlt_vcf.GT)]): raise ValueError("The monoallelic vcf was not properly generated")

        # check that the ploidies are consistent with the allele freq
        print_if_verbose("getting distribution of AF across several ploidies")
        for GT in GT_to_fractionReads:
            if GT in set(noMultiAlt_vcf.GT):

                # get the AF distribution
                df_GT = noMultiAlt_vcf[noMultiAlt_vcf.GT==GT]
                AF = df_GT["AF"].apply(float)

                print_if_verbose("This is the distribution of AF for GT=%s: range=(%.4f, %.4f), median=%.4f, pct_1=%.4f, pct_99=%.4f"%(GT, min(AF), max(AF), np.median(AF), np.percentile(AF, 1), np.percentile(AF, 99)))

        ############################

        # get the 
        vcf_fields = ["#CHROM",  "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

        # redefne the sample
        vcf_df["SAMPLE"] = vcf_df.GT + ":" + vcf_df.AF + ":" + vcf_df.DP + ":" + vcf_df.AD

        # write vcf lines
        vcf_lines_file = "%s.vcf_files.txt"%output_vcf_tmp
        vcf_df[vcf_fields].to_csv(vcf_lines_file, sep="\t", header=False, index=False)

        # define the header lines
        header_lines_file = "%s.header_lines.txt"%output_vcf_tmp
        header_lines_file_stderr = "%s.generating.stderr"%header_lines_file
        print_if_verbose("generating header. The stderr is in %s"%header_lines_file_stderr)
        run_cmd("egrep '^#' %s > %s 2>%s"%(output_vcf_tmp, header_lines_file, header_lines_file_stderr))

        # get the final vcf
        merging_files_stderr = "%s.generating.stderr"%output_vcf_tmp
        print_if_verbose("merging files. The stderr is in %s"%merging_files_stderr)

        run_cmd("cat %s %s > %s 2>%s"%(header_lines_file, vcf_lines_file, output_vcf_tmp, merging_files_stderr))

        # remove packages
        for f in [input_vcf_only_knownGT, joining_std, vcf_lines_file, header_lines_file, generating_knownGT_stderr, header_lines_file_stderr, merging_files_stderr]: remove_file(f)

        # rename
        os.rename(output_vcf_tmp, output_vcf)

def get_normed_bgzip_and_tabix_vcf_file(file, reference_genome, replace=False, threads=4, multiallelics_cmd="-any"):

    """This function takes a file and gzips it, creating a tabix index file.
    multiallelics_cmd is passed to bcftools norm --multiallelics. By default it joins all the multiallelics."""

    file_gz = "%s.gz"%file
    file_tmp_gz = "%s.tmp.gz"%file
    file_gz_tbi = "%s.gz.tbi"%file
    file_tmp_gz_tbi = "%s.tmp.gz.tbi"%file
    sorted_vcf = "%s.sorted.vcf"%file
    normed_vcf = "%s.norm.vcf"%file
    normed_vcf_tmp = "%s.norm.tmp.vcf"%file

    if file_is_empty(file_gz) or replace is True:

        # remove previous files
        for f in [file_tmp_gz, file_gz_tbi, file_tmp_gz_tbi, sorted_vcf]: remove_file(f)

        # sort with bedtools
        sorting_vcf_stderr = "%s.generating.stderr"%sorted_vcf
        print_if_verbose("sorting vcf. The stderr is in %s"%sorting_vcf_stderr)
        run_cmd("%s sort -header -i %s > %s 2>%s"%(bedtools, file, sorted_vcf, sorting_vcf_stderr))

        # normalise with bcftools
        normalising_vcf_std = "%s.std"%normed_vcf_tmp
        print_if_verbose("normalising vcf. STD can be found in %s"%normalising_vcf_std)
        run_bcftools_latest("norm --check-ref ws --fasta-ref %s --multiallelics %s -o %s --output-type v --threads %i %s > %s 2>&1"%( reference_genome, multiallelics_cmd, normed_vcf_tmp, threads, file, normalising_vcf_std))
        os.rename(normed_vcf_tmp, normed_vcf)

        # bgzip
        bgzip_stderr = "%s.generating.stderr"%file_tmp_gz
        print_if_verbose("bgzipping. The stderr is in %s"%bgzip_stderr)
        run_cmd("%s -c %s > %s 2>%s"%(bgzip, normed_vcf, file_tmp_gz, bgzip_stderr))

        tabix_std = "%s.tabixing.std"%file_tmp_gz
        print_if_verbose("tabix-ing. The std is in %s"%tabix_std)
        run_cmd("%s -p vcf %s > %s 2>&1"%(tabix, file_tmp_gz, tabix_std))

        # delete intermediate unnecessary files
        for f in [sorted_vcf, normed_vcf, sorting_vcf_stderr, normalising_vcf_std, bgzip_stderr, tabix_std]: remove_file(f)

        # rename files
        os.rename(file_tmp_gz_tbi, file_gz_tbi)
        os.rename(file_tmp_gz, file_gz)

    return file_gz


def get_altAllele_freq_noMultiAllele_fromAD(ad):

    """Takes AD and returns the alternative allele freq"""
    if ad==".": return 0.0
    else:  
        ad_split = ad.split(",")
        if len(ad_split)!=2: raise ValueError("AD %s is not valid"%ad)
        reads_ref = int(ad_split[0])
        reads_alt = int(ad_split[1])

        if reads_alt==0: return 0.0
        else: return reads_alt / (reads_ref + reads_alt)

def get_readsCoveringVariant(ad):

    """Takes an AD and returns the reads covering the variant"""

    # empty
    if ad==".": return 0
    else:  
        ad_split = ad.split(",")
        if len(ad_split)!=2: raise ValueError("AD %s is not valid"%ad)

        return int(ad_split[1])


def get_series_variant_in_repeats(vcf_df, repeats_table, replace=False):

    """Takes a df that has a vcf and returns a series that contains a boolean array with wether the variant intersects the repeats"""

    print_if_verbose("getting the intersection with repeats")

    # if the repeats_table is None, override
    if repeats_table is None: return [False]*len(vcf_df)

    # copy the df
    vcf_df = cp.deepcopy(vcf_df)

    # generate the repeats_positions
    repeats_positions_file = "%s.repeats_positions.py"%repeats_table
    repeats_positions_file_tmp = "%s.repeats_positions.tmp.py"%repeats_table

    if file_is_empty(repeats_positions_file) or replace is True:

        # define the bed_fields
        bed_fields = ["#chrom", "start", "end", "length"]

        # define a df with the repeats positions
        repeats_df = pd.read_csv(repeats_table, sep="\t").rename(columns={"chromosome":"#chrom", "begin_repeat":"start", "end_repeat":"end"})

        repeats_df["length"] = repeats_df.end - repeats_df.start
        repeats_df["pos_range"] = repeats_df[bed_fields].apply(lambda r: range(r["start"], r["end"]+1), axis=1)
        repeats_df["chr_range"] = repeats_df[bed_fields].apply(lambda r: [r["#chrom"]]*(r["length"]+1), axis=1)

        repeats_positions_series = pd.Series(make_flat_listOflists(repeats_df.pos_range))
        chromosomes_series = pd.Series(make_flat_listOflists(repeats_df.chr_range))

        # check that the length is the same
        if len(repeats_positions_series)!=len(chromosomes_series): raise ValueError("repeats and chromosomes are expected to be the same. This is a bug.")

        # check that the length is the same as expected_length_positions
        expected_length_positions = sum(repeats_df.length + 1)
        if len(repeats_positions_series)>expected_length_positions: raise ValueError("The length of repeats_positions_series can't exceed expected_length_positions")

        if len(repeats_positions_series)!=expected_length_positions: print_if_verbose("WARNING: There are some overlapping repeats in the genome")

        repeats_positions = set(chromosomes_series + "_" + repeats_positions_series.apply(str))
        print_if_verbose("There are %i bp with repeats"%(len(repeats_positions)))

        # save
        save_object(repeats_positions, repeats_positions_file_tmp)
        os.rename(repeats_positions_file_tmp, repeats_positions_file)

    else: repeats_positions = load_object(repeats_positions_file)

    # define a series in vcf_df that has the variant as string
    vcf_df["position_as_str"] = vcf_df["#CHROM"] + "_" + vcf_df["POS"].apply(str)

    # get the overlaps
    vcf_df["overlaps_repeats"] = vcf_df.position_as_str.isin(repeats_positions)

    print_if_verbose("There are %i/%i variants overlapping repeats"%(sum(vcf_df["overlaps_repeats"]), len(vcf_df)))

    return vcf_df["overlaps_repeats"]

def merge_several_vcfsSameSample_into_oneMultiSample_vcf(vcf_iterable, reference_genome, outdir, ploidy,  replace=False, threads=4, repeats_table=None):

    """This function takes an iterable of vcf files and gets the merged output. It writes a vcf into outdir. It only considers PASS vars"""

    # map each vcf to it's program
    program_to_vcf = {get_program_that_called_vcf(vcf) : vcf for vcf in vcf_iterable}

    # get the vcfs into a df
    program_to_vcf_df = {p : get_vcf_as_df_simple_oneSample(vcf) for p,vcf in program_to_vcf.items()}
    program_to_vcf_df = {p : df for p, df in program_to_vcf_df.items() if len(df)>0}

    # define the common 'FORMAT' fields
    common_format_fields = sorted(set.intersection(*[set(df.FORMAT.iloc[0].split(":")) for df in program_to_vcf_df.values()]), key=get_GTto0)
    print_if_verbose("These are the common FORMAT fields:", common_format_fields)
    if len(common_format_fields)==0: raise ValueError("There are no common FORMAT fields")

    # get the sampleID
    sampleIDs = {df.columns[9] for df in program_to_vcf_df.values()}
    if len(sampleIDs)!=1: raise ValueError("You are not trying to merge vcfs from the same sample")
    sampleID = next(iter(sampleIDs))

    # define the vcf fields (it is missing the sample)
    backbone_vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    # map each caller to an abbrebiations
    program_to_abbreviation = {"HaplotypeCaller":"HC", "freebayes":"fb", "bcftools":"bt"}
    
    # go through the types of filters
    for type_filters in ["all"]:
        print_if_verbose(type_filters)

        # deepcopy the df
        p_to_df = cp.deepcopy(program_to_vcf_df)

        # define the outfile
        merged_vcf = "%s/merged_vcfs_%sVars_ploidy%i.vcf"%(outdir, type_filters, ploidy)
        merged_vcf_tmp = "%s/merged_vcfs_%sVars_ploidy%i.tmp.vcf"%(outdir, type_filters, ploidy)

        if file_is_empty(merged_vcf) or replace is True:

            # initialize a list of the important vcfs
            vcfs_to_merge = []

            # go through each vcf and format the fields
            for program, vcf_df in p_to_df.items():
                print_if_verbose(program)

                # define the formatted vcf
                formatted_vcf = "%s.formatted.%sVars.vcf"%(program_to_vcf[program], type_filters)
                formatted_vcf_gz = "%s.gz"%formatted_vcf

                if file_is_empty(formatted_vcf_gz) or replace is True:

                    print_if_verbose("formatting vcf")

                    ########## FORMAT VCF ##########

                    # keep only the PASS variants if necessary
                    if type_filters=="onlyPASS": vcf_df = vcf_df[vcf_df.FILTER=="PASS"]

                    # format the FORMAT to include only 
                    vcf_df["FORMAT"] = ":".join(common_format_fields)
                    vcf_df[program] = vcf_df.apply(lambda r: ":".join([r[f] for f in common_format_fields]), axis=1)

                    # make upper vcfs
                    vcf_df["REF"]  = vcf_df["REF"].apply(lambda x: x.upper())
                    vcf_df["ALT"]  = vcf_df["ALT"].apply(lambda x: x.upper())

                    # format the INFO to include the INFO, FILTER, QUAL and common_format_fields
                    abb = program_to_abbreviation[program]
                    vcf_df["INFO_with_abb"] = vcf_df.INFO.apply(lambda x: ";".join(["%s_%s"%(abb, I) for I in x.split(";")]))
                    vcf_df["FILTER_and_QUAL_with_abb"] = vcf_df.apply(lambda r: "%s_FILTER=%s;%s_QUAL=%.2f"%(abb, r["FILTER"], abb, r["QUAL"]) , axis=1)
                    vcf_df["INFO"] = vcf_df.FILTER_and_QUAL_with_abb + ";" + vcf_df.INFO_with_abb + ";%s_DATA="%abb + vcf_df[program]

                    #################################

                    # write the vcf lines
                    vcf_lines = vcf_df[backbone_vcf_fields + [program]].to_csv(sep="\t", header=True, index=False)
                    
                    # get the header lines 
                    header_lines = [l.strip() for l in open(program_to_vcf[program], "r").readlines() if l.startswith("##")]

                    ##### ADD HEADER LINES #####

                    # map each INFO id to the description
                    infoID_to_Header = {l.split("ID=")[1].split(",")[0] :l for l in header_lines if l.startswith("##INFO=<")}

                    # map each FORMAT id to the description
                    formatID_to_Header = {l.split("ID=")[1].split(",")[0] :l for l in header_lines if l.startswith("##FORMAT=<")}

                    # initialize without filter and info
                    #edited_header_lines = [l for l in header_lines if not any([l.startswith("##%s="%x) for x in ["INFO", "FORMAT"]])]
                    edited_header_lines = [l for l in header_lines if not any([l.startswith("##%s="%x) for x in ["INFO", "FORMAT", "ALT"]])]

                    # INFO

                    # from info
                    edited_header_lines += [l.replace("ID=%s"%infoID, "ID=%s_%s"%(abb, infoID)) for infoID, l in infoID_to_Header.items()]

                    # add filter, qual and data
                    edited_header_lines += ['##INFO=<ID=%s_FILTER,Number=1,Type=String,Description="The FILTER field by %s">'%(abb, program)] # FILTER
                    edited_header_lines += ['##INFO=<ID=%s_QUAL,Number=1,Type=Float,Description="The QUAL field by %s">'%(abb, program)] # QUAL
                    edited_header_lines += ['##INFO=<ID=%s_DATA,Number=.,Type=String,Description="The DATA field by %s">'%(abb, program)] # DATA

                    # FORMAT

                    # only keep format lines for the common ones
                    edited_header_lines += [formatID_to_Header[formatID] for formatID in common_format_fields]

                    # write vcf
                    open(formatted_vcf, "w").write("\n".join(edited_header_lines) + "\n" + vcf_lines)
                    
                    ####################

                    # bgzip, tabix and split mulrialleles into sepparate records
                    get_normed_bgzip_and_tabix_vcf_file(formatted_vcf, reference_genome, replace=True, threads=threads)

                # keep
                vcfs_to_merge.append(formatted_vcf_gz)

            # define the bcftools_merge_std
            bcftools_merge_std = "%s.generating.std"%merged_vcf_tmp
            print_if_verbose("generating merged vcf. The std is in %s"%bcftools_merge_std)

            # run bcftools merge only if there are more than 1 vcf
            if len(vcfs_to_merge)>1:

                run_cmd("%s merge --merge none -o %s -Ov --threads %i %s > %s 2>&1"%(bcftools, merged_vcf_tmp, threads, " ".join(vcfs_to_merge), bcftools_merge_std))

            elif len(vcfs_to_merge)==1: 

                run_cmd("%s view -o %s -Ov --threads %i %s > %s 2>&1"%(bcftools, merged_vcf_tmp, threads, vcfs_to_merge[0], bcftools_merge_std))

            else: raise ValueError("there are no vcfs to merge")   

                
            ######## ADD EXTRA FILEDS TO INFO ######## 
            print_if_verbose("editing INFO")

            # load into df and add the number of PASS vars and also the PASS programs
            header_lines = [line.strip() for line in open(merged_vcf_tmp, "r", encoding='utf-8', errors='ignore') if line.startswith("##")]
            vcf_df = pd.read_csv(merged_vcf_tmp, skiprows=list(range(len(header_lines))), sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)

            # remove the duplicates
            vcf_df = vcf_df.drop_duplicates(subset=["#CHROM", "POS", "REF", "ALT"])


            # replace by a '.' if empty
            def get_point_if_empty(x):
                if len(x)>0: return x
                else: return ["."]


            # add the called and PASS programs
            fields = list(vcf_df.columns)
            all_programs = list(p_to_df)
            vcf_df["called_programs_list"] = vcf_df[all_programs].apply(lambda r: [program_to_abbreviation[p] for p in all_programs if not r[p].endswith(".")], axis=1)
            vcf_df["PASS_programs_list"] = vcf_df.INFO.apply(lambda x: [program_to_abbreviation[p] for p in all_programs if "%s_FILTER=PASS"%program_to_abbreviation[p] in x]).apply(get_point_if_empty)

            # get whether the variant overlaps repeats
            vcf_df["overlaps_repeats"] = get_series_variant_in_repeats(vcf_df, repeats_table, replace=replace)

            # define the 'ID' in  a way that resembles VEP
            vcf_df["ID"] = vcf_df["#CHROM"] + "_" + vcf_df.POS.apply(str) + "_" + vcf_df.REF + "/" + vcf_df.ALT

            # check that the ID is unique (check that the drop_duplicates worked)
            if len(set(vcf_df.ID))!=len(vcf_df): 

                #duplicated_ID = vcf_df[vcf_df.duplicate]
                #print(.duplicate())

                raise ValueError("The ID has to be unique")

            # check if there are empty programs
            for f in ["called_programs_list", "PASS_programs_list"]: 
                if any(vcf_df[f].apply(len)==0): raise ValueError("There are empty programs")

            # debug multialleles
            if any(vcf_df.ALT.apply(lambda x: "," in x)): raise ValueError("There are multiallelic records")

            # define SNPs
            strange_nucs = {"*", "-"}
            vcf_df["is_snp"] = (vcf_df.REF.apply(len)==1) &  (vcf_df.ALT.apply(len)==1) & ~(vcf_df.REF.isin(strange_nucs)==1) & ~(vcf_df.ALT.isin(strange_nucs)==1)

            # remove the FILTER
            vcf_df["FILTER"] = "."

            # add to info
            info_series = vcf_df.INFO + ";CALLEDALGS=" + vcf_df.called_programs_list.apply(lambda x: "-".join(x)) + ";PASSALGS=" + vcf_df.PASS_programs_list.apply(lambda x: "-".join(x)) + ";NCALLED=" + vcf_df.called_programs_list.apply(len).apply(str) + ";NPASS=" + vcf_df.PASS_programs_list.apply(lambda x: len([y for y in x if y!="."])).apply(str) + ";ISSNP=" + vcf_df.is_snp.apply(str)


            # initialize description and type lines for INFO fields
            f_to_description = {"CALLEDALGS":"The algorithms that called this var, sepparated by '-'",
                                "PASSALGS":"The algorithms where this var PASSed the filters, sepparated by '-'",
                                "NCALLED":"The number of algorithms that called this var",
                                "NPASS":"The number of algorithms where this var PASSed the filters",
                                "ISSNP":"Whether it is a SNP"}

            f_to_type = {"CALLEDALGS":"String",
                         "PASSALGS":"String",
                         "NCALLED":"Integer",
                         "NPASS":"Integer",
                         "ISSNP":"String"}

            # add the sets of algorithms
            vcf_df["called_algs_set"] = vcf_df.called_programs_list.apply(lambda x: set(x).difference({"."}))
            vcf_df["PASS_algs_set"] = vcf_df.PASS_programs_list.apply(lambda x: set(x).difference({"."}))

            # define the interesting_algs as those that are either PASS or called if none are PASS
            def get_interesting_algs(r):

                if len(r["PASS_algs_set"])>0: interesting_algs = r["PASS_algs_set"]
                else : interesting_algs = r["called_algs_set"]

                return interesting_algs

            vcf_df["interesting_algs"] = vcf_df.apply(get_interesting_algs, axis=1)

            # add some fields of each program
            print_if_verbose("adding program specific info")
            ADidx = [I for I, field in enumerate(vcf_df.FORMAT.iloc[0].split(":")) if field=="AD"][0]
            DPidx = [I for I, field in enumerate(vcf_df.FORMAT.iloc[0].split(":")) if field=="DP"][0]
            GTidx = [I for I, field in enumerate(vcf_df.FORMAT.iloc[0].split(":")) if field=="GT"][0]

            for p in all_programs:
                abb = program_to_abbreviation[p]

                # get the reads covering this variant
                vcf_df["%s_readsCovVar"%abb] = vcf_df[p].apply(lambda x: x.split(":")[ADidx]).apply(get_readsCoveringVariant)

                # get the AD
                vcf_df["%s_AD"%abb] = vcf_df[p].apply(lambda x: x.split(":")[ADidx])

                # get the total depth at the locus
                vcf_df["%s_DP"%abb] = vcf_df[p].apply(lambda x: x.split(":")[DPidx]).apply(get_int_or0)

                # add the AF by the program
                vcf_df["%s_AF"%abb] = (vcf_df["%s_readsCovVar"%abb] / vcf_df["%s_DP"%abb]).apply(getNaN_to_0)

                # get the genotype
                vcf_df["%s_GT"%abb] = vcf_df[p].apply(lambda x: x.split(":")[GTidx].replace("|", "/")) 

                # check that the GTs only include 0s and 1s, so that tehre are no multiallelic records
                if any(vcf_df["%s_GT"%abb].apply(lambda x: any([gt not in {"0", "1", "."} for gt in x.split("/") ]))):
                    raise ValueError("There may be some multiallelic records")

                # add the reorderedGT
                vcf_df["%s_GTreordered"%abb] = vcf_df["%s_GT"%abb].apply(lambda x: "/".join(sorted(x.split("/"))) )

                # add the PASS and called algorithsm
                vcf_df["%s_called"%abb] = vcf_df.called_algs_set.apply(lambda x: abb in x)
                vcf_df["%s_PASS"%abb] = vcf_df.PASS_algs_set.apply(lambda x: abb in x)
                
                # add to info
                info_series += ";%s_fractionReadsCov="%abb + vcf_df["%s_AF"%abb].apply(lambda x: "%.4f"%x) + ";%s_GT="%abb + vcf_df["%s_GT"%abb] + ";%s_called="%abb + vcf_df["%s_called"%abb].apply(str) + ";%s_PASS="%abb + vcf_df["%s_PASS"%abb].apply(str) + ";%s_readsCovVar="%abb + vcf_df["%s_readsCovVar"%abb].apply(str)

                # add the headers
                f_to_description["%s_fractionReadsCov"%abb] = "The fraction of reads covering this var by %s"%p
                f_to_type["%s_fractionReadsCov"%abb] = "Float"

                f_to_description["%s_readsCovVar"%abb] = "The number of reads covering this var by %s"%p
                f_to_type["%s_readsCovVar"%abb] = "Integer"

                f_to_description["%s_GT"%abb] = "The GT by %s"%p
                f_to_type["%s_GT"%abb] = "String"

                f_to_description["%s_called"%abb] = "Whether the variant was called by %s"%p
                f_to_description["%s_PASS"%abb] = "Whether the variant PASSed the filters by %s"%p

                f_to_type["%s_called"%abb] = "String"
                f_to_type["%s_PASS"%abb] = "String"


            # get the mean AD
            def get_mean_AD(list_ADs):

                # get the mean for each of the ref and alt
                mean_ref_reads = int(np.mean([int(AD.split(",")[0]) for AD in list_ADs]))
                mean_alr_reads = int(np.mean([int(AD.split(",")[1]) for AD in list_ADs]))

                return "%i,%i"%(mean_ref_reads, mean_alr_reads)

            vcf_df["mean_AD"] = vcf_df.apply(lambda r: get_mean_AD([r["%s_AD"%p] for p in r["interesting_algs"]]), axis=1)
            f_to_description["mean_AD"] = "The mean AD across  all programs that are PASS or all called programs if none are PASS"
            f_to_type["mean_AD"] = "String"
            info_series += ";mean_AD=" + vcf_df["mean_AD"]

            print_if_verbose("getting common genotype")
            # get the common GT
            def get_commonGT(r):

                # get all the GTs
                all_GTs = [r["%s_GTreordered"%p] for p in r["interesting_algs"]]
                all_GTs_set = set(all_GTs)

                if len(r["called_algs_set"])>3: raise ValueError("This function does not work for >3 programs")

                # if there is only one GT, return it
                if len(all_GTs_set)==1: commonGT = next(iter(all_GTs_set))

                # if there are 3 programs calling, you may have 2vs1 GT. If so, keep the most common ones
                elif len(r["called_algs_set"])==3:

                    # define the called GTs
                    all_GTs_called = [r["%s_GTreordered"%p] for p in r["called_algs_set"]]
                    all_GTs_called_set = set(all_GTs_called)
                    
                    # if each program calls a different GT, define as 'no-consensus'
                    if len(all_GTs_called_set)==3: commonGT = "."

                    # if not, you can the GT that is most commonly-called
                    else:

                        # map each GT to the number of programs that call it
                        nPrograms_to_GT = {nPrograms : GT for GT, nPrograms in Counter(all_GTs_called).items()}
                        if len(nPrograms_to_GT)<=1: raise ValueError("something went wrong with the parsing")
                        commonGT = nPrograms_to_GT[max(nPrograms_to_GT)]

                else: commonGT = "." 

                return commonGT

            vcf_df["common_GT"] = vcf_df.apply(get_commonGT, axis=1)

            # keep
            f_to_description["common_GT"] = "The GT if it is common by all the PASS algorithms (or called if there are none). If there is no agreement between these algorithms it is '.'"
            f_to_type["common_GT"] = "String"
            info_series += ";common_GT=" + vcf_df["common_GT"]

            # add whether it overlaps repeats
            if repeats_table is not None: 

                f_to_description["INREPEATS"] = "A boolean stating whether the variants overlap any repeat as annotated with RepeatModeler and RepeatMasker"
                f_to_type["INREPEATS"] = "String"
                info_series += ";INREPEATS=" + vcf_df["overlaps_repeats"].apply(str)

            # set the QUAL to be the mean by the interesting_algs
            print_if_verbose("getting QUAL")
            vcf_df["QUAL"] = vcf_df.apply(lambda r: np.mean([float(r["INFO"].split("%s_QUAL="%p)[1].split(";")[0]) for p in r["interesting_algs"]]), axis=1)

            # get the mean DP 
            print_if_verbose("getting mean DP")
            vcf_df["mean_DP"] = vcf_df.apply(lambda r: np.mean([r["%s_DP"%p] for p in r["interesting_algs"]]), axis=1)
            
            f_to_description["mean_DP"] = "The mean read depth by all programs that are PASS or all called programs if none are PASS"
            f_to_type["mean_DP"] = "Float"
            info_series += ";mean_DP=" + vcf_df["mean_DP"].apply(lambda x: "%.4f"%x)

            # add the mean AF for PASS and CALLED programs
            for targetField in ["called_algs_set", "PASS_algs_set"]: 

                # define the allele frequency  for the given set of algs
                type_algs = targetField.split("_")[0]
                new_field = "mean_fractionReadsCov_%s_algs"%type_algs
                vcf_df[new_field] = vcf_df.apply(lambda r: np.mean([r["%s_AF"%alg] for alg in r[targetField] if alg!="."]), axis=1).apply(getNaN_to_0)

                # keep
                f_to_description[new_field] = "The mean fraction of reads covering this variant by the %s algorithms"%type_algs

                f_to_type[new_field] = "Float"
                info_series += ";%s="%new_field + vcf_df[new_field].apply(lambda x: "%.4f"%x)

          
            vcf_df["INFO"] = info_series

            # add the programs' description
            for f, description in f_to_description.items(): header_lines += ['##INFO=<ID=%s,Number=1,Type=%s,Description="%s">'%(f, f_to_type[f], description)]

            # test that some fields have non NaNs
            for field in ["mean_fractionReadsCov_PASS_algs", "mean_DP", "common_GT"]: 
                if any(pd.isna(vcf_df[field])): raise ValueError("There are NaNs in %s"%field)


            # remove the std
            remove_file(bcftools_merge_std)

            # write vcf
            vcf_lines = vcf_df[fields].to_csv(sep="\t", header=True, index=False)
            open(merged_vcf_tmp, "w").write("\n".join(header_lines) +  "\n" + vcf_lines)
            os.rename(merged_vcf_tmp, merged_vcf)

            ################################ 
            
        # define variables to return 
        if type_filters=="all": merged_vcf_all = merged_vcf
        elif type_filters=="onlyPASS": merged_vcf_onlyPASS = merged_vcf


    return merged_vcf_all

def get_float_from_string(x):

    """Gets a float from a string, or empty"""

    if x in {"", ".", " "}: return np.nan
    else: return float(x)

def get_int_from_string(x):

    """Gets a float from a string, or -1"""

    if x in {"", ".", " "}: return np.nan
    else: return int(x)


def get_vep_df_for_vcf_df(vcf_df, outdir, reference_genome, gff_with_biotype, mitochondrial_chromosome, mitochondrial_code, gDNA_code, replace):

    """This function takes a vcf_df and runs vep for this chunk, while keeping. vcf_df index should be unique from 0 ti nvars  """

    # define the prefix
    prefix = "%s/vep_%i_to_%i"%(outdir, vcf_df.index[0], vcf_df.index[-1])

    # define the vcf file
    vcf_file = "%s_variants.vcf"%prefix

    # define the annotated_vcf 
    annotated_vcf = "%s_annotated.tab"%vcf_file; annotated_vcf_tmp = "%s.tmp"%annotated_vcf

    if file_is_empty(annotated_vcf) or replace is True:
        print_if_verbose("running vep for %s"%vcf_file)

        # clean previous files
        for f in os.listdir(outdir):
            path = "%s/%s"%(outdir, f)
            if path.startswith(prefix) and path!=annotated_vcf: remove_file(path)

        # generate the raw vcf
        vcf_df.to_csv(vcf_file, sep="\t", index=False, header=True)

        # run vep for this vcf
        vep_std = "%s_annotating_vep.std"%prefix
        run_cmd("%s --input_vcf %s --outfile %s --ref %s --gff %s --mitochondrial_chromosome %s --mito_code %i --gDNA_code %i > %s 2>&1"%(run_vep, vcf_file, annotated_vcf_tmp, reference_genome, gff_with_biotype, mitochondrial_chromosome, mitochondrial_code, gDNA_code, vep_std))

        # check that the std contains no signs of compressing the gff
        if any(["compressing gff before running vep" in l for l in open(vep_std, "r").readlines()]): raise ValueError("There was a compression of the gff before running vep. This is not acceptable when running in parallel")   

        remove_file(vep_std)

        # keep
        os.rename(annotated_vcf_tmp, annotated_vcf)

    # remove all the files that are related to this prefix
    for f in os.listdir(outdir):
        path = "%s/%s"%(outdir, f)
        if path.startswith(prefix) and path!=annotated_vcf: remove_file(path)

    # load the vep df
    df_vep = pd.read_csv(annotated_vcf, sep="\t")

    return df_vep

def run_vep_parallel(vcf, reference_genome, gff_with_biotype, mitochondrial_chromosome, mitochondrial_code, gDNA_code, threads=4, replace=False):

    """This function runs vep in parallel and returns an annotated_vcf"""

    # define an output file for VEP
    output_vcf = "%s_annotated.tab"%vcf; output_vcf_tmp = "%s.tmp"%output_vcf

    # run annotation by VEP
    if file_is_empty(output_vcf) or replace is True:

        print_if_verbose("Annotating with VEP %s"%vcf)

        # get the input vcf as a df
        print_if_verbose("loading vcf_df")
        df_vcf, header = get_df_and_header_from_vcf(vcf)
        df_vcf.index = list(range(len(df_vcf)))

        # check that the uploaded variation is unique
        if len(set(df_vcf.ID))!=len(df_vcf): raise ValueError("The IDs are not unique")

        # get chunks of the df
        chunk_size = 2000
        chunks_vcf_df = [df_vcf.loc[chunk_idxs] for chunk_idxs in chunks(df_vcf.index, chunk_size)]
        print_if_verbose("There are %i chunks of %i variants. Running on %i threads"%(len(chunks_vcf_df), chunk_size, threads))

        # define an outdir where to write the 
        outdir_intermediate_files = "%s.chunks_vep_annotation"%vcf
        delete_folder(outdir_intermediate_files)
        make_folder(outdir_intermediate_files)

        ####### get the gff tabixed and sorted #######

        # this is a necessary step to run vep, and you don't want each individual run to repeat it

        gff_clean = "%s_clean.gff"%gff_with_biotype
        gff_clean_compressed = "%s_clean.gz"%gff_with_biotype
        gff_clean_compressed_tbi = "%s.tbi"%gff_clean_compressed

        # remove previous files
        remove_file(gff_clean)
        remove_file(gff_clean_compressed)
        remove_file(gff_clean_compressed_tbi)

        if file_is_empty(gff_clean_compressed_tbi):
            print_if_verbose("compressing gff before running vep in parallel")

            # eliminate strange lines,chromosomes and compress
            cleaning_gff_stderr = "%s.generating.stderr"%gff_clean
            print_if_verbose("cleaning gff. The stderr is in %s"%cleaning_gff_stderr)
            run_cmd("%s sort -i %s | egrep -v '^#' | egrep -v '\tchromosome\t' > %s 2>%s"%(bedtools, gff_with_biotype, gff_clean, cleaning_gff_stderr))
            
            compresing_gff_stderr = "%s.generating.stderr"%gff_clean_compressed
            print_if_verbose("compressing gff. The stderr is in %s"%compresing_gff_stderr)
            run_cmd("%s -c %s > %s 2>%s"%(bgzip, gff_clean, gff_clean_compressed, compresing_gff_stderr))

            # index with tabix
            tabixing_std = "%s.tabixing.std"%gff_clean_compressed
            print_if_verbose("tabixing. The std is in %s"%tabixing_std)
            run_cmd("%s %s > %s 2>&1"%(tabix, gff_clean_compressed, tabixing_std))

            for f in [cleaning_gff_stderr, compresing_gff_stderr, tabixing_std]: remove_file(f)

        ################################################

        # get the inputs of the function
        inputs_fn = [(c, outdir_intermediate_files, reference_genome, gff_with_biotype, mitochondrial_chromosome, mitochondrial_code, gDNA_code, replace) for c in chunks_vcf_df]

        # run in parallel
        with multiproc.Pool(threads) as pool:
            list_vep_dfs = pool.starmap(get_vep_df_for_vcf_df, inputs_fn) 
                
            pool.close()
            pool.terminate()

        # run one after the other, this is to test
        #list_vep_dfs = list(map(lambda x: get_vep_df_for_vcf_df(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]), inputs_fn))

        # concatenate
        print_if_verbose("concatenating")
        df_vep = pd.concat(list_vep_dfs).sort_values(by="#Uploaded_variation").drop_duplicates()

        # check that not all of the variants are intergenic (this would be suspicious)
        if (sum(df_vep.Consequence=="intergenic_variant")/len(df_vep)) > 0.9: print_if_verbose("WARNING: There are >90 percent of the variants that are intergenic. Maybe your gff is wrong. ")
        # write
        print_if_verbose("saving")
        df_vep.to_csv(output_vcf_tmp, sep="\t", header=True, index=False)

        # remove the intermediate files
        delete_folder(outdir_intermediate_files)
        
        # keep
        os.rename(output_vcf_tmp, output_vcf)

    return output_vcf

def get_INFOdict(info, infoField_to_typeFN):

    """Takes an INFO string of a vcf and returns a dict with the transformed types"""

    # get the info that is sepparated by '='
    infoField_to_content = {x.split("=")[0] : x.split("=")[1] for x in info.split(";") if "=" in x}

    # get the info_dict
    info_dict = {}

    for infoField, typeFN in infoField_to_typeFN.items():

        # parse the booleans in a special way
        if typeFN==bool:

            # other Bools are Flags. Keep them if they are in info
            info_dict[infoField] = (infoField in info)

        else:

            # if it is there, keep it
            if infoField in infoField_to_content: info_dict[infoField] = typeFN(infoField_to_content[infoField])

            # if not, get the missing value
            else: info_dict[infoField] = typeFN(".")

    return info_dict

def write_variantInfo_table(vcf, variantInfo_table, replace=False):

    """This function takes a vcf from merge_several_vcfsSameSample_into_oneMultiSample_vcf and writes it in tabular format into variantInfo_table"""

    # define the tmp
    variantInfo_table_tmp = "%s.tmp"%variantInfo_table

    if file_is_empty(variantInfo_table) or replace is True:
        print_if_verbose("converting %s to tab"%vcf)

        # load the vcf table
        df, header_lines = get_df_and_header_from_vcf(vcf)

        # map each info field to the type
        typeStr_to_typeFN = {"Float":get_float_from_string, "String":str, "Integer":get_int_from_string, "Flag":bool}
        infoField_to_typeFN = {h.split("<ID=")[1].split(",")[0] : typeStr_to_typeFN[h.split("Type=")[1].split(",")[0]]  for h in header_lines if h.startswith("##INFO=<ID=")}

        # change a couple of types that are incorrectly parsed
        infoField_to_typeFN["bt_DP4"] = str

        # add the INFO formated as in the df
        print_if_verbose("getting INFO_dict")
        df["INFO_dict"] = df.INFO.apply(lambda info: get_INFOdict(info, infoField_to_typeFN))
        all_INFO_tags = sorted(set.union(*df.INFO_dict.apply(set)))

        # create a field that is uploaded variation
        df["#Uploaded_variation"] = df.ID

        # initalize the important fields
        important_fields = ['#Uploaded_variation', '#CHROM', 'POS', 'REF', 'ALT', 'QUAL'] + all_INFO_tags

        # add as columns of df
        print_if_verbose("adding cols to df")
        for f in all_INFO_tags: df[f] = df.INFO_dict.apply(lambda x: x[f])

        # write the important fields
        print_if_verbose("saving to %s"%variantInfo_table)
        df[important_fields].to_csv(variantInfo_table_tmp, sep="\t", index=False, header=True)
        os.rename(variantInfo_table_tmp, variantInfo_table)


    # return as df
    print_if_verbose("returning df")
    return pd.read_csv(variantInfo_table, sep="\t")



def report_variant_calling_statistics(df, variantCallingStats_tablePrefix, programs):

    """This function takes a df with the tabular vcf and generates a report of how the variant calling went. """

    print_if_verbose("generating report of the called variants")

    # define some vars
    p_to_abb = {"HaplotypeCaller":"HC", "freebayes":"fb", "bcftools":"bt"}

    # define the program combinations
    prog_combinations = [("HaplotypeCaller",),
                         ("freebayes",),
                         ("bcftools",),
                         ("HaplotypeCaller","freebayes"),
                         ("HaplotypeCaller","bcftools"),
                         ("freebayes","bcftools"),
                         ("HaplotypeCaller","freebayes","bcftools")]

    # redefine the programs. Sometimes there are no vars called by a program
    programs = [p for p in programs if "%s_called"%(p_to_abb[p]) in df.keys()]

    # go through different types of filters
    for type_filter in ["PASS", "called"]:
        print_if_verbose(type_filter)

        # initialize dict
        stats_dict = {}

        # go through each of the SNPs and indels
        for type_var in ["SNP", "INDEL"]:

            # get the df of the type var
            if type_var=="SNP": df_type_var = df[df.ISSNP]
            elif type_var=="INDEL": df_type_var = df[~df.ISSNP]



            # get the vars of each type
            p_to_vars = {p : set(df_type_var[df_type_var["%s_%s"%(p_to_abb[p], type_filter)]]["#Uploaded_variation"]) for p in programs}

            for prog_comb in prog_combinations: 

                # get if all the programs are in programs 
                if all([x in programs for x in prog_comb]): 

                    # get as string
                    prog_comb_name = "-".join(prog_comb)

                    # define the number of vars that are in the intersection of the programs
                    nVars_in_intersection = len(set.intersection(*[p_to_vars[p] for p in prog_comb]))
                    nVars_in_Programs = len(set.union(*[p_to_vars[p] for p in prog_comb]))
                    if nVars_in_Programs==0: pctVars_in_intersection = 0
                    else: pctVars_in_intersection = (nVars_in_intersection/nVars_in_Programs)*100

                    # keep 
                    stats_dict.setdefault("%s_nVars"%type_var, {}).setdefault(prog_comb_name, nVars_in_intersection)
                    stats_dict.setdefault("%s_pctVars"%type_var, {}).setdefault(prog_comb_name, pctVars_in_intersection)


        # get df
        df_stats = pd.DataFrame(stats_dict)
        print_if_verbose("\n", df_stats)
        df_stats.to_csv("%s_%s.tab"%(variantCallingStats_tablePrefix, type_filter), sep="\t")


def run_perSVade_severalSamples(paths_df, cwd, common_args, time_greasy="48:00:00", threads=4, sampleID_to_parentIDs={}, samples_to_run=set(), repeat=False, mn_queue="bsc_ls", job_array_mode="greasy", max_ncores_queue=48, ploidy=1, variant_calling_fields=["#Uploaded_variation", "QUAL", "fb_DP", "fb_MQM", "fb_MQMR", "fb_PQA", "fb_PQR", "fb_QA", "fb_QR", "fb_fractionReadsCov", "fb_readsCovVar"]):

 
    """
    This function inputs a paths_df, which contains an index as 0-N rows and columns "reads", "sampleID", "readID"  and runs the perSVade pipeline without repeating steps (unless indicated). pths_df can also be a tab-sepparated file were there are this 3 fields. The trimmed_reads_dir has to be the full path to the .fastq file. The sampleID has to be the unique sample identifier and the readID has to be R1 or R2 for paired-end sequencing. The p

    - cwd is the current working directory, where files will be written
    - repeat is a boolean that indicates whether to repeat all the steps of this function
    - threads are the number of cores per task allocated. In mn, you can use 48 cores per node. It has not been tested whether more cores can be used per task
    - samples_to_run is a set of samples for which we want to run all the pipeline
    - time_greasy is the time greasy gets in DD:HH:MM
    - mn_queue is the queue were you want to throw the pipeline. You can use "bsc_queues" to see which you have. debug is a faster one, but you can only specify 2h 
    - job_array_mode can be 'greasy' or 'local'. If local each job will be run after the other
    - sampleID_to_parentIDs is a dictionary that maps each sampleID to the parent sampleIDs (a set), in a way that there will be a col called parentIDs_with_var, which is a string of ||-joined parent IDs where the variant is also found
    - common_args is a string with all the perSVade args except the reads. The arguments added will be -o, -f1, -f2
    - max_ncores_queue is the total number of cores that will be assigned to the job.
    - ploidy is the ploidy with which to run the varcall
    - variant_calling_fields are the fields in variant_calling_ploidy<N>.tab to keep in the concatenated data
    """

    print_if_verbose("Running VarCall pipeline...")

    # if it is a path
    if type(paths_df)==str: paths_df = pd.read_csv(paths_df, sep="\t")

    # create files that are necessary
    VarCallOutdirs = "%s/VarCallOutdirs"%cwd; make_folder(VarCallOutdirs)
    
    # define the samples_to_run
    if len(samples_to_run)==0: samples_to_run = set(paths_df.sampleID)

    # get the info of all the reads and samples
    all_cmds = []

    for sampleID in samples_to_run:

        # define the df for this sample
        df = paths_df[paths_df.sampleID==sampleID]
        df1 = df[df.readID=="R1"]
        df2 = df[df.readID=="R2"]

        # define the reads of interest and keep
        reads1 = df1.reads.values[0]
        reads2 = df2.reads.values[0]

        # create an outdir
        outdir = "%s/%s_VarCallresults"%(VarCallOutdirs, sampleID); make_folder(outdir)

        # define the files that shoud be not empty in order not to run this code
        success_files = ["%s/smallVars_CNV_output/variant_annotation_ploidy%i.tab"%(outdir, ploidy)]
                   
        # define the cmd          
        cmd = "%s -f1 %s -f2 %s -o %s --ploidy %i %s"%(perSVade_py, reads1, reads2, outdir, ploidy, common_args)

        # add cmd if necessary
        if any([file_is_empty(x) for x in success_files]) or repeat is True: all_cmds.append(cmd)

    # submit to cluster or return True
    if len(all_cmds)>0:

        if job_array_mode=="local":

            for Icmd, cmd in enumerate(all_cmds):
                print_if_verbose("running cmd %i/%i"%(Icmd+1, len(all_cmds)))
                run_cmd(cmd)

        elif job_array_mode=="greasy":

            print_if_verbose("Submitting %i jobs to cluster ..."%len(all_cmds))
            jobs_filename = "%s/jobs.run_SNPs_CNV"%cwd
            open(jobs_filename, "w").write("\n".join(all_cmds))

            generate_jobarray_file_greasy(jobs_filename, walltime=time_greasy,  name="perSVade_many_samples", queue=mn_queue, sbatch=True, ncores_per_task=threads, constraint="", number_tasks_to_run_at_once="all", max_ncores_queue=max_ncores_queue )

        else: raise ValueError("%s is not a valid job_array_mode"%job_array_mode)

        return False

    print_if_verbose("Integrating all variants and CNV into one......")

    ###### INTEGRATE VARIANT CALLING ######

    # define the file
    variant_calling_df_file = "%s/integrated_variant_calling_ploidy%i.tab"%(cwd, ploidy)

    if file_is_empty(variant_calling_df_file) or repeat is True:
        print_if_verbose("generating integrated vars")

        # define the columns related to variant_calling_fields
        df_example = pd.read_csv("%s/%s_VarCallresults/smallVars_CNV_output/variant_calling_ploidy%i.tab"%(VarCallOutdirs, next(iter(samples_to_run)), ploidy), sep="\t")
        variant_calling_colNames = ",".join([str(I+1) for I, field in enumerate(df_example.keys()) if field in variant_calling_fields])

        del df_example

        # initialize df
        df_variant_calling = pd.DataFrame()


        for Is, sampleID in enumerate(samples_to_run):
            print_if_verbose("%i/%i: %s"%(Is+1, len(samples_to_run), sampleID))

            # get the partial file
            target_varcall_file = "%s/%s_VarCallresults/smallVars_CNV_output/variant_calling_ploidy%i.tab"%(VarCallOutdirs, sampleID, ploidy)
            partial_varcall_file = "%s/partial_variant_calling.tab"%cwd

            cutting_cols_stderr = "%s.generating.stderr"%partial_varcall_file
            print_if_verbose("getting the important cols. The stderr is in %s"%cutting_cols_stderr)
            run_cmd("cut -f%s %s > %s 2>%s"%(variant_calling_colNames, target_varcall_file, partial_varcall_file, cutting_cols_stderr))
            remove_file(cutting_cols_stderr)

            # load df
            df = pd.read_csv(partial_varcall_file, sep="\t")[variant_calling_fields]
            remove_file(partial_varcall_file)

            # append the sample ID 
            df["sampleID"] = sampleID

            # keep
            df_variant_calling = df_variant_calling.append(df)

            # print the size
            print_if_verbose("Size of df_variant_calling: %.2f MB"%(sys.getsizeof(df_variant_calling)/1000000))

        # save
        variant_calling_df_file_tmp = "%s.tmp"%variant_calling_df_file
        df_variant_calling.to_csv(variant_calling_df_file_tmp, sep="\t", header=True, index=False)
        os.rename(variant_calling_df_file_tmp, variant_calling_df_file)

    else: variant_calling_df = pd.read_csv(variant_calling_df_file, sep="\t")


    ######################################

    ###### INTEGRATE VARIANT ANNOTATION ######

    # define the file
    variant_annotation_df_file = "%s/integrated_variant_annotation_ploidy%i.tab"%(cwd, ploidy)

    if file_is_empty(variant_annotation_df_file) or repeat is True:
        print_if_verbose("generating integrated variant annotation")

        # initialize df
        df_variant_annotation = pd.DataFrame()

        # initialize the previous vars
        already_saved_vars = set()

        for Is, sampleID in enumerate(samples_to_run):
            print_if_verbose("%i/%i: %s"%(Is+1, len(samples_to_run), sampleID))

            # load df
            df = pd.read_csv("%s/%s_VarCallresults/smallVars_CNV_output/variant_annotation_ploidy%i.tab"%(VarCallOutdirs, sampleID, ploidy), sep="\t")

            # get only the new vars
            df_new = df[~df["#Uploaded_variation"].isin(already_saved_vars)]

            # keep 
            if len(df_new)>0: df_variant_annotation = df_variant_annotation.append(df_new)

            # define the already existing vars
            already_saved_vars = set(df_variant_annotation["#Uploaded_variation"])

            # print the size
            print_if_verbose("Size of df_variant_annotation: %.2f MB"%(sys.getsizeof(df_variant_annotation)/1000000))


        # sort
        df_variant_annotation = df_variant_annotation.sort_values(by="#Uploaded_variation").drop_duplicates()

        # add some fields
        """
        df_variant_annotation["chromosome"] = df_variant_annotation["#Uploaded_variation"].apply(lambda x: "_".join(x.split("_")[0:-2]))
        df_variant_annotation["position"] =  df_variant_annotation["#Uploaded_variation"].apply(lambda x: x.split("_")[-2]).apply(int)
        df_variant_annotation["ref"] = df_variant_annotation["#Uploaded_variation"].apply(lambda x: x.split("_")[-1].split("/")[0])
        df_variant_annotation["alt"] = df_variant_annotation["#Uploaded_variation"].apply(lambda x: x.split("_")[-1].split("/")[1])
        """

        # save
        variant_annotation_df_file_tmp = "%s.tmp"%variant_annotation_df_file
        df_variant_annotation.to_csv(variant_annotation_df_file_tmp, sep="\t", header=True, index=False)
        os.rename(variant_annotation_df_file_tmp, variant_annotation_df_file)


    ######################################




    return variant_calling_df



def run_repeat_modeller(reference_genome, threads=4, replace=False):

    """Runs repeat modeller to get a fasta file were all the repeats are. 

    (not tested) The part of integrating the results does not work very well, so that we have to integrate the results into one with RepeatClassifier afterwards.

    I had to improve the databases.
    In the RepeatMasker/Libraries I ran ' makeblastdb -in RepeatPeps.lib -input_type fasta -dbtype prot' to get the proteins, and then copied to RepeatModeler/Libraries (cp -r * ../../RepeatModeler/Libraries/)

    and also for the lib  makeblastdb -in RepeatMasker.lib -input_type fasta -dbtype nucl

    and makeblastdb -in simple.lib  -input_type fasta -dbtype nucl in RepeatModeler

    """

    # get the genome into outdir
    outdir = get_fullpath("%s.repeat_modeler_outdir"%reference_genome)

    # define the final files
    genome_dir = "%s/reference_genome.fasta"%outdir
    repeat_modeler_outfile = "%s-families.fa"%genome_dir

    if file_is_empty(repeat_modeler_outfile) or replace is True:
        print_if_verbose("running repeat modeler")

        # delete the outdir
        delete_folder(outdir)

        # create it
        make_folder(outdir)

        # put the genome under this outdir
        shutil.copy2(reference_genome, genome_dir)

        # run the database
        name_database = get_file(genome_dir)

        bulding_repModeler_db_std = "%s.genearting_db.std"%genome_dir
        print_if_verbose("getting repeat modeler db. The std is in %s"%bulding_repModeler_db_std)

        run_cmd("cd %s && %s -name %s %s > %s 2>&1"%(outdir, repeat_modeller_BuildDatabase, name_database, genome_dir, bulding_repModeler_db_std))

        remove_file(bulding_repModeler_db_std)

        # run repeatmodeller
        njobs = int(threads/4) # Specify the number of parallel search jobs to run. RMBlast jobs wil use 4 cores each and ABBlast jobs will use a single core each. i.e. on a machine with 12 cores and running with RMBlast you would use -pa 3 to fully utilize the machine
        print_if_verbose("Running repeat modeller in %s on %i jobs"%(outdir, njobs))

        #raise ValueError("This has to be fixed!!!!")
        cmd = "export PERL5LIB=%s && cd %s && %s -database %s -pa %i -LTRStruct -debug"%(repeatmoder_dir, outdir, repeat_modeller, name_database, njobs)

        # add the location were eveything is installed and run
        repeatmodeler_std = "%s/repeatmodeler.std"%outdir
        print_if_verbose("running repeatmodeler. The std is in %s"%repeatmodeler_std)        
        cmd += " -abblast_dir %s -cdhit_dir %s -genometools_dir %s -ltr_retriever_dir %s -mafft_dir %s -ninja_dir %s -recon_dir %s -repeatmasker_dir %s -rmblast_dir %s -rscout_dir %s -trf_prgm %s > %s 2>&1"%(abblast_dir, cdhit_dir, genometools_dir, ltr_retriever_dir, mafft_dir, ninja_dir, recon_dir, repeatmasker_dir, rmblast_dir, rscout_dir, trf_prgm_dir, repeatmodeler_std)

        run_cmd(cmd)

        if file_is_empty(repeat_modeler_outfile): 

            # test that there are no families identified
            no_families_identified = any([l.startswith("No families identified") for l in open(repeatmodeler_std, "r").readlines()[-3:]])

            # test that there are errors
            errors_in_repeatModeler = any(["ERROR" in l.upper() for l in open(repeatmodeler_std, "r").readlines()])

            # check if there are no families
            if no_families_identified and not errors_in_repeatModeler: open(repeat_modeler_outfile, "w").write("no_families_identified")

            else: raise ValueError("RepeatModeler did not end properly. Check %s for the std"%repeatmodeler_std)

        remove_file(repeatmodeler_std)

    # check if any new families were identified
    new_families_identified = open(repeat_modeler_outfile, "r").readlines()[0].strip()!="no_families_identified"

    # remove the folder
    for f in os.listdir(outdir):
        path = "%s/%s"%(outdir, f)
        if os.path.isdir(path) and f.startswith("RM_"): delete_folder(path)

    return repeat_modeler_outfile, new_families_identified


def run_repeat_masker(reference_genome, threads=4, replace=False, use_repeat_modeller=True):

    """
    It runs repeat masker for a reference genome, writing the results under a folder where the ref genome is
    """

    # get the reference genome as a full path
    reference_genome = get_fullpath(reference_genome)

    # get the library from repeat_modeller
    if use_repeat_modeller is True: library_repeats_repeatModeller, new_families_identified =  run_repeat_modeller(reference_genome, threads=threads, replace=replace)

    # define the repear masker outdir
    genome_dir = "/".join(reference_genome.split("/")[0:-1])
    genome_name = reference_genome.split("/")[-1]

    # define the outdirs for each 
    repeat_masker_outdir = "%s/%s_repeat_masker_outdir"%(genome_dir, genome_name.split(".")[0]); make_folder(repeat_masker_outdir)
    repeat_masker_outdir_default = get_fullpath("%s/default"%repeat_masker_outdir); make_folder(repeat_masker_outdir_default)
    repeat_masker_outdir_personal = get_fullpath("%s/personal"%repeat_masker_outdir); make_folder(repeat_masker_outdir_personal)

    # run in the default configuration
    repeat_masker_outfile_default = "%s/%s.out"%(repeat_masker_outdir_default, genome_name)
    repeat_masker_std_default = "%s/%s.std.out"%(repeat_masker_outdir_default, genome_name)

    if file_is_empty(repeat_masker_outfile_default) or replace is True:
        print_if_verbose("running repeatmasker to get the repeats of the genome in the default configuration. The std is in %s"%repeat_masker_std_default)
        run_cmd("cd %s && %s -pa %i -dir %s -poly -html -gff %s > %s 2>&1"%(repeat_masker_outdir_default, repeat_masker, threads, repeat_masker_outdir_default, reference_genome, repeat_masker_std_default))
        remove_file(repeat_masker_std_default)

    # run in the personal configuration
    repeat_masker_outfile_personal = "%s/%s.out"%(repeat_masker_outdir_personal, genome_name)
    repeat_masker_std_personal = "%s/%s.std.out"%(repeat_masker_outdir_personal, genome_name)

    if use_repeat_modeller is True and new_families_identified is True:
        
        if file_is_empty(repeat_masker_outfile_personal) or replace is True:
            print_if_verbose("running repeatmasker to get the repeats of the genome with the lib obtained with RepeatModeler. The std is in %s"%repeat_masker_std_personal)
            run_cmd("cd %s && %s -pa %i -dir %s -poly -html -gff -lib %s %s > %s 2>&1"%(repeat_masker_outdir_personal, repeat_masker, threads, repeat_masker_outdir_personal, library_repeats_repeatModeller, reference_genome, repeat_masker_std_personal))
            
    else: 

        # empty file
        print_if_verbose("avoiding the generation of the personal repeat masker. The stderr can be found in %s"%repeat_masker_std_personal)
        run_cmd("head -n 3 %s > %s 2>%s"%(repeat_masker_outfile_default, repeat_masker_outfile_personal, repeat_masker_std_personal))

    remove_file(repeat_masker_std_personal)

       
    return repeat_masker_outfile_personal, repeat_masker_outfile_default

def get_repeat_maskerDF(reference_genome, threads=4, replace=False):

    """gets the repeat masker outfile as a pandas df. The repeatmasker locations are 1-based (https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/rmsk2bed.html)"""

    # define the table 
    repeats_table_file = "%s.repeats.tab"%reference_genome
    repeats_table_file_tmp = "%s.tmp"%repeats_table_file

    if file_is_empty(repeats_table_file) or replace is True:
        print_if_verbose("running RepeatModeler and RepeatMasker into %s"%repeats_table_file)

        # get the file
        repeat_masker_outfile_personal, repeat_masker_outfile_default = run_repeat_masker(reference_genome, threads=threads, replace=replace, use_repeat_modeller=True)

        # load both dfs
        df_list = []
        for repeat_masker_outfile in [repeat_masker_outfile_personal, repeat_masker_outfile_default]:

            # define the header
            header = ["SW_score", "perc_div", "perc_del", "perc_ins", "chromosome", "begin_repeat", "end_repeat", "left_repeat", "strand", "repeat", "type", "position_inRepeat_begin", "position_inRepeat_end", "left_positionINrepeat",  "IDrepeat"]
            function_map_type = [int, float, float, float, str, int, int, str, str, str, str, str, str, str, int]

            # initialize header
            dict_repeat_masker = {}

            # go line by line
            for IDline, line in enumerate(open(repeat_masker_outfile, "r").readlines()):

                # debug lines
                split_line = line.strip().split()
                if len(split_line)==0 or split_line[0] in {"SW", "score"}: continue
                if split_line[-1]=="*": del split_line[-1]

                # keep
                line_content = [function_map_type[I](content) for I, content in enumerate(split_line)]
                dict_repeat_masker[IDline] = dict(zip(header, line_content))

            # get as df
            df = pd.DataFrame(dict_repeat_masker).transpose()
            df_list.append(df)

        # get df
        print_if_verbose("getting both repeats df")
        df = pd.concat(df_list).sort_values(by=["chromosome", "begin_repeat", "end_repeat"])
        df = df.drop_duplicates(subset=[x for x in df.keys() if x not in {"IDrepeat"}], keep="first")
        df["IDrepeat"] = list(range(1, len(df)+1))
        df.index = list(range(1, len(df)+1))

        # write
        df.to_csv(repeats_table_file_tmp, sep="\t", header=True, index=False)
        os.rename(repeats_table_file_tmp, repeats_table_file)

        print_if_verbose("repeat inference finished")

    else: df = pd.read_csv(repeats_table_file, sep="\t")

    return df, repeats_table_file
