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
from collections import OrderedDict
import warnings
from Bio import SeqIO
from scipy import linalg
from math import ceil
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle
import cylowess 
import itertools
import ast
import copy as cp
import re
import logging
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
from sklearn.model_selection import KFold
import collections
import scipy.interpolate as scipy_interpolate
#from ete3 import Tree, NCBITaxa
from shutil import copyfile
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
import psutil
from sklearn.utils import resample

import plotly.plotly as py
import plotly.figure_factory as ff
import plotly.offline as off_py
import plotly.graph_objs as go
from plotly import tools
from plotly.offline import init_notebook_mode, plot, iplot # download_plotlyjs
import cufflinks as cf

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
# import inspect
#import statsmodels.api as statsmodels_api
#import igraph


# chnage the errors to report
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

# import modules
import pyloess_Loess as pyloess_fun

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# define the conda base dir and the env name
CondaDir =  "/".join(sys.executable.split("/")[0:-4])
EnvName = EnvDir.split("/")[-1]

# define the R env name
EnvName_R = "%s_R_env"%EnvName
EnvName_gridss = "%s_gridss_env"%EnvName
EnvName_picard = "%s_picard_env"%EnvName
EnvName_ete3 = "%s_ete3_3.0.0_env"%EnvName
EnvName_RepeatMasker = "%s_RepeatMasker_env"%EnvName
EnvName_CONY = "%s_CONY_env"%EnvName
EnvName_HMMcopy = "%s_HMMcopy_env"%EnvName
EnvName_AneuFinder = "%s_AneuFinder_env"%EnvName

# define other envDirs
picard_EnvDir = "%s/envs/%s"%(CondaDir, EnvName_picard)
RepeatMasker_EnvDir = "%s/envs/%s"%(CondaDir, EnvName_RepeatMasker)

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
vcffilter = "%s/bin/vcffilter"%EnvDir
freebayes_parallel = "%s/bin/freebayes-parallel"%EnvDir
fasta_generate_regions_py = "%s/bin/fasta_generate_regions.py"%EnvDir
wgsim = "%s/bin/wgsim"%EnvDir
svim = "%s/bin/svim"%EnvDir
mosdepth = "%s/bin/mosdepth"%EnvDir
perl = "%s/bin/perl"%EnvDir
esearch = "%s/bin/esearch"%EnvDir
efetch = "%s/bin/efetch"%EnvDir
prefetch = "%s/bin/prefetch"%EnvDir
fastqdump = "%s/bin/fastq-dump"%EnvDir
parallel_fastq_dump = "%s/bin/parallel-fastq-dump"%EnvDir
FASTQC = "%s/bin/fastqc"%EnvDir
seqtk = "%s/bin/seqtk"%EnvDir
fasta_generate_regions_py = "%s/bin/fasta_generate_regions.py"%EnvDir
pigz = "%s/bin/pigz"%EnvDir
unpigz = "%s/bin/unpigz"%EnvDir
bedmap = "%s/bin/bedmap"%EnvDir
cnvpytor_exec = "%s/bin/cnvpytor"%EnvDir
genmap = "%s/bin/genmap"%EnvDir
porechop = "%s/bin/porechop"%EnvDir
sniffles = "%s/bin/sniffles"%EnvDir

# executables that are provided in the repository
external_software = "%s/../installation/external_software"%CWD
gridss_run = "%s/gridss.sh"%external_software
gridss_jar = "%s/gridss-2.9.2-gridss-jar-with-dependencies.jar"%external_software
clove = "%s/clove-0.17-jar-with-dependencies.jar"%external_software
gztool = "%s/gztool"%external_software

# executables that are in other environments
picard_exec = "%s/bin/picard"%picard_EnvDir

repeatmoder_dir = "%s/share/RepeatModeler"%RepeatMasker_EnvDir
repeat_modeller = "%s/bin/RepeatModeler"%RepeatMasker_EnvDir
repeatmasker_dir = "%s/share/RepeatMasker"%RepeatMasker_EnvDir
makeblastdb = "%s/bin/makeblastdb"%RepeatMasker_EnvDir
repeat_modeller_BuildDatabase = "%s/bin/BuildDatabase"%RepeatMasker_EnvDir
abblast_dir = "%s/bin"%RepeatMasker_EnvDir
cdhit_dir = "%s/bin"%RepeatMasker_EnvDir
trf_prgm_dir = "%s/bin/trf"%RepeatMasker_EnvDir
ltr_retriever_dir = "%s/bin"%RepeatMasker_EnvDir
genometools_dir = "%s/bin"%RepeatMasker_EnvDir
repeat_masker = "%s/bin/RepeatMasker"%RepeatMasker_EnvDir
mafft_dir = "%s/bin"%RepeatMasker_EnvDir
recon_dir = "%s/bin"%RepeatMasker_EnvDir
rmblast_dir = "%s/bin"%RepeatMasker_EnvDir
rscout_dir = "%s/bin"%RepeatMasker_EnvDir
ninja_dir = "%s/bin"%RepeatMasker_EnvDir

# define the bcftools=1.10 by activating the conda env
#SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
# CONDA_ACTIVATING_CMD = "conda activate %s;"%EnvName
#bcftools_latest = "%s && conda activate %s_bcftools_1.10.2_env && bcftools"%(SOURCE_CONDA_CMD, EnvName)
#bcftools_latest_cmd = "%s && conda activate %s_bcftools_1.10.2_env && bcftools"%(SOURCE_CONDA_CMD, EnvName)


# scripts that are of this pipeline
create_random_simulatedSVgenome_R = "%s/create_random_simulatedSVgenome.R"%CWD
create_targeted_simulatedSVgenome_R = "%s/create_targeted_simulatedSVgenome.R"%CWD
annotate_simpleEvents_gridssVCF_R = "%s/annotate_simpleEvents_gridssVCF.R"%CWD
run_CONY_R = "%s/run_CONY.R"%CWD
run_HMMcopy_R = "%s/run_HMMcopy.R"%CWD
run_AneuFinder_R = "%s/run_ANEUFINDER.R"%CWD
analyze_svVCF = "%s/generate_files_from_svVCF.R"%CWD
analyze_svVCF_simple = "%s/generate_files_from_svVCF_simple.R"%CWD
TRIMMOMATIC = "%s/run_trimmomatic.py"%CWD 
perSVade_py = "%s/perSVade.py"%CWD 
run_trimmomatic_and_fastqc_py = "%s/run_trimmomatic_and_fastqc.py"%CWD
get_trimmed_reads_for_srr_py = "%s/get_trimmed_reads_for_srr.py"%CWD
run_vep = "%s/run_vep.py"%CWD
calculate_memory_py = "%s/calculate_memory.py"%CWD
get_interestingTaxIDs_distanceToTarget_taxID_to_sciName_py = "%s/get_interestingTaxIDs_distanceToTarget_taxID_to_sciName.py"%CWD
libraries_CONY = "%s/CONY_package_debugged.R "%CWD
run_svim_and_sniffles_py = "%s/run_svim_and_sniffles.py"%CWD

# old code
#create_simulatedSVgenome_R = "%s/create_simulatedSVgenome.R"%CWD
#bbmap_reformat_sh = "%s/bin/reformat.sh"%EnvDir


######################################################
######################################################

####################################
######## DEFINE VARIABLES ##########
####################################

# types of small variants
ALL_MUTATIONS = {'stop_gained', 'intron_variant', 'upstream_gene_variant', '5_prime_UTR_variant', 'inframe_insertion', 'synonymous_variant', 'non_coding_transcript_exon_variant', 'intergenic_variant', 'protein_altering_variant', 'coding_sequence_variant', 'downstream_gene_variant', '3_prime_UTR_variant', 'missense_variant', 'splice_region_variant', 'splice_acceptor_variant', 'inframe_deletion', 'stop_lost', 'non_coding_transcript_variant', 'start_retained_variant', 'frameshift_variant', 'stop_retained_variant', 'start_lost', 'incomplete_terminal_codon_variant', 'splice_donor_variant', '-'}

PROT_ALTERRING_MUTATIONS = {'missense_variant', 'start_lost', 'inframe_deletion', 'protein_altering_variant', 'stop_gained', 'inframe_insertion', 'frameshift_variant', 'stop_lost', 'splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant'}

NON_PROT_ALTERRING_MUTATIONS = ALL_MUTATIONS.difference(PROT_ALTERRING_MUTATIONS)


# types of SVs
SVs_ALL_MUTATIONS = {'coding_sequence_variant_BND', 'upstream_gene_variant_BND', '3_prime_UTR_variant', 'feature_elongation', 'feature_truncation', 'coding_sequence_variant', 'intergenic_variant', 'upstream_gene_variant', '5_prime_UTR_variant', 'transcript_amplification', '5_prime_UTR_variant_BND', 'downstream_gene_variant', 'intron_variant_BND', 'intron_variant', 'intergenic_variant_BND', '3_prime_UTR_variant_BND', 'downstream_gene_variant_BND', 'non_coding_transcript_exon_variant_BND', 'non_coding_transcript_exon_variant', 'transcript_ablation', '-', "inframe_insertion", "frameshift_variant"}

SVs_TRANSCRIPT_DISRUPTING_MUTATIONS = {'coding_sequence_variant_BND', 'feature_elongation', 'feature_truncation', 'coding_sequence_variant', 'transcript_amplification', 'transcript_ablation', 'intron_variant_BND', 'non_coding_transcript_exon_variant_BND', 'intron_variant', 'non_coding_transcript_exon_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', '5_prime_UTR_variant_BND', '3_prime_UTR_variant_BND', "inframe_insertion", "frameshift_variant"}

SVs_NON_TRANSCRIPT_DISRUPTING_MUTATIONS = SVs_ALL_MUTATIONS.difference(SVs_TRANSCRIPT_DISRUPTING_MUTATIONS)

# variant representation data
sorted_consequences = ["downstream_gene_variant", "upstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "stop_retained_variant", "synonymous_variant", "missense_variant", "stop_gained", "frameshift_variant"]

consequence_to_abbreviation  = {'frameshift_variant':"FS",
                                 'inframe_deletion':"del",
                                 'inframe_insertion':"ins",
                                 'missense_variant':"mis",
                                 'protein_altering_variant':"FS",
                                 'splice_acceptor_variant':"spliceAcc",
                                 'splice_donor_variant':"spliceDon",
                                 'splice_region_variant':"spliceReg",
                                 'start_lost':"lostATG",
                                 'stop_gained':"PTC",
                                 'stop_lost':"lostSTOP",
                                 'non_coding_transcript_exon_variant':"nonCodExon",
                                 '3_prime_UTR_variant':"UTR3",
                                 '5_prime_UTR_variant':"UTR5",
                                 'downstream_gene_variant':"down",
                                 'intron_variant':"intr",
                                 'non_coding_transcript_variant':"nonCodTrans",
                                 'start_retained_variant':"retATG",
                                 'stop_retained_variant':"retSTOP",
                                 'synonymous_variant':"syn",
                                 'upstream_gene_variant':"up",
                                 'intergenic_variant':"ig", 
                                 'incomplete_terminal_codon_variant':"incTermCodonVar",
                                 "coding_sequence_variant":"codSeqVar",
                                 '':"empty"}

# define a set of priorities of the variants                                                   
var_to_IDX  =  {'frameshift_variant':0,
                'stop_gained':1,
                'start_lost':2,
                'protein_altering_variant':3,
                'inframe_deletion':4,
                'inframe_insertion':5,
                'stop_lost':6,
                'missense_variant':7,
                'splice_acceptor_variant':8,
                'splice_donor_variant':8,
                'splice_region_variant':9,
                '3_prime_UTR_variant':10,
                'downstream_gene_variant':11,
                'upstream_gene_variant':11,
                '5_prime_UTR_variant':11,
                'intron_variant':12,
                'synonymous_variant':13,
                'coding_sequence_variant':13,
                'incomplete_terminal_codon_variant':14,
                'start_retained_variant':14,
                'stop_retained_variant':15,
                'non_coding_transcript_exon_variant':16,
                'non_coding_transcript_variant':17,
                'intergenic_variant':18,
                '':19}

# map each type of variant to either CDS position, Codons or Protein_position, Amino_acids
info_codons = ("c", "CDS_position", "Codons")
info_aa = ("p", "Protein_position", "Amino_acids")
info_spliciingAndNoncoding = ("g", "CDS_position", "#Uploaded_variation")

protVar_to_info =  {'frameshift_variant':info_aa,
                             'inframe_deletion': info_aa,
                             'inframe_insertion': info_aa,
                             'missense_variant': info_aa,
                             'protein_altering_variant': info_aa,
                             'splice_acceptor_variant': info_spliciingAndNoncoding,
                             'splice_region_variant': info_spliciingAndNoncoding,
                             'start_lost': info_codons,
                             'stop_gained': info_aa,
                             'stop_lost': info_codons,
                             'non_coding_transcript_exon_variant': info_spliciingAndNoncoding,
                             'splice_donor_variant':info_spliciingAndNoncoding,
                             '3_prime_UTR_variant': info_spliciingAndNoncoding,
                             '5_prime_UTR_variant': info_spliciingAndNoncoding,
                             'downstream_gene_variant': info_spliciingAndNoncoding,
                             'non_coding_transcript_variant': info_spliciingAndNoncoding,
                             'start_retained_variant': info_spliciingAndNoncoding,
                             'stop_retained_variant': info_spliciingAndNoncoding,
                             'synonymous_variant': info_codons,
                             'upstream_gene_variant': info_spliciingAndNoncoding,
                             'intron_variant':info_spliciingAndNoncoding,
                             'intergenic_variant':info_spliciingAndNoncoding,
                             "coding_sequence_variant":info_codons,
                             "incomplete_terminal_codon_variant":info_spliciingAndNoncoding,
                             "":info_spliciingAndNoncoding}

# define the strings that have to be considered as NaN in the VCF parsing
vcf_strings_as_NaNs = ['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'n/a', 'nan', 'null']

# define default parameters for gridss filtering. This has changed from v0
default_filtersDict_gridss = {"min_Nfragments":5, "min_af":0.25, "wrong_FILTERtags":("NO_ASSEMBLY",), "filter_polyGC":True, "filter_noSplitReads":False, "filter_noReadPairs":False, "maximum_strand_bias":0.99, "maximum_microhomology":50, "maximum_lenght_inexactHomology":50, "range_filt_DEL_breakpoints":(0, 1), "min_length_inversions":40, "dif_between_insert_and_del":5, "max_to_be_considered_small_event":1000, "wrong_INFOtags":('IMPRECISE',), "min_size":50, "min_af_EitherSmallOrLargeEvent":0.25, "min_QUAL":0, "filter_overlappingRepeats":False} # the minimum af is 0.25 to include both heterozygous and homozygous vars as default

# define other default parameters
default_gridss_blacklisted_regions = ""
default_gridss_maxcoverage = 50000
default_max_rel_coverage_to_consider_del = 0.1
default_min_rel_coverage_to_consider_dup = 1.8

default_parms_dict_HMMcopy = dict(e=0.9999999, mu=[-0.458558247, -0.215877601, -0.002665686, 0.191051578, 0.347816046, 1.664333241], lambda_val=20, nu=2.1, kappa=[50, 50, 700, 100, 50, 50], m=[-0.458558247, -0.215877601, -0.002665686, 0.191051578, 0.347816046, 1.664333241], eta=50000, gamma=3, S=0.02930164, strength="1e7")

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
fractionRAM_to_dedicate = 0.5

# define the fraction of memory available
fraction_available_mem = None

# the minimum size of CNVs to be considered
min_CNVsize_coverageBased = 300

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

def run_cmd(cmd, env=EnvName):

    """This function runs a cmd with a given env"""

    # define the cmds
    SOURCE_CONDA_CMD = "source %s/etc/profile.d/conda.sh"%CondaDir
    cmd_prefix = "%s && conda activate %s &&"%(SOURCE_CONDA_CMD, env)

    # define the running
    cmd_to_run = "%s %s"%(cmd_prefix, cmd)

    # run
    out_stat = os.system(cmd_to_run) 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd_to_run, out_stat))

def get_weigthed_median(df, field, weight_field, type_algorithm="fast"):

    """This functon takes a df with a field. It returns the median weighted by weight_field (which should be a number, not a fraction)"""

    # set to int
    df[weight_field] = df[weight_field].apply(int)

    # check
    if max(df[weight_field])==1 or min(df[weight_field])<0: raise ValueError("this function can't be applied here")

    # unefficient calc
    if type_algorithm=="long":

        # get an array by weight
        weighted_array = make_flat_listOflists(df.apply(lambda r: [r[field]]*r[weight_field], axis=1))
        
        # get the median
        median = np.median(weighted_array)
       
    elif type_algorithm=="fast":

        # get the algorithm from SO
        df_sorted = df.sort_values(field)
        cumsum = df_sorted[weight_field].cumsum()
        cutoff = df_sorted[weight_field].sum() / 2.
        median = df_sorted[cumsum >= cutoff][field].iloc[0]

    else: raise ValueError("%s is not valid"%type_algorithm)

    return median


def get_median_coverage(coverage_df, mitochondrial_chromosome, coverage_field="mediancov_1"):

    """This function takes a coverage df and calculates the median for those non-0 coverage windows. It will return the median coverage weighted by length """

    # keep
    if len(coverage_df)==0: raise ValueError("there should be some data in the coverage_df")
    coverage_df = cp.deepcopy(coverage_df)

    # define the chrom field
    if "#chrom" in coverage_df.keys(): chrom_f = "#chrom"
    elif "chromosome" in coverage_df.keys(): chrom_f = "chromosome"

    # define the set
    mitochondrial_chromosomes = set(mitochondrial_chromosome.split(","))

    # get the nuclear and not 0 idx
    if set(coverage_df[chrom_f])==mitochondrial_chromosomes: idx = (coverage_df[coverage_field]>0)
    else: idx = (~coverage_df[chrom_f].isin(mitochondrial_chromosomes)) & (coverage_df[coverage_field]>0)

    df = coverage_df[idx]

    # get the len
    df["length"] = df.end - df.start

    # get the median wwighted by len
    median = get_weigthed_median(df, coverage_field, "length", type_algorithm="fast")

    return median

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

def get_chr_to_len(genome, replace=False):

    chr_to_len_file = "%s.chr_to_len.py"%genome
    chr_to_len_file_tmp = "%s.tmp"%chr_to_len_file

    if file_is_empty(chr_to_len_file) or replace is True:

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

def get_availableGbRAM(outdir):

    """This function returns a float with the available memory in your system. psutil.virtual_memory().available/1e9 would yield the theoretically available memory.

    fraction_available_mem is the fraction of the total memory available in the node (computer) that is dedicated to the run of perSVade."""

    # get the memory with psutil
    available_mem = psutil.virtual_memory().available/1e9

    # define the fraction_total_running 
    if fraction_available_mem is None:

        # MN4
        if "BSC_MACHINE" in os.environ and os.environ["BSC_MACHINE"]=="mn4":

            availableGbRAM = available_mem*(int(os.environ["SLURM_CPUS_PER_TASK"])/48)

        # nord3
        elif "BSC_MACHINE" in os.environ and os.environ["BSC_MACHINE"]=="nord3":

            # define the available threads
            real_available_threads = get_available_threads(outdir)

            # get thr ram considering that 1 node has 16 threads
            availableGbRAM = available_mem*(int(real_available_threads)/16)

        # BSC machine
        elif str(subprocess.check_output("uname -a", shell=True)).startswith("b'Linux bscls063 4.12.14-lp150.12.48-default"): 

            availableGbRAM = available_mem

        # others. Calculate from filling the memory
        else:
            print_if_verbose("WARNING: calculating the available RAM by filling it. This may be dangerous.")

            # define a file where to calculate memory
            calc_memory_file = "%s/calculating_memory.txt"%outdir
            remove_file(calc_memory_file)

            # run the calculate memory
            cmd_std = "%s.std"%calc_memory_file
            try: run_cmd("%s --outfile %s > %s 2>&1"%(calculate_memory_py, calc_memory_file, cmd_std))
            except: print_if_verbose("memory was calculated")

            # get from file
            availableGbRAM = float(open(calc_memory_file, "r").readlines()[-1])

            # at the end remove the file
            remove_file(calc_memory_file)
            remove_file(cmd_std)

    else: availableGbRAM = available_mem*float(fraction_available_mem)

    return availableGbRAM


def simulate_testing_reads_on_genome(genome, window_l=2000, npairs=50000, read_length=150, median_insert_size=250, median_insert_size_sd=50, threads=4, replace=False):

    """ 
    Takes a genome and simulates reads for it, saving them under <genome>_simulating_reads 
    """

    # define the outdir
    outdir = "%s_simulating_reads"%genome; 
    outdir_reads = "%s/getting_reads"%outdir; 
    
    # remove the outdirs if replace is True
    if replace is True: 
        delete_folder(outdir)
        delete_folder(outdir_reads)

    # make folders 
    make_folder(outdir)
    make_folder(outdir_reads)

    # define the expected reads
    reads1 = "%s/all_reads1.correct.fq.gz"%outdir_reads
    reads2 = "%s/all_reads2.correct.fq.gz"%outdir_reads

    if any([file_is_empty(f) for f in [reads1, reads2]]):

        # run index the genome
        run_cmd("%s faidx %s"%(samtools, genome))

        # get the windows df
        windows_bed = "%s/windows_file.bed"%outdir
        run_cmd("%s makewindows -g %s.fai -w %i > %s"%(bedtools, genome, window_l, windows_bed))
        df_windows = pd.read_csv(windows_bed, sep="\t", header=-1, names=["chromosome", "start", "end"])
        df_windows["predicted_relative_coverage"] = random.sample(list(np.linspace(0.5, 2, 10000)), len(df_windows))

        # simulate reads
        simulate_readPairs_per_window(df_windows, genome, npairs, outdir_reads, read_length, median_insert_size, median_insert_size_sd, replace=False, threads=threads) 

    return reads1, reads2


def get_sorted_bam_test(r1, r2, ref_genome, replace=False):

    """Runs bwa mem on the reads and returns the sorted bam with marked duplicates"""

    # define the outdir
    outdir = "%s/aligning_reads_against_%s"%(get_dir(r1), get_file(ref_genome))

    # if replace is True, delete the outdir
    if replace is True: delete_folder(outdir)

    # make de outdir
    make_folder(outdir)

    # define the inputs of bam
    bamfile = "%s/aligned_reads.bam"%outdir
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam
    name_sample = "test_sample"

    # run
    run_bwa_mem(r1, r2, ref_genome, outdir, bamfile, sorted_bam, index_bam, name_sample, threads=4, replace=False, MarkDuplicates=False)

    return sorted_bam


def get_available_threads(outdir):

    """Returns the avilable number of threads by runnning GATK. It runs gatk on a dummy genome and returns the really available number of threads  """

    # MN4
    if "BSC_MACHINE" in os.environ and os.environ["BSC_MACHINE"]=="mn4":

        available_threads = int(os.environ["SLURM_CPUS_PER_TASK"])

    # Nord3 
    elif "BSC_MACHINE" in os.environ and os.environ["BSC_MACHINE"]=="nord3" and not "LSB_MCPU_HOSTS" in os.environ:

        available_threads = 4

    # Nord3 interactive nodes
    elif "BSC_MACHINE" in os.environ and os.environ["BSC_MACHINE"]=="nord3" and "LSB_MCPU_HOSTS" in os.environ:

        available_threads = int(os.environ["LSB_MCPU_HOSTS"].split()[-1])

    # BSC machine
    elif str(subprocess.check_output("uname -a", shell=True)).startswith("b'Linux bscls063 4.12.14-lp150.12.48-default"): 

        available_threads = 4

    # others. Calculate by running GATK
    else:

        # redefnie the outdir under outdir
        outdir = "%s/getting_available_threads"%(outdir)
        delete_folder(outdir)
        make_folder(outdir)

        # define a genome that has one chromosome
        genome = "%s/genome.fasta"%outdir
        genome_obj = SeqRecord(Seq("ACTGCGATCGACTCGATCGATGAGAGAGAGGACTCTCAACAG"*10), id="chromosomeX")
        SeqIO.write([genome_obj], genome, "fasta")

        # get some simulated reads
        reads1, reads2 = simulate_testing_reads_on_genome(genome, window_l=75, npairs=1000, read_length=50, median_insert_size=15, median_insert_size_sd=5, threads=4, replace=False)

        # get a sorted bam
        sorted_bam = get_sorted_bam_test(reads1, reads2, genome, replace=False)

        # create the files
        create_sequence_dict(genome, replace=False)

        # run GATK HC 
        gatk_out = "%s/output_HC.vcf"%outdir
        gatk_std = "%s.running.std"%gatk_out

        gatk_cmd = "%s HaplotypeCaller -R %s -I %s -O %s -ploidy %i --genotyping-mode DISCOVERY --emit-ref-confidence NONE --stand-call-conf 30 --native-pair-hmm-threads %i > %s 2>&1"%(gatk, genome, sorted_bam, gatk_out, 1, 100000000, gatk_std)

        run_cmd(gatk_cmd)

        # get the available threads
        threads_lines = [l for l in open(gatk_std, "r").readlines() if "IntelPairHmm - Using" in l and "available threads, but" in l and "were requested" in l]
        if len(threads_lines)!=1: raise ValueError("the threads were not properly calculated")

        available_threads = int(threads_lines[0].split("IntelPairHmm - Using ")[1].split("available threads")[0])

        # print
        print_if_verbose("there are %i available threads in this run"%available_threads)

        # remove the outdir
        delete_folder(outdir)

    return available_threads

def write_coverage_per_gene_mosdepth_and_parallel(sorted_bam, reference_genome, cnv_outdir, bed, gene_to_coverage_file, replace=False, threads=4):

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
        coverage_file = generate_coverage_per_window_file_parallel(reference_genome, cnv_outdir, sorted_bam, windows_file=cnv_bed, replace=replace, run_in_parallel=True, delete_bams=True, threads=threads)

        # add the ID
        df_bed = pd.read_csv(bed, sep="\t", header=-1, names=["chromosome", "start", "end", "ID"])
        df_coverage = pd.read_csv(coverage_file, sep="\t")
        df_coverage_with_ID = df_coverage.merge(df_bed, how="right", left_on=["#chrom", "start", "end"], right_on=["chromosome", "start", "end"], right_index=False, left_index=False, validate="one_to_many")[['chromosome', 'start', 'end', 'length', 'mediancov_1', 'nocoveragebp_1', 'percentcovered_1', 'ID']].rename(columns={"mediancov_1":"median_reads_per_gene"})

        # add fields
        df_coverage_with_ID["fraction_covered_by_MoreThan1read"] = df_coverage_with_ID.percentcovered_1 / 100
        df_coverage_with_ID["relative_coverage"] = df_coverage_with_ID.median_reads_per_gene / np.median(df_coverage_with_ID[df_coverage_with_ID.median_reads_per_gene>0].median_reads_per_gene)

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

    if file_is_empty(sorted_bam) or replace is True:

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


        if MarkDuplicates is True:
            bamfile_MarkedDuplicates = get_bam_with_duplicatesMarkedSpark(bamfile, threads=threads, replace=replace)

        else: bamfile_MarkedDuplicates = bamfile

        if file_is_empty(sorted_bam) or replace is True:
            print_if_verbose("Sorting bam")

            # remove all temporary files generated previously in samtools sort (they'd make a new sort to be an error)
            for outdir_file in os.listdir(outdir): 
                fullfilepath = "%s/%s"%(outdir, outdir_file)
                if outdir_file.startswith("aligned_reads") and ".tmp." in outdir_file: os.unlink(fullfilepath)

            # sort
            sorted_bam_tmp = "%s.tmp"%sorted_bam
            bam_sort_std = "%s.generating.txt"%sorted_bam
            print_if_verbose("the sorting bam std is in %s"%bam_sort_std)
            cmd_sort = "%s sort --threads %i -o %s %s > %s 2>&1"%(samtools, threads, sorted_bam_tmp, bamfile_MarkedDuplicates, bam_sort_std); run_cmd(cmd_sort)

            # remove files
            for f in [bam_sort_std, bamfile_MarkedDuplicates, bamfile]: remove_file(f)

            # rename
            os.rename(sorted_bam_tmp, sorted_bam)

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
        affected_region_bed_df["end"] = affected_region_bed_df.apply(lambda r: min([r["end"]+add_interval_bp, chr_to_len[r["chromosome"]]]), axis=1)

    return affected_region_bed_df[["chromosome", "start", "end"]], nSVs_in_interesting_chromosomes

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

def get_ChrB_bp_pos_translocations(r, chr_to_len, first_bp_pos=1):

    """Takes a row of a translocations df and returns the breakpoint position"""

    if r["StartB"]==first_bp_pos: return r["EndB"]
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

            # find the rearranged chromosomes A and B
            chrField_to_putativeChromosomes = {}
            for chrField in ["ChrA", "ChrB"]:
                
                # define the 3
                seq_3prime = sv_row["%s_3seq"%chrField]
                
                # go through each chrom
                chrField_to_putativeChromosomes[chrField] = [chrom for chrom, rearrangedSeq in chr_to_rearrangedSeq.items() if rearrangedSeq.count(seq_3prime)==1]

            # check that the 3' seq is unique
            if any([len(chrField_to_putativeChromosomes[chrField])!=1 for chrField in ["ChrA", "ChrB"]]) or chrField_to_putativeChromosomes["ChrA"][0]==chrField_to_putativeChromosomes["ChrB"][0]: 
                
                print_if_verbose("these are the putative chromosomes attributed to each sequence:")
                print_if_verbose(chrField_to_putativeChromosomes)
                print_if_verbose("WARNING: The sequences for %s are not unique enough to find the position of the bp in the rearranged genome"%ID)
                continue

            # define general parameters of the rearranged genome
            ChrA = chrField_to_putativeChromosomes["ChrA"][0]
            ChrB = chrField_to_putativeChromosomes["ChrB"][0]
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

def insert_balanced_translocations_into_rearranged_genome(reference_genome, input_rearranged_genome, output_rearranged_genome, svDF, translocations_file, svtype_to_svDF, replace=False):

    """This function takes a rearranged genome and insert the translocations generating output_rearranged_genome. This substitutes the translocations_file in case that some translocations cannot be inserted. svDF should contain translocations. The svDF should have 1s on it. svDF has to be in 1-based coordinates"""

    print_if_verbose("inserting translocations inHouse")

    # keep
    svtype_to_svDF = cp.deepcopy(svtype_to_svDF)

    # test that all are balanced translocations
    if not all(svDF.Balanced): raise ValueError("This has not been developed for unbalanced translocations")

    # rewrite translocations file to get them in the correct format
    rewrite_translocations_uniformizedFormat_simulateSV(translocations_file, reference_genome)
    svDF = pd.read_csv(translocations_file, sep="\t")

    # initialize the feasible translocation IDs
    feasible_translocation_IDs = set()

    # define a tmp genome
    output_rearranged_genome_tmp = "%s.tmp.fasta"%output_rearranged_genome

    # initialize the output_rearranged_genome as the input_rearranged genome
    shutil.copy2(input_rearranged_genome, output_rearranged_genome_tmp)

    # go through each translocation and insert it into the output_rearranged_genome
    for Itra, rTra in cp.deepcopy(svDF).iterrows():

        # get as df
        df_I = pd.DataFrame({0: rTra}).transpose()

        # add the currently inserted translocations to svtype_to_svDF
        if len(feasible_translocation_IDs)>0: svtype_to_svDF["translocations"] = cp.deepcopy(svDF.loc[feasible_translocation_IDs])
        else: svtype_to_svDF["translocations"] = pd.DataFrame(columns=svtype_to_fieldsDict["translocations"]["all_fields"])

        # get the svDF in coordinates of the rearranged genome
        df_I_rearrangedCoords = get_svDF_in_coords_of_rearranged_genome(df_I, reference_genome, output_rearranged_genome_tmp, "translocations", svtype_to_svDF)
        if len(df_I_rearrangedCoords)==0: continue

        # add fields to the rearrangedCoords df
        chr_to_rearranged_len = get_chr_to_len(output_rearranged_genome_tmp, replace=True)
        df_I_rearrangedCoords["orientation"] = df_I_rearrangedCoords.apply(lambda r: get_orientation_translocation(r, chr_to_rearranged_len), axis=1) 
        df_I_rearrangedCoords["ChrA_bp_pos"] = df_I_rearrangedCoords["EndA"]
        df_I_rearrangedCoords["ChrB_bp_pos"] = df_I_rearrangedCoords.apply(lambda r: get_ChrB_bp_pos_translocations(r, chr_to_rearranged_len), axis=1)

        # rewrite positions so that they are 0-based (so that each position is the real one)
        for pos_f in ["StartA", "EndA", "StartB", "EndB", "ChrA_bp_pos", "ChrB_bp_pos"]: df_I_rearrangedCoords[pos_f] = df_I_rearrangedCoords[pos_f] - 1

        # load the genome as a fasta
        seqID_to_seq = {seq.id : seq for seq in SeqIO.parse(output_rearranged_genome_tmp, "fasta")}

        ###### INSERT TRANSLOCATIONS ######

        # define as a series
        r = df_I_rearrangedCoords.iloc[0]

        # define chromosomes
        chrA_seq = cp.deepcopy(seqID_to_seq[r["ChrA"]])
        chrB_seq = cp.deepcopy(seqID_to_seq[r["ChrB"]])

        # define the new chromosomes
        if r["orientation"]=="5_to_5":

            new_chrA_seq = cp.deepcopy(chrA_seq[0:r["ChrA_bp_pos"]] + chrB_seq[r["ChrB_bp_pos"]:])
            new_chrB_seq = cp.deepcopy(chrB_seq[0:r["ChrB_bp_pos"]] + chrA_seq[r["ChrA_bp_pos"]:])

        elif r["orientation"]=="5_to_3":

            new_chrA_seq = cp.deepcopy(chrA_seq[0:r["ChrA_bp_pos"]] + chrB_seq[0:r["ChrB_bp_pos"]].reverse_complement())
            new_chrB_seq = cp.deepcopy(chrB_seq[r["ChrB_bp_pos"]:].reverse_complement() + chrA_seq[r["ChrA_bp_pos"]:])

        else: raise ValueError("%s is not valid"%r["orientation"])

        # check that the lengths match
        if (len(new_chrA_seq)+len(new_chrB_seq)) != (len(chrA_seq)+len(chrB_seq)): raise ValueError("the lens should match")

        # add the new chromosomes to the iterable
        new_chrA_seq.id = r["ChrA"]
        new_chrB_seq.id = r["ChrB"]

        new_chrA_seq.name = ""; new_chrA_seq.description = ""
        new_chrB_seq.name = ""; new_chrB_seq.description = ""

        # change in the original dict
        seqID_to_seq[r["ChrA"]] = new_chrA_seq
        seqID_to_seq[r["ChrB"]] = new_chrB_seq

        ###################################

        # print the len of each chromosome
        #print_if_verbose("translocation between %s and %s performed. These are the lens of the chroms:"%(r["ChrA"], r["ChrB"]))
        #for seqID, seq in seqID_to_seq.items(): print_if_verbose("%s %i"%(seqID, len(seq)))

        # write genome
        all_chroms = [s for s in seqID_to_seq.values()]
        SeqIO.write(all_chroms, output_rearranged_genome_tmp, "fasta")

        # record that this is feasible
        feasible_translocation_IDs.add(Itra)
        print_if_verbose("inserted translocation %i/%i"%(Itra+1, len(svDF)))

    # at the end rename
    os.rename(output_rearranged_genome_tmp, output_rearranged_genome)

    #### replace translocations file ####

    # write into translocations file (this is important so that in case any of the translocations could not be mapped from the sequence)
    svDF_final = svDF.loc[feasible_translocation_IDs][svtype_to_fieldsDict["translocations"]["all_fields"]]

    # delete the translocations file because it has to be rewritten as the end point of this function
    remove_file(translocations_file)
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

            if std_rearranging_genome!="stdout": run_cmd("%s > %s 2>&1"%(targetSV_cmd, std_rearranging_genome), env=EnvName_R)
            else: run_cmd(targetSV_cmd, env=EnvName_R)

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
        if "translocations" in svtype_to_svDF: 

            insert_balanced_translocations_into_rearranged_genome(reference_genome, rearranged_genome_InsInvDelTan, rearranged_genome, svtype_to_svDF["translocations"], translocations_file, svtype_to_svDF)

        # rewrite the variants so that they are optimal for comparison. This is important to re-sort the chromosomes if necessary
        print_if_verbose("rewriting %s"%translocations_file)
        if not file_is_empty(translocations_file): rewrite_translocations_uniformizedFormat_simulateSV(translocations_file, reference_genome)

        # write a file that indicates that this has finsihed
        open(rearranged_genome_finalFile, "w").write("finsihed")


def get_affected_positions_bed_row(r):

    """TAkes a bed row and returns the affected positions as a set of chrom_pos"""

    # get the positions as a series
    positions_range = range(int(r["start"]), int(r["end"]))
    index = range(0, len(positions_range))
    positions = "%s_"%r["chromosome"] + pd.Series(positions_range, index=index).apply(str)

    return set(positions)

def get_affected_positions_from_bed_df(bed_df):

    """Takes a bed df and returns a set with the affected positions"""

    return set.union(*bed_df.apply(get_affected_positions_bed_row, axis=1))

def get_affected_positions_from_bedpe_r(r, extra_bp=100):

    """Takes a bedpe df and gets the affected positions"""

    # get a bed_df
    bed_df = pd.DataFrame({"chromosome": [r["chrom1"], r["chrom2"]], 
                           "start": [r["start1"], r["start2"]], 
                           "end": [r["end1"], r["end2"]]})


    # add extra bp
    bed_df.start = bed_df.start - extra_bp
    bed_df.end = bed_df.end + extra_bp

    return get_affected_positions_from_bed_df(bed_df)

def get_length_bedpe_r(r):

    """Gets the len of the bedpe row"""

    if r["chrom1"]==r["chrom2"]: 
        length = int(r["medium2"]-r["medium1"])
        if length<=0: raise ValueError("len should be positive")
        return length

    else: return np.nan

def get_df_bedpe_with_nonOverlappingBreakpoints(df_bedpe):

    """Gets a df_bedpe that has affected_positions_arroundBp and returns the set of breakpoints that are not overlapping"""

    # check that the indices are unique
    if len(df_bedpe)!=len(set(df_bedpe.index)): raise ValueError("indices should be unique")

    # initialize the correct indices
    coorect_indices = set()

    # initialize the already affected positions
    already_affected_positions = set()

    # go through each breakpoint
    for I in df_bedpe.index:

        positons = df_bedpe.loc[I, "affected_positions_arroundBp"]

        # if the positions overlap with previous breakpoints, skip
        if positons.intersection(already_affected_positions)!=set(): continue

        #keep 
        already_affected_positions.update(positons)
        coorect_indices.add(I)

    df_bedpe = df_bedpe.loc[coorect_indices]

    return df_bedpe

def get_SVs_arround_breakpoints(genome_file, df_bedpe, nvars, outdir, svtypes, replace=False, max_ninsertions=100):

    """This function takes a df bedpe and defines svtypes arround these breakpoints. interchromosomal breakpoints will be equally distributed to """

    # get the dir
    make_folder(outdir)
    expected_files = {"%s/%s.tab"%(outdir, s) for s in svtypes}


    if any([file_is_empty(f) for f in expected_files]) or replace is True:
        print_if_verbose("calculating random variants among the provided breakpoints")

        # initialize a dict with each svtype and the current number of locations, a
        svtype_to_svDF = {svtype : pd.DataFrame(columns=[x for x in svtype_to_fieldsDict[svtype]["all_fields"] if x!="ID"]) for svtype in svtypes}

        # deifine the genome len
        chrom_to_len = get_chr_to_len(genome_file)

        # define the chroms
        chroms = set(chrom_to_len)

        # define the interesting bedpe, the one compromised by chroms and reshufle and reindex
        df_bedpe = df_bedpe[(df_bedpe.chrom1.isin(chroms)) & (df_bedpe.chrom2.isin(chroms))].sample(frac=1)
        df_bedpe.index = list(range(0, len(df_bedpe)))

        # set the max_n_breakpoints, 
        max_n_breakpoints = nvars*3000
        if len(df_bedpe)>max_n_breakpoints: 
        
            """
            # setting half of them to be interchromosomal
            half_max_n_breakpoints = int(max_n_breakpoints/2)

            df_bedpe["is_interchromosomal_bp"] = df_bedpe.chrom1!=df_bedpe.chrom2
            df_bedpe_same_chrom = df_bedpe[~(df_bedpe.is_interchromosomal_bp)].iloc[0:half_max_n_breakpoints]
            df_bedpe_dif_chrom =  df_bedpe[df_bedpe.is_interchromosomal_bp].iloc[0:half_max_n_breakpoints]

            df_bedpe = df_bedpe_dif_chrom.append(df_bedpe_same_chrom)
            df_bedpe.index = list(range(0, len(df_bedpe)))
            """
            
            df_bedpe = df_bedpe.iloc[0:max_n_breakpoints]


        if len(df_bedpe)>0:

            # define the interval to add for overlaps
            add_interval_bp = 100

            # add things
            df_bedpe["affected_positions_arroundBp"] = df_bedpe.apply(lambda r: get_affected_positions_from_bedpe_r(r, extra_bp=add_interval_bp), axis=1)
            df_bedpe["medium1"] = (df_bedpe.start1 + (df_bedpe.end1 - df_bedpe.start1)/2).apply(int)
            df_bedpe["medium2"] = (df_bedpe.start2 + (df_bedpe.end2 - df_bedpe.start2)/2).apply(int)
            df_bedpe["length"] = df_bedpe.apply(get_length_bedpe_r, axis=1)

            # get breakpoints that are above 50 bp
            df_bedpe = df_bedpe[(df_bedpe.length>=50) | (pd.isna(df_bedpe.length))]

            # keep breakpoints that are non overlapping by add_interval_bp
            print_if_verbose("getting non redundant breakpoints out of %i of them"%len(df_bedpe))
            df_bedpe = get_df_bedpe_with_nonOverlappingBreakpoints(df_bedpe)

            # define the legths of the variants
            lengths_SVs = list(df_bedpe[(~pd.isna(df_bedpe.length)) & (df_bedpe.length>=50)].length)
            random.shuffle(lengths_SVs)

            # if empty, pick a fixed lenght
            if len(lengths_SVs)==0: lengths_SVs = [1000]

            # already_affected_positions contains a set of "chrom_pos" with the already affected positions. It will be updated with each new variant
            already_affected_positions = set()

            ############ DEFINE THE RANDOM SVs INVOLVING df_bedpe ############

            # initialize a set of the svtypes that are already full, meaning that they already have nvars
            already_nvars_svtypes = set()

            # initialize a set that will contain the bad bedpe rows
            bpID_to_tried_svtypes = {bpID : set() for bpID in df_bedpe.index}

            # define the number of breakpoints traversed
            original_n_breakpoints = cp.deepcopy(len(df_bedpe))

            # sort the bedpe so that interchromosomal events are first
            df_bedpe["order_by_chrom"] =  df_bedpe.apply(lambda r: {True:10, False:1}[r["chrom1"]==r["chrom2"]], axis=1)

            print("There are %i/%i interchromosomal breakpoints"%(sum(df_bedpe.order_by_chrom==1), len(df_bedpe)))

            # go through each breakpoint and assign it to a cahegory. Break if a
            while len(df_bedpe)>0:
                print_if_verbose("already traversed %i/%i breakpoints. There are %i affected positions"%(original_n_breakpoints - len(df_bedpe), original_n_breakpoints, len(already_affected_positions)))

                # get the sorted svtypes
                sorted_svtypes = [x for x in ["translocations", "tandemDuplications", "deletions", "inversions", "insertions"] if x in svtypes]

                # print the current numbers
                for svtype in sorted_svtypes: print_if_verbose("%i/%i %s defined"%(len(svtype_to_svDF[svtype]), nvars, svtype))

                # if you already found all the interchromosomal events, keep only df_breakpoints among the same chromosome
                if already_nvars_svtypes==(svtypes.intersection({"translocations", "insertions"})):
                    df_bedpe = df_bedpe[df_bedpe.order_by_chrom==10]

                # if you already found all the intrachromosomal events and there are no interchromosomal breakpoints, skip
                if already_nvars_svtypes==(svtypes.intersection({"inversions", "tandemDuplications", "deletions"})) and all(df_bedpe.order_by_chrom==10): break

                # interate through each svtype
                for svtype in sorted_svtypes:  

                    # get only bedpe rows that are not overlapping the current positions, in terms of breakpoints
                    df_bedpe = df_bedpe[(df_bedpe.affected_positions_arroundBp.apply(lambda positions: positions.intersection(already_affected_positions)).apply(len))==0]

                    # get only bedpe rows that have not already been tried to assign to all svtypes
                    good_bpIDs = {bpID for bpID, tried_svtypes in bpID_to_tried_svtypes.items() if tried_svtypes!=svtypes}.intersection(set(df_bedpe.index))
                    df_bedpe = df_bedpe.loc[good_bpIDs]

                    # if empty, continue
                    if len(df_bedpe)==0: continue

                    # get the first bedpe row. Sorting in a way that the interchromosomal events will happen always first. This is to prioritize translocations
                    df_bedpe = df_bedpe.sort_values(by="order_by_chrom")
                    r = df_bedpe.iloc[0]

                    # record that this was already tried
                    bpID_to_tried_svtypes[r.name].add(svtype)

                    # if there are already nvars, skip
                    svDF = svtype_to_svDF[svtype]
                    if len(svDF)>=nvars:
                        already_nvars_svtypes.add(svtype)
                        continue

                    # assign the breakpoint to the df_sv_r
                    if svtype in {"inversions", "tandemDuplications", "deletions", "translocations"}:

                        # it should be in the same chromosome unless they are translocations
                        if r["chrom1"]!=r["chrom2"] and svtype!="translocations": continue
                        if r["chrom1"]==r["chrom2"] and svtype=="translocations": continue

                        # define the svdf depending on the svtype
                        if svtype in {"inversions", "tandemDuplications", "deletions"}:

                            # define the resulting df_sv_r
                            df_sv_r = pd.DataFrame({"Chr":[r["chrom1"]], "Start":[int(r["medium1"])+1], "End":[int(r["medium2"])]})

                            # define the chroms
                            chroms_calculating_positions =  {r["chrom1"]}

                        # for translocations
                        else:

                            # define the coordinates. These are always the same:
                            chrA = r["chrom1"]
                            startA = 1
                            endA = int(r["medium1"])
                            chrB = r["chrom2"]

                            # the orientation determines chrB start and end
                            orientation = {True: "5_to_5", False: "5_to_3"}[(random.random()>=0.5)]
                            if orientation=="5_to_5": 
                                startB = 1
                                endB = int(r["medium2"])

                            elif orientation=="5_to_3":
                                startB = int(r["medium2"])
                                endB = chrom_to_len[chrB]

                            else: raise ValueError("%s is not valid"%orientation) 
                           
                            # don't consider translocations at the border of the chromosomes
                            if endA<50 or r["medium2"]<50 or (chrom_to_len[chrA]-endA)<50 or (chrom_to_len[chrB]-r["medium2"])<50: continue


                            # define the resulting df_sv_r
                            df_sv_r = pd.DataFrame({"ChrA":[chrA], "StartA":[startA], "EndA":[endA], "ChrB":[chrB], "StartB":[startB], "EndB":[endB], "Balanced":[True]})

                            # define the chroms
                            chroms_calculating_positions =  {chrA, chrB}

                        # get the affected positions of this svtype
                        affected_positions = get_affected_positions_from_bed_df(get_affected_region_bed_for_SVdf(df_sv_r, svtype, chroms_calculating_positions, add_interval_bp=add_interval_bp, first_position_idx=1, translocations_type="breakpoint_pos", chr_to_len=chrom_to_len, insertions_type="only_one_chrB_breakpoint")[0])
                        
                        # only keep if the already affected positions do not overlap this SV
                        if affected_positions.intersection(already_affected_positions)!=set(): continue

                        # keep
                        already_affected_positions.update(affected_positions)
                        svtype_to_svDF[svtype] = svDF.append(df_sv_r)

                    elif svtype=="insertions":  

                        # go through different lengths. Until all of them or 
                        for Iins in range(0, max_ninsertions):

                            # define randomly a length from the lengths_SVs
                            length_ins = random.choice(lengths_SVs)

                            # the insertion will be from medium1 - length_ins ----> medium2. The Copied will be a 50% possibility
                            startA = int(r["medium1"] - length_ins)
                            endA = int(r["medium1"])
                            startB = int(r["medium2"])
                            endB = int(r["medium2"] + length_ins)
                            copied = (random.random()>=0.5)

                            # don't consider insertions at the border of the chromosome
                            if startA<50 or (chrom_to_len[r["chrom2"]]-endB)<50: continue

                            # define the resulting df_sv_r
                            df_sv_r = pd.DataFrame({"ChrA":[r["chrom1"]], "StartA":[startA], "EndA":[endA], "ChrB":[r["chrom2"]], "StartB":[startB], "EndB":[endB], "Copied":[copied]})

                            # get the affected positions of this svtype
                            affected_positions = get_affected_positions_from_bed_df(get_affected_region_bed_for_SVdf(df_sv_r, svtype, {r["chrom1"], r["chrom2"]}, add_interval_bp=add_interval_bp, first_position_idx=1, translocations_type="breakpoint_pos", chr_to_len=chrom_to_len, insertions_type="only_one_chrB_breakpoint")[0])
                        
                            # only keep if the already affected positions do not overlap this SV
                            if affected_positions.intersection(already_affected_positions)!=set(): continue

                            # keep
                            already_affected_positions.update(affected_positions)
                            svtype_to_svDF[svtype] = svDF.append(df_sv_r)

                            # break if you could get an insertion
                            break 

                    else: raise ValueError("%s is not valid"%svtype)

                # if all the svs have been created, exit
                if already_nvars_svtypes==svtypes:
                    print_if_verbose("all svtypes found")
                    break

        ##################################################################

        # add IDs and write
        for svtype, svDF in svtype_to_svDF.items():
            svDF["index_val"] = list(range(0, len(svDF))) 
            svDF["ID"] = svtype + "_" + svDF.index_val.apply(str)
            svDF["Name"] = svDF.ID
            svDF[svtype_to_fieldsDict[svtype]["all_fields"] + ["Name"]].to_csv("%s/%s.tab"%(outdir, svtype), sep="\t", index=False, header=True)

def get_subset_genome_as_fasta(input_genome, output_genome, chroms, replace=False):

    """Takes an input genome and keeps only the chroms, writing to output_genome"""

    if file_is_empty(output_genome) or replace is True:

        # define all seqs
        all_chroms = [seq for seq in SeqIO.parse(input_genome, "fasta") if seq.id in chroms]

        # write
        output_genome_tmp = "%s.tmp"%output_genome
        SeqIO.write(all_chroms, output_genome_tmp, "fasta")
        os.rename(output_genome_tmp, output_genome)


def simulate_SVs_in_genome(reference_genome, mitochondrial_chromosome, outdir, nvars=200, replace=False, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications", "translocations"}, bedpe_breakpoints=None):

    """This function generates nvars into the reference genome splitting by gDNA and mtDNA with files written under outdir"""

    print_if_verbose("generating simulations")

    # initialize a df that will contain the randomly-simulated vars
    final_svtype_to_svDF = {svtype : pd.DataFrame() for svtype in svtypes}

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

        # simulate random SVs into regions without previous SVs 
        random_sim_dir = "%s/random_SVs"%genome_outdir

        # define the chromosomes above window length
        chroms_above_window_l = {chrom for chrom in chroms if chrom_to_len[chrom]>=window_l}
        if len(chroms_above_window_l)==0: chroms_above_window_l = chroms

        # define a bed file with all the regions
        if bedpe_breakpoints is None:
            print_if_verbose("simulating randomly distributed SVs")
         
            # get a genome with the desired chromosomes
            genome_file_correctChroms = "%s.correct_chromosomes.fasta"%genome_file
            get_subset_genome_as_fasta(genome_file, genome_file_correctChroms, chroms_above_window_l, replace=replace)

            # define the length of the shortest chromosome
            len_shortest_chr = min([chrom_to_len[c] for c in chroms_above_window_l])

            #### GET THE RANDOM INS,INV,DEL,TRA ####

            # initialize the cmd
            randomSV_cmd = "%s --input_genome %s --outdir %s --len_shortest_chr %i"%(create_random_simulatedSVgenome_R, genome_file_correctChroms, random_sim_dir,  len_shortest_chr) 

            # add the numbers of SVs depending on if it is random or not
            expected_svtypes = {"insertions", "deletions", "inversions", "tandemDuplications", "translocations", "translocations"}.intersection(svtypes)

            # if there is only one, get the expecyed SVtypes
            if len(chroms_above_window_l)==1: expected_svtypes = {s for s in expected_svtypes if s!="translocations"}

            # define the expected_files
            expected_files = {"%s/%s.tab"%(random_sim_dir, svtype) for svtype in expected_svtypes}

            # add to cmd        
            svtype_to_arg = {"insertions":"number_Ins", "deletions":"number_Del", "inversions":"number_Inv", "tandemDuplications":"number_Dup", "translocations":"number_Tra"}
            for svtype in expected_svtypes: randomSV_cmd += " --%s %i"%(svtype_to_arg[svtype], vars_to_simulate)

            # define the finalisation file
            random_sims_performed_file = "%s/random_sims_performed.txt"%random_sim_dir

            # remove the file if all the expected svtypes are empty
            if any([file_is_empty(f) for f in expected_files]):
                print_if_verbose("none of the files were performed")
                remove_file(random_sims_performed_file)

            # run the cmd if necessary
            if file_is_empty(random_sims_performed_file) or replace is True:
                print_if_verbose("generating random SVs")

                # make and delete the folder
                delete_folder(random_sim_dir); make_folder(random_sim_dir)

                # run the random simulation
                std_rearranging_genome = "%s/simulation_std.txt"%random_sim_dir
                #std_rearranging_genome = "stdout"
                print_if_verbose("getting random SVs. The std is in %s"%std_rearranging_genome)
                if std_rearranging_genome!="stdout": run_cmd("%s > %s 2>&1"%(randomSV_cmd, std_rearranging_genome), env=EnvName_R)
                else: run_cmd(randomSV_cmd, env=EnvName_R)

                # check that everything went fine
                if any([file_is_empty(f) for f in expected_files]): raise ValueError("none of the SVs were simulated")

                remove_file(std_rearranging_genome)

                # edit the insertions 
                insertions_file = "%s/insertions.tab"%random_sim_dir
                if not file_is_empty(insertions_file): rewrite_insertions_uniformizedFormat_simulateSV(insertions_file)

                open(random_sims_performed_file, "w").write("random simulations performed")          

            ########################################

        ##### GET THE RANDOM TAN,DEL,INV,INS,TRA ARROUND PROVIDED BREAKPOINTS #####
        else:

            # load the bedpe
            bedpe_fields = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]
            df_bedpe = pd.read_csv(bedpe_breakpoints, sep="\t", header=-1, names=bedpe_fields)

            # define the translocations 
            get_SVs_arround_breakpoints(genome_file, df_bedpe, vars_to_simulate, random_sim_dir, svtypes, replace=replace)

        #######################################################################

        ######### KEEP VARS #########
       
        for svtype in final_svtype_to_svDF.keys():
            svDF = final_svtype_to_svDF[svtype]

            # skip translocations if they where not feasible
            if svtype=="translocations" and file_is_empty("%s/%s.tab"%(random_sim_dir, svtype)): continue

            # get the new sv
            file = "%s/%s.tab"%(random_sim_dir, svtype)
            new_svDF = pd.read_csv(file, sep="\t")
            new_svDF = new_svDF[[c for c in new_svDF.keys() if "BpSeq" not in c]]

            # add the name
            new_svDF["ID"] = new_svDF.Name + "_sim_%s"%type_genome

            # append 
            final_svtype_to_svDF[svtype] = svDF.append(new_svDF, sort=True)

        ##############################

    ##### REARRANGE THE GENOME WITH THE CALCULATED VARS #####

    # define the final outdirs
    final_simulated_SVs_dir = "%s/final_simulated_SVs"%(outdir); 
    final_rearranged_genome = "%s/rearranged_genome.fasta"%final_simulated_SVs_dir
    final_rearranged_genome_finalFile = "%s.performed"%(final_rearranged_genome)

    # insert the predicted  
    generate_rearranged_genome_from_svtype_to_svDF(reference_genome, final_svtype_to_svDF, final_simulated_SVs_dir, replace=replace)

    # get the final svtype
    final_svtype_to_svfile = {svtype : "%s/%s.tab"%(final_simulated_SVs_dir, svtype) for svtype in svtypes}

    return final_svtype_to_svfile, final_rearranged_genome

    ############################################################

def get_mosdepth_coverage_per_windows_output_likeBamStats(fileprefix, sorted_bam, windows_bed, replace=False, extra_threads=0, chromosome_id=""):

    """This function uses mosdepth to get the coverage for some regions in bed for a sorted_bam """

    # get outfiles
    fileprefix_tmp = "%s.tmp"%fileprefix
    regions_file = "%s.regions.bed.gz"%fileprefix 
    regions_file_tmp = "%s.regions.bed.gz"%fileprefix_tmp 
    thresholds_file = "%s.thresholds.bed.gz"%fileprefix 
    thresholds_file_tmp = "%s.thresholds.bed.gz"%fileprefix_tmp

    if file_is_empty(regions_file) or file_is_empty(thresholds_file) or replace is True:
        
        # change the end, setting it to +1, and also sorting
        windows_1_based = "%s.1_based.bed"%windows_bed
        windows_1_based_stderr = "%s.generating.stderr"%windows_1_based
        #print_if_verbose("getting 1-based bed file. The stderr is in %s"%windows_1_based_stderr)
        run_cmd(""" awk '{print $1 "\t" ($2+1) "\t" ($3)}' %s | sort -k1,1 -k2,2n > %s 2>%s"""%(windows_bed, windows_1_based, windows_1_based_stderr))
        remove_file(windows_1_based_stderr)

        # get the cmd
        mosdepth_std = "%s.generating.std"%fileprefix_tmp
        #print_if_verbose("running mosdepth. The std is in %s "%mosdepth_std)
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


def get_file_size(file, factor_size=1e9):

    """Gets the size of a file in units"""
    
    return (os.stat(file).st_size)/factor_size

def generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file="none", replace=False, run_in_parallel=True, delete_bams=True, threads=4):

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
        df_windows = pd.read_csv(windows_file, sep="\t", header=None, names=["chromosome", "start", "end"]).sort_values(by=["chromosome", "start", "end"])

        # check that everything is OK
        if any(df_windows.start>=df_windows.end): 

            print(df_windows[(df_windows.start>=df_windows.end)])
            raise ValueError("There can't be windows with start>=end")

        # write it
        df_windows.to_csv(windows_file, header=False, index=False, sep="\t")

        # print the number of windows
        print_if_verbose("calculating coverage for %i windows"%len(df_windows))

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

        start_time =  time.time()

        # define the size of the bam and print it
        size_bam = get_file_size
        print_if_verbose("calculating coverage for a bam of %.3f Gb"%(get_file_size(sorted_bam)))

        # define a prefix
        mosdepth_outprefix = "%s.generating"%coverage_file

        # run mosdepth
        all_df =  get_mosdepth_coverage_per_windows_output_likeBamStats(mosdepth_outprefix, sorted_bam, windows_file, replace=replace, extra_threads=threads, chromosome_id="")

        # get only the important fields
        bamstats_fields = ["#chrom", "start", "end", "length", "mediancov_1", "nocoveragebp_1", "percentcovered_1"]
        all_df[bamstats_fields]

        # check that it is not empty
        if len(all_df)==0: raise ValueError("There is no proper coverage calculation for %s on windows %s"%(sorted_bam, windows_file))

        # write
        coverage_file_tmp = "%s.tmp"%coverage_file
        all_df.to_csv(coverage_file_tmp, sep="\t", header=True, index=False)

        nseconds = time.time() - start_time
        print_if_verbose("the running of mosdepth for %i windows took %s seconds"%(len(all_df), nseconds))

        # rename
        os.rename(coverage_file_tmp, coverage_file)

    return coverage_file

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
                allocated_ram = get_availableGbRAM(gridss_tmpdir)*fractionRAM_to_dedicate
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
                run_cmd(gridss_cmd, env=EnvName_gridss)

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
    if (file_is_empty(gridss_VCFoutput_with_simple_event) or replace is True) and n_breakends>0 : run_cmd("%s %s > %s 2>&1"%(annotate_simpleEvents_gridssVCF_R, gridss_VCFoutput, simple_event_std), env=EnvName_R)

    remove_file(simple_event_std)

    return gridss_VCFoutput_with_simple_event

def load_single_sample_VCF(path):

    """Loads a vcf with a single sample into a df."""

    vcf_df_file = "%s.vcf_df.py"%path

    if file_is_empty(vcf_df_file):

        print_if_verbose("running load_single_sample_VCF")

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

        print_if_verbose("getting INFO as dict")

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

        print_if_verbose("load_single_sample_VCF ran")

        # save
        save_object(df, vcf_df_file)

    # load
    df = load_object(vcf_df_file)

    return df

def getNaN_to_0(x):

    if pd.isna(x): return 0.0
    else: return x

def add_info_to_gridssDF(df, reference_genome, expected_fields={"allele_frequency", "allele_frequency_SmallEvent", "other_coordinates", "other_chromosome", "other_position", "other_orientation",  "inserted_sequence", "len_inserted_sequence", "length_event", "has_poly16GC", "length_inexactHomology", "length_microHomology", "overlaps_repeats", "eventID_as_clove", "coordinates"}, median_insert_size=500, median_insert_size_sd=50):

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

        df["coordinates"] = df["#CHROM"] + ":" + df.POS.apply(int).apply(str)

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
            else: return int(IHOMPOS[1]) - int(IHOMPOS[0])

        df["length_inexactHomology"] = df.INFO_IHOMPOS.apply(add_inexact_homology_length)

        # add the length of homology
        def get_homology_length(HOMLEN):
            
            if pd.isna(HOMLEN): return 0
            else: return int(HOMLEN)

        if "INFO_HOMLEN" in df.keys(): df["length_microHomology"] = df.INFO_HOMLEN.apply(get_homology_length)
        else: df["length_microHomology"] = 0

        # add the actual allele_frequency, which depends on the median_insert_size and the median_insert_size_sd. If the breakend is longer than the insert size it is a large event
        maxiumum_insert_size = median_insert_size + median_insert_size_sd
        def get_allele_freq(r):
            
            if r["length_event"]>maxiumum_insert_size: return r["allele_frequency"]
            else: return r["allele_frequency_SmallEvent"]
        
        df["real_AF"] = df.apply(get_allele_freq, axis=1)

        if any(pd.isna(df["real_AF"])): raise ValueError("There are NaNs in the real_AF field")

        # get the eventID
        df["eventID_as_clove"] = df.INFO_EVENT + "o"


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
                run_cmd("%s %s > %s 2>&1"%(analyze_svVCF_simple, svVCF, r_stdout), env=EnvName_R)                

            else:
                
                print_if_verbose("Getting files for svVCF file. Mainly generating a bedpe file for breakpoints with some extra info. The std can be found in %s"%r_stdout)
                run_cmd("%s %s > %s 2>&1"%(analyze_svVCF, svVCF, r_stdout), env=EnvName_R)

            remove_file(r_stdout)

        else: open(bedpe_file, "w").write("no_vcf_records\n")

    return bedpe_file

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
  
def generate_nt_content_file(genome, target_nts="GC", replace=False):

    """Takes a genome and outputs a file with chromosome, position and 1 or 0 regarding if any of the target_nts is the same in the genome. This is 0-based"""

    target_nt_content_file = "%s.%scontent.tab"%(genome, target_nts)

    if file_is_empty(target_nt_content_file) or replace is True:
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


def get_distanceToTelomere_chromosome_GCcontent_to_coverage_fn(df_coverage_train, genome, outdir, expected_coverage_per_bp, mitochondrial_chromosome="mito_C_glabrata_CBS138", replace=False):

    """This function takes a training df_coverage (with windows of a genome) and returns a lambda function that takes GC content, chromosome and  distance to the telomere and returns coverage according to the model.

    Your genome graph can also be a linear genome, you just have to create it without considering breakpoints"""
    print_if_verbose("getting coverage-predictor function")

    # rename the training df
    df = df_coverage_train.rename(columns={"#chrom":"chromosome", "mediancov_1":"coverage"})

    # add the distance to the telomere
    chr_to_len = get_chr_to_len(genome)
    df["middle_position"] = (df.start + (df.end - df.start)/2).apply(int)
    df["distance_to_telomere"] = df.apply(lambda r: min([r["middle_position"], chr_to_len[r["chromosome"]]-r["middle_position"]]), axis=1)

    # add the gc content
    gcontent_outfile = "%s/GCcontent.py"%outdir
    df = get_df_with_GCcontent(df, genome, gcontent_outfile, replace=replace)

    # define the set of each type of chromosomes
    all_chromosomes = set(df.chromosome)
    if mitochondrial_chromosome!="no_mitochondria": mtDNA_chromosomes = set(mitochondrial_chromosome.split(","))
    else: mtDNA_chromosomes = set()
    gDNA_chromosomes = all_chromosomes.difference(mtDNA_chromosomes)

    # load the genome
    chr_to_len = get_chr_to_len(genome)
    print_if_verbose("there are %i/%i chromsomes above window_l"%(sum([l>=window_l for l in chr_to_len.values()]), len(chr_to_len)))

    ######## find the coeficients for each chromosome #########

    # map each chromosome to the coefs of the quadratic fit that explains coverage form the distance to the telomere and also the coefs of the GC content explaining the resiudal of this fit
    chrom_to_coefType_to_coefs = {}

    # go through each type of genome
    for type_genome, chroms in [("mtDNA", mtDNA_chromosomes), ("gDNA", gDNA_chromosomes)]:
        print_if_verbose("investigating %s"%type_genome)

        # define the training df
        df_g = df[df.chromosome.isin(chroms)]

        # define the relative coverage of each window of this genome
        median_coverage = np.median(df_g[df_g.coverage>0].coverage)
        df_g["relative_coverage"] = df_g.coverage / median_coverage

        # go through each chrom and identify that it is duplicated if the quadratic fit from the prediction of the distance to the telomere suggests a minimum of >=1.6. Also discard chromosomes where the size is below window_l
        bad_chroms = set()
        duplicated_chroms = set()
        for chrom in chroms:

            # get df of this chrom
            df_c = df_g[df_g.chromosome==chrom]

            # flag dup chromosomes
            if np.median(df_c.relative_coverage)>=1.7: duplicated_chroms.add(chrom)

            # flag short chromosomes
            if chr_to_len[chrom]<window_l: bad_chroms.add(chrom)

        # define the training set for the modelling
        df_correct = df_g[(df_g.relative_coverage<=5) & (df_g.relative_coverage>0.05) & ~(df_g.chromosome.isin(duplicated_chroms)) & ~(df_g.chromosome.isin(bad_chroms))]

        # if the filtering is useless, use all the df. This is a way to skip the modelling.
        if len(df_correct)==0: df_correct = df_g[~(df_g.chromosome.isin(bad_chroms))]

        # if still empty, use all
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


        ### calculate the r2 based on all windows ###

        # get r2
        r2 = r2_score(df_correct.coverage, df_correct.coverage_from_dist_to_telomere_and_GC_content)
        fraction_windows_considered = len(df_correct)/len(df_g)

        # get the pval of the coverage not being a normal distriburion
        if len(df_g)>10:

            pvalNormDist_rel_coverage = get_pvalue_is_normal_distribution(df_g.relative_coverage)
            print(pvalNormDist_rel_coverage)

        else: pvalNormDist_rel_coverage = 1

        print_if_verbose("The rsquare for %s is %.3f. %.3f pct of the windows are included. The p of not being a normal distribution is %.6f"%(type_genome, r2, fraction_windows_considered*100, pvalNormDist_rel_coverage))

        # change coefs so that there is no modelling if the fit is bad
        if r2<0.2 or pd.isna(r2) or fraction_windows_considered<0.5 or pvalNormDist_rel_coverage>0.05:
            print_if_verbose("not modelling coverage per window")
            coefs_dist_to_telomere = [expected_coverage_per_bp, 0, 0]
            coefs_GCcontent = [0, 0, 0]

        # re-calculate
        df_g["coverage_from_dist_to_telomere"] = poly.polyval(df_g.distance_to_telomere, coefs_dist_to_telomere)
        df_g["residualCoverage_from_dist_to_telomere"] = df_g.coverage - df_g.coverage_from_dist_to_telomere
        df_g["residualCoverage_from_dist_to_telomere_from_GC_content"] = poly.polyval(df_g.GCcontent, coefs_GCcontent)
        df_g["coverage_from_dist_to_telomere_and_GC_content"] = df_g["coverage_from_dist_to_telomere"] + df_g["residualCoverage_from_dist_to_telomere_from_GC_content"]

        #############################################

        # save the coefficients
        for chrom in chroms: chrom_to_coefType_to_coefs[chrom] = {"dist_telomere":coefs_dist_to_telomere, "GCcontent":coefs_GCcontent}

        # plot
        outfile = "%s/coverage_modelling_%s.pdf"%(outdir, type_genome)

        if file_is_empty(outfile) or replace is True:

            # define the chroms to plot
            if len(chroms)<10: chroms_plot = sorted(chroms)
            else: chroms_plot = sorted(chroms.difference(bad_chroms))
            print_if_verbose("plotting coverage modelling for %i chroms"%len(chroms_plot))

            # plot the coverage for each of the chromosomes
            fig = plt.figure(figsize=(7, len(chroms_plot)*5))

            for I, chrom in enumerate(chroms_plot):

                # initialize a subplot, where each row is one chromosome
                ax = plt.subplot(len(chroms_plot), 1, I+1)

                # get df of this chrom
                df_c = df_g[df_g.chromosome==chrom]

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

    ###############################################################

    # define the function that takes a tuple of (distToTelomere, chromosome and GCcontent) and returns the predicted relative coverage
    final_function = (lambda dist_telomere, chrom, GCcontent:  # this is suposed to be the tuple

                        (poly.polyval([dist_telomere], chrom_to_coefType_to_coefs[chrom]["dist_telomere"]) + # from the dist to tel
                        poly.polyval([GCcontent], chrom_to_coefType_to_coefs[chrom]["GCcontent"]))[0] # residual predicted from GC

                     )

    # check that it works
    df_g["cov_predicted_from_final_lambda"] = df_g.apply(lambda r: final_function(r["distance_to_telomere"], r["chromosome"], r["GCcontent"]), axis=1)

    if any(((df_g["coverage_from_dist_to_telomere_and_GC_content"]-df_g["cov_predicted_from_final_lambda"]).apply(abs))>0.01): raise ValueError("error in lambda function generation for coverage")

      
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

    # change the ID so that it ends always with an o
    df["ID"] = df.ID.apply(lambda x: "+".join([y[0:-1]+"o" for y in re.split("\+|\-", x)]))

    return df

def dummy_fun(x): return x

def parallelization_is_possible(threads):

    """Tests whether a multiprocessing can be launched, and returns whtehther it is possible"""

    try: 

        with multiproc.Pool(threads) as pool:

            # run a dummy function
            pool.starmap(dummy_fun, [(1,), (2,), (3,)])

            # close the pool
            pool.close()

        return True

    except AssertionError: return False


def get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, windows_file, replace=False, run_in_parallel=True, delete_bams=False, threads=4):

    """This function takes a windows file and a bam, and it runs generate_coverage_per_window_file_parallel but only for regions that are not previously calculated"""

    # check if it can be run in parallel
    if parallelization_is_possible(threads) is False: run_in_parallel = False
    print_if_verbose("running get_coverage_per_window_df_without_repeating with %i threads"%threads)

    # define the query_windows
    if open(windows_file).readlines()[0].startswith("chromosome"): query_windows_df = pd.read_csv(windows_file, sep="\t")
    else: query_windows_df = pd.read_csv(windows_file, sep="\t", header=None, names=["chromosome", "start", "end"])

    query_windows_df = query_windows_df.set_index(["chromosome", "start", "end"], drop=False)
    if len(query_windows_df)==0: return pd.DataFrame()
    
    # chek the initial length
    query_df_len = len(query_windows_df)

    # define the file were the coverage will be calculated
    calculated_coverage_file = "%s.coverage_per_window.tab"%sorted_bam

    ######### define previously measured regions ##########
    print_if_verbose("getting the previously-calculated regions")

    if not file_is_empty(calculated_coverage_file):

        # get into a df all the previously calculated coverages
        df_previosuly_calculated_coverages = pd.read_csv(calculated_coverage_file, sep="\t").set_index(["chromosome", "start", "end"], drop=False)

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
        print_if_verbose("calculating coverage for uncalculated windows")

        # remove the already calculated windows
        remove_file(calculated_coverage_file)

        # write a bed with the regions to measure
        bed_windows_to_measure = "%s.measuring.bed"%(calculated_coverage_file)
        windows_to_measure_df[["chromosome", "start", "end"]].to_csv(bed_windows_to_measure, sep="\t", header=False, index=False)

        # get the coverage df
        destination_dir = "%s.calculating_windowcoverage"%sorted_bam
        coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file=bed_windows_to_measure, replace=replace, run_in_parallel=run_in_parallel, delete_bams=delete_bams, threads=threads), sep="\t").rename(columns={"#chrom":"chromosome"}).set_index(["chromosome", "start", "end"], drop=False)

        delete_folder(destination_dir)
        remove_file(bed_windows_to_measure)

        delete_folder("%s.coverage_measurement_destination"%bed_windows_to_measure)
        remove_file("%s.coverage_provided_windows.tab"%bed_windows_to_measure)

        # get the merged dfs
        df_coverage_all = df_previosuly_calculated_coverages.append(coverage_df, sort=True).loc[list(query_windows_df.index)]
        
        # save
        calculated_coverage_file_tmp = "%s.tmp"%calculated_coverage_file
        df_coverage_all.to_csv(calculated_coverage_file_tmp, sep="\t", header=True, index=False)
        os.rename(calculated_coverage_file_tmp, calculated_coverage_file)

    #############################################

    # load the df
    df_coverage_all = pd.read_csv(calculated_coverage_file, sep="\t")

    # drop duplicates
    df_coverage_all = df_coverage_all.drop_duplicates(subset=["chromosome", "start", "end"], keep='first')

    # set the index to be the chromosome,start,end
    df_coverage_all = df_coverage_all.set_index(["chromosome", "start", "end"], drop=False)

    # keep as df_coverage_final those that intersect with query_windows_df
    df_coverage_final = df_coverage_all.loc[query_windows_df.index]

    # get the correct index
    df_coverage_final.index = list(range(len(df_coverage_final)))

    # debug
    if query_df_len!=len(df_coverage_final): raise ValueError("the length has changed in the process")

    return df_coverage_final

def get_int(x):

    # def get nans to -1

    try: return int(x)
    except: return -1

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
        print(r, region_name, maxPos, start, end)
        raise ValueError("start after end")

    # return a series of all important fields
    return pd.Series({"chromosome":r["chromosome"], "start":start, "end":end, "region_name":region_name})


def get_chrom_to_bpPositions(df_clove, reference_genome):

    """This function takes a df_clove and maps each chromosome to the breakpoint positions"""

    # all chroms
    all_chromosomes = {seq.id for seq in SeqIO.parse(reference_genome, "fasta")}

    # change to ints
    for f in ["POS", "START", "END"]: df_clove[f] = df_clove[f].apply(int)

    # get each of the chroms
    chrom_to_bpPositions_CHROM = dict(df_clove.groupby("#CHROM").apply(lambda df_c: get_uniqueVals_df(df_c[["POS"]])))
    chrom_to_bpPositions_CHR2 = dict(df_clove.groupby("CHR2").apply(lambda df_c: get_uniqueVals_df(df_c[["START", "END"]]).difference({-1})))

    chrom_to_bpPositions = {chrom:set() for chrom in all_chromosomes}
    for chromDict in [chrom_to_bpPositions_CHROM, chrom_to_bpPositions_CHR2]:
        for chrom, bpPositions in chromDict.items():
            chrom_to_bpPositions[chrom].update(bpPositions)

    return chrom_to_bpPositions


def get_IDwindow_df(r):

    """Takes a row of a df and returns the chromosome_start_end"""

    return "%s_%i_%i"%(r["chromosome"], r["start"], r["end"])

def get_df_with_coverage_per_windows_relative_to_neighbor_regions(df_windows, bed_windows_prefix, reference_genome, sorted_bam, df_clove, median_coverage, replace=True, run_in_parallel=True, delete_bams=True, threads=4):

    """Takes a df with windows of the genome and returns it with the coverage and the relative to the genome. It returns a df with several relative coverage measures."""

    print_if_verbose("getting coverage relative to neighbors")

    # get the initial index
    initial_index = list(df_windows.index)

    # keep
    df_windows = cp.deepcopy(df_windows)

    # get the coverage to len
    chrom_to_maxPos = {seq.id : len(seq.seq)-1 for seq in SeqIO.parse(reference_genome, "fasta")}

    # map the chromosome to the positions with breakpoints
    chrom_to_bpPositions = get_chrom_to_bpPositions(df_clove, reference_genome)

    # initialize a df windows for the 5' and 3' regions
    df_windows = df_windows.sort_values(by=["chromosome", "start", "end"]).drop_duplicates(subset=["chromosome", "start", "end"])

    # get a df with all windows
    all_df_windows = cp.deepcopy(df_windows)
    all_df_windows.index = all_df_windows.apply(get_IDwindow_df, axis=1) + "_originalRegion"
    all_df_windows["IDwindow"] = all_df_windows.index

    # original_indices
    all_df_windows_index = list(all_df_windows.index)

    # go through each region and get a coverage df
    for region in ["target", "5", "3"]: 

        if region=="target":
            df_region = df_windows
            df_region["region_name"] = "target_region"

        else:

            # get a df with the regions
            df_region = df_windows.apply(lambda r: get_target_region_row(r, region, chrom_to_bpPositions[r["chromosome"]], chrom_to_maxPos[r["chromosome"]]), axis=1)

        # add the index
        df_region.index = df_region.apply(get_IDwindow_df, axis=1)

        # get the coverage df
        bed_file = "%s.%s.bed"%(bed_windows_prefix, region)
        df_region.to_csv(bed_file, sep="\t", header=True, index=False)

        coverage_df = get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, bed_file, replace=replace, run_in_parallel=run_in_parallel, delete_bams=delete_bams, threads=threads).drop_duplicates(subset=["chromosome", "start", "end"])
        coverage_df.index = coverage_df.apply(get_IDwindow_df, axis=1)

        # make sure that they are unique
        if len(coverage_df)!=len(set(coverage_df.index)): raise ValueError("coverage_df is not unique")

        # add the coverage to the windows df
        df_region["IDwindow"] = df_region.index
        df_region["coverage"] = df_region.IDwindow.apply(lambda x: coverage_df.loc[x, "mediancov_1"])
        if any(pd.isna(df_region.coverage)): raise ValueError("there should be no NaNs")

        # add the original index
        df_region.index = all_df_windows_index

        # add to all df
        all_df_windows["%s_coverage"%region] = all_df_windows.IDwindow.apply(lambda x: df_region.loc[x, "coverage"])

        # add relative parms
        all_df_windows["relative_coverage_%s"%region] = all_df_windows["%s_coverage"%region]/median_coverage

        # add rge coordinates
        all_df_windows["%s_region_start"%region] = all_df_windows.IDwindow.apply(lambda x: df_region.loc[x, "start"])
        all_df_windows["%s_region_end"%region] = all_df_windows.IDwindow.apply(lambda x: df_region.loc[x, "end"])

    # get the coverage relative to the regions
    for region in ["5", "3"]: all_df_windows["coverage_rel_to_%s"%region] = all_df_windows["target_coverage"]/all_df_windows["%s_coverage"%region]

    # get estimate of both relative coverages
    all_df_windows["mean_rel_coverage_to_neighbor"] = (all_df_windows["coverage_rel_to_5"]+all_df_windows["coverage_rel_to_3"])/2
    all_df_windows["closestTo1_rel_coverage_to_neighbor"] = all_df_windows.apply(lambda r: find_nearest([r["coverage_rel_to_5"], r["coverage_rel_to_3"]], 1), axis=1)

    return all_df_windows

def get_clove_output_with_coverage(outfile_clove, reference_genome, sorted_bam, median_coverage, replace=False, run_in_parallel=True, delete_bams=True, threads=4):

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
            coverage_df = get_df_with_coverage_per_windows_relative_to_neighbor_regions(df_TANDELINS, bed_TANDELINS_regions, reference_genome, sorted_bam, df_clove, median_coverage, replace=replace, run_in_parallel=run_in_parallel, delete_bams=delete_bams, threads=threads)
         
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

    Thesse are the SVTYPE fields that are easily classifiable

    DEL "Deletion"
    TAN "Tandem Duplication"
    INV "Inversion"
    INS "Insertion" 
    DUP "Complex Duplication"
    TRA "Complex Translocation"
    ##ALT=<ID=CIV,Description="Complex Inversion">    
    ##ALT=<ID=CIT,Description="Complex Interchromosomal Translocation">
    ##ALT=<ID=CID,Description="Complex Interchromosomal Duplication">
    ##ALT=<ID=IVD,Description="Complex Inverted Interchromosomal Duplication">

    These are the ones set to remaining:

    ##ALT=<ID=CVT,Description="Complex Inverted Translocation"> This is within the same chromosome. CHR2:START-END is cut, inverted and inserted into CHROM:POS

    ##ALT=<ID=CVD,Description="Complex Inverted Duplication"> This is within the same chromosome. CHR2:START-END is duplicated, inverted and inserted into CHROM:POS

    ##ALT=<ID=IVT,Description="Complex Inverted Interchromosomal Translocation"> CHR2:START-END is duplicated, inverted and inserted into CHROM:POS

    define_insertions_based_on_coverage indicates whether to filter the insertions based on coverage. If so, it will set as copied insertions those that have a coverage above 

    """
    #print_if_verbose("getting SVs from clove")

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

        #print_if_verbose("There are %i translocations"%len(df_tra))

    #############################

    ####### INVERSIONS ##########
    if "CIV" in set(df_clove.SVTYPE) and "inversions" in svtypes_to_consider:

        df_inversions = df_clove[df_clove.SVTYPE=="CIV"][["#CHROM", "POS", "END", "ID"]].rename(columns={"#CHROM":"Chr", "POS":"Start", "END":"End"})
        inversions_bed = "%s.inversions.bed"%fileprefix
        df_inversions.to_csv(inversions_bed, sep="\t", header=True, index=False)

        considered_idxs += list(df_inversions.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]
        svtype_to_svfile["inversions"] = inversions_bed
        series_considered_idxs["inversions"] = list(df_inversions.index)


        #print_if_verbose("There are %i inversions"%len(df_inversions))


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

        #print_if_verbose("There are %i insertions, %i of which are copy-and-paste"%(len(df_ins), sum(df_ins.Copied=="TRUE")))

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

            #print_if_verbose("There are %i %s"%(len(df_svtype), typeSV))

    # keep only df_clove that has not been already used, which should be done after each step
    df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

    ###########################


    # write the remaining events which are not easily assignable
    df_notAssigned = df_clove[~(df_clove.SVTYPE.isin(cloveSVtypes_easy_classification))]
    df_notAssigned_file = "%s.remaining.tab"%(fileprefix)
    df_notAssigned[["ID", "#CHROM", "POS", "CHR2", "START", "END", "SVTYPE"]].to_csv(df_notAssigned_file, sep="\t", header=True, index=False)
    svtype_to_svfile["remaining"] = df_notAssigned_file

    #print_if_verbose("There are %i remaining SVs"%len(df_notAssigned))

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

def run_gridssClove_given_filters(sorted_bam, reference_genome, working_dir, median_coverage, replace=True, threads=4, gridss_blacklisted_regions="", gridss_VCFoutput="", gridss_maxcoverage=50000, median_insert_size=250, median_insert_size_sd=0, gridss_filters_dict=default_filtersDict_gridss, tol_bp=50, run_in_parallel=True, max_rel_coverage_to_consider_del=0.2, min_rel_coverage_to_consider_dup=1.8, replace_FromGridssRun=False, define_insertions_based_on_coverage=False):

    """This function runs gridss and clove with provided filtering and parameters. This can be run at the end of a parameter optimisation process. It returns a dict mapping each SV to a table, and a df with the gridss.

    coverage_field is the field where clove is filtered to detect CNV. It can be relative_coverage or relative_coverage_dist_to_telomere.

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
        coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=run_in_parallel, delete_bams=True, threads=threads), sep="\t")

        median_coverage =  get_median_coverage(coverage_df, "")

    print_if_verbose("running gridss and clove with given parameter with %.2f min_rel_coverage_to_consider_dup"%min_rel_coverage_to_consider_dup)
    make_folder(working_dir)

    # edit the replace, regarding if filtering from the run of GRIDSS
    if replace is True and replace_FromGridssRun is False: replace_FromGridssRun = True

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

    # change the name to be always finishing with 'o', and not 'h'
    df_bedpe["name"] = df_bedpe.name.apply(lambda x: x[0:-1]+"o")

    # write
    df_bedpe.to_csv(raw_bedpe_file, sep="\t", header=False, index=False)
    print_if_verbose("there are %i breakpoints"%len(df_bedpe))

    ###################################

    #################################################
    #################################################
    #################################################

    # run clove without checking filtering
    outfile_clove = "%s.clove.vcf"%(raw_bedpe_file)
    run_clove_filtered_bedpe(raw_bedpe_file, outfile_clove, sorted_bam, replace=replace_FromGridssRun, median_coverage=median_coverage, median_coverage_dev=1, check_coverage=False) #  REPLACE debug

    # add the filter of coverage to the clove output
    df_clove = get_clove_output_with_coverage(outfile_clove, reference_genome, sorted_bam, median_coverage, replace=replace_FromGridssRun, run_in_parallel=run_in_parallel, delete_bams=run_in_parallel, threads=threads)

    if len(df_clove)==0: return {}, df_gridss

    # define the coverage filtering based on the type_coverage_to_filterTANDEL
    df_clove["coverage_FILTER"] = df_clove.apply(lambda r: get_covfilter_cloveDF_row_according_to_SVTYPE(r, max_rel_coverage_to_consider_del=max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup=min_rel_coverage_to_consider_dup, coverage_field="mean_rel_coverage_to_neighbor"), axis=1)

    # annotated clove 
    fileprefix = "%s.structural_variants"%outfile_clove

    remaining_df_clove, svtype_to_SVtable = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, fileprefix, reference_genome, sorted_bam, tol_bp=tol_bp, replace=replace_FromGridssRun, svtypes_to_consider={"insertions", "deletions", "inversions", "translocations", "tandemDuplications", "remaining"}, run_in_parallel=run_in_parallel, define_insertions_based_on_coverage=define_insertions_based_on_coverage)

    # merge the coverage files in one
    #merge_coverage_per_window_files_in_one(sorted_bam)

    return svtype_to_SVtable, df_gridss


###################################################################################################
###################################################################################################
###################################################################################################


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

def run_freebayes_parallel_regions(outdir_freebayes, ref, sorted_bam, ploidy, coverage, threads=4, replace=False, pooled_sequencing=False, window_fb=10000):

    """It parallelizes over the provided threads of the system"""

    # make the dir if not already done
    make_folder(outdir_freebayes)

    #run freebayes
    freebayes_output ="%s/output.raw.vcf"%outdir_freebayes; freebayes_output_tmp = "%s.tmp"%freebayes_output
    if file_is_empty(freebayes_output) or replace is True:

        print_if_verbose("running freebayes in parallel with %i threads"%(threads))

        # define the regions file
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

        # init the df that will have all these vcfs and also the header lines
        all_df = pd.DataFrame()
        all_header_lines = []

        ############ FEED THE ALL_DF #############

        size_chunk = 100 # this has to be constant
        for Ic, chunk_regions in enumerate(chunks(regions, size_chunk)):
            print_if_verbose("working on chunk %i/%i"%(Ic, int(len(regions)/size_chunk)))

            # get the df of this chunk
            all_df_chunk_file = "%s/regions_chunk%i_all_df_concatenated.tab"%(regions_vcfs_dir, Ic)
            all_header_lines_chunk_file = "%s/regions_chunk%i_header_lines.py"%(regions_vcfs_dir, Ic)
            if file_is_empty(all_df_chunk_file) or replace is True:

                # restart the folder
                regions_vcfs_dir_chunk = "%s/regions_chunk%i"%(regions_vcfs_dir, Ic)
                make_folder(regions_vcfs_dir_chunk)

                # define the inputs of the function
                inputs_fn = [(region, regions_vcfs_dir_chunk, ref, sorted_bam, ploidy, coverage, replace, pooled_sequencing) for region in chunk_regions]

                # initialize the pool class with the available CPUs --> this is asyncronous parallelization
                with multiproc.Pool(threads) as pool:

                    # run in parallel the freebayes generation for all the 
                    regions_vcfs = pool.starmap(run_freebayes_for_region, inputs_fn)

                    # close the pool
                    pool.close()
                    pool.terminate()

                # go through each of the chromosomal vcfs and append to a whole df
                print_if_verbose("appending all vcfs of the individual regions together")
                all_header_lines_chunk = []
                all_df_chunk = pd.DataFrame()
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
                    all_header_lines_chunk.append("".join([line for line in header_lines if line.split("=")[0] not in {"##reference", "##commandline", "##fileDate"}]))
                    
                    # append to the previous df
                    all_df_chunk = all_df_chunk.append(df[vcf_header], sort=True)

                # delete unnecessary files
                delete_folder(regions_vcfs_dir_chunk)

                # save
                save_object(all_header_lines_chunk, all_header_lines_chunk_file)
                save_df_as_tab(all_df_chunk, all_df_chunk_file)

            # load files
            all_header_lines_chunk = load_object(all_header_lines_chunk_file)
            all_df_chunk = get_tab_as_df_or_empty_df(all_df_chunk_file)

            # add to the final dfs
            all_header_lines += all_header_lines_chunk
            all_df = all_df.append(all_df_chunk)

        ##########################################

        # check that all headers are the same
        if len(set(all_header_lines))!=1: 
            print_if_verbose("These are the header lines: ", set(all_header_lines))
            print_if_verbose("There are %i unique headers"%len(set(all_header_lines)))
            raise ValueError("Not all headers are the same in the individual chromosomal vcfs. This may indicate a problem with parallelization of freebayes")

        # sort and remove duplicate entries
        all_df = all_df.sort_values(by=["#CHROM", "POS"]).drop_duplicates(subset=["#CHROM", "POS", "REF", "ALT"])

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
    return int(nlines/4)

def get_n_pairs_in_fastqgz_gunzip(file, min_lines=10):

    """Takes a fastqgz file and returns the number of reads"""

    # get the number of lines into file
    file_wc = "%s.gunzip.wc"%file
    file_wc_tmp = "%s.tmp"%file_wc
    if file_is_empty(file_wc):

        unpigz_stderr = "%s.generating.stderr"%file_wc_tmp
        print_if_verbose("calculating # reads for %s. The stderr is in %s"%(file, unpigz_stderr))

        run_cmd("gunzip -c %s | wc -l > %s 2>%s"%(file, file_wc_tmp, unpigz_stderr))
        remove_file(unpigz_stderr)

        # test that it worked
        nlines = int(open(file_wc, "r").readlines()[0].strip())
        if nlines<min_lines: raise ValueError("There should be at least %i lines"%min_lines)

        os.rename(file_wc_tmp, file_wc)

    # get the number
    nlines = int(open(file_wc, "r").readlines()[0].strip())

    # check that it is multiple of 4
    if nlines%4!=0: raise ValueError("nlines %i is not valid"%nlines)

    # get the number of lines between 4
    return nlines/4

def get_approx_n_pairs_in_fastqgz(file, nlines=10000):

    """This function calculates the approximate number of lines in a fastq.gz file"""

    print_if_verbose("calculating approximate npairs")

    # get size of a file
    total_gb = os.path.getsize(file)/1e9

    # get a partial file with nlines
    partial_file = "%s.%ilines.fastq.gz"%(file, nlines)

    if file_is_empty(partial_file):

        # define files
        partial_file_tmp = "%s.tmp"%partial_file
        partial_file_stderr = "%s.generating.stderr"%partial_file

        print_if_verbose("getting partial file. The stderr is in %s"%partial_file_stderr)
        run_cmd("zcat %s | head -n %i | gzip > %s 2>%s"%(file, nlines, partial_file_tmp, partial_file_stderr))

        remove_file(partial_file_stderr)
        os.rename(partial_file_tmp, partial_file)


    # calculate the size of the partial
    partial_gb = os.path.getsize(partial_file)/1e9

    # calculate the real nlines
    real_nlines = int((nlines * total_gb) / partial_gb)

    # get the number of lines between 4
    return int(real_nlines/4)

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

    if len(last_reads)>0:

        # change the reads 
        if last_reads[0].endswith("/1") or last_reads[0].endswith("/2"): last_reads = ["/".join(r.split("/")[0:-1]) for r in last_reads]

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
            any_error = any(["ERROR" in l.upper() for l in std_lines])
            if any_error:
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
        any_error = any(["ERROR" in l.upper() for l in std_lines])
        if any_error:
            raise ValueError("Something went wrong with the fastqdump. Check the log in %s"%stdfile)

        remove_file(stdfile)

        # rename the reads
        tmp_reads = "%s/%s"%(tmpdir, get_file(reads))
        os.rename(tmp_reads, reads)

    # delete the tmpdir
    delete_folder(tmpdir)

    return reads

def run_svim(reads, reference_genome, outdir,  threads=4, replace=False, min_sv_size=50, max_sv_size=100000, aligner="ngmlr", is_nanopore=True, minimum_depth=5):

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
        svim_cmd = "%s reads %s %s %s --min_sv_size %i --max_sv_size %i --cores %i --aligner %s --minimum_depth %s --min_mapq 20 > %s 2>&1"%(svim, outdir, reads, reference_genome, min_sv_size, max_sv_size, threads, aligner, minimum_depth, svim_std)
        if is_nanopore is True: svim_cmd += " --nanopore"
        run_cmd(svim_cmd)
        remove_file(svim_std)
        
        os.rename(sorted_bam_long, sorted_bam_short)
        os.rename(sorted_bam_long_idx, sorted_bam_short_idx)

    return sorted_bam_short

def run_sniffles(sorted_bam, outdir, replace, threads, minimum_depth=5, min_sv_size=50):

    """Takes a sorted bam (from svim) and runs sniffles (generating a vcf)."""

    make_folder(outdir)
    output_vcf = "%s/output.vcf"%outdir

    if file_is_empty(output_vcf) or replace is True:

        # define inputs
        output_vcf_tmp = "%s.tmp.vcf"%output_vcf
        stdfile = "%s.generating.std"%output_vcf
        tmp_file = "%s/tmp_file"%outdir

        # run
        print_if_verbose("running SNIFFLES. The std is in %s"%stdfile)
        run_cmd("%s -m %s -v %s --tmp_file %s -s %i -t %i -l %i -q 20 --genotype > %s 2>&1"%(sniffles, sorted_bam, output_vcf_tmp, tmp_file, minimum_depth, threads, min_sv_size, stdfile))

        # clean
        remove_file(stdfile)
        os.rename(output_vcf_tmp, output_vcf)

    return output_vcf


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

        # get the fraction of reads mapped
        fraction_reads_mapped_list.append(fraction_mapped)

    SRA_runInfo_df["fraction_reads_mapped"] = fraction_reads_mapped_list

    ############################################

    return SRA_runInfo_df

def get_taxID_or_BioSample(r, target_taxID):

    """Takes a row of the all_SRA_runInfo_df and returns the TaxID. If the taxID is the target_taxID it returns the biosample"""

    if r["TaxID"]==target_taxID: return r["BioSample"]
    else: return r["TaxID"]


def get_ancestor_taxID(target_taxID, nancestorNodes, outdir):

    """This function gets the ancestor taxID on nancestorNodes.

    ancestor_taxID = int(ncbi.get_lineage(target_taxID)[-nancestorNodes]) """


    cmds_ancestor_taxID = ["from ete3 import NCBITaxa; ncbi = NCBITaxa()",
                           "ancestor_taxID = int(ncbi.get_lineage(%i)[-%i])"%(target_taxID, nancestorNodes),
                           "print(ancestor_taxID)"]

    ancestor_taxID_std = "%s/ancestor_taxID.std"%outdir              
    run_cmd("python -c '%s' > %s"%("; ".join(cmds_ancestor_taxID), ancestor_taxID_std), env=EnvName_ete3)

    ancestor_taxID = int(open(ancestor_taxID_std, "r").readlines()[0])

    return ancestor_taxID

def get_SRA_runInfo_df(target_taxID, n_close_samples, nruns_per_sample, outdir, reference_genome, min_coverage, replace, threads, coverage_subset_reads, min_fraction_reads_mapped, get_lowest_coverage_possible=False):

    """This function mines the SRA to find n_close_samples and nruns_per_sample, returning the necessary df """

    # check if you have network access
    if connected_to_network() is False: raise ValueError("There is no network connection available, which is necessary to get the get_SRA_runInfo_df working")

    ######## UPDATE NCBI TAXONOMY ########

    # this will update the current database in the computer
    print_if_verbose("Getting genomes for taxID into %s"%(outdir))

    # load the NCBI taxonomy database and upgrade it if not already done
    print_if_verbose("getting NCBI taxonomy database. This may fail if you already had installed an NCBI Taxonomy database from ete3 before.")

    # change the dir and the taxdump
    curdir = get_fullpath(os.getcwd())
    outdir = get_fullpath(outdir)
    
    ncbiTaxa_updated_file = "%s/ncbiTaxa_updated.txt"%outdir
    if file_is_empty(ncbiTaxa_updated_file) or replace is True:

        dir_updating_ete3 = "%s/update_taxonomy_database"%outdir; make_folder(dir_updating_ete3)
        os.chdir(dir_updating_ete3)
        print_if_verbose("updating db into %s"%dir_updating_ete3)

        # update the ncbi taxonomy database
        #os.system("rm -r ~/.etetoolkit/") # remove the previous ncbi tax database. This is not always necessary
        cmd_update = "from ete3 import NCBITaxa; ncbi = NCBITaxa(); ncbi.update_taxonomy_database()"
        run_cmd("python -c '%s'"%cmd_update, env=EnvName_ete3)

        # delete 
        delete_folder(dir_updating_ete3)

        # write file
        open(ncbiTaxa_updated_file, "w").write("NCBItaxa updated\n")

    #######################################

    # set dir to the outdir
    os.chdir(curdir)

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

        # get the ancestor taxID
        ancestor_taxID = get_ancestor_taxID(target_taxID, nancestorNodes, outdir_ancestors)

        # get the runs for this division
        print_if_verbose("Getting WGS info")
        fileprefix = "%s/output"%(outdir_ancestors)
        all_SRA_runInfo_df = get_allWGS_runInfo_fromSRA_forDivision(fileprefix, ancestor_taxID, reference_genome, taxIDs_to_exclude=set(), replace=False, min_coverage=min_coverage).set_index("Run", drop=False)

        # exclude the taxID
        if any(pd.isna(all_SRA_runInfo_df.TaxID)) or any(all_SRA_runInfo_df.TaxID.apply(lambda x: type(x)!=int)): raise ValueError("TaxID is not proerly formated t in all_SRA_runInfo_df")
        all_SRA_runInfo_df = all_SRA_runInfo_df[all_SRA_runInfo_df.TaxID!=target_taxID]

        # if it is empty, continue
        if len(all_SRA_runInfo_df)==0: continue
        print_if_verbose(all_SRA_runInfo_df[["TaxID", "ScientificName"]])

        # exclude the wrong SRRs
        all_SRA_runInfo_df = all_SRA_runInfo_df[~all_SRA_runInfo_df.Run.isin(wrong_SRRs)]

        # if it is empty, continue
        if len(all_SRA_runInfo_df)==0: continue

        # define the runs with target taxID
        runs_target_taxID = set(all_SRA_runInfo_df[all_SRA_runInfo_df.TaxID==target_taxID].Run)

        # get the interesting taxIDs, the distance to the target and the scientific name with an external script
        outfile_interesting_objects = "%s/interestingTaxIDs_distanceToTarget_taxID_to_sciName.py"%outdir_ancestors

        run_cmd("%s %i %i %s"%(get_interestingTaxIDs_distanceToTarget_taxID_to_sciName_py, ancestor_taxID, target_taxID, outfile_interesting_objects), env=EnvName_ete3)

        interesting_taxIDs, taxID_to_distanceToTarget, taxID_to_sciName = load_object(outfile_interesting_objects)

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


def get_median_readLength_fastqgz(fastqgz, nlines=1000, replace=False):

    """Takes a fastqgz file and returns the median read length"""

    # define the file
    read_len_file = "%s.readLengths"%fastqgz
    read_len_file_tmp = "%s.tmp"%read_len_file

    if file_is_empty(read_len_file) or replace is True:

        read_len_file_stderr = "%s.getting.stderr"%read_len_file
        print_if_verbose("getting median read lengths. The stderr is in %s"%read_len_file_stderr)
        run_cmd("zcat %s | head -n %i | sed -n '2~4p' > %s 2>%s"%(fastqgz, nlines, read_len_file_tmp, read_len_file_stderr))

        remove_file(read_len_file_stderr)
        os.rename(read_len_file_tmp, read_len_file)

    return int(np.median([len(l.strip()) for l in open(read_len_file, "r").readlines()]))

def generate_downsampledReads(fastqgz, downsampled_fastqgz, fraction_downsample, replace=False):

    """Takes some reads and downsamples them to fraction_downsample into downsampled_fastqgz."""

    if file_is_empty(downsampled_fastqgz) or replace is True:

        downsampled_fastqgz_tmp = "%s.tmp.fastq.gz"%downsampled_fastqgz
        downsampled_fastqgz_stderr = "%s.generating.stderr"%downsampled_fastqgz
        print_if_verbose("downsampling into %s. The stderr is in %s"%(downsampled_fastqgz, downsampled_fastqgz_stderr))

        try: run_cmd("%s sample -s100 %s %.4f | %s > %s 2>%s"%(seqtk, fastqgz, fraction_downsample, pigz, downsampled_fastqgz_tmp, downsampled_fastqgz_stderr))
        except:
            print_if_verbose("downsampling with pigz failed. downsampling with gzip")
            run_cmd("%s sample -s100 %s %.4f | gzip > %s 2>%s"%(seqtk, fastqgz, fraction_downsample, downsampled_fastqgz_tmp, downsampled_fastqgz_stderr))

        remove_file(downsampled_fastqgz_stderr)
        os.rename(downsampled_fastqgz_tmp, downsampled_fastqgz)

def downsample_close_shortReads_table(close_shortReads_table, close_shortReads_table_df, max_coverage_sra_reads, reference_genome, replace=False, threads=4):

    """This function takes a df with the close_reads downloaded from SRA and replaces them with a downsampled version"""

    # initialize the corresponding file
    new_close_shortReads_table = "%s.max_%ix"%(close_shortReads_table, max_coverage_sra_reads)

    if file_is_empty(new_close_shortReads_table) or replace is True or not close_shortReads_table_is_correct(new_close_shortReads_table):
        print_if_verbose("getting the new fastqc files to coverage %ix..."%max_coverage_sra_reads)

        # initialize a df that has the modified reads (possibly downsampled)
        new_close_shortReads_table_df = pd.DataFrame(columns=close_shortReads_table_df.columns)

        # go through each runID
        for I, row in close_shortReads_table_df.iterrows():
            print_if_verbose(row["runID"])

            # calculate the genome size
            genome_length = sum(get_chr_to_len(reference_genome).values())

            # calculate the coverage
            fraction_downsample_file = "%s.fraction_downsample.py"%row["short_reads1"]

            if file_is_empty(fraction_downsample_file) or replace is True:
                print_if_verbose("calculating the fraction to downsample")

                # calculate the read length
                read_len = get_median_readLength_fastqgz(row["short_reads1"], replace=replace)

                # calculate the number of reads
                npairs = get_approx_n_pairs_in_fastqgz(row["short_reads1"]) # approximate, faster way
                #npairs = get_n_pairs_in_fastqgz(row["short_reads1"])

                # calculate the expected coverage
                expected_coverage = (npairs*read_len)/genome_length
                print_if_verbose("The expected coverage is %.3fx"%expected_coverage)

                # define the maximum number of read pairs and the fraction to downsample
                max_npairs = (max_coverage_sra_reads*genome_length)/read_len
                fraction_downsample = max_npairs/npairs
                
                # save
                save_object(fraction_downsample, fraction_downsample_file)

            else: fraction_downsample = load_object(fraction_downsample_file)

            # downsample if the expected coverage is above max_coverage_sra_reads
            if fraction_downsample < 1:

                # downsample
                new_short_reads1 = "%s.%ix.fastq.gz"%(row["short_reads1"], max_coverage_sra_reads)
                new_short_reads2 = "%s.%ix.fastq.gz"%(row["short_reads2"], max_coverage_sra_reads)

                if file_is_empty(new_short_reads1) or file_is_empty(new_short_reads2):

                    new_short_reads1_tmp = "%s.tmp.fastq.gz"%new_short_reads1
                    new_short_reads2_tmp = "%s.tmp.fastq.gz"%new_short_reads2

                    generate_downsampledReads(row["short_reads1"], new_short_reads1_tmp, fraction_downsample, replace=replace)
                    generate_downsampledReads(row["short_reads2"], new_short_reads2_tmp, fraction_downsample, replace=replace)

                    # check that the reads are correct
                    check_that_paired_reads_are_correct(new_short_reads1_tmp, new_short_reads2_tmp)

                    os.rename(new_short_reads1_tmp, new_short_reads1)
                    os.rename(new_short_reads2_tmp, new_short_reads2)

            else:
                new_short_reads1 = row["short_reads1"]
                new_short_reads2 = row["short_reads2"]

            # add to df
            df_dict = {"sampleID":row["sampleID"], "runID":row["runID"], "short_reads1":new_short_reads1, "short_reads2":new_short_reads2}
            df = pd.DataFrame({I : df_dict}).transpose()
            new_close_shortReads_table_df = new_close_shortReads_table_df.append(df, sort=True)

        # save
        new_close_shortReads_table_tmp = "%s.tmp"%new_close_shortReads_table
        new_close_shortReads_table_df.to_csv(new_close_shortReads_table_tmp, sep="\t", index=False, header=True)
        os.rename(new_close_shortReads_table_tmp, new_close_shortReads_table)


    else: new_close_shortReads_table_df = pd.read_csv(new_close_shortReads_table, sep="\t")

    return new_close_shortReads_table, new_close_shortReads_table_df



def get_close_shortReads_table_close_to_taxID(target_taxID, reference_genome, outdir, ploidy, n_close_samples=3, nruns_per_sample=3, replace=False, threads=4, min_fraction_reads_mapped=0.4, coverage_subset_reads=0.1, min_coverage=30, job_array_mode="local", StopAfter_sampleIndexingFromSRA=False, StopAfterPrefecth_of_reads=False, get_lowest_coverage_possible=False, max_coverage_sra_reads=10000000000000000):

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

    # write a tab
    save_df_as_tab(SRA_runInfo_df, "%s/final_SRA_runInfo_df.tab"%outdir)

    print_if_verbose("these are the samples chosen:\n:", SRA_runInfo_df[["Run", "sampleID", "SampleName", "sci_name", "fraction_reads_mapped"]].sort_values("sampleID"))

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
            print_if_verbose("getting %s reads"%srr)

            # make the outdir for this
            outdir_srr = "%s/%s"%(reads_dir, srr); make_folder(outdir_srr)

            # define files
            trimmed_reads1 = "%s/%s_trimmed_reads_1.fastq.gz"%(outdir_srr, srr)
            trimmed_reads2 = "%s/%s_trimmed_reads_2.fastq.gz"%(outdir_srr, srr)
            SRRfile = "%s/%s.srr"%(outdir_srr, srr)

            # keep the files
            srr_to_readsDict[srr] = {"short_reads1":trimmed_reads1, "short_reads2":trimmed_reads2}

            # define the finalisation file
            if StopAfterPrefecth_of_reads is True: final_file = SRRfile
            else: final_file = trimmed_reads2 # the last one generated
  
            # get the cmd if necessary
            if (file_is_empty(final_file) and file_is_empty(trimmed_reads2)) or replace is True:

                # define the cmd and add it
                cmd = "%s --srr %s --outdir %s --threads %i"%(get_trimmed_reads_for_srr_py, srr, outdir_srr, threads)
                if StopAfterPrefecth_of_reads is True: cmd += " --stop_after_prefetch"

                all_cmds.append(cmd)
                continue

        # if there are all_cmds, run them in a job array
        if len(all_cmds)>0:
            print_if_verbose("Getting reads into %s"%(reads_dir))

            # if you are in local, run them in parallel. Each job in one threads
            if job_array_mode=="local" and StopAfterPrefecth_of_reads is True: 

                print_if_verbose("running prefetch in parallel")
                with multiproc.Pool(threads) as pool:
                    pool.starmap(run_cmd, [(x,) for x in all_cmds])
                    pool.close()

            # local
            elif job_array_mode=="local" and StopAfterPrefecth_of_reads is False: 

                print_if_verbose("running all files one after the other")
                for cmd in all_cmds: run_cmd(cmd)

            # if you are in a greasy environment
            elif job_array_mode=="job_array": 

                read_downloading_dir = "%s/read_downloading_files"%outdir; make_folder(read_downloading_dir)
                print_if_verbose("Submitting %i downloading reads jobs to a job_array. All files will be stored in %s"%(len(all_cmds), read_downloading_dir))

                jobs_filename = "%s/jobs.getting_SRAdatasets"%read_downloading_dir
                open(jobs_filename, "w").write("\n".join(all_cmds))

                generate_jobarray_file(jobs_filename, "gettingCloseShortReads")

                # exit before it starts
                print_if_verbose("You need to run all the jobs from %s"%jobs_filename)
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

    # define the initially important files
    initial_important_files = set(close_shortReads_table_df["short_reads1"]).union(close_shortReads_table_df["short_reads2"])

    # downsample the reads according to max_coverage_sra_reads
    print_if_verbose("downsampling the reads to a max coverage of %ix"%max_coverage_sra_reads)
    close_shortReads_table, close_shortReads_table_df = downsample_close_shortReads_table(close_shortReads_table, close_shortReads_table_df, max_coverage_sra_reads, reference_genome, replace=replace, threads=threads)

    # debug
    if not close_shortReads_table_is_correct(close_shortReads_table): raise ValueError("%s has empty reads files"%close_shortReads_table)
    
    # define the important files
    important_files = set(close_shortReads_table_df["short_reads1"]).union(close_shortReads_table_df["short_reads2"]).union(initial_important_files)

    # remove the files that are not the reads
    for srr in os.listdir(reads_dir):
        srr_dir = "%s/%s"%(reads_dir, srr)
        for file in os.listdir(srr_dir):
            filepath = "%s/%s"%(srr_dir, file)

            if filepath not in important_files: remove_file(filepath)

    # remove all the other unimportant files and folders
    for f in os.listdir(outdir):

        if f not in {"final_SRA_runInfo_df.py", "reads", get_file(close_shortReads_table)}:

            path = "%s/%s"%(outdir, f)
            delete_folder(path)
            remove_file(path)

    # last check
    if not close_shortReads_table_is_correct(close_shortReads_table): raise ValueError("%s has empty reads files"%close_shortReads_table)

    #########################################

    return close_shortReads_table

def get_is_matching_predicted_and_known_rows(rk, rp, equal_fields, approximate_fields, chromField_to_posFields, tol_bp=50, pct_overlap=0.75):

    """Takes a row of a knownID (rk) and a predictedID (rp) and returns a boolean indicating if they match. These rk and rp can be any dict-like structures that have the expected equal_fields and so."""

    # ask if the equal fields match, and return if false
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
        #print_if_verbose("benchmarking %s"%svtype)

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


def get_IDstring_for_svDF_r(r, svtype):

    """This function takes a row of an svDF and returns the ID as a string"""

    if svtype=="insertions":

        bool_to_text = {True:"copyPaste", False:"cutPaste"}
        ID =  "INS|%s:%i-%i|%s:%i|%s"%(r["ChrA"], r["StartA"], r["EndA"], r["ChrB"], r["StartB"], bool_to_text[r["Copied"]])

    elif svtype=="remaining":

        if r["SVTYPE"] in {"ITX1", "ITX2", "INVTX1", "INVTX2", "TAN", "DEL", "INV1", "INV2"}: 

            ID = "%slike|%s:%i-%s:%i"%(r["SVTYPE"], r["#CHROM"], r["POS"], r["CHR2"], r["END"])

        elif r["SVTYPE"] in {"CVT", "CVD", "IVT"}:

            # these are all SVs that are understood but not among the generally classified events. They imply one region copied or translocated into another
            ID = "%s|%s:%i-%i|%s:%i"%(r["SVTYPE"], r["CHR2"], r["START"], r["END"], r["#CHROM"], r["POS"])

        else: raise ValueError("%s is not considered"%r["SVTYPE"])

    elif svtype in {"tandemDuplications", "deletions", "inversions"}:

        sv_to_tag = {"tandemDuplications":"TDUP", "deletions":"DEL", "inversions":"INV"}
        ID = "%s|%s:%i-%i"%(sv_to_tag[svtype], r["Chr"], r["Start"], r["End"])

    elif svtype=="translocations":

        ID = "TRA|%s:%i-%i<>%s:%i-%i"%(r["ChrA"], r["StartA"], r["EndA"], r["ChrB"], r["StartB"], r["EndB"])

    else: 
        print(r)
        raise ValueError("%s has not been considered"%r["SVTYPE"])

    return ID



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
        svDF["real_AF_%s"%estimate_fn_name] = svDF.bends_metadata_dict.apply(get_estimate_AF_for_breakends, estimate_fn=estimate_fn)

    # get the implicated breakpoint IDs
    svDF["bpIDs"] =  svDF.bends_metadata_dict.apply(lambda x: ",".join(sorted(x)))

    # get the quality
    for estimate_fn_name, estimate_fn in [("min", min), ("max", max), ("mean", np.mean)]:
        svDF["QUAL_%s"%estimate_fn_name] = svDF.bends_metadata_dict.apply(get_estimate_AF_for_breakends, AF_field="QUAL", estimate_fn=estimate_fn)


    # get the worst filter tag
    svDF["all_FILTERs"] =  svDF.bends_metadata_dict.apply(lambda x: ",".join(set.union(*[set.union(*[set(bend_info["FILTER"].split(";")) for bend_info in list_breakend_info]) for list_breakend_info in x.values()])))

    # get the breakpoint IDs
    svDF["BREAKPOINTIDs"] = svDF.bends_metadata_dict.apply(lambda x: ",".join(sorted(x.keys())))

    # add ID 
    svDF["IDstring"] = svDF.apply(lambda r: get_IDstring_for_svDF_r(r, svtype), axis=1)

    return svDF

def get_sampleID_to_svtype_to_svDF_filtered(sampleID_to_svtype_to_file, sampleID_to_dfGRIDSS, sampleID_to_parentIDs={}, breakend_info_to_keep=['#CHROM', 'POS', 'other_coordinates', 'allele_frequency', 'allele_frequency_SmallEvent', 'real_AF', 'FILTER', 'inserted_sequence', 'has_poly16GC', 'length_inexactHomology', 'length_microHomology', 'QUAL', 'overlaps_repeats', 'REF']):

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

        # add the breakpoint ID
        df_gridss["BREAKPOINTID"] = df_gridss.eventID

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

                svDF["highAF"] = (svDF.real_AF_min>=min_af)

                # define those that have high confidence as those that have the correct af and are in all replicateIDs
                svDF["high_confidence"] = ((svDF.sv_in_any_replicateIDs) | ( (svDF.highAF) & ((svDF.all_bends_highQUAL) | (svDF.all_bends_PASS)) )) & ~(svDF.sv_in_75pct_samples)

                # keep the df that has high confidence
                ID_to_svtype_to_svDF[ID][svtype] = svDF[svDF.high_confidence]

def get_is_overlapping_query_vs_target_region(q, r):

    """This function takes two 'bed'-like regions and returns whether they are overlapping by some extent """

    return (q["chromosome"]==r["chromosome"]) and ((r["start"]<=q["start"]<=r["end"]) or (r["start"]<=q["end"]<=r["end"]) or (q["start"]<=r["start"]<=q["end"]) or (q["start"]<=r["end"]<=q["end"]))


def get_svtype_to_svfile_from_perSVade_outdir(perSVade_outdir):

    """This function takes from the perSVade outdir the svdict"""

    # define the paths
    outdir = "%s/SVdetection_output/final_gridss_running"%perSVade_outdir
    gridss_vcf = "%s/gridss_output.raw.vcf"%outdir # raw

    # this means that it is a cleaned dir
    if not file_is_empty(gridss_vcf):

         # get the svtype_to_svfile
        svtype_to_svfile = {svtype : "%s/%s.tab"%(outdir, svtype)  for svtype in {"insertions", "deletions", "tandemDuplications", "translocations", "inversions"}}
        svtype_to_svfile["remaining"] = "%s/unclassified_SVs.tab"%outdir 

    else: 

        # assume that it is a not cleaned dir
        gridss_vcf = "%s/gridss_output.vcf.withSimpleEventType.vcf"%outdir # raw
        #gridss_vcf = "%s/gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf"%outdir # filtered
        svtype_to_svfile = {file.split(".structural_variants.")[1].split(".")[0] : "%s/%s"%(outdir, file) for file in os.listdir(outdir) if ".structural_variants." in file}

    # keep only the ones that exist
    svtype_to_svfile = {svtype : file for svtype, file in svtype_to_svfile.items() if not file_is_empty(file)}

    return svtype_to_svfile


def get_svtype_to_svfile_and_df_gridss_from_perSVade_outdir(perSVade_outdir, reference_genome):

    """This function takes from the perSVade outdir the svdict and the df_gridss"""

    # define the paths
    outdir = "%s/SVdetection_output/final_gridss_running"%perSVade_outdir
    gridss_vcf = "%s/gridss_output.raw.vcf"%outdir # raw
    #gridss_vcf = "%s/gridss_output.filt.vcf"%outdir # filt

    # this means that it is a cleaned dir
    if not file_is_empty(gridss_vcf):

         # get the svtype_to_svfile
        svtype_to_svfile = {svtype : "%s/%s.tab"%(outdir, svtype)  for svtype in {"insertions", "deletions", "tandemDuplications", "translocations", "inversions"}}
        svtype_to_svfile["remaining"] = "%s/unclassified_SVs.tab"%outdir 

    else: 

        # assume that it is a not cleaned dir
        gridss_vcf = "%s/gridss_output.vcf.withSimpleEventType.vcf"%outdir # raw
        #gridss_vcf = "%s/gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf"%outdir # filtered
        svtype_to_svfile = {file.split(".structural_variants.")[1].split(".")[0] : "%s/%s"%(outdir, file) for file in os.listdir(outdir) if ".structural_variants." in file}

    # get the df gridss
    print_if_verbose("loading df_gridss")
    df_gridss = add_info_to_gridssDF(load_single_sample_VCF(gridss_vcf), reference_genome)

    # keep only the ones that exist
    svtype_to_svfile = {svtype : file for svtype, file in svtype_to_svfile.items() if not file_is_empty(file)}
    print_if_verbose("There are %i svfiles"%len(svtype_to_svfile))

    return svtype_to_svfile, df_gridss


def get_is_protein_altering_consequence(consequence):

    """Takes tha consequences and returns whether they are protein altering"""

    # get as set
    consequences_set = set(consequence.split(","))

    # ask
    if len(consequences_set.intersection(PROT_ALTERRING_MUTATIONS))>0: return True
    elif len(consequences_set.difference(NON_PROT_ALTERRING_MUTATIONS))==0: return False
    else: raise ValueError("%s contains non-described vars"%consequences_set)

def get_is_transcript_disrupting_consequence_SV(consequence):

    """Takes a consequence of a SVs annotations df and returns whether they are transcript disrupting."""

    # get the consequences
    consequences_set = set(consequence.split(","))

    # define whether the consequences are prot_alterring
    if len(consequences_set.intersection(SVs_TRANSCRIPT_DISRUPTING_MUTATIONS))>0: return True
    elif len(consequences_set.difference(SVs_NON_TRANSCRIPT_DISRUPTING_MUTATIONS))==0: return False
    else: raise ValueError("%s contains non-described vars"%consequences_set)


def copy_file(origin, target):

    """Copy a file with tmp"""

    target_tmp = "%s.tmp"%target

    if file_is_empty(target):

        shutil.copy2(origin, target_tmp)
        os.rename(target_tmp, target)


def rsync_file(origin, target):

    """Copy a file with tmp and rsync"""

    target_tmp = "%s.tmp"%target

    if file_is_empty(target):

        run_cmd("rsync %s %s"%(origin, target_tmp))
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

       # files under CNV calling
       "SVcalling_output/calculating_CNVcoverage"
    ]

    # add all the temporary files
    files_to_remove += [f for f in os.listdir(outdir) if "temporary_file" in f or f.endswith(".tmp") or "coverage_per_window.tab." in f] 

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
                                               "benchmarking_all_filters_for_all_genomes_and_ploidies/plots",
                                               "df_CNV_allKnownRegions.tab"
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

    ######## CNVcalling files ##########

    ####################################

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

def get_compatible_real_bedpe_breakpoints(close_shortReads_table, reference_genome, outdir, replace=False, threads=4, mitochondrial_chromosome="mito_C_glabrata_CBS138", job_array_mode="local", max_nvars=100, parameters_json_file=None):

    """Generates a file under outdir that has the stacked breakpoints arround which to generate SV calls
    realSV_calling_on can be reads or assembly"""

    # load the df
    df_genomes = pd.read_csv(close_shortReads_table, sep="\t").set_index("runID")

    # define an outdir that will store all the real_vars
    make_folder(outdir)
    all_realVars_dir = "%s/all_realVars"%(outdir)
    if replace is True: delete_folder(all_realVars_dir)
    make_folder(all_realVars_dir)

    # initialize a list of cmds to run
    all_cmds = []

    # define the name of the final file name
    final_file_name = "perSVade_finished_file.txt"

    # check if there are some jobs to run SV calling on
    njobs_to_run_SVcalling_on = sum([file_is_empty("%s/shortReads_realVarsDiscovery_%s/%s"%(all_realVars_dir,ID, final_file_name)) for ID, row in df_genomes.iterrows()])

    print_if_verbose("There are %i jobs still to run"%njobs_to_run_SVcalling_on)

    # init a dict with the timing info
    timiming_dict = {}

    # generate all real vars
    for ID, row in df_genomes.iterrows():
        print_if_verbose(ID)

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
            cmd = "python %s -r %s --threads %i --outdir %s  --mitochondrial_chromosome %s --fast_SVcalling --previous_repeats_table %s --min_CNVsize_coverageBased %i --skip_CNV_calling --skip_SV_CNV_calling"%(perSVade_py, reference_genome, threads, outdir_gridssClove, mitochondrial_chromosome, previous_repeats_table, min_CNVsize_coverageBased)

            # add arguments depending on the pipeline
            if replace is True: cmd += " --replace"
            if parameters_json_file is not None: cmd += " --parameters_json_file %s"%parameters_json_file

            # add the input
            all_keys_df = set(df_genomes.keys())

            # reads
            if "short_reads1" in all_keys_df and "short_reads2" in all_keys_df: cmd += " -f1 %s -f2 %s"%(row["short_reads1"], row["short_reads2"])

            # bams
            elif "sorted_bam" in all_keys_df: cmd += " -sbam %s"%(row["sorted_bam"])

            else: raise ValueError("The provided close_shortReads_table is not valid") 

            # if the running in slurm is false, just run the cmd
            if job_array_mode=="local": run_cmd(cmd)
            elif job_array_mode=="job_array": 
                all_cmds.append(cmd)
                continue

            else: raise ValueError("%s is not valid"%job_array_mode)

        else:

            pass
            # get the timings 
            timiming_dict[ID] = {l.split(":")[0].split("time_")[1] : float(l.strip().split(":")[1])/3600 for l in open(final_file, "r").readlines() if l.startswith("time_")}


    # if yoy are running on slurm, get it in a job array
    if job_array_mode=="job_array": 

        if len(all_cmds)>0: 
            print_if_verbose("submitting %i jobs to the cluster for the real data. The files can be monitored from %s"%(len(all_cmds), all_realVars_dir))
            jobs_filename = "%s/jobs.getting_realSVs"%all_realVars_dir
            open(jobs_filename, "w").write("\n".join(all_cmds))

            generate_jobarray_file(jobs_filename, "compatibleRealBedpeObtention")

            print_if_verbose("Exiting... You have to wait until all the jobs in testRealSVs are done. Wait until the jobs are done and rerun this pipeline to continue")
            sys.exit(0)

    # print the time that it took each sample
    timiming_df = pd.DataFrame(timiming_dict).transpose()
    timiming_df["sample_runID"] = timiming_df.index
    save_df_as_tab(timiming_df, "%s/timing_data.tab"%all_realVars_dir)
    #print_if_verbose(timiming_df)

    # get the 
    bedpe_fields = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "ID", "score", "or1", "or2"]
    df_bedpe = pd.DataFrame(columns=bedpe_fields)

    # remove the bam files and reference files, while keeping the bedpe into a df
    print_if_verbose("cleaning files")
    for ID, row in df_genomes.iterrows():

        # run in the gridss and clove with the fast parameters
        outdir_gridssClove = "%s/shortReads_realVarsDiscovery_%s"%(all_realVars_dir,ID)

        # keep the bedpe
        bedpe_file = "%s/SVdetection_output/final_gridss_running/gridss_output.filt.bedpe"%(outdir_gridssClove)
        df_bedpe = df_bedpe.append(pd.read_csv(bedpe_file, sep="\t", names=bedpe_fields, header=-1))

        # clean again
        clean_perSVade_outdir(outdir_gridssClove)

        # remove the file that has the sample
        #remove_file("%s/df_gridss_svtype_to_svfile_tuple_%s.py"%(all_realVars_dir,ID))

    real_bedpe_breakpoints = "%s/integrated_breakpoints.bedpe"%outdir
    df_bedpe = df_bedpe.drop_duplicates(subset=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "or1", "or2"])
    df_bedpe.to_csv(real_bedpe_breakpoints, sep="\t", header=False, index=False)

    return real_bedpe_breakpoints




def add1_unless_it_is_minus1(x):

    """Takes an int and adds 1 unless it is -1"""

    if x==-1: return x
    else: return x+1 

def set_position_to_max(pos, maxPos):

    """Sets a position to a maximum"""

    if pos>maxPos: return maxPos
    else: return pos


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
        run_cmd("%s CollectInsertSizeMetrics HISTOGRAM_FILE=%s INPUT=%s OUTPUT=%s > %s 2>&1"%(picard_exec, hist_file, sampled_bam, outfile_tmp, picard_insertSize_std), env=EnvName_picard)
        remove_file(picard_insertSize_std)
        remove_file(sampled_bam)

        os.rename(outfile_tmp, outfile)

    # get stats
    wrong_foot_lines = [l for l in open(outfile, "r").readlines() if len(l.split("\t"))==2 and not l.startswith("## METRICS CLASS")]
    df = pd.read_csv(outfile, sep="\t", skip_blank_lines=False, header=6, skipfooter=len(wrong_foot_lines)+2, engine='python').iloc[0]

    return (int(float(df["MEDIAN_INSERT_SIZE"])), int(float(df["MEDIAN_ABSOLUTE_DEVIATION"])))



def index_genome(genome, replace=False):

    """Takes a fasta and generates a <genome>.fai file"""

    # index the genome of interest if not already done
    if file_is_empty("%s.fai"%genome) or replace is True: 

        faidx_std = "%s.indexing.std"%genome
        print_if_verbose("running faidx. The std is in %s"%faidx_std)
        run_cmd("%s faidx %s > %s 2>&1"%(samtools, genome, faidx_std))
        remove_file(faidx_std)

def create_sequence_dict(genome, replace=False):

    """Takes a fasta and generates the reference dict"""

    rstrip = genome.split(".")[-1]
    dictionary = "%sdict"%(genome.rstrip(rstrip)); tmp_dictionary = "%s.tmp"%dictionary

    if file_is_empty(dictionary) or replace is True:

        # remove any previously created tmp_file
        remove_file(tmp_dictionary)

        # define the std
        dictionary_std = "%s.generating.std"%dictionary
        print_if_verbose("Creating picard dictionary. The std is in %s"%dictionary_std)

        run_cmd("%s CreateSequenceDictionary R=%s O=%s TRUNCATE_NAMES_AT_WHITESPACE=true > %s 2>&1"%(picard_exec, genome, tmp_dictionary, dictionary_std), env=EnvName_picard) 

        remove_file(dictionary_std)  
        os.rename(tmp_dictionary , dictionary)


def get_windows_infoDF_with_predictedFromFeatures_coverage(genome, distToTel_chrom_GC_to_coverage_fn, expected_coverage_per_bp, replace=False, threads=4, make_plots=True):

    """This function gets a genome and returns a df for windows of the genome and the relative coverage predicted from distToTel_chrom_GC_to_coverage_fn"""

    windows_infoDF_file = "%s_windows_with_predictedRelCov_from_features.tab"%genome

    if file_is_empty(windows_infoDF_file) or replace is True:
        print_if_verbose("getting relCov predicted from feats")

        # index the genome of interest if not already done
        index_genome(genome, replace=replace)

        ##### get the windows df ####

        # get the file
        windows_file = "%s.windows%ibp.bed"%(genome, window_l)
        windows_file_stderr = "%s.generating.stderr"%windows_file
        print_if_verbose("running makewindows. The stderr is in %s"%windows_file_stderr)
        run_cmd("%s makewindows -g %s.fai -w %i > %s 2>%s"%(bedtools, genome, window_l, windows_file, windows_file_stderr)) # debug
        remove_file(windows_file_stderr)

        # get into df 
        df = pd.read_csv(windows_file, sep="\t", header=None, names=["chromosome", "start", "end"])

        #############################

        # get the distance to the telomere in the df
        chr_to_len = get_chr_to_len(genome)
        df["middle_position"] = (df.start + (df.end - df.start)/2).apply(int)
        df["distance_to_telomere"] = df.apply(lambda r: min([r["middle_position"], chr_to_len[r["chromosome"]]-r["middle_position"]]), axis=1)

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
        df.to_csv(windows_infoDF_file, sep="\t", header=True, index=False)

    else: df = pd.read_csv(windows_infoDF_file, sep="\t")

    print_if_verbose("you can test the predicted coverage in %s"%windows_infoDF_file)

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

def run_wgsim_pairedEnd_per_windows_in_parallel(df_windows, genome, outdir, read_length,  median_insert_size, median_insert_size_sd, replace=False, error_rate=0.02, threads=4):

    """Takes a dataframe with windows of ["chromosome", "start", "end", "readPairs"] and writes, under outdir, two fastq.gz files of simulated paired end  reads. The parallel runs are written in subfolder under outir which is finally removed. It returns a tuple of the two fastq.gz files.

    max_n_windows_at_once indicates the number of windows that can be processed at once. """

    # define the max_n_windows_at_once as 2x the threads
    max_n_windows_at_once = threads

    # make the outdir if it does not exist
    make_folder(outdir)
    print_if_verbose("running simulation in a parallelized manner. Each window will be run in a different core. Running on chunks of %i windows"%max_n_windows_at_once)

    # define the final outdirs
    allChunks_fastqgz_1 = "%s/allChunks_reads1.fq.gz"%outdir; allChunks_fastqgz_2 = "%s/allChunks_reads2.fq.gz"%outdir

    if any([file_is_empty(f) for f in [allChunks_fastqgz_1, allChunks_fastqgz_2]]) or replace is True:

        # define the important files
        df_windows = df_windows[["chromosome", "start", "end", "readPairs"]]
        df_windows.index = list(range(len(df_windows)))

        # keep only windows with at more than 1 (pseudocount) red pair
        df_windows = df_windows[df_windows.readPairs>1]

        # define chunks of indices
        chunks_indices = list(chunks(list(df_windows.index), max_n_windows_at_once))

        # initialize chunks files
        chunks_fastq_files = []

        # go through slices of max_n_windows_at_once of the df_windows
        print_if_verbose("there are %i windows"%len(df_windows))
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

                print_if_verbose("opening multiprocessing on %i cores on %i windows"%(threads, len(df_chunk)))

                # initialize the pool
                start_time =  time.time()
                with  multiproc.Pool(threads) as pool:

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
            wholeChr_fastqgz1, wholeChr_fastqgz2 = run_wgsim_pairedEnd_per_windows_in_parallel(df_windows_chromosomes_minCov, genome, outdir_whole_chromosomes, read_length, median_insert_size, median_insert_size_sd, replace=replace, threads=threads)

            # now simulate the reads for the regions, only considering the extra_readPairs
            df_windows_sliding["readPairs"] = df_windows_sliding.extra_readPairs        
            perWindow_fastqgz1, perWindow_fastqgz2 = run_wgsim_pairedEnd_per_windows_in_parallel(df_windows_sliding, genome, outdir_per_windows, read_length, median_insert_size, median_insert_size_sd, replace=replace, threads=threads)

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
        run_bwa_mem(read1_fastqgz, read2_fastqgz, reference_genome, outdir, sim_bamfile, sim_sorted_bam, sim_index_bam, name_sample="simulations_reference_genome", threads=threads, replace=replace, MarkDuplicates=False)

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

def keep_relevant_filters_lists_inparallel(filterName_to_filtersList, df_gridss, reference_genome, type_filtering="keeping_all_filters_that_change",  wrong_INFOtags=("IMPRECISE",), min_size=50, threads=4):

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
    with  multiproc.Pool(threads) as pool:
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

def write_breakpoints_for_parameter_combinations_and_get_filterIDtoBpoints_gridss(df_gridss, df_bedpe, outdir, reference_genome, range_filtering="theoretically_meaningful", expected_AF=1.0, replace=False, threads=4):

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
            wrong_FILTERtags_l = [("NO_ASSEMBLY",), ("NO_ASSEMBLY", "INSUFFICIENT_SUPPORT"), ("NO_ASSEMBLY", "LOW_QUAL"), ("LOW_QUAL", "INSUFFICIENT_SUPPORT"), all_FILTER_tags, meaningful_FILTER_tags] 
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
            min_QUAL_l = [50, 100, 300, 500, 800, 1000, 1000000000]
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
        keep_relevant_filters_lists_inparallel(filterName_to_filtersList, df_gridss_twoBreakEnds, reference_genome, type_filtering="keeping_filters_that_yield_uniqueBPs", wrong_INFOtags=wrong_INFOtags, min_size=min_size, threads=threads) # it can also be keeping_all_filters_that_change or keeping_filters_that_yield_uniqueBPs or none

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
        ncores = threads
        print_if_verbose("Obtaining the list of tuples of breakpoints will take arround %.2f minutes on %i cores"%(((np.mean(times)*I)/ncores)/60, ncores))


        # obtain the list of tuples for each parameter combintaion
        with  multiproc.Pool(threads) as pool:
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

        with multiproc.Pool(threads) as pool:
            pool.starmap(write_bedpeANDfilterdicts_for_breakpoints, inputs_function)
            pool.close()

        # save the map between each filter 
        print_if_verbose("writing files")
        filtersID_to_breakpoints = dict(zip(bpoints_to_ID.values(), bpoints_to_ID.keys()))
        save_object(filtersID_to_breakpoints, filtersID_to_breakpoints_file)

    else: filtersID_to_breakpoints = load_object(filtersID_to_breakpoints_file)

    # return the dataframe with all the parameter combinations and the filter
    return filtersID_to_breakpoints

def benchmark_bedpe_with_knownSVs(bedpe, know_SV_dict, reference_genome, sorted_bam, median_coverage, replace=False, ID_benchmark="defaultID", delete_intermediate_files=True, threads=4, expected_AF=1.0):

    """Takes the full path to a bedpe file and generates files, under the same directory, that indicate the benchmarking.
    expected_AF is the allele frequency expected by the reads. It will determine which parameters are used."""

    # write files under the bedpe outdir
    outdir = "/".join(bedpe.split("/")[0:-1])
    bedpe_filename = bedpe.split("/")[-1]

    # get the benchmark file
    benchmark_df_filename = "%s/df_benchmarking_allParms.py"%outdir
    #remove_file(benchmark_df_filename) # debug

    if file_is_empty(benchmark_df_filename) or replace is True:
        print_if_verbose("benchmarking")

        # get the df clove
        df_clove = get_df_clove_from_bedpe(sorted_bam, bedpe, reference_genome, median_coverage, replace=replace, threads=threads, return_df=True)

        # initialize benchmark_df
        df_benchmark_all = pd.DataFrame()

        # define boundaries of the coverage parameters
        max_max_rel_coverage_to_consider_del = 1.1 - expected_AF  
        min_min_rel_coverage_to_consider_dup = 1 + expected_AF*0.6

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
        df_benchmark_all["clove_max_rel_coverage_to_consider_del"] = [max_max_rel_coverage_to_consider_del]*len(df_benchmark_all)
        df_benchmark_all["clove_min_rel_coverage_to_consider_dup"] = [min_min_rel_coverage_to_consider_dup]*len(df_benchmark_all)


        #########################################################################

        # initialize a benchmark df for tandel
        df_benchmark_TANDEL = pd.DataFrame()

        ##### deletions #################

        # get df
        df_DEL = df_clove[df_clove.SVTYPE=="DEL"]

        # go through different deletion ranges that define true deletion
        max_rel_coverage_to_consider_del_l = [0.0, 0.0001, 0.001, 0.01, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0]
        max_rel_coverage_to_consider_del_l = [x for x in max_rel_coverage_to_consider_del_l if x<=max_max_rel_coverage_to_consider_del]

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
            df_benchmark["clove_min_rel_coverage_to_consider_dup"] = [min_min_rel_coverage_to_consider_dup]*len(df_benchmark)

            # keep
            df_benchmark_TANDEL = df_benchmark_TANDEL.append(df_benchmark, sort=True)

        ################################

        ##### tandem duplications ######

        # get df
        df_TAN = df_clove[df_clove.SVTYPE=="TAN"]
        
        # go through different TAN ranges that define true TAN
        min_rel_coverage_to_consider_dup_l = [0.0, 0.4, 0.5, 1.0, 1.5, 1.8, 2.0, 2.5, 3.5]
        min_rel_coverage_to_consider_dup_l = [x for x in min_rel_coverage_to_consider_dup_l if x>=min_min_rel_coverage_to_consider_dup]

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
            df_benchmark["clove_max_rel_coverage_to_consider_del"] = [max_max_rel_coverage_to_consider_del]*len(df_benchmark)
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

        filenames_to_keep = {bedpe_filename, "%s.clove.vcf.TANDEL.bed.coverage_provided_windows.tab"%bedpe_filename, "unbalanced_translocations_5with5_or_3with3.bed.coverage_provided_windows.tab", "%s.clove.vcf"%bedpe_filename, "uniform_filters_series.py", "variable_filters_df.py", "df_benchmarking_allParms.py", "less_conservative_filtersDict.py", "%s.df_clove.tab"%bedpe_filename}

        for file in os.listdir(outdir):
            if file not in filenames_to_keep and "benchmark_analysis_" not in file: remove_file("%s/%s"%(outdir, file))


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
    """
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

        print(cmap)
        print(df_benchmark)

        sns.scatterplot(x="recall", y="precision", data=df_benchmark, hue="Fvalue", palette=cmap, edgecolors=None, style="svtype", markers=svtype_to_marker)
        #sns.scatterplot(x="recall", y="precision", data=df_benchmark, hue="Fvalue", edgecolors=None, style="svtype", markers=svtype_to_marker)


        fig.tight_layout()  # otherwise the right y-label is slightly 
        filename="%s/scatter_PRvsRCvsFvalue.pdf"%(PlotsDir)
        fig.savefig(filename, bbox_inches='tight');
        #if is_cluster is False: plt.show()
        plt.close(fig)
    """

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

def get_df_clove_from_bedpe(sorted_bam, bedpe, reference_genome, median_coverage, replace=False, threads=4, return_df=True):

    """Takes a bedpe file, runs clove on it and returns the path to the df_clove df """

    df_clove_filename = "%s.df_clove.tab"%bedpe

    if file_is_empty(df_clove_filename) or replace is True:
        print_if_verbose("getting clove df for %s"%get_file(bedpe))


        # first run clove without checking for coverage deviations
        outfile_clove = "%s.clove.vcf"%(bedpe)
        run_clove_filtered_bedpe(bedpe, outfile_clove, sorted_bam, replace=replace, median_coverage=10, median_coverage_dev=1, check_coverage=False)

        # now convert it to a df that has also the coverage for TANDEL REGIONS
        df_clove = get_clove_output_with_coverage(outfile_clove, reference_genome, sorted_bam, median_coverage, replace=replace, run_in_parallel=True, delete_bams=False, threads=threads)

        # write
        df_clove_filename_tmp = "%s.tmp"%df_clove_filename
        df_clove.to_csv(df_clove_filename_tmp, sep="\t", index=False, header=True)
        os.rename(df_clove_filename_tmp, df_clove_filename)

    # return 
    if return_df is True:

        # load the df and return it
        df_clove = pd.read_csv(df_clove_filename, sep="\t")

        return df_clove

    else: return None


def benchmark_GridssClove_for_knownSV(sample_bam, reference_genome, know_SV_dict, outdir, range_filtering="theoretically_meaningful", expected_AF=1.0, replace=False, threads=4, median_insert_size=500, median_insert_size_sd=50, mitochondrial_chromosome="mito_C_glabrata_CBS138", run_in_parallel=True):

    """Runs a benchmarking for several combinations of filters of a GridsssClove pipeline of a given bam file (sample_bam), writing files under outdir. The known SV are provided as a dictionary that maps each type of SV to a path where a table with the SVs are known .

    range_filtering indicates which type of simulation will be performed, it can be "large", "medium", "small", "single" and correlates with the range of parameters to use.
    expected_AF is the expected allele frequency, for a haploid it should be 1.0 

    median_insert_size is used to define small breakends to calculate their allele frequency

    window_l is used to define the coverage for winows of regions"""


    ###### DEFINE GENERAL THINGS

    start_time = time.time()

    # define the median coverage of regions
    print_if_verbose("getting coverage")
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir, sample_bam, windows_file="none", replace=replace, threads=threads), sep="\t")

    median_coverage = get_median_coverage(coverage_df, mitochondrial_chromosome)
    print_if_verbose("The median coverage is %i"%median_coverage)

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
            filtersID_to_breakpoints = write_breakpoints_for_parameter_combinations_and_get_filterIDtoBpoints_gridss(df_gridss, df_bedpe, outdir_parameter_combinations, reference_genome, range_filtering=range_filtering, expected_AF=expected_AF, replace=replace, threads=threads) # this is a dataframe with all the filter combinations and the map between filterID and the actual filtering

            # define the paths to the breakpoints
            paths_to_bedpe_breakpoints = ["%s/%s/filtered_breakpoints.bedpe"%(outdir_parameter_combinations, filterID) for filterID in filtersID_to_breakpoints]

            # define the threads depending on the parallel
            if run_in_parallel is True: parallel_threads = 1
            else: parallel_threads = threads

            # go through each bedpe and get the df clove
            for Ibedpe, bedpe in enumerate(paths_to_bedpe_breakpoints): 
                print_if_verbose("running clove for bedpe %i/%i"%(Ibedpe+1, len(paths_to_bedpe_breakpoints)))
                get_df_clove_from_bedpe(sample_bam, bedpe, reference_genome, median_coverage, replace=replace, threads=threads, return_df=False)

            # define inputs of the benchmarking pipeline
            inputs_benchmarking_pipeline = [(bedpe, know_SV_dict, reference_genome, sample_bam, median_coverage, replace, bedpe.split("/")[-2], True, parallel_threads, expected_AF) for bedpe in paths_to_bedpe_breakpoints]

            if run_in_parallel is True:

                # initialize the list of benchmarking dfs
                all_benchmarking_dfs = []

                # go through each chunk of ncpus
                for Ichunk, chunk_inputs_benchmarking_pipeline in enumerate(chunks(inputs_benchmarking_pipeline, threads)):
                    print_if_verbose("working on chunk %i"%Ichunk)

                    # redefine the threads, so that you have at least 4Gb of RAM per thread
                    available_RAM = get_availableGbRAM(outdir)

                    # define the maxiumum number of threads so that each thread has 8Gb of ram
                    max_threads = max([1, int(available_RAM/6 - 1)]) 
                    if threads>max_threads: threads =  max_threads 

                    # get the parallelized obtention of data
                    print_if_verbose("getting benchmarking for each set of filters in parallel on %i threads"%threads)
                    with multiproc.Pool(threads) as pool:
                        all_benchmarking_dfs += pool.starmap(benchmark_bedpe_with_knownSVs, chunk_inputs_benchmarking_pipeline)
                        pool.close()
                        pool.terminate()
            else:

                all_benchmarking_dfs = list(map(lambda x: benchmark_bedpe_with_knownSVs(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]), inputs_benchmarking_pipeline))

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
            #merge_coverage_per_window_files_in_one(sample_bam)

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

def get_df_accuracy_for_train_filer(r, outdir, test_gridss_info_dict, sorted_bam, reference_genome, median_coverage, replace, median_insert_size, median_insert_size_sd, test_SVdict, threads=4):

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
        predicted_svtype_to_SVtable, df_gridss = run_gridssClove_given_filters(sorted_bam, reference_genome, working_dir, median_coverage, replace=replace, threads=threads, gridss_blacklisted_regions=r["gridss_regionsToIgnoreBed"], gridss_VCFoutput=gridss_VCFoutput_underWorkDir, gridss_maxcoverage=r["gridss_maxcoverage"], median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, gridss_filters_dict=r["filters_dict"], tol_bp=50, run_in_parallel=True, max_rel_coverage_to_consider_del=r["clove_max_rel_coverage_to_consider_del"], min_rel_coverage_to_consider_dup=r["clove_min_rel_coverage_to_consider_dup"], replace_FromGridssRun=False)

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

def get_benchmarking_df_for_testSVs_from_trainSV_filterSets(test_SVdict, outdir, df_filters_train, test_gridss_info_dict, genomeID, ploidy, sorted_bam, reference_genome, median_coverage, median_insert_size, median_insert_size_sd, replace, threads=4):

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

        df_benchmark = pd.concat(list(df_filters_train.apply(lambda r: get_df_accuracy_for_train_filer(r, outdir, test_gridss_info_dict, sorted_bam, reference_genome, median_coverage, replace, median_insert_size, median_insert_size_sd, test_SVdict, threads=threads), axis=1)))

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

    # change the col_cluster
    if len(df.columns)==1: col_cluster = False

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

    # go through each combination
    for interesting_train_plodies in interesting_train_plodies_list:
        for interesting_train_genomeIDs in interesting_train_genomeIDs_list:
            for interesting_train_svtypes in interesting_train_svtypes_list:
                for interesting_test_plodies in interesting_test_plodies_list:
                    for interesting_test_genomeIDs in interesting_test_genomeIDs_list:
                        for interesting_test_svtypes in interesting_test_svtypes_list:        

                            # get the filtered df
                            df = df_cross_benchmark[(df_cross_benchmark.train_ploidy.isin(interesting_train_plodies)) & (df_cross_benchmark.train_genomeID.isin(interesting_train_genomeIDs)) & (df_cross_benchmark.train_svtype.isin(interesting_train_svtypes)) & (df_cross_benchmark.test_ploidy.isin(interesting_test_plodies)) & (df_cross_benchmark.test_genomeID.isin(interesting_test_genomeIDs)) & (df_cross_benchmark.test_svtype.isin(interesting_test_svtypes)) & (df_cross_benchmark.nevents>=5)]

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

    # define the parallel threads, as a function of 
    if run_in_parallel is True: parallel_threads = 1
    else: parallel_threads = threads

    # run a function in parallel that will take a genome and ploidy combination and evaluate the accuracy of all the filters in df_best_filters. first prepare input as list of tuples
    list_inputs = [(d["test_SVdict"], d["outdir"], d["df_filters_train"], d["test_gridss_info_dict"], genomeID, ploidy, d["sorted_bam"], reference_genome, d["median_coverage"], d["median_insert_size"], d["median_insert_size_sd"],  replace, parallel_threads) for (genomeID, ploidy), d in genomeIDandPlody_to_info.items()] # 0 is for debug

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

            list_cross_benchmarking_dfs = list(map(lambda x: get_benchmarking_df_for_testSVs_from_trainSV_filterSets(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12]), list_inputs))
        
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

def get_best_parameters_for_GridssClove_run(sorted_bam, reference_genome, outdir, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", nvars=100, real_bedpe_breakpoints=None, median_insert_size=250, median_insert_size_sd=0):

    """This finds the optimum parameters for running GRIDSS clove and returns them. The parameters are equivalent to the run_GridssClove_optimising_parameters function"""


    # define plots dir
    PlotsDir = "%s/plots"%outdir; make_folder(PlotsDir)

    # map each chromosome to length
    chr_to_len = get_chr_to_len(reference_genome)

    # count the length od the reads
    read_length = get_read_length(sorted_bam, threads=threads, replace=replace)
    print_if_verbose("The median read length is %i"%read_length)

    # count total number of reads
    total_nread_pairs = count_number_read_pairs(sorted_bam, replace=replace, threads=threads)
    #total_nread_pairs  = 100000 # this is to debug the simulation pipeline
    expected_coverage_per_bp = int((total_nread_pairs*read_length) / sum(chr_to_len.values())) +  1 # the expected coverage per position with pseudocount
    print_if_verbose("There are %i read pairs in your library. The expected coverage is %ix."%(total_nread_pairs, expected_coverage_per_bp))

    ###### MODELLING COVERAGE ######
    print_if_verbose("modelling coverage of the sample")

    # get a function that takes the GC content, chromosome and distance to the telomere and returns coverage. This is actually a lambda function
    outdir_coverage_calculation = "%s/coverage_per_regions%ibb"%(outdir, window_l); make_folder(outdir_coverage_calculation)
    df_coverage_train = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_coverage_calculation, sorted_bam, windows_file="none", replace=replace, threads=threads), sep="\t")

    distToTel_chrom_GC_to_coverage_fn = get_distanceToTelomere_chromosome_GCcontent_to_coverage_fn(df_coverage_train, reference_genome, outdir_coverage_calculation, expected_coverage_per_bp, mitochondrial_chromosome=mitochondrial_chromosome, replace=replace)

    print_if_verbose("coverage model obtained")

    ################################

    ############ GENERAL OPERATIONS THAT WILL BE NEEDED FOR ALL THE STEPS #####

    # the dir and genome names
    genome_dir = "/".join(reference_genome.split("/")[0:-1])
    genome_name = reference_genome.split("/")[-1].split(".")[0]


    # simulate reads for the reference if you are not only simulating haploid
    if set(simulation_ploidies)!={"haploid"}: 

        # get the info of the reference genome with predictions of coverage per window
        df_REFgenome_info = get_windows_infoDF_with_predictedFromFeatures_coverage(reference_genome, distToTel_chrom_GC_to_coverage_fn, expected_coverage_per_bp, replace=replace, threads=threads)

        outdir_ref = "%s/simulation_reference_genome_%ibp_windows"%(outdir, window_l)
        simulated_reference_bam_file = simulate_and_align_PairedReads_perWindow(df_REFgenome_info, reference_genome, reference_genome, total_nread_pairs, read_length, outdir_ref, median_insert_size, median_insert_size_sd, replace=replace, threads=threads)

    else: simulated_reference_bam_file = None

    ##############################################################################

    ################ SIMULATION PIPELINE ################ 

    # initialize a df with all the benchmarking data
    df_benchmark_all = pd.DataFrame()
    genomeID_to_knownSVdict = {}

    df_benchmark_all_file = "%s/df_benchmark_all.py"%outdir
    genomeID_to_knownSVdict_file= "%s/genomeID_to_knownSVdict.py"%outdir

    if file_is_empty(df_benchmark_all_file) or file_is_empty(genomeID_to_knownSVdict_file) or replace is True:

        # init all the ploidy dirs
        dirs_to_remove = []

        # go throigh each simulation (these are technical replicates of the pipeline)
        for simulation_ID in range(1, n_simulated_genomes+1):
            print_if_verbose("working on simulation %i"%simulation_ID)

            # get an outdir where all the simulations of this ID will be stored
            simulation_outdir = "%s/simulation_%i"%(outdir, simulation_ID); make_folder(simulation_outdir)

            # get the simulated SVs, which are an integration of 
            sim_svtype_to_svfile, rearranged_genome = simulate_SVs_in_genome(reference_genome, mitochondrial_chromosome, simulation_outdir, nvars=nvars, bedpe_breakpoints=real_bedpe_breakpoints, replace=replace, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"})

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
                dirs_to_remove.append("%s/benchmark_max50000x_ignoreRegionsFalse/several_parameter_combinations_filter_%s_af%.2f"%(ploidy_dir, range_filtering_benchmark, fraction_var))

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

        # clean the too-many files
        print_if_verbose("removing unnecessary files")
        for d in dirs_to_remove: delete_folder(d)

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

def run_GridssClove_optimising_parameters(sorted_bam, reference_genome, outdir, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", nvars=100, fast_SVcalling=False, real_bedpe_breakpoints=None, gridss_VCFoutput="", replace_FromGridssRun_final_perSVade_run=False):

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
    - real_bedpe_breakpoints is a bedpe file with the real breakpoints. We will simulate SVs arround them
    - gridss_VCFoutput is passed to run_gridssClove_given_filters. It can be useful if you don't want to rerun gridss, which is very expensive.

    """

    # make the outdir
    make_folder(outdir)

    # define final dirs
    outdir_gridss_final = "%s/final_gridss_running"%outdir; make_folder(outdir_gridss_final)
    final_file = "%s/gridssClove_finished.txt"%outdir

    if file_is_empty(final_file) or replace is True or replace_FromGridssRun_final_perSVade_run is True:

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

            gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup, df_cross_benchmark_best = get_best_parameters_for_GridssClove_run(sorted_bam, reference_genome, parameter_optimisation_dir, threads=threads, replace=replace, n_simulated_genomes=n_simulated_genomes, mitochondrial_chromosome=mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=range_filtering_benchmark, nvars=nvars, real_bedpe_breakpoints=real_bedpe_breakpoints, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd)

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

        # write the parameters of the running
        json_file = "%s/perSVade_parameters.json"%outdir_gridss_final
        write_gridss_parameters_as_json(gridss_blacklisted_regions, gridss_maxcoverage, gridss_filters_dict, max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup, json_file, replace=replace)

        # define the final vcf dir
        final_gridss_vcf = "%s/output_gridss.vcf"%outdir_gridss_final

        # if there is a provided gridss_VCFoutput, softlink it to the outdir_gridss_final
        if not file_is_empty(gridss_VCFoutput) and file_is_empty(final_gridss_vcf): soft_link_files(gridss_VCFoutput, final_gridss_vcf)

        # define the median coverage across window_l windows of the genome
        coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_gridss_final, sorted_bam, windows_file="none", replace=replace, threads=threads), sep="\t")
        median_coverage = get_median_coverage(coverage_df, mitochondrial_chromosome)
        print_if_verbose("The median coverage is %i"%median_coverage)

        # run the pipeline
        print_if_verbose("running final gridss with parameters...")
        final_sv_dict, df_gridss = run_gridssClove_given_filters(sorted_bam, reference_genome, outdir_gridss_final, median_coverage, replace=replace, threads=threads, gridss_blacklisted_regions=gridss_blacklisted_regions, gridss_VCFoutput=final_gridss_vcf, gridss_maxcoverage=gridss_maxcoverage, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, gridss_filters_dict=gridss_filters_dict, run_in_parallel=True, max_rel_coverage_to_consider_del=max_rel_coverage_to_consider_del, min_rel_coverage_to_consider_dup=min_rel_coverage_to_consider_dup, replace_FromGridssRun=replace_FromGridssRun_final_perSVade_run)

        # not debug: replace_FromGridssRun=replace

        ########################################

        ##### PIPELINE ENDING OPERATIONS ##### 

        # at the end, remove all the mosdepth and windows files under the reference
        clean_reference_genome_windows_files(reference_genome)
        
        print_if_verbose("GRIDSS pipeline finished correctly")

        print_if_verbose("--- the gridss pipeline optimising parameters took %s seconds in %i cores ---"%(time.time() - pipeline_start_time, threads))

        # generate a file that indicates whether the gridss run is finished
        open(final_file, "w").write("gridssClove_finished finished...")

        ######################################

    return outdir_gridss_final


def generate_jobarray_file(jobs_filename, name):
    
    """ 
    This function takes a jobs filename and replaces it creating the necessary STD files. These will be set to be run in a cluster environment
    """

    # define the stddir
    outdir = get_dir(jobs_filename)
    stddir = "%s/STDfiles"%outdir; 
    delete_folder(stddir)
    make_folder(stddir)

    # remove previous rst files
    name_jobs_filename = get_file(jobs_filename)
    for file in os.listdir(get_dir(jobs_filename)): 

        if file.startswith("%s-"%name_jobs_filename) and file.endswith(".rst"): 
            remove_file("%s/%s"%(get_dir(jobs_filename), file))

    # rewrite the jobs_filename so that each std goes to a different file
    std_perJob_prefix = "%s/%s"%(stddir, name)
    jobs_filename_lines = ["%s > %s.%i.out 2>&1"%(l.strip(), std_perJob_prefix, I+1) for I, l in enumerate(open(jobs_filename, "r").readlines())]
    open(jobs_filename, "w").write("\n".join(jobs_filename_lines))

    print("You need to successfully run all jobs in %s to continue"%jobs_filename)

   

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

def report_accuracy_realSVs_perSVadeRuns(close_shortReads_table, reference_genome, outdir, real_bedpe_breakpoints, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", nvars=100, job_array_mode="local", skip_cleaning_simulations_files_and_parameters=False, skip_cleaning_outdir=False, parameters_json_file=None, gff=None, replace_FromGridssRun_final_perSVade_run=False, fraction_available_mem=None, skip_CNV_calling=False, outdir_finding_realVars=None, replace_SV_CNVcalling=False):


    """This function runs the SV pipeline for all the datasets in close_shortReads_table with the fastSV, optimisation based on uniform parameters and optimisation based on realSVs (specified in real_svtype_to_file). The latter is skipped if real_svtype_to_file is empty.

    First, it runs perSVade on all parameters keeping some important files. It returns a dict with the outdir of each configuration."""

    ##### DEFINE INPUTS ######

    # this pipeline requires real data and close_shortReads_table that is not none
    if real_bedpe_breakpoints is None: raise ValueError("You need real data if you want to test accuracy")
    if file_is_empty(close_shortReads_table): raise ValueError("You need real data reads if you want to test accuracy")

    # make the outdir
    make_folder(outdir)

    print_if_verbose("testing the accuracy of perSVade. Running perSVade on each sample with each configuration")

    # get the gff with biotype
    if gff is not None: 
        correct_gff, gff_with_biotype = get_correct_gff_and_gff_with_biotype(gff, replace=replace)
        gff = gff_with_biotype

    # load the real data table
    df_reads = pd.read_csv(close_shortReads_table, sep="\t").set_index("runID", drop=False)

    # init final dict
    final_dict = {}

    # map each runID to the bam file, if it exists, from the previous run
    runID_to_previous_bam = {}
    for runID in set(df_reads.runID):
        if outdir_finding_realVars is None: runID_to_previous_bam[runID] = None
        else:  
            runID_to_previous_bam[runID] = "%s/all_realVars/shortReads_realVarsDiscovery_%s/aligned_reads.bam.sorted"%(outdir_finding_realVars, runID)
            for f in [runID_to_previous_bam[runID], "%s.bai"%runID_to_previous_bam[runID]]:
                if file_is_empty(f): raise ValueError("%s should exist"%f)

    ##########################

    ##### RUN JOBS ######

    # initialize the cmds to run 
    all_cmds = []

    # predefine if some jobs need to be ran
    n_remaining_jobs = sum([sum([file_is_empty("%s/%s/%s/perSVade_finished_file.txt"%(outdir, typeSimulations, runID)) for runID in set(df_reads.runID)]) for typeSimulations in ["uniform", "fast", "realSVs"]])
    print_if_verbose("There are %i remaining jobs"%n_remaining_jobs)

    # go through each run and configuration
    for typeSimulations, bedpe_breakpoints, fast_SVcalling in [("uniform", None, False), ("realSVs", real_bedpe_breakpoints, False), ("fast", None, True)]:

        # define an outdir for this type of simulations
        outdir_typeSimulations = "%s/%s"%(outdir, typeSimulations); make_folder(outdir_typeSimulations)

        # go though each runID
        for runID in set(df_reads.runID):
            print_if_verbose(typeSimulations, runID)

            # define an outdir for this runID
            outdir_runID = "%s/%s"%(outdir_typeSimulations, runID); make_folder(outdir_runID)

            # keep
            final_dict.setdefault(typeSimulations, {}).setdefault(runID, outdir_runID)

            # map the previous bam file to here, to spare running time
            previous_bam = runID_to_previous_bam[runID]
            if previous_bam is not None:
                dest_bam = "%s/aligned_reads.bam.sorted"%outdir_runID
                soft_link_files(previous_bam, dest_bam)
                soft_link_files(previous_bam+".bai", dest_bam+".bai")

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
                cmd = "python %s -r %s --threads %i --outdir %s --nvars %i --nsimulations %i --simulation_ploidies %s --range_filtering_benchmark %s --mitochondrial_chromosome %s -f1 %s -f2 %s --previous_repeats_table %s --min_CNVsize_coverageBased %i --skip_cleaning_outdir --skip_SV_CNV_calling"%(perSVade_py, reference_genome, threads, outdir_runID, nvars, n_simulated_genomes, ",".join(simulation_ploidies), range_filtering_benchmark, mitochondrial_chromosome, r1, r2, previous_repeats_table, min_CNVsize_coverageBased)

                # add arguments depending on the pipeline
                if replace is True: cmd += " --replace"
                if fast_SVcalling is True: cmd += " --fast_SVcalling"
                if bedpe_breakpoints is not None: cmd += " --real_bedpe_breakpoints %s"%bedpe_breakpoints
                if printing_verbose_mode is True: cmd += " --verbose"
                if parameters_json_file is not None: cmd += " --parameters_json_file %s"%parameters_json_file
                if gff is not None: cmd += " --gff %s"%gff
                if replace_FromGridssRun_final_perSVade_run is True: cmd += " --replace_FromGridssRun_final_perSVade_run"
                if fraction_available_mem is not None: cmd += " --fraction_available_mem %.3f"%(float(fraction_available_mem))
                if replace_SV_CNVcalling is True: cmd += " --replace_SV_CNVcalling"
                if skip_CNV_calling is True: cmd += " --skip_CNV_calling"


                # if the running in slurm is false, just run the cmd
                if job_array_mode=="local": run_cmd(cmd)
                elif job_array_mode=="job_array": 
                    all_cmds.append(cmd)
                    continue

                else: raise ValueError("%s is not valid"%job_array_mode)

            # keep the simulation files and clean outdir
            elif n_remaining_jobs==0:

                # keeping simulations and cleaning
                keep_simulation_files_for_perSVade_outdir(outdir_runID, replace=replace, n_simulated_genomes=n_simulated_genomes, simulation_ploidies=simulation_ploidies)

                print(outdir_runID)
                stopbeforecleaningtotestthatsimlationkeepingworked_ThisIsTOCkeckThatAllFIleswwereCorrect

                # clean
                clean_perSVade_outdir(outdir_runID)

    # if you are not running on slurm, just execute one cmd after the other
    if job_array_mode=="job_array":

        if len(all_cmds)>0: 
            print_if_verbose("submitting %i jobs to the cluster for testing accuracy of perSVade on several combinations of parameters. The files of the submission are in %s"%(len(all_cmds), outdir))
            jobs_filename = "%s/jobs.testingRealDataAccuracy"%outdir
            open(jobs_filename, "w").write("\n".join(all_cmds))

            generate_jobarray_file(jobs_filename, "accuracyRealSVs")

            print_if_verbose("You have to wait under all the jobs in testRealSVs are done")
            sys.exit(0)

    #####################

    return final_dict


def keep_simulation_files_for_perSVade_outdir(perSVade_outdir, replace=False, n_simulated_genomes=2, simulation_ploidies=["haploid", "diploid_homo"]):

    """Takes an outdir of perSVade and saves the simulation files for further use"""

    print_if_verbose("keeping simulation files")

    # define the dirs
    svDetection_dir = "%s/SVdetection_output"%perSVade_outdir

    # define the outdir of the simulation files
    outdir = "%s/simulations_files_and_parameters"%perSVade_outdir
    if replace is True: delete_folder(outdir)
    make_folder(outdir)

    # get the parameters json file
    parameters_json_origin = "%s/final_gridss_running/perSVade_parameters.json"%svDetection_dir
    parameters_json_dest = "%s/final_parameters.json"%(outdir)
    rsync_file(parameters_json_origin, parameters_json_dest)

    # check that there are some simulations done, if not, just return
    parameter_optimisation_dir = "%s/parameter_optimisation"%svDetection_dir
    if not os.path.isdir(parameter_optimisation_dir): return

    # init a dict with the metadata of the simulations
    simulations_metadata_dict = {}; I=0

    # keep data for each optimisation
    for Isim in range(n_simulated_genomes):

        # define the name
        simName = "simulation_%i"%(Isim+1)
        sim_outdir = "%s/%s"%(parameter_optimisation_dir, simName)

        # go through each ploidy
        for ploidy in simulation_ploidies:

            # define the destintaion bam 
            destination_bam = "%s/reads_sim%i_%s.bam"%(outdir, Isim+1, ploidy)

            # define the origin bam
            if ploidy=="haploid": suffix_ploidy = "bam.sorted"
            else: suffix_ploidy = "bam.sorted.%s.bam.sorted"%ploidy
            bam_files = ["%s/%s"%(sim_outdir, file) for file in os.listdir(sim_outdir) if file.startswith("aligned_reads") and ".".join(file.split(".")[1:])==suffix_ploidy]
            if len(bam_files)!=1: raise ValueError("There should be only one bam")
            origin_bam = bam_files[0]

            # copy the metrics
            origin_metrics = "%s/aligned_reads.bam.sorted.CollectInsertSizeMetrics.out"%perSVade_outdir
            destination_metrics = "%s.CollectInsertSizeMetrics.out"%destination_bam
            rsync_file(origin_metrics, destination_metrics)

            # change the bai
            rsync_file("%s.bai"%origin_bam, "%s.bai"%destination_bam)

            # change the coverage per window (this is any destination)
            rsync_file("%s.coverage_per_window.tab"%origin_bam, "%s.coverage_per_window.tab"%destination_bam)

            # change the coverage per constant windows
            calculating_windowcoverage_dir = "%s.calculating_windowcoverage"%destination_bam; make_folder(calculating_windowcoverage_dir)
            destination_windowcoverage_file = "%s/coverage_windows_%ibp.tab"%(calculating_windowcoverage_dir, window_l)
            origin_windowcoverage_file = "%s/benchmark_GridssClove_%s/coverage_windows_%ibp.tab"%(sim_outdir, ploidy, window_l)
            rsync_file(origin_windowcoverage_file, destination_windowcoverage_file)

            # rename the bam
            rsync_file(origin_bam, destination_bam)

            # change the gridss vcf
            origin_gridss_vcf = "%s/benchmark_GridssClove_%s/benchmark_max50000x_ignoreRegionsFalse/gridss_output.vcf.withSimpleEventType.vcf"%(sim_outdir, ploidy)
            dest_gridss_vcf = "%s/gridss_vcf_sim%i_%s.vcf"%(outdir, Isim+1, ploidy)
            rsync_file(origin_gridss_vcf, dest_gridss_vcf)

            # change the location of the simulated SVs
            svtables_prefix =  "%s/SVs_sim%i"%(outdir, Isim+1)
            for svtype in {"insertions", "deletions", "translocations", "inversions", "tandemDuplications"}:

                # define the files
                origin_file = "%s/final_simulated_SVs/%s.tab"%(sim_outdir, svtype)
                dest_file = "%s_%s.tab"%(svtables_prefix, svtype)
                rsync_file(origin_file, dest_file)

            # keep metadata
            simulations_metadata_dict[I] = {"simName":simName, "ploidy":ploidy, "sorted_bam":get_file(destination_bam), "gridss_vcf":get_file(dest_gridss_vcf), "svtables_prefix":get_file(svtables_prefix)}; I+=1

    # save the content
    simulations_metadata_df = pd.DataFrame(simulations_metadata_dict).transpose()
    save_df_as_tab(simulations_metadata_df, "%s/directory_content.tab"%outdir)

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


def generate_final_file_report(final_file, start_time_GeneralProcessing, end_time_GeneralProcessing, start_time_alignment, end_time_alignment, start_time_all, end_time_all, start_time_obtentionCloseSVs, end_time_obtentionCloseSVs, start_time_SVcalling, end_time_SVcalling, start_time_SVandCNVcalling, end_time_SVandCNVcalling, start_time_smallVarsCNV, end_time_smallVarsCNV):

    """Generates the final file of perSVade. This is a tab file with many relevant info"""

    # generate the lines
    lines = ["perSVade finished...",
             "This is the time in seconds that each part of the pipeline took:",
             "time_GeneralProcessing:%s"%(end_time_GeneralProcessing-start_time_GeneralProcessing),
             "time_alignment:%s"%(end_time_alignment-start_time_alignment),
             "time_all:%s"%(end_time_all-start_time_all),
             "time_obtentionCloseSVs:%s"%(end_time_obtentionCloseSVs-start_time_obtentionCloseSVs),
             "time_SVcalling:%s"%(end_time_SVcalling-start_time_SVcalling),
             "time_SVandCNVcalling:%s"%(end_time_SVandCNVcalling-start_time_SVandCNVcalling),
             "time_smallVarsCNV:%s"%(end_time_smallVarsCNV-start_time_smallVarsCNV)
             ]

    # write
    open(final_file, "w").write("\n".join(lines))


def report_accuracy_golden_set_runJobs(goldenSet_dir, outdir, reference_genome, real_bedpe_breakpoints, threads=4, replace=False, n_simulated_genomes=2, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], range_filtering_benchmark="theoretically_meaningful", nvars=100, job_array_mode="local", StopAfter_sampleIndexingFromSRA=False, StopAfterPrefecth_of_reads=False, target_taxID=None, parameters_json_file=None, fraction_available_mem=None, StopAfter_goldenSetAnalysis_readObtention=False, verbose=False, StopAfter_goldenSetAnalysis_readTrimming=False, min_coverage=30):


    """This function takes a directory that has the golden set vars and generates plots reporting the accuracy. If auto, it will find them in the SRA and write them under outdir."""

    print_if_verbose("calculating accuracy for golden set SVcalls")
    make_folder(outdir)

    ##########################
    ####### GET READS ########
    ##########################

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

        if StopAfter_sampleIndexingFromSRA:
            print("Golden set analysis. Stop after indexing from SRA")
            sys.exit(0)

        #####################################

        # define the reads
        longReads = "%s/%s/%s.srr.fastq.gz"%(goldenSet_dir, long_reads_srr, long_reads_srr)
        short_reads1 = "%s/%s/%s.srr_1.fastq.gz"%(goldenSet_dir, short_reads_srr, short_reads_srr)
        short_reads2 = "%s/%s/%s.srr_2.fastq.gz"%(goldenSet_dir, short_reads_srr, short_reads_srr)

        # download each of the reads (raw). Stop after fastqdump
        for type_data, srr in [("illumina_paired", short_reads_srr), ("nanopore", long_reads_srr)]:
            print_if_verbose("Getting raw reads for %s"%type_data)

            # define the outdir
            outdir_srr = "%s/%s"%(goldenSet_dir, srr)

            # define the cmd downloading after the fastq-dump
            cmd = "%s --srr %s --outdir %s --threads %i --stop_after_fastqdump --type_data %s"%(get_trimmed_reads_for_srr_py, srr, outdir_srr, threads, type_data)
            if StopAfterPrefecth_of_reads is True: cmd += " --stop_after_prefetch"

            run_cmd(cmd)

    #####################################
    else:

        # define the reads, they are suposed to be called like this
        origin_longReads = "%s/long_reads.fastq.gz"%goldenSet_dir
        origin_short_reads1 = "%s/short_reads_1.fastq.gz"%goldenSet_dir
        origin_short_reads2 = "%s/short_reads_2.fastq.gz"%goldenSet_dir

        # copy under provided_goldenSetReads/
        provided_goldenSetReads_dir = "%s/provided_goldenSetReads"%(outdir); make_folder(provided_goldenSetReads_dir)
        longReads = "%s/long_reads.fastq.gz"%provided_goldenSetReads_dir
        short_reads1 = "%s/short_reads_1.fastq.gz"%provided_goldenSetReads_dir
        short_reads2 = "%s/short_reads_2.fastq.gz"%provided_goldenSetReads_dir

        soft_link_files(origin_longReads, longReads)
        soft_link_files(origin_short_reads1, short_reads1)
        soft_link_files(origin_short_reads2, short_reads2)

    # debug
    if any([StopAfter_sampleIndexingFromSRA, StopAfterPrefecth_of_reads]): 
        print("Golden set analysis. Exiting after sample Indexing or prefetch")
        sys.exit(0)

    if any([file_is_empty(f) for f in [longReads, short_reads1, short_reads2]]): raise ValueError("Your golden dir %s should contain long_reads.fasta, short_reads_1.fastq.gz and short_reads_2.fastq.gz"%goldenSet_dir)

    if StopAfter_goldenSetAnalysis_readObtention:
        print("Golden set analysis. Exiting after sample Indexing or prefetch")
        sys.exit(0)


    ##########################
    ##########################
    ##########################

    #########################
    ####### RUN JOBS ########
    #########################

    # init the cmds
    all_cmds = []

    # add the run svim and sniffles job
    outdir_ONT_calling = "%s/ONT_SV_calling"%outdir; make_folder(outdir_ONT_calling)
    final_file_ONT_calling = "%s/ONT_SV_calling_finished.txt"%outdir_ONT_calling
    if file_is_empty(final_file_ONT_calling) or replace is True: 

        ont_calling_cmd = "%s --ref %s --input_reads %s --outdir %s --aligner ngmlr --threads %i"%(run_svim_and_sniffles_py, reference_genome, longReads, outdir_ONT_calling, threads)
        if replace is True: ont_calling_cmd += " --replace"
        if verbose is True: ont_calling_cmd += " --verbose"

        all_cmds.append(ont_calling_cmd)

    # add the perSVade runs in several combinations
    n_remaining_jobs = sum([file_is_empty("%s/perSVade_calling_%s"%(outdir, typeSimulations)) for typeSimulations in ["uniform", "fast", "realSVs"]])

    for typeSimulations, bedpe_breakpoints, fast_SVcalling in [("uniform", None, False), ("realSVs", real_bedpe_breakpoints, False), ("fast", None, True)]:

        # define an outdir for this type of simulations
        outdir_typeSimulations = "%s/perSVade_calling_%s"%(outdir, typeSimulations); make_folder(outdir_typeSimulations)

        # define the final file 
        final_file = "%s/perSVade_finished_file.txt"%outdir_typeSimulations

        # define the previous repeats file 
        previous_repeats_table = "%s.repeats.tab"%reference_genome
        if file_is_empty(previous_repeats_table): raise ValueError("%s should exist"%previous_repeats_table)
            
        # only contine if the final file is not defined
        if file_is_empty(final_file) or replace is True:

            # define the cmd. This is a normal perSvade.py run with the vars of the previous dir  
            cmd = "python %s -r %s --threads %i --outdir %s --nvars %i --nsimulations %i --simulation_ploidies %s --range_filtering_benchmark %s --mitochondrial_chromosome %s -f1 %s -f2 %s --previous_repeats_table %s --skip_cleaning_outdir --skip_SV_CNV_calling --QC_and_trimming_reads"%(perSVade_py, reference_genome, threads, outdir_typeSimulations, nvars, n_simulated_genomes, ",".join(simulation_ploidies), range_filtering_benchmark, mitochondrial_chromosome, short_reads1, short_reads2, previous_repeats_table)

            # add arguments depending on the pipeline
            if replace is True: cmd += " --replace"
            if fast_SVcalling is True: cmd += " --fast_SVcalling"
            if bedpe_breakpoints is not None: cmd += " --real_bedpe_breakpoints %s"%bedpe_breakpoints
            if printing_verbose_mode is True: cmd += " --verbose"
            if parameters_json_file is not None: cmd += " --parameters_json_file %s"%parameters_json_file
            if fraction_available_mem is not None: cmd += " --fraction_available_mem %.3f"%(float(fraction_available_mem))

            all_cmds.append(cmd)

        # keep the simulation files and clean outdir
        elif n_remaining_jobs==0:

            # keeping simulations and cleaning
            keep_simulation_files_for_perSVade_outdir(outdir_typeSimulations, replace=replace, n_simulated_genomes=n_simulated_genomes, simulation_ploidies=simulation_ploidies)

            print(outdir_typeSimulations)
            stopbeforecleaningtotestthatsimlationkeepingworked_ThisIsTOCkeckThatAllFIleswwereCorrect_GoldenSetTesting

            # clean
            clean_perSVade_outdir(outdir_runID)

    # run jobs
    if job_array_mode=="job_array":

        if len(all_cmds)>0: 
            print_if_verbose("submitting %i jobs to the cluster for testing accuracy of perSVade on short and long reads (Golden set analysis). The files of the submission are in %s"%(len(all_cmds), outdir))
            jobs_filename = "%s/jobs.GoldenSetTesting"%outdir
            open(jobs_filename, "w").write("\n".join(all_cmds))

            generate_jobarray_file(jobs_filename, "accuracyGoldenSet")

            print_if_verbose("You have to wait under all the jobs in testRealSVs are done")
            sys.exit(0)

    elif job_array_mode=="local":
        for cmd in all_cmds: run_cmd(cmd)

    else: raise ValueError("%s is not valid"%job_array_mode)
    
    #########################
    #########################
    #########################

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

                             "variants_atLeast1PASS_ploidy%i_alternative_genome.fasta"%ploidy,
                             "variants_atLeast2PASS_ploidy%i_alternative_genome.fasta"%ploidy,
                             "variants_atLeast3PASS_ploidy%i_alternative_genome.fasta"%ploidy,

                             "variants_atLeast1PASS_ploidy%i.withMultiAlt.vcf"%ploidy,
                             "variants_atLeast2PASS_ploidy%i.withMultiAlt.vcf"%ploidy,
                             "variants_atLeast3PASS_ploidy%i.withMultiAlt.vcf"%ploidy,

                             "variant_calling_stats_ploidy%i_called.tab"%ploidy,
                             "variant_calling_stats_ploidy%i_PASS.tab"%ploidy

                             }

            if f not in files_to_keep: files_to_remove.append(file)

    for f in files_to_remove: remove_file(f)



def remove_smallVarsCNV_nonEssentialFiles_severalPloidies(outdir, ploidies):

    """Removes all the files in outdir that are not essential. The outdir has to be the one of the VarCall outdir."""

    print_if_verbose("cleaning outdir of small vars")

    # initialize the files to remove
    files_to_remove = ["%s/CNV_results/gene_to_coverage_genes.tab"%outdir, # the genes coverage
                       "%s/CNV_results/gene_to_coverage_regions.tab"%outdir # the regions coverage
                       ]

    # init the files to keep
    files_to_keep = set()

    # add the files from the small vars
    for ploidy in ploidies:

        # add the bcftools
        bcftools_dir = "%s/bcftools_ploidy%i_out"%(outdir, ploidy)
        HC_dir = "%s/HaplotypeCaller_ploidy%i_out"%(outdir, ploidy)
        fb_dir = "%s/freebayes_ploidy%i_out"%(outdir, ploidy)

        # go through each dir
        for vcfDir in [bcftools_dir, HC_dir, fb_dir]:

            if os.path.isdir(vcfDir):
                for file in os.listdir(vcfDir):

                    if file not in {"output.raw.vcf", "output.filt.vcf"}: files_to_remove.append("%s/%s"%(vcfDir, file))

        files_to_keep.update({"merged_vcfs_allVars_ploidy%i.vcf"%ploidy,
                              "variant_annotation_ploidy%i.tab"%ploidy,
                              "variant_calling_ploidy%i.tab"%ploidy,

                             "variants_atLeast1PASS_ploidy%i.vcf"%ploidy,
                             "variants_atLeast2PASS_ploidy%i.vcf"%ploidy,
                             "variants_atLeast3PASS_ploidy%i.vcf"%ploidy,

                             "variants_atLeast1PASS_ploidy%i_alternative_genome.fasta"%ploidy,
                             "variants_atLeast2PASS_ploidy%i_alternative_genome.fasta"%ploidy,
                             "variants_atLeast3PASS_ploidy%i_alternative_genome.fasta"%ploidy,

                             "variants_atLeast1PASS_ploidy%i.withMultiAlt.vcf"%ploidy,
                             "variants_atLeast2PASS_ploidy%i.withMultiAlt.vcf"%ploidy,
                             "variants_atLeast3PASS_ploidy%i.withMultiAlt.vcf"%ploidy,

                             "variant_calling_stats_ploidy%i_called.tab"%ploidy,
                             "variant_calling_stats_ploidy%i_PASS.tab"%ploidy})

    # go through the files in the outdir and just keep the essential ones
    for f in os.listdir(outdir):
        file = "%s/%s"%(outdir, f)

        # only set to remove if the file is not in files_to_keep
        if os.path.isfile(file) and f not in files_to_keep: files_to_remove.append(file)

    for f in files_to_remove: remove_file(f)


def get_bam_with_duplicatesMarkedSpark(bam, threads=4, replace=False, remove_duplicates=False):

    """
    This function takes a bam file and runs MarkDuplicatesSpark (most efficient when the input bam is NOT coordinate-sorted) returning the bam with the duplicates sorted. It does not compute metrics to make it faster. Some notes about MarkDuplicatesSpark:

    - it is optimised to run on query-name-sorted bams. It is slower if not
    - by default it does not calculate duplicate metrics
    - for 30x coverage WGS, it is recconnended to have at least 16Gb
    - by default it takes all the cores

    it reurns the bam with duplicates removed

    """

    # define dirs
    bam_dupMarked = "%s.MarkDups.bam"%bam
    bam_dupMarked_tmp = "%s.MarkDups.tmp.bam"%bam

    # calculate availableGbRAM
    availableGbRAM = get_availableGbRAM(get_dir(bam))

    # different memory options
    javaRamGb_common = int(availableGbRAM*fractionRAM_to_dedicate) # At Broad, we run MarkDuplicates with 2GB Java heap (java -Xmx2g) and 10GB hard memory limit
    #javaRamGb_half = int(availableGbRAM*0.5) # At Broad, we run MarkDuplicates with 2GB Java heap (java -Xmx2g) and 10GB hard memory limit
    #javaRamGb_allBut2 = int(availableGbRAM - 2) # rule of thumb from GATK
    #javaRamGb_4 = 4 # this is from a post from 2011, reccommended for a 170Gb RAM

    javaRamGb = javaRamGb_common

    if file_is_empty(bam_dupMarked) or replace is True:
        print_if_verbose("marking duplicate reads")

        # define the java memory
        
        # define the tmpdir
        tmpdir = "%s.runningMarkDups_tmp"%bam
        delete_folder(tmpdir)
        make_folder(tmpdir)

        # remove any files with .tmp.
        for f in os.listdir(get_dir(bam)):
            if f.startswith(get_file(bam_dupMarked_tmp)): 
                path = "%s/%s"%(get_dir(bam), f)
                remove_file(path)
                delete_folder(path)

        # define the remove_duplicates options
        if remove_duplicates is True: remove_duplicates_str = "true"
        else: remove_duplicates_str = "false"

        markduplicates_std = "%s.markingDuplicates.std"%bam
        print_if_verbose("running MarkDuplicates with %iGb of RAM. The std is in %s"%(javaRamGb, markduplicates_std))

        try: 

            run_cmd("%s --java-options '-Xms%ig -Xmx%ig' MarkDuplicatesSpark -I %s -O %s --verbosity INFO --tmp-dir %s  --create-output-variant-index false --create-output-bam-splitting-index false --create-output-bam-index false --remove-all-duplicates %s > %s 2>&1"%(gatk, javaRamGb, javaRamGb, bam, bam_dupMarked_tmp, tmpdir, remove_duplicates_str, markduplicates_std))

            print_if_verbose("MarkDuplicatesSpark worked correctly")

        except: 

            # define the MarkDuplicatesSpark options
            MarkDuplicatesSpark_log = "%s.MarkDuplicatesSpark_failed.log"%bam
            os.rename(markduplicates_std, MarkDuplicatesSpark_log)



            print_if_verbose("MarkDuplicatesSpark did not work on the current memory configuration. Trying with the normal MarkDuplicates. You can check the failed log of MarkDuplicatesSpark_failed in %s"%MarkDuplicatesSpark_log)

            # running with the traditional MarkDuplicates implementation
            bam_dupMarked_metrics = "%s.metrics.txt"%bam_dupMarked_tmp
            run_cmd("%s -Xms%ig -Xmx%ig MarkDuplicates I=%s O=%s M=%s ASSUME_SORT_ORDER=queryname REMOVE_DUPLICATES=%s > %s 2>&1"%(picard_exec, javaRamGb, javaRamGb, bam, bam_dupMarked_tmp, bam_dupMarked_metrics, remove_duplicates_str, markduplicates_std), env=EnvName_picard)
            
            remove_file(bam_dupMarked_metrics)

        remove_file(markduplicates_std)
        delete_folder(tmpdir)

        # keep
        os.rename(bam_dupMarked_tmp, bam_dupMarked)

    return bam_dupMarked

def get_chromosomal_sorted_bam(sorted_bam, sorted_bam_chr, chromosome, replace, threads):

    """Gets the chromosomal bam and index for sorted_bam"""

    if file_is_empty(sorted_bam_chr) or replace is True:

        sorted_bam_chr_stderr = "%s.generating.stderr"%sorted_bam_chr
        sorted_bam_chr_tmp = "%s.tmp"%sorted_bam_chr

        run_cmd("%s view -b %s %s > %s 2>%s"%(samtools, sorted_bam, chromosome, sorted_bam_chr_tmp, sorted_bam_chr_stderr))

        remove_file(sorted_bam_chr_stderr)

        # index the bam
        index_bam(sorted_bam_chr_tmp, threads=threads)

        # renames
        os.rename("%s.bai"%sorted_bam_chr_tmp, "%s.bai"%sorted_bam_chr)
        os.rename(sorted_bam_chr_tmp, sorted_bam_chr)

def get_mpileup_file_one_chrom(sorted_bam, replace, reference_genome, min_basecalling_qual=30, min_map_qual=30):

    """Runs samtools mpileup for the sorted bam, returning the file that has the position and the read depth. The bam should be only one chromosome"""

    mpileup_file = "%s.mpileup.txt"%sorted_bam

    if file_is_empty(mpileup_file) or replace is True:

        mpileup_file_tmp = "%s.tmp"
        mpileup_stderr = "%s.generating.stderr"%mpileup_file

        run_cmd("%s mpileup -f %s --min-MQ %i --min-BQ %i -a %s | cut -f2,4 > %s 2>%s "%(samtools, reference_genome, min_map_qual, min_basecalling_qual, sorted_bam, mpileup_file_tmp, mpileup_stderr))

        remove_file(mpileup_stderr)
        os.rename(mpileup_file_tmp, mpileup_file)

    return mpileup_file




#########################
####### HMMCOPY #########
#########################

def run_HMMcopy(coverage_file, outfile, parms_dict=default_parms_dict_HMMcopy, replace=False):

    """Runs HMMcopy with the different paramters into outfile. It returns the df with the outfile."""

    if file_is_empty(outfile) or replace is True:

        # define tmp
        outfile_tmp = "%s.tmp"%outfile
        std = "%s.std"%outfile

        # redefine parms
        mu = ",".join([str(x) for x in parms_dict["mu"]])
        m = ",".join([str(x) for x in parms_dict["m"]])
        kappa = ",".join([str(x) for x in parms_dict["kappa"]])
        e = parms_dict["e"]
        lambda_val = parms_dict["lambda_val"]
        nu = parms_dict["nu"]
        eta = parms_dict["eta"]
        gamma = parms_dict["gamma"]
        S = parms_dict["S"]
        strength = parms_dict["strength"]

        # run cmd
        cmd = "%s --outfile %s --coverage_table %s --e %s --mu %s --lambda %s --nu %s --kappa %s --m %s --eta %s --gamma %s --S %s --strength %s > %s 2>&1"%(run_HMMcopy_R, outfile_tmp, coverage_file, e, mu, lambda_val,  nu, kappa, m, eta, gamma, S, strength, std)

        run_cmd(cmd, env=EnvName_HMMcopy)

        # clean
        remove_file(std)
        os.rename(outfile_tmp, outfile)

    # get the df and return
    df = get_tab_as_df_or_empty_df(outfile)

    # add the CN as in CONY. Note that HMMcopy gets as NaNs the regions with 0 coverage. I set them as 0s
    state_to_CN = {1:0, 2:1, 3:2, 4:3, 5:4, 6:5}

    def get_CN_as_0(r):
        if r["reads"]==0: return 0
        else: return state_to_CN[r["state_CNsegment"]]

    df["CN"] = df.apply(get_CN_as_0, axis=1)

    # add the relative CN
    df["relative_CN"] = (df.CN / 2)
 
    return df

def get_rsquare_one_chromosome_crossValidated_rsquare_HMMcopy(df_train, df_test, chromosome):

    """This function takes a chromosome and a df train and df test of HMMcopy and returns the rsquare of the fit"""

    # get the dfs of the chromosome
    df_train = df_train[df_train.chr==chromosome]
    df_test = df_test[df_test.chr==chromosome]

    # get the function that interpolates position vs copy number
    interpolation_function = scipy_interpolate.interp1d(df_train.middle_position, df_train.relative_CN, bounds_error=True, kind="linear", assume_sorted=True)

    # get the prediction on the test
    test_relative_CN = interpolation_function(df_test.middle_position)

    # calculate the rsquare between the test CN and the actual coverage
    rsquare = r2_score(df_test.reads.values, test_relative_CN)

    # debug
    if pd.isna(rsquare): raise ValueError("rsquare can't be NaN")

    # get miniumum 0
    return max([0.0, rsquare])


def get_crossValidated_rsquare_HMMcopy_givenParameters(Ip, parameters_dict, outfile_HMMcopy_prefix, training_coverage_files, testing_coverage_dfs):

    """This function takes a parameters dict and returns the R2 of how these parameters work on 5-fold cross validation. coverage_files are a list of matched training and testing files for a balanced cross validation selection. This function will report how good each of the training files works on the testing files. """

    print_if_verbose("checking parameters %i"%Ip)

    # define the rsquares
    rsquares_list = []

    # go through each cvID
    for cvID, training_coverage_file in enumerate(training_coverage_files):

        # define the outfile of this CV 
        outfile_HMMcopy_cvID = "%s.train_cv%i.tab"%(outfile_HMMcopy_prefix, cvID)

        # run hmm copy on the training df
        df_train = run_HMMcopy(training_coverage_file, outfile_HMMcopy_cvID, parms_dict=parameters_dict, replace=True)
        #remove_file(outfile_HMMcopy_cvID)

        # get the test df
        df_test = testing_coverage_dfs[cvID]

        # get the rsquares for each chromosome
        rsquare = np.mean(list(map(lambda c: get_rsquare_one_chromosome_crossValidated_rsquare_HMMcopy(df_train, df_test, c), set(df_train.chr))))

        if rsquare<=0: break

        # keep 
        rsquares_list.append(rsquare)

    # get the mean and std rsquares
    if len(rsquares_list)!=len(training_coverage_files):
        mean_rsquare = 0.0
        std_rsquare = 1.0

    else:
        mean_rsquare = max([0.0, np.mean(rsquares_list)])
        std_rsquare = np.std(rsquares_list)

    return mean_rsquare, std_rsquare


def run_CNV_calling_HMMcopy(outdir, replace, threads, df_coverage, ploidy, reference_genome, mitochondrial_chromosome, min_number_of_regions_CNVcalling=50):

    """This function runs HMMcopy and returns a df with the relative_CN."""

    print_if_verbose("running HMMcopy")
    make_folder(outdir)

    # define chroms
    all_chromosomes = set(get_chr_to_len(reference_genome))
    if mitochondrial_chromosome!="no_mitochondria": mtDNA_chromosomes = set(mitochondrial_chromosome.split(","))
    else: mtDNA_chromosomes = set()
    gDNA_chromosomes = all_chromosomes.difference(mtDNA_chromosomes)

    # initialize the final df
    final_fields = ["chromosome", "start", "end", "corrected_relative_coverage", "relative_CN"]
    df_CNperWindow_HMMcopy = pd.DataFrame(columns=final_fields)

    # iterate through each genome
    for type_genome, chroms in [("gDNA", gDNA_chromosomes), ("mtDNA", mtDNA_chromosomes)]:

        # define the outdir of the genome
        outdir_genome = "%s/%s"%(outdir, type_genome); make_folder(outdir_genome)

        # debug
        if len(chroms)==0: continue

        # get the df coverage for this genome
        df_coverage_genome = df_coverage[df_coverage.chromosome.isin(chroms)]

        # debug the fact that there are not enough regions. If so set everything to relative_CN==1
        if len(df_coverage_genome)<min_number_of_regions_CNVcalling: 

            # get relative_CN as 1
            df_CNperWindow_HMMcopy_genome = df_coverage_genome[["chromosome", "start", "end", "corrected_relative_coverage"]]
            df_CNperWindow_HMMcopy_genome["relative_CN"] = 1.0

        else:

            ######### PREPARE THE INPUT DATA #########

            # define a file that has the coverage per region 
            df_coverage_HMMcopy = df_coverage_genome.rename(columns={"GCcontent":"gc", "median_mappability":"map"})
            df_coverage_HMMcopy["start"] += 1 # it has to be 1-based

            # add fields necessary to run HMMcopy
            df_coverage_HMMcopy["reads"] = df_coverage_HMMcopy.corrected_relative_coverage
            df_coverage_HMMcopy["chr"] = df_coverage_HMMcopy.chromosome
            df_coverage_HMMcopy["cor.gc"] = df_coverage_HMMcopy.reads
            df_coverage_HMMcopy["cor.map"] = df_coverage_HMMcopy.reads

            def get_copy_or_NaN(x):
                if x==0.0: return np.nan
                else: return np.log2(x)
            df_coverage_HMMcopy["copy"] = df_coverage_HMMcopy["cor.map"].apply(get_copy_or_NaN)

            # write as tab
            coverage_file = "%s/coverage_forHMMcopy.tab"%outdir_genome
            HMMcopy_fields = ["chr", "start", "end", "reads", "gc", "map", "cor.gc", "cor.map", "copy"]
            save_df_as_tab(df_coverage_HMMcopy[HMMcopy_fields], coverage_file)

            ##########################################

            ##### DEFINE THE RANGE OF PARAMETERS FOR CROSS-VALIDATION #####

            # these are the parameters:

            """ 
            - e : Probability of extending a segment, increase to lengthen segments, decrase to shorten segments. Range: (0, 1)
            - strength: Strength of initial e suggestion, reducing allows e to change, increasing makes e undefiable. Range: [0, Inf)
            - mu: Suggested median for copy numbers in state, change to readjust classification of states. Range: (-Inf, Inf)
            - lambda: Suggested precision (inversed variance) for copy numbers in state, increase to reduce overlap between states. Range: [0, Inf)
            - nu: Suggested degree of freedom between states, increase to reduce overlap between states. Range: [0, Inf)
            - kappa: Suggested distribution of states. Should sum to 1.
            - m: Optimal value for mu, difference from corresponding mu value determines elasticity of the mu value. i.e. Set to identical value as mu if you don’t want mu to move much.
            - eta Mobility of mu, increase to allow more movement. Range: [0, Inf)
            - gamma Prior shape on lambda, gamma distribution. Effects flexibility of lambda.
            - S Prior scale on lambda, gamma distribution. Effects flexibility of lambda
            """

            print_if_verbose("generating all possible parameters")

            # define the ideal coverage of each state
            ideal_coverage = [0.01, 0.5, 1, 1.5, 2, 2.5]

            # define the default mu
            default_mu = [-0.458558247, -0.215877601, -0.002665686, 0.191051578, 0.347816046, 1.664333241]

            # create a df with training parameters
            parameters_dict = {}; Ip=0
            for e in np.linspace(0.9, 0.9999999, 3):
                for strength in ["1e3", "1e7", "1e30"]:
                    for mu_multiplier in [None, 0.1]: 
                        for lambda_val in [20]:
                            for nu in [2.1]:
                                for kappa in [[50, 50, 700, 100, 50, 50]]:
                                    for m_multiplier in [None, 0.1]: 
                                        for eta in [50000]:
                                            for gamma in [3]:
                                                for S in [0.02930164]:

                                                    # define the mu
                                                    if mu_multiplier is None: mu = default_mu
                                                    else: mu = [np.log2(max([0.01, c+mu_multiplier])) for c in ideal_coverage]

                                                    # define the m
                                                    if m_multiplier is None: m = default_mu
                                                    else: m = [np.log2(max([0.01, c+m_multiplier])) for c in ideal_coverage]

                                                    # keep the data
                                                    parameters_dict[Ip] = dict(e=e, mu=mu, lambda_val=lambda_val, nu=nu, kappa=kappa, m=m, eta=eta, gamma=gamma, S=S, strength=strength)
                                                    Ip+=1

            # get as df
            parameter_fields = ["e", "mu", "lambda_val", "nu", "kappa", "m", "eta", "gamma", "S", "strength"]
            parameters_df = pd.DataFrame(parameters_dict).transpose().sort_index()

            ################################################################

            ####### FIND THE BEST PARAMETERS #######

            print_if_verbose("finding best parameters out of %i"%len(parameters_df))

            # define the cross validation datasets
            typeDF_to_cvID_to_df = {t : {cvID : pd.DataFrame() for cvID in range(5)} for t in {"train", "test"}}

            for chrom in chroms:

                # get the df with the coverage
                df_coverage_HMMcopy_chrom = df_coverage_HMMcopy[df_coverage_HMMcopy.chr==chrom][HMMcopy_fields]

                # add the midle position
                df_coverage_HMMcopy_chrom["middle_position"] = df_coverage_HMMcopy_chrom.start + (df_coverage_HMMcopy_chrom.end-df_coverage_HMMcopy_chrom.start)/2

                # debug too short chromosomes
                if len(df_coverage_HMMcopy_chrom)<10: continue

                # go through CV of each chromosome
                kfold_object = KFold(n_splits=5, random_state=1, shuffle=True)
                for cvID, (numeric_train_index, numeric_test_index) in enumerate(kfold_object.split(df_coverage_HMMcopy_chrom.index)):

                    # get the dfs
                    train_df = df_coverage_HMMcopy_chrom.iloc[numeric_train_index]
                    test_df = df_coverage_HMMcopy_chrom.iloc[numeric_test_index]

                    # get the test df that is within the bounds of train_df
                    test_df = test_df[(test_df.middle_position>=min(train_df.middle_position)) & (test_df.middle_position<=max(train_df.middle_position))]

                    # add to the dfs
                    typeDF_to_cvID_to_df["train"][cvID] = typeDF_to_cvID_to_df["train"][cvID].append(train_df)
                    typeDF_to_cvID_to_df["test"][cvID] = typeDF_to_cvID_to_df["test"][cvID].append(test_df)


            # write the training datasets into outdir outdir_genome
            training_coverage_files = []
            for cvID in range(5): 

                # get file
                df = typeDF_to_cvID_to_df["train"][cvID]
                training_coverage_file = "%s/training_df_cv%i.tab"%(outdir_genome, cvID)
                save_df_as_tab(df, training_coverage_file)

                # keep
                training_coverage_files.append(training_coverage_file)

            # get a list with the testing dfs
            testing_coverage_dfs = [typeDF_to_cvID_to_df["test"][cvID] for cvID in range(5)]

            # get a list with the R2 of each parameter set from cross validation
            list_inputs = [(Ip, parameters_dict[Ip], "%s/HMMcopy_parms%i.tab"%(outdir_genome, Ip), training_coverage_files, testing_coverage_dfs) for Ip in parameters_df.index]

            with multiproc.Pool(threads) as pool:
                list_rsquares_and_stds = pool.starmap(get_crossValidated_rsquare_HMMcopy_givenParameters, list_inputs) 

                pool.close()
                pool.terminate()

            # add the lists in the df
            parameters_df["cv_rsquare_mean"] = [x[0] for x in list_rsquares_and_stds]
            parameters_df["cv_rsquare_std"] = [x[1] for x in list_rsquares_and_stds]
            parameters_df["cv_rsquare_inverse_std"] = 1/parameters_df.cv_rsquare_std

            # get the best parameters
            best_parameters_row = parameters_df.sort_values(by=["cv_rsquare_mean", "cv_rsquare_inverse_std"], ascending=False).iloc[0]

            # if there are some best parameters, pick them
            if best_parameters_row["cv_rsquare_mean"]>0: 

                print_if_verbose("The best parameters have an R2=%.3f +- %.3f (SD)"%(best_parameters_row["cv_rsquare_mean"], best_parameters_row["cv_rsquare_std"]))

                best_parameters_dict = dict(best_parameters_row[parameter_fields])

            else:

                print_if_verbose("There are no optimum parms, set the default ones")
                best_parameters_dict = default_parms_dict_HMMcopy

            ###########################################

            ####### RUN WITH BEST PARAMETERS AND ARRANGE df_CNperWindow_HMMcopy_genome #########

            print_if_verbose("running with best parameters")

            # run HMMcopy with the best parameters
            best_parameters_outfile = "%s/HMMcopy_output_best_parms.tab"%outdir_genome
            df_HMMcopy_best = run_HMMcopy(coverage_file, best_parameters_outfile, replace=replace, parms_dict=best_parameters_dict)

            # change the start to be 0-based
            df_HMMcopy_best["start"] = df_HMMcopy_best.start - 1

            # define the final output
            df_CNperWindow_HMMcopy_genome = df_HMMcopy_best
            df_CNperWindow_HMMcopy_genome["corrected_relative_coverage"] = df_CNperWindow_HMMcopy_genome.reads
            df_CNperWindow_HMMcopy_genome["chromosome"] = df_CNperWindow_HMMcopy_genome.chr

            #####################################################################################

        # keep
        df_CNperWindow_HMMcopy = df_CNperWindow_HMMcopy.append(df_CNperWindow_HMMcopy_genome[final_fields])

        # remove
        delete_folder(outdir_genome)


    return df_CNperWindow_HMMcopy

#########################
#########################
#########################

############################
####### ANEUFINDER #########
############################


def run_AneuFinder(coverage_file, outfile, threads, parms_dict={"R":10, "sig_lvl":0.1}, replace=False):

    """This function runs aneufinder from coverage_file to outfile using the parms dict. It returns a df with coverage_outfile and relative_CN (which is the CN divided by 2) """

    # define the parms
    R = int(parms_dict["R"])
    sig_lvl = float(parms_dict["sig_lvl"])

    # run aneufinder
    if file_is_empty(outfile) or replace is True:

        # get the outfile
        outfile_tmp = "%s.tmp"%outfile
        print_if_verbose("generating %s"%outfile)

        # run
        aneufinder_std = "%s.generating.std"%outfile
        run_cmd("%s --coverage_table %s --outfile %s --threads %i --R %i --sig_lvl %.4f > %s 2>&1"%(run_AneuFinder_R, coverage_file, outfile_tmp, threads, R, sig_lvl, aneufinder_std), env=EnvName_AneuFinder)

        # save
        os.rename(outfile_tmp, outfile)

    # get the out df
    df_out = get_tab_as_df_or_empty_df(outfile).set_index("seqnames")
    df_out["relative_CN"] = df_out["copy.number"] / 2

    # add to df the relative_CN. Everything is called as a diploid
    df = get_tab_as_df_or_empty_df(coverage_file)
    def get_relative_CN(r):

        # get the chromosome df
        df_c = df_out.loc[{r["seqnames"]}]

        # get the window of df_out that includes the CN of r
        all_relative_CNs = list(df_c[(df_c.start<=r["start"]) & (df_c.end>=r["end"])].relative_CN)

        # check that it is just 1
        if len(all_relative_CNs)!=1: raise ValueError("there should only be one relative CN.")

        return all_relative_CNs[0]

    df["relative_CN"] = df.apply(get_relative_CN, axis=1)
    if any(pd.isna(df.relative_CN)): raise ValueError("there can't be NaNs in df.relative_CN")

    return df


def run_CNV_calling_AneuFinder(outdir, replace, threads, df_coverage, ploidy, reference_genome, mitochondrial_chromosome, read_length, min_number_of_regions_CNVcalling=50):

    """This function runs AneuFinder and returns a df with the relative CN"""

    print_if_verbose("running AneuFinder")
    make_folder(outdir)

    # define chroms
    all_chromosomes = set(get_chr_to_len(reference_genome))
    if mitochondrial_chromosome!="no_mitochondria": mtDNA_chromosomes = set(mitochondrial_chromosome.split(","))
    else: mtDNA_chromosomes = set()
    gDNA_chromosomes = all_chromosomes.difference(mtDNA_chromosomes)

    # initialize the final df
    final_fields = ["chromosome", "start", "end", "corrected_relative_coverage", "relative_CN"]
    df_CNperWindow_AneuFinder = pd.DataFrame(columns=final_fields)

    # iterate through each genome
    for type_genome, chroms in [("mtDNA", mtDNA_chromosomes), ("gDNA", gDNA_chromosomes)]:

        # define the outdir of the genome
        outdir_genome = "%s/%s"%(outdir, type_genome); make_folder(outdir_genome)

        # debug
        if len(chroms)==0: continue

        # get the df coverage for this genome
        df_coverage_genome = df_coverage[df_coverage.chromosome.isin(chroms)]

        # debug the fact that there are not enough regions. If so set everything to relative_CN==1
        if len(df_coverage_genome)<min_number_of_regions_CNVcalling: 

            # get relative_CN as 1
            df_CNperWindow_AneuFinder_genome = df_coverage_genome[["chromosome", "start", "end", "corrected_relative_coverage"]]
            df_CNperWindow_AneuFinder_genome["relative_CN"] = 1.0

        else:

            print_if_verbose("running aneufinder")

            ######## DEFINE THE INPUT #########

            # get a df_coverage with 1-based ["seqnames", "ranges", "strand", "counts"], where counts is the corrected_relative_coverage

            # add fields necesary for AneuFinder
            df_coverage_genome["seqnames"] = df_coverage_genome.chromosome
            df_coverage_genome["strand"] = "*" 
            df_coverage_genome["start"] = df_coverage_genome.start+1

            # add the counts for aneufinder, proportional to the median coverage * the corrected relative coverage * length of each window / read length
            median_coverage_positon = get_median_coverage(df_coverage_genome, mitochondrial_chromosome) 
            if median_coverage_positon<=0: raise ValueError("The median coverage per position %i is not valid"%median_coverage_positon)
            df_coverage_genome["counts"] = (df_coverage_genome.corrected_relative_coverage*df_coverage_genome.width*(median_coverage_positon/read_length)).apply(int)

            # define the fields
            AneuFinder_fields = ["seqnames", "start", "end", "strand", "counts"]

            ##################################

            ########## RUN ANEUFINDER  WITH BEST PARMS ########

            best_parameters_dict = {"R":10, "sig_lvl":0.1}
            print_if_verbose("running AneuFinder with best parameters: %s"%best_parameters_dict)

            # write to file the coverage file
            coverage_file = "%s/coverage_file.tab"%outdir_genome
            save_df_as_tab(df_coverage_genome[AneuFinder_fields], coverage_file)

            # run  with the best parameters
            best_parameters_outfile = "%s/AneuFinder_output_best_parms.tab"%outdir_genome
            df_AneuFinder_best = run_AneuFinder(coverage_file, best_parameters_outfile, threads, parms_dict=best_parameters_dict, replace=replace)

            initial_len_df_AneuFinder_best = len(df_AneuFinder_best)

            # add the corrected_relative_coverage
            df_AneuFinder_best = df_AneuFinder_best.merge(df_coverage_genome[["chromosome", "start", "end", "corrected_relative_coverage"]], left_on=["seqnames", "start", "end"], right_on=["chromosome", "start", "end"], how="left", validate="one_to_one")

            # debug merge
            if initial_len_df_AneuFinder_best!=len(df_AneuFinder_best): raise ValueError("There lengths are not the same")

            # change the start to be 0-based
            df_AneuFinder_best["start"] = df_AneuFinder_best.start - 1

            # define the final output
            df_CNperWindow_AneuFinder_genome = df_AneuFinder_best

            ###################################################

        # keep
        df_CNperWindow_AneuFinder = df_CNperWindow_AneuFinder.append(df_CNperWindow_AneuFinder_genome[final_fields])

        # delete folder
        #delete_folder(outdir_genome)

    return df_CNperWindow_AneuFinder


############################
############################
############################

def run_CNV_calling_CONY_one_chromosome(reference_genome, outdir, chromosome, replace, window_size, threads, chrom_len, sample_name, df_coverage, ploidy, min_number_of_regions_CNVcalling=50):

    """ runs CONY on a given chromosome into outdir. It returns a df with the  """

    print_if_verbose("running CONY for %s and sample %s"%(chromosome, sample_name))
    make_folder(outdir)

    # define the curdir
    CurDir = get_fullpath(".")

    # define the outdir as the full path
    outdir = get_fullpath(outdir)

    # change dir to outdir
    os.chdir(outdir)

    ############### PREPARE CONY DATA ###############

    # create a soft link on the CONY library, which is necessary because there are functions in CONY.R which use this
    soft_link_files(libraries_CONY, "%s/CONY.R"%outdir)

    # define a file that has the coverage per chromosome so that CONY can start in the UsedRD part
    df_coverage_chrom = df_coverage[df_coverage.chromosome==chromosome].rename(columns={"chromosome":"seqname", "mediancov_1":"rawRD", "corrected_relative_coverage":"AdjRD", "GCcontent":"GC"})
    df_coverage_chrom["start"] += 1 # it has to be 1-based
    df_coverage_chrom["nonAmb"] = 1 - df_coverage_chrom.fraction_N_bases

    # if there are less than min_number_of_regions_CNVcalling return an empty df, meaning that there are no calls
    if len(df_coverage_chrom)<min_number_of_regions_CNVcalling: return pd.DataFrame()

    # write the file
    coverage_file_chrom = "./CONY.2-TempRegion.%s.%s.SumUp.AdjRD.txt"%(chromosome, sample_name)
    CONY_fields = ["seqname", "start", "end", "width", "nonAmb", "GC", "AdjRD"]
    save_df_as_tab(df_coverage_chrom[CONY_fields], coverage_file_chrom)

    # define the fragment length: the length of analytic fragments for each lane. The suggested default is 500,000 (bps) at a time. RJMCMC would be run with several lanes simultaneously via snow package. The number of lanes (can be thought as threads) is total number of analytic windows/ (fragment length/ window size). If we increase the fragment length the number of lanes (threads) decreases.

    
    fraction_fragment_len = 0.9 # this was the original. It is a very low number, so it implies paralelism
    #fraction_fragment_len = 3 # this is to increase the number of lanes
    fragment_len = int((window_size*(len(df_coverage_chrom)/threads))*fraction_fragment_len) + 1

    # set to a minimum, for avoiding CONY errors
    min_fragment_len = int(min_number_of_regions_CNVcalling*window_size)
    fragment_len = max([fragment_len, min_fragment_len])

    #################################################

    ########## RUN CONY ##########

    # define the final file
    final_file = "%s/CONY_finished.txt"%outdir
    results_file = "%s/CONY.Result.%s.%s.SumUp.Single.Window.txt"%(outdir, chromosome, sample_name)

    if file_is_empty(final_file) or replace is True:

        # run CONY
        cony_std = "%s/running_cony.std"%outdir
        print_if_verbose("Running CONY. The std is in %s. Running on a fragment length of %i"%(cony_std, fragment_len))

        cmd = "%s --libraries_CONY %s --chromosome %s --sample_name %s --outdir %s --ploidy %i --fragment_len %i > %s 2>&1"%(run_CONY_R, libraries_CONY, chromosome, sample_name, outdir, ploidy, fragment_len, cony_std)

        run_cmd(cmd, env=EnvName_CONY)

        # debug
        if file_is_empty(results_file): raise ValueError("there were no results created for chromosome %s"%chromosome)
        
        # remove std
        remove_file(cony_std)

        # make the final file
        open(final_file, "w").write("CONY finished")

    # remove unnecesary files
    for f in os.listdir(outdir):
        if not f.startswith("CONY.Result.") and f!="CONY_finished.txt": remove_file(f)

    # at the end return to the CurDir
    os.chdir(CurDir)

    # load dfs
    df_perWindow = pd.read_csv(results_file, sep=" ")

    # define the chromsomal copy number (indicated by the coverage of the regions with CN=2 (no changes))
    CN1_relative_coverage = float(np.median(df_perWindow[df_perWindow.CN==2].AdjRD))
    if ploidy==1: possible_CN1s = [1, 2]
    elif ploidy==2: possible_CN1s = [0.5, 1, 1.5, 2]
    else: raise ValueError("CNV calling is not possible for ploidy %i"%ploidy)
    chromosomal_CN = find_nearest(possible_CN1s, CN1_relative_coverage)

    # add the relative CN
    df_perWindow["relative_CN"] = (df_perWindow.CN / 2)*chromosomal_CN
    df_perWindow["chromosomal_relative_CN"] = chromosomal_CN


    #################################

    ###### PLOT THE CN ACROSS THE CHROMOSOME ######
    print_if_verbose("plotting CNV...")

    fig = plt.figure(figsize=(20,5))
    nrows = 1

    # plot each type of coverage
    plt.subplot(nrows, 1, 1)
    plt.plot(df_perWindow.start, df_perWindow.AdjRD, color="blue", label="Adjusted Read Depth")
    plt.plot(df_perWindow.start, df_perWindow.relative_CN, color="red", label="Copy Number CONY", linestyle="--")
    plt.axhline(1, color="gray", linestyle="--", linewidth=.9)
    plt.xlabel("start window")
    plt.ylabel("read depth / Copy Number")
    plt.title("CNV for %s"%chromosome)
    plt.legend(bbox_to_anchor=(1, 1))

    # save
    fig.savefig("%s/coverage_%s.pdf"%(outdir, chromosome))
    plt.close(fig)

    ################################################

    ######## PLOT INTERACTIVE ########

    # init fig
    fig = tools.make_subplots(rows=1, cols=1, specs=[[{}]], vertical_spacing=0.0, horizontal_spacing=0.0, subplot_titles=(""), shared_yaxes=True, shared_xaxes=True, print_grid=True)

    # get the relative coverage
    fig.append_trace(go.Scatter(x=list(df_perWindow.start), y=list(df_perWindow.AdjRD), showlegend=True, mode="lines+markers", marker=dict(symbol="circle", color="blue", size=4), opacity=1, hoveron="points+fills", name="Adjusted read depth", line=dict(color="blue", width=2, dash=None)) , 1, 1) 

    # get the CN
    fig.append_trace(go.Scatter(x=list(df_perWindow.start), y=list(df_perWindow.relative_CN), showlegend=True, mode="lines", line=dict(color="red", width=2, dash="dash"), opacity=1, hoveron="points+fills", name="CONY CN"), 1, 1) 

    # get figure
    fig['layout'].update(title="relative coverage and CN %s"%chromosome, margin=go.Margin(l=250, r=250, b=50, t=50), showlegend=True)
    config={'editable': False}
    off_py.plot(fig, filename="%s/coverage_%s_interactive.html"%(outdir, chromosome), auto_open=False, config=config)

    ##################################

    # change the start to be back to 0
    df_perWindow["start"] =  df_perWindow["start"] - 1

    # return the df_CNV
    return df_perWindow

def get_sample_name_from_bam(bam):

    """Gets the sample name of a bam file"""

    sample_name = str(subprocess.check_output("%s view -H %s | egrep '^@RG'"%(samtools, bam), shell=True)).split("ID:")[1].split("\\t")[0]
    
    return sample_name

def get_coverage_df_for_windows_of_genome(sorted_bam, reference_genome, outdir, replace, threads, window_size):

    """Returns a coverage df for windows of window_size, saving files under outdir"""

    final_coverage_file = "%s/coverage_per_windows_%ibp.tab"%(outdir, window_size)
    if file_is_empty(final_coverage_file) or replace is True:

        # first generate the windows file
        windows_file = "%s.windows%ibp.bed"%(reference_genome, window_size)
        if file_is_empty(windows_file) or replace is True:

            windows_file_stderr = "%s.generating.stderr"%windows_file
            windows_file_tmp = "%s.tmp"%windows_file

            print_if_verbose("genearting windows_file. The stderr is in %s"%windows_file_stderr)
            run_cmd("%s makewindows -g %s.fai -w %i > %s 2>%s"%(bedtools, reference_genome, window_size, windows_file_tmp, windows_file_stderr))

            remove_file(windows_file_stderr)
            os.rename(windows_file_tmp, windows_file)

        # get coverage calculated
        coverage_df = get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, windows_file, replace=replace, run_in_parallel=True, delete_bams=True, threads=threads)

        # save
        save_df_as_tab(coverage_df, final_coverage_file)

    # load
    coverage_df = get_tab_as_df_or_empty_df(final_coverage_file)

    return coverage_df

def generate_genome_mappability_file(genome, replace=False, threads=4):

    """Takes a genome in fasta and generates a file indicating the mappability of each position of the genome.
    This writes a file (output) that contains, for each position (i) in the genome, the value of 1/ (number of 30-mers equal to the one that starts in i, allowing 2 missmatches). If it is 1 it is a uniquely mapped position """

    # first create the index
    genome_dir = "/".join(genome.split("/")[0:-1])
    genome_name = genome.split("/")[-1]
    idx_folder = "%s/%s_genmap_idx"%(genome_dir, genome_name.split(".")[0])

    # define expected files
    expected_idx_files = ["index.info.concat", "index.lf.drv.sbl", "index.sa.val", "index.txt.limits", "index.lf.drv", "index.info.limits", "index.rev.lf.drp"]
    if any([file_is_empty("%s/%s"%(idx_folder, x)) for x in expected_idx_files]) or replace is True:

        # remove previously generated indices
        if os.path.isdir(idx_folder): shutil.rmtree(idx_folder)

        print("Generating index")
        run_cmd("%s index -F %s -I %s -v"%(genmap, genome, idx_folder))

    # generate map
    map_folder = "%s/%s_genmap_map"%(genome_dir, genome_name.split(".")[0])
    map_outfile = "%s/%s.genmap.bed"%(map_folder, genome_name.split(".")[0])
    if file_is_empty(map_outfile) or replace is True:

        if not os.path.isdir(map_folder): os.mkdir(map_folder)

        print("Generating map")
        run_cmd("%s map -E 2 -K 30 -I %s -O %s -b --threads %i -v "%(genmap, idx_folder, map_folder, threads))

    # deine the long file
    map_outfile_long = "%s.long.bed"%map_outfile

    if file_is_empty(map_outfile_long) or replace is True:

        print("getting df in the long format")

        # convert to a file where each position in the genome appears. This is important because genmap generates a file that has only ranges
        df_map = pd.read_csv(map_outfile, sep="\t", header=None, names=["chromosome", "start", "end", "strand", "map_idx"])
        df_map["chromosome_real"] = df_map.chromosome.apply(lambda x: x.split()[0])

        # define a list with the positions and the scores for that window
        df_map["positions_list"] = df_map[["start", "end"]].apply(lambda r: list(range(r["start"], r["end"])), axis=1)
        df_map["length_range"] = df_map.positions_list.apply(len)
        df_map["map_idx_list"] = df_map[["length_range", "map_idx"]].apply(lambda r: [r["map_idx"]]*int(r["length_range"]), axis=1)
        df_map["chromosome_list"] = df_map[["length_range", "chromosome_real"]].apply(lambda r: [r["chromosome_real"]]*int(r["length_range"]), axis=1)

        # get chr_to_len
        chr_to_len = get_chr_to_len(genome)

        # initialize a dictionary that will store chromosome, position and mappability_score as lists
        expanded_data_dict = {"chromosome":[], "position":[], "unique_map_score":[]}

        # go through each row of the dataframe and append the lists
        for chromosome_list, positions_list, map_idx_list in df_map[["chromosome_list", "positions_list", "map_idx_list"]].values:

            expanded_data_dict["chromosome"] += chromosome_list
            expanded_data_dict["position"] += positions_list
            expanded_data_dict["unique_map_score"] += map_idx_list

        df_long = pd.DataFrame(expanded_data_dict)

        # add the missing positions with the last windows score
        for chrom, length_chrom in get_chr_to_len(genome).items():

            # get the df_chrom
            df_chrom = df_long[df_long.chromosome==chrom]

            # define the missing positions
            all_positions = set(range(0, length_chrom))
            missing_positions = sorted(all_positions.difference(set(df_chrom.position)))
            n_missing_positions = len(missing_positions)

            # add to df_long
            if n_missing_positions>0:

                # add them with 0 mappability
                df_missing = pd.DataFrame({"chromosome":[chrom]*n_missing_positions, "position":missing_positions, "unique_map_score":[0.0]*n_missing_positions})

                df_long = df_long.append(df_missing)
        
        # sort by chromosome and position
        df_long = df_long.sort_values(by=["chromosome", "position"])

        # add whether it is uniquely mappable
        df_long["is_uniquely_mappable"] = (df_long.unique_map_score>=1.0).apply(int)

        # save
        save_df_as_tab(df_long, map_outfile_long)

    return map_outfile_long

def get_df_with_GCcontent(df_windows, genome, gcontent_outfile, replace=False):

    """This function takes a df with windows of the genome and adds the gc content for each window, writing a file under gcontent_outfile. It will only do those that have been already measured"""

    print_if_verbose("Getting GC content")

    if file_is_empty(gcontent_outfile) or replace is True:

        # define the initial index
        initial_index = list(df_windows.index)
        initial_cols = list(df_windows.columns)

        # resort
        df_windows = df_windows.sort_values(by=["chromosome", "start", "end"])

        print_if_verbose("getting GC content for %i new windows"%len(df_windows))

        # get the GC content file for each position
        gc_content_outfile_perPosition = generate_nt_content_file(genome, replace=replace, target_nts="GC")
        gc_df = pd.read_csv(gc_content_outfile_perPosition, sep="\t")[["chromosome", "position", "is_in_GC"]].sort_values(by=["chromosome", "position"])

        # add the ID
        gc_df["ID"] =  list(range(0, len(gc_df)))

        # add the end
        gc_df["position_plus1"] = gc_df.position + 1

        # generate a bed with the gc positions
        gc_positions_bed = "%s.gc_positions.bed"%gcontent_outfile
        gc_df[gc_df.is_in_GC==1].sort_values(by=["chromosome", "position"])[["chromosome", "position", "position_plus1", "ID"]].to_csv(gc_positions_bed, sep="\t", header=None, index=False)

        # generate the bed with the windows
        target_windows_bed = "%s.target_windows.bed"%gcontent_outfile
        df_windows["IDwindow"] = list(range(0, len(df_windows)))
        df_windows[["chromosome", "start", "end", "IDwindow"]].to_csv(target_windows_bed, sep="\t", header=None, index=False)

        # run bedmap to get a file where each line corresponds to the regions to which each target_windows_bed
        bedmap_outfile = "%s.bedmap.target_windows_overlappingGC.txt"%gcontent_outfile
        bedmap_stderr = "%s.stderr"%bedmap_outfile

        print_if_verbose("running bedmap. The stderr is in %s"%bedmap_stderr)
        run_cmd("%s --fraction-map 1.0 --echo-map-id  --delim '\t' %s %s > %s 2>%s"%(bedmap, target_windows_bed, gc_positions_bed, bedmap_outfile, bedmap_stderr))

        # load bedmap df into df
        df_windows["overlapping_GC_IDs"] = open(bedmap_outfile, "r").readlines()

        # add the n_GC_positions
        def get_NaN_to_empty_str(x):
            if pd.isna(x): return ""
            else: return x

        all_possibleIDs = set(gc_df.ID.apply(str))
        def get_set_from_string(x): 

            # define set
            set_x = set(x.strip().split(";")).difference({""})

            # debug 
            if len(set_x.difference(all_possibleIDs))>0: raise ValueError("there are unexpected IDs in %s"%set_x)

            return set_x

        df_windows["n_GC_positions"] = df_windows.overlapping_GC_IDs.apply(get_NaN_to_empty_str).apply(get_set_from_string).apply(len)

        # add the GC content
        df_windows["length"] = df_windows.end - df_windows.start
        df_windows["GCcontent"] = df_windows.n_GC_positions / df_windows.length

        # debug
        if any(pd.isna(df_windows.GCcontent)) or any(pd.isna(df_windows.n_GC_positions)) or any(df_windows.GCcontent>1): raise ValueError("Something went went wrong with the GC content calcilation")

        for f in [gc_positions_bed, target_windows_bed, bedmap_outfile, bedmap_stderr]: remove_file(f)

        # at the end save the df windows
        df_windows.index = initial_index
        df_windows = df_windows[initial_cols + ["GCcontent"]]
        save_object(df_windows, gcontent_outfile)

    # load
    df_windows = load_object(gcontent_outfile)

    return df_windows

def get_df_windows_with_median_mappability(df_windows, reference_genome, mappability_outfile, replace, threads):

    """Takes a df windws and returns it with a field called median_mappability, which is the median mappability for each position of the region. Mappability for a position i is 1 / (how often a 30-mer starting at i occurs in the genome), with up to 2 mismatches."""
    
    print_if_verbose("Getting mappability")

    if file_is_empty(mappability_outfile) or replace is True:

        # define the initial index
        initial_index = list(df_windows.index)
        initial_cols = list(df_windows.columns)

        # resort
        df_windows = df_windows.sort_values(by=["chromosome", "start", "end"])

        print_if_verbose("getting mappability for %i new windows"%len(df_windows))

        # get the mappability per position into map_df
        mappability_per_position_file = generate_genome_mappability_file(reference_genome, replace=replace, threads=threads)
        map_df = pd.read_csv(mappability_per_position_file, sep="\t")

        # add the ID
        map_df["ID"] =  list(range(0, len(map_df)))

        # add the end for bedmap
        map_df["position_plus1"] = map_df.position + 1

        # generate a bed with the map positions and the unique_map_score
        mappability_bed = "%s.mappability.bed"%mappability_outfile
        map_df.sort_values(by=["chromosome", "position"])[["chromosome", "position", "position_plus1", "ID", "unique_map_score"]].to_csv(mappability_bed, sep="\t", header=None, index=False)

        # generate the bed with the windows
        target_windows_bed = "%s.target_windows.bed"%mappability_outfile
        df_windows["IDwindow"] = list(range(0, len(df_windows)))
        df_windows[["chromosome", "start", "end", "IDwindow"]].to_csv(target_windows_bed, sep="\t", header=None, index=False)

        # run bedmap to get a file where each line corresponds to the regions to which each target_windows_bed
        bedmap_outfile = "%s.bedmap.target_windows_overlappingMappability.txt"%mappability_outfile
        bedmap_stderr = "%s.stderr"%bedmap_outfile

        print_if_verbose("running bedmap. The stderr is in %s"%bedmap_stderr)
        run_cmd("%s --fraction-map 1.0 --median --delim '\t' %s %s > %s 2>%s"%(bedmap, target_windows_bed, mappability_bed, bedmap_outfile, bedmap_stderr))

        # load bedmap df into df
        df_windows["median_mappability"] = open(bedmap_outfile, "r").readlines()
        df_windows["median_mappability"] = df_windows.median_mappability.apply(float)

        # debug
        if any(pd.isna(df_windows.median_mappability)) or any(df_windows.median_mappability>1): 

            print("These are the NaN regions")
            print(df_windows[pd.isna(df_windows.median_mappability)][["chromosome", "start", "end"]])

            raise ValueError("Something went went wrong with the mappability calculation")

        # remove files
        for f in [mappability_bed, target_windows_bed, bedmap_outfile, bedmap_stderr]: remove_file(f)

        # at the end save the df windows
        df_windows.index = initial_index
        df_windows = df_windows[initial_cols + ["median_mappability"]]
        save_df_as_tab(df_windows, mappability_outfile)

    # load
    df_windows = get_tab_as_df_or_empty_df(mappability_outfile)

    return df_windows

def make_plots_coverage_parameters(df_cov, plots_dir):

    """This function makes some plots useful."""

    print_if_verbose("making plots coverage")

    #for xfield, yfield in [("start", "coverage"), ("GCcontent", "coverage"), ("median_mappability", "coverage"), ("start", "median_mappability"), ("raw_distance_to_telomere", "coverage"), ("raw_distance_to_telomere", "median_mappability")]:
    for xfield, yfield in [("raw_distance_to_telomere", "relative_coverage"), ("raw_distance_to_telomere", "median_mappability")]:

        # get fig
        fig = plt.figure(figsize=(7, 5))
        ax = sns.scatterplot(x=xfield, y=yfield, hue="chromosome", data=df_cov, alpha=.2)
        ax.legend(bbox_to_anchor=(1, 1))
        
        # save
        fig.savefig("%s/perChrom_%s_vs_%s.pdf"%(plots_dir, xfield, yfield), bbox_inches='tight')
        plt.close(fig)


def get_compromised_ref_breakpoint(r):

    """Takes a row of a df_gridss and returns the reference breakpoint, in terms of graph posiition that is compromised by the breakpoint in case of being homozygous. So the returned value is a tuple of the adjacent positions.

    Take into consideration that this is the gridss notation:

    s t[p[ piece extending to the right of p is joined after t
    s t]p] reverse comp piece extending left of p is joined after t
    s ]p]t piece extending to the left of p is joined before t
    s [p[t reverse comp piece extending right of p is joined before t
    """

    graph_node = r["graph_node"]

    # ]p]t piece extending to the left of p is joined before t
    # ]mito_C_glabrata_CBS138:20054]TTA
    if r["ALT"].count("]")==2  and r["ALT"].startswith("]"): return (graph_node-1 , graph_node)

    # t]p] reverse comp piece extending left of p is joined after t
    # TTA]mito_C_glabrata_CBS138:20054]
    elif r["ALT"].count("]")==2  and r["ALT"].endswith("]"): return (graph_node , graph_node+1)

    # [p[t reverse comp piece extending right of p is joined before t
    # [mito_C_glabrata_CBS138:1841[T
    elif r["ALT"].count("[")==2  and r["ALT"].startswith("["): return (graph_node-1 , graph_node)

    # t[p[ piece extending to the right of p is joined after t
    # T[mito_C_glabrata_CBS138:1841[
    elif r["ALT"].count("[")==2  and r["ALT"].endswith("["): return (graph_node , graph_node+1)

    else: raise ValueError("The alt %s is  not properly formatted"%r["ALT"])



def rgb_to_hex(rgb):

    # Helper function to convert colour as RGB tuple to hex string

    rgb = tuple([int(255*val) for val in rgb])
    return '#' + ''.join([hex(val)[2:] for val in rgb]).upper()


def get_value_to_color(values, palette="mako", n=100, type_color="rgb", center=None):

    """TAkes an array and returns the color that each array has. Checj http://seaborn.pydata.org/tutorial/color_palettes.html"""

    # get the colors
    colors = sns.color_palette(palette, n)

    # change the colors
    if type_color=="rgb": colors = colors
    elif type_color=="hex": colors = [rgb_to_hex(c) for c in colors]
    else: raise ValueError("%s is not valid"%palette)

    # change

    # if they are strings
    if type(list(values)[0])==str:

        palette_dict = dict(zip(values, colors))
        value_to_color = palette_dict

    # if they are numbers
    else:

        # map eaqually distant numbers to colors
        if center==None:
            min_palette = min(values)
            max_palette = max(values)
        else: 
            max_deviation = max([abs(fn(values)-center) for fn in [min, max]])
            min_palette = center - max_deviation
            max_palette = center + max_deviation

        all_values_palette = list(np.linspace(min_palette, max_palette, n))
        palette_dict = dict(zip(all_values_palette, colors))

        # get value to color
        value_to_color = {v : palette_dict[find_nearest(all_values_palette, v)] for v in values}

    return value_to_color, palette_dict

def get_lowess_fit_y(x, y, xtarget, frac, iterations, fill_value_interpolation):

    """This function takes an x and a y. It returns the prediction on the targetx.

    fill_value_interpolation should be a tuple of the returned if x<min("""

    # get unique x and ys
    x_to_y_unique = pd.DataFrame({"x":x, "y":y}).groupby("x").apply(lambda df_x: np.mean(df_x.y)).sort_index()
    xfit = np.array(x_to_y_unique.index, dtype=float)
    yfit = np.array(x_to_y_unique.values, dtype=float)

    # get the loess fit
    lowess_results = cylowess.lowess(endog=yfit, exog=xfit, frac=frac, it=iterations)
    # delta is the: Distance within which to use linear-interpolation instead of weighted regression.

    # unpack
    sorted_lowess_x = lowess_results[:,0]
    lowess_y = lowess_results[:,1]

    # if there are NaNs, return them
    if any(pd.isna(lowess_y)): return [np.nan]*len(xtarget)

    # get the interploation function
    if fill_value_interpolation is None: bounds_error = True
    else: bounds_error = False

    # get the interpolation function
    interpolation_function = scipy_interpolate.interp1d(xfit, lowess_y, bounds_error=bounds_error, kind="linear", assume_sorted=True, fill_value=fill_value_interpolation)

    # get the prediction on the xfit
    ypredicted_on_xtarget = interpolation_function(xtarget)

    return ypredicted_on_xtarget


def get_LOWESS_benchmarking_series_CV(kfold, frac, it, df, xfield, yfield, min_test_points_CV, outdir):

    """This function takes a df with xfield and yfield. It runs kfold cross validation and returns a series with the accuracies """

    # define all the indices
    all_idx = set(df.index)

    # init rsquares
    rsquares_cv = []

    # this will only work if the unique df is long enough
    if len(df)>kfold and (frac*len(df))>=3: 

        # iterate through 10-fold chross validation 
        kfold_object = KFold(n_splits=kfold, random_state=1, shuffle=True)
        for numeric_train_index, numeric_test_index in kfold_object.split(df.index):

            print_if_verbose("processing data")

            # get the idx test as the index of df
            test_idx = set(df.iloc[numeric_test_index].index)
            train_idx = all_idx.difference(test_idx)

            # if there are not enough data points, break
            if len(test_idx)<min_test_points_CV: break

            # get dfs
            df_train = df.loc[train_idx].sort_values(by=[xfield, yfield])
            df_test = df.loc[test_idx].sort_values(by=[xfield, yfield])
           
            # get train alues
            xtrain = df_train[xfield].values
            ytrain = df_train[yfield].values

            # get the test values
            xtest = df_test[xfield].values
            ytest = df_test[yfield].values

            idx_correct_test = (xtest>min(xtrain)) & (xtest<max(xtrain))
            xtest = xtest[idx_correct_test]
            ytest = ytest[idx_correct_test]

            # if there are not enough points, skip
            if sum(idx_correct_test)<min_test_points_CV: break

            # get the lowess fit on the test data
            print_if_verbose("running loess")
            fill_value_interpolation = None
            ytest_predicted = get_lowess_fit_y(xtrain, ytrain, xtest, frac, it, fill_value_interpolation)

            # if there are NaNs, break the cross-validation
            if any(pd.isna(ytest_predicted)): break

            # if there are any 0 predicted values, break the cross-validation
            if any(ytest_predicted<=0): break

            # debug
            if len(ytest_predicted)!=len(ytest): raise ValueError("xtest and ytest are not the same")
            if any(pd.isna(ytest_predicted)): raise ValueError("There can't be NaNs")

            # calculate the rsquare, making sure it is a float
            rsquare = r2_score(ytest, ytest_predicted)

            # debug
            if pd.isna(rsquare): raise ValueError("The rsquare can't be nan")

            # break trying if there is a 0 rsquare
            if rsquare<=0: break 

            # keep
            rsquares_cv.append(rsquare)

    # discard if any rsquares are 0
    if len(rsquares_cv)!=kfold: 

        mean_rsquare = 0
        std = 1
        inverse_std_rsquare = 0

    else:

        mean_rsquare = np.mean(rsquares_cv)
        std = np.std(rsquares_cv)
        inverse_std_rsquare = 1/std

    # get the final series
    benchmarking_series = pd.Series({"frac":frac, "it":it, "mean_rsquare":mean_rsquare, "inverse_std_rsquare":inverse_std_rsquare, "std_rsquare":std, "kfold":kfold})

    return benchmarking_series
    

def get_y_corrected_by_x_LOWESS_crossValidation(df, xfield, yfield, outdir, threads, replace, plots_prefix, max_y, min_test_points_CV=10):

    """This function takes an x and a y series, returning the y corrected by x. This y corrected is y/(y predicted from LOWESS from x). The index must be unique. The best parameters are taken with 10 fold cross validation"""

    make_folder(outdir)

    # keep
    initial_index = list(df.index)
    df = cp.deepcopy(df)[[xfield, yfield]].sort_values(by=[xfield, yfield])

    # check that the index is unique
    if len(df.index)!=len(set(df.index)): raise ValueError("The index should be unique")


    # define the df_fitting as the one where the yfield is not 0
    df_fitting = df[(df[yfield]>0) & (df[yfield]<=max_y)]

    # sort by the x
    df_fitting = df_fitting.sort_values(by=[xfield, yfield])

    ########## GET THE DF BENCHMARKED DF 10xCV ########## 

    # define the df_benckmarking file
    df_benchmarking_file = "%s/df_benckmarking.tab"%outdir

    if file_is_empty(df_benchmarking_file) or replace is True:
        print_if_verbose("getting benchmarking for %s vs %s"%(xfield, yfield))

        # define parms
        n_frac = 8
        kfold = 4

        # define all the fractions
        min_frac = min([1/len(df_fitting), 0.05])
        all_fractions = list(np.linspace(min_frac, 0.1, n_frac))

        # define several robustifying iteration
        all_its = range(1, 3) 

        # debug
        if any(pd.isna(df_fitting[xfield])) or any(pd.isna(df_fitting[yfield])): raise ValueError("There are NaNs")

        # define the inputs of the benchmarking function
        inputs_fn = make_flat_listOflists([[(kfold, frac, it, df_fitting, xfield, yfield, min_test_points_CV, outdir) for frac in all_fractions] for it in all_its])

        # get a list of the benchmarking series in parallel
        with multiproc.Pool(threads) as pool:

            list_benchmarking_series = pool.starmap(get_LOWESS_benchmarking_series_CV, inputs_fn) 
            
            pool.close()
            pool.terminate()

        # get as df
        df_benchmarking = pd.DataFrame(list_benchmarking_series)

        # save
        save_df_as_tab(df_benchmarking, df_benchmarking_file)

    # load
    df_benchmarking  = get_tab_as_df_or_empty_df(df_benchmarking_file)

    ##################################################### 

    if len(df_benchmarking)==0 or max(df_benchmarking.mean_rsquare)<=0: 

        print("WARNING: There is not enough variability or data points to perform a correction of %s on %s. There will be no correction applied"%(yfield, xfield))
        y_corrected = df[yfield]
        final_rsquare = 0.0
        df["predicted_yvalues"] = np.median(df[yfield])

    else:

        # get the fit data
        print_if_verbose("performing LOWESS regression with best parameters for %s vs %s"%(xfield, yfield))

        # get sorted df
        df_fitting = df_fitting.sort_values(by=[xfield, yfield])

        # sort df benchmarking to get the max rsquare and minimum std
        max_kfold = max(df_benchmarking.kfold)
        df_benchmarking = df_benchmarking[df_benchmarking.kfold==max_kfold].sort_values(by=["mean_rsquare", "inverse_std_rsquare"], ascending=False)

        # get the best parameters
        best_parms_series = df_benchmarking.iloc[0]

        # get the y predicted with the best parms
        best_frac = best_parms_series["frac"]
        best_it = int(best_parms_series["it"])
        #outprefix = "%s/final_loess_fitting"%(outdir)

        # get the final fitting from training based on the df_fitting, but testing on the real df. The interpolation 
        fill_value_interpolation = "extrapolate"
        df["predicted_yvalues"] = get_lowess_fit_y(df_fitting[xfield].values, df_fitting[yfield].values, df[xfield].values, best_frac, best_it, fill_value_interpolation)

        # correct the predicted_yvalues so that if they are negative thet'd be set to 0 given that the input is also negative
        def get_predicted_yvalues(r):
            if r["predicted_yvalues"]<=0.0 and r[yfield]==0.0: return 0.0
            else: return r["predicted_yvalues"]

        df["predicted_yvalues"] = df.apply(get_predicted_yvalues, axis=1)

        # debug 
        if any(pd.isna(df.predicted_yvalues)): raise ValueError("there should be no NaNs in the final prediction")

        # debug if any of the predicted_yvalues is <=0
        if any(df.predicted_yvalues<0): raise ValueError("There can't be any negative values or less predicted yvalues")

        # calculate the final rsquare
        final_rsquare = r2_score(df[yfield], df.predicted_yvalues)
        if pd.isna(final_rsquare): raise ValueError("rsquare can't be NaN")

        ##############################

        ######### MAKE PLOTS #########

        filename = "%s_%s_from_%s_coverage.pdf"%(plots_prefix, yfield, xfield)
        if file_is_empty(filename) or replace is True:

            filename_tmp = "%s/coverage.tmp.pdf"%(outdir)

            # get the plot
            df_plot = df.sort_values(by=[xfield, yfield])

            fig = plt.figure(figsize=(5,5))

            #plt.plot(df_plot[xfield], df_plot[yfield], "o", alpha=0.2, color="gray", label="raw data")
            sns.kdeplot(df_plot[[xfield, yfield]], cmap="gist_gray", shade=True)
            plt.plot(df_plot[xfield], df_plot.predicted_yvalues, "-", color="red", label="LOWESS fit")

            plt.title("Fitting LOWESS with frac=%.3f it=%i. final R2=%.3f. %ix CV R2=%.3f +- %.3f (SD)\n"%(best_frac, best_it, final_rsquare, best_parms_series["kfold"], best_parms_series["mean_rsquare"], best_parms_series["std_rsquare"]))
            plt.legend(bbox_to_anchor=(1, 1))
            plt.xlabel(xfield)
            plt.ylim([0, np.percentile(df_plot[yfield], 95)])
            plt.ylabel(yfield)

            fig.savefig(filename_tmp, bbox_inches='tight')
            plt.close(fig)

            os.rename(filename_tmp, filename)

        filename = "%s_%s_from_%s_coverage_corrected.pdf"%(plots_prefix, yfield, xfield)
        if file_is_empty(filename) or replace is True:

            filename_tmp = "%s/coverage_corrected.tmp.pdf"%(outdir)

            # get the plot
            df_plot = df.sort_values(by=[xfield, yfield])

            # add the correction
            df_plot["y_corrected"] = df_plot[yfield]/df_plot.predicted_yvalues

            fig = plt.figure(figsize=(5,5))

            #plt.plot(df_plot[xfield], df_plot[yfield], "o", alpha=0.2, color="gray", label="raw data")
            sns.kdeplot(df_plot[[xfield, "y_corrected"]], cmap="gist_gray", shade=True)

            plt.title("Fitting LOWESS with frac=%.3f and it=%i. final R2=%.3f. %ix CV R2=%.3f +- %.3f (SD)\n"%(best_frac, best_it, final_rsquare, best_parms_series["kfold"], best_parms_series["mean_rsquare"], best_parms_series["std_rsquare"]))
            plt.legend(bbox_to_anchor=(1, 1))
            plt.xlabel(xfield)
            plt.ylim([0, np.percentile(df_plot[yfield], 95)])
            plt.ylabel("%s corrected by %s"%(yfield, xfield))

            fig.savefig(filename_tmp, bbox_inches='tight')
            plt.close(fig)

            os.rename(filename_tmp, filename)

        ############################    

        # get the corrected vals. If there is no prediction just return the raw vals
        def divide_with_noNaN_correction(r):

            # if the yfield is 0, return it as it is
            if r[yfield]==0 and r["predicted_yvalues"]==0: return 0.0

            # predicted yvalues can't be 0 unless yfield is also
            elif r["predicted_yvalues"]==0: raise ValueError("predicted_yvalues can't be 0 if yfield is not as well") 
            
            # normal division
            else: return r[yfield]/r["predicted_yvalues"]

        if final_rsquare>0: y_corrected = df.apply(divide_with_noNaN_correction, axis=1)
        else: y_corrected = df[yfield]

        # debug
        if any(pd.isna(y_corrected)): raise ValueError("there should be no NaNs in y_corrected")

    # get in the order of the initial index
    df = df.loc[initial_index]
    y_corrected = y_corrected[initial_index]

    # return in the initial order
    return y_corrected, final_rsquare, df.predicted_yvalues


def verify_no_NaNs(series):

    """It throws an error if there are NaNs"""

    if any(pd.isna(series)): raise ValueError("There should be no NaNs")



def get_df_coverage_with_corrected_coverage(df_coverage, reference_genome, outdir, replace, threads, mitochondrial_chromosome, df_gridss):

    """This function will take a df_coverage that has coverage_field as a proxy for coverage. It will add <coverage_field> which is a value that will be a ratio between the coverage_field and the coverage_field predicted from a loess regression taking into account mappability, GC content and distance to the telomere across the windows. The resulting value will be centered arround 1.  """

    # define the initial cols
    initial_cols = list(df_coverage.columns)

    # define the outfile_final
    outfile_final = "%s/df_coverage_with_corrected_relative_coverage.tab"%outdir

    # define the working dir
    working_outdir = "%s/working_dir"%outdir; make_folder(working_outdir)

    # define the results and plots dir, where all the plots of the fits will be saved
    results_dir = "%s/calculating_corrected_coverage"%outdir; make_folder(results_dir)
    plots_dir = "%s/plots"%results_dir; make_folder(plots_dir)

    # define the rsquares
    outfile_rsquares = "%s/rsquares_tables.tab"%results_dir

    if file_is_empty(outfile_final) or file_is_empty(outfile_rsquares) or replace is True:

        # check content
        if any(df_coverage.start>=df_coverage.end): raise ValueError("start can't be after end")

        # add "relative_coverage", which will include coverage_field
        if "relative_coverage" in set(df_coverage.keys()): raise ValueError("coverage can't be in the df keys")
        median_coverage = get_median_coverage(df_coverage, mitochondrial_chromosome)
        df_coverage["relative_coverage"] = df_coverage["mediancov_1"]/median_coverage

        # add the GC content
        gcontent_outfile = "%s/df_coverage_with_gccontent.py"%working_outdir
        df_coverage = get_df_with_GCcontent(df_coverage, reference_genome, gcontent_outfile, replace=replace)

        # add the median mappability
        mappability_outfile = "%s/df_coverage_with_mappability.tab"%working_outdir
        df_coverage = get_df_windows_with_median_mappability(df_coverage, reference_genome, mappability_outfile, replace, threads)

        # add the raw distance to the telomere, in linear space
        chr_to_len = get_chr_to_len(reference_genome)
        df_coverage["middle_position"] = (df_coverage.start + (df_coverage.end - df_coverage.start)/2).apply(int)
        df_coverage["raw_distance_to_telomere"] = df_coverage.apply(lambda r: min([r["middle_position"], chr_to_len[r["chromosome"]]-r["middle_position"]-1]), axis=1)

        # define chroms
        all_chromosomes = set(get_chr_to_len(reference_genome))
        if mitochondrial_chromosome!="no_mitochondria": mtDNA_chromosomes = set(mitochondrial_chromosome.split(","))
        else: mtDNA_chromosomes = set()
        gDNA_chromosomes = all_chromosomes.difference(mtDNA_chromosomes)

        # define a unique index
        df_coverage.index = list(range(0, len(df_coverage)))

        # init the final df_coverage
        final_df_coverage = pd.DataFrame()

        # init the final df_rsquares
        final_df_rsquares = pd.DataFrame()

        # iterate through each genome
        for type_genome, chroms in [("gDNA", gDNA_chromosomes), ("mtDNA", mtDNA_chromosomes)]:
            print_if_verbose("investigating %s"%type_genome)

            # define an outdir for this type of genome
            outdir_type_genome = "%s/%s"%(working_outdir, type_genome); make_folder(outdir_type_genome)

            # get the df coverage of this genome
            df_cov = df_coverage[df_coverage.chromosome.isin(chroms)]
            if len(df_cov)==0: continue

            # make some plots with all these values
            #plots_dir = "%s/plots_coverage_%s"%(working_outdir, type_genome); make_folder(plots_dir)
            #make_plots_coverage_parameters(df_cov, plots_dir)

            # init the predictor fields of coverage
            predictor_fields = ["GCcontent", "median_mappability"]

            # add the distance to the telomere if necessary

            # first check whether there is a correlation between distance to the telomere and only calculate correlation if so
            r_spearman, p_spearman = stats.spearmanr(df_cov.raw_distance_to_telomere, df_cov.relative_coverage, nan_policy="raise")

            if p_spearman<0.05: predictor_fields.append("raw_distance_to_telomere")

            ######## CALCULATE THE RSQUARE OF EACH PREDICTOR ALONE ########

            # define the maximum relative coverage on which to base the fitting
            max_relative_coverage = np.median(df_cov[df_cov.relative_coverage>0].relative_coverage)*4

            # calculate the rsquare for each predictor
            calculate_rsquares_dir = "%s/calculating_rsquares_each_predictor_%s"%(working_outdir, type_genome); make_folder(calculate_rsquares_dir)
            predictor_to_rsquare = {}
            for p in predictor_fields:

                # fit the data
                plots_prefix = "%s/single_predictors"%plots_dir
                lowess_dir_p = "%s/%s"%(calculate_rsquares_dir, p)
                y_corrected, rsquare, y_predicted = get_y_corrected_by_x_LOWESS_crossValidation(df_cov, p, "relative_coverage", lowess_dir_p, threads, replace, plots_prefix, max_relative_coverage)

                # add the correction based on the p
                df_cov["relative_coverage_corrected_by_%s"%p] = y_corrected

                # add the prediction based on p
                df_cov["relative_coverage_predicted_from_%s"%p] = y_predicted

                # add the rsquare
                predictor_to_rsquare[p] = rsquare

            # get as series
            predictor_to_rsquare = pd.Series(predictor_to_rsquare)

            ################################################################

            ######## PERFORM THE FINAL FITTING ########

            # init the corrected coverage  as divided by the median coverage
            df_cov["median_relative_coverage"] = np.median(df_cov[df_cov.relative_coverage>0].relative_coverage)
            df_cov["corrected_relative_coverage"] = df_cov.relative_coverage / df_cov.median_relative_coverage

            # add the rsquare without predictiopn
            predictor_to_rsquare["no_prediction"] = r2_score(df_cov.relative_coverage, df_cov.median_relative_coverage)

            # go through each final predictor
            for pID, predictor in enumerate(predictor_fields):

                # correct the coverage 
                plots_prefix = "%s/final_fitting_round%i"%(plots_dir, pID+1)
                outdir_lowess = "%s/final_fitting_%s_round%i"%(calculate_rsquares_dir, predictor, pID+1)
                df_cov["corrected_relative_coverage"], rsquare = get_y_corrected_by_x_LOWESS_crossValidation(df_cov, predictor, "corrected_relative_coverage", outdir_lowess, threads, replace, plots_prefix, max_relative_coverage)[0:2]

                # add the rsquare
                predictor_to_rsquare["final_rsquare_round%i_from_%s"%(pID+1, predictor)] = rsquare

            # save the predictor_to_rsquare into results
            df_rsquares = pd.DataFrame({"rsquare":predictor_to_rsquare})
            df_rsquares["type_fit"] = df_rsquares.index
            df_rsquares["type_genome"] = type_genome

            #############################################

            # at the end add
            final_df_coverage = final_df_coverage.append(df_cov)
            final_df_rsquares = final_df_rsquares.append(df_rsquares)

        # save
        save_df_as_tab(final_df_coverage, outfile_final)
        save_df_as_tab(final_df_rsquares, outfile_rsquares)

    # load the dfs
    df_coverage = get_tab_as_df_or_empty_df(outfile_final)
    df_rsquares = get_tab_as_df_or_empty_df(outfile_rsquares)

    ############# PLOT COVERAGE #############

    # add an offset position
    df_plot = df_coverage.sort_values(by=["chromosome", "start", "end"])
    df_plot["xposition_plot"] = list(range(0, len(df_plot)))

    # init fig
    fig = tools.make_subplots(rows=1, cols=1, specs=[[{}]], vertical_spacing=0.0, horizontal_spacing=0.0, subplot_titles=(""), shared_yaxes=True, shared_xaxes=True, print_grid=True)

    # init
    all_x = []
    all_colors = []

    all_y_rel_cov = []
    all_y_corrected_rel_cov = []
    all_y_GCcontent = []
    all_y_relDistanceToTelomere = []

    # go through each chromosome and create vals
    sorted_chroms = sorted(set(df_plot.chromosome))
    chrom_to_color = get_value_to_color(sorted_chroms, palette="tab10", n=len(sorted_chroms), type_color="hex")[0]
    current_offset = 0

    for chrom in sorted_chroms:

        # get plot
        df_chrom = df_plot[df_plot.chromosome==chrom]

        # add
        all_x += (list(df_chrom.xposition_plot + current_offset) + [None])
        all_y_rel_cov += (list(df_chrom.relative_coverage) + [None])
        all_y_corrected_rel_cov += (list(df_chrom.corrected_relative_coverage) + [None])
        all_y_GCcontent += (list(df_chrom.GCcontent) + [None])
        all_y_relDistanceToTelomere += (list(df_chrom.raw_distance_to_telomere/np.median(df_chrom.raw_distance_to_telomere)) + [None])

        current_offset += int(len(df_chrom)/2)

    # get scatters

    # get the relative coverage
    fig.append_trace(go.Scatter(x=all_x, y=all_y_rel_cov, showlegend=True, mode="lines", opacity=1, hoveron="points+fills", name="Relative coverage", line=dict(color="blue", width=2, dash="dash")) , 1, 1) 

    # add the corrected coverage
    fig.append_trace(go.Scatter(x=all_x, y=all_y_corrected_rel_cov, showlegend=True, mode="lines", opacity=1, hoveron="points+fills", name="Relative coverage corrected", line=dict(color="red", width=2, dash=None)) , 1, 1) 

    # add the GC content 
    fig.append_trace(go.Scatter(x=all_x, y=all_y_GCcontent, showlegend=True, mode="lines", opacity=1, hoveron="points+fills", name="GC content", line=dict(color="green", width=2, dash=None)) , 1, 1) 

    # add the relative distance to the telomere
    fig.append_trace(go.Scatter(x=all_x, y=all_y_relDistanceToTelomere, showlegend=True, mode="lines", opacity=1, hoveron="points+fills", name="Relative (raw) distance to telomere", line=dict(color="cyan", width=2, dash=None)) , 1, 1) 

    # get the CN
    #fig.append_trace(go.Scatter(x=list(df_perWindow.start), y=list(df_perWindow.relative_CN), showlegend=True, mode="lines", line=dict(color="red", width=2, dash="dash"), opacity=1, hoveron="points+fills", name="CONY CN"), 1, 1) 

    # get figure
    fig['layout'].update(title="relative and corrected coverage", margin=go.Margin(l=250, r=250, b=50, t=50), showlegend=True)
    config={'editable': False}
    off_py.plot(fig, filename="%s/coverage_interactive.html"%(plots_dir), auto_open=False, config=config)

    #########################################

    # remove the working dir
    delete_folder(working_outdir)

    # return 
    return df_coverage

     
def run_fraction_Nbases_window_genome(r, chrom_to_seq):

    """Takes a row of a df_windwos and returns the fraction of characters that are undefinable"""

    # get the number of Ns
    n_Ns = chrom_to_seq[r["chromosome"]][r["start"]:r["end"]].count("N")

    # return divided by width
    return n_Ns/r["width"]


def get_chrom_to_Xoffset_plot(reference_genome):

    """Rteurns a dict mapping each chromosome to an X offset"""

    # get all chroms
    chrom_to_len = get_chr_to_len(reference_genome)
    all_chromosomes = sorted(chrom_to_len)


    chrom_to_Xoffset = {}
    current_offset = 0
    for chrom in all_chromosomes:
        chrom_to_Xoffset[chrom] = current_offset
        current_offset += chrom_to_len[chrom] + 15000

    return all_chromosomes, chrom_to_Xoffset

def get_df_CNV_with_metadata_from_df_coverage(df_coverage, df_CNV, relative_CN_fields):

    """This function takes a df coverage and a df_CNV and returns a series with metadata from the df_CNV to the df_coverage"""

    # sort
    df_CNV = df_CNV.sort_values(by=["chromosome", "start", "end"])

    # add IDs
    df_CNV["CNVid"] = list(range(len(df_CNV)))
    df_CNV.index = df_CNV.CNVid

    # define a function that returns all the metadata for one window
    def get_metadata_for_df_CNV_r(r):

        # get the coverage related to this window
        df_cov = df_coverage[(df_coverage.chromosome==r["chromosome"]) & (df_coverage.start>=r["start"]) & (df_coverage.end<=r["end"])]

        # check that it makes sense
        if min(df_cov.start)!=r["start"] or max(df_cov.end)!=r["end"]: raise ValueError("invalid parsing")

        # get several things
        final_dict = {"median_coverage" : np.median(df_cov.relative_coverage),
                      "median_coverage_corrected" : np.median(df_cov.corrected_relative_coverage)}

        for field_CN in relative_CN_fields: final_dict["median_%s"%field_CN] = np.median(df_cov[field_CN])

        # return as a series
        final_series = pd.Series(final_dict)

        return final_series

    fields_to_add = ["median_coverage", "median_coverage_corrected"] + ["median_%s"%x for x in relative_CN_fields]
    df_CNV[fields_to_add] = df_CNV.apply(get_metadata_for_df_CNV_r, axis=1)

    # add whether it is a duplication or a deletion
    def get_type_CNV_from_coverageState(x):
        if x<1.0: return "DEL"
        elif x>1.0: return "DUP"
        else: raise ValueError("invalid coverage")

    df_CNV["SVTYPE"] = df_CNV.merged_relative_CN.apply(get_type_CNV_from_coverageState)

    return df_CNV

def get_df_coverage_with_relative_coverage_and_for_each_typeGenome(df_coverage, reference_genome, mitochondrial_chromosome):

    """This function adds the relative coverage to df_coverage and also the relative to gDNA or mtDNA"""

    # add the relative coverage
    median_coverage = get_median_coverage(df_coverage, mitochondrial_chromosome)
    df_coverage["relative_coverage"] = df_coverage["mediancov_1"]/median_coverage

    # define chroms
    all_chromosomes = set(get_chr_to_len(reference_genome))
    if mitochondrial_chromosome!="no_mitochondria": mtDNA_chromosomes = set(mitochondrial_chromosome.split(","))
    else: mtDNA_chromosomes = set()
    gDNA_chromosomes = all_chromosomes.difference(mtDNA_chromosomes)

    # add the type of genome
    def get_type_genome(c):
        if c in gDNA_chromosomes: return "gDNA"
        elif c in mtDNA_chromosomes: return "mtDNA"
        else: raise ValueError("%s is not in gDNA nor mtDNA"%c)
    df_coverage["type_genome"] = df_coverage.chromosome.apply(get_type_genome)

    # calculate the median of each genome
    type_genome_to_median_relativeCov = {type_genome : get_median_coverage(df_coverage[df_coverage.type_genome==type_genome], mitochondrial_chromosome, coverage_field="relative_coverage") for type_genome in ["gDNA", "mtDNA"] if sum(df_coverage.type_genome==type_genome)>0}

    # add to the df the relative_coverage_to_genome
    df_coverage["median_relative_coverage_genome"] = df_coverage.type_genome.apply(lambda x: type_genome_to_median_relativeCov[x])
    df_coverage["relative_coverage_to_genome"] = df_coverage.relative_coverage / df_coverage.median_relative_coverage_genome

    return df_coverage

def get_df_coverage_with_corrected_coverage_background(df_coverage, df_coverage_bg, reference_genome, outdir, replace, threads, mitochondrial_chromosome):

    """This function takes a df_coverage and a df_coverage of the background. It returns the same df_coverage with the 'corrected_relative_coverage' which is the correction by the background. """

    # debug inputs
    if len(df_coverage)!=len(df_coverage_bg): raise ValueError("the df coverages are not the same")

    # define the initial cols
    initial_cols = list(df_coverage.columns)

    # define the outfile_final
    outfile_final = "%s/df_coverage_with_corrected_relative_coverage.tab"%outdir

    # define the working dir
    working_outdir = "%s/working_dir"%outdir; make_folder(working_outdir)

    if file_is_empty(outfile_final) or replace is True:

        # check content
        if any(df_coverage.start>=df_coverage.end): raise ValueError("start can't be after end")

        # add "relative_coverage", which will include coverage_field
        if "relative_coverage" in set(df_coverage.keys()): raise ValueError("coverage can't be in the df keys")

        # add the GC content
        gcontent_outfile = "%s/df_coverage_with_gccontent.py"%working_outdir
        df_coverage = get_df_with_GCcontent(df_coverage, reference_genome, gcontent_outfile, replace=replace)

        # add the median mappability
        mappability_outfile = "%s/df_coverage_with_mappability.tab"%working_outdir
        df_coverage = get_df_windows_with_median_mappability(df_coverage, reference_genome, mappability_outfile, replace, threads)

        # add the relative coverage and relative coverage to the genome to each df
        df_coverage = get_df_coverage_with_relative_coverage_and_for_each_typeGenome(df_coverage, reference_genome, mitochondrial_chromosome)
        df_coverage_bg = get_df_coverage_with_relative_coverage_and_for_each_typeGenome(df_coverage_bg, reference_genome, mitochondrial_chromosome)

        # define a unique index
        df_coverage.index = list(range(0, len(df_coverage)))
        df_coverage_bg.index = list(range(0, len(df_coverage_bg)))

        # add the realtive coverage to the genome of the background sample
        df_coverage = df_coverage.merge(df_coverage_bg[["chromosome", "start", "end", "relative_coverage_to_genome"]], on=["chromosome", "start", "end"], how="left", validate="one_to_one", suffixes=("", "_background"))

        # add the corrected_relative_coverage so by normalising the coverage relative to the genome by the same one in the background. If the coverage in the bg is 0 we'll also set the coverage of the current to 0. The corrected relative coaverage will be maxed at 10.0
        def get_corrected_relative_coverage_bg(r):

            # get the vals
            cov = r.relative_coverage_to_genome
            cov_bg = r.relative_coverage_to_genome_background

            if cov_bg>0: return cov/cov_bg
            elif cov_bg==0: return 0.0 # if the reference is deleted, set also as deleted the relative coverage
            else: raise ValueError("the coverages %.2f and %.2f are not valid"%(cov, cov_bg))

        def set_to_max_cov(x): return min([10.0, x])
        df_coverage["corrected_relative_coverage"] = df_coverage.apply(get_corrected_relative_coverage_bg, axis=1).apply(set_to_max_cov)

        # save
        save_df_as_tab(df_coverage, outfile_final)

    # load
    df_coverage = get_tab_as_df_or_empty_df(outfile_final)

    # remove the working dir
    delete_folder(working_outdir)

    # return 
    return df_coverage



def run_CNV_calling(sorted_bam, reference_genome, outdir, threads, replace, mitochondrial_chromosome, df_gridss, window_size, ploidy, plot_correlation=True, bg_sorted_bam_CNV=None, cnv_calling_algs={"HMMcopy", "AneuFinder", "CONY"}):

    """This function takes a sorted bam and runs several programs on it to get the copy-number variation results. It is important that the sorted_bam contains no duplicates. It will correct bu GC content, mappability and distance to the telomere. All coverage will be corrected by GC content, mappability and the distance to the telomere, which will be calculated also taking into account breakpoint information. 

    If bg_sorted_bam_CNV is provided the correction will be performed by this sorted bam instead of the position, GC content and mappability."""

    make_folder(outdir)

    # define the final file
    final_CNV_file = "%s/final_CNVcalling.tab"%outdir
    final_df_coverage_file = "%s/final_df_coverage.tab"%outdir

    if file_is_empty(final_CNV_file) or file_is_empty(final_df_coverage_file) or replace is True:

        # init the files to remove
        files_folders_remove = []

        ############ GET A DF WITH CORRECTED COVERAGE ##############

        # make a df with windows of the genome
        df_coverage = get_coverage_df_for_windows_of_genome(sorted_bam, reference_genome, outdir, replace, threads, window_size)

        # add the files to remove
        files_folders_remove.append("%s/coverage_per_windows_%ibp.tab"%(outdir, window_size))

        # add the 'corrected_relative_coverage' by mappability, GC content and distance to the telomere
        if bg_sorted_bam_CNV is None:

            df_coverage = get_df_coverage_with_corrected_coverage(df_coverage, reference_genome, outdir, replace, threads, mitochondrial_chromosome, df_gridss)

        # if bg_sorted_bam_CNV is provided, calculate the "corrected_relative_coverage" as compared to the bg_sorted_bam_CNV
        else:

            # get the df coverage of the background
            outdir_calculating_coverage_bg = "%s/calculating_bg_coverage"%outdir; make_folder(outdir_calculating_coverage_bg)
            df_coverage_bg = get_coverage_df_for_windows_of_genome(bg_sorted_bam_CNV, reference_genome, outdir_calculating_coverage_bg, replace, threads, window_size)

            # get the df coverage with the corrected_relative_coverage (relative to bg), mappability, GC content
            df_coverage = get_df_coverage_with_corrected_coverage_background(df_coverage, df_coverage_bg, reference_genome, outdir, replace, threads, mitochondrial_chromosome)

            files_folders_remove += [outdir_calculating_coverage_bg, "%s/df_coverage_with_corrected_relative_coverage.tab"%outdir]


        # check that there are no NaNs
        if any(pd.isna(df_coverage.corrected_relative_coverage)): raise ValueError("There should be no NaNs in 'corrected_relative_coverage' ")

        # add the fraction of N bases
        chrom_to_seq = {seq.id : str(seq.seq).upper() for seq in SeqIO.parse(reference_genome, "fasta")}
        df_coverage["width"] = (df_coverage.end - df_coverage.start).apply(int)
        df_coverage["fraction_N_bases"] = df_coverage.apply(run_fraction_Nbases_window_genome, chrom_to_seq=chrom_to_seq, axis=1)

        # debug
        if len(set(cnv_calling_algs).difference({"AneuFinder", "HMMcopy", "CONY"}))>0: raise ValueError("the cnv_calling_algs %s are not correct"%cnv_calling_algs)

        ############################################################

        ######## RUN ANEUFINDER ########

        # init the cnv_calling fields
        relative_CN_fields = []
        
        if "AneuFinder" in cnv_calling_algs:

            # this is from a 2020 paper of CNV and antifungal drugs in C. albicans

            AneuFinder_outdir = "%s/AneuFinder_run"%outdir
            df_CNperWindow_AneuFinder_file = "%s/df_CNperWindow_AneuFinder_final.tab"%outdir

            if file_is_empty(df_CNperWindow_AneuFinder_file) or replace is True:

                # define the threads
                threads_aneufinder = 1 # note that the aneufinder threads is 1 to improve performance

                # calculate the read length
                read_length = get_read_length(sorted_bam, threads=threads, replace=replace)

                df_CNperWindow_AneuFinder = run_CNV_calling_AneuFinder(AneuFinder_outdir, replace, threads_aneufinder, df_coverage, ploidy, reference_genome, mitochondrial_chromosome, read_length)

                save_df_as_tab(df_CNperWindow_AneuFinder, df_CNperWindow_AneuFinder_file)

            df_CNperWindow_AneuFinder = get_tab_as_df_or_empty_df(df_CNperWindow_AneuFinder_file)

            print(df_CNperWindow_AneuFinder)

            df_coverage = df_coverage.merge(df_CNperWindow_AneuFinder[["chromosome", "start", "end", "relative_CN"]], on=["chromosome", "start", "end"], left_index=False, right_index=False, how="left", validate="one_to_one").rename(columns={"relative_CN":"relative_CN_AneuFinder"})

            # keep
            relative_CN_fields.append("relative_CN_AneuFinder")
            files_folders_remove += [AneuFinder_outdir, df_CNperWindow_AneuFinder_file]

        ################################

        ######## RUN CONY ########

        if "CONY" in cnv_calling_algs:

            # This is a recent paper that shows very high accuracy

            # get chrom to len
            chrom_to_len = get_chr_to_len(reference_genome)

            # define the sample name
            sample_name = get_sample_name_from_bam(sorted_bam)

            # run CONY for each chromosome
            df_CNperWindow_CONY_file = "%s/df_CNperWindow_CONY.tab"%outdir

            if file_is_empty(df_CNperWindow_CONY_file) or replace is True:

                # start by the shortest chromosomes
                chroms_sorted_by_len = sorted(set(chrom_to_len), key=(lambda c: chrom_to_len[c]))

                df_CNperWindow_CONY = pd.concat([run_CNV_calling_CONY_one_chromosome(reference_genome, "%s/%s_CONYrun"%(outdir, c), c, replace, window_size, threads, chrom_to_len[c], sample_name, df_coverage, ploidy) for c in chroms_sorted_by_len])

                save_df_as_tab(df_CNperWindow_CONY, df_CNperWindow_CONY_file)

            df_CNperWindow_CONY = get_tab_as_df_or_empty_df(df_CNperWindow_CONY_file).rename(columns={"seqname":"chromosome"})

            # add to df coverage
            df_coverage = df_coverage.merge(df_CNperWindow_CONY[["chromosome", "start", "end", "relative_CN"]], on=["chromosome", "start", "end"], left_index=False, right_index=False, how="left", validate="one_to_one").rename(columns={"relative_CN":"relative_CN_CONY"})

            relative_CN_fields.append("relative_CN_CONY")
            files_folders_remove += (["%s/%s_CONYrun"%(outdir, c) for c in chrom_to_len] + [df_CNperWindow_CONY_file])

        ##########################

        ####### RUN HMMCOPY, SIMILARLY TO A RECENT C. glabrata PAPER #######

        if "HMMcopy" in cnv_calling_algs:

            # This is based on "Understand the genomic diversity and evolution of fungal pathogen Candida glabrata by genome-wide analysis of genetic variations"

            HMMcopy_outdir = "%s/HMMcopy_run"%outdir
            df_CNperWindow_HMMcopy_file = "%s/df_CNperWindow_HMMcopy_final.tab"%outdir

            if file_is_empty(df_CNperWindow_HMMcopy_file) or replace is True:

                df_CNperWindow_HMMcopy = run_CNV_calling_HMMcopy(HMMcopy_outdir, replace, threads, df_coverage, ploidy, reference_genome, mitochondrial_chromosome)

                save_df_as_tab(df_CNperWindow_HMMcopy, df_CNperWindow_HMMcopy_file)

            df_CNperWindow_HMMcopy = get_tab_as_df_or_empty_df(df_CNperWindow_HMMcopy_file)

            df_coverage = df_coverage.merge(df_CNperWindow_HMMcopy[["chromosome", "start", "end", "relative_CN"]], on=["chromosome", "start", "end"], left_index=False, right_index=False, how="left", validate="one_to_one").rename(columns={"relative_CN":"relative_CN_HMMcopy"})

            relative_CN_fields.append("relative_CN_HMMcopy")
            files_folders_remove += [HMMcopy_outdir, df_CNperWindow_HMMcopy_file]

        if len(relative_CN_fields)==0: raise ValueError("you should have some values in relative_CN_fields")

        ####################################################################

        ###### PLOT THE CALLING OF ALL PROGRAMS ######

        PlotsDir = "%s/plots"%outdir; make_folder(PlotsDir)
        df_plot = cp.deepcopy(df_coverage)

        # add the plot position to the df_plot
        all_chromosomes, chrom_to_Xoffset = get_chrom_to_Xoffset_plot(reference_genome)
        df_plot["Xoffset_plot"] = df_plot.chromosome.apply(lambda c: chrom_to_Xoffset[c])
        df_plot["plot_positionX"] = df_plot.start + df_plot.Xoffset_plot

        # define graphic properties
        cnvAlg_to_color = {"HMMcopy":"red", "CONY":"green", "AneuFinder":"black"}

        # init fig
        fig = tools.make_subplots(rows=1, cols=1, specs=[[{}]], vertical_spacing=0.0, horizontal_spacing=0.0, subplot_titles=(""), shared_yaxes=True, shared_xaxes=True, print_grid=True)

        # get data
        all_Xs = []
        all_relative_coverage = []
        cnv_alg_to_all_CN = {cnv_alg : [] for cnv_alg in cnv_calling_algs}

        for chrom in all_chromosomes:

            # get the df of the chrom
            df_c = df_plot[df_plot.chromosome==chrom]

            # keep the x and the relative coverage
            all_Xs += (list(df_c.plot_positionX) + [None])
            all_relative_coverage += (list(df_c.corrected_relative_coverage) + [None])

            # one line for each cnv calling
            for cnv_alg in cnv_calling_algs: cnv_alg_to_all_CN[cnv_alg] += (list(df_c["relative_CN_%s"%cnv_alg]) + [None])
        
        # plot the relative coverage
        fig.append_trace(go.Scatter(x=all_Xs, y=all_relative_coverage, showlegend=True, mode="lines+markers", marker=dict(symbol="circle", color="blue", size=4), opacity=1, hoveron="points+fills", name="corrected read depth", line=dict(color="blue", width=2, dash=None)) , 1, 1) 

        # add the CN of each alg
        for cnv_alg, all_CN in cnv_alg_to_all_CN.items():

            fig.append_trace(go.Scatter(x=all_Xs, y=all_CN, showlegend=True, mode="lines", line=dict(color=cnvAlg_to_color[cnv_alg], width=2, dash="dash"), opacity=1, hoveron="points+fills", name="%s CN"%cnv_alg), 1, 1) 

        # get figure
        fig['layout'].update(title="relative coverage and CN", margin=go.Margin(l=250, r=250, b=50, t=50), showlegend=True)
        config={'editable': False}
        off_py.plot(fig, filename="%s/CNcalling_interactive.html"%(PlotsDir), auto_open=False, config=config)
        
        ##############################################

        ######### PLOT CORRELATION #########

        if plot_correlation is True and len(relative_CN_fields)>=2:

            print_if_verbose("plotting predictions")

            # define a plot df
            df_plot = cp.deepcopy(df_coverage)

            for xfield, yfield in itertools.combinations(relative_CN_fields, 2):

                # get only the df_plot which is under some CNV
                df_plot = df_plot[(df_plot[xfield]!=1.0) | (df_plot[yfield]!=1.0)]

                # get the plot and jitter. This sets to a maximum of 3
                df_plot = df_plot.sort_values(by=[xfield, yfield])
                def set_to_max3(x): return min([x, 3])
                df_plot[xfield] = df_plot[xfield].apply(set_to_max3) + np.random.uniform(-1, 1, len(df_plot))*0.2
                df_plot[yfield] = df_plot[yfield].apply(set_to_max3) + np.random.uniform(-1, 1, len(df_plot))*0.2

                fig = plt.figure(figsize=(5,5))

                plt.plot(df_plot[xfield], df_plot[yfield], "o", alpha=0.05, color="gray", label="CN by each program")

                max_val = max([max(df_plot[xfield]), max(df_plot[yfield])])
                lims = [0, max_val]

                plt.legend(bbox_to_anchor=(1, 1))
                plt.xlim(lims)
                plt.ylim(lims)
                plt.xlabel("CN predicted by %s"%xfield)
                plt.ylabel("CN predicted by %s"%yfield)

                # add line
                plt.plot(np.linspace(0, max_val, 3), np.linspace(0, max_val, 3), "--", color="red", linewidth=.9)

                filename = "%s/CNtwoPrograms_%s_vs_%s.pdf"%(PlotsDir, xfield, yfield)
                fig.savefig(filename, bbox_inches='tight')
                plt.close(fig)

        ####################################

        ############ MERGE PREDICTIONS ###########

        print_if_verbose("merging predictions")
        
        # get the most conservative relative_CN
        def get_closest_to_relative_CN(r): return find_nearest(r, 1.0)
        df_coverage["merged_relative_CN"] = df_coverage[relative_CN_fields].apply(get_closest_to_relative_CN, axis=1)

        # correct the relative_CN so that if the relative_coverage is 0 it is also 0. Sometimes the programs don't work if there is few variability in the data (like in the WT).
        def get_corrected_merged_relative_CN_according_to_relative_coverage(r):
            if r["corrected_relative_coverage"] in {0.0}: return r["corrected_relative_coverage"]
            else: return r["merged_relative_CN"]

        df_coverage["merged_relative_CN"] = df_coverage.apply(get_corrected_merged_relative_CN_according_to_relative_coverage, axis=1)

        print_if_verbose("There are %i/%i windows of the genome under CNV"%(sum(df_coverage.merged_relative_CN!=1.0), len(df_coverage)))

        # get the df with the events of CNV
        df_coverage = df_coverage.sort_values(by=["chromosome", "start", "end"])
        df_coverage["windowID"] = list(range(len(df_coverage))) 
        df_coverage = df_coverage.set_index("windowID", drop=False)

        df_CNV = pd.DataFrame()
        for chrom in sorted(set(df_coverage.chromosome)):
            for CNstate in set(df_coverage.merged_relative_CN).difference({1.0}):

                # get the df of this chromosome under this specific CNstate
                df_chrom_CN = df_coverage[(df_coverage.merged_relative_CN==CNstate) & (df_coverage.chromosome==chrom)]
                if len(df_chrom_CN)==0: continue

                # define the CNVid, which is an ID for adjacent windows
                CNVids = []
                previous_wID = df_chrom_CN.index[0] - 1
                previous_CNVid = 0
                for wID, r in df_chrom_CN.iterrows():

                    # the window is the same
                    if wID==(previous_wID+1): CNVid = previous_CNVid                    

                    # it is a new window
                    else: CNVid = previous_CNVid + 1

                    # keep
                    CNVids.append(CNVid)

                    # define the previous things
                    previous_CNVid = CNVid
                    previous_wID = wID

                # add to the df
                df_chrom_CN["CNVid"] = CNVids

                # get a df 
                df_chrom_CN = df_chrom_CN.groupby("CNVid").apply(lambda df_cn: pd.Series({"chromosome":chrom, "merged_relative_CN":CNstate, "start":min(df_cn.start), "end":max(df_cn.end)}))

                # keep
                df_CNV = df_CNV.append(df_chrom_CN)


        # add metadata
        df_CNV = get_df_CNV_with_metadata_from_df_coverage(df_coverage, df_CNV, relative_CN_fields)

        ############################################

        # clean
        for f in files_folders_remove: 
            remove_file(f)
            delete_folder(f)

        # save
        save_df_as_tab(df_coverage, final_df_coverage_file)
        save_df_as_tab(df_CNV, final_CNV_file)

    # load
    df_CNV = get_tab_as_df_or_empty_df(final_CNV_file)

    return df_CNV
    
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


def write_repeats_table_file(repeats_table_file):

    """This function writes an empty repeats table file"""

    repeats_fields = ["IDrepeat", "SW_score", "begin_repeat", "chromosome", "end_repeat", "left_positionINrepeat", "left_repeat", "perc_del", "perc_div", "perc_ins", "position_inRepeat_begin", "position_inRepeat_end", "repeat", "strand", "type"]

    open(repeats_table_file, "w").write("\t".join(repeats_fields) + "\n")


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

    # check that the ID is unique
    if len(vcf_df)!=len(set(vcf_df.ID)): raise ValueError("ID is not unique")

    # get the repeats as a df
    repeats_df = pd.read_csv(repeats_table, sep="\t").rename(columns={"begin_repeat":"start", "end_repeat":"end"})[["chromosome", "start", "end"]]

    # get the IDs of the vcf df that overlap repeats
    outdir_varIDs = "%s_getting_overlap_with_vcf"%repeats_table; make_folder(outdir_varIDs)
    varIDs_overlapping_repeats = get_varIDs_overlapping_target_regions(vcf_df, repeats_df, outdir_varIDs)
    delete_folder(outdir_varIDs)

    # get the overlaps
    vcf_df["overlaps_repeats"] = vcf_df.ID.isin(varIDs_overlapping_repeats)

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

            # define the 'ID' in  a way that resembles VEP
            vcf_df["ID"] = vcf_df["#CHROM"] + "_" + vcf_df.POS.apply(str) + "_" + vcf_df.REF + "/" + vcf_df.ALT

            # get whether the variant overlaps repeats
            vcf_df["overlaps_repeats"] = get_series_variant_in_repeats(vcf_df, repeats_table, replace=replace)

            # check that the ID is unique (check that the drop_duplicates worked)
            if len(set(vcf_df.ID))!=len(vcf_df): 

                #duplicated_ID = vcf_df[vcf_df.duplicate]

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
        if (sum(df_vep.Consequence=="intergenic_variant")/len(df_vep)) > 0.9: print("WARNING: There are >90 percent of the variants that are intergenic. Maybe your gff is wrong. ")
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

        # run the database
        name_database = get_file(genome_dir)

        # setup the folder
        delete_folder(outdir)
        make_folder(outdir)
        shutil.copy2(reference_genome, genome_dir)

        bulding_repModeler_db_std = "%s.genearting_db.std"%genome_dir
        print_if_verbose("getting repeat modeler db. The std is in %s"%bulding_repModeler_db_std)

        # define the cmd
        build_db_cmd = "cd %s && %s -name %s %s > %s 2>&1"%(outdir, repeat_modeller_BuildDatabase, name_database, genome_dir, bulding_repModeler_db_std)
        # run
        run_cmd(build_db_cmd, env=EnvName_RepeatMasker)

        remove_file(bulding_repModeler_db_std)

        # run repeatmodeller
        njobs = int(threads/4) # Specify the number of parallel search jobs to run. RMBlast jobs wil use 4 cores each and ABBlast jobs will use a single core each. i.e. on a machine with 12 cores and running with RMBlast you would use -pa 3 to fully utilize the machine
        print_if_verbose("Running repeat modeller in %s on %i jobs"%(outdir, njobs))

        #raise ValueError("This has to be fixed!!!!")
        cmd = "export PERL5LIB=%s && unset MAFFT_BINARIES && cd %s && %s -database %s -pa %i -LTRStruct -debug"%(repeatmoder_dir, outdir, repeat_modeller, name_database, njobs)

        # add the location were eveything is installed and run
        repeatmodeler_std = "%s/repeatmodeler.std"%outdir
        print_if_verbose("running repeatmodeler. The std is in %s"%repeatmodeler_std)        
        cmd += " -abblast_dir %s -cdhit_dir %s -genometools_dir %s -ltr_retriever_dir %s -mafft_dir %s -ninja_dir %s -recon_dir %s -repeatmasker_dir %s -rmblast_dir %s -rscout_dir %s -trf_prgm %s > %s 2>&1"%(abblast_dir, cdhit_dir, genometools_dir, ltr_retriever_dir, mafft_dir, ninja_dir, recon_dir, repeatmasker_dir, rmblast_dir, rscout_dir, trf_prgm_dir, repeatmodeler_std)

        run_cmd(cmd, env=EnvName_RepeatMasker)

        if file_is_empty(repeat_modeler_outfile): 

            # test that there are no families identified
            no_families_identified = any([l.startswith("No families identified") for l in open(repeatmodeler_std, "r").readlines()[-3:]])

            # test that there are errors
            errors_in_repeatModeler = any(["ERROR" in l.upper() for l in open(repeatmodeler_std, "r").readlines()])

            # test if the LTRpipeline finished
            ltr_pipeline_finished = any([l.startswith("LTRPipeline Time:") for l in open(repeatmodeler_std, "r").readlines()])

            # check if there are no families
            if no_families_identified and (not errors_in_repeatModeler or ltr_pipeline_finished): open(repeat_modeler_outfile, "w").write("no_families_identified")
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
        run_cmd("cd %s && %s -pa %i -dir %s -poly -html -gff %s > %s 2>&1"%(repeat_masker_outdir_default, repeat_masker, threads, repeat_masker_outdir_default, reference_genome, repeat_masker_std_default), env=EnvName_RepeatMasker)
        remove_file(repeat_masker_std_default)

    # run in the personal configuration
    repeat_masker_outfile_personal = "%s/%s.out"%(repeat_masker_outdir_personal, genome_name)
    repeat_masker_std_personal = "%s/%s.std.out"%(repeat_masker_outdir_personal, genome_name)

    if use_repeat_modeller is True and new_families_identified is True:
        
        if file_is_empty(repeat_masker_outfile_personal) or replace is True:
            print_if_verbose("running repeatmasker to get the repeats of the genome with the lib obtained with RepeatModeler. The std is in %s"%repeat_masker_std_personal)
            run_cmd("cd %s && %s -pa %i -dir %s -poly -html -gff -lib %s %s > %s 2>&1"%(repeat_masker_outdir_personal, repeat_masker, threads, repeat_masker_outdir_personal, library_repeats_repeatModeller, reference_genome, repeat_masker_std_personal), env=EnvName_RepeatMasker)
            
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



def get_bgzip_and_and_tabix_vcf_file(file, reference_genome, replace=False):


    """Takes a vcf file and returns a tabixed and gzipped file"""

    # define files
    file_gz = "%s.gz"%file
    file_tmp_gz = "%s.tmp.gz"%file
    file_gz_tbi = "%s.gz.tbi"%file
    file_tmp_gz_tbi = "%s.tmp.gz.tbi"%file

    if file_is_empty(file_gz) or file_is_empty(file_gz_tbi) or replace is True:

        # bgzip
        bgzip_stderr = "%s.generating.stderr"%file_tmp_gz
        print_if_verbose("bgzipping. The stderr is in %s"%bgzip_stderr)
        run_cmd("%s -c %s > %s 2>%s"%(bgzip, file, file_tmp_gz, bgzip_stderr))

        # tabix
        tabix_std = "%s.tabixing.std"%file_tmp_gz
        print_if_verbose("tabix-ing. The std is in %s"%tabix_std)
        run_cmd("%s -p vcf %s > %s 2>&1"%(tabix, file_tmp_gz, tabix_std))

        # remove files
        remove_file(bgzip_stderr)
        remove_file(tabix_std)

        # rename
        os.rename(file_tmp_gz, file_gz)
        os.rename(file_tmp_gz_tbi, file_gz_tbi)

    return file_gz, file_gz_tbi


def get_alternative_genome(reference_genome, vcf, alternative_genome, replace=False, threads=4, only_SNPs=False):

    """This function takes a vcf (with no multiallelics) and generates the alternative_genome"""


    if file_is_empty(alternative_genome) or replace is True:

        # get only SNPs
        if only_SNPs is True:
            print_if_verbose("geting vcf with only SNPs")

            df_vcf, header = get_df_and_header_from_vcf(vcf)

            def get_isSNP(r):

                if "," not in r["ALT"]: return (len(r["REF"])==1 and len(r["ALT"])==1)
                else: return (len(r["REF"])==1 and all([len(x)==1 for x in r["ALT"].split(",") ]))

            # get only SNPs
            df_vcf = df_vcf[df_vcf.apply(get_isSNP, axis=1)]

            # get the lines
            vcf_lines = df_vcf.to_csv(sep="\t", header=True, index=False)  

            # write
            vcf_to_analyze = "%s.onlySNPs.vcf"%vcf 
            vcf_to_analyze_tmp = "%s.tmp"%vcf_to_analyze
            open(vcf_to_analyze_tmp, "w").write("\n".join(header) + "\n" + vcf_lines)
            os.rename(vcf_to_analyze_tmp, vcf_to_analyze)

        else: vcf_to_analyze = vcf
        
        # get the gzipped vcf
        vcf_gz, vcf_tbi = get_bgzip_and_and_tabix_vcf_file(vcf_to_analyze, reference_genome, replace=replace)

        # remove the vcf to analyze if it is not as vcf
        if only_SNPs is True: remove_file(vcf_to_analyze)

        # files
        altgenome_stderr = "%s.generating.stderr"%alternative_genome
        alternative_genome_tmp = "%s.tmp"%alternative_genome

        print_if_verbose("getting alternative genome %s. The stderr is in %s"%(alternative_genome_tmp, altgenome_stderr))
        run_cmd("%s consensus -f %s --haplotype 1 %s > %s 2>%s"%(bcftools, reference_genome, vcf_gz, alternative_genome_tmp, altgenome_stderr))

        remove_file(altgenome_stderr)
        remove_file(vcf_gz)
        remove_file(vcf_tbi)

        os.rename(alternative_genome_tmp, alternative_genome)



def get_correct_gff_and_gff_with_biotype(gff, replace=False):

    """This function takes a gff and returns the correct_gff and gff_with_biotype"""

    correct_gff = "%s_corrected.gff"%(gff); correct_gff_tmp = "%s_corrected_tmp.gff"%(gff)

    if file_is_empty(correct_gff) or replace is True:
        print("correcting gff")
        correct_gff_cmd = "egrep -v '^#' %s > %s"%(gff, correct_gff_tmp); run_cmd(correct_gff_cmd)
        os.rename(correct_gff_tmp, correct_gff)

    # modify gff to add biotype
    gff_with_biotype = "%s_with_biotype.gff"%correct_gff
    if file_is_empty(gff_with_biotype) or replace is True:
        print("adding biotype")

        starting_lines = [line for line in open(correct_gff, "r") if line.startswith("#")]
        df_gff3 = pd.read_csv(correct_gff, skiprows=list(range(len(starting_lines))), sep="\t", names=["chromosome", "source", "type_feature", "start", "end", "score", "strand", "phase", "attributes"])

        def add_biotype(row):
            if "biotype" not in row["attributes"] and "gene_biotype" not in row["attributes"]: row["attributes"] += ";biotype=%s"%row["type_feature"]
            return row["attributes"]

        # add biotype and save
        df_gff3["attributes"] = df_gff3.apply(lambda row: add_biotype(row), axis=1)
        df_gff3.to_csv(gff_with_biotype, sep="\t", index=False, header=False)

    return correct_gff, gff_with_biotype

def get_vcf_df_from_insertion_r(r, gridss_fields):

    """This function takes an insertion row and returns a df with the vcf"""

    # initialize
    df_vcf = pd.DataFrame()

    # define the backbone info string
    backbone_info = ";".join(["%s=%s"%(f, r[f]) for f in gridss_fields])

    # add the insertion site in ChrB as breakpoint
    df_chrB = pd.DataFrame({0 : {"#CHROM":r["ChrB"], "POS":r["StartB"], "ALT":"<BND>"}}).transpose()
    df_chrB["INFO"] = "SVTYPE=BND;%s"%backbone_info
    df_vcf = df_vcf.append(df_chrB)

    # define the chromosome A region
    if r["Copied"] is True: SVTYPE = "DUP"
    else: SVTYPE = "BND"
    df = pd.DataFrame({1 : {"#CHROM":r["ChrA"], "POS":r["StartA"], "ALT":"<%s>"%SVTYPE}}).transpose()
    df["INFO"] = "SVTYPE=%s;END=%i;%s"%(SVTYPE, r["EndA"], backbone_info)
    df_vcf = df_vcf.append(df)

    # add the ID to all of them
    df_vcf["ID"] =  r["IDstring"]

    return df_vcf

def get_vcf_df_from_remaining_r(r, gridss_fields):

    """This is similar to get_vcf_df_from_insertion_r, but with the svtype remaining. You need #CHROM, POS, ALT, INFO, ID"""

    # define the backbone info string
    backbone_info = ";".join(["%s=%s"%(f, r[f]) for f in gridss_fields])

    # interchromosomal breakpoints have 2 rows
    if r["SVTYPE"] in {"ITX1", "ITX2", "INVTX1", "INVTX2", "TAN", "DEL", "INV1", "INV2"}: 

        # get one BND for each breakend
        df1 =  pd.DataFrame({0 : {"#CHROM":r["#CHROM"], "POS":r["POS"], "ALT":"<BND>"}}).transpose()
        df2 =  pd.DataFrame({1 : {"#CHROM":r["CHR2"], "POS":r["END"], "ALT":"<BND>"}}).transpose()

        df_vcf = df1.append(df2)

        # add infoq
        df_vcf["INFO"] = "SVTYPE=BND;%s"%(backbone_info)

    # events with breakpoints
    elif r["SVTYPE"] in {"CVT", "IVT"}:

        # get one BND for each breakend
        df1 =  pd.DataFrame({0 : {"#CHROM":r["#CHROM"], "POS":r["POS"], "ALT":"<BND>"}}).transpose()
        df2 =  pd.DataFrame({1 : {"#CHROM":r["CHR2"], "POS":r["START"], "ALT":"<BND>"}}).transpose()
        df3 =  pd.DataFrame({2 : {"#CHROM":r["CHR2"], "POS":r["END"], "ALT":"<BND>"}}).transpose()

        df_vcf = df1.append(df2).append(df3)

        # add info
        df_vcf["INFO"] = "SVTYPE=BND;%s"%(backbone_info)

    # events with 1 breakend and 1 dup event
    elif r["SVTYPE"]=="CVD":

        # get the duplicated region
        df_dup = pd.DataFrame({0 : {"#CHROM":r["CHR2"], "POS":r["START"], "ALT":"<DUP>", "INFO":"SVTYPE=DUP;END=%i;%s"%(r["END"], backbone_info)}}).transpose()

        # get the breakpoint region
        df_bnd = pd.DataFrame({1 : {"#CHROM":r["#CHROM"], "POS":r["POS"], "ALT":"<BND>", "INFO":"SVTYPE=BND;%s"%(backbone_info)}}).transpose()

        # keep
        df_vcf = df_dup.append(df_bnd)

    else: 
        print(r)
        raise ValueError("%s has not been considered"%r["SVTYPE"])

    # add the ID
    df_vcf["ID"] = r["IDstring"]

    return df_vcf


def get_INFO_vcf_with_breakendMetadata(r, svDF, breakend_fields):

    """Takes a row of a df_vcf and returns the INFO that includes extra breakpoint metadata"""

    # go through each bpID
    if r["ALT"]=="<BND>":

        # get the corresponding event
        bends_metadata_dict = svDF.loc[r["ID"]]["bends_metadata_dict"]

        # calculate the distance between this breakend and the defined ones
        bpID_Ibend_to_distance_to_r = {}
        for bpID, bendDicts in bends_metadata_dict.items():
            for Ibend, bendDict in enumerate(bendDicts):

                # get the distance
                if bendDict["#CHROM"]==r["#CHROM"]: 

                    bpID_Ibend_to_distance_to_r[(bpID, Ibend)] = abs(bendDict["POS"]-r["POS"])

        # best breakend
        closest_bpID_Ibend = [k for k, v in sorted(bpID_Ibend_to_distance_to_r.items(), key=(lambda item: item[1]))][0]

        closest_bpID, Ibend = closest_bpID_Ibend
        bendDict = bends_metadata_dict[closest_bpID][Ibend]
        if bendDict["#CHROM"]!=r["#CHROM"]: raise ValueError("something went wrong with the bend calculation")

        # print the difference
        #print(r["POS"]-bendDict["POS"])

        # get the info with adds
        INFO = "%s;%s"%(r["INFO"], ";".join(["%s=%s"%(f, bendDict[f]) for f in breakend_fields]))

        return  INFO
        
    else: return r["INFO"]

def get_vcf_df_withInsertedSequence_from_svDF(svDF, gridss_fields, breakend_fields):

    """Takes a df with the svDF and retuns it with current_vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", , "INFO"] and gridss_fields in the INFO. SVTYPE is set to BND"""

    # initialize
    dict_all = {}
    Ir = 0

    for ID, r in svDF.iterrows():

        # initialize INFO
        backbone_info = "SVTYPE=insertionBND;" + ";".join(["%s=%s"%(f, r[f]) for f in gridss_fields]) 

        # go through each BPdict
        for bpID, bendDict_list in r["bends_metadata_dict"].items():

            # go through each breakend
            for bendDict in bendDict_list:

                # define the info
                INFO = "%s;%s"%(backbone_info, ";".join(["%s=%s"%(f, bendDict[f]) for f in breakend_fields]))

                # redefine ALT
                REF = bendDict["REF"]
                ALT = bendDict["inserted_sequence"]
                if len(ALT)==0 and len(REF)==0: continue
                elif ALT==REF: continue
                elif len(ALT)==0 and len(REF)>0: ALT = "-"

                # define fields of each insertion
                dict_all[Ir] = {"#CHROM":bendDict["#CHROM"],
                                "POS":bendDict["POS"],
                                "ID":ID,
                                "REF":REF,
                                "ALT":ALT,
                                "INFO":INFO}

                Ir += 1

    # get as df
    df_vcf_insertions = pd.DataFrame(dict_all).transpose()
    interesting_fields = list(df_vcf_insertions.keys())

    # add all the INFO fields
    df_vcf_insertions = get_vcf_df_with_INFO_as_single_fields(df_vcf_insertions).sort_values(by=["ID", "#CHROM", "POS", "INFO_QUAL", "INFO_real_AF"], ascending=False)

    # discard duplicates
    df_vcf_insertions = df_vcf_insertions.drop_duplicates(subset=["ID", "#CHROM", "POS"], keep="first")

    return df_vcf_insertions[interesting_fields]

def get_vcf_df_for_svDF(svDF, svtype, reference_genome, df_gridss):

    """This function takes an SVdf and rewrites it as a vcf, so that it can be loaded to clove"""

    # add to svDF
    svDF = cp.deepcopy(svDF)
    svDF["BPS_TYPE"] = "GRIDSS-CLOVE"

    gridss_fields = ['real_AF_min', 'real_AF_max', 'real_AF_mean', 'bpIDs', 'QUAL_min', 'QUAL_max', 'QUAL_mean', 'all_FILTERs', 'BPS_TYPE', 'BREAKPOINTIDs']

    # define vcf fields
    vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    # if the svDF is empty, return an empty df
    if len(svDF)==0: return pd.DataFrame(columns=vcf_fields)

    # set the floats to only 3 decimals
    for f in ['real_AF_min', 'real_AF_max', 'real_AF_mean', 'QUAL_min', 'QUAL_max', 'QUAL_mean']: svDF[f] = svDF[f].apply(lambda x: "%.3f"%x)

    # initialize thd df_vcf
    df_vcf = pd.DataFrame()

    # get the simple SVs
    svtype_to_SVTYPE = {"tandemDuplications":"TDUP", "deletions":"DEL", "inversions":"BND"}
    if svtype in {"tandemDuplications", "deletions", "inversions"}:

        # add obvious fields
        df_vcf["#CHROM"] = svDF["Chr"]
        df_vcf["POS"] = svDF["Start"]

        # add the SVTYPE and END
        df_vcf["ALT"] = "<%s>"%(svtype_to_SVTYPE[svtype])
        svDF["SVTYPE"] = svtype_to_SVTYPE[svtype]
        svDF["END"] = svDF.End

        # define the ID
        df_vcf["ID"] = svDF.IDstring 

        # add the info
        info_fields = ["SVTYPE", "END"] + gridss_fields
        df_vcf["INFO"] = svDF.apply(lambda r: ";".join(["%s=%s"%(f, r[f]) for f in info_fields]), axis=1)

    elif svtype=="translocations": 

        # they are all balanced. We record chrA and chrB breakpoints

        # add the chrB bp pos
        chr_to_len = get_chr_to_len(reference_genome)
        svDF["ChrB_bp_pos"] = svDF.apply(lambda r: get_ChrB_bp_pos_translocations(r, chr_to_len, first_bp_pos=0), axis=1)

        # get both dfs
        df_A = svDF[["ChrA", "EndA", "IDstring"] + gridss_fields].rename(columns={"ChrA":"#CHROM", "EndA":"POS"})
        df_B = svDF[["ChrB", "ChrB_bp_pos", "IDstring"] + gridss_fields].rename(columns={"ChrB":"#CHROM", "ChrB_bp_pos":"POS"})
        df_vcf = df_A.append(df_B)

        # add direct things
        df_vcf["ALT"] = "<BND>"
        df_vcf["ID"] = df_vcf.IDstring
        df_vcf["SVTYPE"] = "BND"

        # add the info
        info_fields = ["SVTYPE"] + gridss_fields
        df_vcf["INFO"] = df_vcf.apply(lambda r: ";".join(["%s=%s"%(f, r[f]) for f in info_fields]), axis=1)

    elif svtype=="insertions": df_vcf = pd.concat([get_vcf_df_from_insertion_r(r, gridss_fields) for I, r in svDF.iterrows()])

    elif svtype=="remaining": df_vcf = pd.concat([get_vcf_df_from_remaining_r(r, gridss_fields) for I, r in svDF.iterrows()]) 

    else: raise ValueError("%s is not valid"%svtype) 

    # add the REF 
    df_vcf["REF"] = "."

    # deifine the breakend_fields to add to the vcf
    breakend_fields = ['allele_frequency', 'allele_frequency_SmallEvent', 'real_AF', 'FILTER', 'has_poly16GC', 'length_inexactHomology', 'length_microHomology', 'QUAL', 'overlaps_repeats', 'BREAKPOINTID']

    # add the inserted sequences
    svDF = svDF.set_index("IDstring", drop=False)
    current_vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "INFO"]
    df_vcf = df_vcf[current_vcf_fields].append(get_vcf_df_withInsertedSequence_from_svDF(svDF, gridss_fields, breakend_fields)[current_vcf_fields])

    # add breakpoint metadata to breakend-like elements
    df_vcf["INFO"] = df_vcf.apply(lambda r: get_INFO_vcf_with_breakendMetadata(r, svDF, breakend_fields), axis=1)

    # add general thingd
    
    df_vcf["FORMAT"] = "."
    df_vcf["QUAL"] = "."
    df_vcf["FILTER"] = "."

    return df_vcf[vcf_fields]

def get_corrected_Consequence_for_vep_r(r):

    """Takes a consequence (comma sepparated) and returns it with BND added if it is a BND breakpoint"""

    if r["Allele"]=="BND": return ",".join(["%s_BND"%x for x in r["Consequence"].split(",")])
    else: return r["Consequence"]


def get_df_subwindows_window_r(r, n_subwindows=20):

    """ This function takes a row of a df_windows and returns a df with the n_suwindows within this one"""

    # define the length of each suwindow
    len_subw = int((r["end"] - r["start"])/n_subwindows)

    # go through each subwindow
    data_dict = {}
    for Iw, start in enumerate(np.linspace(r["start"], r["end"]-len_subw, n_subwindows)):

        # define the end
        end = min([r["end"], start+len_subw])

        # keep
        data_dict[Iw] = {"start":int(start), "end":int(end)}

    # get as df
    df = pd.DataFrame(data_dict).transpose()
    df["chromosome"] = r["chromosome"]
    df["ID"] = r["ID"]

    return df


def get_confidence_interval_bootsrap(array, central_val_estimator, confidence=0.95, n_boots=100):

    """Calculates a ci of confidence for a nd array, using nboots"""

    # get an idea of all estimates
    n_size = len(array) - 1
    all_estimates = [central_val_estimator(resample(array, n_samples=n_size)) for i in range(n_boots)]

    # get CI
    tail_p = ((1-confidence)/2)*100
    return (np.percentile(all_estimates, tail_p), np.percentile(all_estimates, 100-tail_p))


def get_pvalue_is_normal_distribution(array):

    """This function returns the pvalue of the array being normally dsitributed"""

    return stats.normaltest(array, nan_policy="raise")[1]

def get_summary_statistics_series_df_subwindows(df):

    """This function takes a df of subwindows and returns a series with the summary statistics of them"""

    # check that it is a df
    if type(df)!=pd.core.frame.DataFrame: raise ValueError("df is not a df")

    # get the data, the relative coverage
    relative_coverage = df.relative_coverage
    position = df.start

    # get CIs
    mean_CI = get_confidence_interval_bootsrap(relative_coverage, np.mean)
    median_CI = get_confidence_interval_bootsrap(relative_coverage, np.median)

    # get correlations, depending on if there is something to correlate
    if len(set(position))==1 or len(set(relative_coverage))==1:

        spearman_r = 0
        spearman_p = 1
        pearson_r = 0
        pearson_p = 1

    else:
        
        spearman_r, spearman_p = stats.spearmanr(position, relative_coverage, nan_policy="raise")
        pearson_r, pearson_p = stats.pearsonr(position, relative_coverage)

        spearman_r = abs(spearman_r)
        pearson_r = abs(pearson_r)

    # debug
    if pd.isna(pearson_r) or pd.isna(spearman_r): 

        print(pearson_r, spearman_r)
        print(df[["start", "relative_coverage"]])
        raise ValueError("correlations should not be NaN")

    # return the series
    data_dict = {"mean_rel_coverage" : np.mean(relative_coverage),
                 "median_rel_coverage" : np.median(relative_coverage),
                 "mean95CI_lower_rel_coverage" : mean_CI[0],
                 "mean95CI_higher_rel_coverage" : mean_CI[1],
                 "median95CI_lower_rel_coverage" : median_CI[0],
                 "median95CI_higher_rel_coverage" : median_CI[1],
                 "abs_spearman_r" : spearman_r,
                 "abs_pearson_r" : pearson_r,
                 "spearman_p" : spearman_p,
                 "pearson_p" : pearson_p}
    # "pvalNormDist_rel_coverage" : get_pvalue_is_normal_distribution(data)

    return pd.Series(data_dict)


def get_df_subwindows_from_df_windows(df_windows, n_subwindows=20):

    """Takes a df_windows and returns an equivalent one with subwindows of n_subwindows"""

    # get df
    rows_iterator = [r for I,r in df_windows.iterrows()]
    df_subwindows = pd.concat(map(lambda r: get_df_subwindows_window_r(r, n_subwindows=n_subwindows), rows_iterator))

    # get only those where the start!=end
    df_subwindows = df_subwindows[df_subwindows.start!=df_subwindows.end]

    # check
    if any(df_subwindows.end<=df_subwindows.start): 
        print(df_subwindows)
        print(df_subwindows[df_subwindows.end<=df_subwindows.start])
        raise ValueError("df_subwindows: start should be < end")

    return df_subwindows


def get_coverage_df_windows_with_within_windows_statistics(df_windows, outdir, sorted_bam, reference_genome, median_coverage, replace=False, threads=4):

    """This function takes a df with windows and returns a similar df with the coverage of n_subwindows calculated """

    # check that the ID is unique
    initial_IDs = set(df_windows.ID)
    if len(df_windows)!=len(initial_IDs): raise ValueError("The ID should be unique")

    # if the df is empty, return it
    if len(df_windows)==0: return df_windows

    # build a df with windows, where each window comes from dividing 
    df_subwindows = get_df_subwindows_from_df_windows(df_windows, n_subwindows=20)
    initial_len_subwindows = cp.deepcopy(len(df_subwindows))

    # debug
    if any(df_subwindows.start>=df_subwindows.end): raise ValueError("there can't be any windows with start>=end")

    # get the  coverage of the subwindows
    windows_file = "%s/subwindows.bed"%outdir
    df_subwindows[["chromosome", "start", "end"]].to_csv(windows_file, sep="\t", header=True, index=False)
    df_subwindows_coverage = get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, windows_file, replace=replace, run_in_parallel=True, delete_bams=False, threads=threads).drop_duplicates(subset=["chromosome", "start", "end"])

    print_if_verbose("merging with subwindows")
    df_subwindows = df_subwindows.merge(df_subwindows_coverage, on=["chromosome", "start", "end"], how="left", validate="many_to_one")
    if len(df_subwindows)!=initial_len_subwindows or any(pd.isna(df_subwindows.mediancov_1)): 
        print(df_subwindows)
        print(initial_len_subwindows)
        print(df_subwindows.mediancov_1)
        print(any(pd.isna(df_subwindows.mediancov_1)))
        raise ValueError("something went wrong with the subwindows coverage calculation")

    # add the relative coverage
    print_if_verbose("calculating relative_coverage")
    df_subwindows["relative_coverage"]  = df_subwindows.mediancov_1 / median_coverage

    # get the final df with summary statistics 
    print_if_verbose("getting summary stats in parallel on %i threads"%threads) 
    df_subwindows = df_subwindows.set_index("ID", drop=False)
    all_IDs = set(df_subwindows.ID)
    all_dfs_tuples = [(df_subwindows.loc[ID],) for ID in all_IDs]

    # run in parallel
    with multiproc.Pool(threads) as pool:
        list_summary_statistics_series = pool.starmap(get_summary_statistics_series_df_subwindows, all_dfs_tuples) 
            
        pool.close()
        pool.terminate()

    # get a df
    df_final = pd.DataFrame(dict(zip(all_IDs, list_summary_statistics_series))).transpose()

    # slow way
    #df_final = df_subwindows.groupby("ID").apply(get_summary_statistics_series_df_subwindows)

    df_final["ID"] = df_final.index
    df_final = df_final.merge(df_windows.set_index("ID", drop=True), left_index=True, right_index=True, validate="one_to_one") 

    # check that the IDs are correct
    missing_IDs_fromInitial = initial_IDs.difference(set(df_final.ID))
    if initial_IDs!=set(df_final.ID): raise ValueError("something went wrong from the IDs. These are missing IDs from the initial: %s"%missing_IDs_fromInitial)

    #### ADD FIELDS ####

    # get chrom_to_len
    chrom_to_len = get_chr_to_len(reference_genome)

    # add fields
    df_final["chr_len"] = df_final.chromosome.map(chrom_to_len)
    if any(pd.isna(df_final.chr_len)): raise ValueError("error in len calculation")
    df_final["fraction_chromosome_covered"] = (df_final.end - df_final.start) / df_final.chr_len

    # add the coverage to which to compare
    df_final["rel_coverage_to_compare"] = df_final.apply(get_rel_coverage_to_compare_region_r , axis=1)

    # add the relative fields
    df_final["mean95CI_higher_rel_coverage_relative"] = df_final.mean95CI_higher_rel_coverage / df_final.rel_coverage_to_compare
    df_final["mean95CI_lower_rel_coverage_relative"] = df_final.mean95CI_lower_rel_coverage / df_final.rel_coverage_to_compare
    df_final["median95CI_higher_rel_coverage_relative"] = df_final.median95CI_higher_rel_coverage / df_final.rel_coverage_to_compare
    df_final["median95CI_lower_rel_coverage_relative"] = df_final.median95CI_lower_rel_coverage / df_final.rel_coverage_to_compare

    ####################

    return df_final

def get_rel_coverage_to_compare_region_r(r):

    """This function returns a coverage value for a region row. This should be the closest coverage of either 5' / 3' regios or 1, in case that the variant covers the whole chromosome"""

    if r["fraction_chromosome_covered"]>=0.9: return 1.0
    else: return find_nearest([r["relative_coverage_3"], r["relative_coverage_5"]], r["relative_coverage_target"])

def clean_chromosomal_bam_files(sorted_bam, reference_genome):

    """Takes a sorted bam and removes the chromosomal bam files"""

    # get the chrom to len
    chrom_to_len = get_chr_to_len(reference_genome)

    # clean each
    for chrom in chrom_to_len: 
        sorted_bam_chr = "%s.%s.bam"%(sorted_bam, chrom)
        remove_file(sorted_bam_chr); remove_file("%s.bai"%sorted_bam_chr)

def check_that_df_index_is_unique(df):

    """checks that df index is unique"""

    if len(set(df.index))!=len(df): raise ValueError("df index should be unique")



def get_list_clusters_from_dict(key_to_setVals):

    """Takes a dictionary mapping strings to strings and returns a list of sets, each of them having uniquely clusters of connected data"""

    # initialize a list_clusters
    list_clusters = []

    # initialize a list of mapped regions
    list_partial_clusters = [cp.deepcopy({key}.union(setVals)) for key, setVals in key_to_setVals.items()]

    # go through each element in dict
    for all_items in list_partial_clusters:

        # if there is any overlap with any of the clusters, add it there
        cluster_found = False
        for cluster in list_clusters:
            if len(cluster.intersection(all_items))>0: 
                cluster.update(all_items)
                cluster_found = True
                break

        # else initialize with all elements of this cluster
        if cluster_found is False: list_clusters.append(all_items)

    # make sure that all clusters are nonoverlapping with each other
    for I in range(len(list_clusters)):
        for J in range(I+1, len(list_clusters)):

            if I==J: raise ValueError("We don't want to compare the same I and J")

            # define clusters
            clusterI = list_clusters[I]
            clusterJ = list_clusters[J]

            if len(clusterI.intersection(clusterJ))!=0: 
                pass
                #print(I, clusterI, "\n", J, clusterJ)
                #raise ValueError("There are some overlapping clusters")
                #print("There are some overlapping clusters")


    return list_clusters

def get_bestID_from_df_CNV_cluster(clustered_nuericIDs, df_CNV):

    """This function takes a df_CNV and the clusteredIDs. It returns the best clusterID"""

    # get the df with the lcuster IDs
    df = df_CNV.loc[clustered_nuericIDs].sort_values(by=["type_CNVcall_int", "length"], ascending=False)

    return df.ID.iloc[0]


def get_list_clusters_overlapping_df_CNV(outdir, df_CNV, pct_overlap, threads):

    """Takes a df_CNV with chromosome, start, and end, where the ID is the index. It returns a list of sets, each set containing IDs of CNVs that overlap by >=pct_overlap and are from the same type.

    It will run bedtools intersect between the df_CNV against itself."""

    print_if_verbose("getting list_clusters_overlapping_df_CNV")

    # checks
    if len(set(df_CNV.SVTYPE).difference({"DUP", "DEL"}))>0: raise ValueError("SVTYPE is not properly formated")
    if len(set(df_CNV.type_CNVcall).difference({"gridssClove", "coverage"}))>0: raise ValueError("type_CNVcall is not properly formated")

    # get the ID_to_overlappingIDs

    # add the combination of chromosome and SVtype. This will be used in the bedtools intersect to get the correct IDs
    df_CNV["chromosome_SVTYPE"] = df_CNV.chromosome + "_" + df_CNV.SVTYPE

    # write a bed with the intersecting IDs
    bed_cnv_regions = "%s/df_CNV.bed"%outdir
    df_CNV[["chromosome_SVTYPE", "start", "end", "numericID"]].to_csv(bed_cnv_regions, sep="\t", header=False, index=False)

    # run bedmap to get a file where each line corresponds to the regions to which each 
    bedmap_outfile = "%s/overlapping_IDs.txt"%outdir
    bedmap_stderr = "%s.stderr"%bedmap_outfile

    print_if_verbose("running bedmap. The stderr is in %s"%bedmap_stderr)
    run_cmd("%s --fraction-both %.2f --echo-map-id  --delim '\t' %s > %s 2>%s"%(bedmap, pct_overlap, bed_cnv_regions, bedmap_outfile, bedmap_stderr))

    remove_file(bedmap_stderr)

    # load as df, which already has the same order as df_CNV
    df_bedmap = pd.read_csv(bedmap_outfile, sep="\t", header=None, names=["overlapping_IDs"])
    df_bedmap["overlapping_IDs_set"] = df_bedmap.overlapping_IDs.apply(str).apply(lambda x: {int(y) for y in x.split(";")})

    # define
    ID_to_overlappingIDs = dict(df_bedmap["overlapping_IDs_set"])

    # get the list of clusters
    print_if_verbose("getting lists of clusters")
    list_clusters = get_list_clusters_from_dict(ID_to_overlappingIDs)

    # check
    all_IDs = set(df_CNV.numericID)
    all_IDs_in_cluster = set.union(*list_clusters)
    if all_IDs!=all_IDs_in_cluster: raise ValueError("all IDs should be in clusters")

    print_if_verbose("list_clusters_overlapping_df_CNV already ran. There are %i clusters"%len(list_clusters))

    return list_clusters


def get_nonRedundant_CNVcalls_coverage(outdir, df_CNV, df_vcf_forCNV, threads, replace, pct_overlap=0.8):

    """Gets a df_CNV with no redudnant calls (those that overlap by more than 80% with other rows in df_CNV or df_vcf_forCNV)"""

    if len(df_CNV)==0: return df_CNV

    print_if_verbose("getting non-redundant CNV calls")

    # get the index
    initial_index = set(df_CNV.index)
    if len(initial_index)!=len(df_CNV): raise ValueError("index should be unique")

    # get the initial fields
    initial_fields = list(df_CNV.columns)

    # keep
    df_CNV = cp.deepcopy(df_CNV)
    df_vcf_forCNV = cp.deepcopy(df_vcf_forCNV)

    # define all called SVs
    fields = ["ID", "chromosome", "start", "end", "SVTYPE", "type_CNVcall"]
    all_df_CNV = df_CNV[fields].append(df_vcf_forCNV[fields]).sort_values(by=["SVTYPE", "chromosome", "start", "end"])

    # sort as bed
    all_df_CNV = all_df_CNV.sort_values(by=["SVTYPE", "chromosome", "start", "end"])

    # add the ID and the numeric ID (the numeric ID takes less space in the bedtools intersect). This will be the index
    all_df_CNV["numericID"] = list(range(0, len(all_df_CNV)))
    all_df_CNV = all_df_CNV.set_index("numericID", drop=False)

    # make sure that the ID is unique
    check_that_df_index_is_unique(all_df_CNV)

    # define the initial length
    initial_len_all_df_CNV = len(all_df_CNV)

    # define the clusters of CNVs that are overlapping by >=80% of their extension
    list_clusters = get_list_clusters_overlapping_df_CNV(outdir, all_df_CNV, pct_overlap, threads)

    # add fields for sorting of redundant variants according to their type and the quality mesurements
    type_CNVcall_to_int = {"gridssClove":1, "coverage":0}
    all_df_CNV["type_CNVcall_int"] = all_df_CNV.type_CNVcall.apply(lambda x: type_CNVcall_to_int[x])

    all_df_CNV["length"] = all_df_CNV.end - all_df_CNV.start

    # get the best IDs from each cluster
    best_NR_IDs = set(map( (lambda x: get_bestID_from_df_CNV_cluster(x, all_df_CNV) ), list_clusters))

    # get the df with these IDs
    df_CNV_NR = df_CNV[df_CNV.ID.isin(best_NR_IDs)]
    
    # at the end set the quality to a '.'
    df_CNV_NR["QUAL"] = "."

    return df_CNV_NR[initial_fields]

def get_df_gridss_chrom_to_positions(df_gridss, chrom_to_len):

    """Gets a df_gridss and returns a dict mapping each chromosome to the breakpoints"""

    # intit breakopints positions with all chroms
    chrom_to_bpPositions = {c:set() for c in set(chrom_to_len)}

    # debug 
    if len(df_gridss)==0: return chrom_to_bpPositions

    # get the bpPOs
    df_gridss["POS_0based"] = df_gridss.POS - 1
    chrom_to_bpPositions_df_gridss = dict(df_gridss.groupby("#CHROM").apply(lambda df_c: set(df_c.POS_0based)))

    for chrom, bpPositions_df_gridss in chrom_to_bpPositions_df_gridss.items(): chrom_to_bpPositions[chrom].update(bpPositions_df_gridss)


    return chrom_to_bpPositions


def get_df_windows_arround_pos_r(r, chrom_to_len, min_sv_size):

    """This function returns a df with chromosome, start, end for the provided r. It checks to not be outside of the boundaries """

    # init dict
    data_dict = {}

    # go through each region
    for region in ["5", "3"]:

        # define the start and end
        if region=="5":

            start = max([r["POS"]-min_sv_size, 0])
            end = r["POS"]

        elif region=="3":

            start = r["POS"]
            end = min([r["POS"]+min_sv_size, chrom_to_len[r["#CHROM"]]])


        if start>=end: raise ValueError("start can't be after end")

        data_dict[region] = {"start":start, "end":end, "chromosome":r["#CHROM"], "ID":r["ID"], "region": "%s_region"%region}


    return pd.DataFrame(data_dict).transpose()

def get_df_gridss_low_confidence_reasonable_breakpoints(df_gridss, ploidy, min_coverage_duplication, max_coverage_deletion, min_sv_size, chrom_to_len, reference_genome, sorted_bam, outdir, replace, threads, median_coverage):

    """This function gets the df_gridss and returns the same one with the breakpoints that are low confidence and may include potential breakpoints."""

    print_if_verbose("getting gridss df low confidence breakpoints")


    # define an AF threshold
    min_AF_gridss = (1/ploidy)*0.1

    # deinfe the initial fields
    initial_fields = list(df_gridss.columns)
           
    # keep low confidence vars 
    df_gridss = df_gridss[(df_gridss.type_BEND=="lowConfidence") & ((df_gridss.allele_frequency_SmallEvent>=min_AF_gridss) | (df_gridss.allele_frequency>=min_AF_gridss))]

    if len(df_gridss)==0: return pd.DataFrame(columns=initial_fields)

    # consider only the df gridss where the POS is different than the ends of the chromosomes
    df_gridss["chromosome_end"] = df_gridss["#CHROM"].apply(lambda c: chrom_to_len[c])
    df_gridss = df_gridss[(df_gridss.POS!=0) & (df_gridss.POS!=1) & (df_gridss.POS!=df_gridss.chromosome_end)]

    # get windows of +-min_sv_size
    df_windows = pd.concat([df for df in df_gridss[["#CHROM", "POS", "ID"]].apply(lambda r: get_df_windows_arround_pos_r(r, chrom_to_len, min_sv_size), axis=1)])

    # get the coverage of these windows
    windows_file = "%s/df_gridss_lowConf_surrounding.bed"%outdir
    df_windows[["chromosome", "start", "end"]].to_csv(windows_file, sep="\t", index=False, header=True)
    df_coverage = get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, windows_file, replace=replace, run_in_parallel=True, delete_bams=False, threads=threads)

    # merge    
    df_coverage["IDwindow"] = df_coverage.chromosome + "_" + df_coverage.start.apply(str) + "_" + df_coverage.end.apply(str)
    df_coverage = df_coverage.set_index("IDwindow")

    df_windows["IDwindow"] = df_windows.chromosome + "_" + df_windows.start.apply(str) + "_" + df_windows.end.apply(str)

    # add the coverage
    df_windows["coverage"] = df_windows.IDwindow.apply(lambda x: df_coverage.loc[x, "mediancov_1"])
    df_windows["relative_coverage"] = df_windows.coverage / median_coverage

    # add the coverage regions
    df_windows = df_windows.set_index(["ID", "region"])

    for region in ["5_region", "3_region"]: df_gridss["relative_coverage_%s"%region] = df_gridss.ID.apply(lambda x: df_windows.loc[(x, region), "relative_coverage"])

    # get the putatively interesting regions
    df_gridss["region5_duplicated"] = df_gridss.relative_coverage_5_region>=min_coverage_duplication
    df_gridss["region3_duplicated"] = df_gridss.relative_coverage_3_region>=min_coverage_duplication
    df_gridss["region5_deleted"] = df_gridss.relative_coverage_5_region<=max_coverage_deletion
    df_gridss["region3_deleted"] = df_gridss.relative_coverage_3_region<=max_coverage_deletion

    df_gridss["regions_5and3_are_different"] = (df_gridss.region5_duplicated!=df_gridss.region3_duplicated) | (df_gridss.region5_deleted!=df_gridss.region3_deleted)

    # get the regoins where there is a change in deletion or duplication status
    df_gridss = df_gridss[df_gridss.regions_5and3_are_different]

    if len(df_gridss)==0: return pd.DataFrame(columns=initial_fields)
    else: return df_gridss[initial_fields]


def get_setIDs_query_df_NR(IDs):

    """Takes an IDs set and returns it as set"""

    return {int(x) for x in IDs.split(";")}

def get_bedmap_df_per_chromosome_NR_query_regions(input_bed, chromosome, replace):

    """This function takes an input bed and runs bedmap against itself. It does it only per chromosome. It returns a df with the important features"""

    # define the outfile
    outfile = "%s.bedmap_run.%s.bed"%(input_bed, chromosome)

    if file_is_empty(outfile) or replace is True:

        outfile_tmp = "%s.tmp"%outfile
        stderr = "%s.generating.stderr"%outfile

        # define the 2 cmds
        cmd_normal = "%s --delim '\t' --fraction-both 0.95 --echo-map-range --echo-map-id --chrom %s %s | sort -u -t '\t' -k4,4 > %s 2>%s"%(bedmap, chromosome, input_bed, outfile_tmp, stderr)

        cmd_low_memory = cmd_normal.replace("--echo-map-range", "--min-memory --echo-map-range")

        # try with more or less memory consumption
        try: run_cmd(cmd_normal)
        except: run_cmd(cmd_low_memory)

        remove_file(stderr)

        os.rename(outfile_tmp, outfile)

    # load as df
    df_NR_bed = pd.read_csv(outfile, sep="\t", header=None, names=["chromosome", "start", "end", "IDs"])

    # get the IDs as set
    df_NR_bed["IDs_set"] = df_NR_bed.IDs.apply(get_setIDs_query_df_NR)

    return df_NR_bed


def get_non_redundant_regions_df_from_query_df(query_df, outdir, threads, replace):

    """Takes a query_df and returns the df where there are only unique regions"""

    print_if_verbose("running get_non_redundant_regions_df_from_query_df")

    # write a sorted bed into outdir
    query_df = query_df.sort_values(by=["chromosome", "start", "end"])
    regions_bed = "%s/query_df_regions.bed"%outdir
    query_df[["chromosome", "start", "end", "numericID"]].to_csv(regions_bed, sep="\t", header=False, index=False)

    # define the inputs of get_bedmap_output_per_chromosome
    inputs_fn = [(regions_bed, chromosome, replace) for chromosome in set(query_df.chromosome)]

    # run in parallel
    """ 
    with multiproc.Pool(threads) as pool:
        list_NR_dfs = pool.starmap(get_bedmap_df_per_chromosome_NR_query_regions, inputs_fn) 
            
        pool.close()
        pool.terminate()

    """

    # not in parallel
    list_NR_dfs = map(lambda x: get_bedmap_df_per_chromosome_NR_query_regions(x[0], x[1], x[2]), inputs_fn)

    # get the concatenated 
    df_NR_regions = pd.concat(list_NR_dfs)

    # define the missing IDs
    missing_IDs = set(query_df.numericID).difference(set.union(*df_NR_regions.IDs_set))
    if len(missing_IDs)>0: raise ValueError("%s are missing in df_NR_regions"%missing_IDs)

    # add a numeric ID
    df_NR_regions["numericID"] = list(range(0, len(df_NR_regions)))

    return df_NR_regions[["chromosome", "start", "end", "IDs_set", "numericID"]]

def get_potentially_CNV_query_df(query_df, sorted_bam, reference_genome, outdir, threads, replace, max_coverage_deletion, min_coverage_duplication, mitochondrial_chromosome, median_coverage):

    """This function takes a query_df (query, start, end). It will filter out regions with no possibility of CNV. The idea is to downsample the bam to 5x and keep those regions 

        -  """

    print_if_verbose("running get_potentially_CNV_query_df")

    # define the initial fields
    initial_fields = list(query_df.keys())

    # add the ID
    query_df["numericID"] = list(range(0, len(query_df)))

    ######### DOWNSAMPLE BAM TO MAKE QUICK CALCULATIONS #########

    # downasmple the bam so that the coverage is keep to 5
    fraction_bam = 5/median_coverage
    if fraction_bam<1: 

        # downsample
        downsampled_bam = "%s/sorted_bam_downsampled_%.3f.bam"%(outdir, fraction_bam)

        if file_is_empty(downsampled_bam) or replace is True:

            downsampled_bam_tmp = "%s.tmp"%downsampled_bam
            downsample_bamfile_keeping_pairs(sorted_bam, fraction_reads=fraction_bam, threads=threads, name="sampleX", sampled_bamfile=downsampled_bam_tmp)

            # index
            index_bam(downsampled_bam_tmp, threads=threads)

            # renames
            os.rename(downsampled_bam_tmp+".bai", downsampled_bam+".bai")
            os.rename(downsampled_bam_tmp, downsampled_bam)
    else: 

        downsampled_bam = sorted_bam
        print("WARNING: The median coverage of your bam is %.3f"%median_coverage)

    # calculate the median_coverage for the downsampled bam
    destination_dir = "%s.calculating_windowcoverage"%downsampled_bam
    downsampled_coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, downsampled_bam, windows_file="none", replace=replace, run_in_parallel=True, delete_bams=True), sep="\t")

    downsampled_median_coverage = get_median_coverage(downsampled_coverage_df, mitochondrial_chromosome)

    #############################################################

    # get a df with the unique regions that overlap 0.95 between them
    df_NR_regions = get_non_redundant_regions_df_from_query_df(query_df, outdir, threads, replace)

    ######### ADD COVERAGE TO THE df_NR_regions #########

    # write as a bed
    NR_regions_file = "%s/NR_query_df_regions.bed"%outdir
    df_NR_regions[["chromosome", "start", "end"]].to_csv(NR_regions_file, sep="\t", index=False, header=True)

    # get df coverage
    df_coverage = get_coverage_per_window_df_without_repeating(reference_genome, downsampled_bam, NR_regions_file, replace=replace, run_in_parallel=True, delete_bams=True, threads=threads)[["chromosome", "start", "end", "mediancov_1"]]

    # merge 
    initial_len_df_NR_regions = len(df_NR_regions)
    df_NR_regions = df_NR_regions.merge(df_coverage, on=["chromosome", "start", "end"], validate="many_to_one", how="left")
    if len(df_NR_regions)!=initial_len_df_NR_regions or any(pd.isna(df_NR_regions.mediancov_1)): raise ValueError("the merge did not work")

    # add the relative coverage
    df_NR_regions["relative_coverage"] = df_NR_regions.mediancov_1/downsampled_median_coverage

    #####################################################


    ###### DEFINE THE INTERESTING REGIONS ######

    # define the CNV IDs
    thshd_max_coverage_deletion = max_coverage_deletion*1.5
    thshd_min_coverage_duplication = min_coverage_duplication*0.9

    # define the series of set IDs that might be under CNV
    CNV_IDs_set = df_NR_regions[(df_NR_regions.relative_coverage>=thshd_min_coverage_duplication) | (df_NR_regions.relative_coverage<=thshd_max_coverage_deletion)]["IDs_set"]

    if len(CNV_IDs_set)>0: CNV_IDs = set.union(*CNV_IDs_set)
    else: CNV_IDs = set()

    query_df = query_df[query_df.numericID.isin(CNV_IDs)]

    ############################################

    return query_df[initial_fields]

def get_bestID_from_df_gridss_lowConf_cluster(clustered_numericIDs, df_lowConf):

    """This function takes a df_gridss_lowConf and the clusteredIDs. It returns the best clusterID"""

    # get the df with the lcuster IDs
    df = df_lowConf.loc[clustered_numericIDs].sort_values(by=["QUAL", "real_AF"], ascending=False)

    return df.numericID.iloc[0]


def get_correct_POS_in1based(r):

    """Takes a row of df_vcf SV calling and returns the 1- based position. Only the ACTG ones will not be considered"""

    if r["ALT"] in {"<TDUP>", "<DUP>", "<DEL>", "<BND>"}: return (r["POS"]+1)
    elif r["INFO"].startswith("SVTYPE=insertionBND"): return r["POS"]
    else: 
        print(r)
        raise ValueError("r is not properly formatted")


def get_correctID_and_INFO_df_vcf_SV_CNV(r):

    """Takes a row of df_vcf and returns the ID and the INFO"""

    # new info
    newINFO = "%s;variantID=%s"%(r["INFO"], r["ID"])

    # new ID depends on the type of variant
    if r["ALT"] in {"<TDUP>", "<DUP>", "<DEL>"}: extraID = "CNV"
    elif r["ALT"]=="<BND>": extraID = "BND-%s-%i"%(r["#CHROM"], r["POS"])
    elif r["INFO"].startswith("SVTYPE=insertionBND"): extraID = "insertion-%s-%i"%(r["#CHROM"], r["POS"])
    else: raise ValueError("r is not properly formatted")

    type_var = r["ID"].split("|")[0]
    other_ID = "|".join(r["ID"].split("|")[1:])
    newID = "%s|%s|%s"%(type_var, extraID, other_ID)

    return pd.Series({"ID" : newID, "INFO" : newINFO})


def get_correct_INFO_withEND_in1based(r, chr_to_len):

    """Adds +1 to INFO_END"""

    # get the INFO_ENDstring
    INFO_ENDstring = [x for x in r["INFO"].split(";") if x.startswith("END=")]
    if len(INFO_ENDstring)==0: return r["INFO"]
    elif len(INFO_ENDstring)==1:

        # get the string
        INFO_ENDstring = INFO_ENDstring[0]

        # get the string with end summed with 1
        end = int(INFO_ENDstring.split("=")[1])
        if end!=chr_to_len[r["#CHROM"]]: end+=1
        new_INFO_ENDstring = "END=%i"%end

        return r["INFO"].replace(INFO_ENDstring, new_INFO_ENDstring)

    else: raise ValueError("your vcf is not properly formatted")

def get_x_as_string(x):

    """gets an x and returns the string of it"""

    if type(x) in [str, float, int]: return str(x)
    elif type(x) in [tuple, list, set]: return ",".join([str(y) for y in x])
    else: return str(x)

def get_correct_INFO_with_bendIDs_and_bendStats(r, df_gridss):

    """Takes a row of the final df_CNV and returns the INFO with the breakend information and the best breakend IDs"""

    # copy dfs
    r = cp.deepcopy(r)

    # set the ID as index
    df_gridss = df_gridss.set_index("ID", drop=False)
    check_that_df_index_is_unique(df_gridss)

    # get the info dict
    info = get_INFO_dict_from_INFO_string(r["INFO"])
    if any({not k.startswith("INFO_") for k in info}): raise ValueError("info is not correct")
    info = {k.split("INFO_")[1] : v for k,v in info.items()}

    ######### GET THE LIST OF BREAKENDS #########

    # get the breakend list
    if "BREAKENDIDs" in info.keys():

        if info["BREAKENDIDs"]=="wholeChrom": breakend_IDs = ["wholeChrom"]
        else: breakend_IDs = info["BREAKENDIDs"].split(",")

    elif "BREAKPOINTIDs" in info.keys():

        # define the interesting df_gridss
        breakpoint_IDs = set(info["BREAKPOINTIDs"].split(","))
        df_gridss = df_gridss[(df_gridss.eventID_as_clove.isin(breakpoint_IDs)) & (df_gridss["#CHROM"]==r["#CHROM"])]
        if len(df_gridss)==0: raise ValueError("there should only be one ID")

        # define the positions where the breakend should be found
        if "END" in info: 

            # sort df_gridss by pos
            df_gridss = df_gridss.sort_values(by="POS")

            # if there are two positions, find the breakpoint that best matches the POS-to-end regime
            def get_score_breakpoint_matchingVariant(df_bp):

                # if there is only one breakend, return a negative score
                if len(df_bp)==1: return -1

                # if there are two breakends return the inverted mean distance to the start and end
                elif len(df_bp)==2: 
                    first_pos = df_bp.iloc[0].POS
                    second_pos =  df_bp.iloc[1].POS
                    return 1 /(np.mean([abs(r["POS"]-first_pos), abs(info["END"]-second_pos)]) + 1)

                else: raise ValueError("df_bp should have 1 or 2 breakends")

            bpointID_to_score = df_gridss.groupby("eventID_as_clove").apply(get_score_breakpoint_matchingVariant)

            # if there is some breakpoint with two breakends, it should be the one that is about this breakpoint
            if any(bpointID_to_score>0):

                best_breakpoint = bpointID_to_score[bpointID_to_score==max(bpointID_to_score)].index[0]

                # keep the breakend IDs of the best breakpoints
                best_bp_df_gridss = df_gridss[df_gridss.eventID_as_clove==best_breakpoint].sort_values(by="POS")
                if len(best_bp_df_gridss)!=2: 
                    raise ValueError("There should be 2 breakpoints in %s"%best_bp_df_gridss)

                # define the breakend IDs
                breakend_IDs = [best_bp_df_gridss.ID.iloc[0], best_bp_df_gridss.ID.iloc[1]]

            # any breakend, even if it is not from the same breakpoint can be interesting
            else:

                # get the breakends that are closest to the positions
                breakend_IDs = [sorted([(bendID, abs(r_bend["POS"]-target_pos)) for bendID, r_bend in df_gridss.iterrows()], key=(lambda x: x[1]))[0][0] for target_pos in [r["POS"], info["END"]]]

                # check that these are two different breakends
                if len(set(breakend_IDs))!=2: raise ValueError("the bend IDs are not unique")

        else: 

            # if there is only one position, find the closest in df_gridss and with the highest QUAL
            bendID_df = df_gridss[df_gridss.POS==find_nearest(df_gridss.POS, r["POS"])].sort_values(by=["QUAL"], ascending=False)

            breakend_IDs = [bendID_df.ID.iloc[0]]

        if len(set(breakend_IDs))!=len(breakend_IDs): raise ValueError("there should be one breakend per position")

    # in the CNV ones there are no breakend IDs, so you should just without breakendIDs
    elif r["ID"].split("|")[0] in  {"coverageDEL", "coverageDUP"}: breakend_IDs = ["."]

    else: raise ValueError("info is not valid. This is it:\n-------\n %s\n-------\n"%info)

    # add them to info
    info["BREAKENDIDs"] = ",".join(breakend_IDs)

    #############################################

    ######### ADD THE BREKEND STATS BASED ON THE BREAKPOINTS #########

    # for whole chromosomes, just keep the best
    if breakend_IDs!=["wholeChrom"] and breakend_IDs!=["."]: 

        # get the quantitative fields
        gridss_quantitative_fields=["allele_frequency", "allele_frequency_SmallEvent", "real_AF", "length_inexactHomology", "length_microHomology", "QUAL", "length_event", "len_inserted_sequence"]

        for estimate_fn_name, estimate_fn in [("min", min), ("max", max), ("mean", np.mean)]:

            # go throug each field
            for quant_field in gridss_quantitative_fields:

                # add to info
                field = "%s_%s"%(quant_field, estimate_fn_name)
                info[field] = estimate_fn([df_gridss.loc[bendID, quant_field] for bendID in breakend_IDs])

        # get the qualitative fields
        filter_to_int = {"LOW_QUAL":0, "REF":1, "INSUFFICIENT_SUPPORT":2, "NO_ASSEMBLY":3, "ASSEMBLY_TOO_SHORT":4, "ASSEMBLY_TOO_FEW_READ":5, "SINGLE_ASSEMBLY":6, "ASSEMBLY_ONLY":7, "PASS":8}

        # get only the worst filter
        df_gridss["FILTER"] = df_gridss.FILTER.apply(lambda x: sorted(x.split(";"), key=(lambda x: filter_to_int[x]))[0])

        # get the best and worse filters
        sorted_filters = sorted([df_gridss.loc[bendID, "FILTER"] for bendID in breakend_IDs], key=(lambda x: filter_to_int[x]))

        info["best_FILTER"] = sorted_filters[-1]
        info["worse_FILTER"] = sorted_filters[0]

        # get the boolean fields
        for f in ["has_poly16GC", "overlaps_repeats"]: info["any_%s"%f] = any([df_gridss.loc[bendID, f] for bendID in breakend_IDs])

        # add all the breakend fields
        for f in gridss_quantitative_fields + ["FILTER", "has_poly16GC", "overlaps_repeats", "coordinates"]:

            info["BREAKEND_%s"%f] = ",".join([str(df_gridss.loc[bendID, f]) for bendID in breakend_IDs])

    ##################################################################

    # get the INFO as a string
    return ";".join(["%s=%s"%(k, get_x_as_string(v)) for k, v in info.items()])


def clean_sorted_bam_coverage_per_window_files(sorted_bam):

    """Removes all files that start with sorted bam and  followed by coverage_per_window"""

    print_if_verbose("cleaning coverage files")

    sorted_bam_dir = get_dir(sorted_bam)

    for f in os.listdir(sorted_bam_dir):

        path = "%s/%s"%(sorted_bam_dir, f)

        if path.startswith("%s.coverage_per_window."%sorted_bam): remove_file(path)


def get_df_vcf_with_df_CNV_coverage_added_nonRedundant(sorted_bam, reference_genome, mitochondrial_chromosome, df_vcf, df_CNV, outdir, df_gridss, df_clove, threads, replace, window_size_CNVcalling, cnv_calling_algs):

    """This function merges the df_vcf with the coverage-based prediction, removing redudnant events."""

    # define the final file
    df_vcf_final_file = "%s/vcf_merged_CNVcalling_SVcalling.vcf"%outdir

    if file_is_empty(df_vcf_final_file) or replace is True:

        # define fields
        data_fields = ["chromosome", "start", "end", "ID", "SVTYPE", "INFO", "median95CI_lower_rel_coverage", "median95CI_higher_rel_coverage", "median95CI_lower_rel_coverage_relative", "median95CI_higher_rel_coverage_relative", "abs_spearman_r", "abs_pearson_r", "spearman_p", "pearson_p"]

        vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

        # calculate median cov
        destination_dir = "%s.calculating_windowcoverage"%sorted_bam
        coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=True, delete_bams=True), sep="\t")
        median_coverage = get_median_coverage(coverage_df, mitochondrial_chromosome)

        ########### GET RID OF REDUNDANT EVENTS AND ADD FIELDS ###########

        # add the ID
        df_CNV["ID"] = "coverage" + df_CNV.SVTYPE + "|" + df_CNV.chromosome + ":" + df_CNV.start.apply(str) + "-" + df_CNV.end.apply(str)

        # get the df_vcf related to CNV
        df_vcf_forCNV = df_vcf[df_vcf.ALT.isin({"<DUP>", "<TDUP>", "<DEL>"})].rename(columns={"POS":"start", "#CHROM":"chromosome"}).set_index("ID", drop=False)
        df_vcf_forCNV["end"] = df_vcf_forCNV.INFO.apply(lambda x: [int(y.split("END=")[1]) for y in x.split(";") if y.startswith("END")][0])

        # add the svtype
        svtype_to_DUPDEL = {"TDUP":"DUP", "DUP":"DUP", "DEL":"DEL"}
        df_vcf_forCNV["SVTYPE"] = df_vcf_forCNV.INFO.apply(lambda x: [svtype_to_DUPDEL[y.split("SVTYPE=")[1]] for y in x.split(";") if y.startswith("SVTYPE")][0])

        # add the type of SVcall
        df_vcf_forCNV["type_CNVcall"] = "gridssClove"
        df_CNV["type_CNVcall"] = "coverage"

        # get only non-redundant CNVs
        df_CNV.index = list(range(0, len(df_CNV)))
        df_CNV = get_nonRedundant_CNVcalls_coverage(outdir, df_CNV, df_vcf_forCNV, threads, replace, pct_overlap=0.8)

        ################################################################

        ###### FORMAT AS VCF ######

        # get the coverage calculation for the input vcf TAN,DUP,DEL
        if len(df_vcf_forCNV)==0: df_vcf_forCNV_final = pd.DataFrame(columns=data_fields)
        
        else:   

            df_vcf_forCNV_final  = df_vcf_forCNV.set_index("ID", drop=False)

            bed_windows_prefix = "%s/calculating_cov_neighbors_SV-based_vcf"%outdir
            df_vcf_forCNV_final = get_df_with_coverage_per_windows_relative_to_neighbor_regions(df_vcf_forCNV_final, bed_windows_prefix, reference_genome, sorted_bam, df_clove, median_coverage, replace=replace, run_in_parallel=True, delete_bams=True, threads=threads)
            df_vcf_forCNV_final = get_coverage_df_windows_with_within_windows_statistics(df_vcf_forCNV_final, outdir, sorted_bam, reference_genome, median_coverage, replace=replace, threads=threads)

            # change the SVTYPE to follow INFO. This is important to get TDUPs back in place
            df_vcf_forCNV_final["SVTYPE"] = df_vcf_forCNV_final.INFO.apply(lambda x: [y.split("SVTYPE=")[1] for y in x.split(";") if y.startswith("SVTYPE")][0])

        # add the INFO to and remaining data_fields to df_CNV
        if len(df_CNV)==0: df_CNV = pd.DataFrame(columns=data_fields)
        else:

            # add the field
            def get_INFO_from_df_CNV_r(r):

                # add the info
                info = "END=%i;SVTYPE=%s;merged_relative_CN=%.3f;median_coverage_corrected=%.3f"%(r["end"], r["SVTYPE"], r["merged_relative_CN"], r["median_coverage_corrected"])

                # add the calling of cnvs
                cnv_calling_algs_fields = ["median_relative_CN_%s"%alg for alg in cnv_calling_algs]
                info += ";%s"%(";".join(["%s=%.3f"%(f, r[f]) for f in cnv_calling_algs_fields]))

                return info 

            df_CNV["INFO"] = df_CNV.apply(get_INFO_from_df_CNV_r, axis=1)

            # filter out SVs that have a size below min_CNVsize_coverageBased
            df_CNV["length_CNV"] = df_CNV.end - df_CNV.start
            df_CNV = df_CNV[df_CNV.length_CNV>=min_CNVsize_coverageBased]

            # add the coverage fields
            bed_windows_prefix = "%s/calculating_cov_neighbors_CNV_vcf"%outdir
            df_CNV = get_df_with_coverage_per_windows_relative_to_neighbor_regions(df_CNV, bed_windows_prefix, reference_genome, sorted_bam, df_clove, median_coverage, replace=replace, run_in_parallel=True, delete_bams=True, threads=threads)
            df_CNV = get_coverage_df_windows_with_within_windows_statistics(df_CNV, outdir, sorted_bam, reference_genome, median_coverage, replace=replace, threads=threads)

        # initialize the final df
        df_vcf_final = df_CNV[data_fields].append(df_vcf_forCNV_final[data_fields])

        # add the INFO
        if len(df_vcf_final)==0: df_vcf_final["INFO"] = ""
        else:   

            df_vcf_final["INFO"] = df_vcf_final.apply(lambda r: "%s;RELCOVERAGE=%.4f,%.4f;RELCOVERAGE_NEIGHBOR=%.4f,%.4f;REGION_ABS_SPEARMANR=%.4f;REGION_ABS_PEARSONR=%.4f;REGION_SPEARMANP=%.4f;REGION_PEARSONP=%.4f"%(r["INFO"], r["median95CI_lower_rel_coverage"], r["median95CI_higher_rel_coverage"], r["median95CI_lower_rel_coverage_relative"], r["median95CI_higher_rel_coverage_relative"], r["abs_spearman_r"], r["abs_pearson_r"], r["spearman_p"], r["pearson_p"]), axis=1)

        # add the ALT
        df_vcf_final["ALT"] = "<" + df_vcf_final.SVTYPE + ">"

        # add empty fields
        for f in  ["REF", "QUAL", "FILTER", "FORMAT"]: df_vcf_final[f] = "."

        # rename fields
        df_vcf_final = df_vcf_final.rename(columns={"chromosome":"#CHROM", "start":"POS"})[vcf_fields]

        # append the initial vcf 
        df_vcf_noCNV = df_vcf[~(df_vcf.ALT.isin({"<DUP>", "<TDUP>", "<DEL>"}))]
        df_vcf_final = df_vcf_final[vcf_fields].append(df_vcf_noCNV[vcf_fields])

        ##########################

        # save
        save_df_as_tab(df_vcf_final[vcf_fields], df_vcf_final_file)

    # load
    df_vcf_final = get_tab_as_df_or_empty_df(df_vcf_final_file).sort_values(by=["#CHROM", "POS"])

    return df_vcf_final

def get_vcf_all_SVs_and_CNV(perSVade_outdir, outdir, sorted_bam, reference_genome, ploidy, df_CNV_coverage, window_size_CNVcalling, cnv_calling_algs, replace=False, threads=4, mitochondrial_chromosome="mito_C_glabrata_CBS138"):

    """This function generates a vcf that has all the variants and CNV"""

    # make the folder
    make_folder(outdir)

    # get the vcf SV calling
    vcf_SVcalling = "%s/SV_and_CNV_variant_calling.vcf"%outdir

    if file_is_empty(vcf_SVcalling) or replace is True:
        print_if_verbose("getting all CNV and SVs into one vcf")

        # clean the sorted bam coverage per window
        clean_sorted_bam_coverage_per_window_files(sorted_bam)

        # define the outdir
        outdir_gridss_final = "%s/SVdetection_output/final_gridss_running"%perSVade_outdir

        # define the clove outfile
        outfile_clove = "%s/gridss_output.vcf.withSimpleEventType.vcf.filtered_default.vcf.bedpe.raw.bedpe.clove.vcf"%outdir_gridss_final
        if file_is_empty(outfile_clove): outfile_clove = "%s/clove_output.vcf"%outdir_gridss_final

        # get the clove df
        df_clove = get_clove_output(outfile_clove)

        # get files from output
        svtype_to_svfile, df_gridss = get_svtype_to_svfile_and_df_gridss_from_perSVade_outdir(perSVade_outdir, reference_genome)

        ######## GET THE VCF OF SVs ########

        if len(svtype_to_svfile)==0:  

            vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            df_vcf = pd.DataFrame(columns=vcf_fields)

        else:

            # get the svDF metadata
            print_if_verbose("getting the svtype_to_svDF")
            svtype_to_svDF = get_sampleID_to_svtype_to_svDF_filtered({"x":svtype_to_svfile}, {"x":df_gridss}, sampleID_to_parentIDs={}, breakend_info_to_keep=['#CHROM', 'POS', 'other_coordinates', 'allele_frequency', 'allele_frequency_SmallEvent', 'real_AF', 'FILTER', 'inserted_sequence', 'has_poly16GC', 'length_inexactHomology', 'length_microHomology', 'QUAL', 'overlaps_repeats', 'REF', 'BREAKPOINTID'])["x"]

            print_if_verbose("svtype_to_svDF got")

            # get a vcf df, that comes from all vcfs
            df_vcf = pd.concat([get_vcf_df_for_svDF(svDF, svtype, reference_genome, df_gridss) for svtype, svDF in svtype_to_svDF.items() if svtype in {"tandemDuplications", "deletions", "inversions", "translocations", "insertions", "remaining"}])

        ####################################



        # add the df_CNV_coverage
        df_vcf = get_df_vcf_with_df_CNV_coverage_added_nonRedundant(sorted_bam, reference_genome, mitochondrial_chromosome, df_vcf, df_CNV_coverage, outdir, df_gridss, df_clove, threads, replace, window_size_CNVcalling, cnv_calling_algs)

        # add a tag to the ID, that makes it unique
        df_vcf[["ID", "INFO"]] = df_vcf.apply(get_correctID_and_INFO_df_vcf_SV_CNV, axis=1)

        # check that it is unique
        if len(df_vcf)!=len(set(df_vcf.ID)): raise ValueError("IDs are not unique")

        # add the POS and END that are correct, these should be 1-based. Note that they wont match the ID
        df_vcf["POS"] = df_vcf.apply(get_correct_POS_in1based, axis=1)

        # add to the END + 1
        chr_to_len = get_chr_to_len(reference_genome)
        df_vcf["INFO"] = df_vcf.apply(lambda r: get_correct_INFO_withEND_in1based(r, chr_to_len), axis=1)        
        
        # add the breakend IDs and the metadata info
        df_vcf["INFO"] = df_vcf.apply(lambda r: get_correct_INFO_with_bendIDs_and_bendStats(r, df_gridss), axis=1)

        # write vcf
        vcf_SVcalling_tmp = "%s.tmp"%vcf_SVcalling
        vcf_lines = df_vcf.to_csv(sep="\t", header=False, index=False)
        header_lines = "\n".join([l.strip() for l in open(outfile_clove, "r").readlines() if l.startswith("#CHROM") or l.startswith("##fileformat")])
        open(vcf_SVcalling_tmp, "w").write(header_lines + "\n" + vcf_lines)
        os.rename(vcf_SVcalling_tmp, vcf_SVcalling)

    return vcf_SVcalling
    
def clean_vep_output_files(run_vep_outfile):

    """Takes the --outfile passed to run_vep and cleans it"""

    outdir = get_dir(run_vep_outfile)
    for f in os.listdir(outdir):
        path = "%s/%s"%(outdir, f)
        if path.startswith(run_vep_outfile): remove_file(path)

def annotate_SVs_inHouse(vcf_SVcalling, gff, reference_genome, replace=False, threads=4, mitochondrial_chromosome="mito_C_glabrata_CBS138", mito_code=3, gDNA_code=1):

    """This function annotates the variants from vcf_SVcalling"""

    # run vep
    annotated_vcf = "%s_annotated_VEP.tab"%vcf_SVcalling
    annotated_vcf_tmp = "%s.tmp"%annotated_vcf

    if file_is_empty(annotated_vcf) or replace is True:

        # clea
        clean_vep_output_files(annotated_vcf_tmp)

        # run vep
        run_cmd("%s --input_vcf %s --outfile %s --ref %s --gff %s --mitochondrial_chromosome %s --mito_code %i --gDNA_code %i"%(run_vep, vcf_SVcalling, annotated_vcf_tmp, reference_genome, gff, mitochondrial_chromosome, mito_code, gDNA_code))

        # rename the files (also the warnings and summary)
        os.rename("%s.raw.tbl.tmp_summary.html"%annotated_vcf_tmp, "%s.raw.tbl.tmp_summary.html"%annotated_vcf)
        os.rename(annotated_vcf_tmp, annotated_vcf)

    clean_vep_output_files(annotated_vcf_tmp)

    return annotated_vcf


def get_defaultWay_median_coverage(sorted_bam, reference_genome, mitochondrial_chromosome, replace=False, threads=4):

    """Takes a sorted bam, reference, and mtDNA and returns the median coverage"""

    destination_dir = "%s.calculating_windowcoverage"%sorted_bam
    coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=True, delete_bams=True, threads=threads), sep="\t")


    median_coverage =  get_median_coverage(coverage_df, mitochondrial_chromosome)

    return median_coverage


def check_that_cloveIDs_are_in_gridss_df(df_clove, df_gridss):

    """This function takes a df clove and a df gridss and checks that all the df_clove IDs are in df_gridss"""

    df_clove = cp.deepcopy(df_clove)

    # get IDs set
    df_clove["IDs_set"] = df_clove.ID.apply(lambda x: set(re.split("\+|\-", x)))
    all_cloveIDs = set.union(*df_clove.IDs_set)

    # get the missing IDs
    missing_IDs = all_cloveIDs.difference(set(df_gridss.eventID_as_clove))
    if len(missing_IDs)>0: raise ValueError("%s can't be found in df_gridss"%missing_IDs)

    print_if_verbose("all clove IDs are in df_gridss")

def get_real_value_from_string(string):

    try: return ast.literal_eval(string)
    except: return string

def get_INFO_dict_from_INFO_string(INFOstring):

    """This function takes an INFO string and returns a dict with real values"""

    # get as a list, considering to put appart whatever has an "="
    list_row = INFOstring.split(";")
    rows_with_equal = [x.split("=") for x in list_row if "=" in x]
    rows_without_equal = [x for x in list_row if "=" not in x]

    # add the values with an "=" to the dictionary
    final_dict = {"INFO_%s"%x[0] : get_real_value_from_string(x[1]) for x in rows_with_equal}

    # add the others collapsed
    final_dict["INFO_misc"] = ";".join(rows_without_equal)

    return final_dict

def get_vcf_df_with_INFO_as_single_fields(df):

    """Takes a vcf df and returns a the same one where the INFO content is split acrros extra fields"""

    if len(df)==0: return df

    ### INFO COLUMN ####

    # add a column that has a dictionary with the info fields
    df["INFO_as_dict"] = df.INFO.apply(get_INFO_dict_from_INFO_string)
    all_INFO_fields = sorted(list(set.union(*df.INFO_as_dict.apply(lambda x: set(x)))))

    # add them as sepparated columns
    def get_from_dict_orNaN(value, dictionary):

        if value in dictionary: return dictionary[value]
        else: return np.nan

    for f in all_INFO_fields: df[f] = df.INFO_as_dict.apply(lambda d: get_from_dict_orNaN(f, d))
    df.pop("INFO_as_dict")

    #########################

    return df

def save_df_as_tab(df, file):

    """Takes a df and saves it as tab"""

    file_tmp = "%s.tmp"%file
    df.to_csv(file_tmp, sep="\t", index=False, header=True)
    os.rename(file_tmp, file)

def save_df_as_tab_with_index(df, file):

    """Takes a df and saves it as tab"""

    file_tmp = "%s.tmp"%file
    df.to_csv(file_tmp, sep="\t", index=True, header=True)
    os.rename(file_tmp, file)

def get_tab_as_df_or_empty_df(file):

    """Gets df from file or empty df"""

    nlines = len([l for l in open(file, "r").readlines() if len(l)>1])

    if nlines==0: return pd.DataFrame()
    else: return pd.read_csv(file, sep="\t")


def get_tab_as_df_or_empty_df_with_index(file):

    """Gets df from file or empty df"""

    nlines = len([l for l in open(file, "r").readlines() if len(l)>1])

    if nlines==0: return pd.DataFrame()
    else: return pd.read_csv(file, sep="\t", index_col=0)

def get_END_vcf_df_r_NaN_to_1(r):

    """Takes a row of a vcf df that has END and POS. If the END is NaN, it set's it to POS+1"""

    if pd.isna(r["INFO_END"]): return int(r["POS"]+1)
    else: return int(r["INFO_END"])

def get_varIDs_overlapping_target_regions(df_vcf, target_regions, outdir):

    """This function takes a df_vcf with #CHROM, POS and END. INFO_END can be NaN. It returns a set with the IDs that overlap target_regions"""

    df_vcf = cp.deepcopy(df_vcf)

    if len(df_vcf)==0: raise ValueError("vcf is empty")

    # get the END to be POS+1 if it is NaN
    if "INFO_END" in df_vcf.keys(): df_vcf["INFO_END"] = df_vcf.apply(get_END_vcf_df_r_NaN_to_1, axis=1)
    else: df_vcf["INFO_END"] = df_vcf.POS + 1

    # get the vcf to bed
    vcf_bed = "%s/variants_locations.bed"%outdir
    df_vcf[["#CHROM", "POS", "INFO_END", "ID"]].to_csv(vcf_bed, sep="\t", header=False, index=False)

    # get the target regions to bed
    target_bed = "%s/target_regions.bed"%outdir
    target_regions[["chromosome", "start", "end"]].to_csv(target_bed, sep="\t", header=False, index=False)

    # if the target regions are empty, define None as overlapping IDs
    if len(target_regions)==0: overlapping_IDs = set()

    else:

        # run bedtools to get the intersection
        intersection_vcf_bed = "%s/variant_locations_intersecting_targetRegions.bed"%outdir
        intersection_vcf_bed_stderr = "%s.generating.stderr"%intersection_vcf_bed
        print_if_verbose("running bedtools to get the variants that intersect the provided regions. The stderr is in %s"%intersection_vcf_bed_stderr)

        intersection_vcf_bed_tmp = "%s.tmp"%intersection_vcf_bed
        run_cmd("%s intersect -a %s -b %s -wa > %s 2>%s"%(bedtools, vcf_bed, target_bed, intersection_vcf_bed_tmp, intersection_vcf_bed_stderr))

        remove_file(intersection_vcf_bed_stderr)
        os.rename(intersection_vcf_bed_tmp, intersection_vcf_bed)

        # get into df
        df_vcf_intersection = pd.read_csv(intersection_vcf_bed, sep="\t", header=None, names=["chromosome", "start", "end", "ID"])

        # check that all IDs are in the beginning
        if len(set(df_vcf_intersection.ID).difference(set(df_vcf.ID)))>0: raise ValueError("There are missing IDs")

        # get the IDs that are overlapping
        if len(df_vcf_intersection)>0: overlapping_IDs = set.union(*df_vcf_intersection.ID.apply(lambda x: set(x.split(";"))))

        else: overlapping_IDs = set() 

    return overlapping_IDs



def remove_files_SV_CNVcalling(outdir):

    """Takes the perSVade outdir and removes all the files necessary to not repeat the SV-CNV calling"""

    print_if_verbose("removing the files related to the CNV SV calling")

    # remove the varcalling
    outdir_SV_calling = "%s/SVcalling_output"%outdir
    delete_folder(outdir_SV_calling)

    # get things under parameter optimistation
    parameter_optimisation_dir = "%s/SVdetection_output/parameter_optimisation"%outdir
    if not os.path.isdir(parameter_optimisation_dir): return

    files_to_remove = ["CNV_optimisation_df_benchmarking.tab", 
                       "CNV_optimisation_df_benchmarking.tab.best_train_parameters.tab",
                       "df_CNV_allKnownRegions.tab"
                       ]

    # add folders inside simulations
    for f in os.listdir(parameter_optimisation_dir):

        # get folder
        folder = "%s/%s"%(parameter_optimisation_dir, f)

        if f.startswith("simulation_"):

            for fCNV in os.listdir(folder):

                # get the folder
                folder_CNV = "%s/%s"%(folder, fCNV)
                if fCNV.startswith("CNVoptimisation_"): files_to_remove += ["%s/%s"%(f, fCNV)]

    # remove 
    for f in files_to_remove: delete_file_or_folder("%s/%s"%(parameter_optimisation_dir, f))


def get_small_variant_calling_withCNstate(varcall_file, df_CNV_coverage, replace=False):

    """This function rewrites varcall_file by adding a column called relative_CN. This includes which is the relative CN state of each variant. """

    print_if_verbose("Adding CNstate to small variant calling")

    # define a function that returns if the varcall_file contains the final "relative_CN" field
    def varcall_file_has_relative_CN(varcall_file): return "relative_CN" in str(subprocess.check_output("head -n 1 %s"%varcall_file, shell=True))

    # only repeat if relative_CN is not in varcall_file
    if not varcall_file_has_relative_CN(varcall_file) or replace is True:
        print_if_verbose("Adding CNstate to small variant calling")

        # load the variant calling df
        df_varcall = get_tab_as_df_or_empty_df(varcall_file).sort_values(by=["#CHROM", "POS"])

        # define the initial fields that can't be realted to relative_CN
        def get_field_is_similar_to_relative_CN(x): return (x.split(".")[0]=="relative_CN") and (set(x.split(".")[1:])=={"1"})
        initial_fields_varcall = [x for x in list(df_varcall.keys()) if x!="relative_CN" and not get_field_is_similar_to_relative_CN(x)]

        # define an outdir for the intersection processing
        outdir_adding_relativeCN = "%s.adding_relativeCN"%varcall_file; make_folder(outdir_adding_relativeCN)

        # create a file with the regions under CNV
        df_CNV_coverage["numericID"] = list(range(len(df_CNV_coverage)))
        df_CNV_coverage["bed_ID"] = df_CNV_coverage.numericID.apply(str) + "|" + df_CNV_coverage.merged_relative_CN.apply(str)
        cnv_bed = "%s/cnv_regions.bed"%outdir_adding_relativeCN
        df_CNV_coverage[["chromosome", "start", "end", "bed_ID"]].sort_values(by=["chromosome", "start", "end"]).to_csv(cnv_bed, sep="\t", index=False, header=False)

        # create a file with the variants
        df_varcall["start_bed"] = df_varcall.POS - 1
        df_varcall["end_bed"] = df_varcall.start_bed + 1
        variants_bed = "%s/variants.bed"%outdir_adding_relativeCN
        df_varcall[["#CHROM", "start_bed", "end_bed", "#Uploaded_variation"]].to_csv(variants_bed, sep="\t", index=False, header=False)

        # run bedmap to get the IDs of the CNVs that overlap each of the variants_bed
        bedmap_outfile = "%s/bedmap_outfile.txt"%outdir_adding_relativeCN
        bedmap_stderr = "%s/bedmap_stderr.txt"%outdir_adding_relativeCN

        print_if_verbose("running bedmap. The stderr is in %s"%bedmap_stderr)
        run_cmd("%s --delim '\t' --range %i --echo-map-id %s %s > %s 2>%s"%(bedmap, min_CNVsize_coverageBased, variants_bed, cnv_bed, bedmap_outfile, bedmap_stderr))

        # add to df
        def get_relative_CN_from_mappedID(mapID):
            mapID = mapID.strip()
            if pd.isna(mapID) or mapID=="": return 1.0
            else: return max([float(x.split("|")[1]) for x in mapID.split(";")])

        df_varcall["relative_CN"] = pd.Series(open(bedmap_outfile, "r").readlines(), index=df_varcall.index).apply(get_relative_CN_from_mappedID)

        # debug
        if any(pd.isna(df_varcall.relative_CN)): raise ValueError("relative_CN can't be NaN")

        # delete the outdir
        delete_folder(outdir_adding_relativeCN)

        # save
        varcall_file_tmp = "%s.tmp"%varcall_file
        save_df_as_tab(df_varcall[initial_fields_varcall + ["relative_CN"]], varcall_file_tmp)
        os.rename(varcall_file_tmp, varcall_file)

    else: print_if_verbose("relative_CN was already added")

#######################################################################################
#######################################################################################
#######################################################################################
###################### EXTRA FUNCTIONS, NOT USED BY THE PIPELINE ######################
#######################################################################################
#######################################################################################
#######################################################################################


def run_jobarray_file_Nord3(jobs_filename, name, time="12:00:00", queue="bsc_ls", threads_per_job=4, RAM_per_thread=1800, max_njobs_to_run=1000):

    """Runs jobs_filename in Nord3"""

   # define dirs
    outdir = get_dir(jobs_filename)
    stddir = "%s/STDfiles"%outdir; 

    delete_folder(stddir); make_folder(stddir)

    # define the std files
    stderr_file = "%s/%s_stderr.txt"%(stddir, name)
    stdout_file = "%s/%s_stdout.txt"%(stddir, name)

    # define the job script
    jobs_filename_run = "%s.run"%jobs_filename

    # define the number of jobs in the job array
    njobs = len(open(jobs_filename, "r").readlines())

    # change the tome
    time = ":".join(time.split(":")[0:2])

    # define the number of jobs to run at once, maximum 1000
    njobs_to_run = min([max_njobs_to_run, njobs])

    arguments =  ["#!/bin/sh",
                  "#BSUB -e  %s"%stderr_file,
                  "#BSUB -o %s"%stdout_file,
                  "#BSUB -cwd %s"%outdir,
                  "#BSUB -W %s"%time,
                  "#BSUB -q %s"%queue,
                  "#BSUB -n %i"%threads_per_job, # the number of processes
                  "#BSUB -M %i"%RAM_per_thread, # the ram per thread in Mb
                  "#BSUB -J %s[1-%i]"%(name, njobs_to_run), # the job array
                  "",
                  "cmdFile=%s/command.$LSB_JOBINDEX &&"%stddir, # define the command file
                  "head -n $LSB_JOBINDEX %s | tail -n 1 > $cmdFile &&"%jobs_filename, # create cmdFile
                  'bash $cmdFile && echo "$cmdFile finished well"' # execute it
                  ]

    # '#BSUB -R "span[ptile=16]"', "export OMP_NUM_THREADS=16"
    
    # define and write the run filename
    with open(jobs_filename_run, "w") as fd: fd.write("\n".join(arguments))
    
    # run in cluster if specified
    print("Submiting jobfile to the cluster from stddir %s"%(stddir))
    run_cmd("bsub < %s"%jobs_filename_run)


def run_jobarray_file_MN4_greasy(jobs_filename, name, time="12:00:00", queue="bsc_ls", threads_per_job=4, nodes=1):

    """
    This function takes a jobs filename and creates a jobscript with args (which is a list of lines to be written to the jobs cript). It works in MN4 for greasy    

    """

    # define dirs
    outdir = get_dir(jobs_filename)
    stddir = "%s/STDfiles"%outdir; 

    delete_folder(stddir); make_folder(stddir)

    # define the std files
    greasy_logfile = "%s/%s_greasy.log"%(stddir, name)
    stderr_file = "%s/%s_stderr.txt"%(stddir, name)
    stdout_file = "%s/%s_stdout.txt"%(stddir, name)

    # define the job script
    jobs_filename_run = "%s.run"%jobs_filename

    # define the number of jobs in the job array
    njobs = len(open(jobs_filename, "r").readlines())

    # define the requested nodes. Each node has 48 threads
    max_nodes = max([int((njobs*threads_per_job)/48), 1])
    requested_nodes = min([nodes, max_nodes])

    # define the number of tasks
    max_ntasks = int((requested_nodes*48)/threads_per_job)
    ntasks = min([njobs, max_ntasks])

    # define the arguments
    arguments = [ "#!/bin/sh",
                  "#SBATCH --error=%s"%stderr_file,
                  "#SBATCH --output=%s"%stdout_file,
                  "#SBATCH --job-name=%s"%name, 
                  "#SBATCH --get-user-env",
                  "#SBATCH --workdir=%s"%outdir,
                  "#SBATCH --time=%s"%time,
                  "#SBATCH --qos=%s"%queue,
                  "#SBATCH --cpus-per-task=%i"%threads_per_job,
                  "#SBATCH --ntasks=%i"%ntasks,
                  "",
                  "module load greasy",
                  "export GREASY_LOGFILE=%s;"%(greasy_logfile),
                  "echo 'running pipeline';",
                  "greasy %s"%jobs_filename
                ]

    # define and write the run filename
    with open(jobs_filename_run, "w") as fd: fd.write("\n".join(arguments))
    
    # run in cluster if specified
    run_cmd("sbatch %s"%jobs_filename_run)


def get_current_clusterName_mareNostrum():

    """Returns the cluster name in which you are. It can be MN4, Nord3"""

    # try MN4
    try: 
        run_cmd("squeue > ~/misc_std.txt 2>&1")
        cluster_name = "MN4"

    except:

        # try Nord3
        try:

            run_cmd("bjobs > ~/misc_std.txt 2>&1")
            cluster_name = "Nord3"

        except: 

            raise ValueError("cluster could not be identified")

    return cluster_name


def get_integrated_variants_into_one_df(df, file_prefix, replace=False, remove_files=False):

    """This function takes a df (or a .tab file) with the following fields:
    
       - smallVars_vcf. The vcf of the small vars provided by perSVade
       - smallVars_var_annotation. The annotation of smallVars_vcf

       - SV_CNV_vcf. The vcf of the integrated SV and CNV calling by perSVade.
       - SV_CNV_var_annotation. The annotation of SV_CNV_vcf

       Either of the two should be provided

       The index of the df shouold be the sample ID.

    """

    # define files
    small_vars_file = "%s_small_vars.tab"%file_prefix
    small_vars_annot_file = "%s_small_vars_annot.tab"%file_prefix
    SV_CNV_file = "%s_SV_CNV.tab"%file_prefix
    SV_CNV_annot_file = "%s_SV_CNV_annot.tab"%file_prefix

    # define the expected files
    expected_files = set()
    if "smallVars_vcf" in df.keys(): expected_files.update({small_vars_file, small_vars_annot_file})
    if "SV_CNV_vcf" in df.keys(): expected_files.update({SV_CNV_file, SV_CNV_annot_file})

    if any([file_is_empty(f) for f in expected_files]) or replace is True:

        # initialize dfs
        small_vars = pd.DataFrame()
        small_vars_annot = pd.DataFrame()
        SV_CNV = pd.DataFrame()
        SV_CNV_annot = pd.DataFrame()

        # go through each sample
        for sampleID, r in df.iterrows():

            # small vars
            if "smallVars_vcf" in r.keys():

                # get dfs
                vcf = get_vcf_df_with_INFO_as_single_fields(get_df_and_header_from_vcf(r["smallVars_vcf"])[0])
                annotation = pd.read_csv(r["smallVars_var_annotation"], sep="\t")

                # add to df
                vcf["sampleID"] = sampleID

                small_vars = small_vars.append(vcf)
                small_vars_annot = small_vars_annot.append(annotation)

            # SVs
            if "SV_CNV_vcf" in df.keys():

                # get dfs
                vcf = get_vcf_df_with_INFO_as_single_fields(get_df_and_header_from_vcf(r["SV_CNV_vcf"])[0])
                annotation = pd.read_csv(r["SV_CNV_var_annotation"], sep="\t")

                # add to df
                vcf["sampleID"] = sampleID

                SV_CNV = SV_CNV.append(vcf)
                SV_CNV_annot = SV_CNV_annot.append(annotation)

        # drop duplicates from the annotations
        small_vars_annot = small_vars_annot.drop_duplicates()
        SV_CNV_annot = SV_CNV_annot.drop_duplicates()

        # save each of them
        save_df_as_tab(small_vars, small_vars_file)
        save_df_as_tab(small_vars_annot, small_vars_annot_file)
        save_df_as_tab(SV_CNV, SV_CNV_file)
        save_df_as_tab(SV_CNV_annot, SV_CNV_annot_file)

    # load them
    small_vars = get_tab_as_df_or_empty_df(small_vars_file)
    small_vars_annot = get_tab_as_df_or_empty_df(small_vars_annot_file)
    SV_CNV = get_tab_as_df_or_empty_df(SV_CNV_file)
    SV_CNV_annot = get_tab_as_df_or_empty_df(SV_CNV_annot_file)


    if remove_files is True:
        for f in [small_vars_file, small_vars_annot_file, SV_CNV_file, SV_CNV_annot_file]: remove_file(f)


    return small_vars, small_vars_annot, SV_CNV, SV_CNV_annot

def run_perSVade_severalSamples(paths_df, cwd, common_args, threads=4, sampleID_to_parentIDs={}, samples_to_run=set(), repeat=False, job_array_mode="job_array", ploidy=1, get_integrated_dfs=True, repeat_job_running=False):

 
    """
    This function inputs a paths_df, which contains an index as 0-N rows and columns "reads", "sampleID", "readID"  and runs the perSVade pipeline without repeating steps (unless indicated). pths_df can also be a tab-sepparated file were there are this 3 fields. The trimmed_reads_dir has to be the full path to the .fastq file. The sampleID has to be the unique sample identifier and the readID has to be R1 or R2 for paired-end sequencing. The p

    - cwd is the current working directory, where files will be written
    - repeat is a boolean that indicates whether to repeat all the steps of this function
    - threads are the number of cores per task allocated. In mn, you can use 48 cores per node. It has not been tested whether more cores can be used per task
    - samples_to_run is a set of samples for which we want to run all the pipeline
    - job_array_mode can be 'job_array' or 'local'. If local each job will be run after the other
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
        #success_files = ["%s/smallVars_CNV_output/variant_annotation_ploidy%i.tab"%(outdir, ploidy)]
        success_files = ["%s/perSVade_finished_file.txt"%(outdir)]

        # if repeat job running is true, remove the success_files
        if repeat_job_running is True: 
            for f in success_files: remove_file(f)

                   
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

        elif job_array_mode=="job_array":

            print_if_verbose("Submitting %i jobs to cluster ..."%len(all_cmds))
            jobs_filename = "%s/jobs.run_SNPs_CNV"%cwd
            open(jobs_filename, "w").write("\n".join(all_cmds))

            generate_jobarray_file(jobs_filename, "perSVade_severalSamples")

        else: raise ValueError("%s is not a valid job_array_mode"%job_array_mode)

        return False

    else: return True


def get_bed_df_from_variantID(varID):

    """Takes a variant ID, such as the ones in SV_CNV vcf 'INFO_variantID'. It returns a df with all chromosome-start-end information that should be matched to be considered as the same variant.

    The chromosomes always have extra information, such as the svtype, which makes that the future bedmap run will not match events from the same chromosome but different svtypes."""

    # get the ID svtype
    svtype = varID.split("|")[0]

    # inferred by coverage 
    if svtype in {"coverageDUP", "coverageDEL", "TDUP", "DEL", "INV"} : 

        chrom = "%s_%s"%(svtype, varID.split("|")[1].split(":")[0])
        start = int(varID.split("|")[1].split(":")[1].split("-")[0])
        end = int(varID.split("|")[1].split(":")[1].split("-")[1])
        type_overlap = "both" # this means that both the positions and fraction of overlap should be matching

        dict_bed = {0 : {"chromosome":chrom, "start":start, "end":end, "ID":varID, "type_overlap":type_overlap}}

    elif svtype.endswith("like"):

        posA, posB = varID.split("|")[1].split("-")

        chromA = "%s_%s"%(svtype, posA.split(":")[0])
        chromB = "%s_%s"%(svtype, posB.split(":")[0])

        startA = int(posA.split(":")[1])
        endA = startA + 1

        startB = int(posB.split(":")[1])
        endB = startB + 1

        type_overlap = "pos" # this means that only the position should be mapping

        dict_bed = {0 : {"chromosome":chromA, "start":startA, "end":endA, "ID":varID+"-A", "type_overlap":type_overlap},
                    1 : {"chromosome":chromB, "start":startB, "end":endB, "ID":varID+"-B", "type_overlap":type_overlap}}

    elif svtype in {"INS"}:

        regionA, posB, typeIns = varID.split("|")[1:]

        chromA = "%s_%s_%s"%(svtype, typeIns, regionA.split(":")[0])
        startA, endA = [int(x) for x in regionA.split(":")[1].split("-")]

        chromB = "%s_%s_%s"%(svtype, typeIns, posB.split(":")[0])
        startB = int(posB.split(":")[1])
        endB = startB + 1

        # here the A region is copied or pasted into a B breakend. This means that the type_overlap is different
        dict_bed = {0 : {"chromosome":chromA, "start":startA, "end":endA, "ID":varID+"-A", "type_overlap":"both"},
                    1 : {"chromosome":chromB, "start":startB, "end":endB, "ID":varID+"-B", "type_overlap":"pos"}}

    # complex inverted duplication. A region is duplicated, inverted and inserted into another region of the genome
    elif svtype in {"CVD"}:

        # get the region A nd the posB
        regionA, posB = varID.split("|")[1:]

        chromA = "%s_%s"%(svtype, regionA.split(":")[0])
        startA, endA = [int(x) for x in regionA.split(":")[1].split("-")]

        chromB = "%s_%s"%(svtype, posB.split(":")[0])
        startB = int(posB.split(":")[1])
        endB = startB + 1

        # here the A region is copied, inverted and pasted into a B breakend. This means that the type_overlap is different
        dict_bed = {0 : {"chromosome":chromA, "start":startA, "end":endA, "ID":varID+"-A", "type_overlap":"both"},
                    1 : {"chromosome":chromB, "start":startB, "end":endB, "ID":varID+"-B", "type_overlap":"pos"}}

    elif svtype in {"TRA"}:

        posA, posB = varID.split("|")[1].split("<>")

        chromA = "TRA_%s"%(posA.split(":")[0])
        chromB = "TRA_%s"%(posB.split(":")[0])

        # get positions as the minimum non-0 position 
        locationsA = {int(x) for x in posA.split(":")[1].split("-")}.difference({0})
        locationsB = {int(x) for x in posB.split(":")[1].split("-")}.difference({0})

        startA = min(locationsA)
        endA = startA + 1

        startB = min(locationsB)
        endB = startB + 1

        type_overlap = "pos" # this means that only the position should be mapping

        dict_bed = {0 : {"chromosome":chromA, "start":startA, "end":endA, "ID":varID+"-A", "type_overlap":type_overlap},
                    1 : {"chromosome":chromB, "start":startB, "end":endB, "ID":varID+"-B", "type_overlap":type_overlap}}


    else: raise ValueError("%s has not been parsed"%varID)


    # get as df
    df_bed = pd.DataFrame(dict_bed).transpose()

    # add the variantID, which will be useful to track, and is not necessarily unique
    df_bed["variantID"] = varID

    return df_bed


def get_SV_CNV_df_with_common_variantID_acrossSamples(SV_CNV, outdir, pct_overlap, tol_bp):

    """
    Takes a SV_CNV df and returns it with the field 'variantID_across_samples'. It uses bedmap to be particularly efficient. The basis of this is that if two variants are of the same type and overlap by pct_overlap or tol_bp they are called to be the same.
    """

    make_folder(outdir)

    # get all variantIDs
    all_variantIDs = SV_CNV[["INFO_variantID"]].drop_duplicates()["INFO_variantID"]

    # create a bed with all the regions that need to be matching in order to be considered the same
    df_bed_all = pd.concat(map(get_bed_df_from_variantID, all_variantIDs)).sort_values(by=["chromosome", "start", "end"])
    df_bed_all.index = list(range(0, len(df_bed_all)))

    # write the bed
    variants_bed = "%s/variants_regions.bed"%outdir
    df_bed_all[["chromosome", "start", "end", "ID"]].to_csv(variants_bed, sep="\t", index=False, header=False)

    ######### RUN BEDMAP TO MAP VARIANTS TO EACH OTHER #########

    # define the stderr
    bedmap_stderr = "%s/running_bedmap_.stderr"%outdir

    # run bedmap tol_bp
    bedmap_outfile = "%s/bedmap_outfile_range.txt"%outdir
    run_cmd("%s --delim '\t' --range %i --echo-map-id %s > %s 2>%s"%(bedmap, tol_bp, variants_bed, bedmap_outfile, bedmap_stderr))

    df_overlap_tolBp  = pd.read_csv(bedmap_outfile, sep="\t", header=None, names=["overlapping_IDs"])
    df_overlap_tolBp["ID"] = list(df_bed_all.ID)
    df_overlap_tolBp = df_overlap_tolBp.set_index("ID", drop=False)

    # run bedmap pct overlap
    bedmap_outfile = "%s/bedmap_outfile_pctOverlap.txt"%outdir
    run_cmd("%s --delim '\t' --fraction-both %.2f --echo-map-id %s > %s 2>%s"%(bedmap, pct_overlap, variants_bed, bedmap_outfile, bedmap_stderr))

    df_overlap_pctOverlap  = pd.read_csv(bedmap_outfile, sep="\t", header=None, names=["overlapping_IDs"])
    df_overlap_pctOverlap["ID"] = list(df_bed_all.ID)
    df_overlap_pctOverlap = df_overlap_pctOverlap.set_index("ID", drop=False)

    # get a unique df_overlap, which only considers overlaps that fullfill both the Bp and pct overlap
    df_overlap = pd.DataFrame()
    df_overlap["ID"] = df_overlap_tolBp.ID
    df_overlap = df_overlap.set_index("ID", drop=False)

    # at the IDs that overlap by either measures 
    df_overlap["IDs_overlap_pos"] = df_overlap_tolBp.loc[df_overlap.index, "overlapping_IDs"].apply(lambda x: set(x.split(";")))

    df_overlap["IDs_overlap_fraction"] = df_overlap_pctOverlap.loc[df_overlap.index, "overlapping_IDs"].apply(lambda x: set(x.split(";")))

    df_overlap["IDs_overlap_both"] = df_overlap.apply(lambda r: r["IDs_overlap_pos"].intersection(r["IDs_overlap_fraction"]), axis=1)

    df_overlap["IDs_overlap_any"] = df_overlap.apply(lambda r: r["IDs_overlap_pos"].union(r["IDs_overlap_fraction"]), axis=1)

    ###############################################

    ######## MAP EACH VARID TO THE OTHER IDs ########

    # map each variantID to the bedIDs
    varID_to_bedIDs = dict(df_bed_all.groupby("variantID").apply(lambda df_varID: set(df_varID.ID)))

    # map each bedID to the type overlap field
    bedID_to_typeOverlapF = dict("IDs_overlap_" + df_bed_all.set_index("ID")["type_overlap"])

    # init
    varID_to_overlapping_varIDs = {}

    # go through each variantID
    for varID_q, bedIDs_q in varID_to_bedIDs.items():

        # init 
        varID_to_overlapping_varIDs[varID_q] = set()

        # go through the targets, and ke
        for varID_t, bedIDs_t in varID_to_bedIDs.items():

            # In order to have the variants overlapping, all the query bedIDs should be overlapped by at least one target bedID according to bedID_to_typeOverlapF

            variants_are_overlapping = all([ any([ bedID_t in df_overlap.loc[bedID_q, bedID_to_typeOverlapF[bedID_q]] for bedID_t in bedIDs_t ]) for bedID_q in bedIDs_q ])

            if variants_are_overlapping is True: varID_to_overlapping_varIDs[varID_q].add(varID_t)

    # print
    for varID, o_varIDs in varID_to_overlapping_varIDs.items(): print(varID, o_varIDs)

    #################################################

    ########## ADD THE variantID_across_samples ##########

    # get the list of clusters of the variant IDs
    list_clusters_varIDs = get_list_clusters_from_dict(varID_to_overlapping_varIDs)

    # create a dict that maps each clusterID to the cluster
    clusterID_to_varIDs = {"cluster%i_%s"%(I, next(iter(varIDs)).split("|")[0]) : varIDs for I, varIDs in enumerate(list_clusters_varIDs)}

    # map each varID to the clusterID
    varID_to_clusterID = {}
    for clusterID, varIDs in clusterID_to_varIDs.items():
        for varID in varIDs: varID_to_clusterID[varID] = clusterID

    # add to the final df
    SV_CNV["variantID_across_samples"] = SV_CNV.INFO_variantID.apply(lambda x: varID_to_clusterID[x])

    #######################################################

    return SV_CNV

def load_gff3_intoDF(gff_path, replace=False):

    """ Takes the path to a gff and loads it into a df"""

    gff_df_file = "%s.df.tab"%gff_path

    if file_is_empty(gff_df_file) or replace is True:

        # define the number of rows to skip
        gff_lines = "%s.gff_lines.gff"%gff_path
        run_cmd("egrep -v '^#' %s > %s"%(gff_path, gff_lines))

        # define the gff fields
        gff_fields = ["chromosome", "source", "feature", "start", "end", "blank1", "strand", "blank2", "annotation"]

        # load
        print("loading gff")
        gff = pd.read_csv(gff_lines, header=None, names=gff_fields, sep="\t")

        # set all possible annotations
        all_annotations = set.union(*[set([x.split("=")[0] for x in an.split(";")]) for an in gff.annotation])

        def list_to_str(list_x):

            if list_x==[]: return ""
            else: return list_x[0]
        
        for anno in all_annotations:

            # define the field 
            anno_field = "ANNOTATION_%s"%anno

            # add the normalQ annotation field
            gff[anno_field] = gff.annotation.apply(lambda x: list_to_str([a.split(anno)[1].lstrip("=") for a in x.split(";") if anno in a]))

            # add the Dbxref sepparated
            if anno=="Dbxref": 

                # get all the dbxref fields
                all_Dbxrefs = set.union(*[set([x.split(":")[0] for x in dbxref.split(",")]) for dbxref in gff[anno_field]]).difference({""})

                # go through each dbxref field and add it to the df
                for dbxref in all_Dbxrefs: 

                    gff["%s_%s"%(anno_field, dbxref)] = gff[anno_field].apply(lambda x: list_to_str([d.split(dbxref)[1].lstrip(":") for d in x.split(",") if dbxref in d]))

        # get the ID
        gff["ID"] = gff.ANNOTATION_ID
        gff["Parent"] = gff.ANNOTATION_Parent

        # change the ID so that all of the features are unique, add numbers for non-unique IDs
        gff["duplicated_ID"] = gff.duplicated(subset="ID", keep=False) # marks all the duplicated IDs with a True
        gff["numeric_idx"] = list(range(0, len(gff)))

        def getuniqueIDs_gff(row):

            """Takes a row and changes the IDs if they are not unique"""
            if row["duplicated_ID"] is False: return row["ID"]
            else: return "%s-%i"%(row["ID"], row["numeric_idx"])

        gff["ID"] = gff.apply(getuniqueIDs_gff, axis=1)

        # check that it is correct
        if len(gff)!=len(set(gff["ID"])): raise ValueError("IDs are not unique in the gff")

        # set the id as index
        gff = gff.set_index("ID", drop=False)

        # add the upmost_parent (which should be the geneID)
        print("getting upmost parent")
        def get_utmost_parent(row, gff_df):

            """Takes a row and the gff_df and returns the highest parent. The index has to be the ID of the GFF"""

            # when you have found the parent it has no parent, so the ID is the upmost_parent
            if row["Parent"]=="": return row["ID"]

            # else you go to the parent
            else: return get_utmost_parent(gff_df.loc[row["Parent"]], gff_df)

        gff["upmost_parent"] = gff.apply(lambda row: get_utmost_parent(row, gff), axis=1)
        gff = gff.set_index("upmost_parent", drop=False)

        # add the type of upmost_parent
        df_upmost_parents = gff[gff.ID==gff.upmost_parent]

        # check that the upmost_parents are unique
        if len(df_upmost_parents)!=len(set(df_upmost_parents.ID)): raise ValueError("upmost parents are not unique")

        # map each upmost_parent to the feature
        upmost_parent_to_feature = dict(df_upmost_parents.set_index("upmost_parent")["feature"])

        # add the upmost_parent_feature to gff df
        gff["upmost_parent_feature"] = gff.upmost_parent.apply(lambda x: upmost_parent_to_feature[x])

        # write df
        gff_df_file_tmp = "%s.tmp"%gff_df_file
        gff.to_csv(gff_df_file_tmp, sep="\t", index=False, header=True)
        os.rename(gff_df_file_tmp, gff_df_file)


    # load
    gff = pd.read_csv(gff_df_file, sep="\t")

    return gff


def get_type_object(string):

    """Gets a string and returns the type of obhect"""

    if "." in string and string.replace('.','').isdigit(): return "float"
    elif string.isdigit(): return "int"
    else: return "str"

def get_value_as_str_forIDXmapping(val):

    """Gets a value as a string"""

    # type
    type_object = get_type_object(str(val))

    if pd.isna(val): return np.nan
    elif type_object=="float" and str(val).endswith(".0"): return str(int(val))
    elif type_object in {"float", "int"}: return str(val)
    elif type_object=="str": return str(val)
    else: raise ValueError("%s can't be parsed by this function"%val)

def get_annotation_df_with_GeneFeature_as_gff(annot_df, gff):

    """Takes the variant annotation df and returns them with added fields that are universal to all samples"""

    print_if_verbose("re-writing Gene and Feature")

    # load gff
    gff_df = load_gff3_intoDF(gff, replace=False)

    # keep only gffs where the upmost parent is a gene
    gff_df = gff_df[gff_df.upmost_parent_feature.isin({"gene", "pseudogene"})]

    # get empty 
    annot_fields = ["#Uploaded_variation", "Gene", "Feature"]
    if len(annot_df)==0: annot_df = pd.DataFrame(columns=annot_fields)

    # change the important fields to strings
    for f in [x for x in gff_df.keys() if x not in {"start", "end", "numeric_idx"}]: gff_df[f] = gff_df[f].apply(get_value_as_str_forIDXmapping)

    annot_df["Gene"] = annot_df.Gene.apply(get_value_as_str_forIDXmapping)
    annot_df["Feature"] = annot_df.Feature.apply(get_value_as_str_forIDXmapping)
   
    # find the fields in gff_df that explain the Gene and the feature
    for field in ["Gene", "Feature"]:

        # define all items
        elements_field = set(annot_df[field]).difference({"-"})

        # init macthing_gff_field, which will have to be defined
        gff_field_to_fraction_mapping = {}

        # init the already matched IDs
        already_matched_elements = set()

        ########## FIND ALL THE MATCHING GFF FIELDS ##########

        # go through all the gff fields
        for gff_field in gff_df.keys():

            # debug
            if all(pd.isna(gff_df[gff_field])): continue

            # define the elements
            elements_gff_field = {x for x in set(gff_df[gff_field]) if not pd.isna(x)}

            # define the elements in field and not in the gff
            elements_in_field_and_in_gff = elements_field.intersection(elements_gff_field)

            # get the fraction
            fraction_mapping = len(elements_in_field_and_in_gff)/len(elements_field)

            # if there is something mapping
            if fraction_mapping>0: 

                gff_field_to_fraction_mapping[gff_field] = fraction_mapping
                already_matched_elements.update(elements_in_field_and_in_gff)

        # debug
        if elements_field!=already_matched_elements: raise ValueError("some elements can't be mapped for %s"%field)

        ##########################################################

        ########## FIND THE MINIMUM NECESSARY GFF FIELDS ##########

        print_if_verbose("getting the minimum necessary GFF fields")

        # convert to series sorted so that the first elements go first
        gff_field_to_fraction_mapping = pd.Series(gff_field_to_fraction_mapping).sort_values(ascending=False)

        # init obhects
        final_gff_fields = []
        already_matched_elements = set()

        # go through each of the gff fields
        for gff_field in gff_field_to_fraction_mapping.index:

            # get the elements
            elements_gff_field = set(gff_df[gff_field])

            # define the elements in field and not in the gff
            elements_in_field_and_in_gff = elements_field.intersection(elements_gff_field)

            # new elements
            new_elements = elements_in_field_and_in_gff.difference(already_matched_elements)

            # add the field
            if len(new_elements)>0: final_gff_fields.append(gff_field)

            # keep
            already_matched_elements.update(elements_in_field_and_in_gff)

            # if  you already got everything, break
            if elements_field==already_matched_elements: break

        ###########################################################

        print_if_verbose("integrating")

        # define the     final_gff_field depending on the 
        target_gff_field = {"Gene":"upmost_parent", "Feature":"ANNOTATION_ID"}[field]

        # create a df that contains the important fields
        interesting_fields_gff_df = sorted(set([target_gff_field] + final_gff_fields))
        reduced_df_gff = gff_df[interesting_fields_gff_df].drop_duplicates()

        # add the intersecting_vals, which are the value of field that can be mapped to the gff_df
        reduced_df_gff["all_vals_set"] = reduced_df_gff[final_gff_fields].apply(set, axis=1)
        def get_intersection_sets(setA, setB): return setA.intersection(setB)
        all_field_vals = set(annot_df[field])
        reduced_df_gff["intersecting_set"] = reduced_df_gff.all_vals_set.apply(get_intersection_sets, setB=all_field_vals)

        # keep only df with some intersection
        reduced_df_gff = reduced_df_gff[reduced_df_gff.intersecting_set.apply(len)>0]

        # check that there are no multiple intersecting vals
        if any(reduced_df_gff["intersecting_set"].apply(len)!=1): raise ValueError("there are some values of gff that map to more than 1 feature")

        # get the set as a string
        reduced_df_gff["intesecting_val"] = reduced_df_gff.intersecting_set.apply(iter).apply(next)

        # map each originalID to the finalIDs
        originalID_to_finalID = reduced_df_gff.groupby("intesecting_val").apply(lambda df_or: ",".join(sorted(df_or[target_gff_field])))

        # add the '-'
        originalID_to_finalID["-"] = "-"

        # check that there are no multiple finalIDs of gene
        if field=="Gene" and any(originalID_to_finalID.apply(lambda x: "," in x)): 

            print(set(originalID_to_finalID))
            raise ValueError("Each gene should be only mapped to one feature")

        # add 
        annot_df[field] = annot_df[field].apply(lambda x: originalID_to_finalID[x])

    return annot_df

def get_overlapping_df_bedpe_multiple_samples(df_bedpe_all, outdir, tol_bp, pct_overlap, threads, replace=False):

    """Takes a df bedpe and returns a df with the overlapping breakpoints """

    # add the unique breakpointID
    df_bedpe_all["unique_bpID"] = df_bedpe_all.sampleID.apply(str) + "_" + df_bedpe_all.name

    # define file
    df_overlapping_BPs_file = "%s/df_bedpe_all_with_overlapping_BPs.py"%outdir
    initial_fields = list(df_bedpe_all.keys())

    if file_is_empty(df_overlapping_BPs_file) or replace is True:

        print("getting overlapping breakpoints")

        # this is a df were the rows are some target breakpoints (the ones that PASS the filters) and each column is a different breakpoint. The cell will be True if they are equivalent breakpoints

        if len(df_bedpe_all)!=len(set(df_bedpe_all.unique_bpID)): raise ValueError("The breakpoint IDs are not unique in the bedpe")
        df_bedpe_all = df_bedpe_all.set_index("unique_bpID", drop=False)

        # add the positions
        df_bedpe_all["pos1"] = (df_bedpe_all.start1 + (df_bedpe_all.end1-df_bedpe_all.start1)/2).apply(int)
        df_bedpe_all["pos2"] = (df_bedpe_all.start2 + (df_bedpe_all.end2-df_bedpe_all.start2)/2).apply(int)
        df_bedpe_all["pos1_plus1"] = df_bedpe_all["pos1"] + 1
        df_bedpe_all["pos2_plus1"] = df_bedpe_all["pos2"] + 1

        # add the type of breakpoint
        bool_to_text = {True:"intra_chromosomal", False:"inter_chromosomal"}
        df_bedpe_all["type_breakpoint"] = (df_bedpe_all.chrom1==df_bedpe_all.chrom2).map(bool_to_text)

        # add the merge of chromosome and strand
        df_bedpe_all["chrom_strand_1"] = df_bedpe_all.chrom1 + "_" + df_bedpe_all.strand1
        df_bedpe_all["chrom_strand_2"] = df_bedpe_all.chrom2 + "_" + df_bedpe_all.strand2
        df_bedpe_all["unique_bpID_1"] = df_bedpe_all.unique_bpID + "_1"
        df_bedpe_all["unique_bpID_2"] = df_bedpe_all.unique_bpID + "_2"

        #### GET BREAKPOINTS THAT OVERLAP BY POSITION #####

        # adds to df_bedpe_all an bpoints_overlapping_by_pos, which is a set of the breakpoints were both breakedns overlap by tol_bp

        # get a bed that has the positions of the breakends
        df_positions = pd.concat([df_bedpe_all[["chrom_strand_%i"%I, "pos%i"%I, "pos%i_plus1"%I, "unique_bpID_%i"%I]].rename(columns=dict(zip(["chrom_strand_%i"%I, "pos%i"%I, "pos%i_plus1"%I, "unique_bpID_%i"%I], ["chrom", "start", "end", "ID"]))) for I in [1, 2]]).sort_values(by=["chrom", "start", "end"])

        bed_positions = "%s/breakend_positions.bed"%outdir
        df_positions.to_csv(bed_positions, sep="\t", header=False, index=False)

        # define the stderr
        bedmap_stderr = "%s/running_bedmap.stderr"%outdir

        # run bedmap tol_bp. These are breakpoints that overlap by at leasttol_bp
        bedmap_outfile = "%s/bedmap_outfile_range.txt"%outdir
        run_cmd("%s --delim '\t' --range %i --echo-map-id %s > %s 2>%s"%(bedmap, tol_bp, bed_positions, bedmap_outfile, bedmap_stderr))

        # add to df bed
        df_overlap_tolBp  = pd.read_csv(bedmap_outfile, sep="\t", header=None, names=["overlapping_IDs"])
        df_overlap_tolBp["ID"] = list(df_positions.ID)
        df_positions = df_positions.merge(df_overlap_tolBp, on="ID", how="left", validate="one_to_one")
        if len(df_overlap_tolBp)!=len(df_positions): raise ValueError("the length should be as the bed") # debug

        # add the unique ID
        def get_bpointID(x): return "_".join(x.split("_")[0:-1])
        df_positions["breakpointID"] = df_positions.ID.apply(get_bpointID)

        # add the other bend ID
        number_to_ther_number = {"1":"2", "2":"1"}
        def get_other_bendID_name(bendID): return number_to_ther_number[bendID.split("_")[-1]]
        df_positions["other_bendID"] = df_positions.breakpointID + "_" + df_positions.ID.apply(get_other_bendID_name)

        # add the overlapping breakendIDs as a set
        def get_as_list(x): return x.split(";")
        df_positions["overlapping_bends"] = df_positions.overlapping_IDs.apply(get_as_list).apply(set)

        # map breakends
        bend_to_bpoint = dict(df_positions.set_index("ID").breakpointID)
        bend_to_otherBend = dict(df_positions.set_index("ID").other_bendID)
        bpoint_to_bends = df_positions.groupby("breakpointID").apply(lambda df_bp: set(df_bp.ID))
        bend_to_overlapping_bends = dict(df_positions.set_index("ID")["overlapping_bends"])

        # get the 
        def get_all_overlapping_bends_from_bpID(bpID): return set.union(*[bend_to_overlapping_bends[bend] for bend in bpoint_to_bends[bpID]])
        df_bedpe_all["overlapping_bends"] =  df_bedpe_all.unique_bpID.apply(get_all_overlapping_bends_from_bpID)

        def get_bpoints_overlapping_by_pos(all_overlapping_bends): return {bend_to_bpoint[bend] for bend in all_overlapping_bends if bend_to_otherBend[bend] in all_overlapping_bends}
        df_bedpe_all["bpoints_overlapping_by_pos"] = df_bedpe_all.overlapping_bends.apply(get_bpoints_overlapping_by_pos)

        ###################################################

        ############ DEFINE BREAKPOINTS OVERLAPPING BY PCT_OVERLAP ######
        print("geting pct_overlaps")

        # add to df_bedpe_all a bpoints_overlapping_by_region, just for the intrachromosomal events
 
        # get a bed that thas the regions of interchromosomal breakpoints
        df_intra_chrom = df_bedpe_all[df_bedpe_all.type_breakpoint=="intra_chromosomal"]
        if not all (df_intra_chrom.pos2>df_intra_chrom.pos1): raise ValueError("pos2 should be after pos1 in intrachromosomal breakpoints")

        bed_intrachromosomal_regions = "%s/intrachromosomal_regions.bed"%outdir
        df_intra_chrom["chrom_and_orientations"] = df_intra_chrom.chrom1 + "_" + df_intra_chrom.strand1 + "_" + df_intra_chrom.strand2
        df_intra_chrom = df_intra_chrom.sort_values(by=["chrom_and_orientations", "pos1", "pos2"])
        df_intra_chrom[["chrom_and_orientations", "pos1", "pos2", "unique_bpID"]].to_csv(bed_intrachromosomal_regions, sep="\t", header=False, index=False)

        # run bedmap pct overlap
        bedmap_outfile = "%s/bedmap_outfile_pctOverlap.txt"%outdir
        run_cmd("%s --delim '\t' --fraction-both %.2f --echo-map-id %s > %s 2>%s"%(bedmap, pct_overlap, bed_intrachromosomal_regions, bedmap_outfile, bedmap_stderr))

        df_overlap_pctOverlap  = pd.read_csv(bedmap_outfile, sep="\t", header=None, names=["overlapping_IDs"])
        df_overlap_pctOverlap["ID"] = list(df_intra_chrom.unique_bpID)
        df_overlap_pctOverlap["overlapping_IDs_set"] = df_overlap_pctOverlap.overlapping_IDs.apply(get_as_list).apply(set)

        # add the set to the breakpoint (also with missing vals)
        bpoint_to_overlappingBpoints_pctOverlap = df_overlap_pctOverlap.set_index("ID")["overlapping_IDs_set"]
        inter_chromosomal_bpoints = df_bedpe_all[df_bedpe_all.type_breakpoint=="inter_chromosomal"].unique_bpID
        bpoint_to_overlappingBpoints_pctOverlap = bpoint_to_overlappingBpoints_pctOverlap.append(pd.Series(dict(zip(inter_chromosomal_bpoints, [set()]*len(inter_chromosomal_bpoints)))))

        df_bedpe_all = df_bedpe_all.set_index("unique_bpID", drop=False)
        df_bedpe_all["bpoints_overlapping_by_region"] = bpoint_to_overlappingBpoints_pctOverlap.loc[df_bedpe_all.index]

        ###################################################################

        # get the overlap, which depends on the type
        def get_overlapping_bpoints(r):

            # define things
            by_region_bps = r["bpoints_overlapping_by_region"].difference({r["unique_bpID"]})
            by_pos_bps = r["bpoints_overlapping_by_pos"].difference({r["unique_bpID"]})

            if r.type_breakpoint=="intra_chromosomal": return by_region_bps.intersection(by_pos_bps)
            elif r.type_breakpoint=="inter_chromosomal": return by_pos_bps
            else: raise ValueError("r %s is not valid"%r)

        df_bedpe_all["overlapping_bpoints_final"] = df_bedpe_all.apply(get_overlapping_bpoints, axis=1)

        # save
        save_object(df_bedpe_all[initial_fields + ["overlapping_bpoints_final"]], df_overlapping_BPs_file)

    # load bepe with adds
    df_bedpe_all_withAdds = load_object(df_overlapping_BPs_file)

    return df_bedpe_all_withAdds

def get_SV_CNV_df_with_overlaps_with_all_samples(SV_CNV, outdir, tol_bp, pct_overlap, cwd, df_bedpe_all, df_CN_all, threads, reference_genome):

    """This function takes a SV_CNV df and returns the same df with some added fields about which other samples may also have that variant."""

    # define the initial fields
    initial_fields = list(SV_CNV.keys())

    # keep
    SV_CNV = cp.deepcopy(SV_CNV)

    # make the folder
    make_folder(outdir)

    # define a file
    SC_CNV_with_added_overlaps_file = "%s/SC_CNV_with_added_overlaps.tab"%outdir

    if file_is_empty(SC_CNV_with_added_overlaps_file):

        print("getting overlaps")

        # add things to the SV df
        SV_CNV["sampleID"] = SV_CNV.sampleID.apply(str)
        SV_CNV["is_coverage_SV"] = SV_CNV.INFO_variantID.apply(lambda ID: ID.startswith("coverage"))
        SV_CNV["INFO_BREAKPOINTIDs_set"] = SV_CNV.INFO_BREAKPOINTIDs.apply(lambda x: set(str(x).split(",")).difference({"nan"}))
        if any(SV_CNV[~SV_CNV.is_coverage_SV].INFO_BREAKPOINTIDs_set.apply(len)==0): raise ValueError("The bpID is not always defined ")
        SV_CNV["INFO_BREAKPOINTIDs_set_withSampleID"] = SV_CNV.apply(lambda r: {"%s_%s"%(r.sampleID, x) for x in r.INFO_BREAKPOINTIDs_set}, axis=1)
        SV_CNV["sampleID_and_variantID"] = SV_CNV.sampleID.apply(str) + "_" + SV_CNV.INFO_variantID

        ###### ADD FIELDS TO DF BEDPE ABOUT THE OVERLAP BETWEEN BREAKPOINTS ######

        # change the name so that it is as in clove
        df_bedpe_all["name"] = df_bedpe_all.name.apply(lambda x: x[0:-1]) + "o"

        # get the df bedpe with the breakpoints that are overlapping across several samples
        df_bedpe_all = get_overlapping_df_bedpe_multiple_samples(df_bedpe_all, outdir, tol_bp, pct_overlap, threads).set_index("unique_bpID", drop=False)

        # add the overlapping breakpoints PASS
        pass_breakpoints = set(df_bedpe_all[df_bedpe_all.PASSed_filters].unique_bpID)
        def get_only_PASSbps(all_bps): return all_bps.intersection(pass_breakpoints)
        df_bedpe_all["overlapping_bpoints_final_PASS"] = df_bedpe_all.overlapping_bpoints_final.apply(get_only_PASSbps)

        # add the samples that have that have the breakpoint
        bpID_to_SID = dict(df_bedpe_all.sampleID)
        def get_sampleIDs_from_bpIDs(bpIDs): return {bpID_to_SID[x] for x in bpIDs}
        df_bedpe_all["overlapping_samples_called"] = df_bedpe_all.overlapping_bpoints_final.apply(get_sampleIDs_from_bpIDs)
        df_bedpe_all["overlapping_samples_PASS"] = df_bedpe_all.overlapping_bpoints_final_PASS.apply(get_sampleIDs_from_bpIDs)

        ##########################################################################

        ###########  add the overlapping samples by breakpoint ###########

        # check that each variant ID has the same breakpoints
        for varID, bps in dict(SV_CNV[~SV_CNV.is_coverage_SV].groupby("sampleID_and_variantID").apply(lambda df_v: set(df_v.INFO_BREAKPOINTIDs))).items():
            if len(bps)!=1: 
                for idx, r in SV_CNV[SV_CNV.INFO_variantID==varID][["sampleID_and_variantID", "ID", "INFO_BREAKPOINTIDs", "sampleID"]].iterrows():
                    print(r["sampleID"], r["sampleID_and_variantID"], r["ID"], r["INFO_BREAKPOINTIDs"])
                raise ValueError("%s has more than 1 (%s) bpIDs"%(varID, bps))

        # define a function that returns the string with the overlapping samples
        def overlapping_samples_byBreakPoints(r, overlapping_samples_f):

            # this only applies to SVs, not CNVs
            if  r.is_coverage_SV is True: return ""
            else:

                all_breakponints_SV = r.INFO_BREAKPOINTIDs_set_withSampleID
                overlapping_samples = set(map(str, set.union(*df_bedpe_all.loc[all_breakponints_SV, overlapping_samples_f])))
                return ",".join(sorted(overlapping_samples))

        SV_CNV["overlapping_samples_byBreakPoints_allCalled"] = SV_CNV.apply(overlapping_samples_byBreakPoints, overlapping_samples_f="overlapping_samples_called", axis=1)

        SV_CNV["overlapping_samples_byBreakPoints_PASS"] = SV_CNV.apply(overlapping_samples_byBreakPoints, overlapping_samples_f="overlapping_samples_PASS", axis=1)

        # init the interesting fields
        final_interesting_fields = initial_fields + ["overlapping_samples_byBreakPoints_allCalled", "overlapping_samples_byBreakPoints_PASS"]

        ##################################################################

        ######## ADD WHETHER THE coverage CNVs overlap other samples ###########
        if df_CN_all is not None:

            print("adding overlaps with CNV")

            # redefine df_CN
            df_CN_all["sampleID"] = df_CN_all.sampleID.apply(str)
            df_CN_all = df_CN_all.sort_values(by=["sampleID", "chromosome", "start", "end"])
            df_CN_all["sorted_window_ID"] = list(range(len(df_CN_all)))

            # add the end to fit SV_CNV. Note that it only changes the end of the chromosome
            chrom_to_len = get_chr_to_len(reference_genome)
            def get_end_as_in_SV_CNV(r):
                if r.end==r.chrom_len: return r.end - 1
                else: return r.end

            df_CN_all["chrom_len"] = df_CN_all.chromosome.map(chrom_to_len)
            df_CN_all["end_as_in_SV_CNV"] = df_CN_all.apply(get_end_as_in_SV_CNV, axis=1)

            # create df_CN with several indices
            sampleChromStart_to_sortedWindowID = dict(df_CN_all.set_index(["sampleID", "chromosome", "start"], drop=False).sorted_window_ID)
            sampleChromEnd_to_sortedWindowID = dict(df_CN_all.set_index(["sampleID", "chromosome", "end_as_in_SV_CNV"], drop=False).sorted_window_ID)
            sortedWindowID_to_relativeCN = df_CN_all.set_index("sorted_window_ID", drop=False).merged_relative_CN

            # go through several pct overlaps
            for min_pct_overlap in [0.1, 0.25, 0.5, 0.75, 0.9, 0.95]: 
                print("getting CNVs overlapping by at least %.2f"%min_pct_overlap)

                # define samples
                all_samples = set(df_CN_all.sampleID)
                sample_to_otherSamples = {sampleID : all_samples.difference({sampleID}) for sampleID in all_samples}

                # define a function that takes an r of the SV_CNV and another sample name, it returns the name of the other sample if it is overlapping by at least min_pct_overlap
                def get_overlapping_sample_string(r, other_sample):

                    # define the relative CN for r
                    r_relative_CN = r.INFO_merged_relative_CN
                    r_start = r.POS-1
                    r_end = r.INFO_END-1

                    # define a df that has the region of the chromosome in other sample
                    other_start_windowID = sampleChromStart_to_sortedWindowID[(other_sample, r["#CHROM"], r_start)]
                    other_end_windowID = sampleChromEnd_to_sortedWindowID[(other_sample, r["#CHROM"], r_end)] + 1
                    other_relativeCN = sortedWindowID_to_relativeCN[other_start_windowID:other_end_windowID]

                    # debug
                    if len(other_relativeCN)==0: raise ValueError("The df_other can't be 0")
                   
                    # if it is a duplication, define the number of windows in other_relativeCN that have a coverage equal or above r_relative_CN
                    if r_relative_CN>1: n_similar_CNV_windows = sum(other_relativeCN>=r_relative_CN)

                    # if it is a deletion, define the number of windows in other_relativeCN that have a covereage equal or below r_relative_CN
                    elif r_relative_CN<1: n_similar_CNV_windows = sum(other_relativeCN<=r_relative_CN)
     
                    else: raise ValueError("The relative_CN can't be 1")

                    # if the ratio of similar windows is above the threshold, return the sample name of the other
                    if (n_similar_CNV_windows/len(other_relativeCN))>=min_pct_overlap: return other_sample
                    else: return ""
                
                def get_overlapping_samples_CNV_atLeast_pct_overlap(r):

                    # only work if it is a CNV variant
                    if r.is_coverage_SV is True: 

                        # get all the overlapping samples
                        all_samples_overlapping = set(map(lambda other_sample: get_overlapping_sample_string(r, other_sample), sample_to_otherSamples[r.sampleID])).difference({""})

                        return ",".join(sorted(all_samples_overlapping))

                    else: return ""

                # add the samples that are overlapping at least by min_pct_overlap windows of coverage
                field_overlapping_samples_CNV = "overlapping_samples_CNV_atLeast_%.2f"%min_pct_overlap
                SV_CNV[field_overlapping_samples_CNV] = SV_CNV.apply(get_overlapping_samples_CNV_atLeast_pct_overlap, axis=1)
                
                # keep field
                final_interesting_fields.append(field_overlapping_samples_CNV)

        #########################################################################

        # save
        save_df_as_tab(SV_CNV[final_interesting_fields], SC_CNV_with_added_overlaps_file)

    # load
    SV_CNV = get_tab_as_df_or_empty_df(SC_CNV_with_added_overlaps_file)

    return SV_CNV


def get_whetherCNVcalling_was_performed(VarCallOutdirs, samples_to_run):

    """Checks whether the CNV calling was correctly performed"""

    existing_cnv_files = [not file_is_empty("%s/%s/CNV_calling/final_df_coverage.tab"%(VarCallOutdirs, f)) for f in os.listdir(VarCallOutdirs) if f.split("_VarCallresults")[0] in samples_to_run]

    if sum(existing_cnv_files)!=0 and not all(existing_cnv_files): raise ValueError("not all CNV have been correctly called")

    return all(existing_cnv_files)

def get_integrated_SV_CNV_smallVars_df_from_run_perSVade_severalSamples(paths_df, cwd, ploidy, pct_overlap, tol_bp, gff, reference_genome, run_ploidy2_ifHaploid=False, generate_integrated_gridss_dataset=False, SV_CNV_called=True, threads=4):

    """Takes the same input as run_perSVade_severalSamples and writes the integrated dfs under cwd, adding to the integrated SV_CNV datasets some added fields that indicate that the SVs and CNV are shared among several samples.

    If generate_integrated_gridss_vcf is True it will generate an integrated gridss dataset that has all breakends and a field called "PASSed_filters" indicating if the breakend passed the filters .

    SV_CNV_called indicates whether SV_CNV were called.
    """

    # define the final filesfiles
    coverage_df_file = "%s/merged_coverage.tab"%cwd
    small_vars_df_file = "%s/smallVars_ploidy%i.tab"%(cwd, ploidy)
    SV_CNV_file = "%s/SV_CNV_ploidy%i.tab"%(cwd, ploidy)
    SV_CNV_file_simple = "%s/SV_CNV_ploidy%i.simple.tab"%(cwd, ploidy)
    SV_CNV_annot_file = "%s/SV_CNV_annot_ploidy%i.tab"%(cwd, ploidy)
    integrated_gridss_df = "%s/integrated_gridss_df_ploidy%i.tab"%(cwd, ploidy)
    integrated_bedpe_df = "%s/integrated_bedpe_df_ploidy%i.tab"%(cwd, ploidy)
    integrated_CNperWindow_df = "%s/integrated_CNperWindow.tab"%cwd

    if run_ploidy2_ifHaploid is True: small_var_annot_file = "%s/smallVars_annot_ploidy1and2.tab"%cwd
    else: small_var_annot_file = "%s/smallVars_annot_ploidy%i.tab"%(cwd, ploidy)

    if run_ploidy2_ifHaploid is True: small_vars_df_ploidy2_file = "%s/smallVars_ploidy2.tab"%cwd

    # define the samples to run
    samples_to_run = set(paths_df.sampleID)
    #samples_to_run = {"8", "74"} # debug        

    # define dirs
    VarCallOutdirs = "%s/VarCallOutdirs"%cwd


    # define the expected files
    if SV_CNV_called is True: expected_files = [coverage_df_file, small_vars_df_file, small_var_annot_file, SV_CNV_file_simple, SV_CNV_annot_file]
    else: expected_files = [coverage_df_file, small_vars_df_file, small_var_annot_file]

    # define whether there was some CNV calling
    CNV_calling_performed = get_whetherCNVcalling_was_performed(VarCallOutdirs, samples_to_run)
    if CNV_calling_performed is True: expected_files.append(integrated_CNperWindow_df)

    ######### GET THE SIMPLY MERGED DFS ##########

    # get the simple dataframes
    if any([file_is_empty(f) for f in expected_files]) or True:

        print_if_verbose("merging files")

        # debug+
        if run_ploidy2_ifHaploid is True and ploidy!=1: raise ValueError("run_ploidy2_ifHaploid can't be true if ploidiy!=1")

        # init dfs
        small_vars_df = pd.DataFrame()
        small_var_annot = pd.DataFrame()
        coverage_df = pd.DataFrame()
        if run_ploidy2_ifHaploid is True: small_vars_df_ploidy2 = pd.DataFrame()
        if CNV_calling_performed is True: df_CN = pd.DataFrame()

        # init dict
        sampleID_to_SV_dataDict = {}

        for sampleID in samples_to_run:
            print(sampleID)

            # create an outdir
            outdir = "%s/%s_VarCallresults"%(VarCallOutdirs, sampleID); make_folder(outdir)

            # add the CN df
            if CNV_calling_performed is True: 
                df_CN_sample = get_tab_as_df_or_empty_df("%s/CNV_calling/final_df_coverage.tab"%(outdir))
                df_CN_sample["sampleID"] = sampleID
                df_CN = df_CN.append(df_CN_sample)

            # add the small vars
            s_small_vars_df = pd.read_csv("%s/smallVars_CNV_output/variant_calling_ploidy%i.tab"%(outdir, ploidy), sep="\t")
            s_small_vars_df["sampleID"] = sampleID
            small_vars_df = small_vars_df.append(s_small_vars_df)

            # add the ploidy 2 small vars
            if run_ploidy2_ifHaploid is True:

                s_small_vars_df_ploidy2 = pd.read_csv("%s/smallVars_CNV_output/variant_calling_ploidy2.tab"%outdir, sep="\t")
                s_small_vars_df_ploidy2["sampleID"] = sampleID
                small_vars_df_ploidy2 = small_vars_df_ploidy2.append(s_small_vars_df_ploidy2)

            # add the small variants annotation
            if run_ploidy2_ifHaploid is True: ploidies_small_var_annot = [1,2]
            else: ploidies_small_var_annot = [ploidy]

            for p in ploidies_small_var_annot:

                small_var_annot = small_var_annot.append(pd.read_csv("%s/smallVars_CNV_output/variant_annotation_ploidy%i.tab"%(outdir, p), sep="\t")).drop_duplicates()

            # add the coverage
            s_coverage_df = pd.read_csv("%s/smallVars_CNV_output/CNV_results/genes_and_regions_coverage.tab"%(outdir), sep="\t")
            s_coverage_df["sampleID"] = sampleID
            coverage_df = coverage_df.append(s_coverage_df)

            # add data to 
            sampleID_to_SV_dataDict[sampleID] = {"sampleID":sampleID,
                                                 "SV_CNV_vcf":"%s/SVcalling_output/SV_and_CNV_variant_calling.vcf"%(outdir),
                                                 "SV_CNV_var_annotation":"%s/SVcalling_output/SV_and_CNV_variant_calling.vcf_annotated_VEP.tab"%(outdir)
                                                 }

        # check
        if run_ploidy2_ifHaploid is True:
            if any(pd.isna(small_vars_df.relative_CN)): raise ValueError("there can't be NaNs in small_vars_df")
            if any(pd.isna(small_vars_df_ploidy2.relative_CN)): raise ValueError("there can't be NaNs in small_vars_df_ploidy2")

        if SV_CNV_called is True:

            # get the integrated SV_CNV dfs
            df_data = pd.DataFrame(sampleID_to_SV_dataDict).transpose()
            file_prefix = "%s/integrated_SV_CNV_calling"%cwd

            SV_CNV, SV_CNV_annot = get_integrated_variants_into_one_df(df_data, file_prefix, replace=True, remove_files=True)[2:]


        # add some fields to the annotation dfs
        print_if_verbose("getting Gene and Feature matching the gff")
        if SV_CNV_called is True: SV_CNV_annot = get_annotation_df_with_GeneFeature_as_gff(SV_CNV_annot, gff)
        small_var_annot = get_annotation_df_with_GeneFeature_as_gff(small_var_annot, gff)

        # add whether the SV is protein alterring
        if SV_CNV_called is True:

            gff_df = load_gff3_intoDF(gff, replace=False)
            all_protein_coding_genes = set(gff_df[gff_df.feature.isin({"CDS", "mRNA"})].upmost_parent)
            SV_CNV_annot["is_protein_coding_gene"] = SV_CNV_annot.Gene.isin(all_protein_coding_genes)
            SV_CNV_annot["is_transcript_disrupting"] = SV_CNV_annot.Consequence.apply(get_is_transcript_disrupting_consequence_SV)

        #### ADD THE SAMPLES THAT SHARE VARIANTS ####
        print("adding samples that have the same var")

        # define all small_vars_df
        if run_ploidy2_ifHaploid is False: small_vars_dfs_list = [small_vars_df]
        else: small_vars_dfs_list = [small_vars_df, small_vars_df_ploidy2]

        # add things to the dfs
        for df in small_vars_dfs_list: 

            # add whether it is an heterozygous variant
            if ploidy==1:

                if "relative_CN" in set(df.keys()): df["is_diploid_heterozygous"] = (df.common_GT=="0/1") & (df.relative_CN>=1.5)
                else: df["is_diploid_heterozygous"] = df.common_GT=="0/1"

            elif ploidy==2: df["is_diploid_heterozygous"] = df.common_GT=="0/1"

            # add whether it a haploid
            if ploidy==1: df["is_haploid"] = (df.common_GT=="1") & (df.mean_fractionReadsCov_PASS_algs>=0.9)
            elif ploidy==2: df["is_haploid"] = (df.mean_fractionReadsCov_PASS_algs>=0.9)
            else: raise ValueError("ploidy %i is not acceptable"%ploidy) 

        # merge all the called variants
        merged_small_vars_df = pd.concat(small_vars_dfs_list) 

        # define dfs
        any_called_variants_df = merged_small_vars_df
        PASSatLeast2_variants_df = merged_small_vars_df[(merged_small_vars_df.NPASS>=2) & ((merged_small_vars_df.is_diploid_heterozygous) | (merged_small_vars_df.is_haploid))]

        # map each sample to the variants of each type
        def get_set_variants_from_df_s(df_s): return set(df_s["#Uploaded_variation"])
        sampleID_to_any_called_variants = any_called_variants_df.groupby("sampleID").apply(get_set_variants_from_df_s)
        sampleID_to_PASSatLeast2_variants = PASSatLeast2_variants_df.groupby("sampleID").apply(get_set_variants_from_df_s)

        # create a dict with all sampleID_to_variants_series
        typeVars_to_sampleID_to_variants_series = {"anyCalled":sampleID_to_any_called_variants, "PASSatLeast2":sampleID_to_PASSatLeast2_variants}

        # add the missing samples
        for sampleID_to_variants_series in typeVars_to_sampleID_to_variants_series.values():
            for sampleID in samples_to_run: 
                if sampleID not in set(sampleID_to_variants_series.index): sampleID_to_variants_series[sampleID] = set()

        # add the intersection of each var in each df of small_vars_dfs_list with the types of vars in typeVars_to_sampleID_to_variants_series
        for df in small_vars_dfs_list:
            for typeVars, sampleID_to_variants in typeVars_to_sampleID_to_variants_series.items():

                # init the field of other samples
                df["other_samples_with_%s"%typeVars] = ""

                # map each sample that is mapping 
                for sampleID, variants in sampleID_to_variants.items():
                    bool_to_text = {True: ","+sampleID, False:""}
                    df["other_samples_with_%s"%typeVars] += df["#Uploaded_variation"].isin(variants).map(bool_to_text)

                # get the final set of samples
                def get_set_other_samples_withVar(r): return set(r["other_samples_with_%s"%typeVars].split(",")).difference({"", r["sampleID"]})
                def convert_list_to_string(x): return ",".join(x)

                df["other_samples_with_%s"%typeVars] = df.apply(get_set_other_samples_withVar, axis=1).apply(sorted).apply(convert_list_to_string)

        #############################################

        # write dfs
        save_df_as_tab(small_vars_df, small_vars_df_file)
        save_df_as_tab(small_var_annot, small_var_annot_file)

        if CNV_calling_performed is True: save_df_as_tab(df_CN, integrated_CNperWindow_df)

        if SV_CNV_called is True:
            save_df_as_tab(SV_CNV, SV_CNV_file_simple)
            save_df_as_tab(SV_CNV_annot, SV_CNV_annot_file)

        if run_ploidy2_ifHaploid is True: save_df_as_tab(small_vars_df_ploidy2, small_vars_df_ploidy2_file) 
        save_df_as_tab(coverage_df, coverage_df_file)

    # load the CNV calling df
    if CNV_calling_performed is True: df_CN_all = get_tab_as_df_or_empty_df(integrated_CNperWindow_df)
    else: df_CN_all = None

    ###########################################################

    ####### GET THE INTEGRATED GRIDSS DF AND BEDPE DF ######

    # generates a gridss df that has all the breakpoints

    if (file_is_empty(integrated_gridss_df) or file_is_empty(integrated_bedpe_df)) and SV_CNV_called is True:

        print("generating integrated gridss df")

        samples_to_run = set(paths_df.sampleID)
        #samples_to_run = {"Cg1_CBS138", "Cg2_921192_2"}

        # make a folder to integrate the gridss df
        outdir_integrating_gridss_df = "%s/intergrating_gridss_df"%cwd
        delete_folder(outdir_integrating_gridss_df)
        make_folder(outdir_integrating_gridss_df)

        # init
        df_gridss_all = pd.DataFrame()
        df_bedpe_all = pd.DataFrame()

        for sampleID in samples_to_run:
            print(sampleID)

            # get the outdir
            outdir = "%s/%s_VarCallresults/SVdetection_output/final_gridss_running"%(VarCallOutdirs, sampleID); make_folder(outdir)

            # define the filenames original
            origin_gridss_vcf_raw_file = "%s/gridss_output.raw.vcf"%outdir
            origin_gridss_vcf_filt_file = "%s/gridss_output.filt.vcf"%outdir

            # put them into the outdir_integrating_gridss_df
            gridss_vcf_raw_file = "%s/%s_gridss_output.raw.vcf"%(outdir_integrating_gridss_df, sampleID)
            gridss_vcf_filt_file = "%s/%s_gridss_output.filt.vcf"%(outdir_integrating_gridss_df, sampleID)
            soft_link_files(origin_gridss_vcf_raw_file, gridss_vcf_raw_file)
            soft_link_files(origin_gridss_vcf_filt_file, gridss_vcf_filt_file)

            ## GET THE BEDPE ##

            # get the bedpe files
            bedpe_raw = get_tab_as_df_or_empty_df(get_bedpe_from_svVCF(gridss_vcf_raw_file, outdir_integrating_gridss_df, replace=False, only_simple_conversion=True))
            bedpe_filt = get_tab_as_df_or_empty_df(get_bedpe_from_svVCF(gridss_vcf_filt_file, outdir_integrating_gridss_df, replace=False, only_simple_conversion=True))

            # add whether it is PASS
            pass_breakpoints = set(bedpe_filt.name)
            bedpe_raw["PASSed_filters"] = bedpe_raw.name.isin(pass_breakpoints)
            
            # add name and keep
            bedpe_raw["sampleID"] = sampleID
            df_bedpe_all = df_bedpe_all.append(bedpe_raw)
            ###################

            ## GET THE GRIDSS ##

            # get the gridss vcfs
            gridss_vcf_raw = get_df_and_header_from_vcf(gridss_vcf_raw_file)[0]
            gridss_vcf_filt = get_df_and_header_from_vcf(gridss_vcf_filt_file)[0]

            # change names
            sample_name_vcf = gridss_vcf_raw.columns[-1]
            gridss_vcf_raw = gridss_vcf_raw.rename(columns={sample_name_vcf:"DATA"})
            gridss_vcf_filt = gridss_vcf_filt.rename(columns={sample_name_vcf:"DATA"})

            # add whether it passed the filters
            pass_variants = set(gridss_vcf_filt.ID)
            gridss_vcf_raw["PASSed_filters"] = gridss_vcf_raw.ID.isin(pass_variants)
            gridss_vcf_raw["sampleID"] = sampleID

            # keep
            df_gridss_all = df_gridss_all.append(gridss_vcf_raw)

            ####################

        # clean
        delete_folder(outdir_integrating_gridss_df)

        # save
        save_df_as_tab(df_gridss_all, integrated_gridss_df)
        save_df_as_tab(df_bedpe_all, integrated_bedpe_df)

    # load the df
    df_gridss_all = get_tab_as_df_or_empty_df(integrated_gridss_df)
    df_bedpe_all = get_tab_as_df_or_empty_df(integrated_bedpe_df)

    ###########################################

    ######### GET THE COMMON DF OF SVs #########

    if file_is_empty(SV_CNV_file) and SV_CNV_called is True:

        print("adding the common variant ID")

        # loading SV_CNV simple
        SV_CNV = get_tab_as_df_or_empty_df(SV_CNV_file_simple)

        # add the common variant ID across samples
        outdir_common_variantID_acrossSamples = "%s/getting_common_variantID_acrossSamples"%cwd
        
        # add the overlaps with other samples
        SV_CNV = get_SV_CNV_df_with_overlaps_with_all_samples(SV_CNV, outdir_common_variantID_acrossSamples, tol_bp, pct_overlap, cwd, df_bedpe_all, df_CN_all, threads, reference_genome)

        # get the common variant ID
        SV_CNV = get_SV_CNV_df_with_common_variantID_acrossSamples(SV_CNV, outdir_common_variantID_acrossSamples, pct_overlap, tol_bp)

        # keep relevant files
        SV_CNV = SV_CNV[[k for k in SV_CNV.keys() if k!="INFO"]]

        delete_folder(outdir_common_variantID_acrossSamples)

        # keep
        save_df_as_tab(SV_CNV, SV_CNV_file)

    #############################################



  

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

