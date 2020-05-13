#!/usr/bin/env python

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
import perSVade_graphics_functions as fun

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

# other executables that may be important in the future
"""
VARCALLPIPELINE = "%s/varcall_cnv_pipeline_miki.py"%FunDir
FASTQC = "%s/software/FastQC/fastqc"%mschikoraParentDir 
TRIMMOMATIC = "%s/run_trimmomatic.py"%FunDir 
ensembl_vep = "%s/run_vep.sh"%FunDir
ninja_dir = "%s/software/NINJA-0.95-cluster_only/NINJA"%mschikoraParentDir
"""

### get the executables that should be in the path ###

gridss_run_list = []
gridss_jar_list = []
clove_list = []
ninja_dir_list = []

# executables that should be in the path
for folder in os.environ['PATH'].split(":"):

    # if it does not exist continue
    if not os.path.isdir(folder): continue

    # go through each file in folder
    for file in os.listdir(folder):
        path = "%s/%s"%(folder, file)

        # gridss.sh
        if file=="gridss.sh" or file=="gridss_2.8.1.sh": gridss_run_list.append(path)

        # gridss with dependencies
        if file.startswith("gridss-") and file.endswith("-jar-with-dependencies.jar"): gridss_jar_list.append(path)

        # clove with dependencies
        if file.startswith("clove-") and file.endswith("-jar-with-dependencies.jar"): clove_list.append(path)

        # the ninja directory, it should contain this file
        if file=="NINJA" and os.path.isdir(path) and "ArrayHeapExtMem.cpp" in os.listdir(path): ninja_dir_list.append(path)

# get the latest version
gridss_run = sorted(gridss_run_list)[-1]
gridss_jar = sorted(gridss_jar_list)[-1]
clove = sorted(clove_list)[-1]
ninja_dir = sorted(ninja_dir_list)[-1]

# print the executables
print("%s is the script to use as gridss.sh"%gridss_run)
print("%s is the script to use as gridss .jar"%gridss_jar)
print("%s is the script to use as clove .jar"%clove)
print("%s is the dir of ninja"%ninja_dir)

######################################################


# scripts that are of this pipeline
simulateSV_R = "%s/create_simulatedSVgenome.R"%CWD
annotate_simpleEvents_gridssVCF_R = "%s/annotate_simpleEvents_gridssVCF.R"%CWD
analyze_svVCF = "%s/generate_files_from_svVCF.R"%CWD
analyze_svVCF_simple = "%s/generate_files_from_svVCF_simple.R"%CWD


######################################################
######################################################

# define the strings that have to be considered as NaN in the VCF parsing
vcf_strings_as_NaNs = ['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'n/a', 'nan', 'null']

def chunks(l, n):
    
    """Yield successive n-sized chunks from a list l"""
    
    for i in range(0, len(l), n):
        yield l[i:i + n]

def run_cmd(cmd):

    """Runs a cmd under the VarCall_CNV_env environment, as defined in CONDA_ACTIVATING_CMD"""
    #out_stat = os.system("%s %s"%(CONDA_ACTIVATING_CMD, cmd)); 
    out_stat = os.system(cmd); 
    if out_stat!=0: raise ValueError("\n%s\n did not finish correctly. Out status: %i"%(cmd, out_stat))

def file_is_empty(path): 
    
    """ask if a file is empty or does not exist """
    
    if not os.path.isfile(path):
        return_val = True
    elif os.stat(path).st_size==0:
        return_val = True
    else:
        return_val = False
            
    return return_val

def make_flat_listOflists(LoL):

    return list(itertools.chain.from_iterable(LoL))

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

def get_object_as_str(x):

    """Takes a builtin object and returns the string of it"""

    # get the type
    tx = type(x)

    if tx==str: return x
    elif tx==int: return "%i"%x
    elif tx==float: return "%.2f"%x
    elif tx in {set, list, tuple}: return "-".join(sorted([str(y) for y in x]))
    elif tx==bool: return "%i"%int(x)
    else: raise ValueError("%s has not been considered"%tx)

def get_dir(filename): return "/".join(filename.split("/")[0:-1])
def get_file(filename): return filename.split("/")[-1]


def remove_file(f):

    if os.path.isfile(f): os.unlink(f)

def clean_reference_genome_windows_files(reference_genome):

    """Cleans all the files under reference_genome that are windows files and bed's """

    print("removing windows files")
    ref_dir = get_dir(reference_genome)
    ref_name = get_file(reference_genome)

    for file in os.listdir(ref_dir):
        if file.startswith(ref_name) and "windows" in file and "bp.bed" in file : remove_file("%s/%s"%(ref_dir, file))


def delete_folder(f):

    if os.path.isdir(f): shutil.rmtree(f)

def make_folder(f):

    if not os.path.isdir(f): os.mkdir(f)

def delete_file_or_folder(f):

    """Takes a path and removes it"""

    if os.path.isdir(f): shutil.rmtree(f)
    if os.path.isfile(f): os.unlink(f)


def id_generator(size=10, chars=string.ascii_uppercase + string.digits, already_existing_ids=set()):

    """ already_existing_ids is a set that indicates whihc IDs can't be picked """

    ID = ''.join(random.choice(chars) for _ in range(size))
    while ID in already_existing_ids:
        ID = ''.join(random.choice(chars) for _ in range(size))

    return ID

def save_object(obj, filename):
    
    """ This is for saving python objects """
    
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    
    """ This is for loading python objects """
        
    return pickle.load(open(filename,"rb"))

def extract_BEDofGENES_of_gff3(gff, bed, window_l=10000, replace=False, reference=""):

    """Takes the full path to a gff annotation and writes a bed file to bed, which only has genes' info. 
    The bed will be 1-indexed, which is not usually the case. Bed files are 0 inexed. 

    It also writes a file called 'bed.regions_lengthwindow_l' that has the regions of - and + window_l of the start and end of the gene, respectively.
    The regions bed is returned by this function.

    """

    

    # load the gff
    df_gff3 = pd.read_csv(gff, skiprows=list(range(len([line for line in open(gff, "r") if line.startswith("#")]))), sep="\t", names=["chromosome", "source", "type_feature", "start", "end", "score", "strand", "phase", "attributes"])
    
    # define a function that takes attribues and returns ID
    def get_ID_gff(attributes):

        ID = [x.lstrip("ID=") for x in attributes.split(";") if x.startswith("ID=")][0]
        if ";part=" in attributes: ID += "_part%s"%([x.lstrip("part=").split("/")[0] for x in attributes.split(";") if x.startswith("part=")][0])

        return ID

    df_gff3["ID"] = df_gff3.attributes.apply(get_ID_gff)

    # define the regions you are interested in 
    interesting_features = {"gene"}

    # get the bed of the interesting features
    df_bed = df_gff3[df_gff3.type_feature.isin(interesting_features)][["chromosome", "start", "end", "ID"]]
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

def write_coverage_per_gene(mpileup, bed, gene_to_coverage_file):

    """Takes an mpileup file that has chromosome, position, base, coverage and a bed file with "chromosome", "start", "end", "ID" 
    It writes under gene_to_coverage_file a tab-sepparated with geneID and median_coverage"""

    # load in dataframes
    df_bed = pd.read_csv(bed, sep="\t", names=["chromosome", "start", "end", "ID"]).set_index("ID", drop=False)
    df_mpileup = pd.read_csv(mpileup, sep="\t", names=["chromosome", "position", "base", "n_reads"])

    def get_medianReads_fractionCovered_chr(row, mpileup_df_chr):

        """Takes a row of the df_bed (for a chr) and the mpileup_df_chr. It returns a series with the medianReads, the fraction covered and the """

        # get coords
        start = row["start"]; end = row["end"]

        # get the df of the region
        df_region = mpileup_df_chr[(mpileup_df_chr.position>=start) & (mpileup_df_chr.position<=end)]

        # check that the region is correct and define medianReads:
        if len(df_region)>=(0.95*(end-start+1)): medianReads = np.median(df_region.n_reads)
        else: raise ValueError("The gene %s, between %i-%i is incorrect. Len df is %i. "%(row.name, start, end, len(df_region)))

        # now fractionCovered
        fractionCovered = len(df_region[df_region.n_reads>=2])/len(df_region)

        return pd.Series({"chromosome":row["chromosome"], "start":row["start"], "end":row["end"], "ID":row["ID"], "median_reads_per_gene":medianReads, "fraction_covered_by_MoreThan1read": fractionCovered})

    def get_relativeCN_chr(row, df_coverage_chr):

        """Takes a row of df_coverage_chr and returns the CN: 0, 1 2, 3. df_coverage_chr is necessary to get the nerby genes (-+15). It is iniced equal to row"""

        # first the losses
        if row["fraction_covered_by_MoreThan1read"]<0.1: return 0
        else:
            # get the surrounding 30 genes, unless you are in some of the boundaries

            # get the genes
            genes_5prime = set(list(df_coverage_chr[df_coverage_chr.start<row["start"]].sort_values(by=["start"], ascending=False).ID)[0:15])
            genes_3prime = set(list(df_coverage_chr[df_coverage_chr.start>row["start"]].sort_values(by=["start"], ascending=True).ID)[0:15])

            # get the median of the surrounding genes
            median_surrounding30 = np.median(df_coverage_chr.loc[genes_5prime.union(genes_3prime)].median_reads_per_gene)
            relative_coverage_to_surrounding30 = row["median_reads_per_gene"] / median_surrounding30

            # return
            if relative_coverage_to_surrounding30 > 2.9: return 3
            elif relative_coverage_to_surrounding30 > 1.9: return 2
            else: return 1

    df_coverage = pd.DataFrame()

    # add the median reads and the fraction covered
    for chromosome in set(df_bed.chromosome):
        print("Median reads , fraction covered and relative CN for chromosome %s"%chromosome)

        # get the mpileup of the chromosome
        mpileup_df_chr = df_mpileup[df_mpileup.chromosome==chromosome]
        df_bed_chr = df_bed[df_bed.chromosome==chromosome]

        # initialize with the median coverage and the fraction covered
        df_coverage_chr = df_bed_chr.apply(lambda row: get_medianReads_fractionCovered_chr(row, mpileup_df_chr), axis=1)

        # add the relative copy number
        df_coverage_chr["relative_copy_number"] = df_coverage_chr.apply(lambda row: get_relativeCN_chr(row, df_coverage_chr), axis=1)

        # add 
        df_coverage = df_coverage.append(cp.deepcopy(df_coverage_chr), sort=True)

    # add the relative coverage of each gene
    df_coverage["relative_coverage"] = df_coverage.median_reads_per_gene / np.median(df_coverage.median_reads_per_gene)

    # write
    df_coverage.to_csv(path_or_buf=gene_to_coverage_file, sep="\t", index=False)





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
    print("calculating features for multialleles")

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

    #print(list(series_dicts.apply(lambda d: d["var_to_frequency"])))

    # get the union of all dicts
    print("Geting final dicts")
    var_to_frequency = dict(j for i in series_dicts.apply(lambda d: d["var_to_frequency"]) for j in i.items())
    var_to_filter = dict(j for i in series_dicts.apply(lambda d: d["var_to_filter"]) for j in i.items())
    var_to_GT = dict(j for i in series_dicts.apply(lambda d: d["var_to_GT"]) for j in i.items())
    var_to_filters = dict(j for i in series_dicts.apply(lambda d: d["var_to_filters"]) for j in i.items())

    # get only the columns that you want to keep the real vcf
    return df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",  data_colname]], var_to_frequency, var_to_filter, var_to_GT, var_to_filters

def load_vep_table_intoDF(vep_filename):

    """Takes the tabular output of ensembl vep and returns a table where the index is (chr, pos, ref, alt)"""

    # load normal df
    vep_df = pd.read_csv(vep_filename, sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)

    # keep the orignial uploaded variation
    vep_df["#Uploaded_variation_original"] = vep_df["#Uploaded_variation"]

    # change the uploaded variation so that it includes only the allele position
    vep_df["#Uploaded_variation"] = vep_df["#Uploaded_variation"].apply(lambda x: x.split("/")[0]) + "/" + vep_df.Allele

    # add the index of the alternative allele which corresponds to this var. For example, if the gt is 1|1|1, and this var has a 1, it will mean that it is the gt of the middle.
    print("geting GT index")
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

def generate_jobarray_file(jobs_filename, stderr, stdout, walltime, memory,  name="JobArray", queue="tg-el7", qsub=False, nnodes=1, rmstd=True ):
    
    """ This function takes:
        jobs_filename: a path to a file in which each line is a command that has to be executed in a sepparate cluster node
        name: the name of the jobs array
        stderr and stdout: paths to directories that will contains the STDERR and STDOUT files
        walltime is the time in "hh:mm:ss"
        memory is the RAM: i.e.: 4G
        nnodes is the number of nodes that each job gets
        
        name is the name prefix
        queue is the queue
        rmstd indicates if the previous std has to be removed
        
        It returns a jobs_filename.run file, which can be qsub to the cluster directly if qsub is True
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
    
    # define arguments
    arguments = ["#!/bin/sh",
    "#$ -S /bin/bash",
    "#$ -N %s"%name,
    "#$ -cwd",
    "#$ -q %s"%queue,
    "#$ -o %s"%stdout,
    "#$ -e %s"%stderr,
    "#$ -j n",
    "#$ -l h_rt=%s"%walltime,
    "#$ -l virtual_free=%s"%memory,
    "#$ -pe smp %i"%nnodes,
    "#$ -p 0",
    "#$ -t 1-%i"%n_jobs,
    "",
    '$(sed $SGE_TASK_ID"q;d" %s);'%jobs_filename,
    ""]
    
    # define and write the run filename
    jobs_filename_run = "%s.run"%jobs_filename
    with open(jobs_filename_run, "w") as fd: fd.write("\n".join(arguments))
    
    # run in cluster if specified
    if qsub is True: out_state = os.system("qsub %s"%jobs_filename_run); print("%s Qsub out state: "%name, out_state)

def generate_jobarray_file_slurm(jobs_filename, stderr="./STDERR", stdout="./STDOUT", walltime="02:00:00",  name="JobArray", queue="bsc_ls", sbatch=False, ncores_per_task=1, rmstd=True, constraint="", number_tasks_to_run_at_once="all" ):
    
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
                 "#SBATCH --nodes=1",
                 constraint_line,
                 "",
                 "echo 'sourcing conda to run pipeline...';",
                 "echo 'running pipeline';",
                 'srun $(head -n $SLURM_ARRAY_TASK_ID %s | tail -n 1)'%jobs_filename]


    # this is left SOURCE_CONDA_CMD, CONDA_ACTIVATING_CMD

    """
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
                 "#SBATCH --nodes=1",
                 constraint_line,
                 "",
                 "echo 'sourcing conda to run pipeline...';",
                 "%s;"%SOURCE_CONDA_CMD,                 
                 "echo 'activating conda environment to run pipeline. This will be killed if it takes more than 10 seconds...';",
                 'for i in {1..10}',
                 'do',
                 "echo 'trying to activate the conda env...';",
                 "timeout 10s conda init bash && timeout 10s %s && break"%(CONDA_ACTIVATING_CMD),
                 'done',
                 "echo 'running pipeline';",
                 'srun $(head -n $SLURM_ARRAY_TASK_ID %s | tail -n 1)'%jobs_filename]
    """

    # define and write the run filename
    jobs_filename_run = "%s.run"%jobs_filename
    with open(jobs_filename_run, "w") as fd: fd.write("\n".join(arguments))
    
    # run in cluster if specified
    if sbatch is True: out_state = os.system("sbatch %s"%jobs_filename_run); print("%s sbatch out state: %i"%(name, out_state))

    # get info about the exit status: sacct -j <jobid> --format=JobID,JobName,MaxRSS,Elapsed

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

def string_in_file(string, file): 

    """Return if a string is anywhere in a file"""

    if file_is_empty(file): return False
    else: return string in "".join(open(file, "r").readlines())

def quoute_spaces(string):

    """Takes a string and returns it as '' if it has spaces """
    if " " in string: return "'%s'"%string
    else: return string

def run_qualityControl_and_filtering_reads(paths_df, cwd, repeatFirstFastQC=False, repeatTrimmomatic=False, repeatSecondFastQC=False, n_threads=8):

    """This function inputs a paths_df (or a tab-separated file), which contains an index as 0-N rows and columns "raw_reads_dir", "fastqc_dir", "trimmed_reads_dir", "trimmed_fastqc_dir", "sampleID", "readID"  and runs fastqc trimmomatic and fastqc again in a parallelized manner for each of the rows, without repeating. 

    sampleID is the unique sampleID, which has one FWD and one RV. readID is either R1 or R2.
    trimmed read names should end with .fq.gz

    It assumes that it is paired end sequencing.

    STDERR, STDOUT and the jobs files will be written in cwd. Returns whether the files have all been generated"""

    print("Running quiality control and filtering of reads ...")

    # define the STDERR and STDOUT files
    STDERR = "%s/STDERR"%cwd; STDOUT = "%s/STDOUT"%cwd

    # debug the fact that there are spaces
    for f in ["raw_reads_dir", "fastqc_dir", "trimmed_reads_dir", "trimmed_fastqc_dir", "sampleID"]: paths_df[f] = paths_df[f].apply(quoute_spaces)

    ###### RUN FASTQC FOR THE RAW READS ########
    all_html_files = []
    all_cmds = []

    for I in paths_df.index:

        fastqc_dir = paths_df.loc[I, "fastqc_dir"]
        if not os.path.isdir(fastqc_dir): os.mkdir(fastqc_dir)
        raw_reads_file = paths_df.loc[I, "raw_reads_dir"]

        # define the possinle html reports
        html_files = ["%s/%s"%(fastqc_dir, x) for x in os.listdir(fastqc_dir) if x.endswith(".html")]

        if len(html_files)==0 or repeatFirstFastQC is True:

            # create the outdir of fasqc
            if not os.path.isdir(fastqc_dir): os.mkdir(fastqc_dir)

            # append the cmd
            all_cmds.append("%s -o %s --threads %i --extract --java %s %s"%(FASTQC, fastqc_dir, n_threads, JAVA, raw_reads_file))

        # keep the report
        else: all_html_files.append(html_files[0])

    # write to jobs file and submit to cluster
    if len(all_cmds)>0:
        print("Submitting FastQC to cluster ...")
        jobs_filename = "%s/jobs.run_RawReadsFastqc"%cwd
        open(jobs_filename, "w").write("\n".join(all_cmds))

        generate_jobarray_file_slurm(jobs_filename, stderr=STDERR, stdout=STDOUT, walltime="01:00:00",  name="FirstFastQC", queue="bsc_ls", sbatch=True, ncores_per_task=n_threads, rmstd=True, constraint="", number_tasks_to_run_at_once="all" )
        return False

    # print 
    print("These are the fastQC output results: %s"%("\n".join(all_html_files)))

    # get adapters
    adapters = set.union(*[get_set_adapter_fastqc_report(x) for x in all_html_files])
    print("These are the adapters: \n", adapters)

    # write adapters to fasta in CWD
    all_seqs = []
    existing_ids = set()
    for adapter in adapters:
        ID = id_generator(already_existing_ids=existing_ids); existing_ids.add(ID)
        all_seqs.append(SeqRecord(Seq(adapter), id=ID, name="", description=""))

    adapters_filename = "%s/adapters.fasta"%cwd
    SeqIO.write(all_seqs, adapters_filename, "fasta")

    ###### RUN TRIMMOMATIC FOR THE RAW READS ########

    # go through each sample and save the cmds for trimmomatic
    all_cmds_trim = []
    for sampleID in set(paths_df.sampleID):

        # define the df for this sample
        df = paths_df[paths_df.sampleID==sampleID]
        df1 = df[df.readID=="R1"]
        df2 = df[df.readID=="R2"]

        # define the reads
        raw_reads1 = df1.raw_reads_dir.values[0]; raw_reads2 = df2.raw_reads_dir.values[0]
        trimmed_reads1 = df1.trimmed_reads_dir.values[0]; trimmed_reads2 = df2.trimmed_reads_dir.values[0];

        # check that all trimmed reads end with ".fq.gz"
        for treads in {trimmed_reads1, trimmed_reads2}: 
            if not treads.endswith(".fq.gz"): raise ValueError("Trimmed reads %s is not a valid name, it should end with '.fq.gz' for proper trimmomatic working."%treads)

        # map the repeat to a value
        repeatType_to_val = {True:"--repeat", False:""}

        # run trimmomatic
        if file_is_empty(trimmed_reads1) or file_is_empty(trimmed_reads2) or repeatTrimmomatic is True:
            trim_cmd = "%s --number_threads %i -rr1 %s -rr2 %s -tr1 %s -tr2 %s -ad %s %s"%(TRIMMOMATIC, n_threads, raw_reads1, raw_reads2, trimmed_reads1, trimmed_reads2, adapters_filename, repeatType_to_val[repeatTrimmomatic])
            all_cmds_trim.append(trim_cmd)

    # run trimmomatic in cluster
    if len(all_cmds_trim)>0:
        print("Submitting Trimmomatic to cluster ...")
        jobs_filename = "%s/jobs.run_Trimmomatic"%cwd
        open(jobs_filename, "w").write("\n".join(all_cmds_trim))

        # with 20 cores it takes maximum 45 min in MN
        generate_jobarray_file_slurm(jobs_filename, stderr=STDERR, stdout=STDOUT, walltime="02:00:00",  name="Trimmomatic", queue="debug", sbatch=True, ncores_per_task=n_threads, rmstd=True, constraint="", number_tasks_to_run_at_once="all" )

        return False

    ##### RUN FASTQC FOR THE TRIMMED READS #########

    all_cmds = []
    for I in paths_df.index:

        trimmed_fastqc_dir = paths_df.loc[I, "trimmed_fastqc_dir"]
        if not os.path.isdir(trimmed_fastqc_dir): os.mkdir(trimmed_fastqc_dir)
        trimmed_reads_file = paths_df.loc[I, "trimmed_reads_dir"]

        # create the outdir of fasqc
        if not os.path.isdir(trimmed_fastqc_dir): os.mkdir(trimmed_fastqc_dir)

        # define the possinle html reports
        html_files = ["%s/%s"%(trimmed_fastqc_dir, x) for x in os.listdir(trimmed_fastqc_dir) if x.endswith(".html")]

        if len(html_files)==0 or repeatSecondFastQC is True:

            # append the cmd
            all_cmds.append("%s -o %s --threads %i --extract --java %s %s"%(FASTQC, trimmed_fastqc_dir, n_threads, JAVA, trimmed_reads_file))

    # write to jobs file and submit to cluster
    if len(all_cmds)>0:
        print("Submitting FastQC of the trimmed reads to cluster ...")
        jobs_filename = "%s/jobs.run_TrimmedReadsFastqc"%cwd
        open(jobs_filename, "w").write("\n".join(all_cmds))

        generate_jobarray_file_slurm(jobs_filename, stderr=STDERR, stdout=STDOUT, walltime="01:00:00",  name="TrimmedFastQC", queue="debug", sbatch=True, ncores_per_task=n_threads, rmstd=True, constraint="", number_tasks_to_run_at_once="all" )

        return False

    # only get here if all files have already been generated
    print("You finished all the Read processing and cleaning")
    return True

def convert_NaN_to_empty_strings(serie):
    
    """ Takes a pandas series containing NaN and strings, and converts one that has the NaN into empty "" """
    
    serie[serie.isna()] = [""]*len(serie[serie.isna()])

def convert_to_string(filter_tag):

    if filter_tag=="": return "notPASS"
    else: return filter_tag 



# define a function that takes a sampleID and returns the per-sample df
def get_perSample_df(sampleID, VarCallOutdirs, sampleSpec_fields, ploidy, Uploaded_variation_to_is_snp):

    print(sampleID)

    # get  the outdir
    outdir = "%s/%s_VarCallresults"%(VarCallOutdirs, sampleID)

    # vars
    df_vars = load_object("%s/integrated_variants_norm_vcflib_ploidy%i.py"%(outdir, ploidy))# debug

    # get the interesting fields
    interesting_fields = [x for x in sampleSpec_fields if x in df_vars.columns]
    df_vars = df_vars[interesting_fields].drop_duplicates()
    df_vars["sampleID"] = [sampleID]*len(df_vars)

    # get duplicaed entries
    df_vars = df_vars.set_index("#Uploaded_variation", drop=False)

    # change the bscftools filtertag, which has NaNs instead of empty strings
    df_vars["bcftools_FILTERtag"] = df_vars.bcftools_FILTERtag.apply(lambda x: convert_to_string(x))
    df_vars["freebayes_FILTERtag"] = df_vars.freebayes_FILTERtag.apply(lambda x:convert_to_string(x))
    df_vars["HaplotypeCaller_FILTERtag"] = df_vars.HaplotypeCaller_FILTERtag.apply(lambda x:convert_to_string(x))

    # the same for the GT tag
    df_vars["bcftools_GT"] = df_vars.bcftools_GT.apply(lambda x: convert_to_string(x))
    df_vars["freebayes_GT"] = df_vars.freebayes_GT.apply(lambda x:convert_to_string(x))
    df_vars["HaplotypeCaller_GT"] = df_vars.HaplotypeCaller_GT.apply(lambda x:convert_to_string(x))

    # the sample for the additional_filters
    df_vars["bcftools_additional_filters"] = df_vars.bcftools_additional_filters.apply(lambda x: convert_to_string(x))
    df_vars["freebayes_additional_filters'"] = df_vars.freebayes_additional_filters.apply(lambda x:convert_to_string(x))
    df_vars["HaplotypeCaller_additional_filters'"] = df_vars.HaplotypeCaller_additional_filters.apply(lambda x:convert_to_string(x))

    # add the info from the varSpec_df_important_fields
    df_vars["is_snp"] = Uploaded_variation_to_is_snp.loc[df_vars.index]
    if any(pd.isna(df_vars.is_snp)): 

        print(df_vars.is_snp)


        raise ValueError("There are nans in the snp definition")
    
    # add the other interesting things
    for caller in ["freebayes", "bcftools", "HaplotypeCaller"]: df_vars["%s_PASS_greaterThan0.9reads"%caller] = ((df_vars["%s_PASS"%caller]) & (df_vars["%s_fractionReadsCoveringThisVariant"%caller]>=0.9))

    df_vars["correct_PASS"] = ( ((df_vars.is_snp) & (df_vars["freebayes_PASS_greaterThan0.9reads"]) & (df_vars["bcftools_PASS_greaterThan0.9reads"]) & (df_vars["HaplotypeCaller_PASS_greaterThan0.9reads"])) # snps should be called by the 3 pipelines
                                | (~(df_vars.is_snp) & (df_vars["freebayes_PASS_greaterThan0.9reads"]) & (df_vars["HaplotypeCaller_PASS_greaterThan0.9reads"])) ) # indels should be called by fb and HC

    # add the correct PASS for diploid
    for caller in ["freebayes", "bcftools", "HaplotypeCaller"]: df_vars["%s_PASS_greaterThan0.1reads"%caller] = ((df_vars["%s_PASS"%caller]) & (df_vars["%s_fractionReadsCoveringThisVariant"%caller]>=0.1))

    # add the correct PASS diploid
    df_vars["correct_PASS_diploid"] = ( ((df_vars.is_snp) & (df_vars["freebayes_PASS_greaterThan0.1reads"]) & (df_vars["bcftools_PASS_greaterThan0.1reads"]) & (df_vars["HaplotypeCaller_PASS_greaterThan0.1reads"])) # snps should be called by the 3 pipelines
                                | (~(df_vars.is_snp) & (df_vars["freebayes_PASS_greaterThan0.1reads"]) & (df_vars["HaplotypeCaller_PASS_greaterThan0.1reads"])) ) # indels should be called by fb and HC

    # get the per sample specific df
    sampleSpec_df = cp.deepcopy(df_vars[sampleSpec_fields].set_index("sampleID", drop=False).drop_duplicates(subset=["sampleID", "#Uploaded_variation"], keep="first").sort_values(by=["#Uploaded_variation"]).sort_index())

    print("Size of sampleSpec_df: %.2f MB"%(sys.getsizeof(sampleSpec_df)/1000000))

    del df_vars

    return sampleSpec_df

def run_SNP_CNV_pipeline_various_options(paths_df, cwd, refgenome, gff3, repeat=False, n_threads=16, samples_to_repeat_cnv=set(), samples_to_run=set(), time_for_job="48:00:00", mitochondrial_code=3, mitochondrial_chromosome="mito_C_glabrata_CBS138", gDNA_code=1, ploidy=1, repeat_integration=False, mn_queue="bsc_ls", rerun_all=False, repeat_cnv=False, repeat_gridss=False, run_gridss=False, sampleID_to_parentIDs={}, run_qualimap=True, replace_vep_integration=False, run_in_cluster=False):

    """This function inputs a paths_df, which contains an index as 0-N rows and columns "trimmed_reads_dir", "sampleID", "readID"  and runs the snp pipeline without repeating steps (unless indicated). pths_df can also be a tab-sepparated file were there are this 3 fields. The trimmed_reads_dir has to be the full path to the .fastq file. The sampleID has to be the unique sample identifier and the readID has to be R1 or R2 for paired-end sequencing.

    - cwd is the current working directory, where files will be written
    - refgenome is the full path to the reference genome in .fasta
    - gff3 is a gff with the annotations
    - repeat is a boolean that indicates whether to repeat all the steps of the varcall pipeline
    - n_threads are the number of cores per task allocated. In mn, you can use 48 cores per node. It has not been tested whether more cores can be used per task
    - samples_to_repeat_cnv can be a set of sampleIDs for which we want to specifically rerun the CNV part
    - samples_to_run is a set of samples for which we want to run all the pipeline
    - time_for_job is the time each job gets in DD:HH:MM
    - mitochondrial_code is the NCBI translation code for the mitochondrial genome, which has to be identified through mitochondrial_chromosome (the name of the mitochondrial chromosome). You can find the possible mitochondrial codes in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . For yeast, it is useful http://www.candidagenome.org/help/code_tables.shtml
    - gDNA_code is the NCBU translation code for the nuclear genes. In C. albicans it is 12.
    - ploidy is an integer that indicates the ploidy.
    - repeat_integration is a boolean that indicates whether all the results have to be integrated in one table
    - mn_queue is the queue were you want to throw the pipeline. You can use "bsc_queues" to see which you have. debug is a faster one, but you can only specify 2h 
    - rerun_all is a boolean that indicates whether the pipeline has to be rerun for all smaples (not necessarily repeating any steps)
    - repeat_cnv is a boolean that indicates whether the cnv steps have to be repeated
    - repeat_gridss is a boolean that indicates whether the gridss steps have to be repeated
    - run_gridss is a boolean that indicates whether to run the gridss SV analysis
    - sampleID_to_parentIDs is a dictionary that maps each sampleID to the parent sampleIDs (a set), in a way that there will be a col called parentIDs_with_var, which is a string of ||-joined parent IDs where the variant is also found
    - run_qualimap indicates whether to run bamqc qualimap to assess how good are the the bam files
    - replace_vep_integration will rerun only the parsing and unification of the vep files. This is for debuging purposes.
    - run_in_cluster indicates whether to run in the cluster as a job array. If false it will run each job in a for loop

    This function will generate one MN task for each sampleID and submit an independent job in a job array manner. It will write the STDOUT and STDERR of each job under STDERR and STDOUT of cwd.

    For a paired-end dataset of C. glabrata and ~600x coverage each sample-specific dataframe is arround 100MB

    """

    print("Running VarCall pipeline...")

    # if it is a path
    if type(paths_df)==str: paths_df = pd.read_csv(paths_df, sep="\t")

    # create files that are necessary
    STDERR = "%s/STDERR"%cwd
    STDOUT = "%s/STDOUT"%cwd
    VarCallOutdirs = "%s/VarCallOutdirs"%cwd
    for x in {STDERR, STDOUT, VarCallOutdirs}: 
        if not os.path.isdir(x): os.mkdir(x)
    
    # define the samples_to_run
    if len(samples_to_run)==0: samples_to_run = set(paths_df.sampleID)

    # get the info of all the reads and samples
    all_cmds = set()

    for sampleID in samples_to_run:

        # define the df for this sample
        df = paths_df[paths_df.sampleID==sampleID]
        df1 = df[df.readID=="R1"]
        df2 = df[df.readID=="R2"]

        # define the reads of interest and keep
        trimmed_reads1 = df1.trimmed_reads_dir.values[0]; trimmed_reads2 = df2.trimmed_reads_dir.values[0]

        # create an outdir
        outdir = "%s/%s_VarCallresults"%(VarCallOutdirs, sampleID)
        if not os.path.isdir(outdir): os.mkdir(outdir)

        # add the --replace optin if indicated
        if repeat is True: repeat_cmd = "--replace"
        else: repeat_cmd = ""

        # WARNING: IT HAS TO BE ADDED THAT THE GATK NORMALISATION WILL BE ONLY PERFORMED IF THERE ARE NO AMBIGUOUS NUCLEOTIDES IN THE GENOME. THIS HAS TO BE IMPLEMENTED HERE. ALSO CONDITION THE NEED OF THE GRIDSS OUTPUT TO THE FACT THAT IT HAS BEEN CALLED

        # define the files that shoud be not empty in order not to run this code
        success_files = ["%s/alternative_genome_vcflib_ploidy%i.fasta"%(outdir, ploidy),
                         "%s/finsihed_file_ploidy_good%i.txt"%(outdir, ploidy)]
                   
        # define the cmd          
        cmd = "%s -r %s -f1 %s -f2 %s -thr %i -caller all -o %s -p %i -glm BOTH -gff %s -cnv -mchr %s -mcode %i -gcode %i %s"%(VARCALLPIPELINE, refgenome, trimmed_reads1, trimmed_reads2, n_threads, outdir, ploidy, gff3, mitochondrial_chromosome, mitochondrial_code, gDNA_code, repeat_cmd)

        # specific things to repeat
        if sampleID in samples_to_repeat_cnv or repeat_cnv is True: cmd = "%s --replace_cnv"%cmd

        # run gridss
        if run_gridss is True: 
            success_files.append("%s/gridss_output/gridss_finished_file_final_gridss_final_with_window_size.txt"%outdir)
            success_files.append("%s/gridss_output/benchmarking_all_filters_for_all_genomes_and_ploidies/df_cross_benchmark.py"%outdir)

            cmd = "%s -sv_gridss"%cmd

        # repeat the gridss
        if repeat_gridss is True: cmd = "%s --replace_gridss"%cmd

        # qualimap
        if run_qualimap is True: cmd = "%s --run_qualimap"%cmd

        # repeat vep intergration
        if replace_vep_integration is True: cmd = "%s --replace_vep_integration"%cmd

        # run repeating
        #if repeat_gridss is True or repeat_cnv is True: all_cmds.add(cmd)
        if replace_vep_integration is True: all_cmds.add(cmd) # this would be the appropriate way

        # add cmd if necessary
        if any([file_is_empty(x) for x in success_files]) or repeat is True or rerun_all is True: all_cmds.add(cmd)

    all_cmds = list(all_cmds)

    # submit to cluster or return True
    if len(all_cmds)>0:
        print("Submitting %i jobs to cluster ..."%len(all_cmds))
        jobs_filename = "%s/jobs.run_SNPs_CNV"%cwd
        open(jobs_filename, "w").write("\n".join(all_cmds))

        if run_in_cluster is True:
            generate_jobarray_file_slurm(jobs_filename, stderr=STDERR, stdout=STDOUT, walltime=time_for_job,  name="SNPsCNV", queue=mn_queue, sbatch=True, ncores_per_task=n_threads, rmstd=True, constraint="", number_tasks_to_run_at_once="all" )

        else: 
            print("running %i jobs"%len(all_cmds))
            for Icmd, cmd in enumerate(all_cmds): 
                print("running job %i"%(Icmd+1)); run_cmd("%s > %s.local_running_jobs.std 2>&1"%(cmd, jobs_filename)); print("cmd finished\n-------------\n\n\n\n\n")


        return False

    # clean the ref genome for gridss concurrence
    clean_reference_genome_windows_files(refgenome)

    print("Integrating all variants and CNV into one......")

    # integrate all the results if you are done
    integrated_cnv_filepath = "%s/integrated_CNV_genes_and_regions.tab"%cwd
    integrated_variants_perSample_filepath = "%s/perSample_variants_ploidy%i.tab"%(cwd, ploidy)
    integrated_variants_perSample_filepath_with_parents = "%s.with_parent_info.tab"%(integrated_variants_perSample_filepath)
    integrated_variantInfo_filepath = "%s/variantInfo_ploidy%i.tab"%(cwd, ploidy)

    ######### vars df ##########
    print("checking varInfo")
    if file_is_empty(integrated_variantInfo_filepath) or repeat_integration is True:
        print("getting varInfo")

        # initialize dict
        varSpec_df = pd.DataFrame()

        # initialize already saved vars
        already_saved_variation = set()

        # define fields that are var specific
        varSpec_fields = ['#Uploaded_variation', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'chromosome', 'position', 'ref', 'alt', 'is_snp', 'var', 'is_protein_altering']

        for I, sampleID in enumerate(samples_to_run):
            print("perVar_df", I, sampleID)
            outdir = "%s/%s_VarCallresults"%(VarCallOutdirs, sampleID)

            # load df with only the necessary fields
            df_vars_all = load_object("%s/integrated_variants_norm_vcflib_ploidy%i.py"%(outdir, ploidy))
            df_vars_all = df_vars_all[list(set(df_vars_all.columns).intersection(set(varSpec_fields)))].set_index("#Uploaded_variation", drop=False)
            
            # define vars that were not previously added
            new_vars = set(df_vars_all.index).difference(already_saved_variation)
            df_vars = cp.deepcopy(df_vars_all.loc[new_vars]).drop_duplicates()
            already_saved_variation.update(set(df_vars.index))
            print("There is a varsDF that has %i vars"%len(df_vars))

            # delete big df
            del df_vars_all

            # add the sample ID
            df_vars["sampleID"] = [sampleID]*len(df_vars)

            # add if it's a snp
            df_vars["ref"] = df_vars["#Uploaded_variation"].apply(lambda x: x.split("_")[-1].split("/")[0])
            df_vars["alt"] = df_vars["#Uploaded_variation"].apply(lambda x: x.split("_")[-1].split("/")[1])
            #df_vars['is_snp'] = ((df_vars["ref"].isin({"A", "C", "G", "T" , "N"})) & df_vars["alt"].isin({"A", "C", "G", "T" , "N"}))
            df_vars['is_snp'] = (df_vars["ref"].apply(len)==1) & (df_vars["ref"]!="-") & (df_vars["alt"].apply(len)==1) & (df_vars["alt"]!="-")

            if any(pd.isna(df_vars["is_snp"])): 
                print("WARNING... NaNs detected in is_snp.")
                return df_vars

            # add the var
            df_vars["var"] = df_vars.apply(lambda r: (r["chromosome"], r["position"], r["ref"], r["alt"]), axis=1)

            # get the consequences set
            df_vars["consequences_set"] = df_vars.Consequence.apply(lambda x: set(str(x).split(",")))

            # add whether it is protein altering
            prot_altering_mutations = {'missense_variant', 'start_lost', 'inframe_deletion', 'protein_altering_variant', 'stop_gained', 'inframe_insertion', 'frameshift_variant', 'stop_lost', 'splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant', 'non_coding_transcript_exon_variant'}
            df_vars["is_protein_altering"] = df_vars.consequences_set.apply(lambda x: len(x.intersection(prot_altering_mutations))>0)

            # per var
            varSpec_df = varSpec_df.append(df_vars[varSpec_fields]).sort_index()
            print("Size of varSpec_df: %.2f MB"%(sys.getsizeof(varSpec_df)/1000000))

        # save
        varSpec_df.to_csv(integrated_variantInfo_filepath, sep="\t", index=False , header=True)

    else: 
        pass
        #varSpec_df = pd.read_csv(integrated_variantInfo_filepath, sep="\t").set_index("#Uploaded_variation", drop=False)

    ###########################

    ##### per sampleDF ########
    print("checking per-sample df")
    if file_is_empty(integrated_variants_perSample_filepath) or repeat_integration is True:
        print("generating per-sample df")


        # generate df that has important fields for the sample df
        Uploaded_variation_to_is_snp = varSpec_df[["is_snp", "#Uploaded_variation"]].drop_duplicates(subset = ["#Uploaded_variation"])["is_snp"]


        # fields that are sample specific
        sampleSpec_fields = ['#Uploaded_variation', 'HaplotypeCaller_called', 'HaplotypeCaller_FILTERtag', 'HaplotypeCaller_fractionReadsCoveringThisVariant', 'HaplotypeCaller_PASS', 'freebayes_called', 'freebayes_FILTERtag', 'freebayes_fractionReadsCoveringThisVariant', 'freebayes_PASS', 'bcftools_called', 'bcftools_FILTERtag', 'bcftools_fractionReadsCoveringThisVariant', 'bcftools_PASS','sampleID', 'freebayes_PASS_greaterThan0.9reads', 'bcftools_PASS_greaterThan0.9reads', 'HaplotypeCaller_PASS_greaterThan0.9reads', 'freebayes_PASS_greaterThan0.1reads', 'bcftools_PASS_greaterThan0.1reads', 'HaplotypeCaller_PASS_greaterThan0.1reads' , 'HaplotypeCaller_GT', 'bcftools_GT', 'freebayes_GT', 'HaplotypeCaller_additional_filters', 'bcftools_additional_filters', 'freebayes_additional_filters','is_snp', 'correct_PASS', 'correct_PASS_diploid', "var", "HaplotypeCaller_GT_index", "bcftools_GT_index", "freebayes_GT_index", "HaplotypeCaller_#Uploaded_variation_original", "bcftools_#Uploaded_variation_original", "freebayes_#Uploaded_variation_original"]


        inputs_parallel = cp.deepcopy([(sample, VarCallOutdirs, sampleSpec_fields, ploidy, Uploaded_variation_to_is_snp) for sample in samples_to_run])
              
        # get them in parallel
        
        """
        with multiproc.Pool(multiproc.cpu_count()) as pool:

            list_dfs = pool.starmap(get_perSample_df, inputs_parallel)
            pool.join()
            pool.close()
            pool.terminate()
        """

        # in a map
        list_dfs = map(lambda x: get_perSample_df(x[0], x[1], x[2], x[3], x[4]), inputs_parallel)        
        sampleSpec_df = pd.concat(list_dfs)

        # save
        print("saving")
        sampleSpec_df.to_csv(integrated_variants_perSample_filepath, sep="\t", index=False , header=True)

    else: 
        pass
        #sampleSpec_df = pd.read_csv(integrated_variants_perSample_filepath, sep="\t")

    print("checking integrated_variants_perSample_filepath_with_parents")
    if file_is_empty(integrated_variants_perSample_filepath_with_parents) or repeat_integration is True: 
        print("getting integrated_variants_perSample_filepath_with_parents")

        ##### ADD THE PARENTS THAT HAVE THE SAME VARIATION ######

        if len(sampleID_to_parentIDs)>0:
            print("getting parents with vars")

            sampleSpec_df_withParentInfo = pd.DataFrame()

            # calculate all the vars in the samples
            sample_to_vars = {s : set(sampleSpec_df[sampleSpec_df.sampleID==s]["var"]) for s in samples_to_run}
            print("Size of sample_to_vars: %.2f MB"%(sys.getsizeof(sample_to_vars)/1000000))

            # go through each sample
            for i, sample in enumerate(samples_to_run):
                print(i, sample)
    
                # get df
                df_s = sampleSpec_df[sampleSpec_df.sampleID==sample]

                # get the parent IDs
                df_s["parentIDs_with_var"] = df_s["var"].apply(lambda v:  "||".join([p for p in sampleID_to_parentIDs[sample] if v in sample_to_vars[p]]) )

                # keep 
                sampleSpec_df_withParentInfo = sampleSpec_df_withParentInfo.append(df_s)

            # keep 
            sampleSpec_df = sampleSpec_df_withParentInfo
        
        ##########################################################

        # save tables
        sampleSpec_df.to_csv(integrated_variants_perSample_filepath_with_parents, sep="\t", index=False , header=True)


    ##### get the CNV df #####
    print("checking df_cnv")
    if file_is_empty(integrated_cnv_filepath) or repeat_integration is True:
        print("getting df_cnv")

        # initialize df
        df_cnv_all = pd.DataFrame()

        for I, sampleID in enumerate(samples_to_run):
            print("cnv", I, sampleID)
            outdir = "%s/%s_VarCallresults"%(VarCallOutdirs, sampleID)

            # load df 
            df_cnv = pd.read_csv("%s/CNV_results/genes_and_regions_coverage.tab"%outdir, sep="\t")
            df_cnv["sampleID"] = [sampleID]*len(df_cnv)
            
            # keep all dataframes
            df_cnv_all = df_cnv_all.append(df_cnv)
            print("Size of df_cnv_all: %.2f MB"%(sys.getsizeof(df_cnv_all)/1000000))

        # save
        df_cnv_all.to_csv(integrated_cnv_filepath, sep="\t", index=False , header=True)

    #########################

    ######### INTEGRATE GRIDSS OUTPUT INTO ONE #########

    ##### generate a sampleID_to_svtype_to_svDF, which has the variants that are not in the parents (non of the breakends overlaps) and is saved as a python dict ######   


    # get the sampleID_to_svtype_to_file and sampleID_to_dfGRIDSS
    sampleID_to_svtype_to_file_filepath = "%s/sampleID_to_svtype_to_file.py"%cwd
    sampleID_to_dfGRIDSS_filepath = "%s/sampleID_to_dfGRIDSS.py"%cwd

    print("checking per sample dicts")
    if file_is_empty(sampleID_to_svtype_to_file_filepath) or file_is_empty(sampleID_to_dfGRIDSS_filepath) or repeat_integration is True:
        print("getting per sample dictionaries")

        sampleID_to_svtype_to_file = {}
        sampleID_to_dfGRIDSS = {}

        for sampleID in samples_to_run:
            print(sampleID)

            # define the outdir of the sample
            outdir_gridss = "%s/%s_VarCallresults/gridss_output/final_gridss_running_with_optimum_parameters"%(VarCallOutdirs, sampleID)

            # calculate the insert size metrics
            bamfile = "%s/%s_VarCallresults/aligned_reads.bam.sorted"%(VarCallOutdirs, sampleID)
            median_insert_size, median_insert_size_sd  = get_insert_size_distribution(bamfile, replace=False, threads=4)

            # get the gridss dataframe
            gridss_output_vcf = "%s/gridss_output.vcf.withSimpleEventType.vcf"%outdir_gridss
            gridss_df = add_info_to_gridssDF(load_single_sample_VCF(gridss_output_vcf), median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd)
            #("\n".join(list(gridss_df.keys())))

            # add to df
            sampleID_to_dfGRIDSS[sampleID] = gridss_df

            # get through each svtype
            for svtype in {"deletions", "tandemDuplications", "inversions", "remaining", "translocations", "insertions"}:

                # find the file that correct file
                possible_files = [f for f in os.listdir(outdir_gridss) if "structural_variants.%s"%svtype in f]

                # get
                if len(possible_files)==0: svfile = ""
                elif len(possible_files)==1: svfile = "%s/%s"%(outdir_gridss, possible_files[0])
                else: raise ValueError("there are problems with %s"%svtype) 

                sampleID_to_svtype_to_file.setdefault(sampleID, {}).setdefault(svtype, svfile)

        # save
        save_object(sampleID_to_svtype_to_file, sampleID_to_svtype_to_file_filepath)
        save_object(sampleID_to_dfGRIDSS, sampleID_to_dfGRIDSS_filepath)

    else:
        sampleID_to_svtype_to_file = load_object(sampleID_to_svtype_to_file_filepath)
        sampleID_to_dfGRIDSS = load_object(sampleID_to_dfGRIDSS_filepath)



    # go through parent removal and no-parent removal
    for type_svs, sample_to_parents_for_SV in [("removing_parentsSVs", sampleID_to_parentIDs), ("allSVs", {s:set() for s in samples_to_run})]:
        print(type_svs)
        sampleID_to_svtype_to_svDF_file = "%s/sampleID_to_svtype_to_svDF_%s.py"%(cwd, type_svs)

        if file_is_empty(sampleID_to_svtype_to_svDF_file) or repeat_integration is True:
            print("getting sampleID_to_svtype_to_svDF")
        

            # get the sample to parents that are in samples to run
            #sample_to_parents_for_SV = {s: sample_to_parents_for_SV[s].intersection(samples_to_run) for s in samples_to_run}  # debug only running samples
            sample_to_parents_for_SV = {s: sample_to_parents_for_SV[s] for s in samples_to_run}


            # get the sampleID_to_svtype_to_svDF
            sampleID_to_svtype_to_svDF = graph_fun.get_sampleID_to_svtype_to_svDF_filtered(sampleID_to_svtype_to_file, sampleID_to_dfGRIDSS, sampleID_to_parentIDs=sample_to_parents_for_SV, breakend_info_to_keep=['#CHROM', 'POS', 'other_coordinates', 'allele_frequency', 'allele_frequency_SmallEvent', 'real_AF', 'FILTER', 'inserted_sequence', 'length_inexactHomology', 'length_microHomology'])

            # save
            save_object(sampleID_to_svtype_to_svDF, sampleID_to_svtype_to_svDF_file)


    ###########################################################################

    ######################### WINDOWS COVERAGE DF #############################
    for window_size in [50, 100, 500, 1000, 5000, 10000, 20000]:
        windows_coverage_df_file = "%s/coverage_per_%ibp_windows_df.tab"%(cwd, window_size)

        if file_is_empty(windows_coverage_df_file) or repeat_integration is True:
            print("Integrating windows coverage df for windows %ibp into one single table"%window_size)

            # initialize df
            df_coverage_all = pd.DataFrame()

            for I, sampleID in enumerate(samples_to_run):
                print("cnv", I, sampleID)
                outdir = "%s/%s_VarCallresults/gridss_output/final_gridss_running_with_optimum_parameters"%(VarCallOutdirs, sampleID)

                # load df 
                df_cnv = pd.read_csv("%s/coverage_windows_%ibp.tab"%(outdir, window_size), sep="\t")[['#chrom', 'end', 'mediancov_1','percentcovered_1', 'start']] # only important fields
                df_cnv["sampleID"] = [sampleID]*len(df_cnv)

                # keep all dataframes
                df_coverage_all = df_coverage_all.append(df_cnv)
                print("Size of df_coverage_all: %.2f MB"%(sys.getsizeof(df_coverage_all)/1000000))

            # sort values
            df_coverage_all = df_coverage_all.sort_values(by=["sampleID", "#chrom", "start", "end"])

            # save
            df_coverage_all.to_csv(windows_coverage_df_file, sep="\t", index=False , header=True)
            del df_coverage_all

    ###########################################################################






    print("getting SV sampleID_to_svtype_to_svDF")




    # at the end return true to state that everything was correct
    return True

def run_SNP_CNV_pipeline_various_options_haploid(paths_df, cwd, refgenome, gff3, repeat=False, n_threads=16, samples_to_repeat_cnv=set(), samples_to_run=set(), time_for_job="48:00:00", mitochondrial_code=3, mitochondrial_chromosome="mito_C_glabrata_CBS138", repeat_integration=False):

    # run the haploid version. This is for consistency with previous runs
    run_SNP_CNV_pipeline_various_options(paths_df, cwd, refgenome, gff3, repeat=repeat, n_threads=n_threads, samples_to_repeat_cnv=samples_to_repeat_cnv, samples_to_run=samples_to_run, time_for_job=time_for_job, mitochondrial_code=mitochondrial_code, mitochondrial_chromosome=mitochondrial_chromosome, ploidy=1, repeat_integration=repeat_integration)

def get_VEPoutput_filt(input_f):

    """Gets the VEP variants filtered"""

    real_consequences =  {'stop_gained', 'start_retained_variant', 'missense_variant', 'splice_region_variant', 'protein_altering_variant', 'start_lost', 'inframe_deletion', 'stop_retained_variant', 'stop_lost', 'frameshift_variant', 'inframe_insertion', 'splice_acceptor_variant', 'intergenic_variant', '5_prime_UTR_variant', 'downstream_gene_variant', '3_prime_UTR_variant', 'synonymous_variant', 'upstream_gene_variant', 'intron_variant'}

    # load df
    df = pd.read_csv(input_f, sep="\t", header=len([x for x in open(input_f, "r") if x.startswith("##")]), na_values=vcf_strings_as_NaNs, keep_default_na=False)

    # get the real consequences
    df["consequences_set"] = df.Consequence.apply(lambda x: set(x.split(",")))

    # filter out non important consequences
    return df[df.consequences_set.apply(lambda x: len(x.intersection(real_consequences))>0)].set_index("#Uploaded_variation", drop=False)

def load_gtf(gtf_file):  

    """Loads a gtf into a df"""

    return pd.read_csv(gtf_file, skiprows=list(range(len([line for line in open(gtf_file, "r") if line.startswith("#")]))), sep="\t", names=["chromosome", "source", "feature", "start", "end", "blank1", "strand", "blank2", "attribute"])

def get_gtf_readyForSIFT(gtf_file, originalChrName_to_changedChrName):

    """Takes a gff file, convert it to gff and change the annotations file to include the biotype (the 9th column (attribute column) says gene_biotype "protein_coding;" for rows which are labelled as exon, CDS, stop_codon, and start_codon). It also changes the chromosomal names"""

    # load df
    df = load_gtf(gtf_file)

    # check if all the items are the ones expected for sift (it needs, CDS, exon, start_codon, stop_codon)
    all_feats = set(df.feature)
    necessary_feats = {"exon", "CDS", "stop_codon", "start_codon"}
    missing_feats = necessary_feats.difference(all_feats)

    # when something is missing:
    if len(missing_feats)>0: print("Warning: the gtf is missing: ", missing_feats)

    # add the biotype
    def get_attribute_with_biotype(row):

        #print(row["attribute"])

        if row["feature"] in {"exon", "CDS", "stop_codon", "start_codon"}: return '%s gene_biotype \"protein_coding\";'%str(row["attribute"])
        else: return str(row["attribute"])

    df["attribute"] = df.apply(get_attribute_with_biotype, axis=1)

    # change the mitochondrial chromosome
    df["chromosome"] = df.chromosome.apply(lambda c: originalChrName_to_changedChrName[c])

    # write
    fd = open(gtf_file, "w")
    for chromosome, source, feature, start, end, blank1, strand, blank2, attribute in df.values:
        fd.write("%s\n"%("\t".join([str(chromosome), source, feature, str(start), str(end), blank1, strand, blank2, attribute])))

    fd.close()

def make_flat_listOflists(LoL):

    return list(itertools.chain.from_iterable(LoL))

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

def load_single_sample_gridss_vcf(path):

    """Loads the gridss vcf"""

    # load the vcf
    df = load_single_sample_VCF(path)

    # add some tags
    df["has_pairBDN"] = ~pd.isna(df.INFO_PARID)
    df["has_dotInALT"] = df.ALT.apply(lambda x: "." in x)
    df["has_moreThanOne_ALTallele"] = df.ALT.apply(lambda x: "," in x)

    # debug the fact that there are diploidies
    if sum(df.has_moreThanOne_ALTallele)>0: raise ValueError("There are some breakpoints with more than one position")

    return df

def write_clove_output_filteredByCoverage(output_vcf_clove, output_vcf_clove_filt, df_cov, bedpe_filename, gridss_VCFoutput, mitochondrial_chromosome="mito_C_glabrata_CBS138"):

    """This function takes the output of clove and a destination file and writes it so that the <DEL> amd <TAN> are normalised by the region.
    bedpe_filename is a file that has the bedpe file with all the individual breakpoints """

    print("Normalising CLOVE output")
    # load vcf into df
    df = get_clove_output(output_vcf_clove)

    # rename the filter tag
    df = df.rename(columns={"FILTER":"cloveFILTER"})

    # calculate parms
    chrom_to_medianCov = {c : np.median(df_cov.loc[c, "coverage"]) for c in set(df_cov.index)}
    medianCov_all = np.median(df_cov.coverage)
    medianCov_nuclear = np.median(df_cov.loc[set(df_cov.index).difference({mitochondrial_chromosome}), "coverage"])
    chrom_to_dfCov = {c : df_cov.loc[c] for c in set(df_cov.index)}

    # filter by coverage DEL and TAN.
    def get_relative_coverages(row):

        """Takes a row of a df and returns a series of relative coverages to the windows. ADP has to be greater than 0"""

        ADP = float(row["ADP"])

        # only if there is some coverage to be considered
        if ADP>0:

            # parse
            POS = int(row["POS"]); END = int(row["END"])
            chrom = row["#CHROM"]; start = min(POS, END); end = max(POS, END);

            # get the df of the Cov
            df_chr_cov = chrom_to_dfCov[chrom]

            # get the median coverage of the region
            median_coverage_region = np.median(df_chr_cov[(df_chr_cov.position>=start) & (df_chr_cov.position<=end)].coverage)
            
            # initialize a dictionary that will contain relative coverages
            relCoverages_dict = {}

            # add the relative coverages to all or chromosomes
            relCoverages_dict["cov_relative_to_all"] = median_coverage_region/medianCov_all
            relCoverages_dict["cov_relative_to_chrom"] = median_coverage_region/chrom_to_medianCov[chrom]
            relCoverages_dict["cov_relative_to_nuclear"] = median_coverage_region/medianCov_nuclear

            # map each type of normalisation to a function
            typeNorm_to_function = {"median": np.median} #, "mean": np.mean}

            for typeNorm, norm_fn in typeNorm_to_function.items():

                # add the coverages relative to the window
                for window_kb in [0, 1, 3, 5, 10, 100, 1000, 2000]:

                    # find the df 
                    df_cov_window = df_chr_cov[(df_chr_cov.position>=(start-(window_kb*1000))) & (df_chr_cov.position<=(end+(window_kb*1000)))]

                    # normalize
                    relCoverages_dict["cov_relative_to_%s_+-%ikb"%(typeNorm, window_kb)] = median_coverage_region / norm_fn(df_cov_window.coverage)

            return relCoverages_dict

        else: return {}

    print("Getting relative coverages")

    df["relCoverages_dict"] = df.apply(get_relative_coverages, axis=1)
    types_relative_coverage = set(df[df.relCoverages_dict.apply(lambda x: len(x)>0)].relCoverages_dict.iloc[0])

    # fill empty ones
    def fill_empty_ones(covDict):

        if len(covDict)>0: return covDict
        else: return {x : np.nan for x in types_relative_coverage}

    df["relCoverages_dict"] = df.relCoverages_dict.apply(fill_empty_ones)

    # add as dict:
    for relCovType in types_relative_coverage: df[relCovType] = df.relCoverages_dict.apply(lambda x: x[relCovType])
    del df["relCoverages_dict"]

    print("Getting bedpe info")
    # add the TAGs of the bedpe file
    df_bedpe = pd.read_csv(bedpe_filename, sep="\t").set_index("name")

    # define a function that takes the ID of the bedpe and returns the set of breakpoint IDs
    def get_br_IDs(cloveID): return re.split("\+|\-", cloveID)

    def get_item_for_bedpe_df(ID, field):

        if ID in df_bedpe.index: return df_bedpe.loc[ID, field]
        else: return "NotFound"


    # add the filters
    df["cloveIDs"] = df.ID.apply(get_br_IDs)
    df["gridssFILTERS"] = df.cloveIDs.apply(lambda IDs: "--".join([get_item_for_bedpe_df(ID, "FILTERS") for ID in IDs]))
    df["gridssALTs"] = df.cloveIDs.apply(lambda IDs: "--".join([get_item_for_bedpe_df(ID, "ALTs") for ID in IDs]))

    # add the allele frequency from the vcf raw
    df_vcf = load_single_sample_gridss_vcf(gridss_VCFoutput).set_index("ID", drop=False)
    df_vcf["allele_frequency"] = df_vcf.apply(lambda r: np.divide(r["DATA_VF"] , (r["DATA_VF"] + r["DATA_REF"] + r["DATA_REFPAIR"])), axis=1)
    df_vcf["allele_frequency_SmallEvent"] = df_vcf.apply(lambda r: np.divide(r["DATA_VF"] , (r["DATA_VF"] + r["DATA_REF"])), axis=1)

    # add the allele frequencies for each of the breakends to the clove output
    def get_field_from_vcf(ID, target="allele_frequency", type_var="float"):

        """Takes an ID and returns a string with the target field of the vcf ID and the partner"""

        # get partner
        partnerID = df_vcf.loc[ID, "INFO_PARID"]

        # return frequencies
        if type_var=="float": return "%.4f||%.4f"%(df_vcf.loc[ID, target], df_vcf.loc[partnerID, target])
        else: return "%s||%s"%(df_vcf.loc[ID, target], df_vcf.loc[partnerID, target])

    df["gridss_allele_frequencies"] = df.cloveIDs.apply(lambda IDs: "--".join([get_field_from_vcf(ID, target="allele_frequency") for ID in IDs]))
    df["gridss_allele_frequencies_SmallEvent"] = df.cloveIDs.apply(lambda IDs: "--".join([get_field_from_vcf(ID, target="allele_frequency_SmallEvent") for ID in IDs]))

    df["gridss_mean_allele_frequencies"] = df.cloveIDs.apply(lambda IDs: np.mean([np.mean([float(x) for x in get_field_from_vcf(ID, target="allele_frequency").split("||")]) for ID in IDs]))
    df["gridss_mean_allele_frequencies_SmallEvent"] = df.cloveIDs.apply(lambda IDs: np.mean([np.mean([float(x) for x in get_field_from_vcf(ID, target="allele_frequency_SmallEvent").split("||")]) for ID in IDs]))

    # get the partner IDs
    df["gridss_IDs"] = df.cloveIDs.apply(lambda IDs: "--".join([get_field_from_vcf(ID, target="ID", type_var="string") for ID in IDs]))


    # WRITE 
    print("writing results")
    df.to_csv(output_vcf_clove_filt, sep="\t", index=False, header=True)

    return df

def get_alternative_genome_from_intersetcion_VCFs(cwd, vcf_iterable, reference_genome, alternative_genome, threads=10, replace=False):

    """Takes an iterable of vcfs and takes the intersection between them to be changing the reference_genome into the alternative_genome. Etra files are written to cwd"""

    if file_is_empty(alternative_genome) or replace is True:

        # generate the bgzipped and indexed vcfs only with PASS
        all_correctlyFormatted_vcfs = []
        for vcf in vcf_iterable: 
            print("working on", vcf)
            
            # sort and get the only pass
            onlyPASS_vcf = "%s.onlyPASS.vcf"%vcf
            run_cmd("grep 'PASS\|^#' %s > %s; "%(vcf, onlyPASS_vcf))

            # sort with bedtools
            sorted_vcf = "%s.sorted"%onlyPASS_vcf
            run_cmd("%s sort -header -i %s > %s"%(bedtools, onlyPASS_vcf, sorted_vcf))
            os.rename(sorted_vcf, onlyPASS_vcf)

            # add the header
            onlyPASS_withHeader_vcf = "%s.withHeader.vcf"%onlyPASS_vcf
            open(onlyPASS_withHeader_vcf, "w").write("##fileformat=VCFv4.2\n%s"%("".join(open(onlyPASS_vcf, "r").readlines())))

            # create the .gz files
            onlyPASS_withHeader_vcf_gz = "%s.gz"%onlyPASS_withHeader_vcf
            run_cmd("%s -c %s > %s; %s -p vcf %s"%(bgzip, onlyPASS_withHeader_vcf, onlyPASS_withHeader_vcf_gz, tabix, onlyPASS_withHeader_vcf_gz))

            # remove and unlink
            os.unlink(onlyPASS_vcf); os.unlink(onlyPASS_withHeader_vcf)
            all_correctlyFormatted_vcfs.append(onlyPASS_withHeader_vcf_gz)

        # get the intersection vcf
        nvcfs = len(all_correctlyFormatted_vcfs)
        intersection_vcf = "%s/variantsPASS.vcf"%cwd; intersection_vcf_gz = "%s.gz"%intersection_vcf

        if nvcfs>1:

            run_cmd("%s isec -o %s -O v --threads %i --nfiles=%i %s"%(bcftools, intersection_vcf, threads, nvcfs, " ".join(all_correctlyFormatted_vcfs)))

            # format as a vcf 4.2
            df_vcf = pd.read_csv(intersection_vcf, sep="\t", header=None, names=["#CHROM", "POS", "REF", "ALT", "QUAL"])
            expected_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "DATA1"]
            for f in set(expected_fields).difference(set(df_vcf.keys())): df_vcf[f] = ["."]*len(df_vcf)
            open(intersection_vcf, "w").write("##fileformat=VCFv4.2\n%s"%(df_vcf[expected_fields].to_csv(sep="\t", index=False, header=True)))

            # compress and index
            run_cmd("%s -c %s > %s; %s -p vcf %s"%(bgzip, intersection_vcf, intersection_vcf_gz, tabix, intersection_vcf_gz))

        else: 
            copyfile(all_correctlyFormatted_vcfs[0], intersection_vcf_gz)
            run_cmd("%s -p vcf %s"%(tabix, intersection_vcf_gz))
      
        # get the alternative genome, if it is diploid it takes tha firts allele for heterozygous snps
        print("getting consensus sequence")
        alternative_genome_tmp = "%s.tmp"%alternative_genome
        run_cmd("%s consensus -f %s --haplotype 1 %s > %s"%(bcftools, reference_genome, intersection_vcf_gz, alternative_genome_tmp))
        os.rename(alternative_genome_tmp, alternative_genome)


def get_samples_with_wrong_pairs(paths_df, cwd):

    """Takes a paths df, such as the one that is passed to run_SNP_CNV_pipeline_various_options(), and returns a set of sampleIDs that do not have matched pairs. cwd is a dir to write files"""

    wrong_samples = set()

    for sample in set(paths_df.sampleID):

        # get df and each of the reads
        df = paths_df[paths_df.sampleID==sample]
        df1 = df[df.readID=="R1"]
        df2 = df[df.readID=="R2"]
        r1 = df1.raw_reads_dir.values[0]; r2 = df2.raw_reads_dir.values[0]
        #print(r1, r2)
        r1_firstLine = "%s/read1_firstline.txt"%cwd
        r2_firstLine = "%s/read2_firstline.txt"%cwd
        stderr = "%s/random_stdr.txt"%cwd

        # get the first line into a file
        run_cmd("zcat %s | head -n 1 > %s 2>%s"%(r1, r1_firstLine, stderr))
        run_cmd("zcat %s | head -n 1 > %s 2>%s"%(r2, r2_firstLine, stderr))

        # get into string
        r1_id = open(r1_firstLine, "r").readlines()[0].rstrip().split("length=")[0]
        r2_id = open(r2_firstLine, "r").readlines()[0].rstrip().split("length=")[0]

        if r1_id!=r2_id: wrong_samples.add(sample)

    return wrong_samples
        
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

def run_gatk_HaplotypeCaller(outdir_gatk, ref, sorted_bam, ploidy, threads, coverage, replace=False):

    """Runs haplotype caller under outdir and returns the filename of the filtered results"""

    # make the outdir if not there
    if not os.path.isdir(outdir_gatk): os.mkdir(outdir_gatk)

    # run GATK
    gatk_out = "%s/output.raw.vcf"%outdir_gatk; gatk_out_tmp = "%s.tmp"%gatk_out
    if file_is_empty(gatk_out) or replace is True:

        print("Running GATK HaplotypeCaller...")
        gatk_cmd = "%s HaplotypeCaller -R %s -I %s -O %s -ploidy %i --genotyping-mode DISCOVERY --emit-ref-confidence NONE --stand-call-conf 30 --native-pair-hmm-threads %i > %s.log"%(gatk, ref, sorted_bam, gatk_out_tmp, ploidy, threads, gatk_out); run_cmd(gatk_cmd)
        os.rename(gatk_out_tmp, gatk_out)

        # rename the index as well
        os.rename("%s.tmp.idx"%gatk_out, "%s.idx"%gatk_out)

    # variant filtration. There's a field called filter that has the FILTER argument
    gatk_out_filtered = "%s/output.filt.vcf"%outdir_gatk; gatk_out_filtered_tmp = "%s.tmp"%gatk_out_filtered
    if file_is_empty(gatk_out_filtered) or replace is True:

        print("Running GATK HaplotypeCaller Variant filtration...")

        # this depends on the ploidy. If ploidy is 2 you don't want to filter out heterozygous positions
        if ploidy==1: filterHeterozygous = '-G-filter-name "heterozygous" -G-filter "isHet == 1"'
        else: filterHeterozygous = ''

        gatk_filt_cmd = '%s VariantFiltration -V %s -O %s -cluster 5 -window 20 %s --filter-name "BadDepthofQualityFilter" -filter "DP <= %i || QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" > %s.log'%(gatk, gatk_out, gatk_out_filtered_tmp, filterHeterozygous , coverage, gatk_out_filtered); run_cmd(gatk_filt_cmd)
        os.rename(gatk_out_filtered_tmp, gatk_out_filtered)

        # rename the index as well
        os.rename("%s.tmp.idx"%gatk_out_filtered, "%s.idx"%gatk_out_filtered)

    # return the filtered file
    return gatk_out_filtered

def run_freebayes(outdir_freebayes, ref, sorted_bam, ploidy, threads, coverage, replace=False):

    # make the dir if not already done
    if not os.path.isdir(outdir_freebayes): os.mkdir(outdir_freebayes)

    #run freebayes
    freebayes_output ="%s/output.raw.vcf"%outdir_freebayes; freebayes_output_tmp = "%s.tmp"%freebayes_output
    if file_is_empty(freebayes_output) or replace is True:
        print("running freebayes")
        cmd_freebayes = "%s -f %s -p %i --min-coverage %i -b %s --haplotype-length -1 -v %s"%(freebayes, ref, ploidy, coverage, sorted_bam, freebayes_output_tmp); run_cmd(cmd_freebayes)
        os.rename(freebayes_output_tmp, freebayes_output)

    # filter the freebayes by quality
    freebayes_filtered = "%s/output.filt.vcf"%outdir_freebayes; freebayes_filtered_tmp = "%s.tmp"%freebayes_filtered
    if file_is_empty(freebayes_filtered) or replace is True:
        print("filtering freebayes")
        cmd_filter_fb = '%s -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" --tag-pass PASS %s > %s'%(vcffilter, freebayes_output, freebayes_filtered_tmp); run_cmd(cmd_filter_fb)
        os.rename(freebayes_filtered_tmp, freebayes_filtered)

    return freebayes_filtered

def run_freebayes_for_chromosome(chromosome_id, outvcf_folder, ref, sorted_bam, ploidy, coverage, replace):

    """Takes a chromosome ID and the fasta file and an outvcf and runs freebayes on it"""

    # define the output vcf file
    outvcf = "%s/%s_freebayes.vcf"%(outvcf_folder, chromosome_id); outvcf_tmp = "%s.tmp.vcf"%outvcf
    print("running freebayes for %s"%chromosome_id)

    # remove previously existing files
    if file_is_empty(outvcf) or replace is True:

        # generate the bam file for this chromosome (and index)
        sorted_bam_chr = "%s.%s.bam"%(sorted_bam, chromosome_id)
        run_cmd("%s view -b %s %s > %s"%(samtools, sorted_bam, chromosome_id, sorted_bam_chr))
        run_cmd("%s index -@ 1 %s"%(samtools, sorted_bam_chr))

        # get the fasta for the chromosome
        fasta_chromosome = "%s.%s.fasta"%(ref, chromosome_id)
        SeqIO.write([seq for seq in SeqIO.parse(ref, "fasta") if seq.id==chromosome_id], fasta_chromosome, "fasta")

        # run freebayes
        run_cmd("%s -f %s -p %i --min-coverage %i -b %s --haplotype-length -1 -v %s"%(freebayes, fasta_chromosome, ploidy, coverage, sorted_bam_chr, outvcf_tmp))

        # remove the intermediate files
        print("%s exists %s"%(fasta_chromosome, str(file_is_empty(fasta_chromosome))))
        remove_file(sorted_bam_chr); remove_file("%s.bai"%sorted_bam_chr); remove_file(fasta_chromosome); remove_file("%s.fai"%fasta_chromosome);

        # rename
        os.rename(outvcf_tmp, outvcf)

    # return the vcfs
    return outvcf

def run_freebayes_parallel(outdir_freebayes, ref, sorted_bam, ploidy, coverage, replace=False):

    """It parallelizes over the current CPUs of the system"""

    # make the dir if not already done
    if not os.path.isdir(outdir_freebayes): os.mkdir(outdir_freebayes)

    #run freebayes
    freebayes_output ="%s/output.raw.vcf"%outdir_freebayes; freebayes_output_tmp = "%s.tmp"%freebayes_output
    if file_is_empty(freebayes_output) or replace is True:

        print("running freebayes in parallel with %i threads"%(multiproc.cpu_count()))

        # define the chromosomes
        all_chromosome_IDs = [seq.id for seq in SeqIO.parse(ref, "fasta")]

        # remove the previous tmp file
        if not file_is_empty(freebayes_output_tmp): os.unlink(freebayes_output_tmp)

        # initialize the pool class with the available CPUs --> this is asyncronous parallelization
        pool = multiproc.Pool(multiproc.cpu_count())

        # make a dir to store the vcfs
        chromosome_vcfs_dir = "%s/chromosome_vcfs"%outdir_freebayes; make_folder(chromosome_vcfs_dir)

        # run in parallel the freebayes generation for all the 
        chromosomal_vcfs = pool.starmap(run_freebayes_for_chromosome, [(ID, chromosome_vcfs_dir, ref, sorted_bam, ploidy, coverage, replace) for ID in all_chromosome_IDs])

        # close the pool
        pool.close()

        # go through each of the chromosomal vcfs and append to a whole df
        all_df = pd.DataFrame()
        all_header_lines = []
        for vcf in chromosomal_vcfs:

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

        # check that all headers are the same
        if len(set(all_header_lines))!=1: 
            print("These are the header lines: ", set(all_header_lines))
            print("There are %i unique headers"%len(set(all_header_lines)))
            raise ValueError("Not all headers are the same in the individual chromosomal vcfs. This may indicate a problem with parallelization of freebayes")

        # write the file
        open(freebayes_output_tmp, "w").write(all_header_lines[0] + all_df[vcf_header].to_csv(sep="\t", index=False, header=True))

        # rename
        os.rename(freebayes_output_tmp, freebayes_output)

    # filter the freebayes by quality
    freebayes_filtered = "%s/output.filt.vcf"%outdir_freebayes; freebayes_filtered_tmp = "%s.tmp"%freebayes_filtered
    if file_is_empty(freebayes_filtered) or replace is True:
        print("filtering freebayes")
        cmd_filter_fb = '%s -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" --tag-pass PASS %s > %s'%(vcffilter, freebayes_output, freebayes_filtered_tmp); run_cmd(cmd_filter_fb)
        os.rename(freebayes_filtered_tmp, freebayes_filtered)

    return freebayes_filtered

def get_onlySNPs_from_vcf_filtering(vcf, filterTAG="PASS"):

    """Takes a vcf and writes only the SNPs that pass the filter (which can be also no)"""

    # load the df
    header_list = [l for l in open(vcf, "r") if l.startswith("##")]
    df = pd.read_csv(vcf, sep="\t", header = len(header_list))

    # get only the ones that are SNPs
    bases = {"A", "C", "G", "T", "N", "a", "c", "g", "t", "n"}
    df = df[df.apply(lambda r: r["REF"] in bases and r["ALT"] in bases, axis=1)]

    # filter
    if filterTAG!="no": df = df[df.FILTER==filterTAG]

    # write and return
    filtered_vcf = "%s.SNPs.%sfilt.vcf"%(vcf, filterTAG)
    open(filtered_vcf, "w").write("%s%s"%("".join(header_list), df.to_csv(sep="\t", header=True, index=False)))
    return filtered_vcf



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

        # initialize a dictionary that will store chromosome, position and mappability_score as lists
        expanded_data_dict = {"chromosome":[], "position":[], "unique_map_score":[]}

        # go through each row of the dataframe and append the lists
        for chromosome_list, positions_list, map_idx_list in df_map[["chromosome_list", "positions_list", "map_idx_list"]].values:

            expanded_data_dict["chromosome"] += chromosome_list
            expanded_data_dict["position"] += positions_list
            expanded_data_dict["unique_map_score"] += map_idx_list

        df_long = pd.DataFrame(expanded_data_dict)
        df_long["is_uniquely_mappable"] = (df_long.unique_map_score>=1.0).apply(int)

        # save
        df_long.to_csv(map_outfile_long, sep="\t", header=True, index=False)

    return map_outfile_long

def generate_nt_content_file(genome, target_nts="GC", replace=False):

    """Takes a genome and outputs a file with chromosome, position and 1 or 0 regarding if any of the target_nts is the same in the genome. This is 0-based"""

    target_nt_content_file = "%s.%scontent.tab"%(genome, target_nts)

    if file_is_empty(target_nt_content_file) or replace is True:

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

def run_gridss_raw(gridss_outdir, sorted_bam, ref, threads=4, replace=False):

    """This file runs the raw gridss into gridss_outdir, returning the output vcf"""

    # make the directories to run gridss
    if not os.path.isdir(gridss_outdir): os.mkdir(gridss_outdir)
    gridss_tmpdir = "%s/tmp"%gridss_outdir
    if not os.path.isdir(gridss_tmpdir): os.mkdir(gridss_tmpdir)

    # create the vcf of the structural variants
    gridss_VCFoutput = "%s/output.vcf"%gridss_outdir
    gridss_assemblyBAM = "%s/assembly.bam"%gridss_outdir # this file is intermediate of gridss and will be created in between

    if file_is_empty(gridss_VCFoutput) or replace is True:

        # remove any previously generated files 
        if not file_is_empty(gridss_VCFoutput): os.unlink(gridss_VCFoutput)
        if not file_is_empty(gridss_assemblyBAM): os.unlink(gridss_assemblyBAM)

        # define the threads were gridss works on (it turns out that more than 8 threads is not reccommended):
        gridss_threads = threads
        if gridss_threads>8: gridss_threads = 8

        print("Running gridss")
        gridss_cmd = "%s --jar %s --reference %s -o %s --assembly %s --threads %i --workingdir %s %s"%(gridss_run, gridss_jar, ref, gridss_VCFoutput, gridss_assemblyBAM, gridss_threads, gridss_tmpdir, sorted_bam); run_cmd(gridss_cmd)

    return gridss_VCFoutput

def run_vep(input_vcf, ref, gff_with_biotype, mitochondrial_chromosome, mitochondrial_code, gDNA_code, replace=False):

    """Runs VEP for an input vcf returning the annotated vcf"""

    # define an output file for VEP
    annotated_vcf = "%s_annotated.tab"%input_vcf; annotated_vcf_tmp = "%s.tmp"%annotated_vcf

    # run annotation by VEP
    if file_is_empty(annotated_vcf) or replace is True:

        print("Annotating with VEP %s"%input_vcf)
        vep_cmd = "%s -i %s -o %s -f %s -gf %s -mch %s -mcode %i -gcode %i"%(ensembl_vep, input_vcf, annotated_vcf_tmp, ref, gff_with_biotype, mitochondrial_chromosome, mitochondrial_code, gDNA_code); run_cmd(vep_cmd)
        os.rename(annotated_vcf_tmp, annotated_vcf)

    return annotated_vcf

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

def run_clove_raw(bedpe_file, sorted_bam, df_cov, mitochondrial_chromosome, replace=False, wrong_filters=["NO_ASSEMBLY"]):

    """Takes a bedpe file and a sorted bam and generates the clove output, only for breakends that do not have the wrong_filers. df_cov is a df of the first 4 columns of the mpileup output"""

    # get the output file
    output_vcf_clove = "%s.clove.vcf"%bedpe_file; output_vcf_clove_tmp = "%s.tmp.vcf"%output_vcf_clove

    if file_is_empty(output_vcf_clove) or replace is True:

        print("Generating the raw clove vcf")

        # write a bedpe with the correct filters
        bedpe_interesting_events = "%s.filtered_events"%bedpe_file
        df_bedpe = pd.read_csv(bedpe_file, sep="\t")
        df_filt = df_bedpe[df_bedpe.FILTERS.apply(lambda f: all([wf not in f for wf in wrong_filters]))][["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]]
        df_filt.to_csv(bedpe_interesting_events, sep="\t", header=False, index=False)

        # calculate the median coverage
        median_cov = np.median(df_cov.loc[set(df_cov.index).difference({mitochondrial_chromosome}), "coverage"])

        # run clove
        run_cmd("%s -jar %s -i %s BEDPE -b %s -o %s -c %i %i"%(JAVA, clove, bedpe_interesting_events, sorted_bam, output_vcf_clove_tmp, median_cov, 1))
        os.rename(output_vcf_clove_tmp, output_vcf_clove)

    return output_vcf_clove

def run_gridss_pipeline(gridss_outdir, sorted_bam, ref, gff_with_biotype, mitochondrial_chromosome, mitochondrial_code, gDNA_code, dataframe_mpileup, threads=4, replace=False):

    """This function runs the whole gridss pipeline given a reference genome and a sorted bam"""

    # first run gridss output to generate the raw vcf file
    gridss_VCFoutput = run_gridss_raw(gridss_outdir, sorted_bam, ref, threads=threads, replace=replace)

    # now run VEP to get the genes in the breakpoints
    annotated_vcf = run_vep(gridss_VCFoutput, ref, gff_with_biotype, mitochondrial_chromosome, mitochondrial_code, gDNA_code, replace=replace)

    # generate the bedpe file 
    bedpe_file = get_bedpe_from_svVCF(gridss_VCFoutput, gridss_outdir, replace=replace)

    # load the dataframe that has the coverage info
    df_cov = load_object(dataframe_mpileup)

    # run clove
    output_vcf_clove = run_clove_raw(bedpe_file, sorted_bam, df_cov, mitochondrial_chromosome, replace=replace)

    # filter clove outpt
    output_vcf_clove_filteredByCov = "%s.filtCoverage.tab"%output_vcf_clove
    if file_is_empty(output_vcf_clove_filteredByCov) or replace is True:

        print("Filtering clove output according to the coverage of the region")
        df_filteredClove = write_clove_output_filteredByCoverage(output_vcf_clove, output_vcf_clove_filteredByCov, df_cov, bedpe_file, gridss_VCFoutput, mitochondrial_chromosome=mitochondrial_chromosome)

    print("The final annotated vcf is in %s"%output_vcf_clove_filteredByCov)

def generate_gcprofiles(reference_genome, gcprofiles_path, replace=False, window_l=1000, threads=4, with_header=False, with_mappability=True):

    """Takes a reference genome and generates the GC profiles path.
    This is a file that has, for windowns of window_l, the following fields:
    -chromosome
    -start position (0-based)
    -GC-content
    -percentage of ACGT-letter per window (1-poly(N)%)
    -*** percentage of uniquely mappable positions per window --> calculated with genmap
    """

    if file_is_empty(gcprofiles_path) or replace is True:

        print("writing GC profiles")

        gc_content_profile_allPositions = "%s.all_positions_window%ibp_mappability%i.tab"%(gcprofiles_path, window_l, int(with_mappability))

        if file_is_empty(gc_content_profile_allPositions) or replace is True:

            print("generating GC profiles")

            # get the GC content file
            gc_content_outfile = generate_nt_content_file(reference_genome, replace=replace, target_nts="GC")
            gc_df = pd.read_csv(gc_content_outfile, sep="\t")[["chromosome", "position", "is_in_GC"]]

            # get the ACGT content file
            atgc_content_outfile = generate_nt_content_file(reference_genome, replace=replace, target_nts="ATGC")
            atgc_df = pd.read_csv(atgc_content_outfile, sep="\t")[["chromosome", "position", "is_in_ATGC"]]

            # define intermediate files
            intermediate_files = [gc_content_outfile, atgc_content_outfile]

            # get the mappabiliy file
            if with_mappability is True:
                mappability_out_file = generate_genome_mappability_file(reference_genome, replace=replace, threads=threads)
                map_df = pd.read_csv(mappability_out_file, sep="\t")[["chromosome", "position", "is_uniquely_mappable"]]
                intermediate_files.append(mappability_out_file)

            # merge the dataframes
            df_merged = gc_df.merge(atgc_df, left_on=["chromosome", "position"], right_on=["chromosome", "position"], validate="one_to_one")
            if with_mappability is True: df_merged = df_merged.merge(map_df, left_on=["chromosome", "position"], right_on=["chromosome", "position"], validate="one_to_one")

            df_merged.to_csv(gc_content_profile_allPositions, sep="\t", header=True, index=False)

            # delete interemediate files
            for f in intermediate_files: remove_file(f)

        else: df_merged = pd.read_csv(gc_content_profile_allPositions, sep="\t")

        # add the start of the window of length window_l
        df_merged["window_start"] = ((df_merged.position / window_l).apply(int))*window_l
        df_merged["window_start_str"] = df_merged.window_start.apply(str)
        df_merged["window_idx"] = df_merged.chromosome + "--||--" + df_merged.window_start_str

        # group by the idx
        mean_g_df = df_merged.groupby("window_idx").mean()
        mean_g_df["window_idx"] = mean_g_df.index
        mean_g_df["window_start"] = mean_g_df.window_start.apply(int)
        mean_g_df["chromosome"] = mean_g_df.window_idx.apply(lambda x: x.split("--||--")[0])
        mean_g_df = mean_g_df.sort_values(by=["chromosome", "window_start"])

        # define the expected cols
        col_to_meanCol_name = {"is_in_GC": "meanGC", "is_in_ATGC": "meanATGC"}
        if with_mappability is True: col_to_meanCol_name["is_uniquely_mappable"] = "meanUniquelyMappable"

        finalCols = ["chromosome", "window_start", "meanGC", "meanATGC"]
        if with_mappability is True: finalCols.append("meanUniquelyMappable")

        # rename cols
        mean_g_df = mean_g_df.rename(columns=col_to_meanCol_name)    

        # write to csv without header
        mean_g_df[finalCols].to_csv(gcprofiles_path, sep="\t", header=with_header, index=False)


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

def write_coverage_per_gene_mosdepth_and_parallel(sorted_bam, reference_genome, cnv_outdir, bed, gene_to_coverage_file, replace=False):

    """Takes a bam, a bed (1-indexed) and a file were the gene-to-coverage should be written, and writes a file with the gene_to_coverage info """

    if file_is_empty(gene_to_coverage_file) or replace is True:
        print("generating %s"%gene_to_coverage_file)

        # get the first three cols, important for the coverage calc
        cnv_bed = "%s/%s"%(cnv_outdir, get_file(bed))
        run_cmd("cut -f1-3 %s > %s"%(bed, cnv_bed))

        # get the coverage file for the first three
        coverage_file = generate_coverage_per_window_file_parallel(reference_genome, cnv_outdir, sorted_bam, windows_file=cnv_bed, replace=replace, window_l=1000, run_in_parallel=True, delete_bams=True)

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

        print("writing %s"%gene_to_coverage_file)
        os.rename(coverage_file_with_ID, gene_to_coverage_file)


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

def plot_coverage_correlations(df_cov, PlotsDir, plots={"simple_correlations", "goodness_coverage_predictions"}, mitochondrial_chromosome="mito_C_glabrata_CBS138"):

    """Takes a dataframe such as the output of get_coverageANDco_per_window_df and makes many plots of coverage correlations"""

    # make the folder if not done
    make_folder(PlotsDir)

    # a plot to correlate GC content and coverage for each chromosome
    if "simple_correlations" in plots:

        for field_correlate in {"meanGC", "fraction_covered_by_repeat", "meanUniquelyMappable", "dist_telomere", "rel_dist_telomere", "n_repeats"}:

            # calculate the r values for each chromosome
            chrom_to_corrvals = {c : scipy.stats.spearmanr(df_cov[df_cov.chromosome==c][[field_correlate, "mediancov_1"]], axis=0)  for c in set(df_cov.chromosome)}
            chrom_to_corrvals_str = {c : "% s, r=%.3f, p=%.4f"%(c, sp.correlation, sp.pvalue) for c, sp in chrom_to_corrvals.items()}
            df_cov["label_for_chr"] = df_cov.chromosome.apply(lambda x: chrom_to_corrvals_str[x])

            # define nuclear and mito chrom
            df_mito = df_cov[df_cov.chromosome==mitochondrial_chromosome]
            df_nuclear = df_cov[df_cov.chromosome!=mitochondrial_chromosome]

            fig = plt.figure(figsize=(20, len(df_nuclear.chromosome.unique())))

            # plot for each genome
            ax =  plt.subplot(2, 1, 1)
            sns.scatterplot(x=field_correlate, y="mediancov_1", data=df_nuclear, s=20, hue="label_for_chr")

            ax = plt.subplot(2, 1, 2)
            sns.scatterplot(x=field_correlate, y="mediancov_1", data=df_mito, s=20, hue="label_for_chr")

            plt.show()
            fig.tight_layout()  # otherwise the right y-label is slightly 
            filename="%s/%s_vs_coverage.pdf"%(PlotsDir, field_correlate)
            fig.savefig(filename, bbox_inches='tight');
            plt.close(fig)

    # subplots to see the correlation between each of the predicted coverages and the actual ones
    if "goodness_coverage_predictions" in plots:

        all_chromosomes = df_cov.chromosome.unique()

        fig = plt.figure(figsize=(10, len(all_chromosomes)*4.5))
        for I, chrom in enumerate(all_chromosomes):

            # get df
            df_c = df_cov[df_cov.chromosome==chrom]

            # calculate the r2 of the different fits with the data
            r2_from_loc = r2_score(df_c.rel_coverage_median, df_c.rel_cov_predicted_from_location)
            r2_all = r2_score(df_c.rel_coverage_median, df_c.predictedFromAllFeats_relative_coverage)

            # print each of the predictions
            ax = plt.subplot(len(all_chromosomes), 1, I+1)
            sns.scatterplot(x="start", y="rel_coverage_median", data=df_c, s=20, color="black", label="observed")
            sns.lineplot(x="start", y="rel_cov_predicted_from_location", data=df_c, linewidth=2, color="blue", label="from position, r²=%.3f"%r2_from_loc)
            sns.lineplot(x="start", y="predictedFromAllFeats_relative_coverage", data=df_c, linewidth=1, color="red", label="from all, r²=%.3f"%r2_all)


            plt.ylim([0, 4])

            ax.set_title(chrom)

        fig.tight_layout()  # otherwise the right y-label is slightly 
        filename="%s/relative_coverage_fitting.pdf"%(PlotsDir)
        fig.savefig(filename, bbox_inches='tight');
        plt.close(fig)

def run_repeat_modeller(reference_genome, threads=4, replace=False):

    """Runs repeat modeller to get a fasta file were all the repeats are. The part of integrating the results does not work very well, so that we have to integrate the results into one with RepeatClassifier afterwards.

    I had to improve the databases.
    In the RepeatMasker/Libraries I ran ' makeblastdb -in RepeatPeps.lib -input_type fasta -dbtype prot' to get the proteins, and then copied to RepeatModeler/Libraries (cp -r * ../../RepeatModeler/Libraries/)

    and also for the lib  makeblastdb -in RepeatMasker.lib -input_type fasta -dbtype nucl

    and makeblastdb -in simple.lib  -input_type fasta -dbtype nucl in RepeatModeler

    """

    # get the genome into outdir
    outdir = "%s.repeat_modeler_outdir"%reference_genome; 
    
    # delete everything
    if replace is True: delete_folder(outdir)

    # create again the genome
    make_folder(outdir)
    possible_genome_dirs = ["%s/%s"%(outdir, f) for f in os.listdir(outdir) if f.startswith("reference_genome_") and f.endswith("fasta")]
    if len(possible_genome_dirs)==0 or replace is True: genome_dir = "%s/reference_genome_%s.fasta"%(outdir, id_generator(10))
    else: genome_dir = possible_genome_dirs[0]

    if file_is_empty(genome_dir): run_cmd("cp  %s %s"%(reference_genome, genome_dir))

    # define the outdir
    repeat_modeler_outfile = "%s-families.fa"%genome_dir

    if file_is_empty(repeat_modeler_outfile) or replace is True:

        # config file
        config_file = "%s/configure"%repeatmoder_dir

        # repeat modeler
        repeat_modeller = "%s/RepeatModeler"%repeatmoder_dir

        # repeat classifier
        repeat_classifier = "%s/RepeatClassifier"%repeatmoder_dir

        # change the path to the place where the config file is
        """
        os.chdir(repeatmoder_dir)

        # configure repeat modeler perl
        
        # this gives problems when running in the cluster for many jobs in parallel, so it is skipped. This is sp because the RepModelConfig.pm is deleted and created again!!

        cmd_config = "%s %s -cdhit_dir %s -abblast_dir %s -genometools_dir %s -ltr_retriever_dir %s -ninja_dir %s -rmblast_dir %s -repeatmasker_dir %s -trf_prgm %s -recon_dir %s -mafft_dir %s -rscout_dir %s"%(perl, config_file, cdhit_dir, abblast_dir, genometools_dir, ltr_retriever_dir, ninja_dir, rmblast_dir, repeatmasker_dir, trf_prgm_dir, recon_dir, mafft_dir, rscout_dir)

        run_cmd(cmd_config)
        """
        
        os.chdir(outdir)

        # run the database
        name_database = get_file(genome_dir)
        run_cmd("%s -name %s %s"%(repeat_modeller_BuildDatabase, name_database, genome_dir))

        # run repeatmodeller
        njobs = int(threads/4) # Specify the number of parallel search jobs to run. RMBlast jobs wil use 4 cores each and ABBlast jobs will use a single core each. i.e. on a machine with 12 cores and running with RMBlast you would use -pa 3 to fully utilize the machine

        #raise ValueError("This has to be fixed!!!!")
        cmd = "export PERL5LIB=%s; %s -database %s -pa %i -LTRStruct"%(repeatmoder_dir, repeat_modeller, name_database, njobs)

        # add the location were eveything is installed and run
        print("running repeatmodeler...")
        cmd += " -abblast_dir %s -cdhit_dir %s -genometools_dir %s -ltr_retriever_dir %s -mafft_dir %s -ninja_dir %s -recon_dir %s -repeatmasker_dir %s -rmblast_dir %s -rscout_dir %s -trf_prgm %s"%(abblast_dir, cdhit_dir, genometools_dir, ltr_retriever_dir, mafft_dir, ninja_dir, recon_dir, repeatmasker_dir, rmblast_dir, rscout_dir, trf_prgm_dir)

        run_cmd(cmd)

    return repeat_modeler_outfile


def run_repeat_masker(reference_genome, threads=4, replace=False, use_repeat_modeller=True):

    """
    It runs repeat masker for a reference genome, writing the results under a folder where the ref genome is
    """

    # get the library from repeat_modeller
    if use_repeat_modeller is True: library_repeats_repeatModeller =  run_repeat_modeller(reference_genome, threads=threads, replace=replace)

    # define the repear masker outdir
    genome_dir = "/".join(reference_genome.split("/")[0:-1])
    genome_name = reference_genome.split("/")[-1]

    # define the outdirs for each 
    repeat_masker_outdir = "%s/%s_repeat_masker_outdir"%(genome_dir, genome_name.split(".")[0]); make_folder(repeat_masker_outdir)
    repeat_masker_outdir_default = "%s/default"%repeat_masker_outdir; make_folder(repeat_masker_outdir_default)
    repeat_masker_outdir_personal = "%s/personal"%repeat_masker_outdir; make_folder(repeat_masker_outdir_personal)

    # run in the default configuration
    repeat_masker_outfile_default = "%s/%s.out"%(repeat_masker_outdir_default, genome_name)
    if file_is_empty(repeat_masker_outfile_default) or replace is True:
        print("running repeatmasker to get the repeats of the genome in the default configuration")
        run_cmd("%s -pa %i -dir %s -poly -html -gff %s"%(repeat_masker, threads, repeat_masker_outdir_default, reference_genome))


    # run in the personal configuration
    repeat_masker_outfile_personal = "%s/%s.out"%(repeat_masker_outdir_personal, genome_name)

    if use_repeat_modeller is True:
        
        if file_is_empty(repeat_masker_outfile_personal) or replace is True:
            print("running repeatmasker to get the repeats of the genome with the lib obtained with RepeatModeler")
            run_cmd("%s -pa %i -dir %s -poly -html -gff %s -lib %s"%(repeat_masker, threads, repeat_masker_outdir_personal, reference_genome, library_repeats_repeatModeller))

    else: 

        # empty file
        run_cmd("head -n 3 %s > %s"%(repeat_masker_outfile_default, repeat_masker_outfile_personal))

       
    return repeat_masker_outfile_personal, repeat_masker_outfile_default

def get_repeat_maskerDF(reference_genome, threads=4, replace=False):

    """gets the repeat masker outfile as a pandas df. The repeatmasker locations are 1-based (https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/rmsk2bed.html)"""

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
    print("getting both repeats df")
    df = pd.concat(df_list).sort_values(by=["chromosome", "begin_repeat", "end_repeat"])
    df = df.drop_duplicates(subset=[x for x in df.keys() if x not in {"IDrepeat"}], keep="first")
    df["IDrepeat"] = list(range(1, len(df)+1))
    df.index = list(range(1, len(df)+1))

    return df

def run_samtools_index(genome, replace=False):

    """Runs samtools index for a genome of interest"""

    index_file = "%s.fai"%genome

    if file_is_empty(index_file) or replace is True:

        print ("Indexing the reference...")
        run_cmd("%s faidx %s"%(samtools, genome)) 

def makeblastdb_for_genome(genome, replace=False):

    """Runs makeblastdb genome"""   

    if any([file_is_empty("%s.%s"%(genome, suffix)) for suffix in {"nhr", "nin", "nsq"}]) or replace is True:

        print("running makeblastdb")
        run_cmd("%s -in %s -input_type fasta -dbtype nucl"%(makeblastdb, genome))

def get_repeat_coverageDF(df_with_windows, destination_dir, reference_genome, threads=4, replace=False):

    """Takes a df that has chromosome, start and end and returns the coverage of repeats. write files under destination_dir"""

    # write a file that has the windows
    regions_file = "%s/regions_file_repeatCoverage.bed"%destination_dir
    df_with_windows[["chromosome", "start", "end"]].to_csv(regions_file, sep="\t", header=False, index=False)

    # write a file with the repeats
    repeat_masker_df = get_repeat_maskerDF(reference_genome, threads=threads, replace=replace)
    repeats_file = "%s/repeats_file_repeatCoverage.bed"%destination_dir
    repeat_masker_df[["chromosome", "begin_repeat", "end_repeat"]].to_csv(repeats_file, sep="\t", header=False, index=False)

    # run fasta index and makeblastdb on reference_genome to get the files necessary
    run_samtools_index(reference_genome, replace=replace)
    makeblastdb_for_genome(reference_genome, replace=replace)

    # get the coverage of repeats on regions
    coverage_repeatsONregions_file =  "%s/coverage_repeatsONregions.bed"%destination_dir; coverage_repeatsONregions_file_tmp = "%s.tmp"%coverage_repeatsONregions_file
    if file_is_empty(coverage_repeatsONregions_file) or replace is True:
        stderr="%s/stderr.txt"%destination_dir
        run_cmd("%s coverage -a %s -b %s -g %s.fai > %s 2>%s"%(bedtools, regions_file, repeats_file, reference_genome, coverage_repeatsONregions_file_tmp, stderr))
        os.rename(coverage_repeatsONregions_file_tmp, coverage_repeatsONregions_file)

    # get as df and check that it is good
    df_coverage = pd.read_csv(coverage_repeatsONregions_file, sep="\t", header=None, names=["chromosome", "start", "end", "n_repeats", "n_bp_covered_by_repeat", "length_region", "fraction_covered_by_repeat"])

    # add the relative number of repeats to the length of the region
    df_coverage.loc[df_coverage.index, "relative_n_repeats"] = df_coverage.n_repeats / df_coverage.length_region

    # remove files
    for f in [repeats_file, regions_file, coverage_repeatsONregions_file]: os.unlink(f)

    return df_coverage

def get_kmercontent_region(chromosome, start, end, chromID_to_seq, dest_dir, replace, k):

    """Takes coordinates of a region and returns the kmer couting file, according to a dictionary that maps chromosomeID to sequence and a destination dir where to replace stuff. It returns a dictionary with the kmers"""

    # first generate the sequence
    ID = "%s_%i_%i"%(chromosome, start, end)
    seq_file = "%s/%s.fasta"%(dest_dir, ID)
    if file_is_empty(seq_file) or replace is True: open(seq_file, "w").write(">%s\n%s\n"%(ID, chromID_to_seq[chromosome][start:end]))

    # now run the kmer counting
    outfile_kmer_count = "%s/%s_kmers.fasta"%(dest_dir, ID)
    stdfile = "%s/%s_std.txt"%(dest_dir, ID)
    if file_is_empty(outfile_kmer_count) or replace is True: run_cmd("%s in=%s out=%s k=%i overwrite=true > %s 2>&1"%(kmercountexact, seq_file, outfile_kmer_count, k, stdfile))

    # return a dictionary with the kmer count
    return [(chromosome, start, end) , {str(seq.seq) : int(seq.id) for seq in SeqIO.parse(outfile_kmer_count, "fasta")}]

def get_kmercontent_per_window_df(df_with_windows,  destination_dir, reference_genome, threads=4, replace=False, k=6):

    """Takes a dataframe that has windows and it returns, for each window, the """
    print("working on kmer obtention")

    # define a file that will contain the final kmer df
    kmer_table_file = "%s/kmers_per_region.tab"%destination_dir

    if file_is_empty(kmer_table_file) or replace is True:

        # define a folder to write files
        writing_files_folder = "%s/kmer_counting_files"%destination_dir; 
        delete_folder(writing_files_folder) # this is to debug previous runs
        make_folder(writing_files_folder)

        # write the sequences
        chrom_to_seq = {c.id : c.seq for c in SeqIO.parse(reference_genome, "fasta")}

        # go through each 
        print("counting kmers in %i threads..."%multiproc.cpu_count())

        # run kmercounting in parallel
        pool = multiproc.Pool(multiproc.cpu_count())

        # run in parallel the coverage generation, which returns a list of dataframes, each with one chromosome
        list_region_and_kmerToCount = pool.starmap(get_kmercontent_region, [(c, s, e, chrom_to_seq, writing_files_folder, replace, k) for c, s, e in df_with_windows[["chromosome", "start", "end"]].values])

        # close the pool
        pool.close()

        # make a dataframe containing the coordingates as index and the kmers
        print("creating kmer dataframe")
        df = pd.DataFrame(dict(pd.DataFrame(list_region_and_kmerToCount).set_index(0)[1])).transpose()

        # add the indices
        for I, field in enumerate(["chromosome", "start", "end"]): df[field] = [idx[I] for idx in df.index]

        # change the NaNs by 0s
        def changeNaNs(value):

            if pd.isna(value): return 0
            else: return value

        df = df.applymap(changeNaNs)

        # correct by the length of the region
        df_freqs = df[[x for x in df.keys() if x not in {"chromosome", "start", "end"}]]
        length_region = df.end - df.start
        df_freqs = df_freqs.div(length_region, axis=0)
        df_freqs = df_freqs.merge(df[["chromosome", "start", "end"]], left_index=True, right_index=True, validate="one_to_one")

        # write
        kmer_table_file_tmp = "%s.tmp"%kmer_table_file
        df_freqs.to_csv(kmer_table_file_tmp, sep="\t", index=False, header=True)

        # remove the folder where all files are
        shutil.rmtree(writing_files_folder)

        # rename at the end
        os.rename(kmer_table_file_tmp, kmer_table_file) 

    else: df_freqs = pd.read_csv(kmer_table_file, sep="\t")

    return df_freqs

def function_product_of_two_sigmoidal_functions(x, b1, c1, b2, c2):

    """This function takes an x value and many parameters and returns the product of two sigmoidal functions data."""

    f1 = 1 / (1 + np.exp(-c1 * (x - b1)))
    f2 = 1 / (1 + np.exp(-c2 * (x - b2)))

    return f1*f2
                  
def get_pvalue(iterable, value, type_test="phigher"):

    """Takes an iterable and a value. Returns the fraction of elemenst in iterable that are lower (plower) or higher (phigher) than value. """

    # prepare
    np_it = np.array(iterable)
    len_it = len(np_it)

    if type_test=="phigher": return sum(np_it>value)/len_it
    elif type_test=="plower": return sum(np_it<value)/len_it
    else: raise ValueError("type_test has to be phigher or plower")

def get_df_chunks_of_coverage(df, n_windows_defining_chunk=10):

    """Takes a dataframe of coverage and returns the main chunks observed. each chunk starts new if there is a statistically significant difference between the -n_windows_defining_chunk and the +n_windows_defining_chunk. doing a ks test

    n_following_windows_to_validate_chunk determines the number of consecutive windows necessary to support that there is a chunk

    This is testes with a couple of examples. It will likely produce 
    """

    chrom_to_len = {chrom : max(df[df.chromosome==chrom].end) for chrom in set(df.chromosome)}

    # define a function that will return, given a sorted df_coverage, the previous coverages, the ones starting at a given window and the statistic of the difference
    def get_window_coverage_info(df_cov, Ichunk, n_windows_defining_chunk):

        # get coverages
        previous_coverages = df_cov.rel_coverage_median.iloc[max([0, Ichunk-n_windows_defining_chunk]) : Ichunk]
        from_w_coverages = df_cov.rel_coverage_median.iloc[Ichunk : Ichunk+n_windows_defining_chunk]

        # get statistic
        statistic = abs(np.median(from_w_coverages) - np.median(previous_coverages)) / (np.std(from_w_coverages) + np.std(previous_coverages))

        # get the median dif
        mean_diff = abs(np.median(from_w_coverages) - np.median(previous_coverages))

        return previous_coverages, from_w_coverages, statistic, mean_diff

    df_chunks = pd.DataFrame()

    # go through each chromosome
    for chrom in df.chromosome.unique():
        df_c = df[df.chromosome==chrom].sort_values(by="start")

        # get all statistics of chunks of n_windows_defining_chunk, and also the mean differences
        all_statistics = []
        all_mean_differences = []

        # compute the mean and std differences between samples of adjacent positions rand
        for Ichunk, (start_w, end_w, cov) in enumerate(df_c[["start", "end", "rel_coverage_median"]].values): 

            previous_coverages, from_w_coverages, statistic, mean_diff = get_window_coverage_info(df_c, Ichunk, n_windows_defining_chunk)
            if not pd.isna(statistic): 

                all_statistics.append(statistic) 
                all_mean_differences.append(mean_diff)

        # go through the coverage and make windows
        chunks = [] # save start and end of each chunk
        start_current_chunk = df_c.start.iloc[0] # initialize the start of the current chunk

        for Ichunk, (start_w, end_w, cov) in enumerate(df_c[["start", "end", "rel_coverage_median"]].values):  

            # skip the first chunk
            if Ichunk==0: continue

            # define the the info about the comparison with the previous window
            previous_coverages, from_w_coverages, statistic, mean_diff = get_window_coverage_info(df_c, Ichunk, n_windows_defining_chunk)

            # make kstest
            ks, p_ks = scipy.stats.ks_2samp(previous_coverages, from_w_coverages)

            # calculate pvalue of a randomisation test
            p_random = get_pvalue(all_statistics, statistic, type_test="phigher")

            # calculate the pvalue of the median difference
            p_diff = get_pvalue(all_mean_differences, mean_diff, type_test="phigher")
            
            # if there are significant difs in both tests
            if p_ks<0.01 and p_random<0.01 and p_diff<0.05:

                # save until this window as a chunk
                chunks.append((start_current_chunk, int(end_w)))

                # initialize a new chunk
                start_current_chunk = df_c.start.iloc[Ichunk+1]

        # at the end you should end the chunk if not already done
        chunks.append((start_current_chunk, int(end_w)))

        # keep and sort
        chunks = sorted(set(chunks), key=(lambda x: x[0]))

        # keep as df
        df_c_chunks = pd.DataFrame(chunks).rename(columns={0:"start", 1:"end"})
        df_c_chunks["chromosome"] = [chrom]*len(df_c_chunks)
        
        # annotate the percentage of the chromosome that it covers
        df_c_chunks["fraction_chromosome"] = df_c_chunks.apply(lambda r: (r["end"]-r["start"])/chrom_to_len[r["chromosome"]], axis=1)

        df_chunks = df_chunks.append(df_c_chunks, sort=True)

    return df_chunks

def add_residual_coverageDF_after_chromosomePosition_correction(df, mitochondrial_chromosome="mito_C_glabrata_CBS138", make_plots=False, PlotsDir=None):

    """This function takes a dataframe with windows of a genome and coverage measures. It returns a DF in the same order with the residual of predicting based on distance to the telomere. If the fit is not good (rsq<0.5) it will assume that there is some strange rearrangement that does not allow a proper modelling, such an inversion. If this is the case we will return a fit function that is equivalent to the mean of the other nuclear chromosomes. There will not be any prediction for the mitochondrial chromosome, since it does not have any telomeres. For mito chromosome there will not be any prediction, and the residual will be the residual of a flat line at the median of the coverage of mito. If plots are indicated, it will create a plot for the fit in each chromosome """

    # for each chromosome, calculate the coeficients of the polynomial fit
    chrom_to_coefs = {} # a list of the coeficients of the polynomial fit
    wrong_chromosomes = set()
    for chromosome in df.chromosome.unique():

        # get the df for the chr
        df_c = df[(df.chromosome==chromosome)]

        # for the mitochondrial chromosome, just return a flat line at the median of the chromosome. Here we don't have telomeres
        if chromosome==mitochondrial_chromosome: coefs = [np.median(df_c.rel_coverage_median), 0, 0]
        
        # nuclear chromosomes. I will fit the product of two sigmoidal functions to the data, which should resemble the smiley face
        else:

            # fit the product with a 2nd degree polynmial, only taking into consideration the +- correct vals
            df_correct = df_c[(df_c.rel_coverage_median<=4) & (df_c.rel_coverage_median>0.05)]
            x_correct = df_correct.relative_position_in_chr; y_correct = df_correct.rel_coverage_median

            coefs = poly.polyfit(x_correct, y_correct, 2)
            yfit_correct = poly.polyval(x_correct, coefs)
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(y_correct, yfit_correct)
            rsquare_fit = r_value**2

            # define the central coverage
            if coefs[0]<0: central_coverage = max(yfit_correct) # This is an inverse smiley-face efect
            else: central_coverage = min(yfit_correct) # This is a smiley-face efect

            # define chromosomes that may be scrambled, which are the ones that fit wrong or the minimum is >1.5 (these may be whole chromosome duplications)
            if central_coverage>=1.5 or central_coverage<=0.1:

                wrong_chromosomes.add(chromosome)
                continue



        # keep 
        chrom_to_coefs[chromosome] = coefs

    # for chromosomes that don't have coefs, just append the mean of the others that are not mito
    non_nuclear_coefs_mean = list(pd.DataFrame([[x2, x1, x0] for chrom, (x2, x1, x0) in chrom_to_coefs.items() if chrom!=mitochondrial_chromosome]).mean())

    if len(wrong_chromosomes)>=(len(set(df.chromosome))-(len(set(df.chromosome))*0.5)):
        raise ValueError("There are too many wrong chromosomes according to the coverage fit. This could be a problem for future runs")

    # initialize the fineal chrom_to_coefs
    final_chrom_to_coefs = {}

    # get the contribution of the relative position
    df_predictions = pd.DataFrame()
    for chromosome in df.chromosome.unique():

        df_c = df[(df.chromosome==chromosome)]

        # if it is there
        if chromosome in chrom_to_coefs: coefs = chrom_to_coefs[chromosome]
        else: coefs = non_nuclear_coefs_mean

        # keep 
        final_chrom_to_coefs[chromosome] = coefs

        # get the fit data on all
        x = df_c.relative_position_in_chr; y = df_c.rel_coverage_median

        yfit = poly.polyval(x, coefs)

        # keep
        df_c.loc[df_c.index, "rel_cov_predicted_from_location"] = yfit
        df_c.loc[df_c.index, "residual_cov_vs_predicted_from_loc"] = df_c.rel_coverage_median - df_c.rel_cov_predicted_from_location
        df_predictions = df_predictions.append(df_c, sort=True)

        # make plots if indicated
        if make_plots is True: 
            fig = plt.figure(figsize=(10,5))
            ax = sns.scatterplot(x="start", y="rel_coverage_median", data=df_c, s=20, color="black")
            sns.lineplot(x="start", y="rel_cov_predicted_from_location", data=df_c, linewidth=1, color="black")

            # add the chunks 
            #for s, e in chunks_df[["start", "end"]].values: plt.plot(np.linspace(s, e, 3), np.linspace(1, 1, 3), linewidth=4)

            # get definition of parms
            ax.set_ylim([0,5])
            plt.show()
            plt.close(fig)

    return df_predictions, wrong_chromosomes, final_chrom_to_coefs

def calculate_rsquare_X_Y(x, y):

    """Takes an x and y arrays and calculates the rsquare of a linear fit between a linear fit of x and y. This cannot be run in parallel"""

    # get the prediction
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    pred_y = intercept + slope*x

    # return rsquared
    return r2_score(y, pred_y) # first the true, then the pred

def get_pval_train_test_df_and_field(train_idx, test_idx, df_X_features, y, xfeature, random_rsquares):

    """Takes a df with the X features and predicts the linear regression between xfeature and y for train_idx. It returns the pvalue on the test data """

    # define the x
    x = df_X_features[xfeature]

    # get the fit on train
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x[train_idx], y[train_idx])

    # get the prediction on test
    pred_y_test = intercept + slope*x[test_idx]

    # calculate the r2
    r2 = r2_score(y[test_idx], pred_y_test)

    # calculate the p val of the r2 on the random set, which is expected to be a numpy array
    pval = sum(random_rsquares>r2)/len(random_rsquares)

    return pval

def get_features_predicting_y_chrossvalidating_on_idx(df_X_features, y, min_signifficant_CVblocks=0.75, pval_tshd=0.05):

    """Takes a df with features that may predict y, where the index refers to CV blocks. It will return those fetaures that offer a signifficant corrrelation (p_spearnan<0.05 FDR corrected) in at least 50% of the CV blocks"""

    print("learning which hexamers predict coverage from interchromosomal CV")

    # define all feats
    all_features = list(df_X_features.keys())
    
    # map each idx to the other idxs
    all_idxs = set(df_X_features.index)
    idx_to_otheridxs = {idx : all_idxs.difference({idx}) for idx in all_idxs}

    # define a list of rsquares when fitting reshuffled fields in X on y
    random_df_X_features = df_X_features.iloc[np.random.permutation(len(df_X_features))]

    def calculate_rsquares(x): return calculate_rsquare_X_Y(x, y)
    random_rsquares = np.array(list(map(calculate_rsquares, [random_df_X_features[xfeature] for xfeature in all_features])))

    # initialize a dict
    feature_to_CVblock_to_rawp = {}

    # go through each idx, which represents a cross-validation block where we train on the non idx and test on the idx dataset
    for I, idx in enumerate(all_idxs):
        print("\tcrossvalidation on: ", I, idx)

        # define the train and test 
        train_idx = idx_to_otheridxs[idx]
        test_idx = idx

        # get a list of the pvals in a map
        raw_pvals = list(map(lambda xfeat: get_pval_train_test_df_and_field(train_idx, test_idx, df_X_features, y, xfeat, random_rsquares), all_features))

        # get correction
        corrected_pvals = list(multitest.fdrcorrection(raw_pvals)[1])

        # keep 
        feature_to_CVblock_to_rawp[idx] = dict(zip(all_features, corrected_pvals))

    # get as a df
    df_pvals = pd.DataFrame(feature_to_CVblock_to_rawp) # here the rows are feats and the columns are CV blocks. The values are the pvals

    # add the fraction of CV blocks that have a signifiicant pval
    df_pvals.loc[df_pvals.index, "fraction_sig"] = df_pvals.apply(lambda r: sum(r<pval_tshd)/len(r), axis=1)

    # get the features that have a high CV blocks being sig
    return set(df_pvals[df_pvals.fraction_sig>=min_signifficant_CVblocks].index)

def get_model_to_predict_residualCovFromLocation(df, replace=True):

    """It takes a df with coverage info (see get_sampleSpecific_genomeInfoDF_and_fitInfo) and returns a dictionary where each feature has its contribution to predicting the residual coverage after position correction. It will consider GC, unqueness, repeat coverage and significant hexamers. This df has to be filtered to include only regions considerable for the fitting"""

    print("predicting coverage from sequence features")

    # define the vars that you predict being important for the coverage
    known_predictors = ["meanGC", "fraction_covered_by_repeat", "relative_n_repeats", "meanUniquelyMappable"]
    X = df[known_predictors]

    # define dependent var
    y = df.loc[df.index, "residual_cov_vs_predicted_from_loc"]

    # fit a linear model to predict the residuals of all the known predictors
    lm = linear_model.LinearRegression()
    model = lm.fit(X,y)

    # get the predicted values
    predictions = pd.Series(lm.predict(X), index=y.index)

    # get the rsquare and print coefficienrs
    rsquare = lm.score(X,y)
    print("\tThe r2 of predicting the residual coverage after position correction from %s is %.5f"%(",".join(known_predictors), rsquare))
    print("\t  These are the coefficeints of each of these features")    
    for I, coeff in enumerate(lm.coef_): print("\t", known_predictors[I], "--> %.4f"%coeff)

    # get the residual of the prediction
    df.loc[df.index, "residual_cov_from_locANDknownFeats"] = y - predictions

    # get a df where the index is the chromosome
    df_idx_chr = df.set_index("chromosome")

    # find hexamers that predict this residuals from the sequence
    df_hexamers = df_idx_chr[[field for field in df.keys() if len(field)==6 and all([c in {"A", "C", "T", "G", "N"} for c in field])]]

    predictor_hexamers = get_features_predicting_y_chrossvalidating_on_idx(df_hexamers, df_idx_chr["residual_cov_from_locANDknownFeats"])

    ## predict with all the important predictors are ##

    # get all predictors
    all_predictors = known_predictors + list(predictor_hexamers)

    # get data
    X_all = df[all_predictors]

    # fit a linear model to predict the residuals of all the known predictors
    lm_all = linear_model.LinearRegression()
    model_all = lm_all.fit(X_all, y)
    predictions_all = lm_all.predict(X_all)
    rsquare_all = lm_all.score(X_all, y)

    # print for all
    print("...")
    print("\tThe r2 of predicting the residual coverage after position correction from %s is %.5f"%(",".join(all_predictors), rsquare_all))
    print("\tThese are the coefficeints of each of these features")    

    for I, coeff in enumerate(lm_all.coef_): print("\t", all_predictors[I], "--> %.4f"%coeff)

    return all_predictors, lm_all

def get_genome_info_df(reference_genome, replace=False, window_l=1000, threads=4, get_repeat_info=True, with_mappability=True):

    """Takes a genome and returns info that can be calculated for it from JUST the sequence."""

    # get the dir where to write stuff
    genome_dir = "/".join(reference_genome.split("/")[0:-1]); genome_name = reference_genome.split("/")[-1].split(".")[0]

    # get the df filename, where it will be saved
    df_info_filename = "%s/%s.genome_info_from_seq_windows%ibp_mappability%i_df.py"%(genome_dir, genome_name, window_l, int(with_mappability))

    if file_is_empty(df_info_filename) or replace is True:

        # first generate the GC profiles
        gcprofiles_path = "%s.gcprofiles_window%ibp_mappability%i.tab"%(reference_genome, window_l, int(with_mappability))
        generate_gcprofiles(reference_genome, gcprofiles_path, replace=replace, window_l=window_l, threads=threads, with_header=True, with_mappability=with_mappability)
        df_gcprofile = pd.read_csv(gcprofiles_path, sep="\t")

        # add the start
        df_gcprofile["start"] = df_gcprofile.window_start

        # add the end
        chrom_to_len = {c.id : len(c) for c in SeqIO.parse(reference_genome, "fasta")}
        ends = []
        for I, (start, chromosome) in enumerate(df_gcprofile[["start", "chromosome"]].values):

            if I==(len(df_gcprofile)-1): ends.append(chrom_to_len[chromosome]); continue

            # next chr is the same
            if df_gcprofile.chromosome.iloc[I+1]==chromosome: ends.append(df_gcprofile.start.iloc[I+1])

            # else add the end
            else: ends.append(chrom_to_len[chromosome])

        df_gcprofile["end"] = ends

        # get the repeat masker info
        if get_repeat_info is True:

            repeat_coverage_dest_dir = "%s/%s_repeats_%ibp"%(genome_dir, genome_name, window_l); make_folder(repeat_coverage_dest_dir)
            df_repeat_coverage = get_repeat_coverageDF(df_gcprofile, repeat_coverage_dest_dir, reference_genome, threads=threads, replace=replace)
            df = df_repeat_coverage.merge(df_gcprofile, left_on=["chromosome", "start", "end"],  right_on=["chromosome", "start", "end"], validate="one_to_one")
        
        else: df = df_gcprofile
        
        # get the hexamer content info
        destination_kmer_counting = "%s/%s_kmer_counting_windows%ibp"%(genome_dir, genome_name, window_l); make_folder(destination_kmer_counting)

        # this fails sometimes for unknown reasons. Try 5 times before finally droping
        for Itry in range(5):

            try: 
                df_hexamers = get_kmercontent_per_window_df(df,  destination_kmer_counting, reference_genome, threads=threads, replace=replace)
                break
            except: 
                print("WARNING: get_kmercontent_per_window_df failed in iteration %i...retrying"%(Itry+1))

        # check if it exists
        try: print("df_hexamers was created and has len=%i"%(len(df_hexamers)))
        except: raise ValueError("df_hexamers could not be created after %s tries"%(Itry+1))

        # merge
        df = df.merge(df_hexamers, left_on=["chromosome", "start", "end"],  right_on=["chromosome", "start", "end"], validate="one_to_one") 

        # add info of distance to the telomere
        
        chrom_to_len = {chrom : max(df[df.chromosome==chrom].end) for chrom in set(df.chromosome)}
        df.loc[df.index, "len_chromosome"] = df.chromosome.apply(lambda x: chrom_to_len[x])
        df.loc[df.index, "relative_position_in_chr"] = df.start / df.len_chromosome
        df.loc[df.index, "dist_telomere"] = df.apply(lambda r: min([r["start"] , r["len_chromosome"]-r["end"]]), axis=1)
        df.loc[df.index, "rel_dist_telomere"] = df.dist_telomere / df.len_chromosome

        # save
        save_object(df, df_info_filename)

    else: df = load_object(df_info_filename)
    
    return df

def get_sampleSpecific_genomeInfoDF_and_fitInfo(reference_genome, destination_dir, sorted_bam, mitochondrial_chromosome="mito_C_glabrata_CBS138", replace=False, window_l=1000, threads=4, replace_covModelObtention=False):

    """ Takes a genome and a sorted bam, returning the df for windows that have info about the genome including:
    - GC content
    - uniqueness
    - repeat coverage
    - relative coverage
    - hexamer couting
    - prediction of coverage based on location (quadratic fit) .

    It also retruns interesting info to predict coverage from the sequence features"""


    # first get all the info that depends solely on the genomic sequence
    df_seq_info =  get_genome_info_df(reference_genome, replace=replace, window_l=window_l, threads=threads)

    # generate the coverage file in a paralllel manner
    coverage_file = generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, replace=replace, window_l=window_l)
    df_coverage = pd.read_csv(coverage_file, sep="\t")

    # check that the two dfs gave the same
    if len(df_coverage)!=len(df_seq_info):

        print(df_coverage, df_seq_info)
        raise ValueError("There df_coverage and df_seq_info are not equal in length")

    # merged
    df = df_seq_info.merge(df_coverage, left_on=["chromosome", "start", "end"], right_on=["#chrom", "start", "end"], validate="one_to_one")

    # add the relative coverage to the median of all regions
    df.loc[df.index, "rel_coverage_median"] = df.mediancov_1 / np.median(df.mediancov_1)

    # add the residuals and fits of correcting the position, and also keep info about the scrambled chromosomes and the coefficients of each chromosome to positions
    df, scrambled_chromosomes, chromosome_to_relPosVScovCoefficients = add_residual_coverageDF_after_chromosomePosition_correction(df, mitochondrial_chromosome="mito_C_glabrata_CBS138", make_plots=False, PlotsDir=None)

    # define a set of windows that are useful to train a linear model that will be useful to predict the coverage that a scrabled genome will have depending on the relative position, GC content, uniqueness score, repeat presence and some hexamers

    df_filt = df[~(df.chromosome.isin(scrambled_chromosomes.union({mitochondrial_chromosome}))) & (df.rel_coverage_median<=4) & (df.rel_coverage_median>0.05)]
    
    # get the coefficients and features that predict the residual of predicting from the position
    tuple_predictors_LMobject_predict_fromPos_object_file = "%s/tuple_predictors_LMobject_predict_fromPos_object_windows_%ibp.py"%(destination_dir, window_l)
    if file_is_empty(tuple_predictors_LMobject_predict_fromPos_object_file) or replace_covModelObtention is True or replace is True:
        
        tuple_predictors_LMobject_predict_fromPos = get_model_to_predict_residualCovFromLocation(df_filt)
        save_object(tuple_predictors_LMobject_predict_fromPos, tuple_predictors_LMobject_predict_fromPos_object_file)

    else: tuple_predictors_LMobject_predict_fromPos = load_object(tuple_predictors_LMobject_predict_fromPos_object_file)

    # get the lm objects
    predictors_predictResidualFromPos, lmObject_predictResidualFromPos = tuple_predictors_LMobject_predict_fromPos

    # get the predicted relative coverage of each window
    df.loc[df.index, "predictedFromAllFeats_relative_coverage"] = df.rel_cov_predicted_from_location + lmObject_predictResidualFromPos.predict(df[predictors_predictResidualFromPos])

    return df, chromosome_to_relPosVScovCoefficients, predictors_predictResidualFromPos, lmObject_predictResidualFromPos, scrambled_chromosomes
   
def write_bed_heterozygous_positions(het_vcf, path_het_positions_bed, replace=False):

    """Takes a vcf and writes a bed with the heterozygous SNP positions"""

    # get into df
    df = load_single_sample_VCF(het_vcf)

    # filter PASS, heterozygous SNPs
    bases = {"A", "C", "G", "T", "N", "a", "c", "g", "t", "n"}
    df  = df[(df.FILTER=="PASS") & (df.is_heterozygous_GT)  & (df.is_heterozygous_coverage) & (df.REF.isin(bases)) & (df.ALT.isin(bases))]

    # add the position in a 0-based manner
    df["0_based_start"] = df.POS - 1
    df["0_based_end"] = df.POS

    # write as bed file
    df[["#CHROM", "0_based_start", "0_based_end"]].to_csv(path_het_positions_bed, sep="\t", header=False, index=False)

def make_folder(folder): 

    if not os.path.isdir(folder): os.mkdir(folder) 

def get_date():

    """Gets the date of today"""

    today = date.today()

    return today.strftime("%d/%m/%Y")

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

    # define chromosome_to_length
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

    # load translocations
    df = pd.read_csv(translocations_file, sep="\t")

    # define a function that takes a translocation and formats it
    df_corrected = df.apply(lambda r: format_translocation_row_simulateSV(r, chr_to_len), axis=1)[["Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB", "Balanced", "BpSeqA", "BpSeqB"]]

    # write the translocations
    df_corrected.to_csv("%s.corrected"%translocations_file, sep="\t", header=True, index=False)

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
    for field in ["BpSeqA", "BpSeqB_5prime", "BpSeqB_3prime"]: df[field] = df[field].apply(change_EmptyString_to_X)
    df["Copied"] = df.Copied.apply(lambda x: str(x).upper())

    # write as corrected
    df.to_csv("%s.corrected"%insertions_file, sep="\t", header=True, index=False)


def rearrange_genomes_simulateSV(reference_genome, outdir, repeats_1kbwindow_bed=None, replace=True, simulation_types=["uniform", "biased_towards_repeats"], mitochondrial_chromosome="mito_C_glabrata_CBS138", percBalancedTrans=0.5):

    """Genereates a randomly rearranged genome under outdir, either uniformly distributed (SV_uniformely_distributed) or biased towards +-1kb regions of repeats (SV_biased_towards_repeats, if repeats_1kbwindow_bed is provided). It returns a dict that maps the simulation type to the rearranged genome. This is based on running the simulateSV function for 50 events in the gDNA and 5 events in the mtDNA (10%). The program tries to simulate genomes that have a beta distribution of events size, in a way that  the longest event will have up to 20% of the shortest chromosome. Several genomes may not allow this to be true, so that what is done is to keep reducing the size until you find one that works.

    percBalancedTrans is the fraction of translocations that are biased"""

    # map each simulation type to the regions_bed argument that will be passed to simulateSV_R
    sim_type_to_regionsBed = {"uniform":"all_regions", "biased_towards_repeats":repeats_1kbwindow_bed}

    # make the outdir
    make_folder(outdir)

    # initialize the dict
    simType_to_rearrangedGenome = {}

    # go through each simulation
    for simulation_type in simulation_types:

        print("simulating SV for %s"%simulation_type)

        # define oudir for this simulation
        sim_outdir = "%s/%s"%(outdir, simulation_type); make_folder(sim_outdir)

        # define the expected genome
        rearranged_genome = "%s/rearranged_genome.fasta"%sim_outdir
        corrected_insertions = "%s/insertions.tab.corrected"%sim_outdir # this is the last to be generated

        # run the simulation if not already done
        if file_is_empty(rearranged_genome) or file_is_empty(corrected_insertions) or replace is True: 

            # calculate the maximum time of rearrangement proportional to the genome size
            genome_size = sum([len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")])
            max_time_rearrangement = int(genome_size/100000)

            # run the pipeline
            std_rearranging_genome = "%s/simulation_std.txt"%sim_outdir
            run_cmd("%s --input_genome %s --outdir %s --regions_bed %s --mitochondrial_chromosome %s --max_time_rearrangement %i  --percBalancedTrans %.2f > %s 2>&1"%(simulateSV_R, reference_genome, sim_outdir, sim_type_to_regionsBed[simulation_type], mitochondrial_chromosome, max_time_rearrangement, percBalancedTrans, std_rearranging_genome))

            # change the name of all the files to have the chromosome_Name
            svtype_to_chrFields = {"insertions":["ChrA", "ChrB"], "deletions":["Chr"], "inversions":["Chr"], "tandemDuplications":["Chr"], "translocations":["ChrA", "ChrB"]}

            for svtype in {"insertions", "deletions", "inversions", "tandemDuplications", "translocations"}:
                    
                # get the file as df
                file = "%s/%s.tab"%(sim_outdir, svtype)
                df = pd.read_csv(file, sep="\t")
                df["Name"] = df.apply(lambda r: "%s_%s"%(r["Name"], "-".join([r[field] for field in svtype_to_chrFields[svtype]])), axis=1)

                # write
                file_tmp = "%s.tmp"%file
                df.to_csv(file_tmp, sep="\t", index=False, header=True)
                os.rename(file_tmp, file)

            # edit the translocations so that the balanced ones are sorted
            translocations_file = "%s/translocations.tab"%sim_outdir
            rewrite_translocations_uniformizedFormat_simulateSV(translocations_file, reference_genome)

            # edit the insertions
            insertions_file = "%s/insertions.tab"%sim_outdir
            rewrite_insertions_uniformizedFormat_simulateSV(insertions_file)
            
        # keep
        simType_to_rearrangedGenome[simulation_type] = rearranged_genome

    print("rearrangement simulations done")
    return simType_to_rearrangedGenome



def get_windows_infoDF_with_predictedFromFeatures_coverage(genome, chromosome_to_relPosVScovCoefficients, predictors_ResCovAfterPosCorrection, lm_predictResidualFromPos, replace=False, window_l=1000, threads=4, make_plots=True, plots_prefix=None):

    """Takes a genome and some info to predict the coverage from the genomic sequence and returns a df with this info and the predicted relative coverage"""


    # first get the info that can be calculated from the sequence per se. 
    df_genome_info = get_genome_info_df(genome, replace=replace, window_l=window_l, threads=threads, get_repeat_info=False, with_mappability=False)

    # get metadata
    chrom_to_df = {chrom : df_genome_info[df_genome_info.chromosome==chrom] for chrom in set(df_genome_info.chromosome)}
    chrom_to_len = {s.id : len(s.seq) for s in SeqIO.parse(genome, "fasta")}

    # check that the genome info is correct per window, and if something is missing replace everything
    if any([max(chrom_to_df[chrom]["end"])!=length for chrom, length in chrom_to_len.items()]):
        print("wrong run of get_genome_info_df. This is likely because you interrupted a previous run. Repeating. ...")

        df_genome_info = get_genome_info_df(genome, replace=True, window_l=window_l, threads=threads, get_repeat_info=False, with_mappability=False)
    
    # initialize predicted coverage
    predicted_cov_all = []

    # go through each chromosome and predict coverage from location
    #print("predicting coverage for windows of the genome from location")
    for chromosome in df_genome_info.chromosome.unique():

        # get the df of this chromosome
        df_c = df_genome_info[df_genome_info.chromosome==chromosome]

        # get the predicted cov from location
        predicted_cov = list(poly.polyval(df_c.relative_position_in_chr.values, chromosome_to_relPosVScovCoefficients[chromosome]))
        predicted_cov_all += predicted_cov

    df_genome_info["cov_pred_from_loc"] = predicted_cov_all

    # get to non-0 values
    def get_to_non0_cov(cov, tsh=0.1):

        if cov<tsh: return tsh
        else: return cov

    # adding the residuals from seq features
    residual_pred_from_loc =  lm_predictResidualFromPos.predict(df_genome_info[predictors_ResCovAfterPosCorrection])
    df_genome_info["predicted_relative_coverage"] =  (df_genome_info.cov_pred_from_loc + residual_pred_from_loc).apply(get_to_non0_cov)
    #print("relative_cov added")

    # plot per chrom
    if make_plots is True:

        all_chromosomes = df_genome_info.chromosome.unique()
        #print("making plots")

        fig = plt.figure(figsize=(10, len(all_chromosomes)*4.5))
        for I, chrom in enumerate(all_chromosomes):

            # get df
            df_c = df_genome_info[df_genome_info.chromosome==chrom]

            # print each of the predictions
            ax = plt.subplot(len(all_chromosomes), 1, I+1)
            sns.lineplot(x="start", y="predicted_relative_coverage", data=df_c, linewidth=2, color="blue")

            ax.set_title(chrom)

        fig.tight_layout()  # otherwise the right y-label is slightly 
        filename="%s_predicted_relative_coverage.pdf"%(plots_prefix)
        fig.savefig(filename, bbox_inches='tight');
        plt.close(fig)
        #print("plots finsihed")

    return df_genome_info

def count_number_read_pairs(bamfile, replace=False, threads=4):

    """counts the total number of reads of a bamfile"""

    read_count_file = "%s.flagstat"%bamfile; read_count_file_tmp = "%s.tmp"%read_count_file

    # get the total n reads
    if file_is_empty(read_count_file) or replace is True:

        print("calculating n reads")
        run_cmd("%s flagstat --threads %i %s > %s"%(samtools, threads, bamfile, read_count_file_tmp))
        os.rename(read_count_file_tmp, read_count_file)

    return [int(l.split()[0]) for l in open(read_count_file, "r").readlines() if " read1" in l][0]

def get_read_length(bamfile, threads=4, nreads=5000, replace=False):

    """Calculates the readlength for a bamfile"""

    readlen_dist_file = "%s.read_length_dist_first%ireads.txt"%(bamfile, nreads); readlen_dist_file_tmp = "%s.tmp"%readlen_dist_file
    if file_is_empty(readlen_dist_file) or replace is True:

        print("The following command will throw a warning stating that 'samtools view: writing to standard output failed: Broken pipe'. This is because the output of samtools view is piped, which is expected.")
        cmd = "%s view --threads %i %s | head -n %i | cut -f10 | perl -ne 'chomp;print length($_) . \"\n\"' | sort > %s"%(samtools, threads, bamfile, nreads, readlen_dist_file_tmp)
        run_cmd(cmd)

        os.rename(readlen_dist_file_tmp, readlen_dist_file)

    return int(np.median([int(l.strip()) for l in open(readlen_dist_file, "r").readlines()]))


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

def simulate_and_align_PairedReads_perWindow(df_windows, genome_interest, reference_genome, npairs, read_length, outdir, median_insert_size, median_insert_size_sd, replace=False, threads=4, mitochondrial_chromosome="mito_C_glabrata_CBS138"):

    """Takes a dataframe with windows of the genome, which also has a predicted_relative_coverage (which indicates by how much should the coverage be multiplied in this window). This function generates a fastq (and deletes it at the end), aligns it and returns it. All files are written under outdir. It returns the aligned bamfile. Any circular chromosome can be passed as mitochondrial_chromosome"""

    print(df_windows)

    print("Simulating reads and aliging them for %s. This will concatenate two mtDNAs because it is circular."%genome_interest)
    start_time = time.time()

    # make folder if it does not exist
    make_folder(outdir)

    # define final outputs
    sim_bamfile = "%s/aligned_reads%ipairs_readlen%i_insertsize%i+-%i.bam"%(outdir, npairs, read_length, median_insert_size, median_insert_size_sd)
    sim_sorted_bam = "%s.sorted"%sim_bamfile
    sim_index_bam = "%s.bai"%sim_sorted_bam

    if file_is_empty(sim_sorted_bam) or replace is True:

        ############## CHANGING mtDNA TO SIMULATE IT CIRCULAR ##############

        # edit the df_windows, so that the mitochondrial chromosome gets two concatenated copies, each of them with half of the coverage. This is to simulate that it is circular
        df_windows = df_windows.sort_values(by=["chromosome", "start"])
        df_mito = df_windows[df_windows.chromosome==mitochondrial_chromosome]
        df_nuclear = df_windows[df_windows.chromosome!=mitochondrial_chromosome]

        # change the coverage to half and double
        df_mito["predicted_relative_coverage"] = df_mito.predicted_relative_coverage / 2
        
        # create a mito2 with shifted windows    
        df_mito2 = cp.deepcopy(df_mito)

        print(df_mito2)
        df_mito2["start"] = df_mito2.start + df_mito2.iloc[-1]["end"]
        df_mito2["end"] = df_mito2.end + df_mito2.iloc[-1]["end"]

        # integrate into one df
        df_mito = df_mito.append(df_mito2, sort=True)
        df_windows = df_nuclear.append(df_mito, sort=True)

        # create a genome that has two copies of the mtDNA
        mtDNA_chr = [c for c in SeqIO.parse(genome_interest, "fasta") if c.id==mitochondrial_chromosome][0]
        mtDNA_chr = mtDNA_chr + cp.deepcopy(mtDNA_chr)
        gDNA_chrs = [c for c in SeqIO.parse(genome_interest, "fasta") if c.id!=mitochondrial_chromosome]

        genome_interest_2mtDNA = "%s.2x_mtDNA.fasta"%genome_interest
        SeqIO.write(gDNA_chrs + [mtDNA_chr], genome_interest_2mtDNA, "fasta")

        #######################################################################

        ##### generate the reads for the genome_interest #####
        outdir_reads = "%s/simulated_reads_%ipairs_readlen%i_insertsize%i+-%i"%(outdir, npairs, read_length, median_insert_size, median_insert_size_sd) 
        read1_fastqgz, read2_fastqgz = simulate_readPairs_per_window(df_windows, genome_interest_2mtDNA, npairs, outdir_reads, read_length, median_insert_size, median_insert_size_sd, replace=replace, threads=threads) 
        ######################################################

        ##### align the reads and retun the bam ######
        run_bwa_mem(read1_fastqgz, read2_fastqgz, reference_genome, outdir, sim_bamfile, sim_sorted_bam, sim_index_bam, name_sample="simulations_reference_genome", threads=threads, replace=replace)

        # remove the fastq files
        print("deleting reads")
        delete_folder(outdir_reads)

    # record the time consumed
    print("--- generating %i reads from %s and aligning them  took %s seconds ---"%(npairs, genome_interest, time.time() - start_time))

    return sim_sorted_bam

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

def downsample_bamfile_keeping_pairs(bamfile, fraction_reads=0.1, replace=True, threads=4, name="sampleX"):

    """Takes a sorted and indexed bam and samples a fraction_reads randomly. This is a fast process so that it can be repeated many times"""

    # define the outfile
    seed =  random.choice(list(range(30)))
    sampled_bamfile = "%s.%ipct_reads_seed%i_%s.bam"%(bamfile, int(fraction_reads*100), seed, name); remove_file(sampled_bamfile)

    # count number of reads
    #npairs = count_number_read_pairs(bamfile, replace=replace, threads=threads)

    # get fraction as str
    fraction_str = str(fraction_reads).split(".")[1]

    # run the sampling
    sampled_bamfile_unedited = "%s.unedited"%sampled_bamfile
    run_cmd("%s view -o %s -s %i.%s --threads %i -b %s"%(samtools, sampled_bamfile_unedited, seed, fraction_str, threads, bamfile))

    # now rename the reads, so that the read name gets a prefix called name
    run_cmd("%s view -h --threads %i %s | sed -e 's/^/%s_/' | sed 's/^%s_@/@/' | %s view --threads %i -bSh > %s"%(samtools, threads, sampled_bamfile_unedited, name, name, samtools, threads, sampled_bamfile))

    # at the end remove the unedited
    remove_file(sampled_bamfile_unedited)

    return sampled_bamfile

def merge_2bams(bamA, bamB, merged, threads=4):

    """Merges two bams"""
    remove_file(merged)

    print("merging %s and %s into %s"%(bamA, bamB, merged))
    run_cmd("%s merge --threads %i %s %s %s"%(samtools, threads, merged, bamA, bamB))

def sort_bam(bam, sorted_bam, threads=4):

    """Sorts a bam file into sorted_bam"""

    print("sorting bam")
    run_cmd("%s sort --threads %i -o %s %s"%(samtools, threads, sorted_bam, bam))

def index_bam(bam, threads=4):

    """indexes bam and creates a .bai file"""

    print("indexing bam")
    run_cmd("%s index -@ %i %s"%(samtools, threads, bam))

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

def get_merged_bamfile_for_ploidy(variant_bamfile, reference_bamfile, ploidy, replace=False, threads=4):

    """This function takes two bam files, one that includes aligned reads for a reference genome and another that includes aligned reads for a variat genome. It subsamples each of the bamfiles so that there is a proportion of reads comming from each reference and variant, and the proportion is indicated by ploidy, which can be any of the strings in the default values of target_ploidies of run_GridssClove_optimising_parameters."""

    print("merging %s and %s in ploidy %s"%(variant_bamfile, reference_bamfile, ploidy))

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
    print("--- merging took for ploidy %s reads took %s seconds ---"%(ploidy, time.time() - start_time))

    return merged_sorted_bam

def get_and_write_repeatsDFbed(reference_genome, threads=4, replace=False):

    """Writes a bed with the repeats of a reference genome"""

    repeats_df = get_repeat_maskerDF(reference_genome, threads=threads, replace=replace)

    repeats_df["start_0based"] = repeats_df.begin_repeat - 1
    repeats_df["end_0based"] = repeats_df.end_repeat - 1

    # write
    repeats_bed = "%s.repeatMaskerRegions.bed"%reference_genome
    repeats_df[["chromosome", "start_0based", "end_0based", "strand"]].to_csv(repeats_bed, sep="\t", header=False, index=False)

    return repeats_bed

def get_bed_regions_arround_bed(input_bed, genome, region_size=50, replace=False):

    """Takes a bed and extends each windwo by +- region_size, writing it under input_bed"""

    # define ouptut
    regions_bed = "%s.regions+-%i.bed"%(input_bed, region_size)

    if file_is_empty(regions_bed) or replace is True:


        # define the fields
        fields = ["chromosome", "start", "end", "strand"]
        
        # find the number of fields
        nfields = len(open(input_bed, "r").readlines()[0].strip().split("\t"))
        fields = fields[0:nfields]

        # get bed into df
        df_bed = pd.read_csv(input_bed, sep="\t", header=None, names=fields)

        # two functions to set the regions
        chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(genome, "fasta")}

        # define a function that sets the correct coordinates
        def get_zeros_for_negatives(coord):

            if coord<0: return 0
            else: return coord

        def get_plus_region(r):

            coord = r["end"] + region_size
            len_chr = chr_to_len[r["chromosome"]]
            if coord > len_chr: return len_chr
            else: return coord

        # add the region
        df_bed["start"] = (df_region["start"] - region_size).apply(get_zeros_for_negatives)
        df_bed["end"] = df_bed.apply(get_plus_region, axis=1)

        # write
        df_bed.to_csv(regions_bed, sep="\t", header=False, index=False)

    # write and return
    return regions_bed


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
        print("running gridss")

        # define the out and error of gridss
        gridss_std = "%s/gridss_run_std.txt"%outdir
        
        max_tries = 2
        for Itry in range(max_tries):
            try: 
                # delete previous files
                delete_folder(gridss_tmpdir); make_folder(gridss_tmpdir)
                remove_file(gridss_assemblyBAM)

                # define the heap size, which depends on the cloud or not
                #jvmheap = "27.5g" # this is the default
                jvmheap = "31g" # this works in MN

                # change the jvmheap if you are in the bsc local computer
                try: 
                    if mschikoraParentDir=="/home/mschikora/samba": jvmheap = "20g" # you may want to change this

                except: print("not running in MN")

                # run
                print("running gridss on %s jvmheap"%jvmheap)
                run_cmd("%s --jar %s --reference %s -o %s --assembly %s --threads %i --workingdir %s --maxcoverage %i --blacklist %s --jvmheap %s %s > %s 2>&1"%(gridss_run, gridss_jar, reference_genome, gridss_VCFoutput, gridss_assemblyBAM, threads, gridss_tmpdir, maxcoverage, blacklisted_regions, jvmheap, sorted_bam, gridss_std))

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

def add_info_to_gridssDF(df, expected_fields={"allele_frequency", "allele_frequency_SmallEvent", "other_coordinates", "other_chromosome", "other_position", "other_orientation",  "inserted_sequence", "len_inserted_sequence", "length_event", "has_poly16GC", "length_inexactHomology", "length_microHomology"}, median_insert_size=500, median_insert_size_sd=50):

    """This function takes a gridss df and returns the same adding expected_fields"""

    #print("adding info to gridss df")
    if len(set(df.keys()).intersection(expected_fields))!=len(expected_fields):

        df = cp.deepcopy(df)

        # add the allele frequencies
        #print("adding allele freq")
        df["allele_frequency"] = df.apply(lambda r: np.divide(r["DATA_VF"] , (r["DATA_VF"] + r["DATA_REF"] + r["DATA_REFPAIR"])), axis=1)
        df["allele_frequency_SmallEvent"] = df.apply(lambda r: np.divide(r["DATA_VF"] , (r["DATA_VF"] + r["DATA_REF"])), axis=1)

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

# def get nans to -1
def get_int(x):

    try: return int(x)
    except: return -1


def get_gridssDF_filtered_from_filtersDict(df_gridss, filters_dict):

    """Takes a df gridss and returns the filtered one, according to filters_dict"""

    # debug the fact that there is no min_af_EitherSmallOrLargeEvent
    if "min_af_EitherSmallOrLargeEvent" not in filters_dict: filters_dict["min_af_EitherSmallOrLargeEvent"] = 0.0

    # get the filtered df
    df_filt = get_gridssDF_filtered(df_gridss, min_Nfragments=filters_dict["min_Nfragments"], min_af=filters_dict["min_af"], wrong_INFOtags=filters_dict["wrong_INFOtags"], wrong_FILTERtags=filters_dict["wrong_FILTERtags"], filter_polyGC=filters_dict["filter_polyGC"], filter_noSplitReads=filters_dict["filter_noSplitReads"], filter_noReadPairs=filters_dict["filter_noReadPairs"], maximum_strand_bias=filters_dict["maximum_strand_bias"], maximum_microhomology=filters_dict["maximum_microhomology"], maximum_lenght_inexactHomology=filters_dict["maximum_lenght_inexactHomology"], range_filt_DEL_breakpoints=filters_dict["range_filt_DEL_breakpoints"], min_length_inversions=filters_dict["min_length_inversions"], dif_between_insert_and_del=filters_dict["dif_between_insert_and_del"], max_to_be_considered_small_event=filters_dict["max_to_be_considered_small_event"], min_size=filters_dict["min_size"], add_columns=False, min_af_EitherSmallOrLargeEvent=filters_dict["min_af_EitherSmallOrLargeEvent"] )

    return df_filt

def get_tupleBreakpoints_for_filters_GRIDSS(df_gridss, filters_dict, return_timing=False):

    """ Takes a df_gridss (the output vcf) and a dictionary with filters, returning a tuple of the breakpoints where both breakends have passed the filters."""

    # initialize time
    start_time = time.time()

    # debug the fact that there is no min_af_EitherSmallOrLargeEvent
    if "min_af_EitherSmallOrLargeEvent" not in filters_dict: filters_dict["min_af_EitherSmallOrLargeEvent"] = 0.0

    # get the filtered df
    df_filt = get_gridssDF_filtered(df_gridss, min_Nfragments=filters_dict["min_Nfragments"], min_af=filters_dict["min_af"], wrong_INFOtags=filters_dict["wrong_INFOtags"], wrong_FILTERtags=filters_dict["wrong_FILTERtags"], filter_polyGC=filters_dict["filter_polyGC"], filter_noSplitReads=filters_dict["filter_noSplitReads"], filter_noReadPairs=filters_dict["filter_noReadPairs"], maximum_strand_bias=filters_dict["maximum_strand_bias"], maximum_microhomology=filters_dict["maximum_microhomology"], maximum_lenght_inexactHomology=filters_dict["maximum_lenght_inexactHomology"], range_filt_DEL_breakpoints=filters_dict["range_filt_DEL_breakpoints"], min_length_inversions=filters_dict["min_length_inversions"], dif_between_insert_and_del=filters_dict["dif_between_insert_and_del"], max_to_be_considered_small_event=filters_dict["max_to_be_considered_small_event"], min_size=filters_dict["min_size"], add_columns=False, min_af_EitherSmallOrLargeEvent=filters_dict["min_af_EitherSmallOrLargeEvent"] )

    # get the breakpoints that have both breakends
    correct_breakpoints = tuple(sorted([bp for bp, N in Counter(df_filt.breakpointID).items() if N==2]))

    if return_timing: return  (time.time() - start_time)
    else: return correct_breakpoints


def get_varName(var):
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    return [k for k, v in callers_local_vars if v is var][0]


########## DEFIE GRIDSS FILTERS ORDERED FROM LESS CONSERVATIVE TO MOST CONSERVATIVE ############

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

# define default parameters for gridss filtering
default_filtersDict_gridss = {"min_Nfragments":0, "min_af":0.0, "wrong_FILTERtags":("",), "filter_polyGC":False, "filter_noSplitReads":False, "filter_noReadPairs":False, "maximum_strand_bias":1.0, "maximum_microhomology":1000000000, "maximum_lenght_inexactHomology":100000000, "range_filt_DEL_breakpoints":(0,1), "min_length_inversions":0, "dif_between_insert_and_del":0, "max_to_be_considered_small_event":1, "wrong_INFOtags":('IMPRECISE',), "min_size":50, "min_af_EitherSmallOrLargeEvent":0.0}

################################################################################################

def get_represenative_filtersDict_for_filtersDict_list(filtersDict_list, type_filters="less_conservative"):

    """Takes a lis, each position with a list of filters like passed to get_tupleBreakpoints_for_filters_GRIDSS and returns a representative dict, according to less_conservative"""

    # map a score for each dict
    score_to_dict = {sum([g_filterName_to_filterValue_to_Number[fname][find_nearest(g_filterName_to_filtersList[fname], fvalue)] for fname, fvalue in filtersDict.items()]) : filtersDict for filtersDict in filtersDict_list}

    # get the dict with the min or max score, depedning on the approach
    type_filters_to_getPositionFun = {"less_conservative":min, "most_conservative":max}

    return score_to_dict[type_filters_to_getPositionFun[type_filters](score_to_dict.keys())]


def write_bedpeANDfilterdicts_for_breakpoints(df_bedpe, breakpoints, filterDicts, outdir):

    """Takes a df_bedpe that is already filtered (only has the fields to write) and it writes the breakpoints into outdir. It also writes a series with the most conservative filter set that gace with the filters that have rise to this bedpe"""

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

def keep_relevant_filters_lists_inparallel(filterName_to_filtersList, df_gridss, type_filtering="keeping_all_filters_that_change",  wrong_INFOtags=("IMPRECISE",), min_size=50):

    """Takes a dictionary that maps the filterName to the list of possible filters. It modifies each of the lists in filterName_to_filtersList in a way that only those values that yield a unique set of breakpoints when being applied in the context of a set of unconservative filters are preserved"""

    # define a set of filters that are very unconservative (they take all the breakpoints)
    unconservative_filterName_to_filter = {"min_Nfragments":-1, "min_af":-1, "wrong_FILTERtags":("",), "filter_polyGC":False, "filter_noSplitReads":False, "filter_noReadPairs":False, "maximum_strand_bias":1.1, "maximum_microhomology":1000000000000, "maximum_lenght_inexactHomology":1000000000000, "range_filt_DEL_breakpoints":(0,1), "min_length_inversions":-1, "dif_between_insert_and_del":0, "max_to_be_considered_small_event":1, "wrong_INFOtags":wrong_INFOtags, "min_size":min_size, "min_af_EitherSmallOrLargeEvent":-1}

    # define an unconservative set of breakpoints
    #unconservative_breakpoints = get_tupleBreakpoints_for_filters_GRIDSS(df_gridss, unconservative_filterName_to_filter)
    unconservative_breakpoints = tuple(sorted([bp for bp, N in Counter(df_gridss.breakpointID).items() if N==2]))
    print("There are %i bp in total"%len(unconservative_breakpoints))

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
    inputs_fn = [(df_gridss, fd) for fd in filters_dict_list]

    # pool
    with  multiproc.Pool(multiproc.cpu_count()) as pool:
        bp_tuples_list = pool.starmap(get_tupleBreakpoints_for_filters_GRIDSS, inputs_fn)
        pool.close()

    # map
    #bp_tuples_list = list(map(lambda x: get_tupleBreakpoints_for_filters_GRIDSS(x[0], x[1]), inputs_fn))

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

def write_breakpoints_for_parameter_combinations_and_get_filterIDtoBpoints_gridss(df_gridss, df_bedpe, outdir, range_filtering="theoretically_meaningful", expected_AF=1.0, replace=False):

    """Gets, for a range of filters defined byrange_filtering, a dictionary that maps a string defining all these filters to a df that has the filtered breakpoints (bedpe) . If range_filtering is large, we expect ~13 h to run on 48 cores for the Candida glabrata genome"""

    # define files that will be written at the end of this function
    filtersID_to_breakpoints_file  = "%s/filtersID_to_breakpoints_file.py"%outdir

    print("getting lists of bedpe breakpoints")

    if any([file_is_empty(f) for f in [filtersID_to_breakpoints_file]]) or replace is True:

        ##### WORK WITH DF_BEDPE ########

        # map each brekendID to the breakpointID
        all_bend_IDs = set.union(*df_bedpe.IDs_set)
        bendID_to_bpointID = {bendID : df_bedpe[df_bedpe.IDs_set.apply(lambda IDs_set: bendID in IDs_set)].iloc[0]["name"] for bendID in df_gridss.ID if bendID in all_bend_IDs}

        # get only breakends that are in bendID_to_bpointID, meaning that there are 2 breakends
        df_gridss_twoBreakEnds = df_gridss[df_gridss.ID.isin(bendID_to_bpointID)]
        df_gridss_twoBreakEnds["breakpointID"] = df_gridss_twoBreakEnds.ID.apply(lambda x: bendID_to_bpointID[x])
        df_gridss_twoBreakEnds = df_gridss_twoBreakEnds.set_index("ID")[["INFO_SIMPLE_TYPE", "length_event", "allele_frequency", "allele_frequency_SmallEvent", "DATA_VF", "INFO_misc", "FILTER", "INFO_SB", "length_microHomology", "length_inexactHomology", "len_inserted_sequence", "has_poly16GC", "DATA_SR", "DATA_RP", "breakpointID", "real_AF"]]

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

        else: raise ValueError("%s is not a valid range_filtering parameter, it has to be 'large', 'medium', 'small' or 'single' "%range_filtering)

        # define filters that are always the same
        wrong_INFOtags = ("IMPRECISE",)
        min_size = 50

        # map the filters through a dict
        filterName_to_filtersList = {"min_Nfragments":min_Nfragments_l, "min_af":min_af_l, "wrong_FILTERtags":wrong_FILTERtags_l, "filter_polyGC":filter_polyGC_l, "filter_noSplitReads":filter_noSplitReads_l, "filter_noReadPairs":filter_noReadPairs_l, "maximum_strand_bias":maximum_strand_bias_l, "maximum_microhomology":maximum_microhomology_l, "maximum_lenght_inexactHomology":maximum_lenght_inexactHomology_l, "range_filt_DEL_breakpoints":range_filt_DEL_breakpoints_l, "min_length_inversions":min_length_inversions_l, "dif_between_insert_and_del":dif_between_insert_and_del_l, "max_to_be_considered_small_event":max_to_be_considered_small_event_l, "min_af_EitherSmallOrLargeEvent":min_af_EitherSmallOrLargeEvent_l}

        # edit the filter list, to keep only those that, when applied, change the called breakpoints
        keep_relevant_filters_lists_inparallel(filterName_to_filtersList, df_gridss_twoBreakEnds, type_filtering="keeping_filters_that_yield_uniqueBPs", wrong_INFOtags=wrong_INFOtags, min_size=min_size) # it can also be keeping_all_filters_that_change or keeping_filters_that_yield_uniqueBPs or none

        # initialize objects to store the filtering
        I = 1
        filters_dict_list = []

        # go through each range of filters
        print("generating dictionaries of filters")
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

                                        #print("filter %i"%I)
                                        I+=1

                                        # get the parameters_dict
                                        filters_dict = dict(min_Nfragments=min_Nfragments, min_af=min_af, wrong_INFOtags=wrong_INFOtags, wrong_FILTERtags=wrong_FILTERtags, filter_polyGC=filter_polyGC, filter_noSplitReads=filter_noSplitReads, filter_noReadPairs=filter_noReadPairs, maximum_strand_bias=maximum_strand_bias, maximum_microhomology=maximum_microhomology, maximum_lenght_inexactHomology=maximum_lenght_inexactHomology, range_filt_DEL_breakpoints=range_filt_DEL_breakpoints, min_length_inversions=min_length_inversions, dif_between_insert_and_del=dif_between_insert_and_del, max_to_be_considered_small_event=max_to_be_considered_small_event, min_size=min_size, min_af_EitherSmallOrLargeEvent=min_af_EitherSmallOrLargeEvent)

                                        # keep
                                        filters_dict_list.append(filters_dict)
        
        print("There are %i combinations of parameters"%I)
        
        # first try for some combinations, which will give you the timing 
        times = [get_tupleBreakpoints_for_filters_GRIDSS(df_gridss_twoBreakEnds, filters_dict, return_timing=True) for filters_dict in random.sample(filters_dict_list, min(10, len(filters_dict_list)))]
        ncores = multiproc.cpu_count()
        print("Obtaining the list of tuples of breakpoints will take arround %.2f minutes on %i cores"%(((np.mean(times)*I)/ncores)/60, ncores))

        # obtain the list of tuples for each parameter combintaion
        with  multiproc.Pool(multiproc.cpu_count()) as pool:
            bp_tuples_list = pool.starmap(get_tupleBreakpoints_for_filters_GRIDSS, [(df_gridss_twoBreakEnds, fd) for fd in filters_dict_list])
            pool.close()

        # map each tuple o bp to the dicts of parameters that gave it
        bpointTuple_to_filterDicts = {}
        for bpointTuple, filterDict in zip(bp_tuples_list, filters_dict_list): 
            if len(bpointTuple)>0: bpointTuple_to_filterDicts.setdefault(bpointTuple, []).append(filterDict)
        print("There are %i sets of breakpoints that can be created with %i combinations of parameters"%(len(bpointTuple_to_filterDicts), I))

        # map each tuple pf breakpoints to an ID that will be saved
        bpoints_to_ID = dict(zip(bpointTuple_to_filterDicts.keys(), map(lambda I: "filters_%i"%I, range(len(bpointTuple_to_filterDicts)))))

        # generate under otdir all the breakpoints from df_bedpe
        print("writing bedpefiles")
        bedpe_fields = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]
        df_bedpe = df_bedpe[bedpe_fields]
        
        # run generation
        inputs_function = [(df_bedpe, bpoints, filterDicts, "%s/%s"%(outdir, bpoints_to_ID[bpoints])) for bpoints, filterDicts in bpointTuple_to_filterDicts.items()]

        with multiproc.Pool(multiproc.cpu_count()) as pool:
            pool.starmap(write_bedpeANDfilterdicts_for_breakpoints, inputs_function)
            pool.close()

        # save the map between each filter 
        print("writing files")
        filtersID_to_breakpoints = dict(zip(bpoints_to_ID.values(), bpoints_to_ID.keys()))
        save_object(filtersID_to_breakpoints, filtersID_to_breakpoints_file)

    else: filtersID_to_breakpoints = load_object(filtersID_to_breakpoints_file)

    # return the dataframe with all the parameter combinations and the filter
    return filtersID_to_breakpoints

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

def get_clove_output_with_coverage_forTANDEL(outfile_clove, reference_genome, sorted_bam, replace=False, run_in_parallel=False, delete_bams=False):

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

        else: coverage_df = pd.DataFrame(columns=["chromosome", "end", "length", "mediancov_1", "nocoveragebp_1", "percentcovered_1", "start"])

        # merge
        merged_df = df_clove.merge(coverage_df, how="left", left_on=["#CHROM", "POS", "END"], right_on=["chromosome", "start", "end"], validate="many_to_one")

        # change types of fields
        merged_df["POS"] = merged_df.POS.apply(get_int)
        merged_df["END"] = merged_df.END.apply(get_int)
        merged_df["START"] = merged_df.START.apply(get_int)

        return merged_df 

    else: return pd.DataFrame()


def get_covfilter_cloveDF_row_according_to_SVTYPE(r, maxDELcoverage=2, minDUPcoverage=50):

    # define a function that takes a row of the dataframe and does the filterinf

    if r["SVTYPE"]=="DEL":

        if r["mediancov_1"]<=maxDELcoverage: return "PASS"
        else: return "FAIL" 

    elif r["SVTYPE"]=="TAN":

        if r["mediancov_1"]>=minDUPcoverage: return "PASS"
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

def get_bedpe_for_clovebalTRA_5with3(r, chr_to_len):

    """Takes a row of the df_balTRA_5with3 df and returns a bedpe row, sorted"""

    return pd.Series({"ChrA":r["#CHROM"], "StartA":0, "EndA":r["POS"], "ChrB":r["CHR2"], "StartB":r["END"], "EndB":chr_to_len[r["CHR2"]]-1, "Balanced":True})

def set_float_to_maximum(val, maximum=1.0):

    if val>maximum: return maximum
    else: return val

def calculate_coverage_df_clove_regions_5end_3end(df_clove, fileprefix, sorted_bam, reference_genome, window_l, replace=False, run_in_parallel=True):

    """Takes a df of clove, i.e. the one for ITX and INVITX and returns the regions with coverage as a df, writing under fileprefix"""

    # map chr_to_len
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

    # make a df where you have the 5' and 3' regions of each #CHROM,POS or CHR2,END positions
    regionID_to_field_to_value = {}; I=0
    for CHROM, POS, CHR2, END, ID in df_clove[["#CHROM", "POS", "CHR2", "END", "ID"]].values:

        # go through each chromosome
        chromType_to_info  = {"#CHROM":{"chromosome":CHROM, "position":POS}, "CHR2":{"chromosome":CHR2, "position":END}}
        for chromType, info in chromType_to_info.items():
            chromosome = info["chromosome"]
            bp = info["position"]

            # go through each region
            for region in {"5end", "3end"}:

                if region=="5end": start=0; end=bp
                else: start=bp; end=chr_to_len[chromosome]

                regionID_to_field_to_value[I] = {"chromType":chromType, "chromosome":chromosome, "start":start, "end":end, "region":region, "breakpoint_ID":ID}

                # at the end update I
                I+=1

    # get as df
    df = pd.DataFrame(regionID_to_field_to_value).transpose()

    # write as bed, all the files that have not 0 and 0
    bed_windows = "%s.regions.bed"%fileprefix
    bed_fields = ["chromosome", "start", "end"]
    df[~((df.start==0) & (df.end==0))][bed_fields].sort_values(by=bed_fields).to_csv(bed_windows, sep="\t", header=True, index=False)
    
    # get coverage dataframe
    outdir = "/".join(fileprefix.split("/")[0:-1])
    coverage_df = get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, bed_windows, replace=replace, run_in_parallel=run_in_parallel, delete_bams=run_in_parallel)

    # merge with the df
    df_all = df.merge(coverage_df, left_on=bed_fields, right_on=bed_fields, how="left", left_index=False, right_index=False)
    
    # remove files
    for f in [bed_windows]: remove_file(f)

    return df_all

def get_bedpe_for_clove_unbalTRA(df_clove, chr_to_len, outdir, reference_genome, sorted_bam, replace=False, median_coverage=10, threshold_p=0.1, run_in_parallel=False):

    """Takes a dataframe of CLOVE and returns the unbalanced translocations, where one part replaces another. The idea is that any ITX* or INVITX* marks a possible region with imbalanced translocation. This function uses coverage to decide whether it is  or not an umbalanced translocation."""

    # first define a df with possible umbalanced translocations
    df_INVTX_ITX = df_clove[df_clove.SVTYPE.isin({"ITX1", "IXT2", "INVTX1", "INVTX2"})]
    all_breakoint_IDs = set(df_INVTX_ITX.ID)
    if len(df_INVTX_ITX)==0: return pd.DataFrame()

    # get the df of this regions if not provided
    coverage_df_INVTX_ITX = calculate_coverage_df_clove_regions_5end_3end(df_INVTX_ITX, "%s/INVTX_ITX_coverage"%outdir, sorted_bam, reference_genome, 1000, replace=replace, run_in_parallel=run_in_parallel)
    
    # chech that all the breakpoints are defined
    if len(all_breakoint_IDs)>len(all_breakoint_IDs.intersection(set(coverage_df_INVTX_ITX.breakpoint_ID))): raise ValueError("Not all the breakpoints are in the provided coverage_df_INVTX_ITX")

    # get the IDs of interest
    df_all = coverage_df_INVTX_ITX[coverage_df_INVTX_ITX.breakpoint_ID.isin(set(all_breakoint_IDs))].merge(df_INVTX_ITX[["ID", "SVTYPE"]], left_on="breakpoint_ID", right_on="ID", validate="many_to_one")

    # add the relative coverage
    df_all["relative_coverage"] = df_all.mediancov_1 / median_coverage

    # map each svtype to possibly exchanged regions
    svtype_to_possiblyChangingRegions = {"ITX1": [("5end", "5end"), ("3end", "3end")], # ITX* mean that one end replaces the same end
                                         "ITX2": [("5end", "5end"), ("3end", "3end")], # ITX* mean that one end replaces the same end
                                         "INVTX1": [("5end", "3end"), ("3end", "5end")], # INVTX* means that a 5' replaces a 3' or viceversa
                                         "INVTX2": [("5end", "3end"), ("3end", "5end")], # INVTX* means that a 5' replaces a 3' or viceversa
                                         } 

    # go through each breakpoint, and add to bedpe dict
    bedpe_dict = {}
    for ID in set(df_all.breakpoint_ID):

        # get a df for this event
        df = df_all[df_all.breakpoint_ID==ID].set_index(["region", "chromType"], drop=False)

        # svtype
        svtype = df.SVTYPE.iloc[0]

        # initialize hypothesis dict
        hypothesis_to_p = {} # hypothesis will be (ChrA, StartA, EndA, ChrB, StartB, EndB), and the p is how distant it is from duplicating A, and deleting B

        # go through each region
        for regions_exchanged in svtype_to_possiblyChangingRegions[svtype]:

            # go through each chromosome type
            for ChrA_type, ChrB_type in [("#CHROM", "CHR2"), ("CHR2", "#CHROM")]:

                # get the series of the region
                regionA = df.loc[(regions_exchanged[0], ChrA_type)]
                regionB = df.loc[(regions_exchanged[1], ChrB_type)]

                # calculate the probability of the hypothesis
                p_duplication = set_float_to_maximum(regionA["relative_coverage"] - 1) # the higher the better, up to 1
                p_deletion = 1 - regionB["relative_coverage"] # the lower the better

                # keep the probability
                p = np.mean([p_duplication, p_deletion])
                if p>=threshold_p: hypothesis_to_p[(regionA.chromosome, regionA.start, regionA.end, regionB.chromosome, regionB.start, regionB.end, ID)] = p

        # take the best hypothesis
        if len(hypothesis_to_p)>0: bedpe_dict[ID] = max(hypothesis_to_p.items(), key=(lambda x: x[1]))[0]

    # convert to df
    if len(bedpe_dict)>0:
        df_unbal_TRA = pd.DataFrame(bedpe_dict).transpose().rename(columns={0:"ChrA", 1:"StartA", 2:"EndA", 3:"ChrB", 4:"StartB", 5:"EndB", 6: "breakpoint_ID"})
        df_unbal_TRA["Balanced"] = [False]*len(df_unbal_TRA)

        # get the infex so that it resembles the one of CLOVE
        df_unbal_TRA.index = [df_INVTX_ITX[df_INVTX_ITX.ID==ID].index[0] for ID in df_unbal_TRA.breakpoint_ID]
        return df_unbal_TRA

    else: return pd.DataFrame()
      
def write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, fileprefix, reference_genome, sorted_bam, tol_bp=50, replace=False, median_coverage=10, svtypes_to_consider={"insertions", "deletions", "inversions", "translocations", "tandemDuplications", "remaining"}, threshold_p_unbalTRA=0.1, run_in_parallel=False):

    """Takes a clove dataframe and writes the different SV into several files, all starting with fileprefix. it returns a dict mapping each SVtype to the file with the bed or bedpe containing it. tol_bp indicates the basepairs that are considered as tolerated to be regarded as 'the same event' 

    consider_TANDEL indicates whether to write, it requires the coverage_FILTER to PASS """

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
    cloveSVtypes_easy_classification = {"CID", "CIT", "DUP", "TRA", "CIV", "IVD", "DEL", "TAN"}

    ###############################


    ###### TRANSLOCATIONS:  A segment from the 5’ or 3’ end of one chromosome A is exchanged with the 5’ or 3’ end of another chromosome B #######
    if any([x in set(df_clove.SVTYPE) for x in {"ITX1", "ITX2", "IVD", "INVTX1", "INVTX2"}]) and "translocations" in svtypes_to_consider:

        # balanced translocations 5with5
        df_balTRA_5with5_or_3with3 = get_bedpeDF_for_clovebalTRA_5with5_or_3with3(df_clove, tol_bp=tol_bp) # here the index is not balanced
        considered_idxs += make_flat_listOflists(df_balTRA_5with5_or_3with3.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

        # balanced translocations 5with3 (these are the ones with an IVD field, assigned by clove. This may fail sometimes)
        df_balTRA_5with3 = df_clove[(df_clove.SVTYPE=="IVD") & ((df_clove.START - df_clove.END)<=tol_bp)].apply(lambda r: get_bedpe_for_clovebalTRA_5with3(r, chr_to_len), axis=1)
        considered_idxs += list(df_balTRA_5with3.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

        # unbalanced translocations
        outdir_unbalTRA = "/".join(fileprefix.split("/")[0:-1])
        df_unbal_TRA = get_bedpe_for_clove_unbalTRA(df_clove, chr_to_len, outdir_unbalTRA, reference_genome, sorted_bam, replace=replace, median_coverage=median_coverage, threshold_p=threshold_p_unbalTRA, run_in_parallel=run_in_parallel)
        considered_idxs += list(df_unbal_TRA.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

        # merge together and add some fields
        important_fields = ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Balanced"]
        translocations_dfs =  [df_balTRA_5with5_or_3with3, df_balTRA_5with3, df_unbal_TRA]
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
        df_ins = df_clove[df_clove.SVTYPE.isin({"CID", "CIT", "DUP", "TRA"})][["#CHROM", "POS", "CHR2", "START", "END", "ID", "SVTYPE"]]

        # rename
        df_ins = df_ins.rename(columns={"#CHROM":"ChrB", "POS":"StartB", "CHR2":"ChrA", "START":"StartA", "END":"EndA"})
        
        # add the end. It is formated in a way that the insertion length is equivalent to the inserted fragment
        df_ins["EndB"] = df_ins.StartB + (df_ins.EndA - df_ins.StartA)

        # add whether it is copied
        svtype_to_isCopied = {"CID":"TRUE", "CIT":"FALSE", "DUP":"TRUE", "TRA":"FALSE"}
        df_ins["Copied"] = df_ins.SVTYPE.apply(lambda x: svtype_to_isCopied[x])

        # write as bedpe
        important_fields = ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Copied", "ID"]
        bedpe_insertions = "%s.insertions.bedpe.withCopiedINFO"%fileprefix
        df_ins[important_fields].to_csv(bedpe_insertions, sep="\t", header=True, index=False)

        # keep
        svtype_to_svfile["insertions"] = bedpe_insertions
        considered_idxs += list(df_ins.index); df_clove = df_clove.loc[set(df_clove.index).difference(set(considered_idxs))]

        print("There are %i insertions"%len(df_ins))

    ############################

    # write the remaining events which are not easily assignable
    df_noTANDEL = df_clove[~(df_clove.SVTYPE.isin(cloveSVtypes_easy_classification))]
    remaining_noTANDELfile = "%s.remaining_threshold_p_unbalTRA=%.2f.tab"%(fileprefix, threshold_p_unbalTRA)
    df_noTANDEL[["ID", "#CHROM", "POS", "CHR2", "START", "END", "SVTYPE"]].to_csv(remaining_noTANDELfile, sep="\t", header=True, index=False)
    svtype_to_svfile["remaining"] = remaining_noTANDELfile

    print("There are %i remaining SVs"%len(df_noTANDEL))


    ###### DEL and TAN. This is done at the end because some TAN and DEL are filtered before and included to be   #######

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


    # at the end make sure that the considered idxs are unique
    if len(considered_idxs)!=len(set(considered_idxs)): 
        print(fileprefix, considered_idxs)
        #raise ValueError("ERROR: Some clove events are assigned to more than one cathegory. Check the insertions and translocations calling")
        print("WARNING: Some clove events are assigned to more than one cathegory. Check the insertions and translocations calling")

    # return the df_clove and the remaining SVs
    return df_clove, svtype_to_svfile

def analyse_benchmarking_comparing_knownANDpredicted_inHouse(df_benchmark, svtype_to_predsvfile, know_SV_dict, fileprefix, unassigned_BPs_name="remaining", tolerance_bp=50):

    """This function runs an analysis of benchmarking by merging rows of FP, FN and uncalled variants. It is useful to understand if there is any pattern. Everything is written under fileprefix.
    We add to each of the dataframes (false negatives df, false positives df and unassigned df) three fields:

    # coordinates are sets of tuples of (chr,pos). Any event that is united by one of this location is considered to be related, within a window of tolerance_bp
    # string_rep is something that identifies the row.
    # unique_id

    Which are used for the overlap analysis.
    """
    
    # map cmplex events to a boolean field
    complexEvt_to_boolField = {"translocations":"Balanced", "insertions":"Copied"}

    # initiallize a dict that will map each data type to the corresponding dict
    dataType_to_df = {}

    # a function that makes a string representation of a row
    def get_row_string_rep(r, fields): return ";".join(["%s=%s"%(f, r[f]) for f in fields])

    # get a df of the uninitialized events
    df_unassigned = pd.read_csv(svtype_to_predsvfile[unassigned_BPs_name], sep="\t")
    if len(df_unassigned)>0:
        df_unassigned["coordinates"] = df_unassigned.apply(lambda r: {(r["#CHROM"], r["POS"]), (r["CHR2"], r["START"]), (r["CHR2"], r["END"]),}, axis=1)
        df_unassigned["string_rep"] = df_unassigned.apply(lambda r: get_row_string_rep(r, [x for x in r.keys() if x not in {"coordinates", "string_rep"}]), axis=1)
        dataType_to_df["unassignedSV"] = df_unassigned

    # load dataType_to_df for each svtype's FPs and FNs
    #print("getting string representations of dfs")
    for svtype, predfile in svtype_to_predsvfile.items():
        if svtype not in know_SV_dict: continue

        # get dfs of the variants
        df_predicted = pd.read_csv(predfile, sep="\t")
        df_known = pd.read_csv(know_SV_dict[svtype], sep="\t")

        # define ID fields
        if "Name" in df_known.keys(): known_IDfield = "Name"
        else: known_IDfield = "ID"

        predicted_IDfield = "ID"

        # get the IDs of the wrong dfs
        df_benchmark_svtype = df_benchmark[df_benchmark.svtype==svtype]
        FP_IDs = set.union(*df_benchmark_svtype.false_positives_predictedIDs.apply(lambda x: set(str(x).split("||"))))
        FN_IDs = set.union(*df_benchmark_svtype.false_negatives_knownIDs.apply(lambda x: set(str(x).split("||"))))

        # get the dfs of the missing data
        df_FP = df_predicted[df_predicted[predicted_IDfield].isin(FP_IDs)]
        df_FN = df_known[df_known[known_IDfield].isin(FN_IDs)]

        # define the fields that are represnetative depending on the type of sv
        if svtype in {"inversions", "tandemDuplications", "deletions"}:
            coords_fields = [["Chr", "Start"], ["Chr","End"]]
            string_rep_fields_withoutID = ["Chr", "Start", "End"]

        elif svtype in {"translocations", "insertions"}: 
            coords_fields = [["ChrA", "StartA"], ["ChrA","EndA"], ["ChrB", "StartB"], ["ChrB","EndB"]]
            string_rep_fields_withoutID =  ["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", complexEvt_to_boolField[svtype]]

        elif svtype in {"remaining"}:
            coords_fields = [["#CHROM", "POS"], ["CHR2", "START"], ["CHR2", "END"]]
            string_rep_fields_withoutID = ["#CHROM", "POS", "CHR2", "START", "END", "SVTYPE"]

        else: raise ValueError("%s is not a valid svtype to benchmark"%svtype)

        # go through each df and keep
        for eventType, (df, IDfield) in {"FP":[df_FP, predicted_IDfield], "FN":[df_FN, known_IDfield]}.items():

            # if it is empty
            if len(df)>0:

                # add the coordinates and string rep
                df["coordinates"] = df.apply(lambda r: set([(r[cf[0]], r[cf[1]]) for cf in coords_fields]), axis=1)
                df["string_rep"] = df.apply(lambda r: get_row_string_rep(r, string_rep_fields_withoutID + [IDfield]), axis=1)

                # keep df
                dataType_to_df["%s_%s"%(svtype, eventType)] = df

    #print("writing dfs, each with the string representation of the other dfs")

    # make a copy of the dict
    dataType_to_df_copy = cp.deepcopy(dataType_to_df)

    # define a function that writes the overlapping calls in all the other dfs
    def get_overlapping_calls(coords, query_dataType):

        """Takes the corrds of a df and tries to find overlapping calls in all dataType_to_df_copy. coords is a set of (Chr, position) tuples"""

        # initialize calls
        overlapping_calls = ""

        # find overlapping calls in other dfs
        for dataType, df in dataType_to_df_copy.items(): # this is to find overlapping calls

            # discard if they are from the same df
            if dataType==query_dataType: continue

            # find the row overlapping, unles the coord has 0, 1 or -1
            df_overlapping = df[df.coordinates.apply(lambda x: any([any([c[0]==coord[0] and abs(c[1]-coord[1])<=tolerance_bp for c in x]) and coord[1] not in {-1,1,0} for coord in coords]))]

            # add to call
            if len(df_overlapping)>0: overlapping_calls += " ".join(["%s_%s"%(dataType, x) for x in df_overlapping.string_rep])

        return overlapping_calls

    # iterate through the dict twice, writing under fileprefix, but adding any intersection with the other dfs
    for dataType_write, df_write in dataType_to_df.items(): # The first iteration is to write

        # initialize a field that will store the overlapping calls
        df_write["overlapping_calls"] = df_write.coordinates.apply(lambda x: get_overlapping_calls(x, dataType_write))

        # write
        important_fields = [x for x in df_write.keys() if x not in {"coordinates", "string_rep"}]
        filename = "%s.benchmark_analysis_%s.tab"%(fileprefix, dataType_write)
        df_write[important_fields].to_csv(filename, sep="\t", header=True, index=False)

all_svs = {'translocations', 'insertions', 'deletions', 'inversions', 'tandemDuplications', 'remaining'}
def get_integrated_benchmarking_fields_series_for_setFilters_df(df):

    """This function takes a grouped per-filterSet df and returns a row with all the integrated accuracy measurements. The filters of gridss that are best for each SV may vary. If so we will take the least conservative filters of all of the fiters that are best for each SV."""

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
        integrated_benchmarking_results_dict["threshold_p_unbalTRA"] = df_best_filters.loc["translocations", "threshold_p_unbalTRA"]

        integrated_benchmarking_results_dict["median_insert_size"] = df_best_filters.loc["deletions", "median_insert_size"]
        integrated_benchmarking_results_dict["median_insert_size_sd"] = df_best_filters.loc["deletions", "median_insert_size_sd"]
        integrated_benchmarking_results_dict["sorted_bam"] = df_best_filters.loc["deletions", "sorted_bam"]
        integrated_benchmarking_results_dict["median_coverage"] = df_best_filters.loc["deletions", "median_coverage"]

    return pd.Series(integrated_benchmarking_results_dict)

def get_SVbenchmark_dict(df_predicted, df_known, equal_fields=["Chr"], approximate_fields=["Start", "End"], tolerance_bp=50):

    """Takes dfs for known and predicted SVs and returns a df with the benchmark. approximate_fields are fields that have to overlap at least by tolerance_bp. It returns a dict that maps each of the benchmark fields to the value """

    # define the ID fields 
    if "Name" in df_known.keys(): known_IDfield = "Name"
    else: known_IDfield = "ID"

    predicted_IDfield = "ID"

    # get the predictedIDs as those that have the same equal_fields and overlap in all approximate_fields
    if len(df_predicted)>0: df_known["predictedSV_IDs"] = df_known.apply(lambda rk: set(df_predicted[df_predicted.apply(lambda rp: all([rp[f]==rk[f] for f in equal_fields]) and all([abs(rp[f]-rk[f])<=tolerance_bp for f in approximate_fields]), axis=1)][predicted_IDfield]), axis=1)
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

def benchmark_processedSVs_against_knownSVs_inHouse(svtype_to_predsvfile, know_SV_dict, fileprefix, replace=False, analysis_benchmarking=False, tolerance_bp=50, add_integrated_benchmarking=False):

    """Takes two dictionaries that map some SVfiles. It runs, for all the types in svtype_to_predsvfile, a benchmarking against the known ones, writing a file under fileprefix. It returns a df of this benchmark, created with functions written here. It returns as matching events those that have an overlap of at least 50 bp.

    know_SV_dict marks the expected events. If they do not exist you have 0 accuracy.

    add_integrated_benchmarking indicates whether to perform a global benchmarking (not only per svtype)"""

    # map cmplex events to a boolean field
    complexEvt_to_boolField = {"translocations":"Balanced", "insertions":"Copied"}

    # initialize benchmark dict
    benchmark_dict = {}

    # go through each type of event
    for svtype in know_SV_dict:
        print("benchmarking %s"%svtype)

        # load dataframes
        if svtype in svtype_to_predsvfile.keys(): 

            print(svtype_to_predsvfile[svtype])

            df_predicted = pd.read_csv(svtype_to_predsvfile[svtype], sep="\t")
        else: df_predicted = pd.DataFrame()

        df_known = pd.read_csv(know_SV_dict[svtype], sep="\t")

        # define the fields that have to be compared
        if svtype in {"inversions", "tandemDuplications", "deletions"}:
            equal_fields = ["Chr"]
            approximate_fields = ["Start", "End"]

        elif svtype in {"translocations", "insertions"}: 
            equal_fields = ["ChrA", "ChrB", complexEvt_to_boolField[svtype]]
            approximate_fields = ["StartA", "EndA", "StartB", "EndB"]

        elif svtype in {"remaining"}:
            equal_fields = ["#CHROM", "CHR2", "SVTYPE"]
            approximate_fields = ["POS", "START", "END"]

        else: raise ValueError("%s is not a valid svtype"%svtype)

        # get the dict of the benchmark
        dict_benchmark_svtype = get_SVbenchmark_dict(df_predicted, df_known, equal_fields=equal_fields, approximate_fields=approximate_fields, tolerance_bp=tolerance_bp)
        dict_benchmark_svtype["svtype"] = svtype

        # keep
        benchmark_dict[svtype] = dict_benchmark_svtype

    # get the benchmarking
    df_benchmark = pd.DataFrame(benchmark_dict).transpose()

    # analyze the benchmarking
    if analysis_benchmarking is True: analyse_benchmarking_comparing_knownANDpredicted_inHouse(df_benchmark, svtype_to_predsvfile, know_SV_dict, fileprefix, unassigned_BPs_name="remaining", tolerance_bp=tolerance_bp)

    ###### perform integrated benchmarking #####
    if add_integrated_benchmarking is True:

        # get a per-filt row
        integrated_df_benchmark = pd.DataFrame({"integrated": get_integrated_benchmarking_fields_series_for_setFilters_df(df_benchmark)}).transpose()

        # keep
        df_benchmark = df_benchmark.append(integrated_df_benchmark[list(df_benchmark.keys())])
    
    ###########################################

    return df_benchmark

def benchmark_bedpe_with_knownSVs(bedpe, know_SV_dict, reference_genome, sorted_bam, median_coverage, replace=False, ID_benchmark="defaultID", delete_intermediate_files=True):

    """Takes the full path to a bedpe file and generates files, under the same directory, that indicate the benchmarking."""

    # write files under the bedpe outdir
    outdir = "/".join(bedpe.split("/")[0:-1])
    bedpe_filename = bedpe.split("/")[-1]

    # get the benchmark file
    benchmark_df_filename = "%s/df_benchmarking_allParms.py"%outdir
    #remove_file(benchmark_df_filename) # debug

    if file_is_empty(benchmark_df_filename) or replace is True:
    #if True: # debug

        #print("running clove in %s"%outdir)
        print("benchmarking")

        # first run clove without checking for coverage deviations
        outfile_clove = "%s.clove.vcf"%(bedpe)
        run_clove_filtered_bedpe(bedpe, outfile_clove, sorted_bam, replace=replace, median_coverage=10, median_coverage_dev=1, check_coverage=False)

        # now convert it to a df that has also the coverage for TANDEL REGIONS
        df_clove = get_clove_output_with_coverage_forTANDEL(outfile_clove, reference_genome, sorted_bam, replace=replace, run_in_parallel=False)

        # initialize benchmark_df
        df_benchmark_all = pd.DataFrame()

        ##### BENCHMARK INSERTIONS AND INVERSIONS ##########
        #print("benchmarking insertions and inversions variants")
        
        # get files in a way that it is similar to RSVSim, only for complex variants
        fileprefix = "%s/insertionsANDinversions"%(outdir)
        remaining_df_clove, svtype_to_predsvfile = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, fileprefix, reference_genome, sorted_bam, replace=replace, median_coverage=median_coverage, svtypes_to_consider={"insertions", "inversions"})
        svtype_to_predsvfile_INSINV = {svtype : svtype_to_predsvfile[svtype] for svtype in {"insertions", "inversions", "remaining"} if svtype in svtype_to_predsvfile}
        know_SV_dict_INSINV = {svtype : know_SV_dict[svtype] for svtype in {"insertions", "inversions", "remaining"} if svtype in know_SV_dict}

        # benchmark (and write missing events with overlaps)
        df_benchmark_INSINV = benchmark_processedSVs_against_knownSVs_inHouse(svtype_to_predsvfile_INSINV, know_SV_dict_INSINV, fileprefix, replace=replace, analysis_benchmarking=True, tolerance_bp=50)
        
        # add the field of the TRA benchmark that is necessary
        df_benchmark_INSINV["threshold_p_unbalTRA"] = [-1]*len(df_benchmark_INSINV)

        # keep
        df_benchmark_all = df_benchmark_all.append(df_benchmark_INSINV, sort=True)

        #####################################################

        ####### BENCHMARK TRANSLOCATIONS #########
        #print("benchmarking translocations variants")

        # go through several proabilities
        for threshold_p_unbalTRA in [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:

            # get files in a way that it is similar to RSVSim, only for complex variants
            fileprefix = "%s/translocations"%(outdir)
            remaining_df_clove, svtype_to_predsvfile = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, fileprefix, reference_genome, sorted_bam, replace=replace, median_coverage=median_coverage, svtypes_to_consider={"translocations"}, threshold_p_unbalTRA=threshold_p_unbalTRA)
            svtype_to_predsvfile_TRA = {svtype : svtype_to_predsvfile[svtype] for svtype in {"translocations", "remaining"} if svtype in svtype_to_predsvfile}
            know_SV_dict_TRA = {svtype : know_SV_dict[svtype] for svtype in {"translocations", "remaining"} if svtype in know_SV_dict}

            # benchmark (and write missing events with overlaps)
            df_benchmark_TRA = benchmark_processedSVs_against_knownSVs_inHouse(svtype_to_predsvfile_TRA, know_SV_dict_TRA, fileprefix, replace=replace, analysis_benchmarking=True, tolerance_bp=50)
            
            # add the p used
            df_benchmark_TRA["threshold_p_unbalTRA"] = [threshold_p_unbalTRA]*len(df_benchmark_TRA)

            # keep 
            df_benchmark_all = df_benchmark_all.append(df_benchmark_TRA, sort=True)

        # add fields that are necessary to compare with the TANDEL df
        df_benchmark_all["clove_max_rel_coverage_to_consider_del"] = [-1]*len(df_benchmark_all)
        df_benchmark_all["clove_min_rel_coverage_to_consider_dup"] = [-1]*len(df_benchmark_all)

        #########################################

        # initialize a benchmark df for tandel
        df_benchmark_TANDEL = pd.DataFrame()

        ##### deletions #################
        #print("benchmarking DEL")

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
            if len(df_DEL)>0: df_DEL["coverage_FILTER"] = df_DEL.apply(lambda r: get_covfilter_cloveDF_row_according_to_SVTYPE(r, maxDELcoverage=maxDELcoverage, minDUPcoverage=-1), axis=1)
            else: df_DEL["coverage_FILTER"] = []

            # get a dict svtype_to_svfile
            fileprefix = "%s/%s"%(outdir, coveragefiltID)
            remaining_df_clove, svtype_to_predsvfile = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_DEL, fileprefix, reference_genome, sorted_bam, replace=replace, svtypes_to_consider={"deletions"})

            # benchmark
            know_SV_dict_DEL = {svtype : know_SV_dict[svtype] for svtype in {"deletions"} if svtype in know_SV_dict}

            df_benchmark = benchmark_processedSVs_against_knownSVs_inHouse(svtype_to_predsvfile, know_SV_dict_DEL, fileprefix, replace=replace, analysis_benchmarking=False, tolerance_bp=50)
            df_benchmark["clove_max_rel_coverage_to_consider_del"] = [max_rel_coverage_to_consider_del]*len(df_benchmark)
            df_benchmark["clove_min_rel_coverage_to_consider_dup"] = [-1]*len(df_benchmark)

            # keep
            df_benchmark_TANDEL = df_benchmark_TANDEL.append(df_benchmark, sort=True)

        ################################

        ##### tandem duplications ######
        #print("benchmarking TAN")

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
            if len(df_TAN)>0: df_TAN["coverage_FILTER"] = df_TAN.apply(lambda r: get_covfilter_cloveDF_row_according_to_SVTYPE(r, maxDELcoverage=1000000, minDUPcoverage=minDUPcoverage), axis=1)
            else: df_TAN["coverage_FILTER"] = []

            # get a dict svtype_to_svfile
            fileprefix = "%s/%s"%(outdir, coveragefiltID)
            remaining_df_clove, svtype_to_predsvfile = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_TAN, fileprefix, reference_genome, sorted_bam, replace=replace, svtypes_to_consider={"tandemDuplications"})

            # benchmark
            know_SV_dict_TAN = {svtype : know_SV_dict[svtype] for svtype in {"tandemDuplications"} if svtype in know_SV_dict}

            df_benchmark = benchmark_processedSVs_against_knownSVs_inHouse(svtype_to_predsvfile, know_SV_dict_TAN, fileprefix, replace=replace, analysis_benchmarking=False, tolerance_bp=50)
            df_benchmark["clove_max_rel_coverage_to_consider_del"] = [-1]*len(df_benchmark)
            df_benchmark["clove_min_rel_coverage_to_consider_dup"] = [min_rel_coverage_to_consider_dup]*len(df_benchmark)

            # keep
            df_benchmark_TANDEL = df_benchmark_TANDEL.append(df_benchmark, sort=True)

        ################################

        # add things from the other dfs
        df_benchmark_TANDEL["threshold_p_unbalTRA"] = [-1]*len(df_benchmark_TANDEL)

        # keep
        df_benchmark_all = df_benchmark_all.append(df_benchmark_TANDEL, sort=True)

        # append both and return
        df_benchmark_all["benchmarkID"] = [ID_benchmark]*len(df_benchmark_all)
        df_benchmark_all["bedpe"] = [bedpe]*len(df_benchmark_all) # I changed this line at some point because it was raising integers

        # add the parameters that yielded this dict
        filters_dict = load_object("%s/less_conservative_filtersDict.py"%outdir)
        df_benchmark_all["filters_dict"] = [filters_dict]*len(df_benchmark_all)

        # save in disk
        print("saving into %s"%benchmark_df_filename)
        save_object(df_benchmark_all, benchmark_df_filename)

    else: df_benchmark_all = load_object(benchmark_df_filename)

    # delete intermediate fields
    if delete_intermediate_files is True:

        filenames_to_keep = {bedpe_filename, "%s.clove.vcf.TANDEL.bed.coverage_provided_windows.tab"%bedpe_filename, "unbalanced_translocations_5with5_or_3with3.bed.coverage_provided_windows.tab", "%s.clove.vcf"%bedpe_filename, "uniform_filters_series.py", "variable_filters_df.py", "df_benchmarking_allParms.py", "less_conservative_filtersDict.py"}

        #for file in os.listdir(outdir):
        #    if file not in filenames_to_keep and "benchmark_analysis_" not in file: remove_file("%s/%s"%(outdir, file))

    return df_benchmark_all

def makePlots_gridsss_benchmarking_oneGenome(df_benchmark, PlotsDir, plots={"histogram", "scatter_PRvsRC", "scatter_PRvsRCa_eachSVtype", "Fscore_correlation_scatter", "Fscore_correlation_mat"}):

    """Takes a dataframe such as the output of benchmark_GridssClove_for_knownSV and writes plots under PlotsDir, as specified in plots. These are several """

    print("performing plots into %s"%PlotsDir)

    make_folder(PlotsDir)

    # map each svtype to a marker
    svtype_to_color = {"tandemDuplications": "gray", "deletions": "black", "inversions": "blue", "translocations": "olive", "insertions": "red", "remaining":"magenta"}
    svtype_to_marker = {"tandemDuplications": "P", "deletions": "s", "inversions": "^", "translocations": "D", "insertions": "o", "remaining":"v"}

    # define things
    all_events = {'inversions', 'translocations', 'deletions', 'tandemDuplications', 'insertions'}

    # only consider the benchmarking types if there are events
    all_events = [e for e in {'inversions', 'translocations', 'deletions', 'tandemDuplications', 'insertions', 'remaining'} if len(df_benchmark[df_benchmark.svtype==e])>0]
    #print(all_events)

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
            for svtype in set(df_benchmark.svtype): sns.distplot(list(df_benchmark[df_benchmark.svtype==svtype][field]), hist=False, kde=True, rug=False, color=svtype_to_color[svtype], label=svtype, kde_kws=dict(linewidth=3))

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

        #print(df_benchmark)

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

            sns.scatterplot(x="recall", y="precision", data=df_benchmark[df_benchmark.svtype==svtype], hue="Fvalue", palette=cmap, edgecolors=None, style="svtype", markers=svtype_to_marker)
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

def merge_tables_into_file(list_table_files, outfile):

    """Takes a list of table files and merges them into outfile"""
    
    if len(list_table_files)>0:

        df = pd.concat([pd.read_csv(x, sep="\t") for x in list_table_files if os.path.isfile(x)], sort=True)

        # check that the outilfe is not in the list_table_files
        if outfile in list_table_files: outfile += ".%s"%(id_generator(25))
        df.to_csv(outfile, sep="\t", header=True, index=False)

        # remove previously generated files
        for f in list_table_files: remove_file(f)


def benchmark_GridssClove_for_knownSV(sample_bam, reference_genome, know_SV_dict, outdir, range_filtering="theoretically_meaningful", expected_AF=1.0, replace=False, threads=4, median_insert_size=500, median_insert_size_sd=50, window_l=1000, mitochondrial_chromosome="mito_C_glabrata_CBS138", run_in_parallel=False):

    """Runs a benchmarking for several combinations of filters of a GridsssClove pipeline of a given bam file (sample_bam), writing files under outdir. The known SV are provided as a dictionary that maps each type of SV to a path where a table with the SVs are known .

    range_filtering indicates which type of simulation will be performed, it can be "large", "medium", "small", "single" and correlates with the range of parameters to use.
    expected_AF is the expected allele frequency, for a haploid it should be 1.0 

    median_insert_size is used to define small breakends to calculate their allele frequency

    window_l is used to define the coverage for winows of regions"""


    ###### DEFINE GENERAL THINGS

    start_time = time.time()

    # a bed with repeats
    repeats_bed = get_and_write_repeatsDFbed(reference_genome, threads=threads, replace=replace)

    # define a bed with the regions to ignore (in the future, this may include hetreozygous positions)
    regions_to_ignore_bed = repeats_bed 

    # define the median coverage of regions
    print("getting coverage")
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir, sample_bam, windows_file="none", replace=replace, window_l=window_l), sep="\t")
    median_coverage = np.median(coverage_df[coverage_df["#chrom"]!=mitochondrial_chromosome].mediancov_1); print("The median coverage is %i"%median_coverage)

    ##################
    
    # define the combinations of gridss parameters
    #maxcoverage_list = [int(median_coverage*5), 50000] # regions with coverage above this will be excluded. 50000 is the default
    #maxcoverage_list = [3000, 50000] # regions with coverage above this will be excluded. 50000 is the default. 3000 would be a very high coverage also, and it is the same for all for the sake of simplicity 
    maxcoverage_list = [50000]
    #ignore_regions_list  = [True, False] # whether to run blacklisted regions
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
            if ignore_regions is True: blacklisted_regions = regions_to_ignore_bed
            else: blacklisted_regions = ""

            # get the gridss outputs
            #print("running gridss")
            gridss_VCFoutput = run_gridss_and_annotateSimpleType(sample_bam, reference_genome, gridss_outdir, replace=replace, threads=threads, blacklisted_regions=blacklisted_regions, maxcoverage=maxcoverage)
                  
            bedpe_with_adds = get_bedpe_from_svVCF(gridss_VCFoutput, gridss_outdir, replace=replace)

            # get into dfs with generally interesting info
            #print("adding info")
            df_gridss = add_info_to_gridssDF(load_single_sample_VCF(gridss_VCFoutput), median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd) # this is a dataframe with some extra info
            #print("info added")
            df_bedpe = pd.read_csv(bedpe_with_adds, sep="\t")
            df_bedpe["IDs_set"] = df_bedpe.IDs.apply(lambda x: set(x.split("||")))

            # write the breakpoints. The point of this is that with many parameter combinations we may yield the same breakpoints, so that it's more efficient to create them first
            #print("getting combinations")
            outdir_parameter_combinations = "%s/several_parameter_combinations_filter_%s_af%.2f"%(gridss_outdir, range_filtering, expected_AF)
            #delete_folder(outdir_parameter_combinations) # DEBUG
            make_folder(outdir_parameter_combinations)
            filtersID_to_breakpoints = write_breakpoints_for_parameter_combinations_and_get_filterIDtoBpoints_gridss(df_gridss, df_bedpe, outdir_parameter_combinations, range_filtering=range_filtering, expected_AF=expected_AF, replace=replace) # this is a dataframe with all the filter combinations and the map between filterID and the actual filtering


            # define the paths to the breakpoints
            paths_to_bedpe_breakpoints = ["%s/%s/filtered_breakpoints.bedpe"%(outdir_parameter_combinations, filterID) for filterID in filtersID_to_breakpoints]

            # define inputs of the benchmarking pipeline
            inputs_benchmarking_pipeline = [(bedpe, know_SV_dict, reference_genome, sample_bam, median_coverage, replace, bedpe.split("/")[-2], True) for bedpe in paths_to_bedpe_breakpoints]

            if run_in_parallel is True:

                # initialize the list of benchmarking dfs
                all_benchmarking_dfs = []

                # go through each chunk of ncpus
                for Ichunk, chunk_inputs_benchmarking_pipeline in enumerate(chunks(inputs_benchmarking_pipeline, threads)):
                    print("working on chunk %i"%Ichunk)

                    # get the parallelized obtention of data
                    print("getting benchmarking for each set of filters in parallel")
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
            print("concatenating dfs")
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

    # add accuracy measures
    if "Fvalue" not in set(all_benchmarking_df.keys()):

        def calculate_Fvalue(r):

            pr = r["precision"]; rc = r["recall"]
            if pr<=0 or rc<=0: return 0.0
            else: return (2*pr*rc)/(pr+rc)

        all_benchmarking_df["Fvalue"] = all_benchmarking_df.apply(calculate_Fvalue, axis=1)

    # make plots of the benchmarking
    PlotsDir = "%s/plots_benchmark"%outdir
    makePlots_gridsss_benchmarking_oneGenome(all_benchmarking_df, PlotsDir)

    # get the time
    print("----It took %s seconds to run the whole benchmarking of one set of SV----"%(time.time() - start_time))

    return all_benchmarking_df

def get_parents_Ewa_experiments_withRun(sample, other_samples):
    
    """Takes a sample as (run, strain, replicate, condition) and other_samples, which is a set of (strain, replicate, condition(s)) it returns a set with the parents of sample. All are tuples"""
    
    # debug that they may be strings
    is_string=False
    if type(sample)==str:
        sample = tuple(sample.split("_"))
        other_samples = set([tuple(x.split("_")) for x in other_samples])
        is_string = True

    # add the normal parents of the same replicate
    run, strain, replicate, condition = sample
    condition_to_parentConditions = {'ANI': {"YPD", "WT"}, 'FLZ': {"YPD", "WT"}, 'AinF':  {"YPD", "ANI", "WT"}, 'ANIFLZ':  {"YPD", "WT"}, 'FinA':  {"YPD", "FLZ", "WT"}, "YPD":{"WT", "ANI", "FLZ", "AinF", "FinA"}, "WT":set()}

    parents = set([(run, s, r, c) for (run, s, r, c) in other_samples if s==strain and r==replicate and c in condition_to_parentConditions[condition]])

    # add the YPD parent, regardless of being from the same replicate
    if condition not in {"WT", "YPD"}: parents.update(set([(run, s, r, c) for (run, s, r, c) in other_samples if s==strain and c=="YPD"]))

    # add the WT
    if condition!="WT": parents.update(set([(run, s, r, c) for (run, s, r, c) in other_samples if s==strain and c=="WT"]))


    if is_string is True: parents = set(["_".join(x) for x in parents])
    
    return parents

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


def get_int_or_float_as_text(number):

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


def run_gridssClove_given_filters(sorted_bam, reference_genome, working_dir, median_coverage, replace=True, threads=4, gridss_blacklisted_regions="", gridss_VCFoutput="", gridss_maxcoverage=50000, median_insert_size=500, median_insert_size_sd=0, gridss_filters_dict=default_filtersDict_gridss, tol_bp=50, threshold_p_unbalTRA=0.7, run_in_parallel=True, max_rel_coverage_to_consider_del=0.1, min_rel_coverage_to_consider_dup=1.1, replace_FromGridssRun=False):

    """This function runs gridss and clove with provided filtering and parameters. This can be run at the end of a parameter optimisation process. It returns a dict mapping each SV to a table, and a df with the gridss """

    print("running gridss and clove with given parameter")
    make_folder(working_dir)

    # edit the replace, regarding if filtering from the run of GRIDSS
    if replace is True and replace_FromGridssRun is False: replace_FromGridssRun = True

    # first obtain the gridss output if it is not provided
    if file_is_empty(gridss_VCFoutput) or replace is True: gridss_VCFoutput = run_gridss_and_annotateSimpleType(sorted_bam, reference_genome, working_dir, replace=replace, threads=threads, blacklisted_regions=gridss_blacklisted_regions, maxcoverage=gridss_maxcoverage)

    ##### GET A LIST OF FILTERED BREAKPOINTS ########
    
    # get the output of gridss into a df
    print("getting gridss")
    df_gridss = add_info_to_gridssDF(load_single_sample_VCF(gridss_VCFoutput), median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd) # this is a dataframe with some info

    # filter according to gridss_filters_dict
    #print("filtering gridss")
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

    #################################################

    # run clove without checking filtering
    outfile_clove = "%s.clove.vcf"%(raw_bedpe_file)
    run_clove_filtered_bedpe(raw_bedpe_file, outfile_clove, sorted_bam, replace=replace_FromGridssRun, median_coverage=median_coverage, median_coverage_dev=1, check_coverage=False) #  REPLACE debug

    # add the filter of coverage to the clove output
    df_clove = get_clove_output_with_coverage_forTANDEL(outfile_clove, reference_genome, sorted_bam, replace=replace_FromGridssRun, run_in_parallel=run_in_parallel, delete_bams=run_in_parallel)
    maxDELcoverage = int(max_rel_coverage_to_consider_del*median_coverage)
    minDUPcoverage = int(min_rel_coverage_to_consider_dup*median_coverage) 
    df_clove["coverage_FILTER"] = df_clove.apply(lambda r: get_covfilter_cloveDF_row_according_to_SVTYPE(r, maxDELcoverage=maxDELcoverage, minDUPcoverage=minDUPcoverage), axis=1)

    # annotated clove 
    fileprefix = "%s.structural_variants"%outfile_clove
    remaining_df_clove, svtype_to_SVtable = write_clove_df_into_bedORbedpe_files_like_RSVSim(df_clove, fileprefix, reference_genome, sorted_bam, tol_bp=tol_bp, replace=replace_FromGridssRun, median_coverage=median_coverage, svtypes_to_consider={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}, threshold_p_unbalTRA=threshold_p_unbalTRA, run_in_parallel=run_in_parallel)

    # merge the coverage files in one
    merge_coverage_per_window_files_in_one(sorted_bam)

    return svtype_to_SVtable, df_gridss


def generate_tables_of_SV_between_genomes_gridssClove(query_genome, reference_genome, replace=False, threads=4, coverage=75, insert_size=250, read_lengths=[kb*1000 for kb in [0.5, 0.7, 0.9, 1, 1.3, 1.5, 2]], error_rate=0.0, expected_ploidy=1):

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

def get_dict_as_tuple(dictionary):

    """Takes a dict and converts it to a sorted tuple"""

    return tuple([(k, dictionary[k]) for k in sorted(dictionary.keys())])

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
    print("Warning: there is no single type of filtering that can fullfill all the requirements") 

    #lhfaeljkhfaekjhafekjhgkljhgrs
    #return df_best.iloc[0]

    # if you didn't find a single best, raise error
    print("\nthis is the best df:\n", df_best, "printing the non equal fields across all rows:\n")
    changing_fields = get_changing_fields_in_df_benchmark(df_best)
    for f in changing_fields:
        print("\t", f)
        for Irow in range(len(df_best)): print("\t\t", df_best[f].iloc[Irow])


    raise ValueError("There is not a single best filtering")

def get_df_accuracy_for_train_filer(r, outdir, test_gridss_info_dict, sorted_bam, reference_genome, median_coverage, replace, median_insert_size, median_insert_size_sd, test_SVdict):

    """define a function that takes a row of df_filters_train and returns a series with the accuracy values of each filter"""

    # define outdir
    working_dir = "%s/train_on_%s_%s_%s"%(outdir, r["genomeID"], r["ploidy"], r["svtype"]); make_folder(working_dir)
    #print("testing from %s"%working_dir)

    # define the file
    df_benchmark_filename = "%s/df_benchmark.tab"%working_dir

    if file_is_empty(df_benchmark_filename) or replace is True:

        # define the gridss_VCFoutput based on test_gridss_info_dict_under_outdir
        gridss_VCFoutput = test_gridss_info_dict[r["gridss_regionsToIgnoreBed"]][r["gridss_maxcoverage"]]

        # make a link under working_dir
        gridss_VCFoutput_underWorkDir = "%s/gridss_output.vcf"%(working_dir)
        print("testing...", gridss_VCFoutput_underWorkDir)
        if file_is_empty(gridss_VCFoutput_underWorkDir) or replace is True: run_cmd("ln -s %s %s"%(gridss_VCFoutput , gridss_VCFoutput_underWorkDir))

        # get the svs
        predicted_svtype_to_SVtable, df_gridss = run_gridssClove_given_filters(sorted_bam, reference_genome, working_dir, median_coverage, replace=replace, threads=multiproc.cpu_count(), gridss_blacklisted_regions=r["gridss_regionsToIgnoreBed"], gridss_VCFoutput=gridss_VCFoutput_underWorkDir, gridss_maxcoverage=r["gridss_maxcoverage"], median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, gridss_filters_dict=r["filters_dict"], tol_bp=50, threshold_p_unbalTRA=r["threshold_p_unbalTRA"], run_in_parallel=False, max_rel_coverage_to_consider_del=r["clove_max_rel_coverage_to_consider_del"], min_rel_coverage_to_consider_dup=r["clove_min_rel_coverage_to_consider_dup"], replace_FromGridssRun=False)

        # get the benchmarking df
        fileprefix = "%s.benchmarking"%working_dir
        df_benchmark_filtN = benchmark_processedSVs_against_knownSVs_inHouse(predicted_svtype_to_SVtable, test_SVdict, fileprefix, replace=replace, analysis_benchmarking=False, tolerance_bp=50, add_integrated_benchmarking=True)

        # add the metdadata
        df_benchmark_filtN["train_genomeID"] = r["genomeID"]
        df_benchmark_filtN["train_ploidy"] = r["ploidy"]
        df_benchmark_filtN["train_svtype"] = r["svtype"]

        # save
        df_benchmark_filtN.to_csv(df_benchmark_filename, sep="\t", header=True, index=False)

    else: df_benchmark_filtN = pd.read_csv(df_benchmark_filename, sep="\t")

    return df_benchmark_filtN

def get_benchmarking_df_for_testSVs_from_trainSV_filterSets(test_SVdict, outdir, df_filters_train, test_gridss_info_dict, genomeID, ploidy, sorted_bam, reference_genome, median_coverage, median_insert_size, median_insert_size_sd, replace):

    """This function takes a  set of test SVdict and it tries all the filters in df_filters_train on the gridss vcfs in gridss_info_dict, writing files under outdir. It returns a df with the accuracy for each svtype and filters from df_filters_train tested on SVdict  """

    start_time = time.time()

    # check that the df_filters_train contains unique vals for each genomeID, ploidy and svtype
    if len(df_filters_train)!=len(df_filters_train[["genomeID", "ploidy", "svtype"]].drop_duplicates()): raise ValueError('df_filters_train does not contain unique vals for "genomeID", "ploidy", "svtype"')

    # define the df_benchmark
    df_benchmark_all_filename = "%s/df_benchmark_all.tab"%outdir
    print("working on %s"%df_benchmark_all_filename)

    if file_is_empty(df_benchmark_all_filename) or replace is True:

        df_benchmark = pd.concat(list(df_filters_train.apply(lambda r: get_df_accuracy_for_train_filer(r, outdir, test_gridss_info_dict, sorted_bam, reference_genome, median_coverage, replace, median_insert_size, median_insert_size_sd, test_SVdict), axis=1)))

        # add metadata
        df_benchmark["test_genomeID"] = genomeID
        df_benchmark["test_ploidy"] = ploidy
        df_benchmark["test_svtype"] = df_benchmark.svtype

        # save
        print("saving %s"%df_benchmark_all_filename)
        df_benchmark.to_csv(df_benchmark_all_filename, sep="\t", header=True, index=False)

    else: df_benchmark = pd.read_csv(df_benchmark_all_filename, sep="\t")

    print("----It took %s seconds to run the whole benchmarking of one set of test filters----"%(time.time() - start_time))

    return df_benchmark


def get_and_report_filtering_accuracy_across_genomes_and_ploidies(df_benchmark, genomeID_to_knownSVdict, outdir, PlotsDir, reference_genome, replace=False, consider_integrated_filtering=True, threads=4, mitochondrial_chromosome="mito_C_glabrata_CBS138", run_in_parallel=False):

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
            print("getting all the variants integrated for each set of filters")
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
    print("Getting list of best filters for each genome, ploidy and svtype")
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

        return pd.Series({"%s_simName"%tag : simName, "%s_simType"%tag : simType})

    df_cross_benchmark[["train_simName", "train_simType"]] = df_cross_benchmark.train_genomeID.apply(lambda x: add_simulation_name_and_type(x, "train"))
    df_cross_benchmark[["test_simName", "test_simType"]] = df_cross_benchmark.test_genomeID.apply(lambda x: add_simulation_name_and_type(x, "test"))

    ######### get the plots ###########

    best_filters_dict = graph_fun.getPlots_filtering_accuracy_across_genomes_and_ploidies(df_cross_benchmark, PlotsDir)

    ###################################

    # get the best filters
    best_filters_series = df_best_filters.loc[("%s_simType_%s"%(best_filters_dict["simName"], best_filters_dict["simType"]), best_filters_dict["ploidy"], best_filters_dict["svtype"])]


    return df_cross_benchmark, best_filters_series

def run_GridssClove_optimising_parameters(sorted_bam, reference_genome, outdir, threads=4, replace=False, window_l=1000, n_simulated_genomes=1, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulation_types=["uniform", "biased_towards_repeats"], target_ploidies=["haploid", "diploid_homo", "diploid_hetero", "ref:2_var:1", "ref:3_var:1", "ref:4_var:1", "ref:5_var:1", "ref:9_var:1", "ref:19_var:1", "ref:99_var:1"], replace_covModelObtention=False, range_filtering_benchmark="theoretically_meaningful", coverage=20, known_genomes_withSV_and_shortReads_table=None, check_SVfromALNpipeline_simulated_genomes=True, expected_ploidy=1):

    """Takes some aligned reads and runs the GridssPipeline optimising the parameters of GRIDSS filtering. These are the different parameters of the function:

    - sorted_bam: the path to a sorted and indexed bam where we want to find the SV
    - reference_genome: the fasta of the reference genome
    - outdir: a directory where all the sample-specific files will be written. All the reference_genome-specific files will be written in the same dir as the reference_genome is
    - window_l: the length of the windows to generate the coverage model
    - n_simulated_genomes is the number of genomes that will be simulated as simulation replicates
    - mitochondrial_chromosome is the name of the mtDNA chromosome in the reference genome. This is passed to some parts of the pipeline to consider differently the gDNA and the mtDNA.
    - is_circular_*DNA indicates whethere the chromosomes are circular or not. This is important for the coverage simulation.
    - simulation_types indicates the types of SV performed
    - target_ploidies indicates which poplulations or ploidies have to be simulated. 2ref_1sv means that we will simulate a genome that has 2 reference genomes and 1 genome under structural variation (this can be a population sequencing)
    - replace_covModelObtention indicates whether the process of predicting coverage from seq features has to be replaced
    - range_filtering_benchmark indicates which type of simulation will be performed, it can be "large", "medium", "small", "single", "theoretically_meaningful". This is passed to benchmark_GridssClove_for_knownSV.
    - coverage is the filter of coverage applied when calculating heterozygous positions in the reference_genome
    - known_genomes_withSV_and_shortReads_table is a file with a table that has three fields: assembly,shoort_reads1,short_reads2 . This can be, for example a set of NANOPORE assemblies of Candida glabrata and the corresponding short reads' sequencing in YPD
    - check_SVfromALNpipeline_simulated_genomes indicates whether to becnhmark how the generate_tables_of_SV_between_genomes function works on simulated genomes
    - expected_ploidy is a number that states the expected ploidy. If you are not sure it should be very high.

    The workflow of this pipeline is the following. 

    - get a model that predicts relative coverage of windows of the genome (of length window_l) for this sample
    - generate n_simulated_genomes (as replicates) (and for simulation_types) with structural variation. Then, for each of them:
    
        - Generate ensembles of genomes (ploidies) according to target_ploidies. For each ensemble:
                    
            - simulate reads (wgsim) in a way that resembles the coverage model. In the case of the mtDNA, simulate a number of copies of mtDNA that is proportional to the median of the observed mtDNA relative coverage. For circular chromosomes, simulate also reads for a window that will cover the end and start of the chromosome, to simulate that it is circular
            - align them with bwa mem
            - run gridss on aligned reads
            - for each combination of filters (in parallel):

                - run clove for the breakpoints that pass the filter (both breakends considered)
                - calculate precision, recall and F score. Keeping the best for further analysis
    
    # for simulating reads, we reccommend to use a 

    """

    # initialize the start
    pipeline_start_time = time.time()

    # define plots dir
    make_folder(outdir)
    PlotsDir = "%s/plots"%outdir; make_folder(PlotsDir)

    # put the reference genome in the outdir to avoid cross problems with the other runs
    reference_genome_dir = "%s/reference_genome_dir"%(outdir)
    make_folder(reference_genome_dir)
    new_reference_genome_file = "%s/reference_genome.fasta"%reference_genome_dir
    run_cmd("cp %s %s"%(reference_genome, new_reference_genome_file))
    reference_genome = new_reference_genome_file

    # clean the reference genome windows files
    clean_reference_genome_windows_files(reference_genome)

    ########### MODELING COVERAGE ############

    # get info that models the coverage of the sample
    objects_coverage_modelling_file = "%s/objects_coverage_modelling_windows_%ikb.py"%(outdir, window_l)


    if file_is_empty(objects_coverage_modelling_file) or replace is True:

        print("Getting model for the coverage for this sample...")
        df_cov, chromosome_to_relPosVScovCoefficients, predictors_predictResidualFromPos, lmObject_predictResidualFromPos, scrambled_chromosomes = get_sampleSpecific_genomeInfoDF_and_fitInfo(reference_genome, outdir, sorted_bam, mitochondrial_chromosome=mitochondrial_chromosome, replace=replace, window_l=window_l, threads=threads, replace_covModelObtention=replace_covModelObtention)

        # make the plots of the coverage prediction
        plot_coverage_correlations(df_cov, PlotsDir, plots={"goodness_coverage_predictions"})

        # get the predictors that are related to mapping, not real reads, and define the interesting
        predictors_related_to_mapping = {"fraction_covered_by_repeat", "relative_n_repeats", "meanUniquelyMappable"}
        predictors_ResCovAfterPosCorrection = list(set(predictors_predictResidualFromPos).difference(predictors_related_to_mapping)) 

        # build a model that predicts the residual of correcting by position from the predictors_ResCovAfterPosCorrection alone, using df_cov.residual_cov_vs_predicted_from_loc as Y
        print("Getting model that predicts the residuals of the coverage from seq features. Training on normal chromosomes and coverage btw 0.05 and 4")
        df_filt = df_cov[~(df_cov.chromosome.isin(scrambled_chromosomes.union({mitochondrial_chromosome}))) & (df_cov.rel_coverage_median<=4) & (df_cov.rel_coverage_median>0.05)]
        lm_predictResidualFromPos = linear_model.LinearRegression(n_jobs=threads).fit(df_filt[predictors_ResCovAfterPosCorrection], df_filt.residual_cov_vs_predicted_from_loc)

        # define and save the important objects that define the model
        tuple_important_objects_prediction = (chromosome_to_relPosVScovCoefficients, predictors_ResCovAfterPosCorrection, lm_predictResidualFromPos)

        # save
        save_object(tuple_important_objects_prediction, objects_coverage_modelling_file)

    else: chromosome_to_relPosVScovCoefficients, predictors_ResCovAfterPosCorrection, lm_predictResidualFromPos = load_object(objects_coverage_modelling_file)
    print("COVERAGE MODEL OBTAINED...")

    # this info will be used to simulate reads per window from position and sequence features for each simulated dataset

    ##########################################

    ############ GENERAL OPERATIONS THAT WILL BE NEEDED FOR ALL THE STEPS #####

    # the dir and genome names
    genome_dir = "/".join(reference_genome.split("/")[0:-1])
    genome_name = reference_genome.split("/")[-1].split(".")[0]

    # map each chromosome to length
    chr_to_len = {seq.id: len(seq.seq) for seq in SeqIO.parse(reference_genome, "fasta")}

    # define a function that sets the correct coordinates
    def get_zeros_for_negatives(coord):

        if coord<0: return 0
        else: return coord

    def get_plus_1kb(r):

        coord = r["end_0based"] + 1000
        len_chr = chr_to_len[r["chromosome"]]
        if coord > len_chr: return len_chr
        else: return coord

    # get a bed file that contains regions with repeats
    print("getting repeats")
    repeats_df = get_repeat_maskerDF(reference_genome, threads=threads, replace=replace)
    repeats_df["start_0based-1kb"] = (repeats_df.begin_repeat - 1 - 1000).apply(get_zeros_for_negatives)
    repeats_df["end_0based"] = repeats_df.end_repeat - 1
    repeats_df["end_0based+1kb"] = repeats_df.apply(get_plus_1kb, axis=1)
    repeats_regions_bed = "%s/%s_repeats_regions_+-1kb.bed"%(genome_dir, genome_name)
    repeats_df[["chromosome", "start_0based-1kb", "end_0based+1kb", "strand"]].to_csv(repeats_regions_bed, sep="\t", header=False, index=False)

    # get the info of the reference genome with predictions of coverage per window
    df_REFgenome_info =  get_windows_infoDF_with_predictedFromFeatures_coverage(reference_genome, chromosome_to_relPosVScovCoefficients, predictors_ResCovAfterPosCorrection, lm_predictResidualFromPos, replace=replace, window_l=window_l, threads=threads, plots_prefix="%s/reference_genome"%PlotsDir)

    # count the length od the reads
    read_length = get_read_length(sorted_bam, threads=threads, replace=replace)
    print("The median read length is %i"%read_length)

    # count total number of reads
    total_nread_pairs = count_number_read_pairs(sorted_bam, replace=False, threads=threads)
    #total_nread_pairs  = 1000000 # this is to debug the simulation pipeline
    expected_coverage_per_bp = int((total_nread_pairs*read_length) / sum(chr_to_len.values())) +  1 # the expected coverage per position with pseudocount
    print("There are %i read pairs in your library. The expected coverage is %ix."%(total_nread_pairs, expected_coverage_per_bp))

    # calculate the insert size
    median_insert_size, median_insert_size_sd  = get_insert_size_distribution(sorted_bam, replace=replace, threads=threads)
    print("The median insert size is %i, with an absolute deviation of %i"%(median_insert_size, median_insert_size_sd))

    # simulate reads for the reference if you are not only simulating haploid
    if set(target_ploidies)!={"haploid"}: 
        outdir_ref = "%s/simulation_reference_genome_%ibp_windows"%(outdir, window_l)
        simulated_reference_bam_file = simulate_and_align_PairedReads_perWindow(df_REFgenome_info, reference_genome, reference_genome, total_nread_pairs, read_length, outdir_ref, median_insert_size, median_insert_size_sd, replace=replace, threads=threads, mitochondrial_chromosome=mitochondrial_chromosome)

    ###########################################################################

    # initialize a df with all the benchmarking data
    df_benchmark_all = pd.DataFrame()
    genomeID_to_knownSVdict = {}

    df_benchmark_all_file = "%s/df_benchmark_all.py"%outdir
    genomeID_to_knownSVdict_file= "%s/genomeID_to_knownSVdict.py"%outdir

    if file_is_empty(df_benchmark_all_file) or file_is_empty(genomeID_to_knownSVdict_file) or replace is True:

        ################ SIMULATION PIPELINE ################ 

        # go throigh each simulation (these are technical replicates of the pipeline)
        for simulation_ID in range(1, n_simulated_genomes+1):
            print("working on simulation %i"%simulation_ID)

            # get an outdir where all the simulations of this ID will be stored
            simulation_outdir = "%s/simulation_%i"%(outdir, simulation_ID); make_folder(simulation_outdir)

            # get the simulated rearranged genomes related to simulation_types. This is a dictionary with the simulations
            simType_to_rearrangedGenome = rearrange_genomes_simulateSV(reference_genome, simulation_outdir, repeats_1kbwindow_bed=repeats_regions_bed, replace=replace, simulation_types=simulation_types, mitochondrial_chromosome=mitochondrial_chromosome) # !!!! i have to change the representation of translocations with different 3' and 5' end

            # go through each type of simulation (either randomly distributed or focused on repeats)
            for simType, rearranged_genome in simType_to_rearrangedGenome.items(): 
                print("working on %s"%simType)
                genomeID = "simulation_%i_simType_%s"%(simulation_ID, simType)

                # define the directory of the simulation
                simType_dir = "/".join(rearranged_genome.split("/")[0:-1]); make_folder(simType_dir)

                # get a df that has genome info and predicted coverage from location and features (trained previously) without repeat info
                df_genome_info =  get_windows_infoDF_with_predictedFromFeatures_coverage(rearranged_genome, chromosome_to_relPosVScovCoefficients, predictors_ResCovAfterPosCorrection, lm_predictResidualFromPos, replace=replace, window_l=window_l, threads=threads, plots_prefix="%s/simulation"%simType_dir)

                # get the aligned reads to the reference
                simulation_bam_file = simulate_and_align_PairedReads_perWindow(df_genome_info, rearranged_genome, reference_genome, total_nread_pairs, read_length, simType_dir, median_insert_size, median_insert_size_sd, replace=replace, threads=threads, mitochondrial_chromosome=mitochondrial_chromosome)

                # define a path to the known SVs
                know_SV_dict = {var : "%s/%s.tab"%(simType_dir, var) for var in ["tandemDuplications", "translocations", "insertions", "inversions", "deletions"]}
                know_SV_dict["translocations"] = "%s/translocations.tab.corrected"%(simType_dir)
                know_SV_dict["insertions"] = "%s/insertions.tab.corrected"%(simType_dir)

                # add the "remaining" cathegory, as an empty field
                remaining_file = "%s/remaining_sv.tab"%simType_dir
                open(remaining_file, "w").write("ID\t#CHROM\tPOS\tCHR2\tSTART\tEND\tSVTYPE\n"+"iii\tccc\t0\tyyy\t0\t0\tzzz\n")
                know_SV_dict["remaining"] = remaining_file

                # map each genomeID to the known variants
                genomeID_to_knownSVdict[genomeID] = know_SV_dict

                # go through each of the target ploidies and generate the resulting bam files:
                for ploidy in target_ploidies:
                    print("working on %s"%ploidy)

                    # define the final sorted bam depending on the ploidy (which also includes populations)
                    ploidy_merged_bam = get_merged_bamfile_for_ploidy(variant_bamfile=simulation_bam_file, reference_bamfile=simulated_reference_bam_file, ploidy=ploidy, replace=replace, threads=threads)

                    # calculate the expected fraction of reads comming from each genome
                    fraction_var, fraction_ref = get_fractions_reads_for_ploidy(ploidy)

                    # write a table and some files with the benchmarking of several filtering strategies of the data
                    ploidy_dir = "%s/benchmark_GridssClove_%s"%(simType_dir, ploidy); make_folder(ploidy_dir)

                    # get a df with a benchmark of many different parameters. This will also report some plots with the 
                    benchmarking_df = benchmark_GridssClove_for_knownSV(ploidy_merged_bam, reference_genome, know_SV_dict, ploidy_dir, range_filtering=range_filtering_benchmark, expected_AF=fraction_var, replace=replace, threads=threads, median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, window_l=window_l, mitochondrial_chromosome=mitochondrial_chromosome)

                    # add some parms and keep
                    benchmarking_df["genomeID"] = [genomeID]*len(benchmarking_df)
                    benchmarking_df["ploidy"] = [ploidy]*len(benchmarking_df)
                    benchmarking_df["sorted_bam"] = [ploidy_merged_bam]*len(benchmarking_df)
                    benchmarking_df["median_insert_size"] = [median_insert_size]*len(benchmarking_df)
                    benchmarking_df["median_insert_size_sd"] = [median_insert_size_sd]*len(benchmarking_df)
                    df_benchmark_all = df_benchmark_all.append(benchmarking_df, sort=True)

        #####################################################
        print("GRIDSS simulation finished correctly")

        # save important files
        print("saving important files...")
        save_object(df_benchmark_all, df_benchmark_all_file)
        save_object(genomeID_to_knownSVdict, genomeID_to_knownSVdict_file)

        ########## 

    else:
        print("GRIDSS simulation finished correctly. Loading previous files ...")
        df_benchmark_all = load_object(df_benchmark_all_file)
        genomeID_to_knownSVdict = load_object(genomeID_to_knownSVdict_file)

    ################### REPORT ACCURACY BETWEEN PARAMETERS OF DIFFERENT OPTIMISATIONS ####################

    # we will take the parameters that work best for all the simulations that we input

    # define the outputs
    outdir_benchmarking = "%s/benchmarking_all_filters_for_all_genomes_and_ploidies"%outdir; make_folder(outdir_benchmarking)
    PlotsDir_benchmarking = "%s/plots"%outdir_benchmarking; make_folder(PlotsDir_benchmarking)


    print("getting report of the accuracies between simulations")

    df_best_filters, best_f = get_and_report_filtering_accuracy_across_genomes_and_ploidies(df_benchmark_all, genomeID_to_knownSVdict, outdir_benchmarking, PlotsDir_benchmarking, reference_genome, replace=replace, consider_integrated_filtering=True, threads=threads)

    # save them as an object
    series_best_filters_file = "%s/series_best_filters_across_simulations.py"%outdir
    save_object(best_f, series_best_filters_file)

    ######################################################################################################


    # run gridss with the optimised filters
    outdir_gridss_final = "%s/final_gridss_running_with_optimum_parameters"%outdir; make_folder(outdir_gridss_final)
    final_gridss_vcf = "%s/output_gridss.vcf"%outdir_gridss_final

    # define the median coverage across 1kb windows of the genome
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_gridss_final, sorted_bam, windows_file="none", replace=replace, window_l=window_l), sep="\t")
    median_coverage = np.median(coverage_df[coverage_df["#chrom"]!=mitochondrial_chromosome].mediancov_1)
    print("The median coverage is %i"%median_coverage)

    # get coverage per different window sizes
    for size in [50, 100, 500, 5000, 10000, 20000]: coverage_file_window_s = generate_coverage_per_window_file_parallel(reference_genome, outdir_gridss_final, sorted_bam, windows_file="none", replace=replace, window_l=size)


    print("running final gridss with optimised parameters...")
    final_svtype_to_SVtable, df_gridss = run_gridssClove_given_filters(sorted_bam, reference_genome, outdir_gridss_final, median_coverage, replace=replace, threads=multiproc.cpu_count(), gridss_blacklisted_regions=best_f["gridss_regionsToIgnoreBed"], gridss_VCFoutput=final_gridss_vcf, gridss_maxcoverage=best_f["gridss_maxcoverage"], median_insert_size=median_insert_size, median_insert_size_sd=median_insert_size_sd, gridss_filters_dict=best_f["filters_dict"], tol_bp=50, threshold_p_unbalTRA=best_f["threshold_p_unbalTRA"], run_in_parallel=True, max_rel_coverage_to_consider_del=best_f["clove_max_rel_coverage_to_consider_del"], min_rel_coverage_to_consider_dup=best_f["clove_min_rel_coverage_to_consider_dup"], replace_FromGridssRun=False)

    # at the end, remove all the mosdepth and windows files under the reference
    clean_reference_genome_windows_files(reference_genome)
    
    print("GRIDSS pipeline finished correctly")

    print("--- the gridss pipeline optimising parameters took %s seconds in %i cores ---"%(time.time() - pipeline_start_time, threads))

    # at the end remove all the files that are unnecessary

    # simulations
    for simulation_ID in range(1, n_simulated_genomes+1): delete_folder("%s/simulation_%i"%(outdir, simulation_ID))

    # remove the cross-benchmarking files
    delete_folder("%s/cross_benchmarking_files"%outdir_benchmarking)

    # generate a file that indicates whether the gridss run is finished
    final_file = "%s/gridss_finished_file_final_gridss_final_with_window_size.txt"%outdir
    open(final_file, "w").write("gridss finished...")





