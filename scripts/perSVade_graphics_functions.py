#!/usr/bin/env python

# These are functions to perform plots for the VarCall pipeline

##### DEFINE ENVIRONMENT #######

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
import seaborn as sns
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
import webcolors
from colour import Color
import random as rd
import plotly.plotly as py
import plotly.figure_factory as ff
import plotly.offline as off_py
import plotly.graph_objs as go
from plotly import tools
import cufflinks as cf
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

warnings.simplefilter(action='ignore', category=pd.core.common.SettingWithCopyWarning)
#pd.options.mode.chained_assignment = 'raise'

# load a specific matplotlib library for cluster envs
try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

except: import matplotlib.pyplot as plt

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# import functions
import perSVade_functions as fun

#########################################


#################################################################
##################### FUNCTIONS TO PROCESS SVs ##################
#################################################################

def ask_if_overlapping_breakends_in_parents_withEqualChromosomes(r, parents_gridss_df, tol_bp):

    """Returns a boolean that indicates whether ther is any breakend in parents that overlaps with the breakedn in r by less than tol_bp bp, where the orientation of the brekends is expected to be the same"""

    # get the df in parents_gridss_df that fullfils the overlap conditions
    df_parents_overlap = parents_gridss_df[(parents_gridss_df.other_orientation==r["other_orientation"]) 
                                         & ((parents_gridss_df.POS-r["POS"]).apply(abs)<=tol_bp) 
                                         & ((parents_gridss_df.other_position-r["other_position"]).apply(abs)<=tol_bp) ]

    # return True if there is some overlap
    if len(df_parents_overlap)>0: return True
    else: return False

def get_eventIDs_already_in_parents(sampleID, parentIDs, sampleID_to_dfGRIDSS, tol_bp=200):

    """This function returns a set with the eventIDs that are in sample but also in any of the parents. The idea is that an event is considered to be in a parent if any of the breakends overlaps with a breakend in the parent, where the orientation is the same and the position is less or equal far appart than tol_bp """
    print(sampleID)

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

def get_sampleID_to_svtype_to_svDF_filtered(sampleID_to_svtype_to_file, sampleID_to_dfGRIDSS, sampleID_to_parentIDs={}, breakend_info_to_keep=['#CHROM', 'POS', 'other_coordinates', 'allele_frequency', 'allele_frequency_SmallEvent', 'real_AF', 'FILTER', 'inserted_sequence', 'has_poly16GC', 'length_inexactHomology', 'length_microHomology']):

    """This function takes a dictionary that maps sampleIDs to svtpes and the corresponding files (the ones returned in the SV-calling pipeline) and the corresponding gridss dataframes, returning a sampleID_to_svtype_to_svDF, after removal of the variants that are in sampleID_to_parentIDs. The idea is that if any breakend is shared between the sample and any of the parents it is removed and so does any sv that relies on this breakend. Other arguments:

    - working_dir is the path were eventually files will be created

    All across this function, eventID reflects the INFO_EVENT with an extra o

    """

    # debug empty 
    if len(sampleID_to_parentIDs)==0: sampleID_to_parentIDs = {s:set() for s in sampleID_to_svtype_to_file}

    # map each sampleID to the eventIDs (with a final o, as formated in the sampleID_to_svtype_to_file) that are already in the parents
    all_samples = sorted(sampleID_to_parentIDs.keys())
    get_eventIDs_already_in_parents_inputs = [(s, sampleID_to_parentIDs[s], sampleID_to_dfGRIDSS) for s in all_samples]

    # get the overlapping events with parallelization
    
    """
    print("getting parent events in parallel")
    with multiproc.Pool(multiproc.cpu_count()) as pool:

        map_eventIDs_already_in_parents = pool.starmap(get_eventIDs_already_in_parents, get_eventIDs_already_in_parents_inputs)
        pool.close()
        pool.terminate()

    """

    # get the overlapping events with a map
    map_eventIDs_already_in_parents = list(map(lambda x: get_eventIDs_already_in_parents(x[0], x[1], x[2]), get_eventIDs_already_in_parents_inputs))
    

    sampleID_to_eventIDs_alreadyInParents = dict(zip(all_samples, map_eventIDs_already_in_parents))

    #sampleID_to_eventIDs_alreadyInParents = {sampleID:  get_eventIDs_already_in_parents(sampleID, parentIDs, sampleID_to_dfGRIDSS) for sampleID, parentIDs in sampleID_to_parentIDs.items()}

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

    return sampleID_to_svtype_to_svDF

def make_flat_listOflists(LoL):

    return list(itertools.chain.from_iterable(LoL))

def chunks(l, n):
    
    """Yield successive n-sized chunks from a list l"""
    
    for i in range(0, len(l), n):
        yield l[i:i + n]

def range_any_direction(x, y, step=1):

    """Takes an x and a y value and generates a range between them, regardless of the orientation"""

    # define range
    if y>=x: range_obj = list(range(x, y+1, step))
    elif y<x: 

        # get inverted range
        range_obj = list(range(y, x+1, step))
        range_obj.reverse()

    return range_obj

def map_protIDs_to_aaGenomicLocations(protIDs, gffdf, protID_in_gff="upmost_parent", expected_protID_to_protLen={}):

    """Given a set of proteinIDs and a gff dataframne. Return a dictionary that maps each protein ID to the aa location and the genomic coordinates that code this location. protID_in_gff defines the field in gffdf that corresponds to the protID"""

    # get a df only with the CDS
    gff_cds = gffdf[gffdf.feature=="CDS"]

    # add the coordinates that define the start of the protein
    def get_real_coords(row_gff):

        strand = row_gff["strand"]
        start = row_gff["start"]
        end = row_gff["end"]

        if strand=="+":
            prot_start = start
            prot_end = end

        elif strand=="-":
            prot_start = end
            prot_end = start

        return pd.Series([prot_start, prot_end])   

    gff_cds[["prot_start", "prot_end"]] = gff_cds.apply(get_real_coords, axis=1)

    # get only the cds that are from protIDs
    gff_cds = gff_cds[gff_cds[protID_in_gff].isin(protIDs)]

    # get a series that only includes the proteinIDs
    protIDs_series = pd.Series(list(protIDs))

    # define a function that takes a proteinID and returns a dictionary between all it's protein locations and the genomic coordinates
    def get_genomic_aaLocs(protID):

        # get interesting parameters
        gff = gff_cds[gff_cds[protID_in_gff]==protID]
        strand = gff.iloc[0].strand

        # sort depending on the strand
        if strand=="+": gff = gff.sort_values(by="prot_start")
        elif strand=="-": gff = gff.sort_values(by="prot_start", ascending=False)
        else: raise ValueError("strand %s is invalid"%strand)

        # initialize parms:
        aa_positions = []
        previous_cds_incomplete_aa_positions = []

        # go throug the CDS in a way that is in the orientation of the protein, regardless of the strain
        for prot_start, prot_end in gff[["prot_start", "prot_end"]].values:

            # when the previous cds did not yield any aberrant aa position
            if previous_cds_incomplete_aa_positions==[]: 

                # get the aa positions of this cds
                cds_aa_positions = list(chunks(range_any_direction(prot_start, prot_end), 3))

            # if the previous cds yielded an incomplete aa, add it
            else:

                # get the extra positions
                if prot_start < prot_end: shift_in_start_prevaa_calc = 3 - len(previous_cds_incomplete_aa_positions) - 1
                elif prot_start > prot_end: shift_in_start_prevaa_calc = - (3 - len(previous_cds_incomplete_aa_positions) - 1)
                else: raise ValueError("the start and end of a cds can't be the same")

                # get the complete aa with previous cds
                complete_aa_with_previous_cds = previous_cds_incomplete_aa_positions + range_any_direction(prot_start, prot_start+shift_in_start_prevaa_calc)

                # get the shift in the start to get the correct aa calculation
                if prot_start < prot_end: shift_in_start_newaa = 3 - len(previous_cds_incomplete_aa_positions)
                elif prot_start > prot_end: shift_in_start_newaa = - (3 - len(previous_cds_incomplete_aa_positions))

                # get the other aa positions
                cds_aa_positions = [complete_aa_with_previous_cds] + list(chunks(range_any_direction(prot_start+shift_in_start_newaa, prot_end), 3))
                previous_cds_incomplete_aa_positions = []

            # if there is something incomplete in the end, save, and also change cds_aa_positions
            if len(cds_aa_positions[-1])<3: 
                previous_cds_incomplete_aa_positions = cds_aa_positions[-1]
                cds_aa_positions = cds_aa_positions[0:-1]

            # keep
            aa_positions += cds_aa_positions

        # delete the last codon, which is the stop
        aa_positions = aa_positions[0:-1]

        # chek that the lengths are correct:
        """
        if protID in expected_protID_to_protLen:
            if len(aa_positions)!=expected_protID_to_protLen[protID]: 
                print("warn. Expected len: %i, calculated len: %i"%(expected_protID_to_protLen[protID], len(aa_positions)))
                print((gff.end-gff.start)/3, "\n\n")
        """

        # keep in the dictionary
        return dict(zip(list(range(1, len(aa_positions)+1)), aa_positions))

    # return a mapping between protID and the corresponding genomic locations
    return dict(zip(protIDs_series, protIDs_series.apply(get_genomic_aaLocs)))

def merge_InterProDF_and_gff3_df(IPdf, gffdf, protID_in_gff="upmost_parent"):

    """Takes a df from interproscan annotations and one of gff3, it  appends the InterPro fields to the GFF, so that annotations are handled as features. The protID_in_gff indicates the field in the gff that corresponds to the unqiue identifier of protenis, so that all values in "proteinID" in IPdf are in protID_in_gff """

    # get gff that have CDS
    gff_cds = gffdf[gffdf.feature=="CDS"]

    # get all proteinIDs
    all_protIDs = set(IPdf.proteinID)

    # map them to genomic locations
    print("Getting genomic locations")
    protID_to_protLocation_to_GenomicPositions = map_protIDs_to_aaGenomicLocations(all_protIDs, gffdf, protID_in_gff=protID_in_gff, expected_protID_to_protLen=dict(zip(IPdf.proteinID, IPdf.seq_length)))

    # get all the important features
    def get_gff_recordDicts(row, gff):

        """Takes a row of an IP annotation df and returns the interesting information, which are all the gff fields. 
        This function returns a list of dictionaries, each of which has all the info for a gff record. Note that if an IP annotation is split across several cds it has to be recorded as two sepparate records.
        """

        # define accession
        accession = row["signature_accession"]

        # define description
        description = row["signature_description"]
        if description=="no description": description = row["InterPro_annotation"]
        if pd.isna(description): description = ""

        # define feature
        annotation = "%s/%s/%s"%(row["type_analysis"], accession, description)
        feature = "%s/%s"%(accession, description)

        # define things about the feature
        protID = row["proteinID"]
        aa_start = row["start_location"]
        aa_end = row["end_location"]

        # define the genomic coordinates. These will be in a way that start<end
        aa_loc_dict = protID_to_protLocation_to_GenomicPositions[protID]
        
        # if the last aa is not in the dictionary (truncated protein) append
        if aa_end not in aa_loc_dict: aa_end = max(aa_loc_dict.keys())

        # get coordinates
        nt_start = min([min(aa_loc_dict[aa]) for aa in {aa_start, aa_end}])
        nt_end = max([max(aa_loc_dict[aa]) for aa in {aa_start, aa_end}])

        # get the cds information
        cds_info = gff[(gff[protID_in_gff]==protID)].iloc[0]

        gff_info_dict = {"chromosome":cds_info["chromosome"],
                         "source":"InterPro",
                         "feature":feature,
                         "start":nt_start,
                         "end":nt_end,
                         "blank1":cds_info["blank1"],
                         "strand":cds_info["strand"],
                         "blank2":cds_info["blank2"],
                         "annotation":annotation,
                         "orf_classification":"no",
                         "Gene":cds_info["Gene"],
                         "Parent":protID,
                         "parent_feature_type":"ORF",
                         "Name":protID,
                         "Note":"",
                         "ID":"%s-%s-%i"%(protID, accession, aa_start),
                         "Alias":"",
                         "duplicated_ID":False,
                         "upmost_parent":cds_info["upmost_parent"]}

        return [gff_info_dict]

        # consider splitting
        """
        # calculate cds that encdoe this feature. Considering only one
        cds_p = gff[(gff[protID_in_gff]==protID)]
        cds = cds_p[(cds_p.start<=nt_start) & (cds_p.end>=nt_end)] # coordinates overlap

        # print the case where the cds is >2
        if len(cds)==0:
            print("\n\n", nt_start, nt_end, cds_p[["start", "end"]], row)

        """

    # get a list of the recods list
    print("Getting gff-like structure")
    gff_record_info_list =  make_flat_listOflists(list(IPdf.apply(lambda r: get_gff_recordDicts(r, gff_cds), axis=1)))
    max_idxGFF = max(gffdf.numeric_idx)
    gffNumericIDX_to_recordDict = dict(zip(range(max_idxGFF+1, max_idxGFF+1+len(gff_record_info_list)+1) , gff_record_info_list))
    IP_gff_df =  pd.DataFrame(gffNumericIDX_to_recordDict).transpose()
    IP_gff_df["numeric_idx"] = IP_gff_df.index
    IP_gff_df = IP_gff_df[list(gffdf.keys())].set_index("upmost_parent", drop=False)

    # return the merged dataframes
    return gffdf.append(IP_gff_df)

def load_InterProAnnotation(path_to_InterPro_annotation):

    """This function retuns a dataframe out of the interpro annotation"""
    col_names = ["proteinID", "MD5digest", "seq_length", "type_analysis", "signature_accession", "signature_description", "start_location", "end_location", "score", "status", "date", "InterPro_accession", "InterPro_annotation", "GO_InterPro"]
    return pd.read_table(path_to_InterPro_annotation, names=col_names)

def load_gff3_intoDF(gff_path):

    """ Takes the path to a gff and loads it into a df"""

    # load
    gff = pd.read_csv(gff_path, sep="\t", skiprows=len([l for l in open(gff_path, "r") if l.startswith("#")]), names=["chromosome", "source", "feature", "start", "end", "blank1", "strand", "blank2", "annotation"])

    # set all possible annotations
    all_annotations = set.union(*[set([x.split("=")[0] for x in an.split(";")]) for an in gff.annotation])

    def list_to_str(list_x):

        if list_x==[]: return ""
        else: return list_x[0]
    
    for anno in all_annotations:

        gff[anno] = gff.annotation.apply(lambda x: list_to_str([a.split(anno)[1].lstrip("=") for a in x.split(";") if anno in a]))

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
    def get_utmost_parent(row, gff_df):

        """Takes a row and the gff_df and returns the highest parent. The index has to be the ID of the GFF"""

        # when you have found the parent it has no parent, so the ID is the upmost_parent
        if row["Parent"]=="": return row["ID"]

        # else you go to the parent
        else: return get_utmost_parent(gff_df.loc[row["Parent"]], gff_df)

    print("Getting top parent")
    gff["upmost_parent"] = gff.apply(lambda row: get_utmost_parent(row, gff), axis=1)
    gff = gff.set_index("upmost_parent", drop=False)

    return gff

def get_IDs_str_SV_deppending_on_overlapping(r, df_defined_svs, svtype, equal_fields, approximate_fields, tol_bp=50):

    """Takes a row of a sv df and a df of putatively overlapping sv df (df_defined_svs) and returns a set of IDs, which contains the  svtype, equalFields and approximate_fields """

    # define the ID string
    IDs_str = {"%s_%s"%(svtype, "-".join(["%s:%s"%(f, r[f]) for f in equal_fields+approximate_fields]))}
    
    # Find overlapping dfs
    if len(df_defined_svs)>0:

        # get the df that is overlapping
        overlapping_df = df_defined_svs[df_defined_svs.apply(lambda rd: all([rd[f]==r[f] for f in equal_fields]) and all([abs(rd[f]-r[f])<=tol_bp for f in approximate_fields]), axis=1)]

        if len(overlapping_df)>0: IDs_str = set.union(*overlapping_df.IDs_str)

    return IDs_str

def get_df_overlapping_SV(sampleID_to_SV_dict, tol_bp=50, expected_svs={'inversions', 'remaining', 'tandemDuplications', 'translocations', 'insertions', 'deletions'}):

    """This function takes a dictionary that maps each sample of interest to a dictionary that maps each svtype to a table or df with this info. It returns a df where each row is a sample and each col is a combination of svtype_coords and contains whether it is found in each sample (1 or 0). The idea is that if 2 vars ovarlap in <=50 they are considered the same."""

    # first convert data structure to df if not already done
    def get_df(x):
        if type(x)==str: return pd.read_csv(x, sep="\t")
        else: return x

    sampleID_to_svtype_to_df = {sID : {sv : get_df(filename_or_df) for sv, filename_or_df in SV_dict.items()} for sID, SV_dict in sampleID_to_SV_dict.items()}

    ####### map each sample to a set of unique vars, named consistently if they overlap ######

    # map each svtype to a set of equal and flexible (coords) fields to define overlaps
    svtype_to_fieldsDict = {"inversions": {"equal_fields": ["Chr"], "approximate_fields": ["Start", "End"]}, 
                            "tandemDuplications": {"equal_fields": ["Chr"], "approximate_fields": ["Start", "End"]}, 
                            "deletions": {"equal_fields": ["Chr"], "approximate_fields": ["Start", "End"]}, 
                            "translocations": {"equal_fields": ["ChrA", "ChrB", "Balanced"], "approximate_fields": ["StartA", "EndA", "StartB", "EndB"]}, 
                            "insertions": {"equal_fields": ["ChrA", "ChrB", "Copied"], "approximate_fields": ["StartA", "EndA", "StartB", "EndB"]}, 
                            "remaining": {"equal_fields": ["#CHROM", "CHR2", "SVTYPE"], "approximate_fields": ["POS", "START", "END"]}}

    # initialize a df of vars that are already found, these will be 
    svtype_to_already_found_dfs = {svtype : pd.DataFrame() for svtype in set.union(*[set(x) for x in sampleID_to_svtype_to_df.values()]).intersection(expected_svs)}

    # initialize a dict that will store all the vars that each sample has as text
    sample_ID_to_vars = {s:set() for s in sampleID_to_svtype_to_df.keys()}

    # go through each sample and svtype
    for sampleID, svtype_to_df in sampleID_to_svtype_to_df.items():
        for svtype, df in svtype_to_df.items(): 

            if svtype not in expected_svs: continue

            # define the already defined svtype
            df_defined_svs = svtype_to_already_found_dfs[svtype]

            # define the equal and approximate fields
            equal_fields = svtype_to_fieldsDict[svtype]["equal_fields"]
            approximate_fields = svtype_to_fieldsDict[svtype]["approximate_fields"]

            # add a "ID_str" to each sv in df that might be defined by the sv that overlaps in df_defined_svs, or new
            df["IDs_str"] = df.apply(lambda r: get_IDs_str_SV_deppending_on_overlapping(r, df_defined_svs, svtype, equal_fields, approximate_fields, tol_bp=tol_bp), axis=1)

            # keep all the ID_as_str
            sample_ID_to_vars[sampleID].update(set.union(*df["IDs_str"]))

            # keep df_defined_svs
            svtype_to_already_found_dfs[svtype] = df_defined_svs.append(df)

    # get a df mapping sample to var to present (1) or absent (0)
    all_vars = set.union(*sample_ID_to_vars.values())
    df_overlapping_SVs = pd.DataFrame({sID : {sv : int(sv in sample_svs) for sv in all_vars} for sID, sample_svs in sample_ID_to_vars.items()}).transpose()

    print("There are %i SVs found in more than 1 sample"%(len([col for col in df_overlapping_SVs.columns if sum(df_overlapping_SVs[col]==1) not in {1, 0}])))

    ############################

    return df_overlapping_SVs


#################################################################
#################################################################
#################################################################


#################################################################
##################### FUNCTIONS TO PLOT SVs #####################
#################################################################

    
def get_colors_df(list_samples, list_ColorDicts, list_levels):

    """Takes a list of sample IDs where _ represents several levels and a list of dictionaries that map each of the levels to a color. Returns colors_df, which is a dataframe where the index is the sampleID(s) from VARSdf and the columns are the layers of complexity, each of them with a color. list_levels are the names of the levels"""

    # create a dict that has [level][sampleID] = color
    level_to_sample_to_color = {}
    for I, level_to_color in enumerate(list_ColorDicts):

        # go through each sample
        for sample in list_samples:

            # get the different levels
            levels = sample.split("_")

            # add to dict
            level_to_sample_to_color.setdefault(list_levels[I], {}).setdefault(sample, level_to_color[levels[I]])

    return pd.DataFrame(level_to_sample_to_color).loc[list_samples]

def get_uniqueVals_df(df): return set.union(*[set(df[col]) for col in df.columns])

def get_Labels_heatmapObject_rows(Labels_colors_df, opacity=1, numeric_ticks=True):

    """Takes a df with color row names, where the index are Labels and the columns are the depth of the levels. The values are colors. It returns a "data" object (fig_data) that can be appended to a plotly figure object with figure.append_trace(fig_data) """

    # redefine depending on the orientation
    df = Labels_colors_df.transpose()
    if numeric_ticks: df.columns = reversed(list(range(len(df.columns)))) # this is to set everything to numbers

    # get all the colors
    all_colors = get_uniqueVals_df(df)
    
    # map them to numbers
    color_to_colnumber = dict(zip(all_colors, range(0, len(all_colors))))
    colnumber_to_color = {v:k for k,v in color_to_colnumber.items()}

    # define df as numbers
    df_numbers = df.applymap(lambda x: color_to_colnumber[x])

    # define the hover, from the Labels_colors_df
    hover = [["%s: %s <br>(%s)"%(col, row[Icol], "_".join(row)) for Icol, col in enumerate(Labels_colors_df.columns)] for Irow, row in enumerate(Labels_colors_df.index)]

    # get the figure, taking into account the real values
    fig_data = df_numbers.iplot(kind='heatmap', asFigure=True, colors=None, colorscale=False, colorbar=False, legend=False)["data"][0]

    # redefine the colorscale
    fig_data["colorscale"] = [[n/max(colnumber_to_color), colnumber_to_color[n]] for n in sorted(colnumber_to_color)]

    # change the opacity
    fig_data["opacity"] = opacity

    # set the hover
    fig_data["text"] = hover
    fig_data["hoverinfo"] = "text"

    # hide the colorbar
    fig_data["colorbar"] = dict(tickvals=[0,0], borderwidth=0, len=0, showticklabels=False, nticks=0, thickness=0, bgcolor="white")


    return fig_data


def define_colorbar(min, max, n=4, rgb_pos=0, type="btw_two_colors", color_from="red", color_to="blue"):
    
    """ Returns a dictionary in which the keys are a linspace array between min and max (with n steps) and the values are HEX color codes, progressively 
    rgb_pos indicates which of the rgb codes will change"""
    
    color_dict = {}
    
    if type=="btw_two_colors":
        
        # define the colors
        colors = list(Color(color_from).range_to(Color(color_to),n))
        
        for I, value in enumerate(np.linspace(min, max, n)):
            color_dict[value] = colors[I].get_hex()
        
    else:
        
        # this is a white to gray progression
        
        rgb_array = np.linspace(0, 255, n) # generate a linspace of the rgb color that will change
        rgb = [0, 0, 0] # this is going to be used
        color_dict = {}
        
        for I, value in enumerate(np.linspace(min, max, n)):
            
            # define the color
            rgb[0] = int(rgb_array[I]); rgb[1] = int(rgb_array[I]); rgb[2] = int(rgb_array[I]); # from black to gray
            #rgb[1] = int(rgb_array[I]); # from black to 
            
            hex_color = '#%02x%02x%02x'%(rgb[0], rgb[1], rgb[2])
            
            # add to dictionary
            color_dict[value] = hex_color
            
    return color_dict, list(color_dict.keys())

def get_colorbar_hex(iterable, color_from="red", color_to="blue"): 
    
    """ function to take list of elements and map them to different colors """

    return {k:v for k,v in dict(zip(iterable, define_colorbar(1,2, n=len(iterable), color_from=color_from, color_to=color_to)[0].values())).items()}


def get_rectangles_and_annotations_for_chromosomes(chrom_to_Xoffset, chrom_to_lenSeq, chrom_to_color, chrName_to_shortName, ylevel=0, width_rect=0.2, xref="x2", yref="y2", annotation_offset=0.1):

    """This function takes a dict mapping chromosomes to Xoffsets and returns a list of shape objects for plotly, where each shape is a rectangle.
    - ylevel indicates the center of the rectangle
    - width_rect is the +- width
    - xref and yref indicate in which axes they have to be located"""

    # initialize the lists 
    chromosome_rectanges = []
    chromosome_annotations = []

    print("There are %i chromosomes"%len(chrom_to_Xoffset))

    # go through each chromosome
    for I, (chrom, Xoffset) in enumerate(chrom_to_Xoffset.items()):

        # get metadata
        fillcolor = chrom_to_color[chrom]
        xstart = Xoffset+1
        xend = Xoffset+chrom_to_lenSeq[chrom]

        # get the rectangle dict
        chromosome_rectanges.append(dict(type="rect", x0=xstart, x1=xend, y0=ylevel-width_rect, y1=ylevel+width_rect, fillcolor=fillcolor, xref=xref, yref=yref, opacity=0.99, line=dict(color="black", width=3)))

        # get the annotation
        text = "<b>%s</b>"%chrName_to_shortName[chrom]
        chromosome_annotations.append(dict(x=(xstart + (xend-xstart)/2), y=(ylevel+width_rect+annotation_offset), text=text, textfont=dict(color=fillcolor, size=10)))

    return chromosome_rectanges, chromosome_annotations

def get_rectangles_and_annotations_for_gff_features(gff_df, chrom_to_Xoffset, interesting_genes="all", initial_ylevel=1, width_rect=0.2, xref="x2", yref="y2", interesting_features="all", annotation_offset=0.1, geneID_to_name={}):

    """This function takes a gff dataframe and gets the rectangels and text annotations above them. interesting_genes can be a set with the geneIDs to plot. If all, just take them all. It returns the rectangles, dictionaries of texts and the yLevel_to_GeneCoordinates, which will be used to adjust the figure layout. """

    # get interesting genes' df
    if interesting_genes=="all": interesting_genes = set(gff_df.upmost_parent)
    elif type(interesting_genes)==int: interesting_genes = set(rd.sample(list(set(gff_df.upmost_parent)), min([interesting_genes, len(set(gff_df.upmost_parent))])))

    if interesting_features=="all": interesting_features = set(gff_df.feature)

    print("There are %i genes and %i features"%(len(interesting_genes), len(interesting_features)))

    genes_df = gff_df[(gff_df.upmost_parent.isin(interesting_genes)) & (gff_df.feature.isin(interesting_features))].sort_values(by=["start", "end", "strand"])

    ######### map features to colors  #####

    # these are the gff-like features
    feature_to_color = {'gene': 'black', 'exon': 'darkblue', 'repeat_region': 'silver', 'mRNA': 'crimson', 'pseudogene': 'gray', 'centromere': 'silver', 'chromosome': 'silver', 'rRNA': 'coral', 'CDS': 'green', 'tRNA': 'red', 'ncRNA': 'orange', 'long_terminal_repeat': 'silver'}

    # add the InterPro features
    possible_colors = ["black", "maroon", "darkgreen", "darkslatergray", "darkblue", "indigo", "crimson"]
    for feat in set(genes_df.feature).difference(set(feature_to_color)): feature_to_color[feat] = rd.choice(possible_colors)
    
    #########################################

    # modify the start and the end to get adjusted to the chromosome of procedence
    genes_df["start"] = genes_df.apply(lambda r: r["start"]+chrom_to_Xoffset[r["chromosome"]], axis=1)
    genes_df["end"] = genes_df.apply(lambda r: r["end"]+chrom_to_Xoffset[r["chromosome"]], axis=1)

    # define the geneID_to_name
    if geneID_to_name=={}: geneID_to_name = {g:g for g in set(genes_df.upmost_parent)}

    ####### get the rectangles and annotations ########

    # initialize a dictionary  that will map the ylevel where to draw the genes vs the start and end of these genes
    yLevel_to_GeneCoordinates = {}

    # initialize a dictionary that will map each cds to the level where it is
    cdsID_to_yLevel  = {}

    # initialize the rectangles and annotations list
    features_rectangles = []
    features_annotations = [] 

    # go through all items of the gff without overlapping, as indicated by yLevel_to_GeneCoordinates
    print("adding features")
    for ID, source, feature, start, end, strand, upmost_parent in genes_df[["ID", "source", "feature", "start", "end", "strand", "upmost_parent"]].values:

        # define the level in which it has to be represented, from 0 to 100
        for level in range(initial_ylevel, 100):

            # if ot explored before
            if level not in yLevel_to_GeneCoordinates: good_level = level; break

            # if it does not overlap with anything in this level
            if all([end<=l_start or start>l_end for l_start, l_end in yLevel_to_GeneCoordinates[level]]): good_level = level; break

        # keep the cordinates
        yLevel_to_GeneCoordinates.setdefault(good_level, set())
        yLevel_to_GeneCoordinates[good_level].add((start, end))

        # represent the feature as a rectangle
        color_feature = feature_to_color[feature]
        rect =  dict(type="rect", x0=start, x1=end, y0=good_level-width_rect, y1=good_level+width_rect, fillcolor=color_feature, xref=xref, yref=yref, opacity=0.99, line=dict(color=color_feature, width=2))
        features_rectangles.append(rect)

        # get the annotation
        direction = {"+":">", "-":"<", ".":"?"}[strand]
        direction_text = direction*1
        complete_text = "<b>%s%s-%s%s</b>"%(direction_text, feature, geneID_to_name[upmost_parent], direction_text)
        features_annotations.append(dict(x=(start + (end-start)/2), y=(good_level+width_rect+annotation_offset), text=complete_text, textfont=dict(color=color_feature, size=9)))

        if feature=="CDS": cdsID_to_yLevel[ID] = good_level

    ##########################################################

    return features_rectangles , features_annotations, yLevel_to_GeneCoordinates


def get_scatter_line_changing_on_hover(xstart, xend, ystart, yend, color, hovertext, npoints=5, mode="lines+markers", opacity=0.9, linewidth=1.0, dash=None):

    """returns a go.Scatter object that is a line going from x to y and displaying hovertext"""

    line = go.Scatter(x=np.linspace(xstart, xend, npoints), y=np.linspace(ystart, yend, npoints), showlegend=False, mode=mode, text=hovertext, line=dict(color=color, width=linewidth, dash=dash), opacity=opacity, hoveron="points+fills") # dash="dash"

    return line

def format_info_hoverlabel(info, chrName_to_shortName, max_len_seq=5):

    """This function takes the info of a breakend to be added to a hoverlabel and formats it"""

    # initialize
    str_to_return = str(info)

    # floats
    if type(info)==float: str_to_return = "%.3f"%info

    # sequences
    elif type(info)==str and all([x.upper() in {"A", "C", "G", "T", "N"} for x in info]):

        # edit long seqs
        if len(info)>max_len_seq: str_to_return = "%ibp"%(len(info))
        else: str_to_return = info

    # chromosomes
    elif type(info)==str:

        chromosomes_in_info = set([c for c in chrName_to_shortName if c in info])
        if len(chromosomes_in_info)>0:

            for c in chromosomes_in_info: info = info.replace(c, chrName_to_shortName[c])
            str_to_return  = info

    return str_to_return

def get_middle_point_between_regions(ChrA, StartA, EndA, ChrB, StartB, EndB, chrom_to_Xoffset):

    """Takes two genomic regions and a chromosome to offset value, and finds the intermediate point, with 1-based coords"""

    # get real coords
    real_StartA = StartA + chrom_to_Xoffset[ChrA] + 1 
    real_StartB = StartB + chrom_to_Xoffset[ChrB] + 1 
    real_EndA = EndA + chrom_to_Xoffset[ChrA]
    real_EndB = EndB + chrom_to_Xoffset[ChrB]

    # get the intermediate positions
    intermediate_positions = sorted([real_StartA, real_StartB, real_EndA, real_EndB])[1:3]


    return intermediate_positions[0] + (intermediate_positions[1] - intermediate_positions[0])/2



def get_rearrangement_visualization_data(svDF, svtype, chrom_to_Xoffset, ylevel, chrName_to_shortName, list_clusters, xref="x2", yref="y1", width_rect=0.5, interesting_hover_properties=["real_AF", "FILTER", "inserted_sequence"]):

    """This function returns a list with all the shapes of a given svDF and a list with the scatter line traces that are displayed only on hover.

    interesting_hover_properties can be any in ['#CHROM', 'POS', 'other_coordinates', 'allele_frequency', 'allele_frequency_SmallEvent', 'real_AF', 'FILTER', 'inserted_sequence', 'has_poly16GC', 'length_inexactHomology', 'length_microHomology'].

    list_clusters is a list of sets, each of them containing IDs that are overlapping """


    # general things
    if len(list_clusters)>0: all_IDs_in_cluster = set.union(*list_clusters)
    else: all_IDs_in_cluster = set()
    all_IDs_without_cluster = set(svDF.uniqueID).difference(all_IDs_in_cluster)
    bool_to_text = {True: "", False:"not "}

    ######## map each uniqueID to the ylevel to which it corresponds, according to the cluster ########

    # get the ylevels of clusters, and also the width of the rects
    ID_to_ylevel = {}
    ID_to_width_rect = {}
    for cluster in list_clusters: 

        # define the width_rect
        cluster_width_rect = width_rect/len(cluster)

        # define the minimum and maxiumum ylevels according to width_rect
        min_ylevel = ylevel - (width_rect-cluster_width_rect)
        max_ylevel = ylevel + (width_rect-cluster_width_rect)

        # get the ylevels
        ID_to_width_rect = {**{ID:cluster_width_rect for ID in cluster}, **ID_to_width_rect}

        # get the widths
        ID_to_ylevel = {**dict(zip(list(cluster), np.linspace(min_ylevel, max_ylevel, len(cluster)))), **ID_to_ylevel}


    # get the ylevels of others, that are not in clusters
    ID_to_ylevel = {**{ID:ylevel for ID in all_IDs_without_cluster}, **ID_to_ylevel}
    ID_to_width_rect = {**{ID:width_rect for ID in all_IDs_without_cluster}, **ID_to_width_rect}

    #############################################################################################

    # initialize shapes
    rearrangements_shapes = [] 
    traces_displayed_onhover = []

    # add a hover text based on the bends_metadata_dict:
    svDF["bends_metadata_info_text"] = svDF.bends_metadata_dict.apply(lambda d: "<br>".join(["<br>%s:<br>bend1:<br>%s"%(bpid, "<br>bend2:<br>".join(["<br>".join(["   %s:%s   "%(typeinfo, format_info_hoverlabel(bend_info[typeinfo], chrName_to_shortName)) for typeinfo in interesting_hover_properties]) for Ibend, bend_info in enumerate(bend_info_list)])) for bpid, bend_info_list in d.items()]))

    # get interesting labels:
    svtype_to_color =  {"deletions":"red", "tandemDuplications":"blue", "inversions":"green"}

    ##### simple variations #####
    if svtype in {"deletions", "tandemDuplications", "inversions"}:

        # get a rectangle for each
        for I, (ID, Chr, Start, End, bends_metadata_info_text) in enumerate(svDF[["uniqueID", "Chr", "Start", "End", "bends_metadata_info_text"]].values):

            # get the coordinates (they are 0-based)
            real_start = Start + chrom_to_Xoffset[Chr] + 1
            real_end = End + chrom_to_Xoffset[Chr]
            color = svtype_to_color[svtype]
            real_ylevel = ID_to_ylevel[ID]
            id_width_rect = ID_to_width_rect[ID]
            
            # set to less opaque those IDs overlapping
            if ID in all_IDs_in_cluster: opacity=1
            else: opacity = 1

            # get a single rectangle
            rect =  dict(type="rect", x0=real_start, x1=real_end, y0=real_ylevel-id_width_rect, y1=real_ylevel+id_width_rect, fillcolor=color, xref=xref, yref=yref, opacity=opacity, line=dict(color=color, width=0.5))
            rearrangements_shapes.append(rect)

            # get the scatter object to be displayed only on hover, which is a line from real_start to real_end
            hovertext = "%s: %s:%i-%i<br>breakpoint info:<br>%s"%(ID, Chr, Start, End, bends_metadata_info_text)
            scatter_line_changing_on_hover = get_scatter_line_changing_on_hover(real_start, real_end, real_ylevel, real_ylevel, color, hovertext, mode="lines+markers", opacity=opacity)
            traces_displayed_onhover.append(scatter_line_changing_on_hover)

    #############################

    ##### translocations ########
    if svtype=="translocations":

        # go through each event
        for I, (ID, ChrA, StartA, EndA, ChrB, StartB, EndB, bends_metadata_info_text, Balanced) in enumerate(svDF[["uniqueID", "ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "bends_metadata_info_text", "Balanced"]].values):

            # define the fillcolor of each, depending on if it is balanced
            if Balanced is True: 

                # define whether it is an exchange between 5' and 5', which means that there is no inversion. The outline color reflects this
                if (StartA==0 and StartB==0) or (StartA!=0 and StartB!=0): outlinecolor = "gray"
                else: outlinecolor = "olive"

                # define the fill color
                fillcolorA = fillcolorB = outlinecolor
                outlinecolorA = outlinecolorB = outlinecolor
                opacity_inner = 1

            else: 
                fillcolorA = outlinecolorA = "steelblue" # this is copied
                fillcolorB = outlinecolorB = "magenta" # this one is deleted
                opacity_inner = 1

            # define general things
            id_width_rect = ID_to_width_rect[ID]
            opacity_line = 1
            real_ylevel = ID_to_ylevel[ID]
            hovertext = "%s: %s:%i-%i to %s:%i-%i (%sBalanced)<br>breakpoint info:<br>%s"%(ID, format_info_hoverlabel(ChrA, chrName_to_shortName), StartA, EndA, format_info_hoverlabel(ChrB, chrName_to_shortName), StartB, EndB, bool_to_text[Balanced], bends_metadata_info_text)
            middle_point_two_regions = get_middle_point_between_regions(ChrA, StartA, EndA, ChrB, StartB, EndB, chrom_to_Xoffset)# the middle point between the end 

            # make a rectangle for each location
            for Chr, Start, End, fillcolor, region, outlinecolor in [[ChrA, StartA, EndA, fillcolorA, "A", outlinecolorA], [ChrB, StartB, EndB, fillcolorB, "B", outlinecolorB]]: 

                # change coords
                real_start = Start + chrom_to_Xoffset[Chr] + 1
                real_end = End + chrom_to_Xoffset[Chr]

                # get the rectangle of the the inner
                rect =  dict(type="rect", x0=real_start, x1=real_end, y0=real_ylevel-id_width_rect, y1=real_ylevel+id_width_rect, fillcolor=fillcolor, xref=xref, yref=yref, opacity=opacity_inner, line=dict(color=outlinecolor, width=0.1))
                rearrangements_shapes.append(rect)

                # get the rectangle of the outline
                rect =  dict(type="rect", x0=real_start, x1=real_end, y0=real_ylevel-id_width_rect, y1=real_ylevel+id_width_rect, fillcolor=None, xref=xref, yref=yref, opacity=opacity_line, line=dict(color=outlinecolor, width=1.5))
                rearrangements_shapes.append(rect)

                # get the line with the hovertext
                scatter_line_changing_on_hover = get_scatter_line_changing_on_hover(real_start, real_end, real_ylevel, real_ylevel, fillcolor, hovertext, mode="lines+markers", opacity=opacity_line, npoints=3)
                traces_displayed_onhover.append(scatter_line_changing_on_hover)

                # get a line from the middle to the middle point
                mid_point = real_start + (real_end-real_start)/2
                scatter_line_to_the_mid = get_scatter_line_changing_on_hover(mid_point, middle_point_two_regions, real_ylevel, real_ylevel, fillcolor, hovertext, mode="lines", opacity=1, npoints=2, dash="dot", linewidth=0.6) # real_ylevel+id_width_rect
                traces_displayed_onhover.append(scatter_line_to_the_mid)

    #############################

    ###### insertions ##########
    if svtype=="insertions":

        # go through each event
        for I, (ID, ChrA, StartA, EndA, ChrB, StartB, EndB, bends_metadata_info_text, Copied) in enumerate(svDF[["uniqueID", "ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "bends_metadata_info_text", "Copied"]].values):

            # the destiny is always in black, which reflects that there is not a lot of importance
            if Copied is True: colorA = "dodgerblue"
            else: colorA = "maroon"
            colorB = "black"

            # define general things
            real_ylevel = ID_to_ylevel[ID]
            id_width_rect = ID_to_width_rect[ID]
            hovertext = "%s: %s:%i-%i to %s:%i-%i (%sCopied)<br>breakpoint info:<br>%s"%(ID, format_info_hoverlabel(ChrA, chrName_to_shortName), StartA, EndA, format_info_hoverlabel(ChrB, chrName_to_shortName), StartB, EndB, bool_to_text[Copied], bends_metadata_info_text)
            middle_point_two_regions = get_middle_point_between_regions(ChrA, StartA, EndA, ChrB, StartB, EndB, chrom_to_Xoffset)# the middle point between the end 

            # for each region, plot the rect and the lines
            for Chr, Start, End, color in [[ChrA, StartA, EndA, colorA], [ChrB, StartB, EndB, colorB]]: 

                # change coords
                real_start = Start + chrom_to_Xoffset[Chr] + 1
                real_end = End + chrom_to_Xoffset[Chr]

                # get the rectangle of the the inner
                rect =  dict(type="rect", x0=real_start, x1=real_end, y0=real_ylevel-id_width_rect, y1=real_ylevel+id_width_rect, fillcolor=color, xref=xref, yref=yref, opacity=1, line=dict(color=color, width=0.1))
                rearrangements_shapes.append(rect)

          

                # get the line with the hovertext
                scatter_line_changing_on_hover = get_scatter_line_changing_on_hover(real_start, real_end, real_ylevel, real_ylevel, color, hovertext, mode="lines+markers", opacity=1, npoints=3)
                traces_displayed_onhover.append(scatter_line_changing_on_hover)

                # get a line from the middle to the middle point
                mid_point = real_start + (real_end-real_start)/2
                scatter_line_to_the_mid = get_scatter_line_changing_on_hover(mid_point, middle_point_two_regions, real_ylevel, real_ylevel, color, hovertext, mode="lines", opacity=1, npoints=2, dash="dot", linewidth=0.6) # real_ylevel+id_width_rect
                traces_displayed_onhover.append(scatter_line_to_the_mid)

    #############################

    return rearrangements_shapes, traces_displayed_onhover

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

    print("There are %i clusters"%len(list_clusters))

    # make sure that all clusters are nonoverlapping with each other
    for I in range(len(list_clusters)):
        for J in range(I+1, len(list_clusters)):

            if I==J: raise ValueError("We don't want to compare the same I and J")

            # define clusters
            clusterI = list_clusters[I]
            clusterJ = list_clusters[J]

            if len(clusterI.intersection(clusterJ))!=0: 
                print(I, clusterI, "\n", J, clusterJ)
                #raise ValueError("There are some overlapping clusters")
                print("There are some overlapping clusters")


    return list_clusters


def get_clusters_overlapping_vars_allGenome(svtype_to_svDF, chrom_to_Xoffset):

    """This function takes a dictionary that maps svtype to the df and returns a list of sets. Each of them represents a cluster of variants that are overlapping. For variants where there is more than one position it is considered overlapping the whole set of variants"""

    print("getting clusters of overlapping vars, taking into account the whole chromosome")

    # define a function that changes the START by the END if it is -1
    def arrange_minus1_START(r):

        if r["START"] == -1: return r["END"]
        else: return r["START"]

    # initialize a df that will have Chr, Start, End, ID for all vars
    all_regions_df = pd.DataFrame()

    # go through each svtype
    for svtype, svDF in svtype_to_svDF.items():

        # keep
        svDF = cp.deepcopy(svDF)

        ### define a svDF that will include real_start and real_end relative to chrom_to_Xoffset

        # simple types, only one box
        if svtype in {"deletions", "tandemDuplications", "inversions"}: 
            
            svDF["real_start"] = svDF.apply(lambda r: chrom_to_Xoffset[r["Chr"]] + r["Start"], axis=1)
            svDF["real_end"] = svDF.apply(lambda r: chrom_to_Xoffset[r["Chr"]] + r["End"], axis=1)

        # two boxes
        elif svtype in {"insertions", "translocations"}:

            # get the min and the max of each
            svDF["real_start_A"] = svDF.apply(lambda r: chrom_to_Xoffset[r["ChrA"]] + r["StartA"], axis=1)
            svDF["real_start_B"] = svDF.apply(lambda r: chrom_to_Xoffset[r["ChrB"]] + r["StartB"], axis=1)
            svDF["real_end_A"] = svDF.apply(lambda r: chrom_to_Xoffset[r["ChrA"]] + r["EndA"], axis=1)
            svDF["real_end_B"] = svDF.apply(lambda r: chrom_to_Xoffset[r["ChrB"]] + r["EndB"], axis=1)

            # keep min and max
            svDF["real_start"] = svDF.apply(lambda r: min(r["real_start_A"], r["real_start_B"]), axis=1)
            svDF["real_end"] = svDF.apply(lambda r: max(r["real_end_A"], r["real_end_B"]), axis=1)

        elif svtype=="remaining":

            # get the coords of each
            svDF["corrected_START"] = svDF.apply(arrange_minus1_START, axis=1)
            svDF["real_pos"] = svDF.apply(lambda r: chrom_to_Xoffset[r["#CHROM"]] + r["POS"], axis=1)
            svDF["real_START"] = svDF.apply(lambda r: chrom_to_Xoffset[r["CHR2"]] + r["corrected_START"], axis=1)
            svDF["real_END"] = svDF.apply(lambda r: chrom_to_Xoffset[r["CHR2"]] + r["END"], axis=1)

            # get min and max
            svDF["real_start"] = svDF.apply(lambda r: min(r["real_pos"], r["real_START"], r["real_END"]), axis=1)
            svDF["real_end"] = svDF.apply(lambda r: max(r["real_pos"], r["real_START"], r["real_END"]), axis=1)

        else: raise ValueError("%s is not a valid svtype"%svtype)

        # keep
        all_regions_df = all_regions_df.append(svDF[["real_start", "real_end", "uniqueID"]])

        ######

    # check everything is fine
    if len(all_regions_df)!=len(set(all_regions_df.uniqueID)): raise ValueError("There are some non-unique eventIDs across the svtypes")

    # now identify clusters of overlapping IDs, and save the IDs

    # map each ID to the overlapping IDs
    ID_to_overlappingIDs = {}

    # go through each region
    for qStart, qEnd, qID in all_regions_df[["real_start", "real_end", "uniqueID"]].values:

        # find overlapping SVs
        overlappingIDs = set(all_regions_df[((all_regions_df.real_start>=qStart) & (all_regions_df.real_end<=qEnd)) | 
                                  ((all_regions_df.real_start<=qStart) & (all_regions_df.real_end>=qStart)) |  
                                  ((all_regions_df.real_start<=qEnd) & (all_regions_df.real_end>qEnd))]["uniqueID"])

        # if there are overlaps add them to the cluster_list with itself
        ID_to_overlappingIDs[qID] =  overlappingIDs

    # initialize clusters list
    list_clusters = get_list_clusters_from_dict(ID_to_overlappingIDs)

    # check
    all_SVs = set(all_regions_df.uniqueID)
    if len(list_clusters)>0: all_SVs_in_cluster = set.union(*list_clusters)
    else: all_SVs_in_cluster = set()

    print("There are %i clusters of SVs involving %i of %i SVs"%(len(list_clusters), len(all_SVs_in_cluster), len(all_SVs)))

    return list_clusters



def format_svDF(svDF, svtype="unknown", interesting_chromosomes="all"):

    """Takes an svDF and formats the ID if suspicious"""

    if set(svDF.ID)=={""}: svDF["ID"] = svDF.Name

    # add a uniqueID
    svDF["uniqueID"] = ["%s_%i"%(svtype, I+1) for I in range(len(svDF))]

    # filter and get only chromosomes that are in interesting_chromosomes
    if interesting_chromosomes!="all":

        # define the chromosomal fields
        chromosome_fields = {"Chr", "ChrA", "ChrB", "#CHROM", "CHR2"}.intersection(set(svDF.keys()))

        # get the df where the interesting_chromosomes are involved
        svDF = svDF[svDF.apply(lambda r: any([r[c] in interesting_chromosomes for c in chromosome_fields]), axis=1)]

    return svDF

def plot_SV_browser(sampleID_to_svtype_to_svDF, filename, samples_colors_df, gff_df, genome_path, legendGroup_to_samples, interesting_genes="all", interesting_features="all", title="SV browser", chrName_to_shortName={}, geneID_to_name={}, fraction_y_domain_by_browser=0.8, interesting_svtypes={"deletions", "tandemDuplications", "inversions", "translocations", "insertions"}, interesting_chromosomes="all"):

    """
    This function plots all the structural variations into an interactive html plot, in order to browse among them. These are the arguments:
    - sampleID_to_svtype_to_svDF: is a dictionary that maps each of the samples to the svtypes and a df with them. These are all the variants that will be ploted
    - filename is where the plot will be generated (it has to be an html)
    - samples_colors_df is a dataframe where the index is the sampleID(s) and the columns are the layers of complexity, each of them with a color. This defines the grid of the samples to include in the plot
    - gff_df a df returned by load_gff3_intoDF
    - genome_path a fasta with the genome
    - legendGroup_to_samples is a dictionary that maps each legendgroup to a set of samples, which is useful to plot the mutations in a consistent way.
    - interesting_genes may be a set of the genes for which to plot annotations, and "all would indicate"
    - chrName_to_shortName maps each chromosome to a reduced name
    - interesting_features can be a set of the gff features to plot
    - is how the upmost_parent (which is the gene) of the gff has to be represented
    - fraction_y_domain_by_browser defines which fraction of the ydomain has to be taken by the genomic browser. If it is 0.5 it will take up half of the ydomain.
    - interesting_svtypes is a set of the types of SVs of interest. Can be any of {"deletions", "tandemDuplications", "inversions", "translocations", "insertions"}
    - can be a set with all the interesting chromosomes or "all"
    """

    # load genome for those chromosomes in sampleID_to_svtype_to_svDF and get it's length, which will be used for plotting
    all_chromosomes = set.union(*[set.union(*[get_uniqueVals_df(svDF[[col for col in svDF.keys() if col in {"ChrA", "ChrB", "#CHROM", "Chr", "CHR2"}]]) for svDF in svtype_to_svDF.values()]) for svtype_to_svDF in sampleID_to_svtype_to_svDF.values()])
    if interesting_chromosomes!="all": all_chromosomes = all_chromosomes.intersection(interesting_chromosomes)

    chrom_to_lenSeq = {seq.id : len(str(seq.seq)) for seq in SeqIO.parse(genome_path, "fasta") if seq.id in all_chromosomes}

    # define chromosome offset, which is necessary to keep the xaxis organised
    chrom_to_Xoffset = {}
    current_offset = 0
    for chrom, seqLen in chrom_to_lenSeq.items():
        chrom_to_Xoffset[chrom] = current_offset
        current_offset+=seqLen

    # filter the gff df by chromosome if needed
    gff_df = gff_df[gff_df.chromosome.isin(all_chromosomes)]

    # define the colors of the chromosomes

    # add the InterPro features
    colors = ['black', 'darkblue', 'crimson', 'gray', 'coral', 'green', 'red', 'orange', 'maroon', 'darkgreen', 'purple', 'indigo', 'magenta', 'darkred']
    #colors = list(get_colorbar_hex(list(range(len(chrom_to_Xoffset)))).values()); rd.shuffle(colors) # with a consistent ordering
    #colors = rd.sample(list(webcolors.css3_hex_to_names.values()), len(chrom_to_Xoffset)) # these are just shuffled
    chrom_to_color = dict(zip(chrom_to_Xoffset, colors))

    # change the chrName_to_shortName if not provided
    if chrName_to_shortName=={}: chrName_to_shortName = {c:c.split("_")[0] for c in all_chromosomes}

    # setup the layout (2 rows and 2 cols)
    fig = tools.make_subplots(rows=2, cols=2, specs=[[{}, {}], [{}, {}]], vertical_spacing=0.0, horizontal_spacing=0.0, subplot_titles=("", "", "", ""), shared_yaxes=True, shared_xaxes=True, print_grid=True)

    # initialize shapes list
    shapes = []

    ######## 1,1: sample levels ######
    fig.append_trace(get_Labels_heatmapObject_rows(samples_colors_df), 1, 1)

    # map each ylevel to the corresponding label
    ylevels = list(reversed(range(len(samples_colors_df))))
    ylabels = ["_".join(x) for x in samples_colors_df.index]
    ylabel_to_ylevel = dict(zip(ylabels, ylevels))

    # check that all labels are in the colors_df
    if any([x not in ylabels for x in sampleID_to_svtype_to_svDF]): raise ValueError("Some of the samples are not added to the grid")

    ##################################

    ######## 2,2: gene browser ######

    # initialize the features text
    features_text_dicts = [] # this is a list of dictionaries, each of which has the info to get a go.Scatter text object

    # get the rectangles of the chromosomes
    chromosome_rectanges, chromosome_annotations = get_rectangles_and_annotations_for_chromosomes(chrom_to_Xoffset, chrom_to_lenSeq, chrom_to_color, chrName_to_shortName, ylevel=0, width_rect=0.2, annotation_offset=0.1)

    # keep them 
    shapes += chromosome_rectanges

    # get the gff features
    features_rectangles , features_annotations, yLevel_to_GeneCoordinates = get_rectangles_and_annotations_for_gff_features(gff_df, chrom_to_Xoffset, interesting_genes=interesting_genes, interesting_features=interesting_features)

    # keep 
    shapes += features_rectangles

    # keep the text as trace, one with all
    for label, features_text_dicts in [["chromosomes", chromosome_annotations], ["annotations", features_annotations]]:

        # get the data structures to define the plots
        keys_first_text = features_text_dicts[0].keys()
        dict_all_texts = {key : [text_dict[key] for text_dict in features_text_dicts] for key in keys_first_text}
        keys_first_textfont_dict = dict_all_texts["textfont"][0].keys()
        textfont_dict_all = {key:[textfont_dict[key] for textfont_dict in dict_all_texts["textfont"]] for key in keys_first_textfont_dict}

        trace_texts = go.Scatter(x=dict_all_texts["x"], y=dict_all_texts["y"], text=dict_all_texts["text"], mode="text", showlegend=True, name=label, textfont=textfont_dict_all, textposition="center", visible="legendonly")
        fig.append_trace(trace_texts, 2, 2)
       
    ##################################

    ######## 1,2: rearrangements browser ######

    # add the background, related by the last column of samples_colors_df
    lastCol = samples_colors_df.columns[-1]
    lastChr = sorted(chrom_to_Xoffset.items(), key=(lambda x: x[1]))[-1][0]
    end_lastChr = chrom_to_Xoffset[lastChr] + chrom_to_Xoffset[lastChr]
    background_df = pd.DataFrame({x : {s : samples_colors_df.loc[{s}, lastCol].iloc[0] for s in samples_colors_df.index} for x in [0, end_lastChr] }).loc[samples_colors_df.index]
    fig.append_trace(get_Labels_heatmapObject_rows(background_df, opacity=0.15, numeric_ticks=True), 1, 2)

    # go through each type of data and return the shapes' list and the scatter traces that cover them, which will include the hover
    for Is,(sample, svtype_to_svDF) in enumerate(sampleID_to_svtype_to_svDF.items()):
        print(sample)
        #if Is!=1: continue # debug

        # filter to get vars that are in interesting_svtypes
        svtype_to_svDF_filt = {svtype : format_svDF(svDF, svtype=svtype, interesting_chromosomes=interesting_chromosomes) for svtype, svDF in svtype_to_svDF.items() if svtype in interesting_svtypes}
        svtype_to_svDF_filt = {svtype:svDF for svtype, svDF in svtype_to_svDF_filt.items() if len(svDF)>0}

        # define a set of tuples of variants that are under overlap
        #list_clusters = get_clusters_overlapping_vars(svtype_to_svDF_filt)
        list_clusters = get_clusters_overlapping_vars_allGenome(svtype_to_svDF_filt, chrom_to_Xoffset) # this is including the distance in the Xgrid

        for svtype, svDF in svtype_to_svDF_filt.items():
            print(sample, svtype)

            # get the plotting data
            rearrangements_shapes, traces_displayed_onhover = get_rearrangement_visualization_data(svDF, svtype, chrom_to_Xoffset, ylabel_to_ylevel[sample], chrName_to_shortName, list_clusters)

            # keep shapes
            shapes += rearrangements_shapes

            # keep the lines displayed on hover
            for line in traces_displayed_onhover: fig.append_trace(line, 1, 2)

 
    ###########################################

    # formating the yticks has to be done in the layout like in https://plot.ly/python/v3/tick-formatting/

    ###### FIGURE LAYOUT ######

    # define domains of axes
    width_samplelevels = len(samples_colors_df.columns)*0.02
    xdomain = [0, width_samplelevels] # this is the domain of the samples x axis
    xdomain2 = [width_samplelevels, 1] # domain browser
    witdth_browser = max([1 - len(samples_colors_df)*0.1, 0.1]) # this is depending on the mutations
    #witdth_browser = max(yLevel_to_GeneCoordinates)*0.03 # this is depending on the broser
    #witdth_browser = fraction_y_domain_by_browser # this is using the given fraction
    ydomain = [witdth_browser, 1] # mutations
    ydomain2 = [0, witdth_browser] # exons

    # add the layout
    yaxis = {"domain":ydomain, "tickmode":"array", "tickvals":ylevels, "ticktext":ylabels}

    fig['layout'].update(title=title,  xaxis={"domain":xdomain},  xaxis2={"domain":xdomain2}, yaxis=yaxis,  yaxis2={"domain":ydomain2}, font={"size":12}, shapes=shapes, margin=go.Margin(l=250, r=250, b=50, t=50), showlegend=True)

    # get the fig
    init_notebook_mode()
    config={'editable': False}
    print("Writing %s"%filename)
    off_py.plot(fig, filename=filename, auto_open=False, config=config)

    return fig


# define graphical properties
Cglabrata_Nanopore_Assemblies_strain_to_color = {'realData_EB0911Sto': 'lightgreen', 'realData_M12': 'blueviolet', 'realData_CBS138': 'gray', 'realData_P35-2': 'blue', 'realData_BG2': 'darkorange', 'realData_EF1237Blo1': 'cyan'}
svtype_to_color = {'remaining': 'gray', 'tandemDuplications': 'blue', 'insertions': 'cyan', 'deletions': 'red', 'inversions': 'green', 'translocations': 'black'}

def plot_clustermap_for_SVdict(sampleID_to_SV_dict, fileprefix, svtype_to_color=svtype_to_color, sample_to_color=Cglabrata_Nanopore_Assemblies_strain_to_color, consider_only_varying_SVs=False):

    """This function takes a sampleID_to_SV_dict and plots a dendrogram where the rows are samples and the columns are the events"""

    # define svs
    all_svs = {'inversions', 'remaining', 'tandemDuplications', 'translocations', 'insertions', 'deletions'}
    simple_svs = {'tandemDuplications', 'deletions', 'inversions'}
    understanded_svs = {'inversions', 'tandemDuplications', 'translocations', 'insertions', 'deletions'}
    complex_svs = {'translocations', 'insertions'}

    #list_expected_svs = [all_svs, simple_svs, understanded_svs, complex_svs] + [{x} for x in all_svs]
    list_expected_svs = [understanded_svs]

    all_methods = ["single", "complete", "average", "weighted"]


    # go through different SVs to plot
    for expected_svs in list_expected_svs:

        # get ID
        svsID = "_".join([x[0:3] for x in sorted(expected_svs)])
        print(svsID)

        # get df of overlapping vars
        df_overlapping_SV = get_df_overlapping_SV(sampleID_to_SV_dict, expected_svs=expected_svs, tol_bp=50)

        # only get df that is varying 
        if consider_only_varying_SVs is True:
            varying_cols = [col for col in df_overlapping_SV.columns if sum(df_overlapping_SV[col]==1) not in {1, 0}]
            df_overlapping_SV = df_overlapping_SV[varying_cols]

        # get row and column colors
        row_colors_df = pd.DataFrame({s:{"sample":c} for s,c in sample_to_color.items()}).transpose()
        col_colors_df = pd.DataFrame({sv : {"SVtype" : svtype_to_color[sv.split("_")[0]]} for sv in df_overlapping_SV.columns}).transpose()

        # get as clustermap
        figsize = (10, len(row_colors_df)*0.5)

        # go through several clustering algorithms
        for method in all_methods:

            # get clustermap
            cm = sns.clustermap(df_overlapping_SV, col_cluster=True, row_cluster=True, row_colors=row_colors_df, col_colors=col_colors_df, cbar_kws={'label': "presence/absence", "ticks":[0, 1]}, square=False, figsize=figsize, xticklabels=False, cmap=["darkorange", "indigo"], method=method, metric="hamming") # figsize=figsize, linecolor=linecolor, linewidths=linewidths, yticklabels=yticklabels, xticklabels=False,

            cm.fig.suptitle('%s, method:%s'%(svsID, method)) 

def plot_clustermap_with_annotation(df, row_colors_df, col_colors_df, filename, title="clustermap", col_cluster=False, row_cluster=False, colorbar_label="default label", adjust_position=True, legend=True, idxs_separator_pattern="_", texts_to_strip={"L001"}, default_label_legend="control", df_annotations=None, cmap=sns.color_palette("RdBu_r", 50), ylabels_graphics_df=None, grid_lines=True, add_to_legend_x=0.5):

    """Takes a df were the index is the annotation and the cols are samples. It will be saved under filename. ylabels_graphics_df can be a df containing fontweight and color for each index value in df"""

    # define figsize 
    figsize = (len(df.columns)*0.3, len(df)*0.35)

    print(figsize)

    # define the line color from the last item in the df
    if grid_lines is True: 
        linecolor = list(col_colors_df[list(col_colors_df.keys())[-1]])
        linewidths = 1.5

    else: 
        linecolor = "gray"
        linewidths = 0

    # define the yticklabels
    if ylabels_graphics_df is not None: yticklabels = True
    else: yticklabels = False

    # decide whether to add annotations
    if df_annotations is not None: 

        # reorder according to the clustering
        cm_labels = sns.clustermap(df, col_cluster=col_cluster, row_cluster=row_cluster, yticklabels=True, xticklabels=True)
        plt.close()

        xlabels = [x.get_text() for x in cm_labels.ax_heatmap.get_xticklabels()]
        ylabels = [y.get_text() for y in cm_labels.ax_heatmap.get_yticklabels()]

        annot = df_annotations.loc[ylabels][xlabels]

    else: annot = False

    if adjust_position is True:

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

        # get the clustermap
        cm = sns.clustermap(df, col_cluster=col_cluster, row_cluster=row_cluster, row_colors=row_colors_df, col_colors=col_colors_df, cbar_kws={'label': colorbar_label}, xticklabels=False, square=False, figsize=figsize, cmap=cmap, annot=annot, fmt="", annot_kws={"size": 6.5}, linecolor=linecolor, linewidths=linewidths, yticklabels=yticklabels) # figsize=figsize, linecolor=linecolor, 

        hm_pos = cm.ax_heatmap.get_position()

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
        cm.ax_heatmap.legend(loc="best", bbox_to_anchor=(hm_pos.x0+hm_pos.width+add_to_legend_x, hm_pos.y0, hm_pos.width, hm_pos.height), ncol=2) # loc is the lower left corner bbox_to_anchor
        #cm.ax_heatmap.legend(loc="best", bbox_to_anchor=(2, 1), ncol=2) # loc is the lower left corner bbox_to_anchor

        ###################################



    # add graphics to y labels if provided
    if ylabels_graphics_df is not None:

        for I, (fw, c) in enumerate(ylabels_graphics_df[["fontweight", "color"]].values):

            cm.ax_heatmap.get_yticklabels()[I].set_weight(fw) 
            cm.ax_heatmap.get_yticklabels()[I].set_color(c) 

    # print the positions of everything
    print("cm.ax_heatmap.get_position()", cm.ax_heatmap.get_position())
    print("cm.ax_row_colors.get_position()", cm.ax_row_colors.get_position())
    print("cm.ax_col_colors.get_position()", cm.ax_col_colors.get_position())
    print("cm.ax_row_dendrogram.get_position()", cm.ax_row_dendrogram.get_position())
    print("cm.ax_col_dendrogram.get_position()", cm.ax_col_dendrogram.get_position())

    # SAVE
    print("saving %s"%filename)
    cm.savefig(filename)

    return cm

    
def getPlots_filtering_accuracy_across_genomes_and_ploidies(df_cross_benchmark, PlotsDir, simName_to_color={**Cglabrata_Nanopore_Assemblies_strain_to_color, **{"simulation_1":"black", "simulation_2":"red", "simulation_3":"green"}}, simType_to_color={'biased_towards_repeats':"red", 'realData':"black", 'uniform':"blue"}, ploidy_to_color={'consensus_ref': 'gray', 'haploid': 'black', 'diploid_hetero': 'maroon', 'ref:3_var:1': 'red', 'ref:9_var:1': 'lightsalmon', 'ref:99_var:1': 'white'}, svtype_to_color={"tandemDuplications": "gray", "deletions": "black", "inversions": "blue", "translocations": "olive", "insertions": "red", "remaining":"magenta", "integrated":"c"}):

    """This function takes a df that has several training and testing genomes with accuracy in each, and makes several plots to represent the data under outdir """

    print("WARNING: you have to change the colors !!!")

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

    print(all_test_genomeIDs)
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

                            print(interesting_test_genomeIDs)

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
                                        if row==best_train_idx: label = "%s(B)"%label

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

                                        plots_dir_crossaccuracyHeatmaps = "%s/cross_accuracy_heatmaps"%PlotsDir; fun.make_folder(plots_dir_crossaccuracyHeatmaps)
                                        string_title_train = "train_genomes:%s_ploidies:%s_svtypes:%s"%(train_genomeIDs, train_ploidies, train_svtypes)
                                        string_title_test = "test_genomes:%s_ploidies:%s_svtypes:%s"%(test_genomeIDs, test_ploidies, test_svtypes)

                                        #filename = "%s/accuracyHeatmap_%s_%s_rowClust%s_colClust%s_%s.pdf"%(plots_dir_crossaccuracyHeatmaps, string_title_train, string_title_test, row_cluster, col_cluster, accuracy)
                                        #filename = "%s/accuracyHeatmap_%s_%s_rowClust%s_colClust%s_%s.pdf"%(plots_dir_crossaccuracyHeatmaps, string_title_train, string_title_test, row_cluster, col_cluster, accuracy)

                                        filename = "%s/accuracyHeatmap_%s.pdf"%(plots_dir_crossaccuracyHeatmaps, accuracy)
                                        

                                        filename = filename.replace(":", "_")

                                        #filename = "%s/accuracyHeatmap_shortName.pdf"%(plots_dir_crossaccuracyHeatmaps)


                                        #title = "%s\n%s"%(string_title_train, string_title_test)
                                        title ="checking overfitting in SV calling"

                                        #print(df_square)

                                        print(sorted(get_uniqueVals_df(df_square)))

                                        # get the dendrogram, either adjusting or not
                                        try: plot_clustermap_with_annotation(df_square, row_colors_df, col_colors_df, filename, title=title, col_cluster=col_cluster, row_cluster=row_cluster, colorbar_label=accuracy, adjust_position=True, legend=True, idxs_separator_pattern="||||", texts_to_strip={"L001"}, default_label_legend="control", df_annotations=df_annotations, cmap=sns.color_palette("RdBu_r", 50), ylabels_graphics_df=None, grid_lines=False)

                                        except: plot_clustermap_with_annotation(df_square, row_colors_df, col_colors_df, filename, title=title, col_cluster=col_cluster, row_cluster=row_cluster, colorbar_label=accuracy, adjust_position=False, legend=True, idxs_separator_pattern="||||", texts_to_strip={"L001"}, default_label_legend="control", df_annotations=df_annotations, cmap=sns.color_palette("RdBu_r", 50), ylabels_graphics_df=None, grid_lines=False)
    

    # return the best_train_idx
    return dict(zip(["simName", "simType", "ploidy", "svtype"], best_train_idx.split("||||")))

    ##################################################################



#################################################################
#################################################################
#################################################################