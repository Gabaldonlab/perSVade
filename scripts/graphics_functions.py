#!/usr/bin/env python

# this contains all the functions related to the plotings of the genome variation browser

######################################################
###############  DEFINE ENVIRONMENT ##################
######################################################

# module imports
import os
import sys
import seaborn as sns
import pandas as pd
import numpy as np

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import sv_functions as fun



######################################################
######################################################
######################################################


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
        value_to_color = {v : palette_dict[fun.find_nearest(all_values_palette, v)] for v in values}

    return value_to_color, palette_dict


def load_gff3_intoDF(gff_path):

    """ Takes the path to a gff and loads it into a df"""

    gff_df_file = "%s.df.tab"%gff_path

    if fun.file_is_empty(gff_df_file):

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

        gff["upmost_parent"] = gff.apply(lambda row: get_utmost_parent(row, gff), axis=1)
        gff = gff.set_index("upmost_parent", drop=False)

        # write df
        gff_df_file_tmp = "%s.tmp"%gff_df_file
        gff.to_csv(gff_df_file_tmp, sep="\t", index=False, header=True)
        os.rename(gff_df_file_tmp, gff_df_file)


    # load
    gff = pd.read_csv(gff_df_file, sep="\t")

    return gff


def get_integrated_vars_df(df, cache_dir, target_regions, target_genes):

    """Takes the df with the data and returns an integrated dataframe with small_vars, small_vars_annot, SV_CNV, SV_CNV_annot. It will discard any variants that don't overlap target_regions and are not related to target_genes"""


    ######## GET INTEGRATED DFS ########

    # define a prefix for these samples
    file_prefix = "%s/%s"%(cache_dir, "-".join(df.sampleID))

    # define files
    small_vars_file = "%s_small_vars.tab"%file_prefix
    small_vars_annot_file = "%s_small_vars_annot.tab"%file_prefix
    SV_CNV_file = "%s_SV_CNV.tab"%file_prefix
    SV_CNV_annot_file = "%s_SV_CNV_annot.tab"%file_prefix

    # define the expected files
    expected_files = set()
    if "smallVars_vcf" in df.keys(): expected_files.update({small_vars_file, small_vars_annot_file})
    if "SV_CNV_vcf" in df.keys(): expected_files.update({SV_CNV_file, SV_CNV_annot_file})

    if any([fun.file_is_empty(f) for f in expected_files]):

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
                vcf = fun.get_vcf_df_with_INFO_as_single_fields(fun.get_df_and_header_from_vcf(r["smallVars_vcf"])[0])
                annotation = pd.read_csv(r["smallVars_var_annotation"], sep="\t")

                # add to df
                vcf["sampleID"] = sampleID

                small_vars = small_vars.append(vcf)
                small_vars_annot = small_vars_annot.append(annotation)

            # SVs
            if "SV_CNV_vcf" in df.keys():

                # get dfs
                vcf = fun.get_vcf_df_with_INFO_as_single_fields(fun.get_df_and_header_from_vcf(r["SV_CNV_vcf"])[0])
                annotation = pd.read_csv(r["SV_CNV_var_annotation"], sep="\t")

                # add to df
                vcf["sampleID"] = sampleID

                SV_CNV = SV_CNV.append(vcf)
                SV_CNV_annot = SV_CNV_annot.append(annotation)

        # drop duplicates from the annotations
        small_vars_annot = small_vars_annot.drop_duplicates()
        SV_CNV_annot = SV_CNV_annot.drop_duplicates()

        # save each of them
        fun.save_df_as_tab(small_vars, small_vars_file)
        fun.save_df_as_tab(small_vars_annot, small_vars_annot_file)
        fun.save_df_as_tab(SV_CNV, SV_CNV_file)
        fun.save_df_as_tab(SV_CNV_annot, SV_CNV_annot_file)

    # load them
    small_vars = fun.get_tab_as_df_or_empty_df(small_vars_file)
    small_vars_annot = fun.get_tab_as_df_or_empty_df(small_vars_annot_file)
    SV_CNV = fun.get_tab_as_df_or_empty_df(SV_CNV_file)
    SV_CNV_annot = fun.get_tab_as_df_or_empty_df(SV_CNV_annot_file)

    #######################################


    ###### KEEP ONLY VARIANTS THAT ARE EITHER OVERLAPPING TARGET REGIONS OR CAN BE MAPPED TO GENES ######

    # get variants overlapping target regions
    vcf_fields = ["ID", "#CHROM", "POS", "INFO_END"]
    if len(small_vars)==0: small_vars = pd.DataFrame(columns=vcf_fields)
    else: small_vars["INFO_END"] = np.nan

    if len(SV_CNV)==0: SV_CNV = pd.DataFrame(columns=vcf_fields)

    all_vcf_df = small_vars[vcf_fields].append(SV_CNV[vcf_fields])
    variant_IDs_overlapping_target_regions = fun.get_varIDs_overlapping_target_regions(all_vcf_df, target_regions, cache_dir)

    # get variants that overlap target genes
    annot_fields = ["#Uploaded_variation", "Gene"]
    if len(small_vars_annot)==0: small_vars_annot = pd.DataFrame(columns=annot_fields)
    if len(SV_CNV_annot)==0: SV_CNV_annot = pd.DataFrame(columns=annot_fields)

    all_var_annot = small_vars_annot[annot_fields].append(SV_CNV_annot[annot_fields])
    if any(all_var_annot["#Uploaded_variation"].apply(lambda x: ";" in x)): raise ValueError("There should be no multiallelics in variant annotation")

    variant_IDs_in_target_genes = set(all_var_annot[all_var_annot.Gene.isin(target_genes)]["#Uploaded_variation"])

    # get the interesting varIDs
    interesting_variantIDs = variant_IDs_in_target_genes.union(variant_IDs_overlapping_target_regions)

    # filter the annotation DFs
    small_vars_annot = small_vars_annot[small_vars_annot["#Uploaded_variation"].isin(interesting_variantIDs)]
    SV_CNV_annot = SV_CNV_annot[SV_CNV_annot["#Uploaded_variation"].isin(interesting_variantIDs)]

    # filter the small vars
    small_vars["IDset"] = small_vars.ID.apply(lambda x: set(x.split(";")))
    small_vars = small_vars[small_vars.IDset.apply(lambda ids: ids.intersection(interesting_variantIDs)).apply(len)>0]

    # filter the CNVvars
    SV_CNV["IDset"] = SV_CNV.ID.apply(lambda x: set(x.split(";")))
    SV_CNV = SV_CNV[SV_CNV.IDset.apply(lambda ids: ids.intersection(interesting_variantIDs)).apply(len)>0]

    ######################################################################################################

    return small_vars, small_vars_annot, SV_CNV, SV_CNV_annot


def get_genome_variation_browser(df_data, target_regions, target_genes, gff_df, filename, cache_dir, reference_genome, threads=4, sample_group_labels=["default_group"], title="SVbrowser", chrName_to_shortName={}, geneID_to_name={}, skip_coverage=False, fraction_y_domain_by_browser=0.8, min_cov=0, max_cov=2, center_cov=1):

    """This function will draw a genome variation browser for each sample in df_data (each row)."""

    # get index integrated
    df_data = df_data.set_index("sampleID", drop=False)

    # get all the variant dfs, already filtering out by gebe and 
    small_vars, small_vars_annot, SV_CNV, SV_CNV_annot = get_integrated_vars_df(df_data, cache_dir, target_regions, target_genes)

    # check that there are some vars to plot
    if len(small_vars)==0 and len(SV_CNV)==0: raise ValueError("there are no variants to represent.")

    # get the dfs for plotting
    

    print(SV_CNV)
    print(SV_CNV_annot)

    adkjhkjdaadh





