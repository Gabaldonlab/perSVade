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
import plotly.plotly as py
import plotly.figure_factory as ff
import plotly.offline as off_py
import plotly.graph_objs as go
from plotly import tools
from plotly.offline import init_notebook_mode, plot, iplot # download_plotlyjs
import cufflinks as cf
import copy as cp
import random as rd
import ast

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import sv_functions as fun

######################################################
######################################################
######################################################


########### DEFINE VARS ##################


#########################################

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

def get_descriptions_affected_features(r, gff_df, all_ANNOTATION_IDs):

    """Taakes a row of an annotation df and returns all the product descriptons of the proteins"""

    # debug the fact that there is no product in the gff
    if "ANNOTATION_product" not in set(gff_df.keys()): return "no_product_field"

    # return for intergenic regions
    if r["Gene"]=="-": return "no_gene"

    # get the features as a set
    features = set(r["Feature"].split(",")).intersection(all_ANNOTATION_IDs)
    if len(features)==0: raise ValueError("there are no properly intersecting features")    

    # get all the proteins
    all_products = set(gff_df.loc[features, "ANNOTATION_product"])
    if len(all_products)!=1: raise ValueError("there are more than 1 product: %s"%features)

    product = next(iter(all_products))

    return product

def get_integrated_vars_df(df, cache_dir, target_regions, target_genes, gff_df):

    """Takes the df with the data and returns an integrated dataframe with small_vars, small_vars_annot, SV_CNV, SV_CNV_annot. It will discard any variants that don't overlap target_regions and are not related to target_genes"""

    # get the gff
    gff_df = cp.deepcopy(gff_df).set_index("ANNOTATION_ID")

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

    # return if empty
    if len(expected_files)==0: 

        small_vars = pd.DataFrame(columns=["ID"])
        small_vars_annot = pd.DataFrame(columns=["Gene"])
        SV_CNV = pd.DataFrame(columns=["ID"])
        SV_CNV_annot = pd.DataFrame(columns=["Gene"])

        return small_vars, small_vars_annot, SV_CNV, SV_CNV_annot

    # load them
    small_vars = fun.get_tab_as_df_or_empty_df(small_vars_file)
    small_vars_annot = fun.get_tab_as_df_or_empty_df(small_vars_annot_file)
    SV_CNV = fun.get_tab_as_df_or_empty_df(SV_CNV_file)
    SV_CNV_annot = fun.get_tab_as_df_or_empty_df(SV_CNV_annot_file)

    #######################################



    ########## QUALITY CONTROL ##########

    annot_fields = ["#Uploaded_variation", "Gene", "Feature"]

    # make sure that there are no multiallelic records
    if any(small_vars["ALT"].apply(lambda x: "," in x)): raise ValueError("There should not be multiallelic records in the vcf")

    # multiallelics doible check
    all_var_annot = small_vars_annot[annot_fields].append(SV_CNV_annot[annot_fields])
    if any(all_var_annot["#Uploaded_variation"].apply(lambda x: ";" in x)): raise ValueError("There should be no multiallelics in variant annotation")

    # check that all genes are in the gff_df.uploaded variation
    affected_genes = set(SV_CNV_annot.Gene).union(set(small_vars_annot.Gene))
    missing_genes = affected_genes.difference(set(gff_df.upmost_parent)).difference({"-"})
    if len(missing_genes)>0: raise ValueError("%s are not found in the gff"%missing_genes)

    #####################################

    ###### ADD FIELDS ######

    # add the short vars representation to the annotation
    small_vars_annot["consequences_set"] = small_vars_annot.Consequence.apply(lambda x: set(x.split(",")))
    small_vars_annot["string_rep_variant"] = small_vars_annot.apply(get_string_representation_of_var, axis=1)

    # add the gene product description
    print("adding protein description")

    all_ANNOTATION_IDs = set(gff_df.index)
    small_vars_annot["Feature_productDescription"] = small_vars_annot.apply(lambda r: get_descriptions_affected_features(r, gff_df, all_ANNOTATION_IDs), axis=1)

    SV_CNV_annot["Feature_productDescription"] = SV_CNV_annot.apply(lambda r: get_descriptions_affected_features(r, gff_df, all_ANNOTATION_IDs), axis=1)

    ########################

    ###### KEEP ONLY VARIANTS THAT ARE EITHER OVERLAPPING TARGET REGIONS OR CAN BE MAPPED TO GENES ######

    # get variants overlapping target regions
    vcf_fields = ["ID", "#CHROM", "POS", "INFO_END"]
    if len(small_vars)==0: small_vars = pd.DataFrame(columns=vcf_fields)
    else: small_vars["INFO_END"] = np.nan

    if len(SV_CNV)==0: SV_CNV = pd.DataFrame(columns=vcf_fields)

    all_vcf_df = small_vars[vcf_fields].append(SV_CNV[vcf_fields])
    variant_IDs_overlapping_target_regions = fun.get_varIDs_overlapping_target_regions(all_vcf_df, target_regions, cache_dir)

    # define the variants that are in the targets
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



def get_Labels_heatmapObject_rows(Labels_colors_df, opacity=1, numeric_ticks=True):

    """Takes a df with color row names, where the index are Labels and the columns are the depth of the levels. The values are colors. It returns a "data" object (fig_data) that can be appended to a plotly figure object with figure.append_trace(fig_data) """

    # redefine depending on the orientation
    df = Labels_colors_df.transpose()
    if numeric_ticks: df.columns = reversed(list(range(len(df.columns)))) # this is to set everything to numbers

    # get all the colors
    all_colors = fun.get_uniqueVals_df(df)
    
    # map them to numbers
    color_to_colnumber = dict(zip(all_colors, range(0, len(all_colors))))
    colnumber_to_color = {v:k for k,v in color_to_colnumber.items()}

    # define df as numbers
    df_numbers = df.applymap(lambda x: color_to_colnumber[x])
    hover = [[row for Icol, col in enumerate(Labels_colors_df.columns)] for Irow, row in enumerate(Labels_colors_df.index)]

    # make it double
    df_numbers = df_numbers.rename(columns={c : 2*c for c in df_numbers.columns}) 

    # get the figure, taking into account the real values
    fig_data = df_numbers.iplot(kind='heatmap', asFigure=True, colors=None, colorscale=False, colorbar=False, legend=False)["data"][0]

    # redefine the colorscale
    if len(all_colors)>1: 
        fig_data = df_numbers.iplot(kind='heatmap', asFigure=True, colors=None, colorscale=False, colorbar=False, legend=False)["data"][0]
        fig_data["colorscale"] = [[n/max(colnumber_to_color), colnumber_to_color[n]] for n in sorted(colnumber_to_color)]

    else: fig_data = df_numbers.iplot(kind='heatmap', asFigure=True, colors=next(iter(all_colors)), colorscale=False, colorbar=False, legend=False)["data"][0]

    # change the opacity
    fig_data["opacity"] = opacity

    # set the hover
    fig_data["text"] = hover
    fig_data["hoverinfo"] = "text"

    # hide the colorbar
    fig_data["colorbar"] = dict(tickvals=[0,0], borderwidth=0, len=0, showticklabels=False, nticks=0, thickness=0, bgcolor="white")


    return fig_data


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



def add_gff_traces_as_scatters(fig, gff_df, chrom_to_Xoffset, gff_annotation_fields, initial_ylevel=1, width_rect=0.2, xref="x2", yref="y2", interesting_features="all", annotation_offset=0.1, geneID_to_name={}, fig_location=(2,2)):

    """This function adds one trace for each feture in interesting_features. It already modifies fig."""

    # define the interesting genes
    interesting_genes = set(gff_df.upmost_parent)

    # define the interesting features
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

    # debug
    if len(genes_df)==0: return None

    # modify the start and the end to get adjusted to the chromosome of procedence
    genes_df["start"] = genes_df.apply(lambda r: r["start"]+chrom_to_Xoffset[r["chromosome"]], axis=1)
    genes_df["end"] = genes_df.apply(lambda r: r["end"]+chrom_to_Xoffset[r["chromosome"]], axis=1)

    # define the geneID_to_name
    if geneID_to_name=={}: geneID_to_name = {g:g for g in set(genes_df.upmost_parent)}


    ####### GET THE ANNOTATIONS SCATTERS ########

    # initialize a dictionary  that will map the ylevel where to draw the genes vs the start and end of these genes
    yLevel_to_GeneCoordinates = {}

    # go through all items of the gff without overlapping, as indicated by yLevel_to_GeneCoordinates
    
    good_levels_list = []
    for ID, source, feature, start, end, strand, upmost_parent in genes_df[["ID", "source", "feature", "start", "end", "strand", "upmost_parent"]].values:

        ####### GET THE CORRECT LEVEL FOR THIS FEATURE #######

        # define the level in which it has to be represented, from 0 to 100
        for level in range(initial_ylevel, 100):

            # if ot explored before
            if level not in yLevel_to_GeneCoordinates: good_level = level; break

            # if it does not overlap with anything in this level
            if all([end<=l_start or start>l_end for l_start, l_end in yLevel_to_GeneCoordinates[level]]): good_level = level; break

        # keep the cordinates
        yLevel_to_GeneCoordinates.setdefault(good_level, set())
        yLevel_to_GeneCoordinates[good_level].add((start, end))

        # keep the good level
        good_levels_list.append(good_level)


        #######################################################

    # keep the good level
    genes_df["good_level"] = good_levels_list
    fun.print_if_verbose("adding features")


    #  go through each feature
    for feat in set(genes_df.feature):

        # get the df
        df_feat = genes_df[genes_df.feature==feat]
        df_feat.index = list(range(0, len(df_feat)))

        # get the color
        color_feature = feature_to_color[feat]

        # init all the positions
        all_x = []
        all_y = []
        all_hovers = []
        all_symbols = []

        # go through each line
        for Ifeat, r in df_feat.iterrows():

            # define hovertext
            direction = {"+":">", "-":"<", ".":"?"}[r["strand"]]
            direction_text = direction*1

            # define the hover text
            hovertext = "<b>%s%s"%(direction_text, r["ID"])
            for f in gff_annotation_fields: 
                if not pd.isna(r[f]): hovertext += "-%s:%s"%(f, r[f])

            hovertext += "%s</b>"%direction_text


            #hovertext = "<b>%s%s-%s%s</b>"%(direction_text, geneID_to_name[r["upmost_parent"]], r["ID"], direction_text)

            # get the symbol
            symbol = {"+":"triangle-right", "-":"triangle-left", ".":"circle"}[r["strand"]]

            # get the line
            x = list(np.linspace(r["start"], r["end"], 2))
            y = list(np.linspace(r["good_level"], r["good_level"], 2))

            # add
            all_x += (x + [x[-1]+1])
            all_y += (y + [None])
            all_hovers += ([hovertext]*3)
            all_symbols += ([symbol]*3)

        # add line of the feat
        line = go.Scatter(x=all_x, y=all_y, showlegend=True, mode="lines+markers", line=dict(color=color_feature, width=2, dash=None), opacity=0.8, hoveron="points+fills", name=feat, text=all_hovers, marker=dict(symbol=all_symbols)) # dash="dash" # supported markers (cross, circle, x, triangle-up)

        fig.append_trace(line, fig_location[0], fig_location[1])

    #####################################################

def get_offset_each_typeData(plot_small_vars, plot_SV, plot_coverage):

    """Returns the yaxis offset of each type of data"""

    if plot_small_vars and plot_SV and plot_coverage: 
        small_vars_offset = (0.33, 1)
        SV_offset = (-0.33, 0.33)
        coverage_offset = (-1, -0.33)

    elif plot_small_vars and plot_SV and plot_coverage is False: 
        small_vars_offset = (0, 1)
        SV_offset = (-1, 0)
        coverage_offset = None

    elif plot_small_vars is False and plot_SV and plot_coverage: 
        small_vars_offset = None
        SV_offset = (0, 1)
        coverage_offset = (-1, 0)

    elif plot_small_vars and plot_SV is False and plot_coverage: 
        small_vars_offset = (0, 1)
        SV_offset = None
        coverage_offset = (-1, 0)

    elif plot_small_vars and plot_SV is False and plot_coverage is False: 
        small_vars_offset = (-1, 1)
        SV_offset = None
        coverage_offset = None

    elif plot_small_vars is False and plot_SV and plot_coverage is False: 
        small_vars_offset = None
        SV_offset = (-1, 1)
        coverage_offset = None

    elif plot_small_vars is False and plot_SV is False and plot_coverage: 
        small_vars_offset = None
        SV_offset = None
        coverage_offset = (-1, 1)

    else: raise ValueError("there was an error with setting the offsets")

    return small_vars_offset, SV_offset, coverage_offset

def get_unique_values_series(series):

    """Takes a series and returns its unique values"""

    if type(series.iloc[0])==set: return set.union(*series)
    else: return set(series)


def get_plot_data_SV_CNV_r(r, vcf_fields_onHover, chrom_to_Xoffset):

    """Takes a row of the SV_CNV df and returns a series with several data related to the plotting."""

    # get the ID svtype
    IDsvtype = r["ID"].split("|")[0]

    # map each understandable IDsvtype to a color and a name
    IDsvtype_to_color = {"TDUP":"blue", "DEL":"red", "TRA":"black", "INV":"green", "CVT":"purple", "CVD":"navy", "IVT":"olive"}
    IDsvtype_to_name = {"TDUP":"tandem duplication", "DEL":"deletion", "TRA":"balanced translocation", "INV":"inversion", "CVT":"inverted translocation", "CVD":"inverted duplication", "IVT":"inverted interchromosomal duplication"}

    # coverage CNVs
    if IDsvtype=="coverageDUP": 
        x = [r["POS"], r["INFO_END"]]
        color = "cyan"
        symbol = "circle"
        name = "coverageDUP"

    elif IDsvtype=="coverageDEL": 
        x = [r["POS"], r["INFO_END"]]
        color = "magenta"
        symbol = "circle"
        name = "coverageDEL"

    # inferred by GRIDSS+CLOVE events
    elif IDsvtype in {"TDUP", "DEL", "INS", "TRA", "INV", "CVT", "CVD", "IVT"} or IDsvtype.endswith("like"):

        # the variant ID defines the color and the name prefix
        if IDsvtype in {"TDUP", "DEL", "TRA", "INV", "CVT", "CVD", "IVT"}:

            color = IDsvtype_to_color[IDsvtype]
            name_prefix = IDsvtype_to_name[IDsvtype]

        elif IDsvtype=="INS":

            type_ins = r["ID"].split("|")[-1]

            if type_ins=="copyPaste":

                color = "steelblue"
                name_prefix = "copy-paste insertion"

            elif type_ins=="cutPaste":

                color = "gray"
                name_prefix = "cut-paste insertion"

        elif IDsvtype.endswith("like"):

            color = "darkorange"
            name_prefix = "unclassified SV"

        # the actual SVTYPE defines the symbol and the suffix
        if r["INFO_SVTYPE"]=="BND": 

            x = [r["POS"]]
            symbol = "x"
            name_suffix = " BNDs"

        elif r["INFO_SVTYPE"]=="insertionBND": 

            x = [r["POS"]]
            symbol = "cross"
            name_suffix = " insertions"

        elif r["INFO_SVTYPE"] in {"TDUP", "DUP", "DEL"}:

            x = [r["POS"], r["INFO_END"]]
            symbol = "circle"
            name_suffix = ""


        # define the final name
        name = name_prefix+name_suffix


        # understandable vars

    else: raise ValueError("%s could not be parsed"%IDsvtype)

    # reformat x
    x = [chrom_to_Xoffset[r["#CHROM"]] + int(pos) for pos in x]

    # get the hover text
    non_INFO_fields = [f for f in vcf_fields_onHover if not f.startswith("INFO_")]
    INFO_fields = [f for f in vcf_fields_onHover if f.startswith("INFO_")]

    hover_text = "<br>".join(["%s: %s"%(f, r[f]) for f in non_INFO_fields])
    hover_text += "<br><br>INFO:<br>"
    hover_text += "<br>".join(["%s: %s"%(f.split("INFO_")[1], r[f]) for f in INFO_fields if not pd.isna(r[f])])

    return pd.Series({"x":x, "color":color, "marker_symbol":symbol, "type_sv":name, "hover_text":hover_text, "chromosome":r["#CHROM"], "sampleID":r["sampleID"], "ID":r["ID"]})


def get_start_and_end_df_sample_r(xlist):

    """Takes a list and returns the first and last position, or just the first repeated if there is only one"""

    if len(xlist)==1: 
        start = xlist[0]
        end = xlist[0]

    else: 
        start = xlist[0]
        end = xlist[-1]

    return pd.Series({"start":start, "end":end})

def get_clusters_overlapping_vars_allGenome(df_sample, type_cluster="overlapping_any"):

    """This function takes a df_plot for one sample and returns the list of IDs that are clustered"""

    # keep
    df_sample = cp.deepcopy(df_sample)

    # add a start and an end, if not there
    if "start" not in df_sample.keys() and "end" not in df_sample.keys(): df_sample[["start", "end"]] = df_sample.x.apply(get_start_and_end_df_sample_r)

    # map each ID to the overlapping IDs
    ID_to_overlappingIDs = {}

    # go through each region
    for qStart, qEnd, qID in df_sample[["start", "end", "ID"]].values:

        # find overlapping SVs
        if type_cluster=="overlapping_any":

            # get IDs that are somehow overlapping
            overlappingIDs = set(df_sample[((df_sample.start>=qStart) & (df_sample.end<=qEnd)) | 
                                           ((df_sample.start<=qStart) & (df_sample.end>=qStart)) |  
                                           ((df_sample.start<=qEnd) & (df_sample.end>qEnd))]["ID"])


        elif type_cluster=="overlapping_but_not_under_query":


            # gets the IDs that are overlapping, but not included under query
            idx_ov = (((df_sample.start<=qStart) & (df_sample.end>=qStart) & (df_sample.end<=qEnd)) |
                      ((df_sample.start<=qEnd) & (df_sample.end>=qEnd) & (df_sample.start>=qStart)))


            overlappingIDs = set(df_sample[idx_ov]["ID"])


        # if there are overlaps add them to the cluster_list with itself
        ID_to_overlappingIDs[qID] =  overlappingIDs

    # initialize clusters list
    list_clusters = fun.get_list_clusters_from_dict(ID_to_overlappingIDs)

    # check
    all_SVs = set(df_sample.ID)
    all_SVs_in_cluster = set.union(*list_clusters)
    if all_SVs!=all_SVs_in_cluster: raise ValueError("all SVs should be in clusters")
  
    return list_clusters

def get_y_position_variants_df_plot(df_plot, sampleID_to_ylevel, offset):

    """This function returns a list with the same order as df_plot with the y locations of each variant, according to sampleID_to_ylevel and offset. The idea is that variants overlapping on the x should be spaced along the y. This should work for both small and SV vars."""


    # initialize ypositions this will map each variant in df_plot to the y position
    dfPlotID_to_ypositon = pd.Series()
    initial_index = cp.deepcopy(df_plot.index)

    # go through each sampleID
    for sampleID in set(df_plot.sampleID):

        # get the df
        df_s = df_plot[df_plot.sampleID==sampleID]

        # check that all IDs are unique
        if len(set(df_s["ID"]))!=len(df_s): raise ValueError("IDs should be unique")

        # get the clusters of IDs 
        list_clusters = get_clusters_overlapping_vars_allGenome(df_s)

        # define the ylevel of the cluster
        min_ylevel = sampleID_to_ylevel[sampleID] + offset[0]
        max_ylevel = sampleID_to_ylevel[sampleID] + offset[1]
        ylevel_range = max_ylevel - min_ylevel

        ######## map each uniqueID to the ylevel to which it corresponds, according to the cluster ########

        # get the ylevels of clusters, and also the width of the rects
        ID_to_ylevel = {}
        for cluster in list_clusters: 

            # define the ylevel range of the cluster
            cluster_ylevel_range = ylevel_range/len(cluster)

            # go through each ID in the cluster and add the ylevel
            for icluster, ID in enumerate(cluster): 

                ID_to_ylevel[ID] = min_ylevel + cluster_ylevel_range*icluster + (cluster_ylevel_range/2)

        # keep
        df_s["ylevel_position"] = df_s.ID.apply(lambda x: ID_to_ylevel[x])
        dfPlotID_to_ypositon = dfPlotID_to_ypositon.append(df_s["ylevel_position"])

        #############################################################################################

    return dfPlotID_to_ypositon[initial_index]

def add_SV_CNV_to_fig(fig, SV_CNV, chrom_to_Xoffset, vcf_fields_onHover, sampleID_to_ylevel, SV_offset, fig_location=(1,2)):

    """This figure adss a scatter for each type of SV or CNV. The hover will include all the vcf files from interesting_vcf_fields."""

    # reindex
    SV_CNV = cp.deepcopy(SV_CNV)
    SV_CNV.index = list(range(0, len(SV_CNV)))

    # define all fields
    all_vcf_fields = set(SV_CNV.keys())

    # define the interesting vcf_fields_onHover
    if vcf_fields_onHover=="all": 

        # define uninteresting fields
        uninteresting_vcf_fields = {f for f in all_vcf_fields if len(get_unique_values_series(SV_CNV[f]))==1 or all(pd.isna(SV_CNV[f]))}
        uninteresting_vcf_fields.update({"INFO", "sampleID", "REF", "IDset"})

        # keep the others
        vcf_fields_onHover = all_vcf_fields.difference(uninteresting_vcf_fields)

    # sort
    vcf_fields_onHover = sorted(vcf_fields_onHover.intersection(all_vcf_fields))

    # get the plotting data
    df_plot = SV_CNV.apply(lambda r: get_plot_data_SV_CNV_r(r, vcf_fields_onHover, chrom_to_Xoffset), axis=1)

    # add the y positions NEEDS TO BE REFACTORED TO AVOID OVERLAPS
    df_plot["y"] = get_y_position_variants_df_plot(df_plot, sampleID_to_ylevel, SV_offset)

    # add one trace for each name
    for type_sv in sorted(set(df_plot.type_sv)):

        # get df
        df_sv = df_plot[df_plot.type_sv==type_sv]

        # get the color and symbol
        color = df_sv.color.iloc[0]
        marker_symbol = df_sv.marker_symbol.iloc[0]

        # init all the positions
        all_x = []
        all_y = []
        all_hovers = []

        # go through each line
        for I, r in df_sv.iterrows():

            # define the number of points
            npoints = len(r["x"])

            # add 
            all_x += (r["x"] + [r["x"][-1]+1])
            all_y += ([r["y"]]*npoints + [None])
            all_hovers += ([r["hover_text"]]*(npoints+1))

        # add line of the feat
        line = go.Scatter(x=all_x, y=all_y, showlegend=True, mode="lines+markers", line=dict(color=color, width=2, dash=None), opacity=1.0, hoveron="points+fills", name=type_sv, text=all_hovers, marker=dict(symbol=marker_symbol)) # dash="dash" # supported markers (cross, circle, x, triangle-up)

        fig.append_trace(line, fig_location[0], fig_location[1])

def get_plot_data_small_vars_r(r, vcf_fields_onHover, chrom_to_Xoffset):


    """Gets a row of the small vars df and returns the plot data"""

    # define the graphics based on the type of variant
    if r["is_protein_altering"] is True and r["overlaps_repeats"] is True:

        color = "teal"
        symbol = "square"
        name = "smallVars - protein altering - overlapping repeats"

    elif r["is_protein_altering"] is True and r["overlaps_repeats"] is False:

        color = "dodgerblue"
        symbol = "square"
        name = "smallVars - protein altering"

    elif r["is_protein_altering"] is False and r["overlaps_repeats"] is True:

        color = "gray"
        symbol = "square"
        name = "smallVars - non protein altering - overlapping repeats"

    elif r["is_protein_altering"] is False and r["overlaps_repeats"] is False:

        color = "black"
        symbol = "square"
        name = "smallVars - non protein altering"

    else: raise ValueError("wrong is_protein_altering or overlaps_repeats")

    # reformat x to ad
    x = [chrom_to_Xoffset[r["#CHROM"]] + int(r["POS"])]

    # get the hover text from the vcf fields
    non_INFO_fields = [f for f in vcf_fields_onHover if not f.startswith("INFO_")]
    INFO_fields = [f for f in vcf_fields_onHover if f.startswith("INFO_")]

    hover_text = "<br>".join(["%s: %s"%(f, r[f]) for f in non_INFO_fields])
    hover_text += "<br><br>INFO:<br>"
    hover_text += "<br>".join(["%s: %s"%(f.split("INFO_")[1], r[f]) for f in INFO_fields if not pd.isna(r[f])])

    # add the variant annotation
    hover_text += "<br>variant effects:<br>%s"%(r["var_annotation_str"])

    return pd.Series({"x":x, "color":color, "marker_symbol":symbol, "type_var":name, "hover_text":hover_text, "chromosome":r["#CHROM"], "sampleID":r["sampleID"], "ID":r["ID"]})



def get_annotation_info_small_vars_r(r, small_vars_annot):

    """Returns several information line of the annotation df into small_vars df"""

    # get the annotation for this var
    df_annot = small_vars_annot.loc[{r["ID"]}]

    # define the easy fields
    is_protein_altering = any(df_annot.is_protein_altering)
    overlaps_repeats = any(df_annot.overlaps_repeats)

    # get the var_annotation_str, which should include all the annotations for the small vars
    var_annotation_str = "<br>".join(df_annot.var_annotation_str)
    
    return pd.Series({"is_protein_altering":is_protein_altering, "var_annotation_str":var_annotation_str, "overlaps_repeats":overlaps_repeats})


def get_string_representation_of_var(row):
    
    """Takes a row of a vars df such as the one provided in get_genesToSamplesToPrivateMutations_df. It returns a string with the meaning of the var."""
    
    # define the most important consequence according to how you ranked the consequences in var_to_IDX
    cons = sorted(row["consequences_set"], key=lambda y: fun.var_to_IDX[y])[0]
    
    # get the meaning
    return "%s|%s.%s|%s"%(fun.consequence_to_abbreviation[cons], fun.protVar_to_info[cons][0], row[fun.protVar_to_info[cons][1]], row[fun.protVar_to_info[cons][2]]) 
     

def add_smallVars_to_fig(fig, small_vars, small_vars_annot, chrom_to_Xoffset, sampleID_to_ylevel, small_vars_offset, threads, replace, vcf_fields_onHover, fig_location=(1,2)):

    """Adds the small-variants data to the figure"""

    # reindex
    small_vars = cp.deepcopy(small_vars)
    small_vars.index = list(range(0, len(small_vars)))

    # get the annot
    small_vars_annot = cp.deepcopy(small_vars_annot).set_index("#Uploaded_variation", drop=False)

    # define all vars
    all_vcf_fields = set(small_vars.keys())

    # define the interesting vcf_fields_onHover
    if vcf_fields_onHover=="all": 

        # define uninteresting fields
        uninteresting_vcf_fields = {f for f in all_vcf_fields if len(get_unique_values_series(small_vars[f]))==1 or all(pd.isna(small_vars[f]))}
        uninteresting_vcf_fields.update({"INFO", "sampleID", "REF", "IDset"})

        # keep the others
        vcf_fields_onHover = all_vcf_fields.difference(uninteresting_vcf_fields)

    # sort fields
    vcf_fields_onHover = sorted(vcf_fields_onHover.intersection(all_vcf_fields))

    # add the var_annotation_str
    small_vars_annot["var_annotation_str"] = small_vars_annot.Gene.apply(str) + ": " + small_vars_annot.string_rep_variant + ": " +  small_vars_annot.Feature_productDescription

    # add some interesting fields that are related to collapsing the small_vars_annot
    small_vars[["is_protein_altering", "var_annotation_str", "overlaps_repeats"]] = small_vars.apply(get_annotation_info_small_vars_r, small_vars_annot=small_vars_annot, axis=1)

    # get plotting data, with the x
    df_plot = small_vars.apply(get_plot_data_small_vars_r, vcf_fields_onHover=vcf_fields_onHover, chrom_to_Xoffset=chrom_to_Xoffset, axis=1)

    # add the y positions NEEDS TO BE REFACTORED TO AVOID OVERLAPS
    df_plot["y"] = get_y_position_variants_df_plot(df_plot, sampleID_to_ylevel, small_vars_offset)

    # add one trace for each name
    for type_var in sorted(set(df_plot.type_var)):

        # get df
        df_var = df_plot[df_plot.type_var==type_var]

        # get the color and symbol
        color = df_var.color.iloc[0]
        marker_symbol = df_var.marker_symbol.iloc[0]

        # init all the positions
        all_x = []
        all_y = []
        all_hovers = []

        # go through each line
        for I, r in df_var.iterrows():

            # add 
            all_x += (r["x"]+[None])
            all_y += ([r["y"]]+[None])
            all_hovers += ([r["hover_text"]]+[None])

        # add line of the feat
        line = go.Scatter(x=all_x, y=all_y, showlegend=True, mode="markers", line=dict(color=color, width=2, dash=None), opacity=1.0, hoveron="points+fills", name=type_var, text=all_hovers, marker=dict(symbol=marker_symbol)) # dash="dash" # supported markers (cross, circle, x, triangle-up)

        fig.append_trace(line, fig_location[0], fig_location[1])



def get_closest_3region_bendPOS(r, SV_CNV):

    """Takes a row of a SV_CNV df and returns the location of the nearest BP 3'"""


    # get a df with the options
    df_options = SV_CNV[(SV_CNV["#CHROM"]==r["#CHROM"]) & (SV_CNV["POS"]>r["POS"])]

    if len(df_options)==0: return r["POS"] + 1
    else: return fun.find_nearest(df_options.POS, r["POS"])

def get_closest_5region_bendPOS(r, SV_CNV):

    """Takes a row of a SV_CNV df and returns the location of the nearest BP 5'"""


    # get a df with the options
    df_options = SV_CNV[(SV_CNV["#CHROM"]==r["#CHROM"]) & (SV_CNV["POS"]<r["POS"])]

    if len(df_options)==0: return r["POS"]
    else: return fun.find_nearest(df_options.POS, r["POS"])

def get_list_clusters_overlapping_df_windows(df_windows):

    """Takes a df_windows with chromosome, start, and end, where the ID is the index. It returns a list of sets, each set containing IDs of CNVs that overlap by >=80% and are from the same type """

    # init dict
    ID_to_overlappingIDs = {}

    # define the fields for comparison
    equal_fields = ["chromosome"]
    approximate_fields = ["start", "end"]
    chromField_to_posFields = {"chromosome":{"start": "start", "end": "end"}}

    # go through each variant
    for ID, r in df_windows.iterrows(): 

        # define the overlapping df
        df_overlapping = df_windows[(df_windows.apply(lambda r_cnv: fun.get_is_matching_predicted_and_known_rows(r, r_cnv, equal_fields, approximate_fields, chromField_to_posFields, tol_bp=1000000000000, pct_overlap=0.8), axis=1))]

        # keep IDs
        ID_to_overlappingIDs[ID] = set(df_overlapping.index)

    # initialize clusters list
    list_clusters = fun.get_list_clusters_from_dict(ID_to_overlappingIDs)

    # check
    all_IDs = set(df_windows.index)
    all_IDs_in_cluster = set.union(*list_clusters)
    if all_IDs!=all_IDs_in_cluster: raise ValueError("all IDs should be in clusters")
  
    return list_clusters


def get_df_windows_draw_coverage(SV_CNV, chrom_to_Xoffset, reference_genome, threads):

    """This function selects a set of regions for which we should draw coverage.  Those arrownd breakpoints, covering the START-END variants or whole chromosomes.

    This function already merges the clustered regions and splits each of them into slices of 20. The final df can be measured for coverage."""

    ############# GET WINDOWS COVERAGE #############

    # define objs
    SV_CNV = cp.deepcopy(SV_CNV)
    windows_fields = ["chromosome", "start", "end"]

    # get chrom len
    chrom_to_len = fun.get_chr_to_len(reference_genome)

    # initialize with the chromosomal regions (0-based)
    df_windows_chroms = pd.DataFrame({chrom : {"chromosome":chrom, "start":0, "end":chrom_to_len[chrom]} for chrom in chrom_to_Xoffset}).transpose()

    if len(SV_CNV)>0:

        # get a df winows of the SVs (0-based)
        SV_CNV_END = SV_CNV[~pd.isna(SV_CNV.INFO_END)]
        SV_CNV_noEND = SV_CNV[pd.isna(SV_CNV.INFO_END)]
        SV_CNV_noEND["INFO_END"] = SV_CNV_noEND.apply(lambda r: get_closest_3region_bendPOS(r, SV_CNV), axis=1)
        SV_CNV_noEND["POS"] = SV_CNV_noEND.apply(lambda r: get_closest_5region_bendPOS(r, SV_CNV), axis=1)
        SV_CNV_all = SV_CNV_END.append(SV_CNV_noEND)

        df_windows_SVs = SV_CNV_all[["#CHROM", "POS", "INFO_END"]].rename(columns={"#CHROM":"chromosome", "POS":"start", "INFO_END":"end"})
        df_windows_SVs["start"] -= 1

        # set to ints
        for f in ["start", "end"]: df_windows_SVs[f] = df_windows_SVs[f].apply(int)

        # map each chromosome to the positions that are breakpoints
        chrom_to_bpPositions = dict(df_windows_SVs.groupby("chromosome").apply(lambda df_c: set(df_c.start).union(df_c.end)))

        # map each chromosome to the max pos
        chrom_to_maxPos = {chrom : length-1 for chrom, length in chrom_to_len.items()}

        # init dict
        final_df_windows_SVs = pd.DataFrame()

        # get a df_windows of the df_windows_SVs +- 5' and 3' regions
        for region in ["5", "3"]: 

            # get a df with the regions
            df_region = df_windows_SVs.apply(lambda r: fun.get_target_region_row(r, region, chrom_to_bpPositions[r["chromosome"]], chrom_to_maxPos[r["chromosome"]], max_relative_len_neighbors=1), axis=1)

            final_df_windows_SVs = final_df_windows_SVs.append(df_region[windows_fields])


        # add the target
        final_df_windows_SVs = final_df_windows_SVs.append(df_windows_SVs[windows_fields])
        if len(final_df_windows_SVs)!=(len(df_windows_SVs)*3): raise ValueError("something went wrong with the neighbors calculus ")

        # keep all
        df_windows = df_windows_chroms[windows_fields].append(final_df_windows_SVs[windows_fields]).drop_duplicates()

    else: df_windows = df_windows_chroms[windows_fields]

    # get index
    df_windows.index = list(range(0, len(df_windows)))
    df_windows["ID"] = df_windows.index


    #######################################################

    ######## MERGE REDUNDANT WINDOWS ########

    # skip this
    #df_windows_NR = df_windows

    # init
    df_windows_NR = pd.DataFrame()

    # go through each chromosome adding windows
    for chrom in chrom_to_Xoffset:

        # get df 
        df_c = df_windows[df_windows.chromosome==chrom]

        # get the clusters
        list_clusters = get_list_clusters_overlapping_df_windows(df_c)

        # get the new window_c, with the chromosome metrged
        new_df_c = pd.DataFrame({Ic : {"chromosome":chrom, "start":min(df_c.loc[cluster, "start"]), "end":max(df_c.loc[cluster, "end"])} for Ic, cluster in enumerate(list_clusters)}).transpose()

        # keep
        df_windows_NR = df_windows_NR.append(new_df_c[windows_fields])


    df_windows_NR.index = list(range(0, len(df_windows_NR)))
    df_windows_NR["ID"] = df_windows_NR.index

    #########################################

    # get the subwindows
    df_windows_final = fun.get_df_subwindows_from_df_windows(df_windows_NR, n_subwindows=10)
    df_windows_final = df_windows_final.drop_duplicates(subset=["chromosome", "start", "end"])

    return df_windows_final

def get_relative_coverage_to_bg(r):

    """Takes a row of the df_s and returns the coverage relative to the bg"""

    # add pseudocounds
    rel_cov = r["relative_coverage"] + 0.01
    rel_cov_bg = r["relative_coverage_bg"] + 0.01

    return rel_cov/rel_cov_bg

def get_df_coverage_with_coverage_relative_to_bg(df_coverage, sampleID_to_backgroundSamples):

    """Takes a df coverage and returns it with an added col: relative_coverage_to_bg"""

    # define the samples
    all_samples = set(df_coverage.sampleID)

    # get as index
    df_coverage["window_ID"] = df_coverage.chromosome + "_" + df_coverage.start.apply(str) + "_" + df_coverage.end.apply(str)
    df_coverage = df_coverage.set_index("window_ID")

    # init the df_coverage final
    df_coverage_final = pd.DataFrame()

    # add the coverage relative to the 'compare samples'
    for sampleID in all_samples:

        # get the df
        df_s = df_coverage[df_coverage.sampleID==sampleID]

        # check that the idx is unique
        if len(df_s)!=len(set(df_s.index)): raise ValueError("idx should be unique")

        # get the other samples
        if sampleID in sampleID_to_backgroundSamples: bg_samples = sampleID_to_backgroundSamples[sampleID]
        else: bg_samples = all_samples.difference({sampleID})

        # if empty, set them all
        if len(all_samples)==1: bg_samples = all_samples

        # get a df with the relative coverage of each of the regions
        df_bg = df_s[["chromosome"]]
        for bg_sample in bg_samples: df_bg[bg_sample] = df_coverage[df_coverage.sampleID==bg_sample].loc[df_bg.index]["relative_coverage"]
        df_bg.pop("chromosome")

        # check
        if list(df_bg.index)!=list(df_s.index): raise ValueError("error in df_s preparation")

        # keep the coverage relative
        df_s["relative_coverage_bg"] = df_bg.apply(np.median, axis=1)
        df_s["relative_coverage_to_bg"] = df_s.apply(get_relative_coverage_to_bg, axis=1)

        # check 
        if any(pd.isna(df_s.relative_coverage_to_bg)): raise ValueError("there are nans in coverage")

        # keep
        df_coverage_final = df_coverage_final.append(df_s)

    return df_coverage_final


def get_coverage_y_plot_for_df_coverage_r(r, rel_cov_f, sampleID_to_ylevel, coverage_offset, min_cov, max_cov):

    """This function gets the actual position in the plot of a coverage window r"""

    # get the relative coverage
    rel_cov = r[rel_cov_f]

    # get the min_y and max_y
    min_y = sampleID_to_ylevel[r["sampleID"]] + coverage_offset[0]
    max_y = sampleID_to_ylevel[r["sampleID"]] + coverage_offset[1]

    # get the ypos in the relative space
    relative_value = min([abs(rel_cov-min_cov)/abs(max_cov-min_cov), 1])
    ypos = min_y + relative_value*(abs(max_y-min_y))

    return ypos

def add_coverage_to_fig(df_data, fig, SV_CNV, chrom_to_Xoffset, sampleID_to_ylevel, coverage_offset, threads, replace, reference_genome, cache_dir, mitochondrial_chromosome, sampleID_to_backgroundSamples, min_cov, max_cov, center_cov, fig_location=(1,2)):

    """This function adds a coverage trace for each sampleID into fig."""

    ######### GET COVERAGE FOR EACH WINDOW IN EACH SAMPLE #########
    
    # define the file
    coverage_all_file = "%s/coverage_df_%s.tab"%(cache_dir, "-".join(df_data.sampleID))

    # get the chrom to len
    chr_to_len = fun.get_chr_to_len(reference_genome)

    #if fun.file_is_empty(coverage_all_file) or replace is True:

    if True:

        fun.print_if_verbose("getting coverage calculated")

        # define the regions for which to draw coverage (these are the non-redundant, fragmented regions)
        df_windows = get_df_windows_draw_coverage(SV_CNV, chrom_to_Xoffset, reference_genome, threads) 

        # write them as a bed
        windows_bed = "%s/subwindows.bed"%cache_dir
        df_windows[["chromosome", "start", "end"]].to_csv(windows_bed, sep="\t", header=True, index=False)

        # check that all the start is before end
        if any(df_windows.end<=df_windows.start): 
            print(df_windows)
            print(df_windows[df_windows.end<=df_windows.start])
            raise ValueError("there should not be any start, end")

        # init coverage
        df_coverage = pd.DataFrame()

        # go through each sample
        for sampleID, r_data in df_data.iterrows():

            # get the coverage-per-window
            df_s = fun.get_coverage_per_window_df_without_repeating(reference_genome, r_data["sorted_bam"], windows_bed, replace=replace, run_in_parallel=True, delete_bams=True, threads=threads)

            # add the type of data
            df_s["chromosome_length"] = df_s.chromosome.apply(lambda c: chr_to_len[c])
            df_s["fraction_chrom_covered"] = df_s.length / df_s.chromosome_length
            fractionCovBool_to_typeRegion = {True:"chromosome_fraction", False:"region"}
            df_s["type_region"] = (df_s.fraction_chrom_covered>=0.04).apply(lambda f: fractionCovBool_to_typeRegion[f])

            # add the relative coverage of the sample
            median_coverage = fun.get_median_coverage(df_s[df_s.type_region=="chromosome_fraction"], mitochondrial_chromosome)
            df_s["relative_coverage"] = df_s.mediancov_1 / median_coverage

            # add the sample ID
            df_s["sampleID"] = sampleID

            # add
            df_coverage = df_coverage.append(df_s)

        # write
        fun.save_df_as_tab(df_coverage, coverage_all_file)

    # load 
    df_coverage = fun.get_tab_as_df_or_empty_df(coverage_all_file)

    ####################################################################

    # get the coverage relative to the background samples
    df_coverage = get_df_coverage_with_coverage_relative_to_bg(df_coverage, sampleID_to_backgroundSamples)


    ######### GET THE FIGURE DATA #########

    # add the medium position, where the points will be drawn
    df_coverage["medium_pos"] = (df_coverage.start + (df_coverage.length/2)).apply(int)

    # add the plot x position, depending on chrom_to_Xoffset
    df_coverage["plot_xpos"] = df_coverage.chromosome.apply(lambda x: chrom_to_Xoffset[x]) + df_coverage.medium_pos 

    # go through each type of data to plot
    for field, color in [("relative_coverage", "blue"), ("relative_coverage_to_bg", "red")]: 

        # add the real ypos
        df_coverage["plot_ypos"] = df_coverage.apply(lambda r: get_coverage_y_plot_for_df_coverage_r(r, field, sampleID_to_ylevel, coverage_offset, min_cov, max_cov), axis=1)

        # init data
        all_x = []
        all_y = []
        all_hovers = []

        # go through each sample and chromosome
        for sampleID in set(df_coverage.sampleID):
            for chrom in set(df_coverage.chromosome):

                # get data
                df_region = df_coverage[(df_coverage.sampleID==sampleID) & (df_coverage.chromosome==chrom)].sort_values(by=["plot_xpos"])

                # add to lists with nones
                all_x += list(df_region.plot_xpos) + [None]
                all_y += list(df_region.plot_ypos) + [None]
                all_hovers += list(df_region[field].apply(lambda x: "%.2f"%x)) + [None]

        # add line of the feat
        #line = go.Scatter(x=x, y=y, showlegend=True, mode="markers", opacity=1.0, hoveron="points+fills", name=field, text=hover, marker=dict(symbol="circle", color=color, size=4)) 

        line = go.Scatter(x=all_x, y=all_y, showlegend=True, mode="lines+markers", opacity=1.0, hoveron="points+fills", name=field, text=all_hovers, line=dict(color=color, width=1), marker=dict(symbol="circle", color=color, size=3), visible="legendonly") 

        fig.append_trace(line, fig_location[0], fig_location[1])


    ####################################

def get_genome_variation_browser(df_data, samples_colors_df, target_regions, target_genes, gff_df, filename, cache_dir, reference_genome, threads=4, sample_group_labels=["default_group"], title="SVbrowser", chrName_to_shortName={}, geneID_to_name={}, fraction_y_domain_by_gene_browser=0.5, min_cov=0, max_cov=2, center_cov=1, only_affected_genes=False, interesting_features="all", vcf_fields_onHover="all", replace=False, mitochondrial_chromosome="mito_C_glabrata_CBS138", sampleID_to_backgroundSamples={}, gff_annotation_fields=set()):

    """This function will draw a genome variation browser for each sample in df_data (each row). Some clarifications:

    -interesting_features can be a set of the interesting features in the gff, or 'all'.
    -vcf_fields_onHover are the fields from the vcf to be included. It can be a set or all. If you want fields from INFO, they should be provided as 'INFO_<field>'.
    -sampleID_to_backgroundSamples is a dict that maps each sample ID to the samples the coverage should be compared against. For those samples in which it is empty, all the other samples will be considered."""

    # get index integrated
    df_data = df_data.set_index("sampleID", drop=False)

    # get only the gff features where the parentID is in gene and pseudogene
    #gff_df = gff_df[gff_df.upmost_parent_feature.isin({"gene", "pseudogene"})]

    # get all the variant dfs, already filtering out
    small_vars, small_vars_annot, SV_CNV, SV_CNV_annot = get_integrated_vars_df(df_data, cache_dir, target_regions, target_genes, gff_df)

    # check that there are some vars to plot
    #if len(small_vars)==0 and len(SV_CNV)==0: raise ValueError("there are no variants to represent.")

  
    # get the chrom to len
    chrom_to_len = fun.get_chr_to_len(reference_genome)

    # get the chroms for plotting
    all_chromosomes = sorted(set(chrom_to_len))

    # define chromosome offset, which is necessary to keep the xaxis organised
    chrom_to_Xoffset = {}
    current_offset = 0
    for chrom in all_chromosomes:
        chrom_to_Xoffset[chrom] = current_offset
        current_offset += chrom_to_len[chrom] + 15000

    # map each chromosome to a color
    chrom_to_color, palette_chromosome = get_value_to_color(all_chromosomes, palette="tab10", type_color="hex")

    # change the chrName_to_shortName if not provided
    #if chrName_to_shortName=={}: chrName_to_shortName = {c:c.split("_")[0] for c in all_chromosomes}
    if chrName_to_shortName=={}: chrName_to_shortName = {c:c for c in all_chromosomes}

    # keep only df that are not chromosomes
    gff_df = gff_df[gff_df.upmost_parent_feature!="chromosome"]

    # get the gff only with the affected genes
    affected_genes = set(SV_CNV_annot.Gene).union(set(small_vars_annot.Gene))

    # get the gff of these genes
    fun.print_if_verbose("There are %i affected genes"%len(affected_genes))
    if only_affected_genes is True: 

        #gff_df = gff_df[~(gff_df.upmost_parent_feature.isin({"gene", "pseudogene"})) | (gff_df.upmost_parent.isin(affected_genes))]

        gff_df = gff_df[(gff_df.upmost_parent.isin(affected_genes))]

    # setup the layout (2 rows and 2 cols)
    fig = tools.make_subplots(rows=2, cols=2, specs=[[{}, {}], [{}, {}]], vertical_spacing=0.0, horizontal_spacing=0.0, subplot_titles=("", "", "", ""), shared_yaxes=True, shared_xaxes=True, print_grid=True)

    # init the shapes
    shapes = []
    
    ######## 1,1: sample levels ######

    # get the heatmap
    fig.append_trace(get_Labels_heatmapObject_rows(samples_colors_df), 1, 1)

    # map each ylevel to the corresponding label
    ylevels = [x*2 for x in list(reversed(range(len(samples_colors_df))))] # it is always 2 scale
    ylabels = samples_colors_df.index
    sampleID_to_ylevel = dict(zip(ylabels, ylevels))

    ##################################

    ######## 1,2: variants browser ######

    # get the y offset of each type of vars. This depends on what is provided
    plot_small_vars = len(small_vars)>0
    plot_SV = len(SV_CNV)>0
    plot_coverage = "sorted_bam" in df_data.keys()

    small_vars_offset, SV_offset, coverage_offset = get_offset_each_typeData(plot_small_vars, plot_SV, plot_coverage)

    # add the background, related by the last column of samples_colors_df
    lastCol = samples_colors_df.columns[-1]
    lastChr = sorted(chrom_to_Xoffset.items(), key=(lambda x: x[1]))[-1][0]
    end_lastChr = chrom_to_Xoffset[lastChr] + chrom_to_Xoffset[lastChr]
    background_df = pd.DataFrame({x : {s : samples_colors_df.loc[{s}, lastCol].iloc[0] for s in samples_colors_df.index} for x in [0, end_lastChr] }).loc[samples_colors_df.index]
    fig.append_trace(get_Labels_heatmapObject_rows(background_df, opacity=0.15, numeric_ticks=True), 1, 2)

    # draw a line sepparating each data visualization cathegory
    all_y_grid = set.union(*[set(x) for x in [small_vars_offset, SV_offset, coverage_offset] if x is not None])
    for sampleID, ylevel in sampleID_to_ylevel.items():
        for y in all_y_grid: 

            # get the real y
            real_y = ylevel + y

            line = go.Scatter(x=np.linspace(0, end_lastChr, 3), y=np.linspace(real_y, real_y, 3), showlegend=False, mode="lines", line=dict(color="gray", width=2, dash="dot"), opacity=0.8, hoveron="points+fills", text="") # dash="dash" # supported markers (cross, circle, x, triangle-up)

            fig.append_trace(line, 1, 2)


    # add the small vars
    if plot_small_vars: add_smallVars_to_fig(fig, small_vars, small_vars_annot, chrom_to_Xoffset, sampleID_to_ylevel, small_vars_offset, threads, replace, vcf_fields_onHover, fig_location=(1,2))

    # add the coverage
    if plot_coverage: add_coverage_to_fig(df_data, fig, SV_CNV, chrom_to_Xoffset, sampleID_to_ylevel, coverage_offset, threads, replace, reference_genome, cache_dir, mitochondrial_chromosome, sampleID_to_backgroundSamples, min_cov, max_cov, center_cov, fig_location=(1,2))
  
    # add the structural vars
    if plot_SV: add_SV_CNV_to_fig(fig, SV_CNV, chrom_to_Xoffset, vcf_fields_onHover, sampleID_to_ylevel, SV_offset, fig_location=(1,2))

    #####################################

    ######## 2,2: gene browser ######

    # initialize the features text
    features_text_dicts = [] # this is a list of dictionaries, each of which has the info to get a go.Scatter text object

    # get the rectangles of the chromosomes
    chromosome_rectanges, chromosome_annotations = get_rectangles_and_annotations_for_chromosomes(chrom_to_Xoffset, chrom_to_len, chrom_to_color, chrName_to_shortName, ylevel=0, width_rect=0.2, annotation_offset=0.1)

    # keep them 
    shapes += chromosome_rectanges

    # get the gff features
    add_gff_traces_as_scatters(fig, gff_df, chrom_to_Xoffset, gff_annotation_fields, interesting_features=interesting_features, geneID_to_name=geneID_to_name, fig_location=(2,2))


    # keep the text as trace, one with all
    for label, features_text_dicts in [["chromosomes", chromosome_annotations]]:

        # get the data structures to define the plots
        keys_first_text = features_text_dicts[0].keys()
        dict_all_texts = {key : [text_dict[key] for text_dict in features_text_dicts] for key in keys_first_text}
        keys_first_textfont_dict = dict_all_texts["textfont"][0].keys()
        textfont_dict_all = {key:[textfont_dict[key] for textfont_dict in dict_all_texts["textfont"]] for key in keys_first_textfont_dict}

        trace_texts = go.Scatter(x=dict_all_texts["x"], y=dict_all_texts["y"], text=dict_all_texts["text"], mode="text", showlegend=True, name=label, textfont=textfont_dict_all, textposition="center", visible="legendonly")
        fig.append_trace(trace_texts, 2, 2)
       
    ##################################

    ###### FIGURE LAYOUT ######

    # define domains of axes
    width_samplelevels = len(samples_colors_df.columns)*0.02
    xdomain = [0, width_samplelevels] # this is the domain of the samples x axis
    xdomain2 = [width_samplelevels, 1] # domain browser
    #witdth_browser = max([1 - len(samples_colors_df)*0.1, 0.1]) # this is depending on the mutations
    #witdth_browser = max(yLevel_to_GeneCoordinates)*0.03 # this is depending on the broser
    witdth_browser = fraction_y_domain_by_gene_browser # this is using the given fraction
    ydomain = [witdth_browser, 1] # mutations
    ydomain2 = [0, witdth_browser] # exons

    # add the layout
    yaxis = {"domain":ydomain, "tickmode":"array", "tickvals":ylevels, "ticktext":ylabels, "visible":True, "showgrid":True, "zeroline":False}

    fig['layout'].update(title=title,  xaxis={"domain":xdomain, "visible":True, "showgrid":True, "zeroline":False},  xaxis2={"domain":xdomain2, "visible":True, "showgrid":True, "zeroline":False}, yaxis=yaxis,  yaxis2={"domain":ydomain2, "visible":True, "showgrid":True, "zeroline":False}, font={"size":12}, shapes=shapes, margin=go.Margin(l=250, r=250, b=50, t=50), showlegend=True)

    # get the fig
    #init_notebook_mode()
    config={'editable': False}
    print("Writing %s"%filename)
    off_py.plot(fig, filename=filename, auto_open=False, config=config)

    return fig

    

#######################################################
############### FUNCTIONS PLOTS TESTING ###############
#######################################################

# this section includes the functions of the plots to test


#######################################################
#######################################################
#######################################################






