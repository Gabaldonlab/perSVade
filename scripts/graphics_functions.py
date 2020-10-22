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


def get_integrated_vars_df(df, cache_dir, target_regions, target_genes):

    """Takes the df with the data and returns an integrated dataframe with small_vars, small_vars_annot, SV_CNV, SV_CNV_annot. It will discard any variants that don't overlap target_regions and are not related to target_genes"""

    # skip the blanks
    df = df[df.sampleID!="blank_to_skip"]

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



def add_gff_traces_as_scatters(fig, gff_df, chrom_to_Xoffset, initial_ylevel=1, width_rect=0.2, xref="x2", yref="y2", interesting_features="all", annotation_offset=0.1, geneID_to_name={}, fig_location=(2,2)):

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
            hovertext = "<b>%s%s-%s%s</b>"%(direction_text, geneID_to_name[r["upmost_parent"]], r["ID"], direction_text)

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

    # depending on several combinations plot one thing or the other

    # inferred by coverage 
    if IDsvtype=="coverageDUP" and r["INFO_BPS_TYPE"] in {"RealBPs", "wholeChrom"}: 
        x = [r["POS"], r["INFO_END"]]
        color = "cyan"
        symbol = "circle"
        name = "coverageDUP_highConf"

    elif IDsvtype=="coverageDUP" and r["INFO_BPS_TYPE"] in {"InferredBPs", "oneRealBP"}: 
        x = [r["POS"], r["INFO_END"]]
        color = "cyan"
        symbol = "square"
        name = "coverageDUP_lowConf"


    elif IDsvtype=="coverageDEL" and r["INFO_BPS_TYPE"] in {"RealBPs", "wholeChrom"}: 
        x = [r["POS"], r["INFO_END"]]
        color = "magenta"
        symbol = "circle"
        name = "coverageDEL_highConf"


    elif IDsvtype=="coverageDEL" and r["INFO_BPS_TYPE"] in {"InferredBPs", "oneRealBP"}: 
        x = [r["POS"], r["INFO_END"]]
        color = "magenta"
        symbol = "square"
        name = "coverageDEL_lowConf"

    # inferred by GRIDSS+CLOVE
    elif IDsvtype=="TDUP" and r["INFO_SVTYPE"]=="TDUP": 
        x = [r["POS"], r["INFO_END"]]
        color = "blue"
        symbol = "circle"
        name = "tandem duplications"

    elif IDsvtype=="DUP" and r["INFO_SVTYPE"]=="DUP": 
        x = [r["POS"], r["INFO_END"]]
        color = "navy"
        symbol = "circle"
        name = "copy-paste insertions"

    elif IDsvtype=="DEL" and r["INFO_SVTYPE"]=="DEL": 
        x = [r["POS"], r["INFO_END"]]
        color = "red"
        symbol = "circle"
        name = "deletions"

    elif IDsvtype.endswith("like") and r["INFO_SVTYPE"]=="BND": 
        x = [r["POS"]]
        color = "gray"
        symbol = "x"
        name = "unclassified SVs BND"

    elif IDsvtype.endswith("like") and r["INFO_SVTYPE"]=="insertionBND": 
        x = [r["POS"]]
        color = "gray"
        symbol = "cross"
        name = "unclassified SVs insertions"

    elif not IDsvtype.endswith("like") and r["INFO_SVTYPE"]=="BND": 
        x = [r["POS"]]
        color = "black"
        symbol = "x"
        name = "classified SVs BND"


    elif not IDsvtype.endswith("like") and r["INFO_SVTYPE"]=="insertionBND": 
        x = [r["POS"]]
        color = "black"
        symbol = "cross"
        name = "classified SVs insertions"

    else:
        print(r["ID"])
        print(r["INFO_SVTYPE"])
        print(r["INFO_BPS_TYPE"])
        raise ValueError("not considered case")

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

def get_clusters_overlapping_vars_allGenome(df_sample):

    """This function takes a df_plot for one sample and returns the list of IDs that are clustered"""

    # keep
    df_sample = cp.deepcopy(df_sample)

    # add a start and an end
    df_sample[["start", "end"]] = df_sample.x.apply(get_start_and_end_df_sample_r)

    # map each ID to the overlapping IDs
    ID_to_overlappingIDs = {}

    # go through each region
    for qStart, qEnd, qID in df_sample[["start", "end", "ID"]].values:

        # find overlapping SVs
        overlappingIDs = set(df_sample[((df_sample.start>=qStart) & (df_sample.end<=qEnd)) | 
                                       ((df_sample.start<=qStart) & (df_sample.end>=qStart)) |  
                                       ((df_sample.start<=qEnd) & (df_sample.end>qEnd))]["ID"])

        # if there are overlaps add them to the cluster_list with itself
        ID_to_overlappingIDs[qID] =  overlappingIDs


    # initialize clusters list
    list_clusters = get_list_clusters_from_dict(ID_to_overlappingIDs)

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

    # define the interesting vcf_fields_onHover
    if vcf_fields_onHover=="all": 

        # define uninteresting fields
        all_vcf_fields = set(SV_CNV.keys())
        uninteresting_vcf_fields = {f for f in all_vcf_fields if len(get_unique_values_series(SV_CNV[f]))==1 or all(pd.isna(SV_CNV[f]))}
        uninteresting_vcf_fields.update({"INFO", "sampleID", "REF", "IDset"})

        # keep the others
        vcf_fields_onHover = all_vcf_fields.difference(uninteresting_vcf_fields)

    # sort
    vcf_fields_onHover = sorted(vcf_fields_onHover)

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

def get_genome_variation_browser(df_data, samples_colors_df, target_regions, target_genes, gff_df, filename, cache_dir, reference_genome, threads=4, sample_group_labels=["default_group"], title="SVbrowser", chrName_to_shortName={}, geneID_to_name={}, fraction_y_domain_by_gene_browser=0.5, min_cov=0, max_cov=2, center_cov=1, only_affected_genes=False, interesting_features="all", vcf_fields_onHover="all"):

    """This function will draw a genome variation browser for each sample in df_data (each row). Some clarifications:

    -interesting_features can be a set of the interesting features in the gff, or 'all'.
    -interesting_vcf_fields are the fields from the vcf to be included. It can be a set or all. If you want fields from INFO, they should be provided as 'INFO_<field>'"""

    # get index integrated
    df_data = df_data.set_index("sampleID", drop=False)

    # get all the variant dfs, already filtering out
    small_vars, small_vars_annot, SV_CNV, SV_CNV_annot = get_integrated_vars_df(df_data, cache_dir, target_regions, target_genes)

    # check that there are some vars to plot
    if len(small_vars)==0 and len(SV_CNV)==0: raise ValueError("there are no variants to represent.")

    # get the dfs for plotting
    all_chromosomes = sorted(set(SV_CNV["#CHROM"]).union(set(small_vars["#CHROM"])))

    # get the chrom to len
    chrom_to_len = fun.get_chr_to_len(reference_genome)

    # define chromosome offset, which is necessary to keep the xaxis organised
    chrom_to_Xoffset = {}
    current_offset = 0
    for chrom in all_chromosomes:
        chrom_to_Xoffset[chrom] = current_offset
        current_offset += chrom_to_len[chrom] + 15000

    # map each chromosome to a color
    chrom_to_color, palette_chromosome = get_value_to_color(all_chromosomes, palette="tab10", type_color="hex")

    # change the chrName_to_shortName if not provided
    if chrName_to_shortName=={}: chrName_to_shortName = {c:c.split("_")[0] for c in all_chromosomes}

    # keep only df that are not chromosomes
    gff_df = gff_df[gff_df.upmost_parent_feature!="chromosome"]

    # get the gff only with the affected genes
    affected_genes = set(SV_CNV_annot.Gene).union(set(small_vars_annot.Gene))
    fun.print_if_verbose("There are %i affected genes"%len(affected_genes))
    if only_affected_genes is True: 

        gff_df = gff_df[~(gff_df.upmost_parent_feature.isin({"gene", "pseudogene"})) | (gff_df.upmost_parent.isin(affected_genes))]

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
    if plot_small_vars: thishastobedeveloped

    # add the structural vars
    if plot_SV: add_SV_CNV_to_fig(fig, SV_CNV, chrom_to_Xoffset, vcf_fields_onHover, sampleID_to_ylevel, SV_offset, fig_location=(1,2))

    # add the coverage
    #if plot_coverage: thishastobedeveloped
  


    #####################################

    ######## 2,2: gene browser ######

    # initialize the features text
    features_text_dicts = [] # this is a list of dictionaries, each of which has the info to get a go.Scatter text object

    # get the rectangles of the chromosomes
    chromosome_rectanges, chromosome_annotations = get_rectangles_and_annotations_for_chromosomes(chrom_to_Xoffset, chrom_to_len, chrom_to_color, chrName_to_shortName, ylevel=0, width_rect=0.2, annotation_offset=0.1)

    # keep them 
    shapes += chromosome_rectanges

    # get the gff features
    add_gff_traces_as_scatters(fig, gff_df, chrom_to_Xoffset, interesting_features=interesting_features, geneID_to_name=geneID_to_name, fig_location=(2,2))


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

    



