# this has a lot of code that is not used


### a way to get tandem duplications in non random simulations

#### ADD THE TANDEM DUPLICATIONS #### 
"""

# deifine the currently generated svtypes
random_svtypes_generated = {"insertions", "deletions", "translocations", "inversions"}

# get the randomly simulated svs into a dict, together with the randomly generated ones.
random_svtype_to_svfile = {svtype : "%s/%s.tab"%(random_sim_dir, svtype) for svtype in random_svtypes_generated}

# get the remaining bed regions
bed_regions_prefix = "%s/bed_regions_rmaining_afterInsInvDelTra"%genome_outdir
all_regions_bed_df =  pd.read_csv(regions_without_SV_bed, sep="\t", header=None, names=["chromosome", "start", "end"]) # these are the remining after removing the already present ones
regions_without_SV_bed, svtype_to_nSVs = get_bed_df_not_overlapping_with_SVs(all_regions_bed_df, random_svtypes_generated, random_svtype_to_svfile, bed_regions_prefix)

# get the tandem duplications for the 
tandemDuplications_file = "%s/tandemDuplications.tab"%random_sim_dir
sizes = list(pd.read_csv("%s/deletions.tab"%random_sim_dir, sep="\t")["Size"])
get_tandemDuplications_file_in_target_regions(regions_without_SV_bed, sizes, tandemDuplications_file)
"""

####################################


def get_tandemDuplications_file_in_target_regions(target_regions_bed, sizes, tanDup_file, max_n_copies=5):

    """This function writes a file that contains tandem duplications (format like RSVsim) of 'sizes' into tanDupl_file. max_n_copies would be the maxium number of tanDups inserted.

    We will try all the sizes down to 50, and if it is not possible to find them we'll drop them """

    # load the regions bed and sort by size
    regions_df = pd.read_csv(target_regions_bed, sep="\t", header=None, names=["chromosome", "start", "end"])
    regions_df["size"] = regions_df.end - regions_df.start
    regions_df = regions_df.sort_values(by="size", ascending=False)

    # define the total number of tanDels
    maximum_tan = len(sizes); n_tan_simulated = 0

    # initialize a dict
    tanDict = {}

    # define the order of the sizes to try
    if min(sizes)<=50: tried_sizes = list(reversed(sorted(sizes)))
    else: tried_sizes = list(reversed(sorted(sizes))) + list(range(min(sizes), 50))

    for size in tried_sizes:

        # break if you have already simulated all the desired tandels
        if n_tan_simulated>=maximum_tan: break

        # find if there is a region in regions_df where it fits
        df = regions_df[regions_df["size"]>=size]

        # try another size if not found
        if len(df)==0: continue

        # get the largest window (upmost) to insert the tandem duplication and calculate the randomly inserted tanDup
        target_window = df.iloc[0]
        start_tanDup = random.randrange(target_window.start, (target_window.end-size-1))
        end_tanDup = start_tanDup + size

        # keep them
        n_tan_simulated += 1
        tanDict[n_tan_simulated] = {"Chr": target_window.chromosome, "Start":start_tanDup, "End":end_tanDup, "Name":"tandemDuplication%i_%s"%(n_tan_simulated, target_window.chromosome), "Duplications":random.randrange(2, max_n_copies), "Size":(end_tanDup-start_tanDup)}

        # get a df with the left and right parts that would be kept after the tanDup is inserted
        maxI = max(regions_df.index) # get the maxium current idx of the SVs
        df_surrounding_regions = pd.DataFrame(
            {maxI+1: {"chromosome":target_window.chromosome, "start":target_window.start, "end":start_tanDup-1000}, # 
            maxI+2: {"chromosome":target_window.chromosome, "start":end_tanDup+1000, "end":target_window.end}}
            ).transpose()

        df_surrounding_regions["size"] = df_surrounding_regions.end - df_surrounding_regions.start

        # add this df and remove the previous one
        regions_df = regions_df.drop(target_window.name, axis=0) # remove the window were the dup was inserted
        regions_df = regions_df.append(df_surrounding_regions, sort=True).sort_values(by="size", ascending=False) # add the newly broken regions 

    # get the tan dict as numbers
    tan_cols = ["Name", "Chr", "Start", "End", "Size", "Duplications"]
    if len(tanDict)==0: tan_df = pd.DataFrame(columns=tan_cols)
    else: tan_df = pd.DataFrame(tanDict).transpose()

    # write
    tan_df[tan_cols].to_csv(tanDup_file, sep="\t", header=True, index=False)


############# distance to the telomere #################


def get_distance_to_telomere_through_breakpoints(genome, df_bedpe, df_gridss_filt, genomeGraph_outfileprefix):

    """This function returns a df where each genomic position is mappend to the real distance to the telomere through all the breakpoints"""

    adlkadjhdlajthisdoesnotwork


    # get the graoh as a genome
    genome_graph, chrom_to_offset, chrom_to_lenSeq, all_positions, chromosome_end_nodes, chromosome_start_nodes = get_genomeGraph_object(genome, df_bedpe, df_gridss_filt)

    # initialize a df where each col is a telomeric position and the rows are the the positions of the genome graph
    telomericNode_to_pos_to_dist = {}

    # get the shortest paths mapping each node to the 
    for I, telomeric_node in enumerate(chromosome_end_nodes.union(chromosome_start_nodes)):
        print("working on telomeric_node %i"%telomeric_node)

        # get the connected positions list
        connected_positions = genome_graph.subcomponent(telomeric_node)
        print("There are %i connected positions"%len(connected_positions))

        # get the shortest paths to all the others
        all_shorest_paths_len = map(len, genome_graph.get_all_shortest_paths(v=telomeric_node, to=connected_positions, mode="OUT"))
        pos_to_shortestPathLen = dict(zip(connected_positions, all_shorest_paths_len))

        telomericNode_to_pos_to_dist[telomeric_node] = pos_to_shortestPathLen

        if I==2: break

    df_telomeric_dist = pd.DataFrame(telomericNode_to_pos_to_dist)

    print(df_telomeric_dist)


    ldjhdljjhlsfd

    return df_telomeric_dist


##### GRAPH GENOME OPERATIONS #####
if type_coverage_to_filterTANDEL=="coverage_rel_to_predFromFeats":

    # only take into account if you want to correct by coverage_rel_to_predFromFeats

    # get a graph of the genome
    genomeGraph_outfileprefix = "%s.genomeGraph_incluingBPs%s"%(raw_bedpe_file, include_breakpoints_in_genomeGraph)
    if include_breakpoints_in_genomeGraph is True:
        df_bedpe_arg = df_bedpe
        df_gridss_filt_arg = df_gridss_filt
    
    else: df_bedpe_arg = df_gridss_filt_arg = None

    genome_graph, df_positions_graph = get_genomeGraph_object(reference_genome, df_bedpe_arg, df_gridss_filt_arg, genomeGraph_outfileprefix, replace=replace_FromGridssRun)

    # get a function that takes the GC content, chromosome and distance to the telomere and returns coverage. This is actually a lambda function
    outdir_coverage_calculation = "%s/coverage_per_regions2kb_incluingBPs%s"%(working_dir, include_breakpoints_in_genomeGraph); make_folder(outdir_coverage_calculation)
    df_coverage_train = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir_coverage_calculation, sorted_bam, windows_file="none", replace=replace_FromGridssRun, window_l=2000), sep="\t")

    distToTel_chrom_GC_to_coverage_fn = get_distanceToTelomere_chromosome_GCcontent_to_coverage_fn(df_coverage_train, reference_genome, genome_graph, df_positions_graph, outdir_coverage_calculation, mitochondrial_chromosome=mitochondrial_chromosome, replace=replace_FromGridssRun)


elif type_coverage_to_filterTANDEL=="mediancov_1": distToTel_chrom_GC_to_coverage_fn = genome_graph = df_positions_graph = None


####### monophyly of vars

    # add to the tree if it is there
            for l in species_tree.get_leaves():

                if svID in ID_to_svIDs[l.name]: l.svID = svID
                else: l.svID = "none"

            # ask if it is monophyletic
            svID_is_monophyletic, type_phyly, nodes_monophyletic = species_tree.check_monophyly(values=[svID], target_attr="svID")

            # length monophyletic nodes
            max_len_monophyletic_nodes = {len(mn) for mn in nodes_monophyletic}



            print(max_len_monophyletic_nodes)


    adlkdalkad


def get_clove_output_with_coverage_forTANDEL(outfile_clove, reference_genome, sorted_bam, distToTel_chrom_GC_to_coverage_fn, genome_graph, df_positions_graph, replace=False, run_in_parallel=False, delete_bams=False):

    """Takes the output of clove and adds the coverage of the TAN and DEL, or -1.

    If you set genome_graph to None it will skip adding the coverage relative to sequence features"""

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
            if genome_graph is not None:

                outdir_rel_cov_calculation = "%s_calculatig_rel_coverage"%outfile_clove
                coverage_df["coverage_rel_to_predFromFeats"] = get_coverage_list_relative_to_predictedFromTelomereAndGCcontent(coverage_df, reference_genome, distToTel_chrom_GC_to_coverage_fn, genome_graph, df_positions_graph, outdir_rel_cov_calculation, real_coverage_field="mediancov_1", replace=replace)

            else: coverage_df["coverage_rel_to_predFromFeats"] = -1

        else: coverage_df = pd.DataFrame(columns=["chromosome", "end", "length", "mediancov_1", "nocoveragebp_1", "percentcovered_1", "start"])

        # merge
        merged_df = df_clove.merge(coverage_df, how="left", left_on=["#CHROM", "POS", "END"], right_on=["chromosome", "start", "end"], validate="many_to_one")

        # change types of fields
        merged_df["POS"] = merged_df.POS.apply(get_int)
        merged_df["END"] = merged_df.END.apply(get_int)
        merged_df["START"] = merged_df.START.apply(get_int)

        return merged_df 

    else: return pd.DataFrame()






