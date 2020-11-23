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

bbmap_reformat_sh = "%s/bin/reformat.sh"%EnvDir

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





########### CALCULATE GC CONTENT BASED ON PREVIOUSLY MEASURED REGIONS  ###########

def get_df_with_GCcontent(df_windows, genome, gcontent_outfile, replace=False):

    """This function takes a df with windows of the genome and adds the gc content for each window, writing a file under gcontent_outfile. It will only do those that have been already measured"""

    print("Getting GC content")

    if file_is_empty(gcontent_outfile) or replace is True:

        # define the initial index
        initial_index = list(df_windows.index)

        # define a file that has all the GC content per windows of the genome
        all_gcontent_outfile = "%s.GCcontent_all_windows.py"%genome

        # load the previously generated windows
        GCcontent_cols = ["chromosome", "start", "end", "GCcontent"]
        if file_is_empty(all_gcontent_outfile) or replace is True: 
        #if True:
            remove_file(all_gcontent_outfile)
            previous_df_windows = pd.DataFrame(columns=GCcontent_cols)

        else: previous_df_windows = load_object(all_gcontent_outfile)

        # set as index the combination of each chrom, start and end
        previous_df_windows = previous_df_windows.set_index(["chromosome", "start", "end"], drop=False)
        all_previous_windows = set(previous_df_windows.index)

        # resort
        df_windows = df_windows.sort_values(by=["chromosome", "start", "end"]).set_index(["chromosome", "start", "end"], drop=False)

        # get the new windows
        new_windows = list(set(df_windows.index).difference(all_previous_windows))
        df_windows_new = df_windows.loc[new_windows]

        if len(df_windows_new)>0:
            print("adding %i new windows"%len(df_windows_new))

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
            startWindow_to_gc = gc_df[["chromosome", "start_window", "end_window", "is_in_GC"]].groupby(["chromosome", "start_window", "end_window"]).mean()["is_in_GC"]

        else: startWindow_to_gc = pd.Series()

        # add the previously generated windows
        startWindow_to_gc = startWindow_to_gc.append(previous_df_windows["GCcontent"])

        # get into df_windows
        print(df_windows, startWindow_to_gc)
        df_windows["GCcontent"] = list(startWindow_to_gc.loc[df_windows.index])

        # save the ultimate windows file
        all_df_windows = previous_df_windows.append(df_windows[GCcontent_cols])
        all_gcontent_outfile_tmp = "%s.%s"%(all_gcontent_outfile, id_generator(15))
        save_object(all_df_windows, all_gcontent_outfile_tmp)
        os.rename(all_gcontent_outfile_tmp, all_gcontent_outfile)

        # at the end save the df windows
        df_windows.index = initial_index
        save_object(df_windows, gcontent_outfile)

    else: df_windows = load_object(gcontent_outfile)

    return df_windows




###### get compatible translocations ####


def get_compatible_translocations_df(df, chr_to_len):

    """Takes a df with translocations and returns those that are compatible with each other, so that two arms of a chromosome not have a translocation"""


    # change the ID
    df = df.set_index("ID", drop=False)

    # define the vars that are real or simulated
    real_df = df[df.ID.apply(lambda x: x.endswith("_realSV"))]
    sim_df = df[~df.ID.apply(lambda x: x.endswith("_realSV"))]

    # define all the chromosomes
    all_chromosomes = set(df["ChrA"]).union(set(df["ChrB"]))

    # initalize a bed with the interesting bed regions
    df_bed_allRegions, nSVs = get_affected_region_bed_for_SVdf(real_df, "translocations", all_chromosomes) 

    # initialize the interesting IDs
    interesting_IDs = list(real_df.ID)

    # go through each sv simulated and keep it if it does not overlap any of the regions in df_bed_allRegions
    for ID, sv_series in sim_df.iterrows():

        # define series as df
        sv_series_df = pd.DataFrame({0 : sv_series}).transpose()

        # get the affected regions
        sv_bed, nSVs = get_affected_region_bed_for_SVdf(sv_series_df, "translocations", all_chromosomes)

        # get if they are overlapinmg
        any_regions_overlapping = any(df_bed_allRegions.apply(lambda rprevious: any(sv_bed.apply(lambda rnew: get_is_overlapping_query_vs_target_region(rprevious, rnew), axis=1)), axis=1))

        if not any_regions_overlapping: interesting_IDs.append(ID)

    # get the final df
    df = df.loc[interesting_IDs]
    #df = real_df
    #df = sim_df


    print(df[svtype_to_fieldsDict["translocations"]["all_fields"]])

    # change the positions so that they do not exceed chr boundaries
    posF_to_chrF = svtype_to_fieldsDict["translocations"]["positionField_to_chromosome"]

    for f in svtype_to_fieldsDict["translocations"]["position_fields"]: df[f] = df.apply(lambda r: set_position_to_max(r[f], chr_to_len[r[posF_to_chrF[f]]]), axis=1)

    return df 






def rearrange_genomes_simulateSV(reference_genome, outdir, replace=False, nvars=50, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulated_svtype_to_svfile={}, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}):

    """Runs a simulation of nvars SVs of each type into the reference genome. It will be sepparated by gDNA and mtDNA. mtDNA will only include 5% of the gDNA vars. Everything is written to outdir. simulated_svtype_to_svfile is a dictionary that maps each svtype to a file with it. This function will insert these SVs plus the remaining ones up to nvars (which will be placed randomly in the remaining spots on the genome). Note that only balanced translocations will be simulated, as unbalanced translocations are hard to bechmark based on coverage. 

    Keep in mind that all RSVSim are 1-based coordinates"""


    # change the simulated simulated_svtype_to_svfile to not include translocations
    #simulated_svtype_to_svfile = {svtype : svfile for svtype, svfile in simulated_svtype_to_svfile.items() if svtype!="translocations"} # debug    

    # only random var
    #simulated_svtype_to_svfile = {}

    # map each chrom to a len
    chr_to_len = get_chr_to_len(reference_genome)

    # check that the vars are consistent
    check_consistency_of_svtype_to_svDF(simulated_svtype_to_svfile, set(chr_to_len))


    sakjhakahd

    # define the final outdirs
    final_simulated_SVs_dir = "%s/final_simulated_SVs"%(outdir); 
    final_rearranged_genome = "%s/rearranged_genome.fasta"%final_simulated_SVs_dir
    final_rearranged_genome_finalFile = "%s.performed"%(final_rearranged_genome)


    if file_is_empty(final_rearranged_genome_finalFile) or replace is True:

        # make the folder again
        #delete_folder(outdir) # debug
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
                svDF["ID"] = svDF.ID + "_realSV"
                svDF["Name"] = svDF.ID

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
            all_regions_bed_df = pd.DataFrame({chrom: {"start":1, "end":chrom_to_len[chrom]} for chrom in chroms}).transpose()
            all_regions_bed_df["chromosome"] = all_regions_bed_df.index
            all_regions_bed_df = all_regions_bed_df[["chromosome", "start", "end"]]

            # get the regions without SV where simulations should be placed
            bed_regions_prefix = "%s/bed_regions"%genome_outdir
            regions_without_SV_bed, svtype_to_nSVs = get_bed_df_not_overlapping_with_SVs(all_regions_bed_df, svtypes, simulated_svtype_to_svfile, bed_regions_prefix)

            # simulate random SVs into regions without previous SVs 
            random_sim_dir = "%s/random_SVs"%genome_outdir

            #### GET THE RANDOM INS,INV,DEL,TRA ####

            #if any([file_is_empty("%s/%s.tab"%(random_sim_dir, svtype)) for svtype in {"insertions", "deletions", "translocations", "inversions", "tandemDuplications"}]) or replace is True:
            if True:

                print("generating random SVs")

                # make and delete the folder
                delete_folder(random_sim_dir); make_folder(random_sim_dir)

                # get the cmd of the simulation
                randomSV_cmd = "%s --input_genome %s --outdir %s --regions_bed %s"%(create_random_simulatedSVgenome_R, genome_file, random_sim_dir, regions_without_SV_bed)

                # add the number of each SV that should be added
                svtype_to_arg = {"insertions":"number_Ins", "deletions":"number_Del", "inversions":"number_Inv", "translocations":"number_Tra", "tandemDuplications":"number_Dup"}
                #svtype_to_arg = {"insertions":"number_Ins", "deletions":"number_Del", "inversions":"number_Inv", "translocations":"number_Tra"}
            
                for svtype, number_alreadyGeneratedSVs in svtype_to_nSVs.items(): 
                    if svtype not in svtype_to_arg: continue

                    # define the number of vars to simulate depending on the type
                    if svtype=="translocations": real_vars_to_simulate = len(chroms)-1
                    else: real_vars_to_simulate = vars_to_simulate

                    randomSV_cmd += " --%s %i"%(svtype_to_arg[svtype], max([0, (real_vars_to_simulate-svtype_to_nSVs[svtype])]))

                # define the regions where to simulate translocations, as these are not 
                if "translocations" in simulated_svtype_to_svfile:

                    regions_without_translocations_or_otherSV_bed = get_bed_df_not_overlapping_with_translocations_allChromARM(regions_without_SV_bed, simulated_svtype_to_svfile["translocations"], genome_outdir, chr_to_len)

                else: regions_without_translocations_or_otherSV_bed = regions_without_SV_bed

                randomSV_cmd += " --regions_tra_bed %s"%regions_without_translocations_or_otherSV_bed

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

            # add the simulations into simulated_svtype_to_svDF
            for svtype in final_svtype_to_svDF.keys():
                svDF = final_svtype_to_svDF[svtype]

                # get the new sv
                new_svDF = pd.read_csv("%s/%s.tab"%(random_sim_dir, svtype), sep="\t")
                new_svDF = new_svDF[[c for c in new_svDF.keys() if "BpSeq" not in c]]

                # add the name
                new_svDF["ID"] = new_svDF.Name + "_sim_%s"%type_genome

                # append 
                final_svtype_to_svDF[svtype] = svDF.append(new_svDF, sort=True)

        ###############################################################################################################

        # debug to check that everything but translocations works
        #final_svtype_to_svDF = {svtype : svDF for svtype, svDF in final_svtype_to_svDF.items() if svtype!="translocations"}

        # check that the final_svtype_to_svDF is self consistent

        adkhghafgd


        ####### generate a rearranged genome with all the simulations in final_svtype_to_svDF #########
        print("inserting these random simulations")
        if file_is_empty(final_rearranged_genome) or replace is True:

            # initialize a cmd to create the simulated genome
            targetSV_cmd = "%s --input_genome %s --output_genome %s"%(create_targeted_simulatedSVgenome_R, reference_genome, final_rearranged_genome)

            # write the SVs into files and add to the cmd
            for svtype, svDF in final_svtype_to_svDF.items(): 

                # shift the insertions by 15 bp so that they are not at the beginning of the chrom
                if svtype=="insertions": svDF["StartA"] = svDF["StartA"] + 15

                # keep only the interesting svs
                svDF = svDF[svtype_to_fieldsDict[svtype]["all_fields"]]

                # write file
                svfile = "%s/%s.tab"%(final_simulated_SVs_dir, svtype)
                svDF.to_csv(svfile, sep="\t", header=True, index=False)

                # get cmd
                targetSV_cmd += " --%s_file %s"%(svtype, svfile)

            # run the cmd
            #std_rearranging_genome = "%s/simulation_std.txt"%final_simulated_SVs_dir
            std_rearranging_genome = "stdout"

            if std_rearranging_genome!="stdout": run_cmd("%s > %s 2>&1"%(targetSV_cmd, std_rearranging_genome))
            else: run_cmd(targetSV_cmd)

        jladnjkdjkdha

        # transform the cut-and-paste insertions to copy-and-paste, whenever necessary
        insertions_file = "%s/insertions.tab"%final_simulated_SVs_dir
        transform_cut_and_paste_to_copy_and_paste_insertions(reference_genome, final_rearranged_genome, insertions_file, final_svtype_to_svDF)

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
    print("simulations correctly generated")

    print(final_svtype_to_svfile)

    ljbadjkadjbad

    return final_svtype_to_svfile, final_rearranged_genome




######## random generation of balanced translocatipns #####




            # get the translocations randomly placed
            get_translocations_randomly_placed_in_target_regions(noInsInvDelTan_bed, translocations_file, chrom_to_len, nvars=nvars)


def get_translocations_randomly_placed_in_target_regions(target_regions_bed, translocations_file, chr_to_len, nvars=100):

    """Writes nvars randomly placed translocations in target_regions_bed, and writes them to translocations_file. It will draw as maximum number of translocations as possbile. Half of them will be inverted and half in the same orientation. All of them are balanced."""


    thisdoesnotwork

    print("getting randomly inserted translocations")

    # get the bed into a df
    target_regions_df = pd.read_csv(target_regions_bed, sep="\t", names=["chromosome", "start", "end"], header=None).drop_duplicates()
    target_regions_df.index = list(range(len(target_regions_df)))

    # add the breakpoint region in the middle
    target_regions_df["bp_pos"] = (target_regions_df.start + (target_regions_df.end-target_regions_df.start)/2).apply(int)

    # initialize a dict that will be used for the tra df
    varID_to_colName_to_value = {}

    # initialize a dict that will take each chromosome and simulate it
    chrom_to_nTRA = {chrom:0 for chrom in chr_to_len}

    # keep simulating translocations until you have nvars
    nvars_simulated = 0

    while nvars_simulated<nvars:

        # define the available chromosomes and their target regions
        available_chroms = {chrom for chrom, nTRA in chrom_to_nTRA.items() if nTRA<2}
        target_regions_df = target_regions_df[target_regions_df.chromosome.isin(available_chroms)]

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
        is_inverted = bool(random.randrange(0, 2))
        #is_inverted = False

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

        # update the number of times the chromosome is affected
        for chrom in [tra_dict["ChrA"], tra_dict["ChrB"]]: chrom_to_nTRA[chrom] += 1

        # update the number of simulated
        nvars_simulated+=1

    # get df
    tra_df = pd.DataFrame(varID_to_colName_to_value).transpose()[["ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Balanced"]]
    tra_df["Name"] = ["translocation_%i"%(I+1) for I in range(len(tra_df))]
    tra_df["ID"] = tra_df["Name"]

    print(tra_df)


    # write a subset of the vars (with 5 it works)
    #tra_df = tra_df.iloc[0:1] 




    tra_df.to_csv(translocations_file, sep="\t")

    print(chr_to_len)


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


   #### get the correct translocations together with all the simulations ####
        print("inserting final list of translocations...")

        # initialize a cmd to create the simulated genome
        #targetSV_cmd = "%s --input_genome %s --output_genome %s"%(create_targeted_simulatedSVgenome_R, reference_genome, rearranged_genome)

        # add all the positions of the different CMDs
        #for svtype in svtype_to_svDF.keys(): targetSV_cmd += " --%s_file %s/%s.tab"%(svtype, outdir, svtype)

        # run
        #run_cmd(targetSV_cmd)

        # not working

        ## #########################################################################


            # get the variants from simulating reads from an assembly. Always ploidy 1 to get homozygous SVs
            if type_data=="assembly":
                print("getting SVs from assemblies")

                predicted_svtype_to_svfile, df_gridss = generate_tables_of_SV_between_genomes_gridssClove(rearranged_genome, reference_genome, replace=replace, threads=threads)


# old way of considering assemblies

    """
    # OLD WAY
    # get the assembly
    assembly_name = "%s/%s.renamed.fasta"%(assemblies_dir, lrPrefix)
    copied_assembly_name = "%s/%s.nanopore_assembly.fasta"%(outdir, lrPrefix)

    # copy
    if file_is_empty(copied_assembly_name): run_cmd("cp %s %s"%(assembly_name , copied_assembly_name))

    # generate table
    dict_for_df[lrPrefix] = {"assembly":copied_assembly_name, "short_reads1":R1, "short_reads2":R2, "ID":lrPrefix, "short_reads_real1":realR1, "short_reads_real2":realR2}
    """


### OLD get_ID_to_svtype_to_svDF_for_setOfGenomes_highConfidence


  # get the real vars on a different way depending on the type of 
            if realSV_calling_on=="assembly":
                print("getting real variants based on the genome assemblies")

                # softlink the genome
                dest_genomeFile = "%s/genome_%s.fasta"%(all_realVars_dir, ID)
                if file_is_empty(dest_genomeFile): run_cmd("ln -s %s %s"%(row["assembly"], dest_genomeFile))

                # find the real vars
                svtype_to_svfile, df_gridss = generate_tables_of_SV_between_genomes_gridssClove(dest_genomeFile, reference_genome, replace=replace, threads=threads)



# if you stated auto in the reads, generate a 30x coverage bam file
if any([x=="auto" for x in {opt.fastq1, opt.fastq2}]):

    # define the number of reads as a function of the coverage
    read_length = 150
    total_nread_pairs = int((genome_length*30)/read_length)
    print("Simulating %.2fM reads"%(total_nread_pairs/1000000))

    sorted_bam, index_bam = fun.get_simulated_bamFile(opt.outdir, opt.ref, replace=opt.replace, threads=opt.threads, total_nread_pairs=total_nread_pairs, read_length=read_length)
    print("using simulated bam file from %s"%sorted_bam)



   #SNPthreshold_sameSample = get_fractionGenome_different_samplings_from_sorted_bam(sorted_bam, reference_genome, outdir_resamplingBam, replace=replace, threads=threads, coverage_subset_reads=coverage_subset_reads)

    print("We will say that if two samples differ by less than %.4f pct of the genome they are from the same sample. This has been calculated by resampling the input sorted bam with %ix coverage many times. We take twice the value of maxium divergence observed from these value."%(SNPthreshold_sameSample*100, coverage_subset_reads))



    # download the reference genome from GenBank given the taxID and also the gff annotation
    if opt.ref=="auto": opt.ref, opt.gff = fun.get_reference_genome_from_GenBank(opt.target_taxID, new_reference_genome_file, replace=opt.replace)

    # get by GenBank annotation
    elif opt.ref.startswith("GCA_"): opt.ref, opt.gff = fun.get_reference_genome_from_GenBank(opt.ref, new_reference_genome_file, replace=opt.replace)

    # just move the ref genome in the outdir
    else:





def get_SRA_runInfo_df_with_sampleID(SRA_runInfo_df, reference_genome, outdir, replace=False, threads=4, SNPthreshold=0.0001, coverage_subset_reads=10):

    """This function takes an SRA_runInfo_df and adds the sampleID. samples with the sample sampleID are those that have less tha SNPthreshold fraction of positions of the reference genome with SNPs. By default it is 0.01%. """

    make_folder(outdir)

    # change index
    SRA_runInfo_df = SRA_runInfo_df.set_index("Run", drop=False)

    # calculate the length of the genome
    length_genome = sum(get_chr_to_len(reference_genome).values())

    # add the subset_n_reads depending on the coverage and the read length
    SRA_runInfo_df["subset_n_reads"] = (length_genome*coverage_subset_reads / SRA_runInfo_df.avgLength).apply(int)

    ##### GET THE SNPS AND COVERAGE #####

    # get the SNPs for each run
    inputs_getSNPs_for_SRR = [(srr, reference_genome, "%s/%s"%(outdir, srr), SRA_runInfo_df.loc[srr, "subset_n_reads"], 1, replace) for srr in SRA_runInfo_df.Run]

    with multiproc.Pool(threads) as pool:
        list_vars_and_coverage = pool.starmap(getSNPs_for_SRR, inputs_getSNPs_for_SRR)
        pool.close()


    # change the vars
    list_vars = [x[0] for x in list_vars_and_coverage]
    list_mean_coverage = [x[1] for x in list_vars_and_coverage]
    list_fraction_genome_covered = [x[2] for x in list_vars_and_coverage]

    # add to the df
    SRA_runInfo_df["vars_set"] = list_vars
    SRA_runInfo_df["mean_coverage"] = list_mean_coverage
    SRA_runInfo_df["fraction_genome_covered"] = list_fraction_genome_covered

    ###################################

    # add the fraction of mapping reads
    print("calculating fraction of mapped reads")
    SRA_runInfo_df["fraction_reads_mapped"] = SRA_runInfo_df.Run.apply(lambda run: get_fraction_readPairsMapped("%s/%s/aligned_reads.bam.sorted"%(outdir, run), replace=replace, threads=threads))

    # add the divergence from the reference genome
    SRA_runInfo_df["fraction_genome_different_than_reference"] = SRA_runInfo_df.vars_set.apply(lambda x: len(x)/length_genome)

    # add the fraction of reads that were maintained 

    # initialize vars
    runA_to_runB_to_fraction_different_positions = {}

    # assign the  ID based on the comparison of SNPs 
    for runA in SRA_runInfo_df.Run:

        # get the snps
        snpsA = SRA_runInfo_df.loc[runA, "vars_set"]

        for runB in SRA_runInfo_df.Run:

            # get the snps
            snpsB = SRA_runInfo_df.loc[runB, "vars_set"]

            # calculate the fraction of positions of the genome that are different
            fraction_different_positions = get_fraction_differing_positions_AvsB(snpsA, snpsB, length_genome)
            runA_to_runB_to_fraction_different_positions.setdefault(runA, {}).setdefault(runB, fraction_different_positions)


    # get a df
    df_divergence = pd.DataFrame(runA_to_runB_to_fraction_different_positions)

    # get the clusters of IDs through a graph
    g =  igraph.Graph(directed=False)
    list_IDs = list(range(len(SRA_runInfo_df.Run)))
    srr_to_ID = dict(zip(SRA_runInfo_df.Run, list_IDs))
    ID_to_srr = dict(zip(list_IDs, SRA_runInfo_df.Run))
    g.add_vertices(list_IDs)
    pairs_equal_IDs = set.union(*[{tuple(sorted((srr_to_ID[runA], srr_to_ID[runB]))) for runB in SRA_runInfo_df.Run if df_divergence.loc[runA, runB]<SNPthreshold} for runA in SRA_runInfo_df.Run])
    g.add_edges(pairs_equal_IDs)
    clustersIDs = list(get_graph_subcomponents(g))

    # map each cluster to the SRRs
    cluster_names = list(range(len(clustersIDs)))
    clusterName_to_srrs = {name : {ID_to_srr[ID] for ID in IDs}  for name, IDs in dict(zip(cluster_names, clustersIDs)).items()}
    run_to_sampleID = {}
    for clusterName, srrs in clusterName_to_srrs.items():
        for srr in srrs: run_to_sampleID[srr] = clusterName+1

    # add to the df
    SRA_runInfo_df["sampleID"] = SRA_runInfo_df.Run.apply(lambda x: run_to_sampleID[x])

    # go through each sampleID and print the sample names
   # for s in set(SRA_runInfo_df.sampleID): print("Sample %i has these names: "%s, set(SRA_runInfo_df[SRA_runInfo_df.sampleID==s].SampleName))

    # drop the vars
    SRA_runInfo_df = SRA_runInfo_df.drop("vars_set", axis=1)

    return SRA_runInfo_df, df_divergence



original function
def get_close_shortReads_table_close_to_taxID(target_taxID, reference_genome, outdir, n_close_samples=3, nruns_per_sample=3, replace=False, threads=4, max_fraction_genome_different_than_reference=0.15, min_fraction_reads_mapped=0.8, min_fraction_genome_covered=0.8, coverage_subset_reads=5, min_coverage=30, min_fraction_coverage_subset_reads=0.5, run_in_slurm=False, walltime="02:00:00", queue="debug", StopAfter_sampleIndexingFromSRA=False, SNPthreshold_sameSample=0.001):

    """
    This function takes a taxID and returns the close_shortReads_table that is required to do optimisation of parameters
    """

    print("Getting genomes for taxID into %s"%(outdir))

    # load the NCBI taxonomy database and upgrade it if not already done
    print("getting NCBI taxonomy database")

    ncbi = NCBITaxa()
    ncbiTaxa_updated_file = "%s/ncbiTaxa_updated.txt"%outdir
    if file_is_empty(ncbiTaxa_updated_file) or replace is True: 

        # update
        ncbi.update_taxonomy_database()

        # write file
        open(ncbiTaxa_updated_file, "w").write("NCBItaxa updated\n")


    # calculate the expected difference between two runs of the same sample from the given sample
    outdir_resamplingBam = "%s/resampling_bam_andGetting_fractionDifPositions"%outdir

    print("We will say that if two samples differ by less than %.4f pct of the genome they are from the same sample. This has been provided by this function "%(SNPthreshold_sameSample*100))

    # get sampleID 
    outdir_gettingID = "%s/getting_sample_IDs"%outdir; make_folder(outdir_gettingID)

    # define the total number of runs that you need
    total_nruns = n_close_samples*nruns_per_sample
    print("Looking for %i runs"%total_nruns)

    SRA_runInfo_df_file = "%s/final_SRA_runInfo_df.py"%outdir

    if file_is_empty(SRA_runInfo_df_file) or replace is True:
        print("getting SRRs")

        # initialize a set that defines the runs of the previous node
        runs_previous_nodes = set()

        # define all potentially interesting taxIDs close to the target_taxIDs
        for nancestorNodes in range(1, 100): # one would mean to consider only IDs that are under the current species
            print("Considering %i ancestor nodes"%nancestorNodes)

            # create a folder for this number of ancestors
            outdir_ancestors = "%s/all_runsWithWGS_arround_target_taxID_%i_considering%iAncestors"%(outdir, target_taxID, nancestorNodes); make_folder(outdir_ancestors)

            # get the ancestor
            ancestor_taxID = ncbi.get_lineage(target_taxID)[-nancestorNodes]

            # get the tree
            tree = ncbi.get_descendant_taxa(ancestor_taxID, collapse_subspecies=False, return_tree=True, intermediate_nodes=True)

            # define interesting taxIDs (the leafs and the species names that may be intermediate)
            interesting_taxIDs = {int(x) for x in set(tree.get_leaf_names()).union({n.name for n in tree.traverse() if n.rank=="species"})}

            # map the distance between each leave and the target
            taxID_to_distanceToTarget = {taxID : tree.get_distance(str(target_taxID), str(taxID)) for taxID in interesting_taxIDs.difference({target_taxID})}
            taxID_to_distanceToTarget[target_taxID] = 0.0

            # get the taxIDs sorted by the distance (so that the closest )
            interesting_taxIDs_sorted = sorted(interesting_taxIDs, key=(lambda x: taxID_to_distanceToTarget[x]))

            # get the run info of all WGS datasets from SRA
            print("Getting WGS info")
            fileprefix = "%s/output"%(outdir_ancestors)
            all_SRA_runInfo_df = get_allWGS_runInfo_fromSRA_forTaxIDs(fileprefix, interesting_taxIDs_sorted, reference_genome, replace=True, min_coverage=min_coverage).set_index("Run", drop=False)

            # if you did not find anything get to farther ancestors
            if len(all_SRA_runInfo_df)<total_nruns: continue

            # reorder runs by coverage
            all_SRA_runInfo_df = all_SRA_runInfo_df.sort_values(by="expected_coverage", ascending=False)
            print("Looking for %i runs"%total_nruns)

            # go throough several fractions of all_SRA_runInfo_df
            size_chunk = max([(total_nruns*1), 1])

            # define a list that represents the idx of the runs
            idx_runs = list(range(len(all_SRA_runInfo_df)))

            # keep growing the all_SRA_runInfo_df until you find the desired number of runs.
            for chunk_run_idx in chunks(idx_runs, size_chunk):

                # get the runs
                SRA_runInfo_df = all_SRA_runInfo_df.iloc[0:chunk_run_idx[-1]]
                nruns = len(SRA_runInfo_df)
                print("Getting %i/%i runs"%(nruns, len(all_SRA_runInfo_df)))

                # get the sra info of this chunk
                SRA_runInfo_df, df_divergence = get_SRA_runInfo_df_with_sampleID(SRA_runInfo_df, reference_genome, outdir_gettingID, replace=replace, threads=threads, coverage_subset_reads=coverage_subset_reads, SNPthreshold=SNPthreshold_sameSample)


                ljdahjkhdadkjhaad

                # define the minimum coverage that the sample should have in order to pass, as afraction of the expected coverage
                min_mean_coverage = min_fraction_coverage_subset_reads*coverage_subset_reads

                print("These are the stats of this chunk of data")
                print(SRA_runInfo_df[["fraction_genome_different_than_reference", "fraction_reads_mapped", "fraction_genome_covered", "mean_coverage", "expected_coverage"]])

                # apply all the filters
                idx = ((SRA_runInfo_df.fraction_genome_different_than_reference<=max_fraction_genome_different_than_reference) &
                      (SRA_runInfo_df.fraction_genome_different_than_reference>=SNPthreshold_sameSample) & 
                      (SRA_runInfo_df.fraction_reads_mapped>=min_fraction_reads_mapped) &
                      (SRA_runInfo_df.fraction_genome_covered>=min_fraction_genome_covered) &
                      (SRA_runInfo_df.mean_coverage>=min_mean_coverage))

                SRA_runInfo_df = SRA_runInfo_df[idx]

                print("There are %i/%i runs that pass the filters"%(sum(idx), len(idx)))

                # add the number of runs that each sample has
                sampleID_to_nRuns = Counter(SRA_runInfo_df.sampleID)
                SRA_runInfo_df["nRuns_with_sampleID"] = SRA_runInfo_df.sampleID.apply(lambda x: sampleID_to_nRuns[x])

                # keep only the SRA_runInfo_df that has above the desired nruns per sample
                SRA_runInfo_df = SRA_runInfo_df[SRA_runInfo_df.nRuns_with_sampleID>=nruns_per_sample]

                # debug
                if len(set(SRA_runInfo_df.sampleID))<n_close_samples: continue

                # map each sampleID to the runs
                sampleID_to_runs = dict(SRA_runInfo_df.groupby("sampleID").apply(lambda df_s: set(df_s.Run)))

                # for each run, add the divergence to the other runs of the same sample
                SRA_runInfo_df["median_divergence_from_run_to_runsSameSample"] = SRA_runInfo_df.Run.apply(lambda run: np.median([df_divergence.loc[run][other_run] for other_run in sampleID_to_runs[SRA_runInfo_df.loc[run, "sampleID"]] if run!=other_run]) )

                # keep the nruns_per_sample that have the lowest median_divergence_from_run_to_runsSameSample nruns_per_sample
                interesting_runs = set.union(*[set(SRA_runInfo_df.loc[runs, "median_divergence_from_run_to_runsSameSample"].sort_values().iloc[0:nruns_per_sample].index) for sampleID, runs in sampleID_to_runs.items()])

                SRA_runInfo_df = SRA_runInfo_df.loc[interesting_runs]

                # redefine sampleID_to_runs with the closest divergence
                sampleID_to_runs = dict(SRA_runInfo_df.groupby("sampleID").apply(lambda df_s: set(df_s.Run)))

                # keep the n_close_samples that have the highest divergence between each other
                samplesDivergenceDict = {}; I = 0
                for sampleA, runsA in sampleID_to_runs.items():
                    for sampleB, runsB in sampleID_to_runs.items():

                        divergence = np.mean([np.mean(df_divergence.loc[runA][list(runsB)]) for runA in runsA])
                        samplesDivergenceDict[I] = {"sampleA":sampleA, "sampleB":sampleB, "divergence":divergence}

                        I+=1

                df_divergence_samples = pd.DataFrame(samplesDivergenceDict).transpose()
                for f in ["sampleA", "sampleB"]: df_divergence_samples[f] = df_divergence_samples[f].apply(int)

                # delete the comparisons that are the same
                df_divergence_samples["sampleA_and_sampleB"] = df_divergence_samples.apply(lambda r: tuple(sorted([r["sampleA"], r["sampleB"]])), axis=1)
                df_divergence_samples = df_divergence_samples.drop_duplicates(subset="sampleA_and_sampleB")
                df_divergence_samples = df_divergence_samples[df_divergence_samples.sampleA!=df_divergence_samples.sampleB].sort_values(by="divergence", ascending=False)

                # iterate through the df_divergence from the most divergent comparisons until you have n_close_samples
                final_samples = set()
                for sampleA, sampleB in df_divergence_samples[["sampleA", "sampleB"]].values:

                    # once you have all the samples break
                    if len(final_samples)>=n_close_samples: break

                    # add the samples
                    final_samples.update({sampleA, sampleB})

                # get the n_close_samples
                final_samples = list(final_samples)[0:n_close_samples]
                SRA_runInfo_df = SRA_runInfo_df[SRA_runInfo_df.sampleID.isin(final_samples)]

                # break the loop
                break

            # if there are no new nodes, break
            runs_in_this_node = set(SRA_runInfo_df.Run)
            if len(runs_in_this_node.difference(runs_previous_nodes))==0: break
            runs_previous_nodes.update(runs_in_this_node)

            # if you already found the IDs, break
            if len(SRA_runInfo_df)==total_nruns: break

        # debug
        if len(SRA_runInfo_df)!=total_nruns: raise ValueError("You could not find any datasets in SRA that would be useful")

        # load df
        save_object(SRA_runInfo_df, SRA_runInfo_df_file)

    else: SRA_runInfo_df = load_object(SRA_runInfo_df_file)

    print("these are the samples chosen:\n:", SRA_runInfo_df[["Run", "sampleID", "SampleName"]].sort_values("sampleID"))

    if StopAfter_sampleIndexingFromSRA is True: 
        print("stopping after generation of SRA_runInfo_df into %s"%SRA_runInfo_df_file)
        exit(0)

    ###### GETTING THE FINAL DATASETS ######

    # define the path to the final table
    close_shortReads_table = "%s/close_shortReads_table.tbl"%outdir

    if file_is_empty(close_shortReads_table) or replace is True:

        # change things
        SRA_runInfo_df["sampleID"] = "sample" + SRA_runInfo_df.sampleID.apply(str)
        SRA_runInfo_df["runID"] = SRA_runInfo_df.sampleID + "_" + SRA_runInfo_df.Run

        # define dirs
        downloads_dir = "%s/final_downloading_SRRs"%outdir; make_folder(downloads_dir)
        final_reads_dir = "%s/final_trimmed_reads_SRRs"%outdir; make_folder(final_reads_dir)

        # download the desired runs (this may be done locally)
        print("downloading each SRR")
        for srr in SRA_runInfo_df.Run: download_srr_parallelFastqDump(srr, "%s/%s"%(downloads_dir,srr), threads=threads, replace=replace)
        
        # define a dict that maps each srr to the reads of reads 
        srr_to_readsDict = {}

        # initialize the cmds to submit to the cluster
        all_cmds = [] 

        # replace each of the downloaded SRRs by the trimmed ones
        for srr in SRA_runInfo_df.Run:
            print("trimming %s reads"%srr)

            # define the raw reads
            reads1 = "%s/%s/%s_1.fastq.gz"%(downloads_dir, srr, srr)
            reads2 = "%s/%s/%s_2.fastq.gz"%(downloads_dir, srr, srr)

            # define the trimmed reads
            trimmed_reads1 = "%s.trimmed.fastq.gz"%reads1
            trimmed_reads2 = "%s.trimmed.fastq.gz"%reads2

            # define the final trimmed reads
            final_trimmed_reads1 = "%s/%s_1.fastq.gz"%(final_reads_dir, srr)
            final_trimmed_reads2 = "%s/%s_2.fastq.gz"%(final_reads_dir, srr)

            if file_is_empty(trimmed_reads1) or file_is_empty(trimmed_reads2) or replace is True:

                # get the trimmed reads
                cmd = "%s -f1 %s -f2 %s --threads %i"%(run_trimmomatic_and_fastqc_py, reads1, reads2, threads)
                if replace is True: cmd += " --replace"

                # add to cmds
                if run_in_slurm is True: 
                    all_cmds.append(cmd)
                    continue
                else: 

                    # get the trimmed reads
                    run_cmd(cmd)

            if file_is_empty(final_trimmed_reads1): os.rename(trimmed_reads1, final_trimmed_reads1)
            if file_is_empty(final_trimmed_reads2): os.rename(trimmed_reads2, final_trimmed_reads2)

            # add to the dict
            srr_to_readsDict[srr] = {"short_reads1":final_trimmed_reads1, "short_reads2":final_trimmed_reads2}

        # if there are all_cmds, run them in a job array
        if len(all_cmds)>0:
            print("Submitting %i trimmomatic jobs to the MN"%(len(all_cmds)))


            STDOUT = "%s/STDOUT"%final_reads_dir
            STDERR = "%s/STDERR"%final_reads_dir

            jobs_filename = "%s/jobs.trimming_SRAdatasets"%final_reads_dir
            open(jobs_filename, "w").write("\n".join(all_cmds))

          
            generate_jobarray_file_slurm(jobs_filename, stderr=STDERR, stdout=STDOUT, walltime=walltime,  name="trimming_SRAreads", queue=queue, sbatch=True, ncores_per_task=threads, rmstd=True, constraint="", number_tasks_to_run_at_once="all" )

            raise ValueError("You need to wait until the trimming of downloaded reads is finsihed")

        # add to the df
        for f in ["short_reads1", "short_reads2"]: SRA_runInfo_df[f] = SRA_runInfo_df.Run.apply(lambda srr: srr_to_readsDict[srr][f])

        # remove the downloads dir
        delete_folder(downloads_dir)

        # get the final df
        final_df = SRA_runInfo_df[["sampleID", "runID", "short_reads1", "short_reads2"]]

        # write
        final_df.to_csv(close_shortReads_table, sep="\t", header=True, index=False)


    print(final_df)

    dlkajhdkhadh

    #########################################

    return close_shortReads_table






            idx = ((SRA_runInfo_df.fraction_genome_different_than_reference<=max_fraction_genome_different_than_reference) &
                  (SRA_runInfo_df.fraction_genome_different_than_reference>=SNPthreshold_sameSample) & 

        # reorder runs by coverage
        all_SRA_runInfo_df = all_SRA_runInfo_df.sort_values(by="expected_coverage", ascending=False)
        print("Looking for %i runs"%total_nruns)

        # go throough several fractions of all_SRA_runInfo_df
        size_chunk = max([(total_nruns*2), 1])

        # define a list that represents the idx of the runs
        idx_runs = list(range(len(all_SRA_runInfo_df)))

        # keep growing the all_SRA_runInfo_df until you find the desired number of runs.
        for chunk_run_idx in chunks(idx_runs, size_chunk):

            # get the runs
            SRA_runInfo_df = all_SRA_runInfo_df.iloc[0:chunk_run_idx[-1]]
            nruns = len(SRA_runInfo_df)
            print("Getting %i/%i runs"%(nruns, len(all_SRA_runInfo_df)))

            # get the sra info of this chunk
            SRA_runInfo_df, df_divergence = get_SRA_runInfo_df_with_sampleID_popStructure(SRA_runInfo_df, reference_genome, outdir_gettingID, ploidy, replace=replace, threads=threads, coverage_subset_reads=coverage_subset_reads)

            #SRA_runInfo_df, df_divergence = get_SRA_runInfo_df_with_sampleID(SRA_runInfo_df, reference_genome, outdir_gettingID, replace=replace, threads=threads, coverage_subset_reads=coverage_subset_reads)


            ljdahjkhdadkjhaad

            # define the minimum coverage that the sample should have in order to pass, as afraction of the expected coverage
            min_mean_coverage = min_fraction_coverage_subset_reads*coverage_subset_reads

            print("These are the stats of this chunk of data")
            print(SRA_runInfo_df[["fraction_genome_different_than_reference", "fraction_reads_mapped", "fraction_genome_covered", "mean_coverage", "expected_coverage"]])

            # apply all the filters
            idx = ((SRA_runInfo_df.fraction_genome_different_than_reference<=max_fraction_genome_different_than_reference) &
                  (SRA_runInfo_df.fraction_genome_different_than_reference>=SNPthreshold_sameSample) & 
                  (SRA_runInfo_df.fraction_reads_mapped>=min_fraction_reads_mapped) &
                  (SRA_runInfo_df.fraction_genome_covered>=min_fraction_genome_covered) &
                  (SRA_runInfo_df.mean_coverage>=min_mean_coverage))

            SRA_runInfo_df = SRA_runInfo_df[idx]

            print("There are %i/%i runs that pass the filters"%(sum(idx), len(idx)))

            # add the number of runs that each sample has
            sampleID_to_nRuns = Counter(SRA_runInfo_df.sampleID)
            SRA_runInfo_df["nRuns_with_sampleID"] = SRA_runInfo_df.sampleID.apply(lambda x: sampleID_to_nRuns[x])

            # keep only the SRA_runInfo_df that has above the desired nruns per sample
            SRA_runInfo_df = SRA_runInfo_df[SRA_runInfo_df.nRuns_with_sampleID>=nruns_per_sample]

            # debug
            if len(set(SRA_runInfo_df.sampleID))<n_close_samples: continue

            # map each sampleID to the runs
            sampleID_to_runs = dict(SRA_runInfo_df.groupby("sampleID").apply(lambda df_s: set(df_s.Run)))

            # for each run, add the divergence to the other runs of the same sample
            SRA_runInfo_df["median_divergence_from_run_to_runsSameSample"] = SRA_runInfo_df.Run.apply(lambda run: np.median([df_divergence.loc[run][other_run] for other_run in sampleID_to_runs[SRA_runInfo_df.loc[run, "sampleID"]] if run!=other_run]) )

            # keep the nruns_per_sample that have the lowest median_divergence_from_run_to_runsSameSample nruns_per_sample
            interesting_runs = set.union(*[set(SRA_runInfo_df.loc[runs, "median_divergence_from_run_to_runsSameSample"].sort_values().iloc[0:nruns_per_sample].index) for sampleID, runs in sampleID_to_runs.items()])

            SRA_runInfo_df = SRA_runInfo_df.loc[interesting_runs]

            # redefine sampleID_to_runs with the closest divergence
            sampleID_to_runs = dict(SRA_runInfo_df.groupby("sampleID").apply(lambda df_s: set(df_s.Run)))

            # keep the n_close_samples that have the highest divergence between each other
            samplesDivergenceDict = {}; I = 0
            for sampleA, runsA in sampleID_to_runs.items():
                for sampleB, runsB in sampleID_to_runs.items():

                    divergence = np.mean([np.mean(df_divergence.loc[runA][list(runsB)]) for runA in runsA])
                    samplesDivergenceDict[I] = {"sampleA":sampleA, "sampleB":sampleB, "divergence":divergence}

                    I+=1

            df_divergence_samples = pd.DataFrame(samplesDivergenceDict).transpose()
            for f in ["sampleA", "sampleB"]: df_divergence_samples[f] = df_divergence_samples[f].apply(int)

            # delete the comparisons that are the same
            df_divergence_samples["sampleA_and_sampleB"] = df_divergence_samples.apply(lambda r: tuple(sorted([r["sampleA"], r["sampleB"]])), axis=1)
            df_divergence_samples = df_divergence_samples.drop_duplicates(subset="sampleA_and_sampleB")
            df_divergence_samples = df_divergence_samples[df_divergence_samples.sampleA!=df_divergence_samples.sampleB].sort_values(by="divergence", ascending=False)

            # iterate through the df_divergence from the most divergent comparisons until you have n_close_samples
            final_samples = set()
            for sampleA, sampleB in df_divergence_samples[["sampleA", "sampleB"]].values:

                # once you have all the samples break
                if len(final_samples)>=n_close_samples: break

                # add the samples
                final_samples.update({sampleA, sampleB})

            # get the n_close_samples
            final_samples = list(final_samples)[0:n_close_samples]
            SRA_runInfo_df = SRA_runInfo_df[SRA_runInfo_df.sampleID.isin(final_samples)]

            # break the loop
            break

        # if there are no new nodes, break
        runs_in_this_node = set(SRA_runInfo_df.Run)
        if len(runs_in_this_node.difference(runs_previous_nodes))==0: break
        runs_previous_nodes.update(runs_in_this_node)

        # if you already found the IDs, break
        if len(SRA_runInfo_df)==total_nruns: break

    # debug
    if len(SRA_runInfo_df)!=total_nruns: raise ValueError("You could not find any datasets in SRA that would be useful")


    ##### DOWNLOAD SRR FILES ####

    if run_in_cluster is False:

        threads_available = multiproc.cpu_count()

        # define all SRR files
        srr_to_SRRfile = {srr : "%s/%s/%s.srr"%(readsDir, srr, srr) for srr in all_SRA_runInfo_df.Run}

        # remove previous files
        for srr in srrs_to_remove:
            if srr in srr_to_SRRfile: 
                print("removing %s"%srr)
                fun.remove_file(srr_to_SRRfile[srr])

        # define the inputs of the downloads
        inputs_downloads = [(srr, SRRfile) for srr, SRRfile in srr_to_SRRfile.items() if fun.file_is_empty(SRRfile)]

        if len(inputs_downloads)>0:
            print("Downloading %i SRR files if not already done on %i threads"%(len(inputs_downloads), threads_available))

            with multiproc.Pool(threads_available) as pool:
                list_srr_files = pool.starmap(fun.download_srr_with_prefetch, inputs_downloads)
                pool.close()

        continue

    #############################


parser.add_argument("--queue_jobs", dest="queue_jobs", type=str, default="debug", help="The name of the queue were to submit the jobs when running with greasy")
parser.add_argument("--max_ncores_queue", dest="max_ncores_queue", type=int, default=768, help="The maximum number of cores that the queue can handle in a single job")

# timings of queues
parser.add_argument("--wallclock_read_obtention", dest="wallclock_read_obtention", type=str, default="02:00:00", help="The time that the fastqdumping of reads will take to perform this task")







def run_freebayes_pooledSeq(outdir_freebayes, ref, sorted_bam, ploidy, coverage, replace=False, threads=4):

    """Runs freebayes for pooled sequencing data, it writes the output as output.filt.vcf, although it is not filtered"""

    # define the output
    outvcf = "%s/output.filt.vcf"%outdir_freebayes; outvcf_tmp = "%s.tmp.vcf"%outvcf

    make_folder(outdir_freebayes)

    if file_is_empty(outvcf) or replace is True:
        print("running freebayes for pooled data into %s"%get_dir(outvcf))

        # not parallel
        #run_cmd("%s -f %s --haplotype-length -1 --min-alternate-count %i --min-alternate-fraction 0 --pooled-continuous -b %s -v %s"%(freebayes, ref, coverage, sorted_bam, outvcf_tmp))

        # generate a regions file
        regions_file = "%s.regions.tab"%outvcf_tmp
        run_cmd("%s %s.fai 100000 > %s"%(fasta_generate_regions_py, ref, regions_file))

        run_cmd("%s -f %s --haplotype-length -1 --min-alternate-count %i --min-alternate-fraction 0 --pooled-continuous -b %s -v %s --region ChrA_C_glabrata_CBS138:0-100000"%(freebayes, ref, coverage, sorted_bam, outvcf_tmp))

        akdjhjkhda

        """
        vcffirstheader \
        |   vcfstreamsort -w 1000 | vcfuniq # remove duplicates at region edges
        """


        run_cmd("%s %s %i -f %s --haplotype-length -1 --min-alternate-count %i --min-alternate-fraction 0 --pooled-continuous -use-best-n-alleles 20 %s > %s"%(freebayes_parallel, regions_file, threads, ref, coverage, sorted_bam, outvcf_tmp))

        remove_file(regions_file)


        os.rename(outvcf_tmp, outvcf)

    return outvcf






    # define the outdirs
    repeats_bed = "%s/repeats_table.bed"%outdir
    vcf_bed = "%s/variant_positions.bed"%outdir


    # make the bed file for the repeats
    repeats_df = pd.read_csv(repeats_table, sep="\t").rename(columns={"chromosome":"#chrom", "begin_repeat":"start", "end_repeat":"end"})
    repeats_df[["#chrom", "start", "end"]].to_csv(repeats_bed, sep="\t", header=True, index=False)

    # make a bed for the variants
    vcf_df["end"] = vcf_df.POS
    vcf_df.rename(columns={"#CHROM":"#chrom", "POS":"start"})[["#chrom", "start", "end"]].to_csv(vcf_bed, sep="\t", header=True, index=False)

    # run bedtools get the intersecting positions in vcf_df
    intersection_bed = "%s/intersection_vcf_and_repeats.bed"%outdir
    run_cmd("%s intersect -a %s -b %s -header > %s"%(bedtools, repeats_bed, vcf_bed, intersection_bed))

    # load the df and define the repeats variants
    intersecting_df = pd.read_csv(intersection_bed, sep="\t")
    variants_in_repeats = set(intersecting_df["#chrom"] + "_" + intersecting_df["start"].apply(str)).union(set(intersecting_df["#chrom"] + "_" + intersecting_df["end"].apply(str)))

    # debug the fact that there is no intersection
    if len(intersecting_df)==0: return [False]*len(vcf_df)

    print(intersecting_df)

    # define a series in vcf_df that has the variant as string
    vcf_df["var_as_str"] = vcf_df["#CHROM"] + "_" + vcf_df["POS"].apply(str)

    # check that all the variants_in_repeats are in vcf_df["var_as_str"]

    variants_in_repeats_not_in_vcf_df = variants_in_repeats.difference(set(vcf_df["var_as_str"]))
    if len(variants_in_repeats_not_in_vcf_df)>0:

        raise ValueError("There are %i/%i  variants in the intersection that can't be found in vcf_df"%(len(variants_in_repeats_not_in_vcf_df), len(variants_in_repeats)))



    # get the overlaps
    vcf_df["overlaps_repeats"] = vcf_df.var_as_str.isin(variants_in_repeats)

    print("There are %i/%i variants overlapping repeats"%(sum(vcf_df["overlaps_repeats"]), len(vcf_df)))

    ljadhjkadhkjadhk

    # clean
    for f in [repeats_bed, vcf_bed, intersection_bed]: remove_file(f)

    return vcf_df["overlaps_repeats"]




  # define an output file for VEP
    annotated_vcf = "%s_annotated.tab"%merged_vcf_all; annotated_vcf_tmp = "%s.tmp"%annotated_vcf

    # run annotation by VEP
    if fun.file_is_empty(annotated_vcf) or opt.replace is True:

        print("Annotating with VEP %s"%merged_vcf_all)

        # clean previous files
        fun.remove_file(annotated_vcf)
        fun.remove_file(annotated_vcf_tmp)
        for f in os.listdir(fun.get_dir(annotated_vcf)): 
            if ".tmp.raw." in f: fun.remove_file("%s/%s"%(fun.get_dir(annotated_vcf), f))


        vep_cmd = "%s --input_vcf %s --outfile %s --ref %s --gff %s --mitochondrial_chromosome %s --mito_code %i --gDNA_code %i "%(run_vep, merged_vcf_all, annotated_vcf_tmp, opt.ref, gff_with_biotype, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code)

        fun.run_cmd(vep_cmd)

        os.rename(annotated_vcf_tmp, annotated_vcf)
        
# test the accuracy on each of the simulations types
if opt.testSimulationsAccuracy is True: fun.report_accuracy_simulations(sorted_bam, opt.ref, "%s/testing_SimulationsAccuracy"%opt.outdir, real_svtype_to_file, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, job_array_mode=opt.job_array_mode)

# test accuracy on real data
if opt.testRealDataAccuracy is True:  fun.report_accuracy_realSVs(opt.close_shortReads_table, opt.ref, "%s/testing_RealSVsAccuracy"%opt.outdir, real_svtype_to_file, outdir_finding_realVars, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_ploidies=simulation_ploidies, range_filtering_benchmark=opt.range_filtering_benchmark, nvars=opt.nvars, job_array_mode=opt.job_array_mode)



#### bamqc
if opt.run_qualimap is True:
    
    bamqc_outdir = "%s/bamqc_out"%opt.outdir
    if fun.file_is_empty("%s/qualimapReport.html"%bamqc_outdir) or opt.replace is True:
        print("Running bamqc to analyze the bam alignment")
        qualimap_std = "%s/std.txt"%bamqc_outdir
        try: bamqc_cmd = "%s bamqc -bam %s -outdir %s -nt %i > %s 2>&1"%(qualimap, sorted_bam, bamqc_outdir, opt.threads, qualimap_std); fun.run_cmd(bamqc_cmd)
        except: print("WARNING: qualimap failed likely due to memory errors, check %s"%qualimap_std)

parser.add_argument("--run_qualimap", dest="run_qualimap", action="store_true", default=False, help="Run qualimap for quality assessment of bam files. This may be inefficient sometimes because of the ")


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

wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.4/vcf_validator_linux


vcf_validator = "%s/vcf_validator_linux"%external_software

def validate_vcf(vcf_file, replace=False):

    """This function takes a vcf and reports whether it is valid according to vcf_validator"""

    # define the outdir
    outdir = "%s_vcf_validator_outdir"%vcf_file; 
    delete_folder(outdir); make_folder(outdir)

    print("running vcf_validator into %s"%outdir)

    run_cmd("%s -i %s -l warning -r summary,text -o %s --require-evidence"%(vcf_validator, vcf_file, outdir))



def assembly_row_has_annotation(r, prefix):

    """Takes a row of the get_GenBank_assembly_statistics_df and returns whether it has a gff annotation"""

    print("looking for annotation")

    # define the ftp direction
    ftp_site = r["ftp_path"]
    ftp_last = ftp_site.split("/")[-1]# something that has to be added for the ftp direction
    full_ftp_direction = ftp_site+"/"+ftp_last

    # define the gff name
    origin_file = "%s_genomic.gff.gz"%(full_ftp_direction)
    dest_file = "%s_testing.gff"%prefix


    # print getting file
    try: 
        urllib.request.urlretrieve(origin_file, dest_file)
        remove_file(dest_file)
        return True

    except: return False

def get_GenBank_assembly_statistics_df(file, assembly_summary_genbank_url="ftp://ftp.ncbi.nih.gov/genomes/genbank/assembly_summary_genbank.txt", replace=False):

    """
    Downloads the assembly summary statistics into file and returns a df.

    Index(['# assembly_accession', 'bioproject', 'biosample', 'wgs_master',
       'refseq_category', 'taxid', 'species_taxid', 'organism_name',
       'infraspecific_name', 'isolate', 'version_status', 'assembly_level',
       'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter',
       'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
       'excluded_from_refseq', 'relation_to_type_material'],
      dtype='object')
    """

    df_file = "%s.df.py"%file

    if file_is_empty(df_file) or replace is True:
        print("getting GeneBank genomes")

        # download the file 
        urllib.request.urlretrieve(assembly_summary_genbank_url, file)

        # get into df and save
        df = pd.read_csv(file, header=1, sep="\t").rename(columns={"# assembly_accession":"assembly_accession"}) 

        save_object(df, df_file)

    else: df = load_object(df_file)

    return df



def get_reference_genome_from_GenBank(ID, dest_genome_file, replace=False):

    """Downloads the reference genome for a given ID into dest_genome_file from genBank """

    # get the genBank df
    gb_file = "%s.GenBankAssemblySummary.txt"%dest_genome_file
    df_assemblies = get_GenBank_assembly_statistics_df(gb_file, replace=replace)

    # it is a genbank assembly just get the accession
    if ID.startswith("GCA_"): rep_genomes = df_assemblies[df_assemblies.assembly_accession==ID]

    # it is a taxID
    else:

        # get the taxID genomes
        df_taxID = df_assemblies[(df_assemblies.taxid==ID) | (df_assemblies.species_taxid==ID)]

        # get the reference genomes
        rep_genomes = df_taxID[(df_taxID.refseq_category.isin({"reference genome", "representative genome"}))]

    if len(rep_genomes)==0: 

        df_taxID.to_csv("%s/../testing/failed_assemblies_taxID%i"%(CWD, ID), sep="\t", header=True, index=False)
        raise ValueError("There are no representative genomes in GenBank for this ID")

    if len(rep_genomes)>1: print("WARNING: There are more than 1 reference genomes for this ID")

    # download genome and annotations
    ftp_site = rep_genomes.iloc[0]["ftp_path"]

    # get the full_ftp_direction
    ftp_last = ftp_site.split("/")[-1]# something that has to be added for the ftp direction
    full_ftp_direction = ftp_site+"/"+ftp_last

    # map each file type to the final file
    file_type_to_finalFile = {"report":"%s.assembly_report.txt"%dest_genome_file, "genome":dest_genome_file, "gff":"%s.features.gff"%dest_genome_file}

    # download the genome and annotations 
    for file_type, file_name in [("report", "assembly_report.txt"), ("genome", "genomic.fna.gz"), ("gff", "genomic.gff.gz")]:

        # define the files
        origin_file = "%s_%s"%(full_ftp_direction, file_name)
        dest_file = "%s.%s"%(dest_genome_file, file_name)
        dest_file_unzipped = dest_file.rstrip(".gz")
        final_file = file_type_to_finalFile[file_type]

        if file_is_empty(final_file) or replace is True:
            print("getting", file_type, "from genBank")

        
            try:

                # print getting file
                urllib.request.urlretrieve(origin_file, dest_file)

                # unzip if necessary
                if file_type in {"genome", "gff"}: run_cmd("gunzip --force %s"%dest_file)

                # move to the final file
                os.rename(dest_file_unzipped, final_file)

            except: print("WARNING: %s could not be found"%file_type)

    # define the files to return
    genome_file = file_type_to_finalFile["genome"]
    gff_file = file_type_to_finalFile["gff"]
    if file_is_empty(gff_file): gff_file = None

    return genome_file, gff_file






def get_SRA_runInfo_df_with_sampleID_popStructure(SRA_runInfo_df, reference_genome, outdir, ploidy, replace=False, threads=4, coverage_subset_reads=10):

    """This function takes an SRA_runInfo_df and adds the sampleID. The sampleID is defined as samples that are from the same population."""

    make_folder(outdir)

    # define the outdir that will store the seq data
    seq_data_dir = "%s/seq_data"%outdir; make_folder(seq_data_dir)

    # change index
    SRA_runInfo_df = SRA_runInfo_df.set_index("Run", drop=False)

    # calculate the length of the genome
    length_genome = sum(get_chr_to_len(reference_genome).values())

    # add the subset_n_reads depending on the coverage and the read length
    SRA_runInfo_df["subset_n_reads"] = (length_genome*coverage_subset_reads / SRA_runInfo_df.avgLength).apply(int)

    ###### DOWNLOAD THE SUBSET OF READS WITH FASTQDUMP ####

    inputs_download_srr = [(srr, "%s/%s/reads_dir"%(seq_data_dir, srr), SRA_runInfo_df.loc[srr, "subset_n_reads"]) for srr in SRA_runInfo_df.Run]

    with multiproc.Pool(threads) as pool:
        list_reads_tuples = pool.starmap(download_srr_subsetReads_onlyFastqDump, inputs_download_srr)
        pool.close()

    #######################################################

    ##### GET THE SNPs CALLED AND COVERAGE ##### 

    # initialize lists
    mean_coverage_list = []
    fraction_genome_covered_list = []
    vcf_list = []

    # get the files
    for I, (reads1, reads2) in enumerate(list_reads_tuples):
        #print("working on sample %i"%I)

        # define the outdir
        srr = SRA_runInfo_df.iloc[I]["Run"]
        outdir_srr = "%s/%s"%(seq_data_dir, srr)

        vcf, mean_coverage, fraction_genome_covered = get_vcf_and_coverage_for_reads(reads1, reads2, reference_genome, outdir_srr, ploidy, threads=threads, replace=replace)

        # keep
        mean_coverage_list.append(mean_coverage)
        fraction_genome_covered_list.append(fraction_genome_covered)
        vcf_list.append(vcf)


    SRA_runInfo_df["mean_coverage"] = mean_coverage_list
    SRA_runInfo_df["fraction_genome_covered"] = fraction_genome_covered_list
    SRA_runInfo_df["vcf"] = vcf_list

    ############################################

    # add the sampleID by popStructure
    pop_structure_dir = "%s/pop_structure"%outdir; make_folder(pop_structure_dir)
    SRA_runInfo_df = get_SRA_runInfo_df_with_populationID(SRA_runInfo_df, pop_structure_dir, replace=replace)


    # add the fraction of mapping reads
    print("calculating fraction of mapped reads")
    SRA_runInfo_df["fraction_reads_mapped"] = SRA_runInfo_df.Run.apply(lambda run: get_fraction_readPairsMapped("%s/%s/aligned_reads.bam.sorted"%(seq_data_dir, run), replace=replace, threads=threads))

    print(SRA_runInfo_df["fraction_reads_mapped"])


def get_vcf_and_coverage_for_reads(reads1, reads2, reference_genome, outdir, ploidy, threads=4, replace=False):

    """This function runs GATK without coverage filtering for some reads and mosdepth to return a vcf with SNPs   """

    # get the trimmed reads
    #print("running trimmomatic")
    trimmed_reads1, trimmed_reads2 = run_trimmomatic(reads1, reads2, replace=replace, threads=threads)

    # get the aligned reads
    #print("running bwa mem")
    bamfile = "%s/aligned_reads.bam"%outdir
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam

    run_bwa_mem(trimmed_reads1, trimmed_reads2, reference_genome, outdir, bamfile, sorted_bam, index_bam, "nameSample", threads=threads, replace=replace)

    # get the variants vcf
    
    # run freebayes
    outdir_freebayes = "%s/freebayes_ploidy%i_out"%(outdir, ploidy)
    vcf =  run_freebayes_parallel(outdir_freebayes, reference_genome, sorted_bam, ploidy, 1, replace=replace) 

    # define the parallel running of mosdepth 
    if threads==1: run_in_parallel=False
    else: run_in_parallel = True

    # get the coverage df
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=True), sep="\t")

    # define stats
    mean_coverage = np.mean(coverage_df.mediancov_1)
    fraction_genome_covered = np.mean(coverage_df.percentcovered_1)/100

    #print("The mean coverage is %.3f for windows of %ibp"%(mean_coverage, window_l))
    #print("The mean fraction coverage for is %.3f for windows of %ibp"%(fraction_genome_covered, window_l))

    return vcf, mean_coverage, fraction_genome_covered

def run_freebayes_for_chromosome(chromosome_id, outvcf_folder, ref, sorted_bam, ploidy, coverage, replace, pooled_sequencing):

    """Takes a chromosome ID and the fasta file and an outvcf and runs freebayes on it"""

    # define the output vcf file
    outvcf = "%s/%s_freebayes.vcf"%(outvcf_folder, chromosome_id); outvcf_tmp = "%s.tmp.vcf"%outvcf
    #print("running freebayes for %s"%chromosome_id)

    # remove previously existing files
    if file_is_empty(outvcf) or replace is True:
        remove_file(outvcf_tmp)

        """
        optional removal of bam file
        # generate the bam file for this chromosome (and index)
        sorted_bam_chr = "%s.%s.bam"%(sorted_bam, chromosome_id)
        run_cmd("%s view -b %s %s > %s"%(samtools, sorted_bam, chromosome_id, sorted_bam_chr))
        run_cmd("%s index -@ 1 %s"%(samtools, sorted_bam_chr))

        """

        # map each chromosome to a length
        chr_to_len = get_chr_to_len(ref)

        # get the fasta for the chromosome
        fasta_chromosome = "%s.%s.fasta"%(ref, chromosome_id)
        SeqIO.write([seq for seq in SeqIO.parse(ref, "fasta") if seq.id==chromosome_id], fasta_chromosome, "fasta")

        # run freebayes
        freebayes_std = "%s.std"%outvcf_tmp
        print("running freebayes with STD %s"%freebayes_std)

        # define the region
        region = "%s:0-%i"%(chromosome_id, chr_to_len[chromosome_id])

        if pooled_sequencing is True:
            print("running for pooled data")
            run_cmd("%s -f %s --haplotype-length -1 --use-best-n-alleles 20 --min-alternate-count %i --min-alternate-fraction 0 --pooled-continuous -b %s -v %s --region %s > %s 2>&1"%(freebayes, fasta_chromosome, coverage, sorted_bam, outvcf_tmp, region, freebayes_std))
        else: 
            print("running unpooled sequencing")
            run_cmd("%s -f %s -p %i --min-coverage %i -b %s --haplotype-length -1 -v %s --region %s > %s 2>&1"%(freebayes, fasta_chromosome, ploidy, coverage, sorted_bam, outvcf_tmp, region, freebayes_std))

        # remove the intermediate files
        #print("%s exists %s"%(fasta_chromosome, str(file_is_empty(fasta_chromosome))))
        remove_file(fasta_chromosome); remove_file("%s.fai"%fasta_chromosome); remove_file(freebayes_std)

        #remove_file(sorted_bam_chr); remove_file("%s.bai"%sorted_bam_chr); remove_file(fasta_chromosome); remove_file("%s.fai"%fasta_chromosome); remove_file(freebayes_std)

        # rename
        os.rename(outvcf_tmp, outvcf)

    # return the vcfs
    return outvcf



def run_freebayes_parallel(outdir_freebayes, ref, sorted_bam, ploidy, coverage, threads=4, max_threads=8, replace=False, pooled_sequencing=False):

    """It parallelizes over the current CPUs of the system"""

    # make the dir if not already done
    if not os.path.isdir(outdir_freebayes): os.mkdir(outdir_freebayes)

    #run freebayes
    freebayes_output ="%s/output.raw.vcf"%outdir_freebayes; freebayes_output_tmp = "%s.tmp"%freebayes_output
    if file_is_empty(freebayes_output) or replace is True:

        #print("running freebayes in parallel with %i threads"%(multiproc.cpu_count()))

        # define the chromosomes
        all_chromosome_IDs = [seq.id for seq in SeqIO.parse(ref, "fasta")]

        # remove the previous tmp file
        if not file_is_empty(freebayes_output_tmp): os.unlink(freebayes_output_tmp)

        # initialize the pool class with the available CPUs --> this is asyncronous parallelization
        threads = min([max_threads, threads])
        with multiproc.Pool(threads) as pool:

            # make a dir to store the vcfs
            chromosome_vcfs_dir = "%s/chromosome_vcfs"%outdir_freebayes; make_folder(chromosome_vcfs_dir)

            # run in parallel the freebayes generation for all the 
            chromosomal_vcfs = pool.starmap(run_freebayes_for_chromosome, [(ID, chromosome_vcfs_dir, ref, sorted_bam, ploidy, coverage, replace, pooled_sequencing) for ID in all_chromosome_IDs])

            # close the pool
            pool.close()
            pool.terminate()


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

        # remove tmp vcfs
        delete_folder(chromosome_vcfs_dir)

        # rename
        os.rename(freebayes_output_tmp, freebayes_output)

    # filter the freebayes by quality
    freebayes_filtered = "%s/output.filt.vcf"%outdir_freebayes; freebayes_filtered_tmp = "%s.tmp"%freebayes_filtered
    if file_is_empty(freebayes_filtered) or replace is True:

        if pooled_sequencing is True: soft_link_files(freebayes_output, freebayes_filtered)
        else:

            #print("filtering freebayes")
            cmd_filter_fb = '%s -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" --tag-pass PASS %s > %s'%(vcffilter, freebayes_output, freebayes_filtered_tmp); run_cmd(cmd_filter_fb)
            os.rename(freebayes_filtered_tmp, freebayes_filtered)

    return freebayes_filtered




def get_SRA_runInfo_df_with_sampleID(SRA_runInfo_df, reference_genome, outdir, replace=False, threads=4, SNPthreshold=0.0001, coverage_subset_reads=10):

    """This function takes an SRA_runInfo_df and adds the sampleID. samples with the sample sampleID are those that have less tha SNPthreshold fraction of positions of the reference genome with SNPs. By default it is 0.01%. """

    make_folder(outdir)

    # change index
    SRA_runInfo_df = SRA_runInfo_df.set_index("Run", drop=False)

    # calculate the length of the genome
    length_genome = sum(get_chr_to_len(reference_genome).values())

    # add the subset_n_reads depending on the coverage and the read length
    SRA_runInfo_df["subset_n_reads"] = (length_genome*coverage_subset_reads / SRA_runInfo_df.avgLength).apply(int)

    ##### GET THE SNPS AND COVERAGE #####

    # get the SNPs for each run
    inputs_getSNPs_for_SRR = [(srr, reference_genome, "%s/%s"%(outdir, srr), SRA_runInfo_df.loc[srr, "subset_n_reads"], 1, replace) for srr in SRA_runInfo_df.Run]

    with multiproc.Pool(threads) as pool:
        list_vars_and_coverage = pool.starmap(getSNPs_for_SRR, inputs_getSNPs_for_SRR)
        pool.close()


    # change the vars
    list_vars = [x[0] for x in list_vars_and_coverage]
    list_mean_coverage = [x[1] for x in list_vars_and_coverage]
    list_fraction_genome_covered = [x[2] for x in list_vars_and_coverage]

    # add to the df
    SRA_runInfo_df["vars_set"] = list_vars
    SRA_runInfo_df["mean_coverage"] = list_mean_coverage
    SRA_runInfo_df["fraction_genome_covered"] = list_fraction_genome_covered

    ###################################

    # add the fraction of mapping reads
    print("calculating fraction of mapped reads")
    SRA_runInfo_df["fraction_reads_mapped"] = SRA_runInfo_df.Run.apply(lambda run: get_fraction_readPairsMapped("%s/%s/aligned_reads.bam.sorted"%(outdir, run), replace=replace, threads=threads))

    # add the divergence from the reference genome
    SRA_runInfo_df["fraction_genome_different_than_reference"] = SRA_runInfo_df.vars_set.apply(lambda x: len(x)/length_genome)

    # add the fraction of reads that were maintained 

    # initialize vars
    runA_to_runB_to_fraction_different_positions = {}

    # assign the  ID based on the comparison of SNPs 
    for runA in SRA_runInfo_df.Run:

        # get the snps
        snpsA = SRA_runInfo_df.loc[runA, "vars_set"]

        for runB in SRA_runInfo_df.Run:

            # get the snps
            snpsB = SRA_runInfo_df.loc[runB, "vars_set"]

            # calculate the fraction of positions of the genome that are different
            fraction_different_positions = get_fraction_differing_positions_AvsB(snpsA, snpsB, length_genome)
            runA_to_runB_to_fraction_different_positions.setdefault(runA, {}).setdefault(runB, fraction_different_positions)


    # get a df
    df_divergence = pd.DataFrame(runA_to_runB_to_fraction_different_positions)

    # get the clusters of IDs through a graph
    g =  igraph.Graph(directed=False)
    list_IDs = list(range(len(SRA_runInfo_df.Run)))
    srr_to_ID = dict(zip(SRA_runInfo_df.Run, list_IDs))
    ID_to_srr = dict(zip(list_IDs, SRA_runInfo_df.Run))
    g.add_vertices(list_IDs)
    pairs_equal_IDs = set.union(*[{tuple(sorted((srr_to_ID[runA], srr_to_ID[runB]))) for runB in SRA_runInfo_df.Run if df_divergence.loc[runA, runB]<SNPthreshold} for runA in SRA_runInfo_df.Run])
    g.add_edges(pairs_equal_IDs)
    clustersIDs = list(get_graph_subcomponents(g))

    # map each cluster to the SRRs
    cluster_names = list(range(len(clustersIDs)))
    clusterName_to_srrs = {name : {ID_to_srr[ID] for ID in IDs}  for name, IDs in dict(zip(cluster_names, clustersIDs)).items()}
    run_to_sampleID = {}
    for clusterName, srrs in clusterName_to_srrs.items():
        for srr in srrs: run_to_sampleID[srr] = clusterName+1

    # add to the df
    SRA_runInfo_df["sampleID"] = SRA_runInfo_df.Run.apply(lambda x: run_to_sampleID[x])

    # go through each sampleID and print the sample names
   # for s in set(SRA_runInfo_df.sampleID): print("Sample %i has these names: "%s, set(SRA_runInfo_df[SRA_runInfo_df.sampleID==s].SampleName))

    # drop the vars
    SRA_runInfo_df = SRA_runInfo_df.drop("vars_set", axis=1)

    return SRA_runInfo_df, df_divergence



def getSNPs_for_SRR(srr, reference_genome, outdir, subset_n_reads=100000, threads=1, replace=False):

    """This function runs fast SNP calling for an SRR, saving data into outdir an srr. It returns the path to the VCF file. By default it runs on one core"""


    start_time = time.time()

    print("running getSNPs_for_SRR for %i reads "%subset_n_reads)

    # make the outdir 
    make_folder(outdir)

    # first get the reads into a downloading dir
    reads_dir = "%s/reads_dir"%outdir
    #print("downloading fastq files")
    reads1, reads2 = download_srr_subsetReads_onlyFastqDump(srr, reads_dir, subset_n_reads=subset_n_reads)

    # get the trimmed reads
    #print("running trimmomatic")
    trimmed_reads1, trimmed_reads2 = run_trimmomatic(reads1, reads2, replace=replace, threads=threads)

    # get the aligned reads
    #print("running bwa mem")
    bamfile = "%s/aligned_reads.bam"%outdir
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam

    run_bwa_mem(trimmed_reads1, trimmed_reads2, reference_genome, outdir, bamfile, sorted_bam, index_bam, srr, threads=threads, replace=replace)

    # get the SNPs
    #print("getting SNPs ")
    snps_set = get_SNPs_from_bam(sorted_bam, outdir, reference_genome, replace=replace, threads=threads)

    # define the parallel running of mosdepth 
    if threads==1: run_in_parallel=False
    else: run_in_parallel = True

    # get the coverage df
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=run_in_parallel), sep="\t")

    # define stats
    mean_coverage = np.mean(coverage_df.mediancov_1)
    fraction_genome_covered = np.mean(coverage_df.percentcovered_1)/100

    #print("The mean coverage for %s is %.3f for windows of %ibp"%(srr, mean_coverage, window_l))
    #print("The mean fraction coverage for %s is %.3f for windows of %ibp"%(srr, fraction_genome_covered, window_l))
    #print(coverage_df)


    #print("--- the running of getSNPs_for_SRR took %s seconds in %i cores for a mean_coverage=%.3f ---"%(time.time() - start_time, threads, mean_coverage))


    return snps_set, mean_coverage, fraction_genome_covered



def get_fractionGenome_different_samplings_from_sorted_bam(sorted_bam, reference_genome, outdir, replace=False, threads=4, coverage_subset_reads=5, nsamples=3):

    """ This function makes nsamples (of subset_n_reads each) of a sorted bam and calculates the fraction of positions of the genome that are different between the different samples. It returns the mean of all the measurements. """

    make_folder(outdir)

    # count number of reads
    npairs = count_number_read_pairs(sorted_bam, replace=replace, threads=threads)

    # count the read length
    read_length = get_read_length(sorted_bam, threads=threads, replace=replace)

    # get the length of the reference genome
    length_genome = sum(get_chr_to_len(reference_genome).values())

    # expected coverage
    expected_coverage = (npairs*read_length)/length_genome

    # calculate the fraction of reads that each sample should have
    fraction_reads_per_sample = coverage_subset_reads/expected_coverage
    print("Subsampling %.4f of reads to get %ix coverage"%(fraction_reads_per_sample, coverage_subset_reads))

    # run for each sample in parallel
    inputs_get_SNPs_for_a_sample_of_a_bam = [(sorted_bam, "%s/sample%i"%(outdir, sample), reference_genome, fraction_reads_per_sample, replace, 1) for sample in range(nsamples)]

    # run in parallel
    with multiproc.Pool(threads) as pool:
        list_vars = pool.starmap(get_SNPs_for_a_sample_of_a_bam, inputs_get_SNPs_for_a_sample_of_a_bam)
        pool.close()

    # get the fraction of different positions
    all_fraction_different_positions = []

    # go through each combination of vars
    for snpsA in list_vars:
        for snpsB in list_vars:
            if snpsA==snpsB: continue

            # get the fraction of SNPs that are different
            fraction_different_positions = get_fraction_differing_positions_AvsB(snpsA, snpsB, length_genome)
            all_fraction_different_positions.append(fraction_different_positions) 

    return max(all_fraction_different_positions)



def get_SNPs_for_a_sample_of_a_bam(sorted_bam, outdir, reference_genome, fraction_reads=0.1, replace=False, threads=4):

    """This function samples a bam and gets a subsample and the SNPs on this Subsample """

    print("getting SNPs")

    make_folder(outdir)

    # define the sampled bamfile
    sampled_bamfile = "%s/sampled_bamfile.bam"%outdir
    sampled_bamfile_tmp = "%s.tmp"%sampled_bamfile

    if file_is_empty(sampled_bamfile) or replace is True:

        remove_file(sampled_bamfile_tmp)

        # get the subsampled bam
        downsample_bamfile_keeping_pairs(sorted_bam, fraction_reads=fraction_reads, replace=replace, threads=threads, sampled_bamfile=sampled_bamfile_tmp)

        os.rename(sampled_bamfile_tmp, sampled_bamfile)

    # index the bam
    index_bam = "%s.bai"%sampled_bamfile
    if file_is_empty(index_bam) or replace is True:
        #print("Indexing bam")
        run_cmd("%s index -@ %i %s"%(samtools, threads, sampled_bamfile))

    # get the snps
    snps = get_SNPs_from_bam(sampled_bamfile, outdir, reference_genome, replace=replace, threads=threads)


    #### FIND COVERAGE ####

    """
    THISISNECESSARYFORCHECKINGCOVERAGE

    # define the parallel running of mosdepth 
    if threads==1: run_in_parallel=False
    else: run_in_parallel = True

    # get the coverage df
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir, sampled_bamfile, windows_file="none", replace=replace, run_in_parallel=run_in_parallel), sep="\t")

    # define stats
    mean_coverage = np.mean(coverage_df.mediancov_1)
    fraction_genome_covered = np.mean(coverage_df.percentcovered_1)/100

    print("The mean coverage for this bam sample is %.3f for windows of 10Kb"%(mean_coverage))
    print("The mean fraction coverage for this bam sample is %.3f for windows of 10Kb"%(fraction_genome_covered))

    """

    #######################

    return snps



def get_SNPs_from_bam(sorted_bam, outdir, reference_genome, replace=False, threads=4):

    """Runs freebayes on the sorted bam and returns a set with the SNPs found with a fast freebayes for ploidy 1"""

    make_folder(outdir)

    # read calling with freebayes
    vcf = run_freebayes_withoutFiltering(outdir, reference_genome, sorted_bam, 1, threads, 2, replace=replace)

    # get the vcf as a df
    df = pd.read_csv(vcf, skiprows=list(range(len([line for line in open(vcf, "r", encoding='utf-8', errors='ignore') if line.startswith("##")]))), sep="\t", na_values=vcf_strings_as_NaNs, keep_default_na=False)


    # kepp only snps
    df["is_snp"] = (df.REF.isin({"A", "C", "T", "G"})) & (df.ALT.isin({"A", "C", "T", "G"}))
    df  = df[df["is_snp"]]

    # add the var
    df["variant"] = df["#CHROM"] + "_" + df.POS.apply(str) + "_" + df.REF + "_" + df.ALT

    return set(df.variant)


def run_freebayes_withoutFiltering(outdir_freebayes, ref, sorted_bam, ploidy, threads, coverage, replace=False):

    # make the dir if not already done
    make_folder(outdir_freebayes)

    #run freebayes
    freebayes_output ="%s/output.raw.vcf"%outdir_freebayes; freebayes_output_tmp = "%s.tmp"%freebayes_output
    if file_is_empty(freebayes_output) or replace is True:
        #print("running freebayes")
        cmd_freebayes = "%s -f %s -p %i --min-coverage %i -b %s --haplotype-length -1 -v %s"%(freebayes, ref, ploidy, coverage, sorted_bam, freebayes_output_tmp); run_cmd(cmd_freebayes)
        os.rename(freebayes_output_tmp, freebayes_output)

    return freebayes_output


def get_fraction_differing_positions_AvsB(snpsA, snpsB, length_genome):

    """Takes two sets of SNPs and returns the fraction of variable positions in the genome"""

    fraction_different_positions = max([len(snpsA.difference(snpsB)), len(snpsB.difference(snpsA))]) / length_genome

    return fraction_different_positions



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



def insert_translocations_into_rearranged_genome_graphBased(reference_genome, input_rearranged_genome, output_rearranged_genome, svDF, translocations_file, svtype_to_svDF, replace=False):

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

   

def generate_rearranged_genome_from_svtype_to_svDF_oldBAdMemory(reference_genome, svtype_to_svDF, outdir, replace=False):

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



def get_svDF_in_coords_of_rearranged_genome_noPreviousTranslocations(svDF, reference_genome, rearranged_genome, svtype, svtype_to_svDF):

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


def simulate_SVs_in_genome_randomTranslocation(reference_genome, mitochondrial_chromosome, outdir, nvars=200, replace=False, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications", "translocations"}, bedpe_breakpoints=None):

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

        # define a bed file with all the regions
        if bedpe_breakpoints is None:
            print_if_verbose("simulating randomly distributed SVs")

            all_regions_bed_df = pd.DataFrame({chrom: {"start":1, "end":chrom_to_len[chrom]} for chrom in chroms}).transpose()
            all_regions_bed_df["chromosome"] = all_regions_bed_df.index
            all_regions_bed_df = all_regions_bed_df[["chromosome", "start", "end"]]
            all_regions_bed = "%s/all_regions_index1.bed"%genome_outdir

        else:
            print_if_verbose("simulating SVs arround %s"%bedpe_breakpoints)
            df_bedpe = pd.read_csv(bedpe_breakpoints, sep="\t")
            df_bedpe.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2']
            bed_1 = df_bedpe[['chrom1', 'start1', 'end1']].rename(columns={"chrom1":"chromosome", "start1":"start", "end1":"end"})
            bed_2 = df_bedpe[['chrom2', 'start2', 'end2']].rename(columns={"chrom2":"chromosome", "start2":"start", "end2":"end"})
            all_regions_bed_df = bed_1.append(bed_2)
            all_regions_bed_df["medium"] = (all_regions_bed_df.start + ((all_regions_bed_df.end - all_regions_bed_df.start)/2) ).apply(int)
            all_regions_bed_df["start"] = all_regions_bed_df.medium
            all_regions_bed_df["end"] = all_regions_bed_df.start + 1
            all_regions_bed_df = all_regions_bed_df[all_regions_bed_df.chromosome.isin(chroms)]

            all_regions_bed = "%s/provided_regions.bed"%genome_outdir

        # write
        all_regions_bed_df[["chromosome", "start", "end"]].to_csv(all_regions_bed, sep="\t", header=False, index=False)

        # simulate random SVs into regions without previous SVs 
        random_sim_dir = "%s/random_SVs"%genome_outdir

        # define the len_shortest_chr
        len_shortest_chr = min([lenChrom for lenChrom in chrom_to_len.values() if lenChrom>=window_l])

        #### GET THE RANDOM INS,INV,DEL ####

        # initialize the cmd
        randomSV_cmd = "%s --input_genome %s --outdir %s --regions_bed %s --len_shortest_chr %i"%(create_random_simulatedSVgenome_R, genome_file, random_sim_dir, all_regions_bed,  len_shortest_chr) 

        # add the numbers of SVs depending on if it is random or not
        if bedpe_breakpoints is None:  expected_svtypes = {"insertions", "deletions", "inversions", "tandemDuplications", "translocations"}.intersection(svtypes)
        else: expected_svtypes = {"translocations"}.intersection(svtypes)

        # if there is only one, get the expecyed SVtypes
        if len(chroms)==1: expected_svtypes = {s for s in expected_svtypes if s!="translocations"}

        # add to cmd        
        svtype_to_arg = {"insertions":"number_Ins", "deletions":"number_Del", "inversions":"number_Inv", "tandemDuplications":"number_Dup", "translocations":"number_Tra"}
        for svtype in expected_svtypes: randomSV_cmd += " --%s %i"%(svtype_to_arg[svtype], vars_to_simulate)

        # define the finalisation file
        random_sims_performed_file = "%s/random_sims_performed.txt"%random_sim_dir

        # run the cmd if necessary
        if file_is_empty(random_sims_performed_file) or replace is True:
            print_if_verbose("generating random SVs")

            # make and delete the folder
            delete_folder(random_sim_dir); make_folder(random_sim_dir)

            # run the random simulation
            #std_rearranging_genome = "%s/simulation_std.txt"%random_sim_dir
            std_rearranging_genome = "stdout"
            print_if_verbose("getting random SVs. The std is in %s"%std_rearranging_genome)
            if std_rearranging_genome!="stdout": run_cmd("%s > %s 2>&1"%(randomSV_cmd, std_rearranging_genome))
            else: run_cmd(randomSV_cmd)
            remove_file(std_rearranging_genome)

            # edit the insertions 
            insertions_file = "%s/insertions.tab"%random_sim_dir
            if not file_is_empty(insertions_file): rewrite_insertions_uniformizedFormat_simulateSV(insertions_file)

            open(random_sims_performed_file, "w").write("random simulations performed")          

        ########################################

        ##### GET THE RANDOM TAN,DEL,INV,INS ARROUND PROVIDED BREAKPOINTS #####
        if bedpe_breakpoints is not None:

            # define the translocations 
            translocations_file = "%s/translocations.tab"%random_sim_dir
            get_SVs_arround_breakpoints_notInTranslocations(genome_file, df_bedpe, vars_to_simulate, translocations_file, svtypes.difference({"translocations"}), replace=replace)

        #######################################################################

        ######### KEEP VARS #########
       
        for svtype in final_svtype_to_svDF.keys():
            svDF = final_svtype_to_svDF[svtype]

            # skip translocations of 1 chromosome
            if len(chroms)==1 and svtype=="translocations": continue

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



def get_SVs_arround_breakpoints_notInTranslocations(genome_file, df_bedpe, nvars, translocations_file, svtypes, replace=False):

    """This function takes a df bedpe and defines svtypes arround these breakpoints."""

    # get the dir
    outdir = get_dir(translocations_file)
    expected_files = {"%s/%s.tab"%(outdir, s) for s in svtypes}

    if any([file_is_empty(f) for f in expected_files]) or replace is True:
        print_if_verbose("calculating random variants among the provided breakpoints")
        # deifine the genome
        chrom_to_len = get_chr_to_len(genome_file)

        # define the chroms
        chroms = set(chrom_to_len)

        # define the interesting bedpe, the one compromised by chroms and reshufle and reindex
        df_bedpe = df_bedpe[(df_bedpe.chrom1.isin(chroms)) & (df_bedpe.chrom2.isin(chroms))].sample(frac=1)
        df_bedpe.index = list(range(0, len(df_bedpe)))

        # define the interval to add for overlaps
        add_interval_bp = 100

        # add things
        df_bedpe["affected_positions_arroundBp"] = df_bedpe.apply(lambda r: get_affected_positions_from_bedpe_r(r, extra_bp=add_interval_bp), axis=1)
        df_bedpe["medium1"] = (df_bedpe.start1 + (df_bedpe.end1 - df_bedpe.start1)/2).apply(int)
        df_bedpe["medium2"] = (df_bedpe.start2 + (df_bedpe.end2 - df_bedpe.start2)/2).apply(int)
        df_bedpe["length"] = df_bedpe.apply(get_length_bedpe_r, axis=1)

        # define the legths of the variants
        lengths_SVs = list(df_bedpe[(~pd.isna(df_bedpe.length)) & (df_bedpe.length>=50)].length)
        random.shuffle(lengths_SVs)

        # get the bed regions affected by translocations (arround the bp). Initialize the already_affected_positions by those that are around the breakpoit
        if not file_is_empty(translocations_file):

            # get the bed
            svDF_tra  = pd.read_csv(translocations_file, sep="\t")
            regions_bed_translocations, nsvs = get_affected_region_bed_for_SVdf(svDF_tra, "translocations", chroms, add_interval_bp=add_interval_bp, first_position_idx=0, translocations_type="breakpoint_pos", chr_to_len=chrom_to_len, insertions_type="only_one_chrB_breakpoint")

            # get the affected positions
            already_affected_positions = get_affected_positions_from_bed_df(regions_bed_translocations)

        else: already_affected_positions = set()

        # already_affected_positions contains a set of "chrom_pos" with the already affected positions. It will be updated with each 

        ############ DEFINE THE RANDOM SVs INVOLVING df_bedpe ############


        # initialize a dict with each svtype and the current number of locations, a
        svtype_to_svDF = {svtype : pd.DataFrame(columns=[x for x in svtype_to_fieldsDict[svtype]["all_fields"] if x!="ID"]) for svtype in svtypes}

        # initialize a set of the svtypes that are already full, meaning that they already have nvars
        already_nvars_svtypes = set()

        # initialize a set that will contain the bad bedpe rows
        bpID_to_tried_svtypes = {bpID : set() for bpID in df_bedpe.index}

        # define the number of breakpoints traversed
        original_n_breakpoints = len(df_bedpe)

        # go through each breakpoint and assign it to a cahegory. Break if a
        while len(df_bedpe)>0:
            print("already traversed %i/%i breakpoints"%(original_n_breakpoints - len(df_bedpe), original_n_breakpoints))

            # interate through each svtype
            for svtype in sorted(svtypes):
                svDF = svtype_to_svDF[svtype]

                # if there are already nvars, skip
                if len(svDF)>=nvars:
                    already_nvars_svtypes.add(svtype)
                    continue

                # get only bedpe rows that are not overlapping the current positions, in terms of breakpoints
                df_bedpe = df_bedpe[(df_bedpe.affected_positions_arroundBp.apply(lambda positions: positions.intersection(already_affected_positions)).apply(len))==0]
                #print("There are %i breakpoints corresponding to translocations"%(original_n_breakpoints - len(df_bedpe)))

                # get only bedpe rows that have not already been tried to assign to all svtypes
                good_bpIDs = {bpID for bpID, tried_svtypes in bpID_to_tried_svtypes.items() if tried_svtypes!=svtypes}.intersection(set(df_bedpe.index))
                df_bedpe = df_bedpe.loc[good_bpIDs]

                # if empty, continue
                if len(df_bedpe)==0: continue

                # get the first bedpe row
                r = df_bedpe.iloc[0]

                # record that this was already tried
                bpID_to_tried_svtypes[r.name].add(svtype)

                # assign the breakpoint to the df_sv_r
                if svtype in {"inversions", "tandemDuplications", "deletions"}:

                    # it should be in the same chromosome
                    if r["chrom1"]!=r["chrom2"]: continue

                    # define the resulting df_sv_r
                    df_sv_r = pd.DataFrame({"Chr":[r["chrom1"]], "Start":[int(r["medium1"])], "End":[int(r["medium2"])]})

                    # get the affected positions of this svtype
                    affected_positions = get_affected_positions_from_bed_df(get_affected_region_bed_for_SVdf(df_sv_r, svtype, {r["chrom1"]}, add_interval_bp=add_interval_bp, first_position_idx=0, translocations_type="breakpoint_pos", chr_to_len=chrom_to_len, insertions_type="only_one_chrB_breakpoint")[0])
                    
                    # only keep if the already affected positions do not overlap this SV
                    if affected_positions.intersection(already_affected_positions)!=set(): continue

                    # keep
                    already_affected_positions.update(affected_positions)
                    svtype_to_svDF[svtype] = svDF.append(df_sv_r)

                elif svtype=="insertions":  

                    # go through different lengths
                    for length_ins in lengths_SVs:

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
                        affected_positions = get_affected_positions_from_bed_df(get_affected_region_bed_for_SVdf(df_sv_r, svtype, {r["chrom1"], r["chrom2"]}, add_interval_bp=add_interval_bp, first_position_idx=0, translocations_type="breakpoint_pos", chr_to_len=chrom_to_len, insertions_type="only_one_chrB_breakpoint")[0])
                    
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



def run_CNV_calling(sorted_bam, reference_genome, outdir, threads, replace, mitochondrial_chromosome, df_bedpe, df_gridss, window_size=100):

    """This function takes a sorted bam and runs several programs on it to get the copy-number variation results. It is important that the sorted_bam contains no duplicates. It will correct bu GC content, mappability and distance to the telomere. All coverage will be corrected by GC content, mappability and the distance to the telomere, which will be calculated also taking into account breakpoint information. """

    make_folder(outdir)


    # define the chrom to len
    chrom_to_len = get_chr_to_len(reference_genome)

    

    # define the breakend locations
    chrom_to_bpPositions = get_chrom_to_bpPositions(df_clove, reference_genome)
    for chrom in chrom_to_bpPositions: chrom_to_bpPositions[chrom].update({0, chrom_to_len[chrom]})



    # define the sample name
    sample_name = get_sample_name_from_bam(sorted_bam)


    thewholeCNVnator_doesntWork


    # run cnvnator
    outdir_cnvnator = "%s/cnvnator_outdir"%outdir
    cnvnator_outfile = run_CNVNATOR(sorted_bam, reference_genome, outdir_cnvnator, threads, replace, mitochondrial_chromosome, df_clove, window_size, sample_name, chrom_to_len)

    ################## RUN CONY ##################


    # go through each chrom
    for c, chrom_len in chrom_to_len.items():

        # define the outdir 
        outdir_chrom = "%s/%s_CONYrun"%(outdir, c)

        # run CONY for this chromosome, this is parallelized
        run_CNV_calling_CONY_one_chromosome(sorted_bam, reference_genome, outdir_chrom, c, replace, window_size, threads, chrom_len, sample_name, chrom_to_bpPositions[c])


    #############################################

    finished_CONYallChroms


    kdajkajdhd




    #run_CNV_calling_CONY_one_chromosome(sorted_bam, reference_genome, outdir, chromosome, replace, window_size, threads):

    ########## PREPARE DATA FOR CONY RUNNING ########## 

    # create an mpileup file for all chromosomes with mosdepth





    available_RAM = get_availableGbRAM()

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



    gridss_run

    # prepare a data table that contains the 1-based positions for each chromosome
    chrom_to_len = get_chr_to_len(reference_genome)
    target_df = pd.DataFrame({c : {"seqname":c, "start":1, "end":length} for c, length in chrom_to_len.items() if c!=mitochondrial_chromosome}).transpose()
    regions_file = "%s/target_regions.bed1based"%outdir
    save_df_as_tab(target_df[["seqname", "start", "end"]], regions_file)

    ###################################################

    # define the final file
    final_file = "%s/CONVY_finished.txt"%outdir

    if file_is_empty(final_file) or replace is True:

        # define the curdir
        CurDir = get_fullpath(".")

        # change the dir to the outdir
        os.chdir(outdir)

        # run CONY
        cony_std = "%s/running_cony.std"%outdir
        print_if_verbose("Running CONY. The std is in %s"%cony_std)

        cmd = "%s --reference_genome %s --sorted_bam %s --regions_file %s --libraries_CONY %s --window_size %i > %s 2>&1"%(run_CONY_R, reference_genome, sorted_bam, regions_file, libraries_CONY, window_size, cony_std)
        run_cmd(cmd, env=EnvName_CONY)

        remove_file(cony_std)

        # return to the initial dir
        os.chdir(CurDir)

        # make the final file
        open(final_file, "w").write("CONY finished")



    #Error in .Call2("solve_user_SEW", refwidths, start, end, width, translate.negative.coord,  : 
    #solving row 14: 'allow.nonnarrowing' is FALSE and the supplied end (20100) is > refwidth
    #Calls: WindowInfo ... narrow -> narrow -> solveUserSEW -> .Call2 -> .Call
    #In addition: There were 42 warnings (use warnings() to see them)



    adkjhkaj


    libraries_CONY


    # run samtools mpileup to get the 
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



def generate_rearranged_genome_from_svtype_to_svDF_tryiningToInsertTranslocationsButFailing(reference_genome, svtype_to_svDF, outdir, replace=False):

    """This function generates a rearranged genome for the provided svDFs, writing their fileds under outdir. It simulates all of them, and then changes the insertions from cut-and-paste to copy-and-paste"""    

    # define the rearranged genome, and generated if not already done
    rearranged_genome = "%s/rearranged_genome.fasta"%outdir
    rearranged_genome_tmp = "%s/rearranged_genome.tmp.fasta"%outdir
    
    # generate the rearranged genome with all variants. Just that all insertions are cut-and-paste
    if file_is_empty(rearranged_genome) or replace is True:
        print_if_verbose("rearranging genome")

        # remove the outdir
        delete_folder(outdir); make_folder(outdir)

        # define the outdir
        outdir_sim = "%s/rearranging_genomeRSVsim"%outdir; 
        delete_folder(outdir_sim); make_folder(outdir_sim)

        # initialize a cmd to create the simulated genome
        targetSV_cmd = "%s --input_genome %s --output_genome %s --output_dir %s"%(create_targeted_simulatedSVgenome_R, reference_genome, rearranged_genome_tmp, outdir_sim)

        for svtype, svDF in svtype_to_svDF.items():

            # shift the insertions by 15 bp so that they are not at the beginning of the chrom
            #if svtype=="insertions": svDF["StartA"] = svDF["StartA"] + 10

            # keep only the interesting svs
            svDF = svDF[svtype_to_fieldsDict[svtype]["all_fields"]]

            # for translocations, write a bed with the breakpoint positions
            if svtype=="translocations":

                translocations_file = "%s/translocations.tab"%outdir
                svDF.to_csv(translocations_file, sep="\t", header=True, index=False)
                rewrite_translocations_uniformizedFormat_simulateSV(translocations_file, reference_genome)
                svDF = pd.read_csv(translocations_file, sep="\t")
                remove_file(translocations_file)


                # get the bed df
                chroms_all = set(svDF.ChrA).union(set(svDF.ChrB))
                add_interval_bp = 100
                chrom_to_len = get_chr_to_len(reference_genome)
                bed_df = get_affected_region_bed_for_SVdf(svDF, "translocations", chroms_all, add_interval_bp=add_interval_bp, first_position_idx=1, translocations_type="breakpoint_pos", chr_to_len=chrom_to_len, insertions_type="only_one_chrB_breakpoint")[0]

                # write
                svfile = "%s/translocations_regions.bed"%outdir
                bed_df.to_csv(svfile, sep="\t", header=False, index=False)

                continue

                # add the number of transloactions
                targetSV_cmd += " --number_Tra %i"%len(svDF)

            else:
                # write file
                svfile = "%s/%s.tab"%(outdir, svtype)
                svDF.to_csv(svfile, sep="\t", header=True, index=False)

            # add the generation of SVs into targetSV_cmd
            targetSV_cmd += " --%s_file %s"%(svtype, svfile)

        # run the cmd
        #std_rearranging_genome = "%s/simulation_std.txt"%outdir
        std_rearranging_genome = "stdout" # debug
        print_if_verbose("rearranging genome. The std is in %s"%std_rearranging_genome)

        if std_rearranging_genome!="stdout": run_cmd("%s > %s 2>&1"%(targetSV_cmd, std_rearranging_genome))
        else: run_cmd(targetSV_cmd)

        sakjsakhj
        
        # transform the cut-and-paste insertions to copy-and-paste, whenever necessary
        insertions_file = "%s/insertions.tab"%outdir
        if not file_is_empty(insertions_file):

            # transform the insertions
            transform_cut_and_paste_to_copy_and_paste_insertions(reference_genome, rearranged_genome_tmp, insertions_file, svtype_to_svDF)

            # edit the insertions so that they are in the correct format
            rewrite_insertions_uniformizedFormat_simulateSV(insertions_file)


        # rewrite translocations
        translocations_file = "%s/translocations.tab"%outdir
        if not file_is_empty(translocations_file): rewrite_translocations_uniformizedFormat_simulateSV(translocations_file, reference_genome)

  
        # write a file that indicates that this has finsihed
        rearranged_genome_finalFile = "%s.performed"%(rearranged_genome)
        open(rearranged_genome_finalFile, "w").write("finsihed")

        # keep
        remove_file(std_rearranging_genome)
        os.rename(rearranged_genome_tmp, rearranged_genome)



def rearrange_genomes_simulateSV_oldVersionBadMemory(reference_genome, outdir, replace=False, nvars=50, mitochondrial_chromosome="mito_C_glabrata_CBS138", simulated_svtype_to_svfile={}, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications"}, check_consistency=False):

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


#### test how well the finding of SVs in an assembly works ####
if opt.testSVgen_from_DefaulReads:

    outdir_test_FindSVinAssembly = "%s/test_FindSVfromDefaultSimulations"%opt.outdir
    if __name__ == '__main__': fun.test_SVgeneration_from_DefaultParms(opt.ref, outdir_test_FindSVinAssembly, sorted_bam, threads=opt.threads, replace=opt.replace, n_simulated_genomes=opt.nsimulations, mitochondrial_chromosome=opt.mitochondrial_chromosome, nvars=opt.nvars)

###############################################################


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

def get_compatible_real_bedpe_breakpoints_old_realSVfiles(close_shortReads_table, reference_genome, outdir, replace=False, threads=4, max_nvars=100, mitochondrial_chromosome="mito_C_glabrata_CBS138", job_array_mode="local", max_ncores_queue=768, time_perSVade_running="48:00:00", queue_jobs="bsc_ls", name_job_array="getRealSVs"):

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



def plot_coverage_across_genome_pairedEndReads(sorted_bam, reference_genome, threads=4, replace=False):

    """Takes a sorted_bam and plots the coverage for windows of the genome"""

    print_if_verbose("plotting coverage across genome")

    # get coverage df  
    calculate_coverage_dir = "%s.calculating_windowcoverage"%sorted_bam; make_folder(calculate_coverage_dir)
    coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, calculate_coverage_dir, sorted_bam, replace=replace, threads=threads), "\t")
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





def get_availableGbRAM_withFiles():

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


def get_simpleSVtype_from_bedpeRow(r):

    """gets the clove notation from a bedpe row"""

    if r["chrom1"]==r["chrom2"]: 

        if r["strand1"]=="+" and r["strand2"]=="+": return "INV1"
        if r["strand1"]=="-" and r["strand2"]=="-": return "INV2"
        if r["strand1"]=="+" and r["strand2"]=="-": return "DEL"
        if r["strand1"]=="-" and r["strand2"]=="+": return "TAN"

    if r["chrom1"]!=r["chrom2"]: 

        if r["strand1"]=="+" and r["strand2"]=="+": return "INVTX1"
        if r["strand1"]=="-" and r["strand2"]=="-": return "INVTX2"
        if r["strand1"]=="+" and r["strand2"]=="-": return "ITX1"
        if r["strand1"]=="-" and r["strand2"]=="+": return "ITX2"

    raise ValueError("It is impossible to establish the clove sig")


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
    coverage_df =  pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, outdir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=True, threads=threads), sep="\t")

    # define stats
    mean_coverage = np.mean(coverage_df.mediancov_1)
    fraction_genome_covered = np.mean(coverage_df.percentcovered_1)/100

    return mean_coverage, fraction_genome_covered



def get_distanceToTelomere_chromosome_GCcontent_to_coverage_fn(df_coverage_train, genome, outdir, mitochondrial_chromosome="mito_C_glabrata_CBS138", replace=False):

    """This function takes a training df_coverage (with windows of a genome) and returns a lambda function that takes GC content, chromosome and  distance to the telomere and returns coverage according to the model.

    This function is edited so that it will not consider the effect of the distance"""
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
    good_chroms = {chrom for chrom, l in chr_to_len.items() if l>=window_l}
    print_if_verbose("there are %i/%i chromsomes above window_l. These will be used to train the GC-content model"%(len(good_chroms), len(chr_to_len)))

    ######## find the coeficients for each chromosome #########

    # map each chromosome to the coefs of the quadratic fit that explains coverage form the distance to the telomere and also the coefs of the GC content explaining the resiudal of this fit
    chrom_to_coefType_to_coefs = {}

    # define the relative coverage of each window of this genome
    median_coverage = np.median(df[(df.coverage>0)  & (df.coverage<10)].coverage)
    df["relative_coverage"] = df.coverage / median_coverage

    # plot the coverage-per-non-0-window distribution
    outfile = "%s/coverage_distribution.pdf"%(outdir)
    if file_is_empty(outfile) or replace is True:

        fig = plt.figure(figsize=(7, 7))
        ax = sns.distplot(df[df.coverage>0].coverage, rug=True)

        ax.set_ylabel("n windows")
        ax.set_xlabel("coverage")
        fig.savefig(outfile, bbox_inches="tight")

    # define the correct regions for modelling. THey have a reasonable coverage
    df_correct = df[(df.relative_coverage<=10) & (df.relative_coverage>0.1) & (df.chromosome.isin(good_chroms))]

    # if the filtering is useless, use all the df. This is a way to skip the modelling.
    if len(df_correct)==0: raise ValueError("There should be some regions for coverage modelling")

    # add the distance to the telomere
    df_correct["coverage_from_dist_to_telomere"] = median_coverage

    # get the residual variation in coverage
    df_correct["residualCoverage_from_dist_to_telomere"] = df_correct.coverage - df_correct.coverage_from_dist_to_telomere

    # get a quadratic fit that predicts coverage from GC content
    coefs_GCcontent = poly.polyfit(df_correct.GCcontent, df_correct.residualCoverage_from_dist_to_telomere, 2)
    df_correct["residualCoverage_from_dist_to_telomere_from_GC_content"] = poly.polyval(df_correct.GCcontent, coefs_GCcontent)

    df_correct["coverage_from_dist_to_telomere_and_GC_content"] = df_correct["coverage_from_dist_to_telomere"] + df_correct["residualCoverage_from_dist_to_telomere_from_GC_content"]

    # get the rsquare of the model
    r2 = r2_score(df_correct.coverage, df_correct.coverage_from_dist_to_telomere_and_GC_content)
    print_if_verbose("The rsquare is %.3f"%(r2))

    # keep coefs
    no_model_coefs =  [median_coverage, 0, 0]
    all_chroms = set(chr_to_len)
    for chrom in all_chroms: chrom_to_coefType_to_coefs[chrom] = {"dist_telomere":no_model_coefs, "GCcontent":coefs_GCcontent}

    # plot
    outfile = "%s/coverage_modelling.pdf"%(outdir)

    if file_is_empty(outfile) or replace is True:

        # define the chroms to plot
        chroms_plot = sorted(good_chroms)
        print_if_verbose("plotting coverage modelling for %i chroms"%len(chroms_plot))

        # plot the coverage for each of the chromosomes
        fig = plt.figure(figsize=(7, len(chroms_plot)*5))

        for I, chrom in enumerate(chroms_plot):

            # initialize a subplot, where each row is one chromosome
            ax = plt.subplot(len(chroms_plot), 1, I+1)

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





# add fake insertion # debug
df_insertions = pd.read_csv("/home/mschikora/samba/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/5478_Candida_glabrata/testing_Accuracy/realSVs/M12_ANI/SVdetection_output/parameter_optimisation/simulation_1/final_simulated_SVs/insertions.tab", sep="\t")
df_insertions["ID"] = "gridss1bf_41o" # this is a fake of the subsampled C. glabrata sample
df_insertions.to_csv("~/Desktop/random_insertions.tab", sep="\t")
svtype_to_svfile["insertions"] = "~/Desktop/random_insertions.tab"




def get_best_coverage_thresholds_knownSVs_df_train(df_train, threads, replace, reference_genome, mitochondrial_chromosome):

    """Try an array of parameters on several training dfs and return a series with the least conservative filters (and the resulting accuracy measures). It will return the filters that are best for both CNV calling and breakpoint identification"""

    # get several parameters of df_train, that are always the same
    r = df_train.iloc[0]
    ploidy = r["ploidy"]
    unique_simID = r["unique_simID"]
    outdir_cnv = r["outdir_cnv"]
    sorted_bam = r["sorted_bam"]

    # define the file
    df_benchmark_CNVcalling_file = "%s/training_df_benchmark_CNVcalling.tab"%outdir_cnv
    df_benchmark_BPident_file = "%s/training_df_benchmark_BPidentification.tab"%outdir_cnv

    if file_is_empty(df_benchmark_CNVcalling_file) or file_is_empty(df_benchmark_BPident_file) or replace is True:
        print_if_verbose("getting the best parameters from training %s. On %i threads"%(unique_simID, threads))

        ########### GET PARAMETERS ###########

        # check that the index is unique
        df_train.index = list(range(0, len(df_train)))

        # make general calculations
        fraction_var, fraction_ref = get_fractions_reads_for_ploidy(ploidy)
        expected_AF = fraction_var

        # define boundaries of the parameters
        max_max_coverage_deletion = 1.1 - expected_AF  
        min_min_coverage_duplication = 1 + expected_AF*0.6

        # define the parms
        n_coverage_filts = 10
        n_r_filts = 10

        # get the combinations of parameters
        inputs_fn_CNVcalling_accuracy = []
        for min_coverage_duplication in np.linspace(min_min_coverage_duplication, 3, n_coverage_filts):
            for max_coverage_deletion in np.linspace(0, max_max_coverage_deletion, n_coverage_filts):

                # get the inputs for the breakpoint identification
                inputs_fn_BPident_accuracy.append([min_coverage_duplication, max_coverage_deletion])


                for min_r_pearson_noFlatRegions in np.linspace(0, 0.5, n_r_filts):
                    for min_r_spearman_noFlatRegions in np.linspace(0, 0.5, n_r_filts):

                        # keep the inputs to optimise
                        inputs_fn_CNVcalling_accuracy.append([df_train, min_coverage_duplication, max_coverage_deletion, min_r_pearson_noFlatRegions, min_r_spearman_noFlatRegions])

        #####################################


        ######## GET ACCURACY BP identification
        if file_is_empty(df_benchmark_BPident_file) or replace is True:

            # define the real chromosome to breakpoint
            df_cnv_known = df_train[df_train.CN!=1]
            known_chrom_to_bpPos = dict(df_cnv_known.groupby("chromosome").apply(lambda df_c: set(df_c.start).union(set(df_c.end))))

            # go through each combination of coverages and return the accuracy of identification of the BPs
            outdir_BPident = "%s/BPidentification"%outdir_cnv
            list_accuracy_measures = [get_accuracy_measures_BPidentification(known_chrom_to_bpPos, sorted_bam, outdir_BPident, reference_genome, threads, replace, min_coverage_duplication, max_coverage_deletion, mitochondrial_chromosome) for min_coverage_duplication, max_coverage_deletion in inputs_fn_BPident_accuracy]

            # get as df
            df_benchmark_BPident = pd.DataFrame(list_accuracy_measures)

            # save
            df_benchmark_BPident_file_tmp = "%s.tmp"%df_benchmark_BPident_file
            df_benchmark_BPident.to_csv(df_benchmark_BPident_file_tmp, sep="\t", index=False, header=True)
            os.rename(df_benchmark_BPident_file_tmp, df_benchmark_BPident_file)

        ######## GET ACCURACY CNV CALLING ########
        if file_is_empty(df_benchmark_CNVcalling_file) or replace is True:

            # get the accuracy measures of CNV calling
            with multiproc.Pool(threads) as pool:
                list_accuracy_measures = pool.starmap(get_accuracy_measures_df_CNV, inputs_fn_CNVcalling_accuracy) 
                    
                pool.close()
                pool.terminate()

            # no parallel
            #list_accuracy_measures = list(map(lambda x: inputs_fn_CNVcalling_accuracy(x[0], x[1], x[2], x[3], x[4]), inputs_fn_accuracy_measures))

            # get as df
            df_benchmark_CNVcalling = pd.DataFrame(list_accuracy_measures)

            # save
            df_benchmark_CNVcalling_file_tmp = "%s.tmp"%df_benchmark_CNVcalling_file
            df_benchmark_CNVcalling.to_csv(df_benchmark_CNVcalling_file_tmp, sep="\t", index=False, header=True)
            os.rename(df_benchmark_CNVcalling_file_tmp, df_benchmark_CNVcalling_file)

        ##########################################

    # load dfs
    df_benchmark_CNVcalling = pd.read_csv(df_benchmark_CNVcalling_file, sep="\t")
    df_benchmark_BPident = pd.read_csv(df_benchmark_BPident_file, sep="\t")

    # get the best combination of filters. The least conservative among CNV calling and BP identification
    list_field_best_fn_CNVcalling = [("Fvalue", max), ("max_coverage_deletion", max), ("min_coverage_duplication", min), ("min_r_pearson_noFlatRegions", max), ("min_r_spearman_noFlatRegions", max)]
    best_filters_series_CNVcalling = get_best_row_in_df_benchmark(df_benchmark_CNVcalling,list_field_best_fn_CNVcalling)

    list_field_best_fn_BPident= [("Fvalue", max), ("max_coverage_deletion", max), ("min_coverage_duplication", min)]
    best_filters_series_BPident = get_best_row_in_df_benchmark(df_benchmark_BPident, list_field_best_fn_BPident)


    ######### GET INTEGRATED BEST FILTERS #########

    # get floats
    for f in ["min_coverage_duplication", "max_coverage_deletion"]:
        df_benchmark_CNVcalling[f] =  df_benchmark_CNVcalling[f].apply(float)
        df_benchmark_BPident[f] =  df_benchmark_BPident[f].apply(float)

    # get the best_filters_series.
    final_min_r_pearson_noFlatRegions = best_filters_series_CNVcalling["min_r_pearson_noFlatRegions"]
    final_min_r_spearman_noFlatRegions = best_filters_series_CNVcalling["min_r_spearman_noFlatRegions"]

    """
    mean_min_coverage_duplication = float(np.mean([best_filters_series_CNVcalling["min_coverage_duplication"], best_filters_series_BPident["min_coverage_duplication"]]))
    final_min_coverage_duplication = find_nearest(df_benchmark_CNVcalling.min_coverage_duplication, mean_min_coverage_duplication)

    mean_max_coverage_deletion = float(np.mean([best_filters_series_CNVcalling["max_coverage_deletion"], best_filters_series_BPident["max_coverage_deletion"]]))
    final_max_coverage_deletion = find_nearest(df_benchmark_CNVcalling.max_coverage_deletion, mean_max_coverage_deletion)
    """

    # get the least conservative calls
    final_min_coverage_duplication = float(min([best_filters_series_CNVcalling["min_coverage_duplication"], best_filters_series_BPident["min_coverage_duplication"]]))

    final_max_coverage_deletion = float(max([best_filters_series_CNVcalling["max_coverage_deletion"], best_filters_series_BPident["max_coverage_deletion"]]))

    # get the best filters
    r_final_CNVcalling = df_benchmark_CNVcalling[(df_benchmark_CNVcalling.min_r_pearson_noFlatRegions==final_min_r_pearson_noFlatRegions) & (df_benchmark_CNVcalling.min_r_spearman_noFlatRegions==final_min_r_spearman_noFlatRegions) & (df_benchmark_CNVcalling.min_coverage_duplication==final_min_coverage_duplication) & (df_benchmark_CNVcalling.max_coverage_deletion==final_max_coverage_deletion)].iloc[0]

    r_final_BPident = df_benchmark_BPident[(df_benchmark_BPident.min_coverage_duplication==final_min_coverage_duplication) & (df_benchmark_BPident.max_coverage_deletion==final_max_coverage_deletion)].iloc[0]


    best_filters_series = pd.Series({"min_r_pearson_noFlatRegions":final_min_r_pearson_noFlatRegions,
                                     "min_r_spearman_noFlatRegions":final_min_r_spearman_noFlatRegions,
                                     "min_coverage_duplication":final_min_coverage_duplication,
                                     "max_coverage_deletion":final_max_coverage_deletion,
                                     "CNVcalling_precision":r_final_CNVcalling["precision"],
                                     "CNVcalling_recall":r_final_CNVcalling["recall"],
                                     "CNVcalling_Fvalue":r_final_CNVcalling["Fvalue"],
                                     "BPident_precision":r_final_BPident["precision"],
                                     "BPident_recall":r_final_BPident["recall"],
                                     "BPident_Fvalue":r_final_BPident["Fvalue"]
                                    })


    #############################################

    return best_filters_series




def get_accuracy_measures_BPidentification(known_chrom_to_bpPos, sorted_bam, outdir, reference_genome, threads, replace, min_coverage_duplication, max_coverage_deletion, mitochondrial_chromosome):

    """This function takes a series that maps each chromosome to the known Bp positions. It runs get_chrom_to_bpPos_inferredByCoverage and returns the accuracy measurements of the given filters. The outdir should be the same for all"""

    print_if_verbose("min_coverage_duplication", min_coverage_duplication, "max_coverage_deletion", max_coverage_deletion)

    # deepcopy
    known_chrom_to_bpPos = cp.deepcopy(known_chrom_to_bpPos)

    # make the outdir
    make_folder(outdir)

    # get the one inferred by coverage
    chrom_to_bpPos_inferredByCoverage, window_size = get_chrom_to_bpPos_inferredByCoverage(sorted_bam, outdir, reference_genome, threads=threads, replace=replace, min_coverage_duplication=min_coverage_duplication, max_coverage_deletion=max_coverage_deletion, mitochondrial_chromosome=mitochondrial_chromosome)

    # initialize the accuracy counters, with those chromosomes that are missing
    FP_chromosomes = set(chrom_to_bpPos_inferredByCoverage).difference(set(known_chrom_to_bpPos))
    FN_chromosomes = set(known_chrom_to_bpPos).difference(set(chrom_to_bpPos_inferredByCoverage))
    FP = sum({len(chrom_to_bpPos_inferredByCoverage[chrom]) for chrom in FP_chromosomes})
    FN = sum({len(known_chrom_to_bpPos[chrom]) for chrom in FN_chromosomes})
    TP = 0

    # go through each chromosome
    for chrom, positions_known in known_chrom_to_bpPos.items():

        # get the inferred positions
        if chrom in chrom_to_bpPos_inferredByCoverage: inferred_positions = chrom_to_bpPos_inferredByCoverage[chrom]
        else: inferred_positions = set()

        # map each position to the inferred ones
        known_to_approx_positions = pd.Series({known_pos : {inf_pos for inf_pos in inferred_positions if abs(inf_pos-known_pos)<=window_size} for known_pos in positions_known})

        # define the matched inferred_positions
        matched_inferred_positions = set.union(*known_to_approx_positions)
        unmatched_inferred_positions = inferred_positions.difference(matched_inferred_positions)

        # define the unmatched known positions
        unmatched_known_positions = set(known_to_approx_positions[known_to_approx_positions.apply(len)==0].index)

        # debug
        if len(unmatched_known_positions.intersection(positions_known))!=len(unmatched_known_positions): raise ValueError("something went wrong with new positions")
        if len(matched_inferred_positions.intersection(inferred_positions))!=len(matched_inferred_positions): raise ValueError("something went wrong with inferred_positions")

        # keep 
        TP += len(matched_inferred_positions)
        FP += len(unmatched_inferred_positions)
        FN += len(unmatched_known_positions)

    # get accuracy
    nevents = TP + FP + FN
    precision, recall, Fvalue = get_precision_recall_and_Fvalue(nevents, TP, FP, FN)

    return pd.Series({"min_coverage_duplication":min_coverage_duplication, "max_coverage_deletion":max_coverage_deletion, "precision":precision, "recall":recall, "Fvalue":Fvalue})



############ overlaps repeats 





    # get the variant IDs that are in repeats
    bed_fields = ["#chrom", "start", "end", "length"]


    


   
    bedtools intersect

    print(vcf_df, vcf_df.keys())





    adkjhakjadha

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



def get_chrom_to_bpPos_inferredByCoverage(sorted_bam, outdir, reference_genome, threads=4, replace=False, fraction_window_l=0.005, min_coverage_duplication=1.8, max_coverage_deletion=0.01, mitochondrial_chromosome="mito_C_glabrata_CBS138"):

    """This function calculates the coverage for windows of a fraction of the window_l, and returns a dict that maps each chromosome to potential breakpoint. """

    # define a window size as a fraction of window_l
    window_size = max([int(window_l*fraction_window_l), 100])

    # get the windows
    windows_file = "%s/genome_windows_%ibp.bed"%(outdir, window_size)
    if file_is_empty(windows_file) or replace is True:

        windows_file_stderr = "%s.generating.stderr"%windows_file
        windows_file_tmp = "%s.tmp"%windows_file

        print_if_verbose("running make windows. The stderr is in %s"%windows_file_stderr)
        run_cmd("%s makewindows -g %s.fai -w %i > %s 2>%s"%(bedtools, reference_genome, window_size, windows_file_tmp, windows_file_stderr))

        # add the header
        df = pd.read_csv(windows_file_tmp, sep="\t", header=None, names=["chromosome", "start", "end"])
        df.to_csv(windows_file_tmp, sep="\t", index=False, header=True)

        remove_file(windows_file_stderr)
        os.rename(windows_file_tmp, windows_file)

    # get the coverage_df of the small regions
    coverage_df = get_coverage_per_window_df_without_repeating(reference_genome, sorted_bam, windows_file, replace=replace, run_in_parallel=True, delete_bams=False, threads=threads)
    coverage_df["relative_coverage"] = coverage_df.mediancov_1/get_median_coverage(coverage_df, mitochondrial_chromosome)
    coverage_df["medium"] = coverage_df.start + ((coverage_df.end - coverage_df.start)/2).apply(int)

    # add the coverage of this window
    cov_previous_window = np.array(coverage_df.relative_coverage.iloc[0:-1])
    cov_previous_window = np.append([coverage_df.relative_coverage.iloc[0]], cov_previous_window)
    coverage_df["cov_previous_window"] =  cov_previous_window

    # add the coverage of two regions before
    cov_previous_2windows = np.array(coverage_df.relative_coverage.iloc[0:-2])
    cov_previous_2windows = np.append([coverage_df.relative_coverage.iloc[0]]*2, cov_previous_2windows)
    coverage_df["cov_previous_2windows"] =  cov_previous_2windows

    # add the coverage in the next window
    cov_next_window = np.array(coverage_df.relative_coverage.iloc[1:])
    cov_next_window = np.append(cov_next_window, [coverage_df.relative_coverage.iloc[-1]])
    coverage_df["cov_next_window"] =  cov_next_window

    # add whether this region is a new deletion
    coverage_df["potential_new_deletion"] = (coverage_df.relative_coverage<max_coverage_deletion) & (coverage_df.cov_next_window<max_coverage_deletion) & (coverage_df.cov_previous_window>=max_coverage_deletion) & (coverage_df.cov_previous_2windows>=max_coverage_deletion)

    # add whether this region is a new duplication
    coverage_df["potential_new_duplication"] = (coverage_df.relative_coverage>min_coverage_duplication) & (coverage_df.cov_next_window>min_coverage_duplication) & (coverage_df.cov_previous_window<=min_coverage_duplication) & (coverage_df.cov_previous_2windows<=min_coverage_duplication)

    # add whether this is a potential loss-of-deletion
    coverage_df["potential_new_loss_deletion"] = (coverage_df.relative_coverage>max_coverage_deletion) & (coverage_df.cov_next_window>max_coverage_deletion) & (coverage_df.cov_previous_window<=max_coverage_deletion) & (coverage_df.cov_previous_2windows<=max_coverage_deletion)

    # add whether this is a potential loss-of-duplication
    coverage_df["potential_new_loss_duplication"] = (coverage_df.relative_coverage<min_coverage_duplication) & (coverage_df.cov_next_window<min_coverage_duplication) & (coverage_df.cov_previous_window>=min_coverage_duplication) & (coverage_df.cov_previous_2windows>=min_coverage_duplication)

    # keep those regions that have a change in coverage outside the range
    coverage_df = coverage_df[(coverage_df.potential_new_duplication) | (coverage_df.potential_new_deletion) | (coverage_df.potential_new_loss_deletion) | (coverage_df.potential_new_loss_duplication)]

    # get as dict
    chrom_to_breakpoints = dict(coverage_df.groupby("chromosome").apply(lambda df_c: set(df_c.medium)))

    return chrom_to_breakpoints, window_size




def get_CNV_calling_with_coverageBased_added(df_vcf, outdir, sorted_bam, reference_genome, df_clove, df_gridss, replace=False, threads=4, mitochondrial_chromosome="mito_C_glabrata_CBS138", max_coverage_deletion=0.01, min_coverage_duplication=1.8, min_sv_size=50, min_r_pearson_noFlatRegions=0.2, min_r_spearman_noFlatRegions=0.2):

    """This function gets a df_vcf with calls of CNV (<DUP> or <DEL>) based on coverage. It generates a bed with the query regions, which include genes and any region that is surrounded by breakpoints (this also includes whole chromosomes). """

    print_if_verbose("calculating CNV based on coverage and breakpoints")

    # keep objects
    df_clove = cp.deepcopy(df_clove)
    df_vcf = cp.deepcopy(df_vcf)

    # add the IDs set
    df_clove["IDs_set"] = df_clove.ID.apply(lambda x: set(re.split("\+|\-", x)))
    check_that_cloveIDs_are_in_gridss_df(df_clove, df_gridss)

    # make outdir
    make_folder(outdir)

    # define files
    raw_cnv_file = "%s/all_CNV_events.tab"%outdir
    merged_vcf_file = "%s/merged_SV_and_CNV.tab"%outdir

    if file_is_empty(merged_vcf_file) or replace is  True:

        # get the chrom to len
        chrom_to_len = get_chr_to_len(reference_genome)

        # calculate median cov
        destination_dir = "%s.calculating_windowcoverage"%sorted_bam
        coverage_df = pd.read_csv(generate_coverage_per_window_file_parallel(reference_genome, destination_dir, sorted_bam, windows_file="none", replace=replace, run_in_parallel=True, delete_bams=False), sep="\t")
        median_coverage = get_median_coverage(coverage_df, mitochondrial_chromosome)

        if file_is_empty(raw_cnv_file) or replace is True:
            print_if_verbose("getting raw CNV calls")

            ######### DEFINE THE QUERY REGIONS #########

            # define the filename
            query_df_filename = "%s/query_regions_CNV.bed"%outdir
            if file_is_empty(query_df_filename) or replace is True:
                print_if_verbose("calculating potential regions")

                # initialize all IDs
                query_regions_dict = {}

                # get the chrom_to_bpPositions, adding the 0 and end 
                chrom_to_bpPositions_real = get_chrom_to_bpPositions(df_clove, reference_genome)
                for chrom in chrom_to_bpPositions_real: chrom_to_bpPositions_real[chrom].update({0, chrom_to_len[chrom]})

                # get the breakpoint positions that are called by gridss and not by CLOVE
                #chrom_to_bpPositions_filteredOut = get_chrom_to_bpPositions_filteredOut_breakpoints(df_gridss, df_clove)





                # get the one inferred by coverage
                chrom_to_bpPos_inferredByCoverage, window_size = get_chrom_to_bpPos_inferredByCoverage(sorted_bam, outdir, reference_genome, threads=threads, replace=replace, min_coverage_duplication=min_coverage_duplication, max_coverage_deletion=max_coverage_deletion, mitochondrial_chromosome=mitochondrial_chromosome)

                # go through each chromose, keeping the IDs that are together
                for chrom in chrom_to_len:
                    print_if_verbose(chrom)

                    # get the two sets
                    if chrom in chrom_to_bpPositions_real: bpPositions_real = chrom_to_bpPositions_real[chrom]
                    else: bpPositions_real = set()

                    if chrom in chrom_to_bpPos_inferredByCoverage: bpPositions_inferred = chrom_to_bpPos_inferredByCoverage[chrom]
                    else: bpPositions_inferred = set()
                    #bpPositions_inferred = set() # debug

                    # get as a np.array
                    bpPositions_set = bpPositions_real.union(bpPositions_inferred)
                    bpPositions = np.array(sorted(bpPositions_set))

                    # go through each combination of them
                    for bpPosA in bpPositions:
                        for bpPosB in bpPositions[bpPositions>bpPosA]:

                            # define the length
                            len_region = bpPosB - bpPosA

                            # skip if the length is less than min_sv_size
                            if len_region<min_sv_size: continue

                            # define the region_type
                            if bpPosA==0 and bpPosB==chrom_to_len[chrom]: region_type = "wholeChrom"
                            elif bpPosA in bpPositions_real and bpPosB in bpPositions_real: region_type = "RealBPs"
                            elif (bpPosA in bpPositions_real and bpPosB in bpPositions_inferred) or (bpPosB in bpPositions_real and bpPosA in bpPositions_inferred): region_type = "oneRealBP"
                            elif bpPosA in bpPositions_inferred and bpPosB in bpPositions_inferred: region_type = "InferredBPs"
                            else: raise ValueError("there is an error with the breakpoints")

                            # define the ID of the region
                            ID = "%s:%i-%i|%s"%(chrom, bpPosA, bpPosB, region_type)
                            query_regions_dict[ID] = {"chromosome":chrom, "start":bpPosA, "end":bpPosB, "ID":ID}

                # get the df
                print_if_verbose("gettting df")
                query_df = pd.DataFrame(query_regions_dict).transpose()[["chromosome", "start", "end", "ID"]].sort_values(by=["chromosome", "start", "end"])
                query_df = query_df.drop_duplicates()

                # save
                print_if_verbose("saving")
                query_df_filename_tmp = "%s.tmp"%query_df_filename
                query_df.to_csv(query_df_filename_tmp, sep="\t", header=True, index=False)
                os.rename(query_df_filename_tmp, query_df_filename)

            # load
            query_df = pd.read_csv(query_df_filename, sep="\t")
            total_n_regions = len(query_df)
            print_if_verbose("There are %i regions to test CNV for"%(total_n_regions))

            # debug 
            #query_df = query_df.iloc[0:100]

            ############################################

            ######### GET THE REGIONS UNDER CNV #########
            print_if_verbose("getting actual regions under CNV")

            # add the coverage of neighbor regions to query_df
            bed_windows_prefix = "%s/calculating_cov_neighbors_queryRegions"%outdir
            query_df = get_df_with_coverage_per_windows_relative_to_neighbor_regions(query_df, bed_windows_prefix, reference_genome, sorted_bam, df_clove, median_coverage, replace=replace, run_in_parallel=True, delete_bams=False, threads=threads)
            # remember that the coverage relative to neighbors is meaningless for whole chromosomes (or those vars that cover >90% of the chromosome)


            # just keep the regions where the coverage is according
            query_df = query_df[(query_df.relative_coverage_target>=min_coverage_duplication) | (query_df.relative_coverage_target<=max_coverage_deletion)]
            print_if_verbose("There are %i/%i that might have CNV. Testing with CI-95pct..."%(len(query_df), total_n_regions))

            # get a query_df with the coverage statistics of windows within them
            query_df = get_coverage_df_windows_with_within_windows_statistics(query_df, outdir, sorted_bam, reference_genome, median_coverage, replace=replace, threads=threads)

            # save
            all_potential_regions_file = "%s/all_potential_regions_CNV.tab"%outdir 
            query_df.to_csv(all_potential_regions_file, sep="\t", index=False, header=True)

            # the notFlat coverage is indicated by spearman r and pearson r
            query_df["notFlatCoverage"] = ((query_df.spearman_p<0.05) | (query_df.pearson_p<0.05)) & (query_df.abs_pearson_r>=min_r_pearson_noFlatRegions) & (query_df.abs_spearman_r>=min_r_spearman_noFlatRegions)

            query_df["is_duplication"] = (query_df.mean95CI_lower_rel_coverage_relative>=min_coverage_duplication) & (query_df.median95CI_lower_rel_coverage_relative>=min_coverage_duplication) & (query_df.mean95CI_lower_rel_coverage>=min_coverage_duplication) & (query_df.median95CI_lower_rel_coverage>=min_coverage_duplication)
            
            query_df["is_deletion"] = (query_df.mean95CI_higher_rel_coverage_relative<=max_coverage_deletion) & (query_df.median95CI_higher_rel_coverage_relative<=max_coverage_deletion) & (query_df.mean95CI_higher_rel_coverage<=max_coverage_deletion) & (query_df.median95CI_higher_rel_coverage<=max_coverage_deletion)

            # prints
            print_if_verbose("There are %i/%i regions with signs of duplication"%(sum(query_df.is_duplication), len(query_df)))
            print_if_verbose("There are %i/%i regions with signs of deletion"%(sum(query_df.is_deletion), len(query_df)))

            # save
            df_CNV = query_df[~(query_df.notFlatCoverage) & ((query_df.is_duplication) | (query_df.is_deletion))]

            raw_cnv_file_tmp = "%s.tmp"%raw_cnv_file
            df_CNV.to_csv(raw_cnv_file_tmp, sep="\t", index=False, header=True)
            os.rename(raw_cnv_file_tmp, raw_cnv_file)

            #################################################

        ############ GET RID OF REDUNDANT EVENTS ############

        # load df
        df_CNV = pd.read_csv(raw_cnv_file, sep="\t")

        # add the SVTYPE
        bool_to_text = {True:"DUP", False:"DEL"}
        df_CNV["SVTYPE"] = df_CNV.is_duplication.apply(lambda x: bool_to_text[x])

        # get the vcf formatted as df_CNV
        df_vcf_forCNV = df_vcf[df_vcf.ALT.isin({"<DUP>", "<TDUP>", "<DEL>"})].rename(columns={"POS":"start", "#CHROM":"chromosome"}).set_index("ID", drop=False)
        df_vcf_forCNV["end"] = df_vcf_forCNV.INFO.apply(lambda x: [int(y.split("END=")[1]) for y in x.split(";") if y.startswith("END")][0])

        # add the svtype
        svtype_to_DUPDEL = {"TDUP":"DUP", "DUP":"DUP", "DEL":"DEL"}
        df_vcf_forCNV["SVTYPE"] = df_vcf_forCNV.INFO.apply(lambda x: [svtype_to_DUPDEL[y.split("SVTYPE=")[1]] for y in x.split(";") if y.startswith("SVTYPE")][0])

        # add the type of SVcall
        df_vcf_forCNV["type_CNVcall"] = "gridssClove"
        df_CNV["type_CNVcall"] = "coverage"

        # get the type of BPs
        df_vcf_forCNV["typeBPs"] = "RealBPs"
        df_CNV["typeBPs"] = df_CNV.ID.apply(lambda x: x.split("|")[-1])

        # define all called SVs
        fields = ["ID", "chromosome", "start", "end", "SVTYPE", "type_CNVcall", "typeBPs"]
        all_df_CNV = df_CNV[fields].append(df_vcf_forCNV[fields]).set_index("ID")

        # get only those CNVs that don't overlap other CNVs
        df_CNV["overlaps_otherSVs"] = df_CNV.apply(lambda r: is_overlapping_other_SVs_CNVdf_r(r, all_df_CNV), axis=1)
        df_CNV = df_CNV[~(df_CNV.overlaps_otherSVs)]

        nDups = sum(df_CNV.SVTYPE=="DUP")
        nDels = sum(df_CNV.SVTYPE=="DEL")
        print_if_verbose("There are %i DUPs and %i DELs that do not ovelap any variants already inferred by GRIDSS+CLOVE"%(nDups, nDels))

        #####################################################

        ###### FORMAT AS VCF ######

        # define fields
        vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        data_fields = ["chromosome", "start", "end", "ID", "SVTYPE", "INFO", "median95CI_lower_rel_coverage", "median95CI_higher_rel_coverage", "median95CI_lower_rel_coverage_relative", "median95CI_higher_rel_coverage_relative", "abs_spearman_r", "abs_pearson_r", "spearman_p", "pearson_p"]

        # get the coverage calculation for the input vcf TAN,DUP,DEL

        if len(df_vcf_forCNV)==0: df_vcf_forCNV_final = pd.DataFrame(columns=data_fields)
        
        else:   

            df_vcf_forCNV_final  = df_vcf_forCNV.set_index("ID", drop=False)
            bed_windows_prefix = "%s/calculating_cov_neighbors_SV-based_vcf"%outdir
            df_vcf_forCNV_final = get_df_with_coverage_per_windows_relative_to_neighbor_regions(df_vcf_forCNV_final, bed_windows_prefix, reference_genome, sorted_bam, df_clove, median_coverage, replace=replace, run_in_parallel=True, delete_bams=False, threads=threads)
            df_vcf_forCNV_final = get_coverage_df_windows_with_within_windows_statistics(df_vcf_forCNV_final, outdir, sorted_bam, reference_genome, median_coverage, replace=replace, threads=threads)

            # change the SVTYPE to follow INFO. This is important to get TDUPs back in place
            df_vcf_forCNV_final["SVTYPE"] = df_vcf_forCNV_final.INFO.apply(lambda x: [y.split("SVTYPE=")[1] for y in x.split(";") if y.startswith("SVTYPE")][0])

        # initialize the final vcf
        if len(df_CNV)==0: df_CNV = pd.DataFrame(columns=data_fields)

        else:

            # add the fields
            df_CNV["region_type"] = df_CNV.ID.apply(lambda x: x.split("|")[1])
            df_CNV["coords_string"] = df_CNV.ID.apply(lambda x: x.split("|")[0])
            df_CNV["ID"] = "coverage" + df_CNV.SVTYPE +  "|" + df_CNV.coords_string
            df_CNV["INFO"] = df_CNV.apply(lambda r: "END=%i;SVTYPE=%s;BPS_TYPE=%s"%(r["end"], r["SVTYPE"], r["region_type"]), axis=1)

        # initialize merge with the info fields
        df_vcf_final = df_CNV[data_fields].append(df_vcf_forCNV_final[data_fields])

        # add the INFO
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

        # delete bams        
        clean_chromosomal_bam_files(sorted_bam, reference_genome)
     
        # write
        merged_vcf_file_tmp = "%s.tmp"%merged_vcf_file
        df_vcf_final.to_csv(merged_vcf_file_tmp, sep="\t", index=False, header=True)
        os.rename(merged_vcf_file_tmp, merged_vcf_file)

    # load
    df_vcf_final = pd.read_csv(merged_vcf_file, sep="\t")

    return  df_vcf_final



def is_overlapping_other_SVs_CNVdf_r(r, df_CNV_all):

    """This function takes a row of df_CNV and all the CNVs and returns whether the row is below another r in df_CNV_all"""

    # get df with the other CNVs
    df_CNV_test = df_CNV_all.loc[set(df_CNV_all.index).difference({r["ID"]})]
    if (len(df_CNV_all)-len(df_CNV_test))!=1: raise ValueError("r.ID should be in the df_CNV_all")


    # define the df with the same chromosome and SVTYPE
    df_same_chr = df_CNV_test[(df_CNV_test.SVTYPE==r["SVTYPE"]) & (df_CNV_test.chromosome==r["chromosome"])]

    # return false if it is from the same chromosome
    if len(df_same_chr)==0: return False

    # define the df with the SVs that incude this SVtype
    df_overlapping_all = df_same_chr[(df_same_chr.start<=r["start"]) & (df_same_chr.end>=r["end"])]

    # define the df that overlaps fraction_overlap (0.8)
    min_fraction_overlap = 0.8
    equal_fields = ["chromosome"]
    approximate_fields = ["start", "end"]
    chromField_to_posFields = {"chromosome":{"start":"start", "end":"end"}}
    df_overlapping_fraction = df_same_chr[df_same_chr.apply(lambda rtarget: get_is_matching_predicted_and_known_rows(rtarget, r, equal_fields, approximate_fields, chromField_to_posFields, tol_bp=window_l*10, pct_overlap=min_fraction_overlap), axis=1)]

    # define dfs that are overlaping
    df_overlapping_fraction_RealBPs = df_overlapping_fraction[df_overlapping_fraction.typeBPs=="RealBPs"]
    df_overlapping_all_RealBPs = df_overlapping_all[df_overlapping_all.typeBPs=="RealBPs"]

    # if these are high-quality BPs, the only overlapping dfs are those that are also high-quality and involve it
    if r["typeBPs"]=="RealBPs": df_CNV_overlapping = df_overlapping_all_RealBPs

    # if they are low quality, we'll take as overlapping regions those that are embedding them and those real variants that are partially overlapping them 
    elif r["typeBPs"] in {"InferredBPs", "oneRealBP", "wholeChrom"}: df_CNV_overlapping = df_overlapping_all.append(df_overlapping_fraction_RealBPs)

    else: raise ValueError("%s is not a valid type_BPs"%r["typeBPs"])

    # return if it is 0
    if len(df_CNV_overlapping)==0: return False
    else: return True



def get_chrom_to_bpPositions_filteredOut_breakpoints(df_gridss, df_clove, reference_genome):

    """This function returns the mapping between each chromosome and the df_gridss breakpoint positions that are not in clove"""

    # intit breakopints positions with all chroms
    all_chromosomes = set(get_chr_to_len(reference_genome))
    chrom_to_bpPositions = {c:set() for c in all_chromosomes}

    # filter df_gridss
    all_cloveIDs = set.union(*df_clove.IDs_set)
    df_gridss_filt = df_gridss[~(df_gridss.eventID_as_clove.isin(all_cloveIDs))]

    # debug 
    if len(df_gridss_filt)==0: return chrom_to_bpPositions

    # get the bpPOs
    df_gridss_filt["POS_0based"] = df_gridss_filt.POS - 1
    chrom_to_bpPositions_df_gridss = dict(df_gridss_filt.groupby("#CHROM").apply(lambda df_c: set(df_c.POS_0based)))

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

        data_dict[region] = {"region":"%s_region"%region, "start":start, "end":end, "chromosome":r["#CHROM"], "ID":r["ID"]}


    return pd.DataFrame(data_dict).transpose()


def simulate_SVs_in_genome(reference_genome, mitochondrial_chromosome, outdir, nvars=200, replace=False, svtypes={"insertions", "deletions", "inversions", "translocations", "tandemDuplications", "translocations"}, bedpe_breakpoints=None, threads=4):

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
       
        # define the df_bedpe arround which to place the breakpoints
        if bedpe_breakpoints is None:

            # define the expected SV types
            expected_svtypes = {"insertions", "deletions", "inversions", "tandemDuplications", "translocations", "translocations"}.intersection(svtypes)

            # if there is only one, get the expecyed SVtypes
            if len(chroms)==1: expected_svtypes = {s for s in expected_svtypes if s!="translocations"}

            # keep trying to generate n_BPs until you have enough intrachromosomal and interchromosomal breakpoints
            factor_n_BPs = 10
            min_n_BPs = vars_to_simulate * len(expected_svtypes) * factor_n_BPs
            previous_len_df_bedpe = -1

            print_if_verbose("looking for at least %i of each type"%min_n_BPs)

            # enter a loop in which you keep adding 
            while True:

                # define the factor as twice the previous
                factor_n_BPs = factor_n_BPs*3

                # get the number of breakpoints related to factor_n_BPs
                n_BPs = vars_to_simulate * len(expected_svtypes) * factor_n_BPs
                df_bedpe = get_random_df_bedpe(n_BPs, genome_file, replace=replace, threads=threads)

                # define the number of each typeBP
                n_inter_chromosomal = sum(df_bedpe.chrom1!=df_bedpe.chrom2)
                n_intra_chromosomal = sum(df_bedpe.chrom1==df_bedpe.chrom2)

                print_if_verbose("There are %i/%i inter/intra chromosomal events"%(n_inter_chromosomal, n_intra_chromosomal))

                # once you have enough BPs of each type, break
                if "translocations" in expected_svtypes and n_inter_chromosomal>=min_n_BPs and n_intra_chromosomal>=min_n_BPs: break

                elif "translocations" not in expected_svtypes and n_intra_chromosomal>=min_n_BPs: break

                # if you are not finding any new breakpoints
                if len(df_bedpe)<=previous_len_df_bedpe: break

                # keep previous_len_df_bedpe
                previous_len_df_bedpe = len(df_bedpe)

           
        else:

            # load the bedpe
            bedpe_fields = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]
            df_bedpe = pd.read_csv(bedpe_breakpoints, sep="\t", header=-1, names=bedpe_fields)


        # define the SVtypes arround the df_bedpe 
        get_SVs_arround_breakpoints(genome_file, df_bedpe, vars_to_simulate, random_sim_dir, svtypes, replace=replace)

        ######### KEEP VARS #########
       
        for svtype in final_svtype_to_svDF.keys():
            svDF = final_svtype_to_svDF[svtype]

            # skip translocations of 1 chromosome
            if len(chroms)==1 and svtype=="translocations": continue

            # get the new sv
            file = "%s/%s.tab"%(random_sim_dir, svtype)
            new_svDF = pd.read_csv(file, sep="\t")
            new_svDF = new_svDF[[c for c in new_svDF.keys() if "BpSeq" not in c]]

            # add the name
            new_svDF["ID"] = new_svDF.Name + "_sim_%s"%type_genome

            # append 
            final_svtype_to_svDF[svtype] = svDF.append(new_svDF, sort=True)

        ##############################


    kadjhgkhajh

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


def get_bedpe_series_brekends(Ileft, Iright, Ibp, df_windows):

    """Takes a df_windows and returns a bedpe dict for windows Ileft and Iright"""

    # map bools to strand
    bool_to_strand = {True:"+", False:"-"}

    # define the positions
    r1 = df_windows.loc[Ileft]
    r2 = df_windows.loc[Iright]

    # get the dict
    data_dict =  {"chrom1" : r1["chromosome"], 
                  "start1" : r1["pos"],
                  "end1" : r1["pos"]+5,

                  "chrom2" : r2["chromosome"], 
                  "start2" : r2["pos"],
                  "end2" : r2["pos"]+5,

                  "name": "bp_%i"%Ibp,
                  "score": 1000,

                  "strand1": bool_to_strand[random.uniform(0, 1)>=0.5],
                  "strand2": bool_to_strand[random.uniform(0, 1)>=0.5]}

    return pd.Series(data_dict)



def get_random_df_bedpe(n_BPs, reference_genome, replace=False, threads=4):

    """This function generates random bBPs breakpoints for a reference genome, and returns a bedpe with the reference genome"""

    print_if_verbose("getting random df_bedpe")

    ###### GET GENOME WINDOWS INTO DF ######

    # index the genome
    index_genome(reference_genome, replace=replace)

    # get df_windows of 1000 bp
    chrom_to_len = get_chr_to_len(reference_genome)
    length_genome = sum(chrom_to_len.values())
    window_size = max([200, int(length_genome/n_BPs)])
    windows_file = "%s.windows%ibp.bed"%(reference_genome, window_size)
    windows_file_stderr = "%s.generating.stderr"%windows_file

    if file_is_empty(windows_file) or replace is True:

        windows_file_tmp = "%s.tmp"%windows_file

        print_if_verbose("generating windows_file. The stderr is in %s"%windows_file_stderr)
        run_cmd("%s makewindows -g %s.fai -w %i > %s 2>%s"%(bedtools, reference_genome, window_size, windows_file_tmp, windows_file_stderr))
        remove_file(windows_file_stderr)

        os.rename(windows_file_tmp, windows_file)

    # get df, the index is already unique
    df_windows = pd.read_csv(windows_file, sep="\t", header=None, names=["chromosome", "start", "end"])
    df_windows["pos"] = (df_windows.start + (df_windows.end - df_windows.start)/2).apply(int)

    ##########################################

    ##### define SVsizes with the estimateSVSizes R function #####

    # define the max SVsize
    all_lens = [lenChrom for lenChrom in chrom_to_len.values() if lenChrom>=window_l]
    if len(all_lens)==0: all_lens = chrom_to_len.values()

    len_shortest_chr = min(all_lens)
    max_size_SV = int(len_shortest_chr*0.05)
    min_size_SV = min([200, max_size_SV-1])

    # define an Rscript that runs the estimateSVSizes
    rscript = "%s.Rscript_calculateSVsizes.R"%reference_genome
    rscript_lines = ["#!%s/bin/Rscript"%EnvDir,
                     "library(RSVSim, quietly=TRUE)",
                     "",
                     "sizes = estimateSVSizes(n=%i, minSize=%i, maxSize=%i, default='inversions', hist=FALSE)"%(n_BPs, min_size_SV, max_size_SV),
                     "print(sizes)"]

    open(rscript, "w").write("\n".join(rscript_lines))
    run_cmd("chmod u+x %s"%rscript)

    # run
    rscript_stdout = "%s.stdout"%rscript 
    rscript_stderr = "%s.stderr"%rscript 
    run_cmd("%s > %s 2>%s"%(rscript, rscript_stdout, rscript_stderr))

    # get the SV sizes
    SVsizes = make_flat_listOflists([[int(x) for x in l.split("]")[1].split()] for l in  open(rscript_stdout, "r").readlines()])

    ##############################################################

    ######## generate random breakpoints ########

    print_if_verbose("placing random breakpoints. Looking for %i of them"%n_BPs)

    # reshuffle df. This will ensure that the selection is random
    df_windows = df_windows.sample(frac=1)[["chromosome", "pos"]]

    # initialize the counter
    Ibp = 0

    # define all bpIDs
    all_BPids = np.array(df_windows.index)

    # define the Ileft, Iright combinations
    get_bedpe_series_brekends_inputs = []

    # iterate through each window, randomly
    for Ileft in all_BPids:

        # break if you already have enough BPs
        if Ibp>=n_BPs: break

        # iterate through the downstream windows
        for Iright in all_BPids[all_BPids>Ileft]:

            # break if you already have enough BPs
            if Ibp>=n_BPs: break

            # get the expected SVlen
            SV_len = random.choice(SVsizes)

            # skip if the length of the event is larger than SV_len
            r1 = df_windows.loc[Ileft]
            r2 = df_windows.loc[Iright]

            if r1["chromosome"]==r2["chromosome"] and (r2["pos"]-r1["pos"]>SV_len): continue

            get_bedpe_series_brekends_inputs.append((Ileft, Iright, Ibp, df_windows))
            Ibp += 1

    # get the df_bedpe from a pool
    with multiproc.Pool(threads) as pool:

        # run a dummy function
        df_bedpe = pd.DataFrame(pool.starmap(get_bedpe_series_brekends, get_bedpe_series_brekends_inputs))

        # close the pool
        pool.close() 
        pool.terminate()

    # get 3x as much intrachromosomal than interchromosomal events
    df_bedpe["intra_chromosome"] = df_bedpe.chrom1==df_bedpe.chrom2
    df_bedpe_intra = df_bedpe[df_bedpe.intra_chromosome]
    df_bedpe_inter = df_bedpe[~(df_bedpe.intra_chromosome)]

    max_len_df_inter = len(df_bedpe_intra)
    df_bedpe = df_bedpe_intra.append(df_bedpe_inter.iloc[0:max_len_df_inter])

    print(df_bedpe)

    # keep non-redundant events
    add_interval_bp = 100
    df_bedpe["affected_positions_arroundBp"] = df_bedpe.apply(lambda r: get_affected_positions_from_bedpe_r(r, extra_bp=add_interval_bp), axis=1)
    df_bedpe = get_df_bedpe_with_nonOverlappingBreakpoints(df_bedpe)

    # define the fields
    bedpe_fields = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"]

    #############################################

    return df_bedpe[bedpe_fields]





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



def run_perSVade_severalSamples(paths_df, cwd, common_args, threads=4, sampleID_to_parentIDs={}, samples_to_run=set(), repeat=False, job_array_mode="job_array", ploidy=1, variant_calling_fields=["#Uploaded_variation", "QUAL", "fb_DP", "fb_MQM", "fb_MQMR", "fb_PQA", "fb_PQR", "fb_QA", "fb_QR", "fb_fractionReadsCov", "fb_readsCovVar"]):

 
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

    print_if_verbose("Integrating all variants and CNV into one......")

    checkthathteintegrationmakessense

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




except:

    # define the error_log
    error_log = "building the db for RepeatModeler did not work. Check %s"%bulding_repModeler_db_std

    # if you are running in Nord3, this can be circumvented by running the makeblastdb sepparately
    if "BSC_MACHINE" in os.environ and os.environ["BSC_MACHINE"]=="nord3":

        # setup the db
        delete_folder(outdir)
        make_folder(outdir)
        shutil.copy2(reference_genome, genome_dir)

        # make the blast db
        makeblastdb_stdout = "%s/makeblastdb_nord3.std"%outdir
        print_if_verbose("runnning makeblastdb in nord3. The std is in %s"%makeblastdb_stdout)

        nord3_makeblastdb = "/apps/BLAST/2.2.28/bin/makeblastdb"
        run_cmd("cd %s && %s -out %s -parse_seqids -dbtype nucl -in %s > %s 2>&1"%(outdir, nord3_makeblastdb, name_database, name_database, makeblastdb_stdout), env=EnvName_RepeatMasker)

        remove_file(makeblastdb_stdout)

    else: raise ValueError(error_log)



  # setup variables if you are in specific machines
        if "BSC_MACHINE" in os.environ and os.environ["BSC_MACHINE"]=="nord3":

            # define paths
            nord3_makeblastdb = "/apps/BLAST/2.2.28/bin/makeblastdb"
            abblast_dir 
            rmblast_dir



            # make the blast db
            makeblastdb_stdout = "%s/makeblastdb_nord3.std"%outdir
            print_if_verbose("runnning makeblastdb in nord3. The std is in %s"%makeblastdb_stdout)

            
            run_cmd("cd %s && %s -out %s -parse_seqids -dbtype nucl -in %s > %s 2>&1"%(outdir, nord3_makeblastdb, name_database, name_database, makeblastdb_stdout), env=EnvName_RepeatMasker)

            remove_file(makeblastdb_stdout)




def get_set_overlappingIDs_df_CNV_r(r, df_CNV, equal_fields, approximate_fields, chromField_to_posFields, tol_bp, pct_overlap):

    """This function takes a row from a df_CNV and the same df, returning the set of IDs that are overlapping with this row according to several fields"""

    # define the df_overlap
    df_overlapping = df_CNV[(df_CNV.apply(lambda r_cnv: get_is_matching_predicted_and_known_rows(r, r_cnv, equal_fields, approximate_fields, chromField_to_posFields, tol_bp=tol_bp, pct_overlap=pct_overlap), axis=1))]

    return set(df_overlapping.index)


def get_nonRedundant_CNVcalls_coverage(df_CNV, df_vcf_forCNV, threads, pct_overlap=0.8, tol_bp=1000000000):

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

    # add the quality and AF, which are useful for the sorting of redundant variants
    df_CNV["QUAL"] = df_CNV.QUAL_mean
    df_CNV["AF"] = df_CNV.real_AF_min
    df_vcf_forCNV["QUAL"] = 1000000000
    df_vcf_forCNV["AF"] = 1.0

    # define all called SVs
    fields = ["ID", "chromosome", "start", "end", "SVTYPE", "type_CNVcall", "typeBPs", "QUAL", "AF"]
    all_df_CNV = df_CNV[fields].append(df_vcf_forCNV[fields])
    all_df_CNV = all_df_CNV.set_index("ID")

    # make sure that the ID is unique
    check_that_df_index_is_unique(all_df_CNV)

    # define the clusters of CNVs that are overlapping by >=80% of their extension
    list_clusters = get_list_clusters_overlapping_df_CNV(all_df_CNV, pct_overlap, tol_bp, threads)

    # add fields for sorting of redundant variants according to their type and the quality mesurements
    type_CNVcall_to_int = {"gridssClove":1, "coverage":0}
    all_df_CNV["type_CNVcall_int"] = all_df_CNV.type_CNVcall.apply(lambda x: type_CNVcall_to_int[x])

    typeBPs_to_int = {"RealBPs":3, "wholeChrom":2, "oneRealBP":1, "FilteredOutBPs":0}
    all_df_CNV["typeBPs_int"] = all_df_CNV.typeBPs.apply(lambda x: typeBPs_to_int[x])

    all_df_CNV["length"] = all_df_CNV.end - all_df_CNV.start

    # get the best IDs from each cluster
    best_NR_IDs = set(map( (lambda x: get_bestID_from_df_CNV_cluster(x, all_df_CNV) ), list_clusters))

    # get the df with these IDs
    df_CNV_NR = df_CNV[df_CNV.ID.isin(best_NR_IDs)]
    
    # at the end set the quality to a '.'
    df_CNV_NR["QUAL"] = "."

    return df_CNV_NR[initial_fields]


def get_list_clusters_overlapping_df_CNV(df_CNV, pct_overlap, tol_bp, threads):

    """Takes a df_CNV with chromosome, start, and end, where the ID is the index. It returns a list of sets, each set containing IDs of CNVs that overlap by >=80% and are from the same type """

    print_if_verbose("getting list_clusters_overlapping_df_CNV")

    # checks
    if len(set(df_CNV.SVTYPE).difference({"DUP", "DEL"}))>0: raise ValueError("SVTYPE is not properly formated")
    if len(set(df_CNV.type_CNVcall).difference({"gridssClove", "coverage"}))>0: raise ValueError("type_CNVcall is not properly formated")
    if len(set(df_CNV.typeBPs).difference({"RealBPs", "wholeChrom", "oneRealBP", "FilteredOutBPs"}))>0: raise ValueError("typeBPs is not properly formated")

    # init dict
    ID_to_overlappingIDs = {}

    # define the fields for comparison
    equal_fields = ["chromosome", "SVTYPE"]
    approximate_fields = ["start", "end"]
    chromField_to_posFields = {"chromosome":{"start": "start", "end": "end"}}

    # define two lists
    print_if_verbose("mapping each row to the cluster")

    all_IDs_list = list(df_CNV.index)
    inputs_fn = [(x[1], df_CNV, equal_fields, approximate_fields, chromField_to_posFields, tol_bp, pct_overlap) for x in df_CNV.iterrows()]

    # get the list of sets in parallel
    with multiproc.Pool(threads) as pool:
        list_overlappingIDs_sets = pool.starmap(get_set_overlappingIDs_df_CNV_r, inputs_fn) 
            
        pool.close()
        pool.terminate()

    # merge
    ID_to_overlappingIDs = dict(zip(all_IDs_list, list_overlappingIDs_sets))

    # initialize clusters list
    print_if_verbose("getting lists of clusters")
    list_clusters = get_list_clusters_from_dict(ID_to_overlappingIDs)

    # check
    all_IDs = set(df_CNV.index)
    all_IDs_in_cluster = set.union(*list_clusters)
    if all_IDs!=all_IDs_in_cluster: raise ValueError("all IDs should be in clusters")

    print_if_verbose("list_clusters_overlapping_df_CNV already ran. There are %i clusters"%len(list_clusters))
  
    return list_clusters




def get_list_clusters_overlapping_df_CNV(outdir, df_CNV, pct_overlap, threads):

    """Takes a df_CNV with chromosome, start, and end, where the ID is the index. It returns a list of sets, each set containing IDs of CNVs that overlap by >=pct_overlap and are from the same type.

    It will run bedtools intersect between the df_CNV against itself."""

    print_if_verbose("getting list_clusters_overlapping_df_CNV")

    # checks
    if len(set(df_CNV.SVTYPE).difference({"DUP", "DEL"}))>0: raise ValueError("SVTYPE is not properly formated")
    if len(set(df_CNV.type_CNVcall).difference({"gridssClove", "coverage"}))>0: raise ValueError("type_CNVcall is not properly formated")
    if len(set(df_CNV.typeBPs).difference({"RealBPs", "wholeChrom", "oneRealBP", "FilteredOutBPs"}))>0: raise ValueError("typeBPs is not properly formated")

    ####### GET A DF WITH THE INTERSECTION #######

    # add the combination of chromosome and SVtype. This will be used in the bedtools intersect to get the correct IDs
    df_CNV["chromosome_SVTYPE"] = df_CNV.chromosome + "_" + df_CNV.SVTYPE

    # add the ID and the numeric ID (the numeric ID takes less space in the bedtools intersect)
    df_CNV["ID"] = df_CNV.index
    df_CNV["numericID"] = list(range(0, len(df_CNV)))

    # define a directory where the bedtools will be run
    outdir_bedtools = "%s/running_get_list_clusters_overlapping_df_CNV"%outdir
    make_folder(outdir_bedtools)

    # define the input of a function that will take data from one chromosome and 
    inputs_fn = [(df_CNV, chromosome, outdir_bedtools) for chrom in set(df_CNV.chromosome_SVTYPE)] 



    # get 2 bed files, to run df_CNV against itself
    bedA = "%s/df_CNV_A.bed"%outdir_bedtools
    bedB = "%s/df_CNV_B.bed"%outdir_bedtools
    intersection_bed = "%s/intersection.bed"%outdir_bedtools

    df_CNV[["chromosome_SVTYPE", "start", "end", "ID"]].to_csv(bedA, sep="\t", header=False, index=False)
    df_CNV[["chromosome_SVTYPE", "start", "end", "ID"]].to_csv(bedB, sep="\t", header=False, index=False)

    df_CNV.pop("ID")

    # run bedtools intersect
    bedtools_stderr = "%s.generating.stderr"%intersection_bed
    print_if_verbose("Running bedtools intersect. The stderr is in %s"%bedtools_stderr)

    run_cmd("%s intersect -a %s -b %s -f %.2f -r -wa -wb | cut -f4,8 > %s 2>%s"%(bedtools, bedA, bedB, pct_overlap, intersection_bed, bedtools_stderr))

    # load as df, and keep only those where the IDs are different
    df_intersection = pd.read_csv(intersection_bed, sep="\t", header=None, names=["IDa", "IDb"])
    df_intersection = df_intersection[df_intersection.IDa!=df_intersection.IDb]

    # clean
    delete_folder(outdir_bedtools)

    ##########################################

    # define the ID_to_overlappingIDs from a groupby
    print_if_verbose("calculating the overlapping IDs")
    ID_to_overlappingIDs = dict(df_intersection.groupby("IDa").apply(lambda df_a: set(df_a["IDb"])))

    # define the IDs that are not overlapping 
    nonoverlapping_IDs = set(df_CNV.index).difference(set(ID_to_overlappingIDs))
    ID_to_overlappingIDs = {**ID_to_overlappingIDs, **dict(zip(nonoverlapping_IDs, [set()]*len(nonoverlapping_IDs)))}

    # initialize clusters list
    print_if_verbose("getting lists of clusters")
    list_clusters = get_list_clusters_from_dict(ID_to_overlappingIDs)

    # check
    all_IDs = set(df_CNV.index)
    all_IDs_in_cluster = set.union(*list_clusters)
    if all_IDs!=all_IDs_in_cluster: raise ValueError("all IDs should be in clusters")

    print_if_verbose("list_clusters_overlapping_df_CNV already ran. There are %i clusters"%len(list_clusters))
  
    return list_clusters



def get_coverage_per_window_for_chromosomeDF(chromosome_id, destination_dir, windows_bed, sorted_bam, replace, window_l):

    """Takes a chromosome id, a destination dir where to write files, a windows file (provided by generate_coverage_per_window_file_parallel) and a sorted bam and it generates a dataframe with the coverage stats"""

    # define the output coverage file
        
    # generate a randomID
    randID = id_generator(25)

    # define a file for the coverage
    windows_bed_chromsome = "%s.%s.%s.bed"%(windows_bed, chromosome_id, randID)
    egrepping_windows_stderr = "%s.generaing.stderr"%windows_bed_chromsome
    #print_if_verbose("Running egrep on chromosome. The stderr is in %s"%egrepping_windows_stderr)
    run_cmd("egrep '%s\t' %s > %s 2>%s"%(chromosome_id, windows_bed, windows_bed_chromsome, egrepping_windows_stderr))
    remove_file(egrepping_windows_stderr)

    # if there is nothing, return an empty df
    bamstats_fields = ["#chrom", "start", "end", "length", "mediancov_1", "nocoveragebp_1", "percentcovered_1"]
    if file_is_empty(windows_bed_chromsome): return pd.DataFrame(columns=bamstats_fields)

    # calculate extra threads
    extra_threads = multiproc.cpu_count() - 1

    # define a file prefix on which to calculate the coverage
    mosdepth_outprefix = "%s.mosdepth_output"%windows_bed_chromsome

    # get a df that has all the info
    df_coverage =  get_mosdepth_coverage_per_windows_output_likeBamStats(mosdepth_outprefix, sorted_bam, windows_bed_chromsome, replace=replace, extra_threads=extra_threads, chromosome_id=chromosome_id)

    remove_file(windows_bed_chromsome)

    return df_coverage[bamstats_fields]

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
                pool = multiproc.Pool(threads)

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


############### NONREDUNDANT DATA ###################


def get_list_clusters_for_chromosome_and_SVTYPE(outdir, chromosome_SVTYPE, df_CNV, pct_overlap):

    """
    This fucntion runs bedtools intersect to get regions that reciprocally intersect by >pct_overlap in chromosome_SVTYPE. It returns a series mapping the ID to overlappingIDs.

    Everything is done on the basis of the numeric ID

    """

    ######## RUN BEDTOOLS INTERSECT #########

    # define the beds
    bedA = "%s/%s_regionsA.bed"%(outdir, chromosome_SVTYPE)
    bedB = "%s/%s_regionsB.bed"%(outdir, chromosome_SVTYPE)
    intersection_bed = "%s/%s_intersection.bed"%(outdir, chromosome_SVTYPE)

    # get the chromosomal df
    df_c = df_CNV[df_CNV.chromosome_SVTYPE==chromosome_SVTYPE]
    if len(df_c)==0: raise ValueError("this can't be 0")

    # get the sorted df
    df_c = df_c[["chromosome_SVTYPE", "start", "end", "numericID"]].sort_values(by="start")

    # write the same df twice
    df_c.to_csv(bedA, sep="\t", header=False, index=False)
    df_c.to_csv(bedB, sep="\t", header=False, index=False)

    # run bedtools intersect
    bedtools_stderr = "%s.generating.stderr"%intersection_bed
    print_if_verbose("Running bedtools intersect. The stderr is in %s"%bedtools_stderr)

    run_cmd("%s intersect -a %s -b %s -f %.2f -r -wa -wb -sorted | cut -f4,8 > %s 2>%s"%(bedtools, bedA, bedB, pct_overlap, intersection_bed, bedtools_stderr))

    # clean files
    for f in [bedA, bedB, bedtools_stderr]: remove_file(f)

    # load as df, and keep only those where the IDs are different
    df_intersection = pd.read_csv(intersection_bed, sep="\t", header=None, names=["IDa", "IDb"])
    df_intersection = df_intersection[df_intersection.IDa!=df_intersection.IDb]

    #########################################

    # define the ID_to_overlappingIDs from a groupby
    ID_to_overlappingIDs = dict(df_intersection.groupby("IDa").apply(lambda df_a: set(df_a["IDb"])))

    # define the IDs that are not overlapping 
    nonoverlapping_IDs = set(df_c.numericID).difference(set(ID_to_overlappingIDs))
    ID_to_overlappingIDs = {**ID_to_overlappingIDs, **dict(zip(nonoverlapping_IDs, [set()]*len(nonoverlapping_IDs)))}

    # get the list of clusters
    list_clusters = get_list_clusters_from_dict(ID_to_overlappingIDs)

    # check that IDs are in a cluster
    all_IDs = set(df_c.numericID)
    all_IDs_in_cluster = set.union(*list_clusters)
    if all_IDs!=all_IDs_in_cluster: raise ValueError("all IDs should be in clusters")

    # clean the intersection bed
    remove_file(intersection_bed)

    del df_intersection

    return list_clusters

def get_list_clusters_overlapping_df_CNV(outdir, df_CNV, pct_overlap, threads):

    """Takes a df_CNV with chromosome, start, and end, where the ID is the index. It returns a list of sets, each set containing IDs of CNVs that overlap by >=pct_overlap and are from the same type.

    It will run bedtools intersect between the df_CNV against itself.

    The returned list of clusters are based on a numeric ID"""

    print_if_verbose("getting list_clusters_overlapping_df_CNV")

    # checks
    if len(set(df_CNV.SVTYPE).difference({"DUP", "DEL"}))>0: raise ValueError("SVTYPE is not properly formated")
    if len(set(df_CNV.type_CNVcall).difference({"gridssClove", "coverage"}))>0: raise ValueError("type_CNVcall is not properly formated")
    if len(set(df_CNV.typeBPs).difference({"RealBPs", "wholeChrom", "oneRealBP", "FilteredOutBPs"}))>0: raise ValueError("typeBPs is not properly formated")

    # add the combination of chromosome and SVtype. This will be used in the bedtools intersect to get the correct IDs
    df_CNV["chromosome_SVTYPE"] = df_CNV.chromosome + "_" + df_CNV.SVTYPE

    # define a directory where the bedtools will be run
    outdir_bedtools = "%s/running_get_list_clusters_overlapping_df_CNV"%outdir
    delete_folder(outdir_bedtools); make_folder(outdir_bedtools)

    # go through each chromosome and create the list of clusters
    list_clusters = []
    for chrom_SVTYPE in set(df_CNV.chromosome_SVTYPE): list_clusters += get_list_clusters_for_chromosome_and_SVTYPE(outdir_bedtools, chrom_SVTYPE, df_CNV, pct_overlap)

    print_if_verbose("list_clusters_overlapping_df_CNV already ran. There are %i clusters and %i regions"%(len(list_clusters), len(df_CNV)))
    
    delete_folder(outdir_bedtools)

    return list_clusters

def get_bestID_from_df_CNV_cluster(clustered_nuericIDs, df_CNV):

    """This function takes a df_CNV and the clusteredIDs. It returns the best clusterID"""

    # get the df with the lcuster IDs
    df = df_CNV.loc[clustered_nuericIDs].sort_values(by=["type_CNVcall_int", "typeBPs_int", "QUAL", "AF", "length"], ascending=False)

    return df.ID.iloc[0]

def get_nonRedundant_CNVcalls_coverage(outdir, df_CNV, df_vcf_forCNV, threads, pct_overlap=0.8):

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

    # add the quality and AF, which are useful for the sorting of redundant variants
    df_CNV["QUAL"] = df_CNV.QUAL_mean
    df_CNV["AF"] = df_CNV.real_AF_min
    df_vcf_forCNV["QUAL"] = 1000000000
    df_vcf_forCNV["AF"] = 1.0

    # define all called SVs
    fields = ["ID", "chromosome", "start", "end", "SVTYPE", "type_CNVcall", "typeBPs", "QUAL", "AF"]
    all_df_CNV = df_CNV[fields].append(df_vcf_forCNV[fields])

    # add the ID and the numeric ID (the numeric ID takes less space in the bedtools intersect). This will be the index
    all_df_CNV["numericID"] = list(range(0, len(all_df_CNV)))
    all_df_CNV = all_df_CNV.set_index("numericID", drop=False)

    # make sure that the ID is unique
    check_that_df_index_is_unique(all_df_CNV)

    # define the clusters of CNVs that are overlapping by >=80% of their extension
    list_clusters = get_list_clusters_overlapping_df_CNV(outdir, all_df_CNV, pct_overlap, threads)

    # add fields for sorting of redundant variants according to their type and the quality mesurements
    type_CNVcall_to_int = {"gridssClove":1, "coverage":0}
    all_df_CNV["type_CNVcall_int"] = all_df_CNV.type_CNVcall.apply(lambda x: type_CNVcall_to_int[x])

    typeBPs_to_int = {"RealBPs":3, "wholeChrom":2, "oneRealBP":1, "FilteredOutBPs":0}
    all_df_CNV["typeBPs_int"] = all_df_CNV.typeBPs.apply(lambda x: typeBPs_to_int[x])

    all_df_CNV["length"] = all_df_CNV.end - all_df_CNV.start

    # get the best IDs from each cluster
    best_NR_IDs = set(map( (lambda x: get_bestID_from_df_CNV_cluster(x, all_df_CNV) ), list_clusters))

    # get the df with these IDs
    df_CNV_NR = df_CNV[df_CNV.ID.isin(best_NR_IDs)]
    
    # at the end set the quality to a '.'
    df_CNV_NR["QUAL"] = "."

    return df_CNV_NR[initial_fields]


def get_list_clusters_for_chromosome_and_SVTYPE(outdir, chromosome_SVTYPE, df_CNV, pct_overlap):

    """
    This fucntion runs bedtools intersect to get regions that reciprocally intersect by >pct_overlap in chromosome_SVTYPE. It returns a series mapping the ID to overlappingIDs.

    Everything is done on the basis of the numeric ID

    """

    ######## RUN BEDTOOLS INTERSECT #########

    # define the beds
    bedA = "%s/%s_regionsA.bed"%(outdir, chromosome_SVTYPE)
    bedB = "%s/%s_regionsB.bed"%(outdir, chromosome_SVTYPE)
    intersection_bed = "%s/%s_intersection.bed"%(outdir, chromosome_SVTYPE)

    # get the chromosomal df
    df_c = df_CNV[df_CNV.chromosome_SVTYPE==chromosome_SVTYPE]
    if len(df_c)==0: raise ValueError("this can't be 0")

    # get the sorted df
    df_c = df_c[["chromosome_SVTYPE", "start", "end", "numericID"]].sort_values(by="start")

    # write the same df twice
    df_c.to_csv(bedA, sep="\t", header=False, index=False)
    df_c.to_csv(bedB, sep="\t", header=False, index=False)

    # run bedtools intersect
    bedtools_stderr = "%s.generating.stderr"%intersection_bed
    print_if_verbose("Running bedtools intersect. The stderr is in %s"%bedtools_stderr)

    run_cmd("%s intersect -a %s -b %s -f %.2f -r -wa -wb -sorted | cut -f4,8 > %s 2>%s"%(bedtools, bedA, bedB, pct_overlap, intersection_bed, bedtools_stderr))


    kjdhjkha

    # clean files
    for f in [bedA, bedB, bedtools_stderr]: remove_file(f)

    # load as df, and keep only those where the IDs are different
    df_intersection = pd.read_csv(intersection_bed, sep="\t", header=None, names=["IDa", "IDb"])
    df_intersection = df_intersection[df_intersection.IDa!=df_intersection.IDb]

    #########################################

    # define the ID_to_overlappingIDs from a groupby
    ID_to_overlappingIDs = dict(df_intersection.groupby("IDa").apply(lambda df_a: set(df_a["IDb"])))

    # define the IDs that are not overlapping 
    nonoverlapping_IDs = set(df_c.numericID).difference(set(ID_to_overlappingIDs))
    ID_to_overlappingIDs = {**ID_to_overlappingIDs, **dict(zip(nonoverlapping_IDs, [set()]*len(nonoverlapping_IDs)))}

    # get the list of clusters
    list_clusters = get_list_clusters_from_dict(ID_to_overlappingIDs)

    # check that IDs are in a cluster
    all_IDs = set(df_c.numericID)
    all_IDs_in_cluster = set.union(*list_clusters)
    if all_IDs!=all_IDs_in_cluster: raise ValueError("all IDs should be in clusters")

    kjadhdjkhakadh

    # clean the intersection bed
    remove_file(intersection_bed)

    del df_intersection

    return list_clusters

def get_list_clusters_overlapping_df_CNV(outdir, df_CNV, pct_overlap, threads):

    """Takes a df_CNV with chromosome, start, and end, where the ID is the index. It returns a list of sets, each set containing IDs of CNVs that overlap by >=pct_overlap and are from the same type.

    It will run bedtools intersect between the df_CNV against itself.

    The returned list of clusters are based on a numeric ID"""

    print_if_verbose("getting list_clusters_overlapping_df_CNV")

    # checks
    if len(set(df_CNV.SVTYPE).difference({"DUP", "DEL"}))>0: raise ValueError("SVTYPE is not properly formated")
    if len(set(df_CNV.type_CNVcall).difference({"gridssClove", "coverage"}))>0: raise ValueError("type_CNVcall is not properly formated")
    if len(set(df_CNV.typeBPs).difference({"RealBPs", "wholeChrom", "oneRealBP", "FilteredOutBPs"}))>0: raise ValueError("typeBPs is not properly formated")

    # add the combination of chromosome and SVtype. This will be used in the bedtools intersect to get the correct IDs
    df_CNV["chromosome_SVTYPE"] = df_CNV.chromosome + "_" + df_CNV.SVTYPE

    # define a directory where the bedtools will be run
    outdir_bedtools = "%s/running_get_list_clusters_overlapping_df_CNV"%outdir
    delete_folder(outdir_bedtools); make_folder(outdir_bedtools)

    # go through each chromosome and create the list of clusters
    list_clusters = []
    for chrom_SVTYPE in set(df_CNV.chromosome_SVTYPE): list_clusters += get_list_clusters_for_chromosome_and_SVTYPE(outdir_bedtools, chrom_SVTYPE, df_CNV, pct_overlap)

    print_if_verbose("list_clusters_overlapping_df_CNV already ran. There are %i clusters and %i regions"%(len(list_clusters), len(df_CNV)))

    kljhadjkhdakjhda
    
    delete_folder(outdir_bedtools)

    return list_clusters



def get_bed_intersecting_numericIDs(df_CNV, outdir, threads, pct_overlap, replace):

    """This function takes a df_CNV and writes a bed under outdir with the intersecting numeric IDs"""

    print_if_verbose("running get_bed_intersecting_numericIDs")

    # checks
    if len(set(df_CNV.SVTYPE).difference({"DUP", "DEL"}))>0: raise ValueError("SVTYPE is not properly formated")
    if len(set(df_CNV.type_CNVcall).difference({"gridssClove", "coverage"}))>0: raise ValueError("type_CNVcall is not properly formated")
    if len(set(df_CNV.typeBPs).difference({"RealBPs", "wholeChrom", "oneRealBP", "FilteredOutBPs"}))>0: raise ValueError("typeBPs is not properly formated")

    # define the final file
    bed_intersecting_numericIDs = "%s/intersecting_numericIDs.bed"%outdir

    if file_is_empty(bed_intersecting_numericIDs) or replace is True:

        # add the combination of chromosome and SVtype. This will be used in the bedtools intersect to get the correct IDs
        df_CNV["chromosome_SVTYPE"] = df_CNV.chromosome + "_" + df_CNV.SVTYPE

        # define a directory where the bedtools will be run
        outdir_bedtools = "%s/running_get_list_clusters_overlapping_df_CNV"%outdir
        delete_folder(outdir_bedtools); make_folder(outdir_bedtools)

        # get all beds on the intersection
        inputs_fn = [(outdir_bedtools, chromosome_SVTYPE, df_CNV, pct_overlap) for chromosome_SVTYPE in set(df_CNV.chromosome_SVTYPE)]

        # generate the files 
        with multiproc.Pool(threads) as pool:
            pool.starmap(get_intersecting_bed_chromosome_SVTYPE, inputs_fn) 
                
            pool.close()
            pool.terminate()

        # cat all the files in one
        bed_intersecting_numericIDs_tmp = "%s.tmp"%bed_intersecting_numericIDs
        bed_intersecting_numericIDs_stderr = "%s.generating.stderr"%bed_intersecting_numericIDs 

        print_if_verbose("catting all files. The stderr is in %s"%bed_intersecting_numericIDs_stderr)
        run_cmd("cat %s/* > %s 2>%s"%(outdir_bedtools, bed_intersecting_numericIDs_tmp, bed_intersecting_numericIDs_stderr))

        # clean
        remove_file(bed_intersecting_numericIDs_stderr)
        delete_folder(outdir_bedtools)

        # rename
        os.rename(bed_intersecting_numericIDs_tmp, bed_intersecting_numericIDs)

    return bed_intersecting_numericIDs




def get_intersecting_bed_chromosome_SVTYPE(outdir, chromosome_SVTYPE, df_CNV, pct_overlap):

    """
    This fucntion runs bedtools intersect to get regions that reciprocally intersect by >pct_overlap in chromosome_SVTYPE. It returns the file where this intersection was written

    Everything is done on the basis of the numeric ID.

    """

    # define the beds
    bedA = "%s/%s_regionsA.bed"%(outdir, chromosome_SVTYPE)
    bedB = "%s/%s_regionsB.bed"%(outdir, chromosome_SVTYPE)
    intersection_bed = "%s/%s_intersection.bed"%(outdir, chromosome_SVTYPE)

    # get the chromosomal df
    df_c = df_CNV[df_CNV.chromosome_SVTYPE==chromosome_SVTYPE]
    if len(df_c)==0: raise ValueError("this can't be 0")

    # get the sorted df
    df_c = df_c[["chromosome_SVTYPE", "start", "end", "numericID"]].sort_values(by="start")

    # write the same df twice
    df_c.to_csv(bedA, sep="\t", header=False, index=False)
    df_c.to_csv(bedB, sep="\t", header=False, index=False)

    # run bedtools intersect
    bedtools_stderr = "%s.generating.stderr"%intersection_bed
    print_if_verbose("Running bedtools intersect. The stderr is in %s"%bedtools_stderr)

    run_cmd("%s intersect -a %s -b %s -f %.2f -r -wa -wb -sorted | cut -f4,8 > %s 2>%s"%(bedtools, bedA, bedB, pct_overlap, intersection_bed, bedtools_stderr))

    # clean files
    for f in [bedA, bedB, bedtools_stderr]: remove_file(f)

    return intersection_bed

dlajhdjkah

# get a file that has the 

# get a bed with the intersecting numericIDs
bed_intersecting_numericIDs = get_bed_intersecting_numericIDs(all_df_CNV, outdir, threads, pct_overlap, replace)


adlhjalhjdahkdahk

# define  a series that will have the number of overlapping features for each numericID
bed_number_overlaps = "%s.number_overlapping_regions.tab"%bed_intersecting_numericIDs
bed_number_overlaps_stderr = "%s.stderr"%bed_number_overlaps

run_cmd("cut -f1 %s | uniq -c | sed 's/^[ \t]*//g'> %s 2>%s"%(bed_intersecting_numericIDs, bed_number_overlaps, bed_number_overlaps_stderr))
remove_file(bed_number_overlaps_stderr)

n_overlaps_series = pd.read_csv(bed_number_overlaps, sep=" ", header=None, names=["n_overlaps", "numericID"]).set_index("numericID")["n_overlaps"]

# add to the df
all_df_CNV["n_overlaps"] = n_overlaps_series[all_df_CNV.index]
if len(all_df_CNV)!=initial_len_all_df_CNV: raise ValueError("There was an error with the number of overlaps")

ljhadlhadlhdakhadkjhkajdh

print(n_overlaps_series)

adkjhdkjha


def get_last_position_from_pipeSepparatedString(x):

    """Takes the last position of a pipe-sepparated string"""

    return x.split("|")[-1]

def  get_query_df_CNVregions_with_metadata_for_redundance(query_df, df_gridss, chrom_to_len):

    """This function takes a CNV of regions that match the df_gridss and it returns it with several added fields that are necessary to get non-redundant regions by the function get_nonRedundant_CNVcalls_coverage """

    if len(query_df)==0: return query_df

    print_if_verbose("adding metadata to query_df")

    # keep
    df_gridss = cp.deepcopy(df_gridss)

    # define fields that are generally applicable
    query_df["type_CNVcall"] = "coverage"
    query_df["typeBPs"] = query_df.ID.apply(get_last_position_from_pipeSepparatedString)
    query_df["SVTYPE"] = "DUP"

 
    ######### REFORMAT df_gridss #########

    # add the 0-based-pos
    df_gridss["POS_0based"] = df_gridss.POS - 1

    # keep important fields
    df_gridss = df_gridss[["#CHROM", "POS_0based", "QUAL", "real_AF"]]
    df_gridss["bend_POS"] = df_gridss["#CHROM"] + "_" + df_gridss.POS_0based.apply(str)
    df_gridss = df_gridss.set_index("bend_POS")
    
    # define the maxQUAL
    maxQUAL = max(df_gridss.QUAL)+1

    # add the first and last positions
    for chrom, length in chrom_to_len.items():

        # add the start
        df_gridss = df_gridss.append(pd.DataFrame({"%s_0"%chrom: {"QUAL":maxQUAL, "real_AF":1.0, "#CHROM":chrom, "POS_0based":0}}).transpose())

        # add the end
        df_gridss = df_gridss.append(pd.DataFrame({"%s_%i"%(chrom, length): {"QUAL":maxQUAL, "real_AF":1.0, "#CHROM":chrom, "POS_0based":length}}).transpose())

    # sort
    df_gridss = df_gridss.sort_values(by=["#CHROM", "POS_0based", "QUAL"], ascending=False)
    df_gridss = df_gridss.drop_duplicates(subset=["#CHROM", "POS_0based"], keep="first")

    check_that_df_index_is_unique(df_gridss)

    ########################################

    # init len
    initial_query_df_len = len(query_df)

    # add the indicators of each breakend
    for field in ["start", "end"]:

        # add the position
        bend_POSitions = np.array(query_df.chromosome + "_" + query_df[field].apply(str))

        # add the metadata
        query_df["QUAL_%s"%field] = df_gridss.loc[bend_POSitions, "QUAL"].values
        query_df["real_AF_%s"%field] = df_gridss.loc[bend_POSitions, "real_AF"].values

    # debug
    if len(query_df)!=initial_query_df_len: raise ValueError("something went wrong with the merging")
    for f in ["QUAL_start", "QUAL_end", "real_AF_start", "real_AF_end"]: 
        if any(pd.isna(query_df[f])): raise ValueError("There are NaNs after the merhing")

    # get the important fields
    query_df["QUAL_mean"] = query_df[["QUAL_start", "QUAL_end"]].apply(np.mean, axis=1)
    query_df["real_AF_min"] = query_df[["real_AF_start", "real_AF_end"]].apply(min, axis=1)

    return query_df



# get the query df with metadata to skip the redundant regions
#query_df = get_query_df_CNVregions_with_metadata_for_redundance(query_df, df_gridss, chrom_to_len)

# skip the redundant regions
#df_empty = pd.DataFrame(columns=list(query_df.keys()))
#query_df = get_nonRedundant_CNVcalls_coverage(outdir, query_df, df_empty, threads, replace, pct_overlap=0.9)


def run_perSVade_severalSamples(paths_df, cwd, common_args, threads=4, sampleID_to_parentIDs={}, samples_to_run=set(), repeat=False, job_array_mode="job_array", ploidy=1, get_integrated_dfs=True):

 
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

    # keep a dict with the paths
    sampleID_to_dataDict = {}

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
                   
        # define the cmd          
        cmd = "%s -f1 %s -f2 %s -o %s --ploidy %i %s"%(perSVade_py, reads1, reads2, outdir, ploidy, common_args)

        # add cmd if necessary
        if any([file_is_empty(x) for x in success_files]) or repeat is True: all_cmds.append(cmd)

        # keep data
        #sampleID_to_dataDict[sampleID] = 

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


    if get_integrated_dfs is True:

        print_if_verbose("Integrating all variants and CNV into one......")


    checkthathteintegrationmakessense

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

    if get_integrated_dfs is True:

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

    e




    return variant_calling_df




def get_bed_df_from_variantID(varID):

    """Takes a variant ID, such as the ones in SV_CNV vcf 'INFO_variantID'. It returns a df with all chromosome-start-end information that should be matched to be considered as the same variant"""

    # get the ID svtype
    svtype = varID.split("|")[0]

    # inferred by coverage 
    if svtype in {"coverageDUP", "coverageDEL", "TDUP", "DEL", "INV"} : 

        chrom = "%s_%s"%(svtype, varID.split("|")[1].split(":")[0])
        start = int(varID.split("|")[1].split(":")[1].split("-")[0])
        end = int(varID.split("|")[1].split(":")[1].split("-")[1])

        dict_bed = {0 : {"chromosome":chrom, "start":start, "end":end, "ID":varID, "type_overlap":"both"}}

    elif svtype.endswith("like"):

        posA, posB = varID.split("|")[1].split("-")


        chromA = "%s_%s"%(svtype, posA.split(":")[0])
        chromB = "%s_%s"%(svtype, posB.split(":")[0])

        startA = int(posA.split(":")[1])
        endA = startA + 1

        startB = endB = int(posB.split(":")[1])
        endB = startB + 1

        dict_bed = {0 : {"chromosome":chromA, "start":startA, "end":endA, "ID":varID+"-A", "type_overlap":"bp_pos"},
                    1 : {"chromosome":chromB, "start":startB, "end":endB, "ID":varID+"-B"}}

    else: raise ValueError("%s has not been parsed"%varID)


    # get as df
    df_bed = pd.DataFrame(dict_bed).transpose()

    # add the variantID, which will be useful to track, and is not necessarily unique
    df_bed["variantID"] = varID

    return df_bed



    # inferred by coverage 
    if IDsvtype=="coverageDUP" and r["INFO_BPS_TYPE"] in {"RealBPs", "wholeChrom"}: 
        x = [r["POS"], r["INFO_END"]]
        color = "cyan"
        symbol = "circle"
        name = "coverageDUP_highConf"

    elif IDsvtype=="coverageDUP" and r["INFO_BPS_TYPE"] in {"FilteredOutBPs", "oneRealBP"}: 
        x = [r["POS"], r["INFO_END"]]
        color = "cyan"
        symbol = "square"
        name = "coverageDUP_lowConf"


    elif IDsvtype=="coverageDEL" and r["INFO_BPS_TYPE"] in {"RealBPs", "wholeChrom"}: 
        x = [r["POS"], r["INFO_END"]]
        color = "magenta"
        symbol = "circle"
        name = "coverageDEL_highConf"


    elif IDsvtype=="coverageDEL" and r["INFO_BPS_TYPE"] in {"FilteredOutBPs", "oneRealBP"}: 
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

    elif IDsvtype=="INS" and r["INFO_SVTYPE"]=="DUP": 
        x = [r["POS"], r["INFO_END"]]
        color = "navy"
        symbol = "circle"
        name = "copy-paste insertions"

    elif IDsvtype=="CVD" and r["INFO_SVTYPE"]=="DUP": 
        x = [r["POS"], r["INFO_END"]]
        color = "green"
        symbol = "circle"
        name = "complex inverted duplications"


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

    elif not IDsvtype.endswith("like"):

        # define an X for all
        x = [r["POS"]]

        needtoworkhere

        # all the classified BND go into the same name 
        if r["INFO_SVTYPE"]=="BND": 
            symbol = "x"
            name = "classified SVs BND"
            color = "black"

        elif r["INFO_SVTYPE"]=="insertionBND": 

            pass
            youshouldworkonthis

            x = [r["POS"]]
            color = "black"
            symbol = "cross"
            name = "classified SVs insertions"







    else:
        print(r["ID"])
        print(r["INFO_SVTYPE"])
        print(r["INFO_BPS_TYPE"])
        raise ValueError("not considered case")



def get_Gene_Feature_varAnnot_from_VepOutput_r(r, gff_df):

    """Takes a row of the vep output and returns the gene and the feature according to upmost_parent of the gff_df (Gene) and ID (Feature) """

    # define the annotation fields
    annotation_fields = [k for k in gff_df.keys() if k.startswith("ANNOTATION_")]

    # get the fields where to look
    fields = ["upmost_parent", "ANNOTATION_ID", "ID", "ANNOTATION_Dbxref_GeneID"] + annotation_fields
    fields = [f for f in fields if f in gff_df.keys() and f not in {"start", "end", "numeric_idx"}]

    # intergenic vars
    if r["Gene"]=="-": gene = feature = "-"

    else:

        ######### FIND THE GENE ########

        # get the genes df
        df_genes = gff_df[gff_df.feature.isin({"gene", "pseudogene"})]

        # find the upmost_parent from the different fields
        for f in fields:

            # see if the df is here
            df = df_genes[df_genes[f]==r["Gene"]]

            # once you find a match, return
            if len(df)>1: 
                gene = df.iloc[0]["upmost_parent"]
                break

        ################################

        ######### FIND THE ANNOTATION ######### 

        # get the all the gff for this gene
        gff_df = gff_df[gff_df.upmost_parent==gene]

        print(gff_df)

        kajdhkhgad




        #######################################

    return pd.Series({"Gene":gene, "Feature":feature})


def get_Gene_varAnnot_from_gene(gene, genes_df):

    """Get the gene so that it is in the gff_df"""

    # set the gene as a str
    gene = str(gene)

    # debug the genes
    if gene=="-": return "-"

    # get gene of is is in upmost parent
    if gene in genes_df["upmost_parent"]: return gene

    # try to assign it by geneID_Dbxref
    elif "geneID_Dbxref" in set(genes_df.keys()) and gene in set(genes_df.geneID_Dbxref): 

        all_upmost_parents = set(genes_df[genes_df.geneID_Dbxref==gene].upmost_parent)

        if len(all_upmost_parents)!=1: raise ValueError("%s can be assigned (through dbxref geneID) to these upmost parents: %s"%(gene, all_upmost_parents))

        return next(iter(all_upmost_parents))

    else:  raise ValueError("%s could not be found in the gff"%gene)

def get_GeneID_from_dbxref(x):

    """returns geneID"""

    all_geneIDs_list = [y.split(":")[1] for y in x.split(",") if y.startswith("GeneID")]

    if len(all_geneIDs_list)==0: return "no_gene"
    elif len(all_geneIDs_list)==1: return all_geneIDs_list[0]
    else: raise ValueError("%s is not valid"%all_geneIDs_list)



def get_sortedBam_with_duplicatesMarked(sorted_bam, threads=4, replace=False, remove_duplicates):

    """This function takes a sorted bam and returns the equivalent with the duplicates marked with picard MarkDuplicates. It also indexes this bam"""

    # define dirs
    sorted_bam_dupMarked = "%s.MarkDups.bam"%sorted_bam
    sorted_bam_dupMarked_tmp = "%s.MarkDups.tmp.bam"%sorted_bam
    sorted_bam_dupMarked_metrics = "%s.MarkDups.metrics"%sorted_bam

    if file_is_empty(sorted_bam_dupMarked) or replace is True:
        print_if_verbose("marking duplicate reads")

        # define the java memory
        #javaRamGb = int(get_availableGbRAM(get_dir(sorted_bam))*fractionRAM_to_dedicate) # At Broad, we run MarkDuplicates with 2GB Java heap (java -Xmx2g) and 10GB hard memory limit
        #javaRamGb = int(get_availableGbRAM(get_dir(sorted_bam))*0.5) # At Broad, we run MarkDuplicates with 2GB Java heap (java -Xmx2g) and 10GB hard memory limit
        javaRamGb = int(get_availableGbRAM(get_dir(sorted_bam)) - 2) # rule of thumb from GATK
        #javaRamGb = 4 # this is from a post from 2011, reccommended for a 170Gb RAM

        # define the MAX_RECORDS_IN_RAM
        MAX_RECORDS_IN_RAM = int(250000*javaRamGb*0.8) # 250,000 reads for each Gb given (for SortSam, I don't know if this will work for Picard tools)

        # define the number of MAX_FILE_HANDLES_FOR_READ_ENDS_MAP
        max_nfilehandles = int(subprocess.check_output("ulimit -n", shell=True))
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = int(max_nfilehandles*0.5) # a little lower than ulimit -n

        # SORTING_COLLECTION_SIZE_RATIO is 0.25 by default. If I have memory issues I can reduce this number.

        markduplicates_std = "%s.markingDuplicates.std"%sorted_bam
        print_if_verbose("running MarkDuplicates with %iGb of RAM and %i MAX_FILE_HANDLES_FOR_READ_ENDS_MAP. The std is in %s"%(javaRamGb, MAX_FILE_HANDLES_FOR_READ_ENDS_MAP, markduplicates_std))

        # running with the traditional MarkDuplicates implementation
        run_cmd("%s -Xmx%ig MarkDuplicates I=%s O=%s M=%s ASSUME_SORT_ORDER=coordinate MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=%i MAX_RECORDS_IN_RAM=%i > %s 2>&1"%(picard_exec, javaRamGb, sorted_bam, sorted_bam_dupMarked_tmp, sorted_bam_dupMarked_metrics, MAX_FILE_HANDLES_FOR_READ_ENDS_MAP, MAX_RECORDS_IN_RAM, markduplicates_std), env=EnvName_picard)
        #REMOVE_DUPLICATES=Boolean

        remove_file(markduplicates_std)

        # keep
        os.rename(sorted_bam_dupMarked_tmp, sorted_bam_dupMarked)

    # index the bam with the duplicate reads
    index_sorted_bam_dupMarked = "%s.bai"%sorted_bam_dupMarked
    if file_is_empty(index_sorted_bam_dupMarked) or replace is True:
        index_bam(sorted_bam_dupMarked, threads=threads)

    return sorted_bam_dupMarked


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

        print(gc_df)

        lkadjladkad


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


def get_y_corrected_by_x_LOWESS_crossValidation(df, xfield, yfield, outdir, threads, replace):

    """This function takes an x and a y series, returning the y corrected by x. This y corrected is y/(y predicted from LOWESS from x). The index must be unique. The best parameters are taken with 10 fold cross validation"""

    make_folder(outdir)

    # keep
    df = cp.deepcopy(df)

    # check that the index is unique
    if len(df.index)!=len(set(df.index)): raise ValueError("The index should be unique")

    # get index as list
    initial_index = list(df.index)

    # sort by the x
    df = df.sort_values(by=[xfield, yfield])

    ########## GET THE DF BENCHMARKED DF 10xCV ########## 

    # define the df_benckmarking file
    df_benchmarking_file = "%s/df_benckmarking.tab"%outdir

    if file_is_empty(df_benchmarking_file) or replace is True:
        print_if_verbose("getting benchmarking for %s vs %s"%(xfield, yfield))

        # init benchmarking df
        dict_benchmarking = {}; Idict = 0

        # define parms
        n_frac = 10
        n_it = 5

        # define all the indices
        all_idx = set(df.idx)

        # iterate through several parameters of the lowess fitting
        for frac in np.linspace(0.001, 0.3, n_frac): # the fraction of data points used to weight the 
            for it in range(1, n_it+1): # the number of residual-based reweightings to perform.
                
                # init rsquares
                rsquares_cv = []

                # iterate through 10-fold chross validation 
                kfold = KFold(n_splits=10, shuffle=True, random_state=1)
                kfold = 10
                for cvID,  in range(1,11):

                    df.iloc[::5, :]


                    print(train_idx)
                    print(test_idx)
                    print(initial_index)

                    kjahdkjhad
                    
                    # get the data, sorted by x
                    df_train = df.iloc[train_idx].sort_values(by=[xfield, yfield])
                    df_test = df.iloc[test_idx].sort_values(by=[xfield, yfield])

                    # get values
                    xtrain = df_train[xfield].values
                    ytrain = df_train[yfield].values

                    xtest = df_test[xfield].values
                    ytest = df_test[yfield].values

                    # get the test yfit based on the train 
                    ytest_predicted = statsmodels_api.nonparametric.lowess(endog=ytrain, exog=xtrain, frac=frac, it=it, xvals=xtest, is_sorted=True, missing="raise") # return sorted returns values storted by exog

                    if len(ytest_predicted)!=len(ytest): raise ValueError("xtest and ytest are not the same")

                    # calculate the rsquare, making sure it is a float
                    if any(pd.isna(ytest_predicted)): rsquare = 0.0
                    else: rsquare = r2_score(ytest, ytest_predicted)
                    if pd.isna(rsquare): rsquare = 0.0

                    # keep
                    rsquares_cv.append(rsquare)

                # keep
                std = np.std(rsquares_cv)
                mean_rsquare = np.mean(rsquares_cv)
                dict_benchmarking[Idict] = {"frac":frac, "it":int(it), "mean_rsquare":mean_rsquare, "inverse_std_rsquare":1/std, "std_rsquare":std}; Idict+=1

                print_if_verbose(frac, it, mean_rsquare)


        # get as df
        df_benchmarking = pd.DataFrame(dict_benchmarking).transpose()

        # save
        save_df_as_tab(df_benchmarking, df_benchmarking_file)

    # load
    df_benchmarking  = get_tab_as_df_or_empty_df(df_benchmarking_file).sort_values(by=["mean_rsquare", "inverse_std_rsquare"], ascending=False)

    # get the best parameters
    best_parms_series = df_benchmarking.iloc[0]

    ##################################################### 

    # get the fit data

    print_if_verbose("performing LOWESS regression with best parameters")

    # get sorted df
    df = df.sort_values(by=[xfield, yfield])

    # when no fit could be obtained, define the corrected_y_values equal to the original ones
    if best_parms_series.mean_rsquare==0: y_corrected = df[yfield]

    else:

        # get the y predicted with the best parms
        best_frac = best_parms_series["frac"]
        best_it = int(best_parms_series["it"])
        lowess_results = statsmodels_api.nonparametric.lowess(endog=df[yfield], exog=df[xfield], frac=best_frac, it=best_it, xvals=None, is_sorted=True, missing="raise", return_sorted=True) # return sorted returns values storted by exog

        predicted_y_values = pd.Series(lowess_results[:,1], index=df.index)
        y_corrected = df[yfield] / predicted_y_values

        # debug 
        if any(pd.isna(predicted_y_values)) or any(pd.isna(y_corrected)): raise ValueError("there should be no NaNs")

        final_rsquare = r2_score(df[yfield], predicted_y_values)

        ##############################

        ######### MAKE PLOTS #########

        fig = plt.figure(figsize=(5,5))
        plt.plot(df[xfield], df[yfield], "o", alpha=0.2, color="gray", label="raw data")
        plt.plot(df[xfield], predicted_y_values, "-", color="red", label="LOWESS fit")

        plt.title("Fitting LOWESS with frac=%.3f and it=%i. final R2=%.3f. 10x CV R2=%.3f +- %.3f (SD)\n"%(best_frac, best_it, final_rsquare, best_parms_series["mean_rsquare"], best_parms_series["std_rsquare"]))
        plt.legend(bbox_to_anchor=(1, 1))
        plt.xlabel(xfield)
        plt.ylabel(yfield)

        fig.savefig("%s/coverage.pdf"%(outdir), bbox_inches='tight')
        plt.close(fig)

        ############################

    # debug
    if any(pd.isna(y_corrected)): raise ValueError("there should be no NaNs")

    # return in the initial order
    return y_corrected.loc[initial_index]



def get_genomeGraph_object(reference_genome, df_gridss, min_AF_homo, min_QUAL, threads, replace):

    """This function takes a reference genome and a df_gridss with the breakpoints. It gererates a graph were each position is one node. 
    """

    # define the files
    prefix = "%s.graph_breakpoints_minQUAL=%.3f_minAFhomo=%.3f"%(reference_genome, min_QUAL, min_AF_homo)
    genome_graph_file = "%s.graph.py"%prefix
    positions_df_file = "%s.df_positions.tab"%prefix

    if any([file_is_empty(x) for x in {genome_graph_file, positions_df_file}]) or replace is True:
        print_if_verbose("getting genome graph object for %s"%reference_genome)

        # map each chromosome to an offset
        chrom_to_lenSeq = get_chr_to_len(reference_genome)

        # deffine an offset for each chromosome, which is necessary to keep all the positions of the genome as independent numbers
        chrom_to_offset = {}
        current_offset = 0
        for chrom, seqLen in chrom_to_lenSeq.items():
            chrom_to_offset[chrom] = current_offset
            current_offset+=seqLen

        # create the undirected graph, meaning that all edges work in both directions
        genome_graph = igraph.Graph(directed=False)

        # add one vertex (node) for each position in the genome
        npositions = sum(chrom_to_lenSeq.values())
        genome_graph.add_vertices(npositions)

        # create a df with the positions
        positions_df = pd.DataFrame()
        for chrom, offset in chrom_to_offset.items():

            df_c = pd.DataFrame({"position" : list(range(0, chrom_to_lenSeq[chrom]))})
            df_c["chromosome"] = chrom
            df_c["offset"] = offset
            df_c["graph_node"] = df_c.position + df_c.offset

            positions_df = positions_df.append(df_c[["chromosome", "position", "graph_node"]])

        # set the index
        positions_df = positions_df.set_index("graph_node", drop=False)

        # add the plus 1
        positions_df["position_1based"] = positions_df.position + 1

        # define the nodes that are the end of chromosomes
        chromosome_start_nodes = {offset for chrom, offset in chrom_to_offset.items()}
        chromosome_end_nodes = {(offset + chrom_to_lenSeq[chrom] - 1) for chrom, offset in chrom_to_offset.items()}
        all_nodes = set(range(npositions))

        # initialize a series with the edges of the adjacent positions
        all_edges = set(map(lambda pos: (pos, pos+1), all_nodes.difference(chromosome_end_nodes)))

        ######### MODIFY EDGES ACCORDING TO df_gridss #########

        if df_gridss is not None:
  
            # filter the df_gridss and check that it is unique
            df_gridss = df_gridss[df_gridss.QUAL>=min_QUAL]
            if len(set(df_gridss.ID))!=len(df_gridss): raise ValueError("df_gridss.ID should be unique")

            if len(df_gridss)>0:

                # keep only bends that have mates
                eventID_to_nBends = df_gridss.groupby("eventID_as_clove").apply(len)
                if any(eventID_to_nBends>2): raise ValueError("there can't be events with more than 2 breakend")
                eventID_to_nBends = dict(eventID_to_nBends)
                df_gridss["n_bends"] = df_gridss.eventID_as_clove.apply(lambda x: eventID_to_nBends[x])
                df_gridss = df_gridss[df_gridss.n_bends==2]

                if len(df_gridss)>0:

                    # add whether the breakend is homozygous. For homozygous breakends the adjacent edges will be removed
                    df_gridss["is_homo_bend"] = df_gridss.real_AF>=min_AF_homo

                    # keep only important fields
                    df_gridss = df_gridss[["#CHROM", "ALT", "POS", "is_homo_bend", "ID", "INFO_MATEID", "real_AF"]]

                    # add the node position 
                    df_gridss = df_gridss.merge(positions_df[["chromosome", "position_1based", "graph_node"]], left_on=["#CHROM", "POS"], right_on=["chromosome", "position_1based"], validate="many_to_one", how="left")

                    # add the node of the mate
                    ID_to_position = dict(df_gridss.set_index("ID")["graph_node"])
                    df_gridss["mate_graph_node"] = df_gridss.INFO_MATEID.apply(lambda x: ID_to_position[x])

                    # add all the edges related to the breakpoints
                    def get_df_gridss_sorted_edge(r): return sorted([r["graph_node"], r["mate_graph_node"]])
                    df_gridss["sorted_edge"] = df_gridss.apply(get_df_gridss_sorted_edge, axis=1).apply(tuple)
                    all_edges.update(set(df_gridss.sorted_edge))

                    # for homozygous breakends, remove the corresponding adjacent edges. Make sure that these are no adjacent breakpoints
                    df_gridss["distance_between_nodes"] = (df_gridss.graph_node - df_gridss.mate_graph_node).apply(abs)
                    df_gridss_homo = df_gridss[(df_gridss.is_homo_bend) & (df_gridss.distance_between_nodes>1)]

                    if len(df_gridss_homo)>0:

                        # add the compromised reference breakpoint
                        df_gridss_homo["compromised_ref_breakpoint"] = df_gridss_homo.apply(get_compromised_ref_breakpoint, axis=1)

                        all_ref_breakpoints_to_remove = set(df_gridss_homo.compromised_ref_breakpoint)

                        # remove the all_ref_breakpoints_to_remove
                        all_edges = all_edges.difference(all_ref_breakpoints_to_remove)

        #######################################################

        # add the edges to the graph
        genome_graph.add_edges(all_edges)
        print_if_verbose("genome graph got")
 
        ########## CALCULATE THE DISTANCE TO THE CLOSEST TELOMERE ##########

        # add the number of neighbors that each position has
        def get_neighbors(node): return genome_graph.neighbors(node, mode="ALL")
        print_if_verbose("calculating number of neighbors")
        positions_df["number_neighbors"] = positions_df.graph_node.apply(get_neighbors).apply(len)

        # add whether it is a telomere
        positions_df["is_telomere"] = positions_df.number_neighbors<=1

        # get the component_ID, which is a unique identifier that all the components that are united
        print_if_verbose("calculating subcomponents")
        subcomponents_list = get_graph_subcomponents_undirected_graph(genome_graph)
        print_if_verbose("There are %i graph subcomponents"%len(subcomponents_list))

        node_to_componentID = pd.concat([pd.Series([Icomp]*len(component_nodes), index=component_nodes) for Icomp, component_nodes in enumerate(subcomponents_list)])
        positions_df["component_ID"] = positions_df.graph_node.map(node_to_componentID)
        if any(pd.isna(positions_df.component_ID)): raise ValueError("all positions should have a component")

        # map each graph subcomponent to the telomeric nodes
        component_ID_to_telomericNodes = dict(positions_df.groupby("component_ID").apply(lambda df_comp: df_comp[df_comp.is_telomere].graph_node.values))

        # calculate the distance from each position to the 
        print_if_verbose("calculating shortest distance to the telomere")
        node_to_distanceToTelomere = pd.Series()
        for compID, telomericNodes in component_ID_to_telomericNodes.items():

            # get the positions of the component
            component_nodes = positions_df[positions_df.component_ID==compID].graph_node.values

            # if there are no telomericNodes, it means that there is a circular chromosome. Pick the first and last parts instead
            if len(telomericNodes)==0: telomericNodes = np.array([component_nodes[0], component_nodes[-1]])

            # get a df with the distance between each node of the component (cols) and the telomeric nodes (index)
            shortestPath_lengths_df = pd.DataFrame(genome_graph.shortest_paths(source=telomericNodes, target=component_nodes, mode="ALL"), columns=component_nodes, index=telomericNodes) # this may include weights

            # get a series that has the minimum distance to the telomere for each node
            distance_telomere_series = shortestPath_lengths_df.apply(min, axis=0).apply(int)

            node_to_distanceToTelomere = node_to_distanceToTelomere.append(distance_telomere_series)

        # add to the df
        positions_df["graph_shortest_distance_to_telomere"] = positions_df.graph_node.map(node_to_distanceToTelomere)
        if any(pd.isna(positions_df.graph_shortest_distance_to_telomere)): raise ValueError("all positions should have a distance to the telomere")

        ######################################################################

        # save 
        save_df_as_tab(positions_df[["chromosome", "position", "graph_node", "number_neighbors", "is_telomere", "component_ID", "graph_shortest_distance_to_telomere"]], positions_df_file)

        # save the graph
        save_object(genome_graph, genome_graph_file)

    print_if_verbose("loading graph genome")

    # load positions df
    positions_df = get_tab_as_df_or_empty_df(positions_df_file)

    # create genome graph
    genome_graph = load_object(genome_graph_file)

    return genome_graph, positions_df

def get_df_positions_from_genomeGraph(reference_genome, df_gridss, min_AF_homo, min_QUAL, threads, replace):

    """This function takes a reference genome and a df_gridss with the breakpoints. It gererates a graph were each position is one node. In C. glabrata, the creation of the graph took ~4Gb of RAM.
    """

    # define the files
    prefix = "%s.graph_breakpoints_minQUAL=%.3f_minAFhomo=%.3f"%(reference_genome, min_QUAL, min_AF_homo)
    genome_graph_file = "%s.graph.py"%prefix
    positions_df_file = "%s.df_positions.tab"%prefix
    componentID_to_nodes_file = "%s.componentID_to_nodes_file.py"%prefix

    if any([file_is_empty(x) for x in {genome_graph_file, positions_df_file, componentID_to_nodes_file}]) or replace is True:
        print_if_verbose("getting genome graph object for %s"%reference_genome)

        # map each chromosome to an offset
        chrom_to_lenSeq = get_chr_to_len(reference_genome)

        # deffine an offset for each chromosome, which is necessary to keep all the positions of the genome as independent numbers
        chrom_to_offset = {}
        current_offset = 0
        for chrom, seqLen in chrom_to_lenSeq.items():
            chrom_to_offset[chrom] = current_offset
            current_offset+=seqLen

        # create the undirected graph, meaning that all edges work in both directions
        genome_graph = igraph.Graph(directed=False)

        # add one vertex (node) for each position in the genome
        npositions = sum(chrom_to_lenSeq.values())
        genome_graph.add_vertices(npositions)

        # create a df with the positions
        positions_df = pd.DataFrame()
        for chrom, offset in chrom_to_offset.items():

            df_c = pd.DataFrame({"position" : list(range(0, chrom_to_lenSeq[chrom]))})
            df_c["chromosome"] = chrom
            df_c["offset"] = offset
            df_c["graph_node"] = df_c.position + df_c.offset

            positions_df = positions_df.append(df_c[["chromosome", "position", "graph_node"]])

        # set the index
        positions_df = positions_df.set_index("graph_node", drop=False)

        # add the plus 1
        positions_df["position_1based"] = positions_df.position + 1

        # define the nodes that are the end of chromosomes
        chromosome_start_nodes = {offset for chrom, offset in chrom_to_offset.items()}
        chromosome_end_nodes = {(offset + chrom_to_lenSeq[chrom] - 1) for chrom, offset in chrom_to_offset.items()}
        all_nodes = set(range(npositions))

        # initialize a series with the edges of the adjacent positions
        all_edges = set(map(lambda pos: (pos, pos+1), all_nodes.difference(chromosome_end_nodes)))

        ######### MODIFY EDGES ACCORDING TO df_gridss #########

        if df_gridss is not None:
  
            # filter the df_gridss and check that it is unique
            df_gridss = df_gridss[df_gridss.QUAL>=min_QUAL]
            if len(set(df_gridss.ID))!=len(df_gridss): raise ValueError("df_gridss.ID should be unique")

            if len(df_gridss)>0:

                # keep only bends that have mates
                eventID_to_nBends = df_gridss.groupby("eventID_as_clove").apply(len)
                if any(eventID_to_nBends>2): raise ValueError("there can't be events with more than 2 breakend")
                eventID_to_nBends = dict(eventID_to_nBends)
                df_gridss["n_bends"] = df_gridss.eventID_as_clove.apply(lambda x: eventID_to_nBends[x])
                df_gridss = df_gridss[df_gridss.n_bends==2]

                if len(df_gridss)>0:

                    # add whether the breakend is homozygous. For homozygous breakends the adjacent edges will be removed
                    df_gridss["is_homo_bend"] = df_gridss.real_AF>=min_AF_homo

                    # keep only important fields
                    df_gridss = df_gridss[["#CHROM", "ALT", "POS", "is_homo_bend", "ID", "INFO_MATEID", "real_AF"]]

                    # add the node position 
                    df_gridss = df_gridss.merge(positions_df[["chromosome", "position_1based", "graph_node"]], left_on=["#CHROM", "POS"], right_on=["chromosome", "position_1based"], validate="many_to_one", how="left")

                    # add the node of the mate
                    ID_to_position = dict(df_gridss.set_index("ID")["graph_node"])
                    df_gridss["mate_graph_node"] = df_gridss.INFO_MATEID.apply(lambda x: ID_to_position[x])

                    # add all the edges related to the breakpoints
                    def get_df_gridss_sorted_edge(r): return sorted([r["graph_node"], r["mate_graph_node"]])
                    df_gridss["sorted_edge"] = df_gridss.apply(get_df_gridss_sorted_edge, axis=1).apply(tuple)
                    all_edges.update(set(df_gridss.sorted_edge))

                    # for homozygous breakends, remove the corresponding adjacent edges. Make sure that these are no adjacent breakpoints
                    df_gridss["distance_between_nodes"] = (df_gridss.graph_node - df_gridss.mate_graph_node).apply(abs)
                    df_gridss_homo = df_gridss[(df_gridss.is_homo_bend) & (df_gridss.distance_between_nodes>1)]

                    if len(df_gridss_homo)>0:

                        # add the compromised reference breakpoint
                        df_gridss_homo["compromised_ref_breakpoint"] = df_gridss_homo.apply(get_compromised_ref_breakpoint, axis=1)

                        all_ref_breakpoints_to_remove = set(df_gridss_homo.compromised_ref_breakpoint)

                        # remove the all_ref_breakpoints_to_remove
                        all_edges = all_edges.difference(all_ref_breakpoints_to_remove)

        #######################################################

        # add the edges to the graph
        genome_graph.add_edges(all_edges)
        print_if_verbose("genome graph got")
 
        ########## ADD METADATA TO THE POSITIONS DF ##########

        # add the number of neighbors that each position has
        def get_neighbors(node): return genome_graph.neighbors(node, mode="ALL")
        print_if_verbose("calculating number of neighbors")
        positions_df["number_neighbors"] = positions_df.graph_node.apply(get_neighbors).apply(len)

        # add whether it is a telomere
        positions_df["is_telomere"] = positions_df.number_neighbors<=1

        # get the component_ID, which is a unique identifier that all the components that are united
        print_if_verbose("calculating subcomponents")
        subcomponents_list = get_graph_subcomponents_undirected_graph(genome_graph)
        print_if_verbose("There are %i graph subcomponents"%len(subcomponents_list))

        node_to_componentID = pd.concat([pd.Series([Icomp]*len(component_nodes), index=component_nodes) for Icomp, component_nodes in enumerate(subcomponents_list)])
        positions_df["component_ID"] = positions_df.graph_node.map(node_to_componentID)
        if any(pd.isna(positions_df.component_ID)): raise ValueError("all positions should have a component")

        # map each componentID to the nodes
        df_components = pd.DataFrame({"componentID":node_to_componentID})
        df_components["node"] = df_components.index

        def get_nodes_from_componentID_df(df_c): return set(df_c["node"])
        componentID_to_nodes = df_components.groupby("componentID").apply(get_nodes_from_componentID_df)

        ######################################################################

        # save 
        save_df_as_tab(positions_df[["chromosome", "position", "graph_node", "number_neighbors", "is_telomere", "component_ID"]], positions_df_file)

        # save the componentID_to_nodes
        save_object(componentID_to_nodes, componentID_to_nodes_file)

        # save the graph
        save_object(genome_graph, genome_graph_file)

    print_if_verbose("loading graph genome")

    # load positions df
    positions_df = get_tab_as_df_or_empty_df(positions_df_file)

    # create genome graph
    genome_graph = load_object(genome_graph_file)

    # load componentID to nodes
    componentID_to_nodes = load_object(componentID_to_nodes_file)

    return genome_graph, positions_df, positions_df_file, componentID_to_nodes


def get_positions_df_with_graph_distance_to_telomere(positions_df, outfile, genome_graph, componentID_to_nodes, replace):

    """This function returns positions_df with a field 'graph_shortest_distance_to_telomere' . All the telomere positions must be in the df in order of this to work properly"""

    if file_is_empty(outfile) or replace is True:

        # map each graph subcomponent to the telomeric nodes
        component_ID_to_telomericNodes = dict(positions_df.groupby("component_ID").apply(lambda df_comp: df_comp[df_comp.is_telomere].graph_node.values))

        # plot the graph
        graph_plotfile = "%s.plot.png"%outfile
        print_if_verbose("plotting graph into %s"%graph_plotfile)

        genome_graph.vs['label'] = ""
        igraph.plot(genome_graph, graph_plotfile, vertex_size=1)

        lkjad


        for comp in get_graph_subcomponents_undirected_graph(genome_graph): print("\n\n", comp)

        print(componentID_to_nodes)


        # calculate the distance from each position to the 
        print_if_verbose("calculating shortest distance to the telomere")
        node_to_distanceToTelomere = pd.Series()
        for compID, telomericNodes in component_ID_to_telomericNodes.items():
            print_if_verbose("component %i/%i"%(compID+1, len(component_ID_to_telomericNodes)))

            print(componentID_to_nodes[compID])

            # get a subgraph with the nodes of the component
            subgraph_vs = list(map(lambda n: genome_graph.vs(name=n)[0].index, componentID_to_nodes[compID]))

            print(subgraph_vs)

            kjasdjkdhjdhjkads

            genome_subgraph = genome_graph.subgraph(subgraph_vs, implementation="auto")

            print(genome_graph)



            # check that he
            all_nodes_subraph = set(genome_subgraph.vs.indices)
            missing_nodes = componentID_to_nodes[compID].difference(all_nodes_subraph)
            if len(missing_nodes)>0: 

                all_nodes_graph = set(genome_graph.vs.indices)
                missing_nodes_in_genome_graph = componentID_to_nodes[compID].difference(all_nodes_graph)

                print(missing_nodes_in_genome_graph)

                raise ValueError("There are missing nodes: %s"%missing_nodes)

            print(len(get_graph_subcomponents_undirected_graph(genome_subgraph)))


            hjdahghgdahjg


            # get the positions of the component
            component_nodes = positions_df[positions_df.component_ID==compID].graph_node.values

            # if there are no telomericNodes, it means that there is a circular chromosome. Pick the first and last parts instead
            if len(telomericNodes)==0: telomericNodes = np.array([component_nodes[0], component_nodes[-1]])

            # get a df with the distance between each node of the component (cols) and the telomeric nodes (index)
            shortestPath_lengths_df = pd.DataFrame(genome_subgraph.shortest_paths(source=telomericNodes, target=component_nodes, mode="ALL"), columns=component_nodes, index=telomericNodes) # this may include weights

            # get a series that has the minimum distance to the telomere for each node
            distance_telomere_series = shortestPath_lengths_df.apply(min, axis=0).apply(int)

            node_to_distanceToTelomere = node_to_distanceToTelomere.append(distance_telomere_series)

        # add to the df
        positions_df["graph_shortest_distance_to_telomere"] = positions_df.graph_node.map(node_to_distanceToTelomere)
        if any(pd.isna(positions_df.graph_shortest_distance_to_telomere)): raise ValueError("all positions should have a distance to the telomere")

        # save
        save_df_as_tab(positions_df, outfile)

    # load
    positions_df = get_tab_as_df_or_empty_df(outfile)

    return positions_df




def get_positions_df_with_graph_distance_to_telomere(positions_df, outfile, genome_graph, componentID_to_nodes, replace):

    """This function returns positions_df with a field 'graph_shortest_distance_to_telomere' . All the telomere positions must be in the df in order of this to work properly"""

    if file_is_empty(outfile) or replace is True:

        # map each graph subcomponent to the telomeric nodes
        component_ID_to_telomericNodes = dict(positions_df.groupby("component_ID").apply(lambda df_comp: df_comp[df_comp.is_telomere].graph_node.values))

        # plot the graph
        graph_plotfile = "%s.plot.png"%outfile
        print_if_verbose("plotting graph into %s"%graph_plotfile)

        genome_graph.vs['label'] = ""
        igraph.plot(genome_graph, graph_plotfile, vertex_size=1)

        lkjad


        for comp in get_graph_subcomponents_undirected_graph(genome_graph): print("\n\n", comp)

        print(componentID_to_nodes)


        # calculate the distance from each position to the 
        print_if_verbose("calculating shortest distance to the telomere")
        node_to_distanceToTelomere = pd.Series()
        for compID, telomericNodes in component_ID_to_telomericNodes.items():
            print_if_verbose("component %i/%i"%(compID+1, len(component_ID_to_telomericNodes)))

            print(componentID_to_nodes[compID])

            # get a subgraph with the nodes of the component
            subgraph_vs = list(map(lambda n: genome_graph.vs(name=n)[0].index, componentID_to_nodes[compID]))

            print(subgraph_vs)

            kjasdjkdhjdhjkads

            genome_subgraph = genome_graph.subgraph(subgraph_vs, implementation="auto")

            print(genome_graph)



            # check that he
            all_nodes_subraph = set(genome_subgraph.vs.indices)
            missing_nodes = componentID_to_nodes[compID].difference(all_nodes_subraph)
            if len(missing_nodes)>0: 

                all_nodes_graph = set(genome_graph.vs.indices)
                missing_nodes_in_genome_graph = componentID_to_nodes[compID].difference(all_nodes_graph)

                print(missing_nodes_in_genome_graph)

                raise ValueError("There are missing nodes: %s"%missing_nodes)

            print(len(get_graph_subcomponents_undirected_graph(genome_subgraph)))


            hjdahghgdahjg


            # get the positions of the component
            component_nodes = positions_df[positions_df.component_ID==compID].graph_node.values

            # if there are no telomericNodes, it means that there is a circular chromosome. Pick the first and last parts instead
            if len(telomericNodes)==0: telomericNodes = np.array([component_nodes[0], component_nodes[-1]])

            # get a df with the distance between each node of the component (cols) and the telomeric nodes (index)
            shortestPath_lengths_df = pd.DataFrame(genome_subgraph.shortest_paths(source=telomericNodes, target=component_nodes, mode="ALL"), columns=component_nodes, index=telomericNodes) # this may include weights

            # get a series that has the minimum distance to the telomere for each node
            distance_telomere_series = shortestPath_lengths_df.apply(min, axis=0).apply(int)

            node_to_distanceToTelomere = node_to_distanceToTelomere.append(distance_telomere_series)

        # add to the df
        positions_df["graph_shortest_distance_to_telomere"] = positions_df.graph_node.map(node_to_distanceToTelomere)
        if any(pd.isna(positions_df.graph_shortest_distance_to_telomere)): raise ValueError("all positions should have a distance to the telomere")

        # save
        save_df_as_tab(positions_df, outfile)

    # load
    positions_df = get_tab_as_df_or_empty_df(outfile)

    return positions_df


def get_LOWESS_benchmarking_series_CV(kfold, frac, it, unique_df, xfield, yfield, min_test_points_CV):

    """This function takes a df with xfield and yfield. It runs kfold cross validation and returns a series with the accuracies """


    print_if_verbose(frac, it)

    # define all the indices
    all_idx = set(unique_df.index)
    
    # init rsquares
    rsquares_cv = []

    # iterate through 10-fold chross validation 
    for cvID in range(kfold):

        # get the idx test as each kfolth position after cvID
        test_idx = set(unique_df.iloc[cvID:].iloc[::kfold, :].index)
        train_idx = all_idx.difference(test_idx)

        # if there are not enough data points, break
        if len(test_idx)<min_test_points_CV: break

        # get dfs
        df_train = unique_df.loc[train_idx].sort_values(by=[xfield, yfield])
        df_test = unique_df.loc[test_idx].sort_values(by=[xfield, yfield])
       
        # get train alues
        xtrain = df_train[xfield].values
        ytrain = df_train[yfield].values

        # get the lowess fit on the train data
        lowess_results_train = statsmodels_api.nonparametric.lowess(endog=ytrain, exog=xtrain, frac=frac, it=it, xvals=None, is_sorted=True, missing="raise", return_sorted=True) # return sorted returns values storted by exog

        # unpack
        lowess_train_x = lowess_results_train[:,0]
        lowess_train_y = lowess_results_train[:,1]

        # debug if there was no LOWESS interpolation
        if any(pd.isna(lowess_train_y)): raise ValueError("the output of lowess can't have NaNs")

        # generate a linear interpolation function between the train results. It will only work for values in the range of the train
        interpolation_function = interp1d(lowess_train_x, lowess_train_y, bounds_error=True, kind="linear", assume_sorted=True)

        # get the test values. Only those where the x is in the range of the train
        xtest = df_test[xfield].values
        ytest = df_test[yfield].values

        idx_correct_test = (xtest>=min(lowess_train_x)) & (xtest<=max(lowess_train_x))
        xtest = xtest[idx_correct_test]
        ytest = ytest[idx_correct_test]

        # get the predicted y test by linear interpolation
        ytest_predicted = interpolation_function(xtest)

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
    benchmarking_series = pd.Series({"frac":frac, "it":int(it), "mean_rsquare":mean_rsquare, "inverse_std_rsquare":inverse_std_rsquare, "std_rsquare":std, "kfold":kfold})


    return benchmarking_series
    

def get_LOWESS_benchmarking_series_CV(kfold, frac, it, unique_df, xfield, yfield, min_test_points_CV):

    """This function takes a df with xfield and yfield. It runs kfold cross validation and returns a series with the accuracies """

    # define all the indices
    all_idx = set(unique_df.index)

    # init rsquares
    rsquares_cv = []

    # this will only work if the unique df is long enough
    if len(unique_df)>kfold: 

        # iterate through 10-fold chross validation 
        kfold_object = KFold(n_splits=kfold, random_state=1, shuffle=True)
        for numeric_train_index, numeric_test_index in kfold_object.split(unique_df.index):

            # get the idx test as the index of unique_df
            test_idx = set(unique_df.iloc[numeric_test_index].index)
            train_idx = all_idx.difference(test_idx)

            # if there are not enough data points, break
            if len(test_idx)<min_test_points_CV: break

            # get dfs
            df_train = unique_df.loc[train_idx].sort_values(by=[xfield, yfield])
            df_test = unique_df.loc[test_idx].sort_values(by=[xfield, yfield])
           
            # get train alues
            xtrain = df_train[xfield].values
            ytrain = df_train[yfield].values

            # get the lowess fit on the train data
            lowess_results_train = statsmodels_api.nonparametric.lowess(endog=ytrain, exog=xtrain, frac=frac, it=it, xvals=None, is_sorted=True, missing="raise", return_sorted=True) # return sorted returns values storted by exog

            # unpack
            lowess_train_x = lowess_results_train[:,0]
            lowess_train_y = lowess_results_train[:,1]

            # debug if there was no LOWESS interpolation
            if any(pd.isna(lowess_train_y)): raise ValueError("the output of lowess can't have NaNs")

            # generate a linear interpolation function between the train results. It will only work for values in the range of the train

            pyloess_fun
            
            interpolation_function = interp1d(lowess_train_x, lowess_train_y, bounds_error=True, kind="linear", assume_sorted=True)

            # get the test values. Only those where the x is in the range of the train
            xtest = df_test[xfield].values
            ytest = df_test[yfield].values

            idx_correct_test = (xtest>=min(lowess_train_x)) & (xtest<=max(lowess_train_x))
            xtest = xtest[idx_correct_test]
            ytest = ytest[idx_correct_test]

            # get the predicted y test by linear interpolation
            ytest_predicted = interpolation_function(xtest)

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


    print_if_verbose(frac, it, mean_rsquare)


    # get the final series
    benchmarking_series = pd.Series({"frac":frac, "it":int(it), "mean_rsquare":mean_rsquare, "inverse_std_rsquare":inverse_std_rsquare, "std_rsquare":std, "kfold":kfold})


    return benchmarking_series
    

def get_y_corrected_by_x_LOWESS_crossValidation(df, xfield, yfield, outdir, threads, replace, min_test_points_CV=10):

    """This function takes an x and a y series, returning the y corrected by x. This y corrected is y/(y predicted from LOWESS from x). The index must be unique. The best parameters are taken with 10 fold cross validation"""

    make_folder(outdir)

    # keep
    df = cp.deepcopy(df)[[xfield, yfield]]

    # check that the index is unique
    if len(df.index)!=len(set(df.index)): raise ValueError("The index should be unique")

    # get the df that has the collapsed xfields. lowess does not work if there are repeated xfields
    def get_median_for_xfield_df(df_x): return np.median(df_x[yfield])
    x_to_y_series = df.groupby(xfield).apply(get_median_for_xfield_df)
    unique_df = pd.DataFrame({xfield : x_to_y_series.index, yfield : x_to_y_series.values })

    # discard yfields that are 0
    unique_df = unique_df[unique_df[yfield]>0]

    # check
    if any(unique_df[yfield]<=0): raise ValueError("There can't be any 0 values in the unique df")

    # sort by the x
    unique_df = unique_df.sort_values(by=[xfield, yfield])

    ########## GET THE DF BENCHMARKED DF 10xCV ########## 

    # define the df_benckmarking file
    df_benchmarking_file = "%s/df_benckmarking.tab"%outdir

    if file_is_empty(df_benchmarking_file) or replace is True:
        print_if_verbose("getting benchmarking for %s vs %s"%(xfield, yfield))

        # define parms
        n_frac = 5
        n_it = 4
        kfold = 10

        # define all the fractions
        min_frac = min([1/len(unique_df), 0.05])
        all_fractions = list(np.linspace(min_frac, 0.1, n_frac)) + list(np.linspace(0.1, 0.5, n_frac))

        # debug
        if any(pd.isna(unique_df[xfield])) or any(pd.isna(unique_df[yfield])): raise ValueError("There are NaNs")

        # define the inputs of the benchmarking function
        inputs_fn = make_flat_listOflists([[(kfold, frac, it, unique_df, xfield, yfield, min_test_points_CV) for it in range(0, n_it+1)] for frac in all_fractions])

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

    else:

        # get the fit data
        print_if_verbose("performing LOWESS regression with best parameters for %s vs %s"%(xfield, yfield))

        # get sorted df
        unique_df = unique_df.sort_values(by=[xfield, yfield])

        # sort df benchmarking to get the max rsquare and minimum std
        max_kfold = max(df_benchmarking.kfold)
        df_benchmarking = df_benchmarking[df_benchmarking.kfold==max_kfold].sort_values(by=["mean_rsquare", "inverse_std_rsquare"], ascending=False)

        # get the best parameters
        best_parms_series = df_benchmarking.iloc[0]

        # get the y predicted with the best parms
        best_frac = best_parms_series["frac"]
        best_it = int(best_parms_series["it"])
        lowess_results = statsmodels_api.nonparametric.lowess(endog=unique_df[yfield], exog=unique_df[xfield], frac=best_frac, it=best_it, xvals=None, is_sorted=True, missing="raise", return_sorted=True) # return sorted returns values storted by exog

        # define the unique_df predicted_y_va
        unique_df["predicted_yvalues"] = lowess_results[:,1]

        # debug 
        if any(pd.isna(unique_df.predicted_yvalues)): raise ValueError("there should be no NaNs in the final prediction")

        # debug if any of the predicted_yvalues is <=0
        if any(unique_df.predicted_yvalues<=0): raise ValueError("There can't be any 0 predicted yvalues")

        # define the predicted_y_values, a series that has the index of df and the values predicted by unique_df["predicted_yvalues"]
        xval_to_predicted_yval = dict(unique_df.set_index(xfield).predicted_yvalues)

        # add the values related to an xvalue that has 0 y
        missing_y0_related_xvals = set(df[df[yfield]==0][xfield]).difference(set(xval_to_predicted_yval.keys()))
        for xval in missing_y0_related_xvals:  xval_to_predicted_yval[xval] = 0.0

        # add the prediction
        df["predicted_yvalues"] = df[xfield].apply(lambda x: xval_to_predicted_yval[x])

        # calculate the final rsquare
        final_rsquare = r2_score(df[yfield], df.predicted_yvalues)
        if pd.isna(final_rsquare): raise ValueError("rsquare can't be NaN")

        ##############################

        ######### MAKE PLOTS #########

        df_plot = df.sort_values(by=[xfield, yfield])

        fig = plt.figure(figsize=(5,5))
        plt.plot(df_plot[xfield], df_plot[yfield], "o", alpha=0.2, color="gray", label="raw data")
        plt.plot(df_plot[xfield], df_plot.predicted_yvalues, "-", color="red", label="LOWESS fit")

        plt.title("Fitting LOWESS with frac=%.3f and it=%i. final R2=%.3f. %ix CV R2=%.3f +- %.3f (SD)\n"%(best_frac, best_it, final_rsquare, best_parms_series["kfold"], best_parms_series["mean_rsquare"], best_parms_series["std_rsquare"]))
        plt.legend(bbox_to_anchor=(1, 1))
        plt.xlabel(xfield)
        plt.ylim([0, np.percentile(df[yfield], 90)])
        plt.ylabel(yfield)

        fig.savefig("%s/coverage.pdf"%(outdir), bbox_inches='tight')
        plt.close(fig)

        ############################

        # get the corrected vals. If there is no prediction just return the raw vals
        def divide_with_noNaN_correction(r):

            # if the yfield is 0, return it as it is
            if r[yfield]==0: return 0.0

            # if the predicted yvalues are equal to 0 because they are in missing_y0_related_xvals, return 0
            elif r["predicted_yvalues"]==0 and r[xfield] in missing_y0_related_xvals: return 0.0

            # else
            elif r["predicted_yvalues"]==0: raise ValueError("predicted_yvalues can't be 0 if yfield is not as well") 
            
            # normal division
            else: return r[yfield]/r["predicted_yvalues"]

        if final_rsquare>0: y_corrected = df.apply(divide_with_noNaN_correction, axis=1)
        else: y_corrected = df[yfield]

        # debug
        if any(pd.isna(y_corrected)): raise ValueError("there should be no NaNs in y_corrected")

    # return in the initial order
    return y_corrected, final_rsquare

def get_lowess_fit_y(x, y, frac):

    """This function takes an x and a y and returns a numpy array that has the lowess fit according to window_size"""

    # define the window_size
    window_size_loess = int(len(x)*frac) + 1

    # this is very slow

    # get the loess fit
    pyloess_object = pyloess_fun.Loess(x, y)
    lowess_y = np.array(list(map(lambda xval: pyloess_object.estimate(xval, window=window_size_loess, use_matrix=False, degree=1), x))) 

    return lowess_y



def get_lowess_fit_y(x, y, frac, outprefix):

    """This function takes an x and a y and returns a numpy array that has the lowess fit according to frac. It runs the R loess function to do so"""

    # write a df that contains x and y
    df = pd.DataFrame({"x":x, "y":y})
    df_file = "%s_df_file.tab"%outprefix
    save_df_as_tab(df, df_file)

    # create an Rfile that contains the rscript to run loess
    Rfile = "%s_Rscript.R"%outprefix
    df_file_with_predicted = "%s_df_file_with_predicted.tab"%outprefix

    Rscript_lines = ["#!/usr/bin/env Rscript",    
                     'df = read.table("%s", sep="\\t", header=TRUE)'%df_file,
                     'loess_object = loess(y ~ x, df, span=%.4f, degree=1, parametric=FALSE, family="gaussian")'%frac,
                     'df$predicted_y = predict(loess_object)',
                     'write.table(df, file="%s", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)'%df_file_with_predicted
                    ] 

    open(Rfile, "w").write("\n".join(Rscript_lines))

    # run this Rfile
    std_loess = "%s.loess.std"%outprefix
    run_cmd("chmod 755 %s && %s > %s 2>&1"%(Rfile, Rfile, std_loess), env=EnvName_R)

    # load df
    df_with_pred = get_tab_as_df_or_empty_df(df_file_with_predicted)

    # remove unnecessary files
    for f in [df_file, Rfile, df_file_with_predicted, std_loess]: remove_file(f)

    return df_with_pred.predicted_y.apply(float).values
