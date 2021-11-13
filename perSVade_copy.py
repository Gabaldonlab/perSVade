#!/usr/bin/env python

# This is the perSVade pipeline main script, which shoul dbe run on the perSVade conda environment


##### DEFINE ENVIRONMENT #######


# import persvade-specific modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np

# packages installed into the conda environment 
samtools = "%s/bin/samtools"%EnvDir
java = "%s/bin/java"%EnvDir

# scripts that are installed under this software
varcall_cnv_pipeline = "%s/varcall_cnv_pipeline.py"%CWD

#######


#####################

########### resource allocation ##############

general_optional_args = parser.add_argument_group("GENERAL OPTIONAL ARGUMENTS")

general_optional_args.add_argument("-p", "--ploidy", dest="ploidy", default=1, type=int, help="Ploidy, can be 1 or 2")

general_optional_args.add_argument("-gff", "--gff-file", dest="gff", default=None, help="path to the GFF3 annotation of the reference genome. Make sure that the IDs are completely unique for each 'gene' tag. This is necessary for both the CNV analysis (it will look at genes there) and the annotation of the variants.")

general_optional_args.add_argument("-mcode", "--mitochondrial_code", dest="mitochondrial_code", default=3, type=int, help="The code of the NCBI mitochondrial genetic code. For yeasts it is 3. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. The information of this website may be wrong, so you may want to double check with the literature.")

general_optional_args.add_argument("-gcode", "--gDNA_code", dest="gDNA_code", default=1, type=int, help="The code of the NCBI gDNA genetic code. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . For C. albicans it is 12. The information of this website may be wrong, so you may want to double check with the literature.")






######### small variant calling #######

smallVars_args = parser.add_argument_group("SMALL VARIANT CALLING AND COVERAGE PER GENE CALCULATION")

smallVars_args.add_argument("--run_smallVarsCNV", dest="run_smallVarsCNV", action="store_true", default=False, help="Will call small variants and CNV.")

smallVars_args.add_argument("-caller", "--caller", dest="caller", required=False, default="all", help="SNP caller option to obtain vcf file. options: no/all/HaplotypeCaller/bcftools/freebayes. It can be a comma-sepparated string, like 'HaplotypeCaller,freebayes'")

smallVars_args.add_argument("-c", "--coverage", dest="coverage", default=20, type=int, help="minimum Coverage (int) of a position to be considered for small variant calling. This parameter should be related to the coverage of your library. You may be careful with setting it too low (i.e. <15) as it will yield many false positive calls. It is reasonable to check how other similar studies set this parameter.")

smallVars_args.add_argument("--minAF_smallVars", dest="minAF_smallVars", default="infer", help="The minimum fraction of reads covering a variant to be called. The default is 'infer', which will set a threshold based on the ploidy. This is only relevant for the final vcfs, where only PASS vars are considered. It can be a number between 0 and 1.")

smallVars_args.add_argument("--window_freebayes_bp", dest="window_freebayes_bp", default=10000, type=int, help="freebayes is run in parallel by chunks of the genome. This cmd specifies the window (in bp) in which freebayes regions are split to. If you increase this number the splitting will be in larger chunks of the genome.")

smallVars_args.add_argument("--pooled_sequencing", dest="pooled_sequencing", action="store_true", default=False, help="It is a pooled sequencing run, which means that the small variant calling is not done based on ploidy. If you are also running SV calling, you will run parameters optimisation on a sample that has 1-10 pooling.")

smallVars_args.add_argument("--run_ploidy2_ifHaploid", dest="run_ploidy2_ifHaploid", action="store_true", default=False, help="If ploidy==1, run also in diploid configuration. This is useful because there may be diploid mutations in duplicated regions.")

smallVars_args.add_argument("--consider_repeats_smallVarCall", dest="consider_repeats_smallVarCall", action="store_true", default=False, help="If --run_smallVarsCNV, this option will imply that each small  variant will have an annotation of whether it overlaps a repeat region.")

smallVars_args.add_argument("--remove_smallVarsCNV_nonEssentialFiles", dest="remove_smallVarsCNV_nonEssentialFiles", action="store_true", default=False, help="Will remove all the varCall files except the integrated final file and the bam file.")

#######################################





########################################
##### GENERAL PROCESSING OF INPUTS #####
########################################





#### REPLACE THE GFF ####
target_gff = "%s/reference_genome_features.gff"%reference_genome_dir

# if you want to repeat the annotation, remove all the files in the reference genome dir that are related to the annotation
if opt.replace_var_annotation is True:
    for f in os.listdir(reference_genome_dir): 

        file_path = "%s/%s"%(reference_genome_dir, f)
        if file_path.startswith(target_gff): fun.remove_file(file_path)

# copy the gff
if opt.gff is None: print("WARNING: gff was not provided. This will be a problem if you want to annotate small variant calls")
else:

    if fun.file_is_empty(target_gff): fun.soft_link_files(opt.gff, target_gff)

    # change the path
    opt.gff = target_gff

#########################


###############################################

#### define misc args ####



# test whether the gff is correct
if opt.gff is not None: fun.check_that_gff_is_correct(opt.gff, opt.ref, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code, opt.threads, opt.replace)


####### GENERATE real_bedpe_breakpoints AROUND HOMOLOGOUS REGIONS #########


############################################################################




#####################################
##### STRUCTURAL VARIATION ##########
#####################################


###################################################################################################

#####################################
#####################################
#####################################

###########################
###### CNV CALLING ########
###########################

# run CNVcalling by getting df_CNV_coverage
minimal_CNV_fields = ["chromosome", "merged_relative_CN", "start", "end", "CNVid", "median_coverage", "median_coverage_corrected", "SVTYPE"] + ["median_relative_CN_%s"%x for x in cnv_calling_algs]

run_CNV_calls = (opt.skip_CNV_calling is False and not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]))
if run_CNV_calls is True: 
    print("RUNNING COVERAGE-BASED CNV CALLING...")

    # make folder
    cnv_calling_outdir = "%s/CNV_calling"%opt.outdir
    fun.make_folder(cnv_calling_outdir)

    # run CNV calling
    df_CNV_coverage = fun.run_CNV_calling(sorted_bam, opt.ref, cnv_calling_outdir, opt.threads, opt.replace, opt.mitochondrial_chromosome, opt.window_size_CNVcalling, opt.ploidy, bg_sorted_bam_CNV=opt.bg_sorted_bam_CNV, cnv_calling_algs=cnv_calling_algs)

else: df_CNV_coverage = pd.DataFrame(columns=minimal_CNV_fields)

###########################
###########################
###########################

#####################################
###### SV and CNV ANNOTATION ########
#####################################

start_time_SVandCNVcalling =  time.time()

run_SV_CNV_calling = (opt.skip_SVcalling is False and not any([x=="skip" for x in {opt.fastq1, opt.fastq2}]) and opt.skip_SV_CNV_calling is False)
if run_SV_CNV_calling is True:

    print("running CNV calling per window and integrating to SV calling")

    # define outdirs
    outdir_var_calling = "%s/SVcalling_output"%opt.outdir

    # remove folders if there is some replacement to be done. Remove
    if opt.replace_SV_CNVcalling is True: fun.delete_folder(outdir_var_calling)

    # stop after the removal
    if opt.StopAfter_replace_SV_CNVcalling is True: 
        print("exitting after the --replace_SV_CNVcalling action")
        sys.exit(0)

    # make folder
    fun.make_folder(outdir_var_calling)
    
    # get the variant calling 
    SV_CNV_vcf = fun.get_vcf_all_SVs_and_CNV(opt.outdir, outdir_var_calling, sorted_bam, opt.ref, opt.ploidy, df_CNV_coverage, opt.window_size_CNVcalling, cnv_calling_algs, replace=opt.replace, threads=opt.threads, mitochondrial_chromosome=opt.mitochondrial_chromosome)

    print("the SV and CNV calling vcf can be found in %s"%SV_CNV_vcf)

    # get variant annotation
    if opt.gff is not None:

        # remove the annotated vcf if needed 
        if opt.replace_var_annotation is True: 
            fun.remove_file("%s_annotated_VEP.tab"%SV_CNV_vcf)
            fun.remove_file("%s_annotated_VEP.tab.raw.tbl.tmp_summary.html"%SV_CNV_vcf)

        print("annotating SV, CNV variants with VEP")
        SV_CNV_vcf_annotated = fun.annotate_SVs_inHouse(SV_CNV_vcf, gff_with_biotype, opt.ref, replace=opt.replace, threads=opt.threads, mitochondrial_chromosome=opt.mitochondrial_chromosome, mito_code=opt.mitochondrial_code, gDNA_code=opt.gDNA_code)

        print("annotated SV vcf can be found in %s"%SV_CNV_vcf_annotated)
    
    else: print("WARNING: Skipping SV annotation because -gff was not provided.")

end_time_SVandCNVcalling =  time.time()

#####################################
#####################################
#####################################

# stop after the generation of SV and CNV calls

#####################################
###### SMALL VARS AND CNV ###########
#####################################

start_time_smallVarsCNV =  time.time()

if opt.run_smallVarsCNV:

    # define an outdir
    outdir_varcall = "%s/smallVars_CNV_output"%opt.outdir

    # define the basic cmd
    varcall_cmd = "%s -r %s --threads %i --outdir %s -sbam %s --caller %s --coverage %i --mitochondrial_chromosome %s --mitochondrial_code %i --gDNA_code %i --minAF_smallVars %s --window_freebayes_bp %i --log_file_all_cmds %s"%(varcall_cnv_pipeline, opt.ref, opt.threads, outdir_varcall, sorted_bam, opt.caller, opt.coverage, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code, opt.minAF_smallVars, opt.window_freebayes_bp, fun.log_file_all_cmds)

    # add options
    if opt.replace is True: varcall_cmd += " --replace"
    if opt.gff is not None: varcall_cmd += " -gff %s"%opt.gff
    if opt.StopAfter_smallVarCallSimpleRunning is True: varcall_cmd += " --StopAfter_smallVarCallSimpleRunning"
    if opt.replace_var_integration is True: varcall_cmd += " --replace_var_integration"
    if opt.pooled_sequencing is True: varcall_cmd += " --pooled_sequencing"
    if opt.consider_repeats_smallVarCall is True: varcall_cmd += " --repeats_table %s"%repeats_table_file
    if opt.generate_alternative_genome is True: varcall_cmd += " --generate_alternative_genome"
    if opt.skip_cnv_per_gene_analysis is True: varcall_cmd += " --skip_cnv_analysis"

    # define which ploidies to run
    if opt.ploidy==1 and opt.run_ploidy2_ifHaploid is False: ploidies_varcall = [1]
    if opt.ploidy==1 and opt.run_ploidy2_ifHaploid is True: ploidies_varcall = [1, 2]
    else: ploidies_varcall = [opt.ploidy]

    # run for each ploidy
    for ploidy_varcall in ploidies_varcall:

        # remove the annotation-derived outputs
        if opt.replace_var_annotation is True: 
            fun.remove_file("%s/variant_annotation_ploidy%i.tab"%(outdir_varcall, ploidy_varcall))
            fun.delete_folder("%s/CNV_results"%outdir_varcall)
            
        # run the variant calling command
        varcall_cmd += " -p %i"%ploidy_varcall
        if __name__ == '__main__': fun.run_cmd(varcall_cmd)

        # regenerate the variant calling file according to run_SV_CNV_calling
        if run_CNV_calls is True: fun.get_small_variant_calling_withCNstate("%s/variant_calling_ploidy%i.tab"%(outdir_varcall, ploidy_varcall), df_CNV_coverage, replace=(opt.replace or opt.replace_addingCNstate_to_smallVars))
  
    # define the small variants vcf
    small_vars_vcf = "%s/variants_atLeast1PASS_ploidy%i.vcf"%(outdir_varcall, opt.ploidy)

    # define the variant annotation
    small_vars_var_annotation = "%s/variant_annotation_ploidy%i.tab"%(outdir_varcall, opt.ploidy)

    # clean the varcall dir if specified
    if opt.remove_smallVarsCNV_nonEssentialFiles is True: fun.remove_smallVarsCNV_nonEssentialFiles_severalPloidies(outdir_varcall, ploidies_varcall)

end_time_smallVarsCNV =  time.time()

if opt.StopAfter_smallVarCall is True:
    print("WARNING: Ending after the running of small variant calling...")
    sys.exit(0)

#####################################
#####################################
#####################################


#####################################
#####################################

# keep the simulation files
if opt.keep_simulation_files is True: fun.keep_simulation_files_for_perSVade_outdir(opt.outdir, replace=opt.replace, n_simulated_genomes=opt.nsimulations, simulation_ploidies=simulation_ploidies)

# at the end you want to clean the outdir to keep only the essential files
if opt.skip_cleaning_outdir is False: fun.clean_perSVade_outdir(opt.outdir)
