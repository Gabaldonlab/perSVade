#!/usr/bin/env python

# This is a pipeline to call small variants and CNV

import argparse, os
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter
import copy as cp
import pickle
import string
import shutil 
from Bio import SeqIO
import random
import sys
from shutil import copyfile

######################################################
###############  DEFINE ENVIRONMENT ##################
######################################################

# get the cwd were all the scripts are 
CWD = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, CWD)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import sv_functions as fun

# packages installed into the conda environment 
samtools = "%s/bin/samtools"%EnvDir
freebayes = "%s/bin/freebayes"%EnvDir
bcftools = "%s/bin/bcftools"%EnvDir
vcffilter = "%s/bin/vcffilter"%EnvDir
gatk = "%s/bin/gatk"%EnvDir # this is gatk4, and it has to be installed like that
vcfallelicprimitives = "%s/bin/vcfallelicprimitives"%EnvDir
qualimap = "%s/bin/qualimap"%EnvDir
sift4g = "%s/bin/sift4g"%EnvDir
java = "%s/bin/java"%EnvDir
bcftools = "%s/bin/bcftools"%EnvDir
bgzip = "%s/bin/bgzip"%EnvDir
tabix = "%s/bin/tabix"%EnvDir
bedtools = "%s/bin/bedtools"%EnvDir
picard = "%s/share/picard-2.18.26-0/picard.jar"%EnvDir


# scripts installed with perSVade
run_vep = "%s/run_vep.py"%CWD


description = """
Run small variant calling and CNV analysis.
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# general args
parser.add_argument("-r", "--ref", dest="ref", required=True, help="Reference genome. Has to end with .fasta. This pipeline will create files under this fasta.")
parser.add_argument("-thr", "--threads", dest="threads", default=2, type=int, help="Number of threads, Default: 16")
parser.add_argument("-o", "--outdir", dest="outdir", action="store", required=True, help="Directory where the data will be stored")
parser.add_argument("--replace", dest="replace", action="store_true", help="Replace existing files")
parser.add_argument("--replace_vep_integration", dest="replace_vep_integration", action="store_true", help="Replace existing files of the merging of the VEP output and the vcf files. This is for debugging purposes")
parser.add_argument("-p", "--ploidy", dest="ploidy", default=1, type=int, help="Ploidy, can be 1 or 2")

# alignment args
parser.add_argument("-sbam", "--sortedbam", dest="sortedbam", required=True, type=str, help="The path to the sorted bam file, which should have a bam.bai file in the same dir. This is mutually exclusive with providing reads")

# variant calling args
parser.add_argument("-caller", "--caller", dest="caller", required=False, default="all", help="SNP caller option to obtain vcf file. options: no/all/HaplotypeCaller/bcftools/freebayes.")
parser.add_argument("-c", "--coverage", dest="coverage", default=20, type=int, help="minimum Coverage (int)")
parser.add_argument("-mchr", "--mitochondrial_chromosome", dest="mitochondrial_chromosome", default="mito_C_glabrata_CBS138", type=str, help="The name of the mitochondrial chromosome. This is important if you have mitochondrial proteins for which to annotate the impact of nonsynonymous variants, as the mitochondrial genetic code is different. This should be the same as in the gff. If there is no mitochondria just put no_mitochondria")
parser.add_argument("-mcode", "--mitochondrial_code", dest="mitochondrial_code", default=3, type=int, help="The code of the NCBI mitochondrial genetic code. For yeasts it is 3. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")
parser.add_argument("-gcode", "--gDNA_code", dest="gDNA_code", default=1, type=int, help="The code of the NCBI gDNA genetic code. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . For C. albicans it is 12. ")

# CNV args
parser.add_argument("--skip_cnv_analysis", dest="skip_cnv_analysis", action="store_true", default=False, help="Skipp the running of the CNV pipeline, which outputs the number of copies that each gene has according to coverage. The gene ID's are based on the GFF3 files that have been provided in -gff")

# perform the vcf integration
parser.add_argument("--get_merged_vcf", dest="get_merged_vcf", action="store_true", default=False, help="Get the integrated vcf")

# avoid marking duplicates
parser.add_argument("--skip_MarkingDuplicates", dest="skip_MarkingDuplicates", action="store_true", default=False, help="Skips the marking of duplicates in the bam.")


# othe args
parser.add_argument("-gff", "--gff-file", dest="gff", default=None, help="path to the GFF3 annotation of the reference genome. Make sure that the IDs are completely unique for each 'gene' tag. This is necessary for both the CNV analysis (it will look at genes there) and the annotation of the variants.")

# removing args
parser.add_argument("--remove_smallVarsCNV_nonEssentialFiles", dest="remove_smallVarsCNV_nonEssentialFiles", action="store_true", default=False, help="Will remove all the varCall files except the integrated final file, the filtered and normalised vcfs, the raw vcf and the CNV files.")

# stopping options
parser.add_argument("--StopAfter_smallVarCallSimpleRunning", dest="StopAfter_smallVarCallSimpleRunning", action="store_true", default=False, help="Stop after obtaining the filtered vcf outputs of each program.")


# get arguments
opt = parser.parse_args()

######################################################
######################################################
######################################################

# debug commands
fun.make_folder(opt.outdir)
if not opt.gff is None and fun.file_is_empty(opt.gff): raise ValueError("%s is not a valid gff"%opt.gff)

print("running small vars and CNV pipeline into %s"%opt.outdir)


# check that the environment is correct
fun.run_cmd("echo 'This is a check of the environment in which the pipeline is running'; which bedtools")

# correct the gff file, so that it doesn't have lines starting with # and also add biotype (important for ensembl VEP)
if not opt.gff is None:
    correct_gff = "%s_corrected.gff"%(opt.gff); correct_gff_tmp = "%s_corrected_tmp.gff"%(opt.gff)

    if fun.file_is_empty(correct_gff) or opt.replace is True:
        print("correcting gff")
        correct_gff_cmd = "grep -v '^#' %s > %s"%(opt.gff, correct_gff_tmp); fun.run_cmd(correct_gff_cmd)
        os.rename(correct_gff_tmp, correct_gff)

    # modify gff to add biotype
    gff_with_biotype = "%s_with_biotype.gff"%correct_gff
    if fun.file_is_empty(gff_with_biotype) or opt.replace is True:
        print("adding biotype")

        starting_lines = [line for line in open(correct_gff, "r") if line.startswith("#")]
        df_gff3 = pd.read_csv(correct_gff, skiprows=list(range(len(starting_lines))), sep="\t", names=["chromosome", "source", "type_feature", "start", "end", "score", "strand", "phase", "attributes"])

        def add_biotype(row):
            if "biotype" not in row["attributes"] and "gene_biotype" not in row["attributes"]: row["attributes"] += ";biotype=%s"%row["type_feature"]
            return row["attributes"]

        # add biotype and save
        df_gff3["attributes"] = df_gff3.apply(lambda row: add_biotype(row), axis=1)
        df_gff3.to_csv(gff_with_biotype, sep="\t", index=False, header=False)

# First create some files that are important for any program

# Create a reference dictionary
rstrip = opt.ref.split(".")[-1]
dictionary = "%sdict"%(opt.ref.rstrip(rstrip)); tmp_dictionary = "%s.tmp"%dictionary; 
if fun.file_is_empty(dictionary) or opt.replace is True:

    # remove any previously created tmp_file
    if not fun.file_is_empty(tmp_dictionary): os.unlink(tmp_dictionary)

    print("Creating picard dictionary")
    cmd_dict = "%s -jar %s CreateSequenceDictionary R=%s O=%s TRUNCATE_NAMES_AT_WHITESPACE=true"%(java, picard, opt.ref, tmp_dictionary); fun.run_cmd(cmd_dict)   
    os.rename(tmp_dictionary , dictionary)

# Index the reference
if fun.file_is_empty("%s.fai"%opt.ref) or opt.replace is True:
    print ("Indexing the reference...")
    cmd_indexRef = "%s faidx %s"%(samtools, opt.ref); fun.run_cmd(cmd_indexRef) # This creates a .bai file of the reference


# marking duplicates or not
if opt.skip_MarkingDuplicates is False: sorted_bam = fun.get_sortedBam_with_duplicatesMarked(opt.sortedbam, threads=opt.threads, replace=opt.replace)

else: sorted_bam = opt.sortedbam
print("running VarCall for %s"%sorted_bam)

#####################################
############### CNV #################
#####################################

if opt.skip_cnv_analysis is False:

    print("Starting CNV analysis")

    # make a folder for the CNV anlysis
    cnv_outdir = "%s/CNV_results"%opt.outdir
    if not os.path.isdir(cnv_outdir): os.mkdir(cnv_outdir)

    # get the bed file, and also the one of the regions surrounding each gene
    print(opt.ref)
    bed_file = "%s.bed_index1"%correct_gff; bed_file_regions = fun.extract_BEDofGENES_of_gff3(correct_gff, bed_file, replace=opt.replace, reference=opt.ref)

    # define the interetsing beds
    gene_to_coverage_file = "%s/gene_to_coverage_genes.tab"%cnv_outdir
    gene_to_coverage_file_regions = "%s/gene_to_coverage_regions.tab"%cnv_outdir

    # go through each region of bed file
    for bed, final_coverge_file in [(bed_file, gene_to_coverage_file), (bed_file_regions, gene_to_coverage_file_regions)]: fun.write_coverage_per_gene_mosdepth_and_parallel(sorted_bam, opt.ref, cnv_outdir, bed, final_coverge_file, replace=opt.replace)

 
    # write the integrated file
    integrated_coverage_file = "%s/genes_and_regions_coverage.tab"%cnv_outdir; integrated_coverage_file_tmp = "%s.tmp"%integrated_coverage_file
    if fun.file_is_empty(integrated_coverage_file) or opt.replace is True: 

       # integrate in one
        df_genes = pd.read_csv(gene_to_coverage_file, sep="\t")
        df_regions = pd.read_csv(gene_to_coverage_file_regions, sep="\t")
        df_integrated = df_genes.merge(df_regions, on="ID", validate="one_to_one", suffixes=("", "_+-10kb_region"))

        # write
        df_integrated.to_csv(integrated_coverage_file_tmp, sep="\t", header=True, index=False)
        os.rename(integrated_coverage_file_tmp, integrated_coverage_file)

    # remove everyhing that is not the coverage file
    for f in os.listdir(cnv_outdir): 
        if f not in {fun.get_file(gene_to_coverage_file), fun.get_file(gene_to_coverage_file_regions), fun.get_file(integrated_coverage_file)}: fun.delete_file_or_folder("%s/%s"%(cnv_outdir, f))
 
    # In Laia's script, she calculates coverage as the median reads per gene (cov per gene) / mean of the cov per gene across all genes

print("CNV analysis finished")

#####################################
######## VARIANTCALLING #############
#####################################

# initialize an array of files that have the VCF results filtered
filtered_vcf_results = []

# Go through the callers, creating in outdir a folder with the results of each
if opt.caller == "no": print("Stop. Doing the variant calling is not necessary.")
    
if opt.caller == "HaplotypeCaller" or opt.caller == "all":

    print("RUNNING GATK: HaplotypeCaller")

    # create a folder that will contain the output of VCF
    outdir_gatk = "%s/HaplotypeCaller_ploidy%i_out"%(opt.outdir, opt.ploidy)

    # run gatk and get the filtered filename
    gatk_out_filtered = fun.run_gatk_HaplotypeCaller(outdir_gatk, opt.ref, sorted_bam, opt.ploidy, opt.threads, opt.coverage, replace=opt.replace)

    # keep
    filtered_vcf_results.append(gatk_out_filtered)
    
    print("HaplotypeCaller is done")

if opt.caller == "bcftools" or opt.caller == "all":

    print("RUNNING bcftools")

    # create a folder that will contain the output of VCF
    outdir_bcftools = "%s/bcftools_ploidy%i_out"%(opt.outdir, opt.ploidy)
    if not os.path.isdir(outdir_bcftools): os.mkdir(outdir_bcftools)

    # only continue if the final file is not done
    filtered_output = "%s/output.filt.vcf"%outdir_bcftools;     
    if fun.file_is_empty(filtered_output) or opt.replace is True:

        # look for the mpileup bcf in sister directories, as it is the same for any other ploidy
        mpileup_output = "%s/output.mpileup.bcf"%outdir_bcftools; mpileup_output_tmp = "%s.tmp"%mpileup_output
        for folder in os.listdir(opt.outdir):
            if folder.startswith("bcftools_ploidy") and folder.endswith("_out"):

                # look for the potential previously calculated mpielup outputs
                potential_previosuly_calculated_mpileup_output = "%s/%s/output.mpileup.bcf"%(opt.outdir, folder)
                if not fun.file_is_empty(potential_previosuly_calculated_mpileup_output): 
                    print("taking %s from previous run"%potential_previosuly_calculated_mpileup_output)
                    mpileup_output = potential_previosuly_calculated_mpileup_output; break

        # if there is no previous run
        if fun.file_is_empty(mpileup_output) or opt.replace is True:

            print("Running mpileup...")
            cmd_bcftools_mpileup = '%s mpileup -a "AD,DP" -O b -f %s -o %s --threads %i %s'%(bcftools, opt.ref, mpileup_output_tmp, opt.threads, sorted_bam); fun.run_cmd(cmd_bcftools_mpileup)
            os.rename(mpileup_output_tmp, mpileup_output)


        # run bcftools call
        call_output = "%s/output.raw.vcf"%outdir_bcftools; call_output_tmp = "%s.tmp"%call_output
        if fun.file_is_empty(call_output) or opt.replace is True:
            print("Running bcftools call ...")

            # define the ploidy specification
            if opt.ploidy==1: ploidy_cmd = "--ploidy %i"%opt.ploidy # This is all haploid
            else:
                # create a ploidy file if ploidy is 2. There's no way to simpli specify ploidy 2
                ploidy_file_bcftools = "%s/ploidy_file.tab"%outdir_bcftools
                open(ploidy_file_bcftools, "w").write("* * * * %i\n"%opt.ploidy) # CHROM, FROM, TO, SEX, PLOIDY

                ploidy_cmd = "--ploidy-file %s"%ploidy_file_bcftools

            cmd_bcftools_call = "%s call -m -f GQ,GP -v -O v --threads %i -o %s %s %s"%(bcftools, opt.threads, call_output_tmp, ploidy_cmd, mpileup_output); fun.run_cmd(cmd_bcftools_call)

            os.rename(call_output_tmp, call_output)
      
        #As there are no recommendations for bcftools, we decided to apply exclusively the filter for coverage. To apply harder filters please edit this command!
        
        # this generates a filtered, vcf, which only has the PASS ones.
        filtered_output_tmp = "%s.tmp"%filtered_output
        if fun.file_is_empty(filtered_output) or opt.replace is True:
            print("Filtering bcftools ... ")
            cmd_filter = "%s filter -m x -e 'INFO/DP <= %i' -O v --threads %i -o %s %s"%(bcftools, opt.coverage, opt.threads, filtered_output_tmp, call_output); fun.run_cmd(cmd_filter)
            os.rename(filtered_output_tmp, filtered_output)

        # keep
    filtered_vcf_results.append(filtered_output)

    print("bcftools is done")

if opt.caller == "freebayes" or opt.caller == "all":

    print("RUNNING freebayes")

    # create a folder that will contain the output of VCF
    outdir_freebayes = "%s/freebayes_ploidy%i_out"%(opt.outdir, opt.ploidy)

    # run freebayes
    freebayes_filtered =  fun.run_freebayes_parallel(outdir_freebayes, opt.ref, sorted_bam, opt.ploidy, opt.coverage, replace=opt.replace) 

    # keep
    filtered_vcf_results.append(freebayes_filtered)
    
    print("freebayes is done")

if opt.StopAfter_smallVarCallSimpleRunning is True:
    print("stopping after generation of each variant calling")
    sys.exit(0)

##########

##############
# NORMALISATION OF THE VARIANTS. IT IS IMPORTANT TO REPRESENT THE VARIANTS IN THE PRIMITIVE FORM, AS IT ALLOWS TO COMPARE INDELS FROM SEVERAL PROGRAMS.
##############

print("Performing variant normalisation. For in/del variants, several programs may yield various variant representations. This is performed here with VCFLIB ")

# initialize an array that will keep the path to the normalised VCFS
all_normalised_vcfs = set()

# normalise all the filtered_vcf_results with vcfallelicprimitives
for unnormalised_vcf in filtered_vcf_results:

    # define the normalised output
    folder = "/".join(unnormalised_vcf.split("/")[0:-1])
    normalised_vcf = "%s/output.filt.norm_vcflib.vcf"%folder; normalised_vcf_tmp = "%s.tmp"%normalised_vcf

    # generate an unifyed representation of the vcfs
    if fun.file_is_empty(normalised_vcf) or opt.replace is True:
    #if True:

        ##### GET AS UPPERCASE #####

        print("Puting all REF and ALT alleles to uppercase")

        # load into df
        initial_lines_list = [line for line in open(unnormalised_vcf, "r", encoding='utf-8', errors='ignore') if line.startswith("##")]
        vcf_df = pd.read_csv(unnormalised_vcf, skiprows=list(range(len(initial_lines_list))), sep="\t", na_values=fun.vcf_strings_as_NaNs, keep_default_na=False)

        # put to uppercase
        vcf_df["REF"]  = vcf_df["REF"].apply(lambda x: x.upper())
        vcf_df["ALT"]  = vcf_df["ALT"].apply(lambda x: x.upper())

        # write to the same file, including the initial lines of the vcf for consistency
        open(unnormalised_vcf, "w").write("".join(initial_lines_list) + vcf_df.to_csv(sep="\t", index=False))

        #############################

        print("Running vcfallelicprimitives for vcf %s"%unnormalised_vcf)
        normalised_vcf_tmp_lines = "%s.lines.vcf"%normalised_vcf_tmp
        cmd_normalise = "%s --keep-geno %s > %s"%(vcfallelicprimitives, unnormalised_vcf, normalised_vcf_tmp_lines); fun.run_cmd(cmd_normalise)

        # get the header. vcfallelicprimitives removes the header info
        header_lines = "%s.header"%unnormalised_vcf
        fun.run_cmd("grep '^##' %s > %s"%(unnormalised_vcf, header_lines))

        # add header
        fun.run_cmd("cat %s %s > %s"%(header_lines, normalised_vcf_tmp_lines, normalised_vcf_tmp))

        # remove unnecessary files
        for f in [header_lines, normalised_vcf_tmp_lines]: fun.remove_file(f)

        os.rename(normalised_vcf_tmp, normalised_vcf)

    # keep
    all_normalised_vcfs.add(normalised_vcf)

print("VCFLIB Normalisation is done")


###################################
##### GET THE INTEGRATED VARS ##### 
################################### 

# merge the variants
if opt.get_merged_vcf is True:

    # get the merged vcf records (these are multiallelic)
    print("getting merged vcf without multialleles")
    # merged_vcf_all, merged_vcf_onlyPASS 
    merged_vcf_all, merged_vcf_onlyPASS = fun.merge_several_vcfsSameSample_into_oneMultiSample_vcf(filtered_vcf_results, opt.ref, opt.outdir, replace=opt.replace, threads=opt.threads)

    jadhjdkah

    # split the multiallelic records for each of them
    print("splitting multialleles and bgzipping")
    merged_vcf_all_noMultialleles_gz = fun.get_normed_bgzip_and_tabix_vcf_file(merged_vcf_all, opt.ref, replace=opt.replace, threads=opt.threads, multiallelics_cmd="-any")
    merged_vcf_onlyPASS_noMultialleles_gz = fun.get_normed_bgzip_and_tabix_vcf_file(merged_vcf_onlyPASS, opt.ref, replace=opt.replace, threads=opt.threads, multiallelics_cmd="-any")

    print(merged_vcf_all_noMultialleles_gz, merged_vcf_onlyPASS_noMultialleles_gz)






    adkghdakhgdag

###################################
###################################
###################################


############################
# ANNOTATE VARIANTS WITH VEP
############################

if opt.gff is None: print("No gff provided. Skipping the annotation AND integration of the variants")

else:
    print("getting vep annotations")

    # create a dictionary with [typeNormalisation][sofware] = vep_df
    normalisation_to_software_to_vepDf = {}

    # go through each of the vcfs
    for normalised_vcf in all_normalised_vcfs:

        # get names
        software = normalised_vcf.split("/")[-2].split("_")[0]
        typeNormalisation = normalised_vcf.split("/")[-1].split(".")[2]
        #if software not in {"freebayes"}: continue # DEBUUUG

        # check if any of the integrated datasets are already created
        fileprefix = "%s/integrated_variants_%s_ploidy%i"%(opt.outdir, typeNormalisation, opt.ploidy)
        if any([fun.file_is_empty("%s.%s"%(fileprefix, x)) for x in ["py", "tab"]]) or opt.replace is True or opt.replace_vep_integration is True:

            print("Annotating variants with vep for %s"%normalised_vcf)

            # filter
            print("Loading vcf to calculate the allele frequencies")
            print(normalised_vcf)
            vcf_df , variant_to_frequency, variant_to_filter, var_to_GT, var_to_filters = fun.load_vcf_intoDF_GettingFreq_AndFilter(normalised_vcf)

            vcf_df.to_csv(normalised_vcf, sep="\t", index=False)

            # define an output file for VEP
            annotated_vcf = "%s_annotated.tab"%normalised_vcf; annotated_vcf_tmp = "%s.tmp"%annotated_vcf

            # run annotation by VEP
            if fun.file_is_empty(annotated_vcf) or opt.replace is True or opt.replace_vep_integration is True:
            #if True: # DEBUG
                print("Annotating with VEP %s"%normalised_vcf)
                fun.remove_file(annotated_vcf)
                fun.remove_file(annotated_vcf_tmp)

                vep_cmd = "%s --input_vcf %s --outfile %s --ref %s --gff %s --mitochondrial_chromosome %s --mito_code %i --gDNA_code %i "%(run_vep, normalised_vcf, annotated_vcf_tmp, opt.ref, gff_with_biotype, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code)

                fun.run_cmd(vep_cmd)

                os.rename(annotated_vcf_tmp, annotated_vcf)


            # add the allele frequency of each variant, as calculated in alternative_allelle_frequencies
            print("getting vep table")
            vep_df = fun.load_vep_table_intoDF(annotated_vcf)
            vep_df["fraction_reads_coveringThisVariant"] = [variant_to_frequency[idx] for idx in vep_df.index]
            vep_df["FILTERtag"] = [variant_to_filter[idx] for idx in vep_df.index]
            vep_df["GT"] = [var_to_GT[idx] for idx in vep_df.index]
            vep_df["additional_filters"] = [var_to_filters[idx] for idx in vep_df.index]

            annotated_vcf_corrected = "%s.corrected"%annotated_vcf
            vep_df.to_csv(annotated_vcf_corrected, sep="\t", index=False)

            # add to the dictionary
            normalisation_to_software_to_vepDf.setdefault(typeNormalisation, {}).setdefault(software, vep_df)


    # generate integrated table with each software and also the FILTERtag
    for norm, software_to_vepDF in normalisation_to_software_to_vepDf.items():
        print("Working on, ", norm)

        # define the fileprefix and continur if it has not been generated
        fileprefix = "%s/integrated_variants_%s_ploidy%i"%(opt.outdir, norm, opt.ploidy)

        if any([fun.file_is_empty("%s.%s"%(fileprefix, x)) for x in ["py", "tab"]]) or opt.replace is True or opt.replace_vep_integration is True:
            print("stacking dfs...")

            # stack all the dfs with the common parts, indicating the software
            df = pd.DataFrame()
            for software, vepDF in software_to_vepDF.items(): df = df.append(vepDF[['#Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'Extra', 'chromosome', 'position', 'ref', 'alt']])

            # remove the duplicated entries
            df["is_duplicate"] = df.duplicated(subset=['chromosome', 'position', 'ref', 'alt', 'Gene'], keep="first") # the first occurrence is False
            df = df[~(df.is_duplicate)] 
            df["var"] = df.index

            # check that the duplication also applies to the variant + gene combintaion
            if sum(df["is_duplicate"])!=sum(df.duplicated(subset=['#Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'Extra', 'chromosome', 'position', 'ref', 'alt'], keep="first")): raise ValueError("The variants are not equaly annotated (with VEP) in all the callers")

            # add which programs called the variants and with the tags
            software_to_setVars = {software : set(vepDF.index) for software, vepDF in software_to_vepDF.items()}
            for software, setVars in software_to_setVars.items(): 
                print("Working on %s"%software)

                # add the sole calling
                df["%s_called"%software] = df["var"].apply(lambda x: x in setVars)

                # map each of the vars to a tag
                var_to_tag = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["FILTERtag"])), 
                              **{var : "" for var in df[~df["%s_called"%software]].index}}

                # map the original uploaded variation and the GT filter
                var_to_GT_index = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["GT_index"])), 
                    **{var : "" for var in df[~df["%s_called"%software]].index}}

                var_to_OriginalUploadedVar = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["#Uploaded_variation_original"])), 
                    **{var : "" for var in df[~df["%s_called"%software]].index}}

                # map each of the vars to the frequency of calling
                var_to_freq = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["fraction_reads_coveringThisVariant"])), 
                               **{var : 0.0 for var in df[~df["%s_called"%software]].index}}

                # map each of the vars to the GT
                var_to_GT = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["GT"])), 
                               **{var : "" for var in df[~df["%s_called"%software]].index}}

                # map each of the vars to the additional_filters
                var_to_filters = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["additional_filters"])), 
                               **{var : "" for var in df[~df["%s_called"%software]].index}}


                # add to df and define if True
                df["%s_FILTERtag"%software] = df["var"].apply(lambda x: var_to_tag[x])
                df["%s_fractionReadsCoveringThisVariant"%software] = df["var"].apply(lambda x: var_to_freq[x])
                df["%s_PASS"%software] = df["%s_FILTERtag"%software].apply(lambda x: x=="PASS")
                df["%s_GT"%software] = df["var"].apply(lambda x: var_to_GT[x])
                df["%s_additional_filters"%software] = df["var"].apply(lambda x: var_to_filters[x])
                df["%s_GT_index"%software] = df["var"].apply(lambda x: var_to_GT_index[x])
                df["%s_#Uploaded_variation_original"%software] = df["var"].apply(lambda x: var_to_OriginalUploadedVar[x])

                # check that all the variants that have been called by this software also have a fractionReadsCoveringThisVariant
                if sum(df["%s_fractionReadsCoveringThisVariant"%software]>0.0) <= sum(df["%s_called"%software])*0.95: 
                    raise ValueError("In %s, %s not almost all the variants that have been called also have a fraction of reads covering them. This may reflect an error in the way how read frequencies are calculated"%(name_sample, software))

            # write and save object
            df.to_csv("%s.tab.tmp"%fileprefix, sep="\t", index=False)
            fun.save_object(df, "%s.py.tmp"%fileprefix)
            os.rename("%s.tab.tmp"%fileprefix, "%s.tab"%fileprefix)
            os.rename("%s.py.tmp"%fileprefix, "%s.py"%fileprefix)

            print("files saved into %s"%("%s.tab"%fileprefix))



    # get the variants that are present in more than x programs
    df = fun.load_object("%s/integrated_variants_norm_vcflib_ploidy%i.py"%(opt.outdir, opt.ploidy))

    ##############################################
    ######## GENERATE THE ANNOTATION FILE ########
    ##############################################

    variantAnnotation_table = "%s/variant_annotation_ploidy%i.tab"%(opt.outdir, opt.ploidy)
    if fun.file_is_empty(variantAnnotation_table) or replace is True:

        print("generating variant annotation table")

        # add fields 
        df['is_snp'] = (df["ref"].apply(len)==1) & (df["ref"]!="-") & (df["alt"].apply(len)==1) & (df["alt"]!="-")

        prot_altering_mutations = {'missense_variant', 'start_lost', 'inframe_deletion', 'protein_altering_variant', 'stop_gained', 'inframe_insertion', 'frameshift_variant', 'stop_lost', 'splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant', 'non_coding_transcript_exon_variant'}
        df["consequences_set"] = df.Consequence.apply(lambda x: set(str(x).split(",")))
        df["is_protein_altering"] = df.consequences_set.apply(lambda x: len(x.intersection(prot_altering_mutations))>0)

        # generate a table that has all the variant annotation info
        varSpec_fields = ['#Uploaded_variation', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'is_snp', 'is_protein_altering']

        # write the df were there are some PASS vars
        programs = {"HaplotypeCaller", "freebayes", "bcftools"}
        df["number_PASS_programs"] = df.apply(lambda r: sum([r["%s_PASS"%p] for p in programs if "%s_PASS"%p in df.keys()]), axis=1)
        variantAnnotation_table_PASS = "%s/variant_annotation_ploidy%i_anyPASS.tab"%(opt.outdir, opt.ploidy)
        df_PASS = df[df.number_PASS_programs>0]
        df_PASS[varSpec_fields].drop_duplicates().to_csv(variantAnnotation_table_PASS, sep="\t", header=True, index=False)

        # write the final vars
        variantAnnotation_table_tmp = "%s.tmp"%variantAnnotation_table
        df[varSpec_fields].drop_duplicates().to_csv(variantAnnotation_table_tmp, sep="\t", header=True, index=False)
        os.rename(variantAnnotation_table_tmp, variantAnnotation_table)

    ##############################################
    ##############################################
    ##############################################

    ##############################################
    ######## GENERATE THE VARIANTS VEP FILE ######
    ##############################################

    # generate a file that has all the information about the samples



    ##############################################
    ##############################################
    ##############################################


    #################################################
    ######## GENERATE THE INTERSECTING FILES ########
    #################################################

    """
    if opt.caller=="all": 
        print("getting intersecting VCFs")
     
        # add the number of programs with PASS
        programs = {"HaplotypeCaller", "freebayes", "bcftools"}
        df["number_PASS_programs"] = df.apply(lambda r: sum([r["%s_PASS"%p] for p in programs]), axis=1)


        for minPrograms in [1, 2, 3]:

            # filter
            df_PASS = df[df.number_PASS_programs>=minPrograms]

            # write a vcf with this df
            intersecting_vcf = "%s/integrated_variants_PASSatLeast%i_ploidy%i.vcf"%(opt.outdir, minPrograms, opt.ploidy)
            #if fun.file_is_empty(intersecting_vcf) or opt.replace is True:
            if True: # debug

                print("getting vcf of samples that are called by at least %i programs"%minPrograms)
                fun.write_integrated_smallVariantsTable_as_vcf(df_PASS, intersecting_vcf, opt.ploidy)
    """

    #################################################
    #################################################
    #################################################

print("VarCall Finished")

# create outfile
open("%s/finsihedVarCall_CNV_file_ploidy%i.txt"%(opt.outdir, opt.ploidy), "w").write("finsihed with pipeline\n")

# at the end remove all the non-essential files
if opt.remove_smallVarsCNV_nonEssentialFiles is True: fun.remove_smallVarsCNV_nonEssentialFiles(opt.outdir, opt.ploidy)

