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
 
# avoid marking duplicates
parser.add_argument("--smallVarsCNV_markDuplicates_inBam", dest="smallVarsCNV_markDuplicates_inBam", action="store_true", default=False, help="Mark duplicates on the input bam")


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
if opt.replace is True: fun.delete_folder(opt.outdir)
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
if opt.smallVarsCNV_markDuplicates_inBam is True: sorted_bam = fun.get_sortedBam_with_duplicatesMarked(opt.sortedbam, threads=opt.threads, replace=opt.replace)

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

###################################
##### GET THE INTEGRATED VARS ##### 
################################### 

# get the merged vcf records (these are multiallelic)
print("getting merged vcf without multialleles")
merged_vcf_all = fun.merge_several_vcfsSameSample_into_oneMultiSample_vcf(filtered_vcf_results, opt.ref, opt.outdir, opt.ploidy, replace=opt.replace, threads=opt.threads)

# get the variants in a tabular format
variantInfo_table = "%s/variant_calling_ploidy%i.tab"%(opt.outdir, opt.ploidy)
df_variants = fun.write_variantInfo_table(merged_vcf_all, variantInfo_table, replace=opt.replace)

# define the used programs
if opt.caller=="all": all_programs = sorted(["HaplotypeCaller", "bcftools", "freebayes"])
else: all_programs = sorted(opt.caller.split(","))

# generate a report of the variant calling
variantCallingStats_tablePrefix = "%s/variant_calling_stats_ploidy%i"%(opt.outdir, opt.ploidy)
fun.report_variant_calling_statistics(df_variants, variantCallingStats_tablePrefix, all_programs)

##### KEEP VCFS THAT PASS some programs #########

for minPASS_algs in [1, 2, 3]:

    simplified_vcf_PASSalgs = "%s/variants_atLeast%iPASS_ploidy%i.vcf"%(opt.outdir, minPASS_algs, opt.ploidy)
    simplified_vcf_PASSalgs_tmp = "%s.tmp"%simplified_vcf_PASSalgs

    if fun.file_is_empty(simplified_vcf_PASSalgs) or opt.replace is True:

        print("getting vcf with vars called by >=%i programs"%minPASS_algs)

        # define the interesting variants
        df_PASS = df_variants[df_variants["NPASS"]>=minPASS_algs]

        # add the FORMAT
        vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

        # rename the ID
        df_PASS = df_PASS.rename(columns={"#Uploaded_variation":"ID"})

        # set the FILTER ad the number of pass programs
        df_PASS["FILTER"] = df_PASS.NPASS.apply(str) + "xPASS"

        # set an empty INFO
        df_PASS["INFO"] = "."

        # set the format
        df_PASS["FORMAT"] = "GT:AF:DP"

        # add the sample according to FORMAT
        df_PASS["SAMPLE"] = df_PASS.common_GT + ":" + df_PASS.mean_fractionReadsCov_PASS_algs.apply(lambda x: "%.4f"%x) + ":" + df_PASS.mean_DP.apply(lambda x: "%.4f"%x)

        # initialize header lines
        valid_header_starts = ["fileformat", "contig", "reference", "phasing"]
        header_lines = [l for l in fun.get_df_and_header_from_vcf(merged_vcf_all)[1] if any([l.startswith("##%s"%x) for x in valid_header_starts])]

        # add headers
        header_lines += ['##FILTER=<ID=1xPASS,Description="The variant PASSed the filters for 1 algorithm">',
                         '##FILTER=<ID=2xPASS,Description="The variant PASSed the filters for 2 algorithms">',
                         '##FILTER=<ID=3xPASS,Description="The variant PASSed the filters for 3 algorithms">',

                         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype. If there are discrepacncies in the GT between the algorithms where this var PASSed the filters GT is set to .">',
                         '##FORMAT=<ID=DP,Number=1,Type=Float,Description="Mean read depth of the locus from algorithms where this variant PASSed the filters">',
                         '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Mean fraction of reads covering the ALT allele from algorithms where this variant PASSed the filters">',

                         '##source=%s'%("_".join(all_programs))]

        # write
        open(simplified_vcf_PASSalgs_tmp, "w").write("\n".join(header_lines) + "\n" + df_PASS[vcf_fields].to_csv(sep="\t", index=False, header=True))
        os.rename(simplified_vcf_PASSalgs_tmp, simplified_vcf_PASSalgs)

#################################################


# stop if there is no GFF provided
if opt.gff is None: 
    print("WARNING: No gff provided. Skipping the annotation of the variants")
    sys.exit(0)

######### RUN VEP AND GENERATE ANNOTATION TABLE #########

variantAnnotation_table = "%s/variant_annotation_ploidy%i.tab"%(opt.outdir, opt.ploidy)
if fun.file_is_empty(variantAnnotation_table) or opt.replace is True:

    # define an output file for VEP
    annotated_vcf = "%s_annotated.tab"%merged_vcf_all; annotated_vcf_tmp = "%s.tmp"%annotated_vcf

    # run annotation by VEP
    if fun.file_is_empty(annotated_vcf) or opt.replace is True or opt.replace_vep_integration is True:

        print("Annotating with VEP %s"%merged_vcf_all)
        fun.remove_file(annotated_vcf)
        fun.remove_file(annotated_vcf_tmp)
        for f in os.listdir(fun.get_dir(annotated_vcf)): 
            if ".tmp.raw." in f: fun.remove_file("%s/%s"%(fun.get_dir(annotated_vcf), f))

        vep_cmd = "%s --input_vcf %s --outfile %s --ref %s --gff %s --mitochondrial_chromosome %s --mito_code %i --gDNA_code %i "%(run_vep, merged_vcf_all, annotated_vcf_tmp, opt.ref, gff_with_biotype, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code)

        fun.run_cmd(vep_cmd)

        os.rename(annotated_vcf_tmp, annotated_vcf)

    # get into df
    df_vep = pd.read_csv(annotated_vcf, sep="\t")

    # check that the relationship between the VEP Uploaded_var and merged_vcf_all is 1:1
    uploaded_variation = set(df_vep["#Uploaded_variation"])
    all_variants = set(fun.get_df_and_header_from_vcf(merged_vcf_all)[0]["ID"])

    if len(uploaded_variation.difference(all_variants))>0: raise ValueError("There are some uploaded variations that can't be found in all_variants")

    # deinfe the unnanotated vars as those that are not in the VEP output and are also not missing 
    missing_vars = all_variants.difference(uploaded_variation)
    unnanotated_vars = {v for v in missing_vars if v.split("/")[-1]!="*"}

    if len(unnanotated_vars)>0: 
        print("WARNING: There are some variants that have not been annotated with VEP:\n%s\n (%i/%i in total)"%("\n".join(unnanotated_vars), len(unnanotated_vars), len(all_variants)))

    # get variant annotation table
    print("generating variant annotation table")

    # add fields 
    df_vep["ref"] = df_vep["#Uploaded_variation"].apply(lambda x: x.split("/")[-2])
    df_vep["alt"] = df_vep["#Uploaded_variation"].apply(lambda x: x.split("/")[-1])

    df_vep['is_snp'] = (df_vep["ref"].apply(len)==1) & (df_vep["ref"]!="-") & (df_vep["alt"].apply(len)==1) & (df_vep["alt"]!="-")

    prot_altering_mutations = {'missense_variant', 'start_lost', 'inframe_deletion', 'protein_altering_variant', 'stop_gained', 'inframe_insertion', 'frameshift_variant', 'stop_lost', 'splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant', 'non_coding_transcript_exon_variant'}



    df_vep["consequences_set"] = df_vep.Consequence.apply(lambda x: set(str(x).split(",")))
    df_vep["is_protein_altering"] = df_vep.consequences_set.apply(lambda x: len(x.intersection(prot_altering_mutations))>0)

    # generate a table that has all the variant annotation info
    varSpec_fields = ['#Uploaded_variation', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'is_snp', 'is_protein_altering']

    # write the final vars
    variantAnnotation_table_tmp = "%s.tmp"%variantAnnotation_table
    df_vep[varSpec_fields].drop_duplicates().to_csv(variantAnnotation_table_tmp, sep="\t", header=True, index=False)
    os.rename(variantAnnotation_table_tmp, variantAnnotation_table)

#############################################


print("VarCall Finished")

# at the end remove all the non-essential files
if opt.remove_smallVarsCNV_nonEssentialFiles is True: fun.remove_smallVarsCNV_nonEssentialFiles(opt.outdir, opt.ploidy)

