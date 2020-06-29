#!/usr/bin/env python

# This is a script that runs the testing of perSVade on several species

##### DEFINE ENVIRONMENT #######

# module imports
import os
import sys

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False    
    threads = 4
else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    threads = 48


# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions
import sv_functions as fun

# define paths
perSVade_py = "%s/perSVade.py"%perSVade_dir

# define dirs
outdir_testing = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies"%ParentDir; fun.make_folder(outdir_testing)
outdir_genomes_and_annotations = "%s/scripts/perSVade/perSVade_repository/testing/genomes_and_annotations"%ParentDir

################################


"""
This is how the genomes were obtained:

C. glabrata: reference genome from CGD: the latest version by 12/03/2019, which is v_s02-m07-r35 

C. albicans: 

    ref genome CGD: http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_version_A22-s07-m01-r110_chromosomes.fasta.gz

    gff from CGD: http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_version_A22-s07-m01-r110_features.gff

    From here I keep 'haplotype A' for both files

C. neoformans: ref genome from GenBank GCA_000149245.3

A. fumigatus: 

    gDNA from reference NCBI:

    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.gff.gz
    
    mtDNA from https://www.ncbi.nlm.nih.gov/nuccore/CM016889.1

A. thaliana: ref genome from GenBank GCA_000001735.2

D. melanogaster: ref genome from GenBank GCA_000001215.4

D. rerio: ref genome from GenBank removing the alternate haplotypes. (this is GCA_000002035.4)

H. sapiens: ref genome from GenBank removing the alternate haplotypes. (this is GCA_000001405.28)

For C. glabrata I got the nanopore reads from ~/../mmarcet/nanopore/GABALDON02/assembly_files/BG2/nanopore.reads.pass.fastq.gz and the short reads from Ewa's experiment in RUN4_BG2_SRA_WT

"""

# define the table for C. glabrata
close_shortReads_table_Cglabrata = "%s/scripts/perSVade/perSVade_repository/testing/Cglabrata_table_short_reads.tab"%ParentDir
goldenSet_dir_Cglabrata = "%s/scripts/perSVade/perSVade_repository/testing/Cglabrata_goldenSetReads_BG2"%ParentDir

# define important info about each species: taxID, spName, ploidy

species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138"),
                ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314"),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1"),
                ("746128", "Aspergillus_fumigatus", 1, "CM016889.1"),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1"),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2")]
                #("7955", "Danio_rerio", 2, "NC_002333.2")]
                #("9606", "Homo_sapiens", 2, "NC_012920.1")]

"""
species_Info = [("7227", "Drosophila_melanogaster", 2, "KJ947872.2")]
                #("7955", "Danio_rerio", 2, "NC_002333.2")]
                #("9606", "Homo_sapiens", 2, "NC_012920.1")]
"""

taxIDs_with_noON_overalpping = {"5476", "746128"}

# define the type of run
type_run = "normalRun" # can be 'normalRun' or 'goldenSet'
StopAfterPrefecth_of_reads = False

# go through each species
for taxID, spName, ploidy, mitochondrial_chromosome in species_Info:
    print(taxID, spName)

    # define  the genome and annotations
    genome = "%s/%s.fasta"%(outdir_genomes_and_annotations, spName)
    gff = "%s/%s.gff"%(outdir_genomes_and_annotations, spName)

    # create an outdir
    outdir_perSVade = "%s/%s_%s"%(outdir_testing, taxID, spName); fun.make_folder(outdir_perSVade)

    if type_run=="normalRun":

        # define the table with short reads
        if spName=="Candida_glabrata": close_shortReads_table = close_shortReads_table_Cglabrata
        else: close_shortReads_table = "auto"

        # get the reads from SRA. 3 samples, 3 runs per sample. Process with the 
        cmd = "%s --ref %s --threads %i -o %s --close_shortReads_table %s --target_taxID %s --n_close_samples 3 --nruns_per_sample 3 -f1 skip -f2 skip --mitochondrial_chromosome %s --gff %s --StopAfter_readObtentionFromSRA"%(perSVade_py, genome, threads, outdir_perSVade, close_shortReads_table, taxID, mitochondrial_chromosome, gff)

    elif type_run=="goldenSet":

        # define the goldenSet_dir
        if spName=="Candida_glabrata": goldenSet_dir = goldenSet_dir_Cglabrata
        else: goldenSet_dir = "auto"

        # get the golden set running 
        if taxID in taxIDs_with_noON_overalpping: continue
        cmd = "%s --ref %s --threads %i -o %s --target_taxID %s --n_close_samples 3 --nruns_per_sample 3 -f1 skip -f2 skip --mitochondrial_chromosome %s --gff %s --goldenSet_dir %s --skip_SVcalling"%(perSVade_py, genome, threads, outdir_perSVade, taxID, mitochondrial_chromosome, gff, goldenSet_dir)

    # add options depending on the machine
    if run_in_cluster is True: cmd += " --job_array_mode greasy --queue_jobs debug --max_ncores_queue 48 --time_read_obtention 02:00:00 "
    else: cmd += " --job_array_mode local"

    if StopAfterPrefecth_of_reads is True: cmd += " --StopAfterPrefecth_of_reads"

    fun.run_cmd(cmd)

    #if taxID=="5478": adkjhdakg











