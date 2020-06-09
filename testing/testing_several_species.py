#!/usr/bin/env python

# This is a script that runs the testing of perSVade on several species

##### DEFINE ENVIRONMENT #######

# module imports
import os
import sys

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_slurm_cmd = ""    
    threads = 4
else:
    run_in_slurm_cmd = " --run_in_slurm "
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

# define important info about each species: taxID, spName, ploidy
species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138", "C_glabrata_CBS138_current_chromosomes.fasta", "C_glabrata_CBS138_current_features.gff"),

                ("5476", "Candia_albicans", 2, "no_mitochondria", "Calbicans_chromosomes.fasta", "Calbicans_features.gff"),

                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1", "GCA_000149245.3", "GCA_000149245.3"),

                ("746128", "Aspergillus_fumigatus", 1, "no_mitochondria", "GCA_000002655.1", "GCA_000002655.1"),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", "GCA_000001735.2", " GCA_000001735.2"),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2", "GCA_000001215.4", "GCA_000001215.4"),
                ("7955", "Danio_rerio", 2, "NC_002333.2", "Drerio_chromosomes.fasta", "Drerio_features.gff"),
                ("9606", "Homo_sapiens", 2, "J01415.2", "GCA_000001405.28", "GCA_000001405.28")]

"""
These are the links:
C. albicans:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.gff.gz

D. rerio:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.gff.gz

"""


# go through each species
for taxID, spName, ploidy, mitochondrial_chromosome, genome, gff in species_Info:
    print(taxID, spName)

    if genome!="auto" and not genome.startswith("GCA_"): 
        genome = "%s/%s"%(outdir_genomes_and_annotations, genome)
        gff = "%s/%s"%(outdir_genomes_and_annotations, gff)

    # create an outdir
    outdir_perSVade = "%s/%s_%s"%(outdir_testing, taxID, spName); fun.make_folder(outdir_perSVade)

    # get the reads from SRA. 5 samples, 3 runs per sample
    fun.run_cmd("%s --ref %s --threads %i -o %s --close_shortReads_table auto --target_taxID %s --n_close_samples 5 --nruns_per_sample 3 --StopAfter_readObtentionFromSRA -f1 auto -f2 auto %s --mitochondrial_chromosome %s --gff %s --StopAfter_bamFileObtention"%(perSVade_py, genome, threads, outdir_perSVade, taxID, run_in_slurm_cmd, mitochondrial_chromosome, gff))









