#!/bin/bash

# This tests the cnv pipeline
parentDir=/home/mschikora/samba 
threads=4
if [ ! -d $parentDir ]; then 
	parentDir=/gpfs/projects/bsc40/mschikora
	threads=48
fi


curDir="$parentDir"/scripts/perSVade/perSVade_repository/testing
pipeline="$parentDir"/scripts/perSVade/perSVade_repository/scripts/perSVade.py

##### Cgalbrata #####


# RUN1_CST34_2G_FLZ 1M reads
reads1="$curDir"/Cglabrata/RUN1_CST34_2G_FLZ_subsampled/sampled_readsR1.fq
reads2="$curDir"/Cglabrata/RUN1_CST34_2G_FLZ_subsampled/sampled_readsR2.fq
outdir="$curDir"/Cglabrata/RUN1_CST34_2G_FLZ_subsampled/outdir_perSVade

# RUN1_CST34_2G_FLZ all reads
#reads1="$curDir"/Cglabrata/RUN1_CST34_2G_FLZ/readsR1.fq.gz
#reads2="$curDir"/Cglabrata/RUN1_CST34_2G_FLZ/readsR2.fq.gz
#outdir="$curDir"/Cglabrata/RUN1_CST34_2G_FLZ/outdir_perSVade

# outdir the auto generation of bam
#outdir="$curDir"/Cglabrata/RUN1_CST34_2G_FLZ/outdir_perSVade_autoBamGeneration

# this is for the whole genome
refgenome="$parentDir"/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_chromosomes.fasta
gff="$parentDir"/Cglabrata_antifungals/data/Cglabrata_genomes_and_annotations/C_glabrata_CBS138_current_features.gff
previous_repeats_table="$parentDir"/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/5478_Candida_glabrata/reference_genome_dir/reference_genome.fasta.repeats.tab
# generate the table with the known genomes and reads

# define the table with real sv from the nanopore assemblies
python "$parentDir"/scripts/perSVade/perSVade_repository/testing/Cglabrata/generate_SVtable_CglabrataStrains.py
close_shortReads_table="$parentDir"/scripts/perSVade/perSVade_repository/testing/Cglabrata/Cglabrata_strains_nanopore_assemblies/table_assemblies_and_shortReads.tab

# define the taxID
taxID=5478

#####################

# test the running of the SV calling
python $pipeline -r $refgenome -o $outdir -p 1 -f1 $reads1 -f2 $reads2 -mchr mito_C_glabrata_CBS138 --threads $threads -gff $gff --verbose --previous_repeats_table $previous_repeats_table --fast_SVcalling --skip_cleaning_outdir

# test VarCall pipeline
#python $pipeline -r $refgenome -o $outdir -p 1 -f1 $reads1 -f2 $reads2 -mchr mito_C_glabrata_CBS138 --threads $threads --skip_SVcalling --run_smallVarsCNV -gff $gff


# test with 'auto' generation of bam file
#python $pipeline -r $refgenome -o $outdir -p 1 -f1 auto -f2 auto -mchr mito_C_glabrata_CBS138 --threads $threads --close_shortReads_table $close_shortReads_table

# This pipeline can be run with
#sbatch --qos=debug --time=02:00:00 --job-name=testing_perSVade --cpus-per-task=48 --error=/gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/testing_perSVade_benchmark_stderr.txt --output=/gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/testing_perSVade_benchmark_stdout.txt --workdir=.  --get-user-env test_pipeline.sh
