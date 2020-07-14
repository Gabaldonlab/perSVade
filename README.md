# perSVade: personalized Structural Variation detection

perSVade is a method that runs structural variation (SV) calling and interpretation for a set of paired end WGS short reads. It is a pipeline to call breakpoints with  GRIDSS (https://github.com/PapenfussLab/gridss) and summarize them into complex structural variants with CLOVE (https://github.com/PapenfussLab/clove), with some added features. perSVade provides an automated benchmarking and parameter selection for these methods in any genome or sequencing run. This is useful for species without reccommended running and filtering parameters. In addition, it provides an automated report of the SV calling accuracy on simulations and real data, useful to assess the confidence of the results. The pipeline has not been extensively tested in several architectures and species, so that it is not intended for widespread usage. The next release will include these.


## Installation:
perSVade is written in python and R for Linux. Most of the dependecies can be installed through conda. It is advisable that you install these dependencies by creating a conda environment (i.e.: perSVade_env) with all of them, which we provide, with the following command:

`conda env create --file installation/perSVade_env.yml --name perSVade_env`

When running the pipeline make sure that the python interpreter is the one of this environment. To test this execute:

`conda activate perSVade_env`

`which python`

It is expected to print: 

`<path_to_conda>/envs/perSVade_env/bin/python`

This is essential so that all the dependencies of the pipeline are met.

In addition, there are some dependencies that are included in the respository "installation/external_software" (only in the "release" packages). These are gridss (tested on version 2.8.1), clove (tested on version 0.17), gztools (installed from https://github.com/circulosmeos/gztool/releases/download/v0.11.5/gztool-linux.x86_64), vcfvalidator (installed from https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.4/vcf_validator_linux) and bcftools (installed from https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2)

## Comments for the installation of extra dependencies
The non-conda dependencies can be installed like this (if you wanted to reinstall them):

1. change to the directory where you installed the perSVade repository:

`cd <perSVade_installation_dir>/installation/external_software`

2. download dependencies

`wget https://github.com/PapenfussLab/gridss/releases/download/v2.8.1/gridss-2.8.1-gridss-jar-with-dependencies.jar`

`wget https://github.com/PapenfussLab/gridss/releases/download/v2.8.1/gridss.sh`

`wget https://github.com/PapenfussLab/clove/releases/download/v0.17/clove-0.17-jar-with-dependencies.jar`

`wget https://github.com/circulosmeos/gztool/releases/download/v0.11.5/gztool-linux.x86_64`

`wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.4/vcf_validator_linux`

`wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2`

3. setup bcftools

`tar -xvf bcftools-1.10.2.tar.bz2`

`rm bcftools-1.10.2.tar.bz2`

`cd bcftools-1.10.2`

`./configure --prefix=$PWD`

`make`

`make install`


4. give execution permssion to all the files:

`chmod u+x *`

You may want to repeat this in case you have problems running any of the programs with the pipeline

The conda environment can be exported to a .yml file with:

`conda env export --no-builds -n perSVade_env --file perSVade_env.yml`

## Running in MareNostrum

If you are working from any cluster that has access to the BSC /gpfs filesystem you can activate the perSVade environment from its location in mschikora. The only risk is that Miki may be editing the environment for the development of perSVade. You should be running this from an interactive node in MN like this:

`salloc -n 1 -c 48 -t 02:00:00 --qos debug` # runs an interactive session. ESSENTIAL

`source /gpfs/projects/bsc40/mschikora/anaconda3/etc/profile.d/conda.sh`  # activate the conda environment of mschikora

`conda activate perSVade_env` # activate the environment

VERY IMPORTANT NOTE:

NEVER activate the enviroment in the login of MN, as this is a large environment and uses a lot of resources. In fact, sometimes the login FAILS for all MN users when you activate this environent, so please don't do it.

## Running
Once you have installed all the dependencies, you can call the perSVade pipeline with:

`conda activate perSVade_env`

`python ./scripts/perSVade.py -r <path to the reference genome (fasta)> -o <output_directory> -p <ploidy, 1 or 2> -f1 <forward_reads.fastq.gz> -f2 <reverse_reads.fastq.gz>`


## Example: running SmallVariant calling and CNV detection for paired-end WGS

perSVade also includes the possibility of running small variant calling. You can do this by skipping SV detection, with a command like:

`./scripts/perSVade.py --ref reference_genome.fasta --threads 4 -o ./output_directory -f1 reads_FWD.fastq.gz -f2 reads_FWD.fastq.gz --mitochondrial_chromosome chr_mito --mitochondrial_code 3 --gDNA_code 12 -gff features.gff  --run_smallVarsCNV --skip_SVcalling --caller all --coverage 20 --ploidy 2 --remove_smallVarsCNV_nonEssentialFiles`

This will align the reads with `bwa mem` and run `GATK HaplotypeCaller`, `freebayes` and `bcftools call` on the provided reads. The variants are filtered with the default parameters and the specified coverage. The resulting variants are be merged and annotated with `Ensembl VEP`. In addition, the read depth of each gene will be calculated with `mosdepth`.

Type `./scripts/perSVade.py -h` to understand wahat is the meaning of these options. Some important remarks:

1. `--mitochondrial_code` and `--gDNA_code` are the NCBI translation code of your species, check them in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . Note that they may be wrong sometimes, so that it is good to double check with the literature.
2. `--mitochondrial_chromosome` is supposed to include the name of the chromosome in the reference genome fasta.
3. `--run_smallVarsCNV` is the cmd that indicates perSVade to also run smallVariant calling.
4. `--skip_SVcalling` tells the pipeline to skip the calling of structural variation
5. `--coverage` is the minimum coverage for a variant to be kept. This parameter should be related to the coverage of your library. You may be careful with setting it too low (i.e. <15) as it will yield many false positive calls. It is reasonable to check how other similar studies set this parameter.
6. `--remove_smallVarsCNV_nonEssentialFiles` is useful to clean the output directory from non-essential output files.
7. If you may want to understand what the programs of this pipeline do. Check freebayes, bcftools, GATK HaplotypeCaller and ENSEMBL Variant Annotation Predictor.

This will output the following files and folders under `./output_directory`:

1. `aligned_reads.sorted.bam` and `aligned_reads.sorted.bam.bai`: the aligned reads sorted by genomic coordinate.
2. `reference_genome_dir` is a directory that contains files related to the provided genome and annotations. This may be removed if necessary.
3. `aligned_reads.bam.sorted.calculating_windowcoverage/coverage_windows_<n_nucleotides>bp.tab` contains a table with the coverage for windows of the genome that are as long as 5% of the median chromosome length (n_nucleotides). This is the output of mosdepth. This table contatins the following fields:

    1. `#chrom` is the chromosome name.
    2. `start` is the 0-based start coordinates of the region
    3. `end` is the 1-based end coordinate of the region
    4. `length` is the length of the region
    5. `mediancov_1` is the median read depth in the region
    6. `percentcovered_1` is the perecentage of the region that is covered with 1 read or more
    7. `nocoveragebp_1` is the number of bases that have no coverage

4. `smallVars_CNV_output/CNV_results/genes_and_regions_coverage.tab` is a table that contains the coverage of all the genes in the provided gff. These are the fields:

    1. `chromosome`, `start`, `end` and `length` are the coordinates of the gene
    2. `median_reads_per_gene` is the median read depth for this gene
    3. `nocoveragebp_1` is the number of bases that have no coverage
    4. `percentcovered_1` is the perecentage of the region that is covered with 1 read or more
    5. `ID` is the gene ID in the gff (parsed from the attributes). If there are some duplicated IDs with parts it corresponds to the union of ID and part.
    6. `fraction_covered_by_MoreThan1read` is equivalent to `percentcovered_1`, but from 0 to 1.
    7. `relative_coverage` is the `median_reads_per_gene` divided by the median of all `median_reads_per_gene`.
    8. All the fields with a `_+-10kb_region` sufix are equivalent to the ones that do not have it, but for a region that starts at -10kb of the gene start and ends at +10kb of the gene end. This can be useful for CNV analysis.

5. `smallVars_CNV_output` contains some folders and files related to the small variant calling:

    1. `bcftools_ploidy2_out`, `HaplotypeCaller_ploidy2_out` and `freebayes_ploidy2_out` are folders that contain the raw and filtered vcf files of each fo the programs.
    2. `merged_vcfs_allVars_ploidy2.vcf` is a vcf file with the merged output of the algorithms used in the variant calling. You can check the header to understand what are all the tags in the `INFO` field. Here, the multiallelic loci are split to ease the analysis. Note that the `ID` field corresponds to the `#Uploaded_variation` field of `variant_annotation_ploidy2.tab`.
    3. `variant_calling_ploidy2.tab` is a tabular version of `merged_vcfs_allVars_ploidy2.vcf`. Each column contains information parsed from the vcf file. This file is easier to manage because it has all the info in tabular format.
    4. `variant_annotation_ploidy2.tab` is the output of VEP, where each line corresponds to a particular alteration in a given gene.
    5. `variants_atLeast<n_PASS_programs>PASS_ploidy2.vcf` are the variants that PASS de filters in at least n_PASS_programs algorithms. The INFO and FORMAT fields are simplified to take less space. It may be useful to take, for example, variants that PASS de filters by at least 2 algorithms.





## Method
Breakpoints are called using gridss and integrated into complex structural variation with clove. The straightforward implementation of these algorithms is challenging for 1) genomes without established parameters and 2) sequencing runs were there is a "smiley-pattern" in read-depth (such as https://www.cell.com/cell/pdf/S0092-8674(16)31071-6.pdf). The latter occurs when read-depth is correlated with the distance to the telomere, which may be an artifact of library preparation and/or sequencing. This impedes the usage of a single read-depth threshold for filtering deletions and tandem duplications (used by clove). perSVade selects the running and filtering parameters from a simulation-based optimization, including the following pipeline for the input sequencing:

1. Simulation of randomly-placed insertions (27 copy-and-paste and 28 cut-and-paste), translocations (one for each gDNA chromosome, where half of them are unbalanced), 55 inversions, 55 deletions and 55 tandem duplications into the reference genome using the RSVSim package (https://www.bioconductor.org/packages/release/bioc/html/RSVSim.html). This corresponds to 50 gDNA and 5 mtDNA variants except for translocations. The variant length is set to follow a decaying beta-distribution function between 20% of the shortest chromosome and 50 bp, so that smaller variants would be more likely to occur. This 20% is progressively reduced up to 1% if the random set of variants does not fit into the genome for being to long. We perform three such replicate genomes for each sample, hereafter referred as SV genomes.

2. For each SV genome, short paired-end reads generation with wgsim (https://github.com/lh3/wgsim). In order to simulate the biases in read depth found in the sample we generate, for each genomic window of 1kb, a different number of reads pairs calculated from sequence features that may influence coverage. We consider both features influencing actual variation in read number (GC content, hexamer content and the “smiley-pattern”) and those related to read mapping errors (coverage by repetitive elements and mappability). To find the quantitative contribution of these we first predict coverage from the position in the chromosome (with a quadratic fit), which is a model of the contribution of the “simley-pattern” to variation in read depth. The residuals of this fit are predicted with a linear model from GC content, coverage by repetitive elements (inferred with RepeatModeller (https://github.com/Dfam-consortium/RepeatModeler) and RepeatMasker (https://github.com/rmhubley/RepeatMasker)) and mappability (calculated with genmap (https://github.com/cpockrandt/genmap)). We also test whether the content of any hexamer can linearly predict (p<0.05) the residual variation of this model (from cross-validation, training the model for each hexamer on all but one chromosomes and testing on the remaining one), which are kept as coverage predictors. Note that hexamer and GC content (contributing to coverage through changes in real read numbers) may be correlated with mappability or repeat coverage (contributing to coverage because of mapping errors), hampering the inference of the separate contribution of each feature. We simplify this by assuming that any coverage variation explainable from GC or hexamer content contributes to changes in the actual number of reads. Accordingly, we build a linear model that infers the coverage from the position in the chromosome (quadratic fit) and the residual variation from GC and hexamer content. The number of reads generated for each genomic window is proportional to the coverage predicted from sequence features of this model. The simulated read length, mean insert size and standard deviation are those of the actual sample (calculated with samtools (https://github.com/samtools/) and picard (https://github.com/broadinstitute/picard/)). In addition, we set an error rate of 0.02 for wgsim. We simulate reads from both haploid (all reads from the SV genome) and diploid heterozygous (a 1:1 mix of reads from the SV and the reference genome) genomes in order to test the pipeline on structural variants in aneuploid chromosomes or heterozygous conformation. This yields 6 simulations (two for each SV genome).

3. Breakpoint calling with gridss from the aligned simulated reads (using bwa mem (https://github.com/lh3/bwa)) for several combinations of filtering parameters. The following fields are considered (used in a recent update of the software (https://github.com/hartwigmedical/gridss-purple-linx)) for filtering: 

    1. Minimum number of fragments (5 to 30).
    2. Minimum allele frequency (0.05 to 0.9).
    3. Maximum to be considered a small deletion or duplication event (100 to 108).
    4. “FILTER” column (all combinations tried).
    5. Presence of a poly-G or poly-C (16bp) sequence around the breakpoint (always filtered).
    6. Split-reads support (present or absent), only applicable for small deletion or deletion events.
    7. Discordant read-pair support  (present or absent), not applicable for small deletion or deletion events.
    8. Maximum strand bias (0.9 to 0.95), only applicable for small deletion or deletion events.
    9. Maximum microhomology length (100 to 108).
    10. Maximum inexact homology length (100 to 108).
    11. Range of lengths of deletion breakpoints with an inexact microhomology of 6bp or greater that should be filtered out (100-800, 50-900, 200-700 and 0-1).
    12. Minimum inversion length (0, 40 and 108)
    13. Maximum difference between the breakpoint-inserted sequence and the length of a deletion (0, 5, 10 and 108)

  This yields 7.74x1e6 combinations of filters to be tried on the gridss output, each of them resulting in an independent set of called breakpoints.

4. For each set of breakpoints, integration into complex events with clove. The output of clove is parsed to output translocations (balanced or not), insertions (copy-and-paste or cut-and-paste) tandem duplications, deletions and inversions to match the “known” simulated variants (in 1) and allow benchmarking of each filter set. In addition, the set of breakpoints that is not classified into any of these types of variants is kept. This parsing was is trivial for the following types of variants:

    1. Unbalanced translocations (the arm the origin chromosome is copied and replaces an arm of the target chromosome), which may result in non-reciprocal interchromosomal breakpoint. However, not all such breakpoints imply unbalanced translocations, which required an additional method to find the true events. We reasoned that only those that imply a deletion of the target chromosome and a duplication of the origin may be real unbalanced translocations, which is tested from read-depth. For each combination of origin-target chromosomes and breakpoints, we calculate “probability of unbalanced translocation” (p_unbalTRA) as the mean of the “probability of origin duplication” (p_origin_dup) and the “probability of target deletion” (p_target_del). p_origin_dup is calculated as the relative coverage of the origin region - 1 (set to 1 as maximum), while p_target_del is calculated as 1 - relative coverage of the target region. We benchmark several thresholds for p_unbalTRA between 0 and 1 to find the optimum ones.

    2. Deletions and tandem duplications, which generate breakpoints referred as DEL and DUP.  Unfortunately, these can also result from complex structural variants, which hampers direct interpretation. Read-depth validation can be used on such variants, and clove calls as true tandem duplications or deletions those events with a coverage above or below (respectively) a symmetric interval across a provided value. However, we consider that this strategy may lead to inaccurate results when a coverage "smiley-pattern" affects the input sample. We thus filter each type of variant based on an independent coverage threshold. Any DEL event with a relative coverage below the threshold (termed max_coverage_for_del) or DUP event with a relative coverage above the threshold (min_coverage_for_tan) is called as a true deletion or tandem duplication, respectively. We benchmark several of these thresholds (0-1 for max_coverage_for_del and 0-3.5 for min_coverage_for_tan) to find the optimum ones.
	
5. For each filter combination (gridss, p_unbalTRA, max_coverage_for_del and min_coverage_for_tan), benchmarking of the calling accuracy as compared to the known simulated variants. This includes up to 3.01x1e8 sets of filters. The chosen optimal parameters for each SV genome are those with the highest F score (a combination of precision and recall). Each event is defined as a true positive if there is a matching “known” variant, with all the breakpoints overlapping by less than 50 bp. For those sets of filters with the same F score we keep the least conservative one. 

6. Choosing the optimal set of filters for structural variant calling on real data. Each of the simulations may yield a different set of optimum filters, where accuracy may be high due to overfitting to the training data. In order to test each of these filter sets on an independent dataset we calculate the accuracy of each of them on the other simulations. We choose the final “optimum set” of filters as the one yielding the maximum minimum F score across all simulations.

7. Running of gridss and clove for the real short-read data, filtering and parsing the output with the chosen optimum set of filters. In addition to the explainable variants (deletions, tandem duplications, insertions and translocations) we keep the unclassified breakpoints for further analysis, which may represent structural variants not included in our simulations. 

We acknowledge that our benchmarking could yield high accuracy because random simulations may not include the sources of false positive calls found in real data. However, if there is no available set of real structural variants in your genome we consider that random simulations are as realistic as possible.





