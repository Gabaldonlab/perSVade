# perSVade: personalized Structural Variation detection

perSVade is a method that runs structural variation (SV) calling and interpretation for a set of paired end WGS short reads. It is a pipeline to call breakpoints with  GRIDSS (https://github.com/PapenfussLab/gridss) and summarize them into complex structural variants with CLOVE (https://github.com/PapenfussLab/clove), with some added features. perSVade provides an automated benchmarking and parameter selection for these methods in any genome or sequencing run. This is useful for species without reccommended running and filtering parameters. In addition, it provides an automated report of the SV calling accuracy on simulations and real data, useful to assess the confidence of the results. 

## Pipeline overview

![alt text](https://github.com/Gabaldonlab/perSVade/blob/master/misc/perSVade_pipeline_cartoon.png)

This is a scheme of the key functions of perSVade. The pipeline runs structural variation (SV), small variant (SNPs and IN/DELs) and read depth-based Copy Number Variation (CNV) calling and annotation with a single bash command. By default, perSVade takes a set of paired-end short reads and a reference genome as inputs. It runs `bwa mem` to align these reads, generating a sorted .bam file that is the core of several downstream analyses. These include:

1.  SV calling. This is the core, most novel function  of perSVade. It uses `gridss` to infer a list of breakpoints (two regions of the genome (two breakends) that are joined in the sample of interest and not in the reference genome) from discordant read pairs, split reads and de novo assembly signatures. The breakpoints are summarized into SVs with `clove`. In order to find the best parameters to run these algorithms, perSVade generates two simulated genomes with 50 SVs of each type and tries several combinations (>13,000,000,000) of filters for the `gridss` and `clove` outputs. It selects the filters that have the highest Fvalue (harmonic mean between precision and recall) for each simulated genome and SV type. In order to reduce overfitting, perSVade selects a final set of "best parameters" that work well for all simulations and SV types. This set of best parameters is used for the final calling on the real data. The accuracy (Fvalue, precision, recall) of these parameters on each simulation and SV type is reported as a heatmap, which is useful to evaluate the expected calling accuracy. By default, the simulations are placed randomly across the genome. However, SVs often appear around repetitive elements or regions of the genome with high similarity (i.e.: transposable elements insertions). This means that random simulations may not be realistic, potentially leading to overestimated calling accuracy and a parameter selection that does not work well for real SVs. perSVade can also generate more realistic simulation occuring around regions with known SVs (i.e. regions with SVs called with perSVade) or homologous regions (inferred from BLAST). See "Output" for more details.

2. Read depth-based Copy Number Variation (CNV) calling. CNVs are one type of SVs where there is an alteration in the genomic content (deletions or duplications). The SV calling feature of perSVade (point 1) identifies some CNVs (insertions, tandem duplications, deletions and complex inverted SVs) but it can miss others (i.e.: whole-chromosome duplications or regions with unknown types of rearrangements yielding CNVs). perSVade also includes a pipeline to call CNVs from read-depth alterations. For example, regions with 0x or 2x read-depth as compared to the mean of the genome can be called duplications, or deletions, respectively. A straight forward implementation of this concept to find CNVs is challenging because many genomic features drive variability in read depth independently of CNV. In order to solve this, perSVade calculates the relative coverage for bins of the genome and corrects the effect of the GC content, mappability and distance to the telomere (using non-parametric regression as in https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-014-0100-8). This corrected coverage is used by `CONY`, `AneuFinder` and/or `HMMcopy` to call CNVs across the genome. perSVade generates consensus CNV calls from the three programs taking always the most conservative copy number for each bin of the genome. For example, if the used programs disagree on the copy number of a region the closest to 1 will be taken as the best estimate.

3. Small variant (SNPs and IN/DELs) and gene CNV calling. perSVade includes an option to run a pipeline that performs small variant calling and calculation of read depth for each gene. It runs any of `freebayes`,  `GATK HaplotypeCaller` and/or `bcftools call` for small variant calling and integrates the results into .tab and .vcf files. It runs `mosdepth` for each gene (it requires a .gff file from the user).

4. Integration of read-depth based CNVs and SVs into a single .vcf file. This is a file that is focused on showing the alteration of SVs on specific genomic regions (see the section "Output" for more details). It also removes redundant calls between the CNVs identified with `gridss`+`clove` and those derived from

## Installation

### 1. Downloading the perSVade source code

Download the perSVade source code from one of the releases and decompress. For example:

`wget https://github.com/Gabaldonlab/perSVade/archive/v0.7.tar.gz`

`tar -xvf v0.7.tar.gz; rm v0.7.tar.gz`

This already contains all the scripts to run the pipeline. Note that the created file (for example `perSVade-v0.7`) will be referred as `<perSVade_dir>`

### 2. Create a conda environment with most dependencies

perSVade is written in python, R and bash for Linux. Most of the dependecies can be installed through conda. It is advisable that you install these dependencies by creating a conda environment (for example called perSVade_env) with all of them, which we provide, with the following commands:

`cd <perSVade_dir>`

`conda env create --file installation/perSVade_env.yml --name <env_name>`

`conda activate <env_name>`

### 3. Automatic installation of additional dependencies 

In addition, you should install some additional dependencies and setup of other environments with the following command:

`./installation/setup_environment.sh`

Make sure that this script ends with the message: `SUCCESS: all perSVade dependencies were properly installed.`

NOTE: This will create the following additional environments:

1. `<env_name>_bcftools_1.10.2_env`
2. `<env_name>_ete3_3.0.0_env`
3. `<env_name>_R_env`
4. `<env_name>_gridss_env`
5. `<env_name>_picard_env`
6. `<env_name>_AneuFinder_env`
7. `<env_name>_CONY_env`
8. `<env_name>_HMMcopy_env`
9. `<env_name>_RepeatMasker_env`
10. `<env_name>_RepeatMasker_env`

Make sure that none of them exist before running this script. You can change `<env_name>` to fullfil this requirement.

We note that this was tested with `conda 4.8.0` on a Linux-x86 64-bit architecture, installed at 03/2019. If you have a different conda version, you may change a bit the perSVade_env.yml file so that it does not create dependency problems.

### 4. Test installation

We highly recommend to test that all dependencies were properly installed with the following commands:

`./installation/test_installation/test_installation.py`

This process should take arround 45 minutes on 4 cores. Verify that it finishes with the following message:

`SUCCESS: perSVade was properly installed`

There is a WARNING message that you should look for after running this script:

`WARNING: The connection to SRA did not work`. perSVade includes the option to query and download from the SRA database for the benchmarking of SV calling. This requires a proper network access and SRA connection, which may not always be available. This warning indicates that this process is not possible on your machine. You can skip this connection by providing the reads on your own through `--close_shortReads_table`.


## Running in MareNostrum

If you are working from any cluster that has access to the BSC /gpfs filesystem you can activate the perSVade environment from its location in mschikora. You don't need to re-install anything if you are working in the BSC.

`source /gpfs/projects/bsc40/mschikora/anaconda3/etc/profile.d/conda.sh`  # activate the conda environment of mschikora

`conda activate perSVade_v0.7_env` # activate the environment of perSVade version 0.7. You can change the version

You can next run perSVade from the releases folder (these are stable versions of the pipeline). For example:

`python /gpfs/projects/bsc40/mschikora/scripts/perSVade/releases/perSVade-0.7/scripts/perSVade.py -r <path to the reference genome (fasta)> -o <output_directory> -p <ploidy, 1 or 2> -f1 <forward_reads.fastq.gz> -f2 <reverse_reads.fastq.gz> `

IMPORTANT NOTE: The activation of the perSVade conda environment works well from version 0.7 on. This means that you can activate from the login of MN or interactive nodes. However, the activation of older versions (v0.4 and below) is costly, and it overloads the login nodes. If you want to use an old version of perSVade, always activate it on an interactive node (i.e.: `salloc`). In addition, you can't run the perSVade pipeline from a login, because it takes too many resources. You can submit perSVade as a job or run from an interactive session with `salloc -n 1 --time=02:00:00 -c 48 --qos debug`.


## Running in other systems

Once you have installed all the dependencies, you can call the perSVade pipeline with:

`conda activate perSVade_env`

`python ./scripts/perSVade.py -r <path to the reference genome (fasta)> -o <output_directory> -p <ploidy, 1 or 2> -f1 <forward_reads.fastq.gz> -f2 <reverse_reads.fastq.gz>`



## Example 1: running SmallVariant calling and CNV detection for paired-end WGS

perSVade also includes the possibility of running small variant calling. You can do this by skipping SV detection, with a command like:

`./scripts/perSVade.py --ref reference_genome.fasta --threads 4 -o ./output_directory -f1 reads_FWD.fastq.gz -f2 reads_FWD.fastq.gz --mitochondrial_chromosome chr_mito --mitochondrial_code 3 --gDNA_code 12 -gff features.gff  --run_smallVarsCNV --skip_SVcalling --caller all --coverage 20 --ploidy 2 --remove_smallVarsCNV_nonEssentialFiles --skip_repeat_analysis`

This will align the reads with `bwa mem` and run `GATK HaplotypeCaller`, `freebayes` and `bcftools call` on the provided reads. The variants are filtered with the default parameters and the specified coverage. The resulting variants are be merged and annotated with `Ensembl VEP`. In addition, the read depth of each gene will be calculated with `mosdepth`.

Type `./scripts/perSVade.py -h` to understand wahat is the meaning of these options. Some important remarks:

1. `--mitochondrial_code` and `--gDNA_code` are the NCBI translation code of your species, check them in https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . Note that they may be wrong sometimes, so that it is good to double check with the literature.
2. `--mitochondrial_chromosome` is supposed to include the name of the chromosome in the reference genome fasta.
3. `--run_smallVarsCNV` is the cmd that indicates perSVade to also run smallVariant calling.
4. `--skip_SVcalling` tells the pipeline to skip the calling of structural variation
5. `--coverage` is the minimum coverage for a variant to be kept. This parameter should be related to the coverage of your library. You may be careful with setting it too low (i.e. <15) as it will yield many false positive calls. It is reasonable to check how other similar studies set this parameter.
6. `--remove_smallVarsCNV_nonEssentialFiles` is useful to clean the output directory from non-essential output files.
7. If you may want to understand what the programs of this pipeline do. Check freebayes, bcftools, GATK HaplotypeCaller and ENSEMBL Variant Annotation Predictor.
8. `--skip_repeat_analysis` does not analyze whether your variants overlap repetitive regions. By default it does, and finding repeats in the genome can slow down the analysis.

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


## Output 

The output of `clove` is processed by custom functions to generate 6 different .tab files, each with a different type of SV. These are the types of SVs(note that text in "" indicates the column names of the corresponding .tab files):

![alt text](https://github.com/Gabaldonlab/perSVade/blob/master/misc/simple_SVs.png)
<img src="https://github.com/Gabaldonlab/perSVade/blob/master/misc/simple_SVs.png" width="50%" height="50%">

1. Simple SVs: deletions, inversions and tandemDuplications (duplication of a region which gets inserted next to the affected region). These are described by a chromosome ("Chr"), "Start" and "End" coordinates of the SV. perSVade outputs one .tab file for each of these SV types.

![alt text](https://github.com/Gabaldonlab/perSVade/blob/master/misc/insertions.png)

2. Insertions: a region of the genome (indicated by "ChrA", "StartA", "EndA") is copied ("Copied" is TRUE) or cut ("Copied" is FALSE) and inserted into another region (indicated by "ChrB", "StartB"). "EndB" comes from adding to "StartB" the length of the inserted region. There is one .tab file for insertions.

![alt text](https://github.com/Gabaldonlab/perSVade/blob/master/misc/translocations.png)

3. Translocations: balanced translocations between two chromosomes ("ChrA" and "ChrB"). "StartA" indicates the start of "ChrA" (position 0). "EndA" indicates the position in "ChrA" where the translocation happens. For non-inverted translocations, "StartB" is 0 and "EndB" is the position in "ChrB" where the translocation happens. For inverted translocations, "StartB" is the position in "ChrB" where the translocation happens and "EndB" is the last position of "ChrB". There is one .tab file for translocations.

![alt text](https://github.com/Gabaldonlab/perSVade/blob/master/misc/unclassified_SVs.png)

4. Unclassified SVs: There is one .tab file ("unclassified_SVs.tab") that reports all the variants that are called by `clove` and cannot be assigned to any of the above SV types. These include unclassified breakpoints (which could be part of unresolved/unkown complex variants) and complex inverted SVs (which are non-standard SVs). These types of SVs are not included in the simulations, so that the accuracy on them is unknown. This is why we group them together into a single file. For unclassified breakpoints, the "SVTYPE" indicates which is the orientation of the two breakends (there are 4 possible orientations, and the "SVTYPE" is different depending on if the two breakends are in the same chromosome or not). "#CHROM" - "POS" indicate one breakend and "#CHR" - "END" the other. "START" is -1 for such unclassified breakpoints. Complex inverted SVs represent variants were a region (indicated by "CHR2", "START", "END") is copied ("SVTYPE" is CVD) or cut ("SVTYPE" is CVT (same chromosome) or IVT (different chromosomes)), inverted and inserted into "#CHROM"-"POS".


## Resource consumption