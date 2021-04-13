# perSVade: personalized Structural Variation detection

perSVade is a method for structural variation (SV), small variant (SNPs and IN/DELs) and read depth-based Copy Number Variation (CNV) calling and annotation with a single bash command using a set of paired-end short reads as input. The novel SV calling pipeline finds breakpoints with  GRIDSS (https://github.com/PapenfussLab/gridss) and summarizes them into complex structural variants with CLOVE (https://github.com/PapenfussLab/clove), with some added features. perSVade provides an automated benchmarking and parameter selection for these methods in any genome or sequencing run. This is useful for species without reccommended running and filtering parameters. In addition, it provides an automated report of the SV calling accuracy on simulations and real data, useful to assess the confidence of the results. 

## Pipeline overview

<img src="https://github.com/Gabaldonlab/perSVade/blob/master/misc/perSVade_pipeline_cartoon.png" width="150%" height="150%">

perSVade runs `bwa mem` to align the short reads, generating a sorted .bam file that is the core of several downstream analyses. These include:

- SV calling. This is the core, most novel function  of perSVade. It uses `gridss` to infer a list of breakpoints (two regions of the genome (two breakends) that are joined in the sample of interest and not in the reference genome) from discordant read pairs, split reads and de novo assembly signatures. The breakpoints are summarized into SVs with `clove`. In order to find the best parameters to run these algorithms, perSVade generates two simulated genomes with 50 SVs of each type and tries several combinations (>13,000,000,000) of filters for the `gridss` and `clove` outputs. It selects the filters that have the highest Fvalue (harmonic mean between precision and recall) for each simulated genome and SV type. In order to reduce overfitting, perSVade selects a final set of "best parameters" that work well for all simulations and SV types. This set of best parameters is used for the final calling on the real data. The accuracy (Fvalue, precision, recall) of these parameters on each simulation and SV type is reported as a heatmap, which is useful to evaluate the expected calling accuracy. By default, the simulations are placed randomly across the genome. However, SVs often appear around repetitive elements or regions of the genome with high similarity (i.e.: transposable elements insertions). This means that random simulations may not be realistic, potentially leading to overestimated calling accuracy and a parameter selection that does not work well for real SVs. perSVade can also generate more realistic simulation occuring around regions with known SVs (i.e. regions with SVs called with perSVade) or homologous regions (inferred from BLAST). See "Output" for more details.

- Read depth-based Copy Number Variation (CNV) calling. CNVs are one type of SVs where there is an alteration in the genomic content (deletions or duplications). The SV calling feature of perSVade (point 1) identifies some CNVs (insertions, tandem duplications, deletions and complex inverted SVs) but it can miss others (i.e.: whole-chromosome duplications or regions with unknown types of rearrangements yielding CNVs). perSVade also includes a pipeline to call CNVs from read-depth alterations. For example, regions with 0x or 2x read-depth as compared to the mean of the genome can be called duplications, or deletions, respectively. A straight forward implementation of this concept to find CNVs is challenging because many genomic features drive variability in read depth independently of CNV. In order to solve this, perSVade calculates the relative coverage for bins of the genome and corrects the effect of the GC content, mappability and distance to the telomere (using non-parametric regression as in https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-014-0100-8). This corrected coverage is used by `CONY`, `AneuFinder` and/or `HMMcopy` to call CNVs across the genome. perSVade generates consensus CNV calls from the three programs taking always the most conservative copy number for each bin of the genome. For example, if the used programs disagree on the copy number of a region the closest to 1 will be taken as the best estimate.

- Small variant (SNPs and IN/DELs) and gene CNV calling. perSVade includes an option to run a pipeline that performs small variant calling and calculation of read depth for each gene. It runs any of `freebayes`,  `GATK HaplotypeCaller` and/or `bcftools call` for small variant calling and integrates the results into .tab and .vcf files. It also runs `mosdepth` to get the coverage for this gene (it requires a .gff file from the user).

- Integration of read-depth based CNVs and SVs into a single .vcf file. This is a file that is focused on showing the alteration of SVs on specific genomic regions (see the section "Output" for more details). It also removes redundant calls between the CNVs identified with `gridss`+`clove` and those derived from


## Installation

These are the steps to install perSVade. If you are running it in the BSC (internal use) you can skip this (see **Running in BSC clusters**).

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

Make sure that none of them exist before running this script. You can change `<env_name>` to fullfil this requirement.

We note that this was tested with `conda 4.8.0` on a Linux-x86 64-bit architecture, installed at 03/2019. If you have a different conda version, you may change a bit the perSVade_env.yml file so that it does not create dependency problems.

### 4. Test installation

We highly recommend to test that all dependencies were properly installed with the following commands:

`./installation/test_installation/test_installation.py`

This process should take arround 45 minutes on 4 cores. Verify that it finishes with the following message:

`SUCCESS: perSVade was properly installed`

There is a WARNING message that you should look for after running this script:

`WARNING: The connection to SRA did not work`. perSVade includes the option to query and download from the SRA database for the benchmarking of SV calling. This requires a proper network access and SRA connection, which may not always be available. This warning indicates that this process is not possible on your machine. You can skip this connection by providing the reads on your own through `--close_shortReads_table`.


## Running in BSC clusters (internal use only)

If you are working from any cluster that has access to the BSC /gpfs filesystem you can activate the perSVade environment from its location in mschikora. YOU DON'T NEED TO INSTALL ANYTHING.

`source /gpfs/projects/bsc40/mschikora/anaconda3/etc/profile.d/conda.sh`  # activate the conda environment of mschikora

`conda activate perSVade_v0.7_env` # activate the environment of perSVade version 0.7. You can change the version

You can next run perSVade from the releases folder (these are stable versions of the pipeline). For example:

`python /gpfs/projects/bsc40/mschikora/scripts/perSVade/releases/perSVade-0.7/scripts/perSVade.py -r <path to the reference genome (fasta)> -o <output_directory> -p <ploidy, 1 or 2> -f1 <forward_reads.fastq.gz> -f2 <reverse_reads.fastq.gz> `

IMPORTANT NOTE: The activation of the perSVade conda environment works well from version 0.7 on. This means that you can activate from the login of MN or interactive nodes. However, the activation of older versions (v0.4 and below) is costly, and it overloads the login nodes. If you want to use an old version of perSVade, always activate it on an interactive node (i.e.: `salloc`). In addition, you can't run the perSVade pipeline from a login, because it takes too many resources. You can submit perSVade as a job or run from an interactive session with `salloc -n 1 --time=02:00:00 -c 48 --qos debug`.

OTHER NOTES:

- All perSVade modules work in MareNostrum, and most of them in Nord3. The generation of "repeat files" and finding of "homologous regions" does not work in Nord3 because the GCLIB versions are old.

## Quick start

You can call the perSVade pipeline with:

`conda activate perSVade_env`

`python ./scripts/perSVade.py -r <path to the reference genome (fasta)> -o <output_directory> -p <ploidy, 1 or 2> -f1 <forward_reads.fastq.gz> -f2 <reverse_reads.fastq.gz>`

By default, perSVade runs SV calling optimised on random simulations, read depth-based CNV calling with AneuFinder/HMMcopy and integration of SVs and CNVs into a single vcf.

However, this is a flexible pipeline with several options. Type `./scripts/perSVade.py -h` to understand what is the meaning of these options. Note that there are some of them which are flagged as "ONLY TO BE USED FOR DEVELOPMENT. USE AT YOUR OWN RISK!". These are arguments that are useful for developing the pipeline, but are not intended for most users.

In addition, you can import the file `./scripts/sv_functions.py` as a python module, which is useful to run some specific functions from a python script/shell. This module can be loaded with these python commands:

`import sys; sys.path.append("./scripts")`

`import sv_functions as fun`


## Examples

Below are some examples of different analyses that can be done with perSVade. We STRONGLY ADVICE that you check the meaning of each of the indicated arguments. In addition, the section **FAQs** will also help why some arguments are specified.

- Traditional variant calling pipeline (small variants and coverage per gene, without SV or CNV calling):

`./scripts/perSVade.py --ref reference_genome.fasta --threads 4 -o ./output_directory -f1 reads_FWD.fastq.gz -f2 reads_FWD.fastq.gz --mitochondrial_chromosome chr_mito --mitochondrial_code 3 --gDNA_code 12 -gff features.gff  --run_smallVarsCNV --skip_SVcalling --skip_CNV_calling --caller all --coverage 20 --ploidy 2 --remove_smallVarsCNV_nonEssentialFiles --skip_repeat_analysis`

- SV and read depth-based CNV calling (on bins of 300 bp) and annotation personalizing the number of simulations. Parameter optimisation will be ran on random SV simulations. There will be a quality control and trimming of the reads:

`./scripts/perSVade.py --ref reference_genome.fasta --threads 4 -o ./output_directory -f1 reads_FWD.fastq.gz -f2 reads_FWD.fastq.gz --mitochondrial_chromosome chr_mito --mitochondrial_code 3 --gDNA_code 12 -gff features.gff  --coverage 20 --ploidy 2 --skip_repeat_analysis --nvars 50 --nsimulations 2 --simulation_ploidies diploid_hetero --range_filtering_benchmark theoretically_meaningful --QC_and_trimming_reads --min_chromosome_len 5000 --window_size_CNVcalling 300 --cnv_calling_algs HMMcopy,AneuFinder`

## Output 

perSVade outputs several files depending on the arguments. Before checking the output make sure that there is a file in the output directory called `perSVade_finished_file.txt`, which indicates that the pipeline worked and includes a report of the timing. 

These are all the important files that can be generated with perSVade (check the section **FAQs** for clarification):

- `aligned_reads.sorted.bam` and `aligned_reads.sorted.bam.bai`: the aligned reads sorted by genomic coordinate.

- `aligned_reads.bam.sorted.flagstat` is a report of the performance of the read alignment.

- `reference_genome_dir` is a directory that contains files related to the provided genome and annotations. This may be removed if necessary.

- `aligned_reads.bam.sorted.calculating_windowcoverage/coverage_windows_<n_nucleotides>bp.tab` contains a table with the coverage for windows of the genome that are as long as 5% of the median chromosome length (n_nucleotides). This is the output of mosdepth. This table contatins the following fields:

    - `#chrom` is the chromosome name.
    - `start` is the 0-based start coordinates of the region
    - `end` is the 1-based end coordinate of the region
    - `length` is the length of the region
    - `mediancov_1` is the median read depth in the region
    - `percentcovered_1` is the perecentage of the region that is covered with 1 read or more
    - `nocoveragebp_1` is the number of bases that have no coverage

- `SVdetection_output` contains the output of the SV calling and parameter optimisation. There can be two folders inside:

    - `parameter_optimisation` contains files that report the accuracy on simulations used for the parameter optimisation. This is only created if such parameter optimisation was performed. These are the relevant files:

        - `SVfiles` includes files of the simulated SVs. They are formated as explained in the FAQ **How are the SVs encoded into single files**.

        - `plots/cross_accuracy_heatmaps/` includes heatmaps that show how the best parameters for one simulation and SV type (training parameters in the rows) work for the other simulations and SV types (testing simulations in the columns). There is one plot for each accuracy measurement (recall, precision and Fvalue). Only testing svtypes that have at least 5 SVs simulated are shown. perSVade takes as a final set of "best (training) parameters" those that work best across all the testing simulations and SVs. The heatmap shows the accuracy (as an annotation on the plot) for the chose set of "best parameters". These plots are useful to infer the SV calling accuracy in the input sample and assess the overfitting associated to parameter optimization.

        - `benchmarking_all_filters_for_all_genomes_and_ploidies/df_cross_benchmark_best.tab` is a .tab file indicating the accuracy (precision, recall, Fvalue) of the final chosen set of "best parameters" when tested on the different simulations and SVs. These are the relevant columns:

            - `test_ploidy`, `test_simName` and `test_svtype` indicate the testing simulation and svtype. For example, if you ran parameter optimisation for 2 simulations you expect `test_simName` to have "simulation_1" and "simulation_2". If you ran for only one ploidy (i.e.: haploid) you expect `ploidy` to be always "haploid".

            - `FN`, `FP`, `TP`, `nevents` are the numbers of false negatives, false positives, true positives and number of all real positives, respectively.

            - `precision`, `recall` and `Fvalue` are the resulting accuracy.

    - `final_gridss_running` includes the output files of the SV calling run with the optimised parameters. These are the relevant files:

        - `deletions.tab`,  `inversions.tab`,  `tandemDuplications.tab`,  `insertions.tab`,  `translocations.tab`,  `unclassified_SVs.tab`: These are the called SVs. Each SV type is saved into a different files.  Check the FAQ **How are the SVs encoded into single files** for more information on these files.

        - `gridss_output.raw.vcf` and `gridss_output.filt.vcf` are the raw and filtered (with optimised parameters) outputs of `gridss`, respectively. These are a set of breakpoints encoded in a .vcf file. Check the `gridss` documentation (https://github.com/PapenfussLab/gridss) for more details.

        -  `gridss_output.filt.bedpe` is a representation in the .bedpe format of the breakpoints from `gridss_output.filt.vcf`. The column names are "chrom1", "start1", "end1", "chrom2", "start2", "end2", "breakpoint ID", "QUAL", "strand 1", "strand 2". Each row indicates a breakpoint between two breakends (1 and 2). Each breakend has a "start" and "end" (for example "start1" and "start2") because the exact position is not always resolved (and hence the range indicates a region where the breakend is found).

        - `clove_output.vcf` is the output of clove without any coverage-based filtering of tandem duplications and deletions. Check https://github.com/PapenfussLab/clove for more details.

        - `perSVade_parameters.json` is a file that indicates the used parameters for this optimised perSVade run. You can pass these to be used as parameters for running perSVade in another sample (see FAQ **perSVAde's parameter optimisation is too slow. I just want to run SV calling with some pre-defined filters. How can I do this?**)

- `CNV_calling` contains the output of the read depth-based CNV calling. These are the relevant files:

    - `final_CNVcalling.tab` includes the regions with a copy number different from 1.  These are the relevant fields (some of them depend on the chosen CNV calling algorithms specified with the option `--cnv_calling_algs`):

        - `chromosome`, `start` and `end` indicate the region under CNV.

        - `median_relative_CN_AneuFinder`, `median_relative_CN_HMMcopy` and `median_relative_CN_CONY` indicate the predicted copy number by each of the three CNV callers. This value is relative to 1.0 (absence of CNV). In a diploid genome, a complete deletion would be 0.0, a monosomy 0.5, a trisomy 1.5 and a tetraploidy 2.0.

        -  `merged_relative_CN` is the output of merging `median_relative_CN_AneuFinder`, `median_relative_CN_HMMcopy` and `median_relative_CN_CONY`. This merging is done in a conservative way, so that the value that is closer to 1.0 (no CNV) is always taken.

        - `median_coverage` is the median read depth of this region.

        - `median_coverage_corrected` is the coverage corrected by mappability, GC content and distance to the telomere. This is the coverage that is used to run CNV calling.

- `SVcalling_output` contains the output of merging the SV and read depth-based CNV into a single .vcf file. These are the relevant files:

    - `SV_and_CNV_variant_calling.vcf` is the merged .vcf file. Check the FAQ **What is in SV_and_CNV_variant_calling.vcf?** for more information.

    - `SV_and_CNV_variant_calling.vcf_annotated_VEP.tab` is the functional annotation with `VEP`, where each line corresponds to a particular alteration in a given gene. You should check https://www.ensembl.org/info/docs/tools/vep/online/results.html for more details on how is the annotation.

- `smallVars_CNV_output` contains some folders and files related to the small variant calling and calculation of coverage per gene:

    - `bcftools_ploidy2_out`, `HaplotypeCaller_ploidy2_out` and `freebayes_ploidy2_out` are folders that contain the raw and filtered vcf files of each fo the programs.
    - `merged_vcfs_allVars_ploidy2.vcf` is a vcf file with the merged output of the algorithms used in the variant calling. You can check the header to understand what are all the tags in the `INFO` field. Here, the multiallelic loci are split to ease the analysis. Note that the `ID` field corresponds to the `#Uploaded_variation` field of `variant_annotation_ploidy2.tab`.
    - `variant_calling_ploidy2.tab` is a tabular version of `merged_vcfs_allVars_ploidy2.vcf`. Each column contains information parsed from the vcf file. This file is easier to manage because it has all the info in tabular format.
    - `variant_annotation_ploidy2.tab` is the output of VEP, where each line corresponds to a particular alteration in a given gene.
    - `variants_atLeast<n_PASS_programs>PASS_ploidy2.vcf` are the variants that PASS de filters in at least n_PASS_programs algorithms. The INFO and FORMAT fields are simplified to take less space. It may be useful to take, for example, variants that PASS de filters by at least 2 algorithms.

    - `CNV_results/genes_and_regions_coverage.tab` is a table that contains the coverage of all the genes in the provided gff. These are the fields:

        - `chromosome`, `start`, `end` and `length` are the coordinates of the gene
        - `median_reads_per_gene` is the median read depth for this gene
        - `nocoveragebp_1` is the number of bases that have no coverage
        - `percentcovered_1` is the perecentage of the region that is covered with 1 read or more
        - `ID` is the gene ID in the gff (parsed from the attributes). If there are some duplicated IDs with parts it corresponds to the union of ID and part.
        - `fraction_covered_by_MoreThan1read` is equivalent to `percentcovered_1`, but from 0 to 1.
        - `relative_coverage` is the `median_reads_per_gene` divided by the median of all `median_reads_per_gene`.
        - All the fields with a `_+-10kb_region` sufix are equivalent to the ones that do not have it, but for a region that starts at -10kb of the gene start and ends at +10kb of the gene end. This can be useful for CNV analysis.

## FAQs

### How are repetitive elements considered in perSVade and why is it important?

Read mapping around repetitive elements of the genome can be innacurate, potentially yielding false positive SVs around them. Conversely,  there are some SVs expected to appear around repeats (i.e.: those derived from transposable elements insertions). This means that the decision on how to handle SVs around such regions is not trivial. By default, perSVade runs `RepeatModeler` and `RepeatMasker` to annotate repetitive regions, which is used to assess whether removing SV calls overlapping these elements increases the overall accuracy. If so, they are removed.

This can be problematic for two reasons. First, the prediction of repeats is time-consuming for large genomes. Second, some repeat families (i.e.: simple repeats, low complexity regions) may yield more false positive calls than others (i.e.: large transposons), so that treating all repetitive elements together may not be always the best option.

perSVade includes two options to circumvent this:

1. Skipping the analysis of repeats (`--skip_repeat_analysis`)

2. Providing a .tab file that has manually-curated regions with repeats to be discarded (`--previous_repeats_table`). Note that this option is also useful if you want to predict only once the repeats and use the resulting table for many samples with the same genome.

It is also possible to annotate whether small variants (when running with `run_smallVarsCNV`) are found around repeats (using `--consider_repeats_smallVarCall`). However, if you just want to get small variants calls with no repeat annotation we recommend using `--skip_repeat_analysis` to avoid excessive, useless, computation time.

### perSVade crashed due to insufficient resources (allocated RAM or time), do I have to repeat all the running?** 

No worries. perSVade is designed to restart from the last correct step if you re-run with the same command in the same output directory. This means that you can simply re-run on a system with higher resources and expect no steps to be repeated.

In addition, you may consider readjusting the balance between allocated threads (with `--threads`) and RAM (depends on the running system). Some programs used by perSVade require a minimum amount of RAM/thread, which depends on the input. We have developed perSVade to balance this automatically, but we can't guarantee that it will always work (specially if input genomes/reads are much larger than those tested). 

Another option that can be useful to deal with memory errors is to specify the fraction of RAM that is being allocated to a given perSVade run. In several steps, the pipeline needs to calculate the available memory (using the python command psutil.virtual_memory()) to allocate resources accordingly. This command returns all the available memory in the computer. Sometimes you may be running on a fraction of the computers' resources (for example if you are using a fraction of a complete node of a computing cluster). If this is the case the psutil calculation will  overestimate the available RAM. If that is the case you can provide the fraction available through the argument `--fraction_available_mem`. By default, it will calculate the available RAM by filling the memory, which may give errors. It is highly reccommended that you provide this option. If you want to use all the allocated memory you should specify `--fraction_available_mem 1.0`. 

IMPORTANT NOTE: If you are running perSVade in MareNostrum4 or Nord3 you should skip the `--fraction_available_mem` argument, as perSVade already calculates it.

### How are the SVs encoded into single files?

perSVade generates a single .tab file for each type of SV (inversions, deletions, tandemDuplications, insertions, translocations and unclassified SVs) into the folder `SVdetection_output/final_gridss_running`. Each of these has a unique set of columns that represent the variants. This is the meaning of these columns (note that text in "" indicates the column names of the corresponding .tab files):

<img src="https://github.com/Gabaldonlab/perSVade/blob/master/misc/simple_SVs.png" width="40%" height="40%">

- Simple SVs: deletions, inversions and tandemDuplications (duplication of a region which gets inserted next to the affected region). These are described by a chromosome ("Chr"), "Start" and "End" coordinates of the SV. perSVade outputs one .tab file for each of these SV types.

<img src="https://github.com/Gabaldonlab/perSVade/blob/master/misc/insertions.png" width="40%" height="40%">

- Insertions: a region of the genome (indicated by "ChrA", "StartA", "EndA") is copied ("Copied" is TRUE) or cut ("Copied" is FALSE) and inserted into another region (indicated by "ChrB", "StartB"). "EndB" comes from adding to "StartB" the length of the inserted region. There is one .tab file for insertions.

<img src="https://github.com/Gabaldonlab/perSVade/blob/master/misc/translocations.png" width="40%" height="40%">

- Translocations: balanced translocations between two chromosomes ("ChrA" and "ChrB"). "StartA" indicates the start of "ChrA" (position 0). "EndA" indicates the position in "ChrA" where the translocation happens. For non-inverted translocations, "StartB" is 0 and "EndB" is the position in "ChrB" where the translocation happens. For inverted translocations, "StartB" is the position in "ChrB" where the translocation happens and "EndB" is the last position of "ChrB". There is one .tab file for translocations.

<img src="https://github.com/Gabaldonlab/perSVade/blob/master/misc/unclassified_SVs.png" width="50%" height="50%">

- Unclassified SVs: There is one .tab file ("unclassified_SVs.tab") that reports all the variants that are called by `clove` and cannot be assigned to any of the above SV types. These include unclassified breakpoints (which could be part of unresolved/unkown complex variants) and complex inverted SVs (which are non-standard SVs). These types of SVs are not included in the simulations, so that the accuracy on them is unknown. This is why we group them together into a single file. For unclassified breakpoints, the "SVTYPE" indicates which is the orientation of the two breakends (there are 4 possible orientations, and the "SVTYPE" is different depending on if the two breakends are in the same chromosome or not). "#CHROM" - "POS" indicate one breakend and "#CHR" - "END" the other. "START" is -1 for such unclassified breakpoints. Complex inverted SVs represent variants were a region (indicated by "CHR2", "START", "END") is copied ("SVTYPE" is CVD) or cut ("SVTYPE" is CVT (same chromosome) or IVT (different chromosomes)), inverted and inserted into "#CHROM"-"POS".

### perSVAde's parameter optimisation is too slow. I just want to run SV calling with some pre-defined filters. How can I do this?

perSVade includes an option to skip the parameter optimisation: `--fast_SVcalling`. This uses some default parameters for the filtering. You can also customize the filters by providing a .json file with them through `--parameters_json_file`. The default parameters encoded as .json are in  `misc/default_perSVade_parameters.json`, and this file can be used as a template to provide custom parameters with `--parameters_json_file`.

### What is in SV_and_CNV_variant_calling.vcf?

This is a vcf that contains all called SVs and CNVs, in a way that is focused on how each SV affects particular regions of the genome. We consider that this is useful for further functional annotation. These are the important fields:

- 

- 

## Resource consumption