
These are all the important files that can be generated with the script perSVade.py (check the section **FAQs** for clarification):

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