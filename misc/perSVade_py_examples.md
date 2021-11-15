Below are some examples of different analyses that can be done with perSVade. We STRONGLY ADVICE that you check the meaning of each of the indicated arguments. We note that these commands are related to the traditional installation. If you want to use the docker or singularity installation you need to do it as in the section **Quick start**. In addition, the section **FAQs** will also help why some arguments are specified.

- Traditional variant calling pipeline (small variants and coverage per gene, without SV or CNV calling), with all callers and filtering out variants with <20 reads supporting them. The consideration of repeats will be avoided:

`./scripts/perSVade.py --ref reference_genome.fasta --threads 4 -o ./output_directory -f1 reads_FWD.fastq.gz -f2 reads_FWD.fastq.gz --mitochondrial_chromosome chr_mito --mitochondrial_code 3 --gDNA_code 12 -gff features.gff --type_variant_calling small_vars --caller all --coverage 20 --ploidy 2 --skip_repeat_analysis`.

- SV and read depth-based CNV calling (on bins of 300 bp) and annotation personalizing the number of simulations. Parameter optimisation will be ran on random SV simulations. There will be a quality control and trimming of the reads. The consideration of repeats will be avoided:

`./scripts/perSVade.py --ref reference_genome.fasta --threads 4 -o ./output_directory -f1 reads_FWD.fastq.gz -f2 reads_FWD.fastq.gz --mitochondrial_chromosome chr_mito --mitochondrial_code 3 --gDNA_code 12 -gff features.gff  --coverage 20 --ploidy 2 --skip_repeat_analysis --nvars 50 --nsimulations 2 --simulation_ploidies diploid_hetero --range_filtering_benchmark theoretically_meaningful --QC_and_trimming_reads --min_chromosome_len 5000 --window_size_CNVcalling 300 --cnv_calling_algs HMMcopy,AneuFinder`


- SV, read depth-based CNV (on bins of 300 bp) and small variant calling. There will be a quality control and trimming of the reads. The consideration of repeats will be avoided:

`./scripts/perSVade.py --ref reference_genome.fasta --threads 16 -o ./output_directory -f1 reads_FWD.fastq.gz  -f2 reads_FWD.fastq.gz --mitochondrial_chromosome  chr_mito -gff annotations.gff --run_smallVarsCNV --caller all --coverage 20 --mitochondrial_code 4  --gDNA_code 12 --ploidy 2 --remove_smallVarsCNV_nonEssentialFiles --min_chromosome_len 100000 --verbose --nvars 50 --nsimulations 2 --simulation_ploidies auto --range_filtering_benchmark theoretically_meaningful --min_CNVsize_coverageBased 600 --window_size_CNVcalling 300 --cnv_calling_algs HMMcopy,AneuFinder --skip_repeat_analysis --QC_and_trimming_reads`


- SV and read depth-based CNV calling (on bins of 300 bp) and annotation without parameter optimisation (with `--fast_SVcalling`), and SV calling will use some user-defined parameters (from `--parameters_json_file perSVade_parameters.json`). Note that [here](https://github.com/Gabaldonlab/perSVade/blob/master/misc/default_perSVade_parameters.json) you can find a template .json file (with the default parameters). There will be a quality control and trimming of the reads. The consideration of repeats will be avoided:

`./scripts/perSVade.py --ref reference_genome.fasta --threads 4 -o ./output_directory -f1 reads_FWD.fastq.gz -f2 reads_FWD.fastq.gz --mitochondrial_chromosome chr_mito --mitochondrial_code 3 --gDNA_code 12 -gff features.gff --ploidy 2 --skip_repeat_analysis --QC_and_trimming_reads --min_chromosome_len 5000 --window_size_CNVcalling 300 --cnv_calling_algs HMMcopy,AneuFinder --fast_SVcalling --parameters_json_file perSVade_parameters.json`


- Run RepeatModeller and RepeatMasker to obtain the repeats annotation in your reference genome (For BSC users: this has only been tested in MareNostrum4):

`./scripts/perSVade.py --ref reference_genome.fasta --threads 16 -o ./output_directory -f1 skip -f2 skip --mitochondrial_chromosome  no_mitochondria --StopAfter_repeatsObtention`