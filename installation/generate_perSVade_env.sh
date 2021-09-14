#!/bin/bash
set -e

# this script has all the steps to create the perSVade_env
echo 'generating perSVade_env from an environment that has already python 3.6...'

# add channels
conda config --add channels conda-forge &&
conda config --add channels bioconda &&
conda config --add channels plotly &&
conda config --add channels anaconda &&

# install packages that should be loaded 
conda install -y pandas=0.24.2 &&
conda install -y biopython=1.73 &&
conda install -y scipy=1.4.1 &&
conda install -y scikit-learn=0.21.3 &&
conda install -y -c plotly plotly=2.7.0 &&
conda install -y -c conda-forge cufflinks-py=0.13.0 &&
conda install -c anaconda -y seaborn=0.9.0 && # updates: ca-certificates    conda-forge/label/cf201901::ca-certif~ --> anaconda::ca-certificates-2020.6.24-0, certifi  conda-forge/label/cf201901::certifi-2~ --> anaconda::certifi-2020.6.20-py36_0
conda install -c conda-forge -y matplotlib=3.3.0 && # by default, seaborn installs 3.3.1, which raises errors
conda install -c anaconda -y psutil=5.7.2 &&
conda install -c anaconda -y pigz=2.4 &&

# install packages related to software
conda install -c bioconda -y bwa=0.7.17 &&
conda install -c bioconda -y bcftools=1.9 &&
conda install -c bioconda -y samtools=1.9 && # this will downgrade: ncurses.2-he1b5a44_1 --> 6.1-hf484d3e_1002 python 3.6.11-h425cb1d_1_cpython --> 3.6.10-h8356626_1011_cpython readline  8.0-he28a2e2_2 --> 8.0-h46ee950_1
conda install -c bioconda -y bedtools=2.29.0 &&
conda install -c bioconda -y wgsim=1.0 &&
conda install -c bioconda -y seqtk=1.3 &&
conda install -c bioconda -y gatk4=4.1.2.0 && # downgrades openjdk  11.0.1-hacce0ff_1021 --> 8.0.192-h516909a_1005 (the old openjdk was 8.0.152)
conda install -c bioconda -y freebayes=1.3.1 &&
conda install -c bioconda -y mosdepth=0.2.6 &&
conda install -c bioconda -y ensembl-vep=100.2 && # downgrades libtiff 4.1.0-hc7e4089_6 --> 4.1.0-hc3755c2_3. libwebp-base-1.1.0-h516909a_3 is removed
conda install -c bioconda -y vcflib=1.0.0_rc2 &&
#conda install -c bioconda -y sra-tools=2.10.0 && # this was originally throught to work
conda install -c bioconda -y sra-tools=2.10.9 && # this solves some problems with prefetch

conda install -c bioconda -y trimmomatic=0.38 &&
conda install -c bioconda -y parallel-fastq-dump=0.6.3 &&
conda install -c bioconda -y fastqc=0.11.9 &&
conda install -c bioconda -y bedops=2.4.39 &&
conda install -c bioconda -y genmap=1.3.0 &&
conda install -c anaconda -y cython=0.29.21 && # necessary for the cylowess running
conda install -c anaconda -y xlrd=1.2.0 &&
conda install -c bioconda -y entrez-direct=13.3 && # would be to install esearch and efetch, but they are already installed with VEP
conda install -c bioconda -y porechop=0.2.4 &&
conda install -c bioconda -y svim=1.4.2 &&
conda install -c bioconda -y sniffles=1.0.12 &&

echo 'SUCCESS: perSVade_env was correctly generated'

