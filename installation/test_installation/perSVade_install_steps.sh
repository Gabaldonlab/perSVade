#!/bin/bash

# This recapitulates all the steps to create the env (on conda 4.8.0)

# conda distribution 'Anaconda3-2019.03-Linux-x86_64.sh', followed by 'conda install conda=4.8.0'

# remove previous env
#conda activate base
#conda remove --name perSVade_env_ --all

# You can download anaconda from 
wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
bash Anaconda3-2020.07-Linux-x86_64.sh

# the tested version is 
# wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh

# you may have to install the following commands to compile binaries: 

#conda install gcc_linux-64
#conda install gxx_linux-64

#################### PERSVADE ENV ####################

# create and activate the env
conda create -n perSVade_env python=3.6 &&
conda activate perSVade_env &&

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
conda install -c bioconda svim=1.4.2 &&
conda install -c bioconda sniffles=1.0.12

conda install -c anaconda spyder # we may skip this one



# export env
conda env export --no-builds --from-history -n perSVade_env --file perSVade_env.yml


# things that don't work
#conda install -c bioconda cnvpytor=1.0 # a package to run CNV calls. It does not install any R pack
#conda install -c conda-forge -y python-igraph=0.8.3 && # this updates many things. It may be necessary to stick to python-igraph=0.7.1.post7. This 
#conda install -c anaconda -y statsmodels=0.12.0 &&
#conda install -c bioconda bioconductor-aneufinder=1.16.0 # it updates a couple things
#conda install -c anaconda -y cython=0.29.21 # this is necessary for the lowess package
#conda install -c conda-forge -y python-igraph=0.7.1.post7 # this have errors with the subgraph
#conda install -c conda-forge -y igraph=0.7.1
#conda install -c conda-forge -y pycairo=1.20.0 # this is to plot igraph. It updates a lot of things so that it may not be worth
#conda install -c conda-forge -y igraph=0.7.1 && # updates openssl pkgs/main::openssl-1.1.1g-h7b6447c_0 --> conda-forge::openssl-1.1.1g-h516909a_1
#conda install -c conda-forge -y python-igraph=0.7.1.post7 &&
#conda create -n perSVade_env_R
# install R packages
#conda install -c conda-forge -y r-base=4.0.2 &&
#conda install -c conda-forge -y r-argparser=0.6 && # at some point this gave problems. The initially installed package was 0.4, which gave inconsistencies
#conda install -c bioconda -y bioconductor-rsvsim=1.28 &&
#conda install -c conda-forge -y r-emdbook=1.3.12 &&
#conda install -c bioconda -y bioconductor-rtracklayer=1.48.0 &&
#conda install -c conda-forge -y r-r.utils=2.9.2 &&
#conda install -c bioconda -y bioconductor-structuralvariantannotation=1.4.0;
# conda install -y gcc_linux-64 &&
#conda install -y gxx_linux-64
#conda install -c bioconda -y repeatmasker=4.0.9_p2 && # downgrades: gmp 6.2.0-he1b5a44_2 --> 6.1.2-hf484d3e_1000
#conda install -c bioconda -y repeatmodeler=2.0.1 &&
# conda install -c asmeurer -y glibc=2.19 # This was necessary to run in Nord3
#conda config --remove channels defaults
#conda config --remove channels r
#conda config --remove channels biocore
#conda install -c bioconda -y picard=2.18.26 &&
#conda install -c conda-forge -y r-base=4.0.2 # originally: conda install -c conda-forge -y r-base=4.0.2
# The following packages will be UPDATED:
#   gsl                                        2.5-h294904e_1 --> 2.6-h294904e_0
#   openssl               anaconda::openssl-1.1.1g-h7b6447c_0 --> conda-forge::openssl-1.1.1g-h516909a_1
# The following packages will be SUPERSEDED by a higher-priority channel:
#   ca-certificates     anaconda::ca-certificates-2020.6.24-0 --> conda-forge::ca-certificates-2020.6.20-hecda079_0
#   certifi                anaconda::certifi-2020.6.20-py36_0 --> conda-forge::certifi-2020.6.20-py36h9f0ad1d_0

# The following packages will be DOWNGRADED:

#   bcftools                                   1.9-h68d8f2e_9 --> 1.9-ha228f0b_4
#   libdeflate                                 1.2-h516909a_1 --> 1.0-h14c3975_1

# The following NEW packages will be INSTALLED:

#   colorlover         conda-forge/noarch::colorlover-0.3.0-py_0
#   cufflinks-py       conda-forge/noarch::cufflinks-py-0.13.0-py_1

# The following packages will be UPDATED:

#   certifi                anaconda::certifi-2020.6.20-py36_0 --> conda-forge::certifi-2020.6.20-py36h9880bd3_2
#   openssl               anaconda::openssl-1.1.1g-h7b6447c_0 --> conda-forge::openssl-1.1.1h-h516909a_0

# The following packages will be SUPERSEDED by a higher-priority channel:

#   ca-certificates     anaconda::ca-certificates-2020.7.22-0 --> conda-forge::ca-certificates-2020.6.20-hecda079_0

# the environment can be exported as:
#conda env create --file perSVade_env.yml --name perSVade_env

# it can be setup with 
# ./setup_environment.sh


# R packages
#conda install -c r rstudio=1.1.456  # this is the last version. This takes too long.

# buggy packages:
#conda install -c bioconda -y qualimap=2.2.2d && (I can instead install 2.2.2, but it fails). I decide to skip this package

# GENERAL NOTES
# conda install -c bioconda entrez-direct=13.3 would be to install esearch and efetch, but they are already installed with VEP
# at the beginning, the installation of the conda packages is not enough to be able to load them from Rscript

# conda install -c r r-stringr=1.4.0 was already installed


# OLD VERSIONS OF THE R PACKAGES
# installed
# conda install -c conda-forge -y r-argparser=0.4 && 
# # installed in r42 (0.4)
# conda install -c conda-forge -y r-emdbook=1.3.11 && 
# # installed in r75 (1.3.11)
# conda install -c bioconda -y bioconductor-rsvsim=1.26.0 && 
# # installed in r74 (1.24.0), updated in r92 (1.26.0)

# conda install -c bioconda -y bioconductor-rtracklayer=1.46.0 && 
# # installed in r9 (1.34.2), updated in r15 (1.44.2), removed in r26 (1.44.2), reinstalled in r27 (1.44.2), removed in r57 (1.44.2), readed in r58 (1.44.2), removed in r59 (1.44.2), added in r60 (1.36.6), updated in r62 (1.44.2), updated in r92 (1.46.0)

# conda install -c conda-forge -y r-r.utils=2.8.0 && 
# # installed in r43 (2.9.0), downgraded in r55 (2.8.0)

# conda install -c conda-forge -y r-base=3.6.1 
# # installed in r9 (3.3.2), updated in r15 (3.6.1)

# conda install -c bioconda -y bioconductor-variantannotation=1.32.0 &&
# # intsalled in r24 (1.30.1), updated in r92 (1.32.0)

# conda install -c bioconda -y bioconductor-structuralvariantannotation=1.2.0 &&
# # same as bioconductor-variantannotation


#conda install -c conda-forge -y igraph=0.7.1 && # updates openssl pkgs/main::openssl-1.1.1g-h7b6447c_0 --> conda-forge::openssl-1.1.1g-h516909a_1
#conda install -c conda-forge -y python-igraph=0.7.1.post7 &&

# the ete installation may give rise to problems
#conda install -c etetoolkit -y ete3=3.1.1 && # at some point 3.1.2 was also installed (I try if ete 3.1.1 does not give segmentation errors). This avoids a previous error that was not working.
# I will install ete3 in a sepparate environment

#conda config --add channels biocore &&

#conda config --add channels etetoolkit &&


# conda config --add channels r &&



# Ninja installation 
#git clone https://github.com/TravisWheelerLab/NINJA
#cd NINJA
#make build

#wget https://github.com/TravisWheelerLab/NINJA/archive/0.95-cluster_only.tar.gz
#tar -xvf 0.95-cluster_only.tar.gz; rm 0.95-cluster_only.tar.gz; cd NINJA-0.95-cluster_only/NINJA
#make all

#wget https://github.com/TravisWheelerLab/NINJA/archive/0.97-cluster_only.tar.gz
#tar -xvf 0.97-cluster_only.tar.gz; rm 0.97-cluster_only.tar.gz; cd NINJA-0.97-cluster_only/NINJA
#make all

#ln -s ~/anaconda3/bin/x86_64-conda-cos6-linux-gnu-gcc ~/anaconda3/bin/gcc
#ln -s ~/anaconda3/bin/x86_64-conda_cos6-linux-gnu-g++ ~/anaconda3/bin/g++

#chmod 755 ~/anaconda3/bin/gcc
#chmod 755 ~/anaconda3/bin/g++



