#!/bin/bash

# This recapitulates all the steps to create the env (on conda 4.8.0)

# conda distribution 'Anaconda3-2019.03-Linux-x86_64.sh', followed by 'conda install conda=4.8.0'

# remove previous
#conda activate base
#conda remove --name perSVade_env_test --all

# create and activate the env
#conda create -n perSVade_env_test python=3.6
#conda activate perSVade_env_test 

# add channels
conda config --add channels conda-forge &&
conda config --add channels biocore &&
conda config --add channels bioconda &&
conda config --add channels etetoolkit &&
conda config --add channels r &&

# install packages that should be loaded 
conda install -y pandas=0.24.2 &&
conda install -y biopython=1.73 &&
conda install -y scipy=1.4.1 &&
conda install -y scikit-learn=0.21.3 &&
conda install -c conda-forge -y igraph=0.7.1 && # updates openssl pkgs/main::openssl-1.1.1g-h7b6447c_0 --> conda-forge::openssl-1.1.1g-h516909a_1
conda install -c conda-forge -y python-igraph=0.7.1.post7 &&
conda install -c etetoolkit -y ete3=3.1.1 &&
conda install -c anaconda -y seaborn=0.9.0 && # updates: ca-certificates    conda-forge/label/cf201901::ca-certif~ --> anaconda::ca-certificates-2020.6.24-0, certifi  conda-forge/label/cf201901::certifi-2~ --> anaconda::certifi-2020.6.20-py36_0

# install packages related to software
conda install -c bioconda -y repeatmasker=4.0.9_p2 && # downgrades: gmp 6.2.0-he1b5a44_2 --> 6.1.2-hf484d3e_1000
conda install -c bioconda -y repeatmodeler=2.0.1 &&
conda install -c bioconda -y bwa=0.7.17 &&
conda install -c bioconda -y picard=2.18.26 &&
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
conda install -c bioconda -y sra-tools=2.10.0 &&
conda install -c bioconda -y trimmomatic=0.38 &&
conda install -c bioconda -y parallel-fastq-dump=0.6.3 &&
conda install -c bioconda -y fastqc=0.11.9 &&


conda install -c r rstudio=1.1.456 

# this is the report:

#The following packages will be UPDATED:

#  libxml2            conda-forge::libxml2-2.9.10-hee79883_0 --> pkgs/main::libxml2-2.9.10-he19cac6_1
#  libxslt            conda-forge::libxslt-1.1.33-h31b3aaa_0 --> pkgs/main::libxslt-1.1.34-hc22bd24_0

# The following packages will be SUPERSEDED by a higher-priority channel:

#  matplotlib-base    conda-forge::matplotlib-base-3.2.2-py~ --> pkgs/main::matplotlib-base-3.2.2-py36hef1b27d_0
#  qt                       conda-forge::qt-5.9.7-h0c104cb_3 --> pkgs/main::qt-5.6.3-h8bf5577_3

#The following packages will be DOWNGRADED:

#  cairo                                1.16.0-hcf35c78_1003 --> 1.16.0-h18b612c_1001
#  fontconfig                           2.13.1-h86ecdb6_1001 --> 2.13.1-he4413a7_1000
#  harfbuzz                                 2.4.0-h9f30f68_3 --> 2.4.0-h37c48d4_1
#  icu                                       64.2-he1b5a44_1 --> 58.2-hf484d3e_1000
#  igraph                                0.7.1-h9e3b1fc_1007 --> 0.7.1-h2166141_1005
#  ld_impl_linux-64                          2.34-hc38a660_9 --> 2.34-h53a641e_7
#  libgd                                 2.2.5-h307a58e_1007 --> 2.2.5-h0d07dcb_1005
#  pyqt                                 5.9.2-py36hcca6a23_4 --> 5.6.0-py36h13b7fb3_1008
#  sip                              4.19.8-py36hf484d3e_1000 --> 4.18.1-py36hf484d3e_1000

# ALL THE NEXT PACKAGES DO NOT WORK IF NOT WITH RSTUDIO (MAYBE)


conda install -c conda-forge -y r-base=3.6.1 && # updates openssl anaconda::openssl-1.1.1g-h7b6447c_0 --> conda-forge::openssl-1.1.1g-h516909a_1; packages supersed to other priority: ca-certificates certifi                anaconda::certifi-2020.6.20-py36_0 --> conda-forge::certifi-2020.6.20-py36h9f0ad1d_0 anaconda::ca-certificates-2020.6.24-0 --> conda-forge::ca-certificates-2020.6.20-hecda079_0; downgrades: curl: 7.71.1-he644dc0_4 --> 7.68.0-hf8cf82a_0 krb5: 1.17.1-hfafb76e_2 --> 1.16.4-h2fd8d38_0, libcurl 7.71.1-hcdd3856_4 --> 7.68.0-hda55be3_0

# R packages
conda install -c bioconda -y bioconductor-rsvsim=1.26.0 &&
conda install -c conda-forge -y r-argparser=0.4 &&
conda install -c conda-forge -y r-emdbook=1.3.11 &&
conda install -c bioconda -y bioconductor-rtracklayer=1.46.0 &&
conda install -c conda-forge -y r-r.utils=2.8.0 &&
conda install -c bioconda -y bioconductor-variantannotation=1.32.0 &&
conda install -c bioconda -y bioconductor-structuralvariantannotation=1.2.0 &&

conda install -c r rstudio=1.1.456 # this will update 


#The following packages will be UPDATED:

#  libxml2            conda-forge::libxml2-2.9.10-hee79883_0 --> pkgs/main::libxml2-2.9.10-he19cac6_1
#  libxslt            conda-forge::libxslt-1.1.33-h31b3aaa_0 --> pkgs/main::libxslt-1.1.34-hc22bd24_0

#The following packages will be SUPERSEDED by a higher-priority channel:

#  matplotlib-base    conda-forge::matplotlib-base-3.2.2-py~ --> pkgs/main::matplotlib-base-3.2.2-py36hef1b27d_0
#  qt                       conda-forge::qt-5.9.7-h0c104cb_3 --> pkgs/main::qt-5.6.3-h8bf5577_3

#The following packages will be DOWNGRADED:

#  cairo                                1.16.0-hcf35c78_1003 --> 1.16.0-h18b612c_1001
#  fontconfig                           2.13.1-h86ecdb6_1001 --> 2.13.1-he4413a7_1000
#  genometools-genom~                   1.6.1-py36h30d060d_2 --> 1.5.10-py36h997e34b_3
#  graphviz                                2.42.3-h0511662_0 --> 2.40.1-h0dab3d1_0
#  harfbuzz                                 2.4.0-h9f30f68_3 --> 2.4.0-h37c48d4_1
#  icu                                       64.2-he1b5a44_1 --> 58.2-hf484d3e_1000
#  igraph                                0.7.1-h9e3b1fc_1007 --> 0.7.1-h2166141_1005
#  libgd                                 2.2.5-h307a58e_1007 --> 2.2.5-h0d07dcb_1005
#  pango                                   1.42.4-h7062337_4 --> 1.40.14-he7ab937_1005
#  pyqt                                 5.9.2-py36hcca6a23_4 --> 5.6.0-py36h13b7fb3_1008
#  r-base                                   3.6.1-h3a67422_6 --> 3.6.1-h8900bf8_2
#  sip                              4.19.8-py36hf484d3e_1000 --> 4.18.1-py36hf484d3e_1000



# buggy packages:
#conda install -c bioconda -y qualimap=2.2.2d && (I can instead install 2.2.2, but it fails). I decide to skip this package

# GENERAL NOTES
# conda install -c bioconda entrez-direct=13.3 would be to install esearch and efetch, but they are already installed with VEP
# at the beginning, the installation of the conda packages is not enough to be able to load them from Rscript

# conda install -c r r-stringr=1.4.0 was already installed









