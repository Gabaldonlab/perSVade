#!/bin/bash

# This recapitulates all the steps to create the env (on conda 4.8.0)

# conda distribution 'Anaconda3-2019.03-Linux-x86_64.sh', followed by 'conda install conda=4.8.0'

# remove previous
#conda activate base
#conda remove --name perSVade_env --all

# create and activate the env
#conda create -n perSVade_env python=3.6 &&
#conda activate perSVade_env &&

# add channels
conda config --add channels conda-forge &&
conda config --add channels biocore &&
conda config --add channels bioconda &&

# install packages that should be loaded 
conda install pandas=0.24.2 &&
conda install biopython=1.73 &&
conda install scipy=1.4.1 &&
conda install scikit-learn=0.21.3 &&
conda install -c conda-forge igraph=0.7.1 && # updates openssl pkgs/main::openssl-1.1.1g-h7b6447c_0 --> conda-forge::openssl-1.1.1g-h516909a_1
conda install -c conda-forge python-igraph=0.7.1.post7 &&
conda install -c conda-forge/label/cf201901 ete3=3.1.1 && # this changes the priority
conda install -c anaconda seaborn=0.9.0 && # updates: ca-certificates    conda-forge/label/cf201901::ca-certif~ --> anaconda::ca-certificates-2020.6.24-0, certifi            conda-forge/label/cf201901::certifi-2~ --> anaconda::certifi-2020.6.20-py36_0

# install packages related to software
conda install -c bioconda repeatmasker=4.0.9_p2 && # downgrades: gmp              6.2.0-he1b5a44_2 --> 6.1.2-hf484d3e_1000
conda install -c bioconda repeatmodeler=2.0.1 &&

