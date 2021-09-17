FROM continuumio/miniconda3:4.8.2

######### SETUP ENVIRONMENT ###########

# define the working directory (creates the folder /perSVade inside the virtual machine)
WORKDIR /perSVade

# copy all the necessary files into /perSVade. This is everything from the github repository.
COPY . .

# remove the scripts and test installation files, in order to be able to mount the folder from outside (only for debugging the container)
#RUN rm -r /perSVade/scripts/*
#RUN rm -r /perSVade/installation/test_installation/*

# give permissions to necessary files
RUN chmod 755 /perSVade/scripts/*
RUN chmod 755 /perSVade/installation/*
RUN chmod -R 777 /perSVade/installation/test_installation/

# install mamba
RUN conda install -y -c conda-forge mamba=0.15.3

# create the perSVade environment from the .yml file. mamba is way faster than conda
RUN mamba env create --file ./installation/perSVade_env.yml --name perSVade_env

#######################################

########### ADD EXTRA DEPENDENCIES #########

# define variables
ARG env_path=/opt/conda/envs/perSVade_env/bin
ARG env_name=perSVade_env
ARG conda_envs_dir=/opt/conda/envs
ARG conda_dir=/opt/conda
ARG installation_dir=/perSVade/installation
ARG external_software_dir=/perSVade/installation/external_software
ARG lowess_dir=$external_software_dir/lowess-work-carljv
ARG repeatMasker_env_name="$env_name"_RepeatMasker_env
ARG repeatMasker_env_path=$conda_envs_dir/$repeatMasker_env_name
ARG ninja_binary=$installation_dir/NINJA_repo_05112020/NINJA/Ninja
ARG new_ninja_binary=$repeatMasker_env_path/bin/Ninja
ARG bcftools_env_name="$env_name"_bcftools_1.10.2_env
ARG ete3_env_name="$env_name"_ete3_3.0.0_env
ARG R_env_name="$env_name"_R_env
ARG CONY_env_name="$env_name"_CONY_env
ARG AneuFinder_env_name="$env_name"_AneuFinder_env
ARG HMMcopy_env_name="$env_name"_HMMcopy_env
ARG gridss_env_name="$env_name"_gridss_env
ARG picard_env_name="$env_name"_picard_env

# Make below RUN commands use the perSVade_env:
SHELL ["conda", "run", "-n", "perSVade_env", "/bin/bash", "-e", "-c"]

# make dirs
RUN mkdir $external_software_dir
RUN mkdir $lowess_dir

# install the lowess package (and necessary conda dependencies) 
RUN mamba install -y -n perSVade_env -c conda-forge gcc_linux-64=11.1.0
RUN mamba install -y -n perSVade_env -c anaconda lxml=4.5.1
RUN git clone https://github.com/histed/lowess-work-carljv $lowess_dir/
RUN cd $lowess_dir && python setup.py install

# test that cylowess works
RUN python -c 'import cylowess'

# download gztool
RUN wget https://github.com/circulosmeos/gztool/releases/download/v0.11.5/gztool-linux.x86_64
RUN mv gztool-linux.x86_64 $external_software_dir/gztool

# download other platform-independent dependencies
RUN cd $external_software_dir && wget https://github.com/PapenfussLab/gridss/releases/download/v2.9.2/gridss-2.9.2-gridss-jar-with-dependencies.jar
RUN cd $external_software_dir && wget https://github.com/PapenfussLab/gridss/releases/download/v2.9.2/gridss.sh
RUN cd $external_software_dir && wget https://github.com/PapenfussLab/clove/releases/download/v0.17/clove-0.17-jar-with-dependencies.jar
RUN cd $external_software_dir && wget https://raw.githubusercontent.com/weiyuchung/CONY/master/CONY.R

# give permissions to all
RUN cd $external_software_dir && chmod 755 *

# add conda channels
RUN conda config --add channels conda-forge
RUN conda config --add channels anaconda
RUN conda config --add channels bioconda

# create an environment to run repeat masker and modeller
RUN mamba create -y --name $repeatMasker_env_name -c bioconda repeatmasker=4.0.9_p2 repeatmodeler=2.0.1 
RUN mamba install --name $repeatMasker_env_name -c anaconda make=4.2.1

# run the commands below on $repeatMasker_env_name
SHELL ["conda", "run", "-n", "perSVade_env_RepeatMasker_env", "/bin/bash", "-e", "-c"]

# configure the Repeat Libraries
RUN bash $installation_dir/configure_repeat_libraries.sh

# install Ninja and important dependencies
RUN mamba install -y -n $repeatMasker_env_name gxx_linux-64
RUN ln -s $repeatMasker_env_path/bin/x86_64-conda_cos6-linux-gnu-g++ $repeatMasker_env_path/bin/g++
RUN chmod 755 $repeatMasker_env_path/bin/g++
RUN tar -xvf $installation_dir/Ninja_data/NINJA_repo_05112020.tar.gz && mv NINJA $installation_dir/NINJA_repo_05112020
RUN cd $installation_dir/NINJA_repo_05112020 && make build
RUN cp $ninja_binary $new_ninja_binary

# test that NINJA was properly installed (skip in debug mode)
RUN $new_ninja_binary --in $installation_dir/test_installation/testing_inputs/MERS_CoV_genome.fasta --out $installation_dir/Ninja_data/ninja_test.out && echo "Ninja installation worked!!!"

# Reactivate the perSVade_env:
SHELL ["conda", "run", "-n", "perSVade_env", "/bin/bash", "-e", "-c"]

# create bcftools environment
RUN mamba create -y --name $bcftools_env_name -c bioconda bcftools=1.10.2

# create a submenvironment that has ete3=3.0.0
RUN conda config --add channels conda-forge
RUN mamba create -y --name $ete3_env_name python=3.6
RUN mamba install -n $ete3_env_name -c etetoolkit -y ete3=3.1.2
RUN source $conda_dir/etc/profile.d/conda.sh && conda activate $ete3_env_name && python -c 'import ete3'

# create a subenvironment with all the R dependencies
RUN mamba create -y --name $R_env_name
RUN mamba install -n $R_env_name -c conda-forge -y r-base=4.0.2
RUN mamba install -n $R_env_name -c conda-forge -y r-argparser=0.6
RUN mamba install -n $R_env_name -c bioconda -y bioconductor-rsvsim=1.28
RUN mamba install -n $R_env_name -c conda-forge -y r-emdbook=1.3.12
RUN mamba install -n $R_env_name -c bioconda -y bioconductor-rtracklayer=1.48.0
RUN mamba install -n $R_env_name -c conda-forge -y r-r.utils=2.9.2
RUN mamba install -n $R_env_name -c bioconda -y bioconductor-structuralvariantannotation=1.4.0

# create environment to run CONY
RUN mamba create -y --name $CONY_env_name -c bioconda bioconductor-iranges=2.8.2 bioconductor-exomecopy=1.22.0 r-snow=0.4
RUN mamba install -y -n $CONY_env_name -c conda-forge r-argparser=0.4

# create an environment to run AneuFinder
RUN mamba create -y --name $AneuFinder_env_name -c bioconda bioconductor-aneufinder=1.18.0
RUN mamba install -y -n $AneuFinder_env_name -c conda-forge r-argparser=0.6

# create a subenvironment to run HMMCOPY
RUN mamba create -y --name $HMMcopy_env_name
RUN conda config --add channels conda-forge
RUN mamba install -y -n $HMMcopy_env_name -c bioconda bioconductor-hmmcopy=1.32.0
RUN mamba install -y -n $HMMcopy_env_name -c conda-forge r-argparser=0.6
RUN mamba install -y -n $HMMcopy_env_name -c conda-forge r-kernsmooth=2.23

# create a subenvironment to run gridss
RUN mamba create -y --name $gridss_env_name -c conda-forge r-base=4.0.2
RUN mamba install -n $gridss_env_name -c bioconda -y samtools=1.10
RUN mamba install -n $gridss_env_name -c bioconda -y bwa=0.7.17
RUN mamba install -n $gridss_env_name -c anaconda -y openjdk=8.0.152

# create a subenvironment to run picard
RUN mamba create -y --name $picard_env_name -c conda-forge r-base=4.0.2
RUN mamba install -n $picard_env_name -c bioconda -y picard=2.18.26

############################################

# activate perSVade's environment on running. "conda run --no-capture-output" would be to skip the generation of STDs
ENTRYPOINT ["conda", "run", "--no-capture-output", "--live-stream", "-n", "perSVade_env"] 

####### GENERAL COMMENTS ########

# this image can be built with 'docker build -t mikischikora/persvade:v1 .'

# This image contains the perSVade environment, and it can be used through 'docker run -i mikischikora/persvade:v1 <commands> '

# we created a .dockerignore file with some folders not to consider

# We tested perSVade on version 4.8.0, but this was not available on dockerhub. 

# we can pass arguments like -v <host_path>:<docker_path>

# 'FROM continuumio/miniconda3:4.8.2' would be to have the directly installed conda, closest to reality.

# run a terminal with 'docker run -it mikischikora/persvade:v1 bash'

# publish the image with 'docker push mikischikora/persvade:v1'

# you can save this image with 'docker save mikischikora/persvade:v1 | gzip > ./perSVade_docker_image.tar.gz'

# and load with 'docker load -i ./perSVade_docker_image.tar'

# testing docker with linked data in /home/mschikora (debug mode):
# docker run --memory "20g" -v /home/mschikora/samba/scripts/perSVade/perSVade_repository/scripts:/perSVade/scripts -v /home/mschikora/samba/scripts/perSVade/perSVade_repository/installation/test_installation:/perSVade/installation/test_installation mikischikora/persvade:v1 python -u ./installation/test_installation/test_installation.py

# testing installation with the actual scripts inside the container
# docker run --memory "4g" -v /home/mschikora/samba/scripts/perSVade/perSVade_repository/installation/test_installation/testing_outputs:/perSVade/installation/test_installation/testing_outputs mikischikora/persvade:v1 python -u ./installation/test_installation/test_installation.py

#################################
