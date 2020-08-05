#!/bin/bash
set -e

# This script sets up the environment to adjust the RepeatMasker libraries

# define the path to the conda environment
env_path=$(which python | sed 's/\/bin\/python//g');
env_name=$(echo $env_path | rev | cut -d '/' -f1 | rev)

# define the path to the installation dir
installation_dir=$(readlink -f $0 | sed 's/\/setup_environment.sh//g')

# go to the path where the RepeatMasker is
cd $env_path/share/RepeatMasker;

# This is just to create the RepeatMasker.lib. Press Ctrl+C when the program asks for any input
echo "generating RepeatMasker.lib"
expected_file=./Libraries/RepeatMasker.lib
t="0"
while [ ! -f $expected_file ]
do
	t=$[$t+10]
	timeout $t ./configure || echo 'exiting the generation on RepeatMasker.lib';

done

# go to the folder where all the pipelines are
cd Libraries;

# generate the database for RepeatPeps
makeblastdb -in RepeatPeps.lib -input_type fasta -dbtype prot;

# make the database for RepeatMasker
makeblastdb -in RepeatMasker.lib -input_type fasta -dbtype nucl;

# go to the installation dir
cd $installation_dir
if [ -d external_software ] 
then 
	rm -r external_software
fi
mkdir external_software
cd external_software

echo 'NOTE: Befeore setting up the installation of further dependencies, you have to make sure that the folder containing the binary of Ninja (https://github.com/TravisWheelerLab/NINJA, preferably release 0.95-cluster_only) can be found in your $PATH. This is an example of how you can do this:'
echo '---'
echo 'cd <path_to_install_Ninja>'
echo 'wget https://github.com/TravisWheelerLab/NINJA/archive/0.95-cluster_only.tar.gz'
echo 'tar -xvf 0.95-cluster_only.tar.gz'
echo 'rm 0.95-cluster_only.tar.gz'
echo 'cd NINJA-0.95-cluster_only/NINJA'
echo 'make'
echo 'export PATH=$PATH:<path_to_install_Ninja>/NINJA-0.95-cluster_only/NINJA'
echo '---'
echo 'You should check that a <path_to_install_Ninja>/NINJA-0.95-cluster_only/NINJA/Ninja binary was created. This is enough for running perSVade'

# copy the Ninja directory under external_software
ninja_binary=$(which Ninja)
ninja_dir=${ninja_binary%/*}
cp -r $ninja_dir ./NINJA
echo 'NINJA is installed successfully in your system'

# download the gztool depending on the architecture
architechture=$(uname -i)

if [ $architechture = "x86_64" ]
then
	echo "installing gztool for architechture $architechture ..."
	wget https://github.com/circulosmeos/gztool/releases/download/v0.11.5/gztool-linux.x86_64
	mv gztool-linux.x86_64 gztool

elif [ $architechture = "aarch64" ]
then
	echo "installing gztool for architechture $architechture"
	wget https://github.com/circulosmeos/gztool/releases/download/v0.11.5/gztool-linux.aarch64
	mv gztool-linux.aarch64 gztool

else
	echo "WARNING: gztool can't be installed in architechture $architechture. perSVade will use more standard and slower tools instead. "
fi

# download extra dependencies (platform-independent)
wget https://github.com/PapenfussLab/gridss/releases/download/v2.9.2/gridss-2.9.2-gridss-jar-with-dependencies.jar
wget https://github.com/PapenfussLab/gridss/releases/download/v2.9.2/gridss.sh
wget https://github.com/PapenfussLab/clove/releases/download/v0.17/clove-0.17-jar-with-dependencies.jar

# give execution permission to all
chmod u+x *;

# create a subenvironment that has bctools=1.10.2
bcftools_env_name="$env_name"_bcftools_1.10.2_env
echo "creating conda env $bcftools_env_name"
conda create -y --name $bcftools_env_name -c bioconda bcftools=1.10.2;

# print that everything went well
echo 'SUCCESS: all perSVade dependencies were properly installed.'
