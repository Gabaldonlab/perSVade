#!/bin/bash
set -e

# This script sets up the environment to adjust the RepeatMasker libraries

# define the path to the conda environment
env_path=$(which python | sed 's/\/bin\/python//g');
env_name=$(echo $env_path | rev | cut -d '/' -f1 | rev)
conda_envs_dir=$(echo $env_path | sed "s/\/$env_name//g")
conda_dir=$(echo $env_path | sed "s/\/envs\/$env_name//g")

# define the path to the installation dir
installation_dir=$(readlink -f $0 | sed 's/\/setup_environment.sh//g')


####### INSTALL EXTERNAL SOFTWARE #######

# go to the installation dir
cd $installation_dir
if [ -d external_software ] 
then 
	rm -rf external_software
fi
mkdir external_software
cd external_software

# install the optimised lowess package
git clone https://github.com/histed/lowess-work-carljv
cd lowess-work-carljv
python setup.py install

# back to the external_software dir
cd $installation_dir/external_software

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
wget https://raw.githubusercontent.com/weiyuchung/CONY/master/CONY.R

# give execution permission to all
chmod u+x *;

###################################

######### CREATE THE REPEATMASKER SUBENVIRONMENT #########

# create a subenvitonment to run repeat modeller. This needs a specific GLIBC version
repeatMasker_env_name="$env_name"_RepeatMasker_env
echo "creating conda env $repeatMasker_env_name"
conda create -y --name $repeatMasker_env_name -c bioconda repeatmasker=4.0.9_p2 repeatmodeler=2.0.1

# activate this environment
source $conda_dir/etc/profile.d/conda.sh && conda activate $repeatMasker_env_name

# define the path to the environment
repeatMasker_env_path=$conda_envs_dir/$repeatMasker_env_name

# go to the path where the RepeatMasker is
cd $repeatMasker_env_path/share/RepeatMasker;

# This is just to create the RepeatMasker.lib. Press Ctrl+C when the program asks for any input
echo "configuring RepeatModeler"
expected_file=./Libraries/RepeatMasker.lib
t="0"
while [ ! -s $expected_file ]
do

	# remove files so that in a next run the repeats library will be re-generated
	rm ./Libraries/RepeatMasker.lib || echo 'already removed' 
	rm ./Libraries/RepeatMaskerLib.embl || echo 'already removed' 

	# get the libraries
	t=$[$t+20]
	timeout $t ./configure || echo 'exiting the generation on RepeatMasker.lib';

done

# go to the folder where all the pipelines are
cd Libraries;

# generate the database for RepeatPeps
makeblastdb -in RepeatPeps.lib -input_type fasta -dbtype prot;

# make the database for RepeatMasker
makeblastdb -in RepeatMasker.lib -input_type fasta -dbtype nucl;


############################################################

######## GET THE NINJA BINARY ########

# go to the installation dir of Ninja
cd $installation_dir/Ninja_data

# define a function that takes the Ninja and the type_installation and the ninja binary
function test_Ninja_installation {

	infile=$installation_dir/test_installation/testing_inputs/MERS_CoV_genome.fasta
	outfile=$installation_dir/Ninja_data/ninja_test.out
               
	$1 --in $infile --out $outfile && echo "Ninja installation worked through $2" && ninja_binary=$1

    }

# define a function that teststs that g++ is properly installed
function test_gxx_installation { 
	
	g++ --help > /dev/null 2>&1 
	
	}

# first try to run as if it was already installed
previously_installed_Ninja=$(which Ninja) || echo 'Ninja was nor previously installed'
test_Ninja_installation $previously_installed_Ninja "previously_installed"  || {

	echo "Ninja was not previously installed. Installing through several options"

	# install g++ with conda, if not installed
	test_gxx_installation || conda install -y gxx_linux-64 # this will install gcc and gxx
	test_gxx_installation || ln -s $repeatMasker_env_path/bin/x86_64-conda_cos6-linux-gnu-g++ $repeatMasker_env_path/bin/g++
	test_gxx_installation || chmod 755 $repeatMasker_env_path/bin/g++

	# test that installation worked
	test_gxx_installation && echo 'g++ is properly installed'

	# try to install NINJA 0.95
	rm  -rf 0.95-cluster_only* || echo 'It is the first time this folder is attempted to be accessed'
	wget https://github.com/TravisWheelerLab/NINJA/archive/0.95-cluster_only.tar.gz
	tar -xvf 0.95-cluster_only.tar.gz && rm 0.95-cluster_only.tar.gz 
	cd NINJA-0.95-cluster_only/NINJA

	make all || echo 'NINJA-0.95-cluster_only could not be compiled'

	# check 
	test_Ninja_installation $installation_dir/Ninja_data/NINJA-0.95-cluster_only/NINJA/Ninja "NINJA-0.95-cluster_only"  || {

		echo 'installation through NINJA-0.95-cluster_only did not work. Trying with the repository'
		cd $installation_dir/Ninja_data

		# decompress
		rm -rf NINJA_repo_05112020 || echo 'Ninja repo was never accessed'
		tar -xvf NINJA_repo_05112020.tar.gz && mv NINJA NINJA_repo_05112020

		# go inside the folder and make the binaries
		cd NINJA_repo_05112020
		make build || echo 'The Ninja binary could not be built from the repository'

		# test
		test_Ninja_installation  $installation_dir/Ninja_data/NINJA_repo_05112020/NINJA/Ninja "NINJA_repo_05112020" || {


			echo 'ERROR: Ninja could not be installed. You should do it on your own for perSVade to run. Check the options in their github repo https://github.com/TravisWheelerLab/NINJA. Another option is to use one of the binaries provided under installation/Ninja_binaries, but this may not work for some systems. Once you have NINJA installed, you should add folder that contains it to $PATH, and re-run this script to finish the setup of perSVade.'
			exit
		}
	}
}
echo "The used Ninja binary is $ninja_binary"

# copy the Ninja binary under the conda environment 
new_ninja_binary=$repeatMasker_env_path/bin/Ninja
if [ $ninja_binary != $new_ninja_binary ]
then

	if [ -f $new_ninja_binary ] 
	then 
		rm $new_ninja_binary
	fi
	cp $ninja_binary $new_ninja_binary

fi
echo 'NINJA is installed successfully in your system'

##########################################################


###### INSTALL OTHER SOFTWARE ######

# activate the normal env
source $conda_dir/etc/profile.d/conda.sh && conda activate $env_name

# create a subenvironment that has bctools=1.10.2
bcftools_env_name="$env_name"_bcftools_1.10.2_env
echo "creating conda env $bcftools_env_name"
conda create -y --name $bcftools_env_name -c bioconda bcftools=1.10.2;

# create a submenvironment that has ete3=3.0.0
ete3_env_name="$env_name"_ete3_3.0.0_env
echo "creating conda env $ete3_env_name"
conda config --add channels conda-forge
rm -rf $conda_dir/envs/$ete3_env_name || echo 'ete3 was not created'
conda create -y --name $ete3_env_name python=3.6 # I install python because it is missing
conda install -n $ete3_env_name -c conda-forge -y ete3=3.0.0;

# fix the ete3 script
ete3_script="$conda_envs_dir"/"$ete3_env_name"/lib/python3.6/site-packages/ete3/ncbi_taxonomy/ncbiquery.py
python $installation_dir/fix_ete3_script.py $ete3_script;

# create a subenvironment with all the R dependencies
R_env_name="$env_name"_R_env
echo "creating conda env $R_env_name"
conda create -y --name $R_env_name;

conda install -n $R_env_name -c conda-forge -y r-base=4.0.2;
conda install -n $R_env_name -c conda-forge -y r-argparser=0.6;
conda install -n $R_env_name -c bioconda -y bioconductor-rsvsim=1.28;
conda install -n $R_env_name -c conda-forge -y r-emdbook=1.3.12;
conda install -n $R_env_name -c bioconda -y bioconductor-rtracklayer=1.48.0;
conda install -n $R_env_name -c conda-forge -y r-r.utils=2.9.2;
conda install -n $R_env_name -c bioconda -y bioconductor-structuralvariantannotation=1.4.0;

# create a subenvironment for running CONY. The last commit of CONY was the 12/06/2017. This is 3 years and 4 months ago., so that I have to install software from a similar age

# packages from bioconda:

# bioconductor-iranges=2.8.2  -> 3 years and 5 months -> R 3.3.2 / R 3.3.1
# bioconductor-exomecopy=1.22.0 -> 3 years and 4 months -> R 3.4.1 / R 3.3.2
# r-snow=0.4 -> 3 years and 5 months -> R 3.3.2 

# conda-forge
# r-argparser=0.4 -> 3 years and 2 months -> R 3.3.2

# 

CONY_env_name="$env_name"_CONY_env
echo "creating conda env $CONY_env_name"
conda create -y --name $CONY_env_name -c bioconda bioconductor-iranges=2.8.2 bioconductor-exomecopy=1.22.0 r-snow=0.4; # this should install r-base=3.3.2, which fits all the dependencies
conda install -y -n $CONY_env_name -c conda-forge r-argparser=0.4;

# create a subenvironment to run gridss
gridss_env_name="$env_name"_gridss_env
echo "creating conda env $gridss_env_name"
conda create -y --name $gridss_env_name -c conda-forge r-base=4.0.2
conda install -n $gridss_env_name -c bioconda -y samtools=1.10
conda install -n $gridss_env_name -c bioconda -y bwa=0.7.17
conda install -n $gridss_env_name -c anaconda -y openjdk=8.0.152

# create a subenvironment to run picard
picard_env_name="$env_name"_picard_env
echo "creating conda env $picard_env_name"
conda create -y --name $picard_env_name -c conda-forge r-base=4.0.2
conda install -n $picard_env_name -c bioconda -y picard=2.18.26

####################################


# print that everything went well
echo 'SUCCESS: all perSVade dependencies were properly installed.'
