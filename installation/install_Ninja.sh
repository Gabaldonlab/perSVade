#!/bin/bash
set -e

# This script installs Ninja

# define the path to the conda environment
env_path=$(which python | sed 's/\/bin\/python//g');
env_name=$(echo $env_path | rev | cut -d '/' -f1 | rev)
conda_envs_dir=$(echo $env_path | sed "s/\/$env_name//g")
conda_dir=$(echo $env_path | sed "s/\/envs\/$env_name//g")

# define the path to the installation dir
installation_dir=$(readlink -f $0 | sed 's/\/install_Ninja.sh//g')

# define the env path of repeat masker
repeatMasker_env_name=$env_name'_RepeatMasker_env'
repeatMasker_env_path=$conda_envs_dir/$repeatMasker_env_name

# activate this environment
source $conda_dir/etc/profile.d/conda.sh && conda activate $repeatMasker_env_name

# go to the installation dir of Ninja
cd $installation_dir/Ninja_data

# define a function that takes the Ninja and the type_installation and the ninja binary and tests if it works
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
	test_gxx_installation || mamba install -y gxx_linux-64=12.1.0 # this will install gcc and gxx
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


			echo 'ERROR: Ninja could not be installed. This means that the repeat masker features will not worl. You should do it on your own for perSVade to run. Check the options in their github repo https://github.com/TravisWheelerLab/NINJA. Another option is to use one of the binaries provided under installation/Ninja_binaries, but this may not work for some systems. Once you have Ninja installed, you should add folder that contains it to $PATH, and re-run this script to finish the setup of perSVade.'
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
echo 'SUCCESS: NINJA is installed successfully in your system'
