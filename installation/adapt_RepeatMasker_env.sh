#!/bin/bash
set -e

# This script sets up the '$PERSVADE_ENV_NAME'_RepeatMasker_env environment to adjust the RepeatMasker libraries

# define the path to the conda environment
env_path=$(which python | sed 's/\/bin\/python//g');
env_name=$(echo $env_path | rev | cut -d '/' -f1 | rev)
conda_envs_dir=$(echo $env_path | sed "s/\/$env_name//g")
conda_dir=$(echo $env_path | sed "s/\/envs\/$env_name//g")

# activate this environment
repeatMasker_env_name=$env_name'_RepeatMasker_env'
source $conda_dir/etc/profile.d/conda.sh && conda activate $repeatMasker_env_name

# define the path to the environment
repeatMasker_env_path=$conda_envs_dir/$repeatMasker_env_name

# go to the path where the RepeatMasker is
cd $repeatMasker_env_path/share/RepeatMasker;

# Create the RepeatMasker.lib
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

echo 'SUCCESS: The RepeatMasker env was properly formated.'