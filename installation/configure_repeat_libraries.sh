#!/bin/bash
set -e

# this script configures the repeat libraries inside the Docker image of perSVade. It has to be ran from perSVade_env_RepeatMasker_env

# define the path to the environment
repeatMasker_env_path=/opt/conda/envs/perSVade_env_RepeatMasker_env

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

echo 'RepeatLibraries Properly configured'

