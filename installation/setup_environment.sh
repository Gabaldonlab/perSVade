#!/bin/bash

# This script sets up the environment to adjust the RepeatMasker libraries

# define the oath to conda
env_path=$(which python | sed 's/\/bin\/python//g');

# go to the path where the RepeatMasker is
cd $env_path/share/RepeatMasker;

# This is just to create the RepeatMasker.lib. Press Ctrl+C when the program asks for any input
echo "generating RepeatMasker.lib"
expected_file=./Libraries/RepeatMasker.lib
t="0"
while [ ! -f $expected_file ]
do
	t=$[$t+10]
	timeout $t ./configure; 

done

# go to the folder where all the pipelines are
cd Libraries;

# generate the database for RepeatPeps
makeblastdb -in RepeatPeps.lib -input_type fasta -dbtype prot;

# make the database for RepeatMasker
makeblastdb -in RepeatMasker.lib -input_type fasta -dbtype nucl 
