#!/usr/bin/env python

######### define environment ##########

# module imports
import sys
import os

# get the cwd were all the scripts are 
test_dir = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, test_dir)
scripts_dir = "%s/../../scripts"%test_dir; sys.path.insert(0, scripts_dir)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# import functions
import test_functions as test_fun

# define the testing inputs (based on the first 2 genes of two chromosomes of Cglabrata)
testing_inputs_dir = "%s/testing_inputs"%test_dir
test_ref_genome = "%s/reduced_genome.fasta"%testing_inputs_dir
test_mutated_genome = "%s/reduced_genome_mutated.fasta"%testing_inputs_dir
test_gff = "%s/reduced_annotation.gff"%testing_inputs_dir

# define the testing inuts dir
testing_outputs_dir = "%s/testing_outputs"%test_dir
test_output_perSVade = "%s/perSVade_output"%testing_outputs_dir

########################################

# test if you can import python packages
import sv_functions as fun
print("loading python packages worked successfully")

# redefine the reference genome location
ref_genome = "%s/reduced_genome.fasta"%testing_outputs_dir
if fun.file_is_empty(ref_genome): fun.run_cmd("cp %s %s"%(test_ref_genome, ref_genome))

# test repeat masker obtention
test_fun.test_get_repeat_maskerDF(ref_genome)

