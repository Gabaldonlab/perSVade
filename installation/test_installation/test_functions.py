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
import sv_functions as fun

######################################

def test_get_repeat_maskerDF(test_genome):

    """Tests the generation of repeats"""

    # define the ref genome
    df_repeats, repeat_masker_outfile_default = fun.get_repeat_maskerDF(test_genome, threads=4, replace=False)

    # test
    expected_fields = {'perc_ins', 'perc_del', 'type', 'end_repeat', 'perc_div', 'position_inRepeat_begin', 'repeat', 'left_repeat', 'IDrepeat', 'strand', 'left_positionINrepeat', 'SW_score', 'chromosome', 'begin_repeat', 'position_inRepeat_end'}

    if set(list(df_repeats.keys()))!=expected_fields: raise ValueError("something went wrong with the repeats generation")

    print("repeats were generated correctly")



