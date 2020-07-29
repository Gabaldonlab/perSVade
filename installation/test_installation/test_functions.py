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
EnvName = EnvDir.split("/")[-1]
AnacondaDir = "/".join(sys.executable.split("/")[0:-4])

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


def test_conda_env_generation(outdir):

    """This function exports the current perSVade_env to a .yml file, and generates a conda file"""

    # define the file that indicates that the enviornment is correct
    correct_env_file = "%s/correct_env.txt"%outdir

    # define a test_env_name
    test_env_name = "%s_test"%EnvName

    if fun.file_is_empty(correct_env_file):

        # remove previous env
        print("removing previous env")
        try: fun.run_cmd("conda remove -y -n %s --all"%test_env_name)
        except: print("%s does not exist"%test_env_name)

        # export file
        print("creating %s yml"%EnvName)
        yml_file = "%s/%s.yml"%(outdir, test_env_name)
        fun.run_cmd("conda env export --no-builds --from-history -n %s --file %s"%(EnvName, yml_file))

        # create environment
        print("re-generating as %s"%test_env_name)
        fun.run_cmd("conda env create --file %s --name %s"%(yml_file, test_env_name))
        
        # test that the activation works
        print("activating %s"%test_env_name)
        cmd_activation = "source %s/etc/profile.d/conda.sh && conda activate %s && python -c 'import sys; sys.path.insert(0, \"%s\"); import sv_functions as fun'"%(AnacondaDir, test_env_name, fun.get_fullpath(scripts_dir))
        fun.run_cmd(cmd_activation)

        # remove file
        print("removing envs")
        fun.remove_file(yml_file)

        # remove env
        fun.run_cmd("conda remove -y -n %s --all"%test_env_name)

        # create file stating that the env is correct
        open(correct_env_file, "w").write("env is correct")

    print("%s can be correctly regenerated"%EnvName)

