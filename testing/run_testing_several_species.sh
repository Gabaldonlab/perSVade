#!/bin/sh

./testing_several_species.py

# this can be run in a node of the cluster with
#sbatch --qos=bsc_ls --time=20:00:00 --job-name=testing_perSVade --cpus-per-task=48 --error=/gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/testing_perSVade_benchmark_stderr.txt --output=/gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/testing_perSVade_benchmark_stdout.txt --workdir=.  --get-user-env run_testing_several_species.sh
