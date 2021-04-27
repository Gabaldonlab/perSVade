#!/bin/sh

./testing_humanGoldenSet.py

# this can be run in a node of the cluster with
#sbatch --qos=bsc_ls --time=48:00:00 --job-name=testing_perSVade --cpus-per-task=48 --error=/gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/testing_perSVade_benchmark_humanGoldenSet_stderr.txt --output=/gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/testing_perSVade_benchmark_humanGoldenSet_stdout.txt --workdir=.  --get-user-env run_testing_humanGoldenSet.sh

#bsub -e /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/testing_perSVade_benchmark_humanGoldenSet_stderr.txt -o /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/testing_perSVade_benchmark_humanGoldenSet_stdout.txt -cwd . -W 48:00 -q bsc_ls -n 16 -M 4000 run_testing_humanGoldenSet.sh

