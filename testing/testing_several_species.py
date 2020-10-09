#!/usr/bin/env python

# This is a script that runs the testing of perSVade on several species

##### DEFINE ENVIRONMENT #######

# module imports
import os
import sys
import pandas as pd

# define the parent dir of the cluster or not
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir):
    run_in_cluster = False    
    threads = 4
else:
    run_in_cluster = True    
    ParentDir = "/gpfs/projects/bsc40/mschikora"
    threads = 48


# define the dir where all perSVade code is
perSVade_dir = "%s/scripts/perSVade/perSVade_repository/scripts"%ParentDir
sys.path.insert(0, perSVade_dir)

# import functions
import sv_functions as fun

# define paths
perSVade_py = "%s/perSVade.py"%perSVade_dir

# define dirs
outdir_testing = "%s/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies"%ParentDir; fun.make_folder(outdir_testing)
CurDir = "%s/scripts/perSVade/perSVade_repository/testing"%ParentDir; fun.make_folder(outdir_testing)
outdir_genomes_and_annotations = "%s/scripts/perSVade/perSVade_repository/testing/genomes_and_annotations"%ParentDir

################################


"""
This is how the genomes were obtained:

C. glabrata: reference genome from CGD: the latest version by 12/03/2019, which is v_s02-m07-r35 

C. albicans: 

    ref genome CGD: http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_version_A22-s07-m01-r110_chromosomes.fasta.gz

    gff from CGD: http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_version_A22-s07-m01-r110_features.gff

    From here I keep 'haplotype A' for both files

C. neoformans: ref genome from GenBank GCA_000149245.3

A. fumigatus: 

    gDNA from reference NCBI:

    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.gff.gz
    
    mtDNA from https://www.ncbi.nlm.nih.gov/nuccore/CM016889.1

A. thaliana: ref genome from GenBank GCA_000001735.2

D. melanogaster: ref genome from GenBank GCA_000001215.4

D. rerio: ref genome from GenBank removing the alternate haplotypes. (this is GCA_000002035.4)

H. sapiens: ref genome from GenBank removing the alternate haplotypes. (this is GCA_000001405.28)

For C. glabrata I got the nanopore reads from ~/../mmarcet/nanopore/GABALDON02/assembly_files/BG2/nanopore.reads.pass.fastq.gz and the short reads from Ewa's experiment in RUN4_BG2_SRA_WT

"""

# define the table for C. glabrata
close_shortReads_table_Cglabrata = "%s/scripts/perSVade/perSVade_repository/testing/Cglabrata_table_short_reads.tab"%ParentDir
goldenSet_dir_Cglabrata = "%s/scripts/perSVade/perSVade_repository/testing/Cglabrata_goldenSetReads_BG2"%ParentDir

# define important info about each species: taxID, spName, ploidy

"""
species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138"),
                ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314"),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1"),
                ("746128", "Aspergillus_fumigatus", 1, "CM016889.1"),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1"),
                ("7227", "Drosophila_melanogaster", 2, "KJ947872.2")]
                #("7955", "Danio_rerio", 2, "NC_002333.2")]
                #("9606", "Homo_sapiens", 2, "NC_012920.1")]
"""

"""
species_Info = [("7227", "Drosophila_melanogaster", 2, "KJ947872.2")]
                #("7955", "Danio_rerio", 2, "NC_002333.2")]
                #("9606", "Homo_sapiens", 2, "NC_012920.1")]
"""

"""
species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138"),
                ("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314"),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1"),
                ("746128", "Aspergillus_fumigatus", 1, "CM016889.1")]
"""
#species_Info = [("7227", "Drosophila_melanogaster", 2, "KJ947872.2")]
#species_Info  = [("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1")]
#species_Info = [("5476", "Candida_albicans", 2, "Ca22chrM_C_albicans_SC5314")]
#species_Info = [("5478", "Candida_glabrata", 1, "mito_C_glabrata_CBS138")]
#species_Info = [("5207", "Cryptococcus_neoformans", 1, "CP003834.1")]

"""
species_Info = [("746128", "Aspergillus_fumigatus", 1, "CM016889.1"),
                ("5207", "Cryptococcus_neoformans", 1, "CP003834.1")]
"""
"""
species_Info = [("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30),
                ("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30)]
"""
species_Info = [("3702", "Arabidopsis_thaliana", 2, "BK010421.1,AP000423.1", 30)]
#species_Info = [("7227", "Drosophila_melanogaster", 2, "KJ947872.2", 30)]


taxIDs_with_noON_overalpping = {"5476", "746128"}

# define the type of run
running_type = "normalRun" # can be 'normalRun' or 'goldenSet'
StopAfterPrefecth_of_reads = False
compute_timimng = False

# initialize the df with the timing information
filename_timing_df = "%s/calculating_resources.tab"%CurDir
header_fields = ["species", "sampleID", "type_run", "threads", "nvars", "nsimulations", "range_filtering_benchmark", "simulation_ploidies",  "run_time", "exit_status", "finishing_greasy_time", "overall_runID", "std_file", "last_error_line", "job_cmd"]

#  generate the file if not already done
if fun.file_is_empty(filename_timing_df): 
    open(filename_timing_df, "w").write("\t".join(header_fields) + "\n")

# go through each species
for taxID, spName, ploidy, mitochondrial_chromosome, max_coverage_sra_reads in species_Info:
    print(taxID, spName)

    #if spName=="Candida_glabrata": continue # debug

    # define  the genome and annotations
    genome = "%s/%s.fasta"%(outdir_genomes_and_annotations, spName)
    gff = "%s/%s.gff"%(outdir_genomes_and_annotations, spName)

    # create an outdir
    outdir_perSVade = "%s/%s_%s"%(outdir_testing, taxID, spName); fun.make_folder(outdir_perSVade)

    if running_type=="normalRun":

        ###### compute the timing of previous runs ######

        # define files
        greasy_log = "%s/testing_Accuracy/STDfiles/testAccuracy_greasy.log"%outdir_perSVade
        jobs_file = "%s/testing_Accuracy/jobs.testingRealDataAccuracy"%outdir_perSVade

        # check that both of these files exists to continue
        if all([not fun.file_is_empty(x) for x in [greasy_log, jobs_file]]) and compute_timimng is True:
            print("timing")

            # define the overall_runID
            overall_runID = "run3"
            needtochangetheID

            # define the expected jobIDs
            expected_jobIDs = set(range(1,28))

            # define the finishing greasy time
            finishing_greasy_time_lines = ["_".join(l.split("]")[0].split("[")[1].split()) for l in open(greasy_log, "r").readlines() if "Finished greasing" in l]

            if len(finishing_greasy_time_lines)!=1: 
                print("greasy did not finish due to unexpected errors. skipping")
                continue

            finishing_greasy_time = finishing_greasy_time_lines[0]
            if len(finishing_greasy_time)!=19: raise ValueError("the greasy log is not correct")

            # if the combination of species and finishing_greasy_time is already in the df, skip. It means that it is an already included measurement
            previous_df = pd.read_csv(filename_timing_df, sep="\t")
            previous_species_finishing_time_combinations = set(previous_df.species + "_" + previous_df.finishing_greasy_time)
            if "%s_%s"%(spName, finishing_greasy_time) in previous_species_finishing_time_combinations:
                print("already completed species. skipping")
                continue
  
            # map each jobID to an exit status
            jobID_to_exit_status = {int(l.split("located in line ")[1].split()[0]) : l.split()[9] for l in open(greasy_log, "r").readlines() if "Elapsed:" in l}

            # check that all are failed or completed
            if len(set(jobID_to_exit_status.values()).difference({"failed", "completed"}))>0: raise ValueError("All the exist status should be failed or completed")

            # map each jobID to the elapsed time
            jobID_to_elapsed_time = {int(l.split("located in line ")[1].split()[0]) : l.strip().split()[-1] for l in open(greasy_log, "r").readlines() if "Elapsed:" in l}

            # add the remaining jobIDs
            for remaining_jobID in expected_jobIDs.difference(set(jobID_to_elapsed_time)): 

                jobID_to_elapsed_time[remaining_jobID] = "00:00:00"
                jobID_to_exit_status[remaining_jobID] = "pending"


            # go through each job
            jobID_to_metadata = {}
            for Ijob, job_cmd in enumerate(open(jobs_file, "r").readlines()):

                # define the things derived from the outdir
                outdir_job = job_cmd.split("--outdir ")[1].split()[0]
                sampleID = outdir_job.split("/")[-1]
                type_run = outdir_job.split("/")[-2]

                # other things
                threads = int(job_cmd.split("--threads ")[1].split()[0])
                nvars = int(job_cmd.split("--nvars ")[1].split()[0])
                nsimulations = int(job_cmd.split("--nsimulations ")[1].split()[0])
                simulation_ploidies = job_cmd.split("--simulation_ploidies ")[1].split()[0]
                range_filtering_benchmark = job_cmd.split("--range_filtering_benchmark ")[1].split()[0]

                # get the log
                std_file_original = job_cmd.split()[-2]
                std_file = "%s/std_files_testingAccuracy/%s_%s_%s_%ithreads_%ivars_%isims_range:%s_ploidies:%s_greasyFinish:%s.std"%(CurDir, spName, sampleID, type_run, threads, nvars, nsimulations, range_filtering_benchmark, simulation_ploidies, finishing_greasy_time)
                fun.run_cmd("cp %s %s"%(std_file_original, std_file))

                # get the last line with an error
                lines_with_error = [l.strip() for l in open(std_file_original, "r").readlines() if "ERROR" in l.upper()]
                if len(lines_with_error)==0: last_error_line = "no_error"
                else: last_error_line = lines_with_error[-1]

                # get into dict
                jobID_to_metadata[Ijob+1] = {'species': spName, 'sampleID': sampleID, 'type_run': type_run, 'threads': threads, 'nvars': nvars, 'nsimulations': nsimulations, 'simulation_ploidies': simulation_ploidies, 'run_time': jobID_to_elapsed_time[Ijob+1], 'range_filtering_benchmark':range_filtering_benchmark, 'exit_status': jobID_to_exit_status[Ijob+1], "finishing_greasy_time":finishing_greasy_time, "overall_runID":overall_runID, "std_file":std_file, "last_error_line":last_error_line, "job_cmd":job_cmd.split(" >")[0]}

                # remove the 'failed'
                #if jobID_to_exit_status[Ijob+1]=="failed": fun.delete_folder(outdir_job)

            # deifine as df
            df_timimg = pd.DataFrame(jobID_to_metadata).transpose()[header_fields]

            # add the previous df to drop duplicates
            df_timimg = df_timimg.append(pd.read_csv(filename_timing_df, sep="\t")).drop_duplicates()

            # append
            df_timimg.to_csv(filename_timing_df, sep="\t", header=True, index=False)

        # skip the running of the cmds
        if compute_timimng is True: continue # debug

        ################################################

        ####### delete the folders that did not complete in any of the previous runs #######
        """

        for typeSim in ["fast", "uniform", "realSVs"]:
            outdir_testAccuracy = "%s/testing_Accuracy/%s"%(outdir_perSVade, typeSim)
            if not os.path.isdir(outdir_testAccuracy): continue

            for f in os.listdir(outdir_testAccuracy): 
                outdir_f = "%s/%s"%(outdir_testAccuracy, f)
                if fun.file_is_empty("%s/perSVade_finished_file.txt"%(outdir_f)):

                    print("deleting %s"%outdir_f)
                    fun.delete_folder(outdir_f)

        """

        ####################################################################################

        # define the table with short reads
        if spName=="Candida_glabrata": close_shortReads_table = close_shortReads_table_Cglabrata
        else: close_shortReads_table = "auto"

        # get the reads from SRA. 3 samples, 3 runs per sample. Process with the. --verbose
        cmd = "%s --ref %s --threads %i -o %s --close_shortReads_table %s --target_taxID %s --n_close_samples 3 --nruns_per_sample 3 -f1 skip -f2 skip --mitochondrial_chromosome %s --gff %s --testAccuracy --skip_SVcalling --verbose --skip_cleaning_simulations_files_and_parameters --StopAfter_testAccuracy_perSVadeRunning --max_coverage_sra_reads %i"%(perSVade_py, genome, threads, outdir_perSVade, close_shortReads_table, taxID, mitochondrial_chromosome, gff, max_coverage_sra_reads)
        # --StopAfter_testAccuracy_perSVadeRunning --slurm_constraint, --StopAfter_obtentionOFcloseSVs

    elif running_type=="goldenSet":

        # define the goldenSet_dir
        if spName=="Candida_glabrata": goldenSet_dir = goldenSet_dir_Cglabrata
        else: goldenSet_dir = "auto"

        # get the golden set running 
        if taxID in taxIDs_with_noON_overalpping: continue
        cmd = "%s --ref %s --threads %i -o %s --target_taxID %s --n_close_samples 3 --nruns_per_sample 3 -f1 skip -f2 skip --mitochondrial_chromosome %s --gff %s --goldenSet_dir %s --skip_SVcalling --verbose"%(perSVade_py, genome, threads, outdir_perSVade, taxID, mitochondrial_chromosome, gff, goldenSet_dir)

    # add options depending on the machine
    if run_in_cluster is True: cmd += " --job_array_mode greasy --queue_jobs bsc_ls --max_ncores_queue 48 --time_read_obtention 48:00:00 --time_perSVade_running 48:00:00"

    else: cmd += " --job_array_mode local"

    if StopAfterPrefecth_of_reads is True: cmd += " --StopAfterPrefecth_of_reads"

    fun.run_cmd(cmd)

    #if taxID=="5476": adkjhdakg # stop after C. albicans



# an example CMD to debug de generation of merged vcfs
"""

cd ~/samba/CandidaMine_data_generation/v1/data/Candida_albicans_5476/varCall_output/SRR6669901/

~/samba/scripts/perSVade/perSVade_repository/scripts/varcall_cnv_pipeline.py -r ~/samba/CandidaMine_data_generation/v1/data/Candida_albicans_5476/genome.fasta -thr 4 -o smallVars_CNV_output -p 2 -sbam aligned_reads.bam.sorted -c 12 -mchr Ca22chrM_C_albicans_SC5314 -mcode 4 -gcode 12 -gff ~/samba/CandidaMine_data_generation/v1/data/Candida_albicans_5476/annotations.gff --get_merged_vcf


"""


# an example of running the pipeline for adding the repeats

"""

python /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/scripts/perSVade.py -r /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/5478_Candida_glabrata/reference_genome_dir/reference_genome.fasta --threads 48 --outdir /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/5478_Candida_glabrata/testing_Accuracy/fast/BG2_ANI --nvars 50 --nsimulations 2 --simulation_ploidies haploid,diploid_hetero --range_filtering_benchmark theoretically_meaningful --mitochondrial_chromosome mito_C_glabrata_CBS138 -f1 /gpfs/projects/bsc40/mschikora/Cglabrata_antifungals/data/trimmed_reads/RUN1_BG2_11B_ANI_R1_trimmed.fq.gz -f2 /gpfs/projects/bsc40/mschikora/Cglabrata_antifungals/data/trimmed_reads/RUN1_BG2_11B_ANI_R2_trimmed.fq.gz --fast_SVcalling --previous_repeats_table /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/5478_Candida_glabrata/reference_genome.fasta.repeats.tab


source ~/.bash_profile
run_interactive_session_debug
conda activate perSVade_env

# testing on drosophila:

cd /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/7227_Drosophila_melanogaster/testing_Accuracy/uniform/sample7240_SRR1210633 

python /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/scripts/perSVade.py -r /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/7227_Drosophila_melanogaster/reference_genome_dir/reference_genome.fasta --threads 48 --outdir /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/7227_Drosophila_melanogaster/testing_Accuracy/uniform/sample7240_SRR1210633 --nvars 50 --nsimulations 2 --simulation_ploidies haploid,diploid_hetero --range_filtering_benchmark theoretically_meaningful --mitochondrial_chromosome KJ947872.2 -f1 /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/7227_Drosophila_melanogaster/findingRealSVs_automaticFindingOfCloseReads/getting_closeReads/reads/SRR1210633/SRR1210633_trimmed_reads_1.fastq.gz.30x.fastq.gz -f2 /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/7227_Drosophila_melanogaster/findingRealSVs_automaticFindingOfCloseReads/getting_closeReads/reads/SRR1210633/SRR1210633_trimmed_reads_2.fastq.gz.30x.fastq.gz --previous_repeats_table /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/7227_Drosophila_melanogaster/reference_genome_dir/reference_genome.fasta.repeats.tab --skip_cleaning_outdir --verbose &

# testing on Arabidopsis

cd /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/3702_Arabidopsis_thaliana/testing_Accuracy/uniform/sample2608267_ERR3514863


python /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/scripts/perSVade.py -r /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/3702_Arabidopsis_thaliana/reference_genome_dir/reference_genome.fasta --threads 48 --outdir /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/3702_Arabidopsis_thaliana/testing_Accuracy/uniform/sample2608267_ERR3514863 --nvars 50 --nsimulations 2 --simulation_ploidies haploid,diploid_hetero --range_filtering_benchmark theoretically_meaningful --mitochondrial_chromosome BK010421.1,AP000423.1 -f1 /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/3702_Arabidopsis_thaliana/findingRealSVs_automaticFindingOfCloseReads/getting_closeReads/reads/ERR3514863/ERR3514863_trimmed_reads_1.fastq.gz.30x.fastq.gz -f2 /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/3702_Arabidopsis_thaliana/findingRealSVs_automaticFindingOfCloseReads/getting_closeReads/reads/ERR3514863/ERR3514863_trimmed_reads_2.fastq.gz.30x.fastq.gz --previous_repeats_table /gpfs/projects/bsc40/mschikora/scripts/perSVade/perSVade_repository/testing/outdirs_testing_severalSpecies/3702_Arabidopsis_thaliana/reference_genome_dir/reference_genome.fasta.repeats.tab --skip_cleaning_outdir --verbose &


"""






