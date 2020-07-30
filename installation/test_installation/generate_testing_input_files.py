#!/usr/bin/env python

######### define environment ##########

# this file is to generate the testing inputs for testing that the installation went well

# module imports
import sys
import os
import Bio.SeqIO as SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 

# get the cwd were all the scripts are 
test_dir = "/".join(__file__.split("/")[0:-1]); sys.path.insert(0, test_dir)

# define the EnvDir where the environment is defined
EnvDir = "/".join(sys.executable.split("/")[0:-2])

# define the testing_inputs dir, where everything should be written
testing_inputs_dir = "%s/testing_inputs"%test_dir

# define the full genome and annotations
full_genome = "%s/Candida_glabrata.fasta"%testing_inputs_dir
full_gff = "%s/Candida_glabrata.gff"%testing_inputs_dir

########################################

# load the gff and only the first 
df_gff3 = pd.read_csv(full_gff, skiprows=list(range(len([line for line in open(full_gff, "r") if line.startswith("#")]))), sep="\t", names=["#chromosome", "source", "type_feature", "start", "end", "score", "strand", "phase", "attributes"])

# write only the first 3 genes of 2 chromosomes
all_chroms = {"ChrA_C_glabrata_CBS138", "mito_C_glabrata_CBS138"}
df_gff_reduced = pd.concat([df_gff3[(df_gff3.type_feature=="gene") & (df_gff3["#chromosome"]==chrom)].sort_values(by=["#chromosome", "start"]).iloc[0:15] for chrom in all_chroms])

# write the reduced gff
df_gff_reduced.to_csv("%s/reduced_annotation.gff"%testing_inputs_dir, sep="\t", index=False, header=True)

# write the fasta
chrom_to_seq = {seq.id : seq for seq in SeqIO.parse(full_genome, "fasta")}
all_records = []
all_records_mutated = []

for chrom in all_chroms:

	# get the maximum coordinate
	max_coord = max(df_gff_reduced[(df_gff_reduced["#chromosome"]==chrom)].end) + 1000

	# keep only the beginning
	record = chrom_to_seq[chrom][0:max_coord]
	record.id = chrom
	record.description = ""
	record.name = ""

	# add some repeats
	for I in range(300): record += SeqRecord(Seq("CACACACA"), description="", name="", id=chrom)

	all_records.append(record)

	# keep the mutated one
	mutated_record = SeqRecord(Seq(str(record.seq).replace("AACTGAG", "AACTGAC")), description="", name="", id="mutated_%s"%chrom)
	all_records_mutated.append(mutated_record)

SeqIO.write(all_records, "%s/reduced_genome.fasta"%testing_inputs_dir,"fasta")
SeqIO.write(all_records_mutated, "%s/reduced_genome_mutated.fasta"%testing_inputs_dir,"fasta")
