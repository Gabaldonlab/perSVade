# This is a script to clean the filters samples

import os

for species_dir_name in os.listdir("./outdirs_testing_severalSpecies"):
	
	species_dir = "./outdirs_testing_severalSpecies/%s"%species_dir_name
	if not os.path.isdir(species_dir) or "reference_genome_dir" not in os.listdir(species_dir): continue

	for typeSim in ["realSVs", "uniform"]:
		typeSim_dir = species_dir+"/testing_Accuracy/"+typeSim
		for sampleID in os.listdir(typeSim_dir):
			sample_dir = typeSim_dir+"/"+sampleID

			if "perSVade_finished_file.txt" in os.listdir(sample_dir):
				
				for simID in [1, 2]:
					for ploidy in ["haploid", "diploid_hetero"]:

						af = {"haploid":1.00, "diploid_hetero":0.50}[ploidy]


						folder_remove = "%s/SVdetection_output/parameter_optimisation/simulation_%i/benchmark_GridssClove_%s/benchmark_max50000x_ignoreRegionsFalse/several_parameter_combinations_filter_theoretically_meaningful_af%.2f"%(sample_dir, simID, ploidy, af)

						print(folder_remove)

						if os.path.isdir(folder_remove): os.system("rm -r %s"%folder_remove)
