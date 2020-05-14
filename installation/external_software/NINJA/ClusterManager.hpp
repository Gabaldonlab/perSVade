/*
 * ClusterManager.hpp
 *
 *  Created on: Jul 12, 2019
 *      Author: robert
 */

#include "ExceptionHandler.hpp"
#include "DistanceReader.hpp"
#include "SequenceFileReader.hpp"
#include <stdlib.h>
//#include "DistanceReaderExtMem.hpp"

class ClusterManager {
	struct Float{
		float* pointer;
		int length;
	};
	public:

		enum InputType {alignment, distance};
		enum AlphabetType {dna, amino, null};
		enum CorrectionType {not_assigned, none, JukesCantor/*DNA*/, Kimura2/*DNA*/, FastTree /*amino*/, MismatchesOneGap /*mismatches + 1 for each gap*/};
		enum OutputType {dist, tree, cluster};

		ClusterManager(std::string method, std::string njTmpDir, std::string inFile, FILE* outfile, InputType
        inType, OutputType outType, AlphabetType alphType, CorrectionType corrType, int threads, bool useSSE, bool
        printTime, float clusterCutoff);

		std::string method;
		std::string njTmpDir;
		std::string inFile;
		FILE* outFile;
		std::string** names;
		InputType inType;
		OutputType outType;
		AlphabetType alphType;
		CorrectionType corrType;
		const static int NUM_RAND_DIR_CHARS = 6;
		std::string chars;
		int threads;
                float clusterCutoff;
        bool printTime;

		bool newDistanceMethod;

		std::string doJob();
		static std::string getUniqueID();
		static bool deleteDir(FILE* dir);

};
