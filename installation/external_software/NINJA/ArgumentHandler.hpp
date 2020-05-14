#ifndef ARGUMENTHANDLER_HPP
#define ARGUMENTHANDLER_HPP

#include "ExceptionHandler.hpp"
#include <stdlib.h>
#include <stdio.h>
#include "TreeBuilder.hpp"

class ArgumentHandler{
	public:
		enum InputType {alignment, distance};
		enum OutputType {dist, tree, cluster};
		enum AlphabetType {dna, amino, null};
		enum CorrectionType {not_assigned, none, JukesCantor/*DNA*/, Kimura2/*DNA*/, FastTree /*amino*/, MismatchesOneGap /*mismatches + 1 for each gap*/};

		std::string method; // default should be "bin" for small, "cand" for large.
		std::string njTmpDir;
		std::string inFile; //works?
		InputType inType;
		OutputType outType;
		AlphabetType alphType;
		CorrectionType corrType;
		FILE* outFile;

        bool printTime;

		int threads;
      
                float clusterCutoff;

		bool SSE;

		bool abort;

		ArgumentHandler (char* argv[],int argc);
		std::string getMethod();
		std::string getInFile();
		std::string getNJTmpDir();
		InputType getInputType();
		OutputType getOutputType();
		FILE* getOutpuFile();
		int getNumThreads();
                float getClusterCutoff();

		AlphabetType getAlphabetType();

		CorrectionType getCorrectionType ();

		bool argumentTest();

		bool useSSE();

        bool getPrintTime();
};

#endif




