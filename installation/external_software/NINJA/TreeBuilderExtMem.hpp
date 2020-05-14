/*
 * ArrayHeapExtMem.hpp
 *
 *  Created on: Mar 15, 2016
 *      Author: michel
 */

#ifndef TREEBUILDEREXTMEM_HPP
#define TREEBUILDEREXTMEM_HPP

#include "Stack.hpp"
#include "TreeNode.hpp"
#include "ArrayHeapExtMem.hpp"
#include <float.h>
#include "TreeBuilder.hpp"

class CandidateHeap;

class TreeBuilderExtMem {

	public:
		TreeBuilderExtMem (std::string** names, int namesSize, float *R, std::string njTmpDir, FILE* diskD, float** memD, int memDFirstSize, int firstMemCol, int rowLength, long maxMemory);
		~TreeBuilderExtMem ();

		int K;
		std::string** names;
		TreeNode** nodes;
		int* redirect;

		static const bool returnCandsToHeaps = false;
		static const bool useCandHeaps = true;
		static const bool useBinPairHeaps = true;
		static const int complexCandidateCount = 2000;
		static const int complexCandidateRatio = 40;
		static constexpr float candHeapDecay = 0.6f;
		static const int clustCnt = 30;
		static const int candHeapThresh = 50000;
		int* nextActiveNode;
		int* prevActiveNode;
		int firstActiveNode;

		ArrayHeapExtMem*** arrayHeaps;

		int* clustAssignments;
		float* clustMaxes;
		float* clustMins;
		int* clustersBySize;
		int* clustSizes;

		int newK;

		FILE* diskD=NULL;
		FILE* candFile=NULL;
		long candFilePos;

		int numCandTriplesToDisk; // 16 * 3 pages (or so)


		float* R;
		int RSize;
		float** memD;
		int firstColInMem;
		int curColInMem;
		int memDSize;

		float* fBuff;
		int* candidateCountPerLoop;
		int* candidateViewsPerLoop;
		int* candidateRowsCountPerLoop;

		std::string njTmpDir;

		int nextInternalNode;

		float* candidatesD;
		int* candidatesI;
		int* candidatesJ;
		bool* candidatesActive;
		Stack* freeCandidates;
		int lastCandidateIndex;

		CandidateHeap* starterCandHeap;

		std::vector<CandidateHeap*>* candHeapList;

		bool usingSimpleCandidates;

		float candidatesSize;

		static const int floatSize = 4;
		int rowLength;

		long maxMemory;

		TreeNode** build ();
		void clear();
	private:
		void clusterAndHeap (int maxIndex );
		void returnCandidates ();
		void appendCandidate (float d, int i, int j);
		int appendCandToSimpleList (float d, int i, int j);
		void removeCandidate (int x);
};
#endif
