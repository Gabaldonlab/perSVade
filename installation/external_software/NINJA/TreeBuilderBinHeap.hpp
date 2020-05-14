/*
 * TreeBuilderBinHeap.hpp
 *
 *  Created on: Feb 13, 2016
 *      Author: michel
 */
#include "Stack.hpp"
#include "BinaryHeap_IntKey_TwoInts.hpp"
#include "TreeBuilder.hpp"
#include <limits.h>

class TreeBuilderBinHeap: public TreeBuilder{
	public:
		TreeBuilderBinHeap(std::string** names, int** distances, int namesSIze, int* clusterEqual);
		TreeBuilderBinHeap(std::string** names, int** distances, int namesSIze);
		~TreeBuilderBinHeap();

		int* candidateCountPerLoop;

		TreeNode** build ();

	private:
		BinaryHeap_IntKey_TwoInts*** heaps;

		int* clustAssignments;
		long* clustPercentiles;
		int* clustersBySize;
		int* clustSizes;

		Int candidatesD;
		Int candidatesI;
		Int candidatesJ;
		bool* candidatesActive;
		Stack *freeCandidates;
		int lastCandidateIndex;

		void clusterAndHeap (int maxIndex);
		int appendCandidate (int d, int i, int j);
};


