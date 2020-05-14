/*
 * TreeBuilder.hpp
 *
 *  Created on: Jan 24, 2016
 *      Author: michel
 */
#ifndef TREEBUILDER_HPP
#define TREEBUILDER_HPP

#include <string>
#include <vector>
#include <algorithm>
#include "TreeNode.hpp"


class TreeBuilder{
	public:
		TreeBuilder (std::string** names, int** distances, int namesSize);
		TreeBuilder (std::string** names, int** distances, int namesSize, int* clusterEqual);
		~TreeBuilder();
		TreeNode* build ();

		static bool distInMem;
		static bool rebuildStepsConstant;
		static float rebuildStepRatio;
		static int rebuildSteps;
		static const int clustCnt = 30;
		static int candidateIters;
		static int verbose;
		int K;

	protected:
		std::string** names;
		int** D;
		long *R;
		TreeNode **nodes;
		int *redirect;
		int *nextActiveNode;
		int *prevActiveNode;
		int firstActiveNode;

		void finishMerging();

};

#endif
