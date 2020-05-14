/*
 * TreeBuilder.cpp
 *
 *  Created on: Jan 24, 2016
 *      Author: michel
 */
#include "TreeBuilder.hpp"
#include <assert.h>

bool TreeBuilder::distInMem = false;
bool TreeBuilder::rebuildStepsConstant = false; // otherwise, decreasing
float TreeBuilder::rebuildStepRatio = (float)0.5;
int TreeBuilder::rebuildSteps  = -1;
//TreeBuilder::clustCnt = 30;
int TreeBuilder::candidateIters = 50;
int TreeBuilder::verbose = 1;

TreeBuilder::TreeBuilder (std::string** names, int** distances, int namesSize){ //standard constructor
	this->firstActiveNode = 0;

	this->nextActiveNode = NULL;
	this->prevActiveNode = NULL;


	this->names = names;
	this->D = distances;
	this->K = namesSize;

	int i,j;

	this->R = new long[this->K]();

	this->redirect = new int[(this->K*2)-1];

	this->nodes = new TreeNode*[(this->K*2)-1];

	for (i=0; i<this->K; i++) {
			this->redirect[i] = i;
			this->nodes[i] = new TreeNode();
	}

	//initialize the internal nodes that will be used in the clustering
	for (i=this->K; i<2*this->K-1; i++) {
		this->redirect[i] = -1;
		this->nodes[i] = new TreeNode();
	}

	for (i=0; i<this->K-1; i++) {
		for (j=i+1; j<this->K; j++) {
			this->R[i] += (long)this->D[i][j-i-1];
			this->R[j] += (long)this->D[i][j-i-1];
		}
	}
}

TreeBuilder::TreeBuilder (std::string** names, int** distances, int namesSize, int* clustersEqual){ //constructor that takes uses equal-clusters to speed-up clustering

	this->firstActiveNode = 0;

	this->nextActiveNode = NULL;
	this->prevActiveNode = NULL;


	this->names = names;
	this->D = distances;
	this->K = namesSize;

	int i,j;

	this->R = new long[this->K]();

	this->redirect = new int[(this->K*2)-1];

	this->nodes = new TreeNode*[(this->K*2)-1];


	//find the equal clusters

//	int clustersEqual[this->K];
//
//	for (j=0; j<this->K; j++) {
//		clustersEqual[j] = -1;
//	}

	int countClustered = 0;
	int numClusters = 0;

	std::vector<int> clusterRoots;

	int numLeft[this->K];
	for (i=0; i<this->K; i++) {
		numLeft[i] = 0;
	}

	for (i=0; i<this->K; i++) {
		if (clustersEqual[i] != -1){
			countClustered++;
			if (numLeft[clustersEqual[i]] == 0){
				numClusters++;
				clusterRoots.push_back(clustersEqual[i]);
			}
			numLeft[clustersEqual[i]]++;
		}
	}

//	for (i=0; i<this->K; i++) {
//		foundEqual = false;
//		numLeft[i] = 0;
////		if (clustersEqual[i] != -1) //if this is linked to some cluster already, ignore it. this is necessary because if A == B and B == C, A ?= C.
////			continue;
//		for (j=i+1; j<this->K; j++) {
//			if (this->D[i][j-i-1] == 0 && clustersEqual[j] == -1){ //equal sequences and cluster not assigned yet
//				if ((*this->seq[i]) == (*this->seq[j])){ //if the two sequences are exactly equal (necessary due to non-transitivity of the relation)
//					clustersEqual[j] = i;
//					countClustered++;
//					if (!foundEqual){
//						numClusters++;
//						clusterRoots.push_back(i);
//					}
//					numLeft[i]++;
//					foundEqual = true;
//				}
//			}
//		}
//	}

	if (numClusters > 0){ //there are equal sequences

		TreeNode** equalNodes = new TreeNode*[this->K];

		//initialize equal sequence root nodes
		for (i=0; i<this->K; i++) {
			equalNodes[i] = new TreeNode();
			//equalNodes[i]->length = 0;
		}

		int nodesAdded = 0;
		int rootsAdded = 0;
		//int countExclusion = 0;
		std::vector<int> nodesPos;

		std::sort(clusterRoots.begin(), clusterRoots.end());

		//initialize sequence nodes
		for (i=0; i<this->K; i++) {
			this->redirect[i] = i;
			TreeNode* newTree = new TreeNode(names[i]);
			if (rootsAdded < numClusters && clusterRoots[rootsAdded] == i){ //if it is the root of an equal subtree I will put the subtree instead
				//this->nodes[i] = equalNodes[nodesAdded];
				this->nodes[nodesAdded] = equalNodes[i];
				newTree->length = 0;
				equalNodes[i]->leftChild = newTree;
				nodesAdded++;
				rootsAdded++;
				nodesPos.push_back(i);
			}else if(clustersEqual[i] == -1){ //should be used in the clustering
				this->nodes[nodesAdded] = newTree;
				//this->nodes[i] = newTree;
				nodesAdded++;
				nodesPos.push_back(i);
			}
//			else{
//				countExclusion++;
//			}
//			this->redirect[i] = i+countExclusion;
		}
//		for (i=0; i<this->K-1; i++) {
//			for (j=0; j<this->K-i-1; j++) {
//				printf("%d ", this->D[i][j]);
//			}
//			printf("\n");
//		}
//
//		printf("---------\n");

		//change distances accordingly

		for (size_t i=0; i<nodesPos.size()-1; i++) {
			for (size_t j=0; j<nodesPos.size()-i-1; j++) {
				//printf("%d,%d = %d, %d\n", i, j, nodesPos[i], nodesPos[i+j+1]-i-1);
				this->D[i][j] = this->D[nodesPos[i]][nodesPos[i+j+1]-nodesPos[i]-1];
			}
		}


//		int shifts;
//
//		for(i=0;i<shiftPos.size();i++){
//			shifts = 0;
//			for(int k=0;k<this->K;k++){
//				//shiftPos[i] -> position of shift
//				//take out i to amend the shifts already done
//				//now transform that to the right index by taking out k and 1
//				memcpy(this->D[k][(shiftPos[i]-i)-k-1], this->D[k][(shiftPos[i]-i)-k], sizeof(int)*((this->K-k-1)-(shiftPos[i]-i)));
//			}
//		}



//		for (i=0; i<this->K-1; i++) {
//			for (j=0; j<this->K-i-1; j++) {
//				this->D[i][j] = this->D[i][this->redirect[j+i+1]-i-1];
//			}
//		}
////		for (i=0; i<this->K-1; i++) {
////			this->D[i] = this->D[this->redirect[i]];
////		}
//
//		for (i=0; i<nodesPos.size()-1; i++) {
//			for (j=0; j<nodesPos.size()-i-1; j++) {
//				printf("%d ", this->D[i][j]);
//			}
//			printf("\n");
//		}


		//initialize internal nodes
		for (i=0; i<this->K; i++) {
			if (clustersEqual[i] != -1){ //if assigned to a cluster, I will append it to the right node.
				TreeNode* newTree = new TreeNode(names[i]);
				TreeNode* newTreeInternal = new TreeNode();
				newTree->length = 0;
				newTreeInternal->length = 0;
				if (numLeft[clustersEqual[i]] != 1){
					equalNodes[clustersEqual[i]]->rightChild = newTreeInternal;
					newTreeInternal->leftChild = newTree;
				}else{ //the last one
					equalNodes[clustersEqual[i]]->rightChild = newTree;
				}
				equalNodes[clustersEqual[i]] = newTreeInternal;
				numLeft[clustersEqual[i]]--;
			}
		}

		//printf("Percentage of nodes taken due to by equality: %d", countClustered*100/this->K);

		this->K -= countClustered;

		assert(nodesAdded == this->K);

	}else{
		for (i=0; i<this->K; i++) {
			this->redirect[i] = i;
			this->nodes[i] = new TreeNode();
		}
	}

	//initialize the internal nodes that will be used in the clustering
	for (i=this->K; i<2*this->K-1; i++) {
		this->redirect[i] = -1;
		this->nodes[i] = new TreeNode();
	}

	for (i=0; i<this->K-1; i++) {
		for (j=i+1; j<this->K; j++) {
			this->R[i] += (long)this->D[i][j-i-1];
			this->R[j] += (long)this->D[i][j-i-1];
		}
	}

}
TreeBuilder::~TreeBuilder(){
	delete[](this->R);
	delete[](this->redirect);
	//delete[](this->nodes); TODO: delete it according to new clustered nodes

}
void TreeBuilder::finishMerging(){ //did not change anything expect the ->`s instead of the .`s
	int last_i=0, last_j=0, last_k=0;
	int ri, rj, rk;
	int i=0;
	while (redirect[i++]==-1);
	ri = redirect[i-1];
	last_i=i-1;

	while (redirect[i++]==-1);
	rj = redirect[i-1];
	last_j=i-1;

	while (redirect[i++]==-1);
	rk = redirect[i-1];
	last_k=i-1;

	int nextNode = 2*K-3;

	nodes[nextNode]->leftChild = nodes[last_i];
	nodes[nextNode]->rightChild = nodes[last_j];

	float d_ij = ri<rj? D[ri][rj] : D[rj][ri];
	nodes[last_i]->length = (d_ij + (R[ri]-R[rj])/2.0) / 20000000.0 ;
	nodes[last_j]->length = (d_ij + (R[rj]-R[ri])/2.0) / 20000000.0 ;

	//if a length is negative, move root of that subtree around to compensate.
	if (nodes[last_i]->length < 0) {
		nodes[last_j]->length -= nodes[last_i]->length;
		nodes[last_i]->length = 0;
	} else if (nodes[last_j]->length < 0) {
		nodes[last_i]->length -= nodes[last_j]->length;
		nodes[last_j]->length = 0;
	}


	float d_ik = ri<rk? D[ri][rk] : D[rk][ri];
	float d_jk = rj<rk? D[rj][rk] : D[rk][rj];
	float d_ijk = (d_ik + d_jk - d_ij) / 2;

	nodes[nextNode+1]->leftChild = nodes[nextNode];
	nodes[nextNode+1]->rightChild = nodes[last_k];

	nodes[nextNode]->length = (d_ijk + (R[ri]-R[rj])/2) / 20000000 ;
	nodes[last_k]->length = (d_ijk + (R[rj]-R[ri])/2) / 20000000 ;

	//if a length is negative, move root of that subtree around to compensate.
	if (nodes[last_i]->length < 0) {
		nodes[last_j]->length -= nodes[last_i]->length;
		nodes[last_i]->length = 0;
	} else if (nodes[last_j]->length < 0) {
		nodes[last_i]->length -= nodes[last_j]->length;
		nodes[last_j]->length = 0;
	}
}

