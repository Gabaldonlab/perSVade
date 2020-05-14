/*
 * CandidateHeap.cpp
 *
 *  Created on: Mar 27, 2016
 *      Author: michel
 */

#include "CandidateHeap.hpp"

CandidateHeap::CandidateHeap(std::string dir, int* activeIJs, int kPrime, TreeBuilderExtMem *tb, long sizeExp) : ArrayHeapExtMem(dir,activeIJs,sizeExp){
	initialize( kPrime, tb);
}

CandidateHeap::CandidateHeap(std::string dir, int* activeIJs, int kPrime, TreeBuilderExtMem *tb) : ArrayHeapExtMem(dir,activeIJs){
	initialize( kPrime, tb);
}
CandidateHeap::~CandidateHeap(){
	clear();
}

void CandidateHeap::initialize (int kPrime, TreeBuilderExtMem *tb) {

	this->firstActiveNode = -1;
	this->origSize = 0;
	this->expired = false;
	this->representedRowCount = 0;


	this->tb = tb;
	this->kPrime = kPrime;
	rowCounts = new int[tb->nextInternalNode+1]();
	rDeltas = new float[tb->nextInternalNode+1]();
	nextActiveNode = new int[tb->nextInternalNode+1]();
	prevActiveNode = new int[tb->nextInternalNode+1]();

	rPrimes = new float[tb->RSize]();
	for (int i=0; i< tb->RSize; i++ ) {
		rPrimes[i] = tb->R[i];
	}
}

void CandidateHeap::insert(int i, int j, float key){
	ArrayHeapExtMem::insert(i, j, key);

	rowCounts[i]++;
	rowCounts[j]++;
}


void CandidateHeap::buildNodeList () {
	origSize = ArrayHeapExtMem::size();

	int prev=-1;
	for (int i=0; i<tb->nextInternalNode+1; i++) {
		if (rowCounts[i] > 0) {
			representedRowCount++;
			if (firstActiveNode==-1) {
				firstActiveNode = i;
			} else {
				prevActiveNode[i] = prev;
				nextActiveNode[prev] = i;
			}
			prev = i;
		}
	}
	prevActiveNode[0] = -1;
	nextActiveNode[prev] = -1;
}


void CandidateHeap::removeMin(){
	HeapReturn x = ArrayHeapExtMem::getBinaryHeapWithMin();
	int i = 0,j = 0;
	if(x.which){
		BinaryHeap_FourInts* H = (BinaryHeap_FourInts*)x.h;
		i = H->heap->front().first;
		j = H->heap->front().second;
	}else{
		BinaryHeap_TwoInts* H = (BinaryHeap_TwoInts*)x.h;
		i = H->heap->front().first;
		j = H->heap->front().second;
	}


	int prev = 0, next = 0;
	if (--rowCounts[i] == 0) { // compact list
		representedRowCount--;
		prev = prevActiveNode[i];
		next = nextActiveNode[i];
		if (next != -1) prevActiveNode[next] = prev;
		if (prev != -1) nextActiveNode[prev] = next;
	}

	if (--rowCounts[j] == 0) { // compact list
		representedRowCount--;
		prev = prevActiveNode[j];
		next = nextActiveNode[j];
		if (next != -1) prevActiveNode[next] = prev;
		if (prev != -1) nextActiveNode[prev] = next;
	}

	ArrayHeapExtMem::removeMin();
}



void CandidateHeap::calcDeltaValues (int newK) {
	//prevActiveNode[0] = -1;
	//nextActiveNode[0] = -1;
	int x=firstActiveNode;

	float minRdelt1, minRdelt2;
	minRdelt1 = FLT_MAX;
	minRdelt2 = FLT_MAX;

	k_over_kprime = (float)(newK-2)/(kPrime-2);

	int rx, prev, next;
	while (x != -1) {
		rx = tb->redirect[x];

		if (rx == -1) {
			prev = prevActiveNode[x];
			next = nextActiveNode[x];
			if (next != -1) prevActiveNode[next] = prev;
			if (prev != -1) nextActiveNode[prev] = next;
			x = next;
		} else {
			rDeltas[x] = k_over_kprime * rPrimes[rx]  - tb->R[rx];

			if (rDeltas[x] < minRdelt1) {
				minRdelt2 = minRdelt1;
				minRdelt1 = rDeltas[x];
			} else if (rDeltas[x] < minRdelt2) {
				minRdelt2 = rDeltas[x];
			}

			x = nextActiveNode[x];
		}
	}
	minDeltaSum =  minRdelt1 + minRdelt2;

}
void CandidateHeap::clear(){

	ArrayHeapExtMem::deleteAll();
	if (rPrimes!=NULL){
		delete[] rPrimes;
		rPrimes = NULL;
	}
	if (rDeltas!=NULL){
		delete[] rDeltas;
		rDeltas = NULL;
	}
	if (rowCounts!=NULL){
		delete[] rowCounts;
		rowCounts = NULL;
	}
	if (nextActiveNode!=NULL){
		delete[] nextActiveNode;
		nextActiveNode = NULL;
	}
	if (prevActiveNode!=NULL){
		delete[] prevActiveNode;
		prevActiveNode = NULL;
	}
	/*
	if (this->tb!=NULL){
		delete this->tb;
		tb = NULL;
	}*/
}
