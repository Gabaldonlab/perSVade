/*
 * CandidateHeap.hpp
 *
 *  Created on: Mar 18, 2016
 *      Author: michel
 */

#ifndef CANDIDATEHEAP_HPP
#define CANDIDATEHEAP_HPP

#include "ArrayHeapExtMem.hpp"
#include "TreeBuilderExtMem.hpp"

class CandidateHeap: public ArrayHeapExtMem{

	public:
		CandidateHeap(std::string dir, int* activeIJs, int kPrime, TreeBuilderExtMem* tb, long sizeExp);
		CandidateHeap(std::string dir, int* activeIJs, int kPrime, TreeBuilderExtMem* tb);
		~CandidateHeap();


		int kPrime;
		float* rPrimes;
		float* rDeltas;
		int* rowCounts;

		int* nextActiveNode;
		int* prevActiveNode;
		int firstActiveNode;

		float k_over_kprime;

		float minDeltaSum;

		int origSize;
		bool expired;
		int representedRowCount;

		TreeBuilderExtMem* tb;

		void insert(int i, int j, float key);
		void buildNodeList();
		void removeMin();
		void calcDeltaValues(int newK);
		void clear();

	private:
		void initialize (int kPrime, TreeBuilderExtMem* tb);
};
#endif
