/*
 * ArrayHeapExtMem.hpp
 *
 *  Created on: Mar 15, 2016
 *      Author: michel
 */

#ifndef ARRAYHEAPEXTMEM_HPP
#define ARRAYHEAPEXTMEM_HPP

#include "BinaryHeap_FourInts.hpp"
#include "BinaryHeap_TwoInts.hpp"
#include <math.h>
#include <list>
#include <float.h>
class Node {
	public:
		int i;
		int j;
		float key;


		Node(int i, int j, float key) {
			this->i = i;
			this->j = j;
			this->key = key;
		}

		Node() {
			this->i = 0;
			this->j = 0;
			this->key = 0.0;
		}
};
class SlotPair{
	public:
		int slot;
		long remaining;
		SlotPair(int s, long r){
			slot = s;
			remaining = r;
		}
		static int compareTo(const void* a, const void* b){
			   return ((SlotPair*)b)->remaining < ((SlotPair*)a)->remaining ? -1 : 1;
		}
};

typedef struct HeapReturn{
	void* h;
	bool which;
}HeapReturn;

class ArrayHeapExtMem {

	public:

		ArrayHeapExtMem(std::string dir, int* activeIJs);
		ArrayHeapExtMem(std::string dir, int* activeIJs, long sizeExp);
		~ArrayHeapExtMem();

		std::string tmpDir; //path to directory
		std::string fileName; //path to directory

		static const bool chopBottom = true;
		static const bool trimInactiveFromHeapArray = true;
		static int numArrays;

		int A;
		int B;

		int maxLevels;
		float c;  // was 1/7 in the paper, but I'm storing other values with the key, so I've changed the ratio
		int blockSize; //4KB blocks = ~1000 ints
		long mem;
		int cM;
		int numSlots;
		int numNodesPerBlock;
		int numFieldsPerBlock;
		int n;
		BinaryHeap_TwoInts* H1; //the one I have with a float key instead of int
		BinaryHeap_FourInts* H2; //same thing as the above, with 4 ints instead of 2


		int loopCntr;
		bool startCheckingValues;


		FILE* file;

		long** slotNodeCnt;
		int*** perSlotIntBuffer;
		int** slotBuffPos;
		int** cntOnHeap;
		long** slotPositions;

		std::vector<std::list<int>*> *freeSlots;

		std::vector<SlotPair*> *slotPairList;

		char* buffB;
		int* buffI;
		char* bigBuffB;
		int* bigBuffI;

		int numFields;
		long* cntMax;

		BinaryHeap_TwoInts* sortingHeap;

		FILE* tempFile;

		void closeAndDeleteFile ();
		void insert(int i, int j, float key);
		void prepare();
		void clearAndInitialize();
		bool isEmpty();
		void describeHeap();
		HeapReturn getBinaryHeapWithMin();
		void removeMin();
		int size();
		bool test(bool verbose);
		void deleteAll();

	private:
		int* active;

		void initialize(std::string dir, int* activeIJs, long sizeExp);
		long calcPos (int level, int slot);
		int mergeLevels (int targetLevel, int* is, int* js, Float keys);
		bool mergeSlots (int lvl);
		bool mergeTwoSlots (int lvl);
		void load(int level, int slot);
		int store(int level, int* is, int* js, float* keys, int cnt);
		void deleteLevelAndSlots (BinaryHeap_FourInts* H, int level, int* slots, int numSlots);
		void deleteLevels (BinaryHeap_FourInts* H, int* levels, int numLevels);

};
#endif
