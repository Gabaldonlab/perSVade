/*
 * BinaryHeap_TwoInts.hpp
 *
 *  Created on: Mar 19, 2016
 *      Author: michel
 */
#ifndef BINARYHEAP_TWO_INTS_HPP
#define BINARYHEAP_TWO_INTS_HPP

#include "ExceptionHandler.hpp"
#include <limits.h>
#include "BinaryHeap.hpp"

struct ints2float{
	int first;
	int second;
	float key;
};



class BinaryHeap_TwoInts {
public:
	//constructors
	BinaryHeap_TwoInts();
	BinaryHeap_TwoInts(int maxCapacity);
	BinaryHeap_TwoInts(int* val1s, const int* val2s, Float keys);
	BinaryHeap_TwoInts(int* val1s, const int* val2s, Float keys, int maxCapacity);
	~BinaryHeap_TwoInts();



	int size(); //return BinaryHeap size

	//functions
	void insert(int val1, int val2, float key);
	void deleteMin(); //Remove the smallest item from the priority queue

	int getMin();
	bool isEmpty(); //check if empty
	void makeEmpty(); //empty heap
	void buildAgain(int* val1s, const int* val2s, Float keys);
	void chopBottomK(int* val1s, int* val2s, Float keys);


	bool binHeapTest(bool verbose);

	std::vector< ints2float > *heap;


	static const int DEFAULT_CAPACITY = 1000;

	bool binHeapTwoTest(bool verbose);

};
#endif
