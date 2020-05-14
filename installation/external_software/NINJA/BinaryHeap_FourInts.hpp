/*
 * BinaryHeap_FourInts.hpp
 *
 *  Created on: Mar 19, 2016
 *      Author: michel
 */

#ifndef BINARYHEAP_FOUR_INTS_HPP
#define BINARYHEAP_FOUR_INTS_HPP

#include "ExceptionHandler.hpp"
#include <limits.h>
#include "BinaryHeap.hpp"

struct ints4float{
	int first;
	int second;
	int third;
	int fourth;
	float key;
};



class BinaryHeap_FourInts {
public:
	//constructors
	BinaryHeap_FourInts();
	BinaryHeap_FourInts(int maxCapacity);
	BinaryHeap_FourInts(int* val1s, const int* val2s, const int* val3s, const int* val4s, Float keys);
	BinaryHeap_FourInts(int* val1s, const int* val2s, const int* val3s, const int* val4s, Float keys, int maxCapacity);
	~BinaryHeap_FourInts();



	int size(); //return BinaryHeap size

	//functions
	void insert(int val1, int val2, int val3, int val4, float key);
	void deleteMin(); //Remove the smallest item from the priority queue

	int getMin();
	bool isEmpty(); //check if empty
	void makeEmpty(); //empty heap

	bool binHeapTest(bool verbose);

	void buildAgain(int* val1s, const int* val2s, const int* val3s, const int* val4s, Float keys);

	std::vector< ints4float > *heap;


	static const int DEFAULT_CAPACITY = 1000;

	bool binHeapFourTest(bool verbose);

};
#endif
