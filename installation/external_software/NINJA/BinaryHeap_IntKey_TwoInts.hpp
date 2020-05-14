/*
 * BinaryHeap_IntKey.hpp
 *
 *  Created on: Feb 14, 2016
 *      Author: michel
 */

#include "BinaryHeap.hpp"
#include "ExceptionHandler.hpp"
#include <limits.h>

struct ints3{
	int first;
	int second;
	int key;
};


class BinaryHeap_IntKey_TwoInts {
public:
	//constructors
	BinaryHeap_IntKey_TwoInts();
	BinaryHeap_IntKey_TwoInts(int maxCapacity);
	BinaryHeap_IntKey_TwoInts(int* val1s, const int* val2s, Int keys);
	BinaryHeap_IntKey_TwoInts(int* val1s, const int* val2s, Int keys, int maxCapacity);
	~BinaryHeap_IntKey_TwoInts();



	int size(); //return BinaryHeap size

	//functions
	void insert(int val1, int val2, int key);
	void deleteMin(); //Remove the smallest item from the priority queue

	int getMin();
	bool isEmpty(); //check if empty
	void makeEmpty(); //empty heap


	bool binHeapTest(bool verbose);

	std::vector< ints3 > *heap;


	static const int DEFAULT_CAPACITY = 1000;

	bool binHeapTwoTest(bool verbose);

};
