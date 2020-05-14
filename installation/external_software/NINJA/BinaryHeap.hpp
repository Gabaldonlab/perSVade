//includes
#ifndef BINARYHEAP_HPP
#define BINARYHEAP_HPP

#include <limits.h>
#include <stdio.h>
#include <assert.h>
#include "ExceptionHandler.hpp"

struct Int{
	int* pointer;
	int length;
};
struct Float{
	float* pointer;
	int length;
};

class BinaryHeap {
public:
	//constructors
	BinaryHeap();
	BinaryHeap(int maxCapacity);
	BinaryHeap(const int* val1s, Int keys);
	BinaryHeap(const int* val1s, Int keys, int maxCapacity);

	~BinaryHeap();


	int size(); //return BinaryHeap size

	//functions
	int insert(int val1, int key); //insert element
	void deleteMin(); //Remove the smallest item from the priority queue

	int getMin();
	bool isEmpty(); //check if empty
	void makeEmpty(); //empty heap

	bool binHeapTest(bool verbose);

	std::vector< std::pair<int,int> > *heap;


	static const int DEFAULT_CAPACITY = 1000;


};
#endif
