/*
 * BinaryHeap_IntKey_TwoInts.cpp
 *
 *  Created on: Feb 14, 2016
 *      Author: michel
 */

#include "BinaryHeap_IntKey_TwoInts.hpp"

bool compare3Ints(ints3 x, ints3 y){
	return (x.key > y.key);
}


BinaryHeap_IntKey_TwoInts::BinaryHeap_IntKey_TwoInts( ) {
	this->heap = new std::vector< ints3 >();
	this->heap->reserve(this->DEFAULT_CAPACITY);
	std::make_heap (this->heap->begin(),this->heap->end(),compare3Ints);
}

BinaryHeap_IntKey_TwoInts::BinaryHeap_IntKey_TwoInts(int maxCapacity ) {
	this->heap = new std::vector< ints3 >();
	this->heap->reserve(maxCapacity);
	std::make_heap (this->heap->begin(),this->heap->end(),compare3Ints);
}
BinaryHeap_IntKey_TwoInts::BinaryHeap_IntKey_TwoInts(int* val1s, const int* val2s, Int keys) {
	this->heap = new std::vector< ints3 >();
	this->heap->reserve(this->DEFAULT_CAPACITY);
	std::make_heap (this->heap->begin(),this->heap->end(),compare3Ints);
	for (int i=0;i<keys.length;i++){
		ints3 aux;
		aux.first =val1s[i];
		aux.second = val2s[i];
		aux.key = keys.pointer[i];
		this->heap->push_back(aux);
		std::push_heap (this->heap->begin(),this->heap->end(),compare3Ints);
	}
}

BinaryHeap_IntKey_TwoInts::BinaryHeap_IntKey_TwoInts(int* val1s, const int* val2s, Int keys, int maxCapacity) {
	this->heap = new std::vector< ints3 >();
	this->heap->reserve(maxCapacity);
	std::make_heap (this->heap->begin(),this->heap->end(),compare3Ints);
	for (int i=0;i<keys.length;i++){
		ints3 aux;
		aux.first =val1s[i];
		aux.second = val2s[i];
		aux.key = keys.pointer[i];
		this->heap->push_back(aux);
		std::push_heap (this->heap->begin(),this->heap->end(),compare3Ints);
	}
}

BinaryHeap_IntKey_TwoInts::~BinaryHeap_IntKey_TwoInts(){
	makeEmpty();
	delete this->heap;
}


void BinaryHeap_IntKey_TwoInts::insert(int val1, int val2, int key){
	ints3 aux;
	aux.first =val1;
	aux.second = val2;
	aux.key = key;
	this->heap->push_back(aux);
	std::push_heap (this->heap->begin(),this->heap->end(),compare3Ints);

}

void BinaryHeap_IntKey_TwoInts::deleteMin(){
	std::pop_heap(this->heap->begin(),this->heap->end(),compare3Ints);
	this->heap->pop_back();
}


bool BinaryHeap_IntKey_TwoInts::isEmpty(){
	return this->heap->empty();
}

int BinaryHeap_IntKey_TwoInts::size(){
	return this->heap->size();
}

void BinaryHeap_IntKey_TwoInts::makeEmpty(){
	if(this->heap == NULL) return;
	this->heap->clear();
}
int compareMyType2 (const void * a, const void * b)
{
  if ( *(int*)a <  *(int*)b ) return -1;
  if ( *(int*)a == *(int*)b ) return 0;
  return 1;
}
bool BinaryHeap_IntKey_TwoInts::binHeapTwoTest(bool verbose){
	printf("Binary Heap Two Ints Test...\n");

	srand (time(NULL));
	int *int_aux = new int[10000];
	for(int l=1;l<6;l++){
		BinaryHeap_IntKey_TwoInts *aux;
		for(int k=1;k<6;k++){

			if (verbose) printf("Default initializing...\n");
			aux = new BinaryHeap_IntKey_TwoInts();


			if (verbose) printf("Inserting...\n");
			for(int i=0;i<10000;i++){
				int_aux[i] = rand()%100000;
				aux->insert(i*k*l,i*k*l*2,int_aux[i]);
			}
			qsort((void*)int_aux,10000,sizeof(int),compareMyType2);

			if (verbose) printf("Deleting and asserting...\n");
			for(int i=0;i<5000;i++){
				int test = aux->heap->front().key;
				int test2 = int_aux[i];
				assert(test2==test);
				aux->deleteMin();
			}

			if (verbose) printf("Assert size.\n");
			assert(aux->size() == 5000);

			if (verbose) printf("Inserting...\n");
			for(int i=0;i<5000;i++){
				int_aux[i] = rand()%100000;
				aux->insert(i*k*l,i*k*l*2,int_aux[i]);
			}

			if (verbose) printf("Assert size.\n");
			assert(aux->size() == 10000);

			qsort((void*)int_aux,10000,sizeof(int),compareMyType2);

			if (verbose) printf("Deleting...\n");
			for(int i=0;i<10000;i++){
				int test = aux->heap->front().key;
				int test2 = int_aux[i];
				assert(test2==test);
				aux->deleteMin();
			}

			if (verbose) printf("Assert size.\n");
			assert(aux->size() == 0);
			assert(aux->isEmpty() == true);

			if (verbose) printf("Emptying...\n");
			aux->makeEmpty();

			if (verbose) printf("Emptying again(check for double free)...\n");
			aux->makeEmpty();

			if (verbose) printf("Binary Heap Two Ints Test %d completed successfully.\n",k);
		}
	}
	printf("Binary Heap Two Ints Test completed successfully.\n");
	return true;
}
