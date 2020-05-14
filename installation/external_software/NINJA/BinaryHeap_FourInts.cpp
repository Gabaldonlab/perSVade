/*
 * BinaryHeap_FourInts.cpp
 *
 *  Created on: Mar 19, 2016
 *      Author: michel
 */

#include "BinaryHeap_FourInts.hpp"

bool compare4Ints(ints4float x, ints4float y){
	return (x.key > y.key);
}


BinaryHeap_FourInts::BinaryHeap_FourInts( ) {
	this->heap = new std::vector< ints4float >();
	this->heap->reserve(this->DEFAULT_CAPACITY);
	std::make_heap (this->heap->begin(),this->heap->end(),compare4Ints);
}

BinaryHeap_FourInts::BinaryHeap_FourInts(int maxCapacity ) {
	this->heap = new std::vector< ints4float >();
	this->heap->reserve(maxCapacity);
	std::make_heap (this->heap->begin(),this->heap->end(),compare4Ints);
}
BinaryHeap_FourInts::BinaryHeap_FourInts(int* val1s, const int* val2s, const int* val3s, const int* val4s, Float keys) {
	this->heap = new std::vector< ints4float >();
	this->heap->reserve(this->DEFAULT_CAPACITY);
	std::make_heap (this->heap->begin(),this->heap->end(),compare4Ints);
	for (int i=0;i<keys.length;i++){
		ints4float aux;
		aux.first =val1s[i];
		aux.second = val2s[i];
		aux.third = val3s[i];
		aux.fourth = val4s[i];
		aux.key = keys.pointer[i];
		this->heap->push_back(aux);
		std::push_heap (this->heap->begin(),this->heap->end(),compare4Ints);
	}
}

BinaryHeap_FourInts::BinaryHeap_FourInts(int* val1s, const int* val2s, const int* val3s, const int* val4s, Float keys, int maxCapacity) {
	this->heap = new std::vector< ints4float >();
	this->heap->reserve(maxCapacity);
	std::make_heap (this->heap->begin(),this->heap->end(),compare4Ints);
	for (int i=0;i<keys.length;i++){
		ints4float aux;
		aux.first =val1s[i];
		aux.second = val2s[i];
		aux.key = keys.pointer[i];
		this->heap->push_back(aux);
		std::push_heap (this->heap->begin(),this->heap->end(),compare4Ints);
	}
}
void BinaryHeap_FourInts::buildAgain(int* val1s, const int* val2s, const int* val3s, const int* val4s, Float keys){
	makeEmpty();
	this->heap = new std::vector< ints4float >();
	std::make_heap (this->heap->begin(),this->heap->end(),compare4Ints);
	for (int i=0;i<keys.length;i++){
		ints4float aux;
		aux.first =val1s[i];
		aux.second = val2s[i];
		aux.third = val3s[i];
		aux.fourth = val4s[i];
		aux.key = keys.pointer[i];
		this->heap->push_back(aux);
		std::push_heap (this->heap->begin(),this->heap->end(),compare4Ints);
	}
}

BinaryHeap_FourInts::~BinaryHeap_FourInts(){
	makeEmpty();
}


void BinaryHeap_FourInts::insert(int val1, int val2, int val3, int val4, float key){
	ints4float aux;
	aux.first =val1;
	aux.second = val2;
	aux.third = val3;
	aux.fourth = val4;
	aux.key = key;
	this->heap->push_back(aux);
	std::push_heap (this->heap->begin(),this->heap->end(),compare4Ints);

}

void BinaryHeap_FourInts::deleteMin(){
	std::pop_heap(this->heap->begin(),this->heap->end(),compare4Ints);
	this->heap->pop_back();
}


bool BinaryHeap_FourInts::isEmpty(){
	return this->heap->empty();
}

int BinaryHeap_FourInts::size(){
	return this->heap->size();
}

void BinaryHeap_FourInts::makeEmpty(){
	if(this->heap == NULL) return;
	this->heap->clear();
	delete this->heap;
	this->heap = new std::vector< ints4float >();
	std::make_heap (this->heap->begin(),this->heap->end(),compare4Ints);
}
int compareMyType4 (const void * a, const void * b)
{
  if ( ((ints4float*)a)->key <  ((ints4float*)b)->key ) return -1;
  if ( ((ints4float*)a)->key == ((ints4float*)b)->key ) return 0;
  return 1;
}
bool BinaryHeap_FourInts::binHeapFourTest(bool verbose){
	printf("Binary Heap Four Ints Test (Float Key)...\n");

	srand (time(NULL));
	ints4float int_aux[10000];
	for(int l=1;l<6;l++){
		BinaryHeap_FourInts *aux;
		for(int k=1;k<6;k++){

			if (verbose) printf("Default initializing...\n");
			aux = new BinaryHeap_FourInts();

			if (verbose) printf("Inserting...\n");
			for(int i=0;i<10000;i++){
				int_aux[i].first = rand()%100000;
				int_aux[i].second = rand()%100000;
				int_aux[i].third = rand()%100000;
				int_aux[i].fourth = rand()%100000;
				int_aux[i].key = rand()%100000;
				aux->insert(int_aux[i].first,int_aux[i].second,int_aux[i].third, int_aux[i].fourth,int_aux[i].key);
			}
			qsort((void*)int_aux,10000,sizeof(ints4float),compareMyType4);

			if (verbose) printf("Deleting and asserting...\n");
			for(int i=0;i<5000;i++){

				assert(aux->heap->front().key==int_aux[i].key);
				aux->deleteMin();
			}

			if (verbose) printf("Assert size.\n");
			assert(aux->size() == 5000);

			if (verbose) printf("Inserting...\n");
			for(int i=0;i<5000;i++){
				int_aux[i].first = rand()%100000;
				int_aux[i].second = rand()%100000;
				int_aux[i].third = rand()%100000;
				int_aux[i].fourth = rand()%100000;
				int_aux[i].key = rand()%100000;
				aux->insert(int_aux[i].first,int_aux[i].second,int_aux[i].third, int_aux[i].fourth,int_aux[i].key);
			}

			if (verbose) printf("Assert size.\n");
			assert(aux->size() == 10000);

			qsort((void*)int_aux,10000,sizeof(int),compareMyType4);

			if (verbose) printf("Deleting...\n");
			for(int i=0;i<10000;i++){
				assert(aux->heap->front().key==int_aux[i].key);
				aux->deleteMin();
			}

			if (verbose) printf("Assert size.\n");
			assert(aux->size() == 0);
			assert(aux->isEmpty() == true);

			if (verbose) printf("Emptying...\n");
			aux->makeEmpty();

			if (verbose) printf("Emptying again(check for double free)...\n");
			aux->makeEmpty();

			if (verbose) printf("Binary Heap Four Ints Test %d completed successfully.\n",k);
		}
	}
	printf("Binary Heap Four Ints Test completed successfully.\n");
	return true;
}


