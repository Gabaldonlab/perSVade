/*
 * BinaryHeap_TwoInts.cpp
 *
 *  Created on: Mar 19, 2016
 *      Author: michel
 */

#include "BinaryHeap_TwoInts.hpp"

bool compare2IntsFloat(ints2float x, ints2float y){
	return (x.key > y.key);
}

BinaryHeap_TwoInts::BinaryHeap_TwoInts( ) {
	this->heap = new std::vector< ints2float >();
	this->heap->reserve(this->DEFAULT_CAPACITY);
	std::make_heap (this->heap->begin(),this->heap->end(),compare2IntsFloat);
}

BinaryHeap_TwoInts::BinaryHeap_TwoInts(int maxCapacity ) {
	this->heap = new std::vector< ints2float >();
	this->heap->reserve(maxCapacity);
	std::make_heap (this->heap->begin(),this->heap->end(),compare2IntsFloat);
}
BinaryHeap_TwoInts::BinaryHeap_TwoInts(int* val1s, const int* val2s, Float keys) {
	this->heap = new std::vector< ints2float >();
	this->heap->reserve(this->DEFAULT_CAPACITY);
	std::make_heap (this->heap->begin(),this->heap->end(),compare2IntsFloat);
	for (int i=0;i<keys.length;i++){
		ints2float aux;
		aux.first =val1s[i];
		aux.second = val2s[i];
		aux.key = keys.pointer[i];
		this->heap->push_back(aux);
		std::push_heap (this->heap->begin(),this->heap->end(),compare2IntsFloat);
	}
}

BinaryHeap_TwoInts::BinaryHeap_TwoInts(int* val1s, const int* val2s, Float keys, int maxCapacity) {
	this->heap = new std::vector< ints2float >();
	this->heap->reserve(maxCapacity);
	std::make_heap (this->heap->begin(),this->heap->end(),compare2IntsFloat);
	for (int i=0;i<keys.length;i++){
		ints2float aux;
		aux.first =val1s[i];
		aux.second = val2s[i];
		aux.key = keys.pointer[i];
		this->heap->push_back(aux);
		std::push_heap (this->heap->begin(),this->heap->end(),compare2IntsFloat);
	}
}
void BinaryHeap_TwoInts::buildAgain(int* val1s, const int* val2s, Float keys){
	makeEmpty();
	for (int i=0;i<keys.length;i++){
		ints2float aux;
		aux.first =val1s[i];
		aux.second = val2s[i];
		aux.key = keys.pointer[i];
		this->heap->push_back(aux);
		std::push_heap (this->heap->begin(),this->heap->end(),compare2IntsFloat);
	}
}
BinaryHeap_TwoInts::~BinaryHeap_TwoInts(){
	makeEmpty();
	delete this->heap;
}


void BinaryHeap_TwoInts::insert(int val1, int val2, float key){
	ints2float aux;
	aux.first =val1;
	aux.second = val2;
	aux.key = key;
	this->heap->push_back(aux);
	std::push_heap (this->heap->begin(),this->heap->end(),compare2IntsFloat);

}

void BinaryHeap_TwoInts::deleteMin(){
	std::pop_heap(this->heap->begin(),this->heap->end(),compare2IntsFloat);
	this->heap->pop_back();
}


bool BinaryHeap_TwoInts::isEmpty(){
	return this->heap->empty();
}

int BinaryHeap_TwoInts::size(){
	return this->heap->size();
}

void BinaryHeap_TwoInts::makeEmpty(){
	if(this->heap == NULL) return;
	this->heap->clear();
}
void BinaryHeap_TwoInts::chopBottomK(int* val1s, int* val2s, Float keys){ //no idea if it works
	int k = keys.length;
	int pos;
	for (int i=0; i<k; i++) {
		pos = this->heap->size() - i -1;
		val1s[i] = this->heap->at(pos).first;
		val2s[i] = this->heap->at(pos).second;
		keys.pointer[i] = this->heap->at(pos).key;
	}
	this->heap->erase(this->heap->end()-k,this->heap->end());
	std::make_heap (this->heap->begin(),this->heap->end(),compare2IntsFloat);
}
int compareMyType3 (const void * a, const void * b)
{
  if ( *(int*)a <  *(int*)b ) return -1;
  if ( *(int*)a == *(int*)b ) return 0;
  return 1;
}
bool BinaryHeap_TwoInts::binHeapTwoTest(bool verbose){
	printf("Binary Heap Two Ints Test (Float Key)...\n");

	srand (time(NULL));
	int *int_aux = new int[10000];
	for(int l=1;l<6;l++){
		BinaryHeap_TwoInts *aux;
		for(int k=1;k<6;k++){

			if (verbose) printf("Default initializing...\n");
			aux = new BinaryHeap_TwoInts();


			if (verbose) printf("Inserting...\n");
			for(int i=0;i<10000;i++){
				int_aux[i] = rand()%100000;
				aux->insert(i*k*l,i*k*l*2,int_aux[i]);
			}
			qsort((void*)int_aux,10000,sizeof(int),compareMyType3);

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

			qsort((void*)int_aux,10000,sizeof(int),compareMyType3);

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

			if (verbose) printf("Binary Heap Two Ints Test (Float Key) Test %d completed successfully.\n",k);
		}
	}
	printf("Binary Heap Two Ints Test (Float Key) Test completed successfully.\n");
	return true;
}

