
#include "BinaryHeap.hpp"

bool comp(std::pair<int,int> x, std::pair<int,int> y){
	return (x.second > y.second);
}
int compareMyType (const void * a, const void * b)
{
  if ( *(int*)a <  *(int*)b ) return -1;
  if ( *(int*)a == *(int*)b ) return 0;
  return 1;
}

BinaryHeap::BinaryHeap( ) {
	this->heap = new std::vector< std::pair<int,int> >();
	this->heap->reserve(this->DEFAULT_CAPACITY);
	std::make_heap (this->heap->begin(),this->heap->end(),comp);
}

BinaryHeap::BinaryHeap(int maxCapacity ) {
	this->heap = new std::vector< std::pair<int,int> >();
	this->heap->reserve(maxCapacity);
	std::make_heap (this->heap->begin(),this->heap->end(),comp);
}
BinaryHeap::BinaryHeap(const int* val1s, Int keys) {
	this->heap = new std::vector< std::pair<int,int> >();
	this->heap->reserve(this->DEFAULT_CAPACITY);
	std::make_heap (this->heap->begin(),this->heap->end(),comp);
	for (int i=0;i<keys.length;i++){
		this->heap->push_back(std::make_pair(val1s[i],keys.pointer[i]));
		std::push_heap (this->heap->begin(),this->heap->end(),comp);
	}
}

BinaryHeap::BinaryHeap(const int* val1s, Int keys, int maxCapacity) {
	this->heap = new std::vector< std::pair<int,int> >();
	this->heap->reserve(maxCapacity);
	std::make_heap (this->heap->begin(),this->heap->end(),comp);
	for (int i=0;i<keys.length;i++){
		this->heap->push_back(std::make_pair(val1s[i],keys.pointer[i]));
		std::push_heap (this->heap->begin(),this->heap->end(),comp);
	}
}

BinaryHeap::~BinaryHeap(){
	makeEmpty();
	delete this->heap;
}


int BinaryHeap::insert(int val1, int key){
	this->heap->push_back(std::make_pair(val1,key));
	std::push_heap (this->heap->begin(),this->heap->end(),comp);
	return 0;

}

void BinaryHeap::deleteMin(){
	std::pop_heap(this->heap->begin(),this->heap->end(),comp);
	this->heap->pop_back();
}


bool BinaryHeap::isEmpty(){
	return this->heap->empty();
}

int BinaryHeap::size(){
	return this->heap->size();
}

void BinaryHeap::makeEmpty(){
	if(this->heap == NULL) return;
	this->heap->clear();
}
bool BinaryHeap::binHeapTest(bool verbose){
	printf("Binary Heap Test...\n");

	srand (time(NULL));
	int *int_aux = new int[10000];
	for(int l=1;l<6;l++){
		BinaryHeap *aux;
		for(int k=1;k<6;k++){

			if (verbose) printf("Default initializing...\n");
			aux = new BinaryHeap();


			if (verbose) printf("Inserting...\n");
			for(int i=0;i<10000;i++){
				int_aux[i] = rand()%100000;
				aux->insert(i*k*l,int_aux[i]);
			}
			qsort((void*)int_aux,10000,sizeof(int),compareMyType);

			if (verbose) printf("Deleting and asserting...\n");
			for(int i=0;i<5000;i++){
				int test = aux->heap->front().second;
				int test2 = int_aux[i];
				assert(test2==test);
				aux->deleteMin();
			}

			if (verbose) printf("Assert size.\n");
			assert(aux->size() == 5000);

			if (verbose) printf("Inserting...\n");
			for(int i=0;i<5000;i++){
				int_aux[i] = rand()%100000;
				aux->insert(i*k*l,int_aux[i]);
			}

			if (verbose) printf("Assert size.\n");
			assert(aux->size() == 10000);

			qsort((void*)int_aux,10000,sizeof(int),compareMyType);

			if (verbose) printf("Deleting...\n");
			for(int i=0;i<10000;i++){
				std::pair<int,int> x = aux->heap->front();
				int test = x.second;
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

			if (verbose) printf("Binary Heap Test %d completed successfully.\n",k);
		}
	}
	printf("Binary Heap Test completed successfully.\n");
	return true;
}
