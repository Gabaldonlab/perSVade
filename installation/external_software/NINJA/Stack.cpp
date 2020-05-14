/*
 * Stack.cpp
 *
 *  Created on: Feb 13, 2016
 *      Author: michel
 */
#include "Stack.hpp"

Stack::Stack(){
	this->mystack = new std::stack<int>();
}

Stack::~Stack(){
	delete this->mystack;
}
void Stack::clear(){
	delete this->mystack;
	this->mystack = new std::stack<int>();

}
int Stack::length(){
	return this->mystack->size();
	this->mystack = new std::stack<int>();
}

void Stack::push(int x){
	this->mystack->push(x);
}
int Stack::pop(){
	int ret = this->mystack->top();
	this->mystack->pop();
	return ret;
}
bool Stack::isEmpty(){
	return this->mystack->empty();
}
//bool Stack::stackTest(bool verbose){
//	printf("Stack Test...\n");
//	if (verbose)
//	printf("Emptying...\n");
//	int *auxInt = new int[10000];
//	for(int l=1;l<6;l++){
//		Stack *aux;
//		for(int k=1;k<6;k++){
//			if (verbose)
//			printf("Default initializing...\n");
//			aux = new Stack();
//			if (verbose)
//			printf("Inserting...\n");
//			for(int i=0;i<10000;i++){
//				auxInt[i] = i+k+l;
//				aux->push(i+k+l);
//			}
//			if (verbose)
//			printf("Deleting and asserting...\n");
//			for(int i=0;i<5000;i++){
//				int x = aux->pop();
//				assert(x==auxInt[aux->length()]);
//			}
//			if (verbose)
//			printf("Assert size.\n");
//			assert(aux->length() == 5000);
//			if (verbose)
//			printf("Inserting...\n");
//			for(int i=0;i<4000;i++){
//				auxInt[aux->length()] = i+k+l;
//				aux->push(i+k+l);
//			}
//			if (verbose)
//			printf("Assert size.\n");
//			assert(aux->length() == 9000);
//			if (verbose)
//			printf("Deleting and asserting...\n");
//			for(int i=0;i<9000;i++){
//				int x = aux->pop();
//				assert(x==auxInt[aux->length()]);
//			}
//			if (verbose)
//			printf("Assert size.\n");
//			assert(aux->length() == 0);
//			assert(aux->isEmpty() == true);
//			if (verbose)
//			printf("Emptying...\n");
//			aux->clear();
//			if (verbose)
//			printf("Emptying again(check for double free)...\n");
//			aux->clear();
//			if (verbose)
//			printf("Stack Test %d completed successfully.\n",k*l);
//		}
//	}
//	delete[] auxInt;
//	printf("Stack Test completed successfully.\n");
//	return true;
//}
