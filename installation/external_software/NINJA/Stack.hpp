/*
 * Stack.hpp
 *
 *  Created on: Feb 13, 2016
 *      Author: michel
 */

#ifndef STACK_HPP
#define STACK_HPP

#include "ExceptionHandler.hpp"
#include <string.h>
#include <stack>

class Stack{
	public:
		Stack();
		~Stack();
		Stack(int size);
		void push(int x);
		int pop();
		int length();
		int maxSize();
		bool isEmpty();
		void clear();

		bool stackTest(bool verbose);
	private:
		std::stack<int> *mystack;

		void resize();
};
#endif

