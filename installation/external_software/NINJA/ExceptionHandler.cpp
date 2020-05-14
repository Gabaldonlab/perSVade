/*
 * ExceptionHandler.cpp
 *
 *  Created on: Jan 24, 2016
 *      Author: michel
 */

#include "ExceptionHandler.hpp"

void Exception::critical(){
		fprintf(stderr,"Critical Error, aborting.\n");
		abort();
}
void Exception::criticalErrno(const char* arg){
		if (arg != NULL){
			perror(arg);
		}else{
			strerror(errno);
		}
		fprintf(stderr,"Critical Error, aborting.\n");
		abort();
}
