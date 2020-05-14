/*
 * Ninja.cpp
 *
 *  Created on: Feb 21, 2016
 *      Author: michel
 */
#include "TreeBuilderManager.hpp"
#include "ClusterManager.hpp"
#include "ArgumentHandler.hpp"
#include "BinaryHeap.hpp"
#include <iostream>


#ifdef LINUX //TODO: include support for multiple plataforms
#include <unistd.h>
#endif

#ifdef TESTCASE

//main unit test
void test_case(int argc, char *argv[]){
	#ifdef LINUX
	long maxMemory = sysconf(_SC_PAGE_SIZE)*sysconf(_SC_AVPHYS_PAGES);
	std::cout<<maxMemory<<std::endl;
	#endif

	ArgumentHandler *argHandler = new ArgumentHandler(argv,argc);
	if (!argHandler->argumentTest())
		printf("Argument test failed");
	Stack *stack = new Stack();
	if (!stack->stackTest(false))
		printf("Stack test failed");
	BinaryHeap* binHeap = new BinaryHeap();
	if (!binHeap->binHeapTest(false))
		printf("Binary Heap test failed");
	BinaryHeap_IntKey_TwoInts* binHeapTwo = new BinaryHeap_IntKey_TwoInts();
	if (!binHeapTwo->binHeapTwoTest(false))
		printf("Binary Heap Two Ints test failed");

/*	BinaryHeap_FourInts* binHeapTwo = new BinaryHeap_FourInts();
	if (!binHeapTwo->binHeapFourTest(true))
		printf("Binary Heap Four Ints test failed");*/

	ArrayHeapExtMem* arrayHeap = new ArrayHeapExtMem("",NULL);

	if (!arrayHeap->test(true))
		printf("Array Heap External Memorytest failed");

}

#endif

int main(int argc, char *argv[]){
#ifdef TEST_CASE
	test_case(argc, argv);
#else

	ArgumentHandler *argHandler = new ArgumentHandler(argv,argc);
	if (argHandler->abort)
		return 0;

	std::string method = argHandler->getMethod();
	std::string inFile = argHandler->getInFile();
	std::string njTmpDir = argHandler->getNJTmpDir();
	ArgumentHandler::InputType inType = argHandler->getInputType();
	ArgumentHandler::OutputType outType = argHandler->getOutputType();
	ArgumentHandler::AlphabetType alphType = argHandler->getAlphabetType();
	ArgumentHandler::CorrectionType corrType = argHandler->getCorrectionType();
	float clusterCutoff = argHandler->getClusterCutoff();
	FILE* out = argHandler->getOutpuFile();
	int threads = argHandler->getNumThreads();
	bool useSSE = argHandler->useSSE();
        bool printTime = argHandler->getPrintTime();

	fprintf(stderr,"Reading...\n");

        if ( outType == ArgumentHandler::cluster ) {
          // RMH: Cluster sequences by nearest neighbor rather than build
          //      a tree.
	  ClusterManager* manager = new ClusterManager(method, njTmpDir, inFile, out, (ClusterManager::InputType)
                                           inType,(ClusterManager::OutputType) outType, (ClusterManager::AlphabetType)
                                           alphType,(ClusterManager::CorrectionType) corrType, threads, useSSE, printTime, 
                                           clusterCutoff);

	  std::string clusterString = manager->doJob();

        }else {
	  TreeBuilderManager* manager = new TreeBuilderManager(method, njTmpDir, inFile, out, (TreeBuilderManager::InputType)
                                                inType,(TreeBuilderManager::OutputType) outType, (TreeBuilderManager::AlphabetType)
                                                alphType,(TreeBuilderManager::CorrectionType) corrType, threads, useSSE, printTime);

	  std::string treeString = manager->doJob();

  	  if(outType == ArgumentHandler::dist){
		fprintf(stderr,"Distances successfully generated.\n");
		fclose(out);
		return 0;
	  }
	  //printf("String: %s \n",treeString.c_str());
	  //return 0;
	  if (!treeString.empty()) {
  		fprintf(out,"%s",treeString.c_str() );
		fprintf(stderr,"Tree successfully generated.\n");
	  } else {
		fprintf(stderr,"\n\nTree string not generated for some unknown reason. Aborting.");
		Exception::critical();
	  }

	  fclose(out);
  	  return 0;
        }
#endif
}

