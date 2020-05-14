/*
 * ArrayHeapExtMem.cpp
 *
 *  Created on: Mar 19, 2016
 *      Author: michel
 */

#include "ArrayHeapExtMem.hpp"

#define LINUX 1

#ifdef LINUX
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

ArrayHeapExtMem::ArrayHeapExtMem(std::string dir, int* activeIJs){
	initialize ( dir, activeIJs, (long)pow(2,21));
}
ArrayHeapExtMem::ArrayHeapExtMem(std::string dir, int* activeIJs, long sizeExp){
	initialize ( dir, activeIJs, sizeExp);
}
void ArrayHeapExtMem::deleteAll(){ // same as the ~ArrayHeapExtMem(), but it does not destroy the object.
	this->n = 0;
	if (this->H1 != NULL) this->H1->makeEmpty();
	if (this->H2 != NULL) this->H2->makeEmpty();

	if (this->freeSlots != NULL){
		freeSlots->clear();
	}

	if (this->perSlotIntBuffer!= NULL){
		for(int i=0;i<this->maxLevels;i++){
			if (this->perSlotIntBuffer[i] != NULL){
				for(int j=0;j<this->numSlots;j++){
					if (this->perSlotIntBuffer[i][j] != NULL)
						delete[] this->perSlotIntBuffer[i][j];
				}
				delete[] this->perSlotIntBuffer[i];
			}
		}
		delete[] this->perSlotIntBuffer;
		this->perSlotIntBuffer = NULL;
	}
	if (bigBuffI != NULL){
		delete[] this->bigBuffI;
		this->bigBuffI = NULL;
	}
	if (bigBuffB != NULL){
		delete[] this->bigBuffB;
		this->bigBuffB = NULL;
	}

	if (cntMax != NULL){
		delete[] this->cntMax;
		this->cntMax = NULL;
	}
	if (buffI != NULL){
		delete[] this->buffI;
		this->buffI = NULL;
	}
	if (buffB != NULL){
		delete[] this->buffB;
		this->buffB = NULL;
	}
	if (this->slotNodeCnt != NULL){
		for(int i=0;i<this->maxLevels;i++)
			if (this->slotNodeCnt[i] != NULL)
				delete[] this->slotNodeCnt[i];
		delete[] this->slotNodeCnt;
		this->slotNodeCnt = NULL;
	}
	if (this->cntOnHeap != NULL){
		for(int i=0;i<this->maxLevels;i++)
			if (this->cntOnHeap[i] != NULL)
				delete[] this->cntOnHeap[i];
		delete[] this->cntOnHeap;
		this->cntOnHeap = NULL;
	}
	if (this->slotPositions != NULL){
		for(int i=0;i<this->maxLevels;i++)
			if (this->slotPositions[i] != NULL)
				delete[] this->slotPositions[i];
		delete[] this->slotPositions;
		this->slotPositions = NULL;
	}
	if (this->slotBuffPos != NULL){
		for(int i=0;i<this->maxLevels;i++)
			if (this->slotBuffPos[i] != NULL)
				delete[] this->slotBuffPos[i];
		delete[] this->slotBuffPos;
		this->slotBuffPos = NULL;
	}

	if (this->file!=NULL){
		//fprintf(stderr,"File deleted: ArrayHeapExtMem delete.\n");
		fclose(this->file);
		this->file = NULL;
	}
	remove(this->fileName.c_str());
}
ArrayHeapExtMem::~ArrayHeapExtMem(){
	deleteAll();
}

void ArrayHeapExtMem::initialize(std::string dir, int* activeIJs, long sizeExp){

	this->H1 = NULL;
	this->H2 = NULL;
	this->slotNodeCnt = NULL;
	this->perSlotIntBuffer = NULL;
	this->slotBuffPos = NULL;
	this->cntOnHeap = NULL;
	this->slotPositions = NULL;
	this->buffB = NULL;
	this->buffI = NULL;
	this->bigBuffB = NULL;
	this->bigBuffI = NULL;
	this->freeSlots = NULL;
	this->file = NULL;
	this->sortingHeap = NULL;

	this->slotPairList = new std::vector<SlotPair*>();
	this->freeSlots = new std::vector<std::list<int>*>();

	this->A = -1;
	this->B = -1;

	this->maxLevels = 4;
	this->blockSize = 1024 ; //4KB blocks = ~1000 ints
	this->n=0;
	this->loopCntr = 0;
	this->startCheckingValues = false;

	this->numFields = 3;
	this->cntMax = new long[maxLevels];
	this->c = (float)1/85;

	this->active = activeIJs;

	char num[15];

	sprintf(num, "%d",rand());
	if (dir == "")
		this->tmpDir = "tmp/";
	else
		this->tmpDir = dir;

    mkdir(this->tmpDir.c_str(), 0700);

	this->fileName = this->tmpDir + "arrayHeap" + num;

	if (sizeExp >  pow(2, 22)  /*4MB*/) {
		this->blockSize = 2048;
		if (sizeExp >  pow(2, 23)  /*8MB*/) {
			this->blockSize = 4096;
		}
	}
	this->mem = sizeExp;

	this->cM = (int)(this->c*this->mem);
	this->numSlots = (int)((this->cM/this->blockSize)-1);
	this->numNodesPerBlock = this->blockSize/this->numFields;
	this->numFieldsPerBlock = this->numNodesPerBlock*this->numFields; // note, this might not be blockSize, because of rounding.

	prepare();
}
void ArrayHeapExtMem::prepare(){
	clearAndInitialize();

	this->cntMax = new long[maxLevels];

	if (this->H1 == NULL) this->H1 = new BinaryHeap_TwoInts(this->cM*2);
	if (this->H2 == NULL) this->H2 = new BinaryHeap_FourInts();
	if (this->perSlotIntBuffer== NULL){
		this->perSlotIntBuffer = new int**[this->maxLevels]();
		for(int i=0;i<this->maxLevels;i++){
			this->perSlotIntBuffer[i] = new int*[this->numSlots]();
			for(int j=0;j<this->numSlots;j++)
				this->perSlotIntBuffer[i][j] = new int[this->numFieldsPerBlock]();
		}
	}
	if (this->slotNodeCnt == NULL){
		this->slotNodeCnt = new long*[this->maxLevels];
		for(int i=0;i<this->maxLevels;i++)
			this->slotNodeCnt[i] = new long[this->numSlots];
	}
	if (this->cntOnHeap == NULL){
		this->cntOnHeap = new int*[this->maxLevels]();
		for(int i=0;i<this->maxLevels;i++)
			this->cntOnHeap[i] = new int[this->numSlots]();
	}
	if (this->slotPositions == NULL){
		this->slotPositions = new long*[this->maxLevels];
		for(int i=0;i<this->maxLevels;i++)
			this->slotPositions[i] = new long[this->numSlots];
	}
	if (this->slotBuffPos== NULL){
		this->slotBuffPos = new int*[this->maxLevels];
		for(int i=0;i<this->maxLevels;i++)
			this->slotBuffPos[i] = new int[this->numSlots];
	}

	this->tempFile = fopen(this->fileName.c_str(), "w+");
	if(this->tempFile == NULL) Exception::criticalErrno(this->fileName.c_str());
	if (this->file == NULL) this->file = this->tempFile;
	fclose(this->file);
	this->file = NULL;
	this->tempFile = NULL;


	this->buffI = new int[this->numFieldsPerBlock];
	this->buffB = new char[this->numFieldsPerBlock*4];
	this->bigBuffI = new int[10*this->numFieldsPerBlock];
	this->bigBuffB = new char[10*this->numFieldsPerBlock*4];


	this->cntMax[0] = (long)this->cM;
	for (int level=1; level<this->maxLevels; level++)  {
		this->cntMax[level] = this->cntMax[level-1] * (this->numSlots + 1);

	}
}
void ArrayHeapExtMem::clearAndInitialize(){

	this->n = 0;
	if (this->H1 != NULL) this->H1->makeEmpty();
	if (this->H2 != NULL) this->H2->makeEmpty();

	if (this->cntOnHeap != NULL) {
		for (int i=0; i<this->maxLevels; i++) {
			for (int j=0; j<this->numSlots; j++) {
				this->cntOnHeap[i][j] = 0;
				this->slotNodeCnt[i][j] = 0;
				this->slotPositions[i][j] = 0;
			}
		}
	}

	if (this->freeSlots == NULL){
		this->freeSlots = new std::vector<std::list<int>*>();
		this->freeSlots->reserve(this->maxLevels);
	}else{
		delete freeSlots;
		this->freeSlots = new std::vector<std::list<int>*>();
	}

	for (int i=0; i<this->maxLevels; i++) {
		std::list<int> *l = new std::list<int>();
		for (int j=0; j<this->numSlots; j++) {
			l->push_back(j);
		}
		this->freeSlots->push_back(l);
	}

	if(this->bigBuffB == NULL){
		delete[] this->bigBuffB;
		this->bigBuffB = NULL;
	}

	if(this->bigBuffI == NULL){
		delete[] this->bigBuffI;
		this->bigBuffI = NULL;
	}


}
void ArrayHeapExtMem::insert(int i, int j, float key){
	//too many problems
	this->H1->insert(i, j, key);
	this->n++;


	if (this->H1->size() == 2* this->cM) {
		//pull half the nodes off heap, put into a new Node array, then put into L1 (or push down to L2 or L3)

		//I implement this by just removing the last half of the binary heap.
		//This seems worthwhile, as these will generally be lower valued entries, which
		//are the ones we'd prefer to see dumped to file anyway.


		int* is = new int[this->cM]();
		int* js = new int[this->cM]();
		Float keys;
		keys.pointer = new float[this->cM]();
		keys.length = this->cM;

		if (this->chopBottom) {
			/*new way - pulling from the bottom of the heap*/


			this->H1->chopBottomK(is, js, keys);


			//need to sort them. Use Heap sort
			if (this->sortingHeap == NULL)
				this->sortingHeap = new BinaryHeap_TwoInts(is, js, keys);
			else
				this->sortingHeap->buildAgain(is, js, keys);


			int K = this->cM;

			for (int ii=0; ii<K; ii++)  {

				is[ii] = this->sortingHeap->heap->front().first;
				js[ii] = this->sortingHeap->heap->front().second;
				keys.pointer[ii] = this->sortingHeap->heap->front().key;

				this->sortingHeap->deleteMin();

			}

			this->sortingHeap = NULL;
		} else {

			/* old way - pulling from top of heap*/

			for (int x=0; x<this->cM; x++) {

				is[x] = this->H1->heap->front().first;
				js[x] = this->H1->heap->front().second;
				keys.pointer[x] = this->H1->heap->front().key;

				this->H1->deleteMin();

			}
		}

		//if there are no free slots at level 0, figure out which level has a free slot
		int targetLevel=0;
		while (this->freeSlots->at(targetLevel)->size() == 0 ) {
			if (mergeSlots(targetLevel) )
				break;  // if two half-empty slots can be merged, then there's now room
			targetLevel++;
		}

		if (targetLevel == 0) {
			int newSlot = store(0, is, js, keys.pointer, keys.length);
			load(0, newSlot);
		} else {
			int slot = mergeLevels (targetLevel,  is, js, keys );

			load(targetLevel, slot);

		}
		delete[] js;
		delete[] is;

	}
}
long ArrayHeapExtMem::calcPos (int level, int slot) {
	long pos = 0;
	for (int i=0; i<level; i++) {
		pos += numSlots * cntMax[i] * numFields * 4;
	}
	pos += slot * cntMax[level] * numFields * 4;  // this is a conservative placement of the new slot, since cntMax is already conservative.

	return pos;
}


int ArrayHeapExtMem::mergeLevels (int targetLevel, int* is, int* js, Float keys){
	if (targetLevel  > maxLevels-1) {
		fprintf(stderr,"unexpected occurance: External Memory Array heap tried to write to a level > %d \n",maxLevels);
		Exception::critical();
	}


	//open file for use
	if (this->tempFile == NULL) this->tempFile = fopen(this->fileName.c_str(), "r+");
	if (this->tempFile == NULL) Exception::criticalErrno(this->fileName.c_str());
	if (this->file == NULL) this->file = this->tempFile;


	int targetSlot = freeSlots->at(targetLevel)->front();
	freeSlots->at(targetLevel)->pop_front();

	long outFilePos = calcPos (targetLevel, targetSlot);


	float** tmpKeys = new float*[targetLevel];
	for(int i=0;i<targetLevel;i++)
		tmpKeys[i] = new float[numSlots];
	int* deletedLevels = new int[targetLevel];
	for (int level=0; level<targetLevel; level++) {
		deletedLevels[level] = level;

		for (int slot=0; slot<numSlots; slot++){
			tmpKeys[level][slot] = FLT_MIN;
			slotPositions[level][slot] -= cntOnHeap[level][slot]; // move back the pointers on the L_i arrays to include positions that are still on the heap
			slotBuffPos[level][slot] = numFieldsPerBlock;
		}
	}
	deleteLevels(this->H2, deletedLevels,targetLevel);


	long basePos=-1;

	int minLevel = -1;
	int minSlot = -1;
	float min;
	float key;
	int inputPos=0;
	long diskPos=-1;
	int buffPos;


	int outBuffPos = 0;

	int nonClearedSlots = targetLevel * numSlots ;

	// we are merging all the levels below targetLevel, plus the submitted values,
	// and putting them in an available slot in targetLevel

	int newCnt = 0;
	bool stillNeedNode;
	int* intBuffer;


	while (nonClearedSlots > 0 || inputPos < keys.length) {


		if (inputPos < keys.length) {
			minLevel = -1;
			min = keys.pointer[inputPos];
		} else {
			min = FLT_MAX;
		}


		for (int level=0; level<targetLevel; level++) {

			for (int slot=0; slot<numSlots; slot++){

				stillNeedNode = true;


				while (stillNeedNode) {

					if (slotPositions[level][slot] == slotNodeCnt[level][slot]) {
						stillNeedNode = false;
						continue;
					}

					if (tmpKeys[level][slot] == FLT_MIN) { // pointer for this slot didn't move; use this so I don't have to do intBitsToFloat calculation for every entry
						// read next value from file buffer ... may need to refill buffer


						if (slotBuffPos[level][slot] == numFieldsPerBlock) { // refill buffer
							basePos = calcPos (level, slot);
							diskPos = (long)numFields * 4 * slotPositions[level][slot]; // 12 bytes for each position: int,int,float
							fseek(file,diskPos + basePos,SEEK_SET);

							fread(this->perSlotIntBuffer[level][slot],sizeof(int),this->numFieldsPerBlock,this->file);

							slotBuffPos[level][slot] = 0;
						}
						intBuffer = perSlotIntBuffer[level][slot];
						buffPos = slotBuffPos[level][slot];


						try {
							if (trimInactiveFromHeapArray && active != NULL && (-1 == active[intBuffer[buffPos]] /*i*/ || -1 == active[intBuffer[buffPos+1]]/*j*/ )) {
								slotBuffPos[level][slot] += numFields;
								slotPositions[level][slot]++;
								n--;
								if (slotPositions[level][slot] == slotNodeCnt[level][slot])
									nonClearedSlots--;

								continue; // don't bother copying dead nodes to the merged slot
							}
						} catch (int e) {
							Exception::critical();

/*							LogWriter.stdErrLogln("level=" + level + ", slot=" + slot + ", targetlevel=" + targetLevel + ", slotPos=" +
									slotPositions[level][slot] + ", slotNodeCnt=" + slotNodeCnt[level][slot] + ", buffPos=" +
									buffPos + ", a=" + intBuffer[buffPos] + ", b=" + intBuffer[buffPos+1] + ", c=" + intBuffer[buffPos+2] +
									", active.len=" + active.length);

							LogWriter.stdErrLogln("And for good measure, the next slot position:  buffPos=" +
									(buffPos+3) + ", a=" + intBuffer[buffPos+3] + ", b=" + intBuffer[buffPos+4] + ", c=" + intBuffer[buffPos+5] +
									", active.len=" + active.length);*/
							throw e;
						}
						memcpy(&tmpKeys[level][slot],&intBuffer[buffPos+2],sizeof(float));
						//tmpKeys[level][slot] = (float)intBuffer[buffPos+2];
					}

					stillNeedNode = false;
					key = tmpKeys[level][slot];
					if ( key < min) {
						minLevel = level;
						minSlot = slot;
						min = key;
					}
				}
			}
		}



		if (minLevel == -1) { // read from the arrays pulled from H1
			bigBuffI[outBuffPos++] = is[inputPos];
			bigBuffI[outBuffPos++] = js[inputPos];
			memcpy(&bigBuffI[outBuffPos++],&keys.pointer[inputPos],sizeof(int));

			inputPos++;
		} else { // read from the slots being merged
			buffPos = slotBuffPos[minLevel][minSlot];
			intBuffer = perSlotIntBuffer[minLevel][minSlot];

			bigBuffI[outBuffPos++] = intBuffer[buffPos];
			bigBuffI[outBuffPos++] = intBuffer[buffPos+1];
			bigBuffI[outBuffPos++] = intBuffer[buffPos+2];


			slotPositions[minLevel][minSlot]++;

			if (slotPositions[minLevel][minSlot] == slotNodeCnt[minLevel][minSlot])
				nonClearedSlots--;

			slotBuffPos[minLevel][minSlot] += numFields;
			tmpKeys[minLevel][minSlot] = FLT_MIN;
		}

		newCnt++;


		if (outBuffPos == 10*this->numFieldsPerBlock) { //changed

			//need to write buffer to file if it's full

			fseek(file,outFilePos,SEEK_SET);

			fwrite(bigBuffI,sizeof(int),10*this->numFieldsPerBlock,this->file);

			outFilePos += 10*this->numFieldsPerBlock*4;
			outBuffPos = 0;

		}

	}

	if (outBuffPos > 0) {
		//need to write buffer to file

		fseek(file,outFilePos,SEEK_SET);
		fwrite(bigBuffI, sizeof(int), outBuffPos,this->file);
	}

	//clean up old slots
	for (int level=0; level<targetLevel; level++) {
		std::list<int> *l = freeSlots->at(level);
		l->clear();

		for (int slot=0; slot<numSlots; slot++) {
			l->push_back(slot);
			cntOnHeap[level][slot] = 0;
			slotNodeCnt[level][slot] = 0;
			slotPositions[level][slot] = 0;
		}
	}

	slotNodeCnt[targetLevel][targetSlot] = newCnt;
	cntOnHeap[targetLevel][targetSlot] = 0;
	slotPositions[targetLevel][targetSlot] = 0;

	//close file
	if (this->file!=NULL){
		//fprintf(stderr,"File deleted: ArrayHeapExtMem delete.\n");
		fclose(this->file);
		this->tempFile = NULL;
		this->file = NULL;
	}

	return targetSlot;
}
bool ArrayHeapExtMem::mergeSlots (int lvl){

	//   If it's levels 1-3, just store the slots in memory and write 'em all out (as below)
	//    but ... if it's lvl 4, then use a temporary file (or just a temporary space at the end of the file),
	//    read blocks and write them to that file/slot,
	//    and when it's all done, reset the file handles to use that as one of the merged slots.
	  //  ... I haven't implemented the later part yet
	if (lvl == maxLevels-1) {
/*		LogWriter.stdErrLogln("Request was made to merge two slots from level " + maxLevels + ". " +
				"  That's not implemented yet!");
		throw new Exception("Request was made to merge two slots from level " + maxLevels + ". " +
				"  That's not implemented yet!");*/
		Exception::critical();
		//return false;
	}


	//scan through the slots at the given level, and see if two of them are less than half full.
	// The Ferragina paper provides for doing all this at deletion time, which allows a
	//  guarantee that no more than one slot in a level will be less than half full -  I delay
	//  it until an insert, so I don't bother merging slots that might be emptied before an insert.  This
	// requires a small amount of overhead to find such slots, but it only happens when a store
	// command is being called anyway, so the overhead is effectively irrelevant.
	long remaining;


	// because this implementation allows removal of expired nodes when slots are merged, the result of a merge
	//might already be small.  No use immediately turning around to merge that one ... so here I merge
	// multiple slots - as many as will fit in one new slot.  Pick slots to merge greedily
	slotPairList->clear();
	for (int slot=0; slot<numSlots; slot ++) {
		if (cntOnHeap[lvl][slot]>0) {
			remaining = slotNodeCnt[lvl][slot] - slotPositions[lvl][slot] + cntOnHeap[lvl][slot];
			if ( remaining <= cntMax[lvl]/2 ) {
				slotPairList->push_back(new SlotPair(slot, remaining));
			}
		}
	}

	if (slotPairList->size()<2) return false;

	//LogWriter.stdErrLogln("merging " + slotPairList.size() + " slots in level " + lvl);

	qsort((void*)(&(*slotPairList)[0]),sizeof(SlotPair*),slotPairList->size(),SlotPair::compareTo);
	//Collections.sort(slotPairList);

	int summedSize = (int)(slotPairList->at(0)->remaining + slotPairList->at(1)->remaining);
	int numSlotsToMerge = 2;
	while (numSlotsToMerge < (int)slotPairList->size() && summedSize + slotPairList->at(numSlotsToMerge)->remaining <= cntMax[lvl]) {
		summedSize += slotPairList->at(numSlotsToMerge++)->remaining;
	}


	int nonClearedSlots = numSlotsToMerge ;
	int* slots = new int[numSlotsToMerge]();
	for (int i=0; i<numSlotsToMerge; i++) {
		slots[i] = slotPairList->at(i)->slot;
	}

	float* tmpKeys = new float[numSlots]();


	for (int i=0;i<numSlotsToMerge;i++){
		tmpKeys[slots[i]] = FLT_MIN;
		slotPositions[lvl][slots[i]] -= cntOnHeap[lvl][slots[i]]; // move back the pointers on the L_i arrays to include positions that are still on the heap
		slotBuffPos[lvl][slots[i]] = numFieldsPerBlock;

	}


	int* is = new int[summedSize]();
	int* js = new int[summedSize]();
	float* keys = new float[summedSize]();




	int newCnt = 0;
	bool stillNeedNode;
	int* intBuffer;
	int minSlot = -1;
	float min;
	float key;
	long basePos, diskPos;
	int buffPos;
	int outBuffPos = 0;

	while (nonClearedSlots > 0 ) {
		min = FLT_MAX;

		for (int i=0;i<numSlotsToMerge;i++){

			stillNeedNode = true;

			while (stillNeedNode) {

				if (slotPositions[lvl][slots[i]] == slotNodeCnt[lvl][slots[i]]) {
					stillNeedNode = false;
					continue;
				}

				if (tmpKeys[slots[i]] == FLT_MIN) { // pointer for this slot moved; use this so I don't have to do intBitsToFloat calculation for every entry
					// read next value from file buffer ... may need to refill buffer
					if (slotBuffPos[lvl][slots[i]] == numFieldsPerBlock) { // refill buffer
						basePos = calcPos (lvl, slots[i]);
						diskPos = (long)numFields * 4 * slotPositions[lvl][slots[i]]; // 12 bytes for each position: int,int,float
						fseek(file,diskPos + basePos,SEEK_SET);
						fread(perSlotIntBuffer[lvl][slots[i]],sizeof(int),this->numFieldsPerBlock,file);
						//Arrays.byteToInt(buffB, perSlotIntBuffer[lvl][slots[i]]);
						slotBuffPos[lvl][slots[i]] = 0;
					}
					intBuffer = perSlotIntBuffer[lvl][slots[i]];
					buffPos = slotBuffPos[lvl][slots[i]];
					if (trimInactiveFromHeapArray && active != NULL && (-1 == active[intBuffer[buffPos]] /*i*/ || -1 == active[intBuffer[buffPos+1]]/*j*/ )) {
						slotBuffPos[lvl][slots[i]] += numFields;
						slotPositions[lvl][slots[i]]++;
					//	skippedCnt++;
						n--;
						if (slotPositions[lvl][slots[i]] == slotNodeCnt[lvl][slots[i]])
							nonClearedSlots--;

						continue; // don't bother copying dead nodes to the merged slot
					}
					memcpy(&tmpKeys[slots[i]],&intBuffer[buffPos+2],sizeof(float));
					//tmpKeys[slots[i]] = (float)(intBuffer[buffPos+2]);
				}

				stillNeedNode = false;
				key = tmpKeys[slots[i]];
				if ( key < min) {
					minSlot = slots[i];
					min = key;
				}
			}
		}

		if (minSlot<0) {

/*			for (int i=0;i<numSlotsToMerge;i++){
				if ( slotPositions[lvl][slots[i]] != slotNodeCnt[lvl][slots[i]] ) { // otherwise, it's fine that we ended up here
					LogWriter.stdErrLogln("Surprising minSlot == -1.  Is it just 'cause we're at the end of the list, and they're all expired?");
					LogWriter.stdErrLog("\t");
					for (int slot2 : slots){
						LogWriter.stdErrLog("[ " + slot2 + " @ " + slotPositions[lvl][slot2] + " of " + slotNodeCnt[lvl][slot2] + "] , ");
					}
					LogWriter.stdErrLogln("");

					throw new Exception ("Surprising minSlot == -1.");
				}
			}*/
			if (nonClearedSlots > 0) {
				//LogWriter.stdErrLogln("Strange.  nonClearedSlots should be 0 but isn't.  Making it 0 now. A=" + A + ", B=" + B );
				nonClearedSlots = 0;
			}

		} else {
			buffPos = slotBuffPos[lvl][minSlot];
			intBuffer = perSlotIntBuffer[lvl][minSlot];

			is[outBuffPos] = intBuffer[buffPos];
			js[outBuffPos] = intBuffer[buffPos+1];
			keys[outBuffPos] = intBuffer[buffPos+2];
			outBuffPos++;

			slotPositions[lvl][minSlot]++;

			if (slotPositions[lvl][minSlot] == slotNodeCnt[lvl][minSlot])
				nonClearedSlots--;

			slotBuffPos[lvl][minSlot] += numFields;
			tmpKeys[minSlot] = FLT_MIN; //

			newCnt++;
		}


	}




	//remove the defunct entries from H2 (the appropriate number will get added back in a moment)
	deleteLevelAndSlots(H2, lvl, slots,numSlotsToMerge);

	for (int i=0;i<numSlotsToMerge;i++){
		freeSlots->at(lvl)->push_back(slots[i]);
		cntOnHeap[lvl][slots[i]] = 0;
		slotPositions[lvl][slots[i]] = 0;
		slotNodeCnt[lvl][slots[i]] = 0;
	}


	int loc = store(lvl, is, js, keys, newCnt);
	load(lvl, loc);


	return true;

}
void ArrayHeapExtMem::load (int level, int slot){

	long cntInSlot = slotNodeCnt[level][slot];
	if (slotPositions[level][slot] == cntInSlot) {
		slotNodeCnt[level][slot] = 0;
		slotPositions[level][slot] = 0;
		freeSlots->at(level)->push_back(slot);
		cntOnHeap[level][slot] = 0;

		return;
	}

	//open file for use
	if (this->tempFile == NULL) this->tempFile = fopen(this->fileName.c_str(), "r+");
	if (this->tempFile == NULL) Exception::criticalErrno(this->fileName.c_str());
	if (this->file == NULL) this->file = this->tempFile;


	//long basePos = cntMax[level] * numFields * 4 * slot;
	long basePos = calcPos (level, slot);


	long slotPos = 4 * numFields * slotPositions[level][slot]; // 12 bytes for each position: int,int,float

	fseek(file,slotPos + basePos,SEEK_SET);
	int sizeRead = fread(buffI,sizeof(int),this->numFieldsPerBlock,file);
	//Arrays.byteToInt(buffB, buffI);
	if (sizeRead == 0){
		Exception::criticalErrno(this->fileName.c_str());
	}

	long endPos = slotPositions[level][slot] + numNodesPerBlock;
	if (endPos > cntInSlot) endPos = cntInSlot;
	cntOnHeap[level][slot] = (int)(endPos - slotPositions[level][slot]);


	int i=0;
	float key = 0.f;
	int val1 = 0, val2 = 0;
	while ( slotPositions[level][slot] < endPos )  {

		val1 = buffI[i++];
		val2 = buffI[i++];
		memcpy(&key,&buffI[i++],sizeof(float));
		//key = (float)(buffI[i++]);



		H2->insert(val1, val2, level, slot, key); // "key" holds the "d" value for the pair. That's what I'll put in "fl" for now
		slotPositions[level][slot]++;
	}

	//close file
	if (this->file!=NULL){
		//fprintf(stderr,"File deleted: ArrayHeapExtMem delete.\n");
		fclose(this->file);
		this->tempFile = NULL;
		this->file = NULL;
	}

}

int ArrayHeapExtMem::store(int level, int* is, int* js, float* keys, int cnt){

	//open file for use
	if (this->tempFile == NULL) this->tempFile = fopen(this->fileName.c_str(), "r+");
	if (this->tempFile == NULL) Exception::criticalErrno(this->fileName.c_str());
	if (this->file == NULL) this->file = this->tempFile;

	int freeSlot = freeSlots->at(level)->front();

	freeSlots->at(level)->pop_front();

	int* bI = new int[cnt * numFields]();

	//convert to byte array
	int i=0;
	for (int j=0; j<cnt; j++) {

		bI[i++] = is[j];
		bI[i++] = js[j];
		memcpy(&bI[i++],&keys[j],sizeof(int));

	}
	//Arrays.intToByte(bI, bB);

	long basePos = calcPos (level, freeSlot);

	fseek(file,basePos,SEEK_SET);
	int sizeWritten = fwrite(bI,sizeof(int),cnt * numFields,file);
	if(sizeWritten==0) Exception::criticalErrno(this->fileName.c_str());

	delete[] bI;

	slotPositions[level][freeSlot] = 0;
	cntOnHeap[level][freeSlot] = 0;

	slotNodeCnt[level][freeSlot] = cnt;

	//close file
	if (this->file!=NULL){
		//fprintf(stderr,"File deleted: ArrayHeapExtMem delete.\n");
		fclose(this->file);
		this->tempFile = NULL;
		this->file = NULL;
	}

	return freeSlot;
}
void ArrayHeapExtMem::deleteLevelAndSlots (BinaryHeap_FourInts* H, int level, int* slots, int numSlotsToMerge) {

	bool ok = false;
	for (int i=0;i<numSlotsToMerge;i++){
		if (cntOnHeap[level][slots[i]] > 0) {
			cntOnHeap[level][slots[i]] = 0;
			ok = true;
		}
	}

	if (!ok) {
		//nothing to do
		return;
	}


	// I can re-build the heap without the to-be-removed entries in linear time,
	// and since it'll take linear time just to identify to to-be-removed entries
	// ... that's what I'll do.


	// get list of nodes without this level and slot
	std::list<int> *list = new std::list<int>();
	for (int i=0; i<H->size(); i++){
		ok = true;
		if (H->heap->at(i).third == level) { // otherwise, it's certainly ok
			for (int j=0;j<numSlotsToMerge;j++)
				if (H->heap->at(i).fourth == slots[j])  ok = false;
		}
		if (ok) list->push_back(i);
	}
	int K = list->size();


	if (K == 0) {
		H->makeEmpty();
	} else if (K < H->size()) { // it's possible that there are no entries to delete ... so no reason to rebuild

		//get the values for those nodes, so we can create a new heap
		int *v1 = new int[K];
		int *v2 = new int[K];
		int *v3 = new int[K];
		int *v4 = new int[K];
		float *key =  new float[K];
		int pos;
		int i=0;
		for (std::list<int>::iterator it=list->begin(); it != list->end(); ++it){
				pos = *it;
				v1[i] = H->heap->at(pos).first;
				v2[i] = H->heap->at(pos).second;
				v3[i] = H->heap->at(pos).third;
				v4[i] = H->heap->at(pos).fourth;
				key[i] = H->heap->at(pos).key;
				i++;
		}
		Float keys;
		keys.pointer = key;
		keys.length = K;
		H->buildAgain(v1, v2, v3, v4, keys);
		delete[] v1;
		delete[] v2;
		delete[] v3;
		delete[] v4;
		delete[] key;

	}
}
void ArrayHeapExtMem::deleteLevels (BinaryHeap_FourInts* H, int* levels, int numLevels) {

	// get list of nodes without this level
	for(int i=0;i<numLevels;i++){
		for (int slot=0; slot<numSlots; slot ++) {
			cntOnHeap[levels[i]][slot] = 0;
		}
	}

	std::list<int> *list = new std::list<int>();
	bool ok = true;
	for (int i=0; i<H->size(); i++){
		ok = true;
		for(int j=0;j<numLevels;j++){
				if (H->heap->at(i).third == levels[j])
					ok = false;
		}
		if (ok)
			list->push_back(i);

	}
	int K = list->size();

	if (K == 0) {
		H->makeEmpty();
	} else if (K < H->size()) { // it's possible that there are no entries to delete ... so no reason to rebuild
		//get the values for those nodes, so we can create a new heap

		int *v1 = new int[K];
		int *v2 = new int[K];
		int *v3 = new int[K];
		int *v4 = new int[K];
		float *key =  new float[K];
		int pos;
		int i=0;
		for (std::list<int>::iterator it=list->begin(); it != list->end(); ++it){
				pos = *it;
				v1[i] = H->heap->at(pos).first;
				v2[i] = H->heap->at(pos).second;
				v3[i] = H->heap->at(pos).third;
				v4[i] = H->heap->at(pos).fourth;
				key[i] = H->heap->at(pos).key;
				i++;
		}
		Float keys;
		keys.pointer = key;
		keys.length = K;
		H->buildAgain(v1, v2, v3, v4, keys);
		delete[] v1;
		delete[] v2;
		delete[] v3;
		delete[] v4;
		delete[] key;
	}
}
int ArrayHeapExtMem::size(){
	return this->n;
}

HeapReturn ArrayHeapExtMem::getBinaryHeapWithMin(){

	HeapReturn x;
	x.which = false;
	if (H1->isEmpty() && H2->isEmpty()) {
		if (n!=0) {
			Exception::critical();
		}
	//if (n==0)
		x.h = NULL;
		return x;
	} else if (H1->isEmpty()){
		x.h = (void*)H2;
		x.which = true;
		return x;
	}
	else if (H2->isEmpty()){
		x.h = (void*)H1;
		return x;
	}

	if (H1->heap->front().key <= H2->heap->front().key){
		x.h = (void*)H1;
		return x;
	}else{
		x.h = (void*)H2;
		x.which = true;
		return x;
	}
}
void ArrayHeapExtMem::removeMin(){


	bool which = false;

	if (H1->isEmpty() && H2->isEmpty())
		return;
	else if (H1->isEmpty()){
		which = true;
	}else if (H2->isEmpty())
		which = false;
	else if (H1->heap->front().key <= H2->heap->front().key)
		which = false;
	else{
		which = true;
	}


	n--;

	if (!which) {
		H1->deleteMin();
	} else{ // may be the last from it's slot - possibly refill

		int lvl = H2->heap->front().third;
		int slot = H2->heap->front().fourth;
		H2->deleteMin();
		// max number of nodes in a slot.

		cntOnHeap[lvl][slot]--;

		if (cntOnHeap[lvl][slot] == 0) {
			//if that was the last representative of one of the slots, then load another block from that slot
			load(lvl, slot);  // this will do nothing if the slot is now exhausted

		}


	}
}
bool ArrayHeapExtMem::isEmpty(){
	return this->n==0;
}
int compareMyTypeExt(const void * a, const void * b)
{
  if ( *(float*)a <  *(float*)b ) return -1;
  if ( *(float*)a == *(float*)b ) return 0;
  return 1;
}
bool ArrayHeapExtMem::test(bool verbose){
	std::string njTmpDir = "/tmp/";

	double valsd[] =
	{0.0834,0.01187,0.10279,0.09835,0.09883,0.1001,0.1129,0.09599,0.09468,0.09063,0.09083,0.08194,0.10182,0.09323,0.08796,0.09972,0.09429,0.08069,0.09008,0.10346,0.10594,0.09416,0.06915,0.08638,0.0886,0.09538,0.08546,0.09271,0.0936,0.09941,0.08026,0.0952,0.09446,0.09309,0.09855,0.08682,0.09464,0.0857,0.09154,0.08024,0.08824,0.09442,0.09495,0.08731,0.08428,0.08959,0.07994,0.08034,0.09095,0.09659,0.10066,0.0821,0.09606,0.12346,0.07866,0.07723,0.08642,0.08076,0.07455,0.07961,0.07364,0.08911,0.06946,0.07509,0.087,0.071,0.08653,0.07899,0.09512,0.09456,0.09161,0.08412,0.09649,0.09994,0.10151,0.09751,0.1019,0.10499,0.0873,0.1085,0.10189,0.09987,0.08912,0.10606,0.09552,0.08902,0.09158,0.08046,0.10687,0.0906,0.09937,0.09737,0.09825,0.10234,0.09926,0.09147,0.09071,0.09659,0.09472,0.09327,0.0949,0.09316,0.09393,0.09328,0.01187,0.00848,0.02284,0.03053,0.08393,0.08167,0.10191,0.06527,0.06613,0.06863,0.0652,0.06848,0.06681,0.07466,0.06444,0.05991,0.07031,0.06612,0.06873,0.06598,0.07283,0.06862,0.06437,0.06599,0.07291,0.06355,0.0685,0.06599,0.06593,0.0869,0.07364,0.08118,0.07693,0.06779,0.06605,0.07286,0.05655,0.06352,0.06105,0.09177,0.08312,0.0978,0.07464,0.07977,0.06241,0.07227,0.06255,0.0675,0.07953,0.07806,0.06702,0.08429,0.08567,0.0933,0.087,0.08809,0.07888,0.06351,0.08651,0.08294,0.07282,0.11102,0.08711,0.06192,0.0652,0.06957,0.06763,0.07123,0.0687,0.06773,0.06338,0.06694,0.09871,0.09221,0.08962,0.0879,0.09625,0.09953,0.09532,0.09903,0.0946,0.09406,0.09704,0.09877,0.07257,0.1001,0.09458,0.10141,0.10581,0.09824,0.10668,0.09835,0.10816,0.09667,0.08962,0.08486,0.08572,0.08324,0.08826,0.08801,0.09744,0.09916,0.09996,0.10054,0.10761,0.105,0.10604,0.10161,0.09155,0.10162,0.08549,0.10342,0.09419,0.11429,0.09764,0.09505,0.09394,0.10411,0.08792,0.08887,0.08648,0.07637,0.08544,0.08034,0.12373,0.12963,0.13817,0.13904,0.12648,0.13207,0.10788,0.09605,0.12674,0.08139,0.08326,0.08835,0.10922,0.103,0.12225,0.09854,0.09326,0.11181,0.089,0.12674,0.11631,0.0879,0.09866,0.11393,0.09839,0.09738,0.09922,0.1145,0.09967,0.1032,0.11624,0.10472,0.09999,0.09762,0.1075,0.11558,0.10482,0.10237,0.10776,0.08781,0.08771,0.09751,0.09025,0.09201,0.08731,0.08537,0.0887,0.0844,0.0804,0.08217,0.10216,0.07789,0.08693,0.0833,0.08542,0.09729,0.0937,0.09886,0.092,0.08392,0.09668,0.09444,0.09401,0.08657,0.09659,0.08553,0.0834,0.0846,0.10167,0.10447,0.09838,0.09545,0.09163,0.10475,0.09761,0.09475,0.09769,0.09873,0.09033,0.09202,0.08637,0.0914,0.09146,0.09437,0.08454,0.09009,0.08888,0.0811,0.12672,0.10517,0.11959,0.10941,0.10319,0.10544,0.10717,0.11218,0.12347,0.10637,0.11558,0.1198,0.10133,0.09795,0.10818,0.11657,0.10836,0.11127,0.09611,0.08462,0.1056,0.09537,0.09815,0.10385,0.10246,0.11299,0.11926,0.104,0.10309,0.09494,0.10078,0.09966,0.08215,0.09136,0.10058,0.10078,0.10121,0.09711,0.10072,0.10881,0.09396,0.09925,0.09221,0.0939,0.08804,0.09234,0.09647,0.07966,0.09939,0.09651,0.10765,0.10154,0.07889,0.10452,0.1023,0.10275,0.08817,0.0923,0.09237,0.09481,0.09309,0.08683,0.09903,0.08784,0.09309,0.08876,0.08442,0.097,0.10054,0.09463,0.10038,0.08208,0.10209,0.10181,0.10416,0.08065,0.09581,0.08961,0.08553,0.10272,0.08432,0.08437,0.08946,0.07594,0.07751,0.07935,0.07751,0.07714,0.09572,0.09626,0.08606,0.08031,0.08196,0.09758,0.0754,0.08671,0.10245,0.07644,0.07965,0.09553,0.08362,0.07587,0.08234,0.08611,0.09835,0.09917,0.09264,0.09656,0.0992,0.10802,0.10905,0.09726,0.09911,0.11056,0.08599,0.09095,0.10547,0.08824,0.09831,0.08445,0.09562,0.09378,0.08482,0.08686,0.09192,0.09617,0.09142,0.1024,0.10415,0.10673,0.08337,0.10091,0.08162,0.08284,0.08472,0.1021,0.09073,0.10521,0.09252,0.08545,0.09849,0.0891,0.10849,0.08897,0.08306,0.10775,0.10054,0.09952,0.10851,0.10823,0.10827,0.11254,0.11344,0.10478,0.11348,0.10646,0.12112,0.10183,0.1197,0.12399,0.11847,0.11572,0.14614,0.13348,0.12449,0.12358,0.12792,0.12525,0.12265,0.1305,0.13037,0.12684,0.12374,0.12907,0.12858,0.1285,0.12857,0.15825,0.15937,0.1467,0.128305,0.118165,0.119619995,0.117565,0.12769,0.11013			};

	int reps = 10000;
	int valLength = 502*reps;

	if(verbose)
	printf("Array Heap ExtMem test started... \n");
	if(verbose)
	printf("%d elements\n",valLength);

	float* vals = new float[valLength]();

	int i = 0;
	for (int j=0; j<reps; j++) {
		for (int l = 0;l<502;l++) {
			vals[i++] = (float)(valsd[l] + (.0001 * j));
		}
	}
	ArrayHeapExtMem *h = NULL;

	h = new ArrayHeapExtMem(njTmpDir, NULL);
	i=0;

	if(verbose)
	printf("Inserting... \n");
	for (int i=0;i<valLength;i++) {
		h->insert(i,0, vals[i]);
	}
	qsort((void*)vals,valLength,sizeof(float),compareMyTypeExt);


	i = 0;

	BinaryHeap_TwoInts *bh;
	BinaryHeap_FourInts *bh2;
	int val = 0;
	float key = 0;
	bool isOk = true;
	bool printB = false;

	if(verbose)
	printf("Removing and asserting... \n");
	while (!h->isEmpty() && isOk) {

		bool which = h->getBinaryHeapWithMin().which;
		if(which) bh2 = (BinaryHeap_FourInts*)(h->getBinaryHeapWithMin().h);
		else bh = (BinaryHeap_TwoInts*)(h->getBinaryHeapWithMin().h);


		if(which){
			val = bh2->heap->front().first;
			key = bh2->heap->front().key;
		}else{
			val = bh->heap->front().first;
			key = bh->heap->front().key;
		}

		if(vals[i]!=key)
			printB = true;
		else
			printB = false;

		if (printB) {

			printf("ACK\n");
			printf("%d : %f =? %f (%d)\n",i,vals[i],key,val);
			isOk = false;
		}
		i++;


		h->removeMin();
	}

	if(isOk && verbose)
		printf("Array Heap ExtMem test finished successfully... \n");
	else if(!isOk)
		printf("Array Heap ExtMem test failed... \n");

	delete[] vals;
	delete h;
	return isOk;
}
