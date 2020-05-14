/*
 * TreeBuilderExtMem.cpp
 *
 *  Created on: Apr 4, 2016
 *      Author: michel
 */

#include "TreeBuilderExtMem.hpp"
#include "CandidateHeap.hpp"

TreeBuilderExtMem::TreeBuilderExtMem (std::string** names, int namesSize, float *R, std::string njTmpDir, FILE* diskD, float** memD, int memDFirstSize, int firstMemCol, int rowLength, long maxMemory){
	this->candFilePos = 0;
	this->numCandTriplesToDisk = 16384;
	this->usingSimpleCandidates = true;

	this->candidatesD = NULL;
	this->candidatesI = NULL;
	this->candidatesJ = NULL;
	this->freeCandidates = NULL;
	this->candidatesActive = NULL;
	this->candHeapList = NULL;
	this->arrayHeaps = NULL;

	if (!useBinPairHeaps && !useCandHeaps) {
		Exception::critical();
		//throw new Exception ("external memory method must use one of the heaps");
	}

	//clustCnt = TreeBuilder::clustCnt;
	this->rowLength = rowLength;
	this->njTmpDir = njTmpDir;
	this->names = names;
	this->diskD = diskD;
	this->candFile = NULL;
	this->memD = memD;
	this->maxMemory = maxMemory;

	this->memDSize = memDFirstSize;
	this->fBuff = new float[this->memDSize];


	this->R = R;
	this->RSize = namesSize;

	this->nextInternalNode = this->K = namesSize;

	this->curColInMem = this->firstColInMem = firstMemCol; //K;

	this->candidateCountPerLoop = new int[this->K-1];
	this->candidateViewsPerLoop = new int[this->K-1];
	this->candidateRowsCountPerLoop = new int[this->K-1];
	this->redirect = new int[2*this->K-1];
	this->nodes = new TreeNode*[2*this->K-1];

	int i;

	for (i=0; i<this->K; i++) {
		this->redirect[i] = i;
		this->nodes[i] = new TreeNode(names[i]);
	}

	for (i=this->K; i<2*this->K-1; i++) {
		this->redirect[i] = -1;
		this->nodes[i] = new TreeNode();
	}

	this->firstActiveNode = 0;
	this->nextActiveNode = new int[2*this->K-1];
	this->prevActiveNode = new int[2*this->K-1];
	for ( i=0; i<2*this->K-1; i++) {
		this->nextActiveNode[i] = i+1;
		this->prevActiveNode[i] = i-1;
	}

	this->newK = this->K;

	clusterAndHeap(this->K);
}
TreeBuilderExtMem::~TreeBuilderExtMem(){
	if (this->nextActiveNode!=NULL){
		delete[] this->nextActiveNode;
		this->nextActiveNode = NULL;
	}
	if (this->prevActiveNode!=NULL){
		delete[] this->prevActiveNode;
		this->prevActiveNode = NULL;
	}
	if (this->candFile!=NULL){
		//fprintf(stderr,"File deleted: TreeBuilderExtMem delete.\n");
		fclose(this->candFile);
	}
	if (this->candHeapList!=NULL){
		for (int i=0;i<(signed)this->candHeapList->size();i++)
			if (this->candHeapList->at(i)!=NULL)
				this->candHeapList->at(i)->clear();
		this->candHeapList->clear();
		this->candHeapList = NULL;
	}
	if (this->arrayHeaps != NULL){
		for(int i=0;i<this->clustCnt;i++){
			for(int j=0;j<this->clustCnt;j++)
				delete this->arrayHeaps[i][j];
			delete[] this->arrayHeaps[i];
		}
		delete[] this->arrayHeaps;
		this->arrayHeaps = NULL;
	}
	if(this->names != NULL){
		delete[] this->names;
	}
	if(this->nodes != NULL){
		delete[] this->nodes;
	}
}
void TreeBuilderExtMem::clusterAndHeap (int maxIndex ){

	if (this->candidatesD == NULL) {
		this->candidatesD = new float[10000];
		this->candidatesI = new int[10000];
		this->candidatesJ = new int[10000];
		this->freeCandidates = new Stack();
		this->candidatesActive = new bool[10000];

		this->candidatesSize = 10000;
	} else { //else - must already be created.  Just keep it the same size
		if (this->useCandHeaps) this->freeCandidates->clear();
	}
	this->lastCandidateIndex = -1;


	if (this->useCandHeaps) {
		if (this->candHeapList == NULL) {
			this->candHeapList = new std::vector<CandidateHeap*>();
		} else {
/*			if (TreeBuilder.verbose >= 2 && candHeapListsize() > 0){
				LogWriter.stdErrLogln("Cleared candidate heap list when K = " + (2 * K - nextInternalNode));
			}*/

			for (int i=0;i<(signed)this->candHeapList->size();i++)
				this->candHeapList->at(i)->clear();
			this->candHeapList->clear();
		}
	}

	int i,j, ri, rj;

	if (!this->useBinPairHeaps) { // just using candidate heaps

		this->starterCandHeap = new CandidateHeap(this->njTmpDir, NULL, this->newK, this, this->maxMemory/4);

	} else { // useBinPairHeapsHeaps

		// pick clusters
		this->clustAssignments = new int[this->K]();

		long maxT = 0;
		long minT = INT_MAX;

		i = this->firstActiveNode;
		while (i<maxIndex) {
			ri = this->redirect[i];
			if (this->R[ri] > maxT) maxT = this->R[ri];
			if (this->R[ri] < minT) minT = this->R[ri];
			i=this->nextActiveNode[i];
		}
		this->clustMins = new float[clustCnt]();
		this->clustMaxes = new float[this->clustCnt];
		float seedTRange = maxT - minT;
		for (i=0; i<this->clustCnt-1; i++) {
			this->clustMaxes[i] = minT + (i+1)*seedTRange/this->clustCnt;
		}
		this->clustMaxes[this->clustCnt-1] = maxT;
		this->clustSizes = new int[this->clustCnt]();
		i=this->firstActiveNode;
		while (i<maxIndex) {
			ri = this->redirect[i];
			for (j=0; j<this->clustCnt; j++) {
				if (this->R[ri]<=this->clustMaxes[j]) {
					this->clustAssignments[ri] = j;
					this->clustSizes[j]++;
					break;
				}
			}
			i=this->nextActiveNode[i];
		}


		//sort using a heap
		BinaryHeap *heap = new BinaryHeap();
		for (i=0; i<this->clustCnt; i++)
			heap->insert(i, this->clustSizes[i]);
		this->clustersBySize = new int[this->clustCnt];
		for (i=0; i<this->clustCnt; i++)  {
			this->clustersBySize[i] = heap->heap->front().first;
			heap->deleteMin();
		}
		delete heap;


		if (this->arrayHeaps == NULL){
			this->arrayHeaps = new ArrayHeapExtMem**[this->clustCnt];
			for(int i=0;i<this->clustCnt;i++){
				this->arrayHeaps[i] = new ArrayHeapExtMem*[this->clustCnt];
				for(int j=0;j<this->clustCnt;j++)
					this->arrayHeaps[i][j] = NULL;
			}
		}


		for (i=0; i<this->clustCnt; i++) {
			for (j=i; j<this->clustCnt; j++) {
				if (this->arrayHeaps[i][j] != NULL) {
					this->arrayHeaps[i][j]->deleteAll();
					this->arrayHeaps[i][j]->prepare();
				} else {
					this->arrayHeaps[i][j] = new ArrayHeapExtMem(this->njTmpDir, this->redirect, this->maxMemory/666 /*that's about 3MB if mem is 2GB*/);
					this->arrayHeaps[i][j]->A = i;
					this->arrayHeaps[i][j]->B = j;
				}
			}
		}
	}

	int cra=-1, crb=-1 ;
	float d = 0.0, q = 0.0;
	long diskPos = 0;
	int buffStart = 0;

	i=this->firstActiveNode;
	int numLeft = 0;


	while (i<maxIndex) {
		ri = this->redirect[i];
		buffStart = -this->memDSize; // get a new buffer

		j = this->nextActiveNode[i];
		while (j<maxIndex) {
			rj = this->redirect[j];

			if (this->useBinPairHeaps) {
				if (this->clustAssignments[ri] < this->clustAssignments[rj]) {
					cra = this->clustAssignments[ri];
					crb = this->clustAssignments[rj];
				} else {
					cra = this->clustAssignments[rj];
					crb = this->clustAssignments[ri];
				}
			}


			if (j>=this->firstColInMem) {
				d = this->memD[ri][j-this->firstColInMem];
			} else {
				if (j >= buffStart + memDSize) {
					//read in next page;
					while ( (buffStart += memDSize) + memDSize <= j); // probably just move to the next page

					diskPos = floatSize * ((long)rowLength * ri + buffStart )  ;
					fseek(diskD,diskPos,SEEK_SET);
					//diskD.seek(diskPos);

					numLeft = maxIndex - buffStart+1;
					if (numLeft >= memDSize) {
						fread(fBuff,sizeof(float),memDSize,diskD);
						//diskD.read(bBuff);
						//Arrays.byteToFloat(bBuff, fBuff);
					} else {
						fread(fBuff,sizeof(float),numLeft,diskD);
						//diskD.read(bBuff,0, 4*numLeft);
						//Arrays.byteToFloat(bBuff, fBuff, numLeft);
					}


				}
				d = fBuff[j - buffStart];

			}


			if (useBinPairHeaps) {
				arrayHeaps[cra][crb]->insert(i, j, d);
			} else {
				q = (newK-2) * d - R[ri] - R[rj];
				starterCandHeap->insert(i, j, q);
			}


			j = nextActiveNode[j];
		}
		i = nextActiveNode[i];
	}

	if (!useBinPairHeaps) {
		starterCandHeap->buildNodeList();
	}
}
TreeNode** TreeBuilderExtMem::build (){

	this->nextInternalNode = this->K;

	int cand_cnt = 0;
	int defunct_cnt = 0;

	int i, j, x, ri, rj, rx, cluster;
	float Dxi, Dxj, Dij, tmp;
	int a,b;
	int prev;
	int next;
	int min_i, min_j;


	int clA, clB;

	float minQ, q, qLimit, minD;

	int* maxT1 = new int[this->clustCnt]();
	int* maxT2 = new int[this->clustCnt]();
	float* maxT1val = new float[this->clustCnt]();
	float* maxT2val = new float[this->clustCnt]();

	int stepsUntilRebuild = TreeBuilder::rebuildSteps;

	if ( stepsUntilRebuild == -1 ) stepsUntilRebuild = (int)(K * TreeBuilder::rebuildStepRatio);
	if ( stepsUntilRebuild < 500 )
		stepsUntilRebuild = this->K; //don't bother rebuilding for tiny trees

	CandidateHeap *cHeap;


	float* fBuff_i = new float[this->memDSize]();
	float* fBuff_j = new float[this->memDSize]();


	//try {
		while (this->nextInternalNode<2*this->K-1) {// until there are 3 left ... at which point the merging is obvious

			this->usingSimpleCandidates = true;

			//get two biggest T values for each cluster maxT1[] and maxT2[]
			if (this->useBinPairHeaps) {
				for (i=0; i<clustCnt; i++) {
					maxT1[i] = maxT2[i] = -1;
					maxT1val[i] = maxT2val[i] = FLT_MIN;
				}
				x=this->firstActiveNode;
				while (x < this->nextInternalNode) {
					rx = this->redirect[x];
					cluster = this->clustAssignments[rx];
					if (this->R[rx] > maxT2val[cluster]){
						if (this->R[rx] > maxT1val[cluster]){
							maxT2val[cluster] = maxT1val[cluster];
							maxT1val[cluster] = this->R[rx];
							maxT2[cluster] = maxT1[cluster];
							maxT1[cluster] = rx;
						} else {
							maxT2val[cluster] = this->R[rx];
							maxT2[cluster] = rx;
						}
					}
					x = this->nextActiveNode[x];
				}
			}

			minQ = FLT_MAX;
			minD = FLT_MIN;
			//Go through current list of candidates, and find best so far.
			min_i = min_j = -1;

			float maxTSum;
			int inactiveCnt = 0;
			if (!returnCandsToHeaps) {
				for (x=lastCandidateIndex; x>=0; x--) {
					if (!candidatesActive[x]) {
						inactiveCnt++;
						continue;
					}

					ri = redirect[candidatesI[x]];
					rj = redirect[candidatesJ[x]];

					if (rj == -1 || ri == -1 /*leftovers from prior seqs redirected to this position*/) {
						candidatesActive[x] = false; // dead node ... can safely remove, 'cause we're going backwards through the list
						defunct_cnt++;
						if (x == lastCandidateIndex) {
							//scan left to find it
							int y = x;
							while (y>0 && !candidatesActive[y]) {y--;}
							lastCandidateIndex = y;
						}

					} else {
						q = candidatesD[x] * (newK-2) - R[ri] - R[rj];


						if (q <= minQ) {
							min_i = candidatesI[x];
							min_j = candidatesJ[x];
							minQ = q;
							minD = candidatesD[x];
						}
					}
				}

				candidateViewsPerLoop[K-newK] = candidateCountPerLoop[K-newK] = lastCandidateIndex-inactiveCnt+1;

				if (useCandHeaps) {
					for (int i=0;i<(int)candHeapList->size();i++) {
						candidateCountPerLoop[K-newK] += candHeapList->at(i)->size();
					}
				}

				if (useBinPairHeaps) {
					/*	frequently (every 50 or 100 iters?), scan through the candidates, and return
					 * entries to the array heap that shouldn't be candidates any more
					 */
					int clI, clJ;
					if ( (K-newK)%TreeBuilder::candidateIters == 0  && stepsUntilRebuild > TreeBuilder::candidateIters/2  ) {
						for (x=lastCandidateIndex; x>=0; x--) {
							if (!candidatesActive[x]) continue;

							ri = redirect[candidatesI[x]];
							rj = redirect[candidatesJ[x]];
							clI = clustAssignments[ri];
							clJ = clustAssignments[rj];

							maxTSum = maxT1val[clI] +  ( clI == clJ ? maxT2val[clI] :  maxT1val[clJ]) ;
							qLimit = candidatesD[x] * (newK-2) - maxTSum;

							if (qLimit > minQ ) {
								// it won't be picked as a candidate in next iter, so stick it back in the cluster heaps
								removeCandidate (x);


								if (clI<=clJ) {
									clA = clI;
									clB = clJ;
								} else {
									clA = clJ;
									clB = clI;
								}

								arrayHeaps[clA][clB]->insert(candidatesI[x], candidatesJ[x], candidatesD[x]);

							}
						}
					}
				}


				//compact the condidate list
				if (lastCandidateIndex> 0 && inactiveCnt > (float)lastCandidateIndex/5) {
					int left = 0;
					int right = lastCandidateIndex;
					while (left < right) {
						while (left < right && candidatesActive[left]) left++;
						while (right > left && !candidatesActive[right]) right--;
						if (left < right) {
							candidatesD[left] = candidatesD[right];
							candidatesI[left] = candidatesI[right];
							candidatesJ[left] = candidatesJ[right];
							candidatesActive[right] = false;
							candidatesActive[left] = true;

							left++;
							right--;
						}
					}
					lastCandidateIndex = right;
					inactiveCnt = 0;
					//freeCandidatesPos = -1;
					freeCandidates->clear();
				}


				if (useCandHeaps) {
					float qPrime, d;
					BinaryHeap_TwoInts *H;
					BinaryHeap_FourInts *H2;
					bool which;
					bool expiredHeapExists = false;
					// now go through the candidate heaps
					for (int c=candHeapList->size()-1; c>=0; c--) { // backwards through the list since more recent ones are likely to have bounded values closer to real values, so I want to use them to reduce the minQ first
						cHeap = candHeapList->at(c);

						cHeap->calcDeltaValues(newK);

						//pluck entries off the heap as long as the delta-based condition allows -> stick them in the cand. list
						while (!cHeap->isEmpty()) {
							which = cHeap->getBinaryHeapWithMin().which;
							if(which) H2 = (BinaryHeap_FourInts*)cHeap->getBinaryHeapWithMin().h;
							else H = (BinaryHeap_TwoInts*)cHeap->getBinaryHeapWithMin().h;

							if(which) qPrime = H2->heap->front().key;
							else qPrime = H->heap->front().key;

							candidateViewsPerLoop[K-newK]++;

							if ( cHeap->k_over_kprime * qPrime + cHeap->minDeltaSum  < minQ) { // "key" hold the q_prime value (q at time of insertion to candidate list)

								if(which){
									i = H2->heap->front().first;
									j = H2->heap->front().second;
								}else{
									i = H->heap->front().first;
									j = H->heap->front().second;
								}
								ri = redirect[i];
								rj = redirect[j];
								cHeap->removeMin(); // remove the one we're looking at
								if (rj == -1 || ri == -1 /*leftovers from prior seqs redirected to this position*/) {
									//don't bother to keep it for return to the minheap
									defunct_cnt++;
								} else {
									//convert q' back to d;
									d = (qPrime + cHeap->rPrimes[ri] + cHeap->rPrimes[rj]) / (cHeap->kPrime - 2);
									q = d * (newK-2) - R[ri] - R[rj];
	//								q = cHeap.k_over_kprime * (qPrime + cHeap.rPrimes[ri] + cHeap.rPrimes[rj]) - R[ri] - R[rj];
									appendCandidate( d, i, j);

									//System.err.println("B: " + i + ", " + j + ": " + q);


									if (q <= minQ) {
										min_i = i;
										min_j = j;

										minQ = q;
										minD = d;
									}
								}
							} else {
								break;
							}

						}

						//if the remaining size of this heap is too small, mark this heap as needing to be cleared.
						//if ( cHeap.size() / cHeap.representedRowCount < complexCandidateRatio || cHeap.size() < cHeap.origSize * .5) {
						//if ( cHeap.size() / cHeap.representedRowCount < (complexCandidateRatio/2) || cHeap.size() < cHeap.origSize * .5) {
						if (  cHeap->size() < cHeap->origSize * candHeapDecay) { // by now, lots of these are likely defunct anyway.  Reset.
							cHeap->expired = true;
							expiredHeapExists = true;
						}

					}


					//empty the heaps needing to be cleared - the entries will be merged into a new heap if appropriate limits are reached
					if (expiredHeapExists) {
						for (int c=candHeapList->size()-1; c>=0; c--) {
							cHeap = candHeapList->at(c);
							if (cHeap->expired) {
								while (!cHeap->isEmpty()) {
									which = cHeap->getBinaryHeapWithMin().which;
									if(which) H2 = (BinaryHeap_FourInts*)cHeap->getBinaryHeapWithMin().h;
									else H = (BinaryHeap_TwoInts*)cHeap->getBinaryHeapWithMin().h;

									if(which){
										i = H2->heap->front().first;
										j = H2->heap->front().second;
									}else{
										i = H->heap->front().first;
										j = H->heap->front().second;
									}

									ri = redirect[i];
									rj = redirect[j];
									if ( ri != -1 && rj != -1) {
										if(which) qPrime = H2->heap->front().key;
										else qPrime = H->heap->front().key;

										d =  (qPrime + cHeap->rPrimes[ri] + cHeap->rPrimes[rj]) / (cHeap->kPrime - 2);
										appendCandidate(d, i, j);
									}
									cHeap->removeMin();
								}
								candHeapList->erase(candHeapList->begin()+c);
								//if (TreeBuilder.verbose >= 2)
									//LogWriter.stdErrLogln("Removed a candidate heap (with K = " + cHeap.kPrime + ") when newK =" + newK );

							}
						}

					}
				}
			}

			if (returnCandsToHeaps) {
				//log the number of candidates grabbed at each iteration.
				//This tells us how the inner loop scales with the # seqs.
				candidateViewsPerLoop[K-newK] = candidateCountPerLoop[K-newK] = 0;
			}


			if (!useBinPairHeaps && !starterCandHeap->isEmpty()) {
				//pull off entries from the primary candidate heap that have some chance of having minQ
				float qPrime, d;
				BinaryHeap_TwoInts *H;
				BinaryHeap_FourInts *H2;
				bool which;


				starterCandHeap->calcDeltaValues(newK);

				//pluck entries off the heap as long as the delta-based condition allows -> stick them in the cand. list
				while (!starterCandHeap->isEmpty()) {
					which = cHeap->getBinaryHeapWithMin().which;
					if(which) H2 = (BinaryHeap_FourInts*)cHeap->getBinaryHeapWithMin().h;
					else H = (BinaryHeap_TwoInts*)cHeap->getBinaryHeapWithMin().h;

					if(which) qPrime = H2->heap->front().key;
					else qPrime = H->heap->front().key;

					candidateViewsPerLoop[K-newK]++;

					if ( starterCandHeap->k_over_kprime * qPrime + starterCandHeap->minDeltaSum  <= minQ) {
						if(which){
							i = H2->heap->front().first;
							j = H2->heap->front().second;
						}else{
							i = H->heap->front().first;
							j = H->heap->front().second;
						}
						ri = redirect[i];
						rj = redirect[j];
						starterCandHeap->removeMin(); // remove the one we're looking at
						if (rj == -1 || ri == -1 /*leftovers from prior seqs redirected to this position*/) {
							//don't bother to keep it for return to the minheap
							defunct_cnt++;
						} else {
							//convert q' back to d;
							d = (qPrime + starterCandHeap->rPrimes[ri] + starterCandHeap->rPrimes[rj]) / (starterCandHeap->kPrime - 2);
							q = d * (newK-2) - R[ri] - R[rj];
	//								q = cHeap.k_over_kprime * (qPrime + cHeap.rPrimes[ri] + cHeap.rPrimes[rj]) - R[ri] - R[rj];
							appendCandidate( d, i, j);

	//								System.err.println("C: " + i + ", " + j + ": " + q);


							if (q <= minQ) {
								min_i = i;
								min_j = j;
								minQ = q;
								minD = d;
							}
						}
					} else {
						break;
					}

				}

			} else if (useBinPairHeaps) {
				//pull off entries for the bin-pair heaps that have some chance of having minQ
				int h_i, h_j;
				float h_d;
				BinaryHeap_TwoInts *bh;
				BinaryHeap_FourInts *bh2;
				bool which;
				ArrayHeapExtMem *h;
				for (a=0; a<clustCnt; a++) {
					for (b=a; b<clustCnt; b++) {

						clA = clustersBySize[a]<clustersBySize[b] ? clustersBySize[a] : clustersBySize[b];
						clB = clustersBySize[a]<clustersBySize[b] ? clustersBySize[b] : clustersBySize[a];

						maxTSum = maxT1val[clA] +  ( clA == clB ? maxT2val[clA] :  maxT1val[clB]) ;

						h = arrayHeaps[clA][clB];

						while (!h->isEmpty()) {

							which = h->getBinaryHeapWithMin().which;
							if(which) bh2 = (BinaryHeap_FourInts*)h->getBinaryHeapWithMin().h;
							else bh = (BinaryHeap_TwoInts*)h->getBinaryHeapWithMin().h;

/*							if (bh == NULL) {
								//LogWriter.stdErrLogln("Surprising: null binary heap, for heap " + clA + ", " + clB);
								//h.describeHeap();
								throw new Exception("null binary heap");
							}*/
							if(which){
								h_i = bh2->heap->front().first;
								h_j = bh2->heap->front().second;
								h_d = bh2->heap->front().key;
							}else{
								h_i = bh->heap->front().first;
								h_j = bh->heap->front().second;
								h_d = bh->heap->front().key;
							}

							ri = redirect[h_i];
							rj = redirect[h_j];

							if (rj==-1 || ri==-1 /*that's an old pointer*/) {

								h->removeMin();//pull it off

								defunct_cnt++;
								continue;
							}
							q = h_d * (newK-2);
							qLimit = q - maxTSum;

	//								System.err.println("D: " + h_i + ", " + h_j + ": " + qLimit + " (" + q + ", " + maxTSum + ")");

							if (qLimit <= minQ) {
								// it's possible that this or one of the following nodes on the
								// heap has Q-value less than the best I've seen so far

								arrayHeaps[clA][clB]->removeMin();//pull it off
								appendCandidate (h_d, h_i, h_j);
								cand_cnt++;
								q -=  R[ri] + R[rj];

								if (q <= minQ) { // this is now the best I've seen
									min_i = h_i;
									min_j = h_j;
									minQ = q;
									minD = h_d;
								}
							} else {
								break; // no use pulling more off the heap ... they can't beat the best so far.
							}
						}

						if (returnCandsToHeaps) {

							//log the number of candidates grabbed at each iteration.
							//This tells us how the inner loop scales with the # seqs.
							candidateViewsPerLoop[K-newK] += lastCandidateIndex+1;
							candidateCountPerLoop[K-newK] += lastCandidateIndex+1;
							returnCandidates();
						}

					}
				}
			}

			if (useCandHeaps && !usingSimpleCandidates && !returnCandsToHeaps ) {
				candHeapList->at(candHeapList->size()-1)->buildNodeList();
			}


			//Now I know the position on the candidates array that has the best Q node.
			//Remove it from the candidates, merge the nodes, update D/T values
			// and possibly reset the candidates (every few iterations? when exceed some size?)


			ri = redirect[min_i];
			rj = redirect[min_j];

			nodes[nextInternalNode]->leftChild = nodes[min_i];
			nodes[nextInternalNode]->rightChild = nodes[min_j];

			//assign branch lengths
/*			if (minD == FLT_MIN) {
				throw new Exception("minD was not assigned correctly");
			}*/
			if (newK==2) {
				nodes[min_i]->length = nodes[min_j]->length = (float)minD / 2;
			} else {
				nodes[min_i]->length = (minD + (R[ri]-R[rj])/(newK-2)) / 2 ;
				nodes[min_j]->length = (minD + (R[rj]-R[ri])/(newK-2)) / 2 ;
			}
			//if a length is negative, move root of that subtree around to compensate.
			if (nodes[min_i]->length < 0) {
				nodes[min_j]->length += nodes[min_i]->length;
				nodes[min_i]->length = 0;
			} else if (nodes[min_j]->length < 0) {
				nodes[min_i]->length += nodes[min_j]->length;
				nodes[min_j]->length = 0;
			}


			// remove i,j from active list
			redirect[min_i] = redirect[min_j] = -1;

			prev = prevActiveNode[min_i];
			next = nextActiveNode[min_i];
			prevActiveNode[next] = prev;
			if (prev == -1)
				firstActiveNode = next;
			else
				nextActiveNode[prev] = next;

			prev = prevActiveNode[min_j];
			next = nextActiveNode[min_j];
			prevActiveNode[next] = prev;
			if (prev == -1)
				firstActiveNode = next;
			else
				nextActiveNode[prev] = next;


			//calculate new D and T values
			//    	for all these D reads, should check if the i/j/x value is greater than the most recently written-to-disk index
			//   	 ... if not, then read from the memD, not diskD
			R[ri] = 0;
			long diskPos;


			// I need to get this before going into the "foreach x" loop, 'cause I need Dij for all new vals.
			if (min_i>=firstColInMem) {
				Dij = memD[rj][min_i-firstColInMem];
			} else if (min_j>=firstColInMem) {
				Dij = memD[ri][min_j-firstColInMem];
			} else {
				diskPos = floatSize * ((long)rowLength * ri + min_j )  ;

	//			    	diskPos = 4 * (numDCols * (ri-(numDRowsPerFile*diskFile_i))  + h_j );
				fseek(diskD,diskPos,SEEK_SET);
				fread(&Dij,sizeof(float),1,diskD);
			}


			int buffStart_i = -memDSize;
			int buffStart_j = -memDSize;

			x=firstActiveNode;

			int numLeft;
			while (x<nextInternalNode) {

				rx = redirect[x];

				if (min_i>=firstColInMem) {
					Dxi = memD[rx][min_i-firstColInMem];
				} else if (x>=firstColInMem) {
					Dxi = memD[ri][x-firstColInMem];
				} else {
					if (x >= buffStart_i + memDSize) {
						//read in next page;
						while ( (buffStart_i += memDSize) + memDSize <= x); // probably just move to the next page

						diskPos = floatSize * ((long)rowLength * ri + buffStart_i );
						fseek(diskD,diskPos,SEEK_SET);

						numLeft = nextInternalNode - buffStart_i + 1;
						if (numLeft >= memDSize) {
							fread(fBuff_i,sizeof(float),memDSize,diskD);
						} else {
							fread(fBuff_i,sizeof(float),numLeft,diskD);
						}
					}
					Dxi = fBuff_i[x - buffStart_i];
				}


				if (min_j>=firstColInMem) {
					Dxj = memD[rx][min_j-firstColInMem];
				} else if (x>=firstColInMem) {
					Dxj = memD[rj][x-firstColInMem];
				} else {
					if (x >= buffStart_j + memDSize) {
						//read in next page;
						while ( (buffStart_j += memDSize) + memDSize <= x); // probably just move to the next page

						diskPos = floatSize * ((long)rowLength * rj + buffStart_j ) ;
						fseek(diskD,diskPos,SEEK_SET);

						numLeft = nextInternalNode - buffStart_j + 1;
						if (numLeft >= memDSize) {
							fread(fBuff_i,sizeof(float),memDSize,diskD);
						} else {
							fread(fBuff_i,sizeof(float),numLeft,diskD);
						}
					}
					Dxj = fBuff_j[x - buffStart_j];

				}

				//tmp =  Math.round(1000000 * (Dxi + Dxj - Dij) / 2)/(float)1000000;	 // this is the level of precision of input distances. No use allowing greater (noisy) precision to dominate decisions
				tmp =   (Dxi + Dxj - Dij) / 2;	 // this is the level of precision of input distances. No use allowing greater (noisy) precision to dominate decisions

				R[ri] += tmp;
				R[rx] += tmp - (Dxi + Dxj);

	//				instead of writing to the file, write to an in-memory D buffer.  if that buffer is full, then write it to disk
				//D[ra][rb-ra-1] = tmp;
	//System.err.println(ri + ", " + rx + ": " + tmp + " (R[" + ri +"] = " + R[ri] + "), R[" + rx + "] = " + R[rx]);
				memD[rx][nextInternalNode - firstColInMem] = tmp;
				if (x>=firstColInMem)
					memD[ri][x - firstColInMem] = tmp;


				x = nextActiveNode[x];
			}



			redirect[nextInternalNode] = ri;


			curColInMem++;
			if (curColInMem == firstColInMem + memDSize) {

				//write memD to diskD  ... there's just one buffer per row.

				// first, we append each row in memD to existing rows in the file (one disk block per row)
				x=firstActiveNode;
				while (x<nextInternalNode) {
					rx = redirect[x];

					diskPos = floatSize * ((long)rowLength * rx + firstColInMem )  ;
	//						diskPos = 4 * (numDCols * (rx-(numDRowsPerFile*diskFile)) + firstColInMem )   ;
					fseek(diskD,diskPos,SEEK_SET);
					fwrite(memD[rx],sizeof(float),memDSize,diskD);

					x = nextActiveNode[x];
				}

				//then write each column as a full new row in the file (which will cover multiple blocks)
				float fBuff_horiz[nextInternalNode];
				int ry;
				for (i=firstColInMem; i<curColInMem; i++ ) { // for each new column
					//write distances from the memD into a buffer,
					//which is then written to the row for this new node (ri)

					ry = redirect[i];
					if (ry==-1) continue; // this node was already merged with another, even before being dumped to file

					x=firstActiveNode;
					while (x<nextInternalNode) {
						//many entries in the buffer are left in default.
						// minor speedup could be had in the conversion routine (but I won't bother)
						rx = redirect[x];

						fBuff_horiz[x] = memD[rx][i-firstColInMem];
						x = nextActiveNode[x];
					}

					diskPos = floatSize * ((long)rowLength * ry  )  ;
					//diskPos = 4 * (numDCols * (ry-(numDRowsPerFile*diskFile))  )   ;
					fseek(diskD,diskPos,SEEK_SET);
					fwrite(fBuff_horiz,sizeof(float),nextInternalNode,diskD);

				}

			}

			newK--;

			if ( stepsUntilRebuild == 0 ) {
/*				if (TreeBuilder.verbose >= 3) {
					LogWriter.stdErrLogln ("Resetting the clusters and corresponding PQs after " + (K-newK-1) + " iterations");
				}*/

				redirect[nextInternalNode++] = ri;
				clusterAndHeap(nextInternalNode);

				if (newK < 200) {
					stepsUntilRebuild = newK; // almost done, quit shrinking the rebuild size
				} else {
					if ( TreeBuilder::rebuildSteps == -1 ) {
						if (TreeBuilder::rebuildStepsConstant)
							stepsUntilRebuild = (int)(K * TreeBuilder::rebuildStepRatio);
						else
							stepsUntilRebuild = (int)(newK * TreeBuilder::rebuildStepRatio);
					} else {
						stepsUntilRebuild = TreeBuilder::rebuildSteps;
					}
				}


			} else {
				stepsUntilRebuild--;

				if (useBinPairHeaps) {
					// re-set the max levels for clusters (based on one-iteration-old data)
					for (j=0; j<clustCnt; j++) {
						//LogWriter.stdErrLogln("updating cluster " + j + " from " + clustPercentiles[j] + " to "  + maxT1val[j]);
						clustMaxes[j] = maxT1val[j];
						clustMins[j] = maxT1val[j];
					}
					//then pick new cluster for ri, based on it's T value;
					for (j=0; j<clustCnt; j++) {
						if (R[ri]<=clustMaxes[j]) {
							clustAssignments[ri] = j;
							break;
						}
					}
				}

				// ... then add all new distances to the appropriate place (bin-pair heap or candidates)
				x = firstActiveNode;
				float d;
				while (x<nextInternalNode) {

					rx = redirect[x];
					d = memD[rx][nextInternalNode - firstColInMem];

					if (useBinPairHeaps) {
						if (clustAssignments[ri]<=clustAssignments[rx]) {
							clA = clustAssignments[ri];
							clB = clustAssignments[rx];
						} else {
							clA = clustAssignments[rx];
							clB = clustAssignments[ri];
						}

						arrayHeaps[clA][clB]->insert(nextInternalNode, x, d);
					} else {
						appendCandidate(d, nextInternalNode, x);
					}

					x = nextActiveNode[x];


				}


				nextInternalNode++;

			}

			if (curColInMem == firstColInMem + memDSize)
		//	if (curColInMem == firstColInMem + memD[0].length)
				firstColInMem = curColInMem;

		}

	return nodes;

}

void TreeBuilderExtMem::returnCandidates (){
	int cra, crb, x, ri, rj, i, j;
	float d;

	for (x=lastCandidateIndex; x>=0; x--) {
		if (!candidatesActive[x]) {
			continue;
		}

		ri = redirect[candidatesI[x]];
		rj = redirect[candidatesJ[x]];

		if (rj == -1 || ri == -1 /*leftovers from prior seqs redirected to this position*/) {
			candidatesActive[x] = false; // dead node ... can safely remove, 'cause we're going backwards through the list
			if (x == lastCandidateIndex) {
				//scan left to find it
				while (x>0 && !candidatesActive[x]) {x--;}
				lastCandidateIndex = x;
			}
		} else {

			if (clustAssignments[ri] < clustAssignments[rj]) {
				cra = clustAssignments[ri];
				crb = clustAssignments[rj];
			} else {
				cra = clustAssignments[rj];
				crb = clustAssignments[ri];
			}
			arrayHeaps[cra][crb]->insert(candidatesI[x], candidatesJ[x], candidatesD[x]);

		}
	}

	lastCandidateIndex = -1;
	freeCandidates->clear();

	if (candFilePos>0) {
		int* iBuff = new int[numCandTriplesToDisk*3];

		while (candFilePos>0) {
			//use bBuff, which is already in memory
			fseek(candFile,(candFilePos-numCandTriplesToDisk) * 3 * 4,SEEK_SET);
			fread(iBuff,sizeof(int),numCandTriplesToDisk*3,candFile);

			for (x=0; x<numCandTriplesToDisk; x++) {
				i = iBuff[x*3];
				j = iBuff[x*3+1];

				d = (float)(iBuff[x*3+2]);

				ri = redirect[i];
				rj = redirect[j];

				if (rj == -1 || ri == -1 /*leftovers from prior seqs redirected to this position*/) {
					//
				} else {
					if (clustAssignments[ri] < clustAssignments[rj]) {
						cra = clustAssignments[ri];
						crb = clustAssignments[rj];
					} else {
						cra = clustAssignments[rj];
						crb = clustAssignments[ri];
					}
					arrayHeaps[cra][crb]->insert(i,j,d);
				}

			}
			candFilePos-=numCandTriplesToDisk;
		}
	}
}
void TreeBuilderExtMem::appendCandidate (float d, int i, int j){

	if (useCandHeaps ) {

		if (usingSimpleCandidates) {

			int candCnt = lastCandidateIndex - freeCandidates->length() + 1 ;

			if ( candCnt >= 2000000 || (candCnt >= candHeapThresh && candCnt/newK > complexCandidateRatio ) ) {
				usingSimpleCandidates = false;
				CandidateHeap *heap = new CandidateHeap(njTmpDir, NULL, newK, this, maxMemory/1000 /* roughly 2MB for 2GB RAM */);
				float q=0;

				//try {
					for (int x=0; x<=lastCandidateIndex; x++ ) {
						if (candidatesActive[x] && redirect[candidatesI[x]] != -1 && redirect[candidatesJ[x]] != -1) {
							q = candidatesD[x] * (newK-2) - R[redirect[candidatesI[x]]] - R[redirect[candidatesJ[x]]];
							heap->insert(candidatesI[x], candidatesJ[x], q);
						}
					}
					q = d * (newK-2) - R[redirect[i]] - R[redirect[j]];
					heap->insert(i, j, q);
/*				} catch (Exception e) {
					LogWriter.stdErrLogln("Death while appending candidate.  nextInternalNode = " + nextInternalNode);
					throw e;
				}*/


				//just cleared candidates ... clean up data structure
				lastCandidateIndex = -1;
				freeCandidates->clear();

				CandidateHeap *expiringHeap;
				BinaryHeap_TwoInts *H;
				BinaryHeap_FourInts *H2;
				bool which;
				int ri, rj;
				float qPrime;
				if (candHeapList->size() == 100 ) { // don't let the number of heaps exceed 100 ... merge the oldest 20 into this one
					for (int c=0; c<20;c++) {
						expiringHeap = candHeapList->at(0);
						while (!expiringHeap->isEmpty()) {
							which = expiringHeap->getBinaryHeapWithMin().which;
							if (which) H2 = (BinaryHeap_FourInts*)expiringHeap->getBinaryHeapWithMin().h;
							else H = (BinaryHeap_TwoInts*)expiringHeap->getBinaryHeapWithMin().h;

							if(which){
								i = H2->heap->front().first;
								j = H2->heap->front().second;
							}else{
								i = H->heap->front().first;
								j = H->heap->front().second;
							}
							ri = redirect[i];
							rj = redirect[j];
							if ( ri != -1 && rj != -1) {
								if (which) qPrime = H2->heap->front().key;
								else qPrime = H->heap->front().key;

								d =  (qPrime + expiringHeap->rPrimes[ri] + expiringHeap->rPrimes[rj]) / (expiringHeap->kPrime - 2);
								q = d * (newK-2) - R[ri] - R[rj];
								heap->insert(i, j, q);
							}
							expiringHeap->removeMin();
						}
						candHeapList->erase(candHeapList->begin());
						//if (TreeBuilder.verbose >= 2)
							//LogWriter.stdErrLogln("Removed a candidate heap (with K = " + expiringHeap.kPrime + ") when newK =" + newK );

					}
				}


				candHeapList->push_back(heap);
				//if (TreeBuilder.verbose >= 2)
					//LogWriter.stdErrLogln("Added a candidate heap when newK =" + newK + " : total heap count is " + candHeapList.size());


			} else {

				appendCandToSimpleList (d, i, j);

			}

		} else { // not using simple candidate list.

			d = d * (newK-2) - R[redirect[i]] - R[redirect[j]]; // really q' now

			candHeapList->at(candHeapList->size()-1)->insert(i, j, d);
		}

	}	else { // not using candHeaps
		appendCandToSimpleList (d, i, j);
	}
}
int TreeBuilderExtMem::appendCandToSimpleList (float d, int i, int j){

	int freePos;
	if (freeCandidates->isEmpty()) {
		freePos = lastCandidateIndex + 1;
	} else {
		freePos = freeCandidates->pop();
	}



	//make sure we don't overflow
	if (freePos == this->candidatesSize) { //is it the right number?
		int newLength = (int)(this->candidatesSize * (this->candidatesSize < 1000000 ? 10 : 2 ));

		int x;

		if (newLength>8000000) {
		//if (newLength> 1040) {
		//	System.err.println("Need to change length test and numTriples");

			//need to write some of the current entries to a file, and add entries to freeCandidates
			//try {
				int* iBuff = new int[numCandTriplesToDisk*3];
				for (x=0; x<numCandTriplesToDisk; x++) {
					iBuff[x*3] = candidatesI[lastCandidateIndex];
					iBuff[x*3+1] = candidatesJ[lastCandidateIndex];
					iBuff[x*3+2] = (float)(candidatesD[lastCandidateIndex]);
					freeCandidates->push(lastCandidateIndex);
					candidatesActive[lastCandidateIndex] = false;
					lastCandidateIndex--;
				}
				freePos = lastCandidateIndex + 1;

				if (candFile == NULL) {
					std::string aux = njTmpDir + "/NINJA/candDisk";
					candFile = fopen(aux.c_str(),"rw");
				}
				fseek(candFile,candFilePos * 3 * 4,SEEK_SET);
				fwrite(iBuff,sizeof(int),numCandTriplesToDisk*3,candFile);
				candFilePos += numCandTriplesToDisk;

/*			} catch (Exception e) {
				throw e;
			}*/
		} else {

			float* D = new float[newLength];
			for (x=0; x<this->candidatesSize; x++ )
				D[x] = candidatesD[x];
			candidatesD = D;

			int* I = new int[newLength];
			for (x=0; x<this->candidatesSize; x++ )
				I[x] = candidatesI[x];
			candidatesI = I;

			int* J = new int[newLength];
			for (x=0; x<this->candidatesSize; x++ )
				J[x] = candidatesJ[x];
			candidatesJ = J;

			bool* act = new bool[newLength];
			for (x=0; x<this->candidatesSize; x++ )
				act[x] = candidatesActive[x];
			candidatesActive = act;

			this->candidatesSize = newLength;
		}
	}
	candidatesD[freePos] = d;
	candidatesI[freePos] = i;
	candidatesJ[freePos] = j;
	candidatesActive[freePos] = true;
	if (freePos > lastCandidateIndex) {
		lastCandidateIndex = freePos;
	}
	return freePos;
}
void TreeBuilderExtMem::removeCandidate (int x) {


	//freeCandidates[++freeCandidatesPos] = x;
	freeCandidates->push(x);

	candidatesActive[x] = false;
	if (x == lastCandidateIndex) {
		//scan left to find it
		while (x>0 && !candidatesActive[x]) {x--;}
		lastCandidateIndex = x;
	}
}
