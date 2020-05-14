/*
 * TreeBuilderBinHeap.cpp
 *
 *  Created on: Feb 14, 2016
 *      Author: michel
 */
#include "TreeBuilderBinHeap.hpp"

TreeBuilderBinHeap::TreeBuilderBinHeap(std::string** names, int** distances, int namesSize, int* clusterEqual) : TreeBuilder(names, distances, namesSize, clusterEqual){

	this->firstActiveNode = 0;
	this->nextActiveNode = new int[(this->K*2)-1];
	this->prevActiveNode = new int[(this->K*2)-1];
	for (int i=0; i<2*this->K-1; i++) {
		this->nextActiveNode[i] = i+1;
		this->prevActiveNode[i] = i-1;
	}

	this->candidateCountPerLoop = new int[(this->K)-1];
	this->heaps = NULL;

	this->candidatesD.pointer = NULL;
	this->candidatesI.pointer = NULL;
	this->candidatesJ.pointer = NULL;
	this->freeCandidates = NULL;
	this->candidatesActive = NULL;


	clusterAndHeap(this->K);

/*	} catch (Exception e){
		LogWriter.stdErrLogln("Error in Binary Heap Tree Builder");
		LogWriter.stdErrLogln(e.getMessage());
	}*/
}


TreeBuilderBinHeap::TreeBuilderBinHeap(std::string** names, int** distances, int namesSize) : TreeBuilder(names, distances, namesSize){

	this->firstActiveNode = 0;
	this->nextActiveNode = new int[(this->K*2)-1];
	this->prevActiveNode = new int[(this->K*2)-1];
	for (int i=0; i<2*this->K-1; i++) {
		this->nextActiveNode[i] = i+1;
		this->prevActiveNode[i] = i-1;
	}

	this->candidateCountPerLoop = new int[(this->K)-1];
	this->heaps = NULL;

	this->candidatesD.pointer = NULL;
	this->candidatesI.pointer = NULL;
	this->candidatesJ.pointer = NULL;
	this->freeCandidates = NULL;
	this->candidatesActive = NULL;


	clusterAndHeap(this->K);

/*	} catch (Exception e){
		LogWriter.stdErrLogln("Error in Binary Heap Tree Builder");
		LogWriter.stdErrLogln(e.getMessage());
	}*/
}
TreeBuilderBinHeap::~TreeBuilderBinHeap(){ //not nearly done, review
	delete[](this->nextActiveNode);
	delete[](this->prevActiveNode);
	delete[](this->candidateCountPerLoop);
	delete this->freeCandidates;
	delete[] this->candidatesD.pointer;
	delete[] this->candidatesI.pointer;
	delete[] this->candidatesJ.pointer;
	delete[] this->candidatesActive;
	for(int i=0;i<this->clustCnt;i++){
		for (int j=i; j<clustCnt; j++) {
			delete this->heaps[i][j];
		}
		delete [] this->heaps[i];
	}
	delete[] this->heaps;
}
void TreeBuilderBinHeap::clusterAndHeap (int maxIndex ){

	// pick clusters
	this->clustAssignments = new int[this->K];

	long maxT = 0;
	long minT = LONG_MAX;
	int i,j, ri, rj;

	i = this->firstActiveNode;
	while (i<maxIndex) {
		ri = this->redirect[i];
		if (this->R[ri] > maxT) maxT = this->R[ri];
		if (this->R[ri] < minT) minT = this->R[ri];
		i=this->nextActiveNode[i];
	}
	this->clustPercentiles = new long[this->clustCnt];
	long seedTRange = maxT - minT;
	for (i=0; i<this->clustCnt-1; i++) {
		this->clustPercentiles[i] = minT + (i+1)*seedTRange/this->clustCnt;
	}
	this->clustPercentiles[this->clustCnt-1] = maxT;
	this->clustSizes = new int[this->clustCnt]();
	i=this->firstActiveNode;
	while (i<maxIndex) {
		ri = this->redirect[i];
		for (j=0; j<this->clustCnt; j++) {
			if (this->R[ri]<=this->clustPercentiles[j]) {
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
	//int minPos;
	for (i=0; i<this->clustCnt; i++)  {
		//minPos = heap->heapArray[1];
		this->clustersBySize[i] = heap->heap->front().first;
		heap->deleteMin();
	}
	delete heap;


	// tell me about cluster sizes
	/*			if (verbose >= 4) {
		LogWriter.stdErrLogln("cluster sizes");
		for (i=0; i<clustCnt; i++)
			LogWriter.stdErrLogln(i + " : " + clustSizes[i] + " ( " + clustPercentiles[i] + " )");
	}*/


	//	make PQs for each cluster pair (make candidate list, too)  ... this will allow old ones to be collected
	if (this->candidatesD.pointer == NULL) {
		this->candidatesD.pointer = new int[10000];
		this->candidatesD.length = 10000;

		this->candidatesI.pointer = new int[10000];
		this->candidatesI.length = 10000;

		this->candidatesJ.pointer = new int[10000];
		this->candidatesJ.length = 10000;

		this->freeCandidates = new Stack();
		this->candidatesActive = new bool[10000];
	} //else - must already be created.  Just keep it the same size
	this->lastCandidateIndex = -1;

	if (this->heaps == NULL)
		this->heaps = new BinaryHeap_IntKey_TwoInts**[this->clustCnt];

	for(i=0;i<this->clustCnt;i++){
		this->heaps[i] = new BinaryHeap_IntKey_TwoInts*[this->clustCnt];
	}

	for (i=0; i<clustCnt; i++) {
		for (j=i; j<clustCnt; j++) {
				heaps[i][j] = new BinaryHeap_IntKey_TwoInts();
		}
	}

	int ra, rb, cra, crb ;
	int d;
	i=this->firstActiveNode;
	while (i<maxIndex) {
		ri = this->redirect[i];
		j = this->nextActiveNode[i];
		while (j<maxIndex) {

			rj = this->redirect[j];
			if (ri<rj) {
				ra = ri;
				rb = rj;
			} else {
				ra = rj;
				rb = ri;
			}
			if (this->clustAssignments[ra] < this->clustAssignments[rb]) {
				cra = this->clustAssignments[ra];
				crb = this->clustAssignments[rb];
			} else {
				cra = this->clustAssignments[rb];
				crb = this->clustAssignments[ra];
			}
			d = this->D[ra][rb-ra-1];
			this->heaps[cra][crb]->insert( i, j, d);   // this can be sped up for really big heaps by inserting all as one big list

			j = this->nextActiveNode[j];
		}
		i = this->nextActiveNode[i];
	}

}
TreeNode** TreeBuilderBinHeap::build (){

	int nextInternalNode = this->K;

//	this->nextActiveNode[3]++;
//	this->prevActiveNode[5]--;
//	this->candidatesActive[4] = false;
//	nextInternalNode++;

	int cand_cnt = 0;
	int defunct_cnt = 0;

	int i, j, x, ri, rj, rx, cluster, newK;
	int min_cand;
	int Dxi, Dxj, Dij, tmp;
	int a,b, ra, rb;
	int prev;
	int next;

	int clA, clB;

	long minQ, q, qLimit, minD;

	int maxT1[this->clustCnt];
	int maxT2[this->clustCnt];
	long maxT1val[this->clustCnt];
	long maxT2val[this->clustCnt];

	int stepsUntilRebuild = TreeBuilder::rebuildSteps;

	if ( stepsUntilRebuild == -1 ) stepsUntilRebuild = (int)(this->K * TreeBuilder::rebuildStepRatio);
	if ( stepsUntilRebuild < 500 )
		stepsUntilRebuild = this->K; //don't bother rebuilding for tiny trees



	newK = this->K;



	while (nextInternalNode<2*this->K-1) {// until there are 3 left ... at which point the merging is obvious


		//get two biggest T values for each cluster maxT1[] and maxT2[]
		for (i=0; i<this->clustCnt; i++) {
			maxT1[i] = maxT2[i] = -1;
			maxT1val[i] = maxT2val[i] = LONG_MIN;
		}
		x=this->firstActiveNode;
		while (x < nextInternalNode) {
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

		minQ = LONG_MAX;
		minD = LONG_MIN;

		//Go through current list of candidates, and find best so far.
		min_cand = -1;
		int clI, clJ;
		long maxTSum;
		int inactiveCnt = 0;
		for (x=this->lastCandidateIndex; x>=0; x--) {
			if (!this->candidatesActive[x]) {
				inactiveCnt++;
				continue;
			}

			ri = this->redirect[this->candidatesI.pointer[x]];
			rj = this->redirect[this->candidatesJ.pointer[x]];

			if (rj == -1 || ri == -1 /*leftovers from prior seqs redirected to this position*/) {
				this->candidatesActive[x] = false; // dead node ... can safely remove, 'cause we're going backwards through the list
				defunct_cnt++;
				if (x == this->lastCandidateIndex) {
					//scan left to find it
					int y = x;
					while (y>0 && !this->candidatesActive[y]) {y--;}
					this->lastCandidateIndex = y;
				}

//						if (nextInternalNode == 1118)
//							System.out.println("defuncted " + x);
			} else {
				q = (long)this->candidatesD.pointer[x] * (newK-2) - R[ri] - R[rj];
				if (q <= minQ) {
//							if (nextInternalNode == 1150)
//								System.out.println("(fib a) next node: " + nextInternalNode + ", cand: " + candidatesI[x] + ", " + candidatesJ[x] + " ( q=" +q + " , minQ=" + minQ +")");

					min_cand = x;
					minQ = q;
					minD = this->candidatesD.pointer[x];
				}
			}
		}

		this->candidateCountPerLoop[this->K-newK] = this->lastCandidateIndex-inactiveCnt+1;


		/*	frequently (every 50 or 100 iters?), scan through the candidates, and return
		 * entries to the array heap that shouldn't be candidates any more
		 */
		if ( (this->K-newK)%TreeBuilder::candidateIters == 0  && stepsUntilRebuild > TreeBuilder::candidateIters/2  ) {
			for (x=this->lastCandidateIndex; x>=0; x--) {
				if (!this->candidatesActive[x]) continue;

				ri = this->redirect[this->candidatesI.pointer[x]];
				rj = this->redirect[this->candidatesJ.pointer[x]];
				clI = this->clustAssignments[ri];
				clJ = this->clustAssignments[rj];

				maxTSum = maxT1val[clI] +  ( clI == clJ ? maxT2val[clI] :  maxT1val[clJ]) ;
				qLimit = (long)this->candidatesD.pointer[x] * (newK-2) - maxTSum;

				if (qLimit > minQ ) {
					// it won't be picked as a candidate in next iter, so stick it back in the cluster heaps
					this->candidatesActive[x] = false;
					if (x == this->lastCandidateIndex) {
						//scan left to find it
						int y = x;
						while (y>0 && !this->candidatesActive[y]) {y--;}
						this->lastCandidateIndex = y;
					}
					if (clI<=clJ)
						this->heaps[clI][clJ]->insert(this->candidatesI.pointer[x], this->candidatesJ.pointer[x], this->candidatesD.pointer[x]);
					else
						this->heaps[clJ][clI]->insert(this->candidatesI.pointer[x], this->candidatesJ.pointer[x], this->candidatesD.pointer[x]);
				}
			}
		}

		//compact the candidate list
		if (this->lastCandidateIndex>0) {
			if (inactiveCnt > (float)this->lastCandidateIndex/5) {
				int left = 0;
				int right = this->lastCandidateIndex;
				while (left < right) {
					while (left < right && this->candidatesActive[left]) left++;
					while (right > left && !this->candidatesActive[right]) right--;
					if (left < right) {
						this->candidatesD.pointer[left] = this->candidatesD.pointer[right];
						this->candidatesI.pointer[left] = this->candidatesI.pointer[right];
						this->candidatesJ.pointer[left] = this->candidatesJ.pointer[right];
						this->candidatesActive[right] = false;
						this->candidatesActive[left] = true;

						if (min_cand == right)
							min_cand = left;
//							System.out.println(right + " --> " + left);

						left++;
						right--;
					}
				}
				this->lastCandidateIndex = right;
				this->freeCandidates->clear();
			}
		}

		//pull off entries for the fib-heaps that have some chance of having minQ
		int h_i, h_j;
		int h_d;
		BinaryHeap_IntKey_TwoInts* h;
		for (a=0; a<this->clustCnt; a++) {
			for (b=a; b<this->clustCnt; b++) {

				clA = this->clustersBySize[a]<this->clustersBySize[b] ? this->clustersBySize[a] : this->clustersBySize[b];
				clB = this->clustersBySize[a]<this->clustersBySize[b] ? this->clustersBySize[b] : this->clustersBySize[a];

				maxTSum = maxT1val[clA] +  ( clA == clB ? maxT2val[clA] :  maxT1val[clB]) ;

				h = this->heaps[clA][clB];
				while (! h->isEmpty() ) {
					//minPos = h->heapArray[1];

					h_d = h->heap->front().key;
					h_i = h->heap->front().first;
					h_j = h->heap->front().second;


					ri = this->redirect[h_i];
					rj = this->redirect[h_j];

					if (rj==-1 || ri==-1 /*that's an old pointer*/) {
						h->deleteMin();//pull it off


						defunct_cnt++;
						continue;
					}
					q = (long)h_d * (newK-2);
					qLimit = q - maxTSum;
					q -=  this->R[ri] + this->R[rj];

					if (qLimit <= minQ) {
						// it's possible that this or one of the following nodes on the
						// heap has Q-value less than the best I've seen so far
						h->deleteMin();//pull it off

						int pos = appendCandidate(h_d, h_i, h_j);
//								LogWriter.stdErrLogln("appended " + heapNode.i + ", " + heapNode.j + " to " + pos + " (" + q + ")");
						cand_cnt++;

						if (q <= minQ) { // this is now the best I've seen
							min_cand = pos;
							minQ = q;
							minD = h_d;
						}
					} else {
						break; // no use pulling more off the heap ... they can't beat the best so far.
					}
				}
			}
		}

		//Now I know the position on the candidates array that has the best Q node.
		//Remove it from the candidates, merge the nodes, update D/T values
		// and possibly reset the candidates (every few iterations? when exceed some size?)
		int bestI = this->candidatesI.pointer[min_cand];
		int bestJ = this->candidatesJ.pointer[min_cand];
		this->candidatesActive[min_cand] = false;
		if (min_cand == this->lastCandidateIndex) {
			//scan left to find it
			int y = min_cand;
			while (y>0 && !this->candidatesActive[y]) {y--;}
			this->lastCandidateIndex = y;
		}
		this->freeCandidates->push(min_cand);

//				LogWriter.stdErrLogln("best was at " + min_cand);

		ri = this->redirect[bestI];
		rj = this->redirect[bestJ];

		this->nodes[nextInternalNode]->leftChild = this->nodes[bestI];
		this->nodes[nextInternalNode]->rightChild = this->nodes[bestJ];


		//assign branch lengths
		if (minD == FLT_MIN) {
			fprintf(stderr,"minD was not assigned correctly");
			Exception::critical();
		}
		if (newK==2) {
			this->nodes[bestI]->length = nodes[bestJ]->length = (float)minD / 200000000;
		} else {
			this->nodes[bestI]->length = ((float)minD + (R[ri]-R[rj])/(newK-2)) / 200000000 ;
			this->nodes[bestJ]->length = ((float)minD + (R[rj]-R[ri])/(newK-2)) / 200000000 ;
		}

		//if a length is negative, move root of that subtree around to compensate.
		if (this->nodes[bestI]->length < 0) {
			this->nodes[bestJ]->length += this->nodes[bestI]->length;
			this->nodes[bestI]->length = 0;
		} else if (nodes[bestJ]->length < 0) {
			this->nodes[bestI]->length += this->nodes[bestJ]->length;
			this->nodes[bestJ]->length = 0;
		}


		this->R[ri] = 0;
		this->redirect[bestI] = this->redirect[bestJ] = -1;

		// remove i,j from active list
		prev = this->prevActiveNode[bestI];
		next = this->nextActiveNode[bestI];
		this->prevActiveNode[next] = prev;
		if (prev == -1)
			this->firstActiveNode = next;
		else
			nextActiveNode[prev] = next;

		prev = this->prevActiveNode[bestJ];
		next = this->nextActiveNode[bestJ];
		this->prevActiveNode[next] = prev;
		if (prev == -1)
			this->firstActiveNode = next;
		else
			this->nextActiveNode[prev] = next;

		//calculate new D and T values
		x=this->firstActiveNode;
//				for (x=0; x<nextInternalNode; x++) {
		while (x<nextInternalNode) {

			rx = this->redirect[x];
			Dxj = rx<rj ? this->D[rx][rj-rx-1] : this->D[rj][rx-rj-1];
			Dij = ri<rj ? this->D[ri][rj-ri-1] : this->D[rj][ri-rj-1];

			if (rx < ri)  {
				Dxi = this->D[rx][ri-rx-1];
				a = x;
				ra = rx;
				b = nextInternalNode;
				rb = ri;
			} else {
				Dxi = this->D[ri][rx-ri-1];
				a = nextInternalNode;
				ra = ri;
				b=x;
				rb = rx;
			}

			tmp = (Dxi + Dxj - Dij) / 2;

			this->R[ri] += tmp;
			this->R[rx] += tmp - (Dxi + Dxj);

			//System.err.println(ri + ", " + rx + ": " + tmp + " (R[" + ri +"] = " + R[ri] + "), R[" + rx + "] = " + R[rx]);

			this->D[ra][rb-ra-1] = tmp;

			x = this->nextActiveNode[x];
		}

/*		if (TreeBuilder.verbose >= 3) {
			LogWriter.stdErrLogln("(inmem) next node: " + nextInternalNode + ": " + bestI + "(" + ri +") , " + bestJ + "(" + rj + ") Q = " + minQ + " (" + cand_cnt + " cands; " + defunct_cnt +  " defunct); ");
			LogWriter.stdErrLogln("     lengths: " + nodes[bestI].length + ", " + nodes[bestJ].length );
		}*/
		newK--;

		if ( stepsUntilRebuild == 0 ) {

/*			if (verbose >= 3) {
				LogWriter.stdErrLogln ("Resetting the clusters and corresponding PQs after " + (K-newK-1) + " iterations");
			}*/

			this->redirect[nextInternalNode++] = ri;
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

			// re-set the max levels for clusters (based on one-iteration-old data)
			for (j=0; j<this->clustCnt; j++) {
				//LogWriter.stdErrLogln("updating cluster " + j + " from " + clustPercentiles[j] + " to "  + maxT1val[j]);
				this->clustPercentiles[j] = maxT1val[j];
			}
			//then pick new cluster for ri, based on it's T value;
			for (j=0; j<this->clustCnt; j++) {
				if (this->R[ri]<=this->clustPercentiles[j]) {
					this->clustAssignments[ri] = j;
					break;
				}
			}

			// ... then add all new distances to the appropriate cluster pair
			x = this->firstActiveNode;
			int cra, crb;
			while (x<nextInternalNode) {

				rx = this->redirect[x];
				if (rx < ri)  {
					a = x;
					ra = rx;
					b = nextInternalNode;
					rb = ri;
				} else {
					a = nextInternalNode;
					ra = ri;
					b=x;
					rb = rx;
				}


				if (this->clustAssignments[ra] < this->clustAssignments[rb]) {
					cra = this->clustAssignments[ra];
					crb = this->clustAssignments[rb];
				} else {
					crb = this->clustAssignments[ra];
					cra = this->clustAssignments[rb];
				}
				this->heaps[cra][crb]->insert(a, b, D[ra][rb-ra-1]);


				x = nextActiveNode[x];
			}
			//prev = prevActiveNode[nextInternalNode]; // last active node
			//nextActiveNode[prev] = nextInternalNode; // now this is the last active node.
			//prevActiveNode[nextInternalNode] = prev;
			this->redirect[nextInternalNode++] = ri;

		}


	}



/*
	if (TreeBuilder.verbose >= 1) {
		LogWriter.stdErrLogln(cand_cnt + " candidates added");
		LogWriter.stdErrLogln(defunct_cnt + " defunct nodes removed");
	}
*/


/*	if (TreeBuilder.verbose >= 2) {
		long max_cnt = 0;
		long max_cnt_iter = 0;
		float max_ratio = 0;
		long sum_possible = 0;
		long sum_cnt = 0;
		for (i=1; i<candidateCountPerLoop.length; i++) {
			long cnt = candidateCountPerLoop[i];
			if (cnt > max_cnt) {
				max_cnt = cnt;
				max_cnt_iter = i;
			}
			sum_cnt += cnt;

			long all_pairs_cnt = ( (K-i) * (K-i-1) / 2 ) ;
			sum_possible += all_pairs_cnt;

			float ratio = (float)cnt / all_pairs_cnt;
			if (ratio>max_ratio) max_ratio = ratio;

		}

		LogWriter.stdErrLogln("max # candidates: " + max_cnt + " at iteration " + max_cnt_iter + " of " + candidateCountPerLoop.length);
		LogWriter.stdErrLogln( "max ratio of candidates to possible pairs: " + String.format("%.7f", (float)max_ratio/100) + "%");

		LogWriter.stdErrLogln("average # candidates: " + (sum_cnt/(K-4)) );
		LogWriter.stdErrLogln("total # candidates: " + sum_cnt + "  ( of " + sum_possible + " possible)");

		LogWriter.stdErrLogln("avg ratio of candidates to possible pairs: " + String.format("%.7f", ((float)sum_cnt/sum_possible)/100 ) + "%");

	}*/

	return this->nodes;
}
int TreeBuilderBinHeap::appendCandidate (int d, int i, int j) {
	int freePos;
	if (this->freeCandidates->isEmpty()) {
		freePos = this->lastCandidateIndex + 1;
	} else {
		freePos = this->freeCandidates->pop();
	}

	if (freePos > this->lastCandidateIndex) {
		this->lastCandidateIndex = freePos;
	}

	//make sure we don't overflow
	if (freePos == this->candidatesD.length) {
		int newLength = this->candidatesD.length * 10;
		int* D = new int[newLength];
		int* I = new int[newLength];
		int* J = new int[newLength];
		bool* act = new bool[newLength];
		for (int x=0; x<this->candidatesD.length; x++ ) {
			D[x] = this->candidatesD.pointer[x];
			I[x] = this->candidatesI.pointer[x];
			J[x] = this->candidatesJ.pointer[x];
			act[x] = this->candidatesActive[x];
		}
		delete[] this->candidatesD.pointer;
		delete[] this->candidatesI.pointer;
		delete[] this->candidatesJ.pointer;
		delete[] this->candidatesActive;

		this->candidatesD.pointer = D;
		this->candidatesI.pointer = I;
		this->candidatesJ.pointer = J;
		this->candidatesActive = act;

		this->candidatesD.length = newLength;
		this->candidatesI.length = newLength;
		this->candidatesJ.length = newLength;

	}
	this->candidatesD.pointer[freePos] = d;
	this->candidatesI.pointer[freePos] = i;
	this->candidatesJ.pointer[freePos] = j;
	this->candidatesActive[freePos] = true;

	return freePos;
}
