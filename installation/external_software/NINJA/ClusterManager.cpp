/*
 * ClusterManager.cpp
 *
 * A naive agglomerative implemenation of single-linkage 
 * clustering ( a.k.a. nearest-neighbor clustering or
 * heirarchical clustering ).
 * 
 *  Created on: Jul 12, 2019
 *      Author: robert
 */

#include "ClusterManager.hpp"

#define LINUX 1

#ifdef LINUX
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

//standard constructor
ClusterManager::ClusterManager(std::string method, std::string njTmpDir, std::string inFile, FILE* outFile,
InputType inType, OutputType outType, AlphabetType alphType, CorrectionType corrType, int threads, bool useSSE, bool
printTime, float clusterCutoff ){
	this->method = method;
	this->njTmpDir = njTmpDir;
	this->inFile = inFile;
	this->outFile = outFile;
	this->inType = inType;
	this->outType = outType;
	this->alphType = alphType;
	this->corrType = corrType;
	this->names = NULL;
	this->chars = "abcdefghijklmonpqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
	this->newDistanceMethod = false;
	this->threads = threads;
        this->clusterCutoff = clusterCutoff;
	this->newDistanceMethod = useSSE;
	this->printTime = printTime;
}

/*
int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p){
return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
        ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
}
*/

std::string ClusterManager::doJob(){
	double** distances = NULL;
	float** memD = NULL;
	float* R  = NULL;
	// all for external memory version
	int pageBlockSize = 1024; //that many ints = 4096 bytes;
	FILE* diskD = NULL;

	int rowLength = 0;
	int firstMemCol = -1;

	int numCols = 0;

	//Runtime runtime = Runtime.getRuntime();
	long maxMemory = -1;

	bool ok = true;
	std::string clusterString = "";

	//NinjaLogWriter.printVersion();

	int K=0;

	int* equalCluster;

        struct timespec start, afterRead, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

	SequenceFileReader* seqReader = NULL;
	if (!this->method.compare("extmem")){
		fprintf(stderr,"External memory clustering is not implemented.\n");
		Exception::critical();
	}else{
		DistanceReader* reader = NULL;

		if (this->inType == alignment) {
			seqReader = new SequenceFileReader(&(this->inFile),(SequenceFileReader::AlphabetType) this->alphType);
			std::string** seqs = seqReader->getSeqs();
			this->names = seqReader->getNames();

            clock_gettime(CLOCK_MONOTONIC, &afterRead); 

			this->alphType = (ClusterManager::AlphabetType) seqReader->getAlphType();
			fprintf(stderr,"Calculating distances....\n");
			DistanceCalculator* distCalc = new DistanceCalculator(seqs,(DistanceCalculator::AlphabetType) alphType,(DistanceCalculator::CorrectionType)  corrType, seqReader->numSeqs,newDistanceMethod);
			K = seqReader->numSeqs;
			reader = new DistanceReader(distCalc, K, this->threads);
		}else{
			reader = new DistanceReader(this->inFile);
			K = reader->K;
			this->names = new std::string*[K];
			for (int i = 0;i<K;i++)
				this->names[i] = new std::string();
		}

                // 
                //  Distance matrix occupancy:
                //   
                //          s1 s2 s3 s4    [
                //      s1   0  #  #  #  =>  [ #, #, # ], 
                //      s2   0  0  #  #  =>  [ #, # ],
                //      s3   0  0  0  #  =>  [ # ]
                //                         ]
                //
		distances = new double*[K];
		for (int i=0; i<K; i++) {
			distances[i] = new double[K - i - 1];
		}
		reader->readDoubles( this->names, distances);

                typedef std::vector<int> ClusterMembers;
                std::vector<ClusterMembers> clusters;

                // Heirarchical

                // The sparsified index to the nearest-then-first seq/cluster.
                //   Sparsified means that for seq/cluster s1, if the nearest
                //   cluster is s3 then min_dist[0] = 1;  To get back to the
                //   true index simply: min_dist[m] + m + 1 
                int min_dist[K];
                for (int i=0; i<K-1; i++) {
                  min_dist[i] = 0;
              
                  // Initialize clusters
                  clusters.push_back( std::vector<int>{ i } );
               
                  //printf("Looking at row %d...\n", i);
                  for ( int j=0; j<(K-i-1); j++ ) {
                    //printf("  col %d distance = %f\n", j, distances[i][j]);
                    if ( distances[i][j] < distances[i][min_dist[i]]  ) {
                      min_dist[i] = j;
                      //printf("min_dist[%d] = %d\n", i, j);
                    }
                  }
                }
                clusters.push_back( std::vector<int>{ K-1 } );
              
                if ( 0 ) {
                  // DEBUGING
                  printf("Initialized %d clusters...\n",clusters.size());
                  for( int j = 0; j < clusters.size(); j++ ) {
                    for( int k = 0; k < clusters[j].size(); k++ ) {
                      printf("%d\t%s\n", j, this->names[clusters[j][k]]->c_str());
                    }
                  }
                  printf("Distance matrix:\n");
                  for (int i=0; i<K-1; i++){
                    for ( int j=0; j<K; j++ )
                      if ( j <= i ) 
                        printf("  -- ");
                      else
                        printf(" %0.2f", distances[i][j-i-1]);
                    printf("\n");
                  }
                }
              
                // Cluster
                //printf("clustering started...\n");
                for ( int m=0; m<K-1; m++ ) {
              
                  int clust1 = 0;
                  for ( int i = 0; i<K-1; i++ ) {
                    //printf("Iterating i = %d, min_dist[%d] = %d, distances[i][min_dist[i]] = %f\n", i, i, min_dist[i], distances[i][min_dist[i]]);
                    if ( distances[i][min_dist[i]] < distances[clust1][min_dist[clust1]] )
                      clust1 = i;
                  }
                  int clust2 = min_dist[clust1] + clust1 + 1;
                  //printf("  Minimum %f @ clust1 = %d, clust2 = %d\n", distances[clust1][min_dist[clust1]], clust1, clust2 );
                  
                  // Are we done?
                  if ( distances[clust1][min_dist[clust1]] > (double)clusterCutoff ) 
                    break;
              
                  // merge minimums from row=clust2 into row=clust1 and
                  // max-out values in row=clust2
                  for ( int c2col = 0; c2col < K-clust2-1; c2col++ ) {
                    // This is iterating over row=clust2
                    int c1col = c2col+(clust2-clust1);
                    if ( distances[clust1][c1col] > distances[clust2][c2col] ){ 
                      //printf("RowUpd: Setting distance [%d][%d] to %f\n", clust1, c1col, distances[clust2][c2col] );
                      distances[clust1][c1col] = distances[clust2][c2col];
                    }
                    distances[clust2][c2col] = (double)3.0;
                  }
                  for ( int i = 0; i < clust1; i++ ){
                    int c1col = clust1 - i - 1;
                    int c2col = clust2 - i - 1;
                    if ( distances[i][c1col] > distances[i][c2col] ) 
                    {
                      //printf("ColUpd: Setting distance [%d][%d] to %f\n", i, c1col, distances[i][c2col] );
                      distances[i][c1col] = distances[i][c2col];
                    }
                  }
                  // max out the column
                  for ( int i = 0; i < clust2 ; i++ ) {
                    distances[i][clust2-i-1] = (double)3.0;
                  }
              
                  if ( 0 ) {
                    printf("Distance matrix:\n");
                    for (int i=0; i<K-1; i++){
                      for ( int j=0; j<K; j++ )
                        if ( j <= i ) 
                          printf("  -- ");
                        else
                          printf(" %0.2f", distances[i][j-i-1]);
                      printf("\n");
                    }
                  }
              
                  //printf("Fixing clusters\n");
                  //printf("  Moving clst %d to clst %d\n", clust2, clust1 );
                  for( int k = 0; k < clusters[clust2].size(); k++ )
                    clusters[clust1].push_back(clusters[clust2][k]);
                  clusters[clust2].clear();
                    
                  //printf("Fixing min_dist clust2-clust1-1 = %d\n", clust2-clust1-1);
                  for ( int i = 0; i < clust1; i++ ) 
                    if ( min_dist[i] == clust2-i-1 ) 
                       min_dist[i] = clust1-i-1;
                  for ( int j = 0; j < K - clust1 - 1; j++ ) 
                    if ( distances[clust1][j] < distances[clust1][min_dist[clust1]] )
                      min_dist[clust1] = j;
              
                  if ( 0 ) {
                    printf("Mindist now:\n");
                    for( int i = 0; i < K-1; i++ ) 
                      printf(" min_dist[%d] = %d\n", i, min_dist[i]);
                    }
                  }

                int cidx = 0;
                for( int j = 0; j < clusters.size(); j++ ) {
                  if ( clusters[j].size() > 0 )  {
                    for( int k = 0; k < clusters[j].size(); k++ ) {
                      fprintf(outFile,"%d\t%s\n", cidx, this->names[clusters[j][k]]->c_str());
                      printf("%d\t%s\n", cidx, this->names[clusters[j][k]]->c_str());
                    }
                    cidx++;
                  }
                }
                if ( cidx == 1 ) 
                  printf("There is 1 cluster\n");
                else
                  printf("There are %d clusters\n", cidx);

        }


	if (seqReader != NULL)
		delete seqReader;

	delete[] distances;
	return ( clusterString );
}
