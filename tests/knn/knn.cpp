#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>
#include "../src/UtilsCUDA/GMEMcommonCUDA.h" // for MAX_REF_POINTS
#include "../src/UtilsCUDA/external/book.h"

#include <algorithm>

int main(void)
{
	// Variables and parameters
	float* ref;                 // Pointer to reference point array
	float* query;               // Pointer to query point array: order is x1,y1,z1,x2,y2,z2... to be cache friendly
	int*   ind;                 // Pointer to index array: size query_nb*maxGaussiansPerVoxel
	int    ref_nb     = 1377;   // Reference point number
	int    query_nb   = 100000;   // Query point number
	//Defined as constants now int    dimsImage        = 3;     // Dimension of points
	//int    maxGaussiansPerVoxel          = 5;     // Nearest neighbors to consider
	int    iterations = 5;     //at each iteration we will upload the query points (to simulate our case of maxGaussiansPerVoxel-NN
	int    i;

	if(MAX_REF_POINTS<ref_nb)
	{
		//TODO allow th epossibility of more ref_points by using global memory instead of constant memory
		printf("ERROR!! Increase MAX_REF_POINTS!\n");
		exit(2);
	}
	if(dimsImage!=3)
	{
		printf("ERROR: dimsImage should be 3\n");
		exit(2);
	}
	
	// Memory allocation
	ref    = (float *) malloc(ref_nb   * dimsImage * sizeof(float));
	query  = (float *) malloc(query_nb * dimsImage * sizeof(float));
	ind    = (int *)   malloc(query_nb * maxGaussiansPerVoxel * sizeof(int));

	// Init 
	srand((unsigned int)time(NULL));
	for (i=0 ; i<ref_nb   * dimsImage ; i++) ref[i]    = (float)rand() / (float)RAND_MAX;
	for (i=0 ; i<query_nb * dimsImage ; i++) query[i]  = (float)rand() / (float)RAND_MAX;

	
	//CUDA variables
	int *indCUDA;
	float* queryCUDA;
	float scale[dimsImage]={1.0,1.0,1.0};
	float *anchorCUDA;
	
	//select GPU and check maximum meory available and check that we have enough memory
	int memoryNeededInBytes=query_nb*maxGaussiansPerVoxel*sizeof(int)+query_nb   * dimsImage * sizeof(float);	
	cudaDeviceProp prop;
	int dev;
	memset( &prop, 0, sizeof( cudaDeviceProp ) );
	prop.totalGlobalMem=memoryNeededInBytes;
	HANDLE_ERROR( cudaChooseDevice( &dev, &prop ) );
	printf( "Memory required: %d;CUDA choosing device:  %d\n", memoryNeededInBytes,dev );
	HANDLE_ERROR( cudaSetDevice( dev ) );

	
	
	// allocate memory on the GPU for the output: it will only be done once in the whole program
	HANDLE_ERROR( cudaMalloc( (void**)&indCUDA, query_nb*maxGaussiansPerVoxel*sizeof(int) ) );//should it be a texture memory?NO. It does not fit in Cuda2Darray but it fits in linear 1Dtexture, although it does not seems to bring benefits
	HANDLE_ERROR( cudaMalloc( (void**)&queryCUDA, query_nb*dimsImage*sizeof(float) ) );
	HANDLE_ERROR( cudaMalloc( (void**)&anchorCUDA, ref_nb*dimsImage*sizeof(float) ) );
	
	// capture the start time
	cudaEvent_t     start, stop;
	HANDLE_ERROR( cudaEventCreate( &start ) );
	HANDLE_ERROR( cudaEventCreate( &stop ) );
	HANDLE_ERROR( cudaEventRecord( start, 0 ) );
	  
	  // Copy image data to array
	   HANDLE_ERROR(cudaMemcpy(queryCUDA,query, dimsImage*query_nb*sizeof(float), cudaMemcpyHostToDevice));
	   HANDLE_ERROR(cudaMemcpyToSymbol(scaleCUDA,scale, dimsImage * sizeof(float)));//constant memory	
	
	
	
	// generate a bitmap from our sphere data
	int numThreads=std::min(MAX_THREADS,query_nb);
	int numGrids=std::min(MAX_BLOCKS,(query_nb+numThreads-1)/numThreads);//TODO: play with these numbers to optimize
	printf("NumThreads=%d;numGrids=%d\n",numThreads,numGrids);
	dim3    grids(numGrids,1);
	dim3    threads(numThreads,1);
	for(int ii=0;ii<iterations;ii++)
	{
		//update ref points (nw blob locations fater EM iterations
		for (i=0 ; i<ref_nb   * dimsImage ; i++) ref[i]    = (float)rand() / (float)RAND_MAX;	
		
		
		HANDLE_ERROR(cudaMemcpyToSymbol(refCUDA,ref,ref_nb* dimsImage * sizeof(float)));//constant memory	
		knnKernel<<<grids,threads>>>(indCUDA,queryCUDA,ref_nb,query_nb);
		

		
		HANDLE_ERROR_KERNEL;
		HANDLE_ERROR(cudaMemcpy(ind,indCUDA,query_nb*maxGaussiansPerVoxel*sizeof(int),cudaMemcpyDeviceToHost));//retrieve indexes: memcopy is synchronous unless stated otherwise

		
		//test results
		//testKnnResults(ref,query,ind,ref_nb,query_nb,dimsImage,maxGaussiansPerVoxel);
	}


	// Display informations
	printf("Number of reference points      : %6d\n", ref_nb  );
	printf("Number of query points          : %6d\n", query_nb);
	printf("Dimension of points             : %4d\n", dimsImage     );
	printf("Number of neighbors to consider : %4d\n", maxGaussiansPerVoxel       );
	printf("Processing kNN search           :\n"                );
	// get stop time, and display the timing results
	HANDLE_ERROR( cudaEventRecord( stop, 0 ) );
	HANDLE_ERROR( cudaEventSynchronize( stop ) );
	float   elapsedTime;
	HANDLE_ERROR( cudaEventElapsedTime( &elapsedTime,start, stop ) );
	printf(" done in %f secs for %d iterations (%f secs per iteration)\n", elapsedTime/1000, iterations, elapsedTime/(iterations*1000));

	HANDLE_ERROR( cudaEventDestroy( start ) );
	HANDLE_ERROR( cudaEventDestroy( stop ) );

	HANDLE_ERROR( cudaFree( indCUDA ) );
	  
	  HANDLE_ERROR( cudaFree( queryCUDA ) );
	  
	
	// Destroy cuda event object and free memory
	free(ind);
	free(query);
	free(ref);
	return 0;
}

