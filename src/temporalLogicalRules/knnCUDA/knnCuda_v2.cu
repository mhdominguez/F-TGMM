/*
 * knnCuda.cu
 *
 *  Created on: Jul 15, 2011
 *      Author: amatf
 */

#include "knnCuda_v2.h"
#include "book.h"
#include <algorithm>

#if defined(_WIN32) || defined(_WIN64)
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif


//static const int MAX_REF_POINTS=5000;//we need to predefined this in order to store reference points as constant memory. Total memory needed is MAX_QUERY_POINTS*3*4 bytes. It can not be more than 5400!!!

#ifndef CUDA_MAX_SIZE_CONST //to protect agains teh same constant define in other places in the code
#define CUDA_MAX_SIZE_CONST
    #ifndef CUDA_CONSTANTS_FA
    #define CUDA_CONSTANTS_FA
        static const int MAX_THREADS=1024;//For Quadro4800->256;//to make sure we don't run out of registers; For TeslaC2070 -> 1024
        //static const int MAX_THREADS_CUDA=1024;//For Quadro4800->512;//certain kernels benefit from maximum number of threads
        static const int MAX_BLOCKS=65535;
    #endif
#endif



#ifndef DIMS_IMAGE_CONST //to protect agains teh same constant define in other places in the code
#define DIMS_IMAGE_CONST
static const int dimsImage = 3;//to be able to precompile code
#endif


__constant__ float knn_v2_scaleCUDA[dimsImage];


__device__ inline void findMaxPosition(float *distArray,float* minDist,int *pos, int KNN)
{
	(*minDist)=distArray[0];
	(*pos)=0;
	for(int ii=1;ii<KNN;ii++) 
	{
		if((*minDist)<distArray[ii])
		{
			(*minDist)=distArray[ii];
			(*pos)=ii;
		}
	}
}


//===========================================================================================
__global__ void __launch_bounds__(MAX_THREADS) knnKernelNoConstantMemory(int *indCUDA,float *distCUDA,float *queryCUDA,float* anchorCUDA,int ref_nb,long long int query_nb, int KNN)
{
	// map from threadIdx/BlockIdx to pixel position
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	//int offset = blockDim.x * gridDim.x;


	if(KNN > knn_v2_maxKNN) return;//code is not ready for this
	if(tid>=query_nb) return;

	
	//int kMinusOne=knn_v2_maxKNN-1;
	float minDist[knn_v2_maxKNN];//to mantain distance for each index: since K is very small instead of a priority queue we keep a sorted array
	int indAux[knn_v2_maxKNN];
	float queryAux[dimsImage];
	float minDistThr;
	
	float dist,distAux;
	int jj2,minPos;


	jj2=tid;

	//global memory: organized as x_1,x_2,x_3,....,y_1,y_2,...,z_1,... to have coalescent access
	queryAux[0]=queryCUDA[jj2];
	jj2+=query_nb;
	queryAux[1]=queryCUDA[jj2];
	jj2+=query_nb;
	queryAux[2]=queryCUDA[jj2];

	int refIdx;
	for(int jj=0;jj<KNN;jj++) minDist[jj]=1e32;//equivalent to infinity. Thus, we know this element has not been assigned
	minDistThr=1e32;
	minPos=0;
	for(int ii=0;ii<ref_nb;ii++)
	{
		//__syncthreads();//to access constant memory coherently (this was effective in CUDA 3.2)
		refIdx = ii;
		

		distAux=(queryAux[0]-anchorCUDA[refIdx])*knn_v2_scaleCUDA[0];
		dist=distAux*distAux;
		refIdx += ref_nb;
		if(dist>minDistThr) continue;
		distAux=(queryAux[1]-anchorCUDA[refIdx])*knn_v2_scaleCUDA[1];
		dist+=distAux*distAux;
		refIdx += ref_nb;
		if(dist>minDistThr) continue;
		distAux=(queryAux[2]-anchorCUDA[refIdx])*knn_v2_scaleCUDA[2];
		dist+=distAux*distAux;
		if(dist>minDistThr) continue;


		//insert element" minimize memory exchanges
		minDist[minPos]=dist;
		indAux[minPos]=ii;
		findMaxPosition(minDist,&minDistThr,&minPos, KNN);
	}

	
	__syncthreads();//I need this to have coalescent memory access to inCUDA: speeds up the code by x4

	//copy indexes to global memory
	jj2=tid;
	for(int jj=0;jj<KNN;jj++)
	{
		//indCUDA[jj+jj2]=indAux[jj];
		indCUDA[jj2]=indAux[jj];
		jj2+=query_nb;
	}
	//copy distance if requested by user
	if(distCUDA != NULL)
	{
		jj2=tid;
		for(int jj=0;jj<KNN;jj++)
		{
			//indCUDA[jj+jj2]=indAux[jj];
			distCUDA[jj2]=minDist[jj];
			jj2+=query_nb;
		}
	}

}


//=============================================================================================================
int knnCUDA_v2(int *ind,float* dist, float *query,float *ref,long long int query_nb,int ref_nb, int KNN, float* scale, int devCUDA)
{
	// Variables and parameters
	//float* ref;                 // Pointer to reference point array: order is cache friednly with the GPU
	//float* query;               // Pointer to query point array: order is x1,y1,z1,x2,y2,z2... to be cache friendly
	//int*   ind;                 // Pointer to index array: size query_nb*knn_v2_maxKNN. Again, order is GPU cache friendly.
	//float*   dist;              // Pointer to distance^2 array: size query_nb*knn_v2_maxKNN. Again, order is GPU cache friendly. If pointer is null, scaled euclidean distance to each nearest neighbor is not returned
	//int    ref_nb       // Reference point number
	//int    query_nb    // Query point number
	//float scale[dimsImage] //
	
	if(dimsImage!=3)
	{
		printf("ERROR: at knnCUDA: code is not ready for dimsImage other than 3\n");//TODO: change this to any dimensionality
		return 2;
	}
	
	if(ref_nb <= 0)//nothing to do. There are no possible assignments
	{
		if(dist != NULL)
		{
			for(long long int ii = 0; ii < query_nb*KNN; ii++)
				dist[ii] = 1e32f;//no assignments
		}
		return 0;
	}
	//CUDA variables
	int *indCUDA;
	float* queryCUDA;
	float *anchorCUDA;
	float *distCUDA = NULL;
	
	//set CUDA device
	HANDLE_ERROR( cudaSetDevice( devCUDA ) );
	
	
	// allocate memory on the GPU for the output: it will only be done once in the whole program
	HANDLE_ERROR( cudaMalloc( (void**)&indCUDA, query_nb*KNN*sizeof(int) ) );//should it be a texture memory?NO. It does not fit in Cuda2Darray but it fits in linear 1Dtexture, although it does not seems to bring benefits
	HANDLE_ERROR( cudaMalloc( (void**)&queryCUDA, query_nb*dimsImage*sizeof(float) ) );
	HANDLE_ERROR( cudaMalloc( (void**)&anchorCUDA, ref_nb*dimsImage*sizeof(float) ) );
	if( dist != NULL)
		HANDLE_ERROR( cudaMalloc( (void**)&distCUDA, query_nb*KNN*sizeof(float) ) );

	// Copy image data to array
	HANDLE_ERROR(cudaMemcpy(queryCUDA,query, dimsImage*query_nb*sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(anchorCUDA,ref, dimsImage*ref_nb*sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpyToSymbol(knn_v2_scaleCUDA,scale, dimsImage * sizeof(float)));//constant memory	



	//prepare to launch kernel
	int numThreads=min(MAX_THREADS,(int) query_nb);
	int numGrids=min(MAX_BLOCKS,(int) (query_nb+numThreads-1)/numThreads);//TODO: play with these numbers to optimize

	knnKernelNoConstantMemory<<<numGrids,numThreads>>>(indCUDA, distCUDA, queryCUDA,anchorCUDA,ref_nb,query_nb, KNN);HANDLE_ERROR_KERNEL;	
	

	//copy results back
	HANDLE_ERROR(cudaMemcpy(ind,indCUDA,query_nb*KNN*sizeof(int),cudaMemcpyDeviceToHost));//retrieve indexes: memcopy is synchronous unless stated otherwise
	if( distCUDA != NULL)
		HANDLE_ERROR(cudaMemcpy(dist,distCUDA,query_nb*KNN*sizeof(float),cudaMemcpyDeviceToHost));


	//free memory
	HANDLE_ERROR( cudaFree( indCUDA ) );
	HANDLE_ERROR( cudaFree( queryCUDA ) );
	HANDLE_ERROR( cudaFree( anchorCUDA ) );
	if( distCUDA != NULL)
		HANDLE_ERROR( cudaFree( distCUDA ) );
	return 0;
}

#if 0
//===================================================================================================
int allocateGPUMemoryForKnnCUDA_(float *queryTemp,float **queryCUDA,int **indCUDA,long long int query_nb,float *scale, int KNN)
{
	
	
	HANDLE_ERROR( cudaMalloc( (void**)&(*indCUDA), query_nb*KNN*sizeof(int) ) );
	HANDLE_ERROR( cudaMalloc( (void**)&(*queryCUDA), query_nb*dimsImage*sizeof(float) ) );
	// Copy image data to array
	HANDLE_ERROR(cudaMemcpy((*queryCUDA),queryTemp, dimsImage*query_nb*sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpyToSymbol(knn_v2_scaleCUDA,scale, dimsImage * sizeof(float)));//constant memory
	return 0;
}

void setDeviceCUDA_(int devCUDA)
{
//WE ASSUME qeuryCUDA AND indCUDA HAVE BEEN ALLOCATED ALREADY AND MEMORY TRANSFERRED TO THE GPU
	HANDLE_ERROR( cudaSetDevice( devCUDA ) );
}




//====================================================================================================
void deallocateGPUMemoryForKnnCUDA_(float **queryCUDA,int **indCUDA)
{
	HANDLE_ERROR( cudaFree( *indCUDA ) );
	(*indCUDA)=NULL;
    HANDLE_ERROR( cudaFree( *queryCUDA ) );
    (*queryCUDA)=NULL;
}
//==============================================================
void uploadknn_v2_ScaleCUDA_(float *scale)
{
	HANDLE_ERROR(cudaMemcpyToSymbol(knn_v2_scaleCUDA,scale, dimsImage * sizeof(float)));//constant memory
}

#endif