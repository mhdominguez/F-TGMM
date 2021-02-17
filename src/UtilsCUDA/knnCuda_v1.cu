/*
 * knnCuda.cu
 *
 *  Created on: Jul 15, 2011
 *      Author: amatf
 */

#include "GMEMcommonCUDA.h"
#include "knnCuda_v1.h"
#include "external/book.h"
#include <iostream>
#include <algorithm>

#if defined(_WIN32) || defined(_WIN64)
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif

__constant__ float knn_v1_refCUDA[MAX_REF_POINTS*3];
__constant__ float knn_v1_scaleCUDA[dimsImage];

__device__ inline void findMaxPosition(float *distArray,float* minDist,int *pos)
{
	(*minDist)=distArray[0];
	(*pos)=0;
	for(int ii=1;ii<maxGaussiansPerVoxel;ii++) 
	{
		if((*minDist)<distArray[ii])
		{
			(*minDist)=distArray[ii];
			(*pos)=ii;
		}
	}
}

__global__ void __launch_bounds__(MAX_THREADS) knnKernel(int *indCUDA,float *queryCUDA,int ref_nb,long long int query_nb)
{
	// map from threadIdx/BlockIdx to pixel position
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	//int offset = blockDim.x * gridDim.x;

	if(tid>=query_nb) return;

	//int kMinusOne=maxGaussiansPerVoxel-1;
	float minDist[maxGaussiansPerVoxel];//to mantain distance for each index: since K is very small instead of a priority queue we keep a sorted array
	int indAux[maxGaussiansPerVoxel];
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


	int refIdx=-3;
	for(int jj=0;jj<maxGaussiansPerVoxel;jj++) minDist[jj]=1e32;//equivalent to infinity
	minDistThr=1e32;
	minPos=0;
	for(int ii=0;ii<ref_nb;ii++)
	{
		__syncthreads();//to access constant memory coherently
		refIdx+=3;

		distAux=(queryAux[0]-knn_v1_refCUDA[refIdx])*knn_v1_scaleCUDA[0];
		dist=distAux*distAux;
		if(dist>minDistThr) continue;
		distAux=(queryAux[1]-knn_v1_refCUDA[refIdx+1])*knn_v1_scaleCUDA[1];
		dist+=distAux*distAux;
		if(dist>minDistThr) continue;
		distAux=(queryAux[2]-knn_v1_refCUDA[refIdx+2])*knn_v1_scaleCUDA[2];
		dist+=distAux*distAux;
		if(dist>minDistThr) continue;

		//insert element" minimize memory exchanges
		minDist[minPos]=dist;
		indAux[minPos]=ii;
		findMaxPosition(minDist,&minDistThr,&minPos);
	}

	__syncthreads();//I need this to have coalescent memory access to inCUDA: speeds up the code by x4

	//copy indexes to global memory
	jj2=tid;
	for(int jj=0;jj<maxGaussiansPerVoxel;jj++)
	{
		//indCUDA[jj+jj2]=indAux[jj];
		indCUDA[jj2]=indAux[jj];
		jj2+=query_nb;
	}
}

//===========================================================================================
__global__ void __launch_bounds__(MAX_THREADS) knnKernelNoConstantMemory(int *indCUDA,float *queryCUDA,float* anchorCUDA,int ref_nb,long long int query_nb)
{
	// map from threadIdx/BlockIdx to pixel position
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	//int offset = blockDim.x * gridDim.x;

	if(tid>=query_nb) return;

	float minDist[maxGaussiansPerVoxel];//to mantain distance for each index: since K is very small instead of a priority queue we keep a sorted array
	int indAux[maxGaussiansPerVoxel];
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

	int refIdx=-3;
	for(int jj=0;jj<maxGaussiansPerVoxel;jj++) minDist[jj]=1e32;//equivalent to infinity
	minDistThr=1e32;
	minPos=0;
	for(int ii=0;ii<ref_nb;ii++)
	{
		__syncthreads();//to access constant memory coherently
		refIdx+=3;

		distAux=(queryAux[0]-anchorCUDA[refIdx])*knn_v1_scaleCUDA[0];
		dist=distAux*distAux;
		if(dist>minDistThr) continue;
		distAux=(queryAux[1]-anchorCUDA[refIdx+1])*knn_v1_scaleCUDA[1];
		dist+=distAux*distAux;
		if(dist>minDistThr) continue;
		distAux=(queryAux[2]-anchorCUDA[refIdx+2])*knn_v1_scaleCUDA[2];
		dist+=distAux*distAux;
		if(dist>minDistThr) continue;

		//insert element" minimize memory exchanges
		minDist[minPos]=dist;
		indAux[minPos]=ii;
		findMaxPosition(minDist,&minDistThr,&minPos);
	}

	__syncthreads();//I need this to have coalescent memory access to inCUDA: speeds up the code by x4

	//copy indexes to global memory
	jj2=tid;
	for(int jj=0;jj<maxGaussiansPerVoxel;jj++)
	{
		//indCUDA[jj+jj2]=indAux[jj];
		indCUDA[jj2]=indAux[jj];
		jj2+=query_nb;
	}
}
//===========================================================================================
__device__ inline void Comparator(float& keyA,int& valA,float& keyB,int& valB,unsigned int dir)
{
    float t;
    int v;
    if( (keyA > keyB) == dir ){
        t = keyA; keyA = keyB; keyB = t;
        v = valA; valA = valB; valB = v;
    }
}
//this kernel needs to be called with MAX_THREADS_CUDA
__global__ void __launch_bounds__(MAX_THREADS_CUDA) knnKernelSorting(int *indCUDA,float *queryCUDA,float* anchorCUDA,int ref_nb,long long int query_nb)
{
	
	//Shared memory storage for one or more short vectors
	__shared__ float s_key[MAX_THREADS_CUDA];//distance
	__shared__ int s_val[MAX_THREADS_CUDA];//index
	__shared__ int indAux[maxGaussiansPerVoxel];
	__shared__ float minDist[maxGaussiansPerVoxel];

	float x_n[dimsImage];
	unsigned int dir=0;//ascending order sorting

	// map from threadIdx/BlockIdx to pixel position
	long long int tid = blockIdx.x;
	long long int pos2;
	float dist,aux;

	int maxOffset=((ref_nb+MAX_THREADS_CUDA-1)/MAX_THREADS_CUDA)*MAX_THREADS_CUDA;

	while(tid<query_nb)
	{
		//load query point
		pos2=tid;
		x_n[0]=queryCUDA[pos2];
		pos2+=query_nb;
		x_n[1]=queryCUDA[pos2];
		pos2+=query_nb;
		x_n[2]=queryCUDA[pos2];

		for(int offset=threadIdx.x;offset<maxOffset;offset+=MAX_THREADS_CUDA)
		{
			//calculate distance
			if(offset<ref_nb)
			{
				aux=x_n[0]-anchorCUDA[offset];
				dist=aux*aux;
				offset+=ref_nb;
				aux=x_n[0]-anchorCUDA[offset];
				dist+=aux*aux;
				offset+=ref_nb;
				aux=x_n[2]-anchorCUDA[offset];
				dist+=aux*aux;
			}else{
				dist=1e32;
			}
			s_val[threadIdx.x]=offset;
			s_key[threadIdx.x]=dist;
			__syncthreads();

			//sort value and key in shared memory: bitonc search from Cuda SDK
			for(unsigned int size = 2; size < MAX_THREADS_CUDA; size <<= 1)
			{
				//Bitonic merge
				unsigned int ddd = dir ^ ( (threadIdx.x & (size / 2)) != 0 );
				for(unsigned int stride = size / 2; stride > 0; stride >>= 1)
				{
					__syncthreads();
					unsigned int pos = 2 * threadIdx.x - (threadIdx.x & (stride - 1));
					Comparator(s_key[pos +      0], s_val[pos +      0],s_key[pos + stride], s_val[pos + stride],ddd);
				}
			}

			//ddd == dir for the last bitonic merge step
			{
				for(unsigned int stride = MAX_THREADS_CUDA / 2; stride > 0; stride >>= 1){
					__syncthreads();
					unsigned int pos = 2 * threadIdx.x - (threadIdx.x & (stride - 1));
					Comparator(s_key[pos +      0], s_val[pos +      0],s_key[pos + stride], s_val[pos + stride],dir);
				}
			}
			__syncthreads();

			//merge this batch of distances with short sorted array
			if(offset<maxGaussiansPerVoxel)//we just need to copy in the first iteration
			{
				indAux[offset]=s_val[offset];
				minDist[offset]=s_key[offset];
			}else
			{

				if(threadIdx.x==0){ //merge two sorted arrays
					int ptr1=0;
					//int	ptr2=0;
					while(ptr1<maxGaussiansPerVoxel)
					{
							ptr1++;
							//TODO: finish this although the kernel is really slow (in comparison) even without this part
					}
				}
				
			}
			__syncthreads();
		}
		//copy indexes to global memory
		if(threadIdx.x<maxGaussiansPerVoxel)
		{
			indCUDA[tid+threadIdx.x*query_nb]=indAux[threadIdx.x];//not coalescence
		}
		//update pointer for next query_point to check
		tid+=gridDim.x;
		__syncthreads();

	}
}

__global__ void __launch_bounds__(MAX_THREADS) knnKernelSortedArray(int *indCUDA,float *queryCUDA,int ref_nb,long long int query_nb)
{
	// map from threadIdx/BlockIdx to pixel position
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	//int offset = blockDim.x * gridDim.x;

	if(tid>=query_nb) return;

	int kMinusOne=maxGaussiansPerVoxel-1;
	float minDist[maxGaussiansPerVoxel];//to mantain distance for each index: since K is very small instead of a priority queue we keep a sorted array
	int indAux[maxGaussiansPerVoxel];
	float queryAux[dimsImage];//TODO: I can probably hardcode dimsImage to improve performance (unroll loops)
	float minDistThr;
	float scaleAux[dimsImage];
	scaleAux[0]=knn_v1_scaleCUDA[0];scaleAux[1]=knn_v1_scaleCUDA[1];scaleAux[2]=knn_v1_scaleCUDA[2];

	float dist,distAux;
	int jj2;

	jj2=tid;
	/*texture mmemory
		queryAux[0]=tex1Dfetch(queryTexture,jj2);//stores query point to compare against all the references
		queryAux[1]=tex1Dfetch(queryTexture,jj2+1);
		queryAux[2]=tex1Dfetch(queryTexture,jj2+2);
	 */

	//global memory: organized as x_1,x_2,x_3,....,y_1,y_2,...,z_1,... to have coalescent access
	queryAux[0]=queryCUDA[jj2];
	jj2+=query_nb;
	queryAux[1]=queryCUDA[jj2];
	jj2+=query_nb;
	queryAux[2]=queryCUDA[jj2];

	int refIdx=-3;
	for(int jj=0;jj<maxGaussiansPerVoxel;jj++) minDist[jj]=1e32;//equivalent to infinity
	minDistThr=minDist[kMinusOne];
	for(int ii=0;ii<ref_nb;ii++)
	{
		__syncthreads();//to access constant memory coherently
		refIdx+=3;
		/*
			dist=0;
			for(int jj=0;jj<dimsImage;jj++) 
			{
				dist+=(queryAux[jj]-knn_v1_refCUDA[refIdx])*(queryAux[jj]-knn_v1_refCUDA[refIdx]);
				refIdx++;
			}
		 */

		distAux=(queryAux[0]-knn_v1_refCUDA[refIdx])*scaleAux[0];
		dist=distAux*distAux;
		if(dist>minDistThr) continue;
		distAux=(queryAux[1]-knn_v1_refCUDA[refIdx+1])*scaleAux[1];
		dist+=distAux*distAux;
		if(dist>minDistThr) continue;
		distAux=(queryAux[2]-knn_v1_refCUDA[refIdx+2])*scaleAux[2];
		dist+=distAux*distAux;
		if(dist>minDistThr) continue;

		//decide weather to insert this index or not
		for(jj2=kMinusOne-1;jj2>=0;jj2--)
		{
			if(dist>=minDist[jj2])
			{
				minDist[jj2+1]=dist;
				indAux[jj2+1]=ii;
				break;
			}
			minDist[jj2+1]=minDist[jj2];
			indAux[jj2+1]=indAux[jj2];
		}
		if(jj2==-1)//we need to insert the element at position zero
		{
			minDist[0]=dist;
			indAux[0]=ii;
		}
		minDistThr=minDist[kMinusOne];
	}

	__syncthreads();//I need this to have coalescent memory access to inCUDA: speeds up the code by x4

	//copy indexes to global memory
	jj2=tid;
	for(int jj=0;jj<maxGaussiansPerVoxel;jj++)
	{
		//indCUDA[jj+jj2]=indAux[jj];
		indCUDA[jj2]=indAux[jj];
		jj2+=query_nb;
	}
	//update pointer for next query_point to check
	//tid+=offset;

}


#if 0
//===================================================================================================
int allocateGPUMemoryForKnnCUDA(float *queryTemp,float **queryCUDA,int **indCUDA,long long int query_nb,float *scale)
{
	
	
	HANDLE_ERROR( cudaMalloc( (void**)&(*indCUDA), query_nb*maxGaussiansPerVoxel*sizeof(int) ) );
	HANDLE_ERROR( cudaMalloc( (void**)&(*queryCUDA), query_nb*dimsImage*sizeof(float) ) );
	// Copy image data to array
	HANDLE_ERROR(cudaMemcpy((*queryCUDA),queryTemp, dimsImage*query_nb*sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpyToSymbol(knn_v1_scaleCUDA,scale, dimsImage * sizeof(float)));//constant memory
	return 0;
}
#endif

#if 0
void setDeviceCUDA(int devCUDA)
{
//WE ASSUME qeuryCUDA AND indCUDA HAVE BEEN ALLOCATED ALREADY AND MEMORY TRANSFERRED TO THE GPU
	HANDLE_ERROR( cudaSetDevice( devCUDA ) );
}
#endif

#if 0
//===================================================================================================
void knnCUDA(int *ind,int *indCUDA,float *queryCUDA,float *refTemp,long long int query_nb,int ref_nb)
{


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
	//calculate number of threads and blocks
long long 	int numThreads=std::min((long long int)MAX_THREADS,query_nb);
long long 	int numGrids=std::min((long long int)MAX_BLOCKS,(query_nb+numThreads-1)/numThreads);
	
	//printf("NumThreads=%d;numGrids=%d\n",numThreads,numGrids);
	
	
	HANDLE_ERROR(cudaMemcpyToSymbol(knn_v1_refCUDA,refTemp,ref_nb* dimsImage * sizeof(float)));//constant memory	
	knnKernel<<<numGrids,numThreads>>>(indCUDA,queryCUDA,ref_nb,query_nb);
	HANDLE_ERROR_KERNEL;


	HANDLE_ERROR(cudaMemcpy(ind,indCUDA,query_nb*maxGaussiansPerVoxel*sizeof(int),cudaMemcpyDeviceToHost));//retrieve indexes: memcopy is synchronous unless stated otherwise

}
#endif 

//===================================================================================================
void knn_v1_CUDAinPlace(int *indCUDA,float *queryCUDA,float *refTemp,long long int query_nb,int ref_nb)
{

	
	if(dimsImage!=3)
	{
		printf("ERROR: dimsImage should be 3\n");
		exit(2);
	}
	//calculate number of threads and blocks
long long 	int numThreads=std::min((long long int)MAX_THREADS,query_nb);
long long 	int numGrids=std::min((long long int)MAX_BLOCKS,(query_nb+numThreads-1)/numThreads);
	
	if(MAX_REF_POINTS<ref_nb)//use global memory for anchorPoints (slower)
	{
		knnKernelNoConstantMemory<<<numGrids,numThreads>>>(indCUDA,queryCUDA,refTemp,ref_nb,query_nb);HANDLE_ERROR_KERNEL;	
	}else{//use constant memory for anchor points
		HANDLE_ERROR(cudaMemcpyToSymbol(knn_v1_refCUDA,refTemp,ref_nb* dimsImage * sizeof(float),0,cudaMemcpyDeviceToDevice));//constant memory	
		knnKernel<<<numGrids,numThreads>>>(indCUDA,queryCUDA,ref_nb,query_nb);HANDLE_ERROR_KERNEL;
	}
}


#if 0
//====================================================================================================
void deallocateGPUMemoryForKnnCUDA(float **queryCUDA,int **indCUDA)
{
	HANDLE_ERROR( cudaFree( *indCUDA ) );
	(*indCUDA)=NULL;
    HANDLE_ERROR( cudaFree( *queryCUDA ) );
    (*queryCUDA)=NULL;
}
#endif 

//==============================================================
void knn_v1_uploadScaleCUDA(float *scale)
{
	HANDLE_ERROR(cudaMemcpyToSymbol(knn_v1_scaleCUDA,scale, dimsImage * sizeof(float)));//constant memory
}
