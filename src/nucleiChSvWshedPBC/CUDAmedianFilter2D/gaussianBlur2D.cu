/*
* Copyright (C) 2013 by  Fernando Amat
* See license.txt for full license and copyright notice.
*
* Authors: Fernando Amat 
*  gaussianBlur2D.cu
*
*  Created on: January 17th, 2013
*      Author: Fernando Amat
*
* \brief Code to calculate 2D median filter in CUDA using templates and different window sizes
*
*/

#include "gaussianBlur2D.h"
#include "book.h"
#include "cuda.h"
#include <iostream>
#include <math.h>
#include <cuda_runtime.h>

//#define PROFILE_CODE // uncomment this to time execution

__constant__ int imDimCUDA[dimsImageSlice];//image dimensions



//=====================================================================================
template<class imgType, int radius>
__global__ void __launch_bounds__(BLOCK_SIDE*BLOCK_SIDE) gaussianBlurCUDAkernel(imgType* imCUDAin, imgType* imCUDAout, unsigned int imSize, float *kernel)
{

	const int radiusSize = ( 1 + 2 * radius) * ( 1 + 2 * radius);
	//shared memory to copy global memory
	__shared__ imgType blockNeigh [BLOCK_SIDE * BLOCK_SIDE];//stores values for a whole block
	imgType imgNeigh[ radiusSize ];//store values for each thread. This is the reason why radius is a template parameter and not a function input variable. Here we have a chance that everything fits in register memory for small radiuses
	
	
	int offset_x = blockIdx.x * (BLOCK_SIDE - 2* radius) - radius + threadIdx.x;//upper left corner of the image to start loading into share memory (counting overlap to accomodate radius)
	int offset_y = blockIdx.y * (BLOCK_SIDE - 2* radius) - radius + threadIdx.y;//upper left corner of the image to start loading into share memory (counting overlap to accomodate radius)
	
	int tid = threadIdx.y * BLOCK_SIDE + threadIdx.x;

	//each thread loads one pixel into share memory (colescent access)
	int pos;
	if( offset_x < 0 || offset_y < 0 || offset_x >= imDimCUDA[0] || offset_y >= imDimCUDA[1] )//out of bounds
	{
		pos = -1;
		blockNeigh[ tid ] = 0;//for now we assume zeros outside image boundaries
	}else{
		pos = offset_x + offset_y  * imDimCUDA[0];
		blockNeigh[ tid ] = imCUDAin[pos];
	}

	__syncthreads();

	if( threadIdx.x < radius || threadIdx.x >= BLOCK_SIDE-radius || threadIdx.y < radius || threadIdx.y >= BLOCK_SIDE-radius)
		return;//these threads are not needed (kind of a waste, but it is OK);

	
	//operate on block: this part could be substituted by any other operation in a blokc if we want to apply a different filter than median		
	int pp, count = 0;
	for( int ii = -radius; ii <= radius; ii++)
	{
		pp = threadIdx.x -radius + BLOCK_SIDE * ( threadIdx.y + ii );//initial position for jj for loop		

		for( int jj = -radius; jj <= radius; jj++)
		{				
			imgNeigh[count++] = blockNeigh[pp++];
		}
	}

	//selection algorithm to find the k-th smallest number (k = (radiusSize - 1) /2 (http://en.wikipedia.org/wiki/Selection_algorithm)
	imgType temp;
	for ( int ii=0; ii < (1 + radiusSize) /2; ii++)
	{
		// Find position of minimum element
		pp = ii;//minIndex
		for ( int jj = ii+1; jj < radiusSize; jj++)
		{
			if (imgNeigh[jj] < imgNeigh[pp])
			{
				pp = jj;
			}
		}
		temp = imgNeigh[pp];
		imgNeigh[pp] = imgNeigh[ii];
		imgNeigh[ii] = temp;
	}	

	if( pos>=0 && pos< imSize )
		imCUDAout[ pos ] = imgNeigh[ (radiusSize - 1) /2 ];
};



//===========================================================================

template<class imgType>
int gaussianBlurCUDA(imgType* im,int* imDim,int radius,int devCUDA)
{
	HANDLE_ERROR( cudaSetDevice( devCUDA ) );

	if( radius > (int)(floor( (BLOCK_SIDE -1) / 2.0f) )  || 2 * radius >= BLOCK_SIDE)
	{
		std::cout<<"ERROR: at gaussianBlurCUDA: code is not ready for such a large radius. Maximum radius allowed is "<<(int)(floor( (BLOCK_SIDE - 1) / 2.0f) )<<std::endl;
		return 2;
	}


	imgType* imCUDAinput = NULL;
	imgType* imCUDAoutput = NULL;


	int imSize = imDim[0];
	for( int ii = 1; ii < dimsImageSlice; ii++)
		imSize *= imDim[ii];

	//allocate memory in CUDA (input and output)
	HANDLE_ERROR( cudaMalloc( (void**)&(imCUDAinput), imSize * sizeof(imgType) ) );
	HANDLE_ERROR( cudaMalloc( (void**)&(imCUDAoutput), imSize * sizeof(imgType) ) );

	//transfer input: image and image dimensions
	HANDLE_ERROR(cudaMemcpy(imCUDAinput, im, imSize * sizeof(imgType), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpyToSymbol(imDimCUDA,imDim, dimsImageSlice * sizeof(int)));//constant memory

	//run kernel	
	dim3 threads( BLOCK_SIDE, BLOCK_SIDE );
	int numBlocks[dimsImageSlice];
	for (int ii = 0 ; ii< dimsImageSlice; ii++)
		numBlocks[ii] = (int) (ceil( (float)(imDim[ii] + radius ) / (float)(BLOCK_SIDE - 2 * radius) ) );
	dim3 blocks(numBlocks[0], numBlocks[1]);//enough to cover all the image

	switch(radius)
	{
	case 0:
		//do nothing
		break;
	case 1:
		gaussianBlurCUDAkernel<imgType, 1> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 2:
		gaussianBlurCUDAkernel<imgType, 2> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 3:
		gaussianBlurCUDAkernel<imgType, 3> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 4:
		gaussianBlurCUDAkernel<imgType, 4> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 5:
		gaussianBlurCUDAkernel<imgType, 5> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 6:
		gaussianBlurCUDAkernel<imgType, 6> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 7:
		gaussianBlurCUDAkernel<imgType, 7> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 8:
		gaussianBlurCUDAkernel<imgType, 8> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 9:
		gaussianBlurCUDAkernel<imgType, 9> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 10:
		gaussianBlurCUDAkernel<imgType, 10> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	case 11:
		gaussianBlurCUDAkernel<imgType, 11> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
		break;
	default:
		std::cout<<"ERROR: at gaussianBlurCUDA: code is not ready for such a large radius." <<std::endl;//If I need it at any point, I could extend this up to (int)(floor( (BLOCK_SIDE -1) / 2.0f) )
		return 4;
	}
	//copy result to host
	HANDLE_ERROR(cudaMemcpy(im, imCUDAoutput, imSize * sizeof(imgType), cudaMemcpyDeviceToHost));

	//deallocate memory
	HANDLE_ERROR( cudaFree( imCUDAinput ) );
	HANDLE_ERROR( cudaFree( imCUDAoutput ) );

	return 0;
}

//declare all the possible types so template compiles properly
template int gaussianBlurCUDA<unsigned char>(unsigned char* im,int* imDim,int radius,int devCUDA);
template int gaussianBlurCUDA<unsigned short int>(unsigned short int* im,int* imDim,int radius,int devCUDA);
template int gaussianBlurCUDA<float>(float* im,int* imDim,int radius,int devCUDA);


//===========================================================================

template<class imgType>
int gaussianBlurCUDASliceBySlice(imgType* im,int* imDim,float sigma,int devCUDA)
{
	HANDLE_ERROR( cudaSetDevice( devCUDA ) );

#ifdef PROFILE_CODE
	cudaEvent_t start, stop;
	HANDLE_ERROR( cudaEventCreate(&start ) );
	HANDLE_ERROR( cudaEventCreate(&stop ) );
	HANDLE_ERROR( cudaEventRecord(start,0 ) );
#endif

	int radius = getKernelRadiusForSigmaCUDA(sigma);
	
	//fill kernel
	float kernel[kradius*2+1];
	gaussianBlurKernel(sigma, kradius*2+1, kernel);	


	imgType* imCUDAinput = NULL;
	imgType* imCUDAoutput = NULL;


	int imSize = imDim[0];
	for( int ii = 1; ii < dimsImageSlice; ii++)
		imSize *= imDim[ii];

	//allocate memory in CUDA (input and output)
	HANDLE_ERROR( cudaMalloc( (void**)&(imCUDAinput), imSize * sizeof(imgType) ) );
	HANDLE_ERROR( cudaMalloc( (void**)&(imCUDAoutput), imSize * sizeof(imgType) ) );

	//kernel parameters
	dim3 threads( BLOCK_SIDE, BLOCK_SIDE );
	int numBlocks[dimsImageSlice];
	for (int ii = 0 ; ii< dimsImageSlice; ii++)
		numBlocks[ii] = (int) (ceil( (float)(imDim[ii] + radius ) / (float)(BLOCK_SIDE - 2*radius) ) );
	dim3 blocks(numBlocks[0], numBlocks[1]);//enough to cover all the image

	//perform median filter slice by slice
	for( int slice = 0; slice < imDim[dimsImageSlice ]; slice++)
	{

		//transfer input: image and image dimensions
		HANDLE_ERROR(cudaMemcpy(imCUDAinput, im, imSize * sizeof(imgType), cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpyToSymbol(imDimCUDA,imDim, dimsImageSlice * sizeof(int)));//constant memory

		//run kernel			
		switch(radius)
		{
		case 0:
			//do nothing
			break;
		case 1:
			gaussianBlurCUDAkernel<imgType, 1> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
			break;
		case 3:
			gaussianBlurCUDAkernel<imgType, 3> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
			break;
		case 7:
			gaussianBlurCUDAkernel<imgType, 7> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
			break;
		case 15:
			gaussianBlurCUDAkernel<imgType, 15> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
			break;
		case 31:
			gaussianBlurCUDAkernel<imgType, 31> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
			break;
		case 63:
			gaussianBlurCUDAkernel<imgType, 63> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
			break;
		case 127:
			gaussianBlurCUDAkernel<imgType, 127> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, imSize);HANDLE_ERROR_KERNEL;
			break;
		default:
			std::cout<<"ERROR: at gaussianBlurCUDA: code is not ready for discrete radius "<<radius <<std::endl;//If I need it at any point, I could extend this up to (int)(floor( (BLOCK_SIDE -1) / 2.0f) )
			return 4;
		}
		//copy result to host
		HANDLE_ERROR(cudaMemcpy(im, imCUDAoutput, imSize * sizeof(imgType), cudaMemcpyDeviceToHost));

		im += imSize;//increment pointer to next slice
	}


	//deallocate memory
	HANDLE_ERROR( cudaFree( imCUDAinput ) );
	HANDLE_ERROR( cudaFree( imCUDAoutput ) );


#ifdef PROFILE_CODE	
	HANDLE_ERROR( cudaEventRecord(stop,0 ) );
	HANDLE_ERROR( cudaEventSynchronize(stop ) );

	float elapsedTime;
	HANDLE_ERROR( cudaEventElapsedTime(&elapsedTime, start,stop ) );
	printf ("Time for the CUDA function gaussianBlurCUDASliceBySlice: %f ms\n", elapsedTime);
	HANDLE_ERROR( cudaEventDestroy(start ) );
	HANDLE_ERROR( cudaEventDestroy(stop ) );
#endif
	return 0;
}

//declare all the possible types so template compiles properly
template int gaussianBlurCUDASliceBySlice<unsigned char>(unsigned char* im,int* imDim,int radius,int devCUDA);
template int gaussianBlurCUDASliceBySlice<unsigned short int>(unsigned short int* im,int* imDim,int radius,int devCUDA);
template int gaussianBlurCUDASliceBySlice<float>(float* im,int* imDim,int radius,int devCUDA);

int getKernelRadiusForSigmaCUDA(float sigma) {
	int size = int(ceilf(sigma * 3)); //3 sigmas plus/minus should be decent
	if (size <= 1) {
		return 1;
	} else if ( size <= 3 ) {
		return 3;  
	} else if ( size <= 7 ) {
		return 7;
	} else if ( size <= 15 ) {
		return 15;
	} else if ( size <= 32 ) {
		return 31;
	} else if ( size <= 63 ) {
		return 63;		
	} else { //255 is max Kernel size
		return 127;
	}
}

void gaussianBlurKernel(float sigma, int size, float* kernel)
{
	float sigma2 = sigma * sigma;
	int middle = size / 2;
	for (int i = 0; i < size; ++i)
	{
		float distance = float (middle - i);
		float distance2 = distance * distance;
		float s = 1.0f / (sigma * SQRT_2PI * expf(-distance2 / (2.0f * sigma2));
		kernel[i] = s;
	}
}
