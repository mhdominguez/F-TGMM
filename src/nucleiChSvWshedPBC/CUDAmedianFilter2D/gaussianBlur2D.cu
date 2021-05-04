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

template<int radius>
__constant__ float kernelCUDA[ radius*2+1 ];

// how many threads per block in x (total num threads: x*y)
#define	ROWS_BLOCKDIM_X 16
// how many threads per block in y
#define	ROWS_BLOCKDIM_Y 16
// how many pixels in x are convolved by each thread
#define	ROWS_RESULT_STEPS 8
// these are the border pixels (loaded to support the kernel width for processing)
// the effective border width is ROWS_HALO_STEPS * ROWS_BLOCKDIM_X, which has to be
// larger or equal to the kernel radius to work
#define	ROWS_HALO_STEPS 1

#define	COLUMNS_BLOCKDIM_X 16
#define	COLUMNS_BLOCKDIM_Y 16
#define	COLUMNS_RESULT_STEPS 8
#define	COLUMNS_HALO_STEPS 1

#define	DEPTH_BLOCKDIM_X 16
#define	DEPTH_BLOCKDIM_Z 16
#define	DEPTH_RESULT_STEPS 8
#define	DEPTH_HALO_STEPS 1

static const float SQRT_2PI = sqrt(2*M_PI);


//=====================================================================================
template<class imgType, int radius>
__global__ void __launch_bounds__(MAX_THREADS_CUDA) gaussianBlurCUDAkernel(imgType* imCUDAin, imgType* imCUDAout, float *kernel, unsigned int imSize, bool convolve_y )
{
	//shared memory to copy global memory
	__shared__ double blockNeigh [MAX_THREADS_CUDA];//stores values for a whole block

	double sum;
	int offset_x;
	int offset_y;
	int tid;
	
	if (convolve_y) //y-dimension convolution
	{
		offset_y = blockIdx.y * (MAX_THREADS_CUDA - (2* radius+1)) -radius + threadIdx.x;
		offset_x = blockIdx.x;
		tid = threadIdx.x;
	}
	else //x-dimension convolution
	{
		offset_x = blockIdx.x * (MAX_THREADS_CUDA - (2* radius+1)) -radius + threadIdx.x;
		offset_y = blockIdx.y;
		tid = threadIdx.x;
	}

	//each thread loads one pixel into share memory (colescent access)
	int pos;
	if( offset_x < 0 || offset_y < 0 || offset_x >= imDimCUDA[0] || offset_y >= imDimCUDA[1] )//out of bounds
	{
		pos = -1;
		blockNeigh[ tid ] = 0;//zeros outside image boundaries
	}else{
		pos = offset_x + offset_y  * imDimCUDA[0];
		blockNeigh[ tid ] = (double)imCUDAin[pos];
	}
	__syncthreads();
	
	if( tid < radius || tid >= MAX_THREADS_CUDA-radius )
		return;//these threads are not needed (kind of a waste, but it is OK);	

	//here, we actually calculate the one-dimension convolution for this pixel
	int d;
	for( int k = -radius; k <= radius; k++)
	{
		d = tid + k;
		if(d >= 0 && d < MAX_THREADS_CUDA)
		{
			sum += (double)(blockNeigh[d] * kernelCUDA<radius>[radius - k]);
			//if (blockNeigh[d]>0)
			//	printf("found non-zero blockNeigh value at %d for {%d,%d}: %d\n",d,offset_x, offset_y, (int)blockNeigh[d] );
		}
	}

	if( pos>=0 && pos< imSize )
		imCUDAout[ pos ] = (imgType)(sum);
		//printf("%d: {%d, %d}/%d: %d/%d.\n",tid, offset_x, offset_y, pos, (int)imCUDAout[ pos ], (int)blockNeigh[tid]);
	
};



//===========================================================================

template<class imgType>
int gaussianBlurCUDA(imgType* im,int* imDim,float sigma,int devCUDA)
{
	HANDLE_ERROR( cudaSetDevice( devCUDA ) );

	int kradius = getKernelRadiusForSigmaCUDA(sigma);
	
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

	//transfer input: image and image dimensions
	HANDLE_ERROR(cudaMemcpy(imCUDAinput, im, imSize * sizeof(imgType), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpyToSymbol(imDimCUDA,imDim, dimsImageSlice * sizeof(int)));//constant memory
	
	//run kernel	
	dim3 threads( BLOCK_SIDE, BLOCK_SIDE );
	int numBlocks[dimsImageSlice];
	for (int ii = 0 ; ii< dimsImageSlice; ii++)
		numBlocks[ii] = (int) (ceil( (float)(imDim[ii] + kradius ) / (float)(BLOCK_SIDE - 2 * kradius) ) );
	dim3 blocks(numBlocks[0], numBlocks[1]);//enough to cover all the image

	switch(kradius)
	{
	case 0:
		//do nothing
		break;
		case 1:
			cudaMemcpyToSymbol(  kernelCUDA<1>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 1> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false );HANDLE_ERROR_KERNEL;
			break;
		case 3:
			cudaMemcpyToSymbol(  kernelCUDA<3>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 3> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			break;
		case 7:
			cudaMemcpyToSymbol(  kernelCUDA<7>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 7> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			break;
		case 15:
			cudaMemcpyToSymbol(  kernelCUDA<15>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 15> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			break;
		case 31:
			cudaMemcpyToSymbol(  kernelCUDA<31>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 31> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			break;
		case 63:
			cudaMemcpyToSymbol(  kernelCUDA<63>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 63> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			break;
		case 127:
			cudaMemcpyToSymbol(  kernelCUDA<127>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 127> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			break;
	default:
		std::cout<<"ERROR: at gaussianBlurCUDA: code is not ready for discrete radius "<<kradius <<std::endl;//If I need it at any point, I could extend this up to (int)(floor( (BLOCK_SIDE -1) / 2.0f) )
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
template int gaussianBlurCUDA<unsigned char>(unsigned char* im,int* imDim,float sigma,int devCUDA);
template int gaussianBlurCUDA<unsigned short int>(unsigned short int* im,int* imDim,float sigma,int devCUDA);
template int gaussianBlurCUDA<float>(float* im,int* imDim,float sigma,int devCUDA);


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

	int kradius = getKernelRadiusForSigmaCUDA(sigma);
	
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

	//copy imDim as constant
	HANDLE_ERROR(cudaMemcpyToSymbol(imDimCUDA,imDim, dimsImageSlice * sizeof(int)));//constant memory	
	
	dim3 threads( MAX_THREADS_CUDA ); //we will go line-by-line or column-by-column, or shorter than that, as dictated by max threads
	dim3 blocks;

	
	//perform separable convolution slice by slice
	for( int slice = 0; slice < imDim[dimsImageSlice ]; slice++)
	{
		//transfer input: image and image dimensions
		HANDLE_ERROR(cudaMemcpy(imCUDAinput, im, imSize * sizeof(imgType), cudaMemcpyHostToDevice));
		
		//run kernel			
		switch(kradius)
		{
		case 0:
			//do nothing
			break;
		case 1:
			blocks.x = (int) (ceil( (float)(imDim[0]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			blocks.y = imDim[1];
			if ( slice == 0 )
				cudaMemcpyToSymbol(  kernelCUDA<1>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 1> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			blocks.x = imDim[0];
			blocks.y = (int) (ceil( (float)(imDim[1]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			gaussianBlurCUDAkernel<imgType, 1> <<<blocks, threads>>>(imCUDAoutput, imCUDAinput, kernel, imSize, true);HANDLE_ERROR_KERNEL;
			break;
		case 3:
			blocks.x = (int) (ceil( (float)(imDim[0]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			blocks.y = imDim[1];
			if ( slice == 0 )
				cudaMemcpyToSymbol(  kernelCUDA<3>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 3> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			blocks.x = imDim[0];
			blocks.y = (int) (ceil( (float)(imDim[1]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			gaussianBlurCUDAkernel<imgType, 3> <<<blocks, threads>>>(imCUDAoutput, imCUDAinput, kernel, imSize, true);HANDLE_ERROR_KERNEL;
			break;
		case 7:
			blocks.x = (int) (ceil( (float)(imDim[0]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			blocks.y = imDim[1];
			if ( slice == 0 )
				cudaMemcpyToSymbol(  kernelCUDA<7>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 7> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			blocks.x = imDim[0];
			blocks.y = (int) (ceil( (float)(imDim[1]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			gaussianBlurCUDAkernel<imgType, 7> <<<blocks, threads>>>(imCUDAoutput, imCUDAinput, kernel, imSize, true);HANDLE_ERROR_KERNEL;
			break;
		case 15:
			blocks.x = (int) (ceil( (float)(imDim[0]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			blocks.y = imDim[1];	
			if ( slice == 0 )
				cudaMemcpyToSymbol(  kernelCUDA<15>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 15> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			blocks.x = imDim[0];
			blocks.y = (int) (ceil( (float)(imDim[1]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			gaussianBlurCUDAkernel<imgType, 15> <<<blocks, threads>>>(imCUDAoutput, imCUDAinput, kernel, imSize, true);HANDLE_ERROR_KERNEL;
			break;
		case 31:
			blocks.x = (int) (ceil( (float)(imDim[0]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			blocks.y = imDim[1];	
			if ( slice == 0 )
				cudaMemcpyToSymbol(  kernelCUDA<31>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 31> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			blocks.x = imDim[0];
			blocks.y = (int) (ceil( (float)(imDim[1]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			gaussianBlurCUDAkernel<imgType, 31> <<<blocks, threads>>>(imCUDAoutput, imCUDAinput, kernel, imSize, true);HANDLE_ERROR_KERNEL;
			break;
		case 63:
			blocks.x = (int) (ceil( (float)(imDim[0]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			blocks.y = imDim[1];			
			if ( slice == 0 )
				cudaMemcpyToSymbol(  kernelCUDA<63>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 63> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			blocks.x = imDim[0];
			blocks.y = (int) (ceil( (float)(imDim[1]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			gaussianBlurCUDAkernel<imgType, 63> <<<blocks, threads>>>(imCUDAoutput, imCUDAinput, kernel, imSize, true);HANDLE_ERROR_KERNEL;
			break;
		case 127:
			blocks.x = (int) (ceil( (float)(imDim[0]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			blocks.y = imDim[1];
			if ( slice == 0 )			
				cudaMemcpyToSymbol(  kernelCUDA<127>, kernel, (kradius*2+1) * sizeof(float) );
			gaussianBlurCUDAkernel<imgType, 127> <<<blocks, threads>>>(imCUDAinput, imCUDAoutput, kernel, imSize, false);HANDLE_ERROR_KERNEL;
			blocks.x = imDim[0];
			blocks.y = (int) (ceil( (float)(imDim[1]+kradius) / (float)(MAX_THREADS_CUDA - (2 * kradius+1)) ) );
			gaussianBlurCUDAkernel<imgType, 127> <<<blocks, threads>>>(imCUDAoutput, imCUDAinput, kernel, imSize, true);HANDLE_ERROR_KERNEL;
			break;
		default:
			std::cout<<"ERROR: at gaussianBlurCUDA: code is not ready for discrete radius "<<kradius <<std::endl;//If I need it at any point, I could extend this up to (int)(floor( (BLOCK_SIDE -1) / 2.0f) )
			return 4;
		}
		//copy result to host
		HANDLE_ERROR(cudaMemcpy(im, imCUDAinput, imSize * sizeof(imgType), cudaMemcpyDeviceToHost));
	
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
template int gaussianBlurCUDASliceBySlice<unsigned char>(unsigned char* im,int* imDim,float sigma,int devCUDA);
template int gaussianBlurCUDASliceBySlice<unsigned short int>(unsigned short int* im,int* imDim,float sigma,int devCUDA);
template int gaussianBlurCUDASliceBySlice<float>(float* im,int* imDim,float sigma,int devCUDA);

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
	float sigma2 = 2.0f * sigma * sigma;
	float sigma1 = 1.0f / sqrt( M_PI * sigma2 ); 
	int middle = size / 2;
	//std::cout << "GB kernel: ";
	for (int i = 0; i < size; ++i)
	{
		float distance = float (middle - i);
		//float distance2 = distance * distance;
		//float s = 1.0f / (sigma * sqrtf(2.0f * PI_F)) * expf(-distance2 / (2.0f * sigma2));
		kernel[i] = sigma1 * exp(-(distance * distance) / sigma2);
		//std::cout << i << ":" << kernel[i] << ",";
	}
	//std::cout << std::endl;
}
