//Martin H. Dominguez
//2021 Gladstone Institutes

#include "gaussianBlur2D.h"
#include <iostream>
#include <math.h>
#include <thread>
#include <mutex>
#include <cstring>
#include <vector>

//#define PROFILE_CODE // uncomment this to time execution

int sliceId;//onlye one thread accesses a single slice
std::mutex				g_lockblockId;//so each worker reads a unique slice

using namespace std;

const float SQRT_2PI = sqrtf(2.0f*PI_F));

//=======================================================================
template<class imgType>
void nvsepgold_gaussianBlur2D(const imgType* im,const int* imDim,int kradius, float *kernel, imgType *imOut)
{
    int x, y, k, d;
    imgType sum;
	
	//x dimension
	for(y = 0; y < imDim[1]; y++)
	{
		for(x = 0; x < imDim[0]; x++)
		{
			sum = 0;
			for(k = -kradius; k <= kradius; k++){
				d = x + k;
				if(d >= 0 && d < imDim[0])
					sum += im[y * imDim[0] + d] * kernel[kradius - k];
			}
			imOut[y * imDim[0] + x] = sum;
		}
	}
	
	for(y = 0; y < imDim[1]; y++)
	{
		for(x = 0; x < imDim[0]; x++)
		{
			sum = 0;
			for(k = -kradius; k <= kradius; k++){
				d = y + k;
				if(d >= 0 && d < imDim[1])
					sum += imOut[d * imDim[0] + x] * kernel[kradius - k];
			}
			imOut[y * imDim[0] + x] = sum;
		}
	}
}

//declare all the possible types so template compiles properly
template void nvsepgold_gaussianBlur2D<unsigned char>(const unsigned char* im, const int* imDim, int kradius, float *kernel, unsigned char *imOut);
template void nvsepgold_gaussianBlur2D<unsigned short int>(const unsigned short int* im, const int* imDim, int kradius, float *kernel, unsigned short int *imOut);
template void nvsepgold_gaussianBlur2D<float>(const float* im, const int* imDim, int kradius, float *kernel, float *imOut);


//===============================================
template<class imgType>
void gaussianBlur2Dworker(imgType* im, const int* imDim, int kradius, float *kernel)
{
	uint64_t sliceSize = imDim[0];
	for (int ii = 1; ii < dimsImageSlice; ii++)
		sliceSize *= imDim[ii];

	int sliceId_t;
	const int numSlices = imDim[dimsImageSlice];
	imgType* slicePtr;
	uint64_t imOffset;

	imgType* imBuffer = new imgType[sliceSize];

	while (1)
	{
		//get the blockId resource
		std::unique_lock<std::mutex> locker(g_lockblockId);//exception safe
		sliceId_t = sliceId;
		sliceId++;
		locker.unlock();

		//check if we have more blocks
		if (sliceId_t >= numSlices)
			break;

		imOffset = sliceSize * (uint64_t)sliceId_t;

		slicePtr = &(im[imOffset]);

		memcpy(imBuffer, slicePtr, sizeof(imgType)* sliceSize);

		nvsepgold_gaussianBlur2D(imBuffer, imDim, kradius, slicePtr);

	}

	delete[] imBuffer;
}

template void gaussianBlur2Dworker<unsigned char>(unsigned char* im, const int* imDim, int kradius, float *kernel);
template void gaussianBlur2Dworker<unsigned short int>(unsigned short int* im, const int* imDim, int kradius, float *kernel);
template void gaussianBlur2Dworker<float>(float* im, const  int* imDim, int kradius, float *kernel);
//===========================================================================

template<class imgType>
void gaussianBlur2DSliceBySlice(imgType* im,const  int* imDim, int kradius, float sigma)
{
	int numThreads = std::thread::hardware_concurrency();//all the cores available

	int sliceSize = imDim[0];
	for (int ii = 1; ii < dimsImageSlice; ii++)
		sliceSize *= imDim[ii];
	
	sliceId = 0;//initialize count
	
	//fill kernel
	float kernel[kradius*2+1];
	gaussianBlurKernel(sigma, kradius*2+1, kernel);

	// start the working threads
	std::vector<std::thread> threads;
	std::vector<int> errFlagVec(numThreads, 0);
	for (int i = 0; i < numThreads; ++i)
	{
		threads.push_back(std::thread(gaussianBlur2Dworker<imgType>, im, imDim, kradius));
	}

	//wait for the workers to finish
	for (auto& t : threads)
		t.join();	
}

//declare all the possible types so template compiles properly
template void gaussianBlur2DSliceBySlice<unsigned char>(unsigned char* im, const  int* imDim, int kradius, float *kernel);
template void gaussianBlur2DSliceBySlice<unsigned short int>(unsigned short int* im, const  int* imDim, int kradius, float *kernel);
template void gaussianBlur2DSliceBySlice<float>(float* im, const int* imDim, int kradius), float *kernel;



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

