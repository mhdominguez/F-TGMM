//Martin H. Dominguez
//2021 Gladstone Institutes


#include <iostream>
#include <math.h>
#include "gaussianBlur2D.h"
#include <thread>
#include <mutex>
#include <cstring>
#include <vector>

//#define PROFILE_CODE // uncomment this to time execution

int sliceId_gB2D;//onlye one thread accesses a single slice
std::mutex				g_lockblockId_gB2D;//so each worker reads a unique slice

using namespace std;

static const float SQRT_2PI = sqrt(2.0f*M_PI);

//=======================================================================
template<class imgType>
void nvsepgold_gaussianBlur2D(imgType* im,const int* imDim,int kradius, float *kernel, imgType *imOut, uint64_t sliceSize)
{
    unsigned int x,y,z,y_pos;
	int k,d;	
    //imgType sum;
	double sum; // can't do blur on integers
	
	//x dimension
	for(y_pos=0; y_pos < sliceSize; ) //do nothing to increment, will increment within nested for loop
	{
		for(x = 0; x < imDim[0]; x++)
		{
			sum = 0;
			for(k = -kradius; k <= kradius; k++){
				d = x + k;
				if(d >= 0 && d < imDim[0])
					sum += im[y_pos + k] * kernel[kradius - k];
			}
			imOut[y_pos] = (imgType)sum;
			y_pos++;
		}
	}
	
	/*for(y = 0; y < imDim[1]; y++)
	{
		for(x = 0; x < imDim[0]; x++)
		{
			sum = 0;
			for(k = -kradius; k <= kradius; k++){
				d = x + k;
				if(d >= 0 && d < imDim[0])
					sum += im[y * imDim[0] + d] * kernel[kradius - k];
			}
			imOut[y * imDim[0] + x] = (imgType)sum;
		}
	}*/
	
	//save our work, then do other dimension
	memcpy(im, imOut, sizeof(imgType)* sliceSize);
	
	//y dimension
	for(y = 0, y_pos=0; y < imDim[1]; y++)
	{
		for(x = 0; x < imDim[0]; x++)
		{
			sum = 0;
			for(k = -kradius; k <= kradius; k++){
				d = y + k;
				if(d >= 0 && d < imDim[1])
					sum += im[d * imDim[0] + x] * kernel[kradius - k];
			}
			imOut[y_pos] = (imgType)sum;
			y_pos++;
		}
	}
	
	/*for(y = 0; y < imDim[1]; y++)
	{
		for(x = 0; x < imDim[0]; x++)
		{
			sum = 0;
			for(k = -kradius; k <= kradius; k++){
				d = y + k;
				if(d >= 0 && d < imDim[1])
					sum += im[d * imDim[0] + x] * kernel[kradius - k];
			}
			imOut[y * imDim[0] + x] = sum;
		}
	}*/
}

//declare all the possible types so template compiles properly
template void nvsepgold_gaussianBlur2D<unsigned char>( unsigned char* im, const int* imDim, int kradius, float *kernel, unsigned char *imOut, uint64_t sliceSize);
template void nvsepgold_gaussianBlur2D<unsigned short int>( unsigned short int* im, const int* imDim, int kradius, float *kernel, unsigned short int *imOut, uint64_t sliceSize);
template void nvsepgold_gaussianBlur2D<float>( float* im, const int* imDim, int kradius, float *kernel, float *imOut, uint64_t sliceSize);


//===============================================
template<class imgType>
void gaussianBlur2Dworker(imgType* im, const int* imDim, int kradius, float *kernel)
{
	uint64_t sliceSize = imDim[0];
	for (int ii = 1; ii < dimsImageSlice; ii++)
		sliceSize *= imDim[ii];

	int sliceId_gB2D_t;
	const int numSlices = imDim[dimsImageSlice];
	imgType* slicePtr;
	uint64_t imOffset;

	imgType* imBuffer = new imgType[sliceSize];

	while (1)
	{
		//get the blockId resource
		std::unique_lock<std::mutex> locker(g_lockblockId_gB2D);//exception safe
		sliceId_gB2D_t = sliceId_gB2D;
		sliceId_gB2D++;
		locker.unlock();

		//check if we have more blocks
		if (sliceId_gB2D_t >= numSlices)
			break;

		imOffset = sliceSize * (uint64_t)sliceId_gB2D_t;

		slicePtr = &(im[imOffset]);

		memcpy(imBuffer, slicePtr, sizeof(imgType)* sliceSize);

		nvsepgold_gaussianBlur2D(imBuffer, imDim, kradius, kernel, slicePtr, sliceSize);

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
	
	sliceId_gB2D = 0;//initialize count
	
	//fill kernel
	float *kernel = new float[kradius*2+1];
	gaussianBlurKernel(sigma, kradius*2+1, kernel);
	
	// start the working threads
	std::vector<std::thread> threads;
	std::vector<int> errFlagVec(numThreads, 0);
	for (int i = 0; i < numThreads; ++i)
	{
		threads.push_back(std::thread(gaussianBlur2Dworker<imgType>, im, imDim, kradius, kernel));
	}

	//wait for the workers to finish
	for (auto& t : threads)
		t.join();
	
	delete kernel;
}

//declare all the possible types so template compiles properly
template void gaussianBlur2DSliceBySlice<unsigned char>(unsigned char* im, const  int* imDim, int kradius, float sigma);
template void gaussianBlur2DSliceBySlice<unsigned short int>(unsigned short int* im, const  int* imDim, int kradius, float sigma);
template void gaussianBlur2DSliceBySlice<float>(float* im, const int* imDim, int kradius, float sigma);


void gaussianBlurKernel(float sigma, int size, float* kernel)
{
	float sigma2 = 2.0f * sigma * sigma;
	float sigma1 = 1.0f / sqrt( M_PI * sigma2 ); 
	int middle = size / 2;
	//cout << "GB kernel: ";
	for (int i = 0; i < size; ++i)
	{
		float distance = float (middle - i);
		//float distance2 = distance * distance;
		//float s = 1.0f / (sigma * sqrtf(2.0f * PI_F)) * expf(-distance2 / (2.0f * sigma2));
		kernel[i] = sigma1 * exp(-(distance * distance) / sigma2);
		//cout << i << ":" << kernel[i] << ",";
	}
	//cout << endl;
}


