/*
* Copyright (C) 2014 by  Fernando Amat
* See license.txt for full license and copyright notice.
*
* Authors: Fernando Amat 
*  medianFilter2D.cpp
*
*  Created on: October 17th, 2014
*      Author: Fernando Amat
*
* \brief Code to calculate 2D median filter CPU using multithread for different slices
*
*/

#include "medianFilter2D.h"
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

//=======================================================================
/*
* The following code is public domain.
* Algorithm by Torben Mogensen, implementation by N. Devillard.
* This code in public domain.
*/

//From http://ndevilla.free.fr/median/median/src/torben.c
//Not the fastest but it does not modify input array
//For the fastest use sorting networks for small predefined size kernels http://ndevilla.free.fr/median/median/src/optmed.c
//What is worst: to read all the block for every pixel? (instead of just updating one column) or to have to use torben?

template<class elem_type>
elem_type medianBuffer(const elem_type* m, const int n)
{
	int         i, less, greater, equal;
	elem_type  min, max, guess, maxltguess, mingtguess;

	min = max = m[0];
	for (i = 1; i<n; i++) {
		if (m[i]<min) min = m[i];
		if (m[i]>max) max = m[i];
	}

	while (1) {
		guess = (min + max) / 2;
		less = 0; greater = 0; equal = 0;
		maxltguess = min;
		mingtguess = max;
		for (i = 0; i<n; i++) 
		{
			if (m[i]<guess) {
				less++;
				if (m[i]>maxltguess) 
					maxltguess = m[i];
			}
			else if (m[i]>guess) {
				greater++;
				if (m[i]<mingtguess) 
					mingtguess = m[i];
			}
			else equal++;
		}
		if (less <= (n + 1) / 2 && greater <= (n + 1) / 2) 
			break;
		else if (less>greater) 
			max = maxltguess;
		else 
			min = mingtguess;
	}
	if (less >= (n + 1) / 2) return maxltguess;
	else if (less + equal >= (n + 1) / 2) return guess;
	else return mingtguess;
}
//===========================================================================

template<class imgType>
void amatf_medianFilter2D(const imgType* im,const int* imDim,int radius, imgType *imOut)
{
		
	int imSize = imDim[0];
	for( int ii = 1; ii < dimsImageSlice; ii++)
		imSize *= imDim[ii];

	const int sliceOffset = imDim[0];
	const int win = 2 * radius + 1;
	const int blockSize = win * win;

	imgType* imBlock = new imgType[blockSize];//to store elements for the current element
	


	int pos = 0;
	int imBlockIdx;//to keep track of the circular index for the image block

	//initial pass for the first radius rows
	for (int jj = 0; jj < radius; jj++)
	{
		//reset imBlock: we load image data one column at a time (so the first radius column are zero		
		imBlockIdx = win * radius;
		memset(imBlock, 0, sizeof(imgType)*imBlockIdx);//for now we assume zeros outsize image borders
		for (int aa = 0; aa < radius; aa++)//the lsat colum will be done later
		{
			for (int bb = -radius; bb < -jj; bb++)
			{
				imBlock[imBlockIdx++] = 0;//for now we assume zeros outsize image borders
			}

			for (int bb = -jj; bb <= radius; bb++)
			{
				imBlock[imBlockIdx++] = im[pos + aa + sliceOffset * bb];
			}
		}

		//advance the iwndow in the x-direction one column at a time
		for (int ii = 0; ii < imDim[0] - radius; ii++)
		{
			//update imBlock with the new column
			for (int kk = -radius; kk < -jj; kk++)
			{
				imBlock[imBlockIdx % blockSize] = 0;//for now we assume zeros outsize image borders
				imBlockIdx++;
			}
			int aux = pos + radius - jj * sliceOffset;
			for (int kk = -jj; kk <= radius; kk++)
			{
				imBlock[imBlockIdx % blockSize] = im[aux];
				aux += sliceOffset;//accessing im[pos + radius + kk * sliceOffset];
				imBlockIdx++;
			}

			//calculate median element
			imOut[pos] = medianBuffer(imBlock, blockSize);
			//next element
			pos++;
		}

		//do the final elements where block edge is in the border
		for (int ii = imDim[0] - radius; ii < imDim[0]; ii++)
		{
			//update imBlock with the new column
			for (int kk = -radius; kk <= radius; kk++)
			{
				imBlock[imBlockIdx % blockSize] = 0; //for now we assume zeros outsize image borders
				imBlockIdx++;
			}

			//calculate median element
			imOut[pos] = medianBuffer(imBlock, blockSize);
			//next element
			pos++;
		}
	}

	//main loop when we are in rows between [radius, end-radius]
	for (int jj = radius; jj < imDim[1]-radius; jj++)
	{	
		//reset imBlock: we load image data one column at a time (so the first radius column are zero		
		imBlockIdx = win * radius;
		memset(imBlock, 0, sizeof(imgType)*imBlockIdx);//for now we assume zeros outsize image borders
		for (int aa = 0; aa < radius; aa++)//the lsat colum will be done later
		{			
			for (int bb = -radius; bb <= radius; bb++)
			{
				imBlock[imBlockIdx++] = im[pos + aa + sliceOffset * bb];
			}
		}

		//advance the iwndow in the x-direction one column at a time
		for (int ii = 0; ii < imDim[0]-radius; ii++)
		{
			//update imBlock with the new column
			int aux = pos + radius - radius * sliceOffset;
			for (int kk = -radius; kk <= radius; kk++)
			{
				imBlock[imBlockIdx % blockSize] = im[aux];
				aux += sliceOffset;//accessing im[pos + radius + kk * sliceOffset];
				imBlockIdx++;
			}

			//calculate median element
			imOut[pos] = medianBuffer(imBlock, blockSize);
			//next element
			pos++;
		}

		//do the final elements where block edge is in the border
		for (int ii = imDim[0] - radius; ii < imDim[0]; ii++)
		{
			//update imBlock with the new column
			for (int kk = -radius; kk <= radius; kk++)
			{
				imBlock[imBlockIdx % blockSize] = 0; //for now we assume zeros outsize image borders
				imBlockIdx++;
			}

			//calculate median element
			imOut[pos] = medianBuffer(imBlock, blockSize);
			//next element
			pos++;
		}
	}


	//last pass for the last radius rows
	for (int jj = imDim[1] - radius; jj < imDim[1]; jj++)
	{
		//reset imBlock: we load image data one column at a time (so the first radius column are zero		
		imBlockIdx = win * radius;
		memset(imBlock, 0, sizeof(imgType)*imBlockIdx);//for now we assume zeros outsize image borders
		for (int aa = 0; aa < radius; aa++)//the lsat colum will be done later
		{
			for (int bb = -radius; bb < imDim[1] - jj; bb++)
			{				
				imBlock[imBlockIdx++] = im[pos + aa + sliceOffset * bb];
			}

			for (int bb = imDim[1] - jj; bb <= radius; bb++)
			{
				imBlock[imBlockIdx++] = 0;//for now we assume zeros outsize image borders
			}
		}

		//advance the iwndow in the x-direction one column at a time
		for (int ii = 0; ii < imDim[0] - radius; ii++)
		{
			//update imBlock with the new column
			int aux = pos + radius - radius * sliceOffset;
			for (int kk = -radius; kk < imDim[1] - jj; kk++)
			{				
				imBlock[imBlockIdx % blockSize] = im[aux];
				aux += sliceOffset;//accessing im[pos + radius + kk * sliceOffset];
				imBlockIdx++;
			}			
			for (int kk = imDim[1] - jj; kk <= radius; kk++)
			{
				imBlock[imBlockIdx % blockSize] = 0;//for now we assume zeros outsize image borders
				imBlockIdx++;
			}

			//calculate median element
			imOut[pos] = medianBuffer(imBlock, blockSize);
			//next element
			pos++;
		}

		//do the final elements where block edge is in the border
		for (int ii = imDim[0] - radius; ii < imDim[0]; ii++)
		{
			//update imBlock with the new column
			for (int kk = -radius; kk <= radius; kk++)
			{
				imBlock[imBlockIdx % blockSize] = 0; //for now we assume zeros outsize image borders
				imBlockIdx++;
			}

			//calculate median element
			imOut[pos] = medianBuffer(imBlock, blockSize);
			//next element
			pos++;
		}
	}

	delete[] imBlock;
}

//declare all the possible types so template compiles properly
template void amatf_medianFilter2D<unsigned char>(const unsigned char* im, const int* imDim, int radius, unsigned char *imOut);
template void amatf_medianFilter2D<unsigned short int>(const unsigned short int* im, const int* imDim, int radius, unsigned short int *imOut);
template void amatf_medianFilter2D<float>(const float* im, const int* imDim, int radius, float *imOut);


//===============================================
template<class imgType>
void medianFilter2Dworker(imgType* im, const int* imDim, int radius)
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

		amatf_medianFilter2D(imBuffer, imDim, radius, slicePtr);

	}

	delete[] imBuffer;
}

template void medianFilter2Dworker<unsigned char>(unsigned char* im, const int* imDim, int radius);
template void medianFilter2Dworker<unsigned short int>(unsigned short int* im, const int* imDim, int radius);
template void medianFilter2Dworker<float>(float* im, const  int* imDim, int radius);
//===========================================================================

template<class imgType>
void medianFilter2DSliceBySlice(imgType* im,const  int* imDim, int radius)
{
	int numThreads = std::thread::hardware_concurrency();//all the cores available

	int sliceSize = imDim[0];
	for (int ii = 1; ii < dimsImageSlice; ii++)
		sliceSize *= imDim[ii];
	
	sliceId = 0;//initialize count

	// start the working threads
	std::vector<std::thread> threads;
	std::vector<int> errFlagVec(numThreads, 0);
	for (int i = 0; i < numThreads; ++i)
	{
		threads.push_back(std::thread(medianFilter2Dworker<imgType>, im, imDim, radius));
	}

	//wait for the workers to finish
	for (auto& t : threads)
		t.join();	
}

//declare all the possible types so template compiles properly
template void medianFilter2DSliceBySlice<unsigned char>(unsigned char* im, const  int* imDim, int radius);
template void medianFilter2DSliceBySlice<unsigned short int>(unsigned short int* im, const  int* imDim, int radius);
template void medianFilter2DSliceBySlice<float>(float* im, const int* imDim, int radius);