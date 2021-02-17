/*
 * Copyright (C) 2011-2015 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat
 *  autoThreshold.cpp
 *
 *  Created on: January 7th, 2015
 *      Author: Fernando Amat
 *
 * \brief Simple functions to try to estimate global background threshold in stacks. Inspired by the Fiji plugin of the same name https://github.com/fiji/Auto_Threshold/blob/master/src/main/java/fiji/threshold/Auto_Threshold.java
 *
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include "autoThreshold.h"

#if defined(_WIN32) || defined(_WIN64)
	#define NOMINMAX
#else //Linux headers
#include <unistd.h>
#endif

using namespace std;


//triangle method
template<class imgType, int ndims>
imgType bckgAutoThr_TriangleMethod(imgType* im, const int* imDim)
{
	// Zack, G. W., Rogers, W. E. and Latt, S. A., 1977,
	// Automatic Measurement of Sister Chromatid Exchange Frequency,
	// Journal of Histochemistry and Cytochemistry 25 (7), pp. 741-753
	//
	// modified from Johannes Schindelin plugin
	//


	//basically calculate histogram
	vector<int64_t> data;
	imgType offset = autothreshold_preprocessing<imgType, ndims>(im, imDim, data);

	if (data.size() < 2)
		return 0;
	// find min and max
	int64_t min = 0, max = 0, min2 = 0, dmax = 0;
	for (int64_t i = 0; i < data.size(); i++) 
	{
		if (data[i]>0)
		{
			min = i;
			break;
		}
	}
	if (min > 0) 
		min--; // line to the (p==0) point, not to data[min]
	// The Triangle algorithm cannot tell whether the data is skewed to one side or another.
	// This causes a problem as there are 2 possible thresholds between the max and the 2 extremes
	// of the histogram.
	// Here I propose to find out to which side of the max point the data is furthest, and use that as
	// the other extreme.
	for (int64_t i = data.size() - 1; i > 0; i--) 
	{
		if (data[i] > 0)
		{
			min2 = i;
			break;
		}
	}
	if (min2 < data.size() - 1) 
		min2++; // line to the (p==0) point, not to data[min]

	for (int64_t i = 0; i < data.size(); i++) 
	{
		if (data[i] >dmax) 
		{
			max = i;
			dmax = data[i];
		}
	}
	// find which is the furthest side
	//IJ.log(""+min+" "+max+" "+min2);
	bool inverted = false;
	if ((max - min) < (min2 - max))
	{
		// reverse the histogram
		//IJ.log("Reversing histogram.");
		inverted = true;
		int64_t left = 0; // index of leftmost element
		int64_t right = data.size() - 1; // index of rightmost element
		while (left < right) 
		{
			// exchange the left and right elements
			int64_t temp = data[left];
			data[left] = data[right];
			data[right] = temp;
			// move the bounds toward the center
			left++;
			right--;
		}
		min = data.size() - 1 - min2;
		max = data.size() - 1 - max;
	}
	if (min == max)
	{
		//IJ.log("Triangle: min == max.");
		return min + offset;
	}
	// describe line by nx * x + ny * y - d = 0
	double nx, ny, d;
	// nx is just the max frequency as the other point has freq=0
	nx = data[max]; //-min; // data[min]; // lowest value bmin = (p=0)% in the image
	ny = min - max;
	d = sqrt(nx * nx + ny * ny);
	nx /= d;
	ny /= d;
	d = nx * min + ny * data[min];
	// find split point
	int64_t split = min;
	double splitDistance = 0;
	for (int64_t i = min + 1; i <= max; i++) 
	{
		double newDistance = nx * i + ny * data[i] - d;
		if (newDistance > splitDistance) 
		{
			split = i;
			splitDistance = newDistance;
		}
	}
	split--;
	if (inverted) {
		// The histogram might be used for something else, so let's reverse it back
		int64_t left = 0;
		int64_t right = data.size() - 1;
		while (left < right) 
		{
			int64_t temp = data[left];
			data[left] = data[right];
			data[right] = temp;
			left++;
			right--;
		}
		return (data.size() - 1 - split) + offset;
	}
	else
		return split + offset;


}


//================================================
template<class imgType, int ndims>
imgType* maximumIntensityProjection(const imgType* im, const int* imDim)
{
	int64_t offsetSlice = imDim[0];
	for (int ii = 1; ii < ndims - 1; ii++)
		offsetSlice *= (int64_t)(imDim[ii]); 

	imgType *imMIP = new imgType[offsetSlice];
	memset(imMIP, 0, sizeof(imgType)* offsetSlice);

	int64_t p = 0;
	for (int ii = 0; ii < imDim[ndims - 1]; ii++)
	{
		for (int64_t jj = 0; jj < offsetSlice; jj++)
		{
			imMIP[jj] = max(imMIP[jj], im[p]);
			p++;
		}
	}

	return imMIP;
}


//================================================
template<class imgType, int ndims>
imgType autothreshold_preprocessing(const imgType* im, const int* imDim, std::vector<int64_t>& hist)
{


	int64_t imSize = imDim[0];
	for (int ii = 1; ii < ndims; ii++)
		imSize *= (int64_t)(imDim[ii]);


	//calculate histogram
	vector<int> data(pow(2,sizeof(imgType) * 8 ), 0);	
	for (int64_t ii = 0; ii < imSize; ii++)
		data[im[ii]]++;

	//ignore black or white in the histogram
	//if (noBlack) data[0] = 0;
	//if (noWhite) data[data.size() - 1] = 0;
	data[0] = 0;//to avoid biasing in case images have been background corrected

	// bracket the histogram to the range that holds data to make it quicker
	imgType minbin = 0, maxbin = 0;
	for (int i = 0; i<data.size(); i++){
		if (data[i]>0) maxbin = i;
	}
	for (int i = data.size() - 1; i >= 0; i--){
		if (data[i]>0) minbin = i;
	}
	//IJ.log (""+minbin+" "+maxbin);
	hist.resize((maxbin - minbin) + 1);
	for (int i = minbin; i <= maxbin; i++)
	{
		hist[i - minbin] = data[i];;
	}

	return minbin;
}

//=========================================================
//template declaration
template uint8_t bckgAutoThr_TriangleMethod<uint8_t, 2>( uint8_t* im, const int* imDim);
template uint8_t bckgAutoThr_TriangleMethod<uint8_t, 3>( uint8_t* im, const int* imDim);
template uint16_t bckgAutoThr_TriangleMethod<uint16_t, 2>( uint16_t* im, const int* imDim);
template uint16_t bckgAutoThr_TriangleMethod<uint16_t, 3>( uint16_t* im, const int* imDim);



template uint8_t* maximumIntensityProjection<uint8_t, 2>(const uint8_t* im, const int* imDim);
template uint8_t* maximumIntensityProjection<uint8_t, 3>(const uint8_t* im, const int* imDim);
template uint16_t* maximumIntensityProjection<uint16_t, 2>(const uint16_t* im, const int* imDim);
template uint16_t* maximumIntensityProjection<uint16_t, 3>(const uint16_t* im, const int* imDim);

template uint8_t autothreshold_preprocessing<uint8_t, 2>(const uint8_t* im, const int* imDim, vector<int64_t>& hist);
template uint8_t autothreshold_preprocessing<uint8_t, 3>(const uint8_t* im, const int* imDim, vector<int64_t>& hist);
template uint16_t autothreshold_preprocessing<uint16_t, 2>(const uint16_t* im, const int* imDim, vector<int64_t>& hist);
template uint16_t autothreshold_preprocessing<uint16_t, 3>(const uint16_t* im, const int* imDim, vector<int64_t>& hist);


