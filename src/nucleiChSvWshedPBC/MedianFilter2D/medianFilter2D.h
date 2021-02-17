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

#ifndef __MEDIAN_FILTER_2D_H__
#define __MEDIAN_FILTER_2D_H__



#ifndef DIMS_IMAGE_SLICE
#define DIMS_IMAGE_SLICE
static const int dimsImageSlice = 2;//so thing can be set at co0mpile time

#endif



/*
	\brief Main function to calculaye median filter with CPU (out of place)
	
	\param[in]			im			pointer to image in host. I
	\param[in]			imDim		array of length dimsImageSlice indicating the image dimensions. imDim[0] is the fastest running index in memory.
	\param[in]			radius		radius of the median filter. Median filter window size is 2*radius +1. So a 3x3 median filter has radius = 1
	\param[out]			imOut		pointer to filtered image in host. 
*/

template<class imgType>
void amatf_medianFilter2D(const imgType* im, const int* imDim, int radius, imgType *imOut);


/*
	\brief apply medianFilter2D to a 3D stack slice by slice (or an RGB stack) using multithreading
	
*/

template<class imgType>
void medianFilter2DSliceBySlice(imgType* im,const int* imDim,int radius);




#endif //__MEDIAN_FILTER_2D_H__