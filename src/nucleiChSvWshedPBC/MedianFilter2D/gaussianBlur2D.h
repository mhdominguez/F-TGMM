//Martin H. Dominguez
//2021 Gladstone Institutes

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
	\param[in]			radiusk		radius of the median filter. Median filter window size is 2*radius +1. So a 3x3 median filter has radius = 1
	\param[out]			imOut		pointer to filtered image in host. 
*/

template<class imgType>
void nvsepgold_gaussianBlur2D(const imgType* im, const int* imDim, int radius, imgType *imOut);


/*
	\brief apply gaussianBlur2D to a 3D stack slice by slice (or an RGB stack) using multithreading
	
*/

template<class imgType>
void gaussianBlur2DSliceBySlice(imgType* im,const int* imDim,int radius);


void gaussianBlurKernel(float sigma, int size, float* kernel);

#endif //__MEDIAN_FILTER_2D_H__
