//Martin H. Dominguez
//2021 Gladstone Institutes

#ifndef __GAUSSIAN_BLUR_2D_CUDA_H__
#define __GAUSSIAN_BLUR_2D_CUDA_H__


//define constants

#ifndef CUDA_CONSTANTS_FA
#define CUDA_CONSTANTS_FA
#ifndef CUDA_MAX_SIZE_CONST //to protect agains teh same constant define in other places in the code
#define CUDA_MAX_SIZE_CONST

static const int MAX_THREADS_CUDA = 1024; //adjust it for your GPU. This is correct for a 2.0 architecture
static const int MAX_BLOCKS_CUDA = 65535;
static const int BLOCK_SIDE = 32; //we load squares into share memory. 32*32 = 1024, which is the maximum number of threads per block for CUDA arch 2.0. 

#endif
#endif

#ifndef DIMS_IMAGE_SLICE
#define DIMS_IMAGE_SLICE
static const int dimsImageSlice = 2;//so thing can be set at co0mpile time

#endif



/*
	\brief Main function to calculaye median filter with CUDA
	
	\param[in/out]		im			pointer to image in host. It is overwritten by the median result
	\param[in]			imDim		array of length dimsImageSlice indicating the image dimensions. imDim[0] is the fastest running index in memory.
	\param[in]			radius		radius of the median filter. Median filter window size is 2*radius +1. So a 3x3 median filter has radius = 1
	\param[in]			devCUDA		in case you have multiple GPU in the same computer
*/

template<class imgType>
int gaussianBlurCUDA(imgType* im,int* imDim,float sigma,int devCUDA);


/*
	\brief apply gaussianBlur2D to a 3D stack slice by slice (or an RGB stack)
	
*/

template<class imgType>
int gaussianBlurCUDASliceBySlice(imgType* im,int* imDim,float sigma,int devCUDA);


int getKernelRadiusForSigmaCUDA(float sigma);
void gaussianBlurKernel(float sigma, int size, float* kernel);

#endif //__CONVOLUTION_3D_FFT_H__
