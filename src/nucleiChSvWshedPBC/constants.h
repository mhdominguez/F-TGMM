#ifndef __CONSTANT_WATERSHED_H__
#define __CONSTANT_WATERSHED_H__

#ifndef DIMS_IMAGE_CONST //to protect agains teh same constant define in other places in the code
#define DIMS_IMAGE_CONST
const static int dimsImage = 3;//to be able to precompile code
#endif

#include <stdint.h>

typedef unsigned long long uint64;
typedef long long int64;

typedef unsigned short int imgVoxelType;//decide teh format of the image
typedef uint32_t imgLabelType;//decide teh format of the segmentation mask
static const  uint64 maxNumLabels = (1ULL<<32)-1;//maximum number of labels. Label = 0 is always reserved for background
static const imgVoxelType imgVoxelMinValue = 0;//min possible value for voxel values
static const imgVoxelType imgVoxelMaxValue = 65535;//max possible value for voxel values




#endif