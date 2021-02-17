/*
* Copyright (C) 2011-2015 by  Fernando Amat
* See license.txt for full license and copyright notice.
*
* Authors: Fernando Amat
*  autoThreshold.h
*
*  Created on: January 7th, 2015
*      Author: Fernando Amat
*
* \brief Simple functions to try to estimate global background threshold in stacks. Inspired by the Fiji plugin of the same name https://github.com/fiji/Auto_Threshold/blob/master/src/main/java/fiji/threshold/Auto_Threshold.java
*
*/


#ifndef __FA_AUTO_THRESHOLD_H__
#define __FA_AUTO_THRESHOLD_H__

#include <vector>
#include "stdint.h"

//triangle method: translated from Java from https://github.com/fiji/Auto_Threshold/blob/master/src/main/java/fiji/threshold/Auto_Threshold.java
template<class imgType, int ndims>
imgType bckgAutoThr_TriangleMethod(imgType* im, const int* imDim);


//============================================================================
//auxiliary functions
template<class imgType, int ndims>
imgType* maximumIntensityProjection(const imgType* im, const int* imDim);


template<class imgType, int ndims>
imgType autothreshold_preprocessing(const imgType* im, const int* imDim, std::vector<int64_t>& hist);

#endif
