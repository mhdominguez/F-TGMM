/*
 * knnCuda.h
 *
 *  Created on: Jul 15, 2011
 *      Author: amatf
 */

#ifndef KNNCUDA_V1_H_
#define KNNCUDA_V1_H_

#include "../constants.h"

void knn_v1_CUDAinPlace(int *indCUDA,float *queryCUDA,float *refTemp,long long int query_nb,int ref_nb);
void knn_v1_uploadScaleCUDA(float *scale);

#endif /* KNNCUDA_V1_H_ */
