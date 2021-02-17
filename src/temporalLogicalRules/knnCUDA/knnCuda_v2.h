/*
 * knnCuda.h
 *
 *  Created on: Jul 15, 2011
 *      Author: amatf
 */

#ifndef KNNCUDA_V2_H_
#define KNNCUDA_V2_H_

const int knn_v2_maxKNN=11; //maximum number of nearest neighbors considered for the spatio-temporal graph used in tracking

/*
\brief: simple kNN calculation. It handles all the GPU memory internally. If you need to perform kNN on the same data over and over this method is not recommeneded, but for single kNN searches it is easy to use.

\param[in] query : query points organized as x_1, x_2, x_3, ...,x_N,y_1, y_2,.....z_N  for GPU cache friendliness. Query are the points that we want to find the NN. So size(dist) = KmaxKNN * query_nb

*/

int knnCUDA_v2(int *ind,float* dist, float *query,float *ref,long long int query_nb,int ref_nb, int KNN, float* scale, int devCUDA);


#endif /* KNNCUDA_V2_H_ */
