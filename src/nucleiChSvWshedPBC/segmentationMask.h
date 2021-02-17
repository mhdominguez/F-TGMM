/*
 * Copyright (C) 2011-2014 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat 
 *  segmentationMask.h
 *
 *  Created on: August 17th, 2014
 *      Author: Fernando Amat
 *
 * \brief Implements a segmentation mask contained with fast find (bottleneck in Watershed code). In this class the most costly operation is to merge two regions (we have to copy elements around)
 *
 */


#ifndef __AF_SEGMENTATION_MASK_H__
#define __AF_SEGMENTATION_MASK_H__

#include <vector>
#include <iostream>
#include <assert.h>
#include <limits>
#include "constants.h"

using namespace std;

template<class voxelType, class labelType>
class AF_segmentationMask
{
public:

	//constructor / desctructor
	AF_segmentationMask(int64 imgSize);
	~AF_segmentationMask(){};

	//easy get/set functions
	size_t getImgSize(){ return L.size();};
	labelType find(int64 pos) {return L[pos];};//returns label for image index pos. This was the bottleneck using set_union structure for watershed. Now it is constant time memory access
	voxelType getFmax(labelType label) {return fMaxVec[label];};
	bool sameComponent(int64 pos1, int64 pos2) {return L[pos1] == L[pos2];};
	labelType getNumLabels(){ return numLabels;};
	size_t getMaxLabels() {return PixelIdxList.size();};
	bool emptyLabel(labelType l){ return PixelIdxList[l].empty();};
	size_t sizeLabel(labelType l){ return PixelIdxList[l].size();};
	vector<uint64> getPixelIdxList(labelType l) {return PixelIdxList[l];};
	vector<uint64> getPixelIdxList(labelType l) const {return PixelIdxList[l];};

	//delete set
	void deleteSet(labelType l)
	{
		if( PixelIdxList[l].empty() == false )
		{			
			for(vector<uint64>::const_iterator iter = PixelIdxList[l].begin(); iter != PixelIdxList[l].end(); ++iter)
				L[*iter] = 0;
			PixelIdxList.clear();
			nextLabelAvail = std::min(nextLabelAvail, l);
			numLabels--;
		}
	};

	//merges two regions
	void unionSets(labelType l1, labelType l2);

	//adds elements to region
	void addElementToSet(int64 pos, labelType l,voxelType fVal);

	//start a new region
	void addNewComponent(int64 pos, voxelType fVal);

	int64 getRootIdx(labelType l)
	{
		if( PixelIdxList[l].empty() )
			return -1;
		else
			return PixelIdxList[l][0];
	};



	//----------------------debugging functions-------------------
	bool debugCheckNumLabels();

protected:

private:
	vector<labelType> L;//contains the label for each voxel in the image
	vector< vector<uint64> > PixelIdxList;//PixelIdxList[ii] contains the list of pixels of i-th label. Label 0 is reserved for background
	vector<voxelType> fMaxVec;//maximum value within each region
	size_t nextLabelAvail;//to keep track during insertions and deletions
	labelType numLabels;
};

//=======================================================================
template<class voxelType, class labelType>
AF_segmentationMask<voxelType, labelType>::AF_segmentationMask(int64 imgSize)
{
	L.resize(imgSize,0);
	nextLabelAvail = 1;//0 is reserved for background
	numLabels = 0;

	//reserve some initial space
	fMaxVec.reserve(1000);//we can expect at least 1000 regions
	PixelIdxList.reserve( 1000 );
}

//================================================
template<class voxelType, class labelType>
void AF_segmentationMask<voxelType, labelType>::unionSets(labelType l1, labelType l2)
{
	
	if( l1 == l2 )
		return;

	if( PixelIdxList[l1].empty() || PixelIdxList[l2].empty() )//nothing to merge
		return;

	//copy the smallest region into the largest to avoid copies
	if( PixelIdxList[l1].size() > PixelIdxList[l2].size() )
	{
		PixelIdxList[l1].insert( PixelIdxList[l1].end(), PixelIdxList[l2].begin(), PixelIdxList[l2].end()); 
		fMaxVec[l1] = fmax(fMaxVec[l1],fMaxVec[l2]);

		//change all the labels
		for(vector<uint64>::const_iterator iter = PixelIdxList[l2].begin(); iter != PixelIdxList[l2].end(); ++iter)
			L[*iter] = l1;

		PixelIdxList[l2].clear();
		nextLabelAvail = std::min((labelType)nextLabelAvail, l2);
	}else{
		PixelIdxList[l2].insert( PixelIdxList[l2].end(), PixelIdxList[l1].begin(), PixelIdxList[l1].end()); 
		fMaxVec[l2] = fmax(fMaxVec[l1],fMaxVec[l2]);

		//change all the labels
		for(vector<uint64>::const_iterator iter = PixelIdxList[l1].begin(); iter != PixelIdxList[l1].end(); ++iter)
			L[*iter] = l2;

		PixelIdxList[l1].clear();
		nextLabelAvail = std::min((labelType)nextLabelAvail, l1);
	}

	numLabels--;
}
//=========================================================
template<class voxelType, class labelType>
void AF_segmentationMask<voxelType, labelType>::addElementToSet(int64 pos, labelType l,voxelType fVal)
{
	assert( L[pos] == 0 );//current element should not have been assigned yet
	L[pos] = l;
	PixelIdxList[l].push_back(pos);
	fMaxVec[l] = fmax(fMaxVec[l],fVal);
}

//=========================================================
template<class voxelType, class labelType>
void AF_segmentationMask<voxelType, labelType>::addNewComponent(int64 pos, voxelType fVal)
{
	assert( L[pos] == 0 );//current element should not have been assigned yet
	L[pos] = nextLabelAvail;


	if( nextLabelAvail >= PixelIdxList.size() )
	{
		PixelIdxList.resize( nextLabelAvail + 1 );
		PixelIdxList[nextLabelAvail].reserve( 200 );
		fMaxVec.resize( nextLabelAvail + 1 );
	}

	assert( PixelIdxList[nextLabelAvail].empty() == true );	
	PixelIdxList[nextLabelAvail].push_back(pos);

	fMaxVec[nextLabelAvail] = fVal;

	//update next avilable label
	nextLabelAvail++;
	while(  nextLabelAvail < PixelIdxList.size() && PixelIdxList[nextLabelAvail].empty() == false )
		nextLabelAvail++;

	if( nextLabelAvail >= std::numeric_limits<labelType>::max() )
	{
		cout<<"ERROR: AF_segmentationMask<voxelType, labelType>::addNewComponent: no more labels available. Change the labelType to be able to segment the image with the required regions"<<endl;
		exit(2);
	}

	numLabels++;
}

//=================================================================
template<class voxelType, class labelType>
bool AF_segmentationMask<voxelType, labelType>::debugCheckNumLabels()
{
	size_t n = 0;
	for(size_t ii = 0; ii < PixelIdxList.size(); ii++)
	{
		if( PixelIdxList[ii].empty() == false )
			n++;
	}

	return n==numLabels;
}

#endif //__AF_SEGMENTATION_MASK_H__