/*
 * Copyright (C) 2011-2013 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat 
 *  supportFunctionsMain.h
 *
 *  Created on: January 21st, 2012
 *      Author: Fernando Amat
 *
 * \brief Contains different routines that are called form the main.cpp part of theprogram to not make the code so cluttr
 *        
 *
 */
#ifndef __SUPPORT_FUNCTIONS_MAIN_H__
#define __SUPPORT_FUNCTIONS_MAIN_H__

#include <vector>
#include <string>

#include "lineageHyperTree.h"
//#include "hierarchicalSegmentation.h"
#include "GaussianMixtureModel.h"



/*
\brief Parse lineage hypertree model to GaussianMixtureModel class
*/
template <class imgTypeC>
int parseNucleiList2TGMM(std::vector<GaussianMixtureModel*> &vecGM, lineageHyperTree &lht,int frame, bool regularizeW, float thrDist2LargeDisplacement,imgTypeC* imgPtr);




/*
\brief To read file patterns using ???? for frame numbers with different precisions
*/
void parseImagePath(string& imgRawPath, int frame);


void transposeStackUINT16(mylib::Array *img);

//---------------------implementing LP to fix issues---------------------------


/*
\brief same as int lineageHyperTree::extendDeadNucleiAtTM(int TM, int& numExtensions, int &numDeaths) but using Hierarchical Segmentaion information to entend death track with Hungarian algorithm (so supervoxels can be split)
*/
int extendDeadNucleiAtTMwithHS(lineageHyperTree &lht, hierarchicalSegmentation* hsForward, int TM, int& numExtensions, int &numDeaths,float* imgPtr);
int extendDeadNucleiWithHS(lineageHyperTree &lht, hierarchicalSegmentation* hsForward, TreeNode<ChildrenTypeLineage>* rootDead,float *imgPtr);

/*
\brief 

\param	return			returns root of the lineage at iniTM
*/
TreeNode< ChildrenTypeLineage >* addSupervoxelsPointersFromLineage(std::vector< vector<supervoxel*> > &svIniVec, int iniTM, int endTM, TreeNode< ChildrenTypeLineage >* node);



int mkDirectory(const std::string& folder);

//=============================================================================
//--------------------------debugging------------------------------

/*
\brief Function to visualize hierarchical segmentation as a lineage. Preallocate lht with enough time points

\warning This function does not trim the supervoxels
*/
int parseHierarchicalSegmentation2LineageHyperTree(hierarchicalSegmentation* hs,lineageHyperTree& lht);


/*

\brief function to test split / merge hypothesis starting from a segmentation with a fix tau
\warning This function does trim supervoxels
*/
int debugMergeSplitHStoLHT(hierarchicalSegmentation* hs,lineageHyperTree& lht, imgVoxelType tau);

#endif //__SUPPORT_FUNCTIONS_MAIN_H__
