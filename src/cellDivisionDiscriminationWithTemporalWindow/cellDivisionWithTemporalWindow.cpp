/*
 * Copyright (C) 2011-2013 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat 
 *  cellDivisionWithTemporalWindow.cpp
 *
 *  Created on: April 21st, 2014
 *      Author: Fernando Amat
 *
 * \brief temporal features (image-based and geometrical and spatial) to discriminate between true and false cell divisions detected by the TGMM pipeline
 *        
 *
 */

#include <assert.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include "cellDivisionWithTemporalWindow.h"
#include "utilsAmatf.h"

using namespace std;

static void breakme() {
    cout << "break me:" << __FILE__ << "(" << __LINE__<<")" << endl;
}

#define CHECK(e) do{if(!(e)) {breakme();throw runtime_error("Expression evaluated as false. " #e);}} while(0)

template<typename T> static inline mylib::Value_Type MylibValueType();
template<> inline mylib::Value_Type MylibValueType<float>(){ return mylib::FLOAT32_TYPE; }
template<> inline mylib::Value_Type MylibValueType<unsigned short>(){ return mylib::UINT16_TYPE; }
template<> inline mylib::Value_Type MylibValueType<unsigned char>(){ return mylib::UINT8_TYPE; }

int cellDivisionWithTemporalWindow::temporalWindowRadius;
bool cellDivisionWithTemporalWindow::useFixPrecisionMatrix;//true->we use the same precision W for all time points within the temporal window. So Haar features capture changes in size (but we lose excentricity calculations)
double cellDivisionWithTemporalWindow::precisionWfix[dimsImage * (1+dimsImage) /2];//if useFixPrecisionMatrix == true, then we use the value in here to calculate 3D Elliptical Haar features
int cellDivisionWithTemporalWindow::temporalUndersampling = 1;//default value

short int cellDivisionWithTemporalWindow::HaarPatternTriplet[39] = {0,0,1,0,1,-1,0,1,0,0,1,1,1,-1,-1,1,-1,0,1,-1,1,1,0,-1,1,0,0,1,0,1,1,1,-1,1,1,0,1,1,1}; //39 = 13x3 elements. Each triplet is a unique (linearly independent) Haar pattern. generated with EllipticalHaarFeatures::generateHaarPattern
int cellDivisionWithTemporalWindow::HaarPatternTripletN = 13;



double logbase(double a, double base)
{
   return log(a) / log(base);
}

cellDivisionWithTemporalWindow::cellDivisionWithTemporalWindow()
{
	//temporalWindowRadius = 0;//you have to setup temporal windiw size yourself since it is a static variable
	cellDivisionPtr = NULL;
};

cellDivisionWithTemporalWindow::cellDivisionWithTemporalWindow(int temporalWindowRadius_)
{
	temporalWindowRadius = temporalWindowRadius_;
	fHaarVec.resize(temporalWindowRadius + 1 + 2*temporalWindowRadius);
	cellDivisionPtr = NULL;
};

cellDivisionWithTemporalWindow::~cellDivisionWithTemporalWindow()
{
	
};

//=================================================================
int cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchForCellDivision(const vector< TreeNode< ChildrenTypeLineage >* >& divisionNodes, const vector<mylib::Array*>& imgVec, vector<cellDivisionWithTemporalWindow>& cdtwVec, int devCUDA, int symmetry)
{
	assert( imgVec.size() == 2*temporalWindowRadius + 1 );

	if( divisionNodes.empty() == true )
		return 0;
	

	int numEllipsoids = divisionNodes.size();
	int sizeW = dimsImage * (1+dimsImage) / 2;
	int err = 0;

	//allocate memory for centroid and precision matrix
	double *m = new double[dimsImage * numEllipsoids];
	double *W = new double[sizeW  * numEllipsoids];
	long long int *dimsVec = new long long int[imgVec[0]->ndims];
	
	//auxliary variables
	TreeNode< ChildrenTypeLineage >* auxNode = NULL;
	vector< TreeNode< ChildrenTypeLineage >* > auxNodeVec(numEllipsoids, NULL);//to store pointer for each lineage


	//setup return elements
	cdtwVec.resize( numEllipsoids );
	for(int ii = 0; ii < numEllipsoids; ii++)
	{
		assert(divisionNodes[ii]->getNumChildren() == 2);
		cdtwVec[ii].fHaarVec.resize(temporalWindowRadius + 1 + 2*temporalWindowRadius);
		auxNodeVec[ii] = divisionNodes[ii];
		cdtwVec[ii].cellDivisionPtr = divisionNodes[ii];//initialize this pointer
	}

	//calculate features backwards (including cell division point)
	for( int tt = temporalWindowRadius; tt >= 0; tt--)
	{

		err = calculateBasicEllipticalHaarFeaturesBatchAtTM( auxNodeVec, imgVec[tt], cdtwVec, tt, devCUDA, symmetry, &m , &W, &dimsVec);
		if( err > 0 )
			return err;
		//move backwards in the lineage
		for(int ii = 0; ii < numEllipsoids; ii++)
		{
			if( auxNodeVec[ii] != NULL )
			{				
				auxNodeVec[ii] = auxNodeVec[ii]->parent;
			}
		}				
	}

	//reset elements to left or right daughter	
	int ttOffset = temporalWindowRadius + 1;//to save in the correct place in fHaarVec
	for(int ch = 0; ch < 2; ch++)
	{
		if( ch == 0 )
		{
			for(int ii = 0; ii < numEllipsoids; ii++)	
				auxNodeVec[ii] = divisionNodes[ii]->left;
		}else{
			for(int ii = 0; ii < numEllipsoids; ii++)	
				auxNodeVec[ii] = divisionNodes[ii]->right;
		}
		//iterate over time points
		for( int tt = temporalWindowRadius + 1; tt < 2* temporalWindowRadius + 1; tt++)
		{
			err = calculateBasicEllipticalHaarFeaturesBatchAtTM( auxNodeVec, imgVec[tt], cdtwVec, ttOffset, devCUDA, symmetry, &m , &W, &dimsVec);
			if( err > 0 )
				return err;
			//move forwards in the lineage
			for(int ii = 0; ii < numEllipsoids; ii++)
			{
				auxNodeVec[ii] = moveForwardInlineageEvenWithCellDivision(auxNodeVec[ii]);				
			}				
			ttOffset++;
		}
	}
	
	//release memory
	delete[] m;
	delete[] W;
	delete[] dimsVec;

	return 0;
}

//=================================================================
int cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchForCellDivisionSingleWindowForDaughters(
	const vector< TreeNode< ChildrenTypeLineage >* >& divisionNodes, 
	const vector<mylib::Array*>& imgVec, 
	vector<cellDivisionWithTemporalWindow>& cdtwVec, 
	int devCUDA, 
	int symmetry)
{
	assert( imgVec.size() == 2*temporalWindowRadius + 1 );

	if( divisionNodes.empty() == true )
		return 0;
	
	int numEllipsoids = divisionNodes.size();
	const int sizeW = dimsImage * (1+dimsImage) / 2;
	int err = 0;

	//allocate memory for centroid and precision matrix
	double *m = new double[dimsImage * numEllipsoids];
	double *W = new double[sizeW  * numEllipsoids];
	long long int *dimsVec = new long long int[dimsImage];
	
	//auxliary variables
	TreeNode< ChildrenTypeLineage >* auxNode = NULL;
	vector< TreeNode< ChildrenTypeLineage >* > auxNodeVec(numEllipsoids, NULL);//to store pointer for each lineage
	vector< TreeNode< ChildrenTypeLineage >* > auxNodeVecDR(numEllipsoids, NULL), auxNodeVecDL(numEllipsoids, NULL);//for left and right daughter

	//setup return elements
	cdtwVec.resize( numEllipsoids );
	for(int ii = 0; ii < numEllipsoids; ii++)
	{
		assert(divisionNodes[ii]->getNumChildren() == 2);
		cdtwVec[ii].fHaarVec.resize(1 + 2*temporalWindowRadius);
		auxNodeVec[ii] = divisionNodes[ii];
		cdtwVec[ii].cellDivisionPtr = divisionNodes[ii];//initialize this pointer
	}

	//calculate features backwards (including cell division point)
	for( int tt = temporalWindowRadius; tt >= 0; tt--)	
	{
		err = calculateBasicEllipticalHaarFeaturesBatchAtTM( auxNodeVec, imgVec[tt], cdtwVec, tt, devCUDA, symmetry, &m , &W, &dimsVec);
		if( err > 0 )
			return err;
		//move backwards in the lineage
		for(int ii = 0; ii < numEllipsoids; ii++)
		{
			if( auxNodeVec[ii] != NULL )
			{				
				auxNodeVec[ii] = auxNodeVec[ii]->parent;
			}
		}				
	}


	
	//calculate features forward (we center the box in the center of the two daughters)
	float d1, d2;

	for(int ii = 0; ii < numEllipsoids; ii++)	
	{
		auxNodeVecDL[ii] = divisionNodes[ii]->left;
		auxNodeVecDR[ii] = divisionNodes[ii]->right;
	}
	//iterate over time points
	vector<bool> isDead(numEllipsoids, false);
	for( int tt = temporalWindowRadius + 1; tt < 2* temporalWindowRadius + 1; tt++)
	{

		
		//update centroids, precision matrix and isDead before calling Haar elliptical features
		if( useFixPrecisionMatrix == true )
		{
			for(int ii = 0; ii < numEllipsoids; ii++)
			{
				if( auxNodeVecDL[ii] == NULL || auxNodeVecDR[ii] == NULL )
				{
					isDead[ii] = true;
					//we just preserve old m and W although we will disregard results
				}else{//we have both daughters
					isDead[ii] = false;
					for(int aa = 0; aa<dimsImage; aa++)
					{
						m[ii + aa * numEllipsoids] = 0.5 * (auxNodeVecDL[ii]->data->centroid[aa] + auxNodeVecDR[ii]->data->centroid[aa]);
					}
					for(int aa = 0; aa<sizeW; aa++)
					{
						W[ii + aa * numEllipsoids] = precisionWfix[aa];
					}
				}
			}
		}else{
			double auxW[sizeW], auxM[dimsImage];

			for(int ii = 0; ii < numEllipsoids; ii++)
			{
				if( auxNodeVecDL[ii] == NULL || auxNodeVecDR[ii] == NULL )
				{
					isDead[ii] = true;
					//we just preserve old m and W although we will disregard results
				}else{//we have both daughters
					isDead[ii] = false;
					//assert( supervoxel::dataSizeInBytes / (supervoxel::dataDims[0] * supervoxel::dataDims[1] * supervoxel::dataDims[2]) == 2 );//uint16 pointers					
					//cellDivisionWithTemporalWindow::calculatePrecisionMatrixForJointDuaghters<mylib::uint16>(*(auxNodeVecDL[ii]->data), *(auxNodeVecDR[ii]->data), auxM, auxW);
					cellDivisionWithTemporalWindow::calculatePrecisionMatrixForJointDuaghters(*(auxNodeVecDL[ii]->data), *(auxNodeVecDR[ii]->data), auxW);		//not using the image information			
					for(int aa = 0; aa<dimsImage; aa++)
					{
						//m[ii + aa * numEllipsoids] = auxM[aa];
						m[ii + aa * numEllipsoids] = 0.5 * (auxNodeVecDL[ii]->data->centroid[aa] + auxNodeVecDR[ii]->data->centroid[aa]);
					}
					for(int aa = 0; aa<sizeW; aa++)
					{						
						W[ii + aa * numEllipsoids] = auxW[aa];
					}
				}
			}
		}

		//calculate features
		err = calculateBasicEllipticalHaarFeaturesBatchAtTM( isDead, imgVec[tt], cdtwVec, tt, devCUDA, symmetry, m , W, &dimsVec);
		if( err > 0 )
			return err;
		//move forwards in the lineage
		for(int ii = 0; ii < numEllipsoids; ii++)
		{
			//left daughter
			auxNodeVecDL[ii] = moveForwardInlineageEvenWithCellDivision(auxNodeVecDL[ii]);
			//right daughter
			auxNodeVecDR[ii] = moveForwardInlineageEvenWithCellDivision(auxNodeVecDR[ii]);			
		}				
	}

	//release memory
	delete[] m;
	delete[] W;
	delete[] dimsVec;

	return 0;
}

//============================================================================================
TreeNode<ChildrenTypeLineage>* cellDivisionWithTemporalWindow::moveForwardInlineageEvenWithCellDivision(TreeNode<ChildrenTypeLineage>* auxNode)
{
	if( auxNode != NULL )			
	{				
		switch( auxNode->getNumChildren() )
		{
		case 0://death
			auxNode = NULL;
			break;
		case 1://normal continuation
			if( auxNode->left != NULL )
				auxNode = auxNode->left;
			else
				auxNode = auxNode->right;
			break;
		case 2://another cell division within temporal (tricky). Right now we apply a simple rule to continue with the closest daughter
			{
				float d1 = auxNode->data->Euclidean2Distance(*(auxNode->left->data), supervoxel::getScale());
				float d2 = auxNode->data->Euclidean2Distance(*(auxNode->right->data), supervoxel::getScale());
				if( d1 < d2 )//chose left
					auxNode = auxNode->left;
				else
					auxNode = auxNode->right;
			}
			break;
		}				
	}
	return auxNode;
}

//========================================================
int cellDivisionWithTemporalWindow::calculateFeatures()
{
	//NOTE: you should reserve space for f before calling this function
	assert(fHaarVec.size() ==  3* temporalWindowRadius + 1);

	//features are added automatimatically to this->f
	f.clear();

	//-------------------------calculate image-based features------------------------
	int numFeaturesHaar = calculateHaarBasedFeatures();

	//------------------------calculate trajectory-based (geometry, shape, etc) features-------------------
	int numFeaturesLineage = calculateLineageBasedFeatures();		




	//check that all features are real numbers
	int count = 0;
	for(vector<float>::iterator iter = f.begin(); iter != f.end(); ++iter, count++)
	{
		#if defined(_WIN32) || defined(_WIN64)
		if( _finite( *iter ) == false )//degenerated cases
		{
			if( _isnan( *iter ) == true )//this should not happen
			{
				cout<<"WARNING: at cellDivisionWithTemporalWindow::calculateHaarBasedFeatures: feature at position "<<count<<" is NAN "<<*iter<<endl;
				*iter = 0.0f;//default value
			}else{ //element is inf->just set a large value
				if( *iter < 0 )
					*iter = -1e32;
				else
					*iter = -1e32;
			}
		}
#else
		if( isfinite( *iter ) )//degenerated cases
		{
			if( isnan( *iter ) == true )//this should not happen
			{
				cout<<"WARNING: at cellDivisionWithTemporalWindow::calculateHaarBasedFeatures: feature at position "<<count<<" is NAN "<<*iter<<endl;
				*iter = 0.0f;//default value
			}else{ //element is inf->just set a large value
				if( *iter < 0 )
					*iter = -1e32;
				else
					*iter = -1e32;
			}
		}
#endif
		
	}

	return 0;

}

//========================================================
int cellDivisionWithTemporalWindow::calculateFeaturesSingleWindowForDaughters()
{
	//NOTE: you should reserve space for f before calling this function
	assert(fHaarVec.size() ==  2 * temporalWindowRadius + 1);

	//features are added automatimatically to this->f
	f.clear();
#ifdef CDWT_SAVE_FEATURE_NAME
	fName.clear();
#endif

	//-------------------------calculate image-based features------------------------
	int numFeaturesHaar = calculateHaarBasedFeaturesSingleWindowForDaughters();

	//------------------------calculate trajectory-based (geometry, shape, etc) features-------------------
	int numFeaturesLineage = calculateLineageBasedFeatures();//we can share these features		




	//check that all features are real numbers
	int count = 0;
	for(vector<float>::iterator iter = f.begin(); iter != f.end(); ++iter, count++)
	{
		#if defined(_WIN32) || defined(_WIN64)
		if( _finite( *iter ) == false )//degenerated cases
		{
			if( _isnan( *iter ) == true )//this should not happen
			{
				cout<<"WARNING: at cellDivisionWithTemporalWindow::calculateHaarBasedFeatures: feature at position "<<count<<" is NAN "<<*iter<<endl;
				*iter = 0.0f;//default value
			}else{ //element is inf->just set a large value
				if( *iter < 0 )
					*iter = -1e32;
				else
					*iter = 1e32;
			}
		}
#else
		if( isfinite( *iter ) )//degenerated cases
		{
			if( isnan( *iter ) == true )//this should not happen
			{
				cout<<"WARNING: at cellDivisionWithTemporalWindow::calculateHaarBasedFeatures: feature at position "<<count<<" is NAN "<<*iter<<endl;
				*iter = 0.0f;//default value
			}else{ //element is inf->just set a large value
				if( *iter < 0 )
					*iter = -1e32;
				else
					*iter = 1e32;
			}
		}
#endif
		
	}

	return 0;

}

//===============================================================================================
//this function is meant for speed, not for inteligibility

int cellDivisionWithTemporalWindow::calculateFeaturesSingleWindowForDaughters_featureSelection_v1()
{
	
	assert(fHaarVec.size() ==  2 * temporalWindowRadius + 1);

	//features are added automatimatically to this->f
	f.clear();			
#ifdef CDWT_SAVE_FEATURE_NAME	
	fName.clear();
#endif
	vector<float> timeSeries( getTotalTemporalWindowSizeSingleWindowForDaughters() );



	//temporal features for ring intensity
	for( int ii = 0; ii <= 1; ii++ )//only first 2 rings seem important
	{		
		//generate time series
		for( size_t jj = 0; jj < timeSeries.size(); jj++)
			timeSeries[jj] = fHaarVec[jj].ringAvgIntensity[ii];

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("ellipHaar3D_ring"),ii);
#endif
		//calculate features
		temporalFeatureExtractionSingleWindowForDaughters(timeSeries);				


		if( basicEllipticalHaarFeatureVector::useDoGfeatures == true ) // do the same for DoG
		{
			//generate time series
			for( size_t jj = 0; jj < timeSeries.size(); jj++)
				timeSeries[jj] = fHaarVec[jj].ringAvgIntensityDoG[ii];

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("ellipHaar3D_ring_DoG"),ii);
#endif
			//calculate features
			temporalFeatureExtractionSingleWindowForDaughters(timeSeries);

		}
	}


	//temporal features for excentricity values
	//for( int ii = 0; ii < dimsImage*(dimsImage-1)/2; ii++ )
	for( int ii = 1; ii <= 1; ii++ )//only excentricity 1 seems relevant
	{
		//generate time series
		for( size_t jj = 0; jj < timeSeries.size(); jj++)
			timeSeries[jj] = fHaarVec[jj].excentricity[ii];
#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("excentricity"),ii);
#endif
		//calculate features
		temporalFeatureExtractionSingleWindowForDaughters(timeSeries);	
	}
	
	//---------------------------------------------------------------------------
	//---------------------------lineage-based features -----------------------------------
	
	int nTS = getTotalTemporalWindowSize();
	float auxF;
	timeSeries.resize( nTS ); //preallocation
	TreeNode<ChildrenTypeLineage>* auxNode;

	//extract a vector of pointers to each point in the lineage around teh cell division (NULL indicates death or birth)
	vector< TreeNode<ChildrenTypeLineage>* > lineageVec;
	lineagePointersWithinTemporalWindowWithCellDivision( lineageVec );

	assert( lineageVec.size() == timeSeries.size());

	
	
	//----------------distance of all cells to the midplane defined by two daughters------------------------------
	float *centroidDL = cellDivisionPtr->left->data->centroid;
	float *centroidDR = cellDivisionPtr->right->data->centroid;

	for(int ii = 0; ii < nTS; ii++)
	{
		if( lineageVec[ii] == NULL )
			timeSeries[ii] = 100.0f;//default large number
		else
			timeSeries[ii] = lineageHyperTree::cellDivisionPlaneDistance(lineageVec[ii]->data->centroid, centroidDL, centroidDR);
	}

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("dist2midplane"),0);
#endif
	//calculate features
	temporalFeatureExtraction(timeSeries);


	//----------------distance between daughters------------------------------	
	int tt;
	for(tt = 0; tt <= temporalWindowRadius; tt++)
	{
		timeSeries[tt] = 0;//at the beginning there are no daughters
	}
	for(; tt <= 2 * temporalWindowRadius; tt++)
	{
		if( lineageVec[tt] == NULL || lineageVec[tt+ temporalWindowRadius] == NULL )
		{
			timeSeries[tt] = 0;
			timeSeries[tt+ temporalWindowRadius] = 0;
		}else{
			auxF = lineageVec[tt]->data->Euclidean2Distance( *(lineageVec[tt + temporalWindowRadius]->data), supervoxel::getScale() );
			timeSeries[tt] = auxF;
			timeSeries[tt+ temporalWindowRadius] = auxF;//I cna use this space for some other metric involving the two daughters
		}
	}
	#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("distBetweenDaughters"),0);
#endif
	//calculate features
	temporalFeatureExtraction(timeSeries, false, false);//no multiscale needed


	//----------------difference in volume (pixels) with respect to parent------------------------------	
	for(int ii = 0; ii < nTS; ii++)
	{
		auxNode = lineageVec[ii];
		if( auxNode != NULL && auxNode->parent != NULL )
		{
			auxF = lineageHyperTree::getNucleusVolume(*(auxNode->data));
			timeSeries[ii] = ( auxF - lineageHyperTree::getNucleusVolume(*(auxNode->parent->data)) ) / auxF;
		}else{
			timeSeries[ii] = 0.0f;
		}
	}

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("delta_volume"),0);
#endif
	//calculate features
	temporalFeatureExtraction(timeSeries);

	//----------------offset Jaccard distance with respect to parent------------------------------	
	for(int ii = 0; ii < nTS; ii++)
	{
		auxNode = lineageVec[ii];
		if( auxNode != NULL && auxNode->parent != NULL )
		{
			timeSeries[ii] = lineageHyperTree::JaccardDistance(auxNode, true);
		}else{
			timeSeries[ii] = 2.0f;//Jaccard distanc eis bounded by 1.0
		}
	}
#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("delta_Jaccard_wOffset"),0);
#endif
	//calculate features
	temporalFeatureExtraction(timeSeries);

	//----------------Jaccard distance with respect to parent------------------------------	
	for(int ii = 0; ii < nTS; ii++)
	{
		auxNode = lineageVec[ii];
		if( auxNode != NULL && auxNode->parent != NULL )
		{
			timeSeries[ii] = lineageHyperTree::JaccardDistance(auxNode, false);
		}else{
			timeSeries[ii] = 2.0f;//Jaccard distanc eis bounded by 1.0
		}
	}

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("delta_Jaccard"),0);
#endif
	//calculate features
	temporalFeatureExtraction(timeSeries, false, false);//no multi-scale needed

	//----------------absolute displacement with respect to parent------------------------------	
	for(int ii = 0; ii < nTS; ii++)
	{
		auxNode = lineageVec[ii];
		if( auxNode != NULL && auxNode->parent != NULL )
		{
			timeSeries[ii] = auxNode->data->Euclidean2Distance( *(auxNode->parent->data), supervoxel::getScale() );
		}else{
			timeSeries[ii] = 100.0f;//deafult last number
		}
	}
	#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("delta_displacement"),0);
#endif
	//calculate features
	temporalFeatureExtraction(timeSeries, false, false);		
	
	return 0;
		
}



//================================================================
int cellDivisionWithTemporalWindow::calculateHaarBasedFeatures()
{
	int nFeatures = 0;
	vector<float> timeSeries( getTotalTemporalWindowSize() );


	//temporal features for ring intensity
	for( int ii = 0; ii < basicEllipticalHaarFeatureVector::numRings; ii++ )
	{		
		//generate time series
		for( size_t jj = 0; jj < timeSeries.size(); jj++)
			timeSeries[jj] = fHaarVec[jj].ringAvgIntensity[ii];
		//calculate features
		nFeatures += temporalFeatureExtraction(timeSeries);				


		if( basicEllipticalHaarFeatureVector::useDoGfeatures == true ) // do the same for DoG
		{
			//generate time series
			for( size_t jj = 0; jj < timeSeries.size(); jj++)
				timeSeries[jj] = fHaarVec[jj].ringAvgIntensityDoG[ii];
			//calculate features
			nFeatures += temporalFeatureExtraction(timeSeries);

		}
	}
	//temporal features for radial cell intensity
	for( int ii = 0; ii < basicEllipticalHaarFeatureVector::numCells; ii++ )
	{
		//generate time series
		for( size_t jj = 0; jj < timeSeries.size(); jj++)
			timeSeries[jj] = fHaarVec[jj].cellAvgIntensity[ii];
		//calculate features
		nFeatures += temporalFeatureExtraction(timeSeries);


		if( basicEllipticalHaarFeatureVector::useDoGfeatures == true ) // do the same for DoG
		{
			//generate time series
			for( size_t jj = 0; jj < timeSeries.size(); jj++)
				timeSeries[jj] = fHaarVec[jj].cellAvgIntensityDoG[ii];
			//calculate features
			nFeatures += temporalFeatureExtraction(timeSeries);

		}
	}
	//temporal features for excentricity values
	for( int ii = 0; ii < dimsImage*(dimsImage-1)/2; ii++ )
	{
		//generate time series
		for( size_t jj = 0; jj < timeSeries.size(); jj++)
			timeSeries[jj] = fHaarVec[jj].excentricity[ii];
		//calculate features
		nFeatures += temporalFeatureExtraction(timeSeries);	
	}
	

	return nFeatures;
}

//=================================
#ifdef CDWT_SAVE_FEATURE_NAME	
void cellDivisionWithTemporalWindow::setFeaturePrefix(string& s, int p)
{
	char buffer[256];
	sprintf(buffer,"%s_%d_",s.c_str(),p);
	fNamePrefix = string(buffer);
}

void cellDivisionWithTemporalWindow::addFeatureName(string& s, int p, int ll)
{
	char buffer[256];
	sprintf(buffer,"%s_%d_tsLength%d",s.c_str(),p,ll);
	fName.push_back(fNamePrefix + string(buffer));
}


void cellDivisionWithTemporalWindow::writeFeaturesNames(string& filenameS)
{
	ofstream fs(filenameS);
	if( fs.is_open() == false )
	{
		cout<<"ERROR: at mainExtractFeatures: file "<<filenameS<<" could not be opened to save results"<<endl;
		exit(3);
	}
	
	

	for(vector<string>::const_iterator iterS = fName.begin(); iterS != fName.end(); ++iterS)
		fs<<*iterS<<endl;

	fs.close();

	cout<<"Saved "<<fName.size()<< "names of features at "<<filenameS<<endl;
}

#endif
//================================================================
int cellDivisionWithTemporalWindow::calculateHaarBasedFeaturesSingleWindowForDaughters()
{
	int nFeatures = 0;
	vector<float> timeSeries( getTotalTemporalWindowSizeSingleWindowForDaughters() );


	//temporal features for ring intensity
	for( int ii = 0; ii < basicEllipticalHaarFeatureVector::numRings; ii++ )
	{		
		//generate time series
		for( size_t jj = 0; jj < timeSeries.size(); jj++)
			timeSeries[jj] = fHaarVec[jj].ringAvgIntensity[ii];

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("ellipHaar3D_ring"),ii);
#endif
		//calculate features
		nFeatures += temporalFeatureExtractionSingleWindowForDaughters(timeSeries);				


		if( basicEllipticalHaarFeatureVector::useDoGfeatures == true ) // do the same for DoG
		{
			//generate time series
			for( size_t jj = 0; jj < timeSeries.size(); jj++)
				timeSeries[jj] = fHaarVec[jj].ringAvgIntensityDoG[ii];

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("ellipHaar3D_ring_DoG"),ii);
#endif
			//calculate features
			nFeatures += temporalFeatureExtractionSingleWindowForDaughters(timeSeries);

		}
	}
	//temporal features for radial cell intensity
	for( int ii = 0; ii < basicEllipticalHaarFeatureVector::numCells; ii++ )
	{
		//generate time series
		for( size_t jj = 0; jj < timeSeries.size(); jj++)
			timeSeries[jj] = fHaarVec[jj].cellAvgIntensity[ii];

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("ellipHaar3D_cell"),ii);
#endif
		//calculate features
		nFeatures += temporalFeatureExtractionSingleWindowForDaughters(timeSeries);


		if( basicEllipticalHaarFeatureVector::useDoGfeatures == true ) // do the same for DoG
		{
			//generate time series
			for( size_t jj = 0; jj < timeSeries.size(); jj++)
				timeSeries[jj] = fHaarVec[jj].cellAvgIntensityDoG[ii];

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("ellipHaar3D_cell_DoG"),ii);
#endif
			//calculate features
			nFeatures += temporalFeatureExtractionSingleWindowForDaughters(timeSeries);

		}
	}
	//temporal features for excentricity values
	for( int ii = 0; ii < dimsImage*(dimsImage-1)/2; ii++ )
	{
		//generate time series
		for( size_t jj = 0; jj < timeSeries.size(); jj++)
			timeSeries[jj] = fHaarVec[jj].excentricity[ii];
#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("excentricity"),ii);
#endif
		//calculate features
		nFeatures += temporalFeatureExtractionSingleWindowForDaughters(timeSeries);	
	}
	

	return nFeatures;
}



//======================================================================
int cellDivisionWithTemporalWindow::calculateLineageBasedFeatures()
{

	int nFeatures = 0;
	int nTS = getTotalTemporalWindowSize();
	float auxF;
	vector<float> timeSeries( nTS ); //preallocation
	TreeNode<ChildrenTypeLineage>* auxNode;

	//extract a vector of pointers to each point in the lineage around teh cell division (NULL indicates death or birth)
	vector< TreeNode<ChildrenTypeLineage>* > lineageVec;
	lineagePointersWithinTemporalWindowWithCellDivision( lineageVec );

	assert( lineageVec.size() == timeSeries.size());

	
	
	//----------------distance of all cells to the midplane defined by two daughters------------------------------
	float *centroidDL = cellDivisionPtr->left->data->centroid;
	float *centroidDR = cellDivisionPtr->right->data->centroid;

	for(int ii = 0; ii < nTS; ii++)
	{
		if( lineageVec[ii] == NULL )
			timeSeries[ii] = 100.0f;//default large number
		else
			timeSeries[ii] = lineageHyperTree::cellDivisionPlaneDistance(lineageVec[ii]->data->centroid, centroidDL, centroidDR);
	}

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("dist2midplane"),0);
#endif
	//calculate features
	nFeatures += temporalFeatureExtraction(timeSeries);


	//----------------distance between daughters------------------------------	
	int tt;
	for(tt = 0; tt <= temporalWindowRadius; tt++)
	{
		timeSeries[tt] = 0;//at the beginning there are no daughters
	}
	for(; tt <= 2 * temporalWindowRadius; tt++)
	{
		if( lineageVec[tt] == NULL || lineageVec[tt+ temporalWindowRadius] == NULL )
		{
			timeSeries[tt] = 0;
			timeSeries[tt+ temporalWindowRadius] = 0;
		}else{
			auxF = lineageVec[tt]->data->Euclidean2Distance( *(lineageVec[tt + temporalWindowRadius]->data), supervoxel::getScale() );
			timeSeries[tt] = auxF;
			timeSeries[tt+ temporalWindowRadius] = auxF;//I cna use this space for some other metric involving the two daughters
		}
	}
	#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("distBetweenDaughters"),0);
#endif
	//calculate features
	nFeatures += temporalFeatureExtraction(timeSeries);


	//----------------difference in volume (pixels) with respect to parent------------------------------	
	for(int ii = 0; ii < nTS; ii++)
	{
		auxNode = lineageVec[ii];
		if( auxNode != NULL && auxNode->parent != NULL )
		{
			auxF = lineageHyperTree::getNucleusVolume(*(auxNode->data));
			timeSeries[ii] = ( auxF - lineageHyperTree::getNucleusVolume(*(auxNode->parent->data)) ) / auxF;
		}else{
			timeSeries[ii] = 0.0f;
		}
	}

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("delta_volume"),0);
#endif
	//calculate features
	nFeatures += temporalFeatureExtraction(timeSeries);

	//----------------offset Jaccard distance with respect to parent------------------------------	
	for(int ii = 0; ii < nTS; ii++)
	{
		auxNode = lineageVec[ii];
		if( auxNode != NULL && auxNode->parent != NULL )
		{
			timeSeries[ii] = lineageHyperTree::JaccardDistance(auxNode, true);
		}else{
			timeSeries[ii] = 2.0f;//Jaccard distanc eis bounded by 1.0
		}
	}
#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("delta_Jaccard_wOffset"),0);
#endif
	//calculate features
	nFeatures += temporalFeatureExtraction(timeSeries);

	//----------------Jaccard distance with respect to parent------------------------------	
	for(int ii = 0; ii < nTS; ii++)
	{
		auxNode = lineageVec[ii];
		if( auxNode != NULL && auxNode->parent != NULL )
		{
			timeSeries[ii] = lineageHyperTree::JaccardDistance(auxNode, false);
		}else{
			timeSeries[ii] = 2.0f;//Jaccard distanc eis bounded by 1.0
		}
	}

#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("delta_Jaccard"),0);
#endif
	//calculate features
	nFeatures += temporalFeatureExtraction(timeSeries);

	//----------------absolute displacement with respect to parent------------------------------	
	for(int ii = 0; ii < nTS; ii++)
	{
		auxNode = lineageVec[ii];
		if( auxNode != NULL && auxNode->parent != NULL )
		{
			timeSeries[ii] = auxNode->data->Euclidean2Distance( *(auxNode->parent->data), supervoxel::getScale() );
		}else{
			timeSeries[ii] = 100.0f;//deafult last number
		}
	}
	#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("delta_displacement"),0);
#endif
	//calculate features
	nFeatures += temporalFeatureExtraction(timeSeries);

	//----------------offset Jaccard distance between daughters------------------------------	
	for(tt = 0; tt <= temporalWindowRadius; tt++)
	{
		if( lineageVec[temporalWindowRadius] == NULL || lineageVec[tt+ temporalWindowRadius + 1] == NULL )
			timeSeries[tt] = 0;
		else
		{
			//timeSeries[tt] = 0;//at the beginning there are no daughters. So we calculate Jaccard with respect to the first daughter
			timeSeries[tt] = lineageHyperTree::JaccardDistance(lineageVec[temporalWindowRadius], lineageVec[tt + temporalWindowRadius + 1], true );
		}
	}
	for(; tt <= 2 * temporalWindowRadius; tt++)
	{
		if( lineageVec[tt] == NULL || lineageVec[tt+ temporalWindowRadius] == NULL )
		{
			timeSeries[tt] = 0;
			timeSeries[tt+ temporalWindowRadius] = 0;
		}else{
			timeSeries[tt] = lineageHyperTree::JaccardDistance(lineageVec[tt], lineageVec[tt + temporalWindowRadius], true );
			timeSeries[tt+ temporalWindowRadius] = lineageHyperTree::JaccardDistance(lineageVec[temporalWindowRadius], lineageVec[tt + temporalWindowRadius], true );//I use this space to calculate Jaccard with respect to second daughter
		}
	}

	#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("Jaccard_offset_betweenDaughters"),0);
#endif
	//calculate features
	nFeatures += temporalFeatureExtraction(timeSeries);

	//----------------volume ratio between daughters------------------------------	
	for(tt = 0; tt <= temporalWindowRadius; tt++)
	{
		if( lineageVec[temporalWindowRadius] == NULL || lineageVec[tt+ temporalWindowRadius + 1] == NULL )
			timeSeries[tt] = 0;
		else
		{
			//timeSeries[tt] = 0;//at the beginning there are no daughters. So we calculate Jaccard with respect to the first daughter
			timeSeries[tt] = lineageHyperTree::getNucleusVolume(*(lineageVec[temporalWindowRadius]->data)) / lineageHyperTree::getNucleusVolume(*(lineageVec[tt + temporalWindowRadius + 1]->data));
		}
	}
	for(; tt <= 2 * temporalWindowRadius; tt++)
	{
		if( lineageVec[tt] == NULL || lineageVec[tt+ temporalWindowRadius] == NULL )
		{
			timeSeries[tt] = 0;
			timeSeries[tt+ temporalWindowRadius] = 0;
		}else{
			auxF = lineageHyperTree::getNucleusVolume(*(lineageVec[tt + temporalWindowRadius]->data));
			timeSeries[tt] = lineageHyperTree::getNucleusVolume(*(lineageVec[tt]->data)) / auxF;
			timeSeries[tt + temporalWindowRadius] = lineageHyperTree::getNucleusVolume(*(lineageVec[temporalWindowRadius]->data)) / auxF;//I use this space to calculate volume with respect to second daughter
			
		}
	}

	#ifdef CDWT_SAVE_FEATURE_NAME	
		setFeaturePrefix(string("volume_ratio_betweenDaughters"),0);
#endif
	//calculate features
	nFeatures += temporalFeatureExtraction(timeSeries);
	

	return nFeatures;
}

//=====================================================================
int cellDivisionWithTemporalWindow::temporalFeatureExtraction(vector<float>& timeSeries, bool isRecursiveCall, bool calculateMultiScale)
{
	int nFeatures = 0;
	const int numStats = 4;//so the code can be expanded later on

	//recalculate tmeporal window radius locally based on timeSeries.size(). This is needed since Imight call it recursively for different time window radius (just using the largest window as a precursor)
	int temporalWindowRadiusLocal = ( timeSeries.size() - 1 ) / 3;
	assert( temporalWindowRadiusLocal*3 + 1 == timeSeries.size() );

	//subsample if necessary
	vector<float> timeSeriesOrig(timeSeries);
	if( temporalUndersampling > 1 && isRecursiveCall == false)//I only subsample at the beginning
	{
		timeSeries.clear();
		for(int ii = 0; ii <= 2 * temporalWindowRadiusLocal ; ii+=temporalUndersampling)
			timeSeries.push_back( timeSeriesOrig[ii] );
		//add second daughter
		for(int ii = 2 * temporalWindowRadiusLocal+1; ii <= 3*temporalWindowRadiusLocal; ii+=temporalUndersampling)
			timeSeries.push_back( timeSeriesOrig[ii] );

		temporalWindowRadiusLocal = ( timeSeries.size() - 1 ) / 3;
	}


	if( temporalWindowRadiusLocal < 2 )//we need at least two, so we have two elements on each of the daughter branches
		return 0;

	double stats[numStats];//0 = mean; 1 = std; 2 = min; 3 = max 

	
	utilsAmatf_meanAndStd<float>(timeSeries, &(stats[0]), &(stats[1]));
	stats[2] = *(std::min_element(timeSeries.begin(), timeSeries.end()) );
	stats[3] = *(std::max_element(timeSeries.begin(), timeSeries.end()) );

	for( int ii = 0; ii < numStats; ii++)
	{
		nFeatures++;
		f.push_back( stats[ii] );
#ifdef CDWT_SAVE_FEATURE_NAME	
		addFeatureName(string("stats"),ii, timeSeries.size());
#endif
	}
	

	//create three subvectors to separate mother, daughter1 and daughter2
	vector<float> timeSeries_MDD[3];
	timeSeries_MDD[0] = vector<float>(timeSeries.begin(), timeSeries.begin()+temporalWindowRadiusLocal+1);
	timeSeries_MDD[1] = vector<float>(timeSeries.begin()+ temporalWindowRadiusLocal+1, timeSeries.begin()+ 2 * temporalWindowRadiusLocal + 1);
	timeSeries_MDD[2] = vector<float>(timeSeries.begin()+ 2*temporalWindowRadiusLocal+1, timeSeries.end());
	
	
	
	double stats_MDD[numStats][3];//0 = mean; 1 = std; 2 = min; 3 = max 
	for( int ii = 0; ii < 3; ii++)
	{
		utilsAmatf_meanAndStd<float>(timeSeries_MDD[ii], &(stats_MDD[0][ii]), &(stats_MDD[1][ii]));
		stats_MDD[2][ii] = *(std::min_element(timeSeries_MDD[ii].begin(), timeSeries_MDD[ii].end()) );
		stats_MDD[3][ii] = *(std::max_element(timeSeries_MDD[ii].begin(), timeSeries_MDD[ii].end()) );
	}

	//calculate all combinations of {+1,0,-1} (3^3 -1 / 2) because half of them are just the same with different sign
	float fAux;
	int countT;
	for( int ss = 0; ss< numStats; ss++)
	{
		countT = 0;
		for( int ii = 0; ii < HaarPatternTripletN; ii++)
		{
			
			fAux = 0;
			for( int jj = 0; jj < 3; jj++)
				fAux += HaarPatternTriplet[countT++] * stats_MDD[ss][jj];
			nFeatures++;
			f.push_back(fAux);
#ifdef CDWT_SAVE_FEATURE_NAME	
		addFeatureName(string("stats_MDD_HaarTriplet"),ss, timeSeries.size());
#endif
		}
	}

	//calculate value at cell division time vs neighborhoods	
	const float fAtCellDivision = timeSeries[temporalWindowRadiusLocal];
	const float epsilonF = 1e-3;//to avoid division by zero in some cases
	for( int ss = 0; ss< numStats; ss++)
	{
		nFeatures++;
		f.push_back( fAtCellDivision / (stats[ss] + epsilonF) );
#ifdef CDWT_SAVE_FEATURE_NAME	
		addFeatureName(string("stats_neighbor"),ss, timeSeries.size());
#endif


		for(int ii = 0; ii < 3; ii++ )//3 segments: D1, D2, Mother
		{
			nFeatures++;
			f.push_back( fAtCellDivision / (stats_MDD[ss][ii] + epsilonF) );
#ifdef CDWT_SAVE_FEATURE_NAME	
		addFeatureName(string("stats_MDD_neighbor"),ss*3 + ii, timeSeries.size());
#endif

		}
	}

	//zscore
	nFeatures++;
	f.push_back( (fAtCellDivision-stats[0]) / (stats[1] + epsilonF) );
	#ifdef CDWT_SAVE_FEATURE_NAME	
		addFeatureName(string("zscore"),0, timeSeries.size());
#endif

	for(int ii = 0; ii < 3; ii++ )//3 segments: D1, D2, Mother
	{
		nFeatures++;
		f.push_back( (fAtCellDivision-stats_MDD[0][ii]) / (stats_MDD[1][ii] + epsilonF) );
		#ifdef CDWT_SAVE_FEATURE_NAME	
		addFeatureName(string("zscore_MDD"),ii, timeSeries.size());
#endif

	}


	//recursive to calculate multi-scale features (in case temporal sampling is different for different recordings)
	if( calculateMultiScale == true && temporalWindowRadiusLocal > 3 )//for example for temporalWindowRadius = 5 -> this generates scales 3
	{
		int newRadius = pow(2,floor(logbase(temporalWindowRadiusLocal,2))) -1;

		vector<float> newTimeSeries(3*newRadius+1);
		int count = 0;
		for(int ii = temporalWindowRadiusLocal-newRadius; ii <= temporalWindowRadiusLocal + newRadius; ii++)
			newTimeSeries[count++] = timeSeries[ii];
		//add second daughter
		for(int ii = 2 * temporalWindowRadiusLocal+1; ii <= 2*temporalWindowRadiusLocal + newRadius; ii++)
			newTimeSeries[count++] = timeSeries[ii];

		//call the function recursively
		nFeatures += temporalFeatureExtraction(newTimeSeries, true);
	}

	//add each value of teh time series itself (per-frame feature) to the vector
	if( isRecursiveCall == false )
	{
		//add the features for each single time point
		f.insert( f.end(), timeSeries.begin(), timeSeries.end());
		nFeatures += timeSeries.size();

#ifdef CDWT_SAVE_FEATURE_NAME	
		for( size_t ii = 0; ii < timeSeries.size(); ii++)
			addFeatureName(string("timeSeries"),ii, timeSeries.size());
#endif


		//add the features for each single time point normalized by stats
		for(int ss = 0; ss < numStats; ss++)
		{
			for( size_t ii = 0; ii < timeSeries.size(); ii++)
			{
				f.push_back(timeSeries[ii] / (stats[ss] + epsilonF));
				nFeatures++;
#ifdef CDWT_SAVE_FEATURE_NAME			
			addFeatureName(string("timeSeries_normalized"),ss*3 + ii, timeSeries.size());
#endif
			}
		}
	}

	//restore from subsampling	
	if( timeSeries.size() != timeSeriesOrig.size() )//I only subsample at the beginning
	{
		timeSeries.swap( timeSeriesOrig );
	}

	return nFeatures;
}


//=====================================================================
int cellDivisionWithTemporalWindow::temporalFeatureExtractionSingleWindowForDaughters(vector<float>& timeSeries, bool isRecursiveCall, bool calculateMultiScale)
{
	int nFeatures = 0;
	const int numStats = 4;//so the code can be expanded later on

	//recalculate tmeporal window radius locally based on timeSeries.size(). This is needed since Imight call it recursively for different time window radius (just using the largest window as a precursor)
	int temporalWindowRadiusLocal = ( timeSeries.size() - 1 ) / 2;
	assert( temporalWindowRadiusLocal*2 + 1 == timeSeries.size() );


	//subsample if necessary
	vector<float> timeSeriesOrig(timeSeries);
	if( temporalUndersampling > 1 && isRecursiveCall == false)//I only subsample at the beginning
	{
		timeSeries.clear();
		for(int ii = 0; ii < timeSeriesOrig.size(); ii+= temporalUndersampling)
		{
			timeSeries.push_back( timeSeriesOrig[ii] );
		}
		temporalWindowRadiusLocal = ( timeSeries.size() - 1 ) / 2;
	}

	if( temporalWindowRadiusLocal < 2 )//we need at least two, so we have two elements on each of the daughter branches
		return 0;

	double stats[numStats];//0 = mean; 1 = std; 2 = min; 3 = max 

	
	utilsAmatf_meanAndStd<float>(timeSeries, &(stats[0]), &(stats[1]));
	stats[2] = *(std::min_element(timeSeries.begin(), timeSeries.end()) );
	stats[3] = *(std::max_element(timeSeries.begin(), timeSeries.end()) );

	for( int ii = 0; ii < numStats; ii++)
	{
		nFeatures++;
		f.push_back( stats[ii] );
#ifdef CDWT_SAVE_FEATURE_NAME			
			addFeatureName(string("SD_stats"),ii, timeSeries.size());
#endif
	}
	

	//create three subvectors to separate the time series in 2 segments: before and after division
	vector<float> timeSeries_MDD[2];
	timeSeries_MDD[0] = vector<float>(timeSeries.begin(), timeSeries.begin()+temporalWindowRadiusLocal+1);
	timeSeries_MDD[1] = vector<float>(timeSeries.begin()+ temporalWindowRadiusLocal+1, timeSeries.end());
	
	
	
	double stats_MDD[numStats][2];//0 = mean; 1 = std; 2 = min; 3 = max 
	for( int ii = 0; ii < 2; ii++)
	{
		utilsAmatf_meanAndStd<float>(timeSeries_MDD[ii], &(stats_MDD[0][ii]), &(stats_MDD[1][ii]));
		stats_MDD[2][ii] = *(std::min_element(timeSeries_MDD[ii].begin(), timeSeries_MDD[ii].end()) );
		stats_MDD[3][ii] = *(std::max_element(timeSeries_MDD[ii].begin(), timeSeries_MDD[ii].end()) );
	}

	//calculate all combinations of {+1,0,-1} (3^2 -1 / 2) because half of them are just the same with different sign
	float fAux;
	int countT;
	for( int ss = 0; ss< numStats; ss++)
	{
		countT = 0;
		for( int ii = -1; ii <= 1; ii++)
		{
			for(int kk = ii; kk <= 1; kk++)//we generate [0,0] and [1,1] which are redundant, but the code is simple
			{
				fAux = stats_MDD[ss][0] * ii + stats_MDD[ss][1] * kk;				
				nFeatures++;
				f.push_back(fAux);
				#ifdef CDWT_SAVE_FEATURE_NAME			
			addFeatureName(string("SD_stats_HaarTriplet"),ss, timeSeries.size());
#endif
			}
		}
	}

	//calculate value at cell division time vs neighborhoods	
	const float fAtCellDivision = timeSeries[temporalWindowRadiusLocal];
	const float epsilonF = 1e-3;//to avoid division by zero in some cases
	for( int ss = 0; ss< numStats; ss++)
	{
		nFeatures++;
		f.push_back( fAtCellDivision / (stats[ss] + epsilonF) );
#ifdef CDWT_SAVE_FEATURE_NAME			
			addFeatureName(string("SD_neighbor"),ss, timeSeries.size());
#endif

		for(int ii = 0; ii < 2; ii++ )//2 segments: D12, Mother
		{
			nFeatures++;
			f.push_back( fAtCellDivision / (stats_MDD[ss][ii] + epsilonF) );

	#ifdef CDWT_SAVE_FEATURE_NAME			
			addFeatureName(string("SD_neighbor_MDD"),ss*3 + ii, timeSeries.size());
#endif
		}
	}

	//zscore
	nFeatures++;
	f.push_back( (fAtCellDivision-stats[0]) / (stats[1] + epsilonF) );
	#ifdef CDWT_SAVE_FEATURE_NAME			
			addFeatureName(string("SD_zscore"),0, timeSeries.size());
#endif
	for(int ii = 0; ii < 2; ii++ )//2 segments: D12, Mother
	{
		nFeatures++;
		f.push_back( (fAtCellDivision-stats_MDD[0][ii]) / (stats_MDD[1][ii] + epsilonF) );

#ifdef CDWT_SAVE_FEATURE_NAME			
			addFeatureName(string("SD_zscore_MDD"),ii, timeSeries.size());
#endif
	}


	//recursive to calculate multi-scale features (in case temporal sampling is different for different recordings)
	if( calculateMultiScale == true && temporalWindowRadiusLocal > 3 )//for example for temporalWindowRadius = 5 -> this generates scales 3
	{
		int newRadius = pow(2,floor(logbase(temporalWindowRadiusLocal,2))) -1;

		vector<float> newTimeSeries(2*newRadius+1);
		int count = 0;
		for(int ii = temporalWindowRadiusLocal-newRadius; ii <= temporalWindowRadiusLocal + newRadius; ii++)
			newTimeSeries[count++] = timeSeries[ii];		

		//call the function recursively
		if( temporalWindowRadiusLocal == 5 )
			nFeatures += temporalFeatureExtraction(newTimeSeries, true);//this was initially a bug, but it splits the new timeSeries into 3,2,2 segments and the features seem to be useful
		else
			nFeatures += temporalFeatureExtractionSingleWindowForDaughters(newTimeSeries, true);

	}

	//add each value of teh time series itself (per-frame feature) to the vector
	if( isRecursiveCall == false )
	{
		//add the features for each single time point
		f.insert( f.end(), timeSeries.begin(), timeSeries.end());
		nFeatures += timeSeries.size();

#ifdef CDWT_SAVE_FEATURE_NAME	
		for( size_t ii = 0; ii < timeSeries.size(); ii++)
			addFeatureName(string("SD_timeSeries"),ii, timeSeries.size());
#endif

		//add the features for each single time point normalized by stats
		for(int ss = 0; ss < numStats; ss++)
		{
			for( size_t ii = 0; ii < timeSeries.size(); ii++)
			{
				f.push_back(timeSeries[ii] / (stats[ss] + epsilonF));
				nFeatures++;
#ifdef CDWT_SAVE_FEATURE_NAME			
			addFeatureName(string("SD_timeSeries_normalized"),ss*3 + ii, timeSeries.size());
#endif
			}
		}
	}

	//restore from subsampling	
	if( timeSeries.size() != timeSeriesOrig.size() )//I only subsample at the beginning
	{
		timeSeries.swap( timeSeriesOrig );
	}

	return nFeatures;
}


//==========================================================
void cellDivisionWithTemporalWindow::averagePrecisionMatrixAyTM(const list<nucleus>& listN)
{
	int sizeW = dimsImage * (1+dimsImage) / 2;
	memset(precisionWfix, 0, sizeof(double) * sizeW);

	for(list<nucleus>::const_iterator iterN = listN.begin(); iterN != listN.end(); ++iterN)
	{
		for( int ii = 0; ii < sizeW; ii++)
			precisionWfix[ii] += iterN->precisionW[ii];
	}

	bool flagErr = false;
	for( int ii = 0; ii < sizeW; ii++)
	{
		precisionWfix[ii] /= listN.size();

#if defined(_WIN32) || defined(_WIN64)
		if( _finite( precisionWfix[ii] ) == false )//degenerated cases
#else
		if( isfinite( precisionWfix[ii] ) == 0 )//degenerated cases
#endif
			flagErr = true;
	}

	if( flagErr == true )//reset to default
	{
		memset(precisionWfix, 0, sizeof(double) * sizeW);
		precisionWfix[0] = 0.01;
		precisionWfix[dimsImage] = 0.01;
		if( dimsImage > 2 )
			precisionWfix[sizeW-1] = 0.01;
	}	
}


//==========================================================
int cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchAtTM(const vector< bool >& isDead, mylib::Array* img, vector<cellDivisionWithTemporalWindow>& cdtwVec, int relativeTimeWithinWindow,int devCUDA, int symmetry, double *m, double *W, long long int **dimsVec)
{
	int numEllipsoids = isDead.size();
	int sizeW = dimsImage * (1+dimsImage) / 2;
	TreeNode< ChildrenTypeLineage >* auxNode = NULL;

	//allocate memory for dimsVec
	if( *dimsVec == NULL )
		*dimsVec = new long long int[dimsImage];



	//----------------in order to minimize memory comsumption during TGMM execution------------------
	//the problem is that I trained in raw data (without local background subtraction using supervoxels + otsu, so I am not using this option anymore because the data does not look like the training data
	bool imgWasNull = false;
    CHECK(img!=NULL);
#if 0
	if( img == NULL )
	{
		imgWasNull = true;
		supervoxel* sv = &(*(cdtwVec[0].cellDivisionPtr->data->treeNode.getChildren().front()));
		mylib::Dimn_Type dimsAux[dimsImage];
		uint64 sizeImg = 1;
		for(int ii = 0; ii < dimsImage; ii++)
		{
			dimsAux[ii] = sv->dataDims[ii];
			sizeImg *= sv->dataDims[ii];
		}
		img = mylib::Make_Array(mylib::PLAIN_KIND, mylib::UINT16_TYPE, dimsImage, dimsAux);
		mylib::uint16* imgPtr = (mylib::uint16*) (img->data);

		uint64 typeImg = sv->dataSizeInBytes / sizeImg;

		switch( typeImg )
		{
		case 1://unsigned char 
			{
				unsigned char* svImg = (unsigned char*) (sv->dataPtr);
				for(uint64 ii = 0; ii < sizeImg; ii++)
					imgPtr[ii] = (mylib::uint16) (svImg[ii]);
			}
			break;

		case 2://uint16
			{
				mylib::uint16* svImg = (mylib::uint16*) (sv->dataPtr);
				for(uint64 ii = 0; ii < sizeImg; ii++)
					imgPtr[ii] = (svImg[ii]);
			}
			break;
		case 4://float32. We assume image is normalized between 0 and 1
			{
				mylib::float32* svImg = (mylib::float32*) (sv->dataPtr);
				for(uint64 ii = 0; ii < sizeImg; ii++)
					imgPtr[ii] = (mylib::uint16)(svImg[ii] * 65535);
			}
			break;
		default:
			cout<<"ERROR:cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchAtTM: code not ready for these kind of images"<<endl;
		}
	}
#endif
	//------------------------------------------------------------------------------
	
    CHECK(img->type == mylib::UINT16_TYPE);

	mylib::uint16* imgPtr = (mylib::uint16*) (img->data);		
	for(int aa = 0; aa < img->ndims; aa++)
		(*dimsVec)[aa] = img->dims[aa];
	
	//right now there is a lot of memory allocation/deallocation and copying being made
	//TODO: pass basicEllipticalHaarFeatureVector **fBasic as an argument so we can reuse it;
	basicEllipticalHaarFeatureVector **fBasic = calculateEllipticalHaarFeatures(m, W,numEllipsoids,imgPtr,*dimsVec,devCUDA,symmetry);
	if(fBasic == NULL)
		return 1;//some error happened

	//copy features and release memory
	for(int ii = 0; ii < numEllipsoids; ii++)
	{
		if( isDead[ii] == false )
		{
			cdtwVec[ii].fHaarVec[relativeTimeWithinWindow] = (*(fBasic[ii]));			
		}else{
			cdtwVec[ii].fHaarVec[relativeTimeWithinWindow].setToDefault();
		}
		delete fBasic[ii];
	}	

	delete[] fBasic;	
	if( imgWasNull == true )
	{
		mylib::Free_Array(img);
		img = NULL;
	}

	return 0;

}



//==========================================================
int cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchAtTM(const vector< TreeNode< ChildrenTypeLineage >* >& auxNodeVec, mylib::Array* img, vector<cellDivisionWithTemporalWindow>& cdtwVec, int relativeTimeWithinWindow,int devCUDA, int symmetry, double **m , double **W, long long int **dimsVec)
{	

	
	int numEllipsoids = auxNodeVec.size();
	int sizeW = dimsImage * (1+dimsImage) / 2;
	TreeNode< ChildrenTypeLineage >* auxNode = NULL;

	if( numEllipsoids == 0 )
		return 0;

	//allocate memory for centroid and precision matrix
	if( *m == NULL )
		*m = new double[dimsImage * numEllipsoids];
	if( *W == NULL )
		*W = new double[sizeW  * numEllipsoids];
	if( *dimsVec == NULL )
		*dimsVec = new long long int[dimsImage];


	//copy centroid and precision matrix
	if( useFixPrecisionMatrix == false )
	{
		for(int ii = 0; ii < numEllipsoids; ii++)
		{			
			auxNode = auxNodeVec[ii];
			if( auxNode != NULL )//there might be irths
			{
				for(int aa = 0; aa<dimsImage; aa++)
				{
					(*m)[ii + aa * numEllipsoids] = auxNode->data->centroid[aa];
				}
				for(int aa = 0; aa<sizeW; aa++)
				{
					(*W)[ii + aa * numEllipsoids] = auxNode->data->precisionW[aa];
				}
			}
		}
	}else{
		for(int ii = 0; ii < numEllipsoids; ii++)
		{			
			auxNode = auxNodeVec[ii];
			if( auxNode != NULL )//there might be irths
			{
				for(int aa = 0; aa<dimsImage; aa++)
				{
					(*m)[ii + aa * numEllipsoids] = auxNode->data->centroid[aa];
				}
				for(int aa = 0; aa<sizeW; aa++)
				{
					(*W)[ii + aa * numEllipsoids] = precisionWfix[aa]; 
				}
			}
		}
	}

	
	//----------------in order to minimize memory comsumption during TGMM execution------------------
	bool imgWasNull = false;
    // (ngc) not sure this code path is ever used
    CHECK(img!=NULL);
#if 0
	if( img == NULL )
	{
		imgWasNull = true;
		supervoxel* sv = &(*(cdtwVec[0].cellDivisionPtr->data->treeNode.getChildren().front()));
		mylib::Dimn_Type dimsAux[dimsImage];
		uint64 sizeImg = 1;
		for(int ii = 0; ii < dimsImage; ii++)
		{
			dimsAux[ii] = sv->dataDims[ii];
			sizeImg *= sv->dataDims[ii];
		}
		img = mylib::Make_Array(mylib::PLAIN_KIND, mylib::UINT16_TYPE, dimsImage, dimsAux);
		mylib::uint16* imgPtr = (mylib::uint16*) (img->data);

		uint64 typeImg = sv->dataSizeInBytes / sizeImg;

		switch( typeImg )
		{
		case 1://unsigned char 
			{
				unsigned char* svImg = (unsigned char*) (sv->dataPtr);
				for(uint64 ii = 0; ii < sizeImg; ii++)
					imgPtr[ii] = (mylib::uint16) (svImg[ii]);
			}
			break;

		case 2://uint16
			{
				mylib::uint16* svImg = (mylib::uint16*) (sv->dataPtr);
				for(uint64 ii = 0; ii < sizeImg; ii++)
					imgPtr[ii] = (svImg[ii]);
			}
			break;
		case 4://float32. We assume image is normalized between 0 and 1
			{
				mylib::float32* svImg = (mylib::float32*) (sv->dataPtr);
				for(uint64 ii = 0; ii < sizeImg; ii++)
					imgPtr[ii] = (mylib::uint16)(svImg[ii] * 65535);
			}
			break;
		default:
			cout<<"ERROR:cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchAtTM: code not ready for these kind of images"<<endl;
		}

	}
#endif
	//------------------------------------------------------------------------------

	
	CHECK( img->type == mylib::UINT16_TYPE );
	mylib::uint16* imgPtr = (mylib::uint16*) (img->data);		
	for(int aa = 0; aa < img->ndims; aa++)
		(*dimsVec)[aa] = img->dims[aa];
	

	//right now there is a lot of memory allocation/deallocation and copying being made
	//TODO: pass basicEllipticalHaarFeatureVector **fBasic as an argument so we can reuse it;
	basicEllipticalHaarFeatureVector **fBasic = calculateEllipticalHaarFeatures(*m, *W,numEllipsoids,imgPtr,*dimsVec,devCUDA,symmetry);
	if(fBasic == NULL)
		return 1;//some error happened

	//copy features and release memory
	for(int ii = 0; ii < numEllipsoids; ii++)
	{
		if( auxNodeVec[ii] != NULL )
		{
			cdtwVec[ii].fHaarVec[relativeTimeWithinWindow] = (*(fBasic[ii]));			
		}else{
			cdtwVec[ii].fHaarVec[relativeTimeWithinWindow].setToDefault();
		}
		delete fBasic[ii];
	}

	delete[] fBasic;
	if( imgWasNull == true )
	{
		mylib::Free_Array(img);
		img = NULL;
	}

	
	//recalculate excentricity with the real W matrix to preserve this shape information
	if( useFixPrecisionMatrix == true )
	{
		int count;
		float dAux[dimsImage];//eigenvalues
		for(int ii = 0; ii < numEllipsoids; ii++)
		{			
			auxNode = auxNodeVec[ii];
			if( auxNode != NULL )//there might be irths
			{
				utilsAmatf_eig3<float>(auxNode->data->precisionW, dAux, NULL);

				//calculate excentricity 
				count = 0;
				for(int jj = 0; jj < dimsImage; jj++)
				{
					for(int kk = jj + 1; kk < dimsImage; kk++)
					{
						if(dAux[kk]<1e-10)
							cdtwVec[ii].fHaarVec[relativeTimeWithinWindow].excentricity[count] = 0.0f;
						else
							cdtwVec[ii].fHaarVec[relativeTimeWithinWindow].excentricity[count] = dAux[jj]/dAux[kk];

						count++;
					}
				}
			}
		}
	}

	return 0;

}

//=====================================================================================
void cellDivisionWithTemporalWindow::lineagePointersWithinTemporalWindowWithCellDivision(vector< TreeNode<ChildrenTypeLineage>* >&  lineageVec)
{
	lineageVec.resize( getTotalTemporalWindowSize() );
	
	assert( cellDivisionPtr->getNumChildren() == 2);

	if( cellDivisionPtr == NULL )
	{
		cout<<"WARNING:cellDivisionWithTemporalWindow::lineagePointersWithinTemporalWindowWithCellDivision: cellDivisionPtr = NULL"<<endl;
		return;
	}

	
	
	//trace backwards (including cell division point)
	TreeNode<ChildrenTypeLineage>* auxNode = cellDivisionPtr;
	for(int tt = temporalWindowRadius; tt >= 0; tt--)
	{
		lineageVec[tt] = auxNode;

		if( auxNode != NULL )
			auxNode = auxNode->parent;		
	}

	//trace forward left daughter
	int ttOffset = temporalWindowRadius + 1;//to save in the correct place in fHaarVec
	for(int ch = 0; ch < 2; ch++)
	{
		if( ch == 0 )
		{
			auxNode = cellDivisionPtr->left;
		}else{
			auxNode = cellDivisionPtr->right;
		}
		//iterate over time points
		for( int tt = temporalWindowRadius + 1; tt < 2* temporalWindowRadius + 1; tt++)
		{
			lineageVec[ttOffset] = auxNode;
			auxNode = moveForwardInlineageEvenWithCellDivision(auxNode);			
			ttOffset++;
		}
	}

}
//===================================================================================================================
void cellDivisionWithTemporalWindow::calculatePrecisionMatrixForJointDuaghters(nucleus& DL, nucleus& DR, double W[dimsImage * (1+dimsImage/2)])
{
	
	if(dimsImage != 3)
	{
		cout<<"ERROR: cellDivisionWithTemporalWindow::calculatePrecisionMatrixForJointDuaghters: code only ready for dimsImage =3 in order to be optimized"<<endl;
		exit(3);
	}


	double N_k = 0.0;
	double m_k[dimsImage];
	double S_k[dimsImage * (1+dimsImage) / 2];
	memset(m_k,0, sizeof(double) * dimsImage);
	//memset(W,0, sizeof(double) * dimsImage * (1 + dimsImage) / 2);
	memset(S_k,0, sizeof(double) * dimsImage * (1 + dimsImage) / 2);
	
	int64 coord[dimsImage];
	
	int count;
	nucleus* auxNuc[2];
	auxNuc[0] = &DL;
	auxNuc[1] = &DR;

	for(int hh = 0; hh < 2; hh++)
	{
		for(vector< ChildrenTypeNucleus >::iterator iterS = auxNuc[hh]->treeNode.getChildren().begin(); iterS != auxNuc[hh]->treeNode.getChildren().end(); ++iterS)
		{						
			int64 coordAux;
			for(vector<uint64>::iterator iterP = (*iterS)->PixelIdxList.begin(); iterP != (*iterS)->PixelIdxList.end(); ++iterP)
			{
				coordAux = (*iterP);
				for(int aa = 0; aa<dimsImage - 1; aa++)
				{
					coord[aa] = coordAux % (supervoxel::dataDims[aa]);
					coordAux -= coord[aa];
					coordAux /= (supervoxel::dataDims[aa]);
				}
				coord[dimsImage - 1] = coordAux;

				N_k ++;
				count =0;
				for(int ii = 0; ii<dimsImage; ii++)
				{
					m_k[ii] += coord[ii];
					for(int jj = ii; jj <dimsImage; jj++)
					{
						S_k[count] += coord[ii] * coord[jj];
						count++;
					}
				}
			}
		}
	}
	//finish calculating sufficient statistics
	for(int ii = 0; ii<dimsImage; ii++)
		m_k[ii] /= N_k;

	count =0;
	for(int ii = 0; ii<dimsImage; ii++)
	{
		for(int jj = ii; jj <dimsImage; jj++)
		{
			S_k[count] = S_k[count] / N_k - m_k[ii] * m_k[jj]; 
			count++;
		}
	}
	utilsAmatf_inverseSymmetricW_3D(S_k, W);
}




//=========================================================================================================================
template<class imgTypeC>
void cellDivisionWithTemporalWindow::calculatePrecisionMatrixForJointDuaghters(const TimeSeriesMapT& time_series_map,nucleus& DL, nucleus& DR, double m_k[dimsImage], double W[dimsImage * (1+dimsImage/2)])
{
	
	if(dimsImage != 3)
	{
		cout<<"ERROR: cellDivisionWithTemporalWindow::calculatePrecisionMatrixForJointDuaghters: code only ready for dimsImage =3 in order to be optimized"<<endl;
		exit(3);
	}


	double N_k = 0.0;
	double S_k[dimsImage * (1+dimsImage) / 2];
	memset(m_k,0, sizeof(double) * dimsImage);
	//memset(W,0, sizeof(double) * dimsImage * (1 + dimsImage) / 2);
	memset(S_k,0, sizeof(double) * dimsImage * (1 + dimsImage) / 2);
	
	int64 coord[dimsImage];
	
	int count;
	nucleus* auxNuc[2];
	auxNuc[0] = &DL;
	auxNuc[1] = &DR;

	for(int hh = 0; hh < 2; hh++)
	{
		for(vector< ChildrenTypeNucleus >::iterator iterS = auxNuc[hh]->treeNode.getChildren().begin(); iterS != auxNuc[hh]->treeNode.getChildren().end(); ++iterS)
		{
            imgTypeC *imgPtr=0;
            {
                const auto TM=(*iterS)->TM;
                const auto im=time_series_map.at(TM);
                CHECK(im->type==MylibValueType<imgTypeC>());
                imgPtr=(imgTypeC*)time_series_map.at(TM)->data;
            }
			mylib::float32 imgVal;
			int64 coordAux;
			for(vector<uint64>::iterator iterP = (*iterS)->PixelIdxList.begin(); iterP != (*iterS)->PixelIdxList.end(); ++iterP)
			{
				coordAux = (*iterP);
				imgVal = imgPtr[coordAux];
				for(int aa = 0; aa<dimsImage - 1; aa++)
				{
					coord[aa] = coordAux % (supervoxel::dataDims[aa]);
					coordAux -= coord[aa];
					coordAux /= (supervoxel::dataDims[aa]);
				}
				coord[dimsImage - 1] = coordAux;

				N_k += imgVal;
				count =0;
				for(int ii = 0; ii<dimsImage; ii++)
				{
					m_k[ii] += imgVal * coord[ii];
					for(int jj = ii; jj <dimsImage; jj++)
					{
						S_k[count] += imgVal * coord[ii] * coord[jj];
						count++;
					}
				}
			}
		}
	}
	//finish calculating sufficient statistics
	for(int ii = 0; ii<dimsImage; ii++)
		m_k[ii] /= N_k;

	count =0;
	for(int ii = 0; ii<dimsImage; ii++)
	{
		for(int jj = ii; jj <dimsImage; jj++)
		{
			S_k[count] = S_k[count] / N_k - m_k[ii] * m_k[jj]; 
			count++;
		}
	}
	utilsAmatf_inverseSymmetricW_3D(S_k, W);

}



//==============================================================================
//predefine templates
template void cellDivisionWithTemporalWindow::calculatePrecisionMatrixForJointDuaghters<float>(const TimeSeriesMapT& time_series_map, nucleus& DL,  nucleus& DR, double m_k[dimsImage], double W[dimsImage * (1+dimsImage/2)]);
template void cellDivisionWithTemporalWindow::calculatePrecisionMatrixForJointDuaghters<unsigned char>(const TimeSeriesMapT& time_series_map, nucleus& DL,  nucleus& DR, double m_k[dimsImage], double W[dimsImage * (1+dimsImage/2)]);
template void cellDivisionWithTemporalWindow::calculatePrecisionMatrixForJointDuaghters<unsigned short int>(const TimeSeriesMapT& time_series_map, nucleus& DL,  nucleus& DR, double m_k[dimsImage], double W[dimsImage * (1+dimsImage/2)]);
