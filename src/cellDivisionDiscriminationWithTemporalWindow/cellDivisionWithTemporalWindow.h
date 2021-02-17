/*
 * Copyright (C) 2011-2013 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat 
 *  cellDivisionWithTemporalWindow.h
 *
 *  Created on: April 21st, 2014
 *      Author: Fernando Amat
 *
 * \brief temporal features (image-based and geometrical and spatial) to discriminate between true and false cell divisions detected by the TGMM pipeline
 *        
 *
 */

#ifndef __CELL_DIVISION_WITH_TEMPORAL_WINDOW_H__
#define __CELL_DIVISION_WITH_TEMPORAL_WINDOW_H__

#include <fstream>
#include <string>
#include "lineageHyperTree.h"
#include "EllipticalHaarFeatures.h"


//#define CDWT_SAVE_FEATURE_NAME //uncomment this if want to write the name of each feature so we can understand which features where selected by boosting

namespace mylib
{
	#include "mylib/image.h"
}

class cellDivisionWithTemporalWindow
{
public:
	vector<float> f;//vector of final features
	//pointer to a vector with the basic Haar features for each time point so further features can be calculated. 
	//The size of teh vector should be equal to temporalWindowRadius + 1 + 2 * temporalWindowRadius (we have cell division) 
	//HOW TO ACCES FEATURES: [par(-t), par(-t+1),...,mother = [par[0]], daughter_left_(+1), daughter_left_(+2), ..., daughter_left_(+t), daughter_right_(+1), daughter_right_(+2), ..., daughter_right_(+t)] 
	vector<basicEllipticalHaarFeatureVector> fHaarVec;
	
	TreeNode<ChildrenTypeLineage >* cellDivisionPtr;//pointer to teh cell division element, so we can compute shape, geometrical, etc features



	//constructor / destructor
	cellDivisionWithTemporalWindow();
	cellDivisionWithTemporalWindow(int temporalWindowRadius_);
	~cellDivisionWithTemporalWindow();

	//set / get function
	static long long int getNumFeatures();
	static void setTemporalWindowRadius(int p){ temporalWindowRadius = p;};
	static int getTemporalWindowRadius(void){ return temporalWindowRadius;};
	static int getTotalTemporalWindowSize(void){ return 1 + 3 * temporalWindowRadius;};//we have to accomodate two daughters
	static int getTotalTemporalWindowSizeSingleWindowForDaughters(void){ return 1 + 2 * temporalWindowRadius;};//we have to accomodate two daughters
	static void setUseFixPrecisionMatrix(bool p){ useFixPrecisionMatrix = p;};
	static bool getUseFixPrecisionMatrix(void){ return useFixPrecisionMatrix;};
	static int getSizePrecisionW(void){return dimsImage*(1+dimsImage)/2;};
	static void setTemporalUndersampling(int p){ temporalUndersampling = p;};
	static int gtTemporalUndersampling(void){ return temporalUndersampling;};
	static double getPrecisionWfix(int pos){ return precisionWfix[pos]; };

#ifdef CDWT_SAVE_FEATURE_NAME
	void setFeaturePrefix(string& s, int p);
	void addFeatureName(string& s, int p, int ll);
#endif
	//main functions to calculate features

	/*
	\brief Given a set of of division nodes to check whther they represent true divisions or not, we calculate the basic haar features for each time point. We do it in batch fasion to optimize performance by minimizing img transfers to GPU

	*/
	static int calculateBasicEllipticalHaarFeaturesBatchForCellDivision(const vector< TreeNode< ChildrenTypeLineage >* >& divisionNodes, const vector<mylib::Array*>& imgVec, vector<cellDivisionWithTemporalWindow>& cdtwVec, int devCUDA, int symmetry);
	

	/*
	\brief The same as calculateBasicEllipticalHaarFeaturesBatchForCellDivision but using a single box for each pair of daughters
			

	*/
	static int calculateBasicEllipticalHaarFeaturesBatchForCellDivisionSingleWindowForDaughters(const vector< TreeNode< ChildrenTypeLineage >* >& divisionNodes, const vector<mylib::Array*>& imgVec, vector<cellDivisionWithTemporalWindow>& cdtwVec, int devCUDA, int symmetry);

	/*
	\brief same calculateBasicEllipticalHaarFeaturesBatchForCellDivision but for a scpefici time point
	*/
	static int calculateBasicEllipticalHaarFeaturesBatchAtTM(const vector< TreeNode< ChildrenTypeLineage >* >& auxNodeVec, mylib::Array* img, vector<cellDivisionWithTemporalWindow>& cdtwVec, int relativeTimeWithinWindow,int devCUDA, int symmetry, double **m, double **W, long long int **dimsVec);

	/*
	\brief We assume m and W are already placed appropiately 
	*/
	static int calculateBasicEllipticalHaarFeaturesBatchAtTM(const vector< bool >& isDead, mylib::Array* img, vector<cellDivisionWithTemporalWindow>& cdtwVec, int relativeTimeWithinWindow,int devCUDA, int symmetry, double *m, double *W, long long int **dimsVec);


	/*
	\brief Main function that computes all the different features (image-based, shape, trajectory, etc...). For image-based fHarrVec needs to be calculated beforehand

	\return			Greater than zero if features could not be calculated
	*/
	int calculateFeatures();

	int calculateFeaturesSingleWindowForDaughters();//assumes that number of temporal elements is 2*temporalWindowRadius+1
	

	int calculateFeaturesSingleWindowForDaughters_featureSelection_v1();//after training with boosting, it only calculates the features used in the tree (~300). It is a subset of all the ones calculated at calculateFeaturesSingleWindowForDaughters


	//support functions
	void  writeToBinary(std::ofstream &fout);
	static void averagePrecisionMatrixAyTM(const list<nucleus>& listN);//listN = lht.listNuclei[TM]
	static void scalePrecisionMatrixW(double s);
#ifdef CDWT_SAVE_FEATURE_NAME
	void  writeFeaturesNames(string& filenameS);
#endif
	

protected:

	/*
	\brief given a feature time-series (basically values of a feature over time) it calculates different features (mean, max, harmonic) to try to extract information
	Inspired by M. Kabra, A. A. Robie, M. Rivera-Alba, S. Branson, and K. Branson, “JAABA: interactive machine learning for automatic annotation of animal behavior,” Nature Methods, vol. 10, no. 1, pp. 64–67, Dec. 2012.
	
	@param[in] timeSeriesF: vector containing the value of the time series	
	
	@returns	number of elements written in this->f
	*/
	int temporalFeatureExtraction(vector<float>& timeSeries, bool isRecursiveCall = false, bool calculateMultiScale = true);


	int temporalFeatureExtractionSingleWindowForDaughters(vector<float>& timeSeries, bool isRecursiveCall = false, bool calculateMultiScale = true);

	/*
	\brief extract a vector with all the pointers in the lineage within the temporal window so we can calculate features efficiently	
	*/
	void lineagePointersWithinTemporalWindowWithCellDivision(vector< TreeNode<ChildrenTypeLineage>* >&  lineageVec);

	/*
	\brief Functions to calculate different kinds of features. They should be called from the "mother function" calculateFeatures.
		They all follow the same pattern in terms of input/ output parameters.

	
	@returns	number of elements written in this->f
	*/
	int calculateHaarBasedFeatures();
	int calculateHaarBasedFeaturesSingleWindowForDaughters();

	int calculateLineageBasedFeatures();
	



	static TreeNode<ChildrenTypeLineage>* moveForwardInlineageEvenWithCellDivision(TreeNode<ChildrenTypeLineage>* auxNode);
	template<class imgTypeC>
	static void calculatePrecisionMatrixForJointDuaghters(const TimeSeriesMapT& time_series_map,nucleus& DL,  nucleus& DR, double m_k[dimsImage], double W[dimsImage * (1+dimsImage/2)]);
	static void calculatePrecisionMatrixForJointDuaghters(nucleus& DL, nucleus& DR, double W[dimsImage * (1+dimsImage/2)]);//without using the image intensity information as weights

private:

	static int temporalWindowRadius;//number of time points includes to calculate features
	static bool useFixPrecisionMatrix;//true->we use the same precision W for all time points within the temporal window. So Haar features capture changes in size (but we lose excentricity calculations)
	static double precisionWfix[dimsImage * (1+dimsImage) /2];//if useFixPrecisionMatrix == true, then we use the value in here to calculate 3D Elliptical Haar features
	static int temporalUndersampling;//in order to subsample in time to obtain different training samples

	static short int HaarPatternTriplet[39];
	static int HaarPatternTripletN;//number of HaarPatternTriplets(39/3);

	#ifdef CDWT_SAVE_FEATURE_NAME
	vector<string> fName;//save feature name
	string fNamePrefix;//to update before calling temporal window features
#endif
};


//===========================================================================================
inline void  cellDivisionWithTemporalWindow::writeToBinary(std::ofstream &fout)
{
	fout.write((char*)(&(f[0])), sizeof(float) * f.size() );
}

//======================================================================

inline long long int cellDivisionWithTemporalWindow::getNumFeatures()
{
        return 0;//make sure this matches with backgroundDetectionFeatures::calculateFeatures (there is an assert)
}

//===================================================================
inline void cellDivisionWithTemporalWindow::scalePrecisionMatrixW(double s)
{
	for(int ww = 0; ww < cellDivisionWithTemporalWindow::getSizePrecisionW(); ww++ )
	{
		precisionWfix[ww] *= s;
	}
}

#endif //__CELL_DIVISION_WITH_TEMPORAL_WINDOW_H__
