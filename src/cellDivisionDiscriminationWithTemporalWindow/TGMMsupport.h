/*
 * Copyright (C) 2011-2013 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat 
 *  TGMMsupport.h
 *
 *  Created on: April 21st, 2014
 *      Author: Fernando Amat
 *
 * \brief class to call the classifier from TGMM code
 *        
 *
 */

#ifndef __CELL_DIVISION_WITH_TEMPORAL_WINDOW_TGMM_SUPPORT_H__
#define __CELL_DIVISION_WITH_TEMPORAL_WINDOW_TGMM_SUPPORT_H__

#include <fstream>
#include <string>
#include <vector>
#include "cellDivisionWithTemporalWindow.h"
#include "lineageHyperTree.h"
#include "gentleBoost/gentleBoost.h"
#include "DivisionClassifierFactory.h"
#include "Utils/parseConfigFile.h"

class cellDivisionTemporalWindow_TGMMsupport : public IDivisionClassifier
{
public:
	using config=struct configOptionsTrackingGaussianMixture::CellDivisionClassifier::Amatf2013;

	cellDivisionTemporalWindow_TGMMsupport(const string debugPath, config &config, int temporalWindowRadius);
	cellDivisionTemporalWindow_TGMMsupport(int temporalWindowRadius, float thrCDWT_);


	//main function to call from TGMM
	//this function will disconnect edges of the elements that are below thrCDWT
	//returns -1 if there is an error. Otherwise returns the number of true cell divisions
	//int classifyCellDivisionTemporalWindow(lineageHyperTree& lht, int frame, std::vector<mylib::Array*>& imgVec, int devCUDA, double thrCellDivisionPlaneDistance) final;
	int classifyCellDivisionTemporalWindow(lineageHyperTree& lht, int frame, std::vector<mylib::Array*>& imgVec, int devCUDA, double thrCellDivisionPlaneDistance, float *im_zero, float *im_plus_one, bool regularize_W4DOF, float scaleOrig[3]) final;
	float cellDivisionPlaneDistance2021(float centroidMM[dimsImage], float centroidDL[dimsImage], float centroidDR[dimsImage], const float scale[dimsImage], float &sqrtnorm );
	int classifyCellDivisionTemporalWindow(lineageHyperTree& lht, int frame, std::vector<mylib::Array*>& imgVec, int devCUDA ) final;

	//I/O functions
	int setClassifierModel(string filename);
	void setThrCDWT(float p) { thrCDWT = p;};
	float getThrCDWT() { return thrCDWT;};
	size_t getNumCellDivisions() { return cdwtVec.size();};

	void setWriteCDWTfeaturesFolder(string p) { writeCDWTfeaturesToTrainClassifierFolder = p; } ;
	bool writeCDWTfeatures() const { return (!writeCDWTfeaturesToTrainClassifierFolder.empty()); };	
	string getWriteCDWTfeaturesFolder() const { return writeCDWTfeaturesToTrainClassifierFolder; };

	std::ostream& print(std::ostream& os);
private:
	vector<cellDivisionWithTemporalWindow> cdwtVec;//we need to handle multiple cell divisions at once
	vector< vector< treeStump> > classifierCDWT;
	float thrCDWT;

	std::string writeCDWTfeaturesToTrainClassifierFolder;//variable to set with config file in case we need to output features to train classifier. Leave empty to not write features

};


#endif //__CELL_DIVISION_WITH_TEMPORAL_WINDOW_TGMM_SUPPORT_H__
