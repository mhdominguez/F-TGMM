/*
 * Copyright (C) 2011-2013 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat 
 *  TGMMsupport.cpp
 *
 *  Created on: April 21st, 2014
 *      Author: Fernando Amat
 *
 * \brief class to call the classifier from TGMM code
 *        
 *
 */

#include "TGMMsupport.h"

using namespace std;

extern int mkDirectory(const std::string& folder);

//constructor
cellDivisionTemporalWindow_TGMMsupport::cellDivisionTemporalWindow_TGMMsupport(const string debugPath, config& config, int temporalWindowRadius)
	: thrCDWT(config.thrCDWT)
	, writeCDWTfeaturesToTrainClassifierFolder("")
{	
	
	cellDivisionWithTemporalWindow::setTemporalWindowRadius(temporalWindowRadius);
	
	if(setClassifierModel(config.classifierFileCDTW)) {
		// got an error code
		throw runtime_error("Could not load division classifier data from "+config.classifierFileCDTW);
	}
	cout<<"Loaded classifier "<<config.classifierFileCDTW<<" for cell division with temporal window"<<endl;

	if(config.writeCDWTfeaturesToTrainClassifier)//set folder to save training features
	{
		const auto folder_cdwt_features(debugPath+"CDTWfeatures");
		const auto error=mkDirectory(folder_cdwt_features);
		if(error>0)
			throw error;
		//create folder to write features
		//also triggeres computation of features 
		setWriteCDWTfeaturesFolder(folder_cdwt_features);
	}
};

cellDivisionTemporalWindow_TGMMsupport::cellDivisionTemporalWindow_TGMMsupport(int temporalWindowRadius, float thrCDWT_)
{
	cellDivisionWithTemporalWindow::setTemporalWindowRadius(temporalWindowRadius);
	thrCDWT = thrCDWT_;
};

std::ostream& cellDivisionTemporalWindow_TGMMsupport::print(std::ostream& os) {
	return os<<"DivisionClassifer(AmatF2013: thrCDWT="<<thrCDWT<<")";
}


//==================================================

int cellDivisionTemporalWindow_TGMMsupport::setClassifierModel(string filename)
{
		
	int errC = loadClassifier(classifierCDWT, filename);
	if(errC > 0 ) 
		return errC;

	return 0;
}

//=====================================================
int cellDivisionTemporalWindow_TGMMsupport::classifyCellDivisionTemporalWindow(lineageHyperTree& lht, int frame, vector<mylib::Array*>& imgVec, int devCUDA)
{


	//this functions follows mainCellDivisionClassifierWithTemporalWindow::mainSingleWindowForDaughters( int argc, const char** argv )
	basicEllipticalHaarFeatureVector::useDoGfeatures = true;
	cellDivisionWithTemporalWindow::setUseFixPrecisionMatrix(true);

	list<nucleus>* nucleiList = &(lht.nucleiList[frame]);

	if (nucleiList->empty() == true)
		return 0;

	//cout << "===========================DEBUGGING: NO THRESHOLD!!!!!!====================" << endl;
	if (thrCDWT <= 0.0f && writeCDWTfeatures() == false)//no point in calculating anything. All cell divisions are accepted as valid
		return 0;

	int numTrueCellDivisions = 0;


	//find all the elements that contain a division
	vector<TreeNode< ChildrenTypeLineage >* > divisionNodes;//contains a pointer to each of the divisions in a specific time point
	divisionNodes.reserve(nucleiList->size() / 10);

	for (list<nucleus>::iterator iterN = nucleiList->begin(); iterN != nucleiList->end(); ++iterN)
	{
		if (iterN->treeNodePtr->getNumChildren() == 2)
			divisionNodes.push_back(iterN->treeNodePtr);
	}

	if (divisionNodes.empty() == true)
		return 0;


	//calculate average precision matrix if necessary
	if (cellDivisionWithTemporalWindow::getUseFixPrecisionMatrix() == true)
	{
		cellDivisionWithTemporalWindow::averagePrecisionMatrixAyTM(*nucleiList);
		//increase W by a constant factor to try to incorporate the cell division
		cellDivisionWithTemporalWindow::scalePrecisionMatrixW(0.1);

		//=======================DEBUGGING===================================
		//cout << "=============================DEBUGGING: average precision matrix  time point. Frame " << frame <<"=== == == == == == == == == == == == == " << endl;
		//cout << "Wfix =[" << cellDivisionWithTemporalWindow::getPrecisionWfix(0) << "," << cellDivisionWithTemporalWindow::getPrecisionWfix(1) << "," << cellDivisionWithTemporalWindow::getPrecisionWfix(2) << "," << cellDivisionWithTemporalWindow::getPrecisionWfix(3) << "," << cellDivisionWithTemporalWindow::getPrecisionWfix(4) << "," << cellDivisionWithTemporalWindow::getPrecisionWfix(5) << "," << "]" << endl;
		//==========================================================
	}


	//extract basic Haar features for each time point
	if (imgVec.empty() == true)
		imgVec.resize(2 * cellDivisionWithTemporalWindow::getTemporalWindowRadius() + 1, NULL);//we initialize to null so the code should get the information from superovxels
	//printf("imgVec.size() is %lu, getTemporalWindowRadius is %d \n", imgVec.size(), cellDivisionWithTemporalWindow::getTemporalWindowRadius() );
	//cout << "imgVec.size() is " << imgVec.size() << ", getTemporalWindowRadius is " << cellDivisionWithTemporalWindow::getTemporalWindowRadius() << endl;
	cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchForCellDivisionSingleWindowForDaughters(divisionNodes, imgVec, cdwtVec, devCUDA, 0);


	//extend 4D Haar features using temporal combinations and classify results
	size_t numFeatures = 1000;//to pre-reserve space
	feature Fx;

	ofstream foutXYZ, foutFeat;
	streampos foutPosNumFeatures, foutPosNumSamples;
	vector<float> FxVec;
	if (writeCDWTfeatures())
	{
		cout << "Writing CDTW features to "<< writeCDWTfeaturesToTrainClassifierFolder << endl;
		char buffer[256];
		sprintf(buffer, "%s/CDWTfeatures_TM%.5d.bin", writeCDWTfeaturesToTrainClassifierFolder.c_str(), frame);
		foutFeat.open(buffer, ios::binary | ios::out);
		sprintf(buffer, "%s/CDWTfeatures_TM%.5d.txt", writeCDWTfeaturesToTrainClassifierFolder.c_str(), frame);
		foutXYZ.open(buffer);


		if (foutFeat.is_open() == false)
		{
			cout << "ERROR: at DEBUG_CDWT_FEATURES: file " << buffer << " could not be opened to save results" << endl;
			exit(3);
		}
		int aux = 0;
		foutPosNumFeatures = foutFeat.tellp();//we need it for later to insert the correct number
		foutFeat.write((char*)(&aux), sizeof(int)); //numFeatures
		//aux = nTotal * symmetryN;
		foutPosNumSamples = foutFeat.tellp();//we need it for later to insert the correct number
		foutFeat.write((char*)(&aux), sizeof(int));//numSamples

		//reserve 
		FxVec.reserve(cdwtVec.size());
	}

	for(vector<cellDivisionWithTemporalWindow>::iterator iterF = cdwtVec.begin(); iterF != cdwtVec.end(); ++iterF)
	{
		iterF->f.reserve(numFeatures);
		//iterF->calculateFeaturesSingleWindowForDaughters();
		iterF->calculateFeaturesSingleWindowForDaughters_featureSelection_v1();						
		//record for next iteration to reserve appropriate ammount of memory
		numFeatures = iterF->f.size();

		//call classifier to make the decision
		boostingTreeClassifierTranspose(&(iterF->f[0]),&Fx ,classifierCDWT , 1, numFeatures);
		iterF->cellDivisionPtr->data->probCellDivision = Fx;//store result for cell division GUI (active learning) and debugging
		if( Fx < thrCDWT )//false cell division->we need to "cut" one of the edges. //We cut linkage between mother and the fursthest daughter
		{
			lht.cutLinkageBetweenMotherAndFurthestDaughter(iterF->cellDivisionPtr);
		}else{//true cell division -> nothing to do in the lineage
			numTrueCellDivisions++;
		}

		if (writeCDWTfeatures())
		{
			iterF->writeToBinary(foutFeat);
			foutXYZ << iterF->cellDivisionPtr->data->centroid[0] << " " << iterF->cellDivisionPtr->data->centroid[1] << " " << iterF->cellDivisionPtr->data->centroid[2] << endl;
			FxVec.push_back(Fx);
		}

	}				


	if (writeCDWTfeatures())
	{
		if (FxVec.empty() == false)
		{
			//write yTrain
			foutFeat.write((char*)(&(FxVec[0])), sizeof(float)* FxVec.size());
		}

		//write the appropriate number of features and samples;
		int aux = numFeatures;
		foutFeat.seekp(foutPosNumFeatures);
		foutFeat.write((char*)(&aux), sizeof(int));
		aux = FxVec.size();
		foutFeat.seekp(foutPosNumSamples);
		foutFeat.write((char*)(&aux), sizeof(int));
		foutFeat.close();

		foutXYZ.close();
	}

	return numTrueCellDivisions;
}



