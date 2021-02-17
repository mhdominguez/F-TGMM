/*
* mainTrainSet.cxx
*
*  Created on: May 5, 2014
*      Author: amatf
*
* \brief generate precision recall curve for a given training set and a model classifier
*/


#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "gentleBoost.h"

using namespace std;

//=============================================================================

int main( int argc, const char** argv )
{

	if(argc != 3)
	{
		std::cout<<"ERROR: input arguments are <model filename> <training binary filename>"<<endl;
		return 2;
	}

	string classifierFilename(argv[1]);
	string filename(argv[2]);
	
	

	//variables ot store training data
	long long int numFeatures = -1, numSamples = 0;		
	vector<float> xTrain, yTrain;
	//read training data
	readBinaryTrainingFile(filename, xTrain, yTrain, numFeatures, numSamples);

	//load classifier
	vector<vector< treeStump > > classifier;
	int err = loadClassifier(classifier, classifierFilename);
	if( err > 0 )
		cout<<"ERROR: loading classifier "<<classifierFilename<<endl;


	//run classifier
	float *Fx = new float[ numSamples ];
	
	transposeXtrainOutOfPlace(&(xTrain[0]), numSamples, numFeatures);//we need to perform transposition from GPU features to gentleBoost classifier				
	boostingTreeClassifier(&(xTrain[0]), Fx ,classifier , numSamples, numFeatures);

	//precision recall curve
	string ROCfilename( filename + "_ROCtraining.txt");
	ofstream out( ROCfilename.c_str() );
	if( !out.is_open())
	{
		cout<<"ERROR: opening file "<<ROCfilename<<" to save ROC curve"<<endl;
	}else{
		cout<<"Saving precision,recall at "<<ROCfilename<<endl;
	}
	precisionRecallAccuracyCurve(&(yTrain[0]), Fx, numSamples, out, 0.01);
	out.close();
	
	string FXfilename(filename + "_FXtraining.txt");
	out.open(FXfilename.c_str());
	if (!out.is_open())
	{
		cout << "ERROR: opening file " << FXfilename << " to save ROC curve" << endl;
	}
	else{
		cout << "Saving classifier scores at " << FXfilename << endl;
	}
	for (long long int ii = 0; ii < numSamples; ii++)
		out << Fx[ii] << endl;	
	out.close();

	//release memory
	delete[] Fx;

	return 0;
}
