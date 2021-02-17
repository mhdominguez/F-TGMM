/*
 * parseConfigFile.cpp
 *
 *  Created on: Oct 19, 2011
 *      Author: amatf
 */

#include "parseConfigFile.h"
#include <cstring>
#include <stdlib.h> 
#include <map>
#include <fstream>

const int MAX_CHARS_PER_LINE = 1024;
const int MAX_TOKENS_PER_LINE = 2;//right now we only expect configVariable:value
const char* DELIMITER = "=";


int configOptionsTrackingGaussianMixture::parseConfigFileTrackingGaussianMixture(const string &filename)
{

	ifstream configFile(filename.c_str());
	if(!configFile.good())
	{
		cout<<"ERROR at parseConfigFileTrackingGaussianMixture: config filename "<<filename<<" could not be opened"<<endl;
		return 1;
	}

	int n;
	while (!configFile.eof())
	{
		// read an entire line into memory
		char buf[MAX_CHARS_PER_LINE];
		configFile.getline(buf, MAX_CHARS_PER_LINE);


		if(strncmp(buf,"#",1)==0) continue;//comment

		// parse the line into DELIMITER-delimited tokens

		// array to store memory addresses of the tokens in buf
		char* token[MAX_TOKENS_PER_LINE];

		// parse the line
		token[0] = strtok(buf, DELIMITER); // first token
		n=0;
		if (token[0]!=NULL) // zero if line is blank
		{
			for (n = 1; n < MAX_TOKENS_PER_LINE; n++)
			{
				token[n] = strtok(NULL, DELIMITER); // subsequent tokens
				if (token[n]==NULL) break; // no more tokens
			}
		}

		if(n!=2) continue;

		if(strcmp("imgFilePattern",token[0])==0)
		{
			imgFilePattern=string(token[1]);
		}else if(strcmp("debugPathPrefix",token[0])==0)
		{
			debugPathPrefix=string(token[1]);
			//for windows the folder separator should be \\. For example E:\\TGMMruns


		}else if(strcmp("offsetFile",token[0])==0)
		{
			offsetFile = string(token[1]);		
		}else if(strcmp("maxIterEM",token[0])==0)
		{
			maxIterEM=atoi(token[1]);
		}else if(strcmp("tolLikelihood",token[0])==0)
		{
			tolLikelihood=atof(token[1]);
		}else if(strcmp("betaPercentageOfN_k",token[0])==0)
		{
			betaPercentageOfN_k=atof(token[1]);
		}else if(strcmp("nuPercentageOfN_k",token[0])==0)
		{
			nuPercentageOfN_k=atof(token[1]);
		}else if(strcmp("alphaPercentage",token[0])==0)
		{
			alphaPercentage=atof(token[1]);
		}else if(strcmp("thrSplitScore",token[0])==0)
		{
			thrSplitScore=atof(token[1]);
		}else if(strcmp("thrCellDivisionPlaneDistance",token[0])==0)
		{
			thrCellDivisionPlaneDistance=atof(token[1]);
		}
		/*
		else if(strcmp("useIlastikBackgroundDetector",token[0])==0)
		{
			if(strcmp("false",token[1])==0)
				useIlastikBackgroundDetector=false;
			else
				useIlastikBackgroundDetector=true;
		}
		else if(strcmp("thrSignal",token[0])==0)
		{
			thrSignal=atof(token[1]);
		}*/
		else if(strcmp("estimateOpticalFlow",token[0])==0)
		{
			estimateOpticalFlow=atoi(token[1]);
		}else if(strcmp("maxDistPartitionNeigh",token[0])==0)
		{
			maxDistPartitionNeigh=atof(token[1]);
		}else if(strcmp("deathThrOpticalFlow",token[0])==0)
		{
			deathThrOpticalFlow=atoi(token[1]);
		}else if(strcmp("persistanceSegmentationTau",token[0])==0)
		{
			tau=atoi(token[1]);
		}else if(strcmp("maxNumKNNsupervoxel",token[0])==0)
		{
			KmaxNumNNsupervoxel = atoi(token[1]);
		}else if(strcmp("maxDistKNNsupervoxel",token[0])==0)
		{
			KmaxDistKNNsupervoxel=atof(token[1]);
		}else if(strcmp("temporalWindowForLogicalRules",token[0])==0)
		{
			temporalWindowRadiusForLogicalRules = atoi(token[1]);
		}else if(strcmp("anisotropyZ",token[0])==0)
		{
			anisotropyZ = atof(token[1]);
		}else if(strcmp("thrBackgroundDetectorHigh",token[0])==0)
		{
			thrBackgroundDetectorHigh = atof(token[1]);
		}else if(strcmp("thrBackgroundDetectorLow",token[0])==0)
		{
			thrBackgroundDetectorLow = atof(token[1]);
		}else if(strcmp("SLD_lengthTMthr",token[0])==0)
		{
			SLD_lengthTMthr = atoi(token[1]);
		}else if(strcmp("radiusMedianFilter",token[0])==0)
		{
			radiusMedianFilter = atoi(token[1]);
		}
		else if (strcmp("useMedianFilterForTracking", token[0]) == 0)
		{
			useMedianFilterForTracking = atoi(token[1]);
		}else if(strcmp("minTau",token[0])==0)
		{
			minTau = atof(token[1]);
		}else if(strcmp("backgroundThreshold",token[0])==0)
		{
			backgroundThreshold = atof(token[1]);
		}else if(strcmp("conn3D",token[0])==0)
		{
			conn3D = atoi(token[1]);
		}else if(strcmp("regularizePrecisionMatrixConstants_lambdaMin",token[0])==0)
		{
			lambdaMin = atof(token[1]);
		}else if(strcmp("regularizePrecisionMatrixConstants_lambdaMax",token[0])==0)
		{
			lambdaMax = atof(token[1]);
		}else if(strcmp("regularizePrecisionMatrixConstants_maxExcentricity",token[0])==0)
		{
			maxExcentricity = atof(token[1]);
		}else if(strcmp("minNucleiSize",token[0])==0)
		{
			minNucleiSize = atoi(token[1]);
		}else if(strcmp("maxNucleiSize",token[0])==0)
		{
			maxNucleiSize = atoi(token[1]);
		}else if(strcmp("maxPercentileTrimSV",token[0])==0)
		{
			maxPercentileTrimSV = atof(token[1]);
		}else if(strcmp("conn3DsvTrim",token[0])==0)
		{
			conn3DsvTrim = atoi(token[1]);
		}
		else if (strcmp("ratioHolesSv", token[0]) == 0)
		{
			ratioHolesSv = atof(token[1]);
			if (ratioHolesSv < 0.0f)
				ratioHolesSv = std::numeric_limits<float>::max();
		}
		else if (strcmp("svKeepWithoutTrimming", token[0]) == 0)
		{
			svKeepWithoutTrimming = (atoi(token[1]) != 0);
		}
		//
		// Division classifier LUT properties
		//

		else if (strcmp("cellDivisionClassifierMethod", token[0]) == 0)
		{
			try {
				map<string,DivisionClassiferKind> kinds{
					{"None",      DivisionClassiferKind_None},
	                {"AmatF2013", DivisionClassiferKind_BoostedElipticalHaarFeatures_AmatF2013},
	                {"LUT2018",   DivisionClassiferKind_LookupTable_2018}
				};
				cellDivisionClassifier.method=kinds.at(token[1]);
			} catch(...) {
				cout <<"ERROR: cellDivisionClassifierMethod must be None, AmatF2013 or LUT2018." 
				     << "\tGot: " << token[1] << endl;
				return 1;
			}
		}
		// Amatf2013
		else if(strcmp("thrCellDivisionWithTemporalWindow",token[0])==0)
		{
			cellDivisionClassifier.Amatf2013.thrCDWT = atof(token[1]);
		}
		else if (strcmp("writeCDWTfeaturesToTrainClassifier", token[0]) == 0)
		{
			cellDivisionClassifier.Amatf2013.writeCDWTfeaturesToTrainClassifier = atoi(token[1]);
		}
		else if (strcmp("classifierFileCDTW", token[0]) == 0)
		{
			cellDivisionClassifier.Amatf2013.classifierFileCDTW = string(token[1]);
		}

		// LUT2018
		else if (strcmp("cellDivisionLUT_filename", token[0]) == 0)
		{ 
			cellDivisionClassifier.LUT.filename=string(token[1]);
		}
		else if (strcmp("cellDivisionLUT_spatial_threshold_px", token[0]) == 0)
		{
			cellDivisionClassifier.LUT.spatial_threshold_px=atof(token[1]);
		}
		else if (strcmp("cellDivisionLUT_temporal_threshold_frames", token[0]) == 0)
		{
			cellDivisionClassifier.LUT.temporal_threshold_frames=atoi(token[1]);
		}
		// fall through
		else{
			cout<<"WARNING at parseConfigFileTrackingGaussianMixture: does not recognize config option "<<token[0]<<endl;
		}				
	}
	configFile.close();

	//
	// Fill in any dependant properties
	//
	cellDivisionClassifier.LUT.anisotropyZ=anisotropyZ;

	//
	// Validation
	// check that the mandatory values are not set to default
	//

	if(strcmp(imgFilePattern.c_str(),"empty")==0)
	{
		cout<<"ERROR at parseConfigFileTrackingGaussianMixture: mandatory variable imgFilePattern was not set at config filename "<<filename<<" could not be opened"<<endl;
		return 1;
	}

	if(strcmp(debugPathPrefix.c_str(),"empty")==0)
	{
		cout<<"ERROR at parseConfigFileTrackingGaussianMixture: mandatory variable debugPathPrefix was not set at config filename "<<filename<<" could not be opened"<<endl;
		return 1;
	}

	if( backgroundThreshold < -1e31 )
	{
		cout<<"ERROR at parseConfigFileTrackingGaussianMixture: mandatory variable backgroundThreshold was not set at config filename "<<filename<<" could not be opened"<<endl;
		return 1;
	}

	return cellDivisionClassifier.validate();
}

static string to_string(const DivisionClassiferKind kind) {
	switch(kind) {
		case DivisionClassiferKind_None: return "None";
		case DivisionClassiferKind_BoostedElipticalHaarFeatures_AmatF2013: return "AmatF2013";
		case DivisionClassiferKind_LookupTable_2018: return "LUT2018";
		default:
			throw runtime_error("Error: Require an string representation of the division classifier method.");
	}
}

void configOptionsTrackingGaussianMixture::printConfigFileTrackingGaussianMixture(ostream &outLog) const
{
	outLog<<"imgFilePattern="<<imgFilePattern<<endl;
	outLog<<"maxIterEM="<<maxIterEM<<endl;
	outLog<<"tolLikelihood="<<tolLikelihood<<endl;
	outLog<<"betaPercentageOfN_k="<<betaPercentageOfN_k<<endl;
	outLog<<"nuPercentageOfN_k="<<nuPercentageOfN_k<<endl;
	outLog<<"alphaPercentage = "<<alphaPercentage<<endl;
	outLog<<"thrSplitScore="<<thrSplitScore<<endl;
	outLog<<"thrCellDivisionPlaneDistance="<<thrCellDivisionPlaneDistance<<endl;	
	//outLog<<"thrSignal="<<thrSignal<<endl;
	//outLog<<"Use Ilastik background/foreground detection="<<useIlastikBackgroundDetector<<endl;
	outLog<<"estimateOpticalFlow="<< estimateOpticalFlow<<endl;
	outLog<<"deathThrOpticalFlow="<< deathThrOpticalFlow<<endl;
	outLog<<"maxDistPartitionNeigh="<< maxDistPartitionNeigh<<endl;
	outLog<<"persistanceSegmentationTau="<< tau<<endl;
	outLog<<"maxNumKNNsupervoxel="<<KmaxNumNNsupervoxel<<endl;
	outLog<<"maxDistKNNsupervoxel="<<KmaxDistKNNsupervoxel<<endl;
	outLog<<"temporalWindowForLogicalRules (radius)="<<temporalWindowRadiusForLogicalRules<<endl;
	outLog<<"anisotropyZ="<<anisotropyZ<<endl;
	outLog<<"offset file = "<<offsetFile<<endl;
	

	outLog<<"thrBackgroundDetectorHigh="<<thrBackgroundDetectorHigh<<endl;
	outLog<<"thrBackgroundDetectorLow="<<thrBackgroundDetectorLow<<endl;
	outLog<<"short-lived duaghter thre = "<<SLD_lengthTMthr<<endl;

	outLog<<"cellDivisionClassifierMethod"<<to_string(cellDivisionClassifier.method);
	
	outLog << "classifier file CDTW = " <<cellDivisionClassifier.Amatf2013.classifierFileCDTW << endl;
	outLog<<"thrCellDivisionWithTemporalWindow="<<cellDivisionClassifier.Amatf2013.thrCDWT<<endl;
	outLog<<"writeCDWTfeaturesToTrainClassifier = " << cellDivisionClassifier.Amatf2013.writeCDWTfeaturesToTrainClassifier << endl;
	
	outLog<<"cellDivisionLUT_filename"<<cellDivisionClassifier.LUT.filename<<endl;
	outLog<<"cellDivisionLUT_spatial_threshold_px"<<cellDivisionClassifier.LUT.spatial_threshold_px<<endl;
	outLog<<"cellDivisionLUT_temporal_threshold_frames"<<cellDivisionClassifier.LUT.temporal_threshold_frames<<endl;

	outLog<<"radiusMedianFilter ="<<radiusMedianFilter<<endl;
	outLog << "useMedianFilterForTracking = " << useMedianFilterForTracking << endl;
	outLog<<"minTau ="<<minTau<<endl;
	outLog<<"backgroundThreshold ="<<backgroundThreshold<<endl;
	outLog<<"conn3D ="<<conn3D<<endl;

	outLog<<"regularizePrecisionMatrixConstants_lambdaMin = "<<lambdaMin<<endl;
	outLog<<"regularizePrecisionMatrixConstants_lambdaMax = "<<lambdaMax<<endl;
	outLog<<"regularizePrecisionMatrixConstants_maxExcentricity = "<<maxExcentricity<<endl;


	outLog<<"Trim supervoxel minNucleiSize = "<<minNucleiSize<<endl;
	outLog<<"Trim supervoxel maxNucleiSize ="<<maxNucleiSize<<endl;
	outLog<<"Trim supervoxel maxPercentileTrimSV = "<<maxPercentileTrimSV<<endl;
	outLog<<"Trim supervoxel conn3DsvTrim = "<<conn3DsvTrim<<endl;
	outLog << "Ratio holes vs number of pixels per supervoxel = " << ratioHolesSv << endl;
	outLog << "output full supervoxels without trimming = " << svKeepWithoutTrimming << endl;

}

/// \returns 0 when valid, 1 otherwise
int configOptionsTrackingGaussianMixture::CellDivisionClassifier::validate() const  {
	switch(method) {
		case DivisionClassiferKind_BoostedElipticalHaarFeatures_AmatF2013: {
			if(Amatf2013.classifierFileCDTW.empty()) {
				cout<<"ERROR: Configuration failed to validate.\n\tWhen using teh AmatF2013 division classifer, the classifierFileCDTW must be set."<<endl;
				return 1;
			}
		} break;
		case DivisionClassiferKind_LookupTable_2018: {
			if(LUT.filename.empty()) {
				cout<<"ERROR: Configuration failed to validate.\n\tWhen using teh LUT2018 division classifer, the cellDivisionLUT_filename must be set."<<endl;
				return 1;
			}
		} break;		
		default: //DivisionClassiferKind_None
			return 0; 
	}
	// anything that falls through is ok
	return 0;
}