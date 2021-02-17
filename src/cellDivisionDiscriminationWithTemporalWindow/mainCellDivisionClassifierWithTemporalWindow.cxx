#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <string.h>
#include <cstring>
#include "EllipticalHaarFeatures.h"
#include "external/Nathan/tictoc.h"
#include "external/xmlParser/xmlParser.h"
#include "gentleBoost/gentleBoost.h"
#include "AnnotationEllipsoid.h"
#include "temporalLogicalRules.h"
#include "cellDivisionWithTemporalWindow.h"

#ifdef PICTOOLS_JP2K
	#include "ioFunctions.h"
#endif
#include "klb_imageIO.h"

namespace mylib
{
	#include "mylib/image.h"
}

#if defined(_WIN32) || defined(_WIN64)

#include <Windows.h>

#endif


using namespace std;

static void eraseTM(TimeSeriesMapT &ts, int TM) {
    mylib::Free_Array(ts.at(TM)); //free because it was allocated using mylib
    ts.erase(TM);
}

int mainIndependentWindowsPerDaughter( int argc, const char** argv );//this function calculates 3 * temporalWindowRadius +1 boxes (each daughter has a unique box centered around her)
int mainSingleWindowForDaughters( int argc, const char** argv );//this function calculates 2 * temporalWindowRadius + 1 boxes (we try to capture daughters using a single window)
int debug_mainSingleWindowForDaughters_writeImageBoxes( int argc, const char** argv );


int main( int argc, const char** argv )
{
	//return mainIndependentWindowsPerDaughter( argc, argv );//this function calculates 3 * temporalWindowRadius +1 boxes (each daughter has a unique box centered around her)
	return mainSingleWindowForDaughters( argc, argv );//this function calculates 2 * temporalWindowRadius + 1 boxes (we try to capture daughters using a single window)

	//--------------for debugging purposos
	//return debug_mainSingleWindowForDaughters_writeImageBoxes( argc, argv );
}


//========================================================================
#if defined(_WIN32) || defined(_WIN64)
//from http://stackoverflow.com/questions/11962554/how-to-get-all-file-names-in-current-directory
vector<string> get_all_files_within_folder(const string& pattern)
{
    vector<string> names;
    char search_path[256];
    sprintf(search_path, "%s", pattern.c_str());
    WIN32_FIND_DATA fd; 
    HANDLE hFind = FindFirstFile(search_path, &fd); 
    if(hFind != INVALID_HANDLE_VALUE) 
    { 
        do 
        { 
            // read all (real) files in current folder, delete '!' read other 2 default folder . and ..
            if(! (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) ) 
            {
                names.push_back(fd.cFileName);
            }
        }while(FindNextFile(hFind, &fd)); 
        FindClose(hFind); 
    } 
    return names;
}
#else
vector<string> get_all_files_within_folder(const string& pattern)
{
	vector<string> names;
	cout<<"ERROR: function get_all_files_within_folder not implemented for non-Windows operating systems"<<endl;
	return names;
}
#endif


//=======================================================================
void parseImagePath(string& imgRawPath, int frame)
{	

	size_t found=imgRawPath.find_first_of("?");
	while(found != string::npos)
	{
		int intPrecision = 0;
		while ((imgRawPath[found] == '?') && found != string::npos)
		{
			intPrecision++;
			found++;
			if( found >= imgRawPath.size() )
				break;

		}

		
		char bufferTM[16];
		switch( intPrecision )
		{
		case 2:
			sprintf(bufferTM,"%.2d",frame);
			break;
		case 3:
			sprintf(bufferTM,"%.3d",frame);
			break;
		case 4:
			sprintf(bufferTM,"%.4d",frame);
			break;
		case 5:
			sprintf(bufferTM,"%.5d",frame);
			break;
		case 6:
			sprintf(bufferTM,"%.6d",frame);
			break;
		}
		string itoaTM(bufferTM);
		
		found=imgRawPath.find_first_of("?");
		imgRawPath.replace(found, intPrecision,itoaTM);
		

		//find next ???
		found=imgRawPath.find_first_of("?");
	}
	
}

//========================================================
mylib::Array* readImage(string imgFilePattern, int frame)
{
	string imgBasename( imgFilePattern );
	parseImagePath(imgBasename, frame);

	string imgFilename;	
	void* data = NULL;
	mylib::Value_Type type;
	mylib::Dimn_Type *dims = NULL;
	int ndims;
#ifdef PICTOOLS_JP2K
	imgFilename = string( imgBasename + ".jp2");	
		
	ifstream fin(imgFilename.c_str());
	if( fin.is_open() == true )//jp2 file exists
	{
		fin.close();
		data = readJP2Kfile(imgFilename, type, ndims, &dims);				
	}
#endif		


	//try to read KLB format
	if (data == NULL)
	{
		imgFilename = string((imgBasename + ".klb").c_str());
		//check if file exists
		FILE* fcheck = fopen(imgFilename.c_str(), "rb");
		if (fcheck != NULL)
		{
			fclose(fcheck);
			klb_imageIO imgKLB(imgFilename);
			int err = imgKLB.readHeader();
			if (err > 0)
				return NULL;

			data = malloc(imgKLB.header.getImageSizeBytes());

			err = imgKLB.readImageFull((char*)data, 0);
			if (err > 0)
				return NULL;
			//update parameters for mylib
			switch (imgKLB.header.dataType)//we follow the same convention as mylib
			{
			case KLB_DATA_TYPE::UINT8_TYPE:
				type = mylib::UINT8_TYPE;
				break;
			case KLB_DATA_TYPE::UINT16_TYPE:
				type = mylib::UINT16_TYPE;
				break;
			case KLB_DATA_TYPE::FLOAT32_TYPE:
				type = mylib::FLOAT32_TYPE;
				break;
			default:
				cout << "ERROR: file type stored in KLB file not supported" << endl;
				free(data);
				return NULL;
			}
			ndims = 0;
			if (dims == NULL)
				dims = new mylib::Dimn_Type[KLB_DATA_DIMS];
			while (ndims < KLB_DATA_DIMS && imgKLB.header.xyzct[ndims] > 1)
			{
				dims[ndims] = imgKLB.header.xyzct[ndims];
				ndims++;
			}
		}
	}


	mylib::Array* img = NULL;
	if(data == NULL)//try top read tif image
	{
		imgFilename = string( imgBasename + ".tif" );
		img = mylib::Read_Image((char*)(imgFilename.c_str()),0);
		if( img == NULL)
		{
			cout<<"ERROR: reading image "<<imgFilename<<endl;
			exit(3);
		}
	}else{//wrap data into mylib::array
		img = mylib::Make_Array_Of_Data(mylib::PLAIN_KIND, type, ndims, dims, data);
	}

	return img;
}

//=============================================================================


//=============================================================================
int mainIndependentWindowsPerDaughter( int argc, const char** argv )
{
	cout<<"Generating samples for cell division discrimination using temporal window of features"<<endl;	
	cout<<"Method using an individual box per daughter"<<endl;	
	

	//default values: small test set for debugging
	/*
	int symmetryN = 1;//calculate features for this image using all possible 8 combinations of eigenvectors directions to artifically increment the traning data	
	string imgFilePattern("C:/Users/Fernando/cppProjects/TrackingGaussianMixtures/NM2013-paperRun/data/dataJP2/TM?????_timeFused_blending/SPC0_CM0_CM1_CHN00_CHN01.fusedStack_?????");
	string TGMMoutputXMLfolder("E:/TGMMruns/GMEMtracking3D_2014_5_2_13_51_14_datasetTest_debugCellDivTempWindowDiscrimination");
	int temporalWindowRadiusForLogicalRules = 5;
	*/
	
	//default values: large sample
	
	int symmetryN = 1;//calculate features for this image using all possible 8 combinations of eigenvectors directions to artifically increment the traning data	
	string imgFilePattern("X:/SiMView2/12-08-28/Dme_E1_His2ARFP_01_20120828_144957.corrected/Results/TimeFused.Blending/Dme_E1_His2ARFP.TM????_timeFused_blending/SPC0_TM????_CM0_CM1_CHN00_CHN01.fusedStack");
	string TGMMoutputXMLfolder("E:/TGMMruns/GMEMtracking3D_2014_4_29_2_27_1_dataset12_08_28_drosophila_simview_temporalLRdeactivatedForCellDivisionTraining");
	int temporalWindowRadiusForLogicalRules = 5;
	

	//fixed parameters
	string positiveClassLabel("cellDivisionCorrect");//the class in the annotations that is considered a positive sample
	int devCUDA = 0;
	basicEllipticalHaarFeatureVector::useDoGfeatures = true;
	cellDivisionWithTemporalWindow::setUseFixPrecisionMatrix(true);

	if(argc == 5)
	{	
		TGMMoutputXMLfolder = string(argv[1]);
		imgFilePattern = string(argv[2]);
		temporalWindowRadiusForLogicalRules = atoi(argv[3]);
		symmetryN = atoi( argv[4] );//calculate features for this image using all possible 8 combinations of eigenvectors directions to artifically increment the traning data	

	}else if(argc == 1){
		cout<<"Using default parameters for debugging purposes"<<endl;
	}else{
		std::cout<<"ERROR: input arguments are <numWeakClassifiers> <J> <symmetryN> <crossValidationPercentage>"<<endl;
		return 2;
	}

	
	
	//basicEllipticalHaarFeatureVector::useDoGfeatures = false;//uncomment this if you want to test results without DoG features
	if( basicEllipticalHaarFeatureVector::useDoGfeatures == false)
		std::cout<<"======================WARNING: TRAINING WITHOUT DoG FEATURES ENHANCEMENT========================="<<endl;
	else
		std::cout<<"======================WARNING: TRAINING WITH DoG FEATURES ENHANCEMENT========================="<<endl;
	if( cellDivisionWithTemporalWindow::getUseFixPrecisionMatrix() == true)
		std::cout<<"======================WARNING: TRAINING WITH FIX PRECISION MATRIX W ========================="<<endl;
	else
		std::cout<<"======================WARNING: TRAINING WITHOUT FIX PRECISION MATRIX W ========================="<<endl;
	//=====================================================================================

	TicTocTimer tt = tic();

	symmetryN = max(symmetryN, 1);//minimum value
	symmetryN = min(symmetryN, 8);//maximum value
	cellDivisionWithTemporalWindow::setTemporalWindowRadius( temporalWindowRadiusForLogicalRules );

	//------------------------------------------------------------------------------
	//store positions of ??? to read time point from annotations
	size_t found = imgFilePattern.find_first_of("?");
	size_t qIniPos = found;
	while ((imgFilePattern[found] == '?') && found != string::npos)
	{
		found++;
		if( found >= imgFilePattern.size() )
			break;
	}
	size_t qEndPos = found;


	//----------------read annotations-----------------------------------------
	string annotationFileFolder(TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin");
	
	vector<string> annotationFiles = get_all_files_within_folder(annotationFileFolder + "/*.xml");


	vector< AnnotationEllipsoid> annotationsVec;
	long long int nPos = 0, nNeg = 0, nTotal = 0;
	for(size_t jj = 0; jj < annotationFiles.size(); jj++)
	{
		string auxFile(TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin/" + annotationFiles[jj]);
		cout<<"Reading annotation file "<<auxFile<<endl;
		XMLNode xMainNode = XMLNode::openFileHelper(auxFile.c_str(),"document");
		int n = xMainNode.nChildNode("Surface");

		annotationsVec.resize( nTotal + n );				
		for(int ii=0;ii<n;ii++)
		{
			annotationsVec[nTotal + ii] = AnnotationEllipsoid(xMainNode,ii);
			//parse class name to labels
			if( annotationsVec[nTotal + ii].className.compare(positiveClassLabel) == 0 ) //cell division
			{
				annotationsVec[nTotal + ii].classVal = 1;
				nPos++;
			}else{
				annotationsVec[nTotal + ii].classVal = 0;
				nNeg++;
			}
			//find out TM from imgFilename and compare it to imgFilenamePattern
			annotationsVec[nTotal + ii].TM = std::stoi(annotationsVec[nTotal + ii].imgFilename.substr(qIniPos, qEndPos-qIniPos) );
		}

		nTotal = annotationsVec.size();
	}
	
	//sort by filename (which is equivalent to sorting by time point)
	sort(annotationsVec.begin(), annotationsVec.end());
	int iniFrame = std::max(annotationsVec.front().TM - temporalWindowRadiusForLogicalRules, 0);
	int endFrame = annotationsVec.back().TM + temporalWindowRadiusForLogicalRules;

	vector<float> yTrain;//to store annotations
	yTrain.reserve( nTotal * symmetryN );
	cout<<"Read "<<nTotal<<" annotations in "<<toc(&tt)<<" secs. Positives = "<<nPos<<"; negatives="<<nNeg<<"; symmetry = "<<symmetryN<<endl;


	//-----------------main loop: read TGMM solution in a sliding window fashion and calculate features for annotations---------------------------------
	string basenameXML(TGMMoutputXMLfolder + "/XML_finalResult_lht/GMEMfinalResult_frame");
	lineageHyperTree lht(endFrame+temporalWindowRadiusForLogicalRules + 2);
	//read initial set of frames: function generates all lht structure (sv, nuclei and lineages)
	int err = parseGMMtrackingFilesToHyperTree(imgFilePattern, basenameXML,iniFrame, iniFrame + temporalWindowRadiusForLogicalRules * 2, lht);	
	if( err > 0 )
		return err;

	//read images
    TimeSeriesMapT time_series_map;  // (ngc) it's awkward to have this and imgVec.  Ideally I'd just choose one, but I'm not sure what the right thing is at the moment.
	vector< mylib::Array* > imgVec(endFrame + temporalWindowRadiusForLogicalRules + 2, NULL);
	for(int frame = iniFrame; frame <= iniFrame + temporalWindowRadiusForLogicalRules * 2; frame++)
	{
        auto im=readImage(imgFilePattern,frame);
		imgVec[frame] = im;// (ngc) shouldn't this be imgVec[frame-iniframe]?
        time_series_map[frame]=im;

		if( imgVec[frame]->type != mylib::UINT16_TYPE )
		{
			cout<<"ERROR: at mainCellDivisionWithTemporalWindow: code is not ready for non UINT16 images"<<endl;
			exit(3);
		}
	}


	//prepare output binary file 
	//header with info (we'll rewrite number of samples at the end)	
	//open file for annotations
	cout<<"=========TODO: add specifics to the name of the file (DoG, numFeatures, temporalWIndowSize, etc)"<<endl;
	string filenameOut (TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin/" + "trainFeaturesCellDivDiscr.bin");
	ofstream fout(filenameOut, ios::binary | ios::out);
	if( fout.is_open() == false )
	{
		cout<<"ERROR: at mainExtractFeatures: file "<<filenameOut<<" could not be opened to save results"<<endl;
		return 3;
	}else{
		cout<<"Saving training samples at "<<filenameOut<<endl;
	}
	int aux = 0;
	streampos foutPosNumFeatures = fout.tellp();//we need it for later to insert the correct number
	fout.write( (char*) (&aux), sizeof(int) ); //numFeatures
	aux = nTotal * symmetryN;
	fout.write( (char*) (&aux), sizeof(int) );//numSamples
	


	//variables needed for main loop
	vector< AnnotationEllipsoid> annotationsVecTM;
	annotationsVecTM.reserve(500);
	size_t annotationsVecPos = 0;
	vector<mylib::Array*> imgVecAux(2*cellDivisionWithTemporalWindow::getTemporalWindowRadius()+1);
	vector<cellDivisionWithTemporalWindow> cdtwVec; //stores features for each lineage where a cell division is detected
	vector<TreeNode< ChildrenTypeLineage >* > divisionNodes;//contains a pointer to each of the divisions in a specific time point

	int numFeatures = 0;

	const int sizeW = dimsImage * (1+dimsImage) / 2;

	for(int frame = iniFrame + temporalWindowRadiusForLogicalRules; frame <= endFrame; frame++)
	{
		//find all annotations for this time point
		annotationsVecTM.clear();
		while(  annotationsVecPos < annotationsVec.size() && annotationsVec[annotationsVecPos].TM == frame )
			annotationsVecTM.push_back( annotationsVec[annotationsVecPos++] );//subset of annotations where the division happens in this time point
		
		if( annotationsVecTM.empty() == false )
		{
			//------------find a match between annotation and lht position---------------------
			divisionNodes.resize(annotationsVecTM.size());			
			for(size_t ii = 0; ii < divisionNodes.size(); ii++)
			{				
				float d;
				double centroid[dimsImage];
				memcpy(centroid, annotationsVecTM[ii].mu, sizeof(double) * dimsImage );				
			
				//currently linear approach (I could improve later with kd-tree or knn-cuda). Not to worried since this code is to generate samples
				divisionNodes[ii] = NULL;
				for( list<nucleus>::const_iterator iterN = lht.nucleiList[frame].begin(); iterN != lht.nucleiList[frame].end(); ++iterN )
				{
					d = iterN->Euclidean2Distance(centroid, supervoxel::getScale() );
					if( d < 1.0f )//we found the match
					{
						divisionNodes[ii] = iterN->treeNodePtr;						
						break;
					}
				}
				//check that we found the right match
				if( divisionNodes[ii] == NULL )
				{
					cout<<"ERROR: mainCellDivisionWithTemporalWindow: we could not find a match between annotation and list of points from XML file. ii = "<<ii<<endl;;
					annotationsVecTM[ii].writeXML(cout);
					return 3;
				}
			}

			//--------------calculate features: image-based, trajectory-based, geometry-based------------------
			for (int ii = 0; ii < 2 * temporalWindowRadiusForLogicalRules + 1; ii++)
			{
				int frameAux = frame - temporalWindowRadiusForLogicalRules + ii;
				//check if we need to load any extra image (annotations can be far apart and we do not want to load unnecessary images)
				if (imgVec[frameAux] == NULL)
				{
					
					imgVec[frameAux] = readImage(imgFilePattern, frameAux);
				}

				//copy pointer to auxiliary vectors
				imgVecAux[ii] = imgVec[frameAux];
			}
			//calculate average precision matrix if necessary
			if( cellDivisionWithTemporalWindow::getUseFixPrecisionMatrix() == true )
				cellDivisionWithTemporalWindow::averagePrecisionMatrixAyTM( lht.nucleiList[frame] );

			//iterate over different symmetries to extent the training set
			for(int ss = 0; ss < symmetryN; ss++)
			{
				//extract basic Haar features for each time point
				cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchForCellDivision(divisionNodes, imgVecAux, cdtwVec, devCUDA, ss);

				//extend Haar features using temporal combinations
				int countY = 0;
				for(vector<cellDivisionWithTemporalWindow>::iterator iterF = cdtwVec.begin(); iterF != cdtwVec.end(); ++iterF, countY++)
				{
					iterF->f.reserve(numFeatures);
					iterF->calculateFeatures();
					//write out features
					iterF->writeToBinary(fout);

					//store label for yTrain
					yTrain.push_back( annotationsVecTM[countY].classVal );

					//record for next iteration to reserve appropriate ammount of memory
					numFeatures = iterF->f.size();
				}

				
			}
		}

		#ifdef CDWT_SAVE_FEATURE_NAME
	//write feature names
	string filenameS (TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin/" + "trainFeaturesCellDivDiscr_FeatName.txt");
	cdtwVec.front().writeFeaturesNames(filenameS);
#endif
		cout<<"Added "<<annotationsVecTM.size() * symmetryN <<" samples for TM "<<frame<<endl;

		//update lineage hyper tree temporal window (free oldest element and add a new one)
		int frameOlder = frame - temporalWindowRadiusForLogicalRules;
		int frameNew = frame + temporalWindowRadiusForLogicalRules + 1;        
		err = parseGMMtrackingFilesToHyperTree(imgFilePattern, basenameXML, frameNew, frameNew, lht, false);
		if( err > 0 )
			return err;	
		
		lht.setFrameAsT_o(frameOlder+1);
        eraseTM(time_series_map, frameOlder);
		
		if (imgVec[frameOlder] != NULL)
			mylib::Free_Array(imgVec[frameOlder]);

		imgVec[frameOlder] = NULL;
	}

	//write yTrain
	fout.write( (char*)(&(yTrain[0])), sizeof(float) * yTrain.size() );
	
	//write the appropriate number of features;
	fout.seekp( foutPosNumFeatures );
	fout.write( (char*) (&numFeatures), sizeof(int) );


	//release final memory in sliding window
	fout.close();
	for(vector<mylib::Array*>::iterator iterIm = imgVec.begin(); iterIm != imgVec.end(); ++iterIm)
	{
		if( *iterIm != NULL )
			mylib::Free_Array( *iterIm );
	}

	return 0;
}


//=====================================================================

//=============================================================================
int mainSingleWindowForDaughters( int argc, const char** argv )
{
	cout<<"Generating samples for cell division discrimination using temporal window of features"<<endl;	
	cout<<"Method using the same box for both daughters"<<endl;	
	

	//default values: small test set for debugging
	
	int symmetryN = 1;//calculate features for this image using all possible 8 combinations of eigenvectors directions to artifically increment the traning data	
	string imgFilePattern("C:/Users/Fernando/cppProjects/TrackingGaussianMixtures/NM2013-paperRun/data/dataJP2/TM?????_timeFused_blending/SPC0_CM0_CM1_CHN00_CHN01.fusedStack_?????");
	string TGMMoutputXMLfolder("E:/TGMMruns/GMEMtracking3D_2014_5_2_13_51_14_datasetTest_debugCellDivTempWindowDiscrimination");
	int temporalWindowRadiusForLogicalRules = 3;

	int temporalUndersampling = 2;//It allows us to generate training samples with lower temporal resolution by skipping frames
	
	
	//default values: large sample
	/*
	int symmetryN = 1;//calculate features for this image using all possible 8 combinations of eigenvectors directions to artifically increment the traning data	
	string imgFilePattern("X:/SiMView2/12-08-28/Dme_E1_His2ARFP_01_20120828_144957.corrected/Results/TimeFused.Blending/Dme_E1_His2ARFP.TM????_timeFused_blending/SPC0_TM????_CM0_CM1_CHN00_CHN01.fusedStack");
	string TGMMoutputXMLfolder("E:/TGMMruns/GMEMtracking3D_2014_4_29_2_27_1_dataset12_08_28_drosophila_simview_temporalLRdeactivatedForCellDivisionTraining");
	int temporalWindowRadiusForLogicalRules = 3;
	*/

	//fixed parameters
	string positiveClassLabel("cellDivisionCorrect");//the class in the annotations that is considered a positive sample
	int devCUDA = 0;
	basicEllipticalHaarFeatureVector::useDoGfeatures = true;
	cellDivisionWithTemporalWindow::setUseFixPrecisionMatrix(true);

	if(argc == 6)
	{	
		TGMMoutputXMLfolder = string(argv[1]);
		imgFilePattern = string(argv[2]);
		temporalWindowRadiusForLogicalRules = atoi(argv[3]);
		symmetryN = atoi( argv[4] );//calculate features for this image using all possible 8 combinations of eigenvectors directions to artifically increment the traning data	
		temporalUndersampling = atoi(argv[5]);

	}else if(argc == 1){
		cout<<"Using default parameters for debugging purposes"<<endl;
	}else{
		std::cout<<"ERROR: input arguments are <TGMMoutputXMLfolder> <temporalWindowRadiusForLogicalRules> <imgFilePattern>  <symmetryN> <temporalUndersampling>"<<endl;
		return 2;
	}

	//a simple way (although more memory) is to expand temporal window and then subsample features
	temporalWindowRadiusForLogicalRules *= temporalUndersampling;
	
	//basicEllipticalHaarFeatureVector::useDoGfeatures = false;//uncomment this if you want to test results without DoG features
	if( basicEllipticalHaarFeatureVector::useDoGfeatures == false)
		std::cout<<"======================WARNING: TRAINING WITHOUT DoG FEATURES ENHANCEMENT========================="<<endl;
	else
		std::cout<<"======================WARNING: TRAINING WITH DoG FEATURES ENHANCEMENT========================="<<endl;
	if( cellDivisionWithTemporalWindow::getUseFixPrecisionMatrix() == true)
		std::cout<<"======================WARNING: TRAINING WITH FIX PRECISION MATRIX W ========================="<<endl;
	else
		std::cout<<"======================WARNING: TRAINING WITHOUT FIX PRECISION MATRIX W ========================="<<endl;
	//=====================================================================================

	TicTocTimer tt = tic();

	symmetryN = max(symmetryN, 1);//minimum value
	symmetryN = min(symmetryN, 8);//maximum value
	cellDivisionWithTemporalWindow::setTemporalWindowRadius( temporalWindowRadiusForLogicalRules );
	cellDivisionWithTemporalWindow::setTemporalUndersampling( temporalUndersampling );

	//------------------------------------------------------------------------------
	//store positions of ??? to read time point from annotations
	size_t found = imgFilePattern.find_first_of("?");
	size_t qIniPos = found;
	while ((imgFilePattern[found] == '?') && found != string::npos)
	{
		found++;
		if( found >= imgFilePattern.size() )
			break;
	}
	size_t qEndPos = found;


	//----------------read annotations-----------------------------------------
	string annotationFileFolder(TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin");
	
	vector<string> annotationFiles = get_all_files_within_folder(annotationFileFolder + "/*.xml");


	vector< AnnotationEllipsoid> annotationsVec;
	long long int nPos = 0, nNeg = 0, nTotal = 0;
	for(size_t jj = 0; jj < annotationFiles.size(); jj++)
	{
		string auxFile(TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin/" + annotationFiles[jj]);
		cout<<"Reading annotation file "<<auxFile<<endl;
		XMLNode xMainNode = XMLNode::openFileHelper(auxFile.c_str(),"document");
		int n = xMainNode.nChildNode("Surface");

		annotationsVec.resize( nTotal + n );				
		for(int ii=0;ii<n;ii++)
		{
			annotationsVec[nTotal + ii] = AnnotationEllipsoid(xMainNode,ii);
			//parse class name to labels
			if( annotationsVec[nTotal + ii].className.compare(positiveClassLabel) == 0 ) //cell division
			{
				annotationsVec[nTotal + ii].classVal = 1;
				nPos++;
			}else{
				annotationsVec[nTotal + ii].classVal = 0;
				nNeg++;
			}
			//find out TM from imgFilename and compare it to imgFilenamePattern
			annotationsVec[nTotal + ii].TM = std::stoi(annotationsVec[nTotal + ii].imgFilename.substr(qIniPos, qEndPos-qIniPos) );
		}

		nTotal = annotationsVec.size();
	}
	
	//sort by filename (which is equivalent to sorting by time point)
	sort(annotationsVec.begin(), annotationsVec.end());
	int iniFrame = std::max(annotationsVec.front().TM - temporalWindowRadiusForLogicalRules, 0);
	int endFrame = annotationsVec.back().TM + temporalWindowRadiusForLogicalRules;

	vector<float> yTrain;//to store annotations
	yTrain.reserve( nTotal * symmetryN );
	cout<<"Read "<<nTotal<<" annotations in "<<toc(&tt)<<" secs. Positives = "<<nPos<<"; negatives="<<nNeg<<"; symmetry = "<<symmetryN<<endl;


	//-----------------main loop: read TGMM solution in a sliding window fashion and calculate features for annotations---------------------------------
	string basenameXML(TGMMoutputXMLfolder + "/XML_finalResult_lht/GMEMfinalResult_frame");
	lineageHyperTree lht(endFrame+temporalWindowRadiusForLogicalRules + 2);
	//read initial set of frames: function generates all lht structure (sv, nuclei and lineages)
	int err = parseGMMtrackingFilesToHyperTree(imgFilePattern, basenameXML,iniFrame, iniFrame + temporalWindowRadiusForLogicalRules * 2, lht);	
	if( err > 0 )
		return err;

	//read images
	vector< mylib::Array* > imgVec(endFrame + temporalWindowRadiusForLogicalRules + 2, NULL);
    TimeSeriesMapT time_series_map;
	for(int frame = iniFrame; frame <= iniFrame + temporalWindowRadiusForLogicalRules * 2; frame++)
	{
        auto im=readImage(imgFilePattern, frame);
        imgVec[frame]=im;
        time_series_map[frame]=im;  // (ngc) awkward FIXME

		if( imgVec[frame]->type != mylib::UINT16_TYPE )
		{
			cout<<"ERROR: at mainCellDivisionWithTemporalWindow: code is not ready for non UINT16 images"<<endl;
			exit(3);
		}
#if 0
		//add pointer to supervoxels		
		for(list<supervoxel>::iterator iterS = lht.supervoxelsList[frame].begin(); iterS != lht.supervoxelsList[frame].end(); ++iterS)
		{
			iterS->dataPtr = imgVec[frame]->data;
		}
#endif
	}

	supervoxel::dataSizeInBytes = 2;//for uint16
	for(int ii =0; ii < dimsImage; ii++)
	{
		supervoxel::dataDims[ii] = imgVec[iniFrame]->dims[ii];
		supervoxel::dataSizeInBytes *= supervoxel::dataDims[ii];
	}

	//prepare output binary file 
	//header with info (we'll rewrite number of samples at the end)	
	//open file for annotations
	cout<<"=========TODO: add specifics to the name of the file (DoG, numFeatures, temporalWIndowSize, etc)"<<endl;
	string filenameOut (TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin/" + "trainFeaturesCellDivDiscr.bin");
	ofstream fout(filenameOut, ios::binary | ios::out);
	if( fout.is_open() == false )
	{
		cout<<"ERROR: at mainExtractFeatures: file "<<filenameOut<<" could not be opened to save results"<<endl;
		return 3;
	}else{
		cout<<"Saving training samples at "<<filenameOut<<endl;
	}
	int aux = 0;
	streampos foutPosNumFeatures = fout.tellp();//we need it for later to insert the correct number
	fout.write( (char*) (&aux), sizeof(int) ); //numFeatures
	//aux = nTotal * symmetryN;
	streampos foutPosNumSamples = fout.tellp();//we need it for later to insert the correct number
	fout.write( (char*) (&aux), sizeof(int) );//numSamples
	


	//variables needed for main loop
	vector< AnnotationEllipsoid> annotationsVecTM;
	annotationsVecTM.reserve(500);
	size_t annotationsVecPos = 0;
	vector<mylib::Array*> imgVecAux(2*cellDivisionWithTemporalWindow::getTemporalWindowRadius()+1);
	vector<cellDivisionWithTemporalWindow> cdtwVec; //stores features for each lineage where a cell division is detected
	vector<TreeNode< ChildrenTypeLineage >* > divisionNodes;//contains a pointer to each of the divisions in a specific time point

	int numFeatures = 0;	

	const int sizeW = dimsImage * (1+dimsImage) / 2;

	for(int frame = iniFrame + temporalWindowRadiusForLogicalRules; frame <= endFrame; frame++)
	{
		//find all annotations for this time point
		annotationsVecTM.clear();
		while (annotationsVecPos < annotationsVec.size() && annotationsVec[annotationsVecPos].TM <= frame)
		{
			if (annotationsVec[annotationsVecPos].TM == frame)
			{
				annotationsVecTM.push_back(annotationsVec[annotationsVecPos]);//subset of annotations where the division happens in this time point
			}
			annotationsVecPos++;
		}		
		
		if( annotationsVecTM.empty() == false )
		{
			//------------find a match between annotation and lht position---------------------
			divisionNodes.resize(annotationsVecTM.size());			
			for(size_t ii = 0; ii < divisionNodes.size(); ii++)
			{				
				float d;
				double centroid[dimsImage];
				memcpy(centroid, annotationsVecTM[ii].mu, sizeof(double) * dimsImage );				
			
				//currently linear approach (I could improve later with kd-tree or knn-cuda). Not to worried since this code is to generate samples
				divisionNodes[ii] = NULL;
				for( list<nucleus>::const_iterator iterN = lht.nucleiList[frame].begin(); iterN != lht.nucleiList[frame].end(); ++iterN )
				{
					d = iterN->Euclidean2Distance(centroid, supervoxel::getScale() );
					if( d < 1.0f )//we found the match
					{
						divisionNodes[ii] = iterN->treeNodePtr;						
						break;
					}
				}
				//check that we found the right match
				if( divisionNodes[ii] == NULL )
				{
					cout<<"ERROR: mainCellDivisionWithTemporalWindow: we could not find a match between annotation and list of points from XML file. ii = "<<ii<<endl;;
					annotationsVecTM[ii].writeXML(cout);
					return 3;
				}else{
					if( divisionNodes[ii]->getNumChildren() != 2 )
					{
						cout<<"ERROR: mainCellDivisionWithTemporalWindow: nucleus "<<*(divisionNodes[ii]->data)<<" is not a division"<<endl;
						return 4;
					}
				}

				
			}

			//--------------calculate features: image-based, trajectory-based, geometry-based------------------
			for(int ii = 0; ii < 2* temporalWindowRadiusForLogicalRules + 1; ii++)
				imgVecAux[ii] = imgVec[frame-temporalWindowRadiusForLogicalRules+ii];

			//calculate average precision matrix if necessary
			if( cellDivisionWithTemporalWindow::getUseFixPrecisionMatrix() == true )
			{
				cellDivisionWithTemporalWindow::averagePrecisionMatrixAyTM( lht.nucleiList[frame] );
				//increase W by a constant factor to try to incorporate the cell division
				cellDivisionWithTemporalWindow::scalePrecisionMatrixW(0.1);
			}
			//iterate over different symmetries to extent the training set
			for(int ss = 0; ss < symmetryN; ss++)
			{
				//extract basic Haar features for each time point
				cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchForCellDivisionSingleWindowForDaughters(divisionNodes, imgVecAux, cdtwVec, devCUDA, ss);

				//extend Haar features using temporal combinations
				int countY = 0;
				for(vector<cellDivisionWithTemporalWindow>::iterator iterF = cdtwVec.begin(); iterF != cdtwVec.end(); ++iterF, countY++)
				{
					iterF->f.reserve(numFeatures);
					//iterF->calculateFeaturesSingleWindowForDaughters();
					iterF->calculateFeaturesSingleWindowForDaughters_featureSelection_v1();
					//write out features
					iterF->writeToBinary(fout);
					
					/*
					cout << iterF->cellDivisionPtr->data->TM << "," << iterF->cellDivisionPtr->data->centroid[0] << "," << iterF->cellDivisionPtr->data->centroid[1] << "," << iterF->cellDivisionPtr->data->centroid[2];
					for (int qq = 0; qq < 6; qq++)
						cout << "," << iterF->cellDivisionPtr->data->precisionW[qq];
					cout << endl;
					*/

					//store label for yTrain
					yTrain.push_back( annotationsVecTM[countY].classVal );

					//record for next iteration to reserve appropriate ammount of memory
					numFeatures = iterF->f.size();
				}				
			}
#ifdef CDWT_SAVE_FEATURE_NAME
	//write feature names
	string filenameS (TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin/" + "trainFeaturesCellDivDiscr_FeatName.txt");
	cdtwVec.front().writeFeaturesNames(filenameS);
#endif
		}
		cout<<"Added "<<annotationsVecTM.size() * symmetryN <<" samples for TM "<<frame<<endl;

		//update lineage hyper tree temporal window (free oldest element and add a new one)
		int frameOlder = frame - temporalWindowRadiusForLogicalRules;
		int frameNew = frame + temporalWindowRadiusForLogicalRules + 1;
		err = parseGMMtrackingFilesToHyperTree(imgFilePattern, basenameXML, frameNew, frameNew, lht, false);
		if (err > 0)
		{
			//it means we have no more time poitns to read. So we cannot generate ground truth
			lht.setFrameAsT_o(frameOlder + 1);
            eraseTM(time_series_map, frameOlder);
			//mylib::Free_Array(imgVec[frameOlder]);//not needed since we free the pointer in lht.setFrameAsT
			imgVec[frameOlder] = NULL;
			break;
		}
		imgVec[frameNew] = readImage(imgFilePattern, frameNew);
		//add pointer to supervoxels		
#if 0
		for(list<supervoxel>::iterator iterS = lht.supervoxelsList[frameNew].begin(); iterS != lht.supervoxelsList[frameNew].end(); ++iterS)
		{
			iterS->dataPtr = imgVec[frameNew]->data;
		}
#endif
		
        lht.setFrameAsT_o(frameOlder+1);
        eraseTM(time_series_map, frameOlder);
		//mylib::Free_Array(imgVec[frameOlder]);//not needed since we free the pointer in lht.setFrameAsT
		imgVec[frameOlder] = NULL;
	}

	if (yTrain.empty() == false)
	{
		//write yTrain
		fout.write((char*)(&(yTrain[0])), sizeof(float)* yTrain.size());
	}

	//write the appropriate number of features;
	fout.seekp( foutPosNumFeatures );
	fout.write( (char*) (&numFeatures), sizeof(int) );
	aux = yTrain.size();
	fout.seekp(foutPosNumSamples);
	fout.write((char*)(&aux), sizeof(int));

	
	//release final memory in sliding window
	fout.close();
	for(vector<mylib::Array*>::iterator iterIm = imgVec.begin(); iterIm != imgVec.end(); ++iterIm)
	{
		if( *iterIm != NULL )
			mylib::Free_Array( *iterIm );
	}

	return 0;
}


//================================================================================
//=============================================================================
int debug_mainSingleWindowForDaughters_writeImageBoxes( int argc, const char** argv )
{
	cout<<"==========DEBUGGING CODE: write out boxes=============="<<endl;
	cout<<"================You need to uncomment EllipticalHaarFeatures.cu : line 1252 to write boxes out"<<endl;
	cout<<"Generating samples for cell division discrimination using temporal window of features"<<endl;	
	cout<<"Method using the same box for both daughters"<<endl;	
	

	//default values: small test set for debugging
	
	int symmetryN = 1;//calculate features for this image using all possible 8 combinations of eigenvectors directions to artifically increment the traning data	
	string imgFilePattern("C:/Users/Fernando/cppProjects/TrackingGaussianMixtures/NM2013-paperRun/data/dataJP2/TM?????_timeFused_blending/SPC0_CM0_CM1_CHN00_CHN01.fusedStack_?????");
	string TGMMoutputXMLfolder("E:/TGMMruns/GMEMtracking3D_2014_5_2_13_51_14_datasetTest_debugCellDivTempWindowDiscrimination");
	int temporalWindowRadiusForLogicalRules = 5;
	
	
	//default values: large sample
	/*
	int symmetryN = 1;//calculate features for this image using all possible 8 combinations of eigenvectors directions to artifically increment the traning data	
	string imgFilePattern("X:/SiMView2/12-08-28/Dme_E1_His2ARFP_01_20120828_144957.corrected/Results/TimeFused.Blending/Dme_E1_His2ARFP.TM????_timeFused_blending/SPC0_TM????_CM0_CM1_CHN00_CHN01.fusedStack");
	string TGMMoutputXMLfolder("E:/TGMMruns/GMEMtracking3D_2014_4_29_2_27_1_dataset12_08_28_drosophila_simview_temporalLRdeactivatedForCellDivisionTraining");
	int temporalWindowRadiusForLogicalRules = 3;
	*/

	//fixed parameters
	string positiveClassLabel("cellDivisionCorrect");//the class in the annotations that is considered a positive sample
	int devCUDA = 0;
	basicEllipticalHaarFeatureVector::useDoGfeatures = true;
	cellDivisionWithTemporalWindow::setUseFixPrecisionMatrix(true);

	if(argc == 5)
	{	
		TGMMoutputXMLfolder = string(argv[1]);
		imgFilePattern = string(argv[2]);
		temporalWindowRadiusForLogicalRules = atoi(argv[3]);
		symmetryN = atoi( argv[4] );//calculate features for this image using all possible 8 combinations of eigenvectors directions to artifically increment the traning data	

	}else if(argc == 1){
		cout<<"Using default parameters for debugging purposes"<<endl;
	}else{
		std::cout<<"ERROR: input arguments are <numWeakClassifiers> <J> <symmetryN> <crossValidationPercentage>"<<endl;
		return 2;
	}

	
	
	//basicEllipticalHaarFeatureVector::useDoGfeatures = false;//uncomment this if you want to test results without DoG features
	if( basicEllipticalHaarFeatureVector::useDoGfeatures == false)
		std::cout<<"======================WARNING: TRAINING WITHOUT DoG FEATURES ENHANCEMENT========================="<<endl;
	else
		std::cout<<"======================WARNING: TRAINING WITH DoG FEATURES ENHANCEMENT========================="<<endl;
	if( cellDivisionWithTemporalWindow::getUseFixPrecisionMatrix() == true)
		std::cout<<"======================WARNING: TRAINING WITH FIX PRECISION MATRIX W ========================="<<endl;
	else
		std::cout<<"======================WARNING: TRAINING WITHOUT FIX PRECISION MATRIX W ========================="<<endl;
	//=====================================================================================

	TicTocTimer tt = tic();

	symmetryN = max(symmetryN, 1);//minimum value
	symmetryN = min(symmetryN, 8);//maximum value
	cellDivisionWithTemporalWindow::setTemporalWindowRadius( temporalWindowRadiusForLogicalRules );

	//------------------------------------------------------------------------------
	//store positions of ??? to read time point from annotations
	size_t found = imgFilePattern.find_first_of("?");
	size_t qIniPos = found;
	while ((imgFilePattern[found] == '?') && found != string::npos)
	{
		found++;
		if( found >= imgFilePattern.size() )
			break;
	}
	size_t qEndPos = found;


	//----------------read annotations-----------------------------------------
	string annotationFileFolder(TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin");
	
	vector<string> annotationFiles = get_all_files_within_folder(annotationFileFolder + "/*.xml");


	vector< AnnotationEllipsoid> annotationsVec;
	long long int nPos = 0, nNeg = 0, nTotal = 0;
	for(size_t jj = 0; jj < annotationFiles.size(); jj++)
	{
		string auxFile(TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin/" + annotationFiles[jj]);
		cout<<"Reading annotation file "<<auxFile<<endl;
		XMLNode xMainNode = XMLNode::openFileHelper(auxFile.c_str(),"document");
		int n = xMainNode.nChildNode("Surface");

		annotationsVec.resize( nTotal + n );				
		for(int ii=0;ii<n;ii++)
		{
			annotationsVec[nTotal + ii] = AnnotationEllipsoid(xMainNode,ii);
			//parse class name to labels
			if( annotationsVec[nTotal + ii].className.compare(positiveClassLabel) == 0 ) //cell division
			{
				annotationsVec[nTotal + ii].classVal = 1;
				nPos++;
			}else{
				annotationsVec[nTotal + ii].classVal = 0;
				nNeg++;
			}
			//find out TM from imgFilename and compare it to imgFilenamePattern
			annotationsVec[nTotal + ii].TM = std::stoi(annotationsVec[nTotal + ii].imgFilename.substr(qIniPos, qEndPos-qIniPos) );
		}

		nTotal = annotationsVec.size();
	}
	
	//sort by filename (which is equivalent to sorting by time point)
	sort(annotationsVec.begin(), annotationsVec.end());
	int iniFrame = std::max(annotationsVec.front().TM - temporalWindowRadiusForLogicalRules, 0);
	int endFrame = annotationsVec.back().TM + temporalWindowRadiusForLogicalRules;

	vector<float> yTrain;//to store annotations
	yTrain.reserve( nTotal * symmetryN );
	cout<<"Read "<<nTotal<<" annotations in "<<toc(&tt)<<" secs. Positives = "<<nPos<<"; negatives="<<nNeg<<"; symmetry = "<<symmetryN<<endl;


	//-----------------main loop: read TGMM solution in a sliding window fashion and calculate features for annotations---------------------------------
	string basenameXML(TGMMoutputXMLfolder + "/XML_finalResult_lht/GMEMfinalResult_frame");
	lineageHyperTree lht(endFrame+temporalWindowRadiusForLogicalRules + 2);
	//read initial set of frames: function generates all lht structure (sv, nuclei and lineages)
	int err = parseGMMtrackingFilesToHyperTree(imgFilePattern, basenameXML,iniFrame, iniFrame + temporalWindowRadiusForLogicalRules * 2, lht);	
	if( err > 0 )
		return err;

	//read images
	vector< mylib::Array* > imgVec(endFrame + temporalWindowRadiusForLogicalRules + 2, NULL);
    TimeSeriesMapT time_series_map;
	for(int frame = iniFrame; frame <= iniFrame + temporalWindowRadiusForLogicalRules * 2; frame++)
	{
        auto im =readImage(imgFilePattern,frame);
        imgVec[frame] = im; 
        time_series_map[frame]=im; // (ngc) awkward FIXME


		if( imgVec[frame]->type != mylib::UINT16_TYPE )
		{
			cout<<"ERROR: at mainCellDivisionWithTemporalWindow: code is not ready for non UINT16 images"<<endl;
			exit(3);
		}
#if 0
		//add pointer to supervoxels		
		for(list<supervoxel>::iterator iterS = lht.supervoxelsList[frame].begin(); iterS != lht.supervoxelsList[frame].end(); ++iterS)
		{
			iterS->dataPtr = imgVec[frame]->data;
		}
#endif
	}

	supervoxel::dataSizeInBytes = 2;//for uint16
	for(int ii =0; ii < dimsImage; ii++)
	{
		supervoxel::dataDims[ii] = imgVec[iniFrame]->dims[ii];
		supervoxel::dataSizeInBytes *= supervoxel::dataDims[ii];
	}

	//prepare output binary file 
	//header with info (we'll rewrite number of samples at the end)	
	//open file for annotations		
	string foutS("E:/temp/3DHaarBoxes/aaa_boxIndex.txt");
	ofstream foutD(foutS.c_str());
	if( foutD.is_open() == false )
	{
		cout<<"ERROR: could not open file "<< foutS<<endl;
		exit(3);
	}


	//variables needed for main loop
	vector< AnnotationEllipsoid> annotationsVecTM;
	annotationsVecTM.reserve(500);
	size_t annotationsVecPos = 0;
	vector<mylib::Array*> imgVecAux(2*cellDivisionWithTemporalWindow::getTemporalWindowRadius()+1);
	vector<cellDivisionWithTemporalWindow> cdtwVec; //stores features for each lineage where a cell division is detected
	vector<TreeNode< ChildrenTypeLineage >* > divisionNodes;//contains a pointer to each of the divisions in a specific time point

	int numFeatures = 0;

	const int sizeW = dimsImage * (1+dimsImage) / 2;
	int countBox = 0, sampleId = 0;
	for(int frame = iniFrame + temporalWindowRadiusForLogicalRules; frame <= endFrame; frame++)
	{
		//find all annotations for this time point
		annotationsVecTM.clear();
		while(  annotationsVecPos < annotationsVec.size() && annotationsVec[annotationsVecPos].TM == frame )
			annotationsVecTM.push_back( annotationsVec[annotationsVecPos++] );//subset of annotations where the division happens in this time point
		
		if( annotationsVecTM.empty() == false )
		{
			//------------find a match between annotation and lht position---------------------
			divisionNodes.resize(annotationsVecTM.size());			
			for(size_t ii = 0; ii < divisionNodes.size(); ii++)
			{				
				float d;
				double centroid[dimsImage];
				memcpy(centroid, annotationsVecTM[ii].mu, sizeof(double) * dimsImage );				
			
				//currently linear approach (I could improve later with kd-tree or knn-cuda). Not to worried since this code is to generate samples
				divisionNodes[ii] = NULL;
				for( list<nucleus>::const_iterator iterN = lht.nucleiList[frame].begin(); iterN != lht.nucleiList[frame].end(); ++iterN )
				{
					d = iterN->Euclidean2Distance(centroid, supervoxel::getScale() );
					if( d < 1.0f )//we found the match
					{
						divisionNodes[ii] = iterN->treeNodePtr;						
						break;
					}
				}
				//check that we found the right match
				if( divisionNodes[ii] == NULL )
				{
					cout<<"ERROR: mainCellDivisionWithTemporalWindow: we could not find a match between annotation and list of points from XML file. ii = "<<ii<<endl;;
					annotationsVecTM[ii].writeXML(cout);
					return 3;
				}else{
					if( divisionNodes[ii]->getNumChildren() != 2 )
					{
						cout<<"ERROR: mainCellDivisionWithTemporalWindow: nucleus "<<*(divisionNodes[ii]->data)<<" is not a division"<<endl;
						return 4;
					}
				}

				
			}

			//--------------calculate features: image-based, trajectory-based, geometry-based------------------
			for(int ii = 0; ii < 2* temporalWindowRadiusForLogicalRules + 1; ii++)
				imgVecAux[ii] = imgVec[frame-temporalWindowRadiusForLogicalRules+ii];

			//calculate average precision matrix if necessary
			if( cellDivisionWithTemporalWindow::getUseFixPrecisionMatrix() == true )
			{
				cellDivisionWithTemporalWindow::averagePrecisionMatrixAyTM( lht.nucleiList[frame] );
				//increase W by a constant factor to try to incorporate the cell division
				cellDivisionWithTemporalWindow::scalePrecisionMatrixW(0.1);
			}
			//iterate over different symmetries to extent the training set
			for(int ss = 0; ss < symmetryN; ss++)
			{
				//extract basic Haar features for each time point
				cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchForCellDivisionSingleWindowForDaughters(divisionNodes, imgVecAux, cdtwVec, devCUDA, ss);

				//to identify each box uniquely
				vector<int> sampleIdVec;
				for(vector<cellDivisionWithTemporalWindow>::iterator iterF = cdtwVec.begin(); iterF != cdtwVec.end(); ++iterF)
				{
					sampleIdVec.push_back(sampleId++);//to be able to tag each box uniquelyu
				}

				//write out in a text file a mapping for each box
				for( int tt = cellDivisionWithTemporalWindow::getTemporalWindowRadius(); tt >= 0; tt--)	
				{
					int countY = 0;
					for(vector<cellDivisionWithTemporalWindow>::iterator iterF = cdtwVec.begin(); iterF != cdtwVec.end(); ++iterF, countY++)
					{
						//write out in a text file a mapping for each box
						foutD<<countBox<<","<<sampleIdVec[countY]<<","<<frame+tt<<","<<annotationsVecTM[countY].classVal<<";"<<endl;//boxId (for bin file), sample id, time point, yTrain value
						countBox++;							
					}				
				}

				for( int tt = cellDivisionWithTemporalWindow::getTemporalWindowRadius() + 1; tt < 2* cellDivisionWithTemporalWindow::getTemporalWindowRadius() + 1; tt++)
				{
					int countY = 0;
					for(vector<cellDivisionWithTemporalWindow>::iterator iterF = cdtwVec.begin(); iterF != cdtwVec.end(); ++iterF, countY++)
					{
						//write out in a text file a mapping for each box
						foutD<<countBox<<","<<sampleIdVec[countY]<<","<<frame+tt<<","<<annotationsVecTM[countY].classVal<<";"<<endl;//boxId (for bin file), sample id, time point, yTrain value
						countBox++;							
					}
				}
			}
#ifdef CDWT_SAVE_FEATURE_NAME
	//write feature names
	string filenameS (TGMMoutputXMLfolder + "/annForCellDivDiscrWithTempWin/" + "trainFeaturesCellDivDiscr_FeatName.txt");
	cdtwVec.front().writeFeaturesNames(filenameS);
#endif
		}
		cout<<"Added "<<annotationsVecTM.size() * symmetryN <<" samples for TM "<<frame<<endl;

		//update lineage hyper tree temporal window (free oldest element and add a new one)
		int frameOlder = frame - temporalWindowRadiusForLogicalRules;
		int frameNew = frame + temporalWindowRadiusForLogicalRules + 1;
		err = parseGMMtrackingFilesToHyperTree(imgFilePattern, basenameXML, frameNew, frameNew, lht, false);
		if( err > 0 )
			return err;
		imgVec[frameNew] = readImage(imgFilePattern, frameNew);
#if 0
		//add pointer to supervoxels		
		for(list<supervoxel>::iterator iterS = lht.supervoxelsList[frameNew].begin(); iterS != lht.supervoxelsList[frameNew].end(); ++iterS)
		{
			iterS->dataPtr = imgVec[frameNew]->data;
		}
#endif
		
		lht.setFrameAsT_o(frameOlder+1);
        eraseTM(time_series_map, frameOlder);
		
		imgVec[frameOlder] = NULL;
	}

	for(vector<mylib::Array*>::iterator iterIm = imgVec.begin(); iterIm != imgVec.end(); ++iterIm)
	{
		if( *iterIm != NULL )
			mylib::Free_Array( *iterIm );
	} 
	foutD.close();

	return 0;
}
