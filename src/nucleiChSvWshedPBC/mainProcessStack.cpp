/*
 * Copyright (C) 2011-2012 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat
 *
 * mainTest.cpp
 *
 *  Created on: September 17th, 2012
 *      Author: Fernando Amat
 *
 * \brief Generates a hierachical segmentation of a 3D stack and saves teh result in binary format
 *
 */


#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <functional>
#include <time.h>
#include "watershedPersistanceAgglomeration.h"
#ifndef DO_NOT_USE_CUDA
#include "CUDAmedianFilter2D/medianFilter2D.h"
#include "cuda_runtime_api.h"
#else
#include "MedianFilter2D/medianFilter2D.h"
#endif
#ifdef PICTOOLS_JP2K
#include "ioFunctions.h"
#endif

#include "klb_imageIO.h"

#include "parseConfigFile.h"
#include "autoThreshold.h"

namespace mylib
{
	#include "../temporalLogicalRules/mylib/array.h"
	#include "../temporalLogicalRules/mylib/image.h"
}


using namespace std;

void generateSegmentationMaskFromHS(string hsFilename, int tau, size_t minSvSize);
void parseImageFilePattern(string& imgRawPath, int frame);

static void writeArrayToKLB(const char* filename,mylib::Array *a);

struct Params {
    imgVoxelType backgroundThr;
    int conn3D;
    imgVoxelType minTau;
    int radiusMedianFilter;
};

static string make_hs_output_filename(const string basename, Params & params)
{
	char suffix[256] = { 0 };
	sprintf(suffix, "_seg_conn%d_rad%d", params.conn3D, params.radiusMedianFilter,(int)params.minTau,(int)params.backgroundThr);
	return basename + string(suffix) + string(".bin");
}

static string make_param_output_filename(const string basename, Params & params) {
	char suffix[256] = { 0 };
	sprintf(suffix, "_seg_conn%d_rad%d", params.conn3D, params.radiusMedianFilter,(int)params.minTau,(int)params.backgroundThr);
	return basename + string(suffix) + string(".txt");
}


#define TEST_OUTPUT_FILE
#ifdef TEST_OUTPUT_FILE
static int test_output_file(const string basename, Params & params )
{
	int ok = 1;
	vector<function<string(const string, Params&)>> gs = {
		make_hs_output_filename,
		make_param_output_filename
	};
	for (auto g : gs) {
		auto fileOutHS = g(basename, params)+string(".tmp");
		ofstream os(fileOutHS.c_str(), ios::binary | ios::out);
		if (!os.is_open())
		{
			cout << "ERROR: could not open file " << fileOutHS << " to save hierarchical segmentation" << endl;
			perror("");
			ok=0;
		}
		os << "test data";
		os.close();
		remove(fileOutHS.c_str());
	}
	if (ok)
		cout << "OPENED TEST OUTPUT FILES JUST FINE.  Removing the test files." << endl;
	return ok;
}
#endif

int main( int argc, const char** argv )
{
#ifndef DO_NOT_USE_CUDA
	int devCUDA = 0;
#endif
	cout << "ProcessStack Version " << GIT_TAG << " " << GIT_HASH << endl;

	if( argc == 4)//we have a .bin file from hierarchical segmentation and we want to output a segmentation for a specific tau
	{
		//Call function with ProcessStack <hsFilename.bin> <tau> <minSupervoxelSize>

		string hsFilename(argv[1]);
		int tau = atoi(argv[2]);
		size_t minSvSize = atoi(argv[3]);

		generateSegmentationMaskFromHS(hsFilename, tau, minSvSize);

		return 0;
	}


	//parse input parameters
	string basename;//filename without extension so we can save our binary hierarchical segmentation
	int radiusMedianFilter = 0;//if radius = 1->3x3 medianFilter
	int minTau = 0;
	int backgroundThr = 0;
	int conn3D = 0;

	if( argc == 3 ) //we call program wiht <configFile> <timePoint>
	{
		configOptionsTrackingGaussianMixture configOptions;
		string configFilename(argv[1]);
		int err = configOptions.parseConfigFileTrackingGaussianMixture(configFilename);
		if( err != 0 )
			return err;

		int frame = atoi( argv[2] );

		#ifndef DO_NOT_USE_CUDA
		{
			int ndevices;
			cudaGetDeviceCount(&ndevices);
			devCUDA = frame%ndevices;
			cout << "Detected " << ndevices << " CUDA devices. Using device " << devCUDA << "." << endl;
		}
		#endif


		basename = configOptions.imgFilePattern;
		parseImageFilePattern(basename, frame);

		radiusMedianFilter = configOptions.radiusMedianFilter;
		minTau = configOptions.minTau;
		backgroundThr = configOptions.backgroundThreshold;
		conn3D = configOptions.conn3D;

	}else if( argc == 6)
	{
		basename = string(argv[1]);//filename without extension so we can save our binary hierarchical segmentation
		radiusMedianFilter = atoi(argv[2]);//if radius = 1->3x3 medianFilter
		minTau = atoi(argv[3]);
		backgroundThr = atoi(argv[4]);
		conn3D = atoi(argv[5]);

	}else{
		cout<<"Wrong number of parameters. Call function with <imgBasename> <radiusMedianFilter> <minTau> <backgroundThr>  <conn3D>"<<endl;
		cout<<"Wrong number of parameters. Call function with <configFile> <frame>"<<endl;
		cout << "Wrong number of parameters. Call function with <hsBinFileName> <tau> <minSuperVoxelSizePx>" << endl;
		return 2;
	}

  Params params;
  params.conn3D=conn3D;
  params.radiusMedianFilter=radiusMedianFilter;
  params.minTau=minTau;
  params.backgroundThr=backgroundThr;

#ifdef TEST_OUTPUT_FILE
	if (!test_output_file(basename, params)) {
		cout << "Exiting" << endl;
		exit(1);
	}
#endif

	//================================================================================

	//read current image
	string imgFilename(basename + ".tif");
	mylib::Array* img = mylib::Read_Image((char*)(imgFilename.c_str()),0);


	mylib::Value_Type type;
	int ndims;
	mylib::Dimn_Type  *dims = NULL;
	//try to read KLB format
	if (img == NULL)
	{
		imgFilename = string((basename + ".klb").c_str());
		//check if file exists
		FILE* fcheck = fopen(imgFilename.c_str(), "rb");
		if (fcheck != NULL)
		{
			fclose(fcheck);
			klb_imageIO imgKLB(imgFilename);
			int err = imgKLB.readHeader();
			if (err > 0)
				return err;

			void* data = malloc(imgKLB.header.getImageSizeBytes());

#ifdef DO_NOT_USE_CUDA
      // Use a single thread for reading when targeting no-GPU
      // Usually, we're doing no-GPU because we're spawning 
      // ProcessStack jobs on a cluster.  There it's more
      // important to rely on process-level parallelism than
      // multithreading: we wan't to target single node machines.
      // So, we limit to a single thread here to be nice for 
      // cluster environments.
			err = imgKLB.readImageFull((char*)data, 1);  // single thread
#else
			err = imgKLB.readImageFull((char*)data, 0);  // all the threads
#endif
			if (err > 0)
				return err;
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
				return 15;
			}
			ndims = 0;
			if (dims == NULL)
				dims = new mylib::Dimn_Type[KLB_DATA_DIMS];
			while (ndims < KLB_DATA_DIMS && imgKLB.header.xyzct[ndims] > 1)
			{
				dims[ndims] = imgKLB.header.xyzct[ndims];
				ndims++;
			}

			img = mylib::Make_Array_Of_Data(mylib::PLAIN_KIND, type, ndims, dims, data);
		}
	}



	if( img == NULL )
	{
#ifdef PICTOOLS_JP2K
		//try to read JP2 image
		imgFilename = string(basename + ".jp2");
		void* data = readJP2Kfile(imgFilename, type, ndims, &dims);


		if( data == NULL)
		{
			cout<<"ERROR: could not open file "<<imgFilename<<" to read input image"<<endl;
			return 5;
		}

		img = mylib::Make_Array_Of_Data(mylib::PLAIN_KIND, type, ndims, dims, data);
#else
		cout<<"ERROR: could not open file "<<imgFilename<<" to read input image"<<endl;
		return 5;
#endif
	}

  cout << "Loaded: " << imgFilename << endl;

	//hack to make the code work for uin8 without changing everything to templates
	//basically, parse the image to uint16, since the code was designed for uint16
	if( img->type == mylib::UINT8_TYPE )
	{
		img = mylib::Convert_Array_Inplace (img, img->kind, mylib::UINT16_TYPE, 16, 0);
	}

	//hack to make the code work for 2D without changing everything to templates
	//basically, add one black slice to the image (you should select conn3D = 4 or 8)
	if( img->ndims == 2 )
	{
		mylib::Dimn_Type dimsAux[dimsImage];
		for(int ii = 0; ii<img->ndims; ii++)
			dimsAux[ii] = img->dims[ii];
		for(int ii = img->ndims; ii<dimsImage; ii++)
			dimsAux[ii] = 2;

		mylib::Array *imgAux = mylib::Make_Array(img->kind, img->type, dimsImage, dimsAux);
		memset(imgAux->data,0, (imgAux->size) * sizeof(mylib::uint16) );
		memcpy(imgAux->data, img->data, img->size * sizeof(mylib::uint16) );

		mylib::Array* imgSwap = imgAux;
		img = imgAux;
		mylib::Free_Array( imgSwap);
	}

	if( img->type != mylib::UINT16_TYPE )
	{
		cout<<"ERROR: code is not ready for this types of images (change imgVoxelType and recompile)"<<endl;
		return 2;
	}

	//calculate background automatically
	if (backgroundThr < 0)
	{
		switch (backgroundThr)
		{
		case -1:
		{
				   //Triangle method on
				   imgVoxelType* imMIP = maximumIntensityProjection<imgVoxelType, dimsImage>((imgVoxelType*)(img->data), img->dims);
				   backgroundThr = (int)(bckgAutoThr_TriangleMethod<imgVoxelType, dimsImage-1>(imMIP, img->dims));
				   backgroundThr -= 20;//be conservative
				   if (backgroundThr < 0)
					   backgroundThr = 1;

				   delete[] imMIP;

				   cout << "Automatic threshold with MIP + triangle method is " << backgroundThr << endl;
		}
			break;
		default:
			cout << "ERROR: code does not recognize option " << backgroundThr << " for automatic background calculation " << endl;
			return 10;
		}
	}


	//calculate median filter
#ifndef DO_NOT_USE_CUDA
	medianFilterCUDASliceBySlice((imgVoxelType*) (img->data), img->dims, radiusMedianFilter,devCUDA);
#else
	medianFilter2DSliceBySlice((imgVoxelType*)(img->data), img->dims, radiusMedianFilter);

#endif
    writeArrayToKLB("median.klb",img);

	//build hierarchical tree
	//cout<<"DEBUGGING: building hierarchical tree"<<endl;
	int64 imgDims[dimsImage];
	for(int ii = 0;ii<dimsImage; ii++)
		imgDims[ii] = img->dims[ii];
	hierarchicalSegmentation* hs =  buildHierarchicalSegmentation((imgVoxelType*) (img->data), imgDims, backgroundThr, conn3D, minTau, 1);

	//save hierarchical histogram
	//cout<<"DEBUGGING: saving hierarchical histogram"<<endl;
	string fileOutHS = make_hs_output_filename(basename, params);
	ofstream os(fileOutHS.c_str(), ios::binary | ios:: out);
	if( !os.is_open() )
	{
		cout<<"ERROR: could not open file "<< fileOutHS<<" to save hierarchical segmentation"<<endl;
		return 3;
	}
  cout<<"Writing to "<<fileOutHS<<endl;
	hs->writeToBinary( os );
	os.close();
  cout<<"Done writing." << endl;

	//save parameters
	fileOutHS = make_param_output_filename(basename, params);
	ofstream osTxt(fileOutHS.c_str());
	if( !osTxt.is_open() )
	{
		cout<<"ERROR: could not open file "<< fileOutHS<<" to save hierarchical segmentation paremeters"<<endl;
		return 3;
	}

	osTxt<<"Image basename = "<<basename<<endl;
	osTxt<<"Radius median filter = "<<radiusMedianFilter<<endl;
	osTxt<<"Min tau ="<<minTau<<endl;
	osTxt<<"Background threshold = "<<backgroundThr<<endl;
	osTxt<<"Conn3D = "<<conn3D<<endl;
	osTxt.close();


	// delete hs;  DON'T release memory!  The tree can get pretty big and we're about to quit anyway


	return 0;
}


static void writeArrayToKLB(const char* filename, mylib::Array *a) {
	uint32_t shape[KLB_DATA_DIMS] = { 1,1,1,1,1 };
    const uint32_t block[KLB_DATA_DIMS] = { 64,64,16,1,1 };
	const float pixelSize[KLB_DATA_DIMS] = { 1, 1, 1, 1, 1 };

	auto ndims = min(a->ndims, KLB_DATA_DIMS);
	for (auto i = 0; i < ndims; ++i)
		shape[i] = a->dims[i];

	klb_imageIO imgIO(filename);
	imgIO.header.setHeader(shape, KLB_DATA_TYPE::UINT16_TYPE, pixelSize, block, BZIP2, "");
	imgIO.writeImage((char*)a->data, -1);

}

//==================================================================================
void generateSegmentationMaskFromHS(string hsFilename, int tau, size_t minSvSize)
{

	//read hierarchical segmentation
	ifstream fin(hsFilename.c_str(), ios::binary | ios::in );

	if( fin.is_open() == false )
	{
		cout<<"ERROR: at generateSegmentationMaskFromHS: opening file "<<hsFilename<<endl;
		exit(4);
	}

	hierarchicalSegmentation* hs = new hierarchicalSegmentation(fin);
	fin.close();


	//generate a segmentation for the given tau
	hs->segmentationAtTau(tau);


	//generate array with segmentation mask
	mylib::Dimn_Type dimsVec[dimsImage];
	for(int ii = 0; ii <dimsImage; ii++)
		dimsVec[ii] = supervoxel::dataDims[ii];

	mylib::Array *imgL = mylib::Make_Array(mylib::PLAIN_KIND,mylib::UINT16_TYPE, dimsImage, dimsVec);
	mylib::uint16* imgLptr = (mylib::uint16*)(imgL->data);

	memset(imgLptr, 0, sizeof(mylib::uint16) * (imgL->size) );
	mylib::uint16 numLabels = 0;
	for(vector<supervoxel>::iterator iterS = hs->currentSegmentatioSupervoxel.begin(); iterS != hs->currentSegmentatioSupervoxel.end(); ++iterS)
	{
		if( iterS->PixelIdxList.size() < minSvSize )
			continue;
        if(numLabels==maxNumLabels /*65535*/)
		{
			cout<<"ERROR: at generateSegmentationMaskFromHS: more labels than permitted in imgLabelType"<<endl;
			exit(4);
		}

		numLabels++;
		for(vector<uint64>::iterator iter = iterS->PixelIdxList.begin(); iter != iterS->PixelIdxList.end(); ++iter)
		{
			imgLptr[ *iter ] = numLabels;
		}
	}

	cout<<"A total of "<<numLabels<<" labels for tau="<<tau<<endl;

	//write tiff file
	char buffer[128];
	sprintf(buffer,"_tau%d",tau);
	string suffix(buffer);
	string imgLfilename( hsFilename +  suffix + ".klb");
	writeArrayToKLB((char*)(imgLfilename.c_str()),imgL);

	//write jp2 file
#ifdef PICTOOLS_JP2K
	imgLfilename = string( hsFilename +  suffix + ".jp2");
	writeJP2Kfile(imgL, imgLfilename);
#endif
	cout<<"KLB file written to "<<imgLfilename<<endl;

	//release memory -- or don't bc the process is exiting anyway
#if 0
	delete hs;
	mylib::Free_Array(imgL);
#endif
}



//==================================================================================
//=======================================================================
void parseImageFilePattern(string& imgRawPath, int frame)
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
