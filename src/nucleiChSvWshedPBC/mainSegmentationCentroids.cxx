/*
 * Copyright (C) 2011-2015 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat 
 *
 * mainTest.cpp
 *
 *  Created on: September 17th, 2014
 *      Author: Fernando Amat
 *
 * \brief Small executable for Burkhard to perform segmentation with a specific tau and save centroids in a txtf file
 *
 */


#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include "watershedPersistanceAgglomeration.h"
#include "CUDAmedianFilter2D/medianFilter2D.h"
#ifdef PICTOOLS_JP2K
#include "ioFunctions.h"
#endif

#include "klb_imageIO.h"

#include "parseConfigFile.h"

namespace mylib
{
	#include "../temporalLogicalRules/mylib/array.h"
	#include "../temporalLogicalRules/mylib/image.h"
}

using namespace std;


int main( int argc, const char** argv )
{

	
	//parse input parameters
	//string basename("C:/Users/Fernando/cppProjects/TrackingGaussianMixtures/NM2013-paperRun/data/data/TM00000_timeFused_blending/SPC0_CM0_CM1_CHN00_CHN01.fusedStack_00000");//filename without extension so we can save our binary hierarchical segmentation
	string basename("E:/SPM00_TM000149_CM00_CM01_CHN00.fusedStack");
	int radiusMedianFilter = 2;//if radius = 1->3x3 medianFilter
	int tau = 10;
	int backgroundThr = 200;
	int conn3D = 74;
	float scale[dimsImage] = {1.0f, 1.0f, 1.0f};	
	

	if( argc == 6)
	{
		basename = string(argv[1]);//filename without extension so we can save our binary hierarchical segmentation
		radiusMedianFilter = atoi(argv[2]);//if radius = 1->3x3 medianFilter
		tau = atoi(argv[3]);
		backgroundThr = atoi(argv[4]);
		conn3D = atoi(argv[5]);

	}else if( argc == 9 )//add 3 floating number for scaling pixel units
	{
		basename = string(argv[1]);//filename without extension so we can save our binary hierarchical segmentation
		radiusMedianFilter = atoi(argv[2]);//if radius = 1->3x3 medianFilter
		tau = atoi(argv[3]);
		backgroundThr = atoi(argv[4]);
		conn3D = atoi(argv[5]);

		scale[0] = atof( argv[6] );
		scale[1] = atof( argv[7] );
		scale[2] = atof( argv[8] );
	}else{
		cout<<"Wrong number of parameters. Call function with <imgBasename> <radiusMedianFilter> <tau> <backgroundThr>  <conn3D>"<<endl;
		//return 2;
	}
	
	size_t minSvSize = 50;//in pixels (hard-coded)
	int devCUDA = 0;
	
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

			err = imgKLB.readImageFull((char*)data, 0);
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

	cout<<"Calculating centroids for image "<<imgFilename<<endl;
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
	
	
	//calculate median filter
	medianFilterCUDASliceBySlice((imgVoxelType*) (img->data), img->dims, radiusMedianFilter,devCUDA);
	

	//build hierarchical tree
	//cout<<"DEBUGGING: building hierarchical tree"<<endl;
	int64 imgDims[dimsImage];
	for(int ii = 0;ii<dimsImage; ii++)
		imgDims[ii] = img->dims[ii];
	hierarchicalSegmentation* hs =  buildHierarchicalSegmentation((imgVoxelType*) (img->data), imgDims, backgroundThr, conn3D, tau, 1);//we use minTau = tau -> basic regions contain the segmentation already
	
	//calculate weighted centroids for each basic region
	for(unsigned int ii = 0; ii < hs->getNumberOfBasicRegions(); ii++)
	{
		//TODO: trim supervoxel if we see that centroid accuracy is not sufficient
		hs->basicRegionsVec[ii].weightedCentroid<mylib::uint16>();//calculated weighted centroid
	}

	//save centroids
	char buffer[256];
	sprintf(buffer, "%s_viewSetupId_0.beads.txt", basename.c_str());
	string filenameOut(buffer);
	ofstream os(filenameOut.c_str());
	if( !os.is_open() )
	{
		cout<<"ERROR: could not open file "<< filenameOut<<" to save centroids"<<endl;
		return 3;
	}
	os<<"id"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"z"<<endl;
	for(unsigned int ii = 0; ii < hs->getNumberOfBasicRegions(); ii++)
	{
		//TODO: trim supervoxel if we see that centroid accuracy is not sufficient
		os<<ii<<"\t"<<hs->basicRegionsVec[ii].centroid[0] * scale[0]<<"\t"<<hs->basicRegionsVec[ii].centroid[1]  * scale[1]<<"\t"<<hs->basicRegionsVec[ii].centroid[2] * scale[2]<<endl;//calculated weighted centroid
	}
	os.close();	

	//release memory
	delete hs;


	return 0;
}






