/*
 * mainTrackingGaussianMixtureModel.cpp
 *
 *  Created on: May 12, 2011
 *      Author: amatf
 */

/*
 * Copyright (C) 2020-2021 Martin Dominguez
 * additional code for FTGMM project
 * 
 * Gladstone Institutes
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <ctime>
#include <list>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <utility>
#include <time.h>
#include <map>
#include <stdexcept>
#include "GaussianMixtureModel.h"
#include "responsibilities.h"
#include "variationalInference.h"
#include "external/xmlParser2/xmlParser.h"
#include "Utils/parseConfigFile.h"
#include "external/Nathan/tictoc.h"
#include "sys/stat.h"
#include "UtilsCUDA/GMEMupdateCUDA.h"
#include "constants.h"
#include "temporalLogicalRules/temporalLogicalRules.h"
#include "UtilsCUDA/3DEllipticalHaarFeatures/EllipticalHaarFeatures.h"
#include "UtilsCUDA/3DEllipticalHaarFeatures/gentleBoost/gentleBoost.h"
#include "trackletCalculation.h"

#include "DivisionClassifierFactory.h"
#include "cellDivision.h"

#ifndef DO_NOT_USE_CUDA
#include "CUDAmedianFilter2D/medianFilter2D.h"
#else
#include "MedianFilter2D/medianFilter2D.h"
#endif

#ifdef PICTOOLS_JP2K
#include "ioFunctions.h"
#endif

#include "klb_imageIO.h"

#include "hierarchicalSegmentation.h"
#include "supportFunctionsMain.h"

#include "backgroundDetectionInterface.h"
#include <cuda_runtime_api.h>
#include <UtilsCUDA/external/book.h>

#if defined(_WIN32) || defined(_WIN64)
#define NOMINMAX
#include <Windows.h>
#include <Psapi.h> // for checking memory usage
#include "Shlwapi.h"
#pragma comment(lib, "Shlwapi.lib")
#else //Linux headers
#include <unistd.h>
#endif

#ifdef CELL_DIVISION_WITH_GPU_3DHAAR
#include "CellDivisionWithGPU3DHaar.h"
#endif //  CELL_DIVISION_WITH_GPU_3DHAAR


#define CHECK(e) do{if(!(e)) {throw runtime_error(string(__FUNCTION__)+string("() : Expression evaluated as false. " #e));}} while(0)

#if defined(_WIN32) || defined(_WIN64)
#define TIME(e) do{TicTocTimer tttttt=tic();(e);printf("Elapsed %10.2f s for %s\n",toc(&tttttt),#e);}while(0)
#else //Linux
//#define TIME(e) do{time_t start = time(NULL);(e);printf("Elapsed %6.0f sec for %s\n",difftime(time(NULL),start),#e);}while(0)
//#define TIME(e) do{clock_t start = clock();(e);printf("Elapsed %6.0f sec for %s\n",(clock() - start) / (CLOCKS_PER_SEC / 1000),#e);}while(0)
#define TIME(e) do{clock_t c_start = clock();time_t t_start = time(NULL);(e);double t_diff = difftime(time(NULL),t_start);if ( t_diff > 1 ) { printf("Elapsed %6.0f sec for %s\n",difftime(time(NULL),t_start),#e); } else { double diff = (clock() - c_start) / (CLOCKS_PER_SEC / 1000); printf("Elapsed %6.2f sec for %s\n",diff,#e); }}while(0)
#endif


template<typename T> static inline mylib::Value_Type MylibValueType();
template<> inline mylib::Value_Type MylibValueType<float>(){ return mylib::FLOAT32_TYPE; }
template<> inline mylib::Value_Type MylibValueType<unsigned short>(){ return mylib::UINT16_TYPE; }
template<> inline mylib::Value_Type MylibValueType<unsigned char>(){ return mylib::UINT8_TYPE; }

static long long GenerateSegmentationAtSpecificTau(mylib::Array *img, hierarchicalSegmentation* hs, configOptionsTrackingGaussianMixture& configOptions, int frame);

struct TimeSeriesMap {
    TimeSeriesMap( configOptionsTrackingGaussianMixture* configOptions, int TM0)
        : _config(configOptions),
          _TM0(TM0),
          _windowRadius(configOptions->temporalWindowRadiusForLogicalRules)
    {
        imgVecUINT16.reserve(2 * _windowRadius + 1);
    }
    ~TimeSeriesMap() {
#if 0 // don't care -- process is exiting
        for (size_t ii = 0; ii < imgVecUINT16.size(); ii++)
        {
            mylib::Array* imgUINT16 = imgVecUINT16[ii];
            if (imgUINT16 != NULL)
                mylib::Free_Array(imgUINT16);
        }
        imgVecUINT16.clear();

        for (auto &e : _cache)
            flush(e.first);
#endif
    }

    //save raw input data before doing any pre-processing with supervoxels
    //copy image before converting it to float (we need it for image features in the GPU)
    void push_back(mylib::Array* img) {
        CHECK(img->type == MylibValueType<uint16_t>());
        if (imgVecUINT16.size() == 2 * _windowRadius + 1)//already full
        {
            mylib::Free_Array(imgVecUINT16.front());
            imgVecUINT16.erase(imgVecUINT16.begin());//even if this has to copy elements it is only 11 pointers
            imgVecUINT16.push_back(img);
            _TM0++;
        }
        else{
            imgVecUINT16.push_back(img);
        }

    }

    //! \returns the index into a window vector given the target time
    unsigned inline window_pos(int TM) {
        return TM - _TM0;
    }

    long long process(int TM,hierarchicalSegmentation* hs) {
        // Supervoxel trimming, local thresholding and
        // Rescale to 0,1
        mylib::Array* img2;
        {
            auto img = imgVecUINT16[TM - _TM0];
            img2 = mylib::Make_Array(img->kind, img->type, img->ndims, img->dims);
            memcpy(img2->data, img->data, sizeof(mylib::uint16) * (img->size));
        }
        cout << "Computing GenerateSegmentationAtSpecificTau for " << TM << endl;
        auto query_nb = GenerateSegmentationAtSpecificTau(img2, hs, *_config, TM);
        mylib::Convert_Array_Inplace(img2, mylib::PLAIN_KIND, mylib::FLOAT32_TYPE, 1, 1.0);
        {
            mylib::Value v1a, v0a;
            v1a.fval = 1.0; v0a.fval = 0.0;
            mylib::Scale_Array_To_Range(img2, v0a, v1a);
        }
        _cache[TM] = img2;
        return query_nb;
    }

    //! if TM isn't cached, processes and caches TM, otherwise does nothing
    //! \returns true if TM was cached before the call, otherwise false.
    bool maybe_process(int TM, hierarchicalSegmentation *hs) {
        if (!is_cached(TM)) {
            process(TM, hs);
            return false;
        }
        return true;
    }

    float *getProcessed(const ChildrenTypeLineage& iterN)  {
        return getTM<float>(iterN);
    }

    float *getProcessed(int TM) {
        return getTM<float>(TM);
    }

    mylib::Array *getProcessedArray(int TM) {
        return _cache.at(TM);
    }

    // removes timepoint from processed data cache
    void flush(int TM) {
        mylib::Free_Array(_cache.at(TM)); //free because it was allocated using mylib
        _cache.erase(TM);
    }
    bool is_cached(int TM) {
        return _cache.find(TM) != _cache.end();
    }
//data
    vector<mylib::Array*> imgVecUINT16;//store original raw data in the temporal window for GPU Haar ceel division features(or any other thing we need in the future)

private:
    configOptionsTrackingGaussianMixture* _config;
    int _TM0; // TM for index 0 of imgVecUINT16
    int _windowRadius;
    map<int, mylib::Array*> _cache; // caches processed float data

    typedef map<unsigned,mylib::Array*> TimeSeriesCacheT;

    template<typename T> T* getTM(int TM) {
        T *imgPtr=0;
        //cout << "\tQUERY TM: " << TM << endl;
        const auto im=_cache.at(TM);
        CHECK(im->type==MylibValueType<T>());
        imgPtr=(T*)im->data;
        return imgPtr;
    }

    template<typename T> inline T* getTM(const ChildrenTypeLineage& iterN) {
        return getTM<T>(iterN->TM);
    }



};

using namespace std;

//#define DEBUG_EM_ITER //uncomment this line to write XML files with all the steps within an iteration of the Variational inference EM
//#define DEBUG_TGMM_XML_FILES //uncomment this line to write xml files with intermediate steps during the sequential GMM fitting process

//TODO: WORK_WITH_UNITS only works with CUDA code so far
//#define WORK_WITH_UNITS //uncomment this define to work with metric units (instead of pixels). So scale is incoporated directly in all the calculation.

//uncomment this to save a list of supervoxels and ratio so a good threshold can be set
//#define DEBUG_SV_RATIO_HOLES

#ifndef ROUND
#define ROUND(x) (floor(x+0.5))
#endif

#include<cmath>
#include<sstream>
#define countof(e) (sizeof(e)/sizeof(*(e)))
string HumanReadable(uint64_t t) {
	double decades = log((double)t)/log(2.0);
	const char* designator[] = { " ", " kB", " MB", " GB", " TB", " PB" };
	int i = (int) floor(decades / 10);
	i = (i < 0) ? 0 : (i>countof(designator)) ? countof(designator) : i;
	stringstream s;
	s << (t >> (10 * i)) << designator[i];
	return s.str();
}

#if 0
void PrintMemoryInUse() {
	PROCESS_MEMORY_COUNTERS_EX info;
    GetProcessMemoryInfo(GetCurrentProcess(),(PROCESS_MEMORY_COUNTERS*)&info,sizeof(info));
    cout << "Memory Usage"
        << "\tPageFaultCount:             " << HumanReadable(info.PageFaultCount) << endl
        << "\tPeakWorkingSetSize:         " << HumanReadable(info.PeakWorkingSetSize) << endl
        << "\tWorkingSetSize:             " << HumanReadable(info.WorkingSetSize) << endl
        << "\tQuotaPeakPagedPoolUsage:    " << HumanReadable(info.QuotaPeakPagedPoolUsage) << endl
        << "\tQuotaPagedPoolUsage:        " << HumanReadable(info.QuotaPagedPoolUsage) << endl
        << "\tQuotaPeakNonPagedPoolUsage: " << HumanReadable(info.QuotaPeakNonPagedPoolUsage) << endl
        << "\tQuotaNonPagedPoolUsage:     " << HumanReadable(info.QuotaNonPagedPoolUsage) << endl
        << "\tPagefileUsage:              " << HumanReadable(info.PagefileUsage) << endl
        << "\tPeakPagefileUsage:          " << HumanReadable(info.PeakPagefileUsage) << endl
        << endl << endl
        << "\n\tmylib::Array_Usage() = " << mylib::Array_Usage() << endl << endl;
    cout << "--- " << endl;
}
#else
#define PrintMemoryInUse(...)
#endif

static inline bool path_exists(const string path) {
    struct stat St;
    if (!stat(path.c_str(), &St)) {
        CHECK(errno != EINVAL);
        return false;
    }
    return true; //path exists
}

#include <iomanip>
static string to_string(int i, int npad, char padchar) {
    stringstream ss;
    ss << setw(npad) << setfill(padchar) << i;
    return ss.str();
};

//-----------------------------------------------
//to separate filename from path
//from http://www.cplusplus.com/reference/string/string/find_last_of/
//returns the
string getPathFromFilename (const std::string& str)
{
  unsigned found = str.find_last_of("/\\");
  return str.substr(0,found);
}
//-------------------------------------------------
//structure to sort foreground values based on supervoxels identity
struct RegionId
{
    unsigned short int labelId;//id for the region
	long long int posId;//position in the original vector (so we can recover sort and apply it to other vectors)
	long long int idx;//position in the original image (so we can store super voxels)

	friend bool operator< ( RegionId& lhs, RegionId& rhs);
};

bool operator< (RegionId& lhs, RegionId& rhs)
{
    return lhs.labelId<rhs.labelId;
}

//------------------------------------------------
//structure to accomodate offsets in the data: if in imageJ offset is (x,y,z) -> the file contains (t,y,x,z) [Matlab style]
struct OffsetCoordinates
{
	float x,y,z;

	OffsetCoordinates() { x = y = z = 0;};
	OffsetCoordinates (const OffsetCoordinates& p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
	}
	OffsetCoordinates& operator=(const OffsetCoordinates & p)
	{
		if (this != &p)
		{
			x = p.x;
			y = p.y;
			z = p.z;
		}
		return *this;
	}
};

static string make_hs_output_filename1(const string basename, const configOptionsTrackingGaussianMixture& configOptions) {
    char suffix[256] = { 0 };
    sprintf(suffix, "_seg_conn%d_rad%d",
        configOptions.conn3D,
        configOptions.radiusMedianFilter);
    return basename + string(suffix) + string(".bin");
}
static string make_hs_output_filename2(const string basename,const configOptionsTrackingGaussianMixture& configOptions) {
    char suffix[256]={0};
    sprintf(suffix,"_seg_conn%d_rad%d_tau%d_thr%d",
            configOptions.conn3D,
            configOptions.radiusMedianFilter,
            (int)configOptions.minTau,
            (int)configOptions.backgroundThreshold);
    return basename+string(suffix)+string(".bin");
}
static string make_hs_output_filename3(const string basename,const configOptionsTrackingGaussianMixture& configOptions) {
	char suffix[256]={0};
	sprintf(suffix,"_hierarchicalSegmentation_conn3D%d_medFilRad%d",
			configOptions.conn3D,
			configOptions.radiusMedianFilter);
	return basename+string(suffix)+string(".bin");
}

#if 0
#define HERE cout<<"HERE: "<<__FILE__<<"("<<__LINE__<<")" << endl
#else
#define HERE
#endif

/// \returns the location of this executable
static string basename() {
    //from http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
    //for platform-dependant options
    string pathS;
#if defined(_WIN32) || defined(_WIN64)
    char path[MAX_PATH];
    //find out folder of the executable
    HMODULE hModule=GetModuleHandleW(NULL);
    GetModuleFileName(hModule,path,MAX_PATH);
    //change to process stack
    PathRemoveFileSpec(path);
    pathS=string(path);
#elif defined(__APPLE__)
    //find out folder of the executable
    cout<<"ERROR: this part of the code has not been ported to MacOSx yet!!!"<<endl;
    return 2;
#else //we assume Linux
    //find out folder of the executable
    char buff[1024];
    string buffP;
    ssize_t len=readlink("/proc/self/exe",buff,sizeof(buff)-1);
    if(len!=-1) {
        buff[len]='\0';
        buffP=std::string(buff);
    } else {
        /* handle error condition */
    }
    //change to process stack
    pathS=getPathFromFilename(buffP);
#endif
    return pathS;
}

static void DebugSvRatioHoles(string debugPath, string itoa, hierarchicalSegmentation* hs) {
    cout << "=============DEBUGGING: writing out info on super-voxels holes ratio feature to setup a good threshold=====================" << endl;
    ofstream fsvr(debugPath + "debugSvHoleRatio_" + itoa + ".txt");
    for (size_t ii = 0; ii < hs->currentSegmentatioSupervoxel.size(); ii++)
    {
        int numHoles = hs->currentSegmentatioSupervoxel[ii].numHoles();
        float ratioHoles = float(numHoles) / float(hs->currentSegmentatioSupervoxel[ii].PixelIdxList.size());
        for (int jj = 0; jj < dimsImage; jj++)
        {
            fsvr << hs->currentSegmentatioSupervoxel[ii].centroid[jj] << ",";//x,y,z
        }
        fsvr << numHoles<<","<<ratioHoles << ";" << endl;
    }
    fsvr.close();
}

static long long GenerateSegmentationAtSpecificTau(mylib::Array *img,hierarchicalSegmentation* hs,configOptionsTrackingGaussianMixture& configOptions,int frame) {
    long long int query_nb = 0; //necessary to initiate reposnsibilities
    const int conn3Dsv = configOptions.conn3DsvTrim;
    const int minNucleiSize = configOptions.minNucleiSize;
        //generate segmentation mask at specific tau
        cout << "Generating segmentation mask from HS (trimming supervoxels)." << endl;
        TicTocTimer ttHS = tic();

        mylib::uint16* imgDataUINT16 = (mylib::uint16*) (img->data);
        hs->segmentationAtTau(configOptions.tau);

        int64 boundarySize[dimsImage];//margin we need to calculate connected components
        int64* neighOffsetSv = supervoxel::buildNeighboorhoodConnectivity(conn3Dsv, boundarySize);
        mylib::uint16 thr;
        bool* imgVisited = new bool[img->size];
        memset(imgVisited, 0, sizeof(bool) * img->size);//if visited->true (so I can set teh rest to zero)

        vector<int> eraseIdx;
        eraseIdx.reserve(hs->currentSegmentatioSupervoxel.size() / 10);
        int numSplits = 0;
        for (size_t ii = 0; ii < hs->currentSegmentatioSupervoxel.size(); ii++) {
            thr = hs->currentSegmentatioSupervoxel[ii].trimSupervoxel<mylib::uint16>(imgDataUINT16);//trimmming supervoxels

                                                                                       //split supervoxels if necessary (obvius oversegmentations)
            TreeNode<nodeHierarchicalSegmentation>* rootSplit[2];
            supervoxel rootSplitSv[2];
            float scoreS;
            queue<size_t> q;
            q.push(ii);
            while (q.empty() == false) {
                size_t aa = q.front();
                q.pop();
                if (hs->currentSegmentatioSupervoxel[aa].getDeltaZ() > supervoxel::pMergeSplit.deltaZthr) {
                    scoreS = hs->suggestSplit<mylib::uint16>(hs->currentSegmentationNodes[aa], hs->currentSegmentatioSupervoxel[aa], rootSplit, rootSplitSv, imgDataUINT16);
                    if (scoreS > 0.0f) {
                        hs->currentSegmentatioSupervoxel[aa] = rootSplitSv[0];//root split has the correct nodeHSptr
                        hs->currentSegmentatioSupervoxel.push_back(rootSplitSv[1]);
                        hs->currentSegmentationNodes[aa] = rootSplit[0];
                        hs->currentSegmentationNodes.push_back(rootSplit[1]);

                        numSplits++;
                        q.push(aa);
                        q.push(hs->currentSegmentationNodes.size() - 1);
                    }
                    //recalculate thr
                    thr = numeric_limits<mylib::uint16>::max();
                    for (vector<uint64>::const_iterator iter = hs->currentSegmentatioSupervoxel[ii].PixelIdxList.begin(); iter != hs->currentSegmentatioSupervoxel[ii].PixelIdxList.end(); ++iter) {
                        thr = min(thr, imgDataUINT16[*iter]);
                    }
                }
            }

            //calculate ratio of holes / num pixels to see if we need to remove supervoxel
            float ratioHoles = float(hs->currentSegmentatioSupervoxel[ii].numHoles()) / float(hs->currentSegmentatioSupervoxel[ii].PixelIdxList.size());

            if (hs->currentSegmentatioSupervoxel[ii].PixelIdxList.size() < minNucleiSize || ratioHoles > configOptions.ratioHolesSv)//delete supervoxel
            {
                eraseIdx.push_back(ii);
            }
            else {//local background subtraction and reference points for local geometrical descriptors
                hs->currentSegmentatioSupervoxel[ii].weightedGaussianStatistics<mylib::uint16>(true, imgDataUINT16);

                localGeometricDescriptor<dimsImage>::addRefPt(frame, hs->currentSegmentatioSupervoxel[ii].centroid);

                for (vector<uint64>::const_iterator iter = hs->currentSegmentatioSupervoxel[ii].PixelIdxList.begin(); iter != hs->currentSegmentatioSupervoxel[ii].PixelIdxList.end(); ++iter) {
                    //local background subtraction: since these are the trimmed supervoxels, they are guaranteed to be above thr
                    imgDataUINT16[*iter] -= thr;
                    imgVisited[*iter] = true;
                }
            }

            query_nb += hs->currentSegmentatioSupervoxel[ii].PixelIdxList.size();
            hs->currentSegmentatioSupervoxel[ii].localBackgroundSubtracted = true;

        }
        //set to true for all basic regions
        for (unsigned int ii = 0; ii < hs->getNumberOfBasicRegions(); ii++)
            hs->basicRegionsVec[ii].localBackgroundSubtracted = true;

        //delete supervoxels that are too small
        if (eraseIdx.size() >= hs->currentSegmentatioSupervoxel.size()) {
            cout << "ERROR: we are going to delete all the supervoxels!!!" << endl;
            exit(5);
        }
        size_t auxSize = hs->currentSegmentatioSupervoxel.size();
        hs->eraseSupervoxelFromCurrentSegmentation(eraseIdx);
        cout << "Deleted " << eraseIdx.size() << " supervoxels out of " << auxSize << " for being of size<" << minNucleiSize << " after trimming. Left " << hs->currentSegmentatioSupervoxel.size() << endl;
        cout << "Splitted " << numSplits << " supervoxels for having deltaZ >" << supervoxel::pMergeSplit.deltaZthr << endl;


#ifdef DEBUG_SV_RATIO_HOLES
        DebugSvRatioHoles(debugPath,itoa,hs);
#endif


        //set all the non-visited elements to zero
        for (mylib::Size_Type ii = 0; ii < img->size; ii++) {
            if (imgVisited[ii] == false) {
                imgDataUINT16[ii] = 0;
            }
        }

        delete[] neighOffsetSv;
        delete[] imgVisited;
        cout << "Parsing " << hs->currentSegmentatioSupervoxel.size() << " supervoxels from HS took " << toc(&ttHS) << " secs" << endl;

        return query_nb;
}

static mylib::Array* LoadTimepoint(string imgBasename) {
    //try to read JP2 image

    string imgFilename((imgBasename+".jp2").c_str());
    mylib::Value_Type type;
    int ndims;
    mylib::Dimn_Type  *dims=NULL;
    void* data=NULL;
#ifdef PICTOOLS_JP2K
    ifstream fin(imgFilename.c_str());
    if(fin.is_open()==true)//jp2 file exists
    {
        fin.close();
        data=readJP2Kfile(imgFilename,type,ndims,&dims);
    }
#endif
    //try to read KLB format
    if(data==NULL){
        imgFilename=string((imgBasename+".klb").c_str());
        //check if file exists
        FILE* fcheck=fopen(imgFilename.c_str(),"rb");
        if(fcheck!=NULL){
            fclose(fcheck);
            klb_imageIO imgKLB(imgFilename);
            int err=imgKLB.readHeader();
            if(err>0)
                exit(err);

            data=malloc(imgKLB.header.getImageSizeBytes());

            err=imgKLB.readImageFull((char*)data,0);
            if(err>0)
                exit(err);
            //update parameters for mylib
            switch(imgKLB.header.dataType)//we follow the same convention as mylib
            {
                case KLB_DATA_TYPE::UINT8_TYPE:
                    type=mylib::UINT8_TYPE;
                    break;
                case KLB_DATA_TYPE::UINT16_TYPE:
                    type=mylib::UINT16_TYPE;
                    break;
                case KLB_DATA_TYPE::FLOAT32_TYPE:
                    type=mylib::FLOAT32_TYPE;
                    break;
                default:
                    cout<<"ERROR: file type stored in KLB file not supported"<<endl;
                    free(data);
                    exit( 15 );
            }
            ndims=0;
            if(dims==NULL)
                dims=new mylib::Dimn_Type[KLB_DATA_DIMS];
            while(ndims<KLB_DATA_DIMS && imgKLB.header.xyzct[ndims] > 1){
                dims[ndims]=imgKLB.header.xyzct[ndims];
                ndims++;
            }
        }
    }

    mylib::Array *img=NULL;

    if(data==NULL)//try top read tif image
    {
        imgFilename=string((imgBasename+".tif").c_str());
        img=mylib::Read_Image((char*)(imgFilename.c_str()),0);
        if(img==NULL){
            cout<<"ERROR: reading image "<<imgBasename<<endl;
            exit(3);
        }
    } else{//wrap data into mylib::array
        img=mylib::Make_Array_Of_Data(mylib::PLAIN_KIND,type,ndims,dims,data);
    }

    return img;
}

void MaybeMedianFilter(const configOptionsTrackingGaussianMixture &configOptions, int devCUDA, mylib::Array* img) {
    if(configOptions.useMedianFilterForTracking!=0){
        cout<<"Calculating median filter with radius "<<configOptions.radiusMedianFilter<<" to be used during tracking"<<endl;
#ifndef DO_NOT_USE_CUDA
        if(img->ndims==2)
            medianFilterCUDA((imgVoxelType*)(img->data),img->dims,configOptions.radiusMedianFilter,devCUDA);
        else
            medianFilterCUDASliceBySlice((imgVoxelType*)(img->data),img->dims,configOptions.radiusMedianFilter,devCUDA);
#else
            medianFilter2DSliceBySlice((imgVoxelType*)(img->data),img->dims,configOptions.radiusMedianFilter);
#endif
    }
}

static hierarchicalSegmentation* LoadHierarchicalSegmentation(const string &imgBasename, const configOptionsTrackingGaussianMixture &configOptions) {
    //to load hierarchcal segmentation
    vector<decltype(&make_hs_output_filename1)> filenameFormatters{
        // (ngc) I foolishly decided to fiddle with hs naming at some point
        //       and now I want to be backwards compatible.
        make_hs_output_filename1,
        make_hs_output_filename2,
		make_hs_output_filename3
	};
    ifstream inHS;
    string imgHS;
    for(auto formatter: filenameFormatters) {
        imgHS = formatter(imgBasename, configOptions);
        inHS.open(imgHS.c_str(), ios::binary | ios::in);
        if (inHS.is_open())
            break;
    }
    if(!inHS.is_open()){
        cout<<"ERROR: could not open binary file "<<imgHS<<" to read hierarchical segmentation"<<endl;
        exit (4);
    }
    auto hs = new hierarchicalSegmentation(inHS);
    //TODO: parse directly segmentation to supervoxels
    inHS.close();
    return hs;
}

void ProcessDeathThrOpticalFlow(
    configOptionsTrackingGaussianMixture &configOptions,
    int &numDeaths,
    std::vector<GaussianMixtureModel *> &vecGM,
    int frame,
    const std::string &debugPath,
    float scaleOrig[3],
    mylib::Array * imgFlowMask,
    int imgDims, mylib::Dimn_Type  imgSize[3], mylib::uint8 * imgFlowMaskPtr)
{
    if (configOptions.estimateOpticalFlow == 2)//it means we already tried to apply optical flow in this frame
    {
        configOptions.estimateOpticalFlow = 0;//we just set to zeroto avoid infinite loop in case optical flow does not fix death issue
    }
    else{
        numDeaths = 0;
        list<int> parentDeathList;
        for (unsigned int kk = 0; kk<vecGM.size(); kk++)
        {
            if (vecGM[kk]->isDead() == true)
            {
                numDeaths++;
                parentDeathList.push_back(vecGM[kk]->parentId);
            }
        }
        //activate module to calculate optical flow and redo this last frame
        if (numDeaths >= configOptions.deathThrOpticalFlow)
        {
            cout << "Number of dead cells in this frame " << numDeaths << " is above threshold. We need to rerun with optical flow" << endl;
            frame--;

            //delete existing solution
            for (unsigned int ii = 0; ii<vecGM.size(); ii++) delete vecGM[ii];
            vecGM.clear();

            //reload previous solution
            char extra[128];
            sprintf(extra, "%.4d", frame);
            string itoaD = string(extra);
            string GMxmlFilename(debugPath + "XML_finalResult/GMEMfinalResult_frame" + itoaD + ".xml");
            XMLNode xMainNode = XMLNode::openFileHelper(GMxmlFilename.c_str(), "document");
            int n = xMainNode.nChildNode("GaussianMixtureModel");
            for (int ii = 0; ii<n; ii++) vecGM.push_back(new GaussianMixtureModel(xMainNode, ii));

            memcpy(scaleOrig, GaussianMixtureModel::scale, sizeof(float)*dimsImage);//to save in case we work with units
#ifdef WORK_WITH_UNITS
                                                                                    //transform all given positions
            for (unsigned int ii = 0; ii<vecGM.size(); ii++) vecGM[ii]->pixels2units(scaleOrig);
#endif

            for (unsigned int kk = 0; kk<vecGM.size(); kk++)
            {
                vecGM[kk]->updateNk();//update Nk just in case
            }

            //setup imgFlowMask to select which areas we need flow (around death cells)
            if (imgFlowMask == NULL)//reallocate memory
            {
                imgFlowMask = mylib::Make_Array(mylib::PLAIN_KIND, mylib::UINT8_TYPE, imgDims, imgSize);
                imgFlowMaskPtr = (mylib::uint8*)(imgFlowMask->data);
            }
            //reset mask
            memset(imgFlowMask->data, 0, sizeof(mylib::uint8)*imgFlowMask->size);


            //---------------set imgFlowMask to 1 everywhere so we apply flow everywhere---------------------
            for (long long int aa = 0; aa<imgFlowMask->size; aa++) imgFlowMaskPtr[aa] = 1;

            //activate flag so we compute optical flow in the next iteration of the for loop
            configOptions.estimateOpticalFlow = 2;

        }
        else{
            //deactivate flag to compute optical flow
            configOptions.estimateOpticalFlow = 0;
        }
    }
    throw runtime_error("ERROR: with the update to lineage hyper tree the code is not ready yet to use vector flow calculations, sicne we need to update lht instead of vecGM");
}


static void PrepareWellFormedArray(mylib::Array *img) {
    //hack to make the code work for uin8 without changing everything to templates
    //basically, parse the image to uint16, since the code was designed for uint16
    if(img->type==mylib::UINT8_TYPE){
        img=mylib::Convert_Array_Inplace(img,img->kind,mylib::UINT16_TYPE,16,0);
    }
    //hack to make the code work for 2D without changing everything to templates
    //basically, add one black slice to the image (you should select conn3D = 4 or 8)
    if(img->ndims==2){
        mylib::Dimn_Type dimsAux[dimsImage];
        for(int ii=0; ii<img->ndims; ii++)
            dimsAux[ii]=img->dims[ii];
        for(int ii=img->ndims; ii<dimsImage; ii++)
            dimsAux[ii]=2;

        mylib::Array *imgAux=mylib::Make_Array(img->kind,img->type,dimsImage,dimsAux);
        memset(imgAux->data,0,(imgAux->size) * sizeof(mylib::uint16));
        memcpy(imgAux->data,img->data,img->size * sizeof(mylib::uint16));

        mylib::Array* imgSwap=img;
        img=imgAux;
        mylib::Free_Array(imgSwap);
    }
    if(img->ndims!=dimsImage){
        cout<<"ERROR: constant dimsImage does not match img->ndims. Please change dimsImage and recompile the code"<<endl;
        exit(3);
    }
    if(img->type!=mylib::UINT16_TYPE){
        cout<<"ERROR: code not ready for images that are not uint16"<<endl;
        exit(7);
    }
}

static void ensureConsistentDimensions(mylib::Array * img){
    for(int ii=0; ii < img->ndims; ii++){
        if(img->dims[ii]!=supervoxel::dataDims[ii]){
            if(img->dims[0]==supervoxel::dataDims[1]&&img->dims[1]==supervoxel::dataDims[0]&&img->dims[2]==supervoxel::dataDims[2]){
                cout<<"WARNING:"<<supervoxel::dataDims[0]<<"x"<<supervoxel::dataDims[1]<<"x"<<supervoxel::dataDims[2]<<endl;
                cout<<"WARNING:"<<img->dims[0]<<"x"<<img->dims[1]<<"x"<<img->dims[2]<<endl;
                cout<<"WARNING: transposing image"<<endl;
                transposeStackUINT16(img);
                break;
            } else{

                cout<<"ERROR!: data dimensions do not match between image and supervoxel record. Most likely a flipped x and y due to tiff to jp2 format with Matlab"<<endl;
                cout<<"ERROR:"<<supervoxel::dataDims[0]<<"x"<<supervoxel::dataDims[1]<<"x"<<supervoxel::dataDims[2]<<endl;
                cout<<"ERROR:"<<img->dims[0]<<"x"<<img->dims[1]<<"x"<<img->dims[2]<<endl;
                exit(3);
            }
        }
    }
}

static mylib::Array* MaybeEstimateOpticalFlow(const configOptionsTrackingGaussianMixture& configOptions,const mylib::Array *img,int frame,const string &imgBasename,string itoa) {
    //---------------------------------------------------------------
    //estimate optical flow
    mylib::Array* imgFlow=0;

    if (configOptions.estimateOpticalFlow == 2) {
        //calculate flow on the fly
        //parameters for optical are set inside that routine. They should be the same for all datasets
        cout << "Calculating optical flow. MAKE SURE ALL PARAMETERS WERE SET APPROPIATELY FOR THIS DATASET." << endl;


        string filenameInTarget(configOptions.imgFilePattern);
        parseImagePath(filenameInTarget, frame - 1);
        string filenameInSource(configOptions.imgFilePattern);
        parseImagePath(filenameInSource, frame);


        string filenameInTargetMask(filenameInTarget + "backgroundPredictionIlastik_");//using Ilastik predictions

                                                                                       //int auxMask=(int)(configOptions.thrSignal*255.0);
        int auxMask = 1;//current code uses supervoxels (no more Ilastik background)
        if (auxMask < 0 || auxMask > 255) {
            cout << "ERROR: thrSignal ooutside bounds" << endl;
            exit(2);
        }
        mylib::uint8 maskThr = (mylib::uint8)(auxMask);

        int err = 0;
        if (frame - 1 < 0)//initial frame: we can not calculate flow
        {
            mylib::Dimn_Type* ndims = new mylib::Dimn_Type[dimsImage + 1];
            for (int ii = 0; ii < dimsImage; ii++) ndims[ii] = img->dims[ii];
            ndims[dimsImage] = dimsImage;
            imgFlow = mylib::Make_Array(mylib::PLAIN_KIND, mylib::FLOAT32_TYPE, dimsImage + 1, ndims);
            memset(imgFlow->data, 0, imgFlow->size * sizeof(mylib::float32));
            delete[] ndims;
        }
        else
            //err=opticalFlow_anisotropy5(filenameInSource,filenameInTarget,filenameInTargetMask,maskThr,imgFlow,configOptions.maxDistPartitionNeigh);
            cout << "WARNING: Option not available in the current code. You need to implement your own optical flow code (or can request ours via email) and paste the function call here" << endl;
        if (err > 0)
            exit( err );
    }
    else if (configOptions.estimateOpticalFlow == 1) {//load precalculated flow

        char bufferDmax[64];
        sprintf(bufferDmax, "%d", int(configOptions.maxDistPartitionNeigh));
        string itoaDmax(bufferDmax);

        string filenameFlow((imgBasename + "flowArray_dMax" + itoaDmax + "_" + itoa + ".bin").c_str());
        ifstream inFlow(filenameFlow.c_str(), ios::binary | ios::in);
        if (!inFlow.is_open()) {
            cout << "ERROR: opening flow array file : " << filenameFlow << endl;
            exit( 2 );
        }
        else {
            cout << "Using file " << filenameFlow << " for flow information" << endl;
        }
        mylib::Dimn_Type* ndims = new mylib::Dimn_Type[dimsImage + 1];
        for (int ii = 0; ii < dimsImage; ii++) ndims[ii] = img->dims[ii];
        ndims[dimsImage] = dimsImage;
        imgFlow = mylib::Make_Array(mylib::PLAIN_KIND, mylib::FLOAT32_TYPE, dimsImage + 1, ndims);
        inFlow.read((char*)(imgFlow->data), imgFlow->size);
        inFlow.close();
        delete[] ndims;
    }
    return imgFlow;
}

static void mkdir(const string path) {
    if (path_exists(path))//check if folder exists
    {
        auto error = mkDirectory(string(path));
        if (error > 0)
            exit(error);
    }
}

void WriteOutputXML(string debugPath,string subdir,string name,int frame,function<void(ofstream& outXML)> writer) {
    mkdir(debugPath + subdir);
    ofstream outXML((debugPath+subdir+"/"+name+to_string(frame,4,'0')+".xml").c_str());
    GaussianMixtureModel::writeXMLheader(outXML);
    writer(outXML);
    GaussianMixtureModel::writeXMLfooter(outXML);
    outXML.close();
}

void DebugEMIteration(const string & debugPath, const vector<GaussianMixtureModel*> &vecGM, int frame) {
    int frameToPrint=40000;//print out all the frames
    stringstream itoaFrame;
    itoaFrame<<frame;
    if(frame==frameToPrint||frameToPrint<0){
        WriteOutputXML(debugPath,"XML_GMEMiterations_frame","debugEMGM_iter0000_",frame,[=](ofstream& outXML){
            for(unsigned int ii=0;ii<vecGM.size();ii++){
#ifdef WORK_WITH_UNITS
                vecGM[ii]->units2pixels(scaleOrig);
                vecGM[ii]->writeXML(outXML);
                vecGM[ii]->pixels2units(scaleOrig);
#else
                vecGM[ii]->writeXML(outXML);
#endif
            }
        });
    }
}

void DebugTGMMXmlFiles(const string & debugPath, const vector<GaussianMixtureModel*> & vecGM, int frame,const string &subdir,const string &name) {
    WriteOutputXML(debugPath,subdir,name,frame,[&](ofstream& outXML){
        for(unsigned int ii=0;ii<vecGM.size();ii++){
#ifdef WORK_WITH_UNITS
            vecGM[ii]->units2pixels(scaleOrig);
            vecGM[ii]->writeXML(outXML);
            vecGM[ii]->pixels2units(scaleOrig);
#else
            vecGM[ii]->writeXML(outXML);
#endif
        }
    });
}

void WriteGMEMResultXML(const string & debugPath, const vector<GaussianMixtureModel*> & vecGM, int frameOffset,const lineageHyperTree &lht) {
    WriteOutputXML(debugPath,"XML_finalResult_lht","GMEMfinalResult_frame",frameOffset,[&](ofstream& outXML)
    {
        int countV = 0;
        for (list<nucleus>::iterator iterN = lht.nucleiList[frameOffset].begin(); iterN != lht.nucleiList[frameOffset].end(); ++iterN , ++countV) {
            //modify parentId
            if (iterN->treeNodePtr->parent == NULL)
                vecGM[countV]->parentId = -1;
            else {
                vecGM[countV]->parentId = (int)(iterN->treeNodePtr->parent->data->tempWilcard);
            }
            //write solution
            vecGM[countV]->writeXML(outXML);
        }
    });
}

struct DriftCorrection {
    DriftCorrection(const configOptionsTrackingGaussianMixture & configOptions,int endFrame){
        driftCorrectionVec.reserve( endFrame);
        if( strcmp(configOptions.offsetFile.c_str(),"empty") != 0 )//read file
        {
            ifstream fin( configOptions.offsetFile.c_str() );

            if( fin.is_open() == false )
            {
                cout<<"ERROR: reading offset file "<<configOptions.offsetFile<<endl;
                throw runtime_error("Error reading drift correction file");
            }else{
                cout<<"Reading offset file "<<configOptions.offsetFile<<endl;
            }

            string line;
            OffsetCoordinates auxOC;
            int tt;
            while ( getline (fin,line) )
            {
                std::istringstream(line) >> tt >> auxOC.x >> auxOC.y >> auxOC.z;
                driftCorrectionVec.resize( tt+1);
                driftCorrectionVec[tt] = auxOC;
            }

            fin.close();
        }
    }

    void apply(const vector<GaussianMixtureModel*> &vecGM, int frame) {
        if( driftCorrectionVec.size() > frame )
        {
            for(unsigned int kk=0;kk<vecGM.size();kk++)
            {
                vecGM[kk]->m_k[0] += driftCorrectionVec[frame].x;
                vecGM[kk]->m_k[1] += driftCorrectionVec[frame].y;
                vecGM[kk]->m_k[2] += driftCorrectionVec[frame].z;

                vecGM[kk]->m_o[0] += driftCorrectionVec[frame].x;
                vecGM[kk]->m_o[1] += driftCorrectionVec[frame].y;
                vecGM[kk]->m_o[2] += driftCorrectionVec[frame].z;
            }
        }
    }
  private:
      vector<OffsetCoordinates> driftCorrectionVec;
};

void ParseHStoLineageHyperTree(lineageHyperTree &lht, int frame, const mylib::Array* img, hierarchicalSegmentation* hs,
                              long long query_nb, unsigned numLabels, float* queryHOST, float* imgDataHOST,
                              float* centroidLabelPositionHOST, float* centroidLabelWeightHOST,
                              long long* labelListPtrHOST) {
    long long int idxQuery = 0;
    long long int offset = 0, offsetL = 0;
    mylib::Coordinate *coord;
    mylib::Dimn_Type *ccAux;
    list<supervoxel> *listPtr = (&(lht.supervoxelsList[frame]));//to avoid refering to it all the time
    mylib::float32* imgData=(float*)img->data;
    listPtr->clear();

    labelListPtrHOST[0]=0;
    for(int ii = 0; ii < hs->currentSegmentatioSupervoxel.size(); ii++)
    {
        for(vector<uint64>::const_iterator iter =  hs->currentSegmentatioSupervoxel[ii].PixelIdxList.begin();
            iter != hs->currentSegmentatioSupervoxel[ii].PixelIdxList.end();
            ++iter)
        {
            imgDataHOST[idxQuery] = imgData[*iter];
            //find closest maxGaussiansPerVoxel Gaussians
            coord=mylib::Idx2CoordA((mylib::Array*)(img),*iter); // not declared as const Array, but it really is
            ccAux=(mylib::Dimn_Type*)coord->data;
            offset = idxQuery;
            for(int ii=0;ii<dimsImage;ii++)
            {
#ifdef WORK_WITH_UNITS
					queryHOSTaux[offset]=scaleOrig[ii]*ccAux[ii];//transform index 2 coordinates
#else
                queryHOST[offset] = (float)ccAux[ii];//transform index 2 coordinates
#endif
                offset += query_nb;//interleave memory to favor coalescenc access in GPU
            }
            idxQuery++;
            mylib::Free_Array(coord);
        }

        //set arrays related to each label
        centroidLabelWeightHOST[ii] = hs->currentSegmentatioSupervoxel[ii].intensity;
        labelListPtrHOST[ ii + 1 ] = labelListPtrHOST[ ii ] + hs->currentSegmentatioSupervoxel[ii].PixelIdxList.size();
        offsetL = ii;
        for(auto jj=0;jj<dimsImage;jj++)
        {
            centroidLabelPositionHOST[offsetL] = hs->currentSegmentatioSupervoxel[ii].centroid[jj];
            offsetL+=numLabels;
        }

        //update lineageHyperTree supervoxel (this is duplication needed for "backwards" compatibility)
        listPtr->push_back( hs->currentSegmentatioSupervoxel[ii] );
    }
}



void DebugWriteRNK_FirstRoundEMnoSplit(const string &debugPath, int frame, string& debugPathCUDAa) {
    char buffer[128];
    debugPathCUDAa = "";
    sprintf(buffer,"%.4d",frame);//to save rnk and ind for each supervoxel
    string itoa(buffer);
    debugPathCUDAa = string(debugPath+"XML_finalFirstRoundEMnoSplit");
    mkdir(debugPathCUDAa);
    debugPathCUDAa=string(debugPathCUDAa+ "/rnk_frame"+ itoa + ".bin");
}

void UpdateLineageHyperTree(bool regularize_W4DOF, int iniFrame, lineageHyperTree &lht,const  vector<GaussianMixtureModel*> & vecGM,
                            float scaleOrig[3], int frame,float *imgPtr) {
    list<nucleus> *listNucleiPtr = (&( lht.nucleiList[frame] ) );
    list<nucleus>::iterator listNucleiIter;
    nucleus nucleusAux;
    nucleusAux.TM = frame;


    vector< list<supervoxel>::iterator > vecSupervoxelsIter;
    lht.getSupervoxelListIteratorsAtTM(vecSupervoxelsIter, frame);

    vector< list<nucleus>::iterator > vecNucleiIter;
    if(frame > iniFrame)
        lht.getNucleiListIteratorsAtTM(vecNucleiIter, frame-1);//we need this to retrieve lineage

    for(size_t ii = 0; ii<vecGM.size(); ii++)
    {
        if(vecGM[ii]->isDead() == true) continue;
        if( vecGM[ii]->supervoxelIdx.empty() ) continue;//there was no clear assignment for this nuclei


        /*
        if( vecGM[ii]->supervoxelIdx.size() == 1 )//TOD0: sometimes, there is larger discrepancy between vecGM centorid and supervoxel centroid (even when there is only one super-voxel). this is a quick hack
        {
            for(int aa = 0;aa<dimsImage;aa++)
                nucleusAux.centroid[aa] = vecSupervoxelsIter[ vecGM[ii]->supervoxelIdx[0] ]->centroid[aa];
        }else{
            for(int aa = 0;aa<dimsImage;aa++)
                nucleusAux.centroid[aa] = vecGM[ii]->m_k(aa);
        }
        */

        nucleusAux.avgIntensity = vecGM[ii]->N_k;

        listNucleiPtr->push_back(nucleusAux);
        listNucleiIter = ((++ ( listNucleiPtr->rbegin() ) ).base());//iterator for the last element in the list

        //add supervoxel-nuclei relationship
        for(size_t aa = 0; aa < vecGM[ii]->supervoxelIdx.size();aa++)
        {
            listNucleiIter->addSupervoxelToNucleus(vecSupervoxelsIter[ vecGM[ii]->supervoxelIdx[aa] ]);
            vecSupervoxelsIter[ vecGM[ii]->supervoxelIdx[aa] ]->treeNode.setParent(listNucleiIter);
        }


        //add lineaging relationships
        if(vecGM[ii]->parentId >= 0 ) //parentId contains the position of the parent nucleus
        {
            listNucleiIter->treeNode.setParent ( vecNucleiIter[ vecGM[ii]->parentId ] -> treeNode.getParent() );
            listNucleiIter->treeNode.getParent()->bt.SetCurrent( vecNucleiIter[ vecGM[ii]->parentId ] ->treeNodePtr );
            listNucleiIter->treeNodePtr = listNucleiIter->treeNode.getParent()->bt.insert(listNucleiIter);
            if(listNucleiIter->treeNodePtr == NULL)
            	exit(3);//error


        }else{//new lineage
            lht.lineagesList.push_back( lineage() );
            list<lineage>::iterator listLineageIter = ((++ ( lht.lineagesList.rbegin() ) ).base());//iterator for the last element in the list

            listNucleiIter->treeNode.setParent(listLineageIter);
            listNucleiIter->treeNodePtr = listLineageIter->bt.insert(listNucleiIter);
            if(listNucleiIter->treeNodePtr == NULL)
            	exit(3);//error
        }

        //calculate centroid and precision matrix
        lht.calculateNucleiGMMParameters<float>(listNucleiIter, regularize_W4DOF, scaleOrig,imgPtr);
    }
}

static void CreateDebugFolder(int argc, const char** argv, configOptionsTrackingGaussianMixture & configOptions, int iniFrame, const string & debugPath) {
    string cmd;

    auto error = mkDirectory(debugPath);
    if (error > 0)
    	exit( error );


    //open file to record parameters
    char logBuffer[32];
    sprintf(logBuffer,"%.4d",iniFrame);//we always read the same image. We could make the code faster by reading it only once but I jts want to test the idea
    string itoaLog=string(logBuffer);
    ofstream outLog(debugPath + "experimentLog_"+ itoaLog + ".txt");
    outLog<<"Program called with the following command:"<<endl;
    for (int ii=0;ii<argc;ii++) outLog<<argv[ii]<<" ";
    outLog<<endl;
    outLog<<"minPi_kForDead="<<minPi_kForDead<<endl;
    outLog<<"maxGaussiansPerVoxel="<<maxGaussiansPerVoxel<<endl;
    configOptions.printConfigFileTrackingGaussianMixture(outLog);
}

static int g_gpuid=0;

static void parse_options(int* argc,const char** argv) {
    for(int i=0;i<*argc;++i) { // scan args for option flags
        // Process: --gpu <id>
        if(strcmp(argv[i],"--gpu")==0) {
            if((i+1)>=*argc)
                throw runtime_error("The --gpu <id> option is missing the <id> field.");
			char *end=0;
			g_gpuid=strtol(argv[i+1],&end,10);
			if(end==argv[i+1])
				throw runtime_error("When processing the --gpu <id>, the <id> field did not seem to be an integer.");
            *argc-=2;
        }
    }
}


struct stat_time_acc_t{ float t,n; stat_time_acc_t():t(0),n(0){} void add(float x){ t+=x;n++; } float avg() const{ return t/n; } };
#define STATTIME_INIT map<string,stat_time_acc_t> stattimes
#define STATTIME(expr) do{ TicTocTimer t=tic(); (expr); stattimes[#expr].add(toc(&t)); }while(0)
#define STATTIME_READ do{ for(auto const &it: stattimes) {cout << "\tTotal Elapsed: " << it.second.t << "s for " << it.first << endl ;} }while(0)

void recalculateThrDist2LargeDisplacement(int iniFrame, lineageHyperTree& lht, float& thrDist2LargeDisplacement, int frame)
{
	if (frame > iniFrame)
	{
		STATTIME_INIT;
		vector<float> dispVec;
		dispVec.reserve(lht.nucleiList[frame].size());

		float auxDist;
		for (list<nucleus>::iterator iterNN = lht.nucleiList[frame].begin(); iterNN != lht.nucleiList[frame].end(); ++iterNN)
		{
			if (iterNN->treeNodePtr->parent != NULL)
			{
				STATTIME(auxDist = iterNN->treeNodePtr->data->Euclidean2Distance(*(iterNN->treeNodePtr->parent->data), supervoxel::getScale()));
				STATTIME(dispVec.push_back(sqrt(auxDist)));


				//if displacement is way too large->flag it as background (it might not be, but we want to clean these results anyways)
				if (auxDist > 6.25 * thrDist2LargeDisplacement)//because we use squares, if we want k times above thrDist2LargeDisplacement->k*k  (for us, k = 2.5 -> 10 std are just not acceptable)
					iterNN->probBackground = 1.0f;
			}
		}		
		STATTIME_READ;
		cout<<"\trecalculateThrDist2LargeDisplacement: Displacement count: "<<dispVec.size()<<endl;
		if (dispVec.size() > 30)//for reliable statistics. Otherwise we do not update
		{
			//calculate median
			TIME(std::sort(dispVec.begin(), dispVec.end()));
			float auxMedian = dispVec[(size_t)(dispVec.size() / 2)];
			//calculate median absolute deviation
			for (size_t ii = 0; ii < dispVec.size(); ii++)
			{
				dispVec[ii] -= auxMedian;
				dispVec[ii] = fabs(dispVec[ii]);
			}
			TIME(std::sort(dispVec.begin(), dispVec.end()));
			float auxMAD = dispVec[(size_t)(dispVec.size() / 2)];
			thrDist2LargeDisplacement = auxMedian + 4.0 *  auxMAD * 1.4826;//sigma = 1.4826 * MAD -> 4 std are considered outliers
			thrDist2LargeDisplacement = thrDist2LargeDisplacement * thrDist2LargeDisplacement;
		}
	}
	cout << "Updated thrDist2LargeDisplacement to = " << sqrt(thrDist2LargeDisplacement) << " * " << sqrt(thrDist2LargeDisplacement) << endl;
}

void recalculateCentroidForMotherAndTwoDaughters(bool regularize_W4DOF, lineageHyperTree& lht, float scaleOrig[3], TimeSeriesMap& time_series_map, int TMaux)
{
	cout<<"\trecalculateCentroidForMotherAndTwoDaughters"<<endl
		<<"\t\t#nuclei: "<<lht.nucleiList[TMaux].size()<<endl;
	int counter=0;
	
	STATTIME_INIT;
	for (list<nucleus>::iterator iterN = lht.nucleiList[TMaux].begin(); iterN != lht.nucleiList[TMaux].end(); ++iterN)
	{
		TreeNode<ChildrenTypeLineage>* aux = iterN->treeNodePtr;

		if (aux->getNumChildren() != 2)
			continue;//not a cell division

		//recalculate centroid for mother and two daughters
		float* im;
		STATTIME(im=time_series_map.getProcessed(iterN));
		STATTIME(lht.calculateNucleiGMMParameters<float>(iterN,regularize_W4DOF,scaleOrig,im));
		STATTIME(im=time_series_map.getProcessed(iterN->treeNodePtr->left->data));
		STATTIME(lht.calculateNucleiGMMParameters<float>(iterN->treeNodePtr->left->data,regularize_W4DOF,scaleOrig,im));
		STATTIME(im=time_series_map.getProcessed(iterN->treeNodePtr->right->data));
		STATTIME(lht.calculateNucleiGMMParameters<float>(iterN->treeNodePtr->right->data,regularize_W4DOF,scaleOrig,im));
		counter++;
	}
	STATTIME_READ;
	cout<<"\t\t#recalculated: "<<counter<<endl;
}

void saveSupervoxels(const string & debugPath, const lineageHyperTree & lht, int frameOffset)
{
	string filenameSv(debugPath + "XML_finalResult_lht/GMEMfinalResult_frame" + to_string(frameOffset, 4, '0') + ".svb");
	cout << "Writing: " << filenameSv << endl;
	if (lht.writeListSupervoxelsToBinaryFile(filenameSv, frameOffset) > 0)
		throw 3;
}


//   ____    ____        __        _____   ____  _____
//  |_   \  /   _|      /  \      |_   _| |_   \|_   _|
//    |   \/   |       / /\ \       | |     |   \ | |
//    | |\  /| |      / ____ \      | |     | |\ \| |
//   _| |_\/_| |_   _/ /    \ \_   _| |_   _| |_\   |_
//  |_____||_____| |____|  |____| |_____| |_____|\____|
//

int main( int argc, const char** argv )
{
    try {
        cout << "FTGMM Version " << GIT_TAG << " " << GIT_HASH << endl;

        parse_options(&argc,argv); // (ngc) options have to be after all the other args.  This modifies argc so things work with the older code below.

    //------------------------------parameters that should be set by user later on (Advanced mode)-----------------------
        bool regularize_W4DOF = true;//whether we allow rotation in Z for Gaussian covariance or not. If regularize_W4DOF==true->W_02 = W_12 = 0

        string pathS = basename();
        string backgroundClassifierFilename(pathS + string("/") + "classifierBackgroundTracks.txt");
        int temporalWindowSizeForBackground = 5;//size of the temporal window to calculate temporal features for background detection. This has to match the classifier!!!

        //--------------------------------parse input arguments from config file---------------------------
        if (argc != 4 && argc != 5)
        {
			cout<<"ERROR: wrong number of input arguments"<<endl
				<<"Usage: TGMM <configFile> <start> <end> [options]"<<endl
				<<"Options:"<<endl
				<<"\t--gpu <id>\tSelect GPU <id> for processing. Default GPU is 0."<<endl;
            exit(2);
        }
        if (argc >= 5)
        {
            cout << "ERROR: with new code using sliding window the option of starting after a crash is not coded yet" << endl;
            exit(3);
        }

        //set main parameters for the algorithm
        configOptionsTrackingGaussianMixture configOptions;
        string configFilename(argv[1]);
        if (configOptions.parseConfigFileTrackingGaussianMixture(configFilename) != 0) 
			exit(2);

        regularizePrecisionMatrixConstants::setConstants(configOptions.lambdaMin, configOptions.lambdaMax, configOptions.maxExcentricity);

        //trim supervoxels
        int minNucleiSize = configOptions.minNucleiSize;//in voxels
        int maxNucleiSize = configOptions.maxNucleiSize;
        imgVoxelType tauMax = 300;//to trim hierarchical segmentation
        float maxPercentileTrimSV = configOptions.maxPercentileTrimSV;//percentile to trim supervoxel
        int conn3Dsv = configOptions.conn3DsvTrim;//connectivity to trim supervoxel

        //for logical temporal rules
        const int lengthTMthr = configOptions.SLD_lengthTMthr;//for delete short lived cell (minimum length of daughter before ending track)
        const int minNeighboringVoxels = 10;//if two supervoxels are have less than minNeighboringVoxels as common border using minNeighboringConn3D, then they must belong to a different nucleus
        const int minNeighboringVoxelsConn3D = 4;

        //for merge split
        const int mergeSplitDeltaZ = 11;
        int64 boundarySizeIsNeigh[dimsImage];
        int conn3DIsNeigh = 74;


        int devCUDA;
        {
            devCUDA=g_gpuid;
            //cout << "Overwriting CUDA device selction to default" << endl;
            if (devCUDA < 0)
            {
                exit(2);
            }
            else{
                cout<<"Selected CUDA device "<<devCUDA<<endl;
            }
        }


        //set init and end frame
        int iniFrame = 0, endFrame = 0;
        //parse input arguments
        iniFrame = atoi(argv[2]);
        endFrame = atoi(argv[3]);

        //read last entry to continue from that frame
        string debugPath;
        char extra[128];
        if (argc >= 5)
        {
            debugPath = string(argv[4]);//we read previous xml file as initialization
            sprintf(extra, "%.4d", iniFrame - 1);//we always read the same image. We could make the code faster by reading it only once but I jts want to test the idea
            string itoaD = string(extra);
            //configOptions.GMxmlIniFilename=string(debugPath+"XML_finalResult/GMEMfinalResult_frame"+ itoaD + ".xml");
            //cout<<"Starting from previous run. Make sure parameters are the same!!!. XMLini="<<configOptions.GMxmlIniFilename<<endl;
            cout << "ERROR: Starting from previous run HAS NOT BEEN IMPLEMENTED YET" << endl;
            return 3;
        }
        else{

#if defined(_WIN32) || defined(_WIN64)
            SYSTEMTIME str_t;
            GetSystemTime(&str_t);
            //sprintf(extra,"C:\\Users\\Fernando\\TrackingNuclei\\TGMMruns\\GMEMtracking3D_%d_%d_%d_%d_%d_%d\\",str_t.wYear,str_t.wMonth,str_t.wDay,str_t.wHour,str_t.wMinute,str_t.wSecond);
            sprintf(extra, "%s/GMEMtracking3D_%d_%d_%d_%d_%d_%d/", configOptions.debugPathPrefix.c_str(), str_t.wYear, str_t.wMonth, str_t.wDay, str_t.wHour, str_t.wMinute, str_t.wSecond);
#else
            sprintf(extra, "%s/GMEMtracking3D_%ld/", configOptions.debugPathPrefix.c_str(), time(NULL));
#endif
            debugPath = string(extra);
        }


        //---------------------------------end of parse input arguments--------------------------------------------
        CreateDebugFolder(argc, argv, configOptions, iniFrame, debugPath);

        //creat elineage hypertree to store lineages
        lineageHyperTree lht(endFrame + 1);
        float thrDist2LargeDisplacement = 2048.0f * 2048.0f;//to set low confidence level if a cell displaces more than this threshold (in pixels with scale). This is adapted at each time point from the data itself.

        DriftCorrection driftCorrection(configOptions, endFrame);  //read file with offsets between time points

        //----------------------------------main for loop for sequential processing--------------------------
        char buffer[128];
        //string itoa,itoa2;
        vector<GaussianMixtureModel*> vecGM;
        vecGM.reserve(10000);
        vector< pair<GaussianMixtureModel*, GaussianMixtureModel*> > splitSet;//stores candidates to be split
        splitSet.reserve(1000);
        vector<double> singleCellSplitScore;//store split scores without divisions to be able to compare afterwards
        singleCellSplitScore.reserve(1000);
        vector<GaussianMixtureModel> backupVecGM;
        backupVecGM.reserve(1000);
        vector<feature> Fx;
        vector<double> splitMergeTest;
        Fx.reserve(10000);
        splitMergeTest.reserve(10000);

#ifdef CELL_DIVISION_WITH_GPU_3DHAAR
        CellDivisionWithGPU3DHaar cellDivisionWithGPU3DHaar(pathS);
#endif



        //load classifier for cell division using temporal window
		DivisionClassifierFactory division_classifier_factory(pathS,configOptions.cellDivisionClassifier,configOptions.temporalWindowRadiusForLogicalRules);
		const auto& cdwtClassifier=division_classifier_factory.get();

        float scaleOrig[dimsImage];//to save original scale
        mylib::Dimn_Type imgSize[dimsImage];//to save size of image
        int imgDims = dimsImage;//to save number of dimensions
        for (int ii = 0; ii < dimsImage - 1; ii++)
            scaleOrig[ii] = 1.0f;
        scaleOrig[dimsImage - 1] = configOptions.anisotropyZ;
        GaussianMixtureModel::setScale(scaleOrig);

        //int numDiv = 10;//number of divisions per frame

        //optical flow pointers
        mylib::Array* imgFlow = NULL;
        mylib::Array* imgFlowMask = NULL;//uint8 containing areas where we need flow estimation
        mylib::uint8* imgFlowMaskPtr = NULL;


        vector<hierarchicalSegmentation*> hsVec;//stores hierarchical segmentation for frames in the temporal window
        hsVec.reserve(configOptions.temporalWindowRadiusForLogicalRules*2);//it shoudl match teh temporal window size

        //vector<mylib::Array*> imgVecUINT16;//store original raw data in the temporal window for GPU Haar ceel division features(or any other thing we need in the future)
        //imgVecUINT16.reserve(2 * configOptions.temporalWindowRadiusForLogicalRules + 1);

		TimeSeriesMap time_series_map(&configOptions, iniFrame);
        for (int frame = iniFrame; frame <= endFrame; frame++){
            TicTocTimer tt = tic();
            cout << "Processing frame " << frame << endl;
            //==============================================================
            //read image
            sprintf(buffer, "%d", configOptions.tau);
            string itoaTau = string(buffer);
            string imgBasename(configOptions.imgFilePattern);
            parseImagePath(imgBasename, frame);



            mylib::Array* img = LoadTimepoint(imgBasename);
            MaybeMedianFilter(configOptions, devCUDA, img);
            PrepareWellFormedArray(img);  // adapt img to handle different dimensionality and types, by mapping to a single kind of array (3D,UINT16)
            hierarchicalSegmentation* hs = LoadHierarchicalSegmentation(imgBasename, configOptions);

            //associate each supervoxel with the correct data pointer

            ////save raw input data before doing any pre-processing with supervoxels
            ////copy image before converting it to float (we need it for image features in the GPU)
            //{
            //    mylib::Array* imgUINT16 = mylib::Make_Array(img->kind, img->type, img->ndims, img->dims);
            //    memcpy(imgUINT16->data, img->data, sizeof(mylib::uint16) * (img->size));
            //    if (imgVecUINT16.size() == 2 * configOptions.temporalWindowRadiusForLogicalRules + 1)//already full
            //    {
            //        mylib::Free_Array(imgVecUINT16.front());
            //        imgVecUINT16.erase(imgVecUINT16.begin());//even if this has to copy elements it is only 11 pointers
            //        imgVecUINT16.push_back(imgUINT16);
            //    }
            //    else{
            //        imgVecUINT16.push_back(imgUINT16);
            //    }
            //}
            time_series_map.push_back(img);

            ensureConsistentDimensions(img);  // checks image dimensions against dimensions set for supervoxel ... (ngc) not sure where the class supervoxel gets it's dims from

            {
                supervoxel::dataSizeInBytes = sizeof(mylib::uint16) * img->size;
                supervoxel::setScale(GaussianMixtureModel::scale);
                for (unsigned int ii = 0; ii < hs->getNumberOfBasicRegions(); ii++){
                    hs->basicRegionsVec[ii].TM = frame;
                }

                hsVec.push_back(hs);
                //setup trimming parameters
                supervoxel::setTrimParameters(maxNucleiSize, maxPercentileTrimSV, conn3Dsv, configOptions.svKeepWithoutTrimming);
                hs->setMaxTau(tauMax);

                //merge / split parameters for HS
                int64* neighOffsetIsNeigh = supervoxel::buildNeighboorhoodConnectivity(conn3DIsNeigh + 1, boundarySizeIsNeigh);//using the special neighborhood for coarse sampling
                supervoxel::pMergeSplit.setParam(conn3DIsNeigh, neighOffsetIsNeigh, mergeSplitDeltaZ);
                delete[] neighOffsetIsNeigh;
            }

            //save image dimensions // (ngc) why?
            for (int aa = 0; aa < img->ndims; aa++) imgSize[aa] = img->dims[aa];
            imgDims = img->ndims;

            // Supervoxel trimming, local thresholding and
            // Rescale to 0,1
            auto query_nb = time_series_map.process(frame, hs);
            HERE; PrintMemoryInUse();

            // Prep for the next part
            // update supervoxel informations: img->data pointer has changed after rescaling
            img = time_series_map.getProcessedArray(frame);
            if (img->type != 8)
                throw runtime_error("ERROR: code is only ready for FLOAT32 images");

            supervoxel::dataSizeInBytes = sizeof(float) * img->size;
            //for(unsigned int ii = 0;ii < hs->getNumberOfBasicRegions(); ii++)
            //	hs->basicRegionsVec[ii].dataPtr = img->data;

            //for(unsigned int ii = 0;ii < hs->currentSegmentatioSupervoxel.size(); ii++)
            //	hs->currentSegmentatioSupervoxel[ii].dataPtr = img->data;

            imgFlow = MaybeEstimateOpticalFlow(configOptions, img, frame, imgBasename, to_string(frame, 4, '0'));

            //------------------------------------------
            //estimate threshold of what is signal and what is not

            //initialize responsibilities
            unsigned int numLabels = hs->currentSegmentatioSupervoxel.size();//keeps track of the maximum number of labels that are assigned
            HANDLE_ERROR( cudaSetDevice( devCUDA ) );
            //host temporary memory needed
            float *queryHOST = new float[query_nb*dimsImage];
            float *imgDataHOST = new float[query_nb];
            float *centroidLabelPositionHOST = new float[numLabels*dimsImage];//stores the centroid for each supervoxel
            float *centroidLabelWeightHOST = new float[numLabels];//stores the wieght (based on intensity) for each supervoxel
            long long int *labelListPtrHOST = new long long int[numLabels + 1];//same concept as column-compressed sparse matrices. When we sort positions by region id (i.e. supervoxels), then all the elements for i-th region are bteween labelListPtr[ii]<=p<labelListPtr[ii+1]  ii=0,...,numLabels-1. Thus, labelsListPtr[0]=0 and labelaListPtr[numLabels]=query_nb

            memset(centroidLabelPositionHOST, 0, sizeof(float)*numLabels*dimsImage);//reset
            memset(centroidLabelWeightHOST, 0, sizeof(float)*numLabels);


            //initialize CUDA arrays to hold 3D points that represent signal
            float* queryCUDA = NULL;//holds 3D locations of each signal pixel in order to find k-NN
            float* imgDataCUDA = NULL;//holds intensity of each pixel          t
            float *rnkCUDA = NULL;
            float *rnkCUDAtr = NULL;//not really needed anymore
            int *indCUDA = NULL;
            int *indCUDAtr = NULL;
            float *centroidLabelPositionCUDA = NULL;
            long long int *labelListPtrCUDA = NULL;

            unsigned short int auxL = 0;

            //main loop to parse information from HS to lineageHyperTree and vecGM: remember supervoxels have been trimmed already
            ParseHStoLineageHyperTree(lht, frame, img, hs,
                query_nb, numLabels,
                queryHOST, imgDataHOST, centroidLabelPositionHOST, centroidLabelWeightHOST, labelListPtrHOST);


            cout << "Calculating nearest neighbors for supervoxels" << endl;
            TicTocTimer ttKNN = tic();
            int err;
            if (frame > 0)
            {
                err = lht.supervoxelNearestNeighborsInTimeForward(frame - 1, configOptions.KmaxNumNNsupervoxel, configOptions.KmaxDistKNNsupervoxel, devCUDA);
                if (err > 0) return err;
            }
            err = lht.supervoxelNearestNeighborsInTimeBackward(frame, configOptions.KmaxNumNNsupervoxel, configOptions.KmaxDistKNNsupervoxel, devCUDA);
            if (err > 0) return err;
            err = lht.supervoxelNearestNeighborsInSpace(frame, configOptions.KmaxNumNNsupervoxel, configOptions.KmaxDistKNNsupervoxel, devCUDA);
            if (err > 0) return err;

            cout << "Nearest neighbors took " << toc(&ttKNN) << " secs" << endl;
            //-----------------------------------------------------------------------

            //--------------------------------------------------------------
            //initialize Gaussians with priors from previous frame
            int numDeaths = 0;
            int Wsize = dimsImage * dimsImage;
            if (frame == iniFrame)//generate initial frame from supervoxel segmentation
            {
                //generate one nuclei per supervoxel
                int countS = 0;
                ParentTypeSupervoxel iterN;
                list<lineage>::iterator iterL;
                for (list<supervoxel>::iterator iterS = lht.supervoxelsList[frame].begin(); iterS != lht.supervoxelsList[frame].end(); ++iterS, ++countS)
                {
                    iterN = lht.addNucleusFromSupervoxel(frame, iterS);//returns pointer to the newly created supervoxel
                    //we do not need to create new lineages: the program will create them automatically after the TGMM iteration since parentId = -1;
                }

                //generate vecGM from lht
                TIME(parseNucleiList2TGMM<float>(vecGM, lht, frame, regularize_W4DOF, thrDist2LargeDisplacement, time_series_map.getProcessed(frame)));
                for (unsigned int kk = 0; kk < vecGM.size(); kk++)//set alpha_o to some default value since many of them can die due to oversegmentation
                {
                    vecGM[kk]->alpha_o = vecGM[kk]->N_k * 0.1;
                    vecGM[kk]->parentId = -1;//so we know it is the beginning of a new lineage
                }
                //update priors based on estimations: since this is the first frame it should be a very loose prior
                for (unsigned int kk = 0; kk < vecGM.size(); kk++)
                {
                    vecGM[kk]->updatePriors(0.1, 0.1, 0.1, -1.0);//use alphaTotal=-1 to not update alpha_o
                    //hack to make it similar to previous code without sliding window
                    vecGM[kk]->nu_k *= 2.0;
                    for (int aa = 0; aa < Wsize; aa++)
                        vecGM[kk]->W_k.data()[aa] /= 2.0;
                }
                //remove nuclei from lht since the final result will be added after GMM
                for (list<supervoxel>::iterator iterS = lht.supervoxelsList[frame].begin(); iterS != lht.supervoxelsList[frame].end(); ++iterS, ++countS)
                {
                    iterS->treeNode.deleteParent();
                }
                lht.nucleiList[frame].clear();

            }
            else//this is not the first time point
            {

                //generate vecGM from lht
                TIME(parseNucleiList2TGMM<float>(vecGM, lht, frame - 1, regularize_W4DOF, thrDist2LargeDisplacement,
                    time_series_map.getProcessed(frame - 1)));//we generate nuclei list from previous frame-> we want to extend this solution to t+1 using TGMM framework
                double alphaTotal = 0.0;
                for (unsigned int kk = 0; kk < vecGM.size(); kk++)
                {
                    alphaTotal += vecGM[kk]->alpha_k;
                }
                //cout<<imgFlow<<" "<<imgFlowMask<<endl;
                if (imgFlow != NULL)
                {
                    //TODO: I should use r_nk (responsibities) to calculate a weighted mean of the flow for each cell
                    int err = applyFlowPredictionToNuclei(imgFlow, imgFlowMask, vecGM, true);//we use forward flow
                    if (err > 0) exit(err);
                }
                //update priors based on estimations
                for (unsigned int kk = 0; kk < vecGM.size(); kk++)
                {
                    vecGM[kk]->updatePriors(configOptions.betaPercentageOfN_k, configOptions.nuPercentageOfN_k, configOptions.alphaPercentage, alphaTotal);//use alphaTotal=-1 to not update alpha_o
                    //hack to make it similar to previous code without sliding window
                    vecGM[kk]->nu_k *= 2.0;
                    for (int aa = 0; aa < Wsize; aa++)
                        vecGM[kk]->W_k.data()[aa] /= 2.0;
                }
            }


            //=====================================================
            driftCorrection.apply(vecGM, frame); //apply offset from txt file

            //--------------------------------------------------
            //---------------debug:write out iteration
#if DEBUG_TGMM_XML_FILES
            DebugTGMMXmlFiles(debugPath, vecGM, frame, "XML_KalmanFilterPrediction", "GMEMiniKalmanFilterPrediction_frame");
#endif
#ifdef DEBUG_EM_ITER
            DebugEMIteration(debugPath, St, error, vecGM, frame);
#endif

    //----------------------------------------------------------------------------------------------------------


        //allocate memory in device
        //GMEMinitializeMemory(&queryCUDA,queryHOST,&imgDataCUDA,imgDataHOST,GaussianMixtureModel::scale,query_nb,&rnkCUDA,&indCUDA,&indCUDAtr);
            GMEMinitializeMemoryWithSupervoxels(&queryCUDA, queryHOST, &imgDataCUDA, imgDataHOST, GaussianMixtureModel::scale, query_nb, &rnkCUDA, &indCUDA, &indCUDAtr, numLabels, &centroidLabelPositionCUDA, centroidLabelPositionHOST, &labelListPtrCUDA, labelListPtrHOST);
            copyScaleToDvice(GaussianMixtureModel::scale);

            //-----------------------------------------------------------------

            //==============================================================

            //string debugPathCUDAa(debugPath+"XML_GMEMiterations_CUDA_noSplit_frame");
            //string debugPathCUDAa("");

            //----------------to write rnk---------------------------------
            string debugPathCUDAa;
#ifdef DEBUG_TGMM_XML_FILES
            DebugWriteRNK_FirstRoundEMnoSplit(debugPath, frame, debugPathCUDAa);
#endif
        //----------------------------------------------------

            GaussianMixtureModelCUDA *vecGMHOST = new GaussianMixtureModelCUDA[vecGM.size()];
            for (unsigned int ii = 0; ii < vecGM.size(); ii++)
                copy2GMEM_CUDA(vecGM[ii], &(vecGMHOST[ii]));
            GMEMvariationalInferenceCUDAWithSupervoxels(
                queryCUDA, imgDataCUDA, rnkCUDA, rnkCUDAtr, indCUDA, indCUDAtr, centroidLabelPositionCUDA, labelListPtrCUDA,
                vecGMHOST,
                query_nb, vecGM.size(), numLabels,
                configOptions.maxIterEM, configOptions.tolLikelihood,
                devCUDA, frame, regularize_W4DOF, debugPathCUDAa);
            for (unsigned int ii = 0; ii < vecGM.size(); ii++)
                copyFromGMEM_CUDA(&(vecGMHOST[ii]), vecGM[ii]);

            size_t vecGM_ref_nb = vecGM.size();//I need to save it for later
            //delete[] vecGMHOST; //we need it to estimate local likelihood before split

            if ( configOptions.cellDivisionClassifier.Amatf2013.use_2021_code )
            {
                //2021 division detection is resource intensive; free all the host and device memory that we can
                GMEMreleaseMemoryWithSupervoxels(&queryCUDA, &imgDataCUDA, &rnkCUDA, &indCUDA, &indCUDAtr, &centroidLabelPositionCUDA, &labelListPtrCUDA);
                delete[]queryHOST;
                delete[]imgDataHOST;
                delete[]centroidLabelPositionHOST;
                delete[]centroidLabelWeightHOST;
                delete[]labelListPtrHOST;
			}
            //----------------debug--------------------------------------------



#ifdef DEBUG_TGMM_XML_FILES
            DebugTGMMXmlFiles(debugPath, vecGM, frame, "XML_finalFirstRoundEMnoSplit", "GMEMfinalFirstRoundEMnoSplit_frame");
#endif

#ifndef CELL_DIVISION_WITH_GPU_3DHAAR //to comment out cell division using GPU 3D Haar features + likelihood ratio score
            delete[] vecGMHOST; //we need it to estimate local likelihood before split
#else
        // (ngc) at some point I moved this code to it's own function....
        //       The code I moved seemed be missing a few things and
        //       I also didn't quite finish the refactor...

            cellDivisionWithGPU3DHaar.apply(vecGM, Fx, img, imgUINT16ptr, frame,
                devCUDA, configOptions, splitSet,
                singleCellSplitScore,
                debugPath,
                lht,
                vecGMHOST);
#endif //CELL_DIVISION_WITH_GPU_3DHAAR

        //---------------------------------------------------------------------------------------------------
        //-------------incorporate temporal logical rules in a time window to improve results----------------
			cout<<"******************************************************************************"<<endl;
            TicTocTimer ttTemporalLogicalRules = tic();

            //-----update lineage hyper tree with the final GMM results from this frame-----------------------
            TIME(UpdateLineageHyperTree(regularize_W4DOF, iniFrame, lht, vecGM, scaleOrig, frame, time_series_map.getProcessed(frame)));

            //-------------------------check if there are possible cell divisions(greedy)--------------------------------------------------
            int numCellDivisions;
            int numBirths;
            {
				float* im;
				TIME(im=time_series_map.getProcessed(frame));
				TIME(cellDivisionMinFlowTouchingSupervoxels(lht,frame,minNeighboringVoxelsConn3D,minNeighboringVoxels,
					 numCellDivisions,numBirths,regularize_W4DOF,im));
            }
            
            cout << "Generated (not touching supervoxels) " << numCellDivisions << " cell divisions out of " << lht.nucleiList[frame].size() - numCellDivisions << " cells in frame " << frame << ".Also added " << numBirths << " new tracks" << endl;

            //--------------------------------------------------------
            if (frame > iniFrame)
            {
                //This is just a patch tomake sure that breakCellDivisionBasedOnCellDivisionPlaneConstraint works correctly.
                //Some of the heuristic spatio-temporal rules change nuclei composition of supervoxels but do not update centroids.
                //Thus, in some cases teh midplane feature was miscalculated. This a quick fix to avoid that
                //TODO: find which heuristic rule changes the centroid and needs to be recalculated

                int TMaux = frame - 1;

                TIME(recalculateCentroidForMotherAndTwoDaughters(regularize_W4DOF, lht, scaleOrig, time_series_map, TMaux));

                // At this point we're done with the procesed float data from the previous frame, so flush it here.
                time_series_map.flush(frame - 1);

                //analyze cell divisions and cut links in the ones that do not satisfy the midplane division constraint
                int numCorrectionsD, numSplitsD;

                //2021 NEW
                if ( configOptions.cellDivisionClassifier.Amatf2013.use_2021_code ) 
                {
					err = 0;
                    //TIME(err = lht.breakCellDivisionBasedOnCellDivisionPlaneConstraint(TMaux, -1, numCorrectionsD, numSplitsD));//break all cell divisions and allow classifier code to find best option for each new track (or none if track is not result of a division)
                }
                else
                {
                    TIME(err = lht.breakCellDivisionBasedOnCellDivisionPlaneConstraint(TMaux, configOptions.thrCellDivisionPlaneDistance, numCorrectionsD, numSplitsD));//delete short living daughter
                    cout << "Cut " << numCorrectionsD << " linkages out of " << numSplitsD << " cell divisions because it did not satisfy the cell division midplane constraint with thr=" << configOptions.thrCellDivisionPlaneDistance << endl;
                }
                if (err > 0)
                    throw err;
                
            }
            //-------perform modifications using temporal logical rules----------------------------------
            //my sliding window is t \in [frame - 2 * configOptions.temporalWindowRadiusForLogicalRules, frame] => size of sliding window is 2 * configOptions.temporalWindowRadiusForLogicalRules + 1
            //lht still contains Tm frame - 2 * configOptions.temporalWindowRadiusForLogicalRules - 1 as an "anchor" time point: it cannot be modified by sliding window => solution should be consistent with it (we cannot merge two lineages present at athe acnhor time point)


            //parameters for logical temporal rules corrections (TODO: I should add them to some advance panel options later)
			if (frame >= iniFrame + 2 * configOptions.temporalWindowRadiusForLogicalRules) // frame is out of the initial window period
			{
				int numCorrections, numSplits;

				if (lengthTMthr > 0) {
#ifdef CELL_DIVISION_WITH_GPU_3DHAAR
					lht.mergeShortLivedAndCloseByDaughtersAll(lengthTMthr, frame - lengthTMthr, minNeighboringVoxels, minNeighboringVoxelsConn3D, numCorrections, numSplits);//merge first cell divisions
					cout << "Merged " << numCorrections << " out of " << numSplits << " splits because of a sibling death before " << lengthTMthr << " time points after cell division and touching tracks" << endl;
#else
					//lht.mergeNonSeparatingDaughtersAll(frame - 2 * configOptions.temporalWindowRadiusForLogicalRules + 1, minNeighboringVoxels, minNeighboringVoxelsConn3D, numCorrections, numSplits);//merge first cell divisions
					//cout<<"Merged "<<numCorrections<<" out of "<<numSplits<<" splits because of touching tracks"<<endl;
					lht.deleteShortLivedDaughtersAll(lengthTMthr, frame - lengthTMthr, numCorrections, numSplits);//delete short living daughter
					cout << "Deleted " << numCorrections << " out of " << numSplits << " splits because of a sibling death before " << lengthTMthr << " time points after cell division" << endl;

#endif
					float *im = time_series_map.getProcessed(frame);
					extendDeadNucleiAtTMwithHS(lht, hs, frame - 1, numCorrections, numSplits,im);//strictly speaking this is not a temporal feature, since it does not require a window, but it is still better to do it here (we can extend later)
					cout << "Extended " << numCorrections << " out of " << numSplits << " dead cells in frame " << frame - 1 << " using a simple local Hungarian algorithm with supervoxels" << endl;


					//redo in case other fixes have incorporated
					lht.deleteShortLivedDaughtersAll(lengthTMthr, frame - lengthTMthr, numCorrections, numSplits);//delete short living daughter
					cout << "Deleted " << numCorrections << " out of " << numSplits << " splits because of a sibling death before " << lengthTMthr << " time points after cell division" << endl;
					
					if ( !(configOptions.cellDivisionClassifier.Amatf2013.use_2021_code) ) 
                    {
						//--------------------------------------------------------
						//This is just a patch tomake sure that breakCellDivisionBasedOnCellDivisionPlaneConstraint works correctly. Som of the heuristic spatio-temporal rules change nuclei composition of supervoxels but do not update centroids.
						//Thus, in some cases teh midplane feature was miscalculated. This a quick fix to avoid that
						//TODO: find which heuristic rule changes the centroid and needs to be recalculated
	                    
	
						int TMaux = frame - 2 * configOptions.temporalWindowRadiusForLogicalRules;
						STATTIME_INIT;
	                    bool was_cached0,was_cached1;
	                    STATTIME(was_cached0=time_series_map.maybe_process(TMaux,hsVec[0]));
	                    STATTIME(was_cached1=time_series_map.maybe_process(TMaux+1,hsVec[1]));
						for (list<nucleus>::iterator iterN = lht.nucleiList[TMaux].begin(); iterN != lht.nucleiList[TMaux].end(); ++iterN)
						{
							TreeNode<ChildrenTypeLineage>* aux = iterN->treeNodePtr;
	
							if (aux->getNumChildren() != 2)
								continue;//not a cell division
	
							//recalculate centroid for mother and two daughters						
						    float *im;
							STATTIME(im = time_series_map.getProcessed(iterN));
							STATTIME(lht.calculateNucleiGMMParameters<float>(iterN, regularize_W4DOF, scaleOrig,im));
							STATTIME(im = time_series_map.getProcessed(iterN->treeNodePtr->left->data));
							STATTIME(lht.calculateNucleiGMMParameters<float>(iterN->treeNodePtr->left->data, regularize_W4DOF, scaleOrig,im));
							STATTIME(im = time_series_map.getProcessed(iterN->treeNodePtr->right->data));
							STATTIME(lht.calculateNucleiGMMParameters<float>(iterN->treeNodePtr->right->data, regularize_W4DOF, scaleOrig,im));
						}
						STATTIME_READ;
	                    if (!was_cached0)
	                        time_series_map.flush(TMaux);
	                    if (!was_cached1)
	                        time_series_map.flush(TMaux+1);
						//-----------------------------------------------------------
					}
					//analyze cell divisions and cut links in the ones that do not satisfy the midplane division constraint

                    //2021 NEW
                    if ( configOptions.cellDivisionClassifier.Amatf2013.use_2021_code ) 
                    {
						err = 0;
                        //err = lht.breakCellDivisionBasedOnCellDivisionPlaneConstraint(frame - 2 * configOptions.temporalWindowRadiusForLogicalRules, -1, numCorrections, numSplits);//delete short living daughter
                    }
                    else
                    {
                        err = lht.breakCellDivisionBasedOnCellDivisionPlaneConstraint(frame - 2 * configOptions.temporalWindowRadiusForLogicalRules, configOptions.thrCellDivisionPlaneDistance, numCorrections, numSplits);//delete short living daughter
                        cout << "Cut " << numCorrections << " linkages out of " << numSplits << " cell divisions because it did not satisfy the cell division midplane constraint with thr=" << configOptions.thrCellDivisionPlaneDistance << endl;
                    }
					if (err > 0)
						return err;
					
				} else {
					cout << "======WARNING: you have deactivated temporal logical rules (short-lived daughter threshold < 1)===================" << endl;
				}

            //run cell division discrimination based on 3D Haar elliptical features on a temporal window
                {
					//2021 NEW
                    if ( configOptions.cellDivisionClassifier.Amatf2013.use_2021_code ) 
                    {
						//set up parameters and cache images as needed
						int frameOffset = frame - configOptions.temporalWindowRadiusForLogicalRules;
						bool was_cached0,was_cached1;
						STATTIME_INIT;
	                    STATTIME(was_cached0=time_series_map.maybe_process(frameOffset,hsVec[configOptions.temporalWindowRadiusForLogicalRules]));
	                    STATTIME(was_cached1=time_series_map.maybe_process(frameOffset+1,hsVec[configOptions.temporalWindowRadiusForLogicalRules+1]));
						STATTIME_READ;
						
						//run division classifier
						float *im = time_series_map.getProcessed(frameOffset);
						float *im1 = time_series_map.getProcessed(frameOffset+1);
						TIME(numCorrections = cdwtClassifier->classifyCellDivisionTemporalWindow(lht, frameOffset, time_series_map.imgVecUINT16, devCUDA, configOptions.thrCellDivisionPlaneDistance,im,im1,regularize_W4DOF, scaleOrig ));

						if (numCorrections < 0)
                        	return numCorrections;
                    	TIME(numSplits = cdwtClassifier->getNumCellDivisions());
                    	cout << "Corrected " << numCorrections << " out of " << numSplits << " initially proposed cell divisions because cell division classifier with temporal window with thr =" << cdwtClassifier->getThrCDWT() << endl;
						
						//since bad divisions at frameOffset are now broken (had been left un-broken until this point), re-run the dead cell extender code to those timepoints
						extendDeadNucleiAtTMwithHS(lht, hs, frameOffset - 1, numCorrections, numSplits,im);//strictly speaking this is not a temporal feature, since it does not require a window, but it is still better to do it here (we can extend later)
						cout << "Extended " << numCorrections << " out of " << numSplits << " dead cells in frame " << frameOffset - 1 << " using a simple local Hungarian algorithm with supervoxels" << endl;
						extendDeadNucleiAtTMwithHS(lht, hs, frameOffset, numCorrections, numSplits,im1);//strictly speaking this is not a temporal feature, since it does not require a window, but it is still better to do it here (we can extend later)
						cout << "Extended " << numCorrections << " out of " << numSplits << " dead cells in frame " << frameOffset << " using a simple local Hungarian algorithm with supervoxels" << endl;
						
						//clean up the mess
						if (!was_cached0)
							time_series_map.flush(frameOffset);
						if (!was_cached1)
							time_series_map.flush(frameOffset+1);
                    }
                    else
                    {
						TIME(numCorrections = cdwtClassifier->classifyCellDivisionTemporalWindow(lht, frame - configOptions.temporalWindowRadiusForLogicalRules, time_series_map.imgVecUINT16, devCUDA));
						
                    	if (numCorrections < 0)
                        	return numCorrections;
                    	TIME(numSplits = cdwtClassifier->getNumCellDivisions());
                    	cout << "Corrected " << numCorrections << " out of " << numSplits << " initially proposed cell divisions because cell division classifier with temporal window with thr =" << cdwtClassifier->getThrCDWT() << endl;
					}
                }

                //background detection is always run (then we can decide what to do with the results)
                if (configOptions.temporalWindowRadiusForLogicalRules < temporalWindowSizeForBackground)
					//throw runtime_error("ERROR: mainTrackingGaussianMixtureModel: temporalWindowRadiusForLogicalRules cannot be smaller than temporalWindowSizeForBackground");
					temporalWindowSizeForBackground = configOptions.temporalWindowRadiusForLogicalRules;

                {
                    time_t t_start, t_end;
                    time(&t_start);
                    //store the probability of being background for each nuclei
                    auto err = setProbBackgroundTracksAtTM(lht, frame - configOptions.temporalWindowRadiusForLogicalRules, temporalWindowSizeForBackground, backgroundClassifierFilename, devCUDA);
                    if (err > 0)
                        throw err;
                    time(&t_end);
                    cout << "Scored background tracks in " << difftime(t_end, t_start) << " secs" << endl;
                }


                //-------delete lineages that are shorter than a certain lenght--------				
                lht.deleteShortLivedLineagesAll(frame - 2 * configOptions.temporalWindowRadiusForLogicalRules + 1, lengthTMthr, numCorrections, numSplits);//delete short lineages
                cout << "Deleted " << numCorrections << " out of " << numSplits << " lineages born at frame " << frame - 2 * configOptions.temporalWindowRadiusForLogicalRules + 1 << " because of length shorter than " << lengthTMthr << " time points" << endl;

            }
            else if (frame == iniFrame + 2 * configOptions.temporalWindowRadiusForLogicalRules - 1) // frame is end of initial window (end of burn in)
            {
                //special case to merge tracks due to initial oversegmentation
                int numCorrections, numSplits, numMerges;


				if(lengthTMthr>0) // short lived daughter threshold active
				{
#ifdef CELL_DIVISION_WITH_GPU_3DHAAR
					lht.mergeShortLivedAndCloseByDaughtersAll(lengthTMthr, frame - lengthTMthr, minNeighboringVoxels, minNeighboringVoxelsConn3D, numCorrections, numSplits);//merge first cell divisions
					cout << "Merged " << numCorrections << " out of " << numSplits << " splits because of a sibling death before " << lengthTMthr << " time points after cell division and touching tracks" << endl;
#else
					//lht.mergeNonSeparatingDaughtersAll(frame - 2 * configOptions.temporalWindowRadiusForLogicalRules + 1, minNeighboringVoxels, minNeighboringVoxelsConn3D, numCorrections, numSplits);//merge first cell divisions
					//cout<<"Merged "<<numCorrections<<" out of "<<numSplits<<" splits because of touching tracks"<<endl;
					lht.deleteShortLivedDaughtersAll(lengthTMthr, frame - lengthTMthr, numCorrections, numSplits);////delete short living daughter
					cout << "Deleted " << numCorrections << " out of " << numSplits << " splits because of a sibling death before " << lengthTMthr << " time points after cell division" << endl;
#endif
				}
            //delete all the tracks that have died during this "burn in" period
                TIME(lht.deleteDeadBranchesAll(frame, numMerges));
                cout << "Deleted " << numMerges << " lineages that died during the burning period" << endl;


            }

            //-------------------------------------------------------------------------------------
            //---------------------re-calculate thrDist2LargeDisplacement based on current data-----

	        TIME(recalculateThrDist2LargeDisplacement(iniFrame, lht, thrDist2LargeDisplacement, frame));
	        //-------------------------------------------------------------------------------------


            //------save last element before being removed of the window----------------------------------------
            if (frame >= iniFrame + 2 * configOptions.temporalWindowRadiusForLogicalRules)
            {
                int frameOffset = frame - 2 * configOptions.temporalWindowRadiusForLogicalRules;

                //the only thing I need to modify is the parentId and the the neighbors
                {
					auto was_cached=time_series_map.maybe_process(frameOffset,hsVec[0]);
                    TIME(parseNucleiList2TGMM<float>(vecGM, lht, frameOffset, regularize_W4DOF, thrDist2LargeDisplacement, time_series_map.getProcessed(frameOffset)));
                    if (!was_cached)
                        time_series_map.flush(frameOffset);
                }


                //set wildcard to idx for frameOffset-1 so we can set parentIdx
                if (frameOffset > 0)
                {
                    int countV = 0;
                    for (list< nucleus >::iterator iterN = lht.nucleiList[frameOffset - 1].begin(); iterN != lht.nucleiList[frameOffset - 1].end(); ++iterN, ++countV)
                        iterN->tempWilcard = (float)countV;
                }
                TIME(WriteGMEMResultXML(debugPath, vecGM, frameOffset, lht));
	            TIME(saveSupervoxels(debugPath, lht, frameOffset)); //save supervoxels so I can visualize them

                //delete time point frameOffset since it is not needed as an anchor point anymore
                TIME(lht.setFrameAsT_o(frameOffset));  // this ends up deleting the bound image data at frameOffset-1 === frame - 2 * window - 1
                // eraseTM(time_series_map, frameOffset - 1);

                //delete hierarchical segmentation for this time point
                TIME(delete hsVec.front());
                TIME(hsVec.erase(hsVec.begin()));
            }
            cout << "Applying all the temporal logical rules took " << toc(&ttTemporalLogicalRules) << " secs" << endl;
			cout<<"******************************************************************************"<<endl;
            //--------------end of temporal logical rules--------------------------------------------------
            //---------------------------------------------------------------------------------------------




            //release memory for each frame
            //mylib::Free_Array(img); // (ngc) <-- was commented out | ---> not true? --> //this memory will be freed by lht.setFrameAsT_o (we need the image for the temporal sliding window)
            img = NULL;
            //imgUINT16 = NULL;//this is freed by imgVecUINT16

            if (imgFlow != NULL)
            {
                mylib::Free_Array(imgFlow);
                imgFlow = NULL;
            }
            if (imgFlowMask != NULL)
            {
                mylib::Free_Array(imgFlowMask);
                imgFlowMask = NULL;
                imgFlowMaskPtr = NULL;
            }

            if ( !(configOptions.cellDivisionClassifier.Amatf2013.use_2021_code) ) 
			{
                GMEMreleaseMemoryWithSupervoxels(&queryCUDA, &imgDataCUDA, &rnkCUDA, &indCUDA, &indCUDAtr, &centroidLabelPositionCUDA, &labelListPtrCUDA);
                delete[]queryHOST;
                delete[]imgDataHOST;
                delete[]centroidLabelPositionHOST;
                delete[]centroidLabelWeightHOST;
                delete[]labelListPtrHOST;
			}
            cout << toc(&tt) << " secs" << endl;
			
			/*DEBUG -- detect potential memory leaks
			HANDLE_ERROR( cudaSetDevice( devCUDA ) );
			HANDLE_ERROR( cudaDeviceReset() );*/


            //------------------------------------------------------------------------------
            //check number of deaths and if we need to rerun the same point with optical flow
            if (configOptions.deathThrOpticalFlow >= 0 && frame != iniFrame) //for iniFrame we can not compute flow
            {
                ProcessDeathThrOpticalFlow(configOptions, numDeaths, vecGM, frame, debugPath, scaleOrig, imgFlowMask, imgDims, imgSize, imgFlowMaskPtr);
            }
            //-----------------------------end of if(configOptions.deathThrOpticalFlow>0)--------------------------------


        }//end of for(frame=...) loop


        //--------------------------------------------------------------------------------
        //-----------flush out the last time points in the lineage that were not saved because of teh sliding window approach-------------
        for (int frame = endFrame + 1 - 2 * configOptions.temporalWindowRadiusForLogicalRules; frame <= endFrame; frame++)
        {
            if(frame<iniFrame) continue;
            cout<<"Flushing last time points.  Frame " << frame << endl;
            //analyze cell divisions and cut links in the ones that do not satisfy the midplane division constraint
            if (frame < endFrame)
            {
                int numCorrections, numSplits;
                int err = lht.breakCellDivisionBasedOnCellDivisionPlaneConstraint(
                    frame, configOptions.thrCellDivisionPlaneDistance, numCorrections, numSplits);//delete short living daughter
                if (err > 0)
                    return err;
                cout << "Cut " << numCorrections << " linkages out of " << numSplits << " cell divisions because it did not satisfy the cell division midplane constraint with thr=" << configOptions.thrCellDivisionPlaneDistance << " (frame" << frame << ")" << endl;
            }
            //save time point frame
            {
				auto was_cached=time_series_map.maybe_process(frame,hsVec[0]);
                parseNucleiList2TGMM<float>(vecGM, lht, frame, regularize_W4DOF, thrDist2LargeDisplacement, time_series_map.getProcessed(frame));
                if (!was_cached)
                    time_series_map.flush(frame);
            }

            //set wildcard to idx for frameOffset-1 so we can set parentIdx
            if (frame > 0)
            {
                int countV = 0;
                for (list< nucleus >::iterator iterN = lht.nucleiList[frame - 1].begin(); iterN != lht.nucleiList[frame - 1].end(); ++iterN, ++countV)
                    iterN->tempWilcard = (float)countV;
            }
            mkdir(debugPath + "XML_finalResult_lht");
            string GMxmlFilename = string(debugPath + "XML_finalResult_lht/GMEMfinalResult_frame" + to_string(frame, 4, '0') + ".xml");

            ofstream outXML(GMxmlFilename.c_str());
            GaussianMixtureModel::writeXMLheader(outXML);
            int countV = 0;
            for (list< nucleus >::iterator iterN = lht.nucleiList[frame].begin(); iterN != lht.nucleiList[frame].end(); ++iterN, ++countV)
            {
                //modify parentId
                if (iterN->treeNodePtr->parent == NULL)
                    vecGM[countV]->parentId = -1;
                else{
                    vecGM[countV]->parentId = (int)(iterN->treeNodePtr->parent->data->tempWilcard);
                }
                //write solution
                vecGM[countV]->writeXML(outXML);
            }
            GaussianMixtureModel::writeXMLfooter(outXML);
            outXML.close();

            //save supervoxels so I can visualize them
            string filenameSv(debugPath + "XML_finalResult_lht/GMEMfinalResult_frame" + to_string(frame, 4, '0') + ".svb");
            if (lht.writeListSupervoxelsToBinaryFile(filenameSv, frame) > 0)
                exit(3);

            //delete hierarchical segmentation for this time point
            delete hsVec.front();
            hsVec.erase(hsVec.begin());
        }
        //-----------------------------------------------------



        //------------------------------------------------------
        //release memory
        for (unsigned int ii = 0; ii < vecGM.size(); ii++)
            delete vecGM[ii];
        vecGM.clear();
        Fx.clear();



        supervoxel::freeTrimParameters();

        //run background "cleaner"
        if (configOptions.thrBackgroundDetectorHigh < 1.0f)
        {
            TicTocTimer tt = tic();
            cout << "Running forward-backward pass with hysteresis threshold background =(" << configOptions.thrBackgroundDetectorLow << "," << configOptions.thrBackgroundDetectorHigh << ")" << " to remove non-cell like objects" << endl;

            string outputFolder(debugPath + "XML_finalResult_lht_bckgRm");
            mkdir(outputFolder);
            int err = applyProbBackgroundHysteresisRulePerBranch(string(debugPath + "XML_finalResult_lht/GMEMfinalResult_frame"), iniFrame, endFrame, string(outputFolder + "/"), configOptions.thrBackgroundDetectorLow, configOptions.thrBackgroundDetectorHigh);
            cout << toc(&tt) << " secs" << endl;
            if (err > 0)
                return err;
        }


        return 0;
    }
    catch (exception e) {
        cerr << e.what() << endl;
        return 1;
    }

}
