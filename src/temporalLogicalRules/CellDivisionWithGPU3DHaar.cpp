#include "CellDivisionWithGPUHaar.h"
#include "UtilsCUDA/3DEllipticalHaarFeatures/EllipticalHaarFeatures.h"
#include "UtilsCUDA/3DEllipticalHaarFeatures/gentleBoost/gentleBoost.h"


//#define WORK_WITH_UNITS //uncomment this define to work with metric units (instead of pixels). So scale is incoporated directly in all the calculation.

CellDivisionWithGPU3DHaar::CellDivisionWithGPU3DHaar(string pathS) {
        //load classifier
    string classifierCellDivisionFilename(pathS + string("/") +"classifierCellDivision.txt");
    auto errC = loadClassifier(classifierCellDivision, classifierCellDivisionFilename);
    if(errC > 0 ) 
        throw errC;
    numFeatures = getNumberOfHaarFeaturesPerllipsoid();
}

 void CellDivisionWithGPU3DHaar::apply(vector<GaussianMixtureModel*> &vecGM,
                                       vector<feature> &Fx,
                                       const mylib::Array* img,
                                       const mylib::uint16* imgUINT16ptr,
                                       int frame,
                                       int devCUDA,
                                       configOptionsTrackingGaussianMixture &configOptions,
                                       vector<pair<GaussianMixtureModel*,GaussianMixtureModel*> > &splitSet,
                                       vector<double> &singleCellSplitScore,
                                       string debugPath,
                                       lineageHyperTree &lht,
                                       GaussianMixtureModelCUDA *vecGMHOST
                                      )
{

    vector<GaussianMixtureModel> backupVecGM;
    //----------- check cell divisions using image features (GPU) and machine learning classifier----------------
    TicTocTimer ttCellDivision = tic();
    cout<<"Checking which ellipsoids out of "<<vecGM.size()<< " might be dividing with trained classifier"<<endl;
    Fx.resize(vecGM.size());
    //allocate memory and copy values
    long long int *dimsVec = new long long int[img->ndims];
    for(int aa = 0; aa < img->ndims; aa++)
        dimsVec[aa] = img->dims[aa];
    int numEllipsoids = vecGM.size();
    int sizeW = dimsImage * (1 + dimsImage) / 2;
    double *m = new double[dimsImage * numEllipsoids];
    double *W = new double[sizeW  * numEllipsoids];
    for(int jj = 0; jj<numEllipsoids; jj++)
    {
        int countW = jj;
        int countM = jj;
        for(int aa = 0; aa<dimsImage; aa++)
        {
            m[countM] = vecGM[jj]->m_k(aa);
            for( int bb = aa; bb < dimsImage; bb++)
            {
                W[ countW ] = vecGM[jj]->W_k(aa,bb) * vecGM[jj]->nu_k;
                countW += numEllipsoids;
            }
            countM += numEllipsoids;
        }           
    }
    //cout<<"Allocating memory took "<<toc(&ttCellDivision)<<" secs"<<endl;

    //calculate features 
    basicEllipticalHaarFeatureVector **fBasic = calculateEllipticalHaarFeatures(m, W,vecGM.size(),imgUINT16ptr,dimsVec,devCUDA,0);
    if(fBasic == NULL)
        exit( 1 );//some error happened

                 //cout<<"Calculating basic features took "<<toc(&ttCellDivision)<<" secs"<<endl;
                 //extend Haar features
    int numHaarFeaturesPerEllipsoid = 0;
    float* xTest = NULL;
    calculateCombinationsOfBasicHaarFeatures(fBasic,vecGM.size(),&numHaarFeaturesPerEllipsoid, &xTest);
    if(xTest == NULL)
        exit( 2 );//some error happened
    if( numHaarFeaturesPerEllipsoid != numFeatures )
    {
        cout<<"ERROR: numFeatures "<<numFeatures<<" is different than numHaarFeaturesPerEllipsoid "<<numHaarFeaturesPerEllipsoid<<endl;
        exit( 3 );
    }
    //cout<<"Calculating extended features took "<<toc(&ttCellDivision)<<" secs"<<endl;
    //calculate classifier results
    //transposeXtrainOutOfPlace(xTest, Fx.size(), numFeatures);//we need to perform transposition from GPU features to gentleBoost classifier               
    boostingTreeClassifierTranspose(xTest,&(Fx[0]) ,classifierCellDivision , Fx.size(), numFeatures);

    if(vecGM.size()!=Fx.size())
    {
        cout<<"ERROR: after calling Matlab C Shared library for split classifier: vecGM and Fx are not the same size!"<<endl;
        exit(3);
    }
    for(unsigned int ii=0;ii<vecGM.size();ii++) vecGM[ii]->splitScore=Fx[ii];

    //release memory
    delete[] m;
    delete[] W;
    delete[] dimsVec;
    for(int ii=0;ii<vecGM.size();ii++) 
        delete fBasic[ii];
    delete[] fBasic;
    delete[] xTest;

    cout<<"Running classifer took "<<toc(&ttCellDivision)<<" secs"<<endl;
    //--------------------------end of cell division detection with image features + machine learning -----------------------------------------------


    //generate set S containing all candidates
    unsigned int sizeGMvec=vecGM.size();
    splitSet.clear();
    GaussianMixtureModel* GM;
    singleCellSplitScore.clear();
    backupVecGM.clear();
    vector<int> dividingNucleiIdx;//stores idx so we know which ones to compute local likelihood
    dividingNucleiIdx.reserve(vecGM.size() / 10);
    for(unsigned int kk=0;kk<sizeGMvec;kk++)
    {
        if(vecGM[kk]->splitScore>configOptions.thrSplitScore && vecGM[kk]->isDead()==false)//cell is a candidates to split
        {
            dividingNucleiIdx.push_back(kk);
            backupVecGM.push_back(*(vecGM[kk]));
            singleCellSplitScore.push_back(vecGM[kk]->splitScore);

            //generate split
            //update Gaussian mixture
            vecGM.push_back(new GaussianMixtureModel());
            vecGM[kk]->splitGaussian(vecGM.back(),vecGM.size()-1,imgData,img->dims);//using k-means
            GM=vecGM.back();//returns pointer of the new Gaussian

            splitSet.push_back(make_pair(vecGM[kk],GM));

            //update priors
            GM->updatePriorsAfterSplitProposal();
            vecGM[kk]->updatePriorsAfterSplitProposal();
        }else{
            vecGM[kk]->fixed=true;//we can not modify this mixture
                                  //if(vecGM[kk]->splitScore<0.0) vecGM[kk]->splitScore=0.0;//to avoid rounding errors
        }
    }


    cout<<"Number of proposed splits="<<splitSet.size()<<endl;
    //debug--------------------------------------------
    char buffer[128];
    sprintf(buffer,"%.4d",frame);
    string itoa(buffer);
    struct stat St;
    if (stat( (debugPath+"XML_splitCandidates").c_str(), &St ) != 0)//check if folder exists
    {
        auto error = mkDirectory(string(debugPath + "XML_splitCandidates"));
        if(error>0)
            exit(error);
    }
    {
        string GMxmlFilename(debugPath+"XML_splitCandidates/"+"GMEMsplitCandidates_frame"+ itoa + ".xml");
        ofstream outXML(GMxmlFilename.c_str());
        GaussianMixtureModel::writeXMLheader(outXML);
        for(unsigned int ii=0;ii<vecGM.size();ii++)
        {
#ifdef WORK_WITH_UNITS
            vecGM[ii]->units2pixels(scaleOrig);
            vecGM[ii]->writeXML(outXML);
            vecGM[ii]->pixels2units(scaleOrig);
#else
            vecGM[ii]->writeXML(outXML);
#endif
        }
        GaussianMixtureModel::writeXMLfooter(outXML);
        outXML.close();
    }
    //-------------------------------------------------------------------------
    if(splitSet.empty()==false){

        //--------------------------------calculate local likelihood around all Gaussians that have divided (likelihood before division)---------------------------
        //vecGMHOST still contains the information before split
        //create a list of supervoxel indexes that are used for each local likelihood
        vector< vector<int> > listSupervoxelIdx(dividingNucleiIdx.size());
        vector< list<supervoxel>::iterator > listSupervoxelIterators;
        lht.getSupervoxelListIteratorsAtTM(listSupervoxelIterators,frame);
        int countSN=0;
        for(list< supervoxel>::iterator iter=lht.supervoxelsList[frame].begin(); iter!=lht.supervoxelsList[frame].end(); ++iter,++countSN)
            iter->tempWildcard=(float)countSN;//so I can retrieve ordering

        countSN=0;
        for(vector<int>::iterator iter=dividingNucleiIdx.begin(); iter!=dividingNucleiIdx.end(); ++iter,++countSN){
            listSupervoxelIdx[countSN].reserve(configOptions.KmaxNumNNsupervoxel+1);
            for(int ii=0; ii<vecGMHOST[*iter].supervoxelNum; ii++){
                int supervoxelIdxAux=vecGMHOST[*iter].supervoxelIdx[ii];
                listSupervoxelIdx[countSN].push_back(supervoxelIdxAux);
                //we only consider the local likelihood in the supervoxels belonging to the dividing cell (very local). Otherwise test is not so meaningful
                //TODO: propose split using supervoxels (not K-means with local intensity). There is finite number of combinations
                //for(vector< SibilingTypeSupervoxel >::iterator iterS = listSupervoxelIterators[ supervoxelIdxAux ]->nearestNeighborsInSpace.begin(); iterS != listSupervoxelIterators[ supervoxelIdxAux ]->nearestNeighborsInSpace.end(); ++iterS)
                //  listSupervoxelIdx[countSN].push_back(( (int) ((*iterS)->tempWildcard) ) );
            }
        }
        //calculate and store likelihood for each of the elements
        vector<double> localLikelihoodBeforeSplit;
        localLikelihoodBeforeSplit.reserve(dividingNucleiIdx.size());

        calculateLocalLikelihood(localLikelihoodBeforeSplit,listSupervoxelIdx,queryCUDA,imgDataCUDA,rnkCUDA,indCUDA,vecGMHOST,labelListPtrCUDA,query_nb,vecGM_ref_nb,numLabels);

        delete[] vecGMHOST;
    }
        //-------------------------------------------------------------------------------------------------------------------------

        //string debugPathCUDAb(debugPath+"XML_GMEMiterations_CUDA_afterSplit_frame");
        string debugPathCUDAb("");
        {
            auto vecGMHOST=new GaussianMixtureModelCUDA[vecGM.size()];
            for(unsigned int ii=0;ii<vecGM.size();ii++) copy2GMEM_CUDA(vecGM[ii],&(vecGMHOST[ii]));
            //GMEMvariationalInferenceCUDA(queryCUDA,imgDataCUDA,rnkCUDA,rnkCUDAtr,indCUDA,indCUDAtr,vecGMHOST,query_nb,vecGM.size(),configOptions.maxIterEM,configOptions.tolLikelihood,devCUDA,frame,true,debugPathCUDAb);
            GMEMvariationalInferenceCUDAWithSupervoxels(queryCUDA,imgDataCUDA,rnkCUDA,rnkCUDAtr,indCUDA,indCUDAtr,centroidLabelPositionCUDA,labelListPtrCUDA,vecGMHOST,query_nb,vecGM.size(),numLabels,configOptions.maxIterEM,configOptions.tolLikelihood,devCUDA,frame,regularize_W4DOF,debugPathCUDAb);
            for(unsigned int ii=0;ii<vecGM.size();ii++) copyFromGMEM_CUDA(&(vecGMHOST[ii]),vecGM[ii]);


            //--------------------------------calculate local likelihood around all Gaussians that have divided---------------------------
            vector<double> localLikelihoodAfterSplit;
            localLikelihoodAfterSplit.reserve(dividingNucleiIdx.size());
            calculateLocalLikelihood(localLikelihoodAfterSplit,listSupervoxelIdx,queryCUDA,imgDataCUDA,rnkCUDA,indCUDA,vecGMHOST,labelListPtrCUDA,query_nb,vecGM.size(),numLabels);

            //calculate final likelihood test ratio
            vector<double> splitScoreVec(dividingNucleiIdx.size());
            for(size_t aa=0; aa<dividingNucleiIdx.size(); aa++)
                splitScoreVec[aa]=vecGM[dividingNucleiIdx[aa]]->splitScore;
            likelihoodRatioTestSplitMerge(localLikelihoodBeforeSplit,localLikelihoodAfterSplit,splitScoreVec,splitMergeTest);
            //-------------------------------------------------------------------------------------------------------------------------

            delete[] vecGMHOST;
        }

        //debug--------------------------------------------
        sprintf(buffer,"%.4d",frame);
        string itoa=string(buffer);
        if (stat( (debugPath+"XML_finalSecondRoundEMSplit").c_str(), &St ) != 0)//check if folder exists
        {
            auto error = mkDirectory(string(debugPath + "XML_finalSecondRoundEMSplit"));
            if(error>0)
                exit(error);
        } 
        {
            string GMxmlFilename(debugPath+ "XML_finalSecondRoundEMSplit/GMEMfinalSecondRoundEMSplit_frame"+ itoa + ".xml");
            string xmlFilenameAfterSplit=GMxmlFilename;
            ofstream outXML(GMxmlFilename.c_str());
            GaussianMixtureModel::writeXMLheader(outXML);
            for(unsigned int ii=0;ii<vecGM.size();ii++)
            {
#ifdef WORK_WITH_UNITS
            vecGM[ii]->units2pixels(scaleOrig);
            vecGM[ii]->writeXML(outXML);
            vecGM[ii]->pixels2units(scaleOrig);
#else
                vecGM[ii]->writeXML(outXML);
#endif
            }
            GaussianMixtureModel::writeXMLfooter(outXML);
            outXML.close();
        }
        //-------------------------------------------------------------------------
        //likelihood test ratio to decide which cell divisions where successful (done in the GPU now)

        int numDiv=0;
        int pos=0;
        for(unsigned int kk=0;kk<splitSet.size();kk++)
        {
            if(splitSet[kk].first->isDead() || splitSet[kk].second->isDead())//if one cell is dead, then whether the other cell is alive or also dead the action is the same: delete second cell and restore previous GM
            {

                pos=splitSet[kk].second->id;
                for(unsigned int ss=pos+1;ss<vecGM.size();ss++) vecGM[ss]->id--;
                delete vecGM[pos];
                vecGM.erase(vecGM.begin()+pos);

                //reset Gaussian parameters to before the split
                pos=splitSet[kk].first->id;
                *(splitSet[kk].first)=backupVecGM[kk];//reset Gaussian parameters
                splitSet[kk].first->id=pos;
            }else{//both cell divisions are alive
                  //if((splitSet[kk].second->splitScore+splitSet[kk].first->splitScore)>singleCellSplitScore[kk])//division proposal was not accepted
                  //if(max(splitSet[kk].second->splitScore,splitSet[kk].first->splitScore)>singleCellSplitScore[kk])//division proposal was not accepted
                if(splitMergeTest[kk]<0)//division proposal rejected
                {
                    pos=splitSet[kk].second->id;
                    for(unsigned int ss=pos+1;ss<vecGM.size();ss++) vecGM[ss]->id--;
                    delete vecGM[pos];
                    vecGM.erase(vecGM.begin()+pos);

                    pos=splitSet[kk].first->id;
                    *(splitSet[kk].first)=backupVecGM[kk];//reset Gaussian parameters
                    splitSet[kk].first->id=pos;
                }else{//division proposal was accepted by the likelihood ratio test
                    numDiv++;
                    //update priors since they were changed to test splits
                    splitSet[kk].first->updatePriors(configOptions.betaPercentageOfN_k,configOptions.nuPercentageOfN_k,configOptions.alphaPercentage, -1.0);//alpha is not altered
                    splitSet[kk].second->updatePriors(configOptions.betaPercentageOfN_k,configOptions.nuPercentageOfN_k, configOptions.alphaPercentage, -1.0);//alpha is not altered
                }
            }
        }

        cout<<"Accepted "<<numDiv<<" cell divisions"<<endl;

        //special case for first frame in order to create a single lineage per initial nuclei after division
        if(frame == iniFrame) // && argc<5)
        {
            //redo all the lineages id
            int llId = 0;
            for(vector<GaussianMixtureModel*>::iterator iter = vecGM.begin(); iter != vecGM.end(); ++iter)
            {
                (*iter)->lineageId = llId;
                llId++;
            }
        }


        //string debugPathCUDAc(debugPath+"XML_GMEMiterations_CUDA_finalRound_frame");
        sprintf(buffer,"%.4d",frame);//to save rnk and ind for each supervoxel
        itoa=string(buffer);
        string debugPathCUDAc(debugPath+"XML_finalResult");
        stat St;
        if (stat( debugPathCUDAc.c_str(), &St ) != 0)//check if folder exists
        {
            int error = mkDirectory(debugPathCUDAc);
            if(error>0)
                exit(error);
        }
        debugPathCUDAc=string(debugPathCUDAc+ "/rnk_frame"+ itoa + ".bin");

        GaussianMixtureModelCUDA *vecGMHOST=new GaussianMixtureModelCUDA[vecGM.size()];
        for(unsigned int ii=0;ii<vecGM.size();ii++) copy2GMEM_CUDA(vecGM[ii],&(vecGMHOST[ii]));
        //GMEMvariationalInferenceCUDA(queryCUDA,imgDataCUDA,rnkCUDA,rnkCUDAtr,indCUDA,indCUDAtr,vecGMHOST,query_nb,vecGM.size(),configOptions.maxIterEM,configOptions.tolLikelihood,devCUDA,frame,true,debugPathCUDAc);
        GMEMvariationalInferenceCUDAWithSupervoxels(queryCUDA,imgDataCUDA,rnkCUDA,rnkCUDAtr,indCUDA,indCUDAtr,centroidLabelPositionCUDA,labelListPtrCUDA,vecGMHOST,query_nb,vecGM.size(),numLabels,configOptions.maxIterEM,configOptions.tolLikelihood,devCUDA,frame,regularize_W4DOF, debugPathCUDAc);
        for(unsigned int ii=0;ii<vecGM.size();ii++) copyFromGMEM_CUDA(&(vecGMHOST[ii]),vecGM[ii]);
        delete[] vecGMHOST;

    }//end of if(splitSet.empty==false)

     //allow all the mixtures to be modified again
    unfixedAllMixtures(vecGM);

    //-------------------write out final result for frame-------------------------
    char buffer[128];
    sprintf(buffer,"%.4d",frame);
    string itoa=string(buffer);
    struct stat St;
    if (stat( (debugPath+"XML_finalResult").c_str(), &St ) != 0)//check if folder exists
    {
        int error = mkDirectory(string(debugPath + "XML_finalResult"));
        if (error > 0)
            exit( error);
    }

    {
        string GMxmlFilename(debugPath+"XML_finalResult/GMEMfinalResult_frame"+ itoa + ".xml");
        ofstream outXML(GMxmlFilename.c_str());
        GaussianMixtureModel::writeXMLheader(outXML);
        for(unsigned int ii=0;ii<vecGM.size();ii++)
        {
#ifdef WORK_WITH_UNITS
        vecGM[ii]->units2pixels(scaleOrig);
        vecGM[ii]->writeXML(outXML);
        vecGM[ii]->pixels2units(scaleOrig);
#else
            vecGM[ii]->writeXML(outXML);
#endif
        }
        GaussianMixtureModel::writeXMLfooter(outXML);
        outXML.close();
    }
}