/*
 * Copyright (C) 2011-2013 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat 
 *  supportFunctionsMain.cpp
 *
 *  Created on: January 21st, 2012
 *      Author: Fernando Amat
 *
 * \brief Contains different routines that are called form the main.cpp part of theprogram to not make the code so cluttr
 *        
 *
 */

#include <iostream>
#include <queue>

#if defined(_WIN32) || defined(_WIN64)
	#define NOMINMAX
	#include <Windows.h>
	
#else //Linux headers
	#include <unistd.h>
	#include "sys/stat.h"
#endif

#include "supportFunctionsMain.h"
#include "temporalLogicalRules/trackletCalculation.h"

using namespace std;



//=========================================================================
template <class imgTypeC>
int parseNucleiList2TGMM(std::vector<GaussianMixtureModel*> &vecGM,lineageHyperTree &lht,int frame, bool regularizeW, float thrDist2LargeDisplacement,imgTypeC *imgPtr)
{
	if(frame >= lht.getMaxTM())
		return 0;

	list<nucleus>* listNucleiPtr = (&(lht.nucleiList[frame]));
	size_t N = listNucleiPtr->size();
	//cout << "     P1." << endl;
	//resize vecGM
	if( N > vecGM.size())//I need to allocate
	{
		//cout << "     P1a." << endl;
		vecGM.reserve( N );
		for(size_t ii = vecGM.size(); ii<N; ii++)
			vecGM.push_back( new GaussianMixtureModel(ii) );
	}else{//I need to deallocate
		//cout << "     P1b." << endl;
		for(size_t ii = N; ii < vecGM.size(); ii++)
			delete vecGM[ii];
		vecGM.resize(N);
	}
	//cout << "     P2." << endl;

	//reset supervoxels id
	int count = 0;
	for(list<supervoxel>::iterator iterS = lht.supervoxelsList[frame].begin(); iterS != lht.supervoxelsList[frame].end(); ++iterS, count++)
		iterS->tempWildcard = (float) count;

	//cout << "     P3." << endl;
	//compute all the stats from centroids	
    float N_k;
    float m_k[dimsImage] = { 0 };
    float S_k[dimsImage * (1 + dimsImage) / 2] = { 0 };
    count = 0;
    int countW;
    cout << "     P4." << endl;
    for (list<nucleus>::iterator iterN = listNucleiPtr->begin(); iterN != listNucleiPtr->end(); ++iterN, ++count)
    {

        lht.calculateGMMparametersFromNuclei<imgTypeC>(iterN, m_k, &N_k, S_k, imgPtr);

        //copy results
        vecGM[count]->N_k = N_k;

        vecGM[count]->alpha_k = N_k;
        vecGM[count]->beta_k = N_k;
        vecGM[count]->nu_k = N_k;
        vecGM[count]->alpha_o = 0;//there are no priors here
        vecGM[count]->beta_o = 0;
        vecGM[count]->nu_o = dimsImage + 1;//minimum degrees of freedom are dimsImage +1
        vecGM[count]->sigmaDist_o = iterN->probCellDivision;//we use this now to store cell division probability

        countW = 0;
        for (int ii = 0; ii < dimsImage; ii++)
        {
            vecGM[count]->m_k(ii) = m_k[ii];
            vecGM[count]->m_o(ii) = m_k[ii];
            for (int jj = ii; jj < dimsImage; jj++)
            {
                vecGM[count]->W_k(ii, jj) = S_k[countW];
                vecGM[count]->W_k(jj, ii) = S_k[countW];
                countW++;
            }
        }

        if (vecGM[count]->W_k(dimsImage - 1, dimsImage - 1) < 1e-8)//it means we have a 2D ellipsoid->just regularize W with any value (we will ignore it in the putput anyway)
        {
            vecGM[count]->W_k(dimsImage - 1, dimsImage - 1) = 0.5 * (vecGM[count]->W_k(0, 0) + vecGM[count]->W_k(1, 1));
        }

        vecGM[count]->W_k = (vecGM[count]->W_k.inverse()) / vecGM[count]->nu_k;




        if (regularizeW == true)
            vecGM[count]->regularizePrecisionMatrix(true);//W_4DOF is set to true here

        vecGM[count]->W_o = vecGM[count]->W_k * vecGM[count]->nu_k / vecGM[count]->nu_o;
        //VIP: we do not recover parent information. We just assume that parentId is the index of the nucleiList[frame]
        vecGM[count]->parentId = count;

        vecGM[count]->lineageId = count;//some pieces of the code need this info to be different than -1;in this case everything is from a different lineage
	}
	count = 0;
	cout << "     P5." << endl;
    for (list<nucleus>::iterator iterN = listNucleiPtr->begin(); iterN != listNucleiPtr->end(); ++iterN, ++count)
    {
        //save supervoxel
        vecGM[count]->supervoxelIdx.resize(iterN->treeNode.getNumChildren());
        int countS = 0;
        for (vector<ChildrenTypeNucleus>::iterator iter = iterN->treeNode.getChildren().begin(); iter != iterN->treeNode.getChildren().end(); ++iter, countS++)
        {
            vecGM[count]->supervoxelIdx[countS] = (int)((*iter)->tempWildcard);
        }

        //save background probability for this nucleus
        vecGM[count]->beta_o = iterN->probBackground;//TODO: use its own variable, although here it was set to zero anyways

        //save split score based on confidence: link to CATMAID editing
        //vecGM[count]->splitScore = iterN->confidence;		
        vecGM[count]->splitScore = lht.confidenceScoreForNucleus(iterN->treeNodePtr, thrDist2LargeDisplacement); //TODO: calculate this before hand and save it in confidence

        //vecGM[count]->splitScore = iterN->debugVisualization;//uncomment this to be able to visualize specific cases
    }
	//cout<<"=============DEBUGGING: confidence split score disconnected to study incoherences========================="<<endl;
	//cout<<"=============REMINDER: confidence split score NOT disconnected to study incoherences========================="<<endl;

	return 0;
}


//=============================================================================================================
int parseHierarchicalSegmentation2LineageHyperTree(hierarchicalSegmentation* hs,lineageHyperTree& lht)
{

	queue< TreeNode<nodeHierarchicalSegmentation>* > q;
	q.push( hs->dendrogram.pointer_mainRoot() );
	TreeNode<nodeHierarchicalSegmentation>* aux = NULL;
	//traverse tree top to bottom to generate lineages
	while( q.empty() == false )
	{
		aux = q.front();
		q.pop();

		if( aux->data.thrTau > hs->getMaxTau() )
		{

			//continue searching down the dendrogram for roots of lineages
			if(aux->left != NULL )
				q.push( aux->left );
			if(aux->right != NULL )
				q.push( aux->right );
			//NOTE: we do not write out objects that are clearly segmented (ie, isolated nodes that between minTau and maxTau have no connections in the dendrogram)
		}else{//generate a new lineage from this root

			lht.lineagesList.push_back( lineage() );
			list<lineage>::iterator listLineageIter = ((++ ( lht.lineagesList.rbegin() ) ).base());//iterator for the last element in the list

			queue< TreeNode<nodeHierarchicalSegmentation>* > qHS;
			queue< TreeNode<ChildrenTypeLineage>* > qLHT;//stores parent in LHT binary tree
			qHS.push( aux );
			qLHT.push( NULL );
			
			TreeNode<ChildrenTypeLineage>* auxLHT = NULL;
			int TM = 0;
			supervoxel sv;
			while( qHS.empty() == false )
			{
				//retrieve information
				aux = qHS.front();
				qHS.pop();
				auxLHT = qLHT.front();
				qLHT.pop();
				
				if( auxLHT == NULL )
					TM = 0;//root of dendrogram
				else
					TM = auxLHT->data->TM + 1;

				//generate supervoxel
				hs->supervoxelAtTreeNode(aux, sv);
				sv.TM = TM;

				//add supervoxel to lht
				lht.supervoxelsList[ TM ].push_back( sv );
				ChildrenTypeNucleus iterS = ((++ ( lht.supervoxelsList[ TM ].rbegin() ) ).base());//iterator for the last element in the list
				ParentTypeSupervoxel iterN =  lht.addNucleusFromSupervoxel(TM, iterS);//creates a nucleus from a single supervoxel and updates the hypergraph properly. Returns an iterator to the created element

				//add nucleus to lineage				
				iterN->treeNode.setParent(listLineageIter);
				if( auxLHT != NULL )
					listLineageIter->bt.SetCurrent( auxLHT );
				else //root
					listLineageIter->bt.SetMainRootToNULL();

				iterN->treeNodePtr = listLineageIter->bt.insert(iterN);
				if(iterN->treeNodePtr == NULL)
				{
					cout<<"ERROR: at parseHierarchicalSegmentation2LineageHyperTree: inerting node into binary tree"<<endl;
					return 3;
				}

				//update queues to keep descending in the lineage
				if( aux->left != NULL )
				{
					qHS.push( aux->left );
					qLHT.push( iterN->treeNodePtr );
				}
				
				if( aux->right != NULL )
				{
					qHS.push( aux->right );
					qLHT.push( iterN->treeNodePtr );
				}
			}
		}
	}

	return 0;
}

//=========================================================================================
//=============================================================================================================
int debugMergeSplitHStoLHT(hierarchicalSegmentation* hs,lineageHyperTree& lht, imgVoxelType tau,unsigned short* dataPtr)
{
	
	//generate segmentation
	hs->segmentationAtTau(tau);

	//keep a list of root nodes so I do not duplicates entries
	vector< TreeNode< nodeHierarchicalSegmentation >* > rootNodes;
	rootNodes.reserve( hs->currentSegmentatioSupervoxel.size() );

	//trim all supervoxels
	for(size_t ii = 0; ii < hs->currentSegmentatioSupervoxel.size(); ii++ )
	{
		hs->currentSegmentatioSupervoxel[ii].trimSupervoxel<unsigned short int>(dataPtr);
	}

	for(size_t ii = 0; ii < hs->currentSegmentatioSupervoxel.size(); ii++ )
	{
		//go all the way up by merging
		TreeNode< nodeHierarchicalSegmentation >* auxNode = hs->currentSegmentationNodes[ii], *rootMerge = hs->currentSegmentationNodes[ii];
		supervoxel auxSv;
		supervoxel rootMergeSv = hs->currentSegmentatioSupervoxel[ii];
		float score = 1.0f;

		while( score > 0.0f )
		{			
			auxSv = rootMergeSv;
			auxNode = rootMerge;
            // (ngc) not super sure dataptr is the right thing here?
            // did this when factoring dataPtr reference out of supervoxel.
            // Originaly dataPtr was passed to debugMergeSplitHStoLHT and I needed
            // to modify the following call.  I'm guessing everything was
            // refering to the same image.
			score = hs->suggestMerge<unsigned short int>( auxNode, auxSv, &rootMerge, rootMergeSv, dataPtr );	
		}

		//check if node was already done
		bool isDone = false;
		for(size_t jj = 0; jj < rootNodes.size(); jj++ )
		{
			if( rootNodes[jj] == auxNode )
			{
				isDone = true;
				break;
			}
		}
		if( isDone == true )
			continue;

		rootNodes.push_back( auxNode );

		//traverse all the ways down by splitting starting at auxSv
		//generate a new lineage from this root
		lht.lineagesList.push_back( lineage() );
		list<lineage>::iterator listLineageIter = ((++ ( lht.lineagesList.rbegin() ) ).base());//iterator for the last element in the list

		queue< TreeNode<nodeHierarchicalSegmentation>* > qHS;
		queue< TreeNode<ChildrenTypeLineage>* > qLHT;//stores parent in LHT binary tree
		qHS.push( auxNode );
		qLHT.push( NULL );
		TreeNode<nodeHierarchicalSegmentation>* rootSplit[2];
		supervoxel rootSplitSv[2];

		TreeNode<ChildrenTypeLineage>* auxLHT = NULL;
		int TM = 0;
		//void* dataPtr = hs->basicRegionsVec[0].dataPtr;
		supervoxel sv;
		while( qHS.empty() == false )
		{
			//retrieve information
			auxNode = qHS.front();
			qHS.pop();
			auxLHT = qLHT.front();
			qLHT.pop();

			if( auxLHT == NULL )
				TM = 0;//root of dendrogram
			else
				TM = auxLHT->data->TM + 1;

			//generate supervoxel and trim it
			hs->supervoxelAtTreeNode(auxNode, sv);
			sv.TM = TM;
			//sv.dataPtr = dataPtr;
			sv.trimSupervoxel<unsigned short int>(dataPtr);			


			//add supervoxel to lht
			lht.supervoxelsList[ TM ].push_back( sv );
			ChildrenTypeNucleus iterS = ((++ ( lht.supervoxelsList[ TM ].rbegin() ) ).base());//iterator for the last element in the list
			ParentTypeSupervoxel iterN =  lht.addNucleusFromSupervoxel(TM, iterS);//creates a nucleus from a single supervoxel and updates the hypergraph properly. Returns an iterator to the created element
			

			//add nucleus to lineage				
			iterN->treeNode.setParent(listLineageIter);
			if( auxLHT != NULL )
				listLineageIter->bt.SetCurrent( auxLHT );
			else //root
				listLineageIter->bt.SetMainRootToNULL();

			iterN->treeNodePtr = listLineageIter->bt.insert(iterN);
			if(iterN->treeNodePtr == NULL)
			{
				cout<<"ERROR: at parseHierarchicalSegmentation2LineageHyperTree: inserting node into binary tree"<<endl;
				return 3;
			}

			/*
			//debugging
			if( TM == 0 && lht.supervoxelsList[ TM ].size() == 2363 )
			{
				iterS->weightedCentroid<unsigned short int>();
				cout<<"ii = "<<ii<<";supervoxel = "<<(*iterS)<<";par tauThr = "<<auxNode->parent->data.thrTau<<endl;
				exit(3);
			}
			*/

			//calculate split and update queues to keep descending in the lineage
			hs->suggestSplit<unsigned short int>(auxNode, sv, rootSplit, rootSplitSv,dataPtr);
			if( rootSplit[0] != NULL )
			{
				qHS.push( rootSplit[0] );
				qLHT.push( iterN->treeNodePtr );
			}
			if( rootSplit[1] != NULL )
			{
				qHS.push( rootSplit[1] );
				qLHT.push( iterN->treeNodePtr );
			}
			
		}		
	}

	
	return 0;
}

//=====================================================================
TreeNode< ChildrenTypeLineage >* addSupervoxelsPointersFromLineage(std::vector< vector<supervoxel*> > &svIniVec, int iniTM, int endTM, TreeNode< ChildrenTypeLineage >* node)
{
	//add elements for node
	TreeNode< ChildrenTypeLineage >* auxNode = node;
	while( auxNode->parent != NULL && auxNode->data->TM > iniTM )
		auxNode = auxNode->parent;
	//traverse the node downstream
	queue< TreeNode< ChildrenTypeLineage >* > q;
	q.push( auxNode );
	TreeNode< ChildrenTypeLineage >* root = auxNode;
	while( q.empty() == false )
	{
		auxNode = q.front();
		q.pop();
		if( auxNode->left != NULL && auxNode->left->data->TM <= endTM )
			q.push( auxNode->left );
		if( auxNode->right != NULL && auxNode->right->data->TM <= endTM )
			q.push( auxNode->right );

		//add all supervoxels
		int offsetTM = auxNode->data->TM - iniTM;
		for( vector< ChildrenTypeNucleus >::iterator iterS = auxNode->data->treeNode.getChildren().begin(); iterS != auxNode->data->treeNode.getChildren().end(); ++iterS )
		{
			svIniVec[ offsetTM ].push_back( &(*(*iterS)) );
		}
	}

	return root;
}


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

//===========================================================
void transposeStackUINT16(mylib::Array *img)
{
	//cout<<"WARNING: transposing each slice to agree with Matlab convention"<<endl;
	if(img->type != mylib::UINT16_TYPE)
	{
		cout<<"ERROR: at transposeStackUINT16: array has to be uint16"<<endl;
		exit(3);
	}
	if(img->ndims != 3)
	{
		cout<<"ERROR: at transposeStackUINT16: code only ready for 3D arrays"<<endl;
		exit(3);
	}

	mylib::Size_Type sliceSize = img->dims[0]*img->dims[1];
	mylib::uint16* imgPtr = (mylib::uint16*) (img->data);
	mylib::uint16* imgCpyPtr = new mylib::uint16[sliceSize]; 

	mylib::Size_Type imgIdx = 0;
	for(mylib::Dimn_Type zz = 0; zz<img->dims[2]; zz++)
	{
		//copy plane
		memcpy(imgCpyPtr,&(imgPtr[imgIdx]),sizeof(mylib::uint16)*sliceSize);

		for(mylib::Dimn_Type xx = 0; xx<img->dims[0]; xx++)		
		{
			mylib::Size_Type offset2 = xx;
			for(mylib::Dimn_Type yy = 0; yy<img->dims[1]; yy++)
			{
				imgPtr[imgIdx++] = imgCpyPtr[offset2];
				offset2 += img->dims[0];
			}
		}
	}

	//change array dimensions
	mylib::Dimn_Type aux = img->dims[0];
	img->dims[0] = img->dims[1];
	img->dims[1] = aux;

	delete[] imgCpyPtr;
}

//==============================================================
int extendDeadNucleiAtTMwithHS(lineageHyperTree &lht, hierarchicalSegmentation* hsForward, int TM, int& numExtensions, int &numDeaths,float* imgPtr)
{
	
	numExtensions = 0;
	numDeaths = 0;
	bool dont_repeat = true;

	if( TM < 0 || TM >= (int) ( lht.getMaxTM()) )
		return 0;

	//TreeNode<ChildrenTypeLineage>* aux;
	for ( int ii=0; ii<2; ii++ ) //only allow up to 2 repetitions of iterating nucleus list
	{
		for(list<nucleus>::iterator iterN = lht.nucleiList[TM].begin(); iterN != lht.nucleiList[TM].end(); ++iterN)
		{
			if( iterN->treeNodePtr->getNumChildren() == 0)
			{
				if( iterN->TM != TM ) //if this is not true we will probably 
				{
					cout << "DEBUG Sv: Nucleus with treeNodePtr " << (void *)(iterN->treeNodePtr) << " has timepoint " << iterN->TM << ", that doesn't match " << TM << "!" << endl;
					continue;
				}
				numDeaths++;
				numExtensions += extendDeadNucleiWithHS(lht, hsForward, iterN->treeNodePtr,imgPtr,dont_repeat);
			}
		}
		
		if (dont_repeat) //stop repeating the time point when we are fully done with lineage slicing/dicing
			break;
	}
	
	return 0;
}

//=============================================================================================================================
int extendDeadNucleiWithHS(lineageHyperTree &lht, hierarchicalSegmentation* hsForward, TreeNode<ChildrenTypeLineage>* rootDead,float*imgPtr,bool &dont_repeat)
{
	dont_repeat = true; //let calling function know we don't need to iterate on the nuclei list from this timepoint, by default
	
	if(rootDead == NULL || rootDead->getNumChildren() != 0)
		return 0;//rootDead is not a dead split so we cannnot do anything
	//try to find the most obvious continuation
	
	//1.-Generate a super-supervoxel by merging all the supervoxel belonging to the nucleus
	vector<ChildrenTypeNucleus>::iterator iterS = rootDead->data->treeNode.getChildren().begin();
	supervoxel supervoxelFromNucleus( *(*iterS) );
	++iterS;
	vector< supervoxel* > auxVecS;
	for(;  iterS != rootDead->data->treeNode.getChildren().end(); ++iterS)
	{
		auxVecS.push_back( &(*(*iterS)) );
	}
	supervoxelFromNucleus.mergeSupervoxels(auxVecS);//we add all the components at once

	//2.-find candidate with largest intersection
	uint64 intersectionSize = 0, auxI;
	SibilingTypeSupervoxel intersectionS;
	for(iterS = rootDead->data->treeNode.getChildren().begin();  iterS != rootDead->data->treeNode.getChildren().end(); ++iterS)
	{
		for(vector< SibilingTypeSupervoxel >::iterator iterS2 = (*iterS)->nearestNeighborsInTimeForward.begin(); iterS2 != (*iterS)->nearestNeighborsInTimeForward.end(); ++iterS2)
		{
			auxI = supervoxelFromNucleus.intersectionSize( *(*(iterS2)) );
			if(  auxI > intersectionSize )
			{
				intersectionSize = auxI;
				intersectionS = (*iterS2);
			}
		}
	}
	
	if( intersectionSize == 0 )
		return 0;//no clear option for extending death
	//3.find the nuclei that "owns" the supervoxel
	list< nucleus >::iterator iterNucOwner, iterNucOwnerDaughterL, iterNucOwnerDaughterR, iterNucNew;
	int intIsCellDivision = 0x000;//0->no children;0x0001->left children;0x0010->right children;0x0011->both children
	if( intersectionS->treeNode.hasParent() == false )//the supervoxel with highest intersection has not been claimed by anybody->just take it
	{
		iterNucNew = lht.addNucleusFromSupervoxel(rootDead->data->TM + 1, intersectionS );//returns iterator to newly created nucleus
		//update lineage-nucleus hypergraph
		iterNucNew->treeNode.setParent( rootDead->data->treeNode.getParent() );
		rootDead->data->treeNode.getParent()->bt.SetCurrent( rootDead );
		iterNucNew->treeNodePtr = rootDead->data->treeNode.getParent()->bt.insert( iterNucNew );
		if( iterNucNew->treeNodePtr == NULL )
			exit(3);

		return 1;//we have added one nucleus
	}else{
		
		iterNucOwner = intersectionS->treeNode.getParent();
		//cout << "    R3: " << (void *)(iterNucOwner->treeNodePtr) << "/" << (void *)(iterNucOwner->treeNodePtr->data->treeNodePtr) << endl;
		 if( iterNucOwner->treeNodePtr->parent != NULL )
		 {
			//cout << "     R3parent: " << (void *)(iterNucOwner->treeNodePtr->parent) << "/" << (void *)(iterNucOwner->treeNodePtr->parent->data->treeNodePtr) << endl;
			iterNucOwner = iterNucOwner->treeNodePtr->parent->data; //we were one time step ahead, so now look from the perspective of the parent of the nucelus who owns the supervoxels
			//cout << "    R4: " << (void *)(iterNucOwner->treeNodePtr) << endl;
			 if( iterNucOwner->treeNodePtr->getNumChildren() > 1)
			 {
				 //cout << "    R5a!" << endl;
				 intIsCellDivision = 0x0011;
				 iterNucOwnerDaughterR = iterNucOwner->treeNodePtr->right->data;
				 //cout << "    R5aa!" << endl;
				 iterNucOwnerDaughterL = iterNucOwner->treeNodePtr->left->data;
				 //cout << "    R5ab!" << endl;
			 }else{
				 //cout << "    R5b!" << endl;
				 if ( iterNucOwner->treeNodePtr->left !=NULL )
				 {
					 //cout << "    R5ba!" << endl;
					 intIsCellDivision = 0x0001;
					 iterNucOwnerDaughterL = iterNucOwner->treeNodePtr->left->data;
					 iterNucOwnerDaughterR = iterNucOwnerDaughterL;
				 }else{
					 //cout << "    R5bb!" << endl;
					 intIsCellDivision = 0x0010;
					 iterNucOwnerDaughterR = iterNucOwner->treeNodePtr->right->data;
					 iterNucOwnerDaughterL = iterNucOwnerDaughterR;
				 }
			 }			 

		 }else{
			 //return 0;//there is no parent
			 //2021 new: if there is no parent, attach dead cell to this one to connect the lineages
			//update lineage-nucleus hypergraph
			//we need to make sure all the the elements in the binary tree point to the correct lineage
			TreeNode<ChildrenTypeLineage> *aux = iterNucOwner->treeNodePtr;
			list< lineage >::iterator iterL1 = iterNucOwner->treeNode.getParent();
			queue< TreeNode<ChildrenTypeLineage>* > q;
			q.push( iterNucOwner->treeNodePtr );
			
			//if( aux->left != NULL) q.push(aux->left);
			//if( aux->right != NULL) q.push(aux->right);
			
			while( !q.empty() )
			{
				aux = q.front();
				q.pop();
				if( aux->left != NULL) q.push(aux->left);
				if( aux->right != NULL) q.push(aux->right);
				aux->data->treeNode.setParent( rootDead->data->treeNode.getParent() );
			}
			//now remove the daughter tree since it has been attached to rootDead
			//iterNucOwner->treeNode.getParent()->bt.SetMainRootToNULL();
			//lht.lineagesList.erase( iterNucOwner->treeNode.getParent() );
			iterL1->bt.SetMainRootToNULL();
			lht.lineagesList.erase( iterL1 );
			
			//last step is set correct binary tree and parent nucleus for first linked (daugther) nucleus
			//iterNucOwner->treeNode.setParent( rootDead->data->treeNode.getParent() );
			//rootDead->data->treeNode.getParent()->bt.reset();
			iterNucOwner->treeNodePtr->parent = rootDead;
			rootDead->left = iterNucOwner->treeNodePtr;
			/*if(iterNucOwner->treeNodePtr == NULL)
				cout << "   PROBLEM: extendDeadNuclei, a new birth used to extend a dead cell has NULL treeNodePtr!" << endl;
			if(iterNucOwner->treeNodePtr->parent == NULL)
				cout << "   PROBLEM: extendDeadNuclei, a new birth used to extend a dead cell has NULL parent!" << endl;
			if(rootDead != rootDead->data->treeNodePtr )
				cout << "   PROBLEM: extendDeadNuclei, rootDead treeNodePtr incorrect!" << endl;
			*/
			//cout<<"   DEBUG: lineageHyperTree::extendDeadNuclei: extended a dead cell at TM " << (void *)(rootDead) << " by attaching to a new birth at TM " << (void *)(iterNucOwner->treeNodePtr) << " with " << iterNucOwner->treeNodePtr->getNumChildren() << " children." << endl;
			
			return 1;
		 }
	}

	if( iterNucOwner->TM != rootDead->data->TM)
	{
		cout<<"ERROR: lineageHyperTree::extendDeadNuclei: TM does not agree between two candidate nucleus"<<endl;
		exit(5);
	}
	//cout << "   S0." << endl;
	//4.-run a small Hungarian algorithm in order to decide what is the best matching to solve this issue.
	//We setup a list of cancidates in time point t
	list< supervoxel > svListT0;	
	for(iterS = rootDead->data->treeNode.getChildren().begin();  iterS != rootDead->data->treeNode.getChildren().end(); ++iterS)
	{
		svListT0.push_back( *(*iterS) );
	}
	for(iterS = iterNucOwner->treeNode.getChildren().begin();  iterS != iterNucOwner->treeNode.getChildren().end(); ++iterS)
	{
		svListT0.push_back( *(*iterS) );
	}


	//4.1-Addition with respect to original code in lineageHypertree class
	//Check in supervoxel cancidates in t+1 if splitting them should help. If it does ->do it before calculating Hungarian algorithm
	//cout << "   S1." << endl;
	float ddNoSplit;
	TreeNode< nodeHierarchicalSegmentation >* rootSplit[2];
	supervoxel rootSplitSv[2];
	float ddSplit[2];
	//items needed to know correspondence between hsForward->currentSegmentation and lht.supervoxelsList[rootDead->data->TM + 1]
	const int TMforward = rootDead->data->TM + 1;
	vector<ChildrenTypeNucleus> svListIterVec;
	lht.getSupervoxelListIteratorsAtTM( svListIterVec, TMforward );
	size_t countNN = 0;//to keep id
	for(vector<ChildrenTypeNucleus>::iterator iterNN = svListIterVec.begin(); iterNN != svListIterVec.end(); ++iterNN, countNN++)
		(*iterNN)->tempWildcard = countNN;
	//cout << "   S2." << endl;
	for(list< supervoxel >::iterator iter = svListT0.begin(); iter != svListT0.end(); ++iter)
	{
		/*size_t sizeVec = iter->nearestNeighborsInTimeForward.size();
		for( size_t cc = 0; cc < sizeVec; cc++ )//we need to access using indexes because we push_back into iter->nearestNeighborsInTimeForward if a split is done
		{
			SibilingTypeSupervoxel iterNeigh = iter->nearestNeighborsInTimeForward[cc];*/
		
		//reverse iterate because using push_back on this very vector
		if ( iter->TM + 1 != TMforward )
		{
			cout << "extendDeadNucleiWithHS: supervoxels time frame not matched, " << TMforward << " != " << iter->TM + 1 << endl;
			continue;
		}
		vector<SibilingTypeSupervoxel>::iterator iterNeighEnd = iter->nearestNeighborsInTimeForward.end(); //set this up before, so we don't re-run elements that have been push_back'd
		for(vector<SibilingTypeSupervoxel>::iterator iterNeighThis = iter->nearestNeighborsInTimeForward.begin(); iterNeighThis != iterNeighEnd; ++iterNeighThis) 
		{
			SibilingTypeSupervoxel iterNeigh = (*iterNeighThis);
			if ( (*iterNeigh).TM != TMforward )
			{
				cout << "extendDeadNucleiWithHS: splitting supervoxels time frame not matched, " << (*iterNeigh).TM << " != " << iter->TM + 1 << endl;
				continue;
			}
			ddNoSplit = supervoxelFromNucleus.JaccardDistance(*iterNeigh);//calculate distance without splitting
			float ddThr = ddNoSplit - 0.1;//it has to be significantly better

			//propose split
			float prob = hsForward->suggestSplit<float>((*iterNeigh).nodeHSptr, (*iterNeigh), rootSplit, rootSplitSv,imgPtr);//probability of the split
			if( prob < 1e-3 )
				continue;
			/* DEBUG
			if ((*iterNeigh).treeNode.hasParent())
				cout << "   S3 parent nucleus: " << (void *)((*iterNeigh).treeNode.getParent()->treeNodePtr) << endl;
			else
				cout << "   S3 no parent nucleus." << endl;
				*/
			//calculate Jaccard distance
			ddSplit[0] = supervoxelFromNucleus.JaccardDistance(rootSplitSv[0]);
			ddSplit[1] = supervoxelFromNucleus.JaccardDistance(rootSplitSv[1]);
			
			if( ddSplit[0] < ddThr || ddSplit[1] < ddThr )//incorporate split into solution space
			{
				
				
				//if matching is significantly better->proceed with to split supervoxels and incorporate them in the lineage
				//calculate Gauss stats
				/*DEBUG
				if (rootSplitSv[0].treeNode.hasParent())
					cout << "   S4a parent nucleus[0]: " << (void *)(rootSplitSv[0].treeNode.getParent()->treeNodePtr) << endl;
				else
					cout << "   S4a no parent nucleus[0]." << endl;
				
				if (rootSplitSv[1].treeNode.hasParent())
					cout << "   S4a parent nucleus[1]: " << (void *)(rootSplitSv[1].treeNode.getParent()->treeNodePtr) << endl;
				else
					cout << "   S4a no parent nucleus[1]." << endl;
				*/
				rootSplitSv[0].weightedGaussianStatistics<float>(true,imgPtr);
				rootSplitSv[1].weightedGaussianStatistics<float>(true,imgPtr);

				//update sv-nucleus info
				if ( (*iterNeigh).treeNode.getParent()->TM != TMforward )
				{
					cout << "extendDeadNucleiWithHS: nucleus-supervoxel time mismatch, " << (*iterNeigh).TM << " != " << (*iterNeigh).treeNode.getParent()->TM << endl;
					continue;
				}				
				rootSplitSv[0].treeNode.setParent((*iterNeigh).treeNode.getParent());
				rootSplitSv[1].treeNode.setParent((*iterNeigh).treeNode.getParent());

				/*DEBUG
				if ( rootSplitSv[1].TM != TMforward )
				{
					cout << "extendDeadNucleiWithHS: Supervoxel being added to HS, time mismatch, " << rootSplitSv[1].TM << " != " << TMforward << ", for nucleus with treeNodePtr " << (void *)((*iterNeigh).treeNode.getParent()->treeNodePtr) << endl;
					continue;
				} */
				//lht->supervoxelsList[TM] was created as a copy of hs->currentSegmentation for backwards compatibility (even if it is redundant). So we just need tokeep that structure								
				int aa = (int) ((*iterNeigh).tempWildcard);
				hsForward->currentSegmentatioSupervoxel[aa] = rootSplitSv[0];//root split has the correct nodeHSptr
				hsForward->currentSegmentatioSupervoxel.push_back( rootSplitSv[1] );//this could dynamically reallocate memory and make all the pointers nodeHSptr->data.svPtr invalid (the crash would be obvious). So I reserve when creating a partition.
				
				hsForward->currentSegmentationNodes[aa] = rootSplit[0];
				hsForward->currentSegmentationNodes.push_back( rootSplit[1] );
			

				//update nucleus-sv info (remember: lht->supervoxelsList[TM] was created as a copy of hs->currentSegmentation for backwards compatibility);
				//I actually only need to change nucleus-supervoxel by adding children to nucleus since it was a split;
				(*(svListIterVec[aa])) = rootSplitSv[0];
				(svListIterVec[aa])->tempWildcard = aa;
				lht.supervoxelsList[TMforward].push_back(rootSplitSv[1]);
				SibilingTypeSupervoxel iterSadded = (++ ( lht.supervoxelsList[TMforward].rbegin() ) ).base();//get iterator to added element
				iterSadded->tempWildcard = lht.supervoxelsList[TMforward].size() - 1;
				(*iterNeigh).treeNode.getParent()->treeNode.addChild( iterSadded ); //CRASH here
				/*if ( iterSadded->TM != TMforward )
				{
					cout << "extendDeadNucleiWithHS: iterSadded being added, time mismatch, " << iterSadded->TM << " != " << TMforward << ", for nucleus with treeNodePtr " << (void *)((*iterNeigh).treeNode.getParent()->treeNodePtr) << endl;
					exit( 5);
				}*/				
				//add new supervoxels in iter->NEARESTneighbors list for Hungarian algorithm;
				//added to FRONT so iterator doesn't return to it
				iter->nearestNeighborsInTimeForward.push_back( iterSadded );
			}
		}
	}

	list< supervoxel> nullAssignmentList;//temporary supervoxel list to simulate null assignment
	nullAssignmentList.push_back(supervoxel());
	SibilingTypeSupervoxel nullAssignment = nullAssignmentList.begin();
	nullAssignment->centroid[0] = -1e32f;//characteristic to find out no assignment
	vector<SibilingTypeSupervoxel> assignmentId;
	int err = calculateTrackletsWithSparseHungarianAlgorithm(svListT0, 0, 0.9, assignmentId, &nullAssignment);
	if( err > 0)
		exit(err);

	
	//5.-Parse results and modify assignment accordingly
	int extendedLineages = 0;//1->we have extended it
	int count =0;
	for(iterS = rootDead->data->treeNode.getChildren().begin();  iterS != rootDead->data->treeNode.getChildren().end(); ++iterS, ++count)
	{
		if( assignmentId[count]->centroid[0] < 0.0f )
			continue;//not assigned to anything
			
		if ( assignmentId[count]->TM != TMforward )
		{
			cout << "extendDeadNucleiWithHS: assignmentId[count]->TM, time mismatch, " << assignmentId[count]->TM << " != " << TMforward << ", for attaching to nucleus with treeNodePtr " << (void *)(rootDead) << endl;
			continue;
		}			

		if( assignmentId[count]->treeNode.hasParent() == false )//the assigned element has no parent-> we can claim it directly
		{
			if( extendedLineages == 0)//we need to create new nucleus in the list
			{
				lht.nucleiList[ TMforward ].push_back( nucleus(TMforward, assignmentId[count]->centroid) );
				iterNucNew = (++ ( lht.nucleiList[ TMforward ].rbegin() ) ).base();//iterator to last added nucleus
				
				//update lineage-nucleus hypergraph
				iterNucNew->treeNode.setParent( rootDead->data->treeNode.getParent() );
				rootDead->data->treeNode.getParent()->bt.SetCurrent( rootDead );
				iterNucNew->treeNodePtr = rootDead->data->treeNode.getParent()->bt.insert( iterNucNew );
				if( iterNucNew->treeNodePtr == NULL )
					exit(3);
				//update supervoxel-nucleus hypergraph
				iterNucNew->addSupervoxelToNucleus( assignmentId[count] );
				assignmentId[count]->treeNode.setParent( iterNucNew );

				extendedLineages++;
				//iterNucNew->confidence = 4;//to analyze extended elements
				
				//cout << "    T1: " << (void *)(iterNucNew->treeNodePtr) << " with parent " << (void *)(iterNucNew->treeNodePtr->parent) << "/" << (void *)(iterNucNew->treeNodePtr->parent->data->treeNodePtr) << endl;
			}
			//cout << "    T2" << endl;
		}
		else if( ( assignmentId[count]->treeNode.getParent() == iterNucOwnerDaughterL ) || ( assignmentId[count]->treeNode.getParent() == iterNucOwnerDaughterR ) )//to confirm it is not null assignment && we are "stealing" a supervoxel from iterNucOwner and not from anothe nuclei
		{
			if( extendedLineages == 0)//we need to create new nucleus in the list
			{
				lht.nucleiList[ TMforward ].push_back( nucleus(TMforward, assignmentId[count]->centroid) );
				iterNucNew = (++ ( lht.nucleiList[ TMforward ].rbegin() ) ).base();//iterator to last added nucleus
				
				//update lineage-nucleus hypergraph
				iterNucNew->treeNode.setParent( rootDead->data->treeNode.getParent() );
				rootDead->data->treeNode.getParent()->bt.SetCurrent( rootDead );
				iterNucNew->treeNodePtr = rootDead->data->treeNode.getParent()->bt.insert( iterNucNew );
				if( iterNucNew->treeNodePtr == NULL )
					exit(3);
				extendedLineages++;
				//iterNucNew->confidence = 4;//to analyze extended elements
				
				//cout << "    T3: " << (void *)(iterNucNew->treeNodePtr) << " with parent " << (void *)(iterNucNew->treeNodePtr->parent) << "/" << (void *)(iterNucNew->treeNodePtr->parent->data->treeNodePtr) << endl;
			}
			//update supervoxel-nucleus hypergraph
			if ( assignmentId[count]->treeNode.getParent() == iterNucOwnerDaughterL )
			{
				iterNucOwnerDaughterL->removeSupervoxelFromNucleus( assignmentId[count] );
				//if( ret > 0 )
				//	cout<<"WARNING: lineageHyperTree::extendDeadNuclei: supervoxel not found to be removed from nucleus"<<endl;
			}
			else{
				iterNucOwnerDaughterR->removeSupervoxelFromNucleus( assignmentId[count] );
				//if( ret > 0 )
				//	cout<<"WARNING: lineageHyperTree::extendDeadNuclei: supervoxel not found to be removed from nucleus"<<endl;
			}
			iterNucNew->addSupervoxelToNucleus( assignmentId[count] );
			assignmentId[count]->treeNode.setParent( iterNucNew );
			//cout << "    T4" << endl;
		}
	}


	//make sure original nuclei still has some supervoxels associated
	if( ( (intIsCellDivision & 0x0001) != 0 ) && ( iterNucOwnerDaughterL->treeNode.getNumChildren() == 0 ) )
	{
		int TMaux = iterNucOwnerDaughterL->TM;
		
		//cout << "    U1: deleting daughter " << (void *)(iterNucOwnerDaughterL->treeNodePtr) << " with parent " << (void *)(iterNucOwner->treeNodePtr) << "/" << (void *)(iterNucOwnerDaughterL->treeNodePtr->parent) << "/" << (void *)(iterNucOwnerDaughterL->treeNodePtr->parent->data->treeNodePtr) << endl;
		
		//before removing nucleus, any children need to be reattached to the new nucleus and rootDead's lineage
		TreeNode<ChildrenTypeLineage> *aux = iterNucOwnerDaughterL->treeNodePtr;
		queue< TreeNode<ChildrenTypeLineage>* > q;
			
		if( aux->left != NULL)
		{
			aux->left->parent = iterNucNew->treeNodePtr;
			iterNucNew->treeNodePtr->left = aux->left;
			aux->left = NULL;
			q.push(iterNucNew->treeNodePtr->left);
		}
		if( aux->right != NULL)
		{
			aux->right->parent = iterNucNew->treeNodePtr;
			iterNucNew->treeNodePtr->right = aux->right;
			aux->right = NULL;
			q.push(iterNucNew->treeNodePtr->right);
		}
			
		while( !q.empty() )
		{
			aux = q.front();
			q.pop();
			if( aux->left != NULL) q.push(aux->left);
			if( aux->right != NULL) q.push(aux->right);
			aux->data->treeNode.setParent( rootDead->data->treeNode.getParent() );
		}

		//remove this nucleus
		delete iterNucOwnerDaughterL->treeNodePtr;
		iterNucOwner->treeNodePtr->left = NULL;
		lht.nucleiList[ TMaux ].erase( iterNucOwnerDaughterL );
		
		dont_repeat = false; //see if parent (now detached from daughter) is now dead and needs to find a path forward
	}
	if(  ( (intIsCellDivision & 0x0010) != 0 ) && ( iterNucOwnerDaughterR->treeNode.getNumChildren() == 0 ) )
	{
		int TMaux = iterNucOwnerDaughterR->TM;
		
		//cout << "    U2: deleting daughter " << (void *)(iterNucOwnerDaughterR->treeNodePtr) << " with parent " << (void *)(iterNucOwner->treeNodePtr) << "/" << (void *)(iterNucOwnerDaughterR->treeNodePtr->parent) << "/" << (void *)(iterNucOwnerDaughterR->treeNodePtr->parent->data->treeNodePtr) << endl;
		
		//before removing nucleus, any children need to be reattached to the new nucleus and rootDead's lineage
		TreeNode<ChildrenTypeLineage> *aux = iterNucOwnerDaughterR->treeNodePtr;
		queue< TreeNode<ChildrenTypeLineage>* > q;
			
		if( aux->left != NULL)
		{
			aux->left->parent = iterNucNew->treeNodePtr;
			iterNucNew->treeNodePtr->left = aux->left;
			aux->left = NULL;
			q.push(iterNucNew->treeNodePtr->left);
		}
		if( aux->right != NULL)
		{
			aux->right->parent = iterNucNew->treeNodePtr;
			iterNucNew->treeNodePtr->right = aux->right;
			aux->right = NULL;
			q.push(iterNucNew->treeNodePtr->right);
		}
			
		while( !q.empty() )
		{
			aux = q.front();
			q.pop();
			if( aux->left != NULL) q.push(aux->left);
			if( aux->right != NULL) q.push(aux->right);
			aux->data->treeNode.setParent( rootDead->data->treeNode.getParent() );
		}		
		
		delete iterNucOwnerDaughterR->treeNodePtr;
		iterNucOwner->treeNodePtr->right = NULL;
		lht.nucleiList[ TMaux ].erase( iterNucOwnerDaughterR );
		
		dont_repeat = false; //see if parent (now detached from daughter) is now dead and needs to find a path forward
	}	

	//cout << "   V! " << extendedLineages << endl;
	return extendedLineages;//returns 1 if extension was achieved
}


//=====================================================================
int mkDirectory(const std::string& folder)
{

#if defined(_WIN32) || defined(_WIN64)
		if (GetFileAttributes(folder.c_str()) == INVALID_FILE_ATTRIBUTES)//check if folder exists
		{

			bool errB = CreateDirectory(folder.c_str(), NULL);
			if (errB == false)
			{
				cout << "ERROR: mkDirectory" << folder << ": generating folder " << folder << endl;
				return 1;
			}
		}

#else
		struct stat St;
		if (stat(folder.c_str(), &St) != 0)//check if folder exists
		{
			string cmd("mkdir " + folder);
			int error = system(cmd.c_str());
			if (error>0)
			{
				cout << "ERROR (" << error << "): generating path " << folder << endl;
				cout << "Wtih command " << cmd << endl;
				return error;
			}
		}
#endif
	return 0;
}

//============================================================================
template int parseNucleiList2TGMM<float>(std::vector<GaussianMixtureModel*> &vecGM, lineageHyperTree &lht,int frame, bool regularizeW, float thrDist2LargeDisplacement,float *imgPtr);
template int parseNucleiList2TGMM<unsigned short int>(std::vector<GaussianMixtureModel*> &vecGM, lineageHyperTree &lht,int frame, bool regularizeW, float thrDist2LargeDisplacement,unsigned short *imgPtr);
template int parseNucleiList2TGMM<unsigned char>(std::vector<GaussianMixtureModel*> &vecGM, lineageHyperTree &lht,int frame, bool regularizeW, float thrDist2LargeDisplacement,unsigned char *imgPtr);
