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
/*
 * Copyright (C) 2020-2021 Martin Dominguez
 * additional code for FTGMM project
 * 
 * Gladstone Institutes
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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

//2021 NEW: add parameter for thrCellDivisionPlaneDistance and use this function to find several division mother/sibling pairs to audition for each new track
int cellDivisionTemporalWindow_TGMMsupport::classifyCellDivisionTemporalWindow(lineageHyperTree& lht, int frame, vector<mylib::Array*>& imgVec, int devCUDA, double thrCellDivisionPlaneDistance, float* im_zero, float* im_plus_one, bool regularize_W4DOF, float scaleOrig[3] )

{
	//cout << "===========================DEBUGGING classifyCellDivisionTemporalWindow, TM: " << frame << "====================" << endl;
	//this functions follows mainCellDivisionClassifierWithTemporalWindow::mainSingleWindowForDaughters( int argc, const char** argv )
	basicEllipticalHaarFeatureVector::useDoGfeatures = true;
	cellDivisionWithTemporalWindow::setUseFixPrecisionMatrix(true);
	
	//float scale[dimsImage];
	//supervoxel::getScale(scale);

	list<nucleus>* nucleiList = &(lht.nucleiList[frame+1]);
	list<nucleus>* nucleiList_minus_1 = &(lht.nucleiList[frame]);

	if (nucleiList->empty() == true)
		return 0;
	if (nucleiList_minus_1->empty() == true)
		return 0;

	
	//2021 NEW: each new track may have several prospective divisions that are tested by the classifier, so don't just return
	//if (thrCDWT <= 0.0f && writeCDWTfeatures() == false)//no point in calculating anything. All cell divisions are accepted as valid
	//	return 0;

	int numTrueCellDivisions = 0;
	const float test_threshold = (float) thrCellDivisionPlaneDistance;
	//only one common divisionnodes and parentnodes list, with divisions being broken here and not in any call to breakCellDivisionBasedOnCellDivisionPlaneConstraint
	//find all the elements that contain a division
	vector<TreeNode< ChildrenTypeLineage >* > divisionNodes;//contains a pointer to each of the divisions in a specific time point -- stores originally-classified divisions, by CellDivison.cpp
	vector<TreeNode< ChildrenTypeLineage >* > daughterNodes;//rather than classifying divisions by parent, will chose furthest daughter and multiple prospective parents to choose from
	vector<TreeNode< ChildrenTypeLineage >* > parentNodes;//contains a pointer to closest prospective parent for each of divisionNodes

	divisionNodes.reserve(nucleiList->size() / 5);
	daughterNodes.reserve(nucleiList->size() / 4);
	parentNodes.reserve(nucleiList->size() / 4);
	
	//2021 NEW: push new tracks, not division parents (there won't be division parents at this point for this timepoint)
	float d, d_this, sv_vol, e, half_daughter_distance; //d stores division midplane distance, and e stores absolute distance between daughters
	float min_1; TreeNode< ChildrenTypeLineage >* par_1;
	float min_2; TreeNode< ChildrenTypeLineage >* par_2;
	float min_3; TreeNode< ChildrenTypeLineage >* par_3;
	TreeNode< ChildrenTypeLineage >* mother;
	TreeNode< ChildrenTypeLineage >* left_daughter;
	TreeNode<ChildrenTypeLineage> *aux, *auxD;
	
	/*for (list<nucleus>::iterator iterN = nucleiList->begin(); iterN != nucleiList->end(); ++iterN)
	{
		if (iterN->treeNodePtr->getNumChildren() == 2)
			divisionNodes.push_back(iterN->treeNodePtr);
	}*/
	
	for (list<nucleus>::iterator iterN = nucleiList_minus_1->begin(); iterN != nucleiList_minus_1->end(); ++iterN)
	{
		if (iterN->treeNodePtr->getNumChildren() == 2) //found pre-called division, now have to break it, identify furthest daugther, then find its connected supervoxels and possible siblings (and their parents)
		//if (iterN->treeNodePtr->parent == NULL && iterN->treeNode.getParent()->bt.pointer_mainRoot() == iterN->treeNodePtr ) //detect new tracks
		{
			aux = iterN->treeNodePtr;
			
			if ( aux->left->data->TM != aux->right->data->TM || iterN->TM + 1 != aux->left->data->TM ) //parents and daughters have to be one timeframe apart for this algorithm; if they aren't don't even break link, just leave it
			{
				//cout << "  Timeframe discrepancy: parent-" << iterN->TM << ", left-" << aux->left->data->TM << ", right-" << aux->right->data->TM << endl;
				continue;
			}
			//check which one is thefurthest daughter
			double dL = 0, dR = 0;
			for(int ii = 0; ii < dimsImage; ii++)
			{
				dL += pow(aux->left->data->centroid[ii] - iterN->centroid[ii], 2) * scaleOrig[ii]; //scale[ii]
				dR += pow(aux->right->data->centroid[ii] - iterN->centroid[ii],2) * scaleOrig[ii]; //scale[ii]
			}

			if( dR > dL )//disconnect right daughter
			{
				auxD = aux->right;
				//remove link in mother
				aux->right = NULL;								
			}else{ //disconnect left daughter
				auxD = aux->left;
				//remove link in mother
				aux->left = NULL;				
			}
			//remove link in daughter
			auxD->parent = NULL;
			
			//add the now-broken division to Santa's naughty list
			divisionNodes.push_back(aux);
			divisionNodes.push_back(auxD);			
			//cout << "auxD: " ;
			
			cout<<"   DEBUG: TGMMsupport at frame " << frame << ", just broke " << (void *)(auxD) << " from mother " << (void *)(aux) << endl;
		}
	}
	
	//get ready for below loop
	bool added_original_trio;
	bool added_any_trio;
	size_t ss = 0;
	
	//revisit broken division nodes and find prospective trios including the disconnected daughter
	for (int dd=0; dd<divisionNodes.size(); dd+=2 )
	{	
		aux = divisionNodes[dd];
		auxD = divisionNodes[dd+1];
		
		/*DEBUG
		if ( aux->data->TM + 1 != auxD->data->TM ) //parents and daughters have to be one timeframe apart for this algorithm; if they aren't don't even break link, just leave it
			cout << "  " << dd << " Timeframe discrepancy: parent-" << aux->data->TM << ", prospective daughter-" << auxD->data->TM << endl;
		if ( aux->getNumChildren() != 1 )
			cout << "  " << dd << " aux->getNumChildren: " << aux->getNumChildren() <<endl;
		if ( auxD->parent != NULL )
			cout << "  " << dd << " auxD->parent not null, in TM: " << auxD->parent->data->TM <<endl;
		//END DEBUG*/		
		
		added_original_trio = false;
		added_any_trio = false;
		
		par_1 = NULL;
		par_2 = NULL;
		par_3 = NULL;
		
		min_1 = test_threshold + 1e-6;
		min_2 = test_threshold + 1e-6;
		min_3 = test_threshold + 1e-6;	

					
		//find the closest three prospective siblings for this new track, then push parents onto running list
		for (list<nucleus>::iterator iterN_2 = nucleiList->begin(); iterN_2 != nucleiList->end(); ++iterN_2)
		{
			mother = iterN_2->treeNodePtr->parent;
			
			if ( mother == NULL || mother->getNumChildren() > 1 )
			{
				//cout << "  iterN_2 parent is NULL" << endl; //DEBUG
				continue;					
			}
			
			if ( mother->getNumChildren() != 1 )
			{
				cout << "   DEBUG: iterN_2 parent has " << mother->getNumChildren() << " children!" << endl; //DEBUG
				continue;
			}
			//cout << "  classifyCellDivisionTemporalWindow iterN_2 parent has " << iterN_2->treeNodePtr->parent->getNumChildren() << " children!" << endl;
			//calculate midplane feature, backpedaling to get the best midplane distance in case we missed the true division point by temporalWindowRadius
			
			left_daughter = iterN_2->treeNodePtr;
			//if ( mother->left == NULL && mother->right != NULL ) 
			//	left_daughter = mother->right;
			//else if ( mother->left != NULL ) 
			//	left_daughter = mother->left;
			//else
			//	continue; //really shouldn't get here
			d = -1; e = -1;
			for ( int iii=cellDivisionWithTemporalWindow::getTemporalWindowRadius(); iii>0; iii-- )
			{
				//d_this = lineageHyperTree::cellDivisionPlaneDistance(mother->data->centroid,left_daughter->data->centroid,auxD->data->centroid);
				//use supervoxels' PixelIdxList.size() for mother nucleus to set maximum absolute radius between the mother and attached daughter auxD and reject connection if above this limit
				//cout << "run cellDivisionPlaneDistance2021...";
				d_this = cellDivisionPlaneDistance2021(mother->data->centroid,left_daughter->data->centroid,auxD->data->centroid,scaleOrig,half_daughter_distance);
				//cout << "cellDivisionPlaneDistance2021 " << distance_return[0] << " " << distance_return[1] << endl;
				//d_this = distance_return[0] / distance_return[1];
				if ( d<0 || d_this < d ) {
					d = d_this;
					e = half_daughter_distance;//distance_return[1];
				}
				/*
				if ( e<0 || distance_return[1] < e ) {
					e = distance_return[1];
					d = distance_return[0] / distance_return[1];
				}*/
				if ( mother->parent != NULL )
				{
					left_daughter = mother;
					mother = mother->parent;						
				}
			}
			mother = iterN_2->treeNodePtr->parent; //reset back to prospective mother in prior frame
			
			//nucleus size for mother
			ss = 0; //sv_vol = 0.0f;
			for(vector< ChildrenTypeNucleus >::iterator iterS = mother->data->treeNode.getChildren().begin(); iterS != mother->data->treeNode.getChildren().end(); ++iterS )
				ss += (*iterS)->PixelIdxList.size();
			sv_vol = (float)ss;
			for(int ii = 0; ii<dimsImage; ii++)
				sv_vol *= scaleOrig[ii];
			
				
			d /= half_daughter_distance; //remember that cellDivisionPlaneDistance2021 returns total midplane distance, not ratio, so have to re-solve for ratio to compare with thrCellDivisionPlaneDistance
			if ( d<test_threshold && cbrt(sv_vol)>e ) //2021 thresholds now are 1. thrCellDivisionPlaneDistance and 2. cutoff for absolute distance between daughter cells
			{
				//maintain a top three list for prospective parents to this division child
				if ( d<min_1 )
				{
					par_3 = par_2;
					min_3 = min_2;
					par_2 = par_1;
					min_2 = min_1;
					par_1 = mother;
					min_1  = d;
				}
				else if ( d<min_2 )
				{
					par_3 = par_2;
					min_3 = min_2;
					par_2 = mother;
					min_2  = d;
				} 
				else if ( d<min_3 )
				{	
					par_3 = mother;
					min_3  = d;
				}
			}
		}
		
		if ( min_1 < test_threshold ) //got at least one hit, grouping each prospective parent/daughter candidates together for subsequent iterators, which will go in begin-end sequence
		{
			daughterNodes.push_back(auxD);
			parentNodes.push_back(par_1);
			added_any_trio = true;
			if ( par_1 == aux )
				added_original_trio = true;
			//cout << "    got one hit at TM " << frame << ", position: " << iterN->centroid[0] << "," << iterN->centroid[1] << "," << iterN->centroid[2] << endl;
			//cout << "       parent 1 position: " << par_1->data->centroid[0] << "," << par_1->data->centroid[1] << "," << par_1->data->centroid[2] << endl;
			if ( min_2 < test_threshold )
			{	
				daughterNodes.push_back(auxD);
				parentNodes.push_back(par_2);
				added_any_trio = true;
				if ( par_2 == aux )
					added_original_trio = true;
				//cout << "       parent 2 position: " << par_2->data->centroid[0] << "," << par_2->data->centroid[1] << "," << par_2->data->centroid[2] << endl;
				if ( min_3 < test_threshold )
				{	
					daughterNodes.push_back(auxD);
					parentNodes.push_back(par_3);
					added_any_trio = true;
					if ( par_3 == aux )
						added_original_trio = true;
					//cout << "       parent 3 position: " << par_3->data->centroid[0] << "," << par_3->data->centroid[1] << "," << par_3->data->centroid[2] << endl;
				}
				else
				{
					//parentNodes3.push_back(NULL);
				}
			}
			else
			{
				//parentNodes3.push_back(NULL);
				//parentNodes2.push_back(NULL);
			}
		}
		
		//test if original parent added as part of trio
		if ( aux->left != NULL )
			left_daughter = aux->left;
		else if ( aux->right != NULL )
			left_daughter = aux->right;
		else
			added_original_trio = true; //if aux has no child (shouldn't ever be the case), avoid segfault in next if statement
		
		//new test original parent, but only if not already added and meets division plane distance cutoff
		if ( !(added_original_trio) &&  lineageHyperTree::cellDivisionPlaneDistance(aux->data->centroid,left_daughter->data->centroid,auxD->data->centroid) < test_threshold )
		{
			daughterNodes.push_back(auxD);
			parentNodes.push_back(aux);
			added_any_trio = true;
		}
		
		if ( !(added_any_trio) ) //if no mother-daughter trios to consider, have to finish making own track for now-disconnected daugther
		{
			lht.lineagesList.push_back( lineage() );
			list<lineage>::iterator listLineageIter = ((++ ( lht.lineagesList.rbegin() ) ).base());//iterator for the last element in the list
			
			//set main root to "paste" sublineage starting at daughter
			listLineageIter->bt.SetMainRoot( auxD );
			
			//we need to make sure all the the elements in the binary tree point to the correct lineage
			queue< TreeNode<ChildrenTypeLineage>* > q;
			q.push( auxD );
			//TreeNode<ChildrenTypeLineage>* aux;
			while( !q.empty() )
			{
				aux = q.front();
				q.pop();
				if( aux->left != NULL) q.push(aux->left);
				if( aux->right != NULL) q.push(aux->right);
				aux->data->treeNode.setParent( listLineageIter );
			}
		}
		else
		{
			//recalculate daughter centroid here, since this will only be run once per daughter
			lht.calculateNucleiGMMParameters<float>(auxD->data,regularize_W4DOF,scaleOrig,im_plus_one);	
		}
	}

	if (daughterNodes.empty() == true)
		return 0;
	
	//recalculate parent and sibling centroids, since these are likely to be unique
	for (int dd=0; dd<parentNodes.size(); dd++ )
	{
		lht.calculateNucleiGMMParameters<float>(parentNodes[dd]->data,regularize_W4DOF,scaleOrig,im_zero);
		
		aux = NULL;
		if ( parentNodes[dd]->left != NULL )
			aux = parentNodes[dd]->left;
		else if ( parentNodes[dd]->right != NULL )
			aux = parentNodes[dd]->right;
		
		if (aux != NULL)
		{	
			lht.calculateNucleiGMMParameters<float>(aux->data,regularize_W4DOF,scaleOrig,im_plus_one);
		}
	}

	//cout << "  continue DEBUG classifyCellDivisionTemporalWindow, TM: " << frame << "====================" << endl;
	//calculate average precision matrix if necessary
	if (cellDivisionWithTemporalWindow::getUseFixPrecisionMatrix() == true)
	{
		cellDivisionWithTemporalWindow::averagePrecisionMatrixAyTM(*nucleiList_minus_1);
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
	cout << "   imgVec.size() is " << imgVec.size() << ", getTemporalWindowRadius is " << cellDivisionWithTemporalWindow::getTemporalWindowRadius() << endl;
	cout << "   daughterNodes.size() is " << daughterNodes.size() << ", parentNodes.size() is " << daughterNodes.size() << endl;
	cout << "   nucleiList->size() is " << nucleiList->size() << ", nucleiList_minus_1->size() is " << nucleiList_minus_1->size() << endl;	
	cellDivisionWithTemporalWindow::calculateBasicEllipticalHaarFeaturesBatchForCellDivisionSingleWindowForDaughters(daughterNodes, parentNodes, imgVec, cdwtVec, devCUDA, 0);


	//extend 4D Haar features using temporal combinations and classify results
	size_t numFeatures = 1000;//to pre-reserve space
	feature Fx[4];
	feature Fx_max;

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
		FxVec.reserve(cdwtVec.size()/3);
	}
	//cout << " cellDivisionWithTemporalWindow point E1." <<endl;
	int iFc = 0; //counter for trio candidates
	int iter_this_max = 0;
	//for(vector<cellDivisionWithTemporalWindow>::iterator iterF = cdwtVec.begin(); iterF != cdwtVec.end(); ++iterF ) //advance(iterF,3))
	vector<cellDivisionWithTemporalWindow>::iterator iterF = cdwtVec.begin();
	while (iterF != cdwtVec.end())
	{
		//cout << " cellDivisionWithTemporalWindow point E2." <<endl;
		auxD = iterF->cellDivisionDaughterPtr;
		//cout << " cellDivisionWithTemporalWindow point E3." <<endl;
		//feature Fx = [3];
		Fx_max = 0;
		iFc = 0; //counter for trio candidates
		iter_this_max = 0;
		
		//cout << "  cellDivisionWithTemporalWindow iterF while loop, daughter-TM:" << auxD->data->TM << "; XYZ:" << auxD->data->centroid[0] << "," << auxD->data->centroid[1] << "," << auxD->data->centroid[2] << endl;
		while ( iterF != cdwtVec.end() && auxD == iterF->cellDivisionDaughterPtr ) //take as many candidate trios as have been pushed above (limit is currently four), but when iterF lands on a different daughter, we are done with candidates
		{
			//cout << " cellDivisionWithTemporalWindow point E4, with iterF->cellDivisionPtr->getNumChildren(): " <<iterF->cellDivisionPtr->getNumChildren() <<endl;
			//iter_this_modulo = iter_count % 3;
			//vector<cellDivisionWithTemporalWindow>::iterator iterF = iterF;
			//Fx_max = 0;
			//iter_this_max = -1;
			//for ( int iFc=0; iFc<3; iFc++ ) { // try each of three division parent/children combinations
			/*DEBUG
			if ( iterF->cellDivisionPtr->data->TM + 1 != iterF->cellDivisionDaughterPtr->data->TM ) //parents and daughters have to be one timeframe apart for this algorithm; if they aren't don't even break link, just leave it
				cout << "  " << iFc << " Timeframe discrepancy: parent-" << iterF->cellDivisionPtr->data->TM << ", prospective daughter-" << iterF->cellDivisionDaughterPtr->data->TM << endl;
			if ( iterF->cellDivisionPtr->getNumChildren() != 1 )
				cout << "  " << iFc << " iterF->cellDivisionPtr->getNumChildren: " << iterF->cellDivisionPtr->getNumChildren() <<endl;
			if ( iterF->cellDivisionDaughterPtr->parent != NULL )
				cout << "  " << iFc << " iterF->cellDivisionDaughterPtr->parent not null, in TM: " << iterF->cellDivisionDaughterPtr->parent->data->TM <<endl;
			//END DEBUG*/
		
			Fx[iFc]=-1; //re-initialize
			
			if ( iterF->cellDivisionPtr == NULL || iterF->cellDivisionDaughterPtr == NULL ) //|| iterF->cellDivisionPtr->left == NULL ) 
			{	
				++iterF;
				iFc++;
				continue;
			}
			//cout << " cellDivisionWithTemporalWindow point E4b." << endl;
			if ( iterF->cellDivisionPtr->getNumChildren() > 1 || iterF->cellDivisionDaughterPtr->parent != NULL ) //test cellDivisionPtr already have two children and cellDivisionDaughterPtr already have parent
			{	
				//cout << " cellDivisionWithTemporalWindow point E4c." << endl;
				++iterF;
				iFc++;
				continue;
			}
			//cout << " cellDivisionWithTemporalWindow point E5." <<endl;
			iterF->f.reserve(numFeatures);
			//iterF->calculateFeaturesSingleWindowForDaughters();
			iterF->calculateFeaturesSingleWindowForDaughters_featureSelection_v1_2021();						
			//record for next iteration to reserve appropriate ammount of memory
			numFeatures = iterF->f.size();
	
			//call classifier to make the decision
			boostingTreeClassifierTranspose(&(iterF->f[0]),&Fx[iFc] ,classifierCDWT , 1, numFeatures);
			iterF->cellDivisionPtr->data->probCellDivision = Fx[iFc];//store result for cell division GUI (active learning) and debugging
			
			//cout << " cellDivisionWithTemporalWindow point E6." <<endl;
			
			if ( Fx[iFc] > Fx_max )
			{
				Fx_max = Fx[iFc];
				iter_this_max = iFc;
			}
			//cout << " cellDivisionWithTemporalWindow point E7." <<endl;
			++iterF;
			iFc++;
		}
		
		//okay now backpedal in iterator to the trio with Fx_max using iFc and iter_this_max
		for ( int iter_c=0; iter_c<iFc; iter_c++ )
			--iterF;
		for ( int iter_c=0; iter_c<iter_this_max; iter_c++ )
			++iterF;		
		//cout << "  About to test iter_this_max and Fx_max: " << iter_this_max << " " << Fx_max << "..." << endl;
		if( iter_this_max >= 0 && Fx_max > thrCDWT )//true cell division, so we have to create new links
		{
			//lht.cutLinkageBetweenMotherAndFurthestDaughter(iterF->cellDivisionPtr);
			/*
			//list< lineage >::iterator iterL2 = iterF->cellDivisionDaughterPtr->data->treeNode.getParent();
			list< lineage >::iterator iterL1 = iterF->cellDivisionPtr->data->treeNode.getParent();
			cout<<"Attach Lineages at Division..."<<endl;
			cout<<" Lineage iterL1..."<<endl;
			//lineageHyperTree::debugPrintLineage(iterL1);
			iterL1->debugPrintLineage();
			*/
			//we need to make sure all the the elements in the binary tree point to the correct lineage
			queue< TreeNode<ChildrenTypeLineage>* > q;
			q.push( iterF->cellDivisionDaughterPtr );
			//TreeNode<ChildrenTypeLineage>* aux;
			while( !q.empty() )
			{
				aux = q.front();
				q.pop();
				if( aux->left != NULL) q.push(aux->left);
				if( aux->right != NULL) q.push(aux->right);
				aux->data->treeNode.setParent( iterF->cellDivisionPtr->data->treeNode.getParent() );
			}
			iterF->cellDivisionDaughterPtr->parent = iterF->cellDivisionPtr;
			if ( iterF->cellDivisionPtr->right == NULL ) 
				//lht.mergeBranches( iterF->cellDivisionPtr->left, iterF->cellDivisionDaughterPtr );
				iterF->cellDivisionPtr->right = iterF->cellDivisionDaughterPtr;
			else if ( iterF->cellDivisionPtr->left == NULL ) // should not go here
				//lht.mergeBranches( iterF->cellDivisionPtr->rig/ht, iterF->cellDivisionDaughterPtr );
				iterF->cellDivisionPtr->left = iterF->cellDivisionDaughterPtr;			
			
			//iterF->cellDivisionPtr->data->treeNode.getParent()->bt.SetCurrent( iterF->cellDivisionPtr );
			//iterF->cellDivisionDaughterPtr = iterF->cellDivisionPtr->data->treeNode.getParent()->bt.insert(iterF->cellDivisionDaughterPtr->data );			
			//assert(iterF->cellDivisionDaughterPtr != NULL);			

			//DEBUG: print lineages to confirm attached
			/*
			cout<<"Attach Lineages at Division..."<<endl;
			cout<<" Lineage iterL1..."<<endl;
			//lineageHyperTree::debugPrintLineage(iterL1);
			iterL1->debugPrintLineage();
			cout<<" Lineage iterL2..."<<endl;
			iterL2->debugPrintLineage();
			//lineageHyperTree::debugPrintLineage(iterL2);
			*/

			//iterL2->bt.SetMainRootToNULL();//we have already deallocated memory, so if we do not set it to null the destructor tries ot do it again
			//lht.lineagesList.erase( iterL2 );
			
			/*cout<<" Lineage iterL1 after..."<<endl;
			iterL1->debugPrintLineage();
			cout<<" Lineage iterL2 after..."<<endl;
			iterL2->debugPrintLineage();
			*/
			cout<<"   DEBUG: TGMMsupport at frame " << frame << ": division parent " << (void *)(iterF->cellDivisionPtr) << ", new attached child " << (void *)(iterF->cellDivisionDaughterPtr) << " with " << iterF->cellDivisionPtr->getNumChildren() << " total siblings." << endl;
			
			numTrueCellDivisions++;
		}else{//false cell division -> daughter starts a new lineage
			lht.lineagesList.push_back( lineage() );
			list<lineage>::iterator listLineageIter = ((++ ( lht.lineagesList.rbegin() ) ).base());//iterator for the last element in the list
			
			//set main root to "paste" sublineage starting at daughter
			listLineageIter->bt.SetMainRoot( iterF->cellDivisionDaughterPtr );
			
			//we need to make sure all the the elements in the binary tree point to the correct lineage
			queue< TreeNode<ChildrenTypeLineage>* > q;
			q.push( iterF->cellDivisionDaughterPtr );
			//TreeNode<ChildrenTypeLineage>* aux;
			while( !q.empty() )
			{
				aux = q.front();
				q.pop();
				if( aux->left != NULL) q.push(aux->left);
				if( aux->right != NULL) q.push(aux->right);
				aux->data->treeNode.setParent( listLineageIter );
			}
			
			cout<<"   DEBUG: TGMMsupport at frame " << frame << ": rejected division with mother " << (void *)(iterF->cellDivisionPtr) << ", new birth " << (void *)(iterF->cellDivisionDaughterPtr) << " with " << iterF->cellDivisionDaughterPtr->getNumChildren() << " children." << endl;
		}

		if (writeCDWTfeatures())
		{
			iterF->writeToBinary(foutFeat);
			foutXYZ << iterF->cellDivisionPtr->data->centroid[0] << " " << iterF->cellDivisionPtr->data->centroid[1] << " " << iterF->cellDivisionPtr->data->centroid[2] << endl;
			FxVec.push_back(Fx_max);
		}
		
		//okay now forward pedal in iterator to next new prospective daughter
		for ( int iter_c=0; iter_c<iFc; iter_c++ )
			++iterF;
		for ( int iter_c=0; iter_c<iter_this_max; iter_c++ )
			--iterF;			
				
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

//2021 NEW: this function returns TOTAL midplane distance, not a ratio, and also updates the absolute distance between each daughter as a float pointer passed as parameter "sqrtnorm"
float cellDivisionTemporalWindow_TGMMsupport::cellDivisionPlaneDistance2021(float centroidMM[dimsImage], float centroidDL[dimsImage], float centroidDR[dimsImage], const float scale[dimsImage], float &sqrtnorm )
{
//calculate midplane feature
		float norm = 0.0f;
		float d = 0.0f;
		float p0, n, m;
		//static float *scale[dimsImage] = supervoxel::getScale();
		
		for(int ii = 0; ii<dimsImage; ii++)
		{
			p0 = 0.5 * (centroidDL[ii] + centroidDR[ii]) ; //midpoint
			n = (centroidDL[ii] - p0) * scale[ii];//normal			
			norm += (n * n);
			//calculate distance of mother cell to division plane
			m = (centroidMM[ii] - p0) * scale[ii];

			d += (n * m);
		}				
		//float return_arr[2];
		//return_arr[0] = fabs(d);//midplane distance * halfway distance
		//return_arr[1] = sqrt(norm); //half of distance between daughters
		//fabsd = fabs(d);
		sqrtnorm = sqrt(norm);
		return fabs(d);
		//return fabs(d) / sqrtnorm;
		//return return_arr;
		//return fabs(d) / sqrt(norm);
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
