/*
 * Copyright (C) 2011-2012 by  Fernando Amat
 * See license.txt for full license and copyright notice.
 *
 * Authors: Fernando Amat 
 *
 * temporalLogicalRules.h
 *
 *  Created on: August 17th, 2012
 *      Author: Fernando Amat
 *
 * \brief Implements methods to inforce temporal logical rules in a cell lineage
 *
 */


#ifndef __TEMPORAL_LOGICAL_RULES_H__
#define __TEMPORAL_LOGICAL_RULES_H__

#include <iostream>
#include "GaussianMixtureModel_Redux.h"
#include "lineageHyperTree.h"

using namespace std;


int mainTestTemporalLogicalRules( int argc, const char** argv );

//functions to parse data from GMM pipeline
int parseGMMtrackingFilesToHyperTreeLoadVols(string imgPrefix, string imgSuffix, string basenameXML,int iniFrame,int endFrame,int tau, lineageHyperTree &lht, TimeSeriesMapT &time_series_map);
//more up-to-date version (not only for debuggging) which also reads supervoxels from binary .svb files
int parseGMMtrackingFilesToHyperTree(const string imgFilePattern, const string basenameXML,const int iniFrame,const int endFrame, lineageHyperTree &lht, bool clearLHT = true);
void getImgPath(string imgPrefix, string imgSuffix, int tau, string& imgPath, string& imgLpath, int intPrecision);
void getImgPath2(string imgPrefix, string imgSuffix, string imgSuffix2, int tau, string& imgPath, string& imgLpath, int intPrecision);

//debugging functions
void testListIteratorProperties();

#endif
