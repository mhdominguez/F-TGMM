/*
 * parseConfigFile.h
 *
 *  Created on: Oct 19, 2011
 *      Author: amatf
 *
 *
 *      /brief Reads and parse options from a config file. Each line in the file should be configOption:value. If the line starts with # it is considered a comment
 */

#ifndef PARSECONFIGFILE_H_
#define PARSECONFIGFILE_H_

#include <string>
#include <iostream>
#include <limits> 

enum DivisionClassiferKind {
	DivisionClassiferKind_None,
	DivisionClassiferKind_BoostedElipticalHaarFeatures_AmatF2013,
	DivisionClassiferKind_LookupTable_2018
};

using namespace std;
//contains all the options to be handled by the program
struct configOptionsTrackingGaussianMixture
{

	//string imgPrefix;//image filename is composed as imgPrefix+ii+imgSufix+ii+'.tif'   where ii=frame number
	//string imgSufix;

	string imgFilePattern;//using ??? to adapt to any file notation and folder structure

	string debugPathPrefix;//folder where to output solutions

	string offsetFile;//optional parameter if time points need drift corrections	

	//string GMxmlIniFilename;//xml file containing initialization for first frame
	int maxIterEM{};
	double tolLikelihood{};//we do not consider the priors so far to decide teh stopping criteria


	//controls the strength of the priors from previous frames: it is a proportion with respect to N_k.
	double betaPercentageOfN_k{};//uncertainty in center of the Gaussian. betaPercentageOfN_k<1.0 indicates data dominates over prior to fit GM. betaPercentageOfN_k>1.0 indicates prior dominates over data to fit GM.
	double nuPercentageOfN_k{};//uncertainty in shape (covariance of the Gaussian)
	double alphaPercentage{};//prior to let Gaussian die (the lower the value, the more likely to die). Relative to alphaTotal

	//before we used super-voxels for local background subtraction and oversegmentation
	double thrSplitScore{};//check September 16th notebook to see reasonable thresholds for adaBoost classifier
	//bool useIlastikBackgroundDetector;//inidicates if we want to use background/foreground detection from Ilastik as intensity map

	double thrCellDivisionPlaneDistance{};//to remove false positive cell division detection
	
	//double thrSignal;//if useIlastikBackgroundDetector=true this value represents the cut-off probability to consider as foreground. Otherwise, it is adjusted based on thrPrcentile


	int KmaxNumNNsupervoxel{};//maximum number of nearest neighbors to consider for supervoxels
	float KmaxDistKNNsupervoxel{};//maximum distance (in pixels with scale) to consider a nearest neighbor between supervoxels

	int temporalWindowRadiusForLogicalRules{};//we will use a sliding window of +- temporalWindowRadiusForLogicalRules to improve trackign precision

	/// Aspect ratio of X-and Y-coordinates versus Z-coordinate
	/// example: anisotropyZ=5.01
	float anisotropyZ{};

	//sporious track removal
	float thrBackgroundDetectorHigh{};//parameters for background detection
	float thrBackgroundDetectorLow{};
	int SLD_lengthTMthr{};//short-lived daughter rul	

	//watershed +PBC segmentation
	int		radiusMedianFilter{};
	int useMedianFilterForTracking{};
	float	minTau{};
	float	backgroundThreshold{};
	int		conn3D{};

	//optical flow parameters
	int estimateOpticalFlow{};//0->do not use flow;1->flow is precacalculted (load file from disk);2->calculate on the fly;
	float maxDistPartitionNeigh{};//controls main regularization parameter for optical flow
	int deathThrOpticalFlow{};//if estimateOpticalFlow=0 we use this to selectively activate if there are too many deaths in the process. SET TO 0 IF YOU DON'T WANT THIS FEATURE
	/*typical choices: [2, 80.0, -20]->optical flow is calculated all the time; 
					   [0, 80.0,  20]->optical flow is only calculated if there are more than 20 deaths
					   [0, *   , -20]->optical flow is not used
					   [1, 80.0, -20] ->pre-calculate image flow is read and used all the time
	*/

	int tau{}; //parameter for persistance segmentation

	//precision matrix regularization boundaries
	double lambdaMin{};//aux=scaleSigma/(maxRadius*maxRadius) with scaleSigma=2.0 and maxRadius=10 (adjust with scale)
	double lambdaMax{};//aux=scaleSigma/(maxRadius*maxRadius) with scaleSigma=2.0 and minRadius=3.0 (adjust with scale)  (when nuclei divide they can be very narrow)
	double maxExcentricity{};//maximum excentricity allowed: sigma[i]=1/sqrt(d[i]). Therefore maxExcentricity needs to be squared to used in terms of radius.


	//-------------- Trimming supervoxels using Otzu's threshold to differentiate background from foreground--------------
	//Minimum size (in voxels) of a supervoxel(smaller supervoxel will be deleted)
	int minNucleiSize{};
	//Maximum size (in voxels) of a supervoxel (considered when we apply Otzu's threshold)
	int	maxNucleiSize{};
	//Maximum percentage of voxels in a supervoxel belonging to foreground (considered when we apply Otzu's threshold)
	float	maxPercentileTrimSV{};
	//Connectivity considered when trimming supervoxels using Otzu's threshold
	int	conn3DsvTrim{};
	//Remove supervoxels that are "garbage" by counting the ratio of holes / number of supervoxels(supervoxels should be filled objects).Set the ratio to a negative number to deactivate it
	float ratioHolesSv{};
	//Set to zero to not keep original supervoxels after trimming(default bhavior).Set to something different than zero to keep original supervoxels and write them out in the svb file
	bool svKeepWithoutTrimming{};

	//
	// Cell division classifer
	//

	struct CellDivisionClassifier{
		enum DivisionClassiferKind method;

		struct Amatf2013 {
			float thrCDWT;//threshold for the classifier to discriminate true cell division from false ones based on 3D Haar Elliptical features along a temporal window. If the classifier score is above thrCDWT, then we consider it a true division
			bool writeCDWTfeaturesToTrainClassifier;
			string classifierFileCDTW;
		} Amatf2013;

		struct LUT {
			string filename;
			float  spatial_threshold_px;        ///< spatial radius in isotropic px. "isotropic" means in units of the x-pixels.
			unsigned temporal_threshold_frames; ///< temporal radius in frames
			float anisotropyZ;                  ///< mirrors anisotropyZ from above
		} LUT;

		int validate() const; ///< returns 0 when valid, 1 otherwise
	} cellDivisionClassifier;

	void setDefaultValues()
	{
		//these values are set to "non-sense" values since they are mandatory in the config file
		//imgPrefix=string("empty");
		//imgSufix=string("empty");
		//GMxmlIniFilename=string("empty");
		imgFilePattern = string("empty");
		debugPathPrefix = string("empty");

		offsetFile = string("empty");

		//most of these values should not be modified but user can decide to do it in the config file
		maxIterEM=400;
		tolLikelihood=1e-6;

		betaPercentageOfN_k=0.01;
		nuPercentageOfN_k=1.0;
		alphaPercentage = 0.8;

		thrSplitScore=-7.0;

		thrCellDivisionPlaneDistance = 3.2;
		
		//useIlastikBackgroundDetector=true; //from previous approach before super-voxels
		//thrSignal=0.2;					//from previous approach before super-voxels

		estimateOpticalFlow=0;//optical flow deactivated by default
		maxDistPartitionNeigh=80;
		deathThrOpticalFlow=-20;

		tau=0;

		KmaxNumNNsupervoxel = 10;
		KmaxDistKNNsupervoxel = 50.0f;
		temporalWindowRadiusForLogicalRules = 1;

		anisotropyZ = 1.0f;

		thrBackgroundDetectorHigh = 1.2f;//above 1.0 background detection is not activated
		thrBackgroundDetectorLow = 0.2f;
		SLD_lengthTMthr = 5;

		// cell division classifer

		cellDivisionClassifier.method=DivisionClassiferKind_None;

		cellDivisionClassifier.Amatf2013.thrCDWT = 0.3f;
		cellDivisionClassifier.Amatf2013.writeCDWTfeaturesToTrainClassifier = false;
		cellDivisionClassifier.Amatf2013.classifierFileCDTW="classifierCellDivisionWithTemporalWindow.txt";

		cellDivisionClassifier.LUT.filename="cellDivisionLookupTable.csv";
		cellDivisionClassifier.LUT.spatial_threshold_px = 10.0f;
		cellDivisionClassifier.LUT.temporal_threshold_frames = 5;

		//elements for watershed segmentation
		radiusMedianFilter = 2;
		useMedianFilterForTracking = 0;
		minTau = 2;
		backgroundThreshold = -1e32f;
		conn3D = 74;

		//precision matrix regularization boundaries for GMM fitting
		lambdaMin = 0.02;//so it can be checked if they ar einitialized or not
		lambdaMax = 0.2222;
		maxExcentricity = 3.0 * 3.0;

		//-------------- Trimming supervoxels using Otzu's threshold to differentiate backgroudn from foreground--------------
		//Minimum size (in voxels) of a supervoxel(smaller supervoxel will be deleted)
		minNucleiSize = 50;
		//Maximum size (in voxels) of a supervoxel (considered when we apply Otzu's threshold)
		maxNucleiSize = 3000;
		//Maximum percentage of voxels in a supervoxel belonging to foreground (considered when we apply Otzu's threshold)
		maxPercentileTrimSV = 0.4f;
		//Connectivity considered when trimming supervoxels using Otzu's threshold
		conn3DsvTrim = 6;

		ratioHolesSv = std::numeric_limits<float>::max();;//by default it is deactivated
		//Set to zero to not keep original supervoxels after trimming(default bhavior).Set to something different than zero to keep original supervoxels and write them out in the svb file
		svKeepWithoutTrimming = false;
	}

	
	/*
	 *\brief returns 0 if everything was fine.returns>0 if there was an error
	 */
	int parseConfigFileTrackingGaussianMixture(const string &filename);
	void printConfigFileTrackingGaussianMixture(ostream &outLog) const;

	configOptionsTrackingGaussianMixture(){setDefaultValues();};//default constructor
};


#endif /* PARSECONFIGFILE_H_ */
