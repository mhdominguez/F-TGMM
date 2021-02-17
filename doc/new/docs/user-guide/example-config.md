# Example TGMM Config File

```
#Configuration file to set parameters for Tracking with Gaussian Mixture Models (TGMM) software. 
#This file is a template for use with the test data set provided with the software package.
#Use the pound symbol (#) for comments.
#Visit the manual or paper for more information on each parameter.


#================== PATH PARAMETERS ==================

#File location. "?" symbols will be substituted by time point index.
imgFilePattern=C:/Users/Fernando/cppProjects/TrackingGaussianMixtures/NM2013-paperRun/data/data/TM?????_timeFused_blending/SPC0_CM0_CM1_CHN00_CHN01.fusedStack_?????

#Folder where results will be saved. This folder needs to exist and the software will generate a subfolder for each run.
debugPathPrefix=E:/TGMMruns


#================== MAIN PARAMETERS ==================

#Aspect ratio of X- and Y-coordinates versus Z-coordinate
anisotropyZ=5.01

#Image background level (set conservatively)
backgroundThreshold=200

#Tau level for generating an initial set of supervoxels from the hierarchical segmentation of watershed and persistence-based clustering
persistanceSegmentationTau=14


#================ ADVANCED PARAMETERS ================

#NOTE: These parameters were not changed between the various runs presented in the associated paper.
#Most of these parameters should not be altered. However, you can modify them if required by your data.

#--- Variational inference for Bayesian GMM

betaPercentageOfN_k=0.01
nuPercentageOfN_k=1.0
alphaPercentage=0.8

#Stopping criteria
maxIterEM=100
tolLikelihood=1e-6

#Boundaries for Gaussian shapes
#aux=scaleSigma/(maxRadius*maxRadius) with scaleSigma=2.0 and maxRadius=10 (adjust with scale)
regularizePrecisionMatrixConstants_lambdaMin=0.02
#aux=scaleSigma/(maxRadius*maxRadius) with scaleSigma=2.0 and minRadius=3.0 (adjust with scale) (the shape of the fluorescently-labeled chromatin can become highly anisotropic when nuclei divide)
regularizePrecisionMatrixConstants_lambdaMax=0.2222
#Maximum excentricity allowed: sigma[i]=1/sqrt(d[i]). maxExcentricity needs to be squared to relate to the radius.
regularizePrecisionMatrixConstants_maxExcentricity=9.0

#--- Removal of tracks belonging to background objects instead of cell nuclei

#Radius of temporal window use to apply combinatorial rules
temporalWindowForLogicalRules=5
#Hysteresis values to apply the background classifier
#If thrBackgroundDetectorHigh >= 1.0, then the background classifier is not used
thrBackgroundDetectorHigh=0.7
thrBackgroundDetectorLow=0.2
#Short-lived daughter parameters
SLD_lengthTMthr=5

#--- Hierarchical segmentation with watershed and persistence-based clustering

radiusMedianFilter=2
minTau=2
#Connectivity for three-dimensional elements. Allowed values are 6, 26, 74.
conn3D=74

#--- Optical flow options

#Decide whether we need to use optical flow to compensate for large movements. 0 -> no optical flow; 1 -> load precalculated optical flow; 2 -> calculate optical flow on the fly
estimateOpticalFlow=0
#Distance (in pixels) to consider nearest neighbor cells for coherent optical flow
maxDistPartitionNeigh=80.0
#If positive, optical flow module will be activated automatically if the number of deaths exceeds this value. Usually, when large movements occur, many Gaussians dissapear, i.e. this circumstance can be used to activate optical flow.
deathThrOpticalFlow=-1

#--- Trimming supervoxels using Otsu's threshold to distinguish background from foreground

#Minimum size (in voxels) of a super-voxel (smaller super-voxels will be deleted)
minNucleiSize=50
#Maximum size (in voxels) of a super-voxel (considered when we apply Otsu's threshold)
maxNucleiSize=3000
#Maximum percentage of voxels in a super-voxel belonging to foreground (considered when we apply Otsu's threshold)
maxPercentileTrimSV=0.4
#Connectivity considered when trimming super-voxels using Otsu's threshold.
conn3DsvTrim=6

#--- Other parameters

#Maximum number of nearest neighbors to consider when building the spatio-temporal graph for tracking
maxNumKNNsupervoxel=10
#Maximum distance (in pixels) of nearest neighbors to consider when building the spatio-temporal graph for tracking
maxDistKNNsupervoxel=41.0
#If using 3D Haar-like features for cell division classification, this is the threshold for the machine learning classifier. You need to recompile the code in order to activate 3D Haar-like features.
thrSplitScore=-1.0
#threshold used to remove false positive cell divisions based on the distance of the mother cell to the midplane (cell division plane) between the two daughters
thrCellDivisionPlaneDistance=10000.5

#--- Division classifier parameters

#Supported method: None, AmatF2013, LUT2018
cellDivisionClassifierMethod=LUT2018

#threshold used to remove false positive cell divisions based on a classifier with temporal 3D Haar ellliptical features
thrCellDivisionWithTemporalWindow=-0.1

cellDivisionLUT_filename=C:\Users\vidrio\Downloads\Data_S5\Data_S5.csv
cellDivisionLUT_spatial_threshold_px=30
cellDivisionLUT_temporal_threshold_frames=2
```