# Key parameters

The threshold for persistence-based agglomeration of watershed regions
(persistenceSegmentationTau) and the intensity threshold for defining the
background level in each recording (backgroundThreshold) both refer to image
properties and are straight-forward to determine by visual inspection of the
image volume at a late time point of the time-lapse recording. In general,
inspecting late time points is more useful, since (depending on the
experiment) intensity levels are often slightly dimmer and cell densities
higher. Measurements in this scenario provide a lower bound constraint for
both values.

### Intensity threshold

To determine, inspect a region in the image volume outside the specimen (for
example, by using the open-source software ImageJ) and determine the mean
intensity level in this background region. It is preferable to be
conservative, i.e. to set a lower value so as not to miss cell nuclei, since
false negatives can alter the coherence between time points. Moreover, we
compute a local threshold for each super-voxel using Otsu’s and, thus, even if
background regions are included in the foreground estimate, this will not
affect the final shape of the super-voxels. The only drawback of a lower
background threshold is a small increase in computation time.

### Persistence-based agglomeration threshold

To determine the threshold for persistence-based agglomeration of watershed
regions (τ), plot the intensity profile across the line connecting two of the
dimmest nuclei centroids in the image stack (for example, by using the
open-source software ImageJ). The profile should have two peaks (nuclei
centroids) and a valley (nuclei borders).  The threshold τ should be set to a
value smaller than the difference between the intensity values of the peaks
and the valley, such that the corresponding nuclei are not merged into a
single super-voxel (under-segmentation). In our experience, a value of τ
between 5 and 20 tends to be sufficient to compensate for the watershed
over-segmentation of noisy regions, without risking merging of dim cell
nuclei.

Note that, although care should be taken to set these parameters
appropriately, one should be able to obtain good results for a fairly wide range
of parameter values. For images with lower signal-to-noise ratio (SNR), such
as confocal microscopy images, the value of τ is more critical than the
background intensity because watershed regions fragment the image into smaller
regions. In contrast, in images with high SNR, such as most light-sheet
microscopy images, the intensity background is more relevant because the
watershed algorithm already produces super-voxels that follow nucleus
morphologies fairly well.

# Advanced parameters

In this section, we provide an overview of all advanced framework parameters.
Note that these parameters were not changed across the computational
reconstructions and data sets presented in this study. To complement the
descriptions below, we also provide the default parameter values in the
configuration file `$ROOT_TGMM\data\TGMM_configFile.txt`, which is included in
`Supplementary_Software_1.zip`.

### `betaPercentageOfN_k`

Non-negative floating point scalar. betaPercentageOfN_k defines the prior
probability for the centroid position of a nucleus based on its position in
the previous time point. No motion model is used, unless the optical flow
module is activated. Thus, if cells are moving fast this parameter should be
set close to zero. If cells are moving very little, and the position at time `t`
is a good prediction for the position at time `t + 1`, this parameter should be
set to 1 or greater.

### `nuPercentageOfN_k`

Non-negative floating point scalar. Follows the same concept as
betaPercentageOfN_k, but for shape variation between consecutive time points.
If objects change shape rapidly between two consecutive time points this
parameter should be set close to zero. If object shapes change very little
between time points this parameter should be set to 1 or greater.

### `alphaPercentage`

Floating point scalar. alphaPercentage controls the prior probability of death
for a track. The more likely nuclei are to disappear from the image or undergo
apoptosis, the lower the value of this parameter should be.

### `maxIterEM`

Integer positive number. maxIterEM defines the maximum number of iterations of
variational inference allowed each time a Gaussian mixture model is fitted. In
general, very few rounds are needed (less than 10), since the model is
initialized with the solution from the previous time point. Thus, this
parameter is implemented as a precaution. The terminal output of TGMM.exe can
be used to obtain an estimate of the typical number of iterations needed and
maxIterEM can then be set accordingly.

### `tolLikelihood`

Floating point positive number. tolLikelihood is used as a stopping criterion
for variational inference of the Gaussian mixture model. Optimization is
stopped if the relative increase in likelihood between two consecutive
iterations is less than the value of this parameter. Thus, the lower the
value, the more iterations of variational inference are run.

### `regularizePrecisionMatrixConstants_lambdaMax`
Floating point positive number. This parameter provides the maximum allowed
value for any of the eigenvalues (in pixels) of the covariance matrix defining
each object in the Gaussian mixture model. Thus, the larger the value of
regularizePrecisionMatrixConstants_lambdaMax, the larger the ellipsoids
fitting each nucleus can grow. This parameter is used during regularization of
the variational inference results.

### `regularizePrecisionMatrixConstants_lambdaMin`
Floating point positive number. This parameter provides the minimum allowed
value for any of the eigenvalues (in pixels) of the covariance matrix defining
each object in the Gaussian mixture model. Thus, the lower the value of
regularizePrecisionMatrixConstants_lambdaMin, the smaller the ellipsoids
fitting each nucleus can shrink. This parameter is used during regularization
of the variational inference results.

### `regularizePrecisionMatrixConstants_maxExcentricity`
Floating point positive number. This parameter provides the maximum
eccentricity between any two principal axes of the ellipsoid defining a
nucleus. This parameter is used during regularization of the variational
inference results.

### `temporalWindowForLogicalRules`
Positive integer number. This parameter provides the radius of the temporal
window (total window length is [2 × temporalWindowForLogicalRules + 1] time
points) used to apply spatio-temporal heuristic rules for fixing tracking
errors. The larger the value, the more memory is required, since the program
will keep the image data of all concerned time points in memory to be able to
calculate the features needed to apply the various heuristics.

### `thrBackgroundDetectorHigh` and `thrBackgroundDetectorLow`
Floating point non-negative numbers. These parameters provide the thresholds
applied to the results of the background track detector for removing
trajectories representing non-nuclei objects. They control the behavior of a
hysteresis filter applied over time to the background probability scores of
multiple data points belonging to the same lineage. When the program detects a
data point with a background probability above thrBackgroundDetectorHigh it
proceeds with deleting its descendants until the probability falls below
thrBackgroundDetectorLow. Thus, the higher the value of
thrBackgroundDetectorHigh, the fewer objects are removed. If
thrBackgroundDetectorHigh is above 1, then no background track removal is
applied.

### `SLD_lengthTMthr`
Non-negative integer number. When SLD_lengthTMthr>0, any daughter branch that
ends within less than SLD_lengthTMthr time points after division is considered
a spurious over-segmentation event and is deleted.  Otherwise, setting this to
0 deactivates the `short lived daughter` detection-and-merge rules.

###`radiusMedianFilter`
Positive integer number. This parameter provides the radius (in pixels) of the
median filter applied before the watershed hierarchical segmentation is
performed. The noisier the images, the larger this value should be.

### `minTau`
Non-negative floating point value. This parameter provides the minimum value
of τ used for the hierarchical segmentation using persistence-based clustering
of watershed regions. The higher minTau, the larger the minimum super-voxel
size that can be generated at the lower level of the hierarchical
segmentation. This value should be kept low so as not to compromise the
framework’s capability to recover from under-segmentation.

### `conn3D`
Values allowed are 6, 28 and 74. This parameter defines the 3D local
neighborhood used to run watershed for generating super-voxels. Values of 6
and 28 define traditional 3D neighborhoods, whereas a value of 74 
generates cubes of 5 x 5 x 3 around each point to address the anisotropy of
the point-spread-function typically encountered in 3D microscopy images.

### `estimateOpticalFlow`
Values allowed are 0, 1 and 2. This parameter activates/deactivates the use of
optical flow calculations between time points for the purpose of compensating
for large object displacements. A value of 0 deactivates optical flow
calculations. A value of 1 indicates that pre-calculated optical flow files
are available and can be used to apply local motion displacements between time
points. A value of 2 indicates that the program will calculate optical flow
on-the-fly, using a routine provided by the user.

### `maxDistPartitionNeigh`
Floating point positive number. It is only used if `estimateOpticalFlow` is
equal to 2 and the routine called is the one described in F. Amat et al.,
`Fast and robust optical flow for time-lapse microscopy using super-voxels`
(Bioinformatics, 2013). The parameter provides the maximum allowed distance
(in pixels) between super-voxels for them to be considered neighbors in the
calculation of the optical flow (coherence constraint).

### `deathThrOpticalFlow`
Integer number. If positive, the optical flow module will be activated
automatically when the number of deaths at a specific time point is larger
than the value of this parameter. Usually, when large motions occur (larger
than one nucleus diameter from one time point to the next), many Gaussians in
the model disappear, since the solution from the previous time point is not
well-suited for initialization of the current time point. Thus, monitoring
deaths can be used as a trigger to activate optical flow only when needed.

### `minNucleiSize`
Positive integer number. If the number of voxels belonging to a super-voxel is
less than minNucleiSize the super-voxel is deleted. This parameter is useful
to delete spurious super-voxels representing background intensity.

### `maxNucleiSize`
Positive integer number. This parameter defines the maximum allowed size (in
voxels) of a super-voxel after applying Otsu's threshold. If Otsu's threshold
generates an object larger than maxNucleiSize the threshold is increased until
the objet size falls below maxNucleiSize.

### `maxPercentileTrimSV`
Floating point number between 0 and 1. This parameter defines the maximum
allowed percentage of voxels in a super-voxel belonging to foreground. If
Otsu's threshold generates an object larger than maxPercentileTrimSV the
threshold is increased until the percentage of foreground voxels falls below
below maxPercentileTrimSV.

### `conn3DsvTrim`
Values allowed are 6, 28 and 74. The final super-voxel generated after
trimming the initial super-voxel partition to detect foreground and background
is guaranteed to have this connectivity.

### `maxNumKNNsupervoxel`
Positive integer number. This parameter defines the maximum number of nearest
neighbors to consider for each super-voxel when building the spatio-temporal
graph for tracking. The shorter the nuclear displacement between time points,
the lower the parameter value can be.

### `maxDistKNNsupervoxel`
Floating point positive number. This parameter defines the maximum distance
(in pixels) to consider for each super-voxel when building the spatio-temporal
graph for tracking. The shorter the nuclear displacement between time points,
the lower the parameter value can be.

### `thrSplitScore`
Floating point number. If 3D Haar features are used for cell division
classification, this parameter sets the threshold for the machine learning
classifier to decide whether a cell is dividing or not. The higher the
threshold, the fewer divisions are going to be called by the classifier. In
order to activate 3D Haar features, the code needs to be recompiled with
preprocessor directive CELL_DIVISION_WITH_GPU_3DHAAR.

### `thrCellDivisionPlaneDistance`
Floating point positive number, defining the threshold of a feature for
disregarding cell division false positives. The feature calculates the
distance (in pixels) between mother cell and the midplane defined by the two
daughter cells. If the value is above thrCellDivisionPlaneDistance, the cell
division is considered a false positive and the linkage between mother and
furthest daughter is removed. The default value of 3.2 was determined
empirically from a small training set to maximize precision while maintaining
a high recall of true cell divisions. The lower the value, the lower the
recall of cell divisions and the higher the precision.

### `thrCellDivisionWithTemporalWindow`
Floating point number between 0 and 1. In order to improve cell division
detection accuracy the 3D Haar features have been combined across the temporal
window. This parameter sets the threshold for the machine learning classifier
to decide whether a cell is dividing or not based on these new set of
features. The higher the threshold, the fewer divisions are going to be called
by the classifier.

