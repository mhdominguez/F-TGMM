# TGMM Output Format

The folder `debugPathPrefix\GMEMtracking3D_%date` contains the output of the
TGMM run.

The final result can be found in one of:

* `%debugPathPrefix\GMEMtracking3D_%date\XML_finalResult_lht` or
* `%debugPathPrefix\GMEMtracking3D_%date\XML_finalResult_lht_bckgRm`.

The latter directory is used if the user applied the background classifier.
The output sub-folder contains one XML file and one `.svb` file per time
point. 

## XML

An XML file contains the main tracking and segmentation information. Each
object is stored under the tag `<GaussianMixtureModel>` with the following
attributes:

*	`id [integer]`: unique id of the object in this particular time point.
*	`lineage [integer]`: unique id of the cell lineage the object belongs to.
*   `parent [integer]`: id of the linked object at the previous time point.
    Following the chain of `parent` objects reconstructs the track. A value of
    -1 indicates the birth of a track.
*   `splitScore [float]`: confidence level for the correct tracking of this
    particular object. A value of 0 indicates very low confidence and a value
    of 5 indicates very high confidence. Sorting elements by confidence level
    can guide the user in the data curation process and facilitate more
    effective editing of the TGMM results (see main text and Fig. 4).
*	`scale [float[3]]`: voxel scaling factors along the x-, y- and z-axis.
*	`nu, beta, alpha [float]`: value of the hyper-parameters for the Bayesian GMM.
*	`m [float[3]]`: mean of the Gaussian Mixture (object centroid, in pixels).
*	`W [float[3][3]]`: precision matrix of the Gaussian Mixture (object shape).
*   `Prior`: same as before, but for prior values obtained from the previous
    time point. These values are used during the inference procedure.
*   `svIdx [integer[]]`: list of indices of the super-voxels clustered by this
    Gaussian. Together with the `.svb` file, this information can be used to
    obtain precise segmentation regions for each object.

## SVB

The `.svb` file is a binary file in a custom format that can be read with the
constructor:
```c++
supervoxel::supervox (istream& is)
```

It contains information about all super-voxels generated at a particular time
point. The precise segmentation mask for each object can be recovered using
the `svIdx` attribute. 

A set of Matlab scripts in the zip file `readTGMM_XMLoutput.zip` hosted on
[sourceforge](http://sourceforge.net/projects/tgmm/files/) in order to import
the data to Matlab for further analysis. The zip file also includes a README
file explaining how to use these scripts.
