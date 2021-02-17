#ifndef CELLDIVISIONLOOKUPTABLE_H
#define CELLDIVISIONLOOKUPTABLE_H
#include "DivisionClassifierFactory.h"
#include "../Utils/parseConfigFile.h"
#include <map>

class CellDivisionLookupTable : public IDivisionClassifier {
public:
	using config=struct configOptionsTrackingGaussianMixture::CellDivisionClassifier::LUT;
	struct float3 { float x,y,z; };

	CellDivisionLookupTable(config& config);

	virtual ~CellDivisionLookupTable() = default;

	
	int classifyCellDivisionTemporalWindow(lineageHyperTree& lht,int frame,std::vector<mylib::Array*>& imgVec,int devCUDA) final;
	

	/**
	* \brief  Used for debugging/progress messages
	* \return Number of divisions that were not rejected from the last classifyCellDivisionTemporalWindow call.
	*/
	size_t getNumCellDivisions();

	/**
	* \brief Used for debugging/progress messages
	* \return value of (spatial) threshold used for rejecting putative divisions
	*/
	float getThrCDWT();

	std::ostream& print(std::ostream& os);

private:
	int load(std::string filename);
	
	float threshold_r2_;      // spatial threshold in radius-squared (units: isotropic pixels^2)
	int temporal_window_radius_;  // spatial threshold in radius-squared (units: frames)
	size_t ndivisions_;
	float anisotropy_z_; ///< Aspect ratio of X-and Y-coordinates versus Z-coordinate. e.g anisotropy_z_=5 implies the z sampling interval is 5 times larger than the x sampling interval.

	// division table
	// These are divisions predicted from some (external) classifier
	// and supplied via lookup table.	
	std::map<int,std::vector<float3>> divisions_; // map TM => set of (x,y,z) locus of detected divisions at time TM
};

#endif // CELLDIVISIONLOOKUPTABLE_H
