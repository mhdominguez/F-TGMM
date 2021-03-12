#ifndef DIVISIONCLASSIFIERFACTORY_H
#define DIVISIONCLASSIFIERFACTORY_H

#ifdef __cplusplus

#include <vector>
#include <memory>
#include "../Utils/parseConfigFile.h"
#include "lineageHyperTree.h"
 namespace mylib
 {
 	#include <mylib/array.h>
 }

class DivisionClassifierFactory;
class lineageHyperTree;

/**
 * \brief Cell division classifier interface used by TGMM.
 * 
 * Abstracts over a couple different division classifier methods.
 * Historically, there was one such method and at some point others
 * were added.  This history has shaped the particular way the
 * interface is defined.
 */
class IDivisionClassifier {
protected:
	~IDivisionClassifier() = default;
public:

	/**
	 * \brief main function to call from TGMM
	 *        this function will disconnect edges of the elements that are below thrCDWT
	 * \param lht     Representation of cell lineages through time
	 * \param frame   time-point at which classification should occur.
	 *                Usually this is the frame corresponding to the beginning of imgVec
	 * \param imgVec  vector of volumes over temporal window
	 * \param devCUDA cuda device to use for any gpu-related computation
	 * \return -1 if there is an error. Otherwise returns the number of true cell divisions
	 */
	virtual int    classifyCellDivisionTemporalWindow(lineageHyperTree& lht, int frame, std::vector<mylib::Array*>& imgVec, int devCUDA, double thrCellDivisionPlaneDistance, float* im_zero, float* im_plus_one, bool regularize_W4DOF, float scaleOrig[3])=0;
	//virtual int    classifyCellDivisionTemporalWindow(lineageHyperTree& lht, int frame, std::vector<mylib::Array*>& imgVec, int devCUDA, double thrCellDivisionPlaneDistance)=0;
	virtual int    classifyCellDivisionTemporalWindow(lineageHyperTree& lht, int frame, std::vector<mylib::Array*>& imgVec, int devCUDA)=0;

	/**
	 * \brief  Used for debugging/progress messages
	 * \return Number of divisions that were not rejected from the last classifyCellDivisionTemporalWindow call.
	 */
	virtual size_t getNumCellDivisions()=0;

	/**
	* \brief Used for debugging/progress messages
	* \return value of threshold used for rejecting putative divisions
	*/
	virtual float  getThrCDWT()=0;

	/**
	 * \brief Print a description of this object to an output stream
	 *        for debugging.
	 * \return pass-thru of the output stream
	 */
	virtual std::ostream& print(std::ostream& os)=0;
};

/**
 * \brief Detects the kind of classifier referenced by the `classifier_file_path`
 *        then configures (@see configure) and instances(@see create) the classifier.
 *        Owns the IDivisionClassifier instance.
 *        
 * \aside There really isn't ever more than one classifier instanced into TGMM at a time.
 *	      If it weren't for the configure() function, it should probably just get made into
 *        a single function that returns the instance given the filename.
 */
class DivisionClassifierFactory {
	struct State;
	State* state;

public:
	
	/**
	 * \brief Detects the kind of classifier referenced by the `classifier_file_path`
	 *        then configures (@see configure) and instances(@see create) the classifier.
	 * \param classifier_file_root_path Path to the file containing classifier data.
	 */
	explicit DivisionClassifierFactory(std::string root_path,configOptionsTrackingGaussianMixture::CellDivisionClassifier& config,int temporalWindowRadius);
	~DivisionClassifierFactory();
    
	const std::shared_ptr<IDivisionClassifier>& get() const;

    friend std::ostream& operator<<(std::ostream& os , const DivisionClassifierFactory& f);
};
#endif // __cplusplus
#endif // DIVISIONCLASSIFIERFACTORY_H


