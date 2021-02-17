#pragma once

#include <vector>
#include <gentleBoost/gentleBoost.h>
#include "../GaussianMixtureModel.h"

struct configOptionsTrackingGaussianMixture;
struct GaussianMixtureModelCUDA;
class lineageHyperTree;

namespace mylib  {
    extern "C"
    {
#include "mylib/mylib.h"
#include "mylib/array.h"
    }
};

struct CellDivisionWithGPU3DHaar {
    CellDivisionWithGPU3DHaar(string pathS);
    void apply(vector<GaussianMixtureModel*> &vecGM,
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
          GaussianMixtureModelCUDA *vecGMHOST);
private:
    vector< vector< treeStump> > classifierCellDivision;
    long long numFeatures;
};