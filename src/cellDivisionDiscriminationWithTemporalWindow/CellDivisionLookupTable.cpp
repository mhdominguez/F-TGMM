#include "CellDivisionLookupTable.h"
#include "binaryTree.h"
#include "lineage.h"
#include "lineageHyperTree.h"
#include <sstream>

using namespace std;

using float3 = CellDivisionLookupTable::float3;

static float3 operator-(const float3& a, const float3& b) {
	return float3{a.x-b.x,a.y-b.y,a.z-b.z};
}

static float3 operator*(const float3& a, const float3& b) {
	return float3{a.x*b.x,a.y*b.y,a.z*b.z};
}

static float sum(const float3& a) {
	return a.x+a.y+a.z;
}

CellDivisionLookupTable::CellDivisionLookupTable(config& config) 
	: threshold_r2_(config.spatial_threshold_px*config.spatial_threshold_px)
    , temporal_window_radius_(config.temporal_threshold_frames) 
	, ndivisions_(0)
	, anisotropy_z_(config.anisotropyZ)
{
	if(load(config.filename)) {
		// got an error code
		throw runtime_error("Could not load division classifier data from "+config.filename);
	}
}

std::ostream& CellDivisionLookupTable::print(std::ostream& os) {
	return os<<"DivisionClassifer(LUT2018: dr2="<<threshold_r2_<<"px^2, dt="<<temporal_window_radius_<<"frames)";
}

int CellDivisionLookupTable::classifyCellDivisionTemporalWindow(
	lineageHyperTree& lht,
	int frame,
	std::vector<mylib::Array*>& imgVec, int devCUDA) 
{	
	unsigned num_true_cell_divisions=0;
	//find all the elements that contain a division
	const auto nuclei_list = &(lht.nucleiList[frame]);
	vector<TreeNode< ChildrenTypeLineage >* > divisionNodes;//contains a pointer to each of the divisions in a specific time point
	divisionNodes.reserve(nuclei_list->size() / 10);

	for (const auto & iterN : *nuclei_list)
	{
		if (iterN.treeNodePtr->getNumChildren() == 2)
			divisionNodes.push_back(iterN.treeNodePtr);
	}

	if (divisionNodes.empty())
		return 0;

	for(const auto division: divisionNodes) {	
		auto r0=*reinterpret_cast<float3*>(&division->data->centroid);
		
		// Handling anisotropy
		//	 Assume: centroid is in isotropic coordinates (can this be confirmed?)
		//           LUT is anisotropic coordinates
		//   Known:  The spatial threshold is defined for isotropic coordinates.
		//
		//   Approach: Need to compute in isotropic space so
		//			   Convert LUT coords to isotropic by r.z*=anisotropy
		//             But, don't want to do that N times.
		//             More efficient to scale r0.z, but need to compensate in dr.
		//					dr=a*z-z0=a(z-z0/a)

		// 2018 Feb - looks like the centroid coordinates are _not_ scaled 
		//            for anisotropy.
		//					dr==a*(z-z0)

		auto any=false;
		const auto
			tmin=max(0,frame-temporal_window_radius_),
			tmax=frame+temporal_window_radius_;

		for(auto t=tmin;t<=tmax;++t) {
			const auto it=divisions_.find(t);
			if(it==divisions_.end())
				continue;
			// it->second is the vector of divisions in the look up table at time `frame`
			for(const auto & r:it->second) {
				auto dr=r-r0;
				dr.z*=anisotropy_z_;
				if(sum(dr*dr)<threshold_r2_) {
					any=true;
					goto EndOfSearch;
				}					
			}			
		}
	EndOfSearch:
		if(any) {
			++num_true_cell_divisions;
		} else {  // reject division						
			lht.cutLinkageBetweenMotherAndFurthestDaughter(division);
		}
	}
	ndivisions_=num_true_cell_divisions;
	return num_true_cell_divisions;
}

size_t CellDivisionLookupTable::getNumCellDivisions() {
	return ndivisions_;	
}

float CellDivisionLookupTable::getThrCDWT() {
	return threshold_r2_;
}

/* See below for excerpt of Example CSV */

static pair<int,float3> parse(const string& line) {
	const char sep=',';
	string tok;
	auto in=istringstream(line);
	int t; 
	float3 r;
	getline(in,tok,sep); istringstream(tok)>>r.x;
	getline(in,tok,sep); istringstream(tok)>>r.y;
	getline(in,tok,sep); istringstream(tok)>>r.z;
	getline(in,tok,sep); istringstream(tok)>>t;
	return make_pair(t,r);
}

int CellDivisionLookupTable::load(std::string filename) {
	auto fin=ifstream(filename);
	if(!fin.good())
		throw runtime_error("Could not open file for reading: "+filename);
	auto rowid=-1; // offset for number of header rows
	for (string line; getline(fin, line); ++rowid) {
		if(rowid>=0) {
			const auto v=parse(line);
			divisions_[v.first].push_back(v.second);
		}
	}
	return 0;
}


/* TABLE FORMAT
 * - comma separated values
 * - values look like floats
 * - one header row
 * - score column is ignored (at the moment)
 */

/* EXAMPLE TABLE
X,Y,Z,T,Score
1157.0,883.0,516.0,11.0,0.095082399999999997
1141.0,1017.0,522.0,11.0,0.1546526
1469.0,927.0,439.0,11.0,0.1219951
1493.0,955.0,441.0,11.0,0.161437
1403.0,717.0,495.0,11.0,0.1027024
1303.0,757.0,514.0,11.0,0.1086265
1411.0,1000.0,418.0,11.0,0.097495819999999997
1470.0,1071.0,423.0,11.0,0.14146919999999999
2041.0,1057.0,524.0,12.0,0.1061257
748.0,31.0,536.0,12.0,0.097755770000000006
1332.0,916.0,535.0,12.0,0.1722504
*/
