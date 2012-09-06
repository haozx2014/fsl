// Declarations of class to make Gaussian-Process
// based predictions about diffusion data.
//
// DiffusionGP.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef DiffusionGP_h
#define DiffusionGP_h

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "EddyHelperClasses.h"
#include "DWIPredictionMaker.h"

namespace EDDY {

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class DiffusionGP
//
// Class used to make predictions about diffusion data.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  class DiffusionGP : public DWIPredictionMaker
{
public:
  DiffusionGP() : _pop(true), _valid(false) {};
  DiffusionGP(const std::string&                 scans_fname,
              const std::string&                 var_mask_fname,
              const std::vector<DiffPara>&       dpars);

  bool IsPopulated() const;                  // Returns true if all data present
  bool IsValid() const { return(_valid); }   // Returns true if ready to make predictions
  void SetNoOfScans(unsigned int n);
  void AddScan(const NEWIMAGE::volume<float>& scan,
               const DiffPara&                dp);
  void SetScan(const NEWIMAGE::volume<float>& scan,
               const DiffPara&                dp,
               unsigned int                   indx);
  void UpdateK(const NEWIMAGE::volume<float>& mask);
  void EvaluateModel(const NEWIMAGE::volume<float>& mask) { if (!_valid) this->UpdateK(mask); }
  NEWIMAGE::volume<float> Predict(double       sigman2,
                                  unsigned int indx) const;
  NEWIMAGE::volume<float> Predict(unsigned int indx) const { return(this->Predict(_s2n,indx)); }
  NEWIMAGE::volume<float> Predict(double            sigman2,
				  const DiffPara&   dpars) const;
  NEWIMAGE::volume<float> Predict(const DiffPara&   dpars) const { return(this->Predict(_s2n,dpars)); }
  NEWIMAGE::volume4D<float> PredictAll(double sigman2) const;
  NEWIMAGE::volume4D<float> PredictAll() const { return(this->PredictAll(_s2n)); }

  // Here starts routines that should be used for debugging only.
  const NEWMAT::Matrix& GetKMatrix() const { if (_valid) return(_K); else throw EddyException("DiffusionGP::GetKMatrix: invalid predictor"); }
  const std::vector<double>& GetWMVariances() const { if (_valid) return(_vars); else throw EddyException("DiffusionGP::GetWMVariances: invalid predictor"); }
  double GetNoiseVariance() const { if (_valid) return(_s2n); else throw EddyException("DiffusionGP::GetNoiseVariance: invalid predictor"); }
  const std::vector<double>& GetbValues() const { if (_valid) return(_grpb); else throw EddyException("DiffusionGP::GetbValues: invalid predictor"); }
  const std::vector<vector<int> >& Get_grps() const { if (_valid) return(_grps); else throw EddyException("DiffusionGP::Get_grps: invalid predictor"); }
  const std::vector<int>& Get_grp() const { if (_valid) return(_grp); else throw EddyException("DiffusionGP::Get_grp: invalid predictor"); }
  void WriteDebugInfo(const std::string& fname) const;

private:
  std::vector<boost::shared_ptr<NEWIMAGE::volume<float> > > _sptrs; // Pointers to the scans
  // NEWIMAGE::volume4D<float>    _scans;  // The scans
  std::vector<DiffPara>                                     _dpars;  // Diffusion parameters
  std::vector<double>                                       _vars;   // The white matter variance of the different shells
  std::vector<vector<int> >                                 _grps;   // Indices (0-offs) to indicate what shell a scan belongs to
  std::vector<int>                                          _grp;    // Group index for each scan. Same info as in _grps, organised differently.
  std::vector<double>                                       _grpb;   // b-values of the shells (bvals within a range of 100 = same)
  double                                                    _s2n;    // Tentative noise variance.
  NEWMAT::Matrix                                            _K;      // The resulting K-matrix (in the order given by scans)
  mutable bool                                              _pop;    // Tells if all data is present
  bool                                                      _valid;  // Tells if it is ready to make predictions

  bool get_y(// Input
	     unsigned int           i,
	     unsigned int           j,
	     unsigned int           k,
	     // Output
	     NEWMAT::ColumnVector&  y) const;
  void subtract_group_means(// Input
			    const std::vector<double>& means,
			    // Input//output
			    NEWMAT::ColumnVector&      y) const;
  void add_group_means(// Input
		       const std::vector<double>& means,
		       // Input//output
		       NEWMAT::ColumnVector&      y) const;
  std::vector<double> get_group_means(const NEWMAT::ColumnVector& y) const;
  std::vector<int>    reorganize_grps() const;
  NEWMAT::CroutMatrix get_LUK(double sigman2) const;
  NEWMAT::Matrix      get_iK(double sigman2) const;
  NEWMAT::Matrix      calculate_K() const;
  NEWMAT::RowVector   extract_K_row(unsigned int indx) const;
  NEWMAT::RowVector   calculate_K_row(const DiffPara&    dpar) const;  
};

// }}} End of fold.

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class VarianceCalculator
//
// Helper Class used to calculate variance (across g-vectors) within
// a shell. Calculations are confined to and averaged over white
// matter voxels.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class VarianceCalculator
{
public:
  static std::vector<double> GetWMVariances(// Input
					    const std::vector<boost::shared_ptr<NEWIMAGE::volume<float> > >&   sptrs,
					    const NEWIMAGE::volume<float>&                                     mask,
					    const std::vector<DiffPara>&                                       dpars,
					    // Output
					    std::vector<vector<int> >&                                         grps,
					    std::vector<double>&                                               grpb);

static double GetNoiseVariance(// Input
			       const std::vector<boost::shared_ptr<NEWIMAGE::volume<float> > >&   sptrs,
			       const NEWIMAGE::volume<float>&                                     mask,
			       const std::vector<DiffPara>&                                       dpars);
private:
  static std::vector<vector<int> > get_groups(const std::vector<DiffPara>&   dpars) { 
    std::vector<double> skrutt; return(get_groups(dpars,skrutt));
  }
  static std::vector<vector<int> > get_groups(// Input 
					      const std::vector<DiffPara>&   dpars,
					      // Output
					      std::vector<double>&           grpb);

};

// }}} End of fold.

} // End namespace EDDY

#endif // End #ifndef DiffusionGP_h
