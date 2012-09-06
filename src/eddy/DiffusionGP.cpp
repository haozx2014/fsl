// Definitions of class to make Gaussian-Process
// based predictions about diffusion data.
//
// DiffusionGP.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "EddyHelperClasses.h"
#include "EddyUtils.h"
#include "DiffusionGP.h"

using namespace EDDY;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class DiffusionGP
//
// 
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

DiffusionGP::DiffusionGP(const std::string&            scans_fname,
			 const std::string&            var_mask_fname,
                         const std::vector<DiffPara>&  dpars)
{
  NEWIMAGE::volume4D<float> scans;
  EddyUtils::read_DWI_volume4D(scans,scans_fname,dpars);
  for (int s=0; s<scans.tsize(); s++) _sptrs.push_back(boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(scans[s])));
  _pop = true;
  _dpars = EddyUtils::GetDWIDiffParas(dpars);
  NEWIMAGE::volume<float> var_mask; NEWIMAGE::read_volume(var_mask,var_mask_fname);
  _vars = VarianceCalculator::GetWMVariances(_sptrs,var_mask,_dpars,_grps,_grpb);
  _grp = reorganize_grps();
  _s2n = VarianceCalculator::GetNoiseVariance(_sptrs,var_mask,_dpars);
  _K = calculate_K();
  _valid = true;
}

bool DiffusionGP::IsPopulated() const
{
  if (_pop) return(true);
  else {
    _pop = true;
    for (unsigned int i=0; i<_sptrs.size(); i++) {
      if (!_sptrs[i]) { _pop = false; break; }
    }
  }
  return(_pop);
}

void DiffusionGP::SetNoOfScans(unsigned int n)
{
  if (n == _sptrs.size()) return; // No change
  else if (n > _sptrs.size()) {   // If increasing size
    _sptrs.resize(n,boost::shared_ptr<NEWIMAGE::volume<float> >()); // New elements populated by NULL
    _dpars.resize(n); // New elements populated with (arbitrary) [1 0 0] bvec
    _pop = false;
    _valid = false;
  }
  else { // If decreasing size
    _sptrs.resize(n);
    _dpars.resize(n);
    _valid = false;
    if (_pop==false) { // _pop may potentially go from false to true
      _pop = IsPopulated();
    }
  }
  return;
}

//
// Ideally one would like to check if array is ready populated
// when adding a scan. But that would potentially make it 
// non-thread reentrant so I have choosen not to.
//
void DiffusionGP::AddScan(const NEWIMAGE::volume<float>& scan,
                          const DiffPara&                dp)
{
  if (_sptrs.size() && !NEWIMAGE::samesize(*_sptrs[0],scan)) throw EddyException("DiffusionGP::AddScan: Wrong image dimension");
  _sptrs.push_back(boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(scan)));
  _dpars.push_back(dp);
  _valid = false;
}

void DiffusionGP::SetScan(const NEWIMAGE::volume<float>& scan,
                          const DiffPara&                dp,
                          unsigned int                   indx)
{
  if (int(indx) > (int(_sptrs.size())-1)) throw EddyException("DiffusionGP::SetScan: Invalid image index");
  if (_sptrs.size() && _sptrs[0] && !NEWIMAGE::samesize(*_sptrs[0],scan)) throw EddyException("DiffusionGP::SetScan: Wrong image dimension");
  if (_sptrs[indx] && dp != _dpars[indx]) throw EddyException("DiffusionGP::SetScan: You cannot change shell or direction of scan");
  _sptrs[indx] = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(scan));
  _dpars[indx] = dp;
}

void DiffusionGP::UpdateK(const NEWIMAGE::volume<float>& var_mask)
{
  if (IsPopulated()) {
    _vars = VarianceCalculator::GetWMVariances(_sptrs,var_mask,_dpars,_grps,_grpb);
  _grp = reorganize_grps();
  _s2n = VarianceCalculator::GetNoiseVariance(_sptrs,var_mask,_dpars);
  _K = calculate_K();
  _valid = true;
  }
  else throw EddyException("DiffusionGP::UpdateK: You cannot update K until predictor is fully populated");
}

NEWIMAGE::volume<float> DiffusionGP::Predict(double       sigman2,
					     unsigned int indx) const
{
  if (!_valid) throw EddyException("DiffusionGP::Predict: invalid predictor");

  NEWMAT::RowVector       K_row = extract_K_row(indx);
  NEWMAT::Matrix          iK = get_iK(sigman2);
  NEWMAT::RowVector       pvec = K_row*iK;
  NEWMAT::ColumnVector    y(_sptrs.size());
  NEWIMAGE::volume<float> pvol = *_sptrs[0];
  pvol = 0;
  for (int k=0; k<_sptrs[0]->zsize(); k++) {
    for (int j=0; j<_sptrs[0]->ysize(); j++) {
      for (int i=0; i<_sptrs[0]->xsize(); i++) {
	if (get_y(i,j,k,y)) { 
	  std::vector<double> ymeans = get_group_means(y);
          subtract_group_means(ymeans,y);
	  pvol(i,j,k) = static_cast<float>((pvec*y).AsScalar()) + ymeans[_grp[indx]];
        }
      }
    }
  }  
  return(pvol);  
}

NEWIMAGE::volume<float> DiffusionGP::Predict(double            sigman2,
					     const DiffPara&   dpars) const
{
  if (!_valid) throw EddyException("DiffusionGP::Predict: invalid predictor");
  printf("DiffusionGP::Predict(DiffPara dp): Not fully implemented yet\n");

  NEWMAT::RowVector       K_row = calculate_K_row(dpars);
  NEWMAT::Matrix          iK = get_iK(sigman2);
  NEWMAT::RowVector       pvec = K_row*iK;
  NEWMAT::ColumnVector    y(_sptrs.size());
  NEWIMAGE::volume<float> pvol = *_sptrs[0];
  pvol = 0.0;
  for (int k=0; k<_sptrs[0]->zsize(); k++) {
    for (int j=0; j<_sptrs[0]->ysize(); j++) {
      for (int i=0; i<_sptrs[0]->xsize(); i++) {
	if (get_y(i,j,k,y)) pvol(i,j,k) = static_cast<float>((pvec*y).AsScalar());
      }
    }
  }  
  return(pvol);  
}

NEWIMAGE::volume4D<float> DiffusionGP::PredictAll(double sigman2) const
{
  if (!_valid) throw EddyException("DiffusionGP::PredictAll: invalid predictor");

  NEWIMAGE::volume4D<float> pvols(_sptrs[0]->xsize(),_sptrs[0]->ysize(),_sptrs[0]->zsize(),_sptrs.size());
  pvols = 0;
  for (unsigned int i=0; i<_sptrs.size(); i++) {
    pvols[i] = this->Predict(sigman2,i);
  }
  return(pvols);
}

void DiffusionGP::WriteDebugInfo(const std::string& fname) const
{
  std::ofstream fout;
  fout.open(fname.c_str(),ios::out|ios::trunc);

  fout << "There are " << GetWMVariances().size() << " groups of b-values" << endl;
  fout << "The bvalues are:";
  for (unsigned i=0; i<GetbValues().size(); i++) { fout << "  " << GetbValues()[i]; }
  fout << endl;
  fout << "Their variances are:";
  for (unsigned int i=0; i<GetWMVariances().size(); i++) { fout << "  " << GetWMVariances()[i]; }
  fout << endl;
  fout << "The noise variance is  " << GetNoiseVariance() << endl;
  fout << "The _grps consists of " << Get_grps().size() << " vectors" << endl;
  for (unsigned int i=0; i<Get_grps().size(); i++) {
    fout << "Vector " << i << " is:" << endl;
    for (unsigned int j=0; j<Get_grps()[i].size(); j++) { fout << "  " << Get_grps()[i][j]; }
    fout << endl;
  }
  fout << "The _grp is:" << endl;
  for (unsigned int i=0; i<Get_grp().size(); i++) { fout << "  " << Get_grp()[i]; }
  fout << endl;
  fout << "The K-matrix is:" << endl;
  fout << GetKMatrix() << endl;

  fout.close();
}

bool DiffusionGP::get_y(// Input
			unsigned int           i,
			unsigned int           j,
			unsigned int           k,
			// Output
			NEWMAT::ColumnVector&  y) const
{
  for (unsigned int t=0; t<_sptrs.size(); t++) {
    if (!(*_sptrs[t])(i,j,k)) return(false);
    else y(t+1) = (*_sptrs[t])(i,j,k);
  }
  return(true);
}

void DiffusionGP::subtract_group_means(// Input
				       const std::vector<double>& means,
				       // Input//output
				       NEWMAT::ColumnVector&      y) const
{
  for (unsigned int g=0; g<_grps.size(); g++) {
    for (unsigned int i=0; i<_grps[g].size(); i++) {
      y(_grps[g][i]+1) -= means[g];
    }
  }
  return;
}

void DiffusionGP::add_group_means(// Input
				  const std::vector<double>& means,
				  // Input//output
				  NEWMAT::ColumnVector&      y) const
{
  for (unsigned int g=0; g<_grps.size(); g++) {
    for (unsigned int i=0; i<_grps[g].size(); i++) {
      y(_grps[g][i]+1) += means[g];
    }
  }
  return;
}

std::vector<double> DiffusionGP::get_group_means(const NEWMAT::ColumnVector& y) const
{
  std::vector<double> means(_grps.size(),0);

  for (unsigned int g=0; g<_grps.size(); g++) {
    for (unsigned int i=0; i<_grps[g].size(); i++) {
      means[g] += y(_grps[g][i]+1);
    }
    means[g] /= _grps[g].size();
  }
  return(means);
}

std::vector<int> DiffusionGP::reorganize_grps() const
{
  std::vector<int> grpi(_sptrs.size());
  for (unsigned g=0; g<_grps.size(); g++) {
    for (unsigned int s=0; s<_grps[g].size(); s++) {
      grpi[_grps[g][s]] = g;
    }
  }
  return(grpi);
}

NEWMAT::CroutMatrix DiffusionGP::get_LUK(double sigman2) const
{
  NEWMAT::IdentityMatrix sn2(_sptrs.size());
  sn2 *= sigman2;
  NEWMAT::CroutMatrix LUK = (_K + sn2);

  return(LUK);
}

NEWMAT::Matrix DiffusionGP::get_iK(double sigman2) const
{
  NEWMAT::IdentityMatrix sn2(_sptrs.size());
  sn2 *= sigman2;
  NEWMAT::Matrix iK = (_K + sn2).i();

  return(iK);
}

NEWMAT::Matrix DiffusionGP::calculate_K() const
{
  NEWMAT::Matrix K(_sptrs.size(),_sptrs.size());

  for (unsigned int i=0; i<_sptrs.size(); i++) {
    for (unsigned int j=i; j<_sptrs.size(); j++) {
      double var = sqrt(_vars[_grp[i]])*sqrt(_vars[_grp[j]]);
      K(i+1,j+1) = 0.54*var + 1.54*var*pow(abs((_dpars[i].bVec().t()*_dpars[j].bVec()).AsScalar()),2.0);
      K(j+1,i+1) = K(i+1,j+1);
    } 
  }
  return(K);
}

NEWMAT::RowVector DiffusionGP::extract_K_row(unsigned int indx) const
{
  NEWMAT::RowVector  K_row(_sptrs.size());
  for (unsigned int i=0; i<_sptrs.size(); i++) K_row(i+1) = _K(indx+1,i+1);

  return(K_row);
}

NEWMAT::RowVector DiffusionGP::calculate_K_row(const DiffPara& dpar) const
{
  // Check that we have a valid b-value
  unsigned int g;
  for (g=0; g<_grps.size(); g++) { if (EddyUtils::AreInSameShell(dpar,_dpars[_grps[g][0]])) break; }
  if (g == _grps.size()) throw EddyException("DiffusionGP::calculate_K_row: bval not in existing shell");

  // Make row
  NEWMAT::RowVector  K_row(_sptrs.size());
  for (unsigned int i=0; i<_sptrs.size(); i++) {
    double var = sqrt(_vars[g])*sqrt(_vars[_grp[i]]);
    K_row(i+1) = 0.54*var + 1.54*var*pow(abs((_dpars[i].bVec().t()*dpar.bVec()).AsScalar()),2.0);
  }
  return(K_row);
} 

// }}} End of fold.


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class VarianceCalculator
//
// 
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

std::vector<double> VarianceCalculator::GetWMVariances(// Input
						       const std::vector<boost::shared_ptr<NEWIMAGE::volume<float> > >&   sptrs,
						       const NEWIMAGE::volume<float>&                                     mask,
						       const std::vector<DiffPara>&                                       dpars,
						       // Output
						       std::vector<vector<int> >&                                         grps,
						       std::vector<double>&                                               grpb)
{
  // Erode mask a little to avoid high variance voxels along edge
  NEWIMAGE::volume<float> kernel = NEWIMAGE::box_kernel(3,3,3);
  NEWIMAGE::volume<float> emask = NEWIMAGE::morphfilter(mask,kernel,std::string("erode"));
  emask = NEWIMAGE::morphfilter(emask,kernel,std::string("erode"));

  grps = get_groups(dpars,grpb);
  std::vector<double>  var(grps.size(),0.0);
  
  for (unsigned int g=0; g<grps.size(); g++) {
    // Calculate volume of voxel-wise variances
    NEWIMAGE::volume4D<float> tmp4D(sptrs[0]->xsize(),sptrs[0]->ysize(),sptrs[0]->zsize(),grps[g].size());
    for (unsigned int s=0; s<grps[g].size(); s++) {
      tmp4D[s] = *sptrs[grps[g][s]];
    }
    NEWIMAGE::volume<float> varvol = NEWIMAGE::variancevol(tmp4D);
    // Chose 20% of intracerebral voxels with highest variance
    float thres = varvol.percentile(0.8,emask);
    NEWIMAGE::volume<float> meanmask = NEWIMAGE::binarise(varvol,thres);
    var[g] = varvol.mean(meanmask);
  } 
  return(var);
}

double VarianceCalculator::GetNoiseVariance(// Input
					    const std::vector<boost::shared_ptr<NEWIMAGE::volume<float> > >&   sptrs,
					    const NEWIMAGE::volume<float>&                                     mask,
					    const std::vector<DiffPara>&                                       dpars)
{
  // Dilate mask a little to make sure we are outside object, then invert.
  NEWIMAGE::volume<float> kernel = NEWIMAGE::box_kernel(3,3,3);
  NEWIMAGE::volume<float> emask = NEWIMAGE::morphfilter(mask,kernel,std::string("dilate"));
  emask = NEWIMAGE::morphfilter(emask,kernel,std::string("dilate"));
  emask *= -1; emask += 1;

  std::vector<vector<int> > grps = get_groups(dpars);
  std::vector<double>       gvar(grps.size(),0.0);
  
  for (unsigned int g=0; g<grps.size(); g++) {
    // Calculate volume of voxel-wise variances
    NEWIMAGE::volume4D<float> tmp4D(sptrs[0]->xsize(),sptrs[0]->ysize(),sptrs[0]->zsize(),grps[g].size());
    for (unsigned int s=0; s<grps[g].size(); s++) {
      tmp4D[s] = *sptrs[grps[g][s]];
    }
    NEWIMAGE::volume<float> varvol = NEWIMAGE::variancevol(tmp4D);
    // Chose voxels from top half of volume to avoid extracerebral tissue.
    int n = 0;
    for (int k=varvol.zsize()/2; k<varvol.zsize()-1; k++) {
      for (int j=2; j<varvol.ysize()-2; j++) {
	for (int i=2; i<varvol.xsize()-2; i++) { if (emask(i,j,k)) { gvar[g] += varvol(i,j,k); n++; } }
      }
    }
    gvar[g] /= n;
  }

  double var = 0.0;
  for (unsigned int g=0; g<grps.size(); g++) var += gvar[g];

  return(var);
}

std::vector<vector<int> > VarianceCalculator::get_groups(// Input 
							 const std::vector<DiffPara>&   dpars,
							 // Output
							 std::vector<double>&           grpb)
{
  std::vector<vector<int> > grps;
  std::vector<DiffPara>     templates;
  // First pass to sort out how many different b-values/shells there are
  templates.push_back(dpars[0]);
  for (unsigned int i=1; i<dpars.size(); i++) {
    unsigned int j;
    for (j=0; j<templates.size(); j++) { if (EddyUtils::AreInSameShell(templates[j],dpars[i])) break; }
    if (j == templates.size()) templates.push_back(dpars[i]);
  }
  // Sort them in ascending order
  std::sort(templates.begin(),templates.end());
  grps.resize(templates.size());
  grpb.resize(templates.size(),0.0);
  // Populate vectors of indicies
  for (unsigned int i=0; i<dpars.size(); i++) {
    for (unsigned int j=0; j<templates.size(); j++) {
      if (EddyUtils::AreInSameShell(templates[j],dpars[i])) { 
	grps[j].push_back(i); 
        grpb[j] += dpars[i].bVal();
	break; 
      }
    }
  }
  for (unsigned int i=0; i<grpb.size(); i++) grpb[i] /= double(grps[i].size());

  return(grps);
}

// }}} End of fold.

