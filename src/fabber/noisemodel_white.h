/*  noisemodel_ar.h - Class declaration for the multiple white noise model

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 4.0 (c) 2007, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/1112. */

#include "noisemodel.h"
#include "dist_gamma.h"

class WhiteNoiseModel : public NoiseModel {
 public:

  virtual WhiteNoiseModel* Clone() const;
  // makes a new identical copy of this object

  // Loading priors & saving results
  virtual void LoadPrior( const string& filename );
  virtual const MVNDist GetResultsAsMVN() const;

  // Debugging stuff
  virtual void Dump(const string indent = "") const;
  virtual void DumpPrior(const string indent = "") const;
  virtual void DumpPosterior(const string indent = "") const;

  // Constructor/destructor
  WhiteNoiseModel(const string& pattern);
  // pattern says which phi distribution to use for each data points; this
  // string is repeated as necessary to make up the data length. e.g. for 
  // dual-echo data, pattern = "12".  after 9, use letters (A=10, ...)
  // Simplest case: pattern = "1".

  virtual ~WhiteNoiseModel() { return; }

  // Do all the calculations
  virtual void UpdateNoise(
  	const MVNDist& theta,
  	const LinearFwdModel& model,
  	const ColumnVector& data);

  virtual void UpdateTheta(
  	MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& model,
  	const ColumnVector& data) const;

  virtual double CalcFreeEnergy(
	const MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& model,
  	const ColumnVector& data) const;

 
 void SaveParams(const MVNDist& theta) ;		 
 void RevertParams(MVNDist& theta) ;

 protected:
  const string phiPattern;

  // Diagonal matrices, indicating which data points use each phi
  mutable vector<DiagonalMatrix> Qis; // mutable because it's used as a cache
  void MakeQis(int dataLen) const;

  vector<GammaDist> phis;
  vector<GammaDist> phisPrior;

  vector<GammaDist> phissave;
  MVNDist thetasave;

};
