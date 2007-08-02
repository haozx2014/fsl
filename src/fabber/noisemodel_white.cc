/*  noisemodel_ar.cc - Class implementation for the multiple white noise model

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

#include "noisemodel_white.h"
#include "noisemodel.h"
#include <stdexcept>
#include "miscmaths/miscmaths.h"
using namespace MISCMATHS;
using namespace Utilities;
#include "easylog.h"

WhiteNoiseModel* WhiteNoiseModel::Clone() const
{
  Tracer_Plus tr("WhiteNoiseModel::Clone");
  WhiteNoiseModel* clone = new WhiteNoiseModel(phiPattern);
  clone->phis = phis;
  clone->phisPrior = phisPrior;

  return clone;
}
  

void WhiteNoiseModel::LoadPrior(const string& filename)
{
  Tracer_Plus tr("WhiteNoiseModel::LoadPrior");
  int nPhis = Qis.size();
  assert(nPhis > 0);
  phis.resize(nPhis);
  phisPrior.resize(nPhis);

  if (filename == "hardcoded")
    for (int i = 1; i <= nPhis; i++)
      {
	phis[i-1].b = phisPrior[i-1].b = 1e6;
	phis[i-1].c = phisPrior[i-1].c = 1e-6;

        // Bit of black magic... tiny initial noise precision seems to help
        phis[i-1].b = 1e-8;
        phis[i-1].c = 50;
      }
  else
    throw Invalid_option("Only --noise-priors=hardcoded is supported at present");

}

const MVNDist WhiteNoiseModel::GetResultsAsMVN() const
{
  Tracer_Plus tr("WhiteNoiseModel::GetResultsAsMVN");
  MVNDist phiMVN( phis.size() );
  SymmetricMatrix phiVars( phis.size() );
  phiVars = 0;
  for (unsigned i = 1; i <= phis.size(); i++)
    {
      phiMVN.means(i) = phis[i-1].CalcMean();
      phiVars(i,i) = phis[i-1].CalcVariance();
    }
  phiMVN.SetCovariance(phiVars);
  return phiMVN;
}

void WhiteNoiseModel::Dump(const string indent) const
{
  Tracer_Plus tr("WhiteNoiseModel::Dump");
  LOG << indent << "White noise distribution, " << phis.size() << " separate distributions." << endl;
//  LOG << indent << "  Prior:" << endl;
  DumpPrior(indent + "  ");
//  LOG << indent << "  Posterior:" << endl;
  DumpPosterior(indent + "  ");
}

void WhiteNoiseModel::DumpPrior(const string indent) const
{
  Tracer_Plus tr("WhiteNoiseModel::DumpPrior");
  for (unsigned i = 0; i < phisPrior.size(); i++)
    {
      LOG << indent << "Phi_" << i+1 << " prior: ";
      phisPrior[i].Dump();
    }
}

void WhiteNoiseModel::DumpPosterior(const string indent) const
{
  Tracer_Plus tr("WhiteNoiseModel::DumpPosterior");
  for (unsigned i = 0; i < phisPrior.size(); i++)
    {
      LOG << indent << "Phi_" << i+1 << ": ";
      phis[i].Dump();
    }
}

WhiteNoiseModel::WhiteNoiseModel(const string& pattern) 
  : phiPattern(pattern)
{ 
  Tracer_Plus tr("WhiteNoiseModel::WhiteNoiseModel");
  assert(phiPattern.length() > 0);
  MakeQis(phiPattern.length()); // a quick way to validate the input
}


void WhiteNoiseModel::MakeQis(int dataLen) const
{
  Tracer_Plus tr("WhiteNoiseModel::MakeQis");
  if (!Qis.empty() && Qis[0].Nrows() == dataLen) 
    return;  // Qis are already up-to-date

  // Read the pattern string into a vector pat
  const int patternLen = phiPattern.length();
  assert(patternLen > 0);
  if (patternLen > dataLen)
    throw Invalid_option(
      "Pattern length exceeds data length... this is probably a mistake");

  vector<int> pat;
  int nPhis = 0;

  for (int i = 1; i <= patternLen; i++)
    {
      char c = phiPattern[i-1];
      int n;
      if (c >= '1' && c <= '9') // Index from 1, so 0 is not allowed
	n = (c - '0');
      else if (c >= 'A' && c <= 'Z')
	n = (c - 'A' + 10);
      else if (c >= 'a' && c <= 'z')
	n = (c - 'a' + 10);
      else
	throw Invalid_option(
	          string("Invalid character in pattern: '") + c + "'");
      pat.push_back(n);
      if (nPhis<n) nPhis = n;
    }


  //  cout << "Pat is " << pat << endl;

  // Extend pat to the the full data length by repeating
  while ((int)pat.size() < dataLen)
    pat.push_back(pat.at(pat.size()-patternLen));

  LOG << "Pattern of phis used is " << pat << endl;

  // Regenerate the Qis
  Qis.clear();
  DiagonalMatrix zeroes(dataLen);
  zeroes = 0.0;
  Qis.resize(nPhis, zeroes); // initialize all to zero

  for (int d = 1; d <= dataLen; d++)
    {
    //    Qis[pat[d-1]-1](d,d) = 1;
    Qis.at(pat.at(d-1)-1)(d,d) = 1;
    }

  // Sanity checking
  for (int i = 1; i <= nPhis; i++)
    if (Qis[i-1].Trace() < 1.0) // this phi is never used
      throw Invalid_option(
	 "At least one Phi was unused! This is probably a bad thing.");
}

void WhiteNoiseModel::UpdateNoise(
  	const MVNDist& theta,
  	const LinearFwdModel& linear,
  	const ColumnVector& data)
{
  Tracer_Plus tr("WhiteNoiseModel::UpdateNoise");
  const Matrix& J = linear.Jacobian();
  ColumnVector k = data - linear.Offset() + J*(linear.Centre() - theta.means);
  
  // check the Qis are valid
  MakeQis(data.Nrows());

  // Update each phi distribution in turn
  for (unsigned i = 1; i <= Qis.size(); i++)
    {
      const DiagonalMatrix& Qi = Qis[i-1];
      double tmp = 
	(k.t() * Qi * k).AsScalar()
	+ (theta.GetCovariance() * J.t() * Qi * J).Trace();
      
      phis[i-1].b =
	1/( tmp*0.5 + 1/phisPrior[i-1].b);
      
      double nTimes = Qi.Trace(); // number of data points for this dist
      assert(nTimes == int(nTimes)); // should be an integer

      phis[i-1].c = 
	(nTimes-1)*0.5 + phisPrior[i-1].c;
    }
}

void WhiteNoiseModel::UpdateTheta(
  	MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& linear,
  	const ColumnVector& data) const
{
  Tracer_Plus tr("WhiteNoiseModel::UpdateTheta");

  const ColumnVector &ml = linear.Centre();
  const ColumnVector &gml = linear.Offset();
  const Matrix &J = linear.Jacobian();

  // Make sure Qis are up-to-date
  MakeQis(data.Nrows());

  // Marginalize over phi distributions
  DiagonalMatrix X(data.Nrows()); 
  X = 0;
  for (unsigned i = 1; i <= Qis.size(); i++)
    X += Qis[i-1] * phis[i-1].CalcMean();

  // Calculate Lambda & Lambda*m (without priors)
  SymmetricMatrix Ltmp;
  Ltmp << J.t() * X * J;  
  // use << instead of = because this is considered a lossy assignment
  // (since NEWMAT isn't smart enough to know J'*X*J is always symmetric)
  ColumnVector mTmp = J.t() * X * (data - gml + J*ml);

  // Update Lambda and m (including priors)
  theta.SetPrecisions(thetaPrior.GetPrecisions() + Ltmp);
  theta.means = theta.GetCovariance()
    * ( mTmp + thetaPrior.GetPrecisions() * thetaPrior.means );

  // Error checking
  LogAndSign chk = theta.GetPrecisions().LogDeterminant();
  if (chk.Sign() <= 0)
    LOG << "Note: In UpdateTheta, theta precisions aren't positive-definite: "
	<< chk.Sign() << ", " << chk.LogValue() << endl;
}

double WhiteNoiseModel::CalcFreeEnergy(
	const MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& linear,
  	const ColumnVector& data) const
{

  int nPhis = Qis.size();

  // Calculate some matrices we will need
  const Matrix &J = linear.Jacobian();
  ColumnVector k = data - linear.Offset()
    + J * (linear.Centre() - theta.means);
  const SymmetricMatrix& Linv = theta.GetCovariance();

  // some values we will need
  int nTimes = data.Nrows(); //*NB assume that each row is an individual time point
  int nTheta = theta.means.Nrows();

  // The following is based on noisemodel_ar::CalcFreeEnergy, modified to remove ar parts - MAC 11-7-2007
  // Some modifications have been made for consistency with (MAC)varbayes2.m - these are noted

  // calcualte individual aprts of the free energy
  double expectedLogThetaDist = //bits arising from the factorised posterior for theta
    +0.5 * theta.GetPrecisions().LogDeterminant().LogValue()
    -0.5 * nTheta * (log(2*M_PI) + 1);

  double expectedLogPhiDist = 0; //bits arising fromt he factorised posterior for phi
  vector<double> expectedLogPosteriorParts(10); //bits arising from the likelihood
   for (int i = 0; i<10; i++) 
     expectedLogPosteriorParts[i] = 0;
  
   for (int i = 0; i < nPhis; i++) 
     {
      double si = phis[i].b;
      double ci = phis[i].c;
      double siPrior = phisPrior[i].b;
      double ciPrior = phisPrior[i].c;

      expectedLogPhiDist +=
	-gammaln(ci) - ci*log(si) - ci
	+(ci-1)*(digamma(ci)+log(si));
      
      expectedLogPosteriorParts[0] += 
	(digamma(ci)+log(si)) * ( (nTimes)*0.5 + ciPrior - 1 ); // nTimes rather than nTimes-1
      
      expectedLogPosteriorParts[9] += 
	-2*gammaln(ciPrior) -2*ciPrior*log(siPrior) - si*ci/siPrior;
     } 

   expectedLogPosteriorParts[1] = 0; //*NB not required
  
  expectedLogPosteriorParts[2] =
    -0.5 * (k.t() *  k).AsScalar() 
    -0.5 * (J.t() * J * Linv).Trace(); //*NB remove Qsum
  
  expectedLogPosteriorParts[3] =
    +0.5 * thetaPrior.GetPrecisions().LogDeterminant().LogValue();
  
  expectedLogPosteriorParts[4] = 
    -0.5 * (
            (theta.means - thetaPrior.means).t() 
            * thetaPrior.GetPrecisions() 
            * (theta.means - thetaPrior.means)
	    ).AsScalar();
  
  expectedLogPosteriorParts[5] =
    -0.5 * (Linv * thetaPrior.GetPrecisions()).Trace();

  expectedLogPosteriorParts[6] = 0; //*NB not required
  
  expectedLogPosteriorParts[7] = 0; //*NB not required  
  
  expectedLogPosteriorParts[8] = 0; //*NB not required
 
  // Assemble the parts into F
  double F =  
    - expectedLogThetaDist 
    - expectedLogPhiDist;
  
  for (int i=0; i<10; i++) 
    F += expectedLogPosteriorParts[i];

  // Error checking
  if (! (F - F == 0) )
  {
    LOG_ERR("expectedLogThetaDist == " << expectedLogThetaDist << endl);
    LOG_ERR("expectedLogPhiDist == " << expectedLogPhiDist << endl);
    //LOG_ERR("expectedLogPosteriorParts == " << expectedLogPosteriorParts << endl);
    throw overflow_error("Non-finite free energy!");
  }
  
  return F;
}

void WhiteNoiseModel::SaveParams(const MVNDist& theta) {
  Tracer_Plus tr("WhiteNoiseModel::SaveParams");
  // save the current values of parameters 
  int nPhis = phis.size();
  assert(nPhis > 0);
  phissave.resize(phis.size());
  for (int i = 1; i <= nPhis; i++)
      {
	phissave[i-1].b = phis[i-1].b ;
	phissave[i-1].c = phis[i-1].c ;
      }
  thetasave = theta;
}

void WhiteNoiseModel::RevertParams(MVNDist& theta) {
  Tracer_Plus tr("WhiteNoiseModel::RevertParams");
  int nPhis = phis.size();
  for (int i = 1; i <= nPhis; i++)
      {
	phis[i-1].b = phissave[i-1].b ;
	phis[i-1].c = phissave[i-1].c ;
      }
  theta = thetasave;
}
