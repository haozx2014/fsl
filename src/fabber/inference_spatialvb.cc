/*  inference_spatialvb.cc - implementation of VB with spatial priors

    Adrian Groves, FMRIB Image Analysis Group

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

#include "easylog.h"
using namespace Utilities;
#include "inference_spatialvb.h"
#include "convergence.h"

void SpatialVariationalBayes::Setup(ArgsType& args)
{
    Tracer_Plus tr("SpatialVariationalBayes::Setup");
    // Call parent to do most of the setup
    VariationalBayesInferenceTechnique::Setup(args);
    
    try {
        spatialDims = convertTo<int>(
            args.ReadWithDefault("spatial-dims", "2"));
    } catch (invalid_argument&){
        throw Invalid_option("--spatial-dims= must have an integer parameter");
    }
    
    if (spatialDims < 0 || spatialDims > 3)
        throw Invalid_option("--spatial-dims= must take 0, 1, 2, or 3");
        
    if (spatialDims == 1)
        LOG_ERR("--spatial-dims=1 is very weird... I hope you're just testing something!\n");
}

void SpatialVariationalBayes::DoCalculations(const DataSet& allData)
{
  Tracer_Plus tr("SpatialVariationalBayes::DoCalculations");
  const Matrix& data = allData.GetVoxelData();
  const int Nvoxels = data.Ncols();
  // Rows are volumes
  // Columns are (time) series
  // num Rows is size of (time) series
  // num Cols is size of volumes       
  const int Nparams = model->NumParams();

  // Added to diagonal to make sure the spatial precision matrix
  // doesn't become singular -- and isolated voxels behave sensibly. 
  const double tiny = 1e-6;

  // Sanity checks:

  if (data.Nrows() != model->NumOutputs())
    throw Invalid_option("Data length (" 
      + stringify(data.Nrows())
      + ") does not match model's output length ("
      + stringify(model->NumOutputs())
      + ")!");

  assert(resultMVNs.empty()); // Only call DoCalculations once

  // Initialization:

  // Make the neighbours[] lists
  CalcNeighbours(allData.GetMask());

  // Make each voxel's distributions
  vector<NoiseModel*> noiseVox; // these change
  vector<MVNDist> fwdPriorVox;
  vector<MVNDist> fwdPosteriorVox;
  vector<LinearizedFwdModel> linearVox;

  { Tracer_Plus tr("SpatialVariationalBayes::DoCalculations - initialization");

  noiseVox.resize(Nvoxels, NULL);
  fwdPriorVox.resize(Nvoxels, *initialFwdPrior);

  fwdPosteriorVox.resize(Nvoxels, *initialFwdPosterior);
  linearVox.resize(Nvoxels, LinearizedFwdModel(model) );
  resultMVNs.resize(Nvoxels, NULL);

  // Note that there's a bunch of wasted memory by using a vector of 
  // noisemodels rather than storing the parameters separately.
  // In particular, the Ar1c cache class is replicated for each voxel,
  // but only one set really needs to be stored at a time.

  for (int v = 1; v <= Nvoxels; v++)
    {
      linearVox[v-1].ReCentre( fwdPosteriorVox[v-1].means );
      noiseVox[v-1] = noise->Clone();
      noiseVox[v-1]->Precalculate( data.Column(v) );
    }
  } // end tracer  

  // Make the spatial normalization parameters
  //akmean = 0*distsMaster.theta.means + 1e-8;
  DiagonalMatrix akmean(Nparams); akmean = 1e-8;

  const double F = 1234.5678; // no sensible updates yet
  if (needF)
    throw Invalid_option("Can't calculate free energy, so stop asking for it!");

  conv->Reset();

  // Allow calculation to continue even with bad voxels
  // Note that this is a bad idea when spatialDims>0, because the bad voxel
  // will drag its neighbours around... but can never recover!  Maybe a more
  // sensible approach is to reset bad voxels to the prior on each iteration.
  //  vector<bool> keepGoing;
  //  keepGoing.resize(Nvoxels, true);

  do {
    conv->DumpTo(LOG);
    conv->DumpTo(cout);
    for (int v = 1; v <= Nvoxels; v++)
      {
	//if (keepGoing.at(v-1)) try{
	//cout << fwdPosteriorVox[v-1].means.t();


	// Apply the spatial priors
	// from simple_do_vb_ar1c_spatial.m
 
	//priors(vox).theta.precisions = diag(20*akmean) ...
	//  + priorsMaster.theta.precisions;
    { Tracer_Plus tr("SpatialVariationalBayes::DoCalculations - spatial priors");

        //n1 = neighbours{vox};
        //n2 = vertcat(neighbours{n1});
        //n2(n2==vox) = []; % avoid self-neighbouring
        //weight8 = 8*length(n1);
        //if weight8 ~= 0
        //    % Main component: 8*immediate neighbours
        //    [distsPrevious(n1).theta]; [ans.means];
        //    contrib8 = 8*sum(ans,2);
        //else
        //    contrib8 = 0;
        //end

	double weight8 = 0; // weighted +8
	ColumnVector contrib8(Nparams); contrib8 = 0.0;
	for (vector<int>::iterator nidIt = neighbours[v-1].begin();
	     nidIt != neighbours[v-1].end(); nidIt++) 
	     // iterate over neighbour ids
	  {
	    int nid = *nidIt;
	    const MVNDist& neighbourPost = fwdPosteriorVox[nid-1];
	    contrib8 += 8 * neighbourPost.means;
	    weight8 += 8;
	  }

	
        //weight12 = -1*length(n2);
        //if weight12 ~= 0
        //    % 2nd-order: -1 or -2 * second neighbours (straight/diagonal)
        //    [distsPrevious(n2).theta]; [ans.means];
        //    contrib12 = -sum(ans,2);
        //else
        //    contrib12 = 0;
        //end

	double weight12 = 0; // weighted -1, may be duplicated
	ColumnVector contrib12(Nparams); contrib12 = 0.0;
	for (vector<int>::iterator nidIt = neighbours2[v-1].begin();
	     nidIt != neighbours2[v-1].end(); nidIt++) 
	     // iterate over neighbour ids
	  {
	    int nid = *nidIt;
	    const MVNDist& neighbourPost = fwdPosteriorVox[nid-1];
	    contrib12 += -neighbourPost.means;
	    weight12 += -1;
	  }

        //if length(n1) > 0 % not isolated
        //    priors(vox).theta.means = (contrib8 + contrib12) / (weight8 + weight12);
        //else
        //    priors(vox).theta.means = 0*priorsMaster.theta.means;
        //end

    // Set prior mean & precisions
     
    int nn = neighbours[v-1].size();
    DiagonalMatrix spatialPrecisions = akmean * ( (nn+tiny)*(nn+tiny) + nn );
    
    fwdPriorVox[v-1].SetPrecisions(
      initialFwdPrior->GetPrecisions() + spatialPrecisions );

    ColumnVector mTmp(Nparams);
    
	if (weight8 != 0)
	  mTmp = (contrib8+contrib12)/(weight8+weight12);
	else
	  mTmp = 0;

// equivalent, when non-spatial priors are very weak:
//    fwdPriorVox[v-1].means = mTmp; 

    fwdPriorVox[v-1].means =
        fwdPriorVox[v-1].GetCovariance() *  
            (spatialPrecisions * mTmp 
            + initialFwdPrior->GetPrecisions() * initialFwdPrior->means);

    } // tracer: spatial priors

        LOG << "Voxel " << v << " of " << Nvoxels << endl;
	noiseVox[v-1]->UpdateTheta( fwdPosteriorVox[v-1], fwdPriorVox[v-1],
				 linearVox[v-1], data.Column(v) );
	
	linearVox[v-1].ReCentre( fwdPosteriorVox[v-1].means );
	
	noiseVox[v-1]->UpdateNoise( fwdPosteriorVox[v-1], linearVox[v-1],
				    data.Column(v) );
	//} catch (...) 
	//{keepGoing[v-1] = false; cout << "Bad Voxel! " << v << endl;} 
      }
   
  { Tracer_Plus tr("SpatialVariationalBayes::DoCalculations - spatial norm update");
  // Update spatial normalization term

  //  for k=1:length(distsMaster.theta.means)
  DiagonalMatrix gk(Nparams); //gk = 0.0/0.0;
  for (int k = 1; k <= Nparams; k++)
  {
    ColumnVector wk(Nvoxels); //wk = 0.0/0.0;
    DiagonalMatrix sigmak(Nvoxels); //sigmak=0.0/0.0;
    for (int v = 1; v <= Nvoxels; v++)
	{
  
//      for i=1:length(dists)
//        wk(i) = dists(i).theta.means(k);
        wk(v) = fwdPosteriorVox[v-1].means(k);
//        inv(dists(i).theta.precisions); sigmak(i) = ans(k,k);
        sigmak(v,v) = fwdPosteriorVox[v-1].GetCovariance()(k,k);
//      end
	}
//    wk = wk(:);
//    tmp1 = trace(diag(sigmak)*S'*S)

//cout << "fwdPosteriorVox[0].GetCovariance():" << fwdPosteriorVox[0].GetCovariance();

for (int v=1; v<=Nvoxels; v++) 


    assert(sigmak==sigmak);

 
    // used to make the matrix non-singular, e.g. for unneighboured voxels
    
    double tmp1 = 0.0;
    for (int v = 1; v <= Nvoxels; v++)
    {
        int nn = neighbours[v-1].size();
//cout << v << ": " << nn << "," << wk(v) << "," << sigmak(v,v) << endl;
        tmp1 += sigmak(v,v) * ( (nn+tiny)*(nn+tiny) + nn );
    }

//    tmp2 = wk'*S'*S*wk
    double tmp2 = 0.0;
    ColumnVector Swk = tiny * wk;
    for (int v = 1; v <= Nvoxels; v++)
    {
	  for (vector<int>::iterator v2It = neighbours[v-1].begin();
	       v2It != neighbours[v-1].end(); v2It++)
	    {
	      Swk(v) += wk(v) - wk(*v2It);
	    }
	}
    tmp2 = (Swk.t() * Swk).AsScalar();

//cout << "tmp1=" << tmp1 << ", tmp2=" << tmp2 << endl;
//cout << Swk.t();

//  gk(k) = 1/(0.5*tmp1 + 0.5*tmp2 + 1/10)
    gk(k,k) = 1/(0.5*tmp1 + 0.5*tmp2 + 0.1);
//  end

  }
//hk = nVoxels/2 + 1
//akmean = gk .* hk % around 6e-6 would be okay for stdev of 100 (Qn)
  
//  cout << gk;
  assert(gk==gk);

  akmean = gk * (Nvoxels*0.5 + 1.0);
  
//  cout << "Using fixed values for akmean!!\n";
//  akmean << 0.0011 << 9.9704e-04 << 2.1570e-07 
//    << 5.2722 << 0.0177 << 0.5041;

//  cout << akmean;
  

  } // tracer: spatial norm update

    // next iteration:
  } while (!conv->Test( F ));

  // Phew!

  for (int v = 1; v <= Nvoxels; v++)
    resultMVNs[v-1] = new MVNDist(
      fwdPosteriorVox[v-1], noiseVox[v-1]->GetResultsAsMVN() );

  //  throw Runtime_error("Spatial VB isn't implemented yet, so use --method=vb instead\n"); 
}


inline int binarySearch(const ColumnVector& data, int num)
{
  int first = 1, last = data.Nrows();

  while (first <= last)
    {
      int test = (first+last)/2;

      if (data(test) < num)
	first = test + 1;
      else if (data(test) > num)
	last = test - 1;
      else if (data(test) == num)
	return test;
      else
	assert(false);
    }
  return -1;
}



void SpatialVariationalBayes::CalcNeighbours(const Volume& mask)
{
  Tracer_Plus tr("SpatialVariationalBayes::CalcNeighbours");

  const ColumnVector& preThresh = mask.getPreThresholdPositions();
  const VolumeInfo& info = mask.getInfo();
  const int nVoxels = preThresh.Nrows();

  vector<int> delta;
  if (spatialDims >= 1)
  { 
    delta.push_back(1); // next row
    delta.push_back(-1); // prev row
  }
  if (spatialDims >= 2)
  {
    delta.push_back(info.x); // next column
    delta.push_back(-info.x); // prev column
  }
  if (spatialDims >= 3)
  {
    delta.push_back(info.x*info.y); // next slice
    delta.push_back(-info.x*info.y); // prev slice
  }
  
  LOG << "Using offsets of " << delta 
      << "\nNote these have not been checked!\n"
      << "Also note that they may wrap around if both edges are in mask\n";
  
  neighbours.resize(nVoxels);
  
  for(int vid = 1; vid <= nVoxels; vid++) // voxel id
    {
      for (unsigned n = 0; n < delta.size(); n++) // neighbour-search number
	{
	  int id = binarySearch(preThresh, int(preThresh(vid)) + delta[n]);
	  if (id > 0)
	    neighbours.at(vid-1).push_back(id);
	}
    }
  
  // Neighbours-of-neighbours, excluding self, and duplicated if there 
  // are two routes to get there (diagonally connected)
  neighbours2.resize(nVoxels);
  for(int vid = 1; vid <= nVoxels; vid++)
    {
      for (unsigned n1 = 0; n1 < neighbours.at(vid-1).size(); n1++)
	{
	  int n1id = neighbours[vid-1].at(n1);
	  int checkNofN = 0;
	  for (unsigned n2 = 0; n2 < neighbours.at(n1id-1).size(); n2++)
	    {
	      int n2id = neighbours[n1id-1].at(n2);
	      if (n2id != vid)
		neighbours2[vid-1].push_back(n2id);
	      else
		checkNofN++;
	    }
	  assert(checkNofN == 1);
	  // Each of this voxel's neighbours must have this voxel 
	  // as a neighbour.
	}
    }    

  LOG << "Neighbours are (for dims==" << spatialDims << ")\n";
  for (int v = 1; v <= nVoxels; v++)
    LOG << v << ": " << neighbours.at(v-1) 
	<< "-" << neighbours2.at(v-1) << endl;
}
