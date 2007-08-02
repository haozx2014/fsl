/*  inference_vb.cc - VB inference technique class declarations

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

#include "inference_vb.h"
#include "convergence.h"

void VariationalBayesInferenceTechnique::Setup(ArgsType& args) 
{ 
  Tracer_Plus tr("VariationalBayesInferenceTechnique::Setup");

  // Call ancestor, which does most of the real work
  InferenceTechnique::Setup(args);

  // Load up initial prior and initial posterior
  MVNDist* loadPrior = new MVNDist( model->NumParams() );
  MVNDist* loadPosterior = new MVNDist( model->NumParams() );  
  
  string filePrior = args.ReadWithDefault("fwd-initial-prior", "modeldefault");
  string filePosterior = args.ReadWithDefault("fwd-initial-posterior", "modeldefault");
  if (filePrior == "modeldefault" || filePosterior == "modeldefault")
      model->HardcodedInitialDists(*loadPrior, *loadPosterior);
  if (filePrior != "modeldefault") loadPrior->Load(filePrior);
  if (filePosterior != "modeldefault") loadPosterior->Load(filePosterior);

  // Make these distributions constant:
  assert(initialFwdPrior == NULL);
  assert(initialFwdPosterior == NULL);    
  initialFwdPrior = loadPrior;
  initialFwdPosterior = loadPosterior;
  loadPrior = loadPosterior = NULL; // now, only accessible as consts.

  // Maximum iterations allowed:
  string maxIterations = args.ReadWithDefault("max-iterations","10");
  if (maxIterations.find_first_not_of("0123456789") != string::npos)
    throw Invalid_option("--convergence=its=?? parameter must be a positive number");    
  int its = atol(maxIterations.c_str());
  if (its<=0)
    throw Invalid_option("--convergence=its=?? paramter must be positive");
  
  // Figure out convergence-testing method:
  string convergence = args.ReadWithDefault("convergence", "maxits");
  if (convergence == "maxits")
    conv = new CountingConvergenceDetector(its);
  else if (convergence == "pointzeroone")
    conv = new FchangeConvergenceDetector(its, 0.01);
  else if (convergence == "freduce")
    conv = new FreduceConvergenceDetector(its, 0.01);
  else if (convergence == "trialmode")
    conv = new TrialModeConvergenceDetector(its, 10, 0.01);
  else
    throw Invalid_option("Unrecognized convergence detector: '" 
                           + convergence + "'");

  // Figure out if F needs to be calculated every iteration
  printF = args.ReadBool("print-free-energy");
  needF = conv->UseF() || printF;

  haltOnBadVoxel = !args.ReadBool("allow-bad-voxels");
  if (haltOnBadVoxel)
    LOG << "Note: numerical errors in voxels will cause the program to halt.\n"
        << "Use --allow-bad-voxels (with caution!) to keep on calculating.\n";
  else
    LOG << "Using --allow-bad-voxels: numerical errors in a voxel will\n"
	<< "simply stop the calculation of that voxel.\n"
	<< "Check log for 'Going on to the next voxel' messages.\n"
	<< "Note that you should get very few (if any) exceptions like this;"
	<< "they are probably due to bugs or a numerically unstable model.";
  
}

void VariationalBayesInferenceTechnique::DoCalculations(const DataSet& allData) 
{
  Tracer_Plus tr("VariationalBayesInferenceTechnique::DoCalculations");

  const Matrix& data = allData.GetVoxelData();
  // Rows are volumes
  // Columns are (time) series
  // num Rows is size of (time) series
  // num Cols is size of volumes       
  int Nvoxels = data.Ncols();
  if (data.Nrows() != model->NumOutputs())
    throw Invalid_option("Data length (" 
      + stringify(data.Nrows())
      + ") does not match model's output length ("
      + stringify(model->NumOutputs())
      + ")!");

  assert(resultMVNs.empty()); // Only call DoCalculations once
  resultMVNs.resize(Nvoxels, NULL);

  // Reverse order (to test that everything is being properly reset):
  //for (int voxel = Nvoxels; voxel >=1; voxel--) 
  
  for (int voxel = 1; voxel <= Nvoxels; voxel++)
    {
      ColumnVector y = data.Column(voxel);
      NoiseModel* noiseVox = noise->Clone();
//      NoiseModel* noiseVox = noise;

// I think this problem has been resolved, but an interesting comment for later...
//
// Problem with using a new noise model every time: much poorer 
// convergence.  Re-using previous voxel's posterior improves things
// enormously.  Can a slightly more informative initial posterior
// on alpha do just as well?? 
// Also: why doesn't the MATLAB one have this problem?
// Increased iterations to match (30) --> problem is considerably reduced
// but still present.

//      noise->Dump();
//      noiseVox->Dump();
      
      LOG_ERR("  Voxel " << voxel << " of " << Nvoxels << endl); 
      //  << " sumsquares = " << (y.t() * y).AsScalar() << endl;
      double F = 1234.5678;

      const MVNDist fwdPrior( *initialFwdPrior );
      MVNDist fwdPosterior( *initialFwdPosterior );
      LinearizedFwdModel linear( model );
      try
	{
	  linear.ReCentre( fwdPosterior.means );
	  conv->Reset();

	  noiseVox->Precalculate( y );
    
	  do 
	    {
	      if (needF) F = noiseVox->CalcFreeEnergy( fwdPosterior, fwdPrior, linear, y );
	      if (printF) LOG << "      Fbefore == " << F << endl;

//	      fwdPosterior.Dump();
	     
              // Save old values if called for
	      if ( conv->NeedSave() ) {
		noiseVox->SaveParams( fwdPosterior );
	      }
 
	      // Theta update
	      noiseVox->UpdateTheta( fwdPosterior, fwdPrior, linear, y );


      
	      if (needF) F = noiseVox->CalcFreeEnergy( fwdPosterior, fwdPrior, linear, y );
	      if (printF) LOG << "      Ftheta == " << F << endl;
	      
	      // Linearization update
	      linear.ReCentre( fwdPosterior.means );
	      if (needF) F = noiseVox->CalcFreeEnergy( fwdPosterior, fwdPrior, linear, y );
	      if (printF) LOG << "      Flin == " << F << endl;
//noiseVox->Dump();	      
	      // Alpha & Phi updates
	      noiseVox->UpdateNoise( fwdPosterior, linear, y );

	      // Test of NoiseModel cloning:
	      // NoiseModel* tmp = noise; noise = tmp->Clone(); delete tmp;
	      
	      if (needF) F = noiseVox->CalcFreeEnergy( fwdPosterior, fwdPrior, linear, y );
	      if (printF) LOG << "      Fnoise == " << F << endl;

	    }           
	  while ( !conv->Test( F ) );

	  // Revert to old values at last stage if required
	  if ( conv-> NeedRevert() ) {
	    noiseVox->RevertParams( fwdPosterior );
	  }
	  
	  conv->DumpTo(LOG, "    ");
	} 
      catch (const overflow_error& e)
	{
	  LOG_ERR("    Went infinite!  Reason:" << endl
		  << "      " << e.what() << endl);
	  //todo: write garbage or best guess to memory/file
	  if (haltOnBadVoxel) throw;
	  LOG_ERR("    Going on to the next voxel." << endl);
	}
      catch (Exception)
	{
	  LOG_ERR("    NEWMAT Exception in this voxel:\n"
		  << Exception::what() << endl);
	  if (haltOnBadVoxel) throw;
	  LOG_ERR("    Going on to the next voxel." << endl);  
	}
      catch (...)
	{
	  LOG_ERR("    Other exception caught in main calculation loop!!\n");
	    //<< "    Use --halt-on-bad-voxel for more details." << endl;
	  if (haltOnBadVoxel) throw;
	  LOG_ERR("    Going on to the next voxel" << endl);
	}
      
      try {

	LOG << "    Final parameter estimates are:" << endl;
	linear.DumpParameters(fwdPosterior.means, "      ");
	
	assert(resultMVNs.at(voxel-1) == NULL);
	resultMVNs.at(voxel-1) = new MVNDist(
	  fwdPosterior, noiseVox->GetResultsAsMVN() );

      } catch (...) {
	// Even that can fail, due to results being singular
	LOG << "    Can't give any sensible answer for this voxel; outputting zero +- identity\n";
	MVNDist* tmp = new MVNDist();
	tmp->SetSize(fwdPosterior.means.Nrows()
		    + noiseVox->GetResultsAsMVN().means.Nrows());
	tmp->SetCovariance(IdentityMatrix(tmp->means.Nrows()));
	resultMVNs.at(voxel-1) = tmp;
      }
      delete noiseVox; noiseVox = NULL;
//      assert(noiseVox == noise);
    }
}

VariationalBayesInferenceTechnique::~VariationalBayesInferenceTechnique() 
{ 
  delete conv;
  delete initialFwdPrior;
  delete initialFwdPosterior;
}




