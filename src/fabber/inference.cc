/*  inference.cc - General inference technique base class

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

#include "inference.h"


void InferenceTechnique::Setup(ArgsType& args)
{
  Tracer_Plus tr("InferenceTechnique::Setup");

  // Pick models
  model = FwdModel::NewFromName(args.Read("model"), args);
  assert( model->NumParams() > 0 );
  LOG_ERR("    Forward Model version:\n      " 
	  << model->ModelVersion() << endl);

  noise = NoiseModel::NewFromName(args.Read("noise"), args);
  noise->LoadPrior(args.ReadWithDefault("noise-prior","hardcoded"));
  noise->Dump("  ");

  saveModelFit = args.ReadBool("save-model-fit");
  saveResiduals = args.ReadBool("save-residuals");
}



void InferenceTechnique::SaveResults(const DataSet& data) const
{
  Tracer_Plus tr("InferenceTechnique::SaveResults");
    LOG << "    Preparing to save results..." << endl;
    
    // Save the resultMVNs as two NIFTI files
    // Note: I should probably use a single NIFTI file with
    // NIFTI_INTENT_NORMAL -- but I can't find the detailed 
    // documentation!  (Ordering for a multivariate norm).

    VolumeSeries means;
    VolumeSeries precisions;

    const Volume& mask = data.GetMask();

    int nVoxels = resultMVNs.size();
    assert(nVoxels > 0 && resultMVNs.at(0) != NULL);   
    int nParams = resultMVNs.at(0)->means.Nrows();
    
    means.ReSize(nParams, nVoxels);
    precisions.ReSize(nParams*(nParams+1)/2, nVoxels);
    for (int vox = 1; vox <= nVoxels; vox++)
      {
        means.Column(vox) = resultMVNs.at(vox-1)->means;
        precisions.Column(vox) = 
	  resultMVNs.at(vox-1)->GetCovariance().AsColumn();
	// AsColumn uses row ordering on the lower triangular part, 
	// as NIFTI_INTENT_SYMMATRIX specifies: (1,1) (2,1) (2,2) (3,1)...
      }

    //LOG << "    Writing means..." << endl;

    // Save means
    VolumeInfo info = mask.getInfo();
    info.v = means.Nrows();
    // Note: this is an abuse of NIFTI_INTENT_VECTOR
    // vector terms should be in 5th dim, not 4th.
    // Also, can I attach names somehow?
    info.intent_code = NIFTI_INTENT_VECTOR;
    means.setInfo(info);
    means.setPreThresholdPositions(mask.getPreThresholdPositions());
    means.unthresholdSeries();
    means.writeAsFloat(outputDir + "/finalMeans");
    means.CleanUp();
    
    //LOG << "    Writing covariance matrices..." << endl;
      
    // Save precision matrices
    info = mask.getInfo();
    info.v = precisions.Nrows();
    // Note: this is an abuse of NIFTI_INTENT_SYMMATRIX becuase
    // matrix terms should be in 5th dim, not 4th.
    info.intent_code = NIFTI_INTENT_SYMMATRIX;
    info.intent_p1 = nParams;
    precisions.setInfo(info);
    precisions.setPreThresholdPositions(mask.getPreThresholdPositions());
    precisions.unthresholdSeries();
    precisions.writeAsFloat(outputDir + "/finalCovariance");
    precisions.CleanUp();

    // Write the parameter names into paramnames.txt
    
    LOG << "    Writing paramnames.txt..." << endl;
    ofstream paramFile((outputDir + "/paramnames.txt").c_str());
    vector<string> paramNames;
    model->NameParams(paramNames);
    for (unsigned i = 0; i < paramNames.size(); i++)
    {
        LOG << "      " << paramNames[i] << endl;
        paramFile << paramNames[i] << endl;
    }
    paramFile.close();

    LOG << "    Same information using DumpParameters:" << endl;
    ColumnVector indices(model->NumParams());
    for (int i = 1; i <= indices.Nrows(); i++)
        indices(i) = i;
    model->DumpParameters(indices, "      ");

    // Create individual files for each parameter's mean and Z-stat

    for (unsigned i = 1; i <= paramNames.size(); i++)
    {
        VolumeSeries paramMean, paramZstat;
    	paramMean.ReSize(1, nVoxels);
	paramZstat.ReSize(1, nVoxels);

        for (int vox = 1; vox <= nVoxels; vox++)
        {
	    paramMean(1,vox) = resultMVNs[vox-1]->means(i);
            paramZstat(1,vox) =
              paramMean(1,vox) / 
              sqrt(resultMVNs[vox-1]->GetCovariance()(i,i));
        }
    	LOG << "    Writing means..." << endl;

    	// Save paramMeans
    	VolumeInfo info = mask.getInfo();
    	info.v = 1;
    	info.intent_code = NIFTI_INTENT_NONE; // just data values
        paramMean.setInfo(info);
        paramMean.setPreThresholdPositions(mask.getPreThresholdPositions());
        paramMean.unthresholdSeries();
        paramMean.writeAsFloat(outputDir + "/mean_" + paramNames.at(i-1));
        paramMean.CleanUp();

        info.intent_code = NIFTI_INTENT_ZSCORE; // = mean/stdev
        paramZstat.setInfo(info);
        paramZstat.setPreThresholdPositions(mask.getPreThresholdPositions());
        paramZstat.unthresholdSeries();
        paramZstat.writeAsFloat(outputDir + "/zstat_" + paramNames.at(i-1));
        paramZstat.CleanUp();
    }        

    if (saveModelFit || saveResiduals)
    {
        LOG << "    Writing model fit/residuals..." << endl;
        // Produce the model fit and residual volumeserieses

        VolumeSeries modelFit, residuals;
        modelFit.ReSize(model->NumOutputs(), nVoxels);
	ColumnVector tmp;
        for (int vox = 1; vox <= nVoxels; vox++)
        {
            model->Evaluate(resultMVNs.at(vox-1)->means.Rows(1,model->NumParams()), tmp);
            modelFit.Column(vox) = tmp;
        }
        info.v = model->NumOutputs();
        info.intent_code = NIFTI_INTENT_NONE; // just data values
        if (saveResiduals)
        {
            residuals = data.GetVoxelData() - modelFit;
            residuals.setInfo(info);
            residuals.setPreThresholdPositions(mask.getPreThresholdPositions());
            residuals.unthresholdSeries();
            residuals.writeAsFloat(outputDir + "/residuals");
            residuals.CleanUp();
        }
        if (saveModelFit)
        {
            modelFit.setInfo(info);
            modelFit.setPreThresholdPositions(mask.getPreThresholdPositions());
            modelFit.unthresholdSeries();
            modelFit.writeAsFloat(outputDir + "/modelfit");
            modelFit.CleanUp();
        }

    }

    LOG << "    Done writing results." << endl;
}


InferenceTechnique::~InferenceTechnique() 
{ 
  delete model;
  delete noise;
  while (!resultMVNs.empty())
    {
      delete resultMVNs.back();
      resultMVNs.pop_back();
    }
}

#include "inference_vb.h"
#include "inference_spatialvb.h"

// Whenever you add a new class to inference.h, update this too.
InferenceTechnique* InferenceTechnique::NewFromName(const string& method)
{
  Tracer_Plus tr("PickInferenceTechnique");

  if (method == "vb")
    {
      return new VariationalBayesInferenceTechnique;
    }
  else if (method == "spatialvb")
    {
      return new SpatialVariationalBayes;
    }
  else if (method == "")
    {
      throw Invalid_option("Must include the --method=vb or --method=spatialvb option");
    }
  else 
    {
      throw Invalid_option("Unrecognized --method: " + method);
    }
}
