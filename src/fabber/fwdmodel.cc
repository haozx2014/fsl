/*  fwdmodel.cc - base class for generic forward models

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

#include "fwdmodel.h"

#include <sstream> 
#include "easylog.h"
 
string FwdModel::ModelVersion() const
{
  return "No version info available.";
  // You should overload this function in your FwdModel class
  
  // Something like this:
  // return " $ I d $ "; // without the spaces
  // CVS will automatically replace this with version information that looks
  // like this: $Id: fwdmodel.cc,v 1.10 2007/08/01 15:55:35 adriang Exp $
  // It should probably go in your .cc file, not the header.
}

int FwdModel::NumOutputs() const
{
    ColumnVector params, result;
    params.ReSize(NumParams()); params = 1; // probably safer than 0
    Evaluate(params, result);
    return result.Nrows();
}

void FwdModel::DumpParameters(const ColumnVector& params, 
                              const string& indent) const
{
    LOG << indent << "Parameters:" << endl;
    vector<string> names;
    NameParams(names);
    
    for (int i = 1; i <= NumParams(); i++)
        LOG << indent << "  " << names[i-1] << " == " << params(i) << endl;
        
    LOG << indent << "Total of " << NumParams() << " parameters." << endl;
}

#include "fwdmodel_simple.h"
#include "fwdmodel_quipss2.h"
#include "fwdmodel_asl_eagle.h"
#include "fwdmodel_asl_grase.h"
#include "fwdmodel_asl_buxton.h"
#include "fwdmodel_q2tips.h"
// Add your models here

FwdModel* FwdModel::NewFromName(const string& name, ArgsType& args)
{
    // Update this to add your own models to the code
    
    if (name == "simple")
    {
        LOG_ERR("WARNING: the 'simple' forward model has several hard-coded filenames and probably won't work outside of FMRIB\n");
        return new SimpleFwdModel(args);
    }
    else if (name == "quipss2")
    {
        return new Quipss2FwdModel(args);
    } 
    else if (name == "q2tips")
    {
        return new Q2tipsFwdModel(args);
    } 
    else if (name == "eagle")
      {
	return new EagleFwdModel(args);
      }
    else if (name == "grase")
      {
	return new GraseFwdModel(args);
      }
    else if (name == "buxton")
      {
	return new BuxtonFwdModel(args);
      }
    // Your models go here!
    else
    {
        throw Invalid_option("Unrecognized forward model '" + name + "'");
    }  
}

// If you want usage information from --help --model=yourmodel, add it below.
void FwdModel::ModelUsageFromName(const string& name, ArgsType& args)
{
    // Update this to add your own models to the code
    
//    if (name == "simple")
//    {
//        return new SimpleFwdModel(args);
//    }
//    else 
    if (name == "quipss2")
    {
        Quipss2FwdModel::ModelUsage();
    } 
    if (name == "q2tips")
    {
        cout << "Note: --model=q2tips uses exactly the same options as --model=quipss2:\n";
	Quipss2FwdModel::ModelUsage();
    } 
    else if (name == "grase")
      {
	GraseFwdModel::ModelUsage();
      }
    else if (name == "buxton")
      {
	BuxtonFwdModel::ModelUsage();
      }
//    else if (name == "eagle")
//      {
//    return new EagleFwdModel(args);
//      }
    // Your models go here!
    else
    {
        cout << "No model usage information available for --model=" 
            << name << endl;
    }   
}

