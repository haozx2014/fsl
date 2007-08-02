/*  fwdmodel_asl_eagle.cc - Implements the EAGLE model

    Michael Chappell, FMRIB Image Analysis Group

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

#include "fwdmodel_asl_eagle.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

void EagleFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("EagleFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Tissue bolus perfusion
    prior.means(1) = 0;
    precisions(1,1) = 1e-6;

    // Tissue bolus transit delay
    prior.means(2) = 0.7;
    precisions(2,2) = 0.5;

    // Tissue bolus length
    prior.means(3) = 2;
    precisions(3,3) = 0.5;

    // Arterial Perfusion
    prior.means(4) = 0;
    precisions(4,4) = 1e-6;

    // Arterial bolus delay
    prior.means(5) = 0.7;
    precisions(5,5) = 0.5;

    // T1 ***1.5T***
    prior.means(6) = 0.9;
    precisions(6,6) = 100;
    
    // T1b ***1.5T***
    prior.means(7) = 1.1;
    precisions(7,7) = 100;

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;
    
}    
    
    

void EagleFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("EagleFwdModel::Evaluate");
    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
    for (int i=1;i<8;i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
    }
     // sensible limit on times
    int tpar []= { 2, 5 };
    for (int i=0;i<2;i++) {
      if (params(tpar[i])>2.5) { paramcpy(tpar[i]) = 2.5; }
      }

  // parameters that are inferred - extract and give sensible names
    float ftiss=paramcpy(1);
    float delttiss=paramcpy(2);
    float tautiss=paramcpy(3);
    float fblood=paramcpy(4);
    float deltblood=paramcpy(5);
    float T_1 = paramcpy(6);
    float T_1b = paramcpy(7);



    float lambda = 0.9;

    float T_1app = 1/( 1/T_1 + 0.01/lambda );
    float R = 1/T_1app - 1/T_1b;

    /* According to EAGLE GRASE sequence bolus length is current TI - 0.1s (assuming infite length 'true' bolus)
       However, also allow here bolus length to be finite as recorded in tautiss
       NB tautiss is the 'true' bolus length, tau is what the tissue actually sees as a result of the sequence
    */

    float tau1 = delttiss;
    float tau; //bolus length as seen by kintic curve
    float tau2; //phase that occurs after end of bolus in kc

    float F=0;
    float kctissue;
    float kcblood;


    // loop over tis
    float ti; int it;
    result.ReSize(14);

    result = -9999;

    for(ti=0.2, it=1; it<=14; ti += 0.2, it++)
      {
	F = 2*ftiss/lambda * exp(-ti/T_1app);

	/* deal with bolus length*/
	if(tautiss < ti - 0.1)
	  { tau = tautiss; }
	else
	{ tau = ti - 0.1; }
	tau2 = tau1 + tau;
	/* tissue contribution */
	if(ti < tau1)
	  { kctissue = 0;}
	else if(ti >= tau1 && ti <= tau2)
	  {
	    kctissue = F/R * (exp(R*ti) - exp(R*tau1));
	  }
	else /*(ti > tau2)*/
	  {kctissue = F/R * (exp(R*tau2) - exp(R*tau1)); }
	
	/* arterial contribution */
	if(ti < deltblood)
	  { 
	    //kcblood = 0;
	    kcblood = fblood * exp(-deltblood/T_1b) * (0.98 * exp( (ti-deltblood)/0.1 ) + 0.02 * ti/deltblood );
	  }
	else
	  { kcblood = fblood * exp(-ti/T_1b); }

	/* output */
	result(it) = kctissue + kcblood;


      }


    assert( result(14) != -9999 );

  return;
}

EagleFwdModel::EagleFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    string tagPattern;
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

void EagleFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void EagleFwdModel::NameParams(vector<string>& names) const
{
    
}
