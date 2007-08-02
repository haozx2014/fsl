/*  fwdmodel_asl_grase.h - Implements the dual echo GRASE resting state ASL model

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

#include "fwdmodel_asl_grase.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string GraseFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_asl_grase.cc,v 1.8 2007/07/27 21:58:30 adriang Exp $";
}

void GraseFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("GraseFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Tissue bolus perfusion
     prior.means(tiss_index()) = 0;
     precisions(tiss_index(),tiss_index()) = 1e-6;

    // Tissue bolus transit delay
     prior.means(tiss_index()+1) = 0.7;
     precisions(tiss_index()+1,tiss_index()+1) = 0.5;
    
    
    // Tissue bolus length
     if (infertau) {
       prior.means(tau_index()) = seqtau;
       precisions(tau_index(),tau_index()) = 0.5;
     }

    // Arterial Perfusion & bolus delay

    if (inferart)
      {
	int aidx = art_index();
	prior.means(aidx) = 0;
	prior.means(aidx) = 0.7;
	precisions(aidx+1,aidx+1) = 0.5;
	precisions(aidx,aidx) = 1e-6;
      }
 
    // T1 & T1b
    if (infert1) {
      int tidx = t1_index();
      prior.means(tidx) = t1;  
      prior.means(tidx+1) = t1b; 
      precisions(tidx,tidx) = 100;
      precisions(tidx+1,tidx+1) = 100;
    }

    if (inferinveff) {
      prior.means(inveff_index()) = 0.3;
      precisions(inveff_index(),inveff_index()) = 10;
    }

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;
    
}    
    
    

void GraseFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("GraseFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
  for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
    }

     // sensible limits on transit times
  if (params(tiss_index()+1)>timax-0.2) { paramcpy(tiss_index()+1) = timax-0.2; }
  if (inferart > 0 ) {
    if (params(art_index()+1)>timax-0.2) { paramcpy(art_index()+1) = timax-0.2; }
  }

  // parameters that are inferred - extract and give sensible names
  float ftiss;
  float delttiss;
  float tautiss;
  float fblood;
  float deltblood;
  float T_1;
  float T_1b;
  float inveffslope;

  ftiss=paramcpy(tiss_index());
  delttiss=paramcpy(tiss_index()+1);

  if (infertau) { tautiss=paramcpy(tau_index()); }
  else { tautiss = seqtau; }
  
  if (inferart) {
    fblood=paramcpy(art_index());
    deltblood=paramcpy(art_index());
  }
  else {
    fblood = 0;
    deltblood = 0;
  }

  if (infert1) {
    T_1 = paramcpy(t1_index());
    T_1b = paramcpy(t1_index()+1);
  }
  else {
    T_1 = t1;
    T_1b = t1b;
  }

  if (inferinveff) { inveffslope = paramcpy(inveff_index()); }
  else { inveffslope = 0; }    



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
    float ti;
    result.ReSize(tis.Nrows()*repeats);

    for(int it=1; it<=tis.Nrows(); it++)
      {
	ti = tis(it);
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
	    kctissue = kctissue * (1 + -inveffslope * R - inveffslope*(ti-delttiss)); //deals with the case where we have a slope in inversion efficency
	    if (kctissue<0) {kctissue = 0; } //dont allow negative values
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
	  { 
	    kcblood = fblood * exp(-ti/T_1b); 
	    kcblood = kcblood * (1 - inveffslope*(ti-deltblood)); //deals with the case where we have a slope in inversion efficency
	    if (kcblood<0) { kcblood = 0; } //dont allow negative values
	  }

	/* output */
	// loop over the repeats
	for (int rpt=1; rpt<=repeats; rpt++)
	  {
	    result( (it-1)*repeats+rpt ) = kctissue + kcblood;
	  }


      }


  return;
}

GraseFwdModel::GraseFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      repeats = convertTo<int>(args.Read("repeats")); // number of repeats in data
      t1 = convertTo<double>(args.ReadWithDefault("t1","1.3"));
      t1b = convertTo<double>(args.ReadWithDefault("t1b","1.5"));
      infertau = args.ReadBool("infertau"); // infer on bolus length?
      infert1 = args.ReadBool("infert1"); //infer on T1 values?
      inferart = args.ReadBool("inferart"); //infer on arterial compartment?
      inferinveff = args.ReadBool("inferinveff"); //infer on a linear decrease in inversion efficiency?
      if (inferinveff) { 
	seqtau = 10; //effectively infinite
	if (infertau) { throw Invalid_option("You may not use --infertau with --inferinveff"); }
      }
      else if (infertau)
	{seqtau = convertTo<double>(args.Read("tau")); } //if inferrring upon tau then an intital vaule MUST be set
      else
	{seqtau = convertTo<double>(args.ReadWithDefault("tau","10"));} //bolus length as set by sequence (default of 10 is effectively infinite
	
 
      // Deal with tis
      tis.ReSize(1); //will add extra values onto end as needed
      tis(1) = atof(args.Read("ti1").c_str());
      
      while (true) //get the rest of the tis
	{
	  int N = tis.Nrows()+1;
	  string tiString = args.ReadWithDefault("ti"+stringify(N), "stop!");
	  if (tiString == "stop!") break; //we have run out of tis
	 
	  // append the new ti onto the end of the list
	  ColumnVector tmp(1);
	  tmp = convertTo<double>(tiString);
	  tis &= tmp; //vertical concatenation

	}
      timax = tis.Maximum(); //dtermine the final TI

      // add information about the parameters to the log
      LOG << "    Data parameters: #repeats = " << repeats << ", t1 = " << t1 << ", t1b = " << t1b;
      LOG << ", bolus length (tau) = " << seqtau << endl ;
      if (infertau) {
	LOG << "Infering on bolus length " << endl; }
      if (inferart) {
	LOG << "Infering on artertial compartment " << endl; }
      if (infert1) {
	LOG << "Infering on T1 values " << endl; }
      if (inferinveff) {
	LOG << "Infering on Inversion Efficency slope (value of bolus length will assumed infinite, even if set by --tau= " << endl; }
      LOG << "TIs: ";
      for (int i=1; i <= tis.Nrows(); i++)
	LOG << tis(i) << " ";
      LOG << endl;
	  
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

void GraseFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=grase:\n"
       << "Required parameters:\n"
       << "--repeats=<no. repeats in data>\n"
       << "--ti1=<first_inversion_time_in_seconds>\n"
       << "--ti2=<second_inversion_time>, etc...\n"
       << "Optional arguments:\n"
       << "--tau=<temporal_bolus_length> (default 10s if --infertau not set)\n"
       << "--t1=<T1_of_tissue> (default 1.3)\n"
       << "--t1b=<T1_of_blood> (default 1.5)\n"
       << "--infertau (to infer on bolus length)\n"
       << "--inferart (to infer on arterial compartment)\n"
       << "--infert1 (to infer on T1 values)\n"
    ;
}

void GraseFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void GraseFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("ftiss");
  names.push_back("delttiss");
  if (infertau)
    names.push_back("tautiss");
  if (inferart) {
    names.push_back("fblood");
    names.push_back("deltblood");
  }
  if (infert1) {
    names.push_back("T_1");
    names.push_back("T_1b");
  }
  if (inferinveff) {
    names.push_back("Inveffslope");
  }
}
