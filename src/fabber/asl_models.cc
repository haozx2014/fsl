/* asl_models.cc Kinetic curve models for ASL

Michael Chappell - IBME & FMRIB Analysis Group

Copyright (C) 2010 University of Oxford */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
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
    innovation@isis.ox.ac.uk quoting reference DE/9564. */

#include "asl_models.h"

namespace OXASL {

// --- Kinetic curve functions ---
//Arterial

  ColumnVector kcblood_nodisp(const ColumnVector& tis, float deltblood, float taub, float T_1b, bool casl=false) {
  //Tracer_Plus tr("OXASL:kcblood_nodisp");
  ColumnVector kcblood(tis.Nrows());
  kcblood=0.0;

  // Non dispersed arterial curve (pASL)
  float ti=0.0;
  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);
      
      if(ti < deltblood)
	{ 
	  kcblood(it) = 2 * exp(-deltblood/T_1b) * (0.98 * exp( (ti-deltblood)/0.05 ) + 0.02 * ti/deltblood );
	  // use a arti§fical lead in period for arterial bolus to improve model fitting
	}
      else if(ti >= deltblood && ti <= (deltblood + taub))
	{ 
	  if (casl) kcblood(it) = 2 * exp(-ti/deltblood);
	  else      kcblood(it) = 2 * exp(-ti/T_1b); 
	}
      else //(ti > deltblood + tau)
	{
	  if (casl) kcblood(it) = 2 * exp(-ti/deltblood);
	  else      kcblood(it) = 2 * exp(-(deltblood+taub)/T_1b);

	  kcblood(it) *= (0.98 * exp( -(ti - deltblood - taub)/0.05) + 0.02 * (1-(ti - deltblood - taub)/5));
	  // artifical lead out period for taub model fitting
	  if (kcblood(it)<0) kcblood(it)=0; //negative values are possible with the lead out period equation
	}

    }

  return kcblood;
}

  ColumnVector kcblood_gammadisp(const ColumnVector& tis, float deltblood, float taub, float T_1b, float s, float p, bool casl=false) {
  //Tracer_Plus tr("OXASL:kcblood_gammadisp");
  ColumnVector kcblood(tis.Nrows());
  kcblood=0.0;

  // Gamma dispersed arterial curve (pASL)

  float k=1+p*s;
  float ti=0.0;

  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);

      if(ti < deltblood)
	{ 
	  kcblood(it) = 0.0;
	}
      else if(ti >= deltblood && ti <= (deltblood + taub))
	{ 
	  if (casl) kcblood(it) = 2 * exp(-deltblood/T_1b);
	  else	    kcblood(it) = 2 * exp(-ti/T_1b);

	    kcblood(it) *= ( 1 - igamc(k,s*(ti-deltblood))  ); 
	}
      else //(ti > deltblood + taub)
	{
	  if (casl) kcblood(it) = 2 * exp(-deltblood/T_1b);
	  else	    kcblood(it) = 2 * exp(-ti/T_1b);
	  kcblood(it) *= ( igamc(k,s*(ti-deltblood-taub)) - igamc(k,s*(ti-deltblood)) ) ; 
	  
	}
      //if (isnan(kcblood(it))) { kcblood(it)=0.0; cout << "Warning NaN in blood KC"; }
    }
  return kcblood;
}

  ColumnVector kcblood_gvf(const ColumnVector& tis, float deltblood, float T_1b, float s, float p, bool casl=false) {
    //Tracer_Plus tr("OXASL:kcblood_gammadisp");
  ColumnVector kcblood(tis.Nrows());
  kcblood=0.0;

  // gamma variate arterial curve
  // NOTES: this model is only suitable for pASL
  //        no taub see below 
  float ti=0.0;

  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);

      if(ti < deltblood)
	{ 
	  kcblood(it) = 0.0;
	}
      else //if(ti >= deltblood) && ti <= (deltblood + taub))
	{ 
	  if (casl) kcblood(it) = 2 * exp(-deltblood/T_1b);
	  else	    kcblood(it) = 2 * exp(-ti/T_1b);

	  kcblood(it) *= gvf(ti-deltblood,s,p); 
	}
      // we do not have bolus duration with a GVF AIF - the duration is 'built' into the function shape
      //else //(ti > deltblood + taub)
      //	{
      //	  kcblood(it) = 0.0 ; 
      //	  
      //	}
    }
  return kcblood;
}

  ColumnVector kcblood_gaussdisp(const ColumnVector& tis, float deltblood, float taub, float T_1b, float sig1, float sig2, bool casl=false) {
    //Tracer_Plus tr("OXASL:kcblood_normdisp");
  ColumnVector kcblood(tis.Nrows());
  kcblood=0.0;

  // Gaussian dispersion arterial curve
  // after Hrabe & Lewis, MRM, 2004 
  float ti=0.0;
  float sqrt2 = sqrt(2);

  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);

      if (casl) kcblood(it) = 2 * exp(-deltblood/T_1b);
      else	kcblood(it) = 2 * exp(-ti/T_1b);

      float erf1 =  (ti-deltblood)/(sqrt2*sig1);
      float erf2 = (ti-deltblood-taub)/(sqrt2*sig2);

      if (erf1>5) erf1=5;
      if (erf2>5) erf2=5;
      if (erf1<-5) erf1=-5;
      if (erf2<-5) erf2 = -5;

      kcblood(it) *= 0.5*( erf(erf1)  - erf(erf2)  );

    }
  return kcblood;
}

ColumnVector kcblood_spatialgaussdisp(const ColumnVector& tis, float deltblood, float taub, float T_1b, float k, bool casl=false) {
    //Tracer_Plus tr("OXASL:kcblood_normdisp");
  ColumnVector kcblood(tis.Nrows());
  kcblood=0.0;

  // Gaussian dispersion arterial curve - in spatial rather than temporal domain
  // after Ozyurt ISMRM 2010 (p4065)
  float ti=0.0;

  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);

      if (casl) kcblood(it) = 2 * exp(-deltblood/T_1b);
      else	kcblood(it) = 2 * exp(-ti/T_1b);

      float erf1 =  (ti-deltblood)/(k*sqrt(ti));
      float erf2 = (ti-deltblood-taub)/(k*sqrt(ti));

      if (erf1>5) erf1=5;
      if (erf2>5) erf2=5;
      if (erf1<-5) erf1=-5;
      if (erf2<-5) erf2 = -5;

      kcblood(it) *= 0.5*( erf(erf1)  - erf(erf2)  );

    }
  return kcblood;
}

ColumnVector kcblood_gallichan(const ColumnVector& tis, float deltblood, float taub, float T_1b, float xdivVm, bool casl=false) {
    //Tracer_Plus tr("OXASL:kcblood_normdisp");
  ColumnVector kcblood(tis.Nrows());
  kcblood=0.0;

  assert(casl==false);

  // Model of dispersion based on a geometrical argument from Gallichan MRM 2008
  // Taking equation [6] (so not including QUIPSSII style saturation)
  // including an 'extra' arrival time term as per the paper
  // bolus duration (taub) takes the place of X/V_m and we let it be a variable

  float ti=0.0;

  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);
      // NOTE: the +xdivVm correction applied to the ti to shift the curve so that the 
      // delay associated with the dispersion parameter has been removed, thus BAT is independent
      // of the dispersion.

      if(ti < deltblood)
	{ 
	  kcblood(it) = 0.0;
	}
      else if((ti >= deltblood) && ti <= (deltblood + taub))
	{ 
	  if (casl) kcblood(it) = 2 * exp(-deltblood/T_1b);
	  else	    kcblood(it) = 2 * exp(-ti/T_1b);

	  kcblood(it) *= 1 - xdivVm/(ti + xdivVm - deltblood); 
	}
      else //(ti > deltblood + taub)
      	{
	  if (casl) kcblood(it) = 2 * exp(-deltblood/T_1b);
	  else	    kcblood(it) = 2 * exp(-ti/T_1b);

      	  kcblood(it) *= taub/(ti + xdivVm - deltblood) ; 
	}
    }
  return kcblood;
}

//Tissue
ColumnVector kctissue_nodisp(const ColumnVector& tis, float delttiss, float tau, float T_1b, float T_1app) {
  //Tracer_Plus tr("OXASL::kctissue_nodisp");
ColumnVector kctissue(tis.Nrows());
 kctissue=0.0;
 float ti=0.0;

  // Tissue kinetic curve no dispersion (pASL)
  // Buxton (1998) model

  float R = 1/T_1app - 1/T_1b; 

  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);
      float F = 2 * exp(-ti/T_1app);

      if(ti < delttiss)
	{ kctissue(it) = 0;}
      else if(ti >= delttiss && ti <= (delttiss + tau))
	{
	  kctissue(it) = F/R * ( (exp(R*ti) - exp(R*delttiss)) ) ;
	}
      else //(ti > delttiss + tau)
	{
	  kctissue(it) = F/R * ( (exp(R*(delttiss+tau)) - exp(R*delttiss))  );
	}
    }
  return kctissue;
}

ColumnVector kctissue_gammadisp(const ColumnVector& tis, float delttiss, float tau, float T_1b, float T_1app, float s, float p) {
  //Tracer_Plus tr("OXASL::kctissue_gammadisp");
  ColumnVector kctissue(tis.Nrows());
  kctissue=0.0;
  float ti=0.0;

  float k=1+p*s;
  float A = T_1app - T_1b;
  float B = A + s*T_1app*T_1b;
  if (B<1e-12) B=1e-12;  //really shouldn't happen, but combination of parameters may arise in artefactual voxels?
  float C = pow(s-1/T_1app+1/T_1b,p*s);
  if (s-1/T_1app+1/T_1b<=0) C=1e-12; //really shouldn't happen, but combination of parameters may arise in artefactual voxels?

  //cout << T_1app << " " << A << " " << B << " "<< C << " " << endl ;

  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);

      if(ti < delttiss)
	{ kctissue(it) = 0;}
      else if(ti >= delttiss && ti <= (delttiss + tau))
	{
	  kctissue(it) = 2* 1/A * exp( -(T_1app*delttiss + (T_1app+T_1b)*ti)/(T_1app*T_1b) )*T_1app*T_1b*pow(B,-k)*
	    (  exp(delttiss/T_1app + ti/T_1b) * pow(s*T_1app*T_1b,k) * ( 1 - igamc(k,B/(T_1app*T_1b)*(ti-delttiss)) ) +
	       exp(delttiss/T_1b + ti/T_1app) * pow(B,k) * ( -1 + igamc(k,s*(ti-delttiss)) )  );  
	  
	}
      else //(ti > delttiss + tau)
	{
	  kctissue(it) = 2* 1/(A*B) *
	    (  exp(-A/(T_1app*T_1b)*(delttiss+tau) - ti/T_1app)*T_1app*T_1b/C*
	       (  pow(s,k)*T_1app*T_1b* ( -1 + exp( (-1/T_1app+1/T_1b)*tau )*( 1 - igamc(k,B/(T_1app*T_1b)*(ti-delttiss)) ) +
					  igamc(k,B/(T_1app*T_1b)*(ti-delttiss-tau)) ) - 
		  exp( -A/(T_1app*T_1b)*(ti-delttiss-tau) )*C*B * ( igamc(k,s*(ti-delttiss-tau)) - igamc(k,s*(ti-delttiss)) ) )  );
	  
	}
      //if (isnan(kctissue(it))) { kctissue(it)=0.0; cout << "Warning NaN in tissue KC"; }
    }
  //cout << kctissue.t() << endl;
  return kctissue;
}

ColumnVector kctissue_gvf(const ColumnVector& tis, float delttiss, float T_1b, float T_1app, float s, float p) {
  //Tracer_Plus tr("OXASL::kctissue_gvf");
  ColumnVector kctissue(tis.Nrows());
  kctissue=0.0;
  float ti=0.0;


  float k=1+p*s;
  float A = T_1app - T_1b;
  float B = A + s*T_1app*T_1b;
  float C = pow(s-1/T_1app+1/T_1b,p*s);
  float sps = pow(s,k);

  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);

      if(ti < delttiss)
	{ kctissue(it) = 0.0;}
      else //if(ti >= delttiss && ti <= (delttiss + tau))
	{
	  kctissue(it) = 2* 1/(B*C) * exp(-(ti-delttiss)/T_1app)*sps*T_1app*T_1b * (1 - igamc(k,(s-1/T_1app-1/T_1b)*(ti-delttiss)));
	}
      // bolus duraiton is specified by the CVF AIF shape and is not an explicit parameter
      //else //(ti > delttiss + tau)
      //	{
      //	  kctissue(it) = exp(-(ti-delttiss-tau)/T_1app) * 2* 1/(B*C) * exp(-(delttiss+tau)/T_1app)*sps*T_1app*T_1b * (1 - igamc(k,(s-1/T_1app-1/T_1b)*(delttiss+tau)));
      //	}
    }
  return kctissue;
}

  ColumnVector kctissue_gaussdisp(const ColumnVector& tis, float delttiss, float tau, float T_1b, float T_1app, float sig1, float sig2) {
    //Tracer_Plus tr("OXASL::kctissue_gaussdisp");
ColumnVector kctissue(tis.Nrows());
 kctissue=0.0;
 float ti=0.0;

  // Tissue kinetic curve gaussian dispersion (pASL)
  // Hrabe & Lewis, MRM, 2004

  float R = 1/T_1app - 1/T_1b; 
  float sqrt2 = sqrt(2);

  for(int it=1; it<=tis.Nrows(); it++)
    {
      ti = tis(it);
      float F = 2 * exp(-ti/T_1app);
      float u1 = (ti-delttiss)/(sqrt2*sig1);
      float u2 = (ti - delttiss - tau)/(sqrt2*sig2);

      kctissue(it) = F/(2*R) * (  (erf(u1) - erf(u2))*exp(R*ti) 
				  - (1 + erf(u1 - (R*sig1)/sqrt2))*exp(R*(delttiss+(R*sig1*sig1)/2)) 
				  + (1 + erf(u2 - (R*sig2)/sqrt2))*exp(R*(delttiss+tau+(R*sig2*sig2)/2)) 
				  );
    
    }
  return kctissue;
}

// --- useful general functions ---
float icgf(float a, float x) {
  //Tracer_Plus tr("OXASL::icgf");

  //incomplete gamma function with a=k, based on the incomplete gamma integral

  return gamma(a)*igamc(a,x);
}

float gvf(float t, float s, float p) {
  //Tracer_Plus tr("OXASL::gvf");

  //The Gamma Variate Function (correctly normalised for area under curve) 
  // Form of Rausch 2000
  // NB this is basically a gamma pdf

  if (t<0)    return 0.0;
  else        return pow(s,1+s*p) / gamma(1+s*p) * pow(t,s*p) * exp(-s*t);
}


} //end namespace
