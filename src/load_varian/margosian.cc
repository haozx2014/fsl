/*
  margosian.c: Half k-space reconstruction (Margosian method)
  
  Copyright Stuart Clare
  FMRIB Centre, University of Oxford.
  
  This program should be considered a beta test version
  and must not be used for any clinical purposes.
  
  Part of...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "margosian.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "image_fftw.h"
#include "hermit.h"
#include "procfunc.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

int hks_margosian(float* data,int ppl,int lpi,int ipv)
{
  float *phase;
  int partial_lpi;
  
  phase=(float *)malloc(ppl*lpi*ipv*2*sizeof(float));

  image_oned_ifftw_r(data,ppl,lpi,ipv);
  hks_centre_kspace(data,ppl,lpi,ipv);
  image_oned_fftw_r(data,ppl,lpi,ipv);
  hks_margosian_phase(data,phase,ppl,lpi,&partial_lpi,ipv);
  hks_margosian_hanning(data,ppl,lpi,partial_lpi,ipv);
  image_oned_fftw_c(data,ppl,lpi,ipv);
  hks_margosian_correct(data,phase,ppl,lpi,ipv);
  image_oned_ifftw_c(data,ppl,lpi,ipv);
  return(0);
}

int hks_margosian_phase(float* data,float* phase,int ppl,int lpi,int *partial_lpi,int ipv)
{
  int i,j,k,centre,extra;
  long l;
  float aux,sum;

  /* Find number of extra lines */
  j=0;
  for(i=0;i<ipv;i++){
    for(j=lpi/2;j<lpi;j++){
      sum=0;
      for(k=0;k<ppl;k++)sum+=data[i*ppl*lpi*2 + j*ppl*2 + k*2];
      if(sum==0)break;
    }
  }
  centre=lpi/2;
  extra = j-centre;
  *partial_lpi=centre+extra;

  for(k=0;k<ipv;k++){
    /* Hanning filter */
    for(i=centre-extra;i<centre+extra;i++){
      aux=0.54-0.46*cos((i-centre+extra)*M_PI/(2*extra));
      for(j=0;j<ppl;j++){
	phase[(k*ppl*lpi+i*ppl+j)*2]=aux*data[(k*ppl*lpi+i*ppl+j)*2];
	phase[(k*ppl*lpi+i*ppl+j)*2+1]=aux*data[(k*ppl*lpi+i*ppl+j)*2+1];
      }
    }
    for(i=0;i<centre-extra;i++){
      for(j=0;j<ppl;j++){
	phase[(k*ppl*lpi+i*ppl+j)*2]=0;
	phase[(k*ppl*lpi+i*ppl+j)*2+1]=0;
      }
    }
    for(i=centre+extra;i<lpi;i++){
      for(j=0;j<ppl;j++){
	phase[(k*ppl*lpi+i*ppl+j)*2]=0;
	phase[(k*ppl*lpi+i*ppl+j)*2+1]=0;
      }
    }
  }
  
  image_oned_fftw_c(phase,ppl,lpi,ipv);

  for(l=0;l<ppl*lpi*ipv;l++){
    phase[l]=atan2(phase[l*2+1],phase[l*2]);
  }

  return(0);
}

int hks_margosian_hanning(float* data,int ppl,int lpi,int partial_lpi,int ipv)
{

  int i,j,k,centre,extra;
  float aux;

  centre = lpi/2;
  extra = partial_lpi - centre;

  for(k=0;k<ipv;k++){
    /* Hanning filter */
    for(i=centre-extra;i<centre+extra;i++){
      aux=0.54-0.46*cos((i-centre+extra)*M_PI/(2*extra));
      for(j=0;j<ppl;j++){
	data[(k*ppl*lpi+i*ppl+j)*2]=aux*data[(k*ppl*lpi+i*ppl+j)*2];
	data[(k*ppl*lpi+i*ppl+j)*2+1]=aux*data[(k*ppl*lpi+i*ppl+j)*2+1];
      }
    }
  }
  
  return(0);
}

int hks_margosian_correct(float* data,float* phase,int ppl,int lpi,int ipv)
{
  int i,j,k;
  float re,im,ph,cs,sn;

   for(i=0;i<ipv;i++){
     for(j=0;j<lpi;j++){
       for(k=0;k<ppl;k++){
	 re = data[i*ppl*lpi*2 + j*ppl*2 + k*2];
	 im = data[i*ppl*lpi*2 + j*ppl*2 + k*2 + 1];
	 ph = phase[i*ppl*lpi + j*ppl + k];
	 cs=cos(ph);
	 sn=sin(ph);
	 data[i*ppl*lpi*2 + j*ppl*2 + k*2] = (re*cs + im*sn);
	 data[i*ppl*lpi*2 + j*ppl*2 + k*2 + 1] = (im*cs - re*sn);
       }
     }
   }
   return(0);
}
