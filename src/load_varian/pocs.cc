/*
  hksfunc.c: Half k-space reconstruction (POCS method)
  
  Copyright Stuart Clare
  FMRIB Centre, University of Oxford.
  
  This program should be considered a beta test version
  and must not be used for any clinical purposes.
  
  Part of...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "pocs.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "image_fftw.h"
#include "hermit.h"
#include "procfunc.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

int hks_pocs(float* data,int ppl,int lpi,int ipv)
{
  int i,iter;
  int partial_lpi;
  float *phase,*pseudo;

  iter = 4;

  phase=(float *)malloc(ppl*lpi*ipv*2*sizeof(float));
  pseudo=(float *)malloc(ppl*lpi*ipv*2*sizeof(float));

  image_oned_ifftw_r(data,ppl,lpi,ipv);
  hks_centre_kspace(data,ppl,lpi,ipv);
  image_oned_fftw_r(data,ppl,lpi,ipv);
  hks_pocs_phase(data,phase,ppl,lpi,&partial_lpi,ipv);
  
  for(i=0;i<iter;i++){
    hks_pocs_imgen(data,phase,pseudo,ppl,lpi,ipv);
    hks_pocs_merge(data,pseudo,ppl,lpi,partial_lpi,ipv);
  }

  free(phase);
  free(pseudo);

  return(0);
}

int hks_pocs_phase(float* data,float* phase,int ppl,int lpi,int *partial_lpi,int ipv)
{
  int i,j,k,centre,extra;
  long l;
  float aux,sum;

  j=0;
  /* Find number of extra lines */
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

int hks_pocs_imgen(float* data,float* phase,float* pseudo,int ppl,int lpi,int ipv)
{
  float m,p;
  long l;

  memcpy(pseudo,data,ppl*lpi*ipv*2*sizeof(float));

  image_oned_fftw_c(pseudo,ppl,lpi,ipv);

  for(l=0;l<ppl*lpi*ipv;l++){
    m=magnitude(pseudo[l*2+1],pseudo[l*2]);
    p=phase[l];
    pseudo[l*2]=m*cos(p);
    pseudo[l*2+1]=m*sin(p);
  }

  image_oned_ifftw_c(pseudo,ppl,lpi,ipv);

  return(0);
}

int hks_pocs_merge(float* data,float* pseudo,int ppl,int lpi,int partial_lpi,int ipv)
{
  int i,j,k,centre,extra,blend,b;
  float fblend;

  centre = lpi/2;
  extra = partial_lpi - centre;
  blend = extra/2;

  for(i=0;i<ipv;i++){
    for(k=0;k<ppl;k++){
      for(j=partial_lpi-blend,b=1;j<partial_lpi;j++,b++){
	fblend = b/(float)blend+1;
	data[(i*ppl*lpi+j*ppl+k)*2] =
	  (1.0-fblend)*data[(i*ppl*lpi+j*ppl+k)*2] + fblend*pseudo[(i*ppl*lpi+j*ppl+k)*2];
	data[(i*ppl*lpi+j*ppl+k)*2+1] =
	  (1.0-fblend)*data[(i*ppl*lpi+j*ppl+k)*2+1] + fblend*pseudo[(i*ppl*lpi+j*ppl+k)*2+1];
      }
      for(j=partial_lpi;j<lpi;j++){
	data[(i*ppl*lpi+j*ppl+k)*2]=pseudo[(i*ppl*lpi+j*ppl+k)*2];
	data[(i*ppl*lpi+j*ppl+k)*2+1]=pseudo[(i*ppl*lpi+j*ppl+k)*2+1];
      }
    }
  }
  return(0);
}
