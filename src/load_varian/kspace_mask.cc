/*
  kspace_mask.c: Masks out part of k-space

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "kspace_mask.h"

#include <cstdlib>
#include <cstdio>

void kspace_mask_border(float* data,int ppl,int lpi,int ipv,int border)
{
  int i,j,k;
  
  for(i=0;i<ipv;i++){
    for(j=0;j<border;j++){
      for(k=0;k<ppl;k++){
	data[i*ppl*lpi*2+j*ppl*2+k*2]=0;
	data[i*ppl*lpi*2+j*ppl*2+k*2+1]=0;
      }
    }
    for(j=lpi-border;j<lpi;j++){
      for(k=0;k<ppl;k++){
	data[i*ppl*lpi*2+j*ppl*2+k*2]=0;
	data[i*ppl*lpi*2+j*ppl*2+k*2+1]=0;
      }
    }
    for(k=0;k<border;k++){
      for(j=0;j<lpi;j++){
	data[i*ppl*lpi*2+j*ppl*2+k*2]=0;
	data[i*ppl*lpi*2+j*ppl*2+k*2+1]=0;
      }
    }
    for(k=ppl-border;k<ppl;k++){
      for(j=0;j<lpi;j++){
	data[i*ppl*lpi*2+j*ppl*2+k*2]=0;
	data[i*ppl*lpi*2+j*ppl*2+k*2+1]=0;
      }
    }
  }
}
