/*
  reslice.c: Reslices data to the axial plane

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "reslice.h"
#include <cstdio>
#include <cstdlib>

int sag2ax(float *data,int ppl,int lpi,int ipv)
{
  int i,j,k;
  long l,in,out,volsize;
  float *tmp;

  volsize = ppl*lpi*ipv*2;
  tmp=(float *)malloc(volsize*sizeof(float));

  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl;k++){
	
	in = i*ppl*lpi + j*ppl + k;
	out = (ipv-1-i) + (ppl-1-k)*ipv + j*ipv*ppl;

	if(out*2>volsize){
	  printf("Reslice error\n");
	  return(1);
	}
	tmp[out*2] = data[in*2];
	tmp[out*2+1] = data[in*2+1];
      }
    }
  }
  for(l=0;l<volsize;l++)data[l]=tmp[l];

  free(tmp);
  return(0);
}
int cor2ax(float *data,int ppl,int lpi,int ipv)
{
  int i,j,k;
  long l,in,out,volsize;
  float *tmp;

  volsize = ppl*lpi*ipv*2;
  tmp=(float *)malloc(volsize*sizeof(float));

  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl;k++){
	
	in = i*ppl*lpi + j*ppl + k;
	out = j*ppl*ipv + (ipv-1-i)*ppl + k;

	if(out*2>volsize){
	  printf("Reslice error\n");
	  return(1);
	}

	tmp[out*2] = data[in*2];
	tmp[out*2+1] = data[in*2+1];
      }
    }
  }
  for(l=0;l<volsize;l++)data[l]=tmp[l];

  free(tmp);
  return(0);
}
