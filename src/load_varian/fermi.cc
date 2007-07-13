/*
  fermi.c: Fermi filter

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.
  Based on filter in fft3d Ravi S. Mennon

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "fermi.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>


#define F_FRAC 0.9

void fermi_filter_3d(float* data,int ppl,int lpi,int ipv)
{
  int i,j,k;
  float cutoff,filter,r1,r2,r3;
  float *slice;
  float **fermi;

  if((slice=(float *)calloc(ipv,sizeof(float)))==NULL){
    printf("Malloc failed (slice)");
    return;
  }
  if((fermi=(float **)malloc(lpi*sizeof(float*)))==NULL){
    printf("Malloc failed (fermi)");
    return;
  }
  for(i=0;i<lpi;i++){
    if((fermi[i]=(float *)calloc(ppl,sizeof(float)))==NULL){
      printf("Malloc failed (fermi)");
      return;
    }
  }

  cutoff = 0.5*F_FRAC;
  for (i = 1; i < (1+ipv/2); i++) {
    r1 = (double)i/(double)(ipv/2);
    filter =  1.0/(1.0+exp((r1-cutoff)*3.0));
    slice[(ipv/2)-i] = filter;
    slice[((ipv/2)-1)+i] = filter;
  }

  cutoff = (double)F_FRAC;
  for (j = 1; j < (1+lpi/2); j++) {
    for (k = 1; k < (1+ppl/2); k++) {
      r2 = (double)j/(double)(lpi/2);
      r3 = (double)k/(double)(ppl/2);
      r1 = r2*r2 + r3*r3;
      filter =  1.0/(1.0+exp((r1-cutoff)*5.0));
      fermi[(lpi/2)-j][((ppl/2)-1)+k] = filter;
      fermi[((lpi/2)-1)+j][(ppl/2)-k] = filter;
      fermi[(lpi/2)-j][(ppl/2)-k] = filter;
      fermi[((lpi/2)-1)+j][((ppl/2)-1)+k] = filter;
    }
  }

  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl;k++){
	data[i*ppl*lpi*2 + j*ppl*2 + k*2]*=(fermi[j][k]*slice[i]);
	data[i*ppl*lpi*2 + j*ppl*2 + k*2 + 1]*=(fermi[j][k]*slice[i]);	
      }
    }
  }

  for(i=0;i<lpi;i++)free(fermi[i]);
  free(fermi);
  free(slice);
}

void fermi_filter_2d(float* data,int ppl,int lpi,int ipv)
{
  int i,j,k;
  float cutoff,filter,r1,r2,r3;
  float **fermi;

  if((fermi=(float **)malloc(lpi*sizeof(float*)))==NULL){
    printf("Malloc failed (fermi)");
    return;
  }
  for(i=0;i<lpi;i++){
    if((fermi[i]=(float *)calloc(ppl,sizeof(float)))==NULL){
      printf("Malloc failed (fermi)");
      return;
    }
  }

  cutoff = (double)F_FRAC;
  for (j = 1; j < (1+lpi/2); j++) {
    for (k = 1; k < (1+ppl/2); k++) {
      r2 = (double)j/(double)(lpi/2);
      r3 = (double)k/(double)(ppl/2);
      r1 = r2*r2 + r3*r3;
      filter =  1.0/(1.0+exp((r1-cutoff)*5.0));
      fermi[(lpi/2)-j][((ppl/2)-1)+k] = filter;
      fermi[((lpi/2)-1)+j][(ppl/2)-k] = filter;
      fermi[(lpi/2)-j][(ppl/2)-k] = filter;
      fermi[((lpi/2)-1)+j][((ppl/2)-1)+k] = filter;
    }
  }

  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl;k++){
	data[i*ppl*lpi*2 + j*ppl*2 + k*2]*=(fermi[j][k]);
	data[i*ppl*lpi*2 + j*ppl*2 + k*2 + 1]*=(fermi[j][k]);	
      }
    }
  }

  for(i=0;i<lpi;i++)free(fermi[i]);
  free(fermi);
}
