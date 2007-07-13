/*
  segment.c: Trims down the matrix size

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "segment.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>

int fsegment(float* data,int lowx,int numx,int totx,
	      int lowy,int numy,int toty,
	      int lowz,int numz,int totz,int cplx)
{
  float* tmp;
  long volsize,in,out;
  int i,j,k;

  if((lowx+numx>totx)||(lowy+numy>toty)||(lowz+numz>totz)){
    printf("Error in dimensions to segment.");
    return(1);    
  }

  volsize = numx*numy*numz*cplx;
  if((tmp=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed.");
    return(1);
  }
  
  for(k=0;k<numz;k++){
    for(j=0;j<numy;j++){
      for(i=0;i<numx;i++){
	in = (k+lowz)*toty*totx + (j+lowy)*totx + i+lowx;
	out = k*numy*numx + j*numx + i;
	if(cplx>1){
	  tmp[out*2]=data[in*2];
	  tmp[out*2+1]=data[in*2+1];
	}
	else tmp[out]=data[in];
      }
    }
  }
  memcpy(data,tmp,volsize*sizeof(float));
  return(0);
}
int ssegment(short* data,int lowx,int numx,int totx,
	      int lowy,int numy,int toty,
	      int lowz,int numz,int totz)
{
  short* tmp;
  long volsize,in,out;
  int i,j,k;

  if((lowx+numx>totx)||(lowy+numy>toty)||(lowz+numz>totz)){
    printf("Error in dimensions to segment.");
    return(1);    
  }

  volsize = numx*numy*numz;
  if((tmp=(short *)malloc(volsize*sizeof(short)))==NULL){
    printf("Malloc failed.");
    return(1);
  }
  
  for(k=0;k<numz;k++){
    for(j=0;j<numy;j++){
      for(i=0;i<numx;i++){
	in = (k+lowz)*toty*totx + (j+lowy)*totx + i+lowx;
	out = k*numy*numx + j*numx + i;
	tmp[out]=data[in];
      }
    }
  }
  memcpy(data,tmp,volsize*sizeof(short));
  return(0);
}
