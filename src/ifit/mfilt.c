/*
  mfilt

  Stuart Clare, FMRIB Physics Group

  Copyright (c) 2000 University of Oxford
*/

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 3.3 (c) 2006, The University of
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

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"fslio/fslio.h"

#define MAXLEN 128

int main(int argc,char *argv[])
{
  short ppl,lpi,ipv,vols,dt,window;
  long volsize,i,j,k,v,size,count,s,t,u;
  float *fimage,*out,*array;
  int bounds(short i,short j,short k,short si,short sj,short sk);
  float calc_med(float *array,int size);

  FSLIO *ifp,*ofp;

  count = 0;

  if(argc<3){
    printf("usage: mfilt infile outfile [window]\n");
    return(1);
  }
  if(argc==4)sscanf(argv[3],"%hd",&window);
  else window=2;
  window/=2;

  array=(float *)malloc((window*2+1)*(window*2+1)*(window*2+1)*sizeof(float));
 
  if((ifp=FslOpen(argv[1],"r"))==NULL){
    printf("Failed to open AVW file: %s\n",argv[1]);
    return(1);
  }
  if((ofp=FslOpen(argv[2],"w"))==NULL){
      printf("Failed to open AVW file: %s\n",argv[2]);
      return(1);
  }
  FslCloneHeader(ofp,ifp);
  FslGetDim(ifp,&ppl,&lpi,&ipv,&vols);
  volsize = ppl*lpi*ipv;

  FslGetDataType(ifp,&dt);
  if(dt!=DT_FLOAT){
    printf("Data type not supported\n");
    return(1);
  }

  /* Malloc away */
  if((fimage=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed\n");
    return(1);
  }
  if((out=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed\n");
    return(1);
  }

  for(v=0;v<vols;v++){

    if(FslReadVolumes(ifp,fimage,1)!=1){
      printf("Read error.\n");
      return(1);
    }
    
    for(i=0;i<ipv;i++){
      for(j=0;j<lpi;j++){
	for(k=0;k<ppl;k++){
	  
	  if((fimage[i*ppl*lpi+j*lpi+k]>1)&&(fimage[i*ppl*lpi+j*lpi+k]<10000)){
	    out[i*ppl*lpi+j*lpi+k]=fimage[i*ppl*lpi+j*lpi+k];
	  }
	  else {
	    count++;
	    size=0;
	    for(s=-window;s<=window;s++){
	      for(t=-window;t<=window;t++){
		for(u=-window;u<=window;u++){
		  if(bounds(i+s,j+t,k+u,ipv,lpi,ppl)){
		    array[size++]=fimage[(i+s)*ppl*lpi+(j+t)*lpi+k+u];
		  }
		}
	      }
	    }
	    out[i*ppl*lpi+j*lpi+k]=calc_med(array,size);
	  }
	}
      }
    }
    FslWriteHeader(ofp);
    FslWriteVolumes(ofp,out,1);
  }

  FslClose(ifp);
  FslClose(ofp);
  free(fimage);
  free(out);
  
  printf("Corrected %ld voxels\n",count);
  
  return(0);
}
int bounds(short i,short j,short k,short si,short sj,short sk)
{
  if((i>=si)||(i<0))return(0);
  if((j>=sj)||(j<0))return(0);
  if((k>=sk)||(k<0))return(0);
  return(1);
}

float calc_med(float *data,int size){
  int l;
  float *tmp,med;
  float nrselect(unsigned long k, unsigned long n, float arr[]);
  
  tmp=(float *)malloc(size*sizeof(float));
  for(l=0;l<size;l++)tmp[l]=data[l];
  med=nrselect(size/2,size,tmp-1);
  free(tmp);
  return(med);
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

float nrselect(unsigned long k, unsigned long n, float arr[])
{
  unsigned long i,ir,j,l,mid;
  float a,temp;

  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir])
          }
      return arr[k];
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
        if (arr[l] > arr[ir]) {
          SWAP(arr[l],arr[ir])
            }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir])
          }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1])
          }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j])
          }
      arr[l+1]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}
#undef SWAP

