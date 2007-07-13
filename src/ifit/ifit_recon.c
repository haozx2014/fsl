/*
  ifit_recon: Reconstruction of data set from output images
  Part of the iFit suite of programs for image fitting to functions

  Stuart Clare, FMRIB Physics Group

  Copyright (C) 1999-2000 University of Oxford.
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fslio/fslio.h"

#define LINEAR 1
#define SATREC 2
#define EXPDEC 3
#define INVREC 4
#define INVRECM 5
#define INVRECI 6
#define INVRECIM 7
#define NOCLIP -1.0

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

#define MAXLEN 128
typedef unsigned char flag;

int main(int argc,char* argv[])
{
  char infile[MAXLEN],outfile[MAXLEN],tmpfile[MAXLEN];
  short ppl,lpi,ipv,vols,dt;
  int func,num,mdim,err;
  flag dtshort,*fimage;
  float *Adata,*Bdata,*Cdata;
  long i,j,volsize;
  short *simage;
  FSLIO *ip,*op;

  void usage_error();
  float linear(float x,float A,float B);
  float expdec(float x,float A,float B);
  float satrec(float x,float A,float B);
  float invrec(float x,float A,float B);
  float invrecm(float x,float A,float B);
  float invreci(float x,float A,float B,float C);
  float invrecim(float x,float A,float B,float C);

  dtshort=0;
  fimage=simage=Cdata=NULL;

  /* Read in the arguments */
  if(argc==1)usage_error();

  func=num=err=0;
  infile[0]=outfile[0]='\0';

  for(i=1;i<argc;i++){
    if((!strcmp(argv[i],"-h"))||(!strcmp(argv[i],"-help")))
      usage_error();
    if(!strcmp(argv[i],"-in")){sprintf(infile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-out")){sprintf(outfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-func")){sscanf(argv[i+1],"%d",&func);continue;}
    if(!strcmp(argv[i],"-err")){err=1;continue;}
    if(!strcmp(argv[i],"-x")){
      sscanf(argv[i+1],"%d",&num);
      if(argc<i+num+2){
	printf("Not enough x values specified\n");
	usage_error();
      }
      if((xvals=(float *)malloc(num*sizeof(float)))==NULL){
	printf("Malloc failed.\n");
	return(1);
      }
      for(j=0;j<num;j++){
	sscanf(argv[i+j+2],"%f",&xvals[j]);
      }
      continue;
    }
    if(argv[i][0]=='-'){
      printf("Command flag %s not recognised\\n",argv[i]);
      usage_error();
    }
  }

  switch(func){
  case LINEAR:
  case SATREC:
  case EXPDEC:
  case INVREC:
  case INVRECM:
    mdim=2;
    break;
  case INVRECI:
  case INVRECIM:
    mdim=3;
    break;
  default:
    printf("Function not specified or not valid for this routine.\n");
    printf("Use -func flag\n");
    return(1);
  }
  
  if(num==0){
    printf("No x values specified. Use -x flag\n");
    usage_error();
  }

  if(infile[0]=='\0'){
    printf("No input file specified. Use -in flag\n");
    usage_error();
  }
  if(outfile[0]=='\0'){
    printf("No output file specified. Use -out flag\n");
    usage_error();
  }

  /* Read in coefficients */
  sprintf(tmpfile,"%s_A",infile);
  if((ip=FslOpen(tmpfile,"r"))==NULL){
    printf("Failed to open AVW file: %s\n",tmpfile);
    return(1);
  }
  FslGetDim(ip,&ppl,&lpi,&ipv,&vols);
  volsize=ppl*lpi*ipv;
  if((Adata=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  FslReadVolumes(ip,Adata,1);
  FslClose(ip);

  sprintf(tmpfile,"%s_B",infile);
  if((ip=FslOpen(tmpfile,"r"))==NULL){
    printf("Failed to open AVW file: %s\n",tmpfile);
    return(1);
  }
  if((Bdata=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  FslReadVolumes(ip,Bdata,1);

  /* Open output file */
  if((op=FslOpen(outfile,"w"))==NULL){
    printf("Can't open : %s\n",outfile);
    return(1);  
  }
  FslCloneHeader(op,ip);
  FslSetDataType(op,DT_FLOAT);
  FslSetDim(op,ppl,lpi,ipv,num);
  FslClose(ip);

  if(mdim==3){
    sprintf(tmpfile,"%s_C",infile);
    if((ip=FslOpen(tmpfile,"r"))==NULL){
      printf("Failed to open AVW file: %s\n",tmpfile);
      return(1);
    }
    if((Cdata=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
    FslReadVolumes(ip,Cdata,1);        
    FslClose(ip);
  }

  if((image=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  if(err){
    if((ip=FslOpen(infile,"r"))==NULL){
      printf("Failed to open AVW file: %s\n",infile);
      return(1);
    }
    if((fimage=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
    FslGetDataType(ip,&dt);
    if((dt=DT_SIGNED_SHORT)){
      dtshort=1;
      if((simage=(short *)malloc(volsize*sizeof(short)))==NULL){
	printf("Malloc failed.\n");
	return(1);
      }
    }
    else dtshort=0;
  }

  for(i=0;i<num;i++){
    
    if(err){
      if(dtshort){
	if(FslReadVolumes(ip,simage,1)!=1){
	  printf("Read error\n");
	  return(1);
	}
      }
      else {
	if(FslReadVolumes(ip,fimage,1)!=1){
	  printf("Read error\n");
	  return(1);
	}
      }
    }

    for(j=0;j<volsize;j++){
      if((Adata[j]!=0)&&(Bdata[j]!=0)){
	switch(func){
	case LINEAR:
	  image[j]=linear(xvals[i],Adata[j],Bdata[j]);
	  break;
	case SATREC:
	  image[j]=satrec(xvals[i],Adata[j],Bdata[j]);	
	  break;
	case EXPDEC:
	  image[j]=expdec(xvals[i],Adata[j],Bdata[j]);	
	  break;
	case INVRECM:
	  image[j]=invrecm(xvals[i],Adata[j],Bdata[j]);	
	  break;
	case INVRECIM:
	  image[j]=invrecim(xvals[i],Adata[j],Bdata[j],Cdata[j]);	
	  break;
	case INVREC:
	  image[j]=invrec(xvals[i],Adata[j],Bdata[j]);	
	  break;
	case INVRECI:
	  image[j]=invreci(xvals[i],Adata[j],Bdata[j],Cdata[j]);	
	  break;
	}
	if(err){
	  if(dtshort)image[j]=(float)simage[j]-image[j];
	  else image[j]=fimage[j]-image[j];
	}
      }
      else image[j]=0;
    }

    FslWriteHeader(op);
    FslWriteVolumes(op,image,1);
  }

  if(err)FslClose(ip);
  FslClose(op);
  return(0);
}
float linear(float x,float A,float B)
{
  return(B*x+A);
}
float expdec(float x,float A,float B)
{
  return(A*exp(-x/B));
}
float satrec(float x,float A,float B)
{
  return(A*(1-exp(-x/B)));
}
float invrec(float x,float A,float B)
{
  return(A*(1-2*exp(-x/B)));
}
float invrecm(float x,float A,float B)
{
  return(fabs(A*(1-2*exp(-x/B))));
}
float invreci(float x,float A,float B,float C)
{
  C=cos(C*M_PI/180);
  return(A*(1-(1-C)*exp(-x/B)));
}
float invrecim(float x,float A,float B,float C)
{
  C=cos(C*M_PI/180);
  return(fabs(A*(1-(1-C)*exp(-x/B))));
  return(B*x+A);
}

void usage_error()
{
  printf("usage: ifit_recon [-options]\n");
  printf("\t -in basefile\n");
  printf("\t -func function_number\n");
  printf("\t -x num x1 x2 x3 ...\n");
  printf("\t [-out outfile]\n");
  printf("\t [-err]\n");
  printf("\t -help\n");
  printf("\n");
  printf("functions:\n");
  printf("\t 1 - Linear  Y = A.X + B\n");
  printf("\t 2 - Saturation Recovery  Y = A.(1 - exp(-X/B))\n");
  printf("\t 3 - Exponential Decay  Y = A.exp(-X/B)\n");
  printf("\t 4 - Inversion Recovery  Y = A.(1 - 2exp(-X/B))\n");
  printf("\t 5 - Modulus Inv Rec  Y = A.|1 - 2exp(-X/B)|\n");
  printf("\t 6 - Incomplete Inv Rec  Y = A.(1 - (1-cosC).exp(-X/B))\n");
  printf("\t 7 - Modulus Incomplete IR : Y = A.|1 - (1-cosC).exp(-X/B)|\n\n");
  printf("\n");
  exit(1);
}
