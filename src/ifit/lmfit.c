/*
  lmfit: Image fitting program based on Levenberg-Marquardt Algorithm
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
#include <string.h>
#include "fslio/fslio.h"
#include "vector.h"
#include "matrix.h"
#include "nr.h"

#define LINEAR 1
#define SATREC 2
#define EXPDEC 3
#define INVREC 4
#define INVRECM 5
#define INVRECI 6
#define INVRECIM 7
#define B1MAP 8

#define NOCLIP -1.0

#define MAXLEN 128
#define MAXIT 1000
typedef unsigned char flag;

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

/* Global variable keeps things simple */
float global_flip;

int main(int argc,char* argv[])
{
  flag dtshort,masked,aim,bim,cim,err,wfit;
  char infile[MAXLEN],maskfile[MAXLEN];
  char outfile[MAXLEN],tmpfile[MAXLEN];
  char aimfile[MAXLEN],bimfile[MAXLEN],cimfile[MAXLEN];
  char *mask;
  /*char orient,morient;*/
  short ppl,lpi,ipv,vols,dt;
  short mppl,mlpi,mipv,mvols,mdt;
  int func,num,nc,ng,i,errcnt;
  short *image;
  long j,volsize;
  float *fimage,*xvals,*yvals,*a;
  float **data,*max,*min;
  float *aimage,*cimage,*bimage;
  float uA,lA,uB,lB,A,B,slope;
  FSLIO *ip,*mp,*op;

  void usage_error();
  int lmfit(float *x,float *y,int n,float *a,int ma,int function);

  /* Read in the arguments */
  if(argc==1)usage_error();

  func=num=nc=err=wfit=0;
  infile[0]=outfile[0]=maskfile[0]='\0';
  aimfile[0]=bimfile[0]=cimfile[0]='\0';
  A=B=0;
  uA=uB=lA=lB=NOCLIP;
  aimage=bimage=cimage=NULL;
  yvals=xvals=NULL;
  mask=NULL;

  for(i=1;i<argc;i++){
    if((!strcmp(argv[i],"-h"))||(!strcmp(argv[i],"-help")))
      usage_error();
    if(!strcmp(argv[i],"-in")){sprintf(infile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-out")){sprintf(outfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-func")){sscanf(argv[i+1],"%d",&func);continue;}
    if(!strcmp(argv[i],"-mask")){sprintf(maskfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-nc")){sscanf(argv[i+1],"%d",&nc);continue;}
    if(!strcmp(argv[i],"-Amax")){sscanf(argv[i+1],"%f\n",&uA);continue;}
    if(!strcmp(argv[i],"-Bmax")){sscanf(argv[i+1],"%f\n",&uB);continue;}
    if(!strcmp(argv[i],"-Amin")){sscanf(argv[i+1],"%f\n",&lA);continue;}
    if(!strcmp(argv[i],"-Bmin")){sscanf(argv[i+1],"%f\n",&lB);continue;}
    if(!strcmp(argv[i],"-Aim")){sprintf(aimfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-Bim")){sprintf(bimfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-Cim")){sprintf(cimfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-A")){sscanf(argv[i+1],"%f\n",&A);continue;}
    if(!strcmp(argv[i],"-B")){sscanf(argv[i+1],"%f\n",&B);continue;}
    if(!strcmp(argv[i],"-err")){err=1;continue;}
    if(!strcmp(argv[i],"-w")){wfit=1;continue;}
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
      if((yvals=(float *)malloc(num*sizeof(float)))==NULL){
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

  if(aimfile[0]=='\0')aim=0;else aim=1;
  if(bimfile[0]=='\0')bim=0;else bim=1;
  if(cimfile[0]=='\0')cim=0;else cim=1;
  if(maskfile[0]=='\0')masked=0;else masked=1;

  switch(func){
  case SATREC:
  case EXPDEC:
  case INVREC:
  case B1MAP:
    break;
  case INVRECI:
    if(!cim){
      printf("Incomplete inversion recovery fit only allowed with Cimage\n");
      return(1);
    }
    else break;
  default:
    printf("Function not specified or not valid for this routine.\n");
    printf("Use -func flag\n");
    return(1);
  }

  if(num==0){
    printf("No x values specified. Use -x flag\n");
    usage_error();
  }
  if(nc==0){
    nc=1;
    ng=2;
  } else {
    ng=nc*2;
  }

  if(infile[0]=='\0'){
    printf("No input file specified. Use -in flag\n");
    usage_error();
  }
  if(outfile[0]=='\0')sprintf(outfile,"%s",infile);

  if(err){
    printf("lmfit cannot output error images yet!\n");
    return(1);
  }

  /* Open input file */
  if((ip=FslOpen(infile,"r"))==NULL){
    printf("Failed to open AVW file: %s\n",infile);
    return(1);
  }
  FslGetDim(ip,&ppl,&lpi,&ipv,&vols);
  FslGetDataType(ip,&dt);
  if(dt==DT_SIGNED_SHORT){
    dtshort=1;
  }
  else {
    dtshort=0;
    if(dt!=DT_FLOAT){
      printf("Data type not short ot float\n");
      return(1);
    }
  }
  volsize=ppl*lpi*ipv;

  if(vols<num){
    printf("Not enough volumes in AVW file\n");
    return(1);
  }

  /* Malloc */
  if((image=(short *)malloc(volsize*num*sizeof(short)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  if((fimage=(float *)malloc(volsize*num*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  if((data=(float **)malloc(ng*sizeof(float*)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  for(i=0;i<ng;i++){
    if((data[i]=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
  }
  if((max=(float *)malloc(ng*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  if((min=(float *)malloc(ng*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  if((a=(float *)malloc(ng*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  

  /* Read in mask */
  if(masked){
    if((mask=(char *)malloc(volsize))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
    if((mp=FslOpen(maskfile,"r"))==NULL){
      printf("Failed to open AVW file: %s\n",maskfile);
      return(1);
    }
    FslGetDim(mp,&mppl,&mlpi,&mipv,&mvols);
    /*FslGetOrientation(mp,&morient);
      FslGetOrientation(ip,&orient);*/
    if((ppl!=mppl)||(lpi!=mlpi)||(mipv!=ipv)){
      printf("Mask not compatible with input file\n");
      return(1);
    }
    FslGetDataType(mp,&mdt);
    switch(mdt){
    case DT_UNSIGNED_CHAR:
      if(FslReadVolumes(mp,mask,1)!=1){
	printf("Read error: %s\n",maskfile);
	return(1);
      }
      break;
    case DT_SIGNED_SHORT:
      if(FslReadVolumes(mp,image,1)!=1){
	printf("Read error: %s\n",maskfile);
	return(1);
      }
      for(j=0;j<volsize;j++){
	if(image[j]>0)mask[j]=1;
	else mask[j]=0;
      }
      break;
    case DT_FLOAT:
      if(FslReadVolumes(mp,fimage,1)!=1){
	printf("Read error: %s\n",maskfile);
	return(1);
      }
      for(j=0;j<volsize;j++){
	if(fimage[j]>0)mask[j]=1;
	else mask[j]=0;
      }
      break;
    default:
      printf("Data type not byte or short\n");
      return(1);
    }
    FslClose(mp);
  }

  /* Read in A image */
  if(aim){
    if((mp=FslOpen(aimfile,"r"))==NULL){
      printf("Failed to open AVW file: %s\n",aimfile);
      return(1);
    }
    FslGetDim(mp,&mppl,&mlpi,&mipv,&mvols);
    /*FslGetOrientation(mp,&morient);*/
    if((ppl!=mppl)||(lpi!=mlpi)||(mipv!=ipv)){
      printf("Aimage not compatible with input file\n");
      return(1);
    }
    FslGetDataType(mp,&mdt);
    if((mdt!=DT_SIGNED_SHORT)&&(mdt!=DT_FLOAT)){
      printf("Aimage data type not valid\n");
      return(1);
    }
    if((aimage=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
    if(mdt==DT_SIGNED_SHORT){
      if(FslReadVolumes(mp,image,1)!=1){
	printf("Read error: %s\n",aimfile);
	return(1);
      }
      for(j=0;j<volsize;j++)aimage[j]=(float)image[j];
    }
    else{
      if(FslReadVolumes(mp,aimage,1)!=1){
	printf("Read error: %s\n",aimfile);
	return(1);
      }
    }
    FslClose(mp);
  }

  /* Read in B image */
  if(bim){
    if((mp=FslOpen(bimfile,"r"))==NULL){
      printf("Failed to open AVW file: %s\n",bimfile);
      return(1);
    }
    FslGetDim(mp,&mppl,&mlpi,&mipv,&mvols);
    /*FslGetOrientation(mp,&morient);*/
    if((ppl!=mppl)||(lpi!=mlpi)||(mipv!=ipv)){
      printf("Bimage not compatible with input file\n");
      return(1);
    }
    FslGetDataType(mp,&mdt);
    if((mdt!=DT_SIGNED_SHORT)&&(mdt!=DT_FLOAT)){
      printf("Bimage data type not valid\n");
      return(1);
    }
    if((bimage=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
    if(mdt==DT_SIGNED_SHORT){
      if(FslReadVolumes(mp,image,1)!=1){
	printf("Read error: %s\n",bimfile);
	return(1);
      }
      for(j=0;j<volsize;j++)bimage[j]=(float)image[j];
    }
    else{
      if(FslReadVolumes(mp,bimage,1)!=1){
	printf("Read error: %s\n",bimfile);
	return(1);
      }
    }
    FslClose(mp);
  }

  /* Read in C image */
  if(cim){
    if((mp=FslOpen(cimfile,"r"))==NULL){
      printf("Failed to open AVW file: %s\n",cimfile);
      return(1);
    }
    FslGetDim(mp,&mppl,&mlpi,&mipv,&mvols);
    /*FslGetOrientation(mp,&morient);*/
    if((ppl!=mppl)||(lpi!=mlpi)||(mipv!=ipv)){
      printf("Cimage not compatible with input file\n");
      return(1);
    }
    FslGetDataType(mp,&mdt);
    if((mdt!=DT_SIGNED_SHORT)&&(mdt!=DT_FLOAT)){
      printf("Cimage data type not valid\n");
      return(1);
    }
    if((cimage=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
    if(mdt==DT_SIGNED_SHORT){
      if(FslReadVolumes(mp,image,1)!=1){
	printf("Read error: %s\n",cimfile);
	return(1);
      }
      for(j=0;j<volsize;j++)cimage[j]=(float)image[j];
    }
    else{
      if(FslReadVolumes(mp,cimage,1)!=1){
	printf("Read error: %s\n",cimfile);
	return(1);
      }
    }
    FslClose(mp);
  }

  /* Read in the data */
  if(dtshort){
    if(FslReadVolumes(ip,image,num)!=num){
      printf("Read error: %s\n",infile);
      return(1);
    }
    for(j=0;j<volsize*num;j++)fimage[j]=(float)image[j];
  }
  else{
    if(FslReadVolumes(ip,fimage,num)!=num){
      printf("Read error: %s\n",infile);
      return(1);
    }
  }
  
  /* Error counter */
  errcnt=0;
  
  /* Do the fitting */
  for(j=0;j<volsize;j++){

    /* Extract the yvalues */
    for(i=0;i<num;i++) yvals[i]=fimage[i*volsize+j];

    if((masked)&&(!mask[j])){
      for(i=0;i<ng;i++)data[i][j]=0;
    }
    else {

      if((func==INVREC)||(func==INVRECI)){
	for(i=0,slope=0;i<num-1;i++){
	  if((xvals[i+1]-xvals[i])>0)
	    slope+=(yvals[i+1]-yvals[i])/(xvals[i+1]-xvals[i]);
	}
	if(slope<0)for(i=0;i<num;i++)yvals[i]*=-1.0;
      }
      
      if(aim)a[0]=aimage[j];
      else a[0]=A;
      if(bim)a[1]=bimage[j];
      else a[1]=B;
      if(cim)global_flip=cimage[j]*M_PI/180;
      for(i=2;i<ng;i++)a[i]=1;

      /* Do the fit */
      errcnt+=lmfit(xvals,yvals,num,a,ng,func);
      for(i=0;i<ng;i++)data[i][j]=a[i];
    }

    /* Clip values */
    for(i=0;i<ng;i+=2){
      if((uA!=NOCLIP)&&(data[i][j]>uA))data[i][j]=uA;
      if((lA!=NOCLIP)&&(data[i][j]<lA))data[i][j]=lA;
      if((uB!=NOCLIP)&&(data[i+1][j]>uB))data[i+1][j]=uB;
      if((lB!=NOCLIP)&&(data[i+1][j]<lB))data[i+1][j]=lB;
    }

    /* Find max and min */
    if(j){
      for(i=0;i<ng;i++){
	if(data[i][j]>max[i])max[i]=data[i][j];
	if(data[i][j]<min[i])min[i]=data[i][j];
      }
    } else {
       for(i=0;i<ng;i++){    
	max[i]=data[i][j];
	min[i]=data[i][j];
       }
    }
  }

  /* Clip max and min */
  for(i=0;i<ng;i++){
    if(max[i]>32767)max[i]=32767;
    if(min[i]<-32767)min[i]=-32767;
  }
  
  /* Write out the data */

  for(i=0;i<nc;i++){

    if(nc>1)sprintf(tmpfile,"%s_A%d",outfile,i);
    else sprintf(tmpfile,"%s_A",outfile);
    if((op=FslOpen(tmpfile,"w"))==NULL){
      printf("Can't open : %s\n",tmpfile);
      return(1);

    }
    FslCloneHeader(op,ip);
    FslSetDataType(op,DT_FLOAT);
    FslSetDim(op,ppl,lpi,ipv,1);
    FslSetCalMinMax(op,(int)min[i*2],(int)max[i*2]);
    FslWriteHeader(op);
    FslWriteVolumes(op,data[i*2],1);
    FslClose(op);

    if(nc>1)sprintf(tmpfile,"%s_B%d",outfile,i);
    else sprintf(tmpfile,"%s_B",outfile);

    if((op=FslOpen(tmpfile,"w"))==NULL){
      printf("Can't open : %s\n",tmpfile);
      return(1);

    }
    FslCloneHeader(op,ip);
    FslSetDataType(op,DT_FLOAT);
    FslSetDim(op,ppl,lpi,ipv,1);
    FslSetCalMinMax(op,(int)min[i*2+1],(int)max[i*2+1]);
    FslWriteHeader(op);
    FslWriteVolumes(op,data[i*2+1],1);
    FslClose(op);

  }
  return(0);
}
int lmfit(float *x,float *y,int n,float *a,int ma,int function)
{
  int i,numit,rtn=0;
  int *ia;

  float	*sig,alamda,**covar,**alpha;
  float oldchisq,chisq;
  void (*func)(float, float [], float *, float [], int);
  
  void invrec(float x,float a[],float *y,float dyda[],int ma);
  void invreci(float x,float a[],float *y,float dyda[],int ma);
  void satrec(float x,float a[],float *y,float dyda[],int na);
  void expdec(float x,float a[],float *y,float dyda[],int na);
  void b1map(float x,float a[],float *y,float dyda[],int na);

  switch(function){
  case SATREC:
    func=satrec;
    break;
  case INVREC:
    func=invrec;
    break;
  case INVRECI:
    func=invreci;
    break;
  case EXPDEC:
    func=expdec;
    break;
  case B1MAP:
    func=b1map;
    break;
  default:
    return(0);
  }

  covar = matrix(1,ma,1,ma);
  alpha = matrix(1,ma,1,ma);
  sig = vector(1,n);
  ia = ivector(1,ma);

  /* Initialisation */
  for(i=1;i<=n;i++) sig[i] = 1;
  for(i=1;i<=ma;i++) ia[i] = 1;
  alamda = -1.0;
  oldchisq = 1.0e30;

  if(mrqmin(x-1,y-1,sig,n,a-1,ia,ma,covar,alpha,&chisq,func,&alamda)){
    for(i=0;i<ma;i++)a[i]=0;
    return(1);
  }

  numit = 1;

  while (((oldchisq==chisq)||(fabs((oldchisq-chisq)/chisq)>=0.0001))&&(alamda!=0.0)) {
    oldchisq = chisq;
    numit++;
    if(numit>MAXIT){
      /*printf("Too many iterations\n");*/
      rtn=1;
      break;
    }
    if (chisq<0.00005) break; 
    if(mrqmin(x-1,y-1,sig,n,a-1,ia,ma,covar,alpha,&chisq,func,&alamda)){
      for(i=0;i<ma;i++)a[i]=0;
      return(1);
    }
  }
  alamda = 0.0;
  if(mrqmin(x-1,y-1,sig,n,a-1,ia,ma,covar,alpha,&chisq,func,&alamda)){
    for(i=0;i<ma;i++)a[i]=0;
    return(1);
  }

  free_matrix(covar,1,ma,1,ma);
  free_matrix(alpha,1,ma,1,ma);
  free_vector(sig,1,n);
  free_ivector(ia,1,ma);

  return(rtn);
}
void invrec(float x,float a[],float *y,float dyda[],int na)
{
  int i;
  float arg,ex,mex;
  
  *y = 0.0;
  for (i=1; i<=na-1; i+=2) {
    arg = x/a[i+1];
    ex = exp(-arg);
    mex = 1-(2*ex);
    *y += a[i]*mex;
    dyda[i] = mex;
    dyda[i+1] = (-2*a[i]*x*ex)/(a[i+1]*a[i+1]);
  }  
}
void invreci(float x,float a[],float *y,float dyda[],int na)
{
  int i;
  float arg,ex,mex;
  
  /* NB This function uses the global variable global_flip */

  *y = 0.0;
  for (i=1; i<=na-1; i+=2) {
    arg = x/a[i+1];
    ex = exp(-arg);
    mex = 1-((1-cos(global_flip))*ex);
    *y += a[i]*mex;
    dyda[i] = mex;
    dyda[i+1] = (-2*a[i]*x*ex)/(a[i+1]*a[i+1]);
  }  
}
void satrec(float x,float a[],float *y,float dyda[],int na)
{
  int i;
  float arg,ex,mex;
  
  *y = 0.0;
  for (i=1; i<=na-1; i+=2) {
    arg = x/a[i+1];
    ex = exp(-arg);
    mex = 1-ex;
    *y += a[i]*mex;
    dyda[i] = mex;
    dyda[i+1] = (-a[i]*x*ex)/(a[i+1]*a[i+1]);
  }
}
void expdec(float x,float a[],float *y,float dyda[],int na)
{
  int i;
  float arg,ex;
  
  *y = 0.0;
  for (i=1; i<=na-1; i+=2) {
      arg = x/a[i+1];
      ex = exp(-arg);
      *y += a[i]*ex;
      dyda[i] = ex;
      dyda[i+1] = (a[i]*x*ex)/(a[i+1]*a[i+1]);
  }
}
void b1map(float x,float a[],float *y,float dyda[],int na)
{
  float sina;
  
  sina=sin((M_PI/2)*a[2]*x);
  *y = a[1]*sina*sina*sina;
  dyda[1]=sina*sina*sina;
  dyda[2]=(M_PI/2)*a[1]*3*sina*sina*cos((M_PI/2)*a[2]*x);
}
void usage_error()
{
  printf("usage: lmfit [-options]\n");
  printf("\t -in infile\n");
  printf("\t -func function_number\n");
  printf("\t -x num x1 x2 x3 ...\n");
  printf("\t [-nc num]\n");
  printf("\t [-Aim filename]\n");
  printf("\t [-Bim filename]\n");
  printf("\t [-Cim filename]\n");
  printf("\t [-out outfile]\n");
  printf("\t [-mask maskfile]\n");
  printf("\t [-Amin n -Amax n]\n");
  printf("\t [-Bmin n -Bmax n]\n");
  printf("\t -help\n");
  printf("\n");
  printf("functions:\n");
  printf("\t 2 - Saturation Recovery\n");
  printf("\t 3 - Exponential Decay\n");
  printf("\t 4 - Inversion Recovery\n");
  printf("\t 6 - Incomplete Inversion Recovery (with Cim)\n");
  printf("\n");
  exit(1);
}
