/*
  linfit: Image fitting program using linear regression
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
#include <string.h>
#include <math.h>
#include <fslio.h>

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
  char infile[MAXLEN],maskfile[MAXLEN],outfile[MAXLEN];
  char Afile[MAXLEN],Bfile[MAXLEN],Efile[MAXLEN];
  short ppl,lpi,ipv,vols,dt;
  short mppl,mlpi,mipv,mvols,mdt;
  int func,num;
  flag dtshort,masked,err;
  float *fimage,*xvals,*yvals;
  float *Adata,*Bdata,*Edata,m,c,em,ec;
  float maxA,minA,maxB,minB,maxE,minE;
  float uA,lA,uB,lB,uE,lE;
  short *image;
  char *mask;
  /*char orient,morient;*/
  long i,j,volsize;
  FSLIO *ip,*op,*mp;

  void usage_error();
  void linfit(float* x,float* y,int n,float* m,float* c,float* em,float* ec);

  /* Read in the arguments */
  if(argc==1)usage_error();

  func=num=err=0;
  infile[0]=outfile[0]=maskfile[0]='\0';
  uA=uB=uE=lA=lB=lE=NOCLIP;
  maxA=minA=maxB=minB=maxE=minE=0;
  Edata=NULL;
  yvals=xvals=NULL;
  mask=NULL;

  for(i=1;i<argc;i++){
    if((!strcmp(argv[i],"-h"))||(!strcmp(argv[i],"-help")))
      usage_error();
    if(!strcmp(argv[i],"-in")){sprintf(infile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-out")){sprintf(outfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-func")){sscanf(argv[i+1],"%d",&func);continue;}
    if(!strcmp(argv[i],"-mask")){sprintf(maskfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-Amax")){sscanf(argv[i+1],"%f\n",&uA);continue;}
    if(!strcmp(argv[i],"-Bmax")){sscanf(argv[i+1],"%f\n",&uB);continue;}
    if(!strcmp(argv[i],"-Emax")){sscanf(argv[i+1],"%f\n",&uE);continue;}
    if(!strcmp(argv[i],"-Amin")){sscanf(argv[i+1],"%f\n",&lA);continue;}
    if(!strcmp(argv[i],"-Bmin")){sscanf(argv[i+1],"%f\n",&lB);continue;}
    if(!strcmp(argv[i],"-Emin")){sscanf(argv[i+1],"%f\n",&lE);continue;}
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

  switch(func){
  case LINEAR:
  case EXPDEC:
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
    sprintf(Afile,"%s_A",infile);
    sprintf(Bfile,"%s_B",infile);
    if(err)sprintf(Efile,"%s_E",infile);
  }
  else {
    sprintf(Afile,"%s_A",outfile);
    sprintf(Bfile,"%s_B",outfile);
    if(err)sprintf(Efile,"%s_E",outfile);
  }

  if(maskfile[0]=='\0')masked=0;
  else masked=1;

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
  if((Adata=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  if((Bdata=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  if(err){
    if((Edata=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
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
      printf("ppl: %d, %d\n",ppl,mppl);
      printf("lpi: %d, %d\n",lpi,mlpi);
      printf("ipv: %d, %d\n",ipv,mipv);
      /*printf("orient: %d, %d\n",orient,morient);*/
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

  /* Do the fitting */
  for(j=0;j<volsize;j++){

    /* Extract the yvalues */
    for(i=0;i<num;i++){
      if(func==EXPDEC){
	yvals[i]=log(fimage[i*volsize+j]);
      } else {
	yvals[i]=fimage[i*volsize+j];
      }
    }

    if((masked)&&(!mask[j])){
      Adata[j]=0;
      Bdata[j]=0;
      if(err)Edata[j]=0;
    }
    else {
     
      linfit(xvals,yvals,num,&m,&c,&em,&ec);

      if(func==LINEAR){
	Adata[j]=m;
	Bdata[j]=c;
	if(err)Edata[j]=em;
      }
      if(func==EXPDEC){
	Adata[j]=exp(c);
	Bdata[j]=-1.0/m;
	if(err)Edata[j]=em/(m*m);
      }
    }

    /* Clip values */
    if((uA!=NOCLIP)&&(Adata[j]>uA))Adata[j]=uA;
    if((lA!=NOCLIP)&&(Adata[j]<lA))Adata[j]=lA;
    if((uB!=NOCLIP)&&(Bdata[j]>uB))Bdata[j]=uB;
    if((lB!=NOCLIP)&&(Bdata[j]<lB))Bdata[j]=lB;
    if(err){
      if((uE!=NOCLIP)&&(Edata[j]>uE))Edata[j]=uE;
      if((lE!=NOCLIP)&&(Edata[j]<lE))Edata[j]=lE;
    }

    /* Find max and min */

    if(j){
      if(Adata[j]>maxA)maxA=Adata[j];
      if(Adata[j]<minA)minA=Adata[j];
      if(Bdata[j]>maxB)maxB=Bdata[j];
      if(Bdata[j]<minB)minB=Bdata[j];
      if(err){
	if(Edata[j]>maxE)maxE=Edata[j];
	if(Edata[j]<minE)minE=Edata[j];
      }
    } else {
      maxA=Adata[j];
      minA=Adata[j];
      maxB=Bdata[j];
      minB=Bdata[j];
      if(err){
	maxE=Edata[j];
	minE=Edata[j];
      }
    }
  }

  /* Clip max and min */
  if(maxA>32767)maxA=32767;
  if(maxB>32767)maxB=32767;
  if(minA<-32767)minA=-32767;
  if(minB<-32767)minB=-32767;
  if(err){
    if(minE<-32767)minE=-32767;
    if(minE<-32767)minE=-32767;
  }

  /* Write out the data */
  if((op=FslOpen(Afile,"w"))==NULL){
    printf("Can't open : %s\n",Afile);
    return(1);
  }
  FslCloneHeader(op,ip);
  FslSetDataType(op,DT_FLOAT);
  FslSetDim(op,ppl,lpi,ipv,1);
  FslSetCalMinMax(op,(int)minA,(int)maxA);
  FslWriteHeader(op);
  FslWriteVolumes(op,Adata,1);
  FslClose(op);

  if((op=FslOpen(Bfile,"w"))==NULL){
    printf("Can't open : %s\n",Bfile);
    return(1);
  }
  FslCloneHeader(op,ip);
  FslSetDataType(op,DT_FLOAT);
  FslSetDim(op,ppl,lpi,ipv,1);
  FslSetCalMinMax(op,(int)minB,(int)maxB);
  FslWriteHeader(op);
  FslWriteVolumes(op,Bdata,1);
  FslClose(op);

  if(err){
    if((op=FslOpen(Efile,"w"))==NULL){
      printf("Can't open : %s\n",Efile);
      return(1);
    }
    FslCloneHeader(op,ip);
    FslSetDataType(op,DT_FLOAT);
    FslSetDim(op,ppl,lpi,ipv,1);
    FslSetCalMinMax(op,(int)minE,(int)maxE);
    FslWriteHeader(op);
    FslWriteVolumes(op,Edata,1);
    FslClose(op);
  
  }
  return(0);
}
void linfit(float* x,float* y,int n,float* m,float* c,float* em,float* ec)
{
  int j;
  float xme,yme,s1,s2,s3;

  for(j=0,yme=xme=0;j<n;j++){
    xme+=x[j];
    yme+=y[j];
  }
  xme/=(float)n;
  yme/=(float)n;
  
  for(j=0,s1=s2=0;j<n;j++){
    s1+=((x[j]-xme)*y[j]);
    s2+=((x[j]-xme)*(x[j]-xme));
  }
  *m=s1/s2;
  *c=yme-((*m)*xme);
		
  for(j=0,s3=0;j<n;j++){
    s3+=((y[j]-((*m)*x[j])-(*c))*(y[j]-((*m)*x[j])-(*c)));
  }
  *em=sqrt(s3/(((float)n-2)*s2));
  *ec=sqrt(((1/n)+(xme*xme/s2))*(s3/(n-2)));
}
void usage_error()
{
  printf("usage: linfit [-options]\n");
  printf("\t -in infile\n");
  printf("\t -func function_number\n");
  printf("\t -x num x1 x2 x3 ...\n");
  printf("\t [-out outfile]\n");
  printf("\t [-mask maskfile]\n");
  printf("\t [-err]\n");
  printf("\t [-Amin n -Amax n]\n");
  printf("\t [-Bmin n -Bmax n]\n");
  printf("\t [-Emin n -Emax n]\n");
  printf("\t -help\n");
  printf("\n");
  printf("functions:\n");
  printf("\t 1 - Linear Y = A.X + B\n");
  printf("\t 3 - Exponential Decay  Y = A.exp(-X/B)\n");
  printf("\n");
  exit(1);
}
