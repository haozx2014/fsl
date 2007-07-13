/*
  minfit: Image fitting program using powells method
  Part of the iFit suite of programs for image fitting to functions

  Stuart Clare, FMRIB Physics Group

  Copyright (C) 1999-2000 University of Oxford
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
#include "vector.h"
#include "matrix.h"
#include "nr.h"
#include"fslio/fslio.h"

#define LINEAR 1
#define SATREC 2
#define EXPDEC 3
#define INVREC 4
#define INVRECM 5
#define INVRECI 6
#define INVRECIM 7

#define MAXLEN 128
#define FTOL 1.0e-20
#define ITMAX 100

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

/* Global variables */
float *xvals,*yvals,cval;
int num,aim,cim;

int main(int argc,char* argv[])
{
  char infile[MAXLEN],maskfile[MAXLEN],outfile[MAXLEN];
  char Afile[MAXLEN],Bfile[MAXLEN],Cfile[MAXLEN],Efile[MAXLEN];
  char Aimfile[MAXLEN],Cimfile[MAXLEN];
  int func,err,mdim,masked,iter;
  float A,B,C,*p,**xi,fret;
  float *fimage,*Adata,*Bdata,*Cdata,*Edata;
  float *Aimage,*Cimage;
  float Amax,Amin,Bmax,Bmin,Cmax,Cmin;
  float uA,uB,uC,uE,lA,lB,lC,lE,slope;
  short *image;
  char *mask;
  /*char orient,morient;*/
  short ppl,lpi,ipv,vols,dt,dtshort;
  short mppl,mlpi,mipv,mvols,mdt;
  long i,j,volsize;
  FSLIO *ip,*op,*mp;

  void usage_error();
  float linear(float p[]);
  float expdec(float p[]);
  float invrec(float p[]);
  float satrec(float p[]);
  float invrecm(float p[]);
  float invreci(float p[]);
  float invrecim(float p[]);

  /* Read in the arguments */
  if(argc==1)usage_error();

  func=err=num=0;
  A=B=C=0.0;
  Amax=Bmax=Cmax=Amin=Bmin=Cmin=-1.0;
  infile[0]=outfile[0]=maskfile[0]=Aimfile[0]=Cimfile[0]='\0';
  Aimage=Cimage=Cdata=NULL;
  mask=NULL;

  for(i=1;i<argc;i++){
    if((!strcmp(argv[i],"-h"))||(!strcmp(argv[i],"-help")))
      usage_error();
    if(!strcmp(argv[i],"-in")){sprintf(infile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-out")){sprintf(outfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-func")){sscanf(argv[i+1],"%d",&func);continue;}
    if(!strcmp(argv[i],"-mask")){sprintf(maskfile,"%s",argv[i+1]);continue;}
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
    if(!strcmp(argv[i],"-Aim")){sprintf(Aimfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-Cim")){sprintf(Cimfile,"%s",argv[i+1]);continue;}
    if(!strcmp(argv[i],"-A")){sscanf(argv[i+1],"%f\n",&A);continue;}
    if(!strcmp(argv[i],"-B")){sscanf(argv[i+1],"%f\n",&B);continue;}
    if(!strcmp(argv[i],"-C")){sscanf(argv[i+1],"%f\n",&C);continue;}
    if(!strcmp(argv[i],"-Amax")){sscanf(argv[i+1],"%f\n",&Amax);continue;}
    if(!strcmp(argv[i],"-Bmax")){sscanf(argv[i+1],"%f\n",&Bmax);continue;}
    if(!strcmp(argv[i],"-Cmax")){sscanf(argv[i+1],"%f\n",&Cmax);continue;}
    if(!strcmp(argv[i],"-Amin")){sscanf(argv[i+1],"%f\n",&Amin);continue;}
    if(!strcmp(argv[i],"-Bmin")){sscanf(argv[i+1],"%f\n",&Bmin);continue;}
    if(!strcmp(argv[i],"-Cmin")){sscanf(argv[i+1],"%f\n",&Cmin);continue;}
    if(!strcmp(argv[i],"-err")){err=1;continue;}
    if(argv[i][0]=='-'){
      printf("Command flag %s not recognised\\n",argv[i]);
      usage_error();
    }
}

  switch(func){
  case(LINEAR):
  case(SATREC):
  case(EXPDEC):
  case(INVREC):
  case(INVRECM):
    mdim=2;
    break;
  case(INVRECI):
  case(INVRECIM):
    mdim=3;
    break;
  default:
    printf("Function number not specified or recognised\n");
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
    if(mdim>2)sprintf(Cfile,"%s_C",infile);
    if(err){
      sprintf(Efile,"%s_err",infile);
    }
  }
  else {
    sprintf(Afile,"%s_A",outfile);
    sprintf(Bfile,"%s_B",outfile);
    if(mdim>2)sprintf(Cfile,"%s_C",outfile);
    if(err){
      sprintf(Efile,"%s_E",outfile);
    }
  }

  if(maskfile[0]=='\0')masked=0;
  else masked=1;
  if(Aimfile[0]=='\0')aim=0;
  else aim=1;
  if(Cimfile[0]=='\0')cim=0;
  else cim=1;

  if(cim)mdim=2;

  p=vector(1,mdim);
  xi=matrix(1,mdim,1,mdim);

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
  if((Edata=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  if(mdim>2){
    if((Cdata=(float *)malloc(volsize*sizeof(float)))==NULL){
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
    if((mp=FslOpen(Aimfile,"r"))==NULL){
      printf("Failed to open AVW file: %s\n",Aimfile);
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
    if((Aimage=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
    if(mdt==DT_SIGNED_SHORT){
      if(FslReadVolumes(mp,image,1)!=1){
	printf("Read error: %s\n",Aimfile);
	return(1);
      }
      for(j=0;j<volsize;j++)Aimage[j]=(float)image[j];
    }
    else{
      if(FslReadVolumes(mp,Aimage,1)!=1){
	printf("Read error: %s\n",Aimfile);
	return(1);
      }
    }
    FslClose(mp);
  }

  /* Read in C image */
  if(cim){
    if((mp=FslOpen(Cimfile,"r"))==NULL){
      printf("Failed to open AVW file: %s\n",Cimfile);
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
    if((Cimage=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed.\n");
      return(1);
    }
    if(mdt==DT_SIGNED_SHORT){
      if(FslReadVolumes(mp,image,1)!=1){
	printf("Read error: %s\n",Cimfile);
	return(1);
      }
      for(j=0;j<volsize;j++)Cimage[j]=(float)image[j];
    }
    else{
      if(FslReadVolumes(mp,Cimage,1)!=1){
	printf("Read error: %s\n",Cimfile);
	return(1);
      }
    }
    FslClose(mp);
  }

  /* Read in the data */
  FslSeekVolume(ip,0);
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

  uA=uB=uC=uE=lA=lB=lC=lE=0.0;

  /* Do the fitting */
  for(j=0;j<volsize;j++){

    /* Extract the yvalues */
    for(i=0;i<num;i++){
      yvals[i]=fimage[i*volsize+j];
    }
    
    if((masked)&&(!mask[j])){
      Adata[j]=0;
      Bdata[j]=0;
      if(mdim>2)Cdata[j]=0;
    }
    else {
      /* Initialise vector and matrix */
      xi[1][1]=1;xi[1][2]=0;xi[2][1]=0;xi[2][2]=1;
      if(mdim>2){xi[1][3]=0;xi[2][3]=0;xi[3][1]=0;xi[3][2]=0;xi[3][3]=1;}

      /* If using A image then order of fitting is swapped */
      if(aim){
	if(mdim>2){p[1]=B;p[2]=C;p[3]=Aimage[j];}
	else {p[1]=B;p[2]=Aimage[j];}
      }
      else {
	if(mdim>2){p[1]=A;p[2]=B;p[3]=C;}
	else {p[1]=A;p[2]=B;}
      }
      
      /* cval is a global variable */
      if(cim){
	cval=Cimage[j];
      }
      
      switch(func){
      case LINEAR:
	powell_lim(p,xi,mdim,FTOL,&iter,&fret,linear,ITMAX);
	break;
      case SATREC:
	powell_lim(p,xi,mdim,FTOL,&iter,&fret,satrec,ITMAX);
	break;
      case EXPDEC:
	powell_lim(p,xi,mdim,FTOL,&iter,&fret,expdec,ITMAX);
	break;
      case INVRECM:
	powell_lim(p,xi,mdim,FTOL,&iter,&fret,invrecm,ITMAX);
	break;
      case INVRECIM:
	powell_lim(p,xi,mdim,FTOL,&iter,&fret,invrecim,ITMAX);
	break;
      case INVREC:
	for(i=0,slope=0;i<num-1;i++){
	  if((xvals[i+1]-xvals[i])>0)
	    slope+=(yvals[i+1]-yvals[i])/(xvals[i+1]-xvals[i]);
	}
	if(slope<0)for(i=0;i<num;i++)yvals[i]*=-1.0;
	powell_lim(p,xi,mdim,FTOL,&iter,&fret,invrec,ITMAX);
	break;
      case INVRECI:
	for(i=0,slope=0;i<num-1;i++){
	  if((xvals[i+1]-xvals[i])>0)
	    slope+=(yvals[i+1]-yvals[i])/(xvals[i+1]-xvals[i]);
	}
	if(slope<0)for(i=0;i<num;i++)yvals[i]*=-1.0;
	powell_lim(p,xi,mdim,FTOL,&iter,&fret,invreci,ITMAX);
	break;
      }
      if(aim){
	if(mdim>2){Bdata[j]=p[1];Cdata[j]=p[2];Adata[j]=p[3];}
	else {Bdata[j]=p[1];Adata[j]=p[2];}
      }
      else {
	if(mdim>2){Adata[j]=p[1];Bdata[j]=p[2];Cdata[j]=p[3];}
	else {Adata[j]=p[1];Bdata[j]=p[2];}
      }

      if(mdim>2){
	Cdata[j]=acos(cos(Cdata[j]*M_PI/180))*180/M_PI;
      }

      if(err){
	Edata[j]=sqrt(fret);
	if(Edata[j]>uE)uE=Edata[j];
	if(Edata[j]<lE)lE=Edata[j];
      }
      if((Amax!=-1.0)&&(Adata[j]>Amax))Adata[j]=Amax;
      if((Bmax!=-1.0)&&(Bdata[j]>Bmax))Bdata[j]=Bmax;
      if((Amin!=-1.0)&&(Adata[j]<Amin))Adata[j]=Amin;
      if((Bmin!=-1.0)&&(Bdata[j]<Bmin))Bdata[j]=Bmin;

      if(Adata[j]>uA)uA=Adata[j];
      if(Bdata[j]>uB)uB=Bdata[j];
      if(Adata[j]<lA)lA=Adata[j];
      if(Bdata[j]<lB)lB=Bdata[j];
      if(mdim>2){
	if((Cmax!=-1.0)&&(Cdata[j]>Cmax))Cdata[j]=Cmax;
	if((Cmin!=-1.0)&&(Cdata[j]<Cmin))Cdata[j]=Cmin;
	if(Cdata[j]>uC)uC=Cdata[j];
	if(Cdata[j]<lC)lC=Cdata[j];
      }      

    }
  }
  
  if(uA>32767)uA=32767;
  if(uB>32767)uB=32767;
  if(uC>32767)uC=32767;
  if(uE>32767)uE=32767;
  if(uA<-32767)uA=-32767;
  if(uB<-32767)uB=-32767;
  if(uC<-32767)uC=-32767;
  if(uE<-32767)uE=-32767;

  /* Write out the data */
  if((op=FslOpen(Afile,"w"))==NULL){
    printf("Can't open : %s\n",Afile);
    return(1);  
  }
  FslCloneHeader(op,ip);
  FslSetDataType(op,DT_FLOAT);
  FslSetDim(op,ppl,lpi,ipv,1);
  FslSetCalMinMax(op,(int)lA,(int)uA);
  FslWriteHeader(op);
  FslWriteVolumes(op,Adata,1);
  FslClose(op);
  
  if((op=FslOpen(Bfile,"w"))==NULL){
    printf("Can't open : %s\n",Afile);
    return(1);
  }
  FslCloneHeader(op,ip);
  FslSetDataType(op,DT_FLOAT);
  FslSetDim(op,ppl,lpi,ipv,1);
  FslSetCalMinMax(op,(int)lB,(int)uB);
  FslWriteHeader(op);
  FslWriteVolumes(op,Bdata,1);
  FslClose(op);
  
  if(mdim>2){
    if((op=FslOpen(Cfile,"w"))==NULL){
      printf("Can't open : %s\n",Afile);
      return(1);  
    }
    FslCloneHeader(op,ip);
    FslSetDataType(op,DT_FLOAT);
    FslSetDim(op,ppl,lpi,ipv,1);
    FslSetCalMinMax(op,(int)lC,(int)uC);
    FslWriteHeader(op);
    FslWriteVolumes(op,Cdata,1);
    FslClose(op);
  }

  if(err){
    if((op=FslOpen(Efile,"w"))==NULL){
      printf("Can't open : %s\n",Efile);
      return(1);  
    }
    FslCloneHeader(op,ip);
    FslSetDataType(op,DT_FLOAT);
    FslSetDim(op,ppl,lpi,ipv,1);
    FslSetCalMinMax(op,(int)lE,(int)uE);
    FslWriteHeader(op);
    FslWriteVolumes(op,Edata,1);
    FslClose(op);
  }

  return(0);
}
float linear(float p[])
{
  int i,n;
  float A,B,Y,ssq=0;

  if(aim){B=p[1];A=p[2];}
  else {A=p[1];B=p[2];}

  for(i=0,n=0;i<num;i++){
    Y=A*xvals[i]+B;
    if(yvals[i]==0.0)continue;
    else n++;
    ssq+=((Y-yvals[i])*(Y-yvals[i]));
  }
  return(ssq/n);
}
float expdec(float p[])
{
  int i,n;
  float A,B,Y,ssq=0;

  if(aim){B=p[1];A=p[2];}
  else {A=p[1];B=p[2];}

  for(i=0,n=0;i<num;i++){
    Y=A*exp((-1)*xvals[i]/B);
    if(yvals[i]==0.0)continue;
    else n++;
    ssq+=((Y-yvals[i])*(Y-yvals[i]));
  }
  return(ssq/n);
}
float invrec(float p[])
{
  int i,n;
  float A,B,Y,x,ssq=0;

  if(aim){B=p[1];A=p[2];}
  else {A=p[1];B=p[2];}

  for(i=0,n=0;i<num;i++){
    x=xvals[i];
    if(yvals[i]==0.0)continue;
    else n++;
    Y=A*(1-2*exp((-1)*x/B));
    ssq+=((Y-yvals[i])*(Y-yvals[i]));
  }
  return(ssq/n);
}
float satrec(float p[])
{
  int i,n;
  float A,B,Y,x,ssq=0;

  if(aim){B=p[1];A=p[2];}
  else {A=p[1];B=p[2];}

  for(i=0,n=0;i<num;i++){
    x=xvals[i];
    Y=A*(1-exp(-1*x/B));
    if(yvals[i]==0.0)continue;
    else n++;
    ssq+=((Y-yvals[i])*(Y-yvals[i]));
  }
  return(ssq/n);
}
float invrecm(float p[])
{
  int i,n;
  float A,B,Y,x,ssq=0;
  
  if(aim){B=p[1];A=p[2];}
  else {A=p[1];B=p[2];}

  for(i=0,n=0;i<num;i++){
    x=xvals[i];
    Y=A*fabs(1-2*exp((-1)*x/B));
    if(yvals[i]==0.0)continue;
    else n++;
    ssq+=((Y-yvals[i])*(Y-yvals[i]));
  }
  
  return(ssq/n);
}
float invreci(float p[])
{
  int i,n;
  float A,B,C,Y,x,ssq=0;

  if(cim){
    if(aim){B=p[1];A=p[2];}
    else {A=p[1];B=p[2];}
    C=cos(cval*M_PI/180);
  }
  else {
    if(aim){A=p[3];B=p[1];C=cos(p[2]*M_PI/180);}
    else {A=p[1];B=p[2];C=cos(p[3]*M_PI/180);} 
  }

  for(i=0,n=0;i<num;i++){
    x=xvals[i];
    Y=A*(1-((1-C)*exp((-1)*x/B)));
    if(yvals[i]==0.0)continue;
    else n++;
    ssq+=((Y-yvals[i])*(Y-yvals[i]));
  }
  return(ssq/n);
}
float invrecim(float p[])
{
  int i,n;
  float A,B,C,Y,x,ssq=0;

  if(cim){
    if(aim){B=p[1];A=p[2];}
    else {A=p[1];B=p[2];}
    C=cos(cval*M_PI/180);
  }
  else {
    if(aim){A=p[3];B=p[1];C=cos(p[2]*M_PI/180);}
    else {A=p[1];B=p[2];C=cos(p[3]*M_PI/180);}
  }

  for(i=0,n=0;i<num;i++){
    x=xvals[i];
    Y=A*fabs(1-((1-C)*exp((-1)*x/B)));
    if(yvals[i]==0.0)continue;
    else n++;
    ssq+=((Y-yvals[i])*(Y-yvals[i]));
  }
  return(ssq/n);
}
void usage_error()
{
  printf("\nusage: fitting [-options]\n");
  printf("\t -in infile\n");
  printf("\t -func function_number(1-7)\n");
  printf("\t -x num x1 x2 x3 ...\n");
  printf("\t [-out outfile]\n");
  printf("\t [-mask maskfile]\n");
  printf("\t [-A aguess -B bguess -C cguess]\n");
  printf("\t [-Aim Afile]\n");
  printf("\t [-Cim Cfile]\n");
  printf("\t [-Amax n -Amin n -Bmax n -Bmin n -Cmax n -Cmin n]\n");
  printf("\t [-err]\n");
  printf("\t -help\n\n");
  printf("functions:\n");
  printf("\t 1 - Linear  Y = A.X + B\n");
  printf("\t 2 - Saturation Recovery  Y = A.(1 - exp(-X/B))\n");
  printf("\t 3 - Exponential Decay  Y = A.exp(-X/B)\n");
  printf("\t 4 - Inversion Recovery  Y = A.(1 - 2exp(-X/B))\n");
  printf("\t 5 - Modulus Inv Rec  Y = A.|1 - 2exp(-X/B)|\n");
  printf("\t 6 - Incomplete Inv Rec  Y = A.(1 - (1-cosC).exp(-X/B))\n");
  printf("\t 7 - Modulus Incomplete IR : Y = A.|1 - (1-cosC).exp(-X/B)|\n\n");
  exit(1);
}
