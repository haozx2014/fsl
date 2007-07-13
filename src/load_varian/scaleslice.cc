/*
  scaleslice.c: Functions to scale slice intensities
                from values in a file

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

/*
  Please note:
  This routine requires a look up table of the form
  [Slice offset] -tab- [Scale factor]

  for example
  -14.0    1.20
  -13.0    1.15
  -12.0    1.00 etc.

  This file must be named after the RF coil it represents
  eg head.scale, surface.scale
  and must be in the directory ${FSLDIR}/tcl
  Also, the environment variable FSLDIR must be set
*/


#include "scaleslice.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>

void scale_slice(float* data,long pts,int ipv,float* scale)
{
  long l;
  int i;
  for(i=0;i<ipv;i++){
    for(l=0;l<pts;l++){
      data[(i*pts)+l]*=scale[i];
    }
  }
}
int scale_slice_array(char* filename,float** scale)
{
  int i,n,gotit=0;
  unsigned l;
  long nv2;
  char procfile[512],proc_string[5000],rfcoil[1024];
  char modfile[1024];
  float pss_array[512];
  float pwr1,pwr2,pss1,pss2;
  FILE *ipp,*mfp;

  strcpy(procfile,filename);
  strcat(procfile,"/procpar");

  if((ipp=fopen(procfile,"rb"))==NULL){
    printf("Can't open input procpar file %s.\n",procfile);
    return(1);
  }
  while(fgets(proc_string,5000,ipp)!=NULL){
    if (!strncmp(proc_string, "rfcoil ", 7)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %s\n", rfcoil);
      for (l = 0; l < strlen(rfcoil) - 2; l++)
	rfcoil[l] = rfcoil[l + 1];
      rfcoil[l] = '\0';
      gotit=1;
    }
    if (!strncmp(proc_string, "pss ", 4)){
      fscanf(ipp,"%d",&n);
      for(i=0;i<n;i++)fscanf(ipp,"%f",&pss_array[i]);
    }
    if (!strncmp(proc_string, "nv2 ",4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %ld\n", &nv2);
    }
  }
  fclose(ipp);
  
  if(!gotit)return(1);
  if(nv2>1)return(1);

  if(((*scale)=(float *)malloc(n*sizeof(float)))==NULL){
    printf("Malloc failed 1");
    return(1);
  }

  strcpy(modfile,getenv("FSLDIR"));
  strcat(modfile,"/tcl/");
  strcat(modfile,rfcoil);
  strcat(modfile,".scale");
  
  if((mfp=fopen(modfile,"rb"))==NULL){
    printf("Can't open file %s\n",modfile);
    return(1);
  }
  pwr1=pwr2=pss1=0.0;

  for(i=0;i<n;i++){
    rewind(mfp);
    while(fgets(proc_string,5000,mfp)!=NULL){
      sscanf(proc_string,"%f\t%f",&pss2,&pwr2);
      if(pss2>pss_array[i]){
	(*scale)[i]=pwr1;
	pwr2-=pwr1;
	pwr2*=((pss_array[i]-pss1)/(pss2-pss1));
	(*scale)[i]+=pwr2;
	break;
      }
      pss1=pss2;
      pwr1=pwr2;
    }
  }
  fclose(mfp);
  return(0);
}
int scale_slice_array_file(char* filename,float** scale,char* modfile)
{
  int i,n;
  long nv2;
  char procfile[512],proc_string[5000];
  float pss_array[512];
  float pwr1,pwr2,pss1,pss2;
  FILE *ipp,*mfp;

  strcpy(procfile,filename);
  strcat(procfile,"/procpar");

  if((ipp=fopen(procfile,"rb"))==NULL){
    printf("Can't open input procpar file %s.\n",procfile);
    return(1);
  }

  while(fgets(proc_string,5000,ipp)!=NULL){
    if (!strncmp(proc_string, "nv2 ",4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %ld\n", &nv2);
    }
  }
  if(nv2>1)return(1);
 
  while(fgets(proc_string,5000,ipp)!=NULL){
    if (!strncmp(proc_string, "pss ", 4)){
      fscanf(ipp,"%d",&n);
      for(i=0;i<n;i++)fscanf(ipp,"%f",&pss_array[i]);
    }
  }
  fclose(ipp);
  
  if(((*scale)=(float *)malloc(n*sizeof(float)))==NULL){
    printf("Malloc failed 2");
    return(1);
  }
  
  if((mfp=fopen(modfile,"rb"))==NULL){
    printf("Can't open file %s\n",modfile);
    return(1);
  }
  pwr1=pwr2=pss1=0.0;

  for(i=0;i<n;i++){
    rewind(mfp);
    while(fgets(proc_string,5000,mfp)!=NULL){
      sscanf(proc_string,"%f\t%f",&pss2,&pwr2);
      if(pss2>pss_array[i]){
	(*scale)[i]=pwr1;
	pwr2-=pwr1;
	pwr2*=((pss_array[i]-pss1)/(pss2-pss1));
	(*scale)[i]+=pwr2;
	break;
      }
      pss1=pss2;
      pwr1=pwr2;
    }
  }
  fclose(mfp);
  return(0);
}
