/*
  pssreorder.c: Reorders the slices based on the pss parameter
                in the Varian procpar file

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "pssreorder.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "maxlen.h"
#include "hfunc.h"

int pss_to_array(char* filename,float** pss_array)
{
  int i,n,order,skip;
  long nv2;
  char procfile[512],proc_string[5000];
  FILE *ipp;

  strcpy(procfile,filename);
  strcat(procfile,"/procpar");

  if((ipp=fopen(procfile,"rb"))==NULL){
    printf("Can't open input procpar file %s.\n",procfile);
    return(0);
  }
  
  while(fgets(proc_string,5000,ipp)!=NULL){
    if (!strncmp(proc_string, "nv2 ",4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %ld\n", &nv2);
      break;
    }
  }
  if(nv2>1)return(0);

  rewind(ipp);
  while(fgets(proc_string,5000,ipp)!=NULL){
    if (!strncmp(proc_string, "pss ", 4)){
      fscanf(ipp,"%d",&n);
      break;
    }
  }

  if((*pss_array=(float *)malloc(n*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(0);
  }
  for(i=0;i<n;i++)fscanf(ipp,"%f",&((*pss_array)[i]));
  fclose(ipp);

  /* Determining the reorder */
  order = 0;
  for(i=1;i<n;i++){
    if((*pss_array)[i-1]<(*pss_array)[i])
      order=1;
  }
  if(!order)return(0);

  order = 0;
  for(i=1;i<n;i++){
    if((*pss_array)[i-1]>(*pss_array)[i])
      order=1;
  }
  if(!order)return(1);

  order = 0; skip = (n+1)/2;
  for(i=0;i<n/2;i++){
    if((*pss_array)[i]>(*pss_array)[i+skip])
      order=1;
    if(i+1>=skip)continue;
    if((*pss_array)[i+skip]>(*pss_array)[i+1])
      order=1;
  }
  if(!order)return(2);

  order = 0; skip = (n+1)/2;
  for(i=0;i<n/2;i++){
    if((*pss_array)[i]<(*pss_array)[i+skip])
      order=1;
    if(i+1>=skip)continue;
    if((*pss_array)[i+skip]<(*pss_array)[i+1])
      order=1;
  }
  if(!order)return(3);

  /* Abnormal slice order */
  return(99);
}

void pss_reorder(float* data,long ppi,int ipv,float* array)
{
  int i,j,p=0;
  float *out,hi;
  int *got;

  out=(float *)malloc(ppi*ipv*sizeof(float));

  got=(int *)malloc(ipv*sizeof(int));
  for(i=0;i<ipv;i++)got[i]=0;

  for(i=0;i<ipv;i++){
    hi=-999;
    for(j=0;j<ipv;j++){
      if((!got[j])&&(array[j]>hi)){
	p=j;
	hi=array[j];
      }
    }
    got[p]=1;
    memcpy(&out[i*ppi],&data[p*ppi],ppi*sizeof(float));
  }
  for(i=0;i<ppi*ipv;i++)data[i]=out[i];

  free(out);
}

int output_pss(char* infile,char* outfile)
{
  int i,n;
  long nv2;
  char tmp[MAXLEN],proc_string[5000];
  float pss;
  FILE *ipp,*ofp;

  sprintf(tmp,"%s/procpar",infile);
  if((ipp=fopen(tmp,"rb"))==NULL){
    printf("Can't open input procpar file %s.\n",tmp);
    return(1);
  }

  if(!strcmp(outfile,"-nosave")){
    ofp=stdout;
  }
  else {
    sprintf(tmp,"%s.pss",outfile);
    if((ofp=fopen(tmp,"wb"))==NULL){
      printf("Can't open output file %s.\n",tmp);
      return(1);
    }
  }
  
  while(fgets(proc_string,5000,ipp)!=NULL){
    if (!strncmp(proc_string, "nv2 ",4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %ld\n", &nv2);
      break;
    }
  }
  if(nv2>1){
    printf("output_pss: Image is acquired as 3d!\n");
    return(1);
  }

  rewind(ipp);
  while(fgets(proc_string,5000,ipp)!=NULL){
    if (!strncmp(proc_string, "pss ", 4)){
      fscanf(ipp,"%d",&n);
      break;
    }
  }

  for(i=0;i<n;i++){
    fscanf(ipp,"%f",&pss);
    fprintf(ofp,"%f\n",pss);
  }
  fclose(ipp);
  if(strcmp(outfile,"-nosave"))fclose(ofp);
  return(0);
}

int mean_pss(char* filename,float* pss)
{
  int i,n;
  char procfile[512],proc_string[5000];
  float tmp,sum;
  FILE *ipp;

  strcpy(procfile,filename);
  strcat(procfile,"/procpar");

  if((ipp=fopen(procfile,"rb"))==NULL){
    printf("Can't open input procpar file %s.\n",procfile);
    return(1);
  }

  while(fgets(proc_string,5000,ipp)!=NULL){
    if (!strncmp(proc_string, "pss ", 4)){
      fscanf(ipp,"%d",&n);
      break;
    }
  }
  sum=0;
  for(i=0;i<n;i++){
    fscanf(ipp,"%f",&tmp);
    sum+=tmp;
  }

  *pss=sum/n;
  fclose(ipp);
  return(0);
}

