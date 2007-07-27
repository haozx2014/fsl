/*
  read_varian.c: Converts Varian fid files to complex AVW files

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "read_varian.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "hfunc.h"
#include "data.h"
#include "byteorder.h"
#include "data_byte.h"
#include "axis_scale.h"
#include "pssreorder.h"
#include "xfm.h"
#include "dbh.h"

#include "maxlen.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

#define GAMMA 4267.6995

int read_varian_header(char *infile1,struct dsr *header,struct xfm *transform)
{
  char procfile[MAXLEN],infile[MAXLEN],datestr[MAXLEN];
  char proc_string[5000],navs[MAXLEN],orient=0;
  long np,nv,ns,ne,nv2,vols;
  int phi,theta,psi,nD,nav,num_ints,epik_kh=0,prescan=0;
  int t_theta,t_psi,fract_ky;
  float lro,t_lro,lpe,t_thk,t_gap,thk,vx,vy,vz;
  float slice1,slice2,tr,TR;
  struct datafilehead fh;
  struct datafilehead_byte fhb;
  struct datablockhead_byte bhb;
  long i,date=0;
  unsigned long bsize;
  float *fdata;
  char *data;
  FILE *fp;
  FILE *ipp;
  int FHEADB,HEADB;
  int reverse,sshort,slong;

  FHEADB=sizeof(fhb);
  HEADB=sizeof(bhb);

  transform->rotate=1;
  transform->reflect=0;

  strcpy(infile,infile1);
  strcpy(procfile,infile);
  strcat(infile,"/fid");
  strcat(procfile,"/procpar");

  /* Open files */
  if((fp=fopen(infile,"rb"))==NULL){
    printf("Cannot open file: %s\n",infile);
    return(1);
  }

  /* Read header */
  if(!fread(&fhb,sizeof(fhb),1,fp)){
    printf("Read error\n");
    return(1);
  }

  if(find_byte_order(&reverse,&sshort,&slong)){
    printf("Cannot cope with byte size\n");
    return(1);
  }

  convert_filehead(&fhb,&fh,reverse,sshort,slong);

  bsize=(fh.bbytes-HEADB)/fh.ebytes;
  data=(char *)malloc(fh.ebytes*bsize);
  fdata=(float *)malloc(sizeof(float)*bsize);

  if((fh.ebytes!=2)&&(fh.ebytes!=4)){
    printf("Data type not supported\n");
    return(1);
  }

  /* Find image dimensions */
  if((ipp=fopen(procfile,"rb"))==NULL){
    printf("Can't open input procpar file %s.\n",procfile);
    return(1);
  }

  nav=0;TR=0.0;fract_ky=0;t_lro=0;
  while(fgets(proc_string,5000,ipp)!=NULL){
    if(!strncmp(proc_string,"np ",3)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %ld\n", &np);
    }
    if(!strncmp(proc_string,"nv ",3)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %ld\n", &nv);
    }
    if(!strncmp(proc_string,"ne ",3)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %ld\n", &ne);
    }
    if(!strncmp(proc_string,"nv2 ",4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %ld\n", &nv2);
    }
    if(!strncmp(proc_string,"nD ",3)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %d\n", &nD);
    }
    if (!strncmp(proc_string, "pss ", 4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "%ld", &ns);
      if(ns>1){
	sscanf(proc_string,"%ld %f %f",&ns,&slice1,&slice2);
      }
    }
    if (!strncmp(proc_string, "lro ", 4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %f\n", &lro);
      lro=lro*10;
    }
    if (!strncmp(proc_string, "actual_lro ", 11)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %f\n", &t_lro);
      t_lro=t_lro*10;
    }
    if (!strncmp(proc_string, "lpe ", 4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %f\n", &lpe);
      lpe*=10;
    }
    if (!strncmp(proc_string, "thk ", 4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %f\n", &thk);
    }
    if (!strncmp(proc_string, "t_thk ", 6)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %f\n", &t_thk);
    }
    if (!strncmp(proc_string, "t_gap ", 6)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %f\n", &t_gap);
    }
    if (!strncmp(proc_string, "phi ", 4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %d\n", &phi);
    }
    if (!strncmp(proc_string, "theta ", 6)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %d\n", &theta);
    }
    if (!strncmp(proc_string, "psi ", 4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %d\n", &psi);
    }
    if (!strncmp(proc_string, "nav ", 4)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %s\n", navs);
      if(navs[1]=='y')nav=1;
    }
    if (!strncmp(proc_string, "num_ints ", 9)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %d\n", &num_ints);
    }
    if (!strncmp(proc_string, "epik_kh ", 8)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %d\n", &epik_kh);
    }
    if (!strncmp(proc_string, "prescan_flag ", 13)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %d\n", &prescan);
    }
    if (!strncmp(proc_string, "tr ", 3)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %f\n", &tr);
    }
    if (!strncmp(proc_string, "TR ", 3)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %f\n", &TR);
    }
    if (!strncmp(proc_string, "fract_ky ", 9)){
      fgets(proc_string, 5000, ipp);
      sscanf(proc_string, "1 %d\n", &fract_ky);
    }
    if (!strncmp(proc_string, "date ", 5)){
      fgets(proc_string, 5000, ipp);
      i=3;
      while(proc_string[i]!='\"'){
	datestr[i-3]=proc_string[i];
	i++;
      }
      datestr[i-3] = '\0';
      date=str_to_date(datestr);
    }
  }
  fclose(ipp);

  if(nv2>1){
    ns=nv2;
  }
  if(nD==3){
    vz=t_thk;
  }
  else {
    vz=thk+t_gap;
  }
  if(nav){
    nv+=num_ints;
  }

  if(t_lro>0) lro=t_lro;
  vx=lro/(np/2);
  
  if(epik_kh){
    nv = nv+((nv/epik_kh)*(num_ints-1));
    if(prescan){
      nv=nv/num_ints;
    }
  }
  if(fract_ky){
    if(fract_ky<nv/2){
      nv=nv/2+fract_ky;
    }
  }
  vy=lpe/nv;
  vols = (bsize*fh.nblocks)/(np*nv*ns*ne);

  orient=3;
  t_theta=t_psi=0;
  if((theta<30)&&(theta>-30))t_theta=0;
  if((theta<120)&&(theta>60))t_theta=90;
  if((theta<-60)&&(theta>-120)){t_theta=90;transform->reflect=1;};

  if((psi<30)&&(psi>-30))t_psi=0;
  if((psi<120)&&(psi>60))t_psi=90;
  if((psi<-60)&&(psi>-120)){t_psi=90;transform->reflect=2;}

  if(phi==90)transform->rotate=0;

  if((t_theta==0)&&(t_psi==0))orient=0;
  if((t_theta==90)&&(t_psi==0))orient=1;
  if((t_theta==90)&&(t_psi==90))orient=2;

  /* Axis scale */
  if((date>AXIS_SCALE_START_DATE)&&(date<AXIS_SCALE_STOP_DATE)){
    switch(orient){
    case(0):
      vx*=AXIS_SCALE_X;
      vy*=AXIS_SCALE_Y;
      vz*=AXIS_SCALE_Z;
      break;
    case(1):
      vx*=AXIS_SCALE_Z;
      vy*=AXIS_SCALE_X;
      vz*=AXIS_SCALE_Y;
      break;
    case(2):
      vx*=AXIS_SCALE_Z;
      vy*=AXIS_SCALE_Y;
      vz*=AXIS_SCALE_X;
      break;
    }
  }

  avw_init(header,np/2,nv,ns,vols*ne,DT_COMPLEX);
  avw_set_vox(header,vx,vy,vz);
  avw_set_orient(header,orient);
  if(TR>0.0) avw_set_tr(header,TR);
  else avw_set_tr(header,tr);

  fclose(fp);
  return(0);
}

int volumes_in_varian_header(char *infile1)
{
  struct dsr header;
  struct xfm transform;
  read_varian_header(infile1,&header,&transform);
  return(avw_vdim(&header));
}

unsigned long sizeof_fidfile(char *infile1)
{
  char infile[MAXLEN];
  struct datafilehead fh;
  struct datafilehead_byte fhb;
  struct datablockhead_byte bhb;
  unsigned long bsize;
  FILE *fp;
  int FHEADB,HEADB;
  int reverse,sshort,slong;

  FHEADB=sizeof(fhb);
  HEADB=sizeof(bhb);

  strcpy(infile,infile1);
  strcat(infile,"/fid");

  if((fp=fopen(infile,"rb"))==NULL){
    printf("Cannot open file: %s\n",infile);
    return(1);
  }
  if(!fread(&fhb,sizeof(fhb),1,fp)){
    printf("Read error\n");
    return(1);
  }
  if(find_byte_order(&reverse,&sshort,&slong)){
    printf("Cannot cope with byte size\n");
    return(1);
  }
  convert_filehead(&fhb,&fh,reverse,sshort,slong);

  bsize=(fh.bbytes-HEADB)/fh.ebytes;

  fclose(fp);

  return(bsize*fh.nblocks);
}

int read_varian_chunks(char *infile1,unsigned long chunksize,unsigned int chunks,unsigned int offset,unsigned int step,float* buffer)
{
  char infile[MAXLEN];
  float dcr,dci;
  struct datafilehead fh;
  struct datafilehead_byte fhb;
  struct datablockhead_byte bhb;
  unsigned long bsize,cpb,bpc,j,k,loc;
  int reverse,sshort,slong;
  int FHEADB,BHEADB;
  char *data;
  FILE *fp;

  FHEADB=sizeof(fhb);
  BHEADB=sizeof(bhb);

  strcpy(infile,infile1);
  strcat(infile,"/fid");

  if((fp=fopen(infile,"rb"))==NULL){
    printf("Cannot open file: %s\n",infile);
    return(1);
  }

  /* Read header */
  if(!fread(&fhb,sizeof(fhb),1,fp)){
    printf("Read error\n");
    return(1);
  }

  if(find_byte_order(&reverse,&sshort,&slong)){
    printf("Cannot cope with byte size\n");
    return(1);
  }
  
  convert_filehead(&fhb,&fh,reverse,sshort,slong);

  bsize=(fh.bbytes-BHEADB)/fh.ebytes;
  if((data=(char *)malloc(fh.ebytes*chunksize))==NULL){
    printf("Malloc failed (data)\n");
    return(1);
  }

  cpb=bsize/chunksize;
  if(cpb<1){

    bpc = chunksize/bsize;
    for(j=0;j<chunks;j++){

      loc = offset + j*step;
      fseek(fp,chunksize*loc*fh.ebytes+(loc*bpc*BHEADB)+FHEADB,SEEK_SET);

      for(k=0;k<bpc;k++){

	if(fread(&bhb,sizeof(bhb),1,fp)!=1){
	  printf("Read error\n");
	  return(1);
	}
	get_dc_from_blockhead(bhb,&dcr,&dci);
	
	if(fread(data,fh.ebytes,bsize,fp)!=bsize){
	  printf("Read error\n");
	  return(1);
	}
	
	if(fh.ebytes==2){
	  two_byte_array_to_float(data,&buffer[(j*bpc+k)*bsize],bsize,reverse,sshort,dcr,dci);
	}
	else {
	  four_byte_array_to_float(data,&buffer[(j*bpc+k)*bsize],bsize,reverse,slong,dcr,dci);
	}
      }
    }
    
  } else {

    if(fread(&bhb,sizeof(bhb),1,fp)!=1){
      printf("Read error\n");
      return(1);
    }
    get_dc_from_blockhead(bhb,&dcr,&dci);
    
    for(j=0;j<chunks;j++){
      loc = offset + j*step;
      fseek(fp,(chunksize*loc*fh.ebytes)+((loc/cpb)*BHEADB)+BHEADB+FHEADB,SEEK_SET);
    
      if(fread(data,fh.ebytes,chunksize,fp)!=chunksize){
	printf("Read error\n");
	return(1);
      }
	
      if(fh.ebytes==2){
	two_byte_array_to_float(data,&buffer[j*chunksize],chunksize,reverse,sshort,dcr,dci);
      }
      else {
	four_byte_array_to_float(data,&buffer[j*chunksize],chunksize,reverse,slong,dcr,dci);
      }
    }
  }  
  fclose(fp);
  free(data);

  return(0);
}

int flags_from_procpar(char *infile,struct flags *proc){

  char tmp[MAXLEN];
  int fract_ky,nv;
  float lpe2,pss;
  unsigned int i;
 
  /* Number of echoes */
  if(!proc->me){
    if(!read_procpar(infile,"ne ",tmp)){
      printf("Failed to find number of echos\n");
      return(1);
    }
    proc->me=atoi(tmp);
  }
  if(proc->me>1){
    if(proc->epi){
      if(!read_procpar(infile,"spinecho ",tmp)){
	printf("Failed to find spinecho flag\n");
	return(1);
      }
      if(!strcmp(tmp,"y"))proc->epi=2;
    }
  }

  /* Multi receivers */
  if(read_procpar(infile,"rcvrs ",tmp)){
    proc->rcvrs=0;
    for(i=0;i<strlen(tmp);i++){
      if(tmp[i]=='y')proc->rcvrs++;
    }
  }

  /* Find number of interleaves */
  proc->seg=1;
  if(read_procpar(infile,"num_ints ",tmp)){
    proc->seg=atoi(tmp);
  }
  if(!(proc->seg%2)){
    proc->buo=0;
    proc->buov=0;
  }

  /* Find if epik */
  proc->epik=0;
  if(read_procpar(infile,"epik_kh ",tmp)){
    proc->epik=atoi(tmp);
  }
  
  /* Find if hks */
  fract_ky=0;nv=1;
  if(read_procpar(infile,"fract_ky ",tmp)){
    fract_ky=atoi(tmp);
  }
  if(read_procpar(infile,"nv ",tmp)){
    nv=atoi(tmp);
  }
  proc->fract_ky=0;
  if(fract_ky&&nv){
    if(fract_ky<nv/2){
      proc->fract_ky=fract_ky;
    }
  }
  
  /*epik prescan*/
  if(proc->epik>0){
    if(read_procpar(infile,"prescan_flag ",tmp)){
      if(atoi(tmp)){
	proc->oneseg=1;
      }
    }
  }

  /* Read oversamp */
  if(proc->rover){
    proc->rover=0;
    if(read_procpar(infile,"read_oversamp ",tmp)){
      if((tmp[0]=='y')||(tmp[0]=='Y')){
	proc->rover=1;
      }
    }
  }
  /* Slice oversamp */
  if(proc->sover){
    proc->sover=0;
    if(read_procpar(infile,"slice_oversamp ",tmp)){
      if((tmp[0]=='y')||(tmp[0]=='Y')){
	proc->sover=1;
      }
    }
  }
  
  /* Phase offset for 3D files */
  if(proc->ft==3){
    if(read_procpar(infile,"lpe2 ",tmp)){
      lpe2=atof(tmp);
      proc->phoff=0;
      if(!mean_pss(infile,&pss)){
	proc->phoff = 2.0 * M_PI / lpe2 * pss;
      }
    }
  }  
  
  /* FSE echo train length */
  if(read_procpar(infile,"etl ",tmp)){
    proc->etl=atoi(tmp);
  } else {
    proc->etl=1;
  }

  /* Seqcon parameter */
  read_procpar(infile,"seqcon ",proc->seqcon);

  /* ppe value */
  if(proc->ppe==TOO_LARGE_PPE){
    if(read_procpar(infile,"ppe ",tmp)){
      proc->ppe=atof(tmp)*10;
    } else {
      proc->ppe=0;
    }
  }

  return(0);
}
float sw_readout(char *infile)
{
  char tmp[MAXLEN];
  if(!read_procpar(infile,"sw ",tmp)){
    printf("Failed to find spinecho flag\n");
    return(1);
  }
  return(atof(tmp));
}

int read_procpar(char *filename,char *parameter,char *value)
{
  int i,n=0;
  char proc_string[5000],tmp[5000];
  FILE *fp;

  strcpy(tmp,filename);
  strcat(tmp,"/procpar");

  if((fp=fopen(tmp,"rb"))==NULL){
    printf("Cannot open file: %s",tmp);
    return(0);
  }

  i=strlen(parameter);

  while(fgets(proc_string,5000,fp)!=NULL){
    if(!strncmp(proc_string,parameter,i)){
      fgets(proc_string,5000,fp);
      sscanf(proc_string,"%d %s\n",&n,value);
    }
  }
  fclose(fp);
  
  if(value[0]=='"'){
    strncpy(tmp,&value[1],strlen(value)-2);
    tmp[strlen(value)-2]='\0';
    strcpy(value,tmp);
  }
  return(n);
}
int read_procpar_strings(char *filename,char *parameter,char *value)
{
  int i,ret=0;
  char proc_string[5000],tmp[5000];
  FILE *fp;
  
  strcpy(tmp,filename);
  strcat(tmp,"/procpar");
  
  if((fp=fopen(tmp,"rb"))==NULL){
    printf("Cannot open file: %s",tmp);
    return(0);
  }
  
  i=strlen(parameter);
  
  while(fgets(proc_string,5000,fp)!=NULL){
    if(!strncmp(proc_string,parameter,i)){
      fgets(tmp,5000,fp);
      ret=1;
    }
  }
  fclose(fp);
  if(ret){
    sprintf(value,"%s",strchr(tmp,'\"')+1);
    strrchr(value,'\"')[0]='\0';
  }
  return(ret);
}
int flash_converted(char *filename)
{
  int ret=0;
  char proc_string[5000],tmp[5000];
  FILE *fp;

  strcpy(tmp,filename);
  strcat(tmp,"/procpar");

  if((fp=fopen(tmp,"rb"))==NULL){
    printf("Cannot open file: %s",tmp);
    return(0);
  }

  while(fgets(proc_string,5000,fp)!=NULL){
    if(!strncmp(proc_string,"flash_converted ",16)){
      ret=1;
      break;
    }
  }
  fclose(fp);
  return(ret);
}
int tab_converted(char *filename)
{
  int ret=0;
  char proc_string[5000],tmp[5000];
  FILE *fp;

  strcpy(tmp,filename);
  strcat(tmp,"/procpar");

  if((fp=fopen(tmp,"rb"))==NULL){
    printf("Cannot open file: %s",tmp);
    return(0);
  }

  while(fgets(proc_string,5000,fp)!=NULL){
    if(!strncmp(proc_string,"tab_converted ",14)){
      ret=1;
      break;
    }
  }
  fclose(fp);
  return(ret);
}
long str_to_date(char *datestr)
{
  char month[10];
  int day=0,year=0;
  long date=0;
  sscanf(datestr,"%s %d %d",month,&day,&year);
  if(!strcmp(month,"Jan"))date+=100;
  if(!strcmp(month,"Feb"))date+=200;
  if(!strcmp(month,"Mar"))date+=300;
  if(!strcmp(month,"Apr"))date+=400;
  if(!strcmp(month,"May"))date+=500;
  if(!strcmp(month,"Jun"))date+=600;
  if(!strcmp(month,"Jul"))date+=700;
  if(!strcmp(month,"Aug"))date+=800;
  if(!strcmp(month,"Sep"))date+=900;
  if(!strcmp(month,"Oct"))date+=1000;
  if(!strcmp(month,"Nov"))date+=1100;
  if(!strcmp(month,"Dec"))date+=1200;
  if(date==0){
    printf("Month not recognised %s\n",datestr);
  }
  date+=year*10000;
  date+=day;
  return(date);
}
