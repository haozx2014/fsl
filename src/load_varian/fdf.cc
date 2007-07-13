#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include "maxlen.h"
#include "fdf.h"
#include "read_varian.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

void fdf_set_vox(struct fdf *header,float vx,float vy,float vz)
{
  header->vx=vx;
  header->vy=vy;
  header->vz=vz;
}
void fdf_set_dim(struct fdf *header,int nx, int ny, int nz, int nv)
{
  header->nx=nx;
  header->ny=ny;
  header->nz=nz;
  header->nv=nv;
}
void fdf_set_type(struct fdf *header,const char *storage,int bytes)
{
  strcpy(header->storage,storage);
  header->bytes=bytes;
}
void fdf_set_maxmin(struct fdf *header,float dmax,float dmin)
{
  header->dmin=dmin;
  header->dmax=dmax;
}
void fdf_set_params(struct fdf *header,char *infile)
{
  char tmp[MAXLEN];
  float phi=0,theta=0,psi=0;
  int nD=2;

  if(read_procpar(infile,"phi ",tmp)){
    header->phi=atof(tmp);
    phi=header->phi*M_PI/180;
  }
  if(read_procpar(infile,"psi ",tmp)){
    header->psi=atof(tmp);
    psi=header->psi*M_PI/180;
  }
  if(read_procpar(infile,"theta ",tmp)){
    header->theta=atof(tmp);
    theta=header->theta*M_PI/180;
  }
  if(read_procpar(infile,"pro ",tmp)){
    header->px=atof(tmp);
  }
  if(read_procpar(infile,"ppe ",tmp)){
    header->py=atof(tmp);
  }
  if(read_procpar(infile,"nD ",tmp)){
    nD=atoi(tmp);
  }
  if(nD==3){
    if(read_procpar(infile,"pss ",tmp)){
      header->pz=atof(tmp);
    }
  }
  if(read_procpar(infile,"seqfil ",tmp)){
    strcpy(header->seqfil,tmp);
  }
  if(read_procpar(infile,"studyid_ ",tmp)){
    strcpy(header->studyid,tmp);
  }
  header->orient[0] = -sin(phi)*sin(psi) - cos(phi)*cos(theta)*cos(psi);
  header->orient[1] = sin(phi)*cos(psi) - cos(phi)*cos(theta)*sin(psi);
  header->orient[2] = sin(theta)*cos(phi);
  header->orient[3] = cos(phi)*sin(psi) - sin(phi)*cos(theta)*cos(psi);
  header->orient[4] = -cos(phi)*cos(psi) - sin(phi)*cos(theta)*sin(psi);
  header->orient[5] = sin(theta)*sin(phi);
  header->orient[6] = cos(psi)*sin(theta);
  header->orient[7] = sin(psi)*sin(theta);
  header->orient[8] = cos(theta);
}
int write_fdf_volume(const void *data,char *fidpath,char *outpath,struct fdf *header)
{
  char tmpstring[MAXLEN];
  short j,k;

  FILE *ofp;

  sprintf(tmpstring,"/bin/rm -rf %s",outpath);
  system(tmpstring);
  sprintf(tmpstring,"/bin/mkdir %s",outpath);
  system(tmpstring);
  sprintf(tmpstring,"/bin/cp %s/procpar %s",fidpath,outpath);
  system(tmpstring);

  /*
  sprintf(tmpstring,"%s/minmax",outpath);
  if((ofp=fopen(tmpstring,"w"))==NULL){
    printf("Failed to open FDF file for writing: %s\n",tmpstring);
    return(1);
  }
  fprintf(ofp,"min %f\nmax %f\n",header->dmin,header->dmax);
  fclose(ofp);
  */
  

  for(j=0;j<header->nv;j++){
    
    for(k=0;k<header->nz;k++){
            
      sprintf(tmpstring,"%s/slice%03dimage%03decho001.fdf",outpath,k+1,j+1);
      
      if((ofp=fopen(tmpstring,"w"))==NULL){
	printf("Failed to open FDF file for writing: %s\n",tmpstring);
	return(1);
      }
      
      /*printf("%dx%dx%d %fx%fx%f\n",header->nx,header->ny,header->nz,header->vx,header->vy,header->vz);*/

      fprintf(ofp,"#!/usr/local/fdf/startup\n");
      fprintf(ofp,"float  rank = 2; \n");
      fprintf(ofp,"char   *spatial_rank = \"2dfov\"; \n");
      fprintf(ofp,"char   *storage = \"%s\"; \n",header->storage);
      fprintf(ofp,"float  bits = %d; \n",header->bytes*8);
      fprintf(ofp,"char   *type = \"absval\"; \n");
      fprintf(ofp,"float  matrix[] = {%d, %d}; \n",header->nx,header->ny);
      fprintf(ofp,"char  *abscissa[] = {\"cm\", \"cm\"}; \n");
      fprintf(ofp,"char  *ordinate[] = { \"intensity\" }; \n");
      fprintf(ofp,"float  span[] = {%f,%f}; \n",header->nx*header->vx/10.0,header->ny*header->vy/10.0);
      fprintf(ofp,"float  origin[] = {0.000000,0.000000}; \n");
      fprintf(ofp,"char  *nucleus[] = {\"H1\",\"H1\"}; \n");
      fprintf(ofp,"float  nucfreq[] = {127.336975,127.336975}; \n");
      fprintf(ofp,"float  location[] = {%f,%f,%f}; \n",header->px,header->py,((k-header->nz/2+0.5)*header->vz)/10+header->pz);
      fprintf(ofp,"float  roi[] = {%f,%f,%f}; \n",header->nx*header->vx/10.0,header->ny*header->vy/10.0,header->vz/10.0);
      fprintf(ofp,"char  *file = \"%s/fid\";\n",fidpath);
      fprintf(ofp,"int    slice_no = %d; \n",k+1);
      fprintf(ofp,"int    slices = %d; \n",header->nz);
      fprintf(ofp,"int    echo_no = 1; \n");
      fprintf(ofp,"int    echoes = 1; \n");
      fprintf(ofp,"float  TE = %f; \n",header->te*1000);
      fprintf(ofp,"float  te = %f; \n",header->te);
      fprintf(ofp,"float  TR = %f; \n",header->tr*1000);
      fprintf(ofp,"float  tr = %f; \n",header->tr);
      fprintf(ofp,"int    ro_size = %d;\n",header->nx);
      fprintf(ofp,"int    pe_size = %d;\n",header->ny);
      fprintf(ofp,"char  *sequence = \"%s\";\n",header->seqfil);
      fprintf(ofp,"char  *studyid = \"%s\";\n",header->studyid);
      fprintf(ofp,"char  *position1 = \"\";\n");
      fprintf(ofp,"char  *position2 = \"\";\n");
      fprintf(ofp,"int    array_index = %d; \n",header->array_index+1);
      fprintf(ofp,"float  array_dim = %f; \n",(float)header->array_dim);
      fprintf(ofp,"float  image = %f; \n",(float)header->array_index+1);
      fprintf(ofp,"int    display_order = %d; \n",(header->array_index*header->nz)+k);
      fprintf(ofp,"float  psi = %f; \n",header->psi);
      fprintf(ofp,"float  phi = %f; \n",header->phi);
      fprintf(ofp,"float  theta = %f; \n",header->theta);
      fprintf(ofp,"float  orientation[] = {");
      fprintf(ofp,"%f,%f,%f,",header->orient[0],header->orient[1],header->orient[2]);
      fprintf(ofp,"%f,%f,%f,",header->orient[3],header->orient[4],header->orient[5]);
      fprintf(ofp,"%f,%f,%f}; \n",header->orient[6],header->orient[7],header->orient[8]);
      fprintf(ofp,"int    checksum=0;");
      fprintf(ofp,"\f\n\n\n\n");
      fprintf(ofp,"%c",0);

      

      if(fwrite((char *)data+(header->nx*header->ny*k*header->bytes),header->bytes,header->nx*header->ny,ofp)
	 !=(unsigned long)(header->nx*header->ny)){
	printf("Write error\n");
	return(1);
      }

      fclose(ofp);
    }
  }
  return(0);
}
