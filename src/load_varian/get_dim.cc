#include <cstdio>
#include <cstring>
#include "data.h"
#include "byteorder.h"
#include "data_byte.h"

int main(int argc,char *argv[])
{
  char filename[256],procfile[256],proc_string[5000];
  long np,nv,ns,nv2,vols;
  long bsize;
  float thk,t_thk,vz;
  struct datafilehead fh;
  struct datafilehead_byte fhb;
  struct datablockhead_byte bhb;
  int reverse,sshort,slong,nD;
  int epik_kh=0,prescan=0,ne=1,nvf,num_ints;
  FILE *fp;

  if(argc<2){
    printf("usage: get_dim filename\n");
    return(1);
  }
  strcpy(filename,argv[1]);
  strcat(filename,"/fid");
  strcpy(procfile,argv[1]);
  strcat(procfile,"/procpar");
  
  if((fp=fopen(filename,"rb"))==NULL){
    printf("Cannot open file: %s\n",filename);
    return(1);
  }
  if(!fread(&fhb,sizeof(fhb),1,fp)){
    printf("Read error");
    return(1);
  }
  fclose(fp);
  if(find_byte_order(&reverse,&sshort,&slong)){
    printf("Cannot cope with byte size\n");
    return(1);
  }
  convert_filehead(&fhb,&fh,reverse,sshort,slong);

  if((fp=fopen(procfile,"rb"))==NULL){
    printf("Cannot open file: %s\n",procfile);
    return(1);
  }
  while(fgets(proc_string,5000,fp)!=NULL){
    if(!strncmp(proc_string,"np ",3)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %ld\n", &np);
    }
    if(!strncmp(proc_string,"nv ",3)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %ld\n", &nv);
    }
    if(!strncmp(proc_string,"nv2 ",4)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %ld\n", &nv2);
    }
    if(!strncmp(proc_string,"ne ",3)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %d\n", &ne);
    }
    if (!strncmp(proc_string, "pss ", 4)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "%ld", &ns);
      if(ns>1){
	sscanf(proc_string,"%ld",&ns);
      }
    }
    if(!strncmp(proc_string,"nD ",3)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %d\n", &nD);
    }
    if (!strncmp(proc_string, "thk ", 4)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %f\n", &thk);
    }
    if (!strncmp(proc_string, "t_thk ", 6)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %f\n", &t_thk);
    }
    if (!strncmp(proc_string, "epik_kh ", 8)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %d\n", &epik_kh);
    }
    if (!strncmp(proc_string, "prescan_flag ", 13)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %d\n", &prescan);
    }
    if (!strncmp(proc_string, "num_ints ", 9)){
      fgets(proc_string, 5000, fp);
      sscanf(proc_string, "1 %d\n", &num_ints);
    }
  }
  fclose(fp);
  if(nv2>1){
    ns=nv2;
  }
  if(nv==0){
    nv=1;
  }
  bsize=(fh.bbytes-sizeof(bhb))/fh.ebytes;

  if(epik_kh){
    nvf = nv+((nv/epik_kh)*(num_ints-1));
    if(prescan){
      nvf=nvf/num_ints;
    }
    vols = (bsize*fh.nblocks)/(np*nvf*ns*ne);
  } else {
    vols = (bsize*fh.nblocks)/(np*nv*ns*ne);
  }

  if(nD==3){
    vz=t_thk;
  }
  else vz=thk;

  printf("%ld\n",np/2);
  printf("%ld\n",nv);
  printf("%ld\n",ns);
  printf("%ld\n",vols);
  printf("%f\n",vz);

  return(0);
}
