/*
  concat_info.c: Concatenates the individual info files from the fid
                 directory to the output avw file

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "concat_info.h"

#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "maxlen.h"

int concat_info_files(char *infile,char *outfile)
{
  
  char inf[MAXLEN],outf[MAXLEN],list[10][15];
  char infileseries[MAXLEN], infilestudy[MAXLEN];
  char line[5000];
  int i;
  FILE *ip;
  FILE *op;

  strcpy(list[0],"comment");
  strcpy(list[1],"patient");
  strcpy(list[2],"study");
  strcpy(list[3],"sequence");
  strcpy(list[4],"options");
  strcpy(list[5],"seq_params");
  strcpy(list[6],"scan_range");
  strcpy(list[7],"matrix");
  strcpy(list[8],"misc");
  strcpy(list[9],"fmri_protocol");


  strcpy(outf,outfile);
  strcat(outf,".info");

  if((op=fopen(outf,"wt"))==NULL){
    printf("Cannot open file: %s\n",outf);
    return(1);
  }

  strcpy(infileseries,infile);
  if(strrchr(infileseries,'.')!=NULL)strrchr(infileseries,'.')[0]='\0';
  if(strrchr(infileseries,'.')!=NULL)strrchr(infileseries,'.')[0]='\0';

  strcpy(infilestudy,infile);
  if(strrchr(infilestudy,'/')!=NULL)strrchr(infilestudy,'/')[0]='\0';
  else strcpy(infilestudy,".");

  for(i=0;i<10;i++){
    sprintf(inf,"%s/info.%s",infile,list[i]);
    if((ip=fopen(inf,"rt"))==NULL){
      sprintf(inf,"%s.1.fid/info.%s",infileseries,list[i]);
      if((ip=fopen(inf,"rt"))==NULL){
	sprintf(inf,"%s/info.%s",infilestudy,list[i]);
	if((ip=fopen(inf,"rt"))==NULL){
	  continue;
	}
      }
    }
    while(fgets(line,5000,ip)!=NULL){
      fprintf(op,"%s",line);
    }
    fprintf(op,"\n");
    fclose(ip);
  }

  fclose(op);
  return(0);
}
