/*
  Helper program for load_varian.tcl
  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "maxlen.h"
#include "hfunc.h"

/*

Flags
  -study
  -descrip
  -date
  -scannum
  -patient

*/

int main(int argc,char *argv[])
{
  struct dsr header;
  char filename[MAXLEN];
  char string[80];
  int i;

  if(argc<3){
    printf("usage: set_hdr_info headerfile flag info\n");
    exit(1);
  }

  strcpy(filename,argv[1]);
  if(strstr(filename,".hdr")==NULL)
    strcat(filename,".hdr");
  
  if(avw_read(filename,&header)<0){
    printf("Cannot read header file: %s.hdr\n",filename);
    return(1);
  }

  /* Concatenate arguments to a single string */
  strcpy(string,argv[3]);
  for(i=4;i<argc;i++){
    strcat(string," ");
    strcat(string,argv[i]);
  }

  if(!strcmp(argv[2],"-study"))avw_set_study(&header,string);
  if(!strcmp(argv[2],"-date"))avw_set_date(&header,string);
  if(!strcmp(argv[2],"-descrip"))avw_set_descrip(&header,string);
  if(!strcmp(argv[2],"-scannum"))avw_set_scannum(&header,string);
  if(!strcmp(argv[2],"-patient"))avw_set_patient(&header,string);

  if(avw_write(filename,&header)<0){
    printf("Cannot write header file: %s.hdr\n",filename);
    return(1);
  }

  return(0);
}
