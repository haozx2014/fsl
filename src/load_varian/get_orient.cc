#include <cstdio>
#include <cstring>
#include "hfunc.h"

int main(int argc,char *argv[])
{
  char filename[256];
  char orient;
  struct dsr header;

  if(argc<2){
    printf("usage: get_orient filename\n");
    return(1);
  }
  strcpy(filename,argv[1]);
  
  if(avw_read(filename,&header)<0){
    printf("Cannot read header file: %s.hdr\n",filename);
    return(1);
  }
  avw_get_orient(&header,&orient);

  switch(orient){
  case 0:
    printf("transverse\n");
    break;
  case 1:
    printf("coronal\n");
    break;
  case 2:
    printf("sagittal\n");
    break;
  default:
    printf("oblique\n");
  }
  return(0);
}
