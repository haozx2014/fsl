/*
  istats

  Stuart Clare, FMRIB Physics Group

  Copyright (c) 2000 University of Oxford
   
  Inputs a 4D AVW file and output the mean/stdev/etc.
  within a mask or above a threshold
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


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"avwio/avwio.h"
#include"statsfunc.h"

#define MAXLEN 128

struct cmd_line_opts {
  int mask,out,thold,sl,list,hbins,quiet;
  int mean,std,var,max,med,medd,pc,hist,mode,sum,ent;
  float hmin,hmax;
};

int main(int argc,char *argv[])
{
  short ppl,lpi,ipv,vols,dt;
  short mppl,mlpi,mipv,mvols,mdt;
  int i,j,ims;
  int *hist;
  long volsize,imsize,l;
  short *simage;
  float *fimage;
  char *cimage;
  char *mask;
  struct cmd_line_opts opts;
  float mean,var,med,medd,pc,max,min,mode,sum,ent;
  AVW *ifp,*mfp;
  FILE *ofp;
  void parse_cmd_line(int argc,char **argv,struct cmd_line_opts *opts);

  parse_cmd_line(argc,&argv[0],&opts);
  
  if((ifp=AvwOpen(argv[1],"r"))==NULL){
    printf("Failed to open AVW file: %s\n",argv[1]);
    return(1);
  }

  AvwGetDim(ifp,&ppl,&lpi,&ipv,&vols);
  volsize = ppl*lpi*ipv;
  if(opts.sl){
    imsize = ppl*lpi;
    ims = ipv;
  } else {
    imsize = ppl*lpi*ipv;
    ims = 1;
  }

  AvwGetDataType(ifp,&dt);
  if((dt!=DT_SIGNED_SHORT)&&(dt!=DT_FLOAT)&&(dt!=DT_UNSIGNED_CHAR)){
    printf("Data type not supported\n");
    return(1);
  }

  if(opts.out){
    if((ofp=fopen(argv[opts.out],"w"))==NULL){
      printf("Can't open : %s",argv[opts.out]);
      return(1);
    }
  } else {
    ofp=stdout;
  }

  if((mask=(char *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }
  if((!opts.mask)&&(!opts.thold)){
    for(l=0;l<volsize;l++)mask[l]=1;
  }
  if(opts.mask){
    if((mfp=AvwOpen(argv[opts.mask],"r"))==NULL){
      printf("Failed to open AVW file: %s\n",argv[opts.mask]);
      return(1);
    }
    AvwGetDim(mfp,&mppl,&mlpi,&mipv,&mvols);
    AvwGetDataType(mfp,&mdt);
    
    if((ppl!=mppl)||(lpi!=mlpi)||(mipv!=ipv)){
      printf("Mask not compatible with input file\n");
      return(1);
    }

    switch(mdt){
    case DT_UNSIGNED_CHAR:
      if(AvwReadVolumes(mfp,mask,1)!=1){
	printf("Read error\n");
	return(1);
      }
      break;
    case DT_SIGNED_SHORT:
      if((simage=(short *)malloc(volsize*sizeof(short)))==NULL){
	printf("Malloc failed");
	return(1);
      }
      if(AvwReadVolumes(mfp,simage,1)!=1){
	printf("Read error\n");
	return(1);
      }
      for(l=0;l<volsize;l++){
	if(simage[l]>0)mask[l]=1;
	else mask[l]=0;
      }
      free(simage);     
      break;
    case DT_FLOAT:
      if((fimage=(float *)malloc(volsize*sizeof(float)))==NULL){
	printf("Malloc failed");
	return(1);
      }
      if(AvwReadVolumes(mfp,fimage,1)!=1){
	printf("Read error\n");
	return(1);
      }
      for(l=0;l<volsize;l++){
	if(fimage[l]>0)mask[l]=1;
	else mask[l]=0;
      }      
      free(fimage);
      break;
    default:
      printf("Mask data type not supported\n");
      return(1);
    }
    AvwClose(mfp);
  }
  
  if((fimage=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }

  if(!opts.quiet){
    if(opts.list){
      fprintf(ofp,"Volume\t");
      if(opts.sl)fprintf(ofp,"Slice\t");
    }
    if(opts.mean)fprintf(ofp,"Mean\t");
    if(opts.var)fprintf(ofp,"Variance\t");
    if(opts.std)fprintf(ofp,"StdDev\t");
    if(opts.max)fprintf(ofp,"Maximum\t");
    if(opts.sum)fprintf(ofp,"Sum\t");
    if(opts.med)fprintf(ofp,"Median\t");
    if(opts.medd)fprintf(ofp,"MedDev\t");
    if(opts.pc)fprintf(ofp,"%d-Percentile\t",opts.pc);
    if(opts.mode)fprintf(ofp,"Mode\t");
    if(opts.ent)fprintf(ofp,"Entropy\t");
    fprintf(ofp,"\n");
  }

  for(i=0;i<vols;i++){
    if(dt==DT_UNSIGNED_CHAR){
      if((cimage=(char *)malloc(volsize))==NULL){
	printf("Malloc failed");
	return(1);
      }
      if(AvwReadVolumes(ifp,cimage,1)!=1){
	printf("Read error.\n");
      }
      for(l=0;l<volsize;l++){
	fimage[l]=(float)cimage[l];
	if(opts.thold){
	  if(cimage[l]>opts.thold)mask[l]=1;
	  else mask[l]=0;
	}
      }
      free(cimage);
    } 
    if(dt==DT_SIGNED_SHORT){
      if((simage=(short *)malloc(volsize*sizeof(short)))==NULL){
	printf("Malloc failed");
	return(1);
      }
      if(AvwReadVolumes(ifp,simage,1)!=1){
	printf("Read error.\n");
      }
      for(l=0;l<volsize;l++){
	fimage[l]=(float)simage[l];
	if(opts.thold){
	  if(simage[l]>opts.thold)mask[l]=1;
	  else mask[l]=0;
	}
      }
      free(simage);      
    } 
    if(dt==DT_FLOAT){
      if(AvwReadVolumes(ifp,fimage,1)!=1){
	printf("Read error.\n");
      }
      if(opts.thold){
	for(l=0;l<volsize;l++){
	  if(fimage[l]>opts.thold)mask[l]=1;
	  else mask[l]=0;
	}
      }
    }
    
    for(j=0;j<ims;j++){
      
      if(opts.list){
	fprintf(ofp,"%d\t",i+1);
	if(opts.sl)fprintf(ofp,"%d\t",j+1);
      }

      /* Calculate mean */
      if(opts.mean||opts.var||opts.std){
	mean=calc_mean(&fimage[j*imsize],imsize,&mask[j*imsize]);
	if(opts.mean)fprintf(ofp,"%.2f\t",mean);

	if(opts.var||opts.std){
	  var=calc_var(&fimage[j*imsize],imsize,&mask[j*imsize],mean);
	  if(opts.var)fprintf(ofp,"%.2f\t",var);
	  if(opts.std)fprintf(ofp,"%.2f\t",sqrt(var));
	}
      }

      /* Calculate maximum */
      if(opts.max){
	max=calc_max(&fimage[j*imsize],imsize,&mask[j*imsize]);
	fprintf(ofp,"%.2f\t",max);
      }

      if(opts.sum){
	sum=calc_sum(&fimage[j*imsize],imsize,&mask[j*imsize]);
	fprintf(ofp,"%.2f\t",sum);
      }
      
      /* Calculate median */
      if(opts.med||opts.medd){
	med=calc_med(&fimage[j*imsize],imsize,&mask[j*imsize]);
	if(opts.med)fprintf(ofp,"%.2f\t",med);
	if(opts.medd){
	  medd=sqrt(calc_medv(&fimage[j*imsize],imsize,&mask[j*imsize],med));
	  fprintf(ofp,"%.2f\t",medd);
	}
      }

      /* Calculate N'th percentile */
      if(opts.pc){
	pc=calc_pc(&fimage[j*imsize],imsize,&mask[j*imsize],opts.pc);
	fprintf(ofp,"%.2f\t",pc);
      }

      /* Calculate mode */
      if(opts.mode){
	max=calc_max(&fimage[j*imsize],imsize,&mask[j*imsize]);
	min=calc_min(&fimage[j*imsize],imsize,&mask[j*imsize]);
	mode=calc_mode(&fimage[j*imsize],imsize,&mask[j*imsize],max,min);
	fprintf(ofp,"%.2f\t",mode);
      }
      /* Calculate entropy */
      if(opts.ent){
	ent=calc_entropy(&fimage[j*imsize],imsize,&mask[j*imsize]);
	fprintf(ofp,"%.6f\t",ent);
      }

      if(opts.hist){
	hist = (int *)malloc(opts.hbins*sizeof(int));
	calc_hist(&fimage[j*imsize],imsize,&mask[j*imsize],
		  opts.hmin,opts.hmax,opts.hbins,hist);
	for(i=0;i<opts.hbins;i++)fprintf(ofp,"%d\t",hist[i]);
	free(hist);
      }

      fprintf(ofp,"\n");
    }
  }
  return(0);
}

void parse_cmd_line(int argc,char **argv,struct cmd_line_opts *opts)
{
  int i;
  void usage(void);

  if(argc>1){
    if((!strcmp(argv[1],"-h"))||(!strcmp(argv[1],"-help")))
      usage();
  }
  if(argc<2)usage();
  if(argv[1][0]=='-'){
    printf("Invalid input file\n\n");
    usage();
  }

  opts->out=0;
  opts->mask=0;
  opts->thold=0;
  opts->sl=0;
  opts->list=0;
  opts->mean=0;
  opts->std=0;
  opts->var=0;
  opts->medd=0;
  opts->med=0;
  opts->pc=0;
  opts->max=0;
  opts->hist=0;
  opts->mode=0;
  opts->sum=0;
  opts->quiet=0;
  opts->ent=0;

  for(i=2;i<argc;i++){
    if(!strcmp(argv[i],"-o")){opts->out=i+1;continue;}
    if(!strcmp(argv[i],"-m")){opts->mask=i+1;continue;}
    if(!strcmp(argv[i],"-t")){opts->thold=atoi(argv[i+1]);continue;}
    if(!strcmp(argv[i],"-2d")){opts->sl=1;continue;}
    if(!strcmp(argv[i],"-l")){opts->list=1;continue;}
    if(!strcmp(argv[i],"-mean")){opts->mean=1;continue;}
    if(!strcmp(argv[i],"-std")){opts->std=1;continue;}
    if(!strcmp(argv[i],"-var")){opts->var=1;continue;}
    if(!strcmp(argv[i],"-med")){opts->med=1;continue;}
    if(!strcmp(argv[i],"-medd")){opts->medd=1;continue;}
    if(!strcmp(argv[i],"-sum")){opts->sum=1;continue;}
    if(!strcmp(argv[i],"-pc")){opts->pc=atoi(argv[i+1]);continue;}
    if(!strcmp(argv[i],"-hist")){
      opts->hist=1;
      i++;
      opts->hmin=atof(argv[i]);
      i++;
      opts->hmax=atof(argv[i]);
      i++;
      opts->hbins=atoi(argv[i]);
      continue;
    }
    if(!strcmp(argv[i],"-entropy")){opts->ent=1;continue;}
    if(!strcmp(argv[i],"-mode")){opts->mode=1;continue;}
    if(!strcmp(argv[i],"-max")){opts->max=1;continue;}
    if(!strcmp(argv[i],"-q")){opts->quiet=1;continue;}
    if(argv[i][0]=='-'){
      printf("flag %s not recognised.  Use -h option for help\n",argv[i]);
      exit(1);
    }
  }
}

void usage(void)
{
  printf("usage: istats filename [opts]\n");
  printf(" options\n");
  printf("        -o filename         Output to file\n");
  printf("        -m filename         Use mask\n");
  printf("        -t thold            Use threshold\n");
  printf("        -2d                 Treat as 2D\n");
  printf("        -l                  Use list lables\n");
  printf("        -mean               Output mean\n");
  printf("        -std                Output standard deviation\n");
  printf("        -var                Output variance\n");
  printf("        -max                Output maximum\n");
  printf("        -med                Output median\n");
  printf("        -medd               Output median deviation\n");
  printf("        -sum                Output sum\n");
  printf("        -pc N               Output N'th percentile\n");
  printf("        -mode               Output mode\n");
  printf("        -entropy            Output entropy\n");
  printf("        -hist min max bins  Output histogram\n");
  printf("        -q                  No column headings\n");
  exit(1);
}
