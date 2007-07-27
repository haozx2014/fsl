/*
  LoadVarian: Transforms time data from the Varian scanner to images
  
  Copyright (c) Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Version 1: July 1998

  Version 1.1: Nov 1998
  Sorted out image orientation for medx
  New median baseline correction method

  Version 1.2: Feb 1999
  RRI stuff removed
  ft routines rewritten
  new max and min calculation
  pss reorder and scale slices
  raw data processing

  Version 1.21: Mar 1999
  Polynomial constrian of phase correction
  Hu phase correction
  Verbose option
  Main split into function process and shorter main

  Version 1.22: Mar 1999
  Polynomial fit added for Buonocore phase method
  FT rows only

  Version 1.23: Mar 1999
  Zeropad to any 2^n matrix square size

  Version 1.24: Mar 1999
  User flag for scale or max

  Version 1.25: Apr 1999
  Load whole series option
  Copy data option

  Version 1.26: May 1999
  Override points/lines/slices/volumes
  Image segment option

  Version 1.27: June 1999
  LLEPI / multiecho reorder added
  Volume reorder added
  Navigator echo processing added

  Version 1.28 July 1999
  Looppss added
  Uses new header functions
  Buonocore method can use a separate image as ref and search for central line
  HKS processing added (by Patricia Figueiredo)
  Reslice option added
  Phase rotate option added
  Int nav rotate added

  Version 1.28a Jan 2000
  Reslice option bug fixed
  Buonocore vol by vol readded
  Added me flag for odd types if multiple spin echo
  Added kspace mask border
  Added smod

  Version 1.28b Feb 2000
  HKS processing fixed
  nrutil replaced by nrvector and nrmatrix

  Version 1.28c Mar 2000
  HKS processing fixed properly
  Image border mask added
  Scale slice fixed

  Version 1.29 May 2000
  Field map reordering added
  Read from fid.orig enabled
  Reverse process added
  Fixed buonocore ref scan bug

  Version 1.29b July 2000
  Spin echo array reorder added
  Reverse load added
  Fixed 'rotate' bug
   
  Version 1.29c August 2000
  Read raw data from AVW added

  Version 1.30 September 2000
  Margosian and POCS hks added
  Negative rotate flag added

  Version 1.31 October 2000
  Fixed find centre kspace

  Version 1.32 November 2000
  3D FT added with new 2dft flag
  Output pss to file added
  Shift matrix (tep correction) added
  Slight changes to main
  Z zero fill added

  Version 1.33 December 2000
  FFTW processing added
  No read-oversamp added
  Direct read of data from fid added

  Version 2.0 January 2001
  Rewrite of main and process
 
  Version 2.01 February 2001
  Load of single slice
  All data now read directly from fid file
  Multiecho Spinecho EPI process broken

  Version 2.02 March 2001
  Turned off zerofill
  Optimum tep estimator/groa calc
  Nosave option added
  Sear fixed

  Version 2.03 June 2001
  ctep added
  output slice order to avw header
  output approx number of dropouts
  process at lowres
  corrected 3D processing

  Version 2.04 July 2001
  Phrot file added
  Auto process added
  Float output option

  Version 2.05 Aug 2001
  Fermi filter added
  
  Version 2.06 Oct 2001
  ctep depricated
  Baseline correction with fermi filter fixed
  Fixed etl 'bug'

  Version 2.07 Oct 2001
  Change HKS to detect how much to fill
  Fixed phase correction for nseg>1
  
  Version 2.08 Dec 2001
  Added entropy ghost minimisation
  Added oneseg option
  
  Version 2.09 Dec 2001
  Added avvol averaging
  Reintroduced simple nav
  IEPI ghost correction

  Version 2.10 Apr 2002
  Fixed rotations properly
  Resample images added

  Version 2.11 Nov 2002
  Added EPIK processing
  Added tabc reordering
  Added multi receiver options
  Added average all option

  Version 2.12 Mar 2003
  Fixed for new code directory structure

  Version 2.13 Aug 2003
  Output info files to File.info

  Version 2.14 Nov 2003
  Added phase encode offset (ppe) code 

  Version 2.15 Mar 2004
  Added gz support
  
  Version 2.16 Jun 2004
  Added multi-coil averaging.

  Version 2.17 Jun 2004
  Fixed entropy iEPI code

  Version 2.18 Nov 2004
  New kspace filter option

  Version 2.19 Sep 2005
  Make self-ref default
  Add fast bias field correction
  Fixed bug with offset and multicoil
  Bug fix to hks process of 3D data
  Slice oversamp cropping added

  Version 2.20 Nov 2005
  FSE processing added

  Version 2.21 Dec 2005
  4 channel coil buo scan correction

  Version 2.22 Jan 2006
  Varian FDF output added

  CURRENT BUGS:
  Nav dosn't work
  Recode IEPI for speed
  HKS processing not optimised for speed

  OTHER NOTES:
  The following flags are set automatically from procpar
  flags.me
  flags.rcvrs
  flags.seg
  flags.epik
  flags.rover
  flags.sover
  flags.phoff
  flags.ppe;
  flags.etl;
*/

#define VERSION "2.22"

/*
  Process data using the following flag conventions

  -h            view help
  -help         view help
  -ft2d         2Dfft
  -ft3d         3Dfft
  -ms           multislice reorder
  -epi          epi reorder
  -mod          output modulus
  -phs          output phase
  -re           output real
  -im           output imaginary
  -16           output shorts
  -bl           baseline correct
  -ref fn       ref scan name
  -buo          Buonocore phase correction
  -buov         Buonocore vol-by-vol phase correction
  -con          constrain phase correction
  -poly         constrian with 5th order polynomial
  -rot          rotate as necessary
  -negrot       negative rotate
  -start n      start from
  -num n        number to process
  -scsl [fn]    scale slices
  -ftr          1Dfft rows only
  -zero n       zeropad to n
  -zzero n      zero fill in the z direction
  -scale n      scale factor of n
  -max n        scale 99% of dist to n
  -lpss n       loop pss reorder
  -hks          half k-space reconstruction
  -pocs         POCS half k-space reconstruction
  -marg         Margosian half k-space reconstruction
  -resl         reslice to axial plane
  -phrot n      phase rotate resetting every n volumes
  -segx low num segment image in x dimension
  -segy low num segment image in x dimension
  -segz low num segment image in x dimension
  -ph0 n        zeroth order alternate echo phase correction
  -ph1 n        first order alternate echo phase correction
  -ne n         multiple echo
  -kmb n        kspace mask border (width n)
  -imb n        image mask border (width n)
  -fmap         fieldmap reordering
  -sear n       spin echo array reorder with step n
  -revproc      process volumes in reverse
  -revload      load volumes in reverse
  -opss         output pss to file
  -fftw         use the fftw functions rather than numerical rec (default)
  -nrft         use the numerical rec functions rather than fftw
  -rover        check for read oversamp correction (default)
  -noro         no read oversamp correction
  -noso         no slice oversamp correction
  -overx n      override x dimension to header file
  -overy n      override y dimension to header file
  -overz n      override z dimension to header file
  -overv n      override v dimension to header file
  -slice n-n    Load slice or slice range
  -otep         Output tep correction
  -ogroa        Output groa value
  -ogrof        Output grof value
  -nosave       As second argument - do not output image
  -ores n       Downsample by factor of n to low resolution matrix
  -odrop        Output number of dropouts
  -phrotf file  Phase rotate based on reference file
  -auto         Automatically determine best processing options
  -verbose      print debugging information
  -short        Output data as singed 16-bit integers
  -float        Output data as floats
  -fermi        Fermi filter
  -egr          Entropy ghost reduction
  -egrv         Entropy ghost reduction (vol-by-vol)
  -ecf          Entropy cost function
  -eiepi        Entropy IEPI phase correction
  -oneseg       Process only the first segment of a segmented EPI data set
  -avvol        Average all volumes before processing
  -lim n1 n2    Limits in self ref correction
  -ctep         Correct tep error
  -nav n/m/p/mp navigator echo correction
  -resamp p l s Resample to p x l x s (can be number or p,l or s)
  -tabc [file]  Perform PE table reordering (filename optional)
  -avall        Output the mean of all volumes
  -ssqall       Output the sum of the squares of all volumes
  -info         Output the info.comment etc. to File.info file
  -ppe          Offset in the phase encode direction
  -ppe2         Offset in the 2nd (slice) phase encode direction
  -gz           Output gzipped data
  -multicoil    SSq average multiple coil data
  -kpsf a b     PSF filter in the PE direction
  -3dstruct     Process using dedicated 3dstruct processing
  -fdf          Output in fdf format

  Depricated in v2.0 or later
  -intnav       phase rotate using internal navigator
  -ft           2Dfft (for backward compatibility)
  -hu           Hu phase correction
  -raw x y z v  process raw data
  -series       process all images in the same series
  -cp fn        copy output file to this filename
  -smod n       Signed modulus of the data
  -phmap fn     Apply phase map (complex AVW file)
  -orig         Reads fid.orig rather than fid
  -shift n      Shift data array (tep correction)
  -avwin        raw data in AVW file
  -hksref file  half k-space ref file
  -nodr         no direct read of data
  -16           Replaced with -short
*/


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "flags.h"
#include "process.h"
#include "maxlen.h"
#include "read_varian.h"
#include "navigator.h"

int verbose;
int buo_limit1,buo_limit2;

int main(int argc,char *argv[])
{
  char infile[MAXLEN],outfile[MAXLEN];
  struct flags proc;
  int fid_exists(char *infile);
  void usage(int,char**,struct flags*);
  void show_flags(struct flags flgs);

  /* Interpret command line flags */
  usage(argc,&argv[0],&proc);
  
  /* Find filenames */
  sprintf(infile,"%s",argv[1]);
  if(argv[2][0]=='-')sprintf(outfile,"%s",argv[1]);
  else sprintf(outfile,"%s",argv[2]);

  if(!fid_exists(infile)){
    strcat(infile,".fid");
    if(!fid_exists(infile)){
      printf("Cannot open fid file %s\n",infile);
      return(1);
    }
  }

  /* Set any other flags from reading the procpar file */
  flags_from_procpar(infile,&proc);

  /* fdf flags */
  if(proc.fdf){
    proc.nosave=1;
  }

  if(verbose)show_flags(proc);

  /* Do the processing */
  if(process(infile,outfile,&proc))return(1);

  return(0);
}
void usage(int argc,char** argv,struct flags *flgs)
{
  int i;
  void help_message(void);
  int fid_exists(char *infile);
  void depricated(char *arg);
  void auto_determine(char*,struct flags*);
  void parse_arguments(int low,int argc,char** argv,struct flags *flgs);
  void zero_all_flags(struct flags *flgs);

  if(argc>1){
    if((!strcmp(argv[1],"-h"))||(!strcmp(argv[1],"-help")))
      help_message();
    if((!strcmp(argv[1],"-v"))||(!strcmp(argv[1],"-version"))){
      printf("%s\n",VERSION);
      exit(0);
    }
  }
  if(argc<3){
    printf("Usage: load_varian infile outfile [processing options]\n");
    exit(1);
  }

  zero_all_flags(flgs);

  /* Filename check */
  if(argv[1][0]=='-'){
    printf("Invalid input file\n");
    exit(0);
  }
  if(argv[2][0]=='-'){
    if(!strcmp(argv[2],"-nosave"))flgs->nosave=1;
    else {
      if (strcmp(argv[2],"-auto")) {
	printf("Invalid output file\n");
	exit(0);
      }
    }
  }
  
  for(i=2;i<argc;i++)
    if(!strcmp(argv[i],"-auto"))auto_determine(argv[1],flgs);

  parse_arguments(3,argc,argv,flgs);
  
  /* Check the integrity of the options selected */
  if(!flgs->epi){
    flgs->ref=0;
    flgs->buo=0;
    flgs->mph=0;
  }
  if(!flgs->out)flgs->sint=0;
  if(flgs->resl)flgs->rot=1;
  if(flgs->buov)flgs->buo=0;
  if(flgs->egrv)flgs->egr=0;
  if(flgs->hks==4)if(flgs->out)flgs->out=OUT_REAL;

  if((flgs->otep)&&(!flgs->out))flgs->out=OUT_MOD;
  if((flgs->ogroa)&&(!flgs->out))flgs->out=OUT_MOD;
  if((flgs->ogrof)&&(!flgs->out))flgs->out=OUT_MOD;
  
  if(flgs->hislice-flgs->loslice){
    if(flgs->ms){
      printf("Slice range from multislice data is not allowed\n");
      exit(1);
    }
    if((flgs->fm)||(flgs->sear)||(flgs->lpss)){
      printf("Single slice load not allowed with special reordering\n");
      exit(1);
    } 
  }
}

void zero_all_flags(struct flags *flgs)
{
  int i;
  flgs->dim=0;
  flgs->ft=0;
  flgs->ms=0;
  flgs->epi=0;
  flgs->epik=0;
  flgs->out=0;
  flgs->seg=0;
  flgs->sint=0;
  flgs->bl=0;
  flgs->ref=0;
  flgs->rot=0;
  flgs->buo=0;
  flgs->buov=0;
  flgs->con=0;
  flgs->start=0;
  flgs->num=0;
  flgs->scsl=0;
  flgs->pss=0;
  flgs->zero=0;
  flgs->zzero=0;
  flgs->max=0.0;
  flgs->scale=0.0;
  flgs->lowx=0;
  flgs->lowy=0;
  flgs->lowz=0;
  flgs->numx=0;
  flgs->numy=0;
  flgs->numz=0;
  flgs->vr=0;
  flgs->me=0;
  flgs->lpss=0;
  flgs->hks=0;
  flgs->resl=0;
  flgs->mph=0;
  flgs->ph0=0.0;
  flgs->ph1=0.0;
  flgs->phrot=0;
  flgs->kmb=0;
  flgs->imb=0;
  flgs->fm=0;
  flgs->revproc=0;
  flgs->revload=0;
  flgs->sear=0;
  flgs->opss=0;
  flgs->fftw=1;
  flgs->rover=1;
  flgs->sover=1;
  flgs->overx=0;
  flgs->overy=0;
  flgs->overz=0;
  flgs->overv=0;
  flgs->loslice=-1;
  flgs->hislice=-1;
  flgs->slices=0;
  flgs->nav=0;
  flgs->otep=0;
  flgs->ctep=0;
  flgs->ogroa=0;
  flgs->ogrof=0;
  flgs->odrop=0;
  flgs->ds=0;
  flgs->phrotf=0;
  flgs->nosave=0;  
  flgs->fermi=0;
  flgs->egr=0;
  flgs->egrv=0;
  flgs->ecf=0;
  flgs->eiepi=0;
  flgs->oneseg=0;
  flgs->avvol=0;
  flgs->tabc=0;
  flgs->avall=0;
  flgs->rcvrs=1;
  flgs->phoff=0;
  flgs->info=0;
  flgs->ppe=TOO_LARGE_PPE;
  flgs->ppe2=0;
  flgs->gz=0;
  flgs->multicoil=0;
  flgs->kpsf=0;
  flgs->fdf=0;
  flgs->fract_ky=0;
  flgs->cusfilt=0;
  flgs->etl=0;
  buo_limit1=buo_limit2=0;
  for(i=0;i<3;i++)flgs->resamp[i]=0;
  strcpy(flgs->petable,"");
}

void parse_arguments(int low,int argc,char** argv,struct flags *flgs)
{
  int i;
  void help_message(void);
  int fid_exists(char *infile);
  void depricated(char *arg);

  for(i=low;i<argc;i++){

    /* Flags with no arguments */
    if(!strcmp(argv[i],"-ms")){flgs->ms=1;continue;}
    if(!strcmp(argv[i],"-epi")){flgs->epi=1;continue;}
    if(!strcmp(argv[i],"-mod")){flgs->out=OUT_MOD;continue;}
    if(!strcmp(argv[i],"-phs")){flgs->out=OUT_PHS;continue;}
    if(!strcmp(argv[i],"-re")){flgs->out=OUT_REAL;continue;}
    if(!strcmp(argv[i],"-im")){flgs->out=OUT_IMAG;continue;}
    if(!strcmp(argv[i],"-16")){flgs->sint=1;continue;}
    if(!strcmp(argv[i],"-short")){flgs->sint=1;continue;}
    if(!strcmp(argv[i],"-float")){flgs->sint=0;continue;}
    if(!strcmp(argv[i],"-bl")){flgs->bl=1;continue;}
    if(!strcmp(argv[i],"-rot")){flgs->rot=1;continue;}
    if(!strcmp(argv[i],"-negrot")){flgs->rot=-1;continue;}
    if(!strcmp(argv[i],"-buo")){flgs->buo=1;continue;}
    if(!strcmp(argv[i],"-buov")){flgs->buov=1;continue;}
    if(!strcmp(argv[i],"-con"))if(!flgs->con){flgs->con=1;continue;}
    if(!strcmp(argv[i],"-pss")){flgs->pss=1;continue;}
    if(!strcmp(argv[i],"-poly")){flgs->con=5;continue;}
    if(!strcmp(argv[i],"-ftr")){flgs->ft=1;continue;}
    if(!strcmp(argv[i],"-hks")){flgs->hks=2;continue;}
    if(!strcmp(argv[i],"-hkspad")){flgs->hks=1;continue;}
    if(!strcmp(argv[i],"-pocs")){flgs->hks=3;continue;}
    if(!strcmp(argv[i],"-marg")){flgs->hks=4;continue;}
    if(!strcmp(argv[i],"-resl")){flgs->resl=1;continue;}
    if(!strcmp(argv[i],"-fmap")){flgs->fm=1;continue;}
    if(!strcmp(argv[i],"-revproc")){flgs->revproc=1;continue;}
    if(!strcmp(argv[i],"-revload")){flgs->revload=1;continue;}
    if(!strcmp(argv[i],"-ft2d")){flgs->ft=2;continue;}
    if(!strcmp(argv[i],"-ft3d")){flgs->ft=3;continue;}
    if(!strcmp(argv[i],"-opss")){flgs->opss=1;continue;}
    if(!strcmp(argv[i],"-fftw")){flgs->fftw=1;continue;}
    if(!strcmp(argv[i],"-nrft")){flgs->fftw=0;continue;}
    if(!strcmp(argv[i],"-noro")){flgs->rover=0;continue;}
    if(!strcmp(argv[i],"-noso")){flgs->sover=0;continue;}
    if(!strcmp(argv[i],"-otep")){flgs->otep=1;continue;}
    if(!strcmp(argv[i],"-ogroa")){flgs->ogroa=1;continue;}
    if(!strcmp(argv[i],"-ogrof")){flgs->ogrof=1;continue;}
    if(!strcmp(argv[i],"-odrop")){flgs->odrop=1;continue;}
    if(!strcmp(argv[i],"-verbose")){verbose=1;continue;}
    if(!strcmp(argv[i],"-short")){flgs->sint=1;continue;}
    if(!strcmp(argv[i],"-float")){flgs->sint=0;continue;}
    if(!strcmp(argv[i],"-fermi")){flgs->fermi=1;continue;}
    if(!strcmp(argv[i],"-noresl")){flgs->resl=0;continue;}
    if(!strcmp(argv[i],"-egr")){flgs->egr=1;continue;}
    if(!strcmp(argv[i],"-eiepi")){flgs->eiepi=1;continue;}
    if(!strcmp(argv[i],"-egrv")){flgs->egrv=1;continue;}
    if(!strcmp(argv[i],"-ecf")){flgs->ecf=1;continue;}
    if(!strcmp(argv[i],"-oneseg")){flgs->oneseg=1;continue;}
    if(!strcmp(argv[i],"-avvol")){flgs->avvol=1;continue;}    
    if(!strcmp(argv[i],"-ctep")){flgs->ctep=1;continue;}
    if(!strcmp(argv[i],"-avall")){flgs->avall=1;continue;}
    if(!strcmp(argv[i],"-ssqall")){flgs->avall=2;continue;}
    if(!strcmp(argv[i],"-info")){flgs->info=1;continue;}
    if(!strcmp(argv[i],"-gz")){flgs->gz=1;continue;}
    if(!strcmp(argv[i],"-multicoil")){flgs->multicoil=1;continue;}
    if(!strcmp(argv[i],"-fdf")){flgs->fdf=1;continue;}

    if(!strcmp(argv[i],"-auto")){continue;}
    if((!strcmp(argv[i],"-h"))||(!strcmp(argv[i],"-help")))help_message();
    
    /* Flags with integer arguments */
    if(!strcmp(argv[i],"-start")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-start option requires an argument\n"); exit(1); }
      else flgs->start=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-num")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-num option requires an argument\n"); exit(1); }
      else flgs->num=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-zero")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-zero option requires an argument\n"); exit(1); }
      else flgs->zero=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-zzero")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-zzero option requires an argument\n"); exit(1); }
      else flgs->zzero=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-vr")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-vr option requires an argument\n"); exit(1); }
      else flgs->vr=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-sear")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-sear option requires an argument\n"); exit(1); }
      else flgs->sear=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-lpss")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-lpss option requires an argument\n"); exit(1); }
      else flgs->lpss=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-phrot")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-phrot option requires an argument\n"); exit(1); }
      else flgs->phrot=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-ne")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-ne option requires an argument\n"); exit(1); }
      else flgs->me=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-kmb")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-kmb option requires an argument\n"); exit(1); }
      else flgs->kmb=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-imb")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-imb option requires an argument\n"); exit(1); }
      else flgs->imb=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-overx")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-overx option requires an argument\n"); exit(1); }
      else flgs->overx=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-overy")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-overy option requires an argument\n"); exit(1); }
      else flgs->overy=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-overz")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-overz option requires an argument\n"); exit(1); }
      else flgs->overz=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-overv")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-overv option requires an argument\n"); exit(1); }
      else flgs->overv=atoi(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-ores")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-ores option requires an argument\n"); exit(1); }
      else flgs->ds=atoi(argv[i+1]);
      i++;continue;
    }

    /* Flags with float arguments */
    if(!strcmp(argv[i],"-max")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-max option requires an argument\n"); exit(1); }
      else flgs->max=atof(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-scale")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-scale option requires an argument\n"); exit(1); }
      else flgs->scale=atof(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-ph0")){
      flgs->ph0=atof(argv[i+1]);
      flgs->mph=1;
      i++;continue;
    }
    if(!strcmp(argv[i],"-ph1")){
      flgs->ph1=atof(argv[i+1]);
      flgs->mph=1;
      i++;continue;
    }
    if(!strcmp(argv[i],"-ppe")){
      flgs->ppe=atof(argv[i+1]);
      i++;continue;
    }
    if(!strcmp(argv[i],"-ppe2")){
      flgs->ppe2=atof(argv[i+1]);
      i++;continue;
    }

    /* Flags with two arguments */
    if(!strcmp(argv[i],"-segx")){
      if((i+2==argc)||(argv[i+1][0]=='-')||(argv[i+2][0]=='-')){
	printf("-segx option requires two arguments\n");
	exit(1);
      }
      else {
	flgs->lowx=atoi(argv[i+1]);
	flgs->numx=atoi(argv[i+2]);
      }
      i+=2;continue;
    }
    if(!strcmp(argv[i],"-segy")){
      if((i+2==argc)||(argv[i+1][0]=='-')||(argv[i+2][0]=='-')){
	printf("-segy option requires two arguments\n");
	exit(1);
      }
      else {
	flgs->lowy=atoi(argv[i+1]);
	flgs->numy=atoi(argv[i+2]);
      }
      i+=2;continue;
    }
    if(!strcmp(argv[i],"-segz")){
      if((i+2==argc)||(argv[i+1][0]=='-')||(argv[i+2][0]=='-')){
	printf("-segz option requires two arguments\n");
	exit(1);
      }
      else {
	flgs->lowz=atoi(argv[i+1]);
	flgs->numz=atoi(argv[i+2]);
      }
      i+=2;continue;
    }
    if(!strcmp(argv[i],"-lim")){
      if((i+2==argc)||(argv[i+1][0]=='-')||(argv[i+2][0]=='-')){
	printf("-lim option requires two arguments\n");
	exit(1);
      }
      else {
	buo_limit1=atoi(argv[i+1]);
	buo_limit2=atoi(argv[i+2]);
      }
      i+=2;continue;
    }

    if(!strcmp(argv[i],"-kpsf")){
      if((i+2==argc)||(argv[i+1][0]=='-')||(argv[i+2][0]=='-')){
	printf("-kpsf option requires two arguments\n");
	exit(1);
      }
      else {
	flgs->kpsf=1;
	flgs->kpsf_a=atof(argv[i+1]);
	flgs->kpsf_b=atof(argv[i+2]);
      }
      i+=2;continue;
    }

    /* Flags with three arguments */
    if(!strcmp(argv[i],"-resamp")){
      if((i+3==argc)||(argv[i+1][0]=='-')
	 ||(argv[i+2][0]=='-')||(argv[i+3][0]=='-')){
	printf("-resamp option requires three arguments\n");
	exit(1);
      }
      else {
	if(strcmp(argv[i+1],"p"))flgs->resamp[0]=atoi(argv[i+1]);
	if(strcmp(argv[i+1],"l"))flgs->resamp[1]=atoi(argv[i+2]);
	if(strcmp(argv[i+1],"s"))flgs->resamp[2]=atoi(argv[i+3]);
      }
      i+=3;continue;
    }

   /* Flags with string arguments */
    if(!strcmp(argv[i],"-ref")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-ref option requires an argument\n"); exit(1); }
      else {
	flgs->ref=1;
	strcpy(flgs->reffile,argv[i+1]);
	if(!fid_exists(flgs->reffile)){
	  strcat(flgs->reffile,".fid");
	  if(!fid_exists(flgs->reffile)){
	    printf("Cannot open reference fid file %s\n",flgs->reffile);
	    exit(1);
	  }
	}
      }
      i++;continue;
    }
    if(!strcmp(argv[i],"-phrotf")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-phrotf option requires an argument\n"); exit(1); }
      else {
	flgs->phrotf=1;
	strcpy(flgs->phrotfile,argv[i+1]);
	if(!fid_exists(flgs->phrotfile)){
	  strcat(flgs->phrotfile,".fid");
	  if(!fid_exists(flgs->phrotfile)){
	    printf("Cannot open fid file %s\n",flgs->phrotfile);
	    exit(1);
	  }
	}
      }
      i++;continue;
    }
    if(!strcmp(argv[i],"-scsl")){
      if((i+1!=argc)&&(argv[i+1][0]!='-')){ 
	flgs->scsl=2;
	strcpy(flgs->scslfile,argv[i+1]);
	i++;continue;
      }
      else {
	flgs->scsl=1;
	continue;
      }
    }
    if(!strcmp(argv[i],"-tabc")){
      if((i+1!=argc)&&(argv[i+1][0]!='-')){ 
	flgs->tabc=2;
	strcpy(flgs->petable,argv[i+1]);
	i++;continue;
      }
      else {
	flgs->tabc=1;
	continue;
      }
    }
    if(!strcmp(argv[i],"-slice")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-slice option requires an argument\n"); exit(1); }
      else {
	if(sscanf(argv[i+1],"%d-%d",&flgs->loslice,&flgs->hislice)!=2){
	  if(sscanf(argv[i+1],"%d",&flgs->loslice)!=1){
	    printf("Unable to interpret loslice-hislice argument of -slice option\n");
	    exit(1);
	  } else flgs->hislice=flgs->loslice;
	}
      }
      i++;continue;
    }

    if(!strcmp(argv[i],"-nav")){
      if((i+1==argc)||(argv[i+1][0]=='-')){ printf("-nav option requires an argument\n"); exit(
1); }
      else {
        if(!strcmp(argv[i+1],"n"))flgs->nav=NAVNONE;
        if(!strcmp(argv[i+1],"m"))flgs->nav=NAVMOD;
        if(!strcmp(argv[i+1],"p"))flgs->nav=NAVPHS;
        if(!strcmp(argv[i+1],"mp"))flgs->nav=NAVMOD+NAVPHS;
        if(!strcmp(argv[i+1],"pm"))flgs->nav=NAVMOD+NAVPHS;
	i++; continue;
      }
    }

    /* Depricated arguments */
    if(!strcmp(argv[i],"-series")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-ft")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-cp")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-phmap")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-intnav")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-shift")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-raw")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-rri")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-avwin")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-hu")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-orig")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-nodr")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-smod")){depricated(argv[i]);continue;}
    if(!strcmp(argv[i],"-hksref")){depricated(argv[i]);continue;}

    printf("flag %s not recognised\n",argv[i]);
    exit(1);
  }
} 

void help_message(void)
{
  printf("\nUsage: load_varian infile outfile [processing options]\n\n");
  printf(" -h            View help\n");
  printf(" -version      Print version\n");
  printf(" -ft2d         2Dfft\n");
  printf(" -ft3d         3Dfft\n");
  printf(" -bl           Baseline correct\n");
  printf(" -ms           Multislice reorder\n");
  printf(" -epi          EPI reorder\n");
  printf(" -ref fn       Ref scan name\n");
  printf(" -buo          Phase correction (Buonocore)\n");
  printf(" -buov         Phase correction (Buonocore vol-by-vol)\n");
  printf(" -con          Constrain phase correction\n");
  printf(" -mod          Output modulus\n");
  printf(" -phs          Output phase\n");
  printf(" -re           Output real\n");
  printf(" -im           Output imaginary\n");
  printf(" -short        Output signed integers (also -16)\n");
  printf(" -float        Output floating point\n");
  printf(" -rot          Rotate as necessary\n");
  printf(" -start n      Start from image\n");
  printf(" -num n        Images to process\n");
  printf(" -scsl [fn]    Scale slices (optional file name)\n");
  printf(" -pss          Reorder slices\n");
  printf(" -ftr          1Dfft rows only\n");
  printf(" -zero n       Zeropad to n square matrix\n");
  printf(" -zzero n      Zeropad to n slices\n");
  printf(" -scale n      Scale to n\n");
  printf(" -max n        Scale maximum to n\n");
  printf(" -vr n         Reorders volumes with step n\n");
  printf(" -sear n       Spin echo array reorder with step n\n");
  printf(" -lpss n       Loop pss reorder\n");
  printf(" -hks          Half k-space reconstruction\n");
  printf(" -hkspad       Half k-space zero pad only\n");
  printf(" -pocs         POCS half k-space reconstruction\n");
  printf(" -marg         Margosian half k-space reconstruction\n");
  printf(" -resl         Reslice to axial plane\n");
  printf(" -ph0 n        Zeroth order alternate echo phase correction\n");
  printf(" -ph1 n        First order alternate echo phase correction\n");
  printf(" -phrot n      Rotate the phase so every n volumes are all real\n");
  printf(" -ne n         Number of echoes\n");
  printf(" -kmb n        Kspace mask border (width n)\n");
  printf(" -imb n        Image mask border (width n)\n");
  printf(" -fmap         Field map reordering\n");
  printf(" -opss         Output pss to file\n");
  printf(" -revproc      Process volumes in reverse\n");
  printf(" -revload      Load volumes in reverse\n");
  printf(" -segx low num Segment image in x dimension\n");
  printf(" -segy low num Segment image in y dimension\n");
  printf(" -segz low num Segment image in z dimension\n");
  printf(" -slice n-n    Load slice or slice range\n");
  printf(" -fftw         Use FFTW functions (default)\n");
  printf(" -nrft         Use Numerical Recipes FT functions\n");
  printf(" -noro         Turns off detection for read-oversampling\n");
  printf(" -otep         Output optimum tep\n");
  printf(" -ctep         Correct tep error\n");
  printf(" -ogroa        Output groa\n");
  printf(" -ores n       Downsample by factor of n to low resolution matrix\n");
  printf(" -odrop        Output number of dropouts\n");
  printf(" -fermi        Apply fermi filter\n");
  printf(" -auto         Automatically determine process options\n");
  printf(" -oneseg       Process only the first segment of an EPI data set\n");
  printf(" -avvol        Average all volumes before processing\n");
  printf(" -lim n1 n2    Limits in self ref correction\n");
  printf(" -nav n/m/p/mp Navigator correction (none,modulus,phase)\n");
  printf(" -egr          Entropy ghost reduction\n");
  printf(" -egrv         Entropy ghost reduction (vol-by-vol)\n");
  printf(" -eiepi        Entropy IEPI correction\n");
  printf(" -resamp p l s Resample to p x l x s (can be number or p,l or s)\n");
  printf(" -tabc [file]  Perform PE table reordering (filename optional)\n");
  printf(" -avall        Output the mean of all volumes\n");
  printf(" -ssqall       Output the sum of the squares of all volumes\n");
  printf(" -info         Save the sequence parameters in File.info\n");
  printf(" -ppe          Phase encode offset (in mm)\n");
  printf(" -ppe2         2nd phase encode (slice) offset (in mm)\n");
  printf(" -gz           Output gzipped data\n");
  printf(" -multicoil    SSq average multiple coil data\n");
  printf(" -kpsf a b     PSF filter in the PE direction\n");
  printf(" \nLoadVarian version : %s\n",VERSION);

  exit(1);
}
int fid_exists(char *infile)
{
  char tmp[MAXLEN];
  sprintf(tmp,"%s/fid",infile);
  if(fopen(tmp,"rb")==NULL)return(0);
  else return(1);
}
void depricated(char *arg)
{
  printf("Flag %s has been removed from the latest version of load_varian.\n",arg);
  printf("Please use an older version or email stuart@fmrib.ox.ac.uk for details.\n");
  exit(1);
}
void auto_determine(char* filename,struct flags *flgs)
{
  int i,j,k,blank,ns,rcvrs;
  char tmp[MAXLEN],infile[MAXLEN];
  char **opts;
  int fid_exists(char *infile);
  void parse_arguments(int low,int argc,char** argv,struct flags *flgs);

  strcpy(infile,filename);
  if(!fid_exists(infile)){
    strcat(infile,".fid");
    if(!fid_exists(infile)){
      printf("Cannot open fid file %s\n",infile);
      return;
    }
  }

  flgs->out=OUT_MOD;
  flgs->sint=1;
  flgs->bl=1;
  flgs->rot=1;
  flgs->pss=1;
  flgs->ms=1;
  flgs->ft=2;

  if(read_procpar(infile,"nD ",tmp)){
    if(atoi(tmp)==3){
      flgs->ft=3;
      flgs->fermi=1;
    }
  }
  
  if(read_procpar(infile,"nv2 ",tmp)){
    if(atoi(tmp)>1) ns = atoi(tmp);
    else ns = read_procpar(infile,"pss ",tmp);
  }
  else ns = read_procpar(infile,"pss ",tmp);
  
  if(flgs->ft==3)read_procpar(infile,"t_thk ",tmp);
  else read_procpar(infile,"thk ",tmp);
  if((ns*atoi(tmp))>=100) flgs->resl=1;
  
  if(read_procpar_strings(infile,"load_varian_opts ",tmp)){
    opts=(char **)malloc(sizeof(char*));
    opts[0]=(char *)malloc(MAXLEN);
    if(tmp[0]!='\0'){
      i=0;j=0;k=0;
      blank=0;
      while(tmp[i]==' ')i++;
      while(1){
	if(tmp[i]=='\0')break;
	if(tmp[i]==' ')blank=1;
	else {
	  if(blank==1){
	    opts[j++][k]='\0';
	    opts=(char **)realloc(opts,(j+1)*sizeof(char*));
	    opts[j]=(char *)malloc(MAXLEN);
	    k=0;
	  }
	  opts[j][k]=tmp[i];
	  k++;
	  blank=0;
	}
	i++;
      }
      opts[j++][k]='\0';
      parse_arguments(0,j,opts,flgs);
      for(i=0;i<j;i++)free(opts[i]);
      free(opts);
      if(flgs->fm){
	flgs->out=OUT_CPLX;
	flgs->sint=0;
      }
    }
  }

  if(read_procpar(infile,"rcvrs ",tmp)){
    rcvrs=0;
    for(i=0;i<(int)strlen(tmp);i++){
      if(tmp[i]=='y')rcvrs++;
    }
    if(rcvrs>1)flgs->multicoil=1;
  }

  if(read_procpar(infile,"seq_type ",tmp)){
    if(!strcmp(tmp,"epi")){
      flgs->epi=1;
      flgs->ms=0;
      flgs->con=1;
      
      strcpy(infile,filename);
      if(strrchr(infile,'.')!=NULL)strrchr(infile,'.')[0]='\0';
      strcat(infile,"_ref.fid");
      if(!fid_exists(infile)){
	strcpy(infile,filename);
	if(strrchr(infile,'.')!=NULL)strrchr(infile,'.')[0]='\0';
	if(strrchr(infile,'.')!=NULL)strrchr(infile,'.')[0]='\0';
	strcat(infile,"_ref.fid");
	if(!fid_exists(infile)){
	  flgs->buo=1;
	}
	else {
	  /* 06sep05 SC: Make self-ref processing default */
	  /*flgs->ref=1;
	  strcpy(flgs->reffile,infile);*/
	  flgs->buo=1;
	}
      } else {
	/*flgs->ref=1;
	strcpy(flgs->reffile,infile);*/
	flgs->buo=1;
      }
    }
  }
  if (!flgs->epi) {
    if(read_procpar(infile,"seqcon ",tmp)){
      if(tmp[2]=='c'){
	if(tmp[3]=='s'){
	  if(!flash_converted(infile))flgs->ms=0;
	}
      }
    }
    if(read_procpar(infile,"petable ",tmp)){
      if(strlen(tmp)){
	if(!tab_converted(infile))flgs->tabc=1;
      }
    }
  }
}



void show_flags(struct flags flgs)
{
  printf("flgs.ft \t=\t%d\n",flgs.ft);
  printf("flgs.ms \t=\t%d\n",flgs.ms);
  printf("flgs.epi\t=\t%d\n",flgs.epi);
  printf("flgs.out\t=\t%d\n",flgs.out);
  printf("flgs.sint\t=\t%d\n",flgs.sint);
  printf("flgs.bl \t=\t%d\n",flgs.bl);
  printf("flgs.rot\t=\t%d\n",flgs.rot);
  printf("flgs.zero\t=\t%d\n",flgs.zero);
  printf("flgs.zzero\t=\t%d\n",flgs.zzero);
  printf("flgs.ref\t=\t%d\n",flgs.ref);
  printf("flgs.buo\t=\t%d\n",flgs.buo);
  printf("flgs.buov\t=\t%d\n",flgs.buov);
  printf("flgs.con\t=\t%d\n",flgs.con);
  printf("flgs.fftw\t=\t%d\n",flgs.fftw);
  printf("flgs.rover\t=\t%d\n",flgs.rover);
  printf("flgs.sover\t=\t%d\n",flgs.sover);
  printf("flgs.seg\t=\t%d\n",flgs.seg);
  printf("flgs.start\t=\t%d\n",flgs.start);
  printf("flgs.num\t=\t%d\n",flgs.num);
  printf("flgs.me \t=\t%d\n",flgs.me);
  printf("flgs.revproc\t=\t%d\n",flgs.revproc);
  printf("flgs.revload\t=\t%d\n",flgs.revload);
  printf("flgs.scsl\t=\t%d\n",flgs.scsl);
  printf("flgs.resl\t=\t%d\n",flgs.resl);
  printf("flgs.pss\t=\t%d\n",flgs.pss);
  printf("flgs.vr \t=\t%d\n",flgs.vr);
  printf("flgs.hks\t=\t%d\n",flgs.hks);
  printf("flgs.mph\t=\t%d\n",flgs.mph);
  printf("flgs.fm \t=\t%d\n",flgs.fm);
  printf("flgs.sear\t=\t%d\n",flgs.sear);
  printf("flgs.lpss\t=\t%d\n",flgs.lpss);
  printf("flgs.phrot\t=\t%d\n",flgs.phrot);
  printf("flgs.opss\t=\t%d\n",flgs.opss);
  printf("flgs.kmb\t=\t%d\n",flgs.kmb);
  printf("flgs.imb\t=\t%d\n",flgs.imb);
  printf("flgs.nosave\t=\t%d\n",flgs.nosave);
  printf("flgs.lowx\t=\t%d\n",flgs.lowx);
  printf("flgs.lowy\t=\t%d\n",flgs.lowy);
  printf("flgs.lowz\t=\t%d\n",flgs.lowz);
  printf("flgs.numx\t=\t%d\n",flgs.numx);
  printf("flgs.numy\t=\t%d\n",flgs.numy);
  printf("flgs.numz\t=\t%d\n",flgs.numz);
  printf("flgs.overx\t=\t%d\n",flgs.overx);
  printf("flgs.overy\t=\t%d\n",flgs.overy);
  printf("flgs.overz\t=\t%d\n",flgs.overz);
  printf("flgs.overv\t=\t%d\n",flgs.overv);
  printf("flgs.loslice\t=\t%d\n",flgs.loslice);
  printf("flgs.hislice\t=\t%d\n",flgs.hislice);
  printf("flgs.slices\t=\t%d\n",flgs.slices);
  printf("flgs.otep\t=\t%d\n",flgs.otep);
  printf("flgs.ogroa\t=\t%d\n",flgs.ogroa);
  printf("flgs.ogrof\t=\t%d\n",flgs.ogrof);
  printf("flgs.ctep\t=\t%d\n",flgs.ctep);
  printf("flgs.ds \t=\t%d\n",flgs.ds);
  printf("flgs.odrop\t=\t%d\n",flgs.odrop);
  printf("flgs.phrotf\t=\t%d\n",flgs.phrotf);
  printf("flgs.fermi\t=\t%d\n",flgs.fermi);
  printf("flgs.egr\t=\t%d\n",flgs.egr);
  printf("flgs.egrv\t=\t%d\n",flgs.egrv);
  printf("flgs.ecf\t=\t%d\n",flgs.ecf);
  printf("flgs.eiepi\t=\t%d\n",flgs.eiepi);
  printf("flgs.oneseg\t=\t%d\n",flgs.oneseg);
  printf("flgs.avvol\t=\t%d\n",flgs.avvol);
  printf("flgs.nav\t=\t%d\n",flgs.nav);
  printf("flgs.epik\t=\t%d\n",flgs.epik);
  printf("flgs.avall\t=\t%d\n",flgs.avall);
  printf("flgs.rcvrs\t=\t%d\n",flgs.rcvrs);
  printf("flgs.tabc\t=\t%d\n",flgs.tabc);
  printf("flgs.resamp[0]\t=\t%d\n",flgs.resamp[0]);
  printf("flgs.resamp[1]\t=\t%d\n",flgs.resamp[1]);
  printf("flgs.resamp[2]\t=\t%d\n",flgs.resamp[2]);
  printf("flgs.info\t=\t%d\n",flgs.info);
  printf("flgs.gz\t=\t%d\n",flgs.gz);
  printf("flgs.multicoil\t=\t%d\n",flgs.multicoil);
  printf("flgs.kpsf\t=\t%d\n",flgs.kpsf);

  printf("flgs.scale\t=\t%f\n",flgs.scale);
  printf("flgs.min\t=\t%f\n",flgs.min);
  printf("flgs.max\t=\t%f\n",flgs.max);
  printf("flgs.ph0\t=\t%f\n",flgs.ph0);
  printf("flgs.ph1\t=\t%f\n",flgs.ph1);
  printf("flgs.phoff\t=\t%f\n",flgs.phoff);
  printf("flgs.ppe\t=\t%f\n",flgs.ppe);
  printf("flgs.ppe2\t=\t%f\n",flgs.ppe2);
  printf("flgs.kpsf_a\t=\t%f\n",flgs.kpsf_a);
  printf("flgs.kpsf_b\t=\t%f\n",flgs.kpsf_b);

  printf("flgs.reffile\t=\t%s\n",flgs.reffile);
  printf("flgs.phrotfile\t=\t%s\n",flgs.phrotfile);
  printf("flgs.scslfile\t=\t%s\n",flgs.scslfile);
  printf("flgs.petable\t=\t%s\n",flgs.petable);
  printf("flgs.seqcon\t=\t%s\n",flgs.seqcon);
}
