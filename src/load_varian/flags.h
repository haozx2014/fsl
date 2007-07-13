#if !defined __FLAGS_H
#define __FLAGS_H

#include "maxlen.h"

struct flags {
  char scslfile[MAXLEN],reffile[MAXLEN],phrotfile[MAXLEN];
  char petable[MAXLEN],seqcon[6];
  int ft,ms,epi,out,sint,bl,rot,zero,zzero;
  int ref,buo,buov,con,fftw,rover,seg,sover;
  int start,num,me,revproc,revload;
  int scsl,resl,pss,vr,hks,mph,fm,sear;
  int lpss,phrot,opss,kmb,imb,nosave;
  int lowx,lowy,lowz,numx,numy,numz;
  int overx,overy,overz,overv;
  int loslice,hislice,slices,otep,ogroa,ogrof,ctep;
  int ds,odrop,phrotf,fermi,egr,egrv,ecf,eiepi;
  int oneseg,avvol,nav,epik,kpsf;
  int avall,rcvrs,tabc,info,gz,multicoil,dim;
  int etl,cusfilt,fdf,fract_ky;
  float scale,min,max,ph0,ph1,phoff;
  float ppe,ppe2;
  float kpsf_a,kpsf_b;
  int resamp[3];

  /* Flags not yet in verbose */
};

#define OUT_CPLX 0
#define OUT_MOD 1
#define OUT_PHS 2
#define OUT_REAL 3
#define OUT_IMAG 4


#endif
