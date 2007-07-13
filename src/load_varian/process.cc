/*
  process.c: Main function for turning Varian fids to images

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "process.h"

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>

#include "procfunc.h"
#include "phasefunc.h"
#include "read_varian.h"
#include "image_ft.h"
#include "pssreorder.h"
#include "scaleslice.h"
#include "segment.h"
#include "maxlen.h"
#include "hfunc.h"
#include "hermit.h"
#include "pocs.h"
#include "margosian.h"
#include "reslice.h"
#include "kspace_mask.h"
#include "image_fftw.h"
#include "tep.h"
#include "dropouts.h"
#include "fermi.h"
#include "entghost.h"
#include "navigator.h"
#include "xfm.h"
#include "resamp.h"
#include "concat_info.h"
#include "zlib.h"
#include "fdf.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

#define SOVER_FRAC 0.8

int process(char *infile,char *outfile,struct flags *proc)
{
  char orient;
  char tmp[MAXLEN];
  short *sdata;
  int ppl,lpi,ipv,vols,type,index,swap,dropi,dropn,dropc,dropl;
  int tppl,tlpi,tipv,tvols,zppl,zlpi,zipv,hlpi,cplx,i,j,fix,refs;
  int oindex,maxset,partial_lpi;
  int *cen,*petable;
  unsigned long datasize,volsize,l;
  unsigned long chunksize,chunks,offset,step,refindex;
  float vx,vy,vz,max,min,ftmp,tep,groa,grof;
  float resampx,resampy,resampz,ppe_phase;
  float *data,*average,*ref,*phrot,*navigator,*ctep,*drop,*adata;
  float *scsl_array,*pss_array,*dropouts;
  float *egr;
  FILE *ofp;
  gzFile ofpz;
  struct fdf fdfheader;
  struct dsr header;
  struct xfm transform;
  int twonup(int);

  zppl = zlpi = zipv = 0;
  tppl = tlpi = tipv = tvols = 0;
  phrot = NULL;
  navigator = scsl_array = pss_array = NULL;
  petable = NULL;
  dropi=dropn=dropc=dropl=0;
  hlpi=0;
  refindex=i=0;
  refs=1;
  tep=0.0;
  resampx=resampy=resampz=0;
  average=ctep=drop=adata=dropouts=egr=NULL;
  ofp=NULL;
  ofpz=NULL;


  /* ==================================================================== */
  /*        R E A D   P R O C P A R   F O R   D A T A   S I Z E           */
  /* ==================================================================== */
  
  /* Do the info file read first */
  if(proc->info){
    concat_info_files(infile,outfile);
  }

  if(read_varian_header(infile,&header,&transform)){
    printf("Read varian header error.\n");
  }
  
  /* Calculate input and output volume sizes */
  avw_get_orient(&header,&orient);
  avw_get_dim(&header,&ppl,&lpi,&ipv,&vols);
 
  if(proc->fdf){
    avw_get_vox(&header,&vx,&vy,&vz);
    fdf_set_vox(&fdfheader,vx,vy,vz);
  }

  if(read_procpar(infile,"nD ",tmp)){
    proc->dim=atoi(tmp);
  }
  if(read_procpar(infile,"fract_ky ",tmp)){
    proc->fract_ky=atoi(tmp);
  }

  if(proc->rcvrs>1){
    if(proc->num!=0)proc->num*=proc->rcvrs;
    if(proc->start!=0)proc->start*=proc->rcvrs;
  }
  if(proc->num==0)proc->num=vols;
  if(proc->start>=vols)proc->start=0;
  if((proc->start+proc->num)>vols)proc->num=vols-proc->start;

  /* Slice load */
  if(proc->loslice>=0){
    if(proc->hislice>=ipv){
      printf("Specified slice is out of range\n");
      return(1);
    }
    proc->slices=ipv;
    ipv=proc->hislice-proc->loslice+1;
  }

  datasize = ppl*lpi*ipv*2;

  /* ==================================================================== */
  /*             S E T   C H U N K S I Z E   A N D   S T E P              */
  /* ==================================================================== */

  chunksize = datasize; chunks = 1; offset = 0; step = 1;
  
  if(proc->me>1){
    if(proc->ms)
      {chunksize = ppl*2; chunks = lpi*ipv; step = proc->me;}
    else 
      {chunksize = ppl*lpi*2; chunks = ipv; step = proc->me;}
  }

  if(proc->fm)
    {chunksize = ppl*ipv*2; chunks = lpi;step = 2;}

  if(proc->sear)
    {chunksize = ppl*ipv*2; chunks = lpi; step = proc->sear;}
  
  if(proc->rcvrs>1){
    if(proc->ms) {
      if(proc->seqcon[1]=='c'){
	if(proc->seqcon[2]=='c'){
	  chunksize = ppl*lpi*ipv*2; chunks = 1; step = proc->rcvrs;
	} else {
	  chunksize = ppl*ipv*2; chunks = lpi; step = proc->rcvrs;
	}
      }
      else {
	if(proc->seqcon[3]=='c'){
	  chunksize = ppl*ipv*2; chunks = lpi; step = proc->rcvrs;
	} else {
	  chunksize = ppl*2; chunks = lpi*ipv; step = proc->rcvrs;
	}
      }
    }
    else {
      if(!proc->epi)
	{chunksize = ppl*lpi*2; chunks = ipv; step = proc->rcvrs;}
    }
  }
  if(proc->loslice>=0){
    if(proc->ms){
      chunksize = ppl*2; chunks = ipv*lpi; step = proc->slices;
      offset = proc->loslice;
    } else {
      chunksize = ppl*lpi*2; chunks = ipv; step = 1;
      offset = proc->loslice;
    }
  }

  if(proc->oneseg){
    if(proc->epik){
      datasize = ppl*lpi*ipv*2;
      chunksize = datasize;
      proc->epik=0;
    } else {
      chunksize /= proc->seg;
      lpi /= proc->seg;
    }
    proc->oneseg = proc->seg;
    proc->seg = 1;
  }

  if(proc->tabc){
    if(!strlen(proc->petable)){
      read_procpar_strings(infile,"petable",tmp);
      /* Look for petable within fid file */
      sprintf(proc->petable,"%s/%s",infile,tmp);
    }
    if(read_petable(proc->petable,&petable)!=lpi){
      /* Look for petable within fsldir */
      sprintf(proc->petable,"%s/etc/petables/%s",getenv("FSLDIR"),tmp);
      if(read_petable(proc->petable,&petable)!=lpi){
	printf("PEtable is unreadable or of wrong size/n");
	return(1);
      }
    }
  }

  /* ==================================================================== */
  /*              Z E R O   F I L L   C A L C U L A T I O N S             */
  /* ==================================================================== */
  
  if(proc->zero){
    zppl=twonup(proc->zero);
    zlpi=twonup(proc->zero);
    if(proc->rover)zppl*=2;
  }
  else {
    if(proc->nav)tlpi=lpi-proc->seg;
    else tlpi=lpi;
    if(proc->hks){
      zppl=twonup(ppl);
      zlpi=twonup(tlpi);
      if(read_procpar(infile,"hks_nv",tmp)){
	hlpi=atoi(tmp);
      } else {
	hlpi=zlpi;
      }
      if(proc->fract_ky){
	hlpi=(lpi-proc->fract_ky)*2;
      }
    }
    else {zppl=ppl;zlpi=tlpi;}
    if(proc->fract_ky){
      zlpi=(lpi-proc->fract_ky)*2;
      hlpi=(lpi-proc->fract_ky)*2;
      if(!proc->hks)proc->hks=1;
    }
  }
  if(proc->zzero){
    zipv=twonup(proc->zzero);
  }
  else zipv=ipv;
  
  if((proc->hks)&&(zppl!=zlpi))zlpi=zppl;

  if(zlpi<lpi) volsize=zppl*lpi*zipv*2;
  else volsize=zppl*zlpi*zipv*2;
  
  /* Voxel dimensions */
  if(volsize>datasize){
    avw_get_vox(&header,&vx,&vy,&vz);
    if(zipv>ipv)proc->zero=2;
    else proc->zero=1;
    vx*=((float)ppl/(float)zppl);
    vy*=((float)lpi/(float)zlpi);
    vz*=((float)ipv/(float)zipv);
    avw_set_vox(&header,vx,vy,vz);
    if(proc->fdf)fdf_set_vox(&fdfheader,vx,vy,vz);
  }
  else proc->zero=0;

  /* Data type */
  avw_get_dt(&header,&type);
  if(proc->out){
    if(proc->sint) type=DT_SIGNED_SHORT;
    else type=DT_FLOAT;
    cplx=1;
  }
  else {
    type=DT_COMPLEX;
    cplx=2;
  }
  avw_set_dt(&header,type);

  /* The scale slice option */
  if(proc->scsl){
    if(proc->scsl>1){
      if(scale_slice_array_file(infile,&scsl_array,proc->scslfile))
	proc->scsl=0;
    } else {
      if(scale_slice_array(infile,&scsl_array))
	proc->scsl=0;
    }
  }
  else scsl_array=NULL;

  /* Find pss to reorder later */
  if(proc->pss){
    proc->pss=pss_to_array(infile,&pss_array);
    if(proc->loslice>0)
      for(i=0;i<ipv;i++)pss_array[i]=pss_array[i+proc->loslice];
  }
  else {
    pss_array=NULL;
  }

  if(proc->opss)output_pss(infile,outfile);
  
  /* ==================================================================== */
  /*     A L L O C A T E   M E M O R Y   A N D   O P E N   F I L E S      */
  /* ==================================================================== */
  
  if((data=(float *)malloc(volsize*sizeof(float)))==NULL){
    printf("Malloc failed (data)");
    return(1);
  }
  if((sdata=(short *)malloc(volsize*sizeof(short)))==NULL){
    printf("Malloc failed (sdata)");
    return(1);
  }
  if((proc->avall)||(proc->multicoil)){
    if((average=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed (data)");
      return(1);
    }
    for(l=0;l<volsize;l++)average[l]=0;
  }
  if((proc->phrot)||(proc->phrotf)){
    if((phrot=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed (phrot)");
      return(1);
    }
  }
  if(!proc->nosave){
    if(proc->gz){
      sprintf(tmp,"%s.img.gz",outfile);
      if((ofpz=gzopen(tmp,"wb"))==NULL){
	printf("Failed to open file %s\n",tmp);
	return(1);
      }  
    } else {
      sprintf(tmp,"%s.img",outfile);
      if((ofp=fopen(tmp,"wb"))==NULL){
	printf("Failed to open file %s\n",tmp);
	return(1);
      }
    }
  }
  if(proc->odrop){
    if((drop=(float *)malloc(volsize*11*sizeof(float)))==NULL){
      printf("Malloc failed (drop)");
      return(1);
    }
    if((dropouts=(float *)malloc(ipv*proc->num*sizeof(float)))==NULL){
      printf("Malloc failed (dropouts)");
      return(1);
    }
    dropi=dropn=dropc=dropl=fix=0;
  }
  if(proc->avvol){
    if((adata=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed (adata)");
      return(1);
    }
    for(l=0;l<volsize;l++)adata[l]=0;
  }
  if(proc->egr){
    if((egr=(float *)malloc(ipv*sizeof(float)))==NULL){
      printf("Malloc failed (egr)");
      return(1);
    }
  }
  if(proc->nav){
    if((navigator=(float *)malloc(ppl*2*proc->seg*ipv*sizeof(float)))==NULL){
      printf("Malloc failed (navigator)");
      return(1);
    }
  }


  /* ==================================================================== */
  /*                      T E P   C O R R E C T I O N                     */
  /* ==================================================================== */
  
  if(proc->ctep){

    if((ctep=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed (ctep)");
      return(1);
    }

    if(proc->revproc){
      index = proc->num+proc->start-1;
      if(proc->loslice>=0)
	offset = (index*proc->slices)+proc->loslice;
      else
	offset = index;
    }

    /* Read in first volume */
    if(read_varian_chunks(infile,chunksize,chunks,offset,step,ctep)){
      printf("Failed to read main data file\n");
      return(1);
    }

    tppl=ppl; tlpi=lpi; tipv=ipv;

    if(proc->ds){
      downsample(ctep,tppl,tlpi,tipv,proc->ds);
      tppl/=proc->ds; tlpi/=proc->ds;
    }

    epireorder(ctep,tppl,tlpi,tipv,proc->seg,1);
    lv_modulus(ctep,tppl*tlpi*tipv*2);
    tep = optimum_tep(ctep,tppl,tlpi,tipv);

    free(ctep);
  }

  /* ==================================================================== */
  /*          R E F   S C A N   P H A S E   C O R R E C T I O N           */
  /* ==================================================================== */
  
  if((proc->ref)&&(!proc->buo)){
    
    refs=volumes_in_varian_header(proc->reffile);
    if(!refs)refs=1;
    if(refs>proc->rcvrs)refs=proc->rcvrs;
    
    if((ref=(float *)malloc(volsize*refs*sizeof(float)))==NULL){
      printf("Malloc failed (ref)");
      return(1);
    }
    
    for(j=0;j<refs;j++){
      
      refindex = volsize*j;
      offset = j;
      
      if(proc->loslice>=0)offset=(i*proc->slices)+proc->loslice;
      
      if(read_varian_chunks(proc->reffile,chunksize/proc->seg,chunks,offset,step,
			    &ref[refindex])){
	printf("Failed to read ref scan\n");
	return(1);
      }
      
      /* Read in image volume to find centre kspace */
      if(read_varian_chunks(infile,chunksize/proc->seg,chunks,offset,step,data)){
	printf("Failed to read main data file\n");
	return(1);
      }
      
      tppl=ppl; tlpi=lpi; tipv=ipv;
      
      if(proc->ds){
	downsample(data,tppl,tlpi/proc->seg,tipv,proc->ds);
	downsample(&ref[refindex],tppl,tlpi/proc->seg,tipv,proc->ds);
	tppl/=proc->ds; tlpi/=proc->ds;
      }
      
      
      epireorder(data,tppl,tlpi/proc->seg,tipv,1,1);
      epireorder(&ref[refindex],tppl,tlpi/proc->seg,tipv,1,1);
      if(proc->ctep)tep_interpolate_complex(data,tppl,tlpi*tipv,tep);
      if(proc->ctep)tep_interpolate_complex(&ref[refindex],tppl,tlpi*tipv,tep);
      
      if(proc->pss)pss_reorder(data,tppl*(tlpi/proc->seg)*2,tipv,pss_array);
      if(proc->pss)pss_reorder(&ref[refindex],tppl*(tlpi/proc->seg)*2,tipv,pss_array);
      
      if(proc->zero){
	zerofill_3d(data,tppl,tlpi/proc->seg,tipv,
		    zppl,zlpi/proc->seg,zipv,proc->hks);
	zerofill_3d(&ref[refindex],tppl,tlpi/proc->seg,tipv,
		    zppl,zlpi/proc->seg,zipv,proc->hks);
	tppl=zppl; tlpi=zlpi; tipv=zipv;
      }
      
      if((cen=(int *)malloc(tipv*sizeof(int)))==NULL){
	printf("Malloc failed (cen)\n");
	return(1);
      }      
      /* Dosn't seem to work as well so removed for the time being */
      /*find_centre_kspace_array(data,tppl,tlpi,tipv,cen);*/
      cen[0]=find_centre_kspace(data,tppl,tlpi/proc->seg,tipv);
      for(i=1;i<tipv;i++)cen[i]=cen[0];
      
      if(proc->fftw)
	image_oned_fftw_r(&ref[refindex],tppl,tlpi/proc->seg,tipv);
      else 
	image_oned_fft_r(&ref[refindex],tppl,tlpi/proc->seg,tipv);
      
      if(proc->con)calc_phase_cor(&ref[refindex],tppl,tlpi,tipv,proc->seg,proc->con,cen);
      
      free(cen);
    }
  }
  else {
    ref=NULL;
  }

  /* ==================================================================== */
  /*         B U O N O C O R E   P H A S E   C O R R E C T I O N          */
  /* ==================================================================== */
  
  if(proc->buo){
    if((ref=(float *)malloc(volsize*proc->rcvrs*sizeof(float)))==NULL){
      printf("Malloc failed (ref)");
      return(1);
    }
  
    refs=proc->rcvrs;

    for(j=0;j<refs;j++){
      
      refindex = volsize*j;
      if(proc->ref){
	offset = j;
	if(proc->loslice>=0)offset=(i*proc->slices)+proc->loslice;
	if(read_varian_chunks(proc->reffile,chunksize,chunks,offset,step,data)){
	  printf("Failed to read ref scan\n");
	  return(1);
	}
      } else {
	if(proc->revproc){
 	  index = proc->num+proc->start-1-j;
	  if(proc->loslice>=0)
	    offset = (index*proc->slices)+proc->loslice;
	  else
	    offset = index;
	} else {
	  offset=j;
	}
	/* Read in volume */
	if(read_varian_chunks(infile,chunksize,chunks,offset,step,data)){
	  printf("Failed to read main data file\n");
	  return(1);
	}
      }
      
      tppl=ppl; tlpi=lpi; tipv=ipv;
      
      if(proc->ds){
	downsample(data,tppl,tlpi,tipv,proc->ds);
	tppl/=proc->ds; tlpi/=proc->ds;
      }
      
      epireorder(data,tppl,tlpi,tipv,proc->seg,1);
      if(proc->ctep)tep_interpolate_complex(data,tppl,tlpi*tipv,tep);
      if(proc->bl)median_baseline(data,tppl,tlpi,tipv,1);
      if(proc->pss)pss_reorder(data,tppl*tlpi*2,tipv,pss_array);
      if(proc->zero){
	zerofill_3d(data,tppl,tlpi,tipv,zppl,zlpi,zipv,proc->hks);
	tppl=zppl; tlpi=zlpi;tipv=zipv;
      }
      
      if(proc->fftw)
	image_oned_fftw_r(data,tppl,tlpi,tipv);
      else
	image_oned_fft_r(data,tppl,tlpi,tipv);
      
      buono_calib(data,&ref[refindex],tppl,tlpi,tipv,proc->con);
    }
  }

  /* ==================================================================== */
  /*             M A N U A L   P H A S E   C O R R E C T I O N            */
  /* ==================================================================== */
  if(proc->mph){
    if((ref=(float *)malloc(volsize*sizeof(float)))==NULL){
      printf("Malloc failed (ref)");
      return(1);
    }
    if(proc->zero){tppl=zppl; tlpi=zlpi; tipv=ipv;}
    else {tppl=ppl; tlpi=lpi; tipv=ipv;}
    create_ref_scan(ref,tppl,tlpi,tipv,proc->seg,proc->ph0,proc->ph1);
    proc->ref=1;
    refs=1;
  }

  /* ==================================================================== */
  /*                  P H A S E   R O T A T I O N   F I L E               */
  /* ==================================================================== */
  if(proc->phrotf){
    tppl=ppl; tlpi=lpi; tipv=ipv;
    
    if(read_varian_chunks(proc->phrotfile,chunksize,chunks,0,step,phrot)){
      printf("Failed to read phase rotate ref file\n");
      return(1);
    }
    
    if(proc->ms)msreorder(phrot,tppl,tlpi,tipv,1);
    if(proc->epi)epireorder(phrot,tppl,tlpi,tipv,proc->seg,1);
    if(proc->pss)pss_reorder(phrot,tppl*tlpi*2,tipv,pss_array);
    if(proc->ft==3){
      reverse_slices(phrot,tppl,tlpi,tipv);
      if(proc->zero){
	zerofill_3d(phrot,tppl,tlpi,tipv,zppl,zlpi,zipv,proc->hks);
	tppl=zppl; tlpi=zlpi; tipv=zipv;
      }
      image_oned_fftw_sp(phrot,tppl,tlpi,tipv,proc->phoff);
    } else {
      if(proc->zero){
	if(proc->zero>1)image_oned_ifftw_s(phrot,tppl,tlpi,tipv);
	zerofill_3d(phrot,tppl,tlpi,tipv,zppl,zlpi,zipv,proc->hks);
	tppl=zppl; tlpi=zlpi; tipv=zipv;
	if(proc->zero>1) image_oned_fftw_s(phrot,tppl,tlpi,tipv);
      }
    }
    image_oned_fftw_r(phrot,tppl,tlpi,tipv);
    
    if((proc->ref)&&(!proc->buo))
      epi_phasecor(phrot,ref,tppl,tlpi,tipv);
    if(proc->buo){buono_apply(phrot,ref,tppl,tlpi,tipv);}
    if(proc->buov)buono_cor(phrot,tppl,tlpi,tipv,proc->con);
    if(proc->hks==2)hks_herm(phrot,tppl,tlpi,tipv);
    if(proc->hks==3)hks_pocs(phrot,tppl,tlpi,tipv);
    if(proc->hks==4)hks_margosian(phrot,tppl,tlpi,tipv);
    if((proc->hks)&&(proc->rover)){
      fsegment(phrot,0,tppl,tppl,(tlpi-hlpi)/2,hlpi,tlpi,0,tipv,tipv,2);
      tlpi=hlpi;
    }
    if(!proc->ft)image_oned_ifftw_r(phrot,tppl,tlpi,tipv);
    if(proc->ft>=2)image_oned_fftw_c(phrot,tppl,tlpi,tipv);
    if((proc->ft)&&(proc->rover)){
      fsegment(phrot,tppl/4,tppl/2,tppl,0,tlpi,tlpi,0,tipv,tipv,2);
      tppl/=2;
    }
    if(proc->rot){
      if(transform.rotate){
	rotate(phrot,tppl,tlpi,tipv,proc->rot);
	swap=tlpi;tlpi=tppl;tppl=swap;
      }
      if(transform.reflect==1){
	reflect_lines(phrot,tppl,tlpi,tipv);
	reverse_slices(phrot,tppl,tlpi,tipv);
      }
      if(transform.reflect==2){
	reflect_points(phrot,tppl,tlpi,tipv);
	reverse_slices(phrot,tppl,tlpi,tipv);
      }
    }
    if(proc->resl){
      switch(orient){
      case 1:
	cor2ax(phrot,tppl,tlpi,tipv);
	swap=tipv;tipv=tlpi;tlpi=swap;
	break;
      case 2:
	sag2ax(phrot,tppl,tlpi,tipv);
	swap=tppl;tppl=tipv;tipv=swap;
	swap=tipv;tipv=tlpi;tlpi=swap;
	break;
      }
    }
  }

  /* ==================================================================== */
  /* ==================================================================== */
  /*     S T A R T   O F    M A I N   P R O C E S S I N G   L O O P       */
  /* ==================================================================== */
  /* ==================================================================== */

  /* Loop over volumes */
  maxset=0;
  for(i=proc->start;i<proc->start+proc->num;i++){
    tppl = ppl; tlpi = lpi; tipv = ipv;
    index=i;
    if((proc->revproc)||(proc->revload))index=proc->num+proc->start-1-i;
    if(proc->vr)index = (i%(vols/proc->vr))*proc->vr + (i/(vols/proc->vr));
    refindex=volsize*((i-proc->start)%refs);
    
    /* Calculate offset */
    offset = index;
    if(proc->me>1){
      if(proc->ms) offset = (index/proc->me)*lpi*ipv*proc->me+index%proc->me;
      else offset = (index/proc->me)*ipv*proc->me+index%proc->me;
    }
    if(proc->rcvrs>1){
      if(proc->ms){
	if(proc->seqcon[1]=='c') {
	  offset = (index/proc->rcvrs)*lpi*proc->rcvrs+index%proc->rcvrs;
	} else {
	  if(proc->seqcon[3]=='c') {
	    offset = (index/proc->rcvrs)*proc->rcvrs+index%proc->rcvrs;
	  } else {
	    offset = (index/proc->rcvrs)*lpi*ipv*proc->rcvrs+index%proc->rcvrs;
	  }
	}
      } else {
	if(!proc->epi)
	  offset = (index/proc->rcvrs)*ipv*proc->rcvrs+index%proc->rcvrs;
      }
    }
    if(proc->fm) offset = (index/2)*lpi*2+index%2;
    if(proc->sear) offset = (index/proc->sear)*lpi*proc->sear+index%proc->sear;
    if(proc->lpss) {
      chunksize = ppl*lpi*2; chunks = ipv; step = ipv/proc->lpss;
      offset = (index/(ipv/proc->lpss))*ipv + index%(ipv/proc->lpss);
    }
    if(proc->loslice>=0){
      if(proc->ms) offset = (index*proc->slices)+proc->loslice;
      else offset = (index*proc->slices)+proc->loslice;
    }
    if(proc->oneseg){
      offset*=proc->oneseg;
    }
    
    /* Read in data */
    if((proc->lpss)&&(proc->rcvrs>1)){
      for(j=0;j<ipv/proc->lpss;j++){
	chunks=proc->lpss;
	offset = ((index/proc->rcvrs)/(ipv/proc->lpss))*ipv 
	  + ((index/proc->rcvrs)%(ipv/proc->lpss)) + (index%proc->rcvrs)*ipv + j*ipv*proc->rcvrs;
	if(read_varian_chunks(infile,chunksize,chunks,offset,step,&data[j*chunksize*chunks])){
	  printf("Failed to read main data file\n");
	  return(1);
	}
      }
    } else {
      if(read_varian_chunks(infile,chunksize,chunks,offset,step,data)){
	printf("Failed to read main data file\n");
	return(1);
      }
    }

    if((proc->revload)||(proc->vr))index=i;

    if(proc->avvol){
      for(l=0;l<volsize;l++)adata[l]+=(data[l]/proc->num);
      if(i<(proc->start+proc->num)-1)continue;
      else{
	for(l=0;l<volsize;l++)data[l]=adata[l];
	index = proc->start;
      }
    }
    oindex=index-proc->start;
    
    if(proc->nav){
      extract_navigators(data,navigator,tppl,tlpi,tipv,proc->seg);
      tlpi-=proc->seg;
    }
    if(proc->ds){
      downsample(data,tppl,tlpi,tipv,proc->ds);
      tppl/=proc->ds; tlpi/=proc->ds;
    }
    
  /* ==================================================================== */
  /*             M A I N   P R O C E S S I N G   S T E P S                */
  /* ==================================================================== */
    
    if(proc->lpss)looppss_reorder(data,tppl,tlpi,tipv,proc->lpss,(index/proc->rcvrs)%(ipv/proc->lpss));
    if(proc->ms)msreorder(data,tppl*proc->etl,tlpi/proc->etl,tipv,1);
    if(proc->tabc)tabc(data,tppl,tlpi,tipv,petable);
    if(proc->epi){
      if(proc->epik){
	epikreorder(data,tppl,tlpi,tipv,proc->seg,proc->epik,1);
	tlpi=tlpi*proc->epik/(proc->epik+(proc->seg-1));
      } else epireorder(data,tppl,tlpi,tipv,proc->seg,1);
    }
    if(proc->ctep)tep_interpolate_complex(data,tppl,tlpi*tipv,tep);
    if(proc->kpsf)kpsf_filter(data,tppl,tlpi,tipv,proc->kpsf_a,proc->kpsf_b);
    /*if(proc->cusfilt)custom_filter(data,tppl,tlpi,tipv,infile);*/

    if(proc->dim==3){
      if(fabs(proc->ppe2)>0.0){
	avw_get_vox(&header,&vx,&vy,&vz);
	ppe_phase=(2*M_PI)*((proc->ppe2)/vz);
	apply_ppe2(data,tppl,tlpi,tipv,ppe_phase);
      }
      if(proc->bl)median_baseline_3d(data,tppl,tlpi,tipv,1,proc->hks);
      if(proc->zero){
	zerofill_3d(data,tppl,tlpi,tipv,zppl,zlpi,zipv,proc->hks);
	tppl=zppl; tlpi=zlpi; tipv=zipv;
      }
      if(proc->hks==2){
	image_oned_fftw_s(data,tppl,tlpi,tipv);
	hks_centre_kspace(data,tppl,tlpi,tipv);
	image_oned_fftw_r(data,tppl,tlpi,tipv);
	hks_phase_images(data,tppl,tlpi,&partial_lpi,tipv);
	image_oned_ifftw_r(data,tppl,tlpi,tipv);
	hks_centre_kspace(data,tppl,tlpi,tipv);
	hks_hermitian_conjugation(data,tppl,tlpi,partial_lpi,tipv);
	image_oned_ifftw_s(data,tppl,tlpi,tipv);
      }
      if(proc->fermi)fermi_filter_3d(data,tppl,tlpi,tipv);
      if(proc->kmb)kspace_mask_border(data,tppl,tlpi,tipv,proc->kmb);
      if(proc->ft){
	if(proc->fftw)
	  image_oned_fftw_sp(data,tppl,tlpi,tipv,proc->phoff);
	else
	  image_oned_fft_sp(data,tppl,tlpi,tipv,proc->phoff);
      }
      reverse_slices(data,tppl,tlpi,tipv);
    } else {
      if(proc->bl)median_baseline(data,tppl,tlpi,tipv,1);
      if(proc->fermi)fermi_filter_2d(data,tppl,tlpi,tipv);
      if(proc->kmb)kspace_mask_border(data,tppl,tlpi,tipv,proc->kmb);
      if(proc->zero){
	if(proc->zero>1){
	  if(proc->fftw)
	    image_oned_ifftw_s(data,tppl,tlpi,tipv);
	  else
	    image_oned_ifft_s(data,tppl,tlpi,tipv);
	}
	zerofill_3d(data,tppl,tlpi,tipv,zppl,zlpi,zipv,proc->hks);
	tppl=zppl; tlpi=zlpi; tipv=zipv;
	if(proc->zero>1){
	  if(proc->fftw)
	    image_oned_fftw_s(data,tppl,tlpi,tipv);
	  else
	    image_oned_fft_s(data,tppl,tlpi,tipv);
	}
      }
    }

    if(fabs(proc->ppe)>0){
      avw_get_vox(&header,&vx,&vy,&vz);
      ppe_phase=(2*M_PI)*((proc->ppe)/vy);
      apply_ppe(data,tppl,tlpi,tipv,ppe_phase);
    }
    if(proc->scsl)scale_slice(data,tppl*tlpi*2,tipv,scsl_array);
    if(proc->pss)pss_reorder(data,tppl*tlpi*2,tipv,pss_array);

    /* Fourier Transform Rows */
    if(proc->fftw)
      image_oned_fftw_r(data,tppl,tlpi,tipv);
    else
      image_oned_fft_r(data,tppl,tlpi,tipv);

    if((proc->ref)&&(!proc->buo))
      epi_phasecor(data,&ref[refindex],tppl,tlpi,tipv);
    if(proc->buo)buono_apply(data,&ref[refindex],tppl,tlpi,tipv);
    if(proc->buov)buono_cor(data,tppl,tlpi,tipv,proc->con);
    if(proc->egrv)entropy_ghost_reduction(data,tppl,tlpi,tipv,proc->seg);
    if(proc->egr){
      if((proc->avvol)||(i==proc->start))
	entropy_ghost_reduction_array(data,tppl,tlpi,tipv,proc->seg,egr);
      entropy_ghost_reduction_apply(data,tppl,tlpi,tipv,proc->seg,egr);
    }
    if(proc->ecf)output_entropy_cost_function(data,tppl,tlpi,tipv,proc->seg);
    if(proc->eiepi)iepi_entropy_ghost_reduction(data,tppl,tlpi,tipv,proc->seg);
    
    /* Half kspace options */
    if(proc->dim<3){
      if(proc->hks==2)hks_herm(data,tppl,tlpi,tipv);
      if(proc->hks==3)hks_pocs(data,tppl,tlpi,tipv);
      if(proc->hks==4)hks_margosian(data,tppl,tlpi,tipv);
    }
    if((proc->hks)&&(proc->rover)){
      fsegment(data,0,tppl,tppl,(tlpi-hlpi)/2,hlpi,tlpi,0,tipv,tipv,2);
      avw_get_vox(&header,&vx,&vy,&vz);
      vy=(vy*tlpi)/hlpi;
      avw_set_vox(&header,vx,vy,vz);
      if(proc->fdf)fdf_set_vox(&fdfheader,vx,vy,vz);
      tlpi=hlpi;      
    }

    /* Fourier Transform Columns */
    if(!proc->ft){
      if(proc->fftw)
	image_oned_ifftw_r(data,tppl,tlpi,tipv);
      else
	image_oned_ifft_r(data,tppl,tlpi,tipv);
    }
    if(proc->ft>=2){
      if(proc->fftw)
	image_oned_fftw_c(data,tppl,tlpi,tipv);
      else
	image_oned_fft_c(data,tppl,tlpi,tipv);
    }

    /* ==================================================================== */
    /*                D R O P O U T S   C A L C U L A T I O N               */
    /* ==================================================================== */
    
    if(proc->odrop){
      if(i>0){
	memcpy(&drop[dropi*tppl*tlpi*tipv*2],data,tppl*tlpi*tipv*2*sizeof(float));
	dropn++;dropi++; 
	if(dropn==11){
	  for(dropc=0;dropc<6;dropc++){
	    calc_dropouts(drop,&dropouts[(dropc+dropl*11)*ipv],tppl,tlpi,tipv,11,dropc);
	  }
	}
	if(dropn>11){
	  calc_dropouts(drop,&dropouts[(dropc+dropl*11)*ipv],tppl,tlpi,tipv,11,dropc);
	  dropc++;
	}
	if(i==proc->start+proc->num-1){
	  for(j=0;j<5;j++){
	    calc_dropouts(drop,&dropouts[(dropc+dropl*11)*ipv],tppl,tlpi,tipv,11,dropc);
	    dropc++; if(dropc==11){dropc=0;dropl++;}
	  }

	  printf("TOTAL NUMBER OF DROPOUTS = %d\n",count_dropouts(dropouts,ipv,dropc+dropl*11));
	}
	if(dropi==11)dropi=0;
	if(dropc==11){dropc=0;dropl++;}
      }
    }

  /* ==================================================================== */
  /*             I M A G E   D O M A I N   P R O C E S S I N G            */
  /* ==================================================================== */

    if((proc->ft)&&(proc->bl)&&(!proc->epi)){
      cencor(data,tppl,tlpi,tipv);
    }
    if(proc->imb)kspace_mask_border(data,tppl,tlpi,tipv,proc->imb);
    if((proc->ft)&&(proc->rover)){
      fsegment(data,tppl/4,tppl/2,tppl,
               0,tlpi,tlpi,0,tipv,tipv,2);
      tppl/=2;
    }
    if((proc->ft)&&(proc->sover)){
      fsegment(data,0,tppl,tppl,
               0,tlpi,tlpi,
	       (int)((tipv-(tipv*SOVER_FRAC))/2),
	       (int)(tipv*SOVER_FRAC),tipv,2);
      tipv=(int)(tipv*SOVER_FRAC);
    }
    if(proc->rot){
      if(transform.rotate){
	rotate(data,tppl,tlpi,tipv,proc->rot);
	swap=tlpi;tlpi=tppl;tppl=swap;
      }
      if(transform.reflect==1){
	reflect_lines(data,tppl,tlpi,tipv);
	reverse_slices(data,tppl,tlpi,tipv);
      }
      if(transform.reflect==2){
	reflect_points(data,tppl,tlpi,tipv);
	reverse_slices(data,tppl,tlpi,tipv);	
      }
    }
    if(proc->fdf){
      reflect_lines(data,tppl,tlpi,tipv);
      reverse_slices(data,tppl,tlpi,tipv);
    }

    /* Reslice data to axial */
    if(proc->resl){
      switch(orient){
      case 1:
	cor2ax(data,tppl,tlpi,tipv);
	swap=tipv;tipv=tlpi;tlpi=swap;
	break;
      case 2:
	sag2ax(data,tppl,tlpi,tipv);
	swap=tppl;tppl=tipv;tipv=swap;
	swap=tipv;tipv=tlpi;tlpi=swap;
	break;
      case 3:
	printf("Warning: Orient is oblique.  Anatomical markers in medx may not be correct\n");
	break;
      }
    }
    
  /* ==================================================================== */
  /*               I M A G E   P H A S E   F U N C T I O N S              */
  /* ==================================================================== */
    if(proc->phrot){
      if((i%proc->phrot)==0){
	phase_rot_init(data,phrot,tppl*tlpi*tipv);
	phase_rot(data,phrot,tppl*tlpi*tipv);
      }
      else {
	phase_rot(data,phrot,tppl*tlpi*tipv);
      }
    }
    if(proc->phrotf)phase_rot(data,phrot,tppl*tlpi*tipv);
    
    /* Modulus, phase, real or imaginary */
    if(proc->out==OUT_MOD)lv_modulus(data,tppl*tlpi*tipv*2);
    if(proc->out==OUT_PHS)phase(data,tppl*tlpi*tipv*2);
    if(proc->out==OUT_REAL)real(data,tppl*tlpi*tipv*2);
    if(proc->out==OUT_IMAG)imag(data,tppl*tlpi*tipv*2); 
    
    /* Resample */
    if((proc->resamp[0]>0)||(proc->resamp[1]>0)||(proc->resamp[2]>0)){
      if(proc->resamp[0]==0)proc->resamp[0]=tppl;
      if(proc->resamp[1]==0)proc->resamp[1]=tlpi;
      if(proc->resamp[2]==0)proc->resamp[2]=tipv;
      resamp(data,tppl,tlpi,tipv,proc->resamp[0],proc->resamp[1],proc->resamp[2]);
      resampx=(float)proc->resamp[0]/(float)tppl;
      resampy=(float)proc->resamp[1]/(float)tlpi;
      resampz=(float)proc->resamp[2]/(float)tipv;
      tppl=proc->resamp[0];
      tlpi=proc->resamp[1];
      tipv=proc->resamp[2];      
    }

    if(proc->otep){
      tep = optimum_tep(data,tppl,tlpi,tipv);
      printf("%f\n",-tep*1e6/sw_readout(infile));
    }
    if(proc->ogroa){
      groa = groa_calc(data,tppl,tlpi,tipv);
      printf("%f\n",groa);
    }
    if(proc->ogrof){
      grof = grof_calc(data,tppl,tlpi,tipv);
      printf("%f\n",grof);
    }

    if(proc->avall){
      for(l=0;l<(unsigned long)(tppl*tlpi*tipv*cplx);l++){
	if(proc->avall==2)
	  average[l]+=(data[l]*data[l]);
	else average[l]+=data[l];
      }
      if(i==proc->start+proc->num-1){
	for(l=0;l<(unsigned long)(tppl*tlpi*tipv*cplx);l++){
	  if(proc->avall==2) data[l]=sqrt(average[l]);
	  else data[l]=average[l]/proc->num;
	}
	oindex=0;
      } else continue;
    }

    if(proc->multicoil){
      if((i+1)%proc->rcvrs==1){
	for(l=0;l<(unsigned long)(tppl*tlpi*tipv*cplx);l++){
	  average[l]=0;
	}
      }
      for(l=0;l<(unsigned long)(tppl*tlpi*tipv*cplx);l++){
	average[l]+=(data[l]*data[l]);
      }
      if((i+1)%proc->rcvrs==0){
	for(l=0;l<(unsigned long)(tppl*tlpi*tipv*cplx);l++){
	  data[l]=sqrt(average[l]);
	}
	oindex=i/proc->rcvrs;
      } else continue;
    }
    
  /* ==================================================================== */
  /*                        S C A L E   I M A G E S                       */
  /* ==================================================================== */

    /* Calculate maximum, minimum */
    if(maxset==0){
      CalcMedian(data,tppl*tlpi*tipv*cplx,&max,&min);
      if(proc->max==0.0){
	if(proc->sint)proc->max=20000;
	else proc->max=max;
      }
      if(proc->scale==0.0){
	proc->scale=proc->max/max;
      }
      proc->max=max*proc->scale;
      proc->min=min*proc->scale;
      avw_set_maxmin(&header,(int)ceil(proc->max*1.275),(int)floor(proc->min));
      if(proc->fdf)fdf_set_maxmin(&fdfheader,proc->max,proc->min);
      maxset=1;
    }

    /* Scale data */
    for(l=0;l<(unsigned long)(tppl*tlpi*tipv*cplx);l++){
      ftmp=(data[l]*proc->scale);
      if(proc->sint){
	if(ftmp>32766)ftmp=32766;
	if(ftmp<-32766)ftmp=-32766;
      }
      sdata[l]=(short)ftmp;
    }

    
  /* ==================================================================== */
  /*   S E T   F D F   P A R A M E T E R S                                */
  /* ==================================================================== */

    if(proc->fdf){
      tvols=proc->num;
      if(proc->avvol)tvols=1;
      if(proc->multicoil)tvols=proc->num/proc->rcvrs;
      if((proc->overx)||(proc->overy)||(proc->overz)||(proc->overv)) {
	if(proc->overx)tppl=proc->overx;
	if(proc->overy)tlpi=proc->overy;
	if(proc->overz)tipv=proc->overz;
	if(proc->overv)tvols=proc->overv;
      }
      if(proc->avall)
	fdf_set_dim(&fdfheader,tppl,tlpi,tipv,1);
      else
	fdf_set_dim(&fdfheader,tppl,tlpi,tipv,tvols);
      if(proc->rot&&transform.rotate){
	avw_get_vox(&header,&vx,&vy,&vz);
	fdf_set_vox(&fdfheader,vy,vx,vz);
      }
      if(proc->resl){
	if(orient==1){
	  avw_get_vox(&header,&vx,&vy,&vz);
	  fdf_set_vox(&fdfheader,vx,vz,vy);
	}
	if(orient==2){
	  avw_get_vox(&header,&vx,&vy,&vz);
	  fdf_set_vox(&fdfheader,vz,vx,vy);
	}
      }
      if((proc->resamp[0]>0)||(proc->resamp[0]>0)||(proc->resamp[0]>0)){
	avw_get_vox(&header,&vx,&vy,&vz);
	vx/=resampx;
	vy/=resampy;
	vz/=resampz;
	fdf_set_vox(&fdfheader,vx,vy,vz);
      }
      fdfheader.array_index=i;
      fdfheader.array_dim=proc->num;
      if(read_procpar(infile,"te ",tmp)){
	fdfheader.te=atof(tmp);
      } else {
	fdfheader.te=0;
      }
      if(read_procpar(infile,"tr ",tmp)){
	fdfheader.tr=atof(tmp);
      } else {
	fdfheader.tr=0;
      }
      fdf_set_params(&fdfheader,infile);
    }
    
    
  /* ==================================================================== */
  /*   S E G M E N T   I M A G E S   A N D   W R I T E   T O   F I L E    */
  /* ==================================================================== */

    /* Short output */
    if(proc->sint){
      /* Segment images */
      if((proc->numx)||(proc->numy)||(proc->numz)){
	if(!proc->numx)proc->numx=tppl;
	if(!proc->numy)proc->numy=tlpi;
	if(!proc->numz)proc->numz=tipv;
	ssegment(sdata,proc->lowx,proc->numx,tppl,
		 proc->lowy,proc->numy,tlpi,
		 proc->lowz,proc->numz,tipv);
	tppl=proc->numx;
	tlpi=proc->numy;
	tipv=proc->numz;
      }
      
      if(!proc->nosave){
	/* Write out volume */
	if(proc->gz){
	  gzseek(ofpz,oindex*tppl*tlpi*tipv*sizeof(short),SEEK_SET);
	  if(gzwrite(ofpz,sdata,sizeof(short)*tppl*tlpi*tipv)
	     !=(long)(tppl*tlpi*tipv*sizeof(short))){
	    printf("Failed to write file %s.img\n",outfile);
	    return(1);
	  } 
	} else {
	  fseek(ofp,oindex*tppl*tlpi*tipv*sizeof(short),SEEK_SET);
	  if(fwrite(sdata,sizeof(short),tppl*tlpi*tipv,ofp)
	     !=(unsigned long)(tppl*tlpi*tipv)){
	    printf("Failed to write file %s.img\n",outfile);
	    return(1);
	  }
	}
      }
      /* Write short FDF file */
      if(proc->fdf){
	fdf_set_type(&fdfheader,"integer",sizeof(short));
	write_fdf_volume(sdata,infile,outfile,&fdfheader);
      }
    }

    /* Float output */
    else {
      /* Segment images */
      if((proc->numx)||(proc->numy)||(proc->numz)){
	if(!proc->numx)proc->numx=tppl;
	if(!proc->numy)proc->numy=tlpi;
	if(!proc->numz)proc->numz=tipv;
	fsegment(data,proc->lowx,proc->numx,tppl,
		 proc->lowy,proc->numy,tlpi,
		 proc->lowz,proc->numz,tipv,cplx);
	tppl=proc->numx;
	tlpi=proc->numy;
	tipv=proc->numz;
      }
      
      /* Write out volume */
      if(!proc->nosave){
	if(proc->gz){
	  gzseek(ofpz,oindex*tppl*tlpi*tipv*cplx*sizeof(float),SEEK_SET);
	  if(gzwrite(ofpz,data,sizeof(float)*tppl*tlpi*tipv*cplx)
	     !=(long)(tppl*tlpi*tipv*cplx*sizeof(float))){
	    printf("Failed to write file %s.img\n",outfile);
	    return(1);
	  }
	} else {
	  fseek(ofp,oindex*tppl*tlpi*tipv*cplx*sizeof(float),SEEK_SET);
	  if(fwrite(data,sizeof(float),tppl*tlpi*tipv*cplx,ofp)
	     !=(unsigned long)(tppl*tlpi*tipv*cplx)){
	    printf("Failed to write file %s.img\n",outfile);
	    return(1);
	  }
	}
      }
      /* Write Float FDF file */
      if(proc->fdf){
	fdf_set_type(&fdfheader,"float",sizeof(float));
	write_fdf_volume(data,infile,outfile,&fdfheader);
      }
    }
  }

  /* ==================================================================== */
  /* ==================================================================== */
  /*       E N D   O F    M A I N   P R O C E S S I N G   L O O P         */
  /* ==================================================================== */
  /* ==================================================================== */
  

  /* ==================================================================== */
  /*   U P D A T E   A N D   W R I T E   A V W   H E A D E R   F I L E    */
  /* ==================================================================== */

  /* Update matrix and voxel sizes */
  tvols=proc->num;
  if(proc->avvol)tvols=1;
  if(proc->multicoil)tvols=proc->num/proc->rcvrs;

  /* Override image matrix size */
  if((proc->overx)||(proc->overy)||(proc->overz)||(proc->overv)) {
    if(proc->overx)tppl=proc->overx;
    if(proc->overy)tlpi=proc->overy;
    if(proc->overz)tipv=proc->overz;
    if(proc->overv)tvols=proc->overv;
  }
  if(proc->avall) {
    avw_set_dim(&header,tppl,tlpi,tipv,1);
  }
  else {
    avw_set_dim(&header,tppl,tlpi,tipv,tvols);
  }
  if(proc->rot&&transform.rotate){
    avw_get_vox(&header,&vx,&vy,&vz);
    avw_set_vox(&header,vy,vx,vz);
  }
  if(proc->resl){
    if(orient==1){
      avw_get_vox(&header,&vx,&vy,&vz);
      avw_set_vox(&header,vx,vz,vy);
      orient=0;
      avw_set_orient(&header,orient);
    }
    if(orient==2){
      avw_get_vox(&header,&vx,&vy,&vz);
      avw_set_vox(&header,vz,vx,vy);
      orient=0;
      avw_set_orient(&header,orient); 
    }
  }
  if((proc->resamp[0]>0)||(proc->resamp[0]>0)||(proc->resamp[0]>0)){
    avw_get_vox(&header,&vx,&vy,&vz);
    vx/=resampx;
    vy/=resampy;
    vz/=resampz;
    avw_set_vox(&header,vx,vy,vz);    
  }

  /* Write AVW header file */
  if(!proc->nosave){
    if(proc->gz){
      if(avw_write_gz(outfile,&header)<0){
	printf("Cannot write header file: %s.hdr\n",outfile);
	return(1);
      }
    } else {
      if(avw_write(outfile,&header)<0){
	printf("Cannot write header file: %s.hdr\n",outfile);
	return(1);
      }
    }
  }


  /* ==================================================================== */
  /*                            C L E A N U P                             */
  /* ==================================================================== */

  /* Close files */
  if(!proc->nosave){
    if(proc->gz)gzclose(ofpz);
    else fclose(ofp);
  }
  /* Free memory */
  free(data);
  free(sdata);
  if(proc->ctep)free(ctep);
  if(proc->nav)free(navigator);
  if(proc->avvol)free(adata);
  if(proc->scsl)free(scsl_array);
  if(proc->pss)free(pss_array);
  if(proc->tabc)free(petable);
  if(proc->phrot)free(phrot);
  if((proc->ref)||(proc->buo))free(ref);
  if(proc->odrop){free(drop);free(dropouts);}
  return(0);
}

int twonup(int n)
{
  int r=2;
  while(r<n)r*=2;
  return(r);
}

