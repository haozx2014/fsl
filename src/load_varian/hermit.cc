/*
  hermit.c: Half k-space reconstruction (Hermitian Conjugation)
  
  Copyright Patricia Figueiredo, Peter Jezzard, Stuart Clare
  FMRIB Centre, University of Oxford.
  
  This program should be considered a beta test version
  and must not be used for any clinical purposes.
  
  Part of...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "hermit.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "image_fftw.h"
#include "phasefunc.h"
#include "procfunc.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

int hks_herm(float *data,int ppl,int lpi,int ipv)
{
  int partial_lpi;
  image_oned_ifftw_r(data,ppl,lpi,ipv);
  hks_centre_kspace(data,ppl,lpi,ipv);
  image_oned_fftw_r(data,ppl,lpi,ipv);
  hks_phase_images(data,ppl,lpi,&partial_lpi,ipv);
  image_oned_ifftw_r(data,ppl,lpi,ipv);
  hks_centre_kspace(data,ppl,lpi,ipv);
  hks_hermitian_conjugation(data,ppl,lpi,partial_lpi,ipv);
  image_oned_fftw_r(data,ppl,lpi,ipv);
  return(0);
}
 
int hks_herm_ref(float *data,float *ref,int ppl,int lpi,int ipv)
{
  int partial_lpi;
  image_oned_ifftw_r(data,ppl,lpi,ipv);
  hks_centre_kspace(data,ppl,lpi,ipv);
  image_oned_fftw_r(data,ppl,lpi,ipv);
  hks_phase_images_ref(data,ref,ppl,lpi,&partial_lpi,ipv);
  image_oned_ifftw_r(data,ppl,lpi,ipv);
  hks_hermitian_conjugation(data,ppl,lpi,partial_lpi,ipv);
  image_oned_fftw_r(data,ppl,lpi,ipv);
  return(0);
}
 
float hks_phase_images(float *data,int ppl,int lpi,int *partial_lpi,int ipv)
{
  int i,j,k,n,offset;
  int centre,extra;
  float *copy,*mod,*phs;
  float aux,cs,sn,re,im,sum;

  if((mod=(float *)malloc(ppl*lpi*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }
  if((phs=(float *)malloc(ppl*lpi*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }
  if((copy=(float *)malloc(2*ppl*lpi*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }

  /* Find number of extra lines */
  j=0;
  for(i=0;i<ipv;i++){
    for(j=lpi/2;j<lpi;j++){
      sum=0;
      for(k=0;k<ppl;k++)sum+=data[i*ppl*lpi*2 + j*ppl*2 + k*2];
      if(sum==0)break;
    }
  }
  centre=lpi/2;
  extra = j-centre;
  *partial_lpi=centre+extra;

  /* Loop over slices */
  for(n=0;n<ipv;n++){

    /* Create offset and copy current image */
    offset=n*ppl*lpi*2;
    for(i=0;i<lpi;i++){
      for(j=0;j<ppl*2;j++){
	copy[i*ppl*2 + j] = data[offset + i*ppl*2 + j];
      }
    }

    /* Hanning filter */
    for(i=centre-extra;i<centre+extra;i++){
      /*aux=0.54-0.46*cos((i-centre+extra)*M_PI/(2*extra));*/
      aux=1;
      for(j=0;j<ppl;j++){
	copy[(i*ppl+j)*2]=aux*copy[(i*ppl+j)*2];
	copy[(i*ppl+j)*2+1]=aux*copy[(i*ppl+j)*2+1];
      }
    }
    for(i=0;i<centre-extra;i++){
      for(j=0;j<ppl;j++){
	copy[(i*ppl+j)*2]=0;
	copy[(i*ppl+j)*2+1]=0;
      }
    }
    for(i=centre+extra;i<lpi;i++){
      for(j=0;j<ppl;j++){
	copy[(i*ppl+j)*2]=0;
	copy[(i*ppl+j)*2+1]=0;
      }
    }

    /* FFTW copy and data */
    image_oned_fftw_c(copy,ppl,lpi,1);
    image_oned_fftw_c(&data[offset],ppl,lpi,1);

    /* Find modulus and phase */
    for(i=0;i<lpi;i++){
      for(j=0;j<ppl;j++){
	mod[i*ppl+j]=magnitude(copy[(i*ppl+j)*2+1],copy[(i*ppl+j)*2]);
	phs[i*ppl+j]=atan2(copy[(i*ppl+j)*2+1],copy[(i*ppl+j)*2]);
      }
    }

    /* Phase rotate */
    for(i=0;i<lpi;i++){
      for(j=0;j<ppl;j++){
	cs=cos(phs[i*ppl+j]);
	sn=sin(phs[i*ppl+j]);
	re=data[offset + i*ppl*2 + j*2];
	im=data[offset + i*ppl*2 + j*2 + 1];
	data[offset + i*ppl*2 + j*2]=(re*cs + im*sn);
	data[offset + i*ppl*2 + j*2 + 1]=(cs*im - re*sn);
      }
    }
    
    /* IFFT data */
    image_oned_ifftw_c(&data[offset],ppl,lpi,1);
  }
  
  /* frees */      
  free(copy);
  free(mod);
  free(phs);
  return(0);
}

float hks_phase_images_ref(float *data,float *ref,int ppl,int lpi,int *partial_lpi,int ipv)
{
  int i,j,k,n,offset;
  int centre,extra;
  float *copy,*mod,*phs;
  float aux,cs,sn,re,im,sum;

  if((mod=(float *)malloc(ppl*lpi*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }
  if((phs=(float *)malloc(ppl*lpi*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }
  if((copy=(float *)malloc(2*ppl*lpi*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }

  /* Find number of extra lines */
  j=0;
  for(i=0;i<ipv;i++){
    for(j=lpi/2;j<lpi;j++){
      sum=0;
      for(k=0;k<ppl;k++)sum+=data[i*ppl*lpi*2 + j*ppl*2 + k*2];
      if(sum==0)break;
    }
  }
  centre=lpi/2;
  extra = j-centre;
  *partial_lpi=centre+extra;

  /* Loop over slices */
  for(n=0;n<ipv;n++){

    /* Create offset and copy current image */
    offset=n*ppl*lpi*2;
    for(i=0;i<lpi;i++){
      for(j=0;j<ppl*2;j++){
	copy[i*ppl*2 + j] = ref[offset + i*ppl*2 + j];
      }
    }

    /* Hanning filter */
    for(i=centre-extra;i<centre+extra;i++){
      aux=0.54-0.46*cos((i-centre+extra)*M_PI/(2*extra));
      for(j=0;j<ppl;j++){
	copy[(i*ppl+j)*2]=aux*copy[(i*ppl+j)*2];
	copy[(i*ppl+j)*2+1]=aux*copy[(i*ppl+j)*2+1];
      }
    }
    for(i=0;i<centre-extra;i++){
      for(j=0;j<ppl;j++){
	copy[(i*ppl+j)*2]=0;
	copy[(i*ppl+j)*2+1]=0;
      }
    }
    for(i=centre+extra;i<lpi;i++){
      for(j=0;j<ppl;j++){
	copy[(i*ppl+j)*2]=0;
	copy[(i*ppl+j)*2+1]=0;
      }
    }

    /* FFT copy and data */
    image_oned_fftw_c(copy,ppl,lpi,1);
    image_oned_fftw_c(&data[offset],ppl,lpi,1);

    /* Find modulus and phase */
    for(i=0;i<lpi;i++){
      for(j=0;j<ppl;j++){
	mod[i*ppl+j]=magnitude(copy[(i*ppl+j)*2+1],copy[(i*ppl+j)*2]);
	phs[i*ppl+j]=atan2(copy[(i*ppl+j)*2+1],copy[(i*ppl+j)*2]);
      }
    }

    /* Phase rotate */
    for(i=0;i<lpi;i++){
      for(j=0;j<ppl;j++){
	cs=cos(phs[i*ppl+j]);
	sn=sin(phs[i*ppl+j]);
	re=data[offset + i*ppl*2 + j*2];
	im=data[offset + i*ppl*2 + j*2 + 1];
	data[offset + i*ppl*2 + j*2]=(re*cs + im*sn);
	data[offset + i*ppl*2 + j*2 + 1]=(cs*im - re*sn);
      }
    }
    
    /* IFFT data */
    image_oned_ifftw_c(&data[offset],ppl,lpi,1);
  }
  
  /* frees */      
  free(copy);
  free(mod);
  free(phs);
  return(0);
}

void hks_hermitian_conjugation(float *data,int ppl,int lpi,int partial_lpi,int ipv)
{
  int n,i,j,p,l,b;
  int pmax,lmax,extra,blend;
  long offset;
  float fblend,real,imag,re,im;
  
  extra = partial_lpi-lpi/2;
  blend=extra/2;
  pmax=(ppl/2);
  lmax=(lpi/2);

  for(n=0;n<ipv;n++){

    offset=n*ppl*lpi*2;

    for(i=0;i<lmax-extra;i++){
      l=(lmax*2)-i;
      if(l>=lpi)continue;
      for(j=0;j<ppl;j++){
	p=(pmax*2)-j;
	if((p>=ppl)||(p<0))continue;
	data[offset + l*ppl*2 + p*2] = data[offset + i*ppl*2 + j*2];
	data[offset + l*ppl*2 + p*2 + 1] = -data[offset + i*ppl*2 + j*2 + 1];
      }
    }

    for(i=lmax-extra,b=1;i<lmax-extra+blend;i++,b++){
      fblend=(float)b/((float)blend+1);
      l=(lmax*2)-i;
      for(j=0;j<ppl;j++){
	p=(pmax*2)-j;
	real=data[offset + i*ppl*2 + j*2];
	imag=-data[offset + i*ppl*2 + j*2 + 1];
	re=data[offset + l*ppl*2 + p*2];
	im=data[offset + l*ppl*2 + p*2 + 1];
	data[offset + l*ppl*2 + p*2] = fblend*re + (1.0-fblend)*real;
	data[offset + l*ppl*2 + p*2 + 1] = fblend*im +  (1.0-fblend)*imag;
      }
    }
  }
}

int hks_centre_kspace(float *data,int ppl,int lpi,int ipv)
{
  int i,j,n,offset;
  int pmax,lmax,p,l;
  long c;
  float maxpro,aux;
  float *projection,*tmp;
  int *shifts;

  if((projection=(float *)malloc(2*lpi*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }
  if((shifts=(int *)malloc(2*ipv*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }
  if((tmp=(float *)malloc(2*ppl*lpi*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  } 

  for(n=0;n<ipv;n++){
    
    /* Calculate offset */
    offset=n*ppl*lpi*2;
    
    /* Find maximum intensity row */
    maxpro=0.0;
    lmax=0;
    for(i=0;i<lpi;i++){
      projection[2*i]=0.0;
      projection[2*i+1]=0.0;
      for(j=0;j<ppl;j++){
	projection[2*i]+=data[offset+i*ppl*2+2*j];
	projection[2*i+1]+=data[offset+i*ppl*2+2*j+1];
      }
      aux=magnitude(projection[2*i],projection[2*i+1]);
      if(aux>maxpro){
	maxpro=aux;
	lmax=i;
      }
    }
    shifts[2*n]=lmax-(lpi/2);

    /* Find maximum intensity column */
    maxpro=0.0;
    pmax=0;
    for(i=0;i<ppl;i++){
      projection[2*i]=0.0;
      projection[2*i+1]=0.0;
      for(j=0;j<lpi;j++){
	projection[2*i]+=data[offset+j*ppl*2+2*i];
	projection[2*i+1]+=data[offset+j*ppl*2+2*i+1];
      }
      aux=magnitude(projection[2*i],projection[2*i+1]);
      if(aux>maxpro){
	maxpro=aux;
	pmax=i;
      }
    }
    shifts[2*n+1]=pmax-(ppl/2);
  }
 
  for(n=0;n<ipv;n++){
    
    /* Calculate offset */
    offset=n*ppl*lpi*2;

    /* Re centre the data  */
    if((shifts[2*n])||(shifts[2*n+1])){
      for(c=0;c<ppl*lpi*2;c++)tmp[c]=data[offset+c];
      for(i=0;i<lpi;i++){
	l=i+shifts[n*2];
	if(l<0)l+=lpi;
	if(l>=lpi)l-=lpi;
	for(j=0;j<ppl;j++){
	  p=j+shifts[n*2+1];
	  if(p<0)p+=ppl;
	  if(p>=ppl)p-=ppl;
	  data[offset+i*ppl*2 + j*2]=tmp[l*ppl*2 + p*2];
	  data[offset+i*ppl*2 + j*2 + 1]=tmp[l*ppl*2 + p*2 + 1];	    
	}
      }
    }
  }
  free(tmp);
  free(projection);
  free(shifts);
  return(0);
}
