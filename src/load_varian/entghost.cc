/*
  entghost.c: phase correction function using maximum entropy

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "entghost.h"

#include "procfunc.h"
#include "phasefunc.h"
#include "image_fftw.h"

#if defined(HAVE_ISNAN)
 #define Xisnan(s) (isnan(s))
#else
 #define Xisnan(s) (0)
#endif

#include "miscmaths/optimise.h"
using namespace MISCMATHS;

float *egr_data,*egr_tmp;
int egr_ppl,egr_lpi,egr_seg;
float egr_max,egr_min;
int egr_calc;
extern int verbose;

int egr_order=1;

/****************************************************************/
/* Main correction routines                                     */
/****************************************************************/
void entropy_ghost_reduction_array(float *data,int ppl,int lpi,int ipv,
				  int nseg,float *array)
{
  egr_ppl = ppl;
  egr_lpi = lpi;
  egr_seg = nseg;
  egr_tmp = (float *)malloc(ppl*lpi*2*sizeof(float));

  ColumnVector result(1);
  ColumnVector dir(1);
  ColumnVector tol(1);

  result=0;
  dir=1;
  tol=0.01;

  int iter;
  for(int i=0;i<ipv;i++){
    egr_data = &data[i*ppl*lpi*2];
    egr_calc = 0;
    float minent = 1;
    float bestfit = 0;
    for(float init=-0.1;init<=0.1;init+=0.01){
      result = init;
      optimise1d(result,dir,tol,iter,entropy_cost_function,1000,0,0.01);
      float ent = entropy_cost_function(result);
      if(ent<minent){
	minent = ent;
	bestfit = result(1);
      }
    }
    array[i]=bestfit;
  }
  fix_phase_file(array,ipv);
  free(egr_tmp);
}

void entropy_ghost_reduction_apply(float *data,int ppl,int lpi,int ipv,
				  int nseg,float *array)
{
  for(int i=0;i<ipv;i++){
    phase_correct(&data[i*ppl*lpi*2], &data[i*ppl*lpi*2],
		  ppl, lpi, nseg, array[i]);
  }
}

void entropy_ghost_reduction(float *data,int ppl,int lpi,int ipv,int nseg)
{
  float* cal = (float *)malloc(ipv*sizeof(float));
  egr_order = 1;
  entropy_ghost_reduction_array(data,ppl,lpi,ipv,nseg,cal);
  entropy_ghost_reduction_apply(data,ppl,lpi,ipv,nseg,cal);
  egr_order = 0;
  entropy_ghost_reduction_array(data,ppl,lpi,ipv,nseg,cal);
  entropy_ghost_reduction_apply(data,ppl,lpi,ipv,nseg,cal);
  free(cal);
}

float entropy_cost_function(const ColumnVector &input)
{
  phase_correct(egr_data, egr_tmp, egr_ppl, egr_lpi, egr_seg, input(1));
  image_fftw_c(egr_tmp, egr_ppl, egr_lpi, 1, 1);
  lv_modulus(egr_tmp,egr_ppl*egr_lpi*2);
  //calc_gradient(egr_tmp,egr_ppl,egr_lpi);
  if(egr_calc==0){
    egr_max = calc_max(egr_tmp,egr_ppl*egr_lpi);
    egr_min = calc_min(egr_tmp,egr_ppl*egr_lpi);
    egr_calc = 1;
  }
  float e = calc_entropy(egr_tmp,egr_ppl*egr_lpi);
  return(e);
}

/****************************************************************/
/* Main correction routines                                     */
/****************************************************************/


/****************************************************************/
/* Fix phase file                                               */
/****************************************************************/
void fix_phase_file(float *cal, int ipv)
{
  int i,j;
  float med,medd,*wt;
  float cof[5];

  wt=(float *)malloc(ipv*sizeof(float));

  for(i=0;i<ipv;i++)wt[i]=cal[i];

  /* Calculate weighting term */
  med=median(0.5,wt,ipv);
  for(i=0;i<ipv;i++)wt[i]=pow(wt[i]-med,2);
  medd=sqrt(median(0.5,wt,ipv));
  for(i=0;i<ipv;i++){
    if((cal[i]>med+3*medd)||(cal[i]<med-3*medd))
      wt[i]=0;
    else wt[i]=1;
  }

  /* Do fit */
  polynomial_fit(cal,wt,ipv,cof,5);
  for(i=0;i<ipv;i++){
    if(verbose) printf("%f\t",cal[i]);
    for(j=3,cal[i]=cof[4];j>=0;j--)cal[i]=cal[i]*i+cof[j];
    if(verbose) printf("%f\n",cal[i]);
  }
  free(wt); 
}

/****************************************************************/
/* Phase corretion routine                                      */
/****************************************************************/
void phase_correct(float *in, float *out, int ppl, int lpi,int nseg,float ph1)
{
  int i,j,k;
  float a,b,re,im,mod;

  for(k=0;k<ppl*2;k+=2){
    if(egr_order==1){
      a = cos(ph1*(k-ppl)/2);
      b = sin(ph1*(k-ppl)/2);
    } else {
      a = cos(ph1);
      b = sin(ph1);
    }
    mod = magnitude(a,b);
    for(i=0;i<lpi;i+=2*nseg){
      for(j=0;j<nseg;j++){
	out[((i+j)*ppl*2)+k]=in[((i+j)*ppl*2)+k];
	out[((i+j)*ppl*2)+k+1]=in[((i+j)*ppl*2)+k+1];	
      }
      for(j=0;j<nseg;j++){
	re = in[((i+nseg+j)*ppl*2)+k];
	im = in[((i+nseg+j)*ppl*2)+k+1];
	if(mod!=0.0){
	  out[((i+nseg+j)*ppl*2)+k] = (re*a + im*b)/mod;
	  out[((i+nseg+j)*ppl*2)+k+1] = (im*a - re*b)/mod;
	}
      }
    }
  }
}

/****************************************************************/
/* Helper functions                                             */
/****************************************************************/
/* Histogram entropy */
double calc_entropy(float *data,long imsize)
{
  long l;
  double s,t;
  int hist[100];

  calc_hist(data,imsize,egr_min,egr_max,100,hist);
  
  s=0;
  for(l=0;l<100;l++){
    t=((double)hist[l]/(double)imsize)*log((double)hist[l]/(double)imsize);
    if(!Xisnan(t))
      s-=t;
  }
  return(s/log(100.0)); 
}

/* David Atkinson's method */
double calc_entropy1(float *data,long imsize)
{
  long l;
  double s=0,b=0;
  for(l=0;l<imsize;l++) b+= data[l]*data[l];
  b = sqrt(b);
  for(l=0;l<imsize;l++) s-=((data[l]/b)*log(data[l]/b));
  return(s);
}

/* Simple mod sum */
double calc_entropy2(float *data,long imsize)
{
  long l;
  double s=0;
  for(l=0;l<imsize;l++){
    s+=data[l];
  }
  return(s); 
}

/* Row sum then fft then sum */
double calc_entropy3(float *data,int ppl,int lpi)
{
  int x,y;
  float *tmp;
  double sum=0;
  tmp = (float *)malloc(2*lpi*sizeof(float));
 
  for(y=0;y<lpi;y++){
    tmp[y*2]=0;
    tmp[y*2+1]=0;
    for(x=1;x<ppl;x++){
      tmp[y*2]+=data[y*ppl+x];
    }
  }
  image_fftw_r(tmp,lpi,1,1,-1);
  lv_modulus(tmp,lpi*2);

  for(y=0;y<lpi;y++)sum+=tmp[y];
  return(sum);
}

void calc_gradient(float *data,int ppl,int lpi)
{
  long l;
  int x,y;
  float *tmp;
  tmp = (float *)malloc(ppl*lpi*sizeof(float));
 
  for(y=0;y<lpi;y++){
    tmp[y*ppl]=0;
    for(x=1;x<ppl;x++){
      tmp[y*ppl+x]=data[y*ppl+x]-data[y*ppl+x-1];
    }
  }

  for(l=0;l<ppl*lpi;l++)tmp[l]=data[l];
}

float calc_max(float *data,long imsize)
{
  long l;
  float max=0;
  for(l=0;l<imsize;l++){
    if(data[l]>max)max=data[l];
  }
  return(max);
}
float calc_min(float *data,long imsize)
{
  long l;
  float min;

  min=data[0];
  for(l=0;l<imsize;l++){
    if(data[l]<min)min=data[l];
  }
  return(min);
}
void calc_hist(float *data,long imsize,float min,float max,int bins,int *hist)
{
  float binsize;
  int i;
  long l;

  binsize=(max-min)/(bins-1);
  for(i=0;i<bins;i++)hist[i]=0;
  for(l=0;l<imsize;l++){
    if((data[l]>=min)&&(data[l]<=max))
      hist[(int)((data[l]-min)/binsize)]++;
  }
}

/****************************************************************/
/* Misc stuff                                                   */
/****************************************************************/
void output_entropy_cost_function(float *data,int ppl,int lpi,int ipv,int nseg)
{
  egr_ppl = ppl;
  egr_lpi = lpi;
  egr_seg = nseg;
  ColumnVector input(1);

  egr_tmp = (float *)malloc(ppl*lpi*2*sizeof(float));
  for(float ph0=-0.5;ph0<=0.5;ph0+=0.001){
    input = ph0;
    for(int i=0;i<ipv;i++){
      egr_data = &data[i*ppl*lpi*2];
      //float e = entropy_cost_function(input);
    }
  }
  free(egr_tmp);
}

/****************************************************************/
/* IEPI correction                                              */
/****************************************************************/

void iepi_entropy_ghost_reduction(float *data,int ppl,int lpi,int ipv,int nseg)
{
  egr_ppl = ppl;
  egr_lpi = lpi;
  egr_seg = nseg;

  ColumnVector result(4*nseg);
  ColumnVector tol(4*nseg);
  ColumnVector bound(4*nseg);

  result=0;
  tol=0.01;
  bound = 10;

  int iter;
  for(int i=0;i<ipv;i++){
    egr_data = &data[i*ppl*lpi*2];
    result=0;
    egr_calc = 0;
    optimise(result,nseg*4,tol,iepi_cost_function,iter,100,bound);
    iepi_rephase(egr_data,egr_data,ppl,lpi,nseg,result);
  }
}

float iepi_cost_function(const ColumnVector& input)
{
  egr_tmp = (float *)calloc(egr_ppl*egr_lpi*2,sizeof(float));
  iepi_rephase(egr_data, egr_tmp, egr_ppl, egr_lpi, egr_seg, input);
  image_fftw_c(egr_tmp, egr_ppl, egr_lpi, 1, 1);
  lv_modulus(egr_tmp,egr_ppl*egr_lpi*2);
  
  if(egr_calc==0){
    egr_max = calc_max(egr_tmp,egr_ppl*egr_lpi);
    egr_min = calc_min(egr_tmp,egr_ppl*egr_lpi);
    egr_calc = 1;
  }
  float e = calc_entropy(egr_tmp,egr_ppl*egr_lpi);
  free(egr_tmp);
  
  return(e);
}


void iepi_rephase(float *in,float *out,int ppl,int lpi,int nseg,const ColumnVector &phase)
{
  int i,j,k;
  float a,b,re,im,mod;

  for(i=0;i<lpi;i+=2*nseg){
    for(j=0;j<2*nseg;j++){
      for(k=0;k<ppl*2;k+=2){

	a = cos((phase(j+2*nseg+1)*(k-ppl))+phase(j+1));
	b = sin((phase(j+2*nseg+1)*(k-ppl))+phase(j+1));

	re = in[((i+j)*ppl*2)+k];
	im = in[((i+j)*ppl*2)+k+1];
	if((mod=magnitude(a,b))!=0.0){
	  out[((i+j)*ppl*2)+k] = (re*a + im*b)/mod;
	  out[((i+j)*ppl*2)+k+1] = (im*a - re*b)/mod;
	}
      }
    }
  }
}
