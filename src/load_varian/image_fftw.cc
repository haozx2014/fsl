#include "image_fftw.h"

#include <cmath>
#include <cstdlib>

#include "fftw3.h"
#include "procfunc.h"

/* Needed for change from fftw2 to fftw3 */
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])

void fftw_baseline(float *data,int points,int step)
{
  int i,bn;
  float *rarray,*iarray,rmed,imed;

  bn=points/4;
  rarray=(float *)malloc(bn*sizeof(float));
  iarray=(float *)malloc(bn*sizeof(float));
  
  for(i=0;i<bn;i++){
    rarray[i]=data[i*2*step];
    iarray[i]=data[i*2*step+1];
  }
  rmed=median(0.5,rarray,bn);
  imed=median(0.5,iarray,bn);
  for(i=0;i<points;i++){
    data[i*2*step]-=rmed;
    data[i*2*step+1]-=imed;
  }
  free(rarray);
  free(iarray);
}

void image_fftw_r(float *data,int ppl,int lpi,int ipv,int dir)
{
  int r,c,n;
  float scale;
  long offset;
  fftw_complex *in, *out;
  fftw_plan plan;

  n = ppl;
  scale = 1/sqrt(n);
  if((in=(fftw_complex *)malloc(n*sizeof(fftw_complex)))==NULL)return;
  if((out=(fftw_complex *)malloc(n*sizeof(fftw_complex)))==NULL)return;
  if(dir>0)
    plan = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  else
    plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  for (r=0;r<lpi*ipv;r++){
    
    offset = r*ppl*2;

    for(c=0;c<n;c++){
      if(c%2){
	c_re(in[c])=-data[offset+2*c];
	c_im(in[c])=-data[offset+2*c+1];
      } else {
	c_re(in[c])=data[offset+2*c];
	c_im(in[c])=data[offset+2*c+1];	
      }
    }

    fftw_execute(plan);

    for(c=0;c<n;c++){
      if(c%2){
	data[offset+2*c]=-c_re(out[c])*scale;
	data[offset+2*c+1]=-c_im(out[c])*scale;
      }
      else {
	data[offset+2*c]=c_re(out[c])*scale;
	data[offset+2*c+1]=c_im(out[c])*scale;
      }
    }
  }

  free(in);
  free(out);
  fftw_destroy_plan(plan);
}

void image_fftw_c(float *data,int ppl,int lpi,int ipv,int dir)
{
  int i,c,r,n;
  float scale;
  long offset;
  fftw_complex *in, *out;
  fftw_plan plan;

  n = lpi;
  scale = 1/sqrt(n);

  if((in=(fftw_complex *)malloc(n*sizeof(fftw_complex)))==NULL)return;
  if((out=(fftw_complex *)malloc(n*sizeof(fftw_complex)))==NULL)return;
  if(dir>0)
    plan = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  else
    plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  for(i=0;i<ipv;i++){
    for(c=0;c<ppl;c++){
      offset = i*lpi*ppl*2+c*2;
      
      for(r=0;r<n;r++){
	if(r%2){
	  c_re(in[r])=-data[offset+2*ppl*r];
	  c_im(in[r])=-data[offset+2*ppl*r+1];
	} else {
	  c_re(in[r])=data[offset+2*ppl*r];
	  c_im(in[r])=data[offset+2*ppl*r+1];	
	}
      }
      
      fftw_execute(plan);

      for(r=0;r<n;r++){
	if(r%2){
	  data[offset+2*ppl*r]=-c_re(out[r])*scale;
	  data[offset+2*ppl*r+1]=-c_im(out[r])*scale;
	}
	else {
	  data[offset+2*ppl*r]=c_re(out[r])*scale;
	  data[offset+2*ppl*r+1]=c_im(out[r])*scale;
	}
      }
    }
  } 
  free(in);
  free(out);
  fftw_destroy_plan(plan);
}

void image_fftw_s(float *data,int ppl,int lpi,int ipv,float phase,int dir)
{
  int r,c,i,n;
  float scale,f,*ph;
  long offset;
  fftw_complex *in, *out;
  fftw_plan plan;

  n = ipv;
  scale = 1/sqrt(n);

  if((in=(fftw_complex *)malloc(n*sizeof(fftw_complex)))==NULL)return;
  if((out=(fftw_complex *)malloc(n*sizeof(fftw_complex)))==NULL)return;
  if((ph=(float *)malloc(ipv*2*sizeof(float)))==NULL)return;

  if(dir>0)
    plan = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  else
    plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  if(phase!=0.0){
    for(i=0;i<ipv;i++){
      ph[i*2] = sin(phase * (i - ipv / 2));
      ph[i*2 + 1] = cos(phase * (i - ipv / 2));
    }
  }
  
  for(r=0;r<lpi;r++){
    for(c=0;c<ppl;c++){
      offset = r*ppl*2+c*2;
      
      for(i=0;i<n;i++){
	if(i%2){
	  c_re(in[i])=-data[offset+2*ppl*lpi*i];
	  c_im(in[i])=-data[offset+2*ppl*lpi*i+1];
	} else {
	  c_re(in[i])=data[offset+2*lpi*ppl*i];
	  c_im(in[i])=data[offset+2*lpi*ppl*i+1];	
	}
      }

      if(phase!=0.0){
	for(i=0;i<ipv;i++){
	  f = c_re(in[i]) * ph[2*i+1] - c_im(in[i]) * ph[2*i];
	  c_im(in[i]) = c_re(in[i]) * ph[2*i] + c_im(in[i]) * ph[2*i+1];
	  c_re(in[i]) = f;
	}
      }

      fftw_execute(plan);

      for(i=0;i<n;i++){
	if(i%2){
	  data[offset+2*lpi*ppl*i]=-c_re(out[i])*scale;
	  data[offset+2*lpi*ppl*i+1]=-c_im(out[i])*scale;
	}
	else {
	  data[offset+2*lpi*ppl*i]=c_re(out[i])*scale;
	  data[offset+2*lpi*ppl*i+1]=c_im(out[i])*scale;
	}
      }
    }
  } 
  free(in);
  free(out);
  fftw_destroy_plan(plan);
}
