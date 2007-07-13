#include "image_ft.h"

#include <cmath>
#include <cstdlib>
#include <cstdlib>

void image_oned_fft_r(float *name,int ppl,int lpi,int ipv)
{
  int r;
  for (r=0;r<lpi*ipv;r++){
    centralise(name+r*ppl*2,ppl);
    four1_scaled(name+r*ppl*2-1,ppl,1);
    centralise(name+r*ppl*2,ppl);
  }
}
void image_oned_ifft_r(float *name,int ppl,int lpi,int ipv)
{
  int r;
  for (r=0;r<lpi*ipv;r++){
    centralise(name+r*ppl*2,ppl);
    four1_scaled(name+r*ppl*2-1,ppl,-1);
    centralise(name+r*ppl*2,ppl);
  }
}
void image_oned_fft_c(float *name,int ppl,int lpi,int ipv)
{
  int c,a,size;
  float *xtmp;

  if((xtmp=(float *)malloc(lpi*2*sizeof(float)))==NULL)return;

  size = lpi*ppl*2;
  for(a=0;a<ipv;a++){
    for (c=0;c<ppl*2;c+=2) {
      xtrt(xtmp,0,name+a*size,c,2,lpi,ppl*2,2);
      centralise(xtmp,lpi);
      four1_scaled(xtmp-1,lpi,1);
      centralise(xtmp,lpi);
      xtrt(name+a*size,c,xtmp,0,2,lpi,2,ppl*2);
    }
  }
}
void image_oned_ifft_c(float *name,int ppl,int lpi,int ipv)
{
  int c,a,size;
  float *xtmp;

  if((xtmp=(float *)malloc(lpi*2*sizeof(float)))==NULL)return;

  size = lpi*ppl*2;
  for(a=0;a<ipv;a++){
    for (c=0;c<ppl*2;c+=2) {
      xtrt(xtmp,0,name+a*size,c,2,lpi,ppl*2,2);
      centralise(xtmp,lpi);
      four1_scaled(xtmp-1,lpi,-1);
      centralise(xtmp,lpi);
      xtrt(name+a*size,c,xtmp,0,2,lpi,2,ppl*2);
    }
  }
}
void image_oned_fft_s(float *name,int ppl,int lpi,int ipv)
{
  int a,b;
  float *xtmp;

  if((xtmp=(float *)malloc(ipv*2*sizeof(float)))==NULL)return;
  
  for(a=0;a<lpi;a++){
    for(b=0;b<ppl*2;b+=2){
      xtrt(xtmp,0,name,a*ppl*2+b,2,ipv,ppl*lpi*2,2);
      centralise(xtmp,ipv);
      four1_scaled(xtmp-1,ipv,1);
      centralise(xtmp,ipv);
      xtrt(name,a*ppl*2+b,xtmp,0,2,ipv,2,ppl*lpi*2);
    }
  }
}
void image_oned_fft_sp(float *name,int ppl,int lpi,int ipv,float p)
{
  int a,b,c;
  float *xtmp,*ph,f;

  if((xtmp=(float *)malloc(ipv*2*sizeof(float)))==NULL)return;
  if((ph=(float *)malloc(ipv*2*sizeof(float)))==NULL)return;
  
  if(p!=0.0){
    for(a=0;a<ipv;a++){
      ph[a*2] = sin(p * (a - ipv / 2));
      ph[a*2 + 1] = cos(p * (a - ipv / 2));
    }
  }

  for(a=0;a<lpi;a++){
    for(b=0;b<ppl*2;b+=2){
      xtrt(xtmp,0,name,a*ppl*2+b,2,ipv,ppl*lpi*2,2);
      /* Phase correct */
      if(p!=0){
	for(c=0;c<ipv;c++){
	  f = xtmp[2*c] * ph[2*c+1] - xtmp[2*c+1] * ph[2*c];
	  xtmp[2*c+1] = xtmp[2*c] * ph[2*c] + xtmp[2*c+1] * ph[2*c+1];
	  xtmp[2*c] = f;
	}
      }
      centralise(xtmp,ipv);
      four1_scaled(xtmp-1,ipv,1);
      centralise(xtmp,ipv);
      xtrt(name,a*ppl*2+b,xtmp,0,2,ipv,2,ppl*lpi*2);
    }
  }
}
void image_oned_ifft_s(float *name,int ppl,int lpi,int ipv)
{
  int a,b;
  float *xtmp;

  if((xtmp=(float *)malloc(ipv*2*sizeof(float)))==NULL)return;
  
  for(a=0;a<lpi;a++){
    for(b=0;b<ppl*2;b+=2){
      xtrt(xtmp,0,name,a*ppl*2+b,2,ipv,ppl*lpi*2,2);
      centralise(xtmp,ipv);
      four1_scaled(xtmp-1,ipv,-1);
      centralise(xtmp,ipv);
      xtrt(name,a*ppl*2+b,xtmp,0,2,ipv,2,ppl*lpi*2);
    }
  }
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* four1_scaled: four1 with scaling by 1/sqrt(n) */
/* Based on numerical recipes in C function  */
void four1_scaled(float *data,int nn,int isign)
{
  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta,invn,fn;
  float tempr,tempi;
    
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=2*mmax;
    theta=6.28318530717959/(isign*mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
  fn = nn;
  invn = 1/sqrt(fn);
  for (i=1;i<= nn*2;i++) data[i] *= invn; 

} 
#undef SWAP

void xtrt(register float *result,register int roffs,register float *source,register int soffs,register int nseq,register int nxt,register int inci,register int inco)
{
  register int a,b,c,d;
  for (a=1,c=0,d=0;a<=nxt;a++) {
    for (b=0;b<nseq;b++) 
      *(result+roffs+b+c) = *(source+soffs+b+d);
    c+=inco;
    d+=inci;
  }
}

void centralise(float* data,int pts)
{
  int i;
  for(i=2;i<pts*2;i+=4){
    data[i]*=-1;
    data[i+1]*=-1;
  }
}
