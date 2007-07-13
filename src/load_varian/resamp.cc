#include "resamp.h"

#include <cstdlib>
#include <cstdio>
#include <cstdlib>
#include <cmath>

void resamp(float *data, int in_ppl, int in_lpi, int in_ipv,
	    int out_ppl, int out_lpi, int out_ipv)
{
  float *tmp;
  long j;
  int s,l,p;
 
  if((out_ppl*out_lpi*out_ipv)>(in_ppl*in_lpi*in_ipv*2)){
    data=(float *)realloc(data,out_ppl*out_lpi*out_ipv*sizeof(float));
  }

  if ((tmp=(float *)malloc(out_ppl*out_lpi*out_ipv*sizeof(float)))==NULL){
    printf("resamp: Not enough memory\n");
    exit(0);
  }
  
  for(s=0;s<out_ipv;s++){
    for(l=0;l<out_lpi;l++){
      for(p=0;p<out_ppl;p++){
	tmp[s*out_ppl*out_lpi + l*out_ppl + p] = 
	  interp(data,in_ppl,in_lpi,in_ipv,out_ppl,out_lpi,out_ipv,p,l,s);
      }
    }
  }
  for(j=0;j<out_ppl*out_lpi*out_ipv;j++)data[j]=tmp[j];
  free(tmp);
} 

float interp(float *data, int in_ppl, int in_lpi, int in_ipv,
	     int out_ppl, int out_lpi, int out_ipv, int p, int l, int s)
{
  float P,L,S,dp,dl,ds;
  int lp,ll,ls,up,ul,us;
  float temp1, temp2, temp3, temp4, temp5, temp6;
  float v000,v001,v010,v011,v100,v101,v110,v111;
 
  if(in_ppl!=out_ppl){
    P=((float)p/(float)out_ppl)*(float)in_ppl;
    lp=(int)floor(P);
    if(lp+1>=in_ppl) up=lp;
    else up=lp+1;
    dp=P-lp;
  } else {
    lp=up=p;
    dp=0;
  }

  if(in_lpi!=out_lpi){
    L=((float)l/(float)out_lpi)*(float)in_lpi;
    ll=(int)floor(L);
    if(ll+1>=in_lpi) ul=ll;
    else ul=ll+1;
    dl=L-ll;
  } else {
    ll=ul=l;
    dl=0;
  }
  
  if(in_ipv!=out_ipv){
    S=((float)s/(float)out_ipv)*(float)in_ipv;
    ls=(int)floor(S);
    if(ls+1>=in_ipv) us=ls;
    else us=ls+1;
    ds=S-ls;
  } else {
    ls=us=s;
    ds=0;
  }  
  
  v000 = data[ ls*in_lpi*in_ppl + ll*in_ppl + lp];
  v001 = data[ ls*in_lpi*in_ppl + ll*in_ppl + up];
  v010 = data[ ls*in_lpi*in_ppl + ul*in_ppl + lp];
  v011 = data[ ls*in_lpi*in_ppl + ul*in_ppl + up];
  v100 = data[ us*in_lpi*in_ppl + ll*in_ppl + lp];
  v101 = data[ us*in_lpi*in_ppl + ll*in_ppl + up];
  v110 = data[ us*in_lpi*in_ppl + ul*in_ppl + lp];
  v111 = data[ us*in_lpi*in_ppl + ul*in_ppl + up];

  temp1 = (v100 - v000)*ds + v000;
  temp2 = (v101 - v001)*ds + v001;
  temp3 = (v110 - v010)*ds + v010;
  temp4 = (v111 - v011)*ds + v011;
  temp5 = (temp3 - temp1)*dl + temp1;
  temp6 = (temp4 - temp2)*dl + temp2;
  return (temp6 - temp5)*dp + temp5;
}
