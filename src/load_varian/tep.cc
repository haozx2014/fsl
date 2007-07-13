#include "tep.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "procfunc.h"

#include "miscmaths/optimise.h"
using namespace MISCMATHS;

float *tep_data,*tep_tmp;
int tep_ppl,tep_lpi,tep_ipv;

float optimum_tep(float *data,int ppl,int lpi,int ipv)
{
  int iter;
  tep_ppl = ppl;
  tep_lpi = lpi;
  tep_ipv = ipv;
  tep_data = data;
  tep_tmp = (float *)malloc(ppl*lpi*ipv*sizeof(float));

  ColumnVector result(1);
  ColumnVector dir(1);
  ColumnVector tol(1);
  
  result=tep_first_shift(data,ppl,lpi,ipv);
  dir=1;
  tol=0.01;

  optimise1d(result,dir,tol,iter,tep_cost_func,1000,0,100);
  
  free(tep_tmp);

  return(result(1));
}

void auto_centre_kspace(float *data,int ppl,int lpi,int ipv)
{
  int i,j,k,oline=0,eline=0,shift;
  long l;
  float *even,*odd,emax,omax;
  
  even=(float *)malloc(ppl*sizeof(float));
  odd=(float *)malloc(ppl*sizeof(float));
  
  for(i=0;i<ppl;i++)even[i]=odd[i]=0;
  
  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j+=2){
      for(k=0;k<ppl;k++){
	even[k]+=magnitude(data[i*ppl*lpi*2 + j*ppl*2 + k*2],data[i*ppl*lpi*2 + j*ppl*2 + k*2 + 1]);
	odd[k]+=magnitude(data[i*ppl*lpi*2 + (j+1)*ppl*2 + k*2],data[i*ppl*lpi*2 + (j+1)*ppl*2 + k*2 + 1]);
      }
    }
  }

  emax=omax=0;
  for(k=0;k<ppl;k++){
    if(even[k]>emax){
      emax=even[k];
      eline=k;
    }
    if(odd[k]>omax){
      omax=odd[k];
      oline=k;
    }
  }

  shift = (ppl-oline)-eline;

  if(shift>0){
    for(l=ppl*lpi*ipv*2-1;l>=0;l--){
      if((l-shift)>=0)data[l]=data[l-shift];
      else data[l]=0;
    }
  } else {
    for(l=0;l<ppl*lpi*ipv*2;l++){
      if((l-shift)<ppl*lpi*ipv*2)data[l]=data[l-shift];
      else data[l]=0;
    }
  }

  free(odd);
  free(even);
}

int tep_first_shift(float *data,int ppl,int lpi,int ipv)
{
  int i,j,k,oline=0,eline=0,shift;
  float *even,*odd,emax,omax;
  
  even=(float *)malloc(ppl*sizeof(float));
  odd=(float *)malloc(ppl*sizeof(float));
  
  for(i=0;i<ppl;i++)even[i]=odd[i]=0;
  
  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j+=2){
      for(k=0;k<ppl;k++){
	even[k]+=magnitude(data[i*ppl*lpi*2 + j*ppl*2 + k*2],data[i*ppl*lpi*2 + j*ppl*2 + k*2 + 1]);
	odd[k]+=magnitude(data[i*ppl*lpi*2 + (j+1)*ppl*2 + k*2],data[i*ppl*lpi*2 + (j+1)*ppl*2 + k*2 + 1]);
      }
    }
  }

  emax=omax=0;
  for(k=0;k<ppl;k++){
    if(even[k]>emax){
      emax=even[k];
      eline=k;
    }
    if(odd[k]>omax){
      omax=odd[k];
      oline=k;
    }
  }

  shift = (ppl-oline)-eline;

  free(odd);
  free(even);

  return(shift);
}

float groa_calc(float *data,int ppl,int lpi,int ipv)
{
  Matrix ref(ppl,lpi/2);
  Real y1,y2,y3;
  ColumnVector max(lpi/2);

  // Use central slice
  long offset = (ipv/2)*ppl*lpi;

  // Calculate max intensity offsets
  for(int i=0;i<lpi;i+=2){
    int pos;
    int index=i/2+1;
    for(int j=0;j<ppl;j++)ref(j+1,index)=data[offset+i*ppl+j];
    y2=ref.Column(index).Maximum1(pos);
    if(pos-1>1)
      y1=ref(pos-1,index);
    else y1=y2;
    if(pos+1<=ppl/2)
      y3=ref(pos+1,index);
    else y3=y2;
    if((y1+y3-(2*y2))>0)
      max(index)=pos+((y1-y3)/(2*(y1+y3-(2*y2))));
    else max(index)=pos;
  }
  
  // Initialise Vandermonde matrix
  Matrix V(lpi/2,2);
  for(int i=1;i<=lpi/2;i++) {
    V(i,1)=1;
    V(i,2)=i-1;
  }
  
  // Create weightings matrix
  Matrix W(lpi/2,lpi/2);
  W=0.0;
  for(int i=1;i<=lpi/2;i++) W(i,i)=ref.Column(i).Maximum();
 
  // Linear fit
  Matrix C = (V.t() * W * V).i() * (V.t() * W * max);

  return((float)C(2,1));
}

float grof_calc(float *data,int ppl,int lpi,int ipv)
{
  Matrix ref(ppl,lpi/2);
  Real y1,y2,y3;
  ColumnVector max(lpi/2);
  float grof1,grof2;

  // Use central slice
  long offset = (ipv/2)*ppl*lpi;

  {
  // Calculate max intensity offsets
  for(int i=0;i<lpi;i+=2){
    int pos;
    int index=i/2+1;
    for(int j=0;j<ppl;j++)ref(j+1,index)=data[offset+i*ppl+j];
    y2=ref.Column(index).Maximum1(pos);
    if(pos-1>1)
      y1=ref(pos-1,index);
    else y1=y2;
    if(pos+1<=ppl/2)
      y3=ref(pos+1,index);
    else y3=y2;
    if((y1+y3-(2*y2))>0)
      max(index)=pos+((y1-y3)/(2*(y1+y3-(2*y2))));
    else max(index)=pos;
  }
  
  // Initialise Vandermonde matrix
  Matrix V(lpi/2,2);
  for(int i=1;i<=lpi/2;i++) {
    V(i,1)=1;
    V(i,2)=i-1;
  }
  
  // Create weightings matrix
  Matrix W(lpi/2,lpi/2);
  W=0.0;
  for(int i=1;i<=lpi/2;i++) W(i,i)=ref.Column(i).Maximum();
 
  // Linear fit
  Matrix C = (V.t() * W * V).i() * (V.t() * W * max);
  grof1=C(1,1);

  }

  {
  // Calculate max intensity offsets
  for(int i=0;i<lpi;i+=2){
    int pos;
    int index=i/2+1;
    for(int j=0;j<ppl;j++)ref(j+1,index)=data[offset+(i+1)*ppl+j];
    y2=ref.Column(index).Maximum1(pos);
    if(pos-1>1)
      y1=ref(pos-1,index);
    else y1=y2;
    if(pos+1<=ppl/2)
      y3=ref(pos+1,index);
    else y3=y2;
    if((y1+y3-(2*y2))>0)
      max(index)=pos+((y1-y3)/(2*(y1+y3-(2*y2))));
    else max(index)=pos;
  }
  
  // Initialise Vandermonde matrix
  Matrix V(lpi/2,2);
  for(int i=1;i<=lpi/2;i++) {
    V(i,1)=1;
    V(i,2)=i;
  }
  
  // Create weightings matrix
  Matrix W(lpi/2,lpi/2);
  W=0.0;
  for(int i=1;i<=lpi/2;i++) W(i,i)=ref.Column(i).Maximum();
 
  // Linear fit
  Matrix C = (V.t() * W * V).i() * (V.t() * W * max);
  grof2=C(1,1);
  
  }
  return(grof1-grof2);
}

float tep_cost_func(const ColumnVector &input)
{
  int i,j,k;
  float sum=0;
  ColumnVector A(tep_ppl/2);
  ColumnVector B(tep_ppl/2);

  tep_interpolate(tep_tmp,tep_data,tep_ppl,tep_lpi*tep_ipv,input(1));

  for(i=0;i<tep_ipv;i++){
    for(j=tep_lpi/4;j<tep_lpi*3/4;j+=2){
      for(k=0;k<tep_ppl/2;k++){
	A(k+1) = tep_tmp[i*tep_ppl*tep_lpi+j*tep_ppl+k+tep_ppl/4];
	B(k+1) = tep_tmp[i*tep_ppl*tep_lpi+(j+1)*tep_ppl+k+tep_ppl/4];
      }
      sum += ((A-B).SumAbsoluteValue()/((float)tep_ipv*(float)tep_lpi));
    }
  }
  return((float)sum);
}

void tep_interpolate(float *out,float *in,long points,long lines,float shift)
{
  long i,offset,ini,outi;
  float whole,part,y1,y2;

  for(i=0;i<lines;i++){

    offset = i*points;

    if(i%2){
      whole = floor(-shift);
      part = -shift-whole;
    } else {
      whole = floor(shift);
      part = shift-whole;
    }
    
    outi = 0;
    ini = (long)whole;
    while(outi<points){
      if(ini<0){
	out[offset+outi]=0;
	outi++;
	ini++;
	continue;
      }
      if(ini+1>=points){
	out[offset+outi]=0;
	outi++;
	ini++;
	continue;
      }
      y1 = in[offset+ini];
      y2 = in[offset+ini+1];
      out[offset+outi] = (part*(y2-y1)) + y1;
      outi++;
      ini++;
    }
  }
}

void tep_interpolate_complex(float *data,int points,int lines,float shift)
{
  long i,j,offset,ini,outi;
  float whole,part,y1,y2;

  float *tmp;
  tmp=(float *)malloc(points*2*sizeof(float));

  for(i=0;i<lines;i++){
    offset = i*points*2;

    if(i%2){
      whole = floor(-shift);
      part = -shift-whole;
    } else {
      whole = floor(shift);
      part = shift-whole;
    }
    
    outi = 0;
    ini = (long)whole;
    while(outi<points){
      if(ini<0){
	tmp[outi*2]=0;
	tmp[outi*2+1]=0;
	outi++;
	ini++;
	continue;
      }
      if(ini+1>=points){
	tmp[outi*2]=0;
	tmp[outi*2+1]=0;
	outi++;
	ini++;
	continue;
      }
      y1 = data[offset+ini*2];
      y2 = data[offset+(ini+1)*2];
      tmp[outi*2] = (part*(y2-y1)) + y1;
      y1 = data[offset+ini*2+1];
      y2 = data[offset+(ini+1)*2+1];
      tmp[outi*2+1] = (part*(y2-y1)) + y1;
      outi++;
      ini++;
    }
    for(j=0;j<points*2;j++)data[offset+j]=tmp[j];
  }
  free(tmp);
}
