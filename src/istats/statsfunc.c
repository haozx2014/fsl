#include "statsfunc.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

float calc_entropy(float *data,long imsize,char *mask){
  long l,count;
  float s=0;
  int hist[100];
  float max,min;

  max = calc_max(data,imsize,mask);
  min = calc_min(data,imsize,mask);
  count = calc_hist(data,imsize,mask,min,max,100,hist);
  
  for(l=0;l<100;l++){
    s-=((double)hist[l]/count)*log((double)hist[l]/count);
  }
  return(s/log(100)); 
}

float calc_mean(float *data,long imsize,char *mask){
  long l,c=0;
  float s=0;
  for(l=0;l<imsize;l++){
    if(mask[l]){
      s+=data[l];
      c++;
    }
  }
  return(s/(float)c);
}

float calc_sum(float *data,long imsize,char *mask){
  long l;
  float s=0;
  for(l=0;l<imsize;l++){
    if(mask[l])s+=data[l];
  }
  return(s);
}

float calc_var(float *data,long imsize,char *mask,float mean){
  long l,c=0;
  float s=0;
  for(l=0;l<imsize;l++){
    if(mask[l]){
      s+=(data[l]-mean)*(data[l]-mean);
      c++;
    }
  }
  return(s/(float)c);
}

float calc_max(float *data,long imsize,char *mask){
  long l;
  float max=0;
  for(l=0;l<imsize;l++){
    if(mask[l]){
      if(data[l]>max)max=data[l];
    }
  }
  return(max);
}

float calc_min(float *data,long imsize,char *mask){
  long l;
  float min=0;
  for(l=0;l<imsize;l++){
    if(mask[l]){
      min=data[l];
      break;
    }
  }
  for(l=0;l<imsize;l++){
    if(mask[l]){
      if(data[l]<min)min=data[l];
    }
  }
  return(min);
}

float calc_med(float *data,long imsize,char *mask){
  long l,c=0;
  float *tmp,med;
  
  tmp=(float *)malloc(imsize*sizeof(float));
  for(l=0;l<imsize;l++){
    if(mask[l]){
      tmp[c]=data[l];
      c++;
    }   
  }
  med=nrselect(c/2,c,tmp-1);
  free(tmp);
  return(med);
}

float calc_medv(float *data,long imsize,char *mask,float med){
  long l,c=0;
  float *tmp,medd;
  
  tmp=(float *)malloc(imsize*sizeof(float));
  for(l=0;l<imsize;l++){
    if(mask[l]){
      tmp[c]=(data[l]-med)*(data[l]-med);
      c++;
    }   
  }
  medd=nrselect(c/2,c,tmp-1);
  free(tmp);
  return(medd);
}

float calc_pc(float *data,long imsize,char *mask,int pc){
  long l,c=0;
  float *tmp,ans;
  
  tmp=(float *)malloc(imsize*sizeof(float));
  for(l=0;l<imsize;l++){
    if(mask[l]){
      tmp[c]=data[l];
      c++;
    }   
  }
  ans=nrselect(c*pc/100,c,tmp-1);
  free(tmp);
  return(ans);
}

float calc_mode(float *data,long imsize,char *mask,float max,float min)
{
  long l,c=0,bins;
  float *tmp,ans,q25,q75,w,*hist,hm;
  
  tmp=(float *)malloc(imsize*sizeof(float));
  for(l=0;l<imsize;l++){
    if(mask[l]){
      tmp[c]=data[l];
      c++;
    }
  }
  
  /* Calculate bin size */
  q25 = nrselect(c*0.25,c,tmp-1);
  q75 = nrselect(c*0.75,c,tmp-1);
  w = 2*fabs(q75-q25)/cbrt(c);
  bins = (int)((max-min)/w)+2;
  /*
    printf("max=%f\n",max);
    printf("min=%f\n",min);
    printf("q25=%f\n",q25);
    printf("q75=%f\n",q75);
    printf("w=%f\n",w);
    printf("bins=%ld\n",bins);
  */
  hist=(float *)malloc(bins*sizeof(float));
  for(l=0;l<bins;l++)hist[l]=0;
  for(l=0;l<c;l++)hist[(int)((tmp[l]-min)/w)]++;
  for(l=0,hm=0,ans=0;l<bins;l++){
    if(hist[l]>hm){
      hm=hist[l];
      ans=w*l;
    }
  }
  free(tmp);
  return(ans);
}

long calc_hist(float *data,long imsize,char *mask,float min,float max,int bins,int *hist)
{
  float binsize;
  int i;
  long l,count=0;
  
  binsize=(max-min)/bins;
  for(i=0;i<bins;i++)hist[i]=0;
  for(l=0;l<imsize;l++){
    if(mask[l]){
      if((data[l]>=min)&&(data[l]<=max))
	hist[(int)((data[l]-min)/binsize)]++;
      count++;
    }
  }
  return(count);
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

float nrselect(unsigned long k, unsigned long n, float arr[])
{
  unsigned long i,ir,j,l,mid;
  float a,temp;
  
  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
	  }
      return arr[k];
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
	if (arr[l] > arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
	  }
      arr[l+1]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}
#undef SWAP

