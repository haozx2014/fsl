/*
  procfunc.c: Basic processing functions used by process.c

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "procfunc.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "maxlen.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

int read_petable(char* petable_filename,int** petable)
{
  FILE *fp;
  int tab,i,n,min;

  if((fp=fopen(petable_filename,"r"))==NULL){
    return(0);
  }

  /* Find the size of petable */
  i=0;
  fscanf(fp,"\tt%d =",&tab);
  while(!feof(fp)){
    fscanf(fp,"%d",&tab);
    i++;
  }

  if((*petable=(int *)malloc(i*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(0);
  }

  rewind(fp);
  i=0;
  fscanf(fp,"\tt%d =",&tab);
  while(!feof(fp)){
    fscanf(fp,"%d",&((*petable)[i]));
    i++;
  }
  fclose(fp);

  n=i-1;
  min=(*petable)[0];
  for(i=0;i<n;i++)if((*petable)[i]<min)min=(*petable)[i];
  for(i=0;i<n;i++)(*petable)[i]-=min;
  
  return(n);
}

int tabc(float *data,int ppl,int lpi,int ipv,int* petable)
{
  int i,j,k;
  float *tmp;
  long isize;
  
  isize=ppl*lpi*2;

  if((tmp=(float *)malloc(isize*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }

  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl*2;k++){
	tmp[petable[j]*ppl*2 + k]=data[i*isize + j*ppl*2 + k];
      }
    }
    for(j=0;j<isize;j++)data[(i*isize)+j]=tmp[j];
  }

  free(tmp);
  return(0);
}

int msreorder(float *data,int ppl,int lpi,int ipv,int vols)
{
  /* Reorders multislice data */

  int i,j,k,l,c;
  float *tmp;
  long vsize;

  vsize=ppl*lpi*ipv*2;

  if ((tmp=(float *)malloc(vsize*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }

  for(i=0;i<vols;i++){
    c=0;
    for(j=0;j<ipv;j++){
      for(k=0;k<lpi;k++){
	for(l=0;l<(ppl*2);l++){
	  tmp[c++]=data[(i*vsize)+(k*ppl*ipv*2)+(j*ppl*2)+l];
	}
      }
    }
    for(j=0;j<vsize;j++)data[(i*vsize)+j]=tmp[j];
  }
  free(tmp);
  return(0);
}
int epireorder(float *data,int ppl,int lpi,int ipv,int nseg,int vols)
{
  /* Reorders epi data */

  long i,j,k,l,c,vsize;
  float *tmp;

  vsize=ppl*lpi*ipv*2;

  if ((tmp=(float *)malloc(vsize*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }

  for(i=0;i<vols;i++){
    c=0;
    for(k=0;k<ipv;k++){
      for(l=(lpi/nseg)-1;l>=0;l--){
	for(j=nseg-1;j>=0;j--){
	  copyrow(&data[(i*vsize)+((j*ipv*lpi/nseg)+(k*lpi/nseg)+l)*ppl*2],&tmp[c*ppl*2],ppl,l%2);
	  c++;
	}
      }
    }
    for(j=0;j<vsize;j++)data[(i*vsize)+j]=tmp[j];
  }
  free(tmp);
  return(0);
}
int epikreorder(float *data,int ppl,int lpi,int ipv,int nseg,int kh,int vols)
{
  /* Reorders epik data */

  int olpi,etl_kh,etl_sp;
  long i,j,k,l,c,m,vsize,ovsize;
  float *tmp;
 
  vsize=ppl*lpi*ipv*2;
  olpi = lpi*kh/(kh+nseg-1);
  ovsize=ppl*olpi*ipv*2;
  etl_kh = olpi/kh;
  /* total number of lines for each sparse region */
  etl_sp = (olpi-etl_kh)/2;
  printf("olpi=%d lpi=%d etl_kh=%d etl_sp=%d\n",olpi,lpi,etl_kh,etl_sp);

  if ((tmp=(float *)malloc(vsize*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }

  for(i=0;i<vols;i++){
    c=0;
    for(k=0;k<ipv;k++){
      for(l=(lpi/nseg)-1;l>=0;l--){
	for(j=nseg-1;j>=0;j--){
	  copyrow(&data[(i*vsize)+((j*ipv*lpi/nseg)+(k*lpi/nseg)+l)*ppl*2],&tmp[c*ppl*2],ppl,l%2);
	  c++;
	}
      }
    }

    for(j=0;j<ipv;j++){
      for(k=0;k<etl_sp;k++){
	for(l=0;l<ppl*2;l++){
	  data[i*ovsize+j*ppl*olpi*2+k*ppl*2+l]=tmp[j*ppl*lpi*2+k*ppl*2+l];
	}
      }
      for(k=0;k<etl_kh;k++){
	for(l=0;l<ppl*2;l++){
	  for(m=0;m<nseg;m++){
	    /*data[i*ovsize+j*ppl*olpi*2+(k+etl_sp)*ppl*2+l] += 
	      tmp[j*ppl*lpi*2+(k*nseg+m+etl_sp)*ppl*2+l];*/
	    data[i*ovsize+j*ppl*olpi*2+(k+etl_sp)*ppl*2+l] = 
	      tmp[j*ppl*lpi*2+(k*nseg+m+etl_sp)*ppl*2+l];
	  }
	  /*data[i*ovsize+j*ppl*olpi*2+(k+etl_sp)*ppl*2+l]/=nseg;*/
	}
      }
      for(k=0;k<etl_sp;k++){
	for(l=0;l<ppl*2;l++){
	  data[i*ovsize+j*ppl*olpi*2+(k+etl_sp+etl_kh)*ppl*2+l]
	    =tmp[j*ppl*lpi*2+(k+etl_sp+etl_kh*nseg)*ppl*2+l];
	}
      }
    }
  }
  free(tmp);
  return(0);
}
void copyrow(float *in,float *out,int ppl,int rev)
{
  int i,j;
  if(rev){
    for(i=ppl-1,j=0;i>=0;i--){
      out[j++]=in[(i*2)];
      out[j++]=in[(i*2)+1];
    }
  }
  else{
    for(i=0,j=0;i<ppl;i++){
      out[j++]=in[i*2];
      out[j++]=in[i*2+1];
    }
  }
}
void baseline(float *data,int ppl,int lpi,int ipv,int vols)
{
  /* Calculates the mean pixel intensity of the extremities of k-space
     and subtracts this value from all data points */

  int i,j,k,l;
  long vsize,np;
  float rav,iav;
  int BN=4;

  vsize=ppl*lpi*ipv*2;
  np=(2*lpi*BN)+(2*ppl*BN)-(2*BN*BN);

  for(i=0;i<vols;i++){
    for(j=0;j<ipv;j++){
      rav=iav=0;
      for(k=0;k<BN;k++){
	for(l=0;l<ppl*2;l+=2){
	  rav+=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iav+=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
      }
      for(k=BN;k<(lpi-BN);k++){
	for(l=0;l<BN*2;l+=2){
	  rav+=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iav+=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
	for(l=(ppl-BN)*2;l<ppl*2;l+=2){
	  rav+=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iav+=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
      }
      for(k=(lpi-BN);k<lpi;k++){
	for(l=0;l<ppl*2;l+=2){
	  rav+=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iav+=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
      }
      rav/=(float)np;
      iav/=(float)np;
      for(k=0;k<ppl*lpi*2;k+=2){
	data[i*vsize+j*ppl*lpi*2+k]-=rav;
	data[i*vsize+j*ppl*lpi*2+k+1]-=iav;
      }
    }
  }
}
void median_baseline(float *data,int ppl,int lpi,int ipv,int vols)
{
  /* Calculates the median pixel intensity of the extremities of k-space
     and subtracts this value from all data points */

  int i,j,k,l,bn;
  long vsize,np,n;
  float *rarray,*iarray,rmed,imed;

  bn=lpi/8;

  vsize=ppl*lpi*ipv*2;
  np=(2*lpi*bn)+(2*ppl*bn)-(2*bn*bn);
  rarray=(float *)malloc(np*sizeof(float));
  iarray=(float *)malloc(np*sizeof(float));

  for(i=0;i<vols;i++){
    for(j=0;j<ipv;j++){
      n=0;
      for(k=0;k<bn;k++){
	for(l=0;l<ppl*2;l+=2){
	  rarray[n]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iarray[n++]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
      }
      for(k=bn;k<(lpi-bn);k++){
	for(l=0;l<bn*2;l+=2){
	  rarray[n]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iarray[n++]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
	for(l=(ppl-bn)*2;l<ppl*2;l+=2){
	  rarray[n]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iarray[n++]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
      }
      for(k=(lpi-bn);k<lpi;k++){
	for(l=0;l<ppl*2;l+=2){
	  rarray[n]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iarray[n++]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
      }
      rmed=median(0.5,rarray,n);
      imed=median(0.5,iarray,n);
      for(k=0;k<ppl*lpi*2;k+=2){
	data[i*vsize+j*ppl*lpi*2+k]-=rmed;
	data[i*vsize+j*ppl*lpi*2+k+1]-=imed;
      }
    }
  }
}
void median_baseline_3d(float *data,int ppl,int lpi,int ipv,int vols,int hks)
{
  /* Calculates the median pixel intensity of the extremities of k-space
     and subtracts this value from all data points */

  int i,j,k,l,bn;
  long vsize,n,np;
  float *rarray,*iarray,rmed,imed;

  bn=lpi/8;

  vsize=ppl*lpi*ipv*2;
  np=((2*lpi*bn)+(2*ppl*bn)-(2*bn*bn))*ipv;
  rarray=(float *)malloc(np*sizeof(float));
  iarray=(float *)malloc(np*sizeof(float));

  for(i=0;i<vols;i++){
    n=0;
    for(j=0;j<ipv;j++){
      if(!hks){
	for(k=0;k<bn;k++){
	  for(l=0;l<ppl*2;l+=2){
	    rarray[n]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	    iarray[n++]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	  }
	}
      }
      for(k=bn;k<(lpi-bn);k++){
	for(l=0;l<bn*2;l+=2){
	  rarray[n]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iarray[n++]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
	for(l=(ppl-bn)*2;l<ppl*2;l+=2){
	  rarray[n]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iarray[n++]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
      }
      for(k=(lpi-bn);k<lpi;k++){
	for(l=0;l<ppl*2;l+=2){
	  rarray[n]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l];
	  iarray[n++]=data[i*vsize+j*ppl*lpi*2+k*ppl*2+l+1];
	}
      }
    }
    rmed=median(0.5,rarray,n);
    imed=median(0.5,iarray,n);
    for(k=0;k<vsize;k+=2){
      data[i*vsize+k]-=rmed;
      data[i*vsize+k+1]-=imed;
    }
  }
}
void lv_modulus(float *data,long volsize)
{
  long i,j;
  for(i=0,j=0;i<volsize;i+=2)
    data[j++]=magnitude(data[i+1],data[i]);
}
void phase(float *data,long volsize)
{
  long i,j;
  for(i=0,j=0;i<volsize;i+=2)
    data[j++]=atan2(data[i+1],data[i]);
}
void real(float *data,long volsize)
{
  long i,j;
  for(i=0,j=0;i<volsize;i+=2)
    data[j++]=data[i];
}
void imag(float *data,long volsize)
{
  long i,j;
  for(i=0,j=0;i<volsize;i+=2)
    data[j++]=data[i+1];
}
void CalcMaxMin(float *data,long volsize,float *max,float *min)
{
  long i;
  float lmax,lmin;
  lmax=0;lmin=0;
  for(i=0;i<volsize;i++){
    if(data[i]>lmax)lmax=data[i];
    if(data[i]<lmin)lmin=data[i];
  }
  *max=lmax;
  *min=lmin;
}
void CalcMedian(float *data,long volsize,float *max,float *min)
{
  float *tmp;

  tmp=(float *)malloc(volsize*sizeof(float));
  memcpy(tmp,data,volsize*sizeof(float));

  *max = median(0.999,tmp,volsize);
  *min = median(0.001,tmp,volsize);
  free(tmp);
}
int rotate(float *data,int ppl,int lpi,int nim,int rot)
{
  /* Rotates the complex image matrix by 90 degrees */

  float *tmp;
  int i,k,pos;
  long volsize,j,os;

  volsize=ppl*lpi;
	
  if ((tmp=(float *)malloc(volsize*2*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }

  if(rot>0){
    for(i=0;i<nim;i++){
      os=i*volsize*2;
      pos=0;
      for(j=ppl-1;j>=0;j--){
	for(k=0;k<lpi;k++){
	  tmp[pos++]=data[os + k*ppl*2 + j*2];
	  tmp[pos++]=data[os + k*ppl*2 + j*2 + 1];
	}
      }
      for(j=0;j<volsize*2;j++)data[os+j]=tmp[j];
    }
  }
  if(rot<0){
    for(i=0;i<nim;i++){
      os=i*volsize*2;
      pos=0;
      for(j=0;j<ppl;j++){
	for(k=lpi-1;k>=0;k--){
	  tmp[pos++]=data[os + k*ppl*2 + j*2];
	  tmp[pos++]=data[os + k*ppl*2 + j*2 + 1];
	}
      }
      for(j=0;j<volsize*2;j++)data[os+j]=tmp[j];
    }
  }
  free(tmp);
  return(0);
}

int zerofill(float *data,int dppl,int dlpi,int ppl,int lpi,int ipv,int edgefill)
{
  float *tmp;
  int i,j,k,lppl,llpi;
  long l,in,out;

  if ((tmp=(float *)malloc(ppl*lpi*ipv*2*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  for(l=0;l<ppl*lpi*ipv*2;l++)tmp[l]=0;

  lppl=(ppl-dppl);
  if(edgefill)llpi=0;
  else llpi=(lpi-dlpi)/2;

  for(i=0;i<ipv;i++){
    for(j=llpi;j<llpi+dlpi;j++){
      for(k=lppl;k<lppl+dppl*2;k++){
	in=i*dlpi*dppl*2 + (j-llpi)*dppl*2 +(k-lppl);
	out=i*lpi*ppl*2 + j*ppl*2 + k;
	tmp[out]=data[in];
      }
    }
  }
  for(l=0;l<ppl*lpi*ipv*2;l++)data[l]=tmp[l];
  return(0);
}

int zerofill_3d(float *data,int dppl,int dlpi,int dipv,int ppl,int lpi,int ipv,int edgefill)
{
  float *tmp,sum,max;
  int i,j,k,lppl,llpi,lipv,maxl=0;
  long l,in,out;
  int hw;
  float hanning1,hanning2;

  if ((tmp=(float *)malloc(ppl*lpi*ipv*2*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  for(l=0;l<ppl*lpi*ipv*2;l++)tmp[l]=0;

  lppl=(ppl-dppl);
  lipv = (ipv-dipv)/2;

  /* Hope that this works! */
  if(edgefill){
    max = 0;
    for(i=0;i<dlpi;i++){
      sum = 0;
      for(j=0;j<dipv;j++){
	for(k=0;k<dppl*2;k+=2){
	  in = j*dppl*dlpi*2 + i*dppl*2 + k;
	  sum += magnitude(data[in],data[in+1]);
	}
      }
      if(sum>max){
	max=sum;
	maxl=i;
      }
    }
    llpi=lpi/2-maxl;
    if(llpi<0)llpi=0;
    if(llpi+dlpi>lpi)llpi-=((dlpi+llpi)-lpi);
  }
  else llpi=(lpi-dlpi)/2;

  hw = dlpi/32;
  for(i=lipv;i<lipv+dipv;i++){
    for(j=llpi;j<llpi+dlpi;j++){
      for(k=lppl;k<lppl+dppl*2;k++){
	
	if((j-llpi)<hw) hanning1 = 0.54-0.46*cos((j-llpi)*M_PI/hw);
	else hanning1 = 1.0;
	
	if((k-lppl)<hw) hanning2 = 0.54-0.46*cos((k-lppl)*M_PI/hw);
	else hanning2 = 1.0;

	in=(i-lipv)*dlpi*dppl*2 + (j-llpi)*dppl*2 +(k-lppl);
	out=i*lpi*ppl*2 + j*ppl*2 + k;
	tmp[out]=data[in]*hanning1*hanning2;
      }
    }
  }
  for(l=0;l<ppl*lpi*ipv*2;l++)data[l]=tmp[l];
  return(0);
}

void cpslice(float *in,float *out,long pts)
{
  long i;
  for(i=0;i<pts;i++){
    out[i]=in[i];
  }
}

int reflect_lines(float *data,int ppl,int lpi,int ipv)
{
  float *tmp;
  int i,j,k;
  long l,in,out;

  if ((tmp=(float *)malloc(ppl*lpi*ipv*2*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl*2;k++){
	out = i*lpi*ppl*2 + ((lpi-j-1)*ppl*2) + k;
	in = i*lpi*ppl*2 + j*ppl*2 + k;				   
	tmp[out] = data[in];
      }
    }
  }
  for(l=0;l<ppl*lpi*ipv*2;l++)data[l]=tmp[l];
  free(tmp);
  return(0);
}  

int reflect_points(float *data,int ppl,int lpi,int ipv)
{
  float *tmp;
  int i,j,k;
  long l,in,out;

  if ((tmp=(float *)malloc(ppl*lpi*ipv*2*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl;k++){
	out = i*lpi*ppl*2 + j*ppl*2 + (ppl-k-1)*2;
	in = i*lpi*ppl*2 + j*ppl*2 + k*2;
	tmp[out] = data[in];
	tmp[out+1] = data[in+1];
      }
    }
  }
  for(l=0;l<ppl*lpi*ipv*2;l++)data[l]=tmp[l];
  free(tmp);
  return(0);
}  

float median_greaterthan_zero(double index,float *in,int inlen)
{
  float *tmp,med;
  int i,j=0;

  if((tmp=(float *)malloc(inlen*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(1);
  }
  
  for(i=0;i<inlen;i++){
    if(in[i]!=0)tmp[j++]=in[i];
  }
  med=median(index,tmp,inlen);
  free(tmp);
  return(med);
}

#define CH(a,b) { float d=(a);(a)=(b);(b)=d; }

float median(double index,float *in,int inlen)
{
  unsigned int e=inlen-1,f,g=(unsigned int)(((double)inlen)*index),h,i=0,j;

  while (1) {
    if (e<i+2)
      {
	if ( (e==i+1) && (in[e]<in[i]) ) CH(in[i],in[e]);
	return in[g];
      }
    else
      {
	float c;

	j=(i+e) >> 1;
	CH(in[j],in[i+1]);
	if (in[i]>in[e])   CH(in[i],in[e]);
	if (in[i+1]>in[e]) CH(in[i+1],in[e]);
	if (in[i]>in[i+1]) CH(in[i],in[i+1]);
	f=i+1; h=e; c=in[i+1];
	while (1) {
	  do f++; while (in[f] < c);
	  do h--; while (in[h] > c);
	  if (h < f) break;
	  CH(in[f],in[h]);
	}
	in[i+1]=in[h]; in[h]=c;
	if (h >= g) e=h-1;
	if (h <= g) i=f;
      }
  }
}

#undef CH

double magnitude(double x, double y)
{
  return sqrt((x*x) + (y*y));
}

int find_centre_kspace(float *data,int ppl,int lpi,int ipv)
{
  int i,j,k,max;
  long loc;
  float *tmp,maxf;

  if ((tmp=(float *)malloc(lpi*sizeof(float)))==NULL){
    printf("Malloc failed");
    return(-1);
  }
  for(i=0;i<lpi;i++)tmp[i]=0;
  
  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl;k++){
	loc = i*ppl*lpi + j*ppl + k;
	tmp[j] += magnitude(data[loc*2 + 1],data[loc*2])/(ipv*ppl);
      }
    }
  }
  maxf=0;
  max=-1;
  for(i=0;i<lpi;i+=2){
    if(tmp[i]>maxf){
      maxf=tmp[i];
      max=i;
    }
  }
  free(tmp);
  return(max);
}
void find_centre_kspace_array(float *data,int ppl,int lpi,int ipv,int* array)
{
  int i;

  for(i=0;i<ipv;i++){
    array[i]=find_centre_kspace(&data[i*ppl*lpi*2],ppl,lpi,1);
    if(array[i]<0)array[i]=lpi/2;
  }
}

int looppss_reorder(float *data,int ppl,int lpi,int ipv,int blocks,int index)
{
  int i,j,l,ipb;
  int *order;
  unsigned long imsize;
  float *tmp;

  imsize = ppl*lpi*2;
  ipb = ipv/blocks;

  if((tmp=(float *)malloc(imsize*ipv*sizeof(float)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  if((order=(int *)malloc(ipv*sizeof(int)))==NULL){
    printf("Malloc failed.\n");
    return(1);
  }
  
  for(i=0;i<ipb;i++){
    for(j=0;j<blocks;j++){
      order[i*blocks+j] = j*ipb + i;
    }
  }

  for(i=0;i<ipv;i++){
    l=i+(index*blocks);
    if(l>=ipv)l=l-ipv;
    memcpy(&tmp[order[l]*imsize],&data[i*imsize],imsize*sizeof(float));
  }

  memcpy(data,tmp,imsize*ipv*sizeof(float));
 
  free(tmp);
  free(order);
  return(0);
}
void downsample(float* cplx_data, int ppl, int lpi, int ipv, int downsamp)
{
  int x,y,z;
  int xoff,yoff,yout,xout;
  float *tmp;
  long l;
  xout = ppl/downsamp;
  yout = lpi/downsamp;
  xoff = (ppl-xout)/2;
  yoff = (lpi-yout)/2;
  
  if ((tmp=(float *)malloc(xout*yout*ipv*2*sizeof(float)))==NULL){
    printf("Malloc failed");
    return;
  } 
  
  for(z=0; z<ipv; z++){
    for(y=0; y<yout; y++){
      for(x=0; x<xout; x++){
	tmp[(z*yout*xout + y*xout + x)*2] = 
	  cplx_data[((z*ppl*lpi) + ((y+yoff)*ppl) + x+xoff)*2];
	
	tmp[(z*yout*xout + y*xout + x)*2+1] = 
	  cplx_data[((z*ppl*lpi) + ((y+yoff)*ppl) + x+xoff)*2+1];
      }
    }
  } 

  for(l=0;l<xout*yout*ipv*2;l++)cplx_data[l]=tmp[l];
  free(tmp);
}
void reverse_slices(float *data,int ppl,int lpi,int ipv)
{
  float *tmp;
  int i;
  long l;

  if ((tmp=(float *)malloc(ppl*lpi*ipv*2*sizeof(float)))==NULL){
    printf("Malloc failed");
    return;
  } 
  
  for(i=0;i<ipv;i++){
    for(l=0;l<ppl*lpi*2;l++){
      tmp[(ipv-1-i)*ppl*lpi*2 + l] = data[i*ppl*lpi*2 + l];
    }
  }
  
  for(l=0;l<ppl*lpi*ipv*2;l++)data[l]=tmp[l];
  free(tmp);
}
void apply_ppe(float *data,int ppl,int lpi,int ipv,float angle)
{
  int i,j,k;
  float re,im;

  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl;k++){
	re=data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)];
	im=data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)+1];
	data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)]=re*cos(j*angle/lpi)-im*sin(j*angle/lpi);
	data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)+1]=re*sin(j*angle/lpi)+im*cos(j*angle/lpi);
      }
    }
  }
}
void apply_ppe2(float *data,int ppl,int lpi,int ipv,float angle)
{
  int i,j,k;
  float re,im;

  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl;k++){
	re=data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)];
	im=data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)+1];
	data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)]=re*cos(i*angle/ipv)-im*sin(i*angle/ipv);
	data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)+1]=re*sin(i*angle/ipv)+im*cos(i*angle/ipv);
      }
    }
  }
}
void kpsf_filter(float *data,int ppl,int lpi,int ipv,float a,float b)
{
  int i,j,k,l;
  float f;
  
  for(j=0;j<lpi;j++){

    f=1/(a+((1-a)*pow(b,j)));

    for(i=0;i<ipv;i++){
      for(k=0;k<ppl;k++){
	for(l=0;l<2;l++){
	  data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)+l]/=f;
	}
      }
    }
  }    
}
void custom_filter(float *data,int ppl,int lpi,int ipv,char* infile1)
{
  char infile[MAXLEN],proc_string[MAXLEN];
  int i,j,k,l;
  float f;
  FILE *fp;

  strcpy(infile,infile1);
  strcat(infile,"/filter");
  
  
  for(j=0;j<lpi;j++){

    if((fp=fopen(infile,"rb"))==NULL){
      printf("custom_filter: Cannot open file: %s\n",infile);
      return;
    }
    if(fgets(proc_string,5000,fp)!=NULL){
      sscanf(proc_string, "%f\n", &f);
      for(i=0;i<ipv;i++){
	for(k=0;k<ppl;k++){
	  for(l=0;l<2;l++){
	    data[(i*ppl*lpi*2)+(j*ppl*2)+(k*2)+l]/=f;
	  }
	}
      }
    }    
    else {
      printf("custom_filter: Unexpected end of file: %s\n",infile);
      return;
    }
  }
}
void cencor(float *data, int ppl,int lpi,int ipv)
{
  int z;
  int cx,cy,cz;
  int i,j;
  float re,im;

  cx=ppl/2;
  cy=lpi/2;
  cz=ipv/2;

  for(z=0;z<ipv;z++){
    re=0;
    im=0;
    for(i=-1;i<=1;i+=2){
      for(j=-1;j<=1;j+=2){
	re+=data[(z*ppl*lpi*2)+((cy+i)*ppl*2)+((cx+j)*2)]; 
	im+=data[(z*ppl*lpi*2)+((cy+i)*ppl*2)+((cx+j)*2)+1]; 
      }
    }
    data[(z*ppl*lpi*2)+(cy*ppl*2)+(cx*2)]=re/4; 
    data[(z*ppl*lpi*2)+(cy*ppl*2)+(cx*2)+1]=im/4; 
  }
}
