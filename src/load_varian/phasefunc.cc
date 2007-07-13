/*
  phasefunc.c: Series of phase correction functions used by process.c
               (mainly for EPI data)

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "phasefunc.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "image_fftw.h"
#include "procfunc.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

#define NCOF 5

extern int buo_limit1, buo_limit2, verbose;

/****************************************************************/
/* Self ref phase correction                                    */
/* Based on method by Buonocore et.al. MRM 38,89-100            */
/****************************************************************/

void buono_calib(float *data,float *cal,int ppl,int lpi,int ipv,int constrain)
{
  int i;
  for(i=0;i<ipv;i++){
    buono(&data[i*ppl*lpi*2],&cal[i*ppl*2],ppl,lpi,constrain);
  }
  if((ipv>10)&&(constrain==1))fix_cal_file(cal,ppl,ipv);
}
void buono_apply(float *data,float *cal,int ppl,int lpi,int ipv)
{
  int i;
  for(i=0;i<ipv;i++){
    buono_correct(&data[i*ppl*lpi*2],&cal[i*ppl*2],ppl,lpi);
  }

}
void buono_cor(float *data,int ppl,int lpi,int ipv,int constrain)
{
  int i;
  float *cal;
  cal=(float *)malloc(ppl*2*sizeof(float));

  for(i=0;i<ipv;i++){
    buono(&data[i*ppl*lpi*2],cal,ppl,lpi,constrain);
    buono_correct(&data[i*ppl*lpi*2],cal,ppl,lpi);
  }
  if(verbose) printf("\n");
  free(cal);
}
void fix_cal_file(float *cal,int ppl,int ipv)
{
  int i;
  float *zero,*first,*cen,med,medd,*zfit,*ffit;
  float cof[NCOF];
  int j;

  zero=(float *)malloc(ipv*sizeof(float));
  first=(float *)malloc(ipv*sizeof(float));
  cen=(float *)malloc(ipv*sizeof(float));
  zfit=(float *)malloc(ipv*sizeof(float));
  ffit=(float *)malloc(ipv*sizeof(float));

  for(i=0;i<ipv;i++){
    cen[i]=cal[i*ppl*2+ppl];
    zero[i]=cal[i*ppl*2];
    first[i]=cal[i*ppl*2+2]-cal[i*ppl*2];
  }

  /* Calculate weighting term */
  med=median(0.5,cen,ipv);
  for(i=0;i<ipv;i++)cen[i]=pow(cen[i]-med,2);
  medd=sqrt(median(0.5,cen,ipv));
  for(i=0;i<ipv;i++){
    if((cal[i*ppl*2+ppl]>med+3*medd)||(cal[i*ppl*2+ppl]<med-3*medd))
      cen[i]=0;
    else cen[i]=1;
  }

  // First order term
  polynomial_fit(first,cen,ipv,cof,NCOF);
  for(i=0;i<ipv;i++){
    for(j=NCOF-2,ffit[i]=cof[NCOF-1];j>=0;j--)ffit[i]=ffit[i]*i+cof[j];
  }
  // Zeroth order term
  polynomial_fit(zero,cen,ipv,cof,NCOF);
  for(i=0;i<ipv;i++){
    for(j=NCOF-2,zfit[i]=cof[NCOF-1];j>=0;j--)zfit[i]=zfit[i]*i+cof[j];
  }

  // Find median error and deviation
  //for(i=0;i<ipv;i++)cen[i]=zero[i]-zfit[i];
  //med=median(0.5,cen,ipv);
  //for(i=0;i<ipv;i++)cen[i]=pow(cen[i]-med,2);
  //medd=sqrt(median(0.5,cen,ipv));
  //for(i=0;i<ipv;i++)cen[i]=zero[i]-zfit[i];

  // Regenerate cal
  for(i=0;i<ipv;i++){
    //if((cen[i]>(med+3*medd))||(cen[i]<(med-3*medd))){
    for(j=0;j<ppl;j++)cal[i*ppl*2+j*2]=zfit[i]+ffit[i]*j;
    //}
  }

  free(zero);
  free(first);
  free(zfit);
  free(ffit);
  free(cen);
}
void buono(float *data,float *cal,int ppl,int lpi,int constrain)
{
  long i;
  int j,k,size,n1,n2;
  float *phase,*mod,*odd,*even,a1,a2,a,b,d,m,n;
  float *phs1,*phs2,*wt;
  float max,x,y,sx,sy,sxy,sx2,w,sum;

  float pp,qq;
	
  size=ppl*lpi;

  phase=(float *)malloc(ppl*2*sizeof(float));
  mod=(float *)malloc(ppl*2*sizeof(float));
  odd=(float *)malloc(size*2*sizeof(float));
  even=(float *)malloc(size*2*sizeof(float));

  if((buo_limit1)&&(buo_limit2)){
    n1=buo_limit1;
    n2=buo_limit2;
  } else {
    /* Find max brightness line */
    for(i=0;i<size*2;i++){
      odd[i]=data[i];
    }

    image_oned_fftw_c(odd,ppl,lpi,1);
    
    n1=n2=0;
    for(i=0,pp=qq=0;i<lpi;i++){
      for(j=0,sum=0;j<ppl*2;j+=2)sum+=magnitude(odd[(i*ppl*2)+j],odd[(i*ppl*2)+j+1]);
      pp+=(sum*i);
      qq+=sum;
    }
    n1=(int)(pp/qq)-4;
    n2=(int)(pp/qq)+4;
    if(n1<0)n1=0;
    if(n2>=lpi)n2=lpi-1;
  }

  /* Split up odd and even echos and ft image */
  for(i=0;i<size*2;i++){
    odd[i]=0;
    even[i]=0;
  }
  for(i=0;i<lpi;i+=2){
    for(j=0;j<ppl*2;j++){
      even[(i*ppl*2)+j]=data[(i*ppl*2)+j];
    }
    for(j=0;j<ppl*2;j++){
      odd[((i+1)*ppl*2)+j]=data[((i+1)*ppl*2)+j];
    }
  }
  image_oned_fftw_c(even,ppl,lpi,1);
  image_oned_fftw_c(odd,ppl,lpi,1);

  for(i=0,max=0;i<size*2;i+=2){
    if((m=magnitude(even[i+1],even[i]))>max)max=m;
  }
  max/=10;

  for(i=0;i<ppl*2;i++)cal[i]=0;
  for(k=0;k<ppl*2;k++)mod[k]=0;

  /* Calculate phase difference between odd and even */
  for(i=n1;i<n2;i++){
    for(k=0;k<ppl*2;k+=2){
      if((m=magnitude(even[(i*ppl*2)+k+1],even[(i*ppl*2)+k]))>max){
	mod[k]+=m;
	mod[k+1]+=1.0;
	a1=atan2(even[(i*ppl*2)+k+1],even[(i*ppl*2)+k]);
	a2=atan2(odd[(i*ppl*2)+k+1],odd[(i*ppl*2)+k]);
	a=a1-a2;
	if(a>M_PI)a-=(2*M_PI);
	if(a<(-1*M_PI))a+=(2*M_PI);
	phase[k]=a;
	cal[k+1]+=1.0;
      }
      else phase[k]=0;
    }
    unwrap2(phase,ppl);
    for(k=0;k<ppl*2;k+=2){
      cal[k]+=phase[k];
    }
  }
  for(i=0;i<ppl*2;i+=2){
    if(cal[i+1]>0)cal[i]/=(2*cal[i+1]);
    if(mod[i+1]>0)mod[i]/=mod[i+1];
  }

  if(constrain){
    switch(constrain){
    case 1:
      sxy=sx=sy=sx2=n=0;
      for(i=0;i<ppl*2;i+=2){
	x=i/2;
	y=cal[i];
	w=mod[i];
	n+=w;
	sxy+=(x*y*w);
	sx+=(x*w);
	sy+=(y*w);
	sx2+=(x*x*w);
      }
      if(n>0){
	d = (sx2-((sx*sx)/n));
	if(d>0){
	  b = (sxy-((sx*sy)/n))/d;
	} else b = 0.0;
	a = (sy/n)-(b*(sx/n));
      } else
	a = b = 0.0;

      if(verbose) printf("%f\t%f\t",a,b);
      
      for(i=0;i<ppl*2;i+=2){
	x=i/2;
	cal[i]=a+(b*x);
      }
      break;
    case 5:
      phs1=(float *)malloc(ppl*sizeof(float));
      phs2=(float *)malloc(ppl*sizeof(float));
      wt=(float *)malloc(ppl*sizeof(float));
      for(i=0;i<ppl*2;i+=2){
	phs1[i/2]=0.0;
	phs2[i/2]=cal[i];
	wt[i/2]=mod[i];
      }
      polyfit(phs1,phs2,wt,ppl);
      for(i=0;i<ppl*2;i+=2){
	cal[i]=phs2[i/2];
      }
      free(phs1);
      free(phs2);
      free(wt);
      break;
    default:
      printf("Only fit of order 1 or 5 available\n");
      exit(1);
    }
  }
  free(phase);
  free(even);
  free(odd);
  free(mod);
}
void buono_correct(float *data,float* cal,int ppl,int lpi)
{
  int i,j;
  float a,r1,i1;
  for(i=0;i<lpi;i+=2){
    for(j=0;j<ppl*2;j+=2){
      r1=data[(i*ppl*2)+j];
      i1=data[(i*ppl*2)+j+1];
      a=cal[j];
      data[(i*ppl*2)+j]=(r1*cos(a))+(i1*sin(a));
      data[(i*ppl*2)+j+1]=(i1*cos(a))-(r1*sin(a));
    }
    for(j=0;j<ppl*2;j+=2){
      r1=data[((i+1)*ppl*2)+j];
      i1=data[((i+1)*ppl*2)+j+1];
      a=cal[j];
      data[((i+1)*ppl*2)+j]=(r1*cos(a))-(i1*sin(a));
      data[((i+1)*ppl*2)+j+1]=(i1*cos(a))+(r1*sin(a));
    }
  }
}
void unwrap(float *data,int points)
{
  int i,j,cen;

  cen=points/2;
  for(i=cen+1;i<points;i++){
    if((data[i]-data[i-1])>M_PI){
      for(j=i;j<points;j++)data[j]-=(2*M_PI);
    }
    if((data[i-1]-data[i])>M_PI){
      for(j=i;j<points;j++)data[j]+=(2*M_PI);
    }
  }
  for(i=cen-1;i>=0;i--){
    if((data[i]-data[i+1])>M_PI){
      for(j=i;j>=0;j--)data[j]-=(2*M_PI);
    }
    if((data[i+1]-data[i])>M_PI){
      for(j=i;j>=0;j--)data[j]+=(2*M_PI);
    }
  }
}
void unwrap2(float *data,int ppl)
{
  int i,j;

  for(i=ppl+2;i<ppl*2;i+=2){
    if((data[i]-data[i-2])>M_PI){
      for(j=i;j<ppl*2;j+=2)data[j]-=(2*M_PI);
    }
    if((data[i-2]-data[i])>M_PI){
      for(j=i;j<ppl*2;j+=2)data[j]+=(2*M_PI);
    }
  }
  for(i=ppl-2;i>=0;i-=2){
    if((data[i]-data[i+2])>M_PI){
      for(j=i;j>=0;j-=2)data[j]-=(2*M_PI);
    }
    if((data[i+2]-data[i])>M_PI){
      for(j=i;j>=0;j-=2)data[j]+=(2*M_PI);
    }
  }
}

/****************************************************************/
/* Ref scan phase correction                                    */
/****************************************************************/

void epi_phasecor(float *data,float *ref,int ppl,int lpi,int ipv)
{
  int j;
  float r,i,r1,i1,mod;

  for(j=0;j<ppl*lpi*ipv;j++){
    r=data[j*2];i=data[j*2+1];
    r1=ref[j*2];i1=ref[j*2+1];
    if((mod=magnitude(r1,i1))!=0.0){
      data[j*2] = (r*r1 + i*i1)/mod;
      data[j*2+1] = (r1*i - r*i1)/mod;
    }
  }
}
int calc_full_ref(float *ref,int ppl,int lpi,int ipv,int nseg)
{
  int i,j,k,l;
  long h,segsize;
  float *tmp;

  segsize=ppl*lpi/nseg;

  if ((tmp=(float *)malloc(segsize*ipv*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  for(h=0;h<segsize*ipv;h++)tmp[h]=ref[h];

  for(i=0;i<ipv;i++){
    for(j=0;j<lpi/nseg;j++){
      for(k=0;k<nseg;k++){
	for(l=0;l<ppl*2;l++){
	  ref[i*segsize + j*nseg*ppl*2 + k*ppl*2 + l]=tmp[i*segsize + j*ppl*2 + l];
	}
      }
    }
  }

  return(0);
}
int calc_phase_cor(float *ref,int ppl,int lpi,int ipv,int nseg,int order,int* cenline)
{
  int i,j,k,l;
  long os;
  float *phs1,*phs2,*mod1,*mod2,*conphs1,*conphs2;

  if ((phs1=(float *)malloc(ppl*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  if ((phs2=(float *)malloc(ppl*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  if ((mod1=(float *)malloc(ppl*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  if ((mod2=(float *)malloc(ppl*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  if ((conphs1=(float *)malloc(ipv*ppl*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }
  if ((conphs2=(float *)malloc(ipv*ppl*sizeof(float)))==NULL){
    printf("Malloc failed (process)");
    return(1);
  }


  for(i=0;i<ipv;i++){
    /* Central line must be even */
    if(cenline[i]%2)cenline[i]--;

    os=i*(lpi/nseg)*ppl*2;

    for(k=0;k<ppl*2;k+=2){
      phs1[k/2]=atan2(ref[os + cenline[i]*ppl*2 + k + 1], ref[os + cenline[i]*ppl*2 + k]);
      mod1[k/2]=magnitude(ref[os + cenline[i]*ppl*2 + k + 1], ref[os + cenline[i]*ppl*2 + k]);
    }
    for(k=0;k<ppl*2;k+=2){
      phs2[k/2]=atan2(ref[os + (cenline[i]+1)*ppl*2 + k + 1], ref[os + (cenline[i]+1)*ppl*2 + k]);
      mod2[k/2]=magnitude(ref[os + (cenline[i]+1)*ppl*2 + k + 1], ref[os + (cenline[i]+1)*ppl*2 + k]);
    }
     
    switch(order){
    case 1:
      linfit(phs1,phs2,mod1,ppl);
      break;
    case 5:
      polyfit(phs1,phs2,mod1,ppl);
      break;
    default:
      printf("Only fit of order 1 or 5 available\n");
      exit(1);
    }
    for(k=0;k<ppl;k++){
      conphs1[i*ppl+k]=phs1[k];
      conphs2[i*ppl+k]=phs2[k];
    }
  }
  
  /* Create new ref scan */
  for(i=0;i<ipv;i++){
    os=i*lpi*ppl*2;
    for(j=0;j<lpi;j+=(nseg*2)){
      for(l=0;l<nseg;l++) {
	for(k=0;k<ppl*2;k+=2){
	  ref[os + (j+l)*ppl*2 + k]=cos(conphs1[i*ppl+k/2]);
	  ref[os + (j+l)*ppl*2 + k + 1]=sin(conphs1[i*ppl+k/2]);
	}
      }
      for(l=nseg;l<nseg*2;l++){
	for(k=0;k<ppl*2;k+=2){
	  ref[os + (j+l)*ppl*2 + k]=cos(conphs2[i*ppl+k/2]);
	  ref[os + (j+l)*ppl*2 + k + 1]=sin(conphs2[i*ppl+k/2]);
	}
      }
    }
  }
  
  free(phs1);
  free(phs2);
  free(mod1);
  free(mod2);
  free(conphs1);
  free(conphs2);
  return(0);
}
void linfit(float *phs1,float *phs2,float *mod,int ppl)
{
  float *dif,*wt,zeroth,first;
  float sum,d,f,w,max;
  int i;
  int pt=0;

  dif=(float *)malloc(ppl*sizeof(float));
  wt=(float *)malloc(ppl*sizeof(float));

  for(i=0;i<ppl;i++){
    dif[i]=phs2[i]-phs1[i];
    wt[i]=mod[i];
  }

  f=0;
  sum=0;
  for(i=1;i<ppl;i++){
    d=dif[i]-dif[i-1];
    if((d>(M_PI/4))||(d<(-M_PI/4)))w=0;
    else w=wt[i];
    f+=(w*d);
    sum+=w;
  }

  if(sum)f/=sum;

  for(i=0,max=0;i<ppl;i++){
    if(wt[i]>max){
      max=wt[i];
      pt=i;
    }
  }
  max=dif[pt];

  for(i=0;i<ppl;i++){
    d=f*(i-pt)+max;
    while((d-dif[i])>M_PI)dif[i]+=(2*M_PI);
    while((d-dif[i])<-M_PI)dif[i]-=(2*M_PI);
  }
 
  linear_regression(dif,wt,ppl,&zeroth,&first);
  /*printf("%f\t%f\n",zeroth,first);*/

  for(i=0;i<ppl;i++){
    phs1[i]=0;
    phs2[i]=zeroth+(i*first);
  }
}
void linear_regression(float *phs,float *mod,int ppl,float *a,float *b)
{
  float sxy,sx,sy,sx2,n;
  float x,y,w;
  int i;

  sxy=sx=sy=sx2=n=0;
  for(i=0;i<ppl;i++){
    x=i;
    y=phs[i];
    w=mod[i];
    n+=w;
    sxy+=(x*y*w);
    sx+=(x*w);
    sy+=(y*w);
    sx2+=(x*x*w);
  }
  if(n){
    *b=(sxy-((sx*sy)/n))/(sx2-((sx*sx)/n));
    *a=(sy/n)-(*b*(sx/n));
  }
  else {
    *b=0;
    *a=0;
  }
}

void polyfit(float* phs1,float* phs2,float* mod,int ppl)
{
  float *dif,*wt,cof[5];
  float sum,d,f,w,max;
  int i,j;
  int pt=0;

  dif=(float *)malloc(ppl*sizeof(float));
  wt=(float *)malloc(ppl*sizeof(float));

  for(i=0;i<ppl;i++){
    dif[i]=phs2[i]-phs1[i];
    wt[i]=mod[i];
  }

  f=0;
  sum=0;
  for(i=1;i<ppl;i++){
    d=dif[i]-dif[i-1];
    if((d>(M_PI/4))||(d<(-M_PI/4)))w=0;
    else w=wt[i];
    f+=(w*d);
    sum+=w;
  }

  if(sum)f/=sum;

  for(i=0,max=0;i<ppl;i++){
    if(wt[i]>max){
      max=wt[i];
      pt=i;
    }
  }
  max=dif[pt];
  
  for(i=0;i<ppl;i++){
    d=f*(i-pt)+max;
    while((d-dif[i])>M_PI)dif[i]+=(2*M_PI);
    while((d-dif[i])<-M_PI)dif[i]-=(2*M_PI);
  }

  cof[0]=cof[1]=cof[2]=cof[3]=cof[4]=0;

  polynomial_fit(dif,wt,ppl,cof,5);

  for(i=0;i<ppl;i++){
    phs1[i]=0;
    phs2[i]=cof[4];
    for(j=3;j>=0;j--)phs2[i]=phs2[i]*i+cof[j];
  }
}


#include "newmat.h"
using namespace NEWMAT;

void polynomial_fit(float *data,float *wt,int n,float* cof,int ncof)
{
  // Initialise Y matrix
  Matrix Y(n,1);
  for(int i=1;i<=n;i++) Y(i,1)=data[i-1];

  // Initialise W matrix
  Matrix W(n,n);
  W=0.0;
  for(int i=1;i<=n;i++) W(i,i) = wt[i-1];

  // Initialise Vandermonde matrix
  Matrix V(n,ncof);
  for(int i=1;i<=n;i++) V(i,1) = 1;
  for(int j=2;j<=ncof;j++){
    for(int i=1;i<=n;i++){
      V(i,j) = V(i,j-1)*(i-1);
    }
  }
  
  Matrix C = (V.t() * W * V).i() * (V.t() * W * Y);

  for(int i=1;i<=ncof;i++) cof[i-1] = C(i,1);
}


void create_ref_scan(float *ref,int ppl,int lpi,int ipv,int nseg,float ph0,float ph1)
{
  int i,j,k,l;
  for(i=0;i<ipv;i++){
    for(j=0;j<lpi;j+=2*nseg){
      for(l=0;l<nseg;l++){
	for(k=0;k<ppl*2;k++){
	  ref[i*ppl*lpi*2+(j+l)*ppl*2+k]=0;
	}
      }
      for(l=0;l<nseg;l++){
	for(k=0;k<ppl*2;k+=2){
	  ref[i*ppl*lpi*2+(j+nseg+l)*ppl*2+k]=cos(ph0+(ph1*(k-ppl)/2));
	  ref[i*ppl*lpi*2+(j+nseg+l)*ppl*2+k+1]=sin(ph0+(ph1*(k-ppl)/2));
	}
      }
    }
  }
}

/****************************************************************/
/* Phase encoded ref scan phase correction                      */
/****************************************************************/

int hu_cal(float *data,float *ref,int ppl,int lpi,int ipv,int constrain)
{
  int x,y,p;
  double r1,r2,i1,i2,m1,m2;
  float *mod,*phs1,*phs2;

  /* Calculate phase correction */
  for(y=1;y<(lpi*ipv);y+=2){
    if(constrain){
      if ((phs1=(float *)malloc(ppl*sizeof(float)))==NULL){
	printf("Malloc failed (process)");
	return(1);
      }
      if ((phs2=(float *)malloc(ppl*sizeof(float)))==NULL){
	printf("Malloc failed (process)");
	return(1);
      }
      if ((mod=(float *)malloc(ppl*sizeof(float)))==NULL){
	printf("Malloc failed (process)");
	return(1);
      }
      for(x=0;x<ppl*2;x+=2){
	p=(y*ppl*2)+x;
	mod[x/2]=magnitude(data[p],data[p+1]);
	phs1[x/2]=atan2(data[p+1],data[p]);
	phs2[x/2]=atan2(ref[p+1],ref[p]);
      }
      linfit(phs1,phs2,mod,ppl);
      for(x=0;x<ppl;x++){
	p=(y*ppl*2)+(x*2);
	/*Store cos(ph2-ph1) in ref[p]   */
	/*      sin(ph2-ph1) in ref[p+1] */
	ref[p]=cos(phs2[x]);
	ref[p+1]=sin(phs2[x]);
      }
    }
    else{
      for(x=0;x<ppl*2;x+=2){
	p=(y*ppl*2)+x;
	r1=(double)data[p];i1=(double)data[p+1];
	r2=(double)ref[p];i2=(double)ref[p+1];
	m1=magnitude(r1,i1);
	m2=magnitude(r2,i2);
			
	/*Store cos(ph2-ph1) in ref[p]   */
	/*      sin(ph2-ph1) in ref[p+1] */
	if((m1*m2)>0){
	  ref[p]=(float)((r1*r2)+(i1*i2))/(m1*m2);
	  ref[p+1]=(float)((r1*i2)-(r2*i1))/(m1*m2);
	}
	else{
	  ref[p]=1.0;
	  ref[p+1]=0.0;
	}
      }
    }
  }		
  return(0);
}
void hu_cor(float *data,float *ref,int ppl,int lpi,int ipv)
{
  int x,y;
  long p;
  float r,i;

  /* Apply correction */
  for(y=1;y<(lpi*ipv);y+=2){
    for(x=0;x<ppl*2;x+=2){
      p=(y*ppl*2)+x;
      /* cos(ph2-ph1) stored in ref[p]   */
      /* sin(ph2-ph1) stored in ref[p+1] */
      r=data[p];i=data[p+1];
      data[p]=(r*ref[p])-(i*ref[p+1]);
      data[p+1]=(r*ref[p+1])+(i*ref[p]);
    }
  }
}
void phase_rot_init(float *data,float *ref,long volsize)
{
  long l;
  for(l=0;l<volsize*2;l++){
    ref[l]=data[l];
  }
}
void phase_rot(float *data,float *ref,long volsize)
{
  long j;
  float r,i,r1,i1,mod;

  for(j=0;j<volsize;j++){
    r=data[j*2];i=data[j*2+1];
    r1=ref[j*2];i1=ref[j*2+1];
    if((mod=magnitude(r1,i1))!=0.0){
      data[j*2] = (r*r1 + i*i1)/mod;
      data[j*2+1] = (r1*i - r*i1)/mod;
    }
  }
}
void smodulus_init(float *data,float *ref,long volsize)
{
  long i,j;
  
  for(i=0,j=0;i<volsize;i+=2,j++)
    ref[j]=atan2(data[i+1],data[i]);
}
void smodulus(float *data,float *ref,long volsize)
{
  long i,j;
  float mod,phs;

  for(i=0,j=0;i<volsize;i+=2,j++){
    mod=magnitude(data[i+1],data[i]);
    phs=atan2(data[i+1],data[i])-ref[i/2];
    if(phs<-M_PI)phs+=2*M_PI;
    if(phs>M_PI)phs-=2*M_PI;
    if((phs<-M_PI/2)||(phs>=M_PI/2))data[j]=-mod;
    else data[j]=mod;
  }
}
