/* Copyright Numerical Recipes in C, Cambridge University Press */
/* Not to be distributed with code */

#include "vector.h"
#include "nr.h"

#include <math.h>
#define NRANSI
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

void powell_lim(float p[], float **xi, int n, float ftol, int *iter,
		float *fret,float (*func)(float []),int itmax)
{
  int i,ibig,j;
  float del,fp,fptt,t,*pt,*ptt,*xit;
  float sqrarg;

  pt=vector(1,n);
  ptt=vector(1,n);
  xit=vector(1,n);
  *fret=(*func)(p);
  for (j=1;j<=n;j++) pt[j]=p[j];
  for (*iter=1;;++(*iter)) {
    fp=(*fret);
    ibig=0;
    del=0.0;
    for (i=1;i<=n;i++) {
      for (j=1;j<=n;j++) xit[j]=xi[j][i];
      fptt=(*fret);
      linmin(p,xit,n,fret,func);
      if (fabs(fptt-(*fret)) > del) {
	del=fabs(fptt-(*fret));
	ibig=i;
      }
    }
    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
      free_vector(xit,1,n);
      free_vector(ptt,1,n);
      free_vector(pt,1,n);
      return;
    }
    /* Bomb out if too many iterations */
    if (*iter == itmax) {
      for(j=1;j<=n;j++)p[j]=0;
      return;
    }
    for (j=1;j<=n;j++) {
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    fptt=(*func)(ptt);
    if (fptt < fp) {
      t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
      if (t < 0.0) {
	linmin(p,xit,n,fret,func);
	for (j=1;j<=n;j++) {
	  xi[j][ibig]=xi[j][n];
	  xi[j][n]=xit[j];
	}
      }
    }
  }
}
#undef NRANSI
#undef SQR
