#include "dropouts.h"

#include <cstdlib>
#include <cmath>

#include "procfunc.h"

void calc_dropouts(float *drop,float *errcount,int ppl,int lpi,int ipv,int nvol,int cvol)
{
  int i,j,k,t;
  float *ts,med,cval;

  ts=(float *)malloc(nvol*sizeof(float));

  for(i=0;i<ipv;i++){
    errcount[i]=0;
    for(j=0;j<lpi;j++){
      for(k=0;k<ppl;k++){
	for(t=0;t<nvol;t++){
	  ts[t] = dropout_xtrt(drop,k,j,i,t,ppl,lpi,ipv,nvol);
	}
	cval=ts[cvol];
	med = median(0.5,ts,nvol);
	errcount[i] += fabs(med-cval);
      }
    }
  }
  free(ts);
}
float dropout_xtrt(float *data,int x,int y,int z,int t,int xsize,int ysize,int zsize,int tsize)
{
  return(magnitude(data[(x+(y*xsize)+(z*xsize*ysize)+(t*xsize*ysize*zsize))*2],
		   data[(x+(y*xsize)+(z*xsize*ysize)+(t*xsize*ysize*zsize))*2+1]));
}
int count_dropouts(float *dropouts,int ipv,int nt)
{
  int j,k,fix=0;
  float *fvalues,*thresholds;
  fvalues = (float *)malloc(nt*sizeof(float));
  thresholds = (float *)malloc(ipv*sizeof(float));

  for(k=0; k<ipv; k++) {
    for(j=0; j<nt; j++)
      fvalues[j] = dropouts[j*ipv+k];
    thresholds[k] = (median(0.5,fvalues,nt)*8);
  }
  
  for(j=0;j<nt; j++){
    for(k=0;k<ipv;k++){
      if(dropouts[j*ipv+k]>thresholds[k])fix++;
    }
  }
  free(thresholds);
  free(fvalues);
  return(fix);
}
