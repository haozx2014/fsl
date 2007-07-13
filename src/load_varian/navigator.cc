/*
  navigator.c: Extracts and processes navigatior echoes

  Copyright Stuart Clare, FMRIB Centre, University of Oxford.

  This program should be considered a beta test version
  and must not be used for any clinical purposes.

  Part of ...
  LoadVarian: Turns time data from the Varian fids to images
  For full version history see main.c
*/

#include "navigator.h"
#include <cstdlib>
#include <cstdio>

void extract_navigators(float *data,float *nav,int ppl,int lpi,int ipv,int seg)
{
  int p,l,s,i;
  int lps;

  lps = (lpi-seg)/seg;

  for(i=0;i<ipv;i++){
    for(s=0;s<seg;s++){
      for(p=0;p<ppl*2;p++){
	nav[(i*seg + s)*ppl*2 + p]=data[((i*lpi) + (s*(lps+1)))*ppl*2 + p];
      }
      for(l=0;l<lps;l++){
	for(p=0;p<ppl*2;p++){
	  data[((i*lps*seg) + (s*lps) + l)*ppl*2 + p]
	    = data[((i*lpi) + (s*(lps+1)) + l + 1)*ppl*2 + p];
	}
      }
    }
  }
}
