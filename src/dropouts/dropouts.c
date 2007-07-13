/* {{{ Copyright etc. */

/*  dropouts - fix slice dropouts

    Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 3.3 (c) 2006, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/1112. */

/* }}} */

#include "libss/libss.h"
#include "libss/libavw.h"

#define THRESHOLD 6 /* was 3 when I was testing for all types of dropout */

void usage()
{
  printf("Usage: dropouts <infile> <num_vols_to_ignore> [threshold]\n");
  printf("<infile_fixed> gets written if any slices are fixed\n");
  printf("[threshold] (default 6): each slice's error is compared against <threshold>*<median error for all slices>\n");
  exit(1);
}

int main(argc, argv)
  int   argc;
  char  *argv [];
{
/* {{{ vars */

image_struct im;
FILE   *fp;
FDT    *in, *out, *values, med;
int    x_size, y_size, z_size, t_size, x, y, z, t, tt, i, 
  junk, fixed_slices=0, orientation=1,
  X_SIZE=0, Y_SIZE=0, Z_SIZE=0;
double *fvalues, *errcounts, errcount, *medianerror, threshold=THRESHOLD;
char   thestring[1000];

/* }}} */

  /* {{{ process arguments and read image */

if (argc<3)
     usage();

if (argc>3)
     threshold=atof(argv[3]);

avw_read(argv[1],&im);

sprintf(thestring,"%s.dropouts",argv[1]);
fp=fopen(thestring,"wb");
fprintf(fp,"\nFix-Slice-Dropouts by Steve Smith, FMRIB\n");

junk = atoi(argv[2]);

x_size = im.x;
y_size = im.y;
z_size = im.z;
t_size = im.t - junk;

if ( (im.yv>im.xv) && (im.yv>im.zv) ) orientation=2; /* probably originally coronals */
if ( (im.xv>im.zv) && (im.xv>im.yv) ) orientation=3; /* probably originally sagittals */
fprintf(fp,"\nOrientation=%d\n\n",orientation);

switch (orientation) {
 case 1: X_SIZE=x_size, Y_SIZE=y_size, Z_SIZE = z_size; break;
 case 2: X_SIZE=z_size, Y_SIZE=x_size, Z_SIZE = y_size; break;
 case 3: X_SIZE=y_size, Y_SIZE=z_size, Z_SIZE = x_size; break;
}

in   = (FDT *) (im.i + x_size*y_size*z_size*junk);

out  = (FDT *) malloc(sizeof(FDT)*x_size*y_size*z_size*t_size);
memset (out,0,sizeof(FDT)*x_size*y_size*z_size*t_size);

values =      (FDT *)    malloc(sizeof(FDT)   *Z_SIZE*t_size*2);
fvalues =     (double *) malloc(sizeof(double)*Z_SIZE*t_size*2);
medianerror = (double *) malloc(sizeof(double)*Z_SIZE*       2);
errcounts =   (double *) malloc(sizeof(double)*Z_SIZE*t_size*2);

/* }}} */
  /* {{{ find error values for each slice */

fprintf(fp,"Error counts for each slice:\n");

for(t=0; t<t_size; t++)
{
  fprintf(fp,"[volume %3d] ",t+1+junk);
  for(z=0; z<Z_SIZE; z++)
    {
      errcount=0;
      for(y=0; y<Y_SIZE; y++)
	for(x=0; x<X_SIZE; x++)
	  {
	    FDT tmpi=0;
	    i=0;
	    for(tt=MAX(t-5,0); tt<=MIN(t+5,t_size-1); tt++)
	      {
		switch (orientation) {
		case 1: values[i]=IT(in,x,y,z,tt); break;
		case 2: values[i]=IT(in,y,z,x,tt); break;
		case 3: values[i]=IT(in,z,x,y,tt); break;
		}
		if (tt==t) tmpi=values[i];
		i++;
	      }
	    med = median(0.5,values,i);
	    errcount += (double)ABS(med-tmpi);
	    switch (orientation) {
	    case 1: IT(out,x,y,z,t)=med; break;
	    case 2: IT(out,y,z,x,t)=med; break;
	    case 3: IT(out,z,x,y,t)=med; break;
	    }
	  }

      errcounts[t*Z_SIZE+z] = errcount;
      fprintf(fp,"%5d ",(int)errcounts[t*Z_SIZE+z]);
    }
  fprintf(fp,"\n");
}

/* }}} */
  /* {{{ find median errors (across t) for each slice */

fprintf(fp,"[thresholds] ");
for(z=0; z<Z_SIZE; z++)
{
  double tmpf;

  for(t=0; t<t_size; t++)
    fvalues[t] = errcounts[t*Z_SIZE+z];

  tmpf=dmedian(0.5,fvalues,t_size)*threshold;
  medianerror[z] = MAX(100,tmpf);
  fprintf(fp,"%5d ",(int)medianerror[z]);
}
fprintf(fp,"\n");

/* }}} */
  /* {{{ for each t { divide slice error by cross-t-median-error and then sum overslices } */

fprintf(fp,"\nTotal normalised error for each volume:\n");

for(t=0; t<t_size; t++)
{
  double tmpf=0;

  for(z=0; z<Z_SIZE; z++)
    tmpf += (float)errcounts[t*Z_SIZE+z] / medianerror[z];

  fprintf(fp,"[Volume %d] %.3f\n",t,tmpf);
}
fprintf(fp,"\n");

/* }}} */
  /* {{{ threshold slice errors and apply to data */

fprintf(fp,"\nCorrected slices:\n");

for(t=0; t<t_size; t++)
{
  char thestring2[1000];
  int printthisline=0;

  sprintf(thestring,"[volume %3d] ",t+1+junk);
  for(z=0; z<Z_SIZE; z++)
    {
      if ( errcounts[t*Z_SIZE+z] < medianerror[z] )
	errcounts[t*Z_SIZE+z]=0.1;
      else
	{
	  errcounts[t*Z_SIZE+z]=1.1;
	  for(y=0; y<Y_SIZE; y++)
	    for(x=0; x<X_SIZE; x++)
	      switch (orientation) {
	      case 1: IT(in,x,y,z,t)=IT(out,x,y,z,t); break;
	      case 2: IT(in,y,z,x,t)=IT(out,y,z,x,t); break;
	      case 3: IT(in,z,x,y,t)=IT(out,z,x,y,t); break;
	      }
	  printthisline=1;
	}
      sprintf(thestring2,"%s %d",thestring,(int)errcounts[t*Z_SIZE+z]);
      strcpy(thestring,thestring2);
      if (errcounts[t*Z_SIZE+z]>0.5)
	fixed_slices++;
    }
  if ( printthisline)
    fprintf(fp,"%s\n",thestring);
}

fprintf(fp,"\nFIXED %d\n",fixed_slices);

/* }}} */
  /* {{{ output and exit */

if (fixed_slices>0)
{
  sprintf(thestring,"%s_fixed",argv[1]);
  avw_write(thestring,im);
}

return(0);

/* }}} */
}
