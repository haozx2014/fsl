/* {{{ Copyright etc. */

/*  tsplot - FMRI time series and model plotting

    Stephen Smith, Mark Woolrich and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2007 University of Oxford  */

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
/* {{{ background theory */

/*

GLM : Y = Xb + e

The "partial model fit" shows, in the case of a contrast which selects
a single EV, the part of the full model fit explained by that EV. In
the case of a more complex contrast, it is basically a plot of the sum
of each EV weighted by the product of that EV's PE and that EV's
contrast weighting.

i.e. the partial model fit is X*diag(c)*b, where X is the design and c
is the contrast vector (renormalised to unit length and then turned
into a diagonal matrix) and b is the parameter estimate vector.

Thus we plot this versus the "reduced data" ie plot

Y - Xb + X*diag(c)*b  vs  X*diag(c)*b
i.e.
residuals + partial fit    vs   partial fit

NOTE: this plot cannot simply be used to generate the t/Z score
associated with this contrast (eg by straight correlation) - this
would not be correct. In order to do that you would need to correlate
the Hansen projection c'*pinv(X) with Y instead.

*/

/* }}} */
/* {{{ defines, includes and typedefs */
#include <iomanip>
#include "featlib.h"
#include "libvis/miscplot.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "utils/options.h"
#include <vector>
 
using namespace MISCPLOT;
using namespace MISCMATHS;
using namespace Utilities;
using namespace std;
 
/*
void setupq (ColumnVector &x,ColumnVector &dx,ColumnVector &y,int npoint,Matrix &v,Matrix &a )
{
double diff,prev;
   v(1,4) = x(2) - x(1);
   for(int i=2;i<npoint;i++)
   {
      v(i,4) = x(i+1) - x(i);
      v(i,1) = dx(i-1)/v(i-1,4);
      v(i,2) = - dx(i)/v(i,4) - dx(i)/v(i-1,4);
      v(i,3) = dx(i+1)/v(i,4);
   }
   v(npoint,1) = 0;
   for(int i=2;i<npoint;i++) v(i,5) = v(i,1)*v(i,1) + v(i,2)*v(i,2) + v(i,3)*v(i,3);

   if ( npoint >= 4 ) for(int i=3;i<npoint;i++) v(i-1,6) = v(i-1,2)*v(i,1) + v(i-1,3)*v(i,2);

    v(npoint-1,6) = 0;

   if (npoint >= 5) for(int i=4;i<npoint;i++) v(i-2,7) = v(i-2,3)*v(i,1);

   v(npoint-2,7) = 0;
   v(npoint-1,7) = 0;
//construct  q-transp. * y  in  qty.
   prev = (y(2) - y(1))/v(1,4);
   for(int i=2;i<npoint;i++)
   {  
      diff = (y(i+1)-y(i))/v(i,4);
      a(i,4) = diff - prev;
      prev = diff;
   }
}

//  from  * a practical guide to splines *  by c. de boor    
//  from  * a practical guide to splines *  by c. de boor    
// to be called in  s m o o t h constructs the upper three diags. in v(i,j), i=2,,npoint-1, j=1,3, of
//  the matrix  6*(1-p)*q-transp.*(d**2)*q + p*r, then computes its
//  l*l-transp. decomposition and stores it also in v, then applies
//  forward and backsubstitution to the right side q-transp.*y in  qty
//  to obtain the solution in  u .
//  a(1,4)=qty a(1,3)=u a(1,1)=qu
void chol1d (double p, Matrix &v, Matrix &a,int &npoint)
{
double prev,ratio,six1mp,twop;
   six1mp = 6*(1-p);
   twop = 2*p;
   for(int i=2;i<npoint;i++)
   {
      v(i,1) = six1mp*v(i,5) + twop*(v(i-1,4)+v(i,4));
      v(i,2) = six1mp*v(i,6) + p*v(i,4);
      v(i,3) = six1mp*v(i,7);
   }
   if (npoint < 4)
   {     
     a(1,3) = 0;
     a(2,3) = a(2,4)/v(2,1);
     a(3,3) = 0;
     goto fortyone;
   }
//  factorization
   for(int i=2;i<npoint-1;i++)
   {
     ratio = v(i,2)/v(i,1);
     v(i+1,1) = v(i+1,1) - ratio*v(i,2);
     v(i+1,2) = v(i+1,2) - ratio*v(i,3);
     v(i,2) = ratio;
     ratio = v(i,3)/v(i,1);
     v(i+2,1) = v(i+2,1) - ratio*v(i,3);
     v(i,3) = ratio;
   }
//  forward substitution
   a(1,3) = 0;
   v(1,3) = 0;
   a(2,3) = a(2,4);
   for(int i=2;i<npoint-1;i++) a(i+1,3) = a(i+1,4) - v(i,2)*a(i,3) - v(i-1,3)*a(i-1,3);
//  back substitution
   a(npoint,3) = 0;
   a(npoint-1,3) = a(npoint-1,3)/v(npoint-1,1);
   int i = npoint-2;
   do
   {
     a(i,3) = a(i,3)/v(i,1)-a(i+1,3)*v(i,2)-a(i+2,3)*v(i,3);
   } while (--i > 1);
//  construct q*u
   fortyone:
   prev = 0;
   for (int i=2;i<=npoint;i++)
   {
     a(i,1) = (a(i,3) - a(i-1,3))/v(i-1,4);
     a(i-1,1) = a(i,1) - prev;
     prev = a(i,1);
   }
   a(npoint,1) = -a(npoint,1);
}

void smooth(ColumnVector &x,ColumnVector &y,ColumnVector &dy,int npoint,double s)
{
  Matrix a(npoint,4);
  Matrix v(npoint,7);
  double change,ooss,oosf,p,prevsf,prevq,q=0,sfq,sixp,six1mp,utru;
  setupq(x,dy,y,npoint,v,a);

  if ( s > 0 )                    
  {
   p = 0;                     
   chol1d(p,v,a,npoint);
   sfq = 0;
   for (int i=1;i<=npoint;i++) sfq = sfq + pow(a(i,1)*dy(i),2.0);
   sfq*=36;
   if (sfq < s) goto sixty;
   utru = 0;
   for (int i=2;i<=npoint;i++) utru+= v(i-1,4)*(a(i-1,3)*(a(i-1,3)+a(i,3))+pow(a(i,3),2.0));
   ooss = 1./sqrt(s);
   oosf = 1./sqrt(sfq);
   q = -(oosf-ooss)*sfq/(6.*utru*oosf);
   prevq = 0;
   prevsf = oosf;

   thirty:
   chol1d(q/(1.+q),v,a,npoint);
   sfq = 0;
   for(int i=1;i<=npoint;i++) sfq = sfq + pow(a(i,1)*dy(i),2.0);
   sfq*=36.0/pow(1+q,2.0);
   if (abs(sfq-s) < 0.01*s) goto fiftynine;
   oosf = 1.0/sqrt(sfq);
   change = (q-prevq)/(oosf-prevsf)*(oosf-ooss);
   prevq = q;
   q-= change;
   prevsf = oosf;
   goto thirty;
  }
  else 
  {
    p = 1;                     
    chol1d(p,v,a,npoint);
    sfq = 0;
    goto sixty;             
   }


   fiftynine: 
   p = q/(1.0+q);
//correct value of p has been found.
//compute pol.coefficients from  Q*u (in a(.,1)).
   sixty: 
   six1mp = 6./(1.+q);
   for(int i=1;i<=npoint;i++) a(i,1) = y(i) - six1mp*pow(dy(i),2.0)*a(i,1);
   sixp = 6*p;
   for(int i=1;i<=npoint;i++) 
   {
     a(i,3)*=sixp;
     y(i)=a(i,1); 
   }
   for(int i=1;i<npoint;i++)  
   {
     a(i,4) = (a(i+1,3)-a(i,3))/v(i,4);
     a(i,2) = (a(i+1,1)-a(i,1))/v(i,4)- (a(i,3)+a(i,4)/3.*v(i,4))/2.*v(i,4);
   }
}
*/

/* }}} */
/* {{{ usage */

void usage(void)
{
  printf("Usage: tsplot <feat_directory.feat> [options]\n");
  printf("[-f <4D_data>] input main filtered data, in case it's not <feat_directory.feat>/filtered_func_data\n");
  printf("[-c <X Y Z>] : use X,Y,Z instead of max Z stat position\n");
  printf("[-C <X Y Z output_file.txt>] : use X,Y,Z to output time series only - no stats or modelling\n");
  printf("[-m <mask>] : use mask image instead of thresholded activation images\n");
  printf("[-o <output_directory>] change output directory from default of input feat directory\n");
  printf("[-n] don't weight cluster averaging with Z stats\n");
  printf("[-p] prewhiten data and model timeseries before plotting\n");
  printf("[-d] don't keep raw data text files\n");
  exit(1);
}

/* }}} */

int main(int argc, char **argv)
{
  /* {{{ variables */

FILE         *ifp, *rofp;
ofstream output_file;
int          argi=1, prewhiten=0, ev, i, j, t, v, X=0, Y=0, Z=0, nevs, npts,
  ncon=1, nftests=0, size, ymin, ymax, coordset=0, dataonly=0, modelfree=0, zweight=1,
  use_triggers, level, custommask=0;
double       *model, *contrasts=NULL, *norm_contrasts=NULL, tsmean, tmpf, maxz;
float        *triggers;
char         rofpM[100000], rofpF[100000], rofpP[100000], 
  fmridata[10000], fsldir[10000], featdir[10000], outputdir[10000], gpname[10000], gprootname[10000],
  thestring[10000], datafile[10000], statname[100], vname[100];
 ColumnVector pwts;
 bool textfiles=true;
volume<float> immask;

 int GRPHSIZE=600;
  int PSSIZE=600;

rofpM[0]='\0';

/* }}} */

  /* {{{ process arguments */

if (argc<2) usage();

strcpy(featdir,argv[argi++]);

strcpy(outputdir,featdir);

sprintf(fmridata,"%s/filtered_func_data",featdir);
sprintf(fsldir,getenv("FSLDIR"));

for (;argi<argc;argi++)
{
  if (!strcmp(argv[argi], "-f"))
    /* {{{ alternative fmri data */

{
  argi++;
  if (argc<argi+1)
    {
      printf("Error: no value given following -f\n");
      usage();
    }
  strcpy(fmridata,argv[argi]);
}

/* }}} */
  else if (!strcmp(argv[argi], "-c"))
    /* {{{ alternative voxel position */

{
  coordset=1;
  argi++;
  if (argc<argi+3) /* options following c haven't been given */
    {
      printf("Error: incomplete values given following -c\n");
      usage();
    }
  X=atoi(argv[argi++]);
  Y=atoi(argv[argi++]);
  Z=atoi(argv[argi]);
}

/* }}} */
  else if (!strcmp(argv[argi], "-C"))
    /* {{{ output data only */

{
  coordset=1;
  dataonly=1;
  argi++;
  if (argc<argi+4)
    {
      printf("Error: incomplete values given following -C\n");
      usage();
    }
  X=atoi(argv[argi++]);
  Y=atoi(argv[argi++]);
  Z=atoi(argv[argi++]);
  strcpy(datafile,argv[argi]);
}

/* }}} */
  else if (!strcmp(argv[argi], "-m"))
    /* {{{ alternative mask image */

{
  custommask=1;

  argi++;
  if (argc<argi+1)
    {
      printf("Error: no mask image given following -m\n");
      usage();
    }

  if ( read_volume(immask,argv[argi]) )
    {
      printf("Error: mask image chosen doesn't exist\n");
      usage();
    }
 }

/* }}} */
  else if (!strcmp(argv[argi], "-o"))
    /* {{{ output dir */

{
  argi++;
  if (argc<argi+1)
    {
      printf("Error: no value given following -o\n");
      usage();
    }
  strcpy(outputdir,argv[argi]);
}

  else if (!strcmp(argv[argi], "-d"))
    /* don't output files */
{
  textfiles=false;
}

/* }}} */
  else if (!strcmp(argv[argi], "-n"))
    /* {{{ zweight clusters? */

{
  zweight=0;
}

/* }}} */
  else if (!strcmp(argv[argi], "-p"))
    /* {{{ turn on prewhitening */

{
  prewhiten=1;
}

/* }}} */
}

/* }}} */
  /* {{{ read filtered_func_data */

volume4D<int> im;
read_volume4D(im, fmridata);

size=im.nvoxels();

if (dataonly && textfiles)
  /* {{{ output raw data and exit */
{
  output_file.open(datafile);
  if(!output_file.is_open())
    {
      fprintf(stderr,"Can't open output data file %s\n",datafile);
      exit(1);
    }

  for(t=0; t<im.tsize(); t++) output_file << scientific << im(X,Y,Z,t) << endl;
  output_file.close();
  return 0;
}

/* }}} */

/* }}} */
  /* {{{ read design.mat */

sprintf(thestring,"%s/design.mat",featdir);
model=read_model(thestring,&nevs,&npts);

if (npts==0)
{
  modelfree=1;
  nftests=1;
  ncon=0;
}

npts=im.tsize();

ColumnVector TS_model(npts),TS_copemodel(npts),TS_pemodel(npts*nevs),TS_data(npts),TS_residuals(npts);


/* }}} */
  /* {{{ read auto correlation estimates for prewhitening */

double* pwmodel=0;
volume4D<float> acs;

if ( prewhiten ) {

  prewhiten=0;

  sprintf(thestring,"%s/stats/threshac1",featdir);
  if (fsl_imageexists(string(thestring)))
    {
      read_volume4D(acs, thestring);

      if (acs[1].max()!=0) /* hacky test for whether prewhitening was actually carried out */
	{
	  pwmodel=new double[nevs*npts];
	  prewhiten=1;
	}
    }

}

/* }}} */
  /* {{{ read design.con and PEs */

vector< volume<float> > impe(nevs);
if (!modelfree)
{
  sprintf(thestring,"%s/design.con",featdir);
  contrasts=read_contrasts(thestring,&nevs,&ncon);

  /* create normalised contrasts */
  norm_contrasts=(double*)malloc(sizeof(double)* nevs * ncon);
  for(i=0; i<ncon; i++)
    {
      double norm_factor=0;

      for(ev=0; ev<nevs; ev++)
	norm_factor += contrasts[i*nevs+ev] * contrasts[i*nevs+ev];

      for(ev=0; ev<nevs; ev++)
	norm_contrasts[i*nevs+ev] = contrasts[i*nevs+ev] / sqrt(norm_factor);
    }

  for(i=1;i<=nevs;i++)
    {
      sprintf(thestring,"%s/stats/pe%d",featdir,i);
      read_volume(impe[i-1], thestring);
    }
}

/* }}} */
  /* {{{ read design.fts */

if (!modelfree)
{
  sprintf(thestring,"%s/design.fts",featdir);
  read_ftests(thestring,&nftests);
}

/* }}} */
  /* {{{ read triggers */

use_triggers=read_triggers(featdir,&triggers,nevs,npts);


/*for(ev=0;ev<nevs;ev++)
{
  for(int i=0;i<=triggers[ev]+1;i++) printf("%f ",triggers[i*nevs+ev]);
  printf("\n");
}
*/

/* }}} */
  /* {{{ check analysis level */

level=1;

sprintf(thestring,"%s/design.lev",featdir);
if((ifp=fopen(thestring,"rb"))!=NULL)
{
  fclose(ifp);
  level=2;
}

/* }}} */
  /* {{{ create plot(s) for each contrast */

for(j=0;j<2;j++)
{
  /* {{{ setup stats type */

int maxi;

if (j==0) { sprintf(statname,"zstat");  maxi=ncon; }
else      { sprintf(statname,"zfstat"); maxi=nftests; }

/* }}} */

  for(i=0; i<maxi; i++)
    {
      /* {{{ vars */

volume<float> imcope, imz, imweight;
bool haveclusters=false;

rofpF[0]='\0';
rofpP[0]='\0';

/* }}} */
      /* {{{ read COPE and derived stats; test for f-test output */

/* load zstat or zfstat */
sprintf(thestring,"%s/stats/%s%d",featdir,statname,i+1);
if (fsl_imageexists(string(thestring))) 
{
  read_volume(imz,thestring);
  imweight=imz;
  if (!zweight)
    imweight=1;
}
else
  continue; /* f-test i wasn't valid - no zfstat image */

/* load cope */
if ( (j==0) && (!modelfree) )
{
  sprintf(thestring,"%s/stats/cope%d",featdir,i+1);
  read_volume(imcope,thestring);
}

/* load cluster mask */
if (!coordset)
{
  if (!custommask)
    {
      sprintf(thestring,"%s/cluster_mask_%s%d",featdir,statname,i+1);
      if (fsl_imageexists(string(thestring)))
	read_volume(immask,thestring);
    }
  haveclusters=(immask.max()>0);
}

/* }}} */
      /* {{{ find max Z and X,Y,Z */

if (!coordset)
{
  X=Y=Z=0;
  maxz=-1000;

  for(int z=0; z<im.zsize(); z++)
    for(int y=0; y<im.ysize(); y++)
      for(int x=0; x<im.xsize(); x++)
	if ( (imz(x,y,z)>maxz) &&
	     ( (!haveclusters) ||   /* make max Z be inside a cluster if we found a cluster map */
	       (immask(x,y,z)>0) && (!prewhiten || acs(x,y,z,1)!=0 || acs(x,y,z,2)!=0) ) )
	  {
	    maxz=imz(x,y,z);
	    X=x; Y=y; Z=z;
	  }
}
else
  maxz=imz(X,Y,Z);

/* }}} */

      /* first do peak voxel plotting then do mask-averaged plotting */
      for(v=0;v<2;v++)
	{
          /* {{{ setup */

double wtotal=0;
int count=0;
volume<float> tmp_imweight=imweight, tmp_immask=immask;

if (v==0) {
  vname[0]='\0';
  tmp_imweight=0;
  tmp_imweight(X,Y,Z)=1;
  tmp_immask=tmp_imweight;
} else {
  sprintf(vname,"c");
}

/* }}} */
          /* {{{ create model and data time series */

TS_model=0;
TS_residuals=0;
TS_copemodel=0;
TS_data=0;
TS_pemodel=0;

for(int x=0; x<im.xsize(); x++) for(int y=0; y<im.ysize(); y++) for(int z=0; z<im.zsize(); z++)
  if (tmp_immask(x,y,z)>0 && (!prewhiten || acs(x,y,z,1)!=0 || acs(x,y,z,2)!=0))
  {
    count++;
    wtotal+=tmp_imweight(x,y,z);

    if(prewhiten)
      prewhiten_timeseries(acs.voxelts(x,y,z), im.voxelts(x,y,z), pwts, npts);
    else
      pwts = im.voxelts(x,y,z);
    for(t=1; t<=npts; t++) TS_data(t)+=pwts(t)*tmp_imweight(x,y,z);
    

    if (!modelfree) {
      if (prewhiten)
	prewhiten_model(acs.voxelts(x,y,z), model, pwmodel, nevs, npts);
      else
	pwmodel=model;
      for(t=1; t<=npts; t++)
	for(ev=0; ev<nevs; ev++)
	  {
	    tmpf=pwmodel[(t-1)*nevs+ev]*impe[ev](x,y,z)*tmp_imweight(x,y,z);
            TS_model(t)           += tmpf;
            TS_copemodel(t)       += tmpf*norm_contrasts[i*nevs+ev];
            TS_pemodel(ev*npts+t) += tmpf;
	  }
    }
  }

tsmean=0;
for(t=1; t<=npts; t++)
{
  TS_data(t)/=wtotal;
  tsmean+=TS_data(t);
}
tsmean/=npts;

if (level==2) tsmean=0;

if (!modelfree)
  for(t=1; t<=npts; t++)
    {
      TS_model(t) = TS_model(t)/wtotal + tsmean;
      TS_copemodel(t) = TS_copemodel(t)/wtotal + tsmean;
      TS_residuals(t)=TS_data(t)-TS_model(t);
      for(ev=0; ev<nevs; ev++)
	TS_pemodel(ev*npts+t) = TS_pemodel(ev*npts+t)/wtotal + tsmean;
    }


/* }}} */
	  /* {{{ output data text files */

sprintf(thestring,"%s/tsplot%s_%s%d.txt",outputdir,vname,statname,i+1);
if (textfiles) output_file.open(thestring);
   ymin=ymax=(int)TS_data(1);
   for(t=1; t<=npts; t++)
   {
     if (textfiles) output_file << scientific << TS_data(t);
     ymin=(int)MISCMATHS::Min(TS_data(t),ymin); 
     ymax=(int)MISCMATHS::Max(TS_data(t),ymax);
     if (!modelfree)
     {
       if (j==0)
       {
          if (textfiles) output_file << " " << TS_copemodel(t); 
          ymin=(int)MISCMATHS::Min(TS_copemodel(t),ymin); 
          ymax=(int)MISCMATHS::Max(TS_copemodel(t),ymax);
       }
       if (textfiles) output_file << " " << TS_model(t); 
       ymin=(int)MISCMATHS::Min(TS_model(t),ymin); 
       ymax=(int)MISCMATHS::Max(TS_model(t),ymax);
       if (j==0) output_file << " " << TS_residuals(t)+TS_copemodel(t);
     }
     if (textfiles) output_file << endl;
   }
if (textfiles) output_file.close();
ymax+=(ymax-ymin)/5;
ymin-=(ymax-ymin)/20;

/* }}} */
	  /* {{{ create graphs */

  sprintf(gprootname,"tsplot%s_%s%d",vname,statname,i+1); 
  sprintf(gpname,"%s/%s",outputdir,gprootname);
  miscplot newplot;
  GRPHSIZE= MISCMATHS::Min(MISCMATHS::Max(npts*4,600),3000);
  newplot.set_minmaxscale(1.001);
  newplot.set_xysize(GRPHSIZE,192);
  newplot.set_yrange(ymin,ymax);
  string title=statname+num2str(i+1);
  if (v==0) 
  {
    if (!coordset) title+= ": max Z stat of "+num2str(maxz)+" at voxel ("+num2str(X)+" "+num2str(Y)+" "+num2str(Z)+")";
    else title+= ": Z stat of "+num2str(maxz)+" at selected voxel ("+num2str(X)+" "+num2str(Y)+" "+num2str(Z)+")";
  } 
  else title+= ": averaged over "+num2str(count)+" voxels";
  Matrix blank=TS_data;
  blank=log(-1.0);
  if (!modelfree)
  {
    if (j==0)
    {
      newplot.add_label("full model fit");
      newplot.add_label("cope partial model fit");
      newplot.add_label("data");
      newplot.timeseries((TS_model | TS_copemodel | TS_data).t(),string(gpname),title,1,GRPHSIZE,4,2,false);
      newplot.remove_labels(3);
      newplot.add_label("");
      newplot.add_label("cope partial model fit");
      newplot.add_label("reduced data");
      newplot.timeseries((blank | TS_copemodel | TS_residuals+TS_copemodel ).t(),string(gpname)+"p",title,1,GRPHSIZE,4,2,false);
      newplot.remove_labels(3);
      sprintf(rofpF,"%sFull model fit - <a href=\"%sp.png\">Partial model fit</a> - <a href=\"%s.txt\">Raw data</a><br>\n<IMG BORDER=0 SRC=\"%s.png\"><br><br>\n",rofpF,gprootname,gprootname,gprootname);
    }
    else
    {
      newplot.add_label("full model fit");
      newplot.add_label("");
      newplot.add_label("data");
      newplot.timeseries((TS_data | blank | TS_model).t(),string(gpname),title,1,GRPHSIZE,4,2,false);
      newplot.remove_labels(3);
      sprintf(rofpF,"%sFull model fit - <a href=\"%s.txt\">Raw data</a><br>\n<IMG BORDER=0 SRC=\"%s.png\"><br><br>\n",rofpF,gprootname,gprootname);
    }
  }
  else
  {    
      newplot.add_label("");
      newplot.add_label("");
      newplot.add_label("data");
      newplot.timeseries((blank | blank | TS_data).t(),string(gpname),title,1,GRPHSIZE,4,2,false);
      newplot.remove_labels(3);
      sprintf(rofpF,"%sData plot - <a href=\"%s.txt\">Raw data</a>\n<IMG BORDER=0 SRC=\"%s.png\"><br><br>\n",rofpF,gprootname,gprootname);
  }
/* picture for main web index page */
  if (v==0) sprintf(rofpM,"%s<a href=\"%s.html\"><IMG BORDER=0 SRC=\"%s.png\"></a><br><br>\n",rofpM,gprootname,gprootname);

  	  /* {{{ peri-stimulus: output text and graphs */
  if (use_triggers)
  {
    if (!modelfree) sprintf(rofpP,"%s<table><tr>\n",rofpP);
    for(ev=0; ev<nevs; ev++) if (triggers[ev]>0.5)
    {
      float ps_period=triggers[((int)triggers[ev]+1)*nevs+ev];
      Matrix ps_compact((int)(10*ps_period)+1,3);
      if (!modelfree) ps_compact.ReSize((int)(10*ps_period)+1,6);

      int size=0;
      for(int which_event=1;which_event<=triggers[ev];which_event++)
	for(t=(int)ceil(triggers[which_event*nevs+ev])+1;t<=MISCMATHS::Min(npts-1,(int)ceil(triggers[which_event*nevs+ev])+(int)ps_period);t++) size++;

      Matrix ps_full(size,ps_compact.Ncols()-1); 
      ps_compact=0;
      ps_full=0;
      int ptr=0;
      for(int which_event=1;which_event<=triggers[ev];which_event++)
      {
        double min_t=triggers[which_event*nevs+ev];
        int int_min_t=(int)ceil(min_t),max_t=MISCMATHS::Min(npts-1,int_min_t+(int)ps_period);
        for(t=int_min_t+1;t<=max_t;t++)
	{
          Matrix input(1,ps_compact.Ncols());
          if (!modelfree) input.Row(1) << ((int)((t-min_t-1)*10))/10.0 << TS_residuals(t)+TS_model(t) << TS_model(t) << TS_pemodel(ev*npts+t) << TS_residuals(t)+TS_pemodel(ev*npts+t) << 1;  //(restricted temporal accuraccy (0.1*TR)
          else input.Row(1) << t-min_t-1 << TS_residuals(t)+TS_model(t) << 1;
          ps_compact.Row(((int)((t-min_t-1)*10))+1)+=input.Row(1); 
          ps_full.Row(++ptr) = input.SubMatrix(1,1,1,input.Ncols()-1); 
	}
      }


      sprintf(gprootname,"ps_tsplot%s_%s%d_ev%d",vname,statname,i+1,ev+1);
      sprintf(thestring,"%s/%s.txt",outputdir,gprootname);
      sprintf(gpname,"%s/%s",outputdir,gprootname);  
      ps_full=ps_full.SubMatrix(1,ptr,1,ps_full.Ncols());
      
      if (textfiles) 
      {
	output_file.open(thestring,ofstream::out);
	for(int k=1;k<=ps_full.Nrows();k++)
	{ 
	  output_file << setprecision(1) << fixed << ps_full(k,1) << setprecision(6) << scientific;
	  for (int j=2;j<=ps_full.Ncols();j++) output_file << " " << ps_full(k,j);
	  output_file << endl;
	}
	output_file.close();
      }
      title=statname+num2str(i+1)+" ev"+num2str(ev+1);

      for(int j=1;j<=ps_compact.Nrows();j++) 
	{
	  if (ps_compact(j,6)) ps_compact.Row(j)/=ps_compact(j,6);
	  else  ps_compact.Row(j)=log10(-1.0); //deliberately set to nan
	}

      Matrix ps_interp=ps_compact.t();

      PSSIZE = MISCMATHS::Min(MISCMATHS::Max(ps_period*3,400),3000);
      newplot.set_minmaxscale(1.001);
      newplot.add_xlabel("peristimulus time (TRs)");
      newplot.set_xysize(PSSIZE,192);
      newplot.set_yrange(ymin,ymax);
      if (v==0) 
      {
        if (!coordset) title+= ": max Z stat of "+num2str(maxz)+" at voxel ("+num2str(X)+" "+num2str(Y)+" "+num2str(Z)+")";
        else title+= ": Z stat of "+num2str(maxz)+" at selected voxel ("+num2str(X)+" "+num2str(Y)+" "+num2str(Z)+")";
      } 
      else title+= ": averaged over "+num2str(count)+" voxels";
      if (!modelfree)
      {
        ps_compact=ps_full.SubMatrix(1,ps_full.Nrows(),1,2);
        ps_compact.Column(1)*=10;

        newplot.setscatter(ps_compact,(int)(10*(ps_period+3)));  //12+2
        newplot.add_label("full model fit");
	newplot.add_label("EV "+num2str(ev+1)+" model fit");
        newplot.add_label("data");
        newplot.timeseries(ps_interp.SubMatrix(3,4,1,ps_interp.Ncols()),string(gpname),title,-0.1,PSSIZE,3,2,false);
        newplot.remove_labels(3);
        ps_compact=ps_full.SubMatrix(1,ps_full.Nrows(),1,1) | ps_full.SubMatrix(1,ps_full.Nrows(),5,5);
        ps_compact.Column(1)*=10;
        newplot.setscatter(ps_compact,(int)(10*(ps_period+3)));
        newplot.add_label("");
        newplot.add_label("EV "+num2str(ev+1)+" model fit");
        newplot.add_label("reduced data");
        ps_interp.Row(3)=log(-1.0);
        newplot.timeseries(ps_interp.SubMatrix(3,4,1,ps_interp.Ncols()),string(gpname)+"p",title,-0.1,PSSIZE,3,2,false);
	newplot.deletescatter();
        newplot.remove_labels(3);
	sprintf(rofpP,"%s<td>Full model fit - <a href=\"%sp.png\">Partial model fit</a> - <a href=\"%s.txt\">Raw data</a><br>\n<IMG BORDER=0 SRC=\"%s.png\">\n",rofpP,gprootname,gprootname,gprootname);
      }
      else
      {
	Matrix blank=ps_full.SubMatrix(2,2,1,ps_compact.Ncols());
	blank=log(-1.0);
        newplot.add_label("");
        newplot.add_label("");
        newplot.add_label("data");
        newplot.timeseries(blank & blank & ps_full.SubMatrix(2,2,1,ps_compact.Ncols()),string(gpname),title,-0.1,PSSIZE,3,2,false);
        newplot.remove_labels(3);
        sprintf(rofpP,"%sData plot - <a href=\"%s.txt\">Raw data</a>\n<IMG BORDER=0 SRC=\"%s.png\"><br><br>\n",rofpP,gprootname,gprootname);
      }
      newplot.remove_xlabel();
    }
    if (!modelfree) sprintf(rofpP,"%s</tr></table><br><br>\n",rofpP);
  }
  if (!haveclusters) v=10;
}
      /* {{{ web output */
sprintf(thestring,"%s/tsplot_%s%d.html",outputdir,statname,i+1);

if((rofp=fopen(thestring,"wb"))==NULL)
{
  fprintf(stderr,"Can't open output report file %s\n",outputdir);
  exit(1);
}

  if (use_triggers) fprintf(rofp,"<HTML>\n<TITLE>%s%d</TITLE>\n<BODY BACKGROUND=\"file:%s/doc/images/fsl-bg.jpg\">\n<hr><CENTER>\n<H1>FEAT Time Series Report - %s%d</H1>\n</CENTER>\n<hr><b>Full plots</b><p>\n%s\n<hr><b>Peristimulus plots</b><p>\n%s\n<HR></BODY></HTML>\n\n",statname,i+1,fsldir,statname,i+1,rofpF,rofpP);
  else fprintf(rofp,"<HTML>\n<TITLE>%s%d</TITLE>\n<BODY BACKGROUND=\"file:%s/doc/images/fsl-bg.jpg\">\n<hr><CENTER>\n<H1>FEAT Time Series Report - %s%d</H1>\n</CENTER>\n<hr><b>Full plots</b><p>\n%s\n</BODY></HTML>\n\n",statname,i+1,fsldir,statname,i+1,rofpF);

fclose(rofp);

/* }}} */
    }
 }

/* {{{ main web index page output */

/* first output full index page (eg for use by featquery) */

sprintf(thestring,"%s/tsplot_index.html",outputdir);

if((rofp=fopen(thestring,"wb"))==NULL)
{
  fprintf(stderr,"Can't open output report file %s\n",outputdir);
  exit(1);
}

fprintf(rofp,"<HTML>\n<TITLE>FEAT Time Series Report</TITLE>\n<BODY BACKGROUND=\"file:%s/doc/images/fsl-bg.jpg\">\n<hr><CENTER>\n<H1>FEAT Time Series Report</H1>\n</CENTER>\n<hr>%s<HR></BODY></HTML>\n\n",fsldir,rofpM);

fclose(rofp);


/* now output same thing without start and end, for inclusion in feat report */

sprintf(thestring,"%s/tsplot_index",outputdir);

if((rofp=fopen(thestring,"wb"))==NULL)
{
  fprintf(stderr,"Can't open output report file %s\n",outputdir);
  exit(1);
}

fprintf(rofp,"%s\n\n",rofpM);

fclose(rofp);

/* }}} */

/* }}} */

  exit(0);
 }

