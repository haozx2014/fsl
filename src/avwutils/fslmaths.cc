//     
//     fslmaths.cc Image processing routines, some basic, some not so basic...
//     Steve Smith, David Flitney, Stuart Clare and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2000-2007 University of Oxford  
//     
/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 4.0 (c) 2007, The University of
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
//     

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"


#if defined ( __CYGWIN__ ) ||  defined (__sun)
extern "C" { 
#include <ieeefp.h> 
}
#endif

using namespace MISCMATHS;
using namespace NEWIMAGE;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

void print_usage(const string& progname) 
{
  cout << "Usage: fslmaths [-dt <datatype>] <first_input> [operations and inputs] <output> [-odt <datatype>]" << endl;

  cout << "\nDatatype information:" << endl;
  cout << " -dt sets the native (math operations) datatype (default float for all except double images)" << endl;
  cout << " -odt sets the output datatype (default as original image)" << endl;
  cout << " Possible datatypes are: char short int float double" << endl;
  cout << " Additionally \"-dt input\" will set the native datatype to that of the original image" << endl;

  cout << "\nBinary operations:" << endl;
  cout << "  (some inputs can be either an image or a number)" << endl;
  cout << " -add   : add following input to current image" << endl;
  cout << " -sub   : subtract following input from current image" << endl;
  cout << " -mul   : multiply current image by following input" << endl;
  cout << " -div   : divide current image by following input" << endl;
  cout << " -rem   : modulus remainder - divide current image by following input and take remainder" << endl;
  cout << " -mas   : use (following image>0) to mask current image" << endl;
  cout << " -thr   : use following number to threshold current image (zero anything below the number)" << endl;
  cout << " -thrp  : use following percentage (0-100) of ROBUST RANGE to threshold current image (zero anything below the number)" << endl;
  cout << " -thrP  : use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold below" << endl;
  cout << " -uthr  : use following number to upper-threshold current image (zero anything above the number)" << endl;
  cout << " -uthrp : use following percentage (0-100) of ROBUST RANGE to upper-threshold current image (zero anything above the number)" << endl;
  cout << " -uthrP : use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold above" << endl;
  cout << " -max   : take maximum of following input and current image" << endl;
  cout << " -min   : take minimum of following input and current image" << endl;

  cout << "\nBasic unary operations:" << endl;
  cout << " -exp   : exponential" << endl;
  cout << " -log   : natural logarithm" << endl;
  cout << " -sqr   : square" << endl;
  cout << " -sqrt  : square root" << endl;
  cout << " -abs   : absolute value" << endl;
  cout << " -bin   : use (current image>0) to binarise" << endl;
  cout << " -index : replace each nonzero voxel with a unique (subject to wrapping) index number" << endl;
  cout << " -grid <value> <spacing> : add a 3D grid of intensity <value> with grid spacing <spacing>" << endl;
  cout << " -edge  : edge strength" << endl;
  cout << " -tfce <H> <E> <connectivity>: enhance with TFCE, e.g. -tfce 2 0.5 6 (maybe change 6 to 26 for skeletons)" << endl;
  cout << " -nan   : replace NaNs (improper numbers) with 0" << endl;
  cout << " -nanm  : make NaN (improper number) mask with 1 for NaN voxels, 0 otherwise" << endl;
  cout << " -inm <mean> :  (-i i ip.c) intensity normalisation (per 3D volume mean)" << endl;
  cout << " -ing <mean> :  (-I i ip.c) intensity normalisation, global 4D mean)" << endl;

  cout << "\nKernel operations:" << endl;
  cout << " -kernel 3D : 3x3x3 box centered on target voxel (set as default kernel)" << endl;
  cout << " -kernel 2D : 3x3x1 box centered on target voxel" << endl;
  cout << " -kernel box    <size>     : all voxels in a box of width <size> centered on target voxel" << endl;
  cout << " -kernel boxv   <size>     : <size>x<size>x<size> box centered on target voxel, CAUTION: size should be an odd number" << endl;
  cout << " -kernel gauss  <sigma>    : gaussian kernel (sigma in mm, not voxels)" << endl;
  cout << " -kernel sphere <size>     : all voxels in a sphere of radius <size> mm centered on target voxel" << endl;
  cout << " -kernel file   <filename> : use external file as kernel" << endl;

  cout << "\nSpatial Filtering operations: N.B. all options apart from -s use the kernel specified by -kernel" << endl;
  cout << " -dilM    : Mean Dilation of zero voxels  (using non-zero voxels in kernel)" << endl;
  cout << " -dilD    : Modal Dilation of zero voxels (using non-zero voxels in kernel)" << endl;
  cout << " -dilF    : Maximum filtering of all voxels" << endl;
  cout << " -ero     : Erode by zeroing non-zero voxels when zero voxels found in kernel" << endl;
  cout << " -eroF    : Minimum filtering of all voxels" << endl;
  cout << " -fmedian : Median Filtering " << endl;
  cout << " -fmean   : Mean filtering, kernel weighted (conventionally used with gauss kernel)" << endl;
  cout << " -fmeanu  : Mean filtering, kernel weighted, un-normalised (gives edge effects)" << endl;
  cout << " -s <sigma> : create a gauss kernel of sigma mm and perform mean filtering" << endl;
  cout << " -subsamp2  : downsamples image by a factor of 2 (keeping new voxels centred on old)" << endl;
  cout << " -subsamp2offc  : downsamples image by a factor of 2 (non-centred)" << endl;

  cout << "\nDimensionality reduction operations:" << endl;
  cout << "  (the \"T\" can be replaced by X, Y or Z to collapse across a different dimension)" << endl;
  cout << " -Tmean   : mean across time" << endl;
  cout << " -Tstd    : standard deviation across time" << endl;
  cout << " -Tmax    : max across time" << endl;
  cout << " -Tmaxn   : time index of max across time" << endl;
  cout << " -Tmin    : min across time" << endl;
  cout << " -Tmedian : median across time" << endl;
  cout << " -Tperc <percentage> : nth percentile (0-100) of FULL RANGE across time" << endl;
  cout << " -Tar1    : temporal AR(1) coefficient (use -odt float and probably demean first)" << endl;

  cout << "\nMulti-argument operations:" << endl;
  cout << " -roi <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize> : zero outside roi" << endl;
  cout << " -bptf  <hp_sigma> <lp_sigma> : (-t in ip.c) Bandpass temporal filtering; nonlinear highpass and Gaussian linear lowpass (with sigmas in volumes, not seconds); set either sigma<0 to skip that filter" << endl;
  cout << " -roc <AROC-thresh> <outfile> [4Dnoiseonly] <truth> : take (normally binary) truth and test current image in ROC analysis against truth. <AROC-thresh> is usually 0.05 and is limit of Area-under-ROC measure FP axis. <outfile> is a text file of the ROC curve (triplets of values: FP TP threshold). If the truth image contains negative voxels these get excluded from all calculations. If <AROC-thresh> is positive then the [4Dnoiseonly] option needs to be set, and the FP rate is determined from this noise-only data, and is set to be the fraction of timepoints where any FP (anywhere) is seen, as found in the noise-only 4d-dataset. This is then controlling the FWE rate. If <AROC-thresh> is negative the FP rate is calculated from the zero-value parts of the <truth> image, this time averaging voxelwise FP rate over all timepoints. In both cases the TP rate is the average fraction of truth=positive voxels correctly found." << endl;

  cout << "\ne.g. fslmaths input_volume -add input_volume2 output_volume" << endl;
  cout << "     fslmaths input_volume -add 2.5 output_volume" << endl;
  cout << "     fslmaths input_volume -add 2.5 -mul input_volume2 output_volume\n" << endl;
}

bool isNumber( const string& x )
{
 bool flag=true;
 char *pend;
 strtod(x.c_str(),&pend);
 if (*pend!='\0') flag=false;
 return flag; 
} 


template <class T>
int fmrib_main(int argc, char *argv[], short output_dt)
{
  volume4D<T> input_volume;
  volumeinfo vinfo;
  //for (int i = 2; i < argc-1; i++) cout << argv[i] << endl;
  volume<float> kernel;
  kernel=box_kernel(3,3,3);
  bool separable=false;
  read_volume4D(input_volume,string(argv[1]),vinfo);

  for (int i = 2; i < argc-1; i++)  //main loop
  {    
    volume4D<T> temp_volume;
    /********************Dimensionality Reduction*******************/
    /***************************************************************/
    if (isupper((int)argv[i][1]) && argv[i][0] == '-')  //if first letters are -capital - dimensionality reduction...
    { 
      int xoff=1,yoff=1,zoff=1,toff=1,nsize;
      if (argv[i][1] == 'T') toff=input_volume.tsize(); 
      if (argv[i][1] == 'Z') zoff=input_volume.zsize();  
      if (argv[i][1] == 'Y') yoff=input_volume.ysize();  
      if (argv[i][1] == 'X') xoff=input_volume.xsize();  
      temp_volume=input_volume;
      input_volume.reinitialize(input_volume.xsize()/xoff,input_volume.ysize()/yoff,input_volume.zsize()/zoff,input_volume.tsize()/toff); 
      input_volume.copyproperties(temp_volume);
      nsize=xoff*yoff*zoff*toff;
      volume<T> column_volume(nsize,1,1); //will be size of appropriate dimension, as only 1 arg is non-unitary
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	    {
              for (int j=0;j<nsize;j++) column_volume.value(j,0,0)=temp_volume(x+j*(xoff!=1),y+j*(yoff!=1),z+j*(zoff!=1),t+j*(toff!=1));
              //This goes along the appropriate axis (non unitary offset variable) and fills a "column" volume with data
	      if (string(argv[i]+2) == "max")    input_volume.value(x,y,z,t)=column_volume.max(); 
	      if (string(argv[i]+2) == "min")    input_volume.value(x,y,z,t)=column_volume.min();
              if (string(argv[i]+2) == "mean")   input_volume.value(x,y,z,t)=(T)column_volume.mean();
              if (string(argv[i]+2) == "std")    input_volume.value(x,y,z,t)=(T)column_volume.stddev();
	      if (string(argv[i]+2) == "maxn")   input_volume.value(x,y,z,t)=column_volume.maxcoordx();
              if (string(argv[i]+2) == "median") input_volume.value(x,y,z,t)=column_volume.percentile(0.5);
	      if (string(argv[i]+2) == "perc")   input_volume.value(x,y,z,t)=column_volume.percentile(atof(argv[i+1])/100.0);
              if (string(argv[i]+2) == "ar1") 
              {
                column_volume-=(T)column_volume.mean();
                double sumsq=column_volume.sumsquares();
                input_volume(x,y,z,t)=0;
		if(sumsq!=0) for (int k=1;k<nsize;k++) input_volume(x,y,z,t)+=(T)(column_volume(k,0,0)*column_volume(k-1,0,0)/sumsq);
	      }
	    }
       if (string(argv[i]+2) == "perc") i++;
       
    }
    /********************Binary Operations**************************/
    /***************************************************************/
    else if (string(argv[i])=="-mas")
    {  
        read_volume4D(temp_volume,string(argv[++i]));
        temp_volume.binarise(0,temp_volume.max()+1,exclusive); // needed to binarise max() value + 1 due to
        for (int t=0;t<input_volume.tsize();t++)                     //potential issue with exclusive binarise
          input_volume[t]*=temp_volume[t%temp_volume.tsize()]; //this gives compatibility with 3 and 4D masks
    }                                                            //without needing to have a specific volume3D variable
    /***************************************************************/
    else if (string(argv[i])=="-thr") input_volume.threshold((T)atof(argv[++i]),input_volume.max()+1,inclusive);
    /***************************************************************/
    else if (string(argv[i])=="-thrp") 
    {
      T lowerlimit =(T)(input_volume.robustmin()+(atof(argv[++i])/100.0)*(input_volume.robustmax()-input_volume.robustmin())); 
      input_volume.threshold(lowerlimit,input_volume.max()+1,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-thrP") 
    {
      volume4D<T> mask(input_volume);
      mask.binarise(0,input_volume.max()+1,exclusive);
      T lowerlimit =(T)(input_volume.robustmin(mask)+(atof(argv[++i])/100.0)*(input_volume.robustmax(mask)-input_volume.robustmin(mask))); 
      input_volume.threshold(lowerlimit,input_volume.max()+1,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-uthr") input_volume.threshold(input_volume.min()-1,(T)atof(argv[++i]),inclusive);
    /***************************************************************/
    else if (string(argv[i])=="-uthrp") 
    {
      T upperlimit = (T)(input_volume.robustmin()+(atof(argv[++i])/100.0)*(input_volume.robustmax()-input_volume.robustmin())); 
       input_volume.threshold(input_volume.min()-1,upperlimit,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-uthrP") 
    {
       volume4D<T> mask(input_volume);
       mask.binarise(0,input_volume.max()+1,exclusive);
       T upperlimit = (T)(input_volume.robustmin(mask)+(atof(argv[++i])/100.0)*(input_volume.robustmax(mask)-input_volume.robustmin(mask))); 
       input_volume.threshold(input_volume.min()-1,upperlimit,inclusive);
    }
    /***************************************************************/
    else if (string(argv[i])=="-kernel") 
    {
       kernel.destroy();
       float xdim=input_volume.xdim();
       float ydim=input_volume.ydim();
       float zdim=input_volume.zdim();
       if(string(argv[i+1])=="2D")      kernel=box_kernel(3,3,1);
       else if(string(argv[i+1])=="3D") kernel=box_kernel(3,3,3);
       else
       {
	 float size=atof(argv[i+2]);
	 if(string(argv[i+1])=="box")      kernel=box_kernel(size,xdim,ydim,zdim);
	 if(string(argv[i+1])=="boxv")   kernel=box_kernel((int)size,(int)size,(int)size);
         else if(string(argv[i+1])=="gauss")  kernel=gaussian_kernel3D(size,xdim,ydim,zdim);
         else if(string(argv[i+1])=="sphere") kernel=spherical_kernel(size,xdim,ydim,zdim);
	 else if(string(argv[i+1])=="file")   read_volume(kernel,string(argv[i+2]));
         if(string(argv[i+1])=="box" || string(argv[i+1])=="gauss") separable=true;       
         else separable=false;  
	  i++;
       }
       i++;
       //save_volume(kernel,"kernel");
       
    }
    /***************************************************************/
    else if (string(argv[i])=="-s") 
    {
      kernel.destroy();
      float xdim=input_volume.xdim();
      float ydim=input_volume.ydim();
      float zdim=input_volume.zdim();
      kernel=gaussian_kernel3D(atof(argv[i+1]),xdim,ydim,zdim);
      separable=true;
      input_volume=generic_convolve(input_volume,kernel,separable,true); 
      i++;
    }
    /***************************************************************/
    else if (string(argv[i])=="-subsamp2"){
      temp_volume.clear();
      temp_volume = subsample_by_2(input_volume,true);
      input_volume = temp_volume;
    }
    /***************************************************************/
    else if (string(argv[i])=="-subsamp2offc"){
      temp_volume.clear();
      temp_volume = subsample_by_2(input_volume,false);
      input_volume = temp_volume;
    }
    /***************************************************************/
    else if (string(argv[i])=="-add"){
      i++;
      if (isNumber(string(argv[i]))) input_volume+=(T)atof(argv[i]); 
      else if (FslFileExists(argv[i])) 
      {  
	read_volume4D(temp_volume,string(argv[i]));
        for (int t=0;t<input_volume.tsize();t++) input_volume[t]+=temp_volume[t%temp_volume.tsize()]; 
      }}
    /***************************************************************/
    else if (string(argv[i])=="-sub"){
      i++;
      if (isNumber(string(argv[i]))) input_volume-=(T)atof(argv[i]);
      else if (FslFileExists(argv[i])) 
      {  
	read_volume4D(temp_volume,string(argv[i]));
        for (int t=0;t<input_volume.tsize();t++) input_volume[t]-=temp_volume[t%temp_volume.tsize()]; 
      }}
    /***************************************************************/
    else if (string(argv[i])=="-mul"){
      i++;
      if (isNumber(string(argv[i]))) input_volume*=(T)atof(argv[i]);
      else if (FslFileExists(argv[i])) 
      {  
	read_volume4D(temp_volume,string(argv[i]));
        for (int t=0;t<input_volume.tsize();t++) input_volume[t]*=temp_volume[t%temp_volume.tsize()]; 
      }}
    /***************************************************************/
    else if (string(argv[i])=="-rem"){
      i++;
      if (isNumber(string(argv[i]))) 
      {
	int denom=(int)atof(argv[i]);
	for(int t=0;t<input_volume.tsize();t++) 
          for(int z=0;z<input_volume.zsize();z++)
	    for(int y=0;y<input_volume.ysize();y++)	    
	      for(int x=0;x<input_volume.xsize();x++)
		input_volume(x,y,z,t)=(int)input_volume(x,y,z,t)%denom;
      }
      else if (FslFileExists(argv[i])) 
      {  
	read_volume4D(temp_volume,string(argv[i]));
        for(int t=0;t<input_volume.tsize();t++) 
	  { 
            int t2=t%temp_volume.tsize();     
            for(int z=0;z<input_volume.zsize();z++)
	      for(int y=0;y<input_volume.ysize();y++)	    
	        for(int x=0;x<input_volume.xsize();x++)
		  if(temp_volume(x,y,z,t2)!=0) input_volume(x,y,z,t)=(int)input_volume(x,y,z,t)%(int)temp_volume(x,y,z,t2); 
	  }
      }
    }
    /***************************************************************/
    else if (string(argv[i])=="-div"){
      i++;
      if (isNumber(string(argv[i]))) {if (atof(argv[i])!=0) input_volume/=(T)atof(argv[i]);}
      else if (FslFileExists(argv[i])) 
      {  
        read_volume4D(temp_volume,string(argv[i]));
        for(int t=0;t<input_volume.tsize();t++)      
	{
          int t2=t%temp_volume.tsize();     
          for(int z=0;z<input_volume.zsize();z++)
            for(int y=0;y<input_volume.ysize();y++)	    
	      for(int x=0;x<input_volume.xsize();x++)
		  if(temp_volume(x,y,z,t2)!=0) input_volume.value(x,y,z,t) /= temp_volume.value(x,y,z,t2);
                  else input_volume.value(x,y,z,t)=(T)0.0;
        }          
      }}
    /***************************************************************/
    else if (string(argv[i])=="-max" || string(argv[i])=="-min")
    {
      T param=0;
      bool max=false;
      bool file=false;
      if (string(argv[i])=="-max") max=true;
      i++;
      if (isNumber(string(argv[i]))) param=(T)atof(argv[i]);
      else if (FslFileExists(argv[i])) 
      {                           
	read_volume4D(temp_volume,string(argv[i]));
        file=true;
      }     
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	    {
              if (max && file) input_volume.value(x,y,z,t)=MAX(input_volume.value(x,y,z,t),temp_volume.value(x,y,z,t%temp_volume.tsize()));
              if (max && !file) input_volume.value(x,y,z,t)=MAX(input_volume.value(x,y,z,t),param);
              if (!max && file) input_volume.value(x,y,z,t)=MIN(input_volume.value(x,y,z,t),temp_volume.value(x,y,z,t%temp_volume.tsize()));
              if (!max && !file) input_volume.value(x,y,z,t)=MIN(input_volume.value(x,y,z,t),param);
            }
    }
    /***************************************************************/
    else if (string(argv[i])=="-roc")
    {
      // {{{ variables

float aroc_thresh = atof(argv[++i]);

ofstream ofs(argv[++i]);

int border=5;

int separatenoise=1;
if (aroc_thresh<0)
  {
    separatenoise=0;
    aroc_thresh*=-1;
  }

// }}}

      // {{{ setup minval, read pure-noise and log-transform it, log-transform input

float minval=input_volume.min();

volume4D<float> noise;
if (separatenoise)
  {
    read_volume4D(noise,string(argv[++i]));
    minval=min(noise.min(),minval);
    for(int t=0;t<input_volume.tsize();t++)
      for(int z=0;z<input_volume.zsize();z++)
	for(int y=0;y<input_volume.ysize();y++)	    
	  for(int x=0;x<input_volume.xsize();x++)
	    noise.value(x,y,z,t)=log(noise.value(x,y,z,t)-minval+1);
  }

volume4D<float> loginput(input_volume.xsize(),input_volume.ysize(),input_volume.zsize(),input_volume.tsize());
for(int t=0;t<input_volume.tsize();t++)
  for(int z=0;z<input_volume.zsize();z++)
    for(int y=0;y<input_volume.ysize();y++)	    
      for(int x=0;x<input_volume.xsize();x++)
	loginput.value(x,y,z,t)=log(input_volume.value(x,y,z,t)-minval+1);

// }}}

#ifdef pooooey
      // {{{ read the truth image

volume<float> truth;
read_volume(truth,string(argv[++i]));

// }}}
      if (separatenoise)
	// {{{ get FP from separate noise image

{
  float aroc=0, FP=0, TP=0, TPprev=0, FPprev=0;

  float maxlogval=max(loginput.max(),noise.max());
  float delta=maxlogval/1000;
  if (delta==0) delta=1;

  // {{{ setup truth and invtruth images

truth.binarise(truth.max() * 0.05);
for(int z=0;z<loginput.zsize();z++)
  for(int y=0;y<loginput.ysize();y++)	    
    for(int x=0;x<loginput.xsize();x++)
      if ((x<border)||(x>=loginput.xsize()-border)||(y<border)||(y>=loginput.ysize()-border)||(z<border)||(z>=loginput.zsize()-border))
	truth.value(x,y,z)=0;
int truecount=truth.sum();

volume4D<float> TPim;
copyconvert(loginput,TPim);
for(int t=0;t<loginput.tsize();t++)           
  TPim[t]*=truth;

// volume<float> invtruth = (truth*-1.0)+1;
// volume4D<float> FPim;
// copyconvert(loginput,FPim);
// for(int t=0;t<loginput.tsize();t++)           
//   FPim[t]*=invtruth;

// }}}
					  
  for (float thresh=maxlogval+delta; thresh>=0 && FP<aroc_thresh; thresh-=delta)
    {
      TP=FP=0;
      for(int t=0;t<loginput.tsize();t++)           
	{
	  int anyFP=0;
	  float sigTP=0, sigFP=0;
	  for(int z=border;z<loginput.zsize()-border;z++)
	    for(int y=border;y<loginput.ysize()-border;y++)	    
	      for(int x=border;x<loginput.xsize()-border;x++)
		{
		  sigTP+=TPim.value(x,y,z,t)>thresh;
		  //sigFP+=FPim.value(x,y,z,t)>thresh;
		  anyFP+=noise.value(x,y,z,t)>thresh;
		}
	  //TP += 2*sigTP / ( truecount + sigFP + sigTP ); // Dice overlap measure
	  TP += sigTP / truecount;
	  FP += anyFP>0;
	}
	      
      TP /= loginput.tsize();
      FP /= loginput.tsize();

      if (FP<aroc_thresh)
	aroc += (TP+TPprev)/2 * (FP-FPprev);
      else  // case for when the latest update straddles the FP threshold
	aroc += (TPprev + (TP-TPprev)*0.5*((aroc_thresh-FPprev)/(FP-FPprev))) * (aroc_thresh-FPprev);
	      
      TPprev=TP;
      FPprev=FP;
      
      ofs << FP << " " << TP << " " << exp(thresh)-1+minval << endl;
    }
  
  if (FP<aroc_thresh) // deal with case of when FP never reached aroc_thresh
    aroc += TP*(aroc_thresh-FP);
  
  ofs.close();
  cout << aroc / aroc_thresh << endl;
}

// }}}
      else
	// {{{ "traditional" ROC

{
  float maxlogval=loginput.max();
  float delta=maxlogval/1000; if (delta==0) delta=1;
  //cout << "minval=" << minval << " maxlogval=" << maxlogval << " delta=" << delta << endl;
  float aroc=0, FP=0, TP=0, TPprev=0, FPprev=0;

  // {{{ setup truth and invtruth images

volume<float> maskim=truth;
maskim.binarise(-0.001);
truth.binarise(truth.max() * 0.05);
int truecount=truth.sum();

volume<float> invtruth = ((truth*-1.0)+1)*maskim;
int invtruecount=invtruth.sum();

volume<float> TPim;
copyconvert(loginput[0],TPim);
TPim*=truth;

volume<float> FPim;
copyconvert(loginput[0],FPim);
FPim*=invtruth;

// }}}

  for (float thresh=maxlogval+delta; thresh>=0 && FP<aroc_thresh; thresh-=delta)
    {
      TP=FP=0;
      for(int z=0;z<loginput.zsize();z++)
	for(int y=0;y<loginput.ysize();y++)	    
	  for(int x=0;x<loginput.xsize();x++)
	    {
	      TP+=TPim.value(x,y,z)>thresh;
	      FP+=FPim.value(x,y,z)>thresh;
	    }
      TP /= truecount;
      FP /= invtruecount;

      if (FP<aroc_thresh)
	aroc += (TP+TPprev)/2 * (FP-FPprev);
      else  // case for when the latest update straddles the FP threshold
	aroc += (TPprev + (TP-TPprev)*0.5*((aroc_thresh-FPprev)/(FP-FPprev))) * (aroc_thresh-FPprev);
      
      TPprev=TP;
      FPprev=FP;
      
      ofs << FP << " " << TP << " " << exp(thresh)-1+minval << endl;
    }

  if (FP<aroc_thresh) // deal with case of when FP never reached aroc_thresh
    aroc += TP*(aroc_thresh-FP);
  
  ofs.close();
  cout << aroc / aroc_thresh << endl;
}

// }}}
#endif

      // {{{ get FP from separate noise image

{
  float aroc=0, FP=0, TP=0, TPprev=0, FPprev=0;
  float maxlogval=loginput.max();
  if (separatenoise)
    maxlogval=max(maxlogval,noise.max());
  float delta=maxlogval/1000;
  if (delta==0) delta=1;
  //cout << "minval=" << minval << " maxlogval=" << maxlogval << " delta=" << delta << endl;

  // {{{ setup truth and invtruth etc, images

// read truth
volume<float> truth;
read_volume(truth,string(argv[++i]));

// mask
volume<float> maskim;
if (!separatenoise)
  {
    maskim=truth;
    maskim.binarise(-0.001);
  }

// truth
for(int z=0;z<loginput.zsize();z++) for(int y=0;y<loginput.ysize();y++) for(int x=0;x<loginput.xsize();x++)
  if ((x<border)||(x>=loginput.xsize()-border)||(y<border)||(y>=loginput.ysize()-border)||(z<border)||(z>=loginput.zsize()-border))
    truth.value(x,y,z)=0;
truth.binarise(truth.max() * 0.05);
int truecount=(int)truth.sum();

// TPim
volume4D<float> TPim;
copyconvert(loginput,TPim);
for(int t=0;t<loginput.tsize();t++)           
  TPim[t]*=truth;

// noise
int invtruecount=1;
if (!separatenoise)
  {
    volume<float> invtruth = ((truth*-1.0)+1)*maskim;
    for(int z=0;z<loginput.zsize();z++) for(int y=0;y<loginput.ysize();y++) for(int x=0;x<loginput.xsize();x++)
      if ((x<border)||(x>=loginput.xsize()-border)||(y<border)||(y>=loginput.ysize()-border)||(z<border)||(z>=loginput.zsize()-border))
	invtruth.value(x,y,z)=0;
    invtruecount=(int)invtruth.sum();
    copyconvert(loginput,noise);
    for(int t=0;t<loginput.tsize();t++)           
      noise[t]*=invtruth;
  }

// }}}
					  
  for (float thresh=maxlogval+delta; thresh>0 && FP<aroc_thresh; thresh-=delta)
    {
      TP=FP=0;
      float FPfwe=0;
      ColumnVector TPvals(loginput.tsize()), FPuncorrected(loginput.tsize());
      TPvals=0; FPuncorrected=0;
      for(int t=0;t<loginput.tsize();t++)           
	{
	  float sigTP=0, sigFP=0;
	  for(int z=border;z<loginput.zsize()-border;z++)
	    for(int y=border;y<loginput.ysize()-border;y++)	    
	      for(int x=border;x<loginput.xsize()-border;x++)
		{
		  sigTP+=TPim.value(x,y,z,t)>=thresh;
		  sigFP+=noise.value(x,y,z,t)>=thresh;
		}
	  TPvals(t+1) = sigTP / truecount;
	  FPuncorrected(t+1) = sigFP / invtruecount;
	  FPfwe += sigFP>0;
	  if (!separatenoise)
	    ofs << FPuncorrected(t+1) << " " << TPvals(t+1) << " " << exp(thresh)-1+minval << endl;
	}
	      
      TP = mean(TPvals).AsScalar();
      if (separatenoise)
	FP = FPfwe / loginput.tsize();
      else
	FP = mean(FPuncorrected).AsScalar();

      if (FP<aroc_thresh)
	aroc += (TP+TPprev)/2 * (FP-FPprev);
      else  // case for when the latest update straddles the FP threshold
	aroc += (TPprev + (TP-TPprev)*0.5*((aroc_thresh-FPprev)/(FP-FPprev))) * (aroc_thresh-FPprev);
	      
      TPprev=TP;
      FPprev=FP;
      
      if (separatenoise)
	ofs << FP << " " << TP << " " << sqrt(var(FPuncorrected).AsScalar()) << " " << sqrt(var(TPvals).AsScalar()) << " " << exp(thresh)-1+minval << endl;
    }
  
  if (FP<aroc_thresh) // deal with case of when FP never reached aroc_thresh
    aroc += TP*(aroc_thresh-FP);
  
  ofs.close();
  cout << aroc / aroc_thresh << endl;
}

// }}}
    }
    /*********************Unary Operations**************************/
    /***************************************************************/
    else if (string(argv[i])=="-sqrt")
    {
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	    {
              if (input_volume.value(x,y,z,t)> 0) input_volume.value(x,y,z,t)=(T)sqrt(input_volume.value(x,y,z,t));
              else  input_volume.value(x,y,z,t)= 0; 
            }
    }
    /***************************************************************/
    else if (string(argv[i])=="-sqr") input_volume*=input_volume;
    /***************************************************************/
    else if (string(argv[i])=="-exp")
    {
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
               input_volume.value(x,y,z,t)=(T)exp((double)input_volume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-log")
    {
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
              if (input_volume.value(x,y,z,t)> 0) input_volume.value(x,y,z,t)=(T)log((double)input_volume.value(x,y,z,t));
    }
    /***************************************************************/
    else if (string(argv[i])=="-abs") input_volume=abs(input_volume);
    /***************************************************************/
    else if (string(argv[i])=="-bin") input_volume.binarise(0,input_volume.max()+1,exclusive); 
    /***************************************************************/
    else if (string(argv[i])=="-index")
    {
      int indexval=0;
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	      if (input_volume.value(x,y,z,t)>0) 
		{
		  input_volume.value(x,y,z,t)=(T)indexval;
		  indexval++;
		}
    }
    /***************************************************************/
    else if (string(argv[i])=="-grid")
    {
      double gridvalue = atof(argv[++i]);
      int gridspacing = atoi(argv[++i]);
      for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
	      if ( x%gridspacing==0 || y%gridspacing==0 || z%gridspacing==0 )
		input_volume.value(x,y,z,t)=(T)gridvalue;
    }
    /*****************SPATIAL FILTERING OPTIONS*********************/
    /***********************Mean Dilation***************************/
    else if (string(argv[i])=="-dilM")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"dilateM");
    /***********************Modal Dilation**************************/
    else if (string(argv[i])=="-dilD")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"dilateD");
    /***********************MJ Dilation**************************/
    else if (string(argv[i])=="-dilF")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"dilate");
    /***********************Steves Erosion**************************/
    else if (string(argv[i])=="-ero")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"erodeS");
    /**************************MJ Erosion**************************/
    else if (string(argv[i])=="-eroF")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"erode");
    /***********************Median Filtering***********************/
    else if (string(argv[i])=="-fmedian")
	for(int t=0;t<input_volume.tsize();t++) input_volume[t]=morphfilter(input_volume[t],kernel,"median");
    /******************Mean Filtering*************************/
    else if (string(argv[i])=="-fmean")
       input_volume=generic_convolve(input_volume,kernel,separable,true);
    /******************Mean Filtering Unnormalised************/
    else if (string(argv[i])=="-fmeanu")
       input_volume=generic_convolve(input_volume,kernel,false,false);
    /*****************END OF FILTERING OPTIONS***************/
    else if (string(argv[i])=="-edge")
       input_volume=edge_strengthen(input_volume);
    else if (string(argv[i])=="-tfce")
      {
	float height_power = atof(argv[++i]);
	float size_power = atof(argv[++i]);
	int connectivity = atoi(argv[++i]);

	for(int t=0;t<input_volume.tsize();t++)           
	  {
	    float maxval=input_volume[t].max();
	    volume<float> clusterenhance;
	    copyconvert(input_volume[0],clusterenhance);
	    clusterenhance=0;

	    float delta=0.1;  // this needs fixing!!
	    for (float thresh=delta; thresh<=maxval; thresh+=delta)
	      {
		volume<float> clusters;
		copyconvert(input_volume[t],clusters);
		clusters.binarise(thresh);

		ColumnVector clustersizes;  
		volume<int>tmpvol=connected_components(clusters,clustersizes,connectivity);
		//int xx=112,yy=144,zz=74;
		//if (tmpvol(xx,yy,zz)>0) cout << clustersizes(tmpvol(xx,yy,zz)) << " ";
		clustersizes = pow(clustersizes,size_power) * pow(thresh,height_power);
		//if (tmpvol(xx,yy,zz)>0) cout << clustersizes(tmpvol(xx,yy,zz)) << endl;
		for(int z=0;z<input_volume.zsize();z++)
		  for(int y=0;y<input_volume.ysize();y++)	    
		    for(int x=0;x<input_volume.xsize();x++)
		      if (tmpvol.value(x,y,z)>0)
			clusterenhance.value(x,y,z) += clustersizes(tmpvol.value(x,y,z));
	      }
	    copyconvert(clusterenhance,input_volume[t]);
	  }
      }
    /******************************************************/
    else if (string(argv[i])=="-nanm")
     {   
       for(int t=0;t<input_volume.tsize();t++)           
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               if ( finite((double)input_volume.value(x,y,z,t))) input_volume.value(x,y,z,t)=0;
	       else input_volume.value(x,y,z,t)=1;
     }
     /******************************************************/
    else if (string(argv[i])=="-nan")
     {   
       for(int t=0;t<input_volume.tsize();t++)           
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               if (!finite((double)input_volume.value(x,y,z,t))) input_volume.value(x,y,z,t)=0;     
     }
     /******************************************************/
    else if (string(argv[i])=="-roi")
     {   
       for(int t=0;t<input_volume.tsize();t++)           
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               if((x<atoi(argv[i+1])) || (x>=atoi(argv[i+1])+atoi(argv[i+2])) || (y<atoi(argv[i+3])) || (y>=atoi(argv[i+3])+atoi(argv[i+4])) || (z<atoi(argv[i+5])) || (z>=atoi(argv[i+5])+atoi(argv[i+6])) || (t<atoi(argv[i+7])) || (t>=atoi(argv[i+7])+atoi(argv[i+8])) )
                 input_volume.value(x,y,z,t)=0;
       i+=8;
     }
     /*******************IP functions***********************/
    else if (string(argv[i])=="-inm")
     { 
       double target,tmean;
       target = atof(argv[++i]);
       volume4D<T> mask(input_volume);   
       mask.binarise(0,mask.max()+1,exclusive); 
       for(int t=0;t<input_volume.tsize();t++)     
       {
         tmean=target/input_volume[t].mean(mask[t]);
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               input_volume.value(x,y,z,t)=(T)(input_volume.value(x,y,z,t)*tmean);
       }
     }
    else if (string(argv[i])=="-ing")
     { 
       double tmean,target;
       target = atof(argv[++i]);
       volume4D<T> mask(input_volume);   
       mask.binarise(0,mask.max()+1,exclusive); 
       tmean=target/input_volume.mean(mask);
       for(int t=0;t<input_volume.tsize();t++)     
         for(int z=0;z<input_volume.zsize();z++)
           for(int y=0;y<input_volume.ysize();y++)	    
	     for(int x=0;x<input_volume.xsize();x++)
               input_volume.value(x,y,z,t)=(T)(input_volume.value(x,y,z,t)*tmean);
     }
    else if(string(argv[i])=="-bptf")
     {
       input_volume=bandpass_temporal_filter(input_volume,atof(argv[i+1]),atof(argv[i+2]));
       i+=2;
     }
    else { cout << "\n Error in command line: unknown option \"" << argv[i] << "\"\n" << endl; print_usage("blah"); return 1; }
     /******************************************************/
  } 

  if (dtype(input_volume)>=DT_FLOAT && output_dt < DT_FLOAT)
  {
    for(int t=0;t<input_volume.tsize();t++)           
        for(int z=0;z<input_volume.zsize();z++)
          for(int y=0;y<input_volume.ysize();y++)	    
	    for(int x=0;x<input_volume.xsize();x++)
              input_volume.value(x,y,z,t)=(T) MISCMATHS::round(input_volume.value(x,y,z,t));
  }

  FslSetCalMinMax(&vinfo,input_volume.min(),input_volume.max());
  save_volume4D_dtype(input_volume,string(argv[argc-1]),output_dt,vinfo,true);
  return 0;
}


int main(int argc,char *argv[])
{
  if (argc < 2) 
  { 
    print_usage(string(argv[0]));
    return 1; 
  }
  short original_dt;
  if(string(argv[1]) =="-datatype" || string(argv[1])== "-dt")  original_dt = dtype(string(argv[3]));
  else original_dt = dtype(string(argv[1]));
  short output_dt=original_dt;
  if(string(argv[argc-2])=="-output_datatype" || string(argv[argc-2])== "-odt") //output datatype
  {
    if(string(argv[argc-1])=="char")        output_dt =  DT_UNSIGNED_CHAR;  
    else if(string(argv[argc-1])=="short")  output_dt =  DT_SIGNED_SHORT;
    else if(string(argv[argc-1])=="int")    output_dt =  DT_SIGNED_INT;
    else if(string(argv[argc-1])=="float")  output_dt =  DT_FLOAT;
    else if(string(argv[argc-1])=="double") output_dt =  DT_DOUBLE;
    else {cout << "Error: Unknown datatype \"" << argv[argc-1] << "\" - Possible datatypes are: char short int float double" << endl; return 1;}
    argc-=2;
  }
  if(string(argv[1])=="-datatype" || string(argv[1])== "-dt") //input datatype
  {
    short input_dt=-1;
     if(string(argv[2])=="input") input_dt=original_dt;
     else if(string(argv[2])=="char" || input_dt == DT_UNSIGNED_CHAR)     return fmrib_main<char>(argc-2, argv+2,output_dt);
     else if(string(argv[2])=="short" || input_dt == DT_SIGNED_SHORT)return fmrib_main<short>(argc-2, argv+2,output_dt);
     else if(string(argv[2])=="int" || input_dt == DT_SIGNED_INT)    return fmrib_main<int>(argc-2, argv+2,output_dt);
     else if(string(argv[2])=="float" || input_dt == DT_FLOAT)   return fmrib_main<float>(argc-2, argv+2,output_dt); 
     else if(string(argv[2])=="double" || input_dt == DT_DOUBLE) return fmrib_main<double>(argc-2, argv+2,output_dt); 
     else {cout << "Error: Unknown datatype \"" << argv[2] <<  "\" - Possible datatypes are: char short int float double input" << endl; return 1;}
  }
  else if (dtype(string(argv[1]))==DT_DOUBLE) return fmrib_main<double>(argc,argv,output_dt);
  else return fmrib_main<float>(argc,argv,output_dt);
}



