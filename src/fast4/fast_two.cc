/*  FAST4 - FMRIB's Automated Segmentation Tool v4

    John Vickers, Mark Jenkinson and Steve Smith
    FMRIB Image Analysis Group

    Copyright (C) 2005-2007 University of Oxford  */

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

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "mriseg_two.h"
#include "multi_mriseg_two.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string title="FAST4 (Version 4.0)\nCopyright(c) 2004-7, University of Oxford (John Vickers)";
string examples="fast4 [options] file(s)";
string examples_multi_channel="fast4 [options] <image> [<image2> ... <imagen>]";
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Each (global) object below specificies as option and can be accessed
// anywhere in this file (since they are global).  The order of the
// arguments needed is: name(s) of option, default value, help message,
// whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
// will not be active.


// following not currently user-settable
// Option<float> fsmall(string("-s, --small"), 0.0,
// 		  string("small positive number for EM convergence; default=0.0"),
// 		  false, requires_argument);

Option<int> inititer(string("-W, --init"), 15,
		  string("number of segmentation-initialisation iterations; default=15"),
		  false, requires_argument);
Option<int> nbiter(string("-I, --iter"), 4,
		  string("number of main-loop iterations during bias-field removal; default=4"),
		  false, requires_argument);
Option<int> initfixity(string("-O, --fixed"), 4,
		  string("number of main-loop iterations after bias-field removal; default=4"),
		  false, requires_argument);


Option<float> fbeta(string("-f, --fHard"), 0.02,
		  string("initial segmentation spatial smoothness (during bias field estimation); default=0.02"),
		  false, requires_argument);
Option<float> Hyp(string("-H, --Hyper"), 0.3,
		  string("segmentation spatial smoothness; default=0.3, set < 0 for automatic estimation"),
		  false, requires_argument);
Option<float> fpveMRFmixeltype(string("-R, --mixel"), 0.1,
		  string("spatial smoothness for mixeltype; default=0.1"),
		  false, requires_argument);
Option<float> nblowpass(string("-l, --lowpass"), 20,
		  string("bias field smoothing extent (FWHM) in mm; default=20"),
		  false, requires_argument);

Option<int> typeofimage(string("-t, --type"), 0,
		  string("type of image 0=T1, 1=T2, 2=PD; default=T1"),
		  false, requires_argument);
Option<int> nclass(string("-n, --class"), 3,
		  string("number of tissue-type classes; default=3"),
		  false, requires_argument);
Option<string> inname1(string("-i, --in1"), string(""),  //not used by user any more but still needed internally
		  string("first input filename"),
		  false, requires_argument);
Option<string> outname(string("-o, --out"), string(""),
		  string("output basename"),
		  false, requires_argument);
Option<int> nchannel(string("-S, --channels"), 1,
		  string("number of input images (channels); default 1"),
		  false, requires_argument);
Option<bool> multichannel(string("-m, --multi"), false, //not used by user any more but still needed internally
		  string("uses multi channel segmentation"),
		  false, no_argument);


Option<bool> nopve(string("--nopve"), false,
		  string("turn off PVE (partial volume estimation)"),
		  false, no_argument);
Option<int> pve(string("--pvestep"), 100,
		  string("discretisation levels of pve values; default=100"),
		  false, requires_argument);

Option<bool> segments(string("-g, --segments"), false,
		  string("outputs a separate binary image for each tissue type"),
		  false, no_argument);

Option<bool> out_probs(string("-r, --outprobs"), false, //not used by user any more but still needed internally
		  string("outputs individual probability maps"),
		  false, no_argument);

Option<bool> biasrem(string("-N, --nobias"), false,
		  string("do not remove bias field"),
		  false, no_argument);
Option<bool> out_bias(string("-b, --outbias"), false,
		  string("output bias field image"),
		  false, no_argument);

Option<string> bapriori(string("-a"), "",
		  string("~<standard2input.mat> initialise using priors; you must supply a FLIRT transform"),
		  false, requires_argument);
Option<bool> talaraichiterations(string("-P, --Prior"), false,
		  string("use priors throughout; you must also set the -a option"),
		  false, no_argument);
Option<string> anotherstdspace(string("-A"), "",
		  string("~<prior1> <prior2> <prior3>    alternative prior images"),
		  false, requires_3_arguments);


Option<string> manualsegmentation(string("-s, --manualseg"), "",
		  string("~<filename> Filename containing intensities"),
		  false, requires_argument);

Option<bool> verbose(string("-v, --verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h, --help"), false,
		  string("display this message"),
		  false, no_argument);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Local functions
int prior_registration(string inname, string main_prior_vol, NEWIMAGE::volume<float>& pCSF, NEWIMAGE::volume<float>& pGM, NEWIMAGE::volume<float>& pWM)
{
  string name1, name2, name3;
  if(anotherstdspace.unset())
    {
      string fname=string(getenv("FSLDIR")) + "/data/standard/tissuepriors/avg152T1_";
      name1 = fname+"csf";
      name2 = fname+"gray";
      name3 = fname+"white";
    } else {
      name1 = anotherstdspace.value(0);
      name2 = anotherstdspace.value(1);
      name3 = anotherstdspace.value(2);
    }

  int bapused=0;
  if((bapriori.value()!=""))
      bapused= 1;
  if((nclass.value()!=2)&&(nclass.value()!=3)&(bapused!=0))
    {
      bapused=0;
      cerr<< "Apriori can only be used for 2 or 3-class segmentation\n";
    }
  if(bapused>0)
    {
      string GMpath, WMpath, CSFpath;

      if(fsl_imageexists(name1))
	read_volume(pCSF, name1);
      else
        {
	  cerr<< "prior image " << name1 << " is not found! priors are not used!\n";
	  bapused = 0;
        }

      if(fsl_imageexists(name2))
	read_volume(pGM, name2);
      else
        {  
	  cerr<< "prior image " << name2 << " is not found! priors are not used!\n";
	  bapused = 0;
        }

      if(fsl_imageexists(name3))
	read_volume(pWM, name3);
      else
        {  
	  cerr<< "prior image " << name3 << " is not found! priors are not used!\n";
	  bapused = 0;
        }
    }
    if(bapused>0)
      {
	char reg[1024];
	sprintf(reg, "%s/bin/flirt -ref %s -in %s -out %s -applyxfm -init %s", getenv("FSLDIR"), inname.c_str(), name1.c_str(), (main_prior_vol+"_csf_stdspace").c_str(),  bapriori.value().c_str());
        if(verbose.value())
	  cout<<reg<<endl;
	system(reg);
	sprintf(reg, "%s/bin/flirt -ref %s -in %s -out %s -applyxfm -init %s", getenv("FSLDIR"), inname.c_str(), name2.c_str(), (main_prior_vol+"_gm_stdspace").c_str(),  bapriori.value().c_str());
	if(verbose.value())
	  cout<<reg<<endl;
        system(reg);
	sprintf(reg, "%s/bin/flirt -ref %s -in %s -out %s -applyxfm -init %s", getenv("FSLDIR"), inname.c_str(), name3.c_str(), (main_prior_vol+"_wm_stdspace").c_str(),  bapriori.value().c_str());
	if(verbose.value())
	  cout<<reg<<endl;
	system(reg);
      }
      if(bapused>0)
	{
	  if(fsl_imageexists((main_prior_vol+"_csf_stdspace")))
	    {
	      read_volume(pCSF, (main_prior_vol+"_csf_stdspace"));
	    }
	  else
	    {  
	      cerr << "csf prior image not transformed correctly! priors are not used!\n";
	      bapused = 0;
	      return -1;
	    }
	  if(fsl_imageexists(main_prior_vol+"_gm_stdspace"))
	    {
	      read_volume(pGM, main_prior_vol+"_gm_stdspace");
	    }
	  else
	    {  
	      cerr << "grey matter prior image not transformed correctly! priors are not used!\n";
	      bapused = 0;
	      return -1;
	    }
	  if(fsl_imageexists(main_prior_vol+"_wm_stdspace"))
	    {
	      read_volume(pWM, main_prior_vol+"_wm_stdspace");
	    }
	  else
	    {  
	      cerr << "white matter prior image not transformed correctly! priors are not used!\n";
	      bapused = 0;
	      return -1;
	    }
	  if(talaraichiterations.value())
	    {
	      bapused=2;
	    }

	}
 
      else
	{
	  pCSF=volume<float>();
	  pGM=volume<float>();
	  pWM=volume<float>();
	}

      return bapused;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Single channel main call

int do_work(int argc, char* argv[]) 
{ 
  float pixdim[3];
  int bapused=0;
  volume<float> pCSF, pGM, pWM;
  int width, height, depth;
  int al=0;

  if (verbose.value()) { cout << "Starting Single Image Segmentation" << endl; }

  inname1.set_value(string(argv[argc-1]));
  if (outname.unset()) {
    outname.set_value(inname1.value());
  }

  ZMRISegmentation	mri;
  volume<float> inputimage;
  volumeinfo inputinfo;
  if(read_volume(inputimage,inname1.value(),inputinfo)!=0)
    {
      cerr<<"Image cannot be found";
      return 1;
    }
  width=inputimage.xsize();
  height=inputimage.ysize();
  depth=inputimage.zsize();
  pixdim[0]=inputimage.xdim();
  pixdim[1]=inputimage.ydim();
  pixdim[2]=inputimage.zdim();
  if(inputimage.min()<0.0)
    inputimage-=inputimage.min();
  if(verbose.value()) 
    { 
      switch(typeofimage.value())
	{
	default:
	case 0:
	  if(verbose.value())
	    cout<< "T1-weighted image" << endl;
	  break;
	case 1:
	  if(verbose.value())
	    cout << "T2-weighted image" << endl;
	  break;
	case 2:
	  if(verbose.value())
	    cout << "PD-weighted image" << endl;
	  break;
	}
      if(verbose.value())
	{
	  cout<< "Imagesize : " << width << " x " << height << " x " << depth << endl;
	  cout<< "Pixelsize : " << pixdim[0] << " x " << pixdim[1] << " x "  << pixdim[2] << endl << endl;
	}
    }
  bapused=prior_registration(inname1.value(),outname.value(), pCSF, pGM, pWM);
  if(bapused==0)
    al=0;
  bool biasf=true;
  if(biasrem.value())
    biasf=false;
  mri.TanakaCreate(inputimage, fbeta.value(), nclass.value(), nblowpass.value(), biasf, pve.value(), fpveMRFmixeltype.value(), nbiter.value(),initfixity.value(), inititer.value(), bapused, Hyp.value(), verbose.value(),manualsegmentation.value(),typeofimage.value());
  if (mri.TanakaMain(pCSF, pGM, pWM)) return -1;
  save_volume(mri.m_Segment,outname.value()+"_seg");
  
  if(segments.value())
    {
      volume<int> ind_segments;
      for(int i=1; i<=nclass.value(); i++)
	{
	  ind_segments=mri.m_Segment;
	  for(int z=0;z<depth;z++)
	    {
	      for(int y=0;y<height;y++)
		{
		  for(int x=0;x<width;x++)
		    {
		      ind_segments.value(x, y, z) = (ind_segments.value(x, y, z)==i);
		    }
		}
	    }
	  save_volume(ind_segments,outname.value()+"_seg_"+num2str(i),inputinfo); 
	}
    }
  
  if(out_probs.value())
    {
      for(int i=1; i<=nclass.value(); i++)
	{
	  save_volume(mri.members[i-1],outname.value()+"_prob_"+num2str(i),inputinfo); 
	}
    }
  
  if(pve.value())
    {
      for(int i=1; i<=nclass.value(); i++)
	{
	  save_volume(mri.m_pve[i],outname.value()+"_pve_"+num2str(i),inputinfo); 
	}
      save_volume(mri.m_pveSegment,outname.value()+"_pveseg",inputinfo);
      save_volume(mri.hardPV, outname.value()+"_mixeltype", inputinfo);	  
    }
  
  if(out_bias.value())
    save_volume(mri.m_Finalbias,outname.value()+"_bias",inputinfo);
  
  return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Multi channel main call

int do_work_multi(int argc, char* argv[]) 
{ 

  float pixdim[3];
  int bapused=0;
  volume<float> pCSF, pGM, pWM;
  int width, height, depth;

  if (verbose.value()) { cout << "Starting Multi Image Segmentation" << endl; }

  if(nchannel.value()>=2)
    {
     int al=0;

      ZMRIMULTISegmentation mri;

      volume<float>* images=new volume<float>[nchannel.value()];
      volumeinfo inputinfo;

      for(int c=0;c<=nchannel.value()-1;c++)
	{
	  if (c==0) {
	    inname1.set_value(string(argv[argc-c-1]));
	    if (outname.unset()) {
	      outname.set_value(inname1.value());
	    }
	  }
	  if(read_volume(images[c], argv[argc-c-1], inputinfo)!=0)
	    {
	      cerr<<"Image cannot be found";
	      return 1;
	    }
	  else
	    {
	      if(images[c].min()<0.0)
		 images[c]-=images[c].min();
	    }
	}
      width=images[0].xsize();
      height=images[0].ysize();
      depth=images[0].zsize();
      pixdim[0]=images[0].xdim();
      pixdim[1]=images[0].ydim();
      pixdim[2]=images[0].zdim();
      bapused=prior_registration(inname1.value(),outname.value(), pCSF, pGM, pWM);
      if(bapused==0)
	al=0;
      bool biasf=true;
      if(biasrem.value())
	biasf=false;

      mri.TanakaCreate(images, nclass.value(), false, nbiter.value(), nblowpass.value(), fbeta.value(), bapused, pve.value(), nchannel.value(), al, biasf,initfixity.value(), verbose.value(), pve.value(), inititer.value(),fpveMRFmixeltype.value(), Hyp.value(),manualsegmentation.value(),typeofimage.value());
      if (mri.TanakaMain(pCSF, pGM, pWM)) return -1;
      


      save_volume(mri.m_Segment, outname.value()+"_seg");

      if(segments.value())
	{
	  volume<int> ind_segments;
	  for(int i=1; i<=nclass.value(); i++)
	    {
	      ind_segments=mri.m_Segment;
	      for(int z=0;z<depth;z++)
		{
		  for(int y=0;y<height;y++)
		    {
		      for(int x=0;x<width;x++)
			{
			  ind_segments.value(x, y, z) = (ind_segments.value(x, y, z)==i);
			}
		    }
		}
	      save_volume(ind_segments,outname.value()+"_seg_"+num2str(i),inputinfo); 
	    }
	}

      if(out_probs.value())
	{
	  for(int i=1; i<=nclass.value(); i++)
	    {
	      save_volume(mri.members[i-1],outname.value()+"_prob_"+num2str(i),inputinfo); 
	    }
	}
      if(pve.value())
      {
	  volume<float> ind_pve;
	  for(int i=1; i<=nclass.value(); i++)
	    {
	      save_volume(mri.m_pve[i],outname.value()+"_pve_"+num2str(i),inputinfo); 
	    }
	  save_volume(mri.m_pveSegment,outname.value()+"_pveseg",inputinfo);
	  save_volume(mri.hardPV, outname.value()+"_mixeltype", inputinfo);
     }
     if(out_bias.value())
       {
	 for(int i=1;i<nchannel.value()+1;i++)
	   {
	     save_volume(mri.m_Finalbias[i],outname.value()+"_bias_"+num2str(i),inputinfo);
	   }
       }


      delete[] images;
      return 0;   
    }
  else cerr<<"At least 2 channels required for Multi Channel Segmentation";
  return 1;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  Tracer tr("main");
  int nonoptarg;
  OptionParser options(title, examples);
  try 
    {
      // must include all wanted options here (the order determines how
      // the help message is printed)
      //      options.add(inname1);
      options.add(nclass);
      options.add(nbiter);;
      options.add(nblowpass);
      options.add(typeofimage);
      options.add(fbeta);
      options.add(segments);
      options.add(bapriori);
      options.add(anotherstdspace);
      options.add(nopve);
      options.add(out_bias);
      options.add(biasrem);
      options.add(nchannel);
      options.add(outname);
      options.add(talaraichiterations);
      options.add(inititer);
      options.add(fpveMRFmixeltype);
      options.add(initfixity);
      options.add(Hyp);
      options.add(verbose);
      options.add(help);
      options.add(manualsegmentation);
      // line below stops the program if the help was requested or 
      // a compulsory option was not set
      if ( argc<2 || (help.value()) || (!options.check_compulsory_arguments(true)) )
	{
	  options.usage();
	  exit(EXIT_FAILURE);
	}

      nonoptarg = options.parse_command_line(argc, argv);
      if (nchannel.value()>1) multichannel.set_value("true");

      if (nopve.value()) pve.set_value("0");
    }
  catch(X_OptionError& e)
    {
      options.usage();
      cerr << endl << e.what() << endl;
      exit(EXIT_FAILURE);
    } 
  catch(std::exception &e)
    {
      cerr << e.what() << endl;
    } 

  // Call the local functions
  if(!multichannel.value())
    return do_work(argc,argv);
  else
    return do_work_multi(argc, argv);
}

