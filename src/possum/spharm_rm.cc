#include <iostream>
#include <string>
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"
  
using namespace NEWIMAGE;
using namespace MISCMATHS;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="spharm_rm (Version 2.0)\nCopyright(c) 2006, University of Oxford (Mark Jenkinson)";
string examples="spharm_rm [options] -i <input_image> -o <output_image>";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> unmasked(string("--unmasked"), false,
		  string("do not mask final removal"),
		  false, no_argument);
Option<int> numterms(string("-n"), 9,
		  string("number of terms to remove (order is 1,x,y,z,z^2+(x^2+y^2)/2,zx,zy,xy,x^2-y^2)"),
		  false, requires_argument);
Option<string> inname(string("-i,--in"), string(""),
		  string("input filename"),
		  true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename"),
		  true, requires_argument);
Option<string> maskname(string("-m"), string(""),
		  string("mask filename"),
		  false, requires_argument);
int nonoptarg;


/////////////////////////////////////////////////////////////////////////////

volume4D<float> generate_spherical_harmonics(const volume<float>& ref, int n=9)
{
  if (n<=0) {
    cerr << "WARNING::generate_spherical_harmonics:: n must be > 0" << endl;
    volume4D<float> rubbish;
    return rubbish;
  }
  volume4D<float> confounds(ref.xsize(),ref.ysize(),ref.zsize(),n);
  float xdim = ref.xdim(), ydim = ref.ydim(),zdim = ref.zdim();
  confounds.setdims(xdim,ydim,zdim,1.0);
  int xc = ref.xsize()/2;
  int yc = ref.ysize()/2;
  int zc = ref.zsize()/2;
  for (int z=0; z<ref.zsize(); z++) {
    for (int y=0; y<ref.ysize(); y++) {
      for (int x=0; x<ref.xsize(); x++) {
	// Set the spherical harmonics: x,y,z,z^2+(x^2+y^2)/2,zx,zy,xy,x^2-y^2
	float x0=(x-xc)*xdim, y0=(y-yc)*ydim, z0=(z-zc)*zdim;
	if (n>0)  confounds(x,y,z,0) = 1.0;
	if (n>1)  confounds(x,y,z,1) = x0;
	if (n>2)  confounds(x,y,z,2) = y0;
	if (n>3)  confounds(x,y,z,3) = z0;
	if (n>4)  confounds(x,y,z,4) = z0*z0 + 0.5*(x0*x0+y0*y0);
	if (n>5)  confounds(x,y,z,5) = z0*x0;
	if (n>6)  confounds(x,y,z,6) = z0*y0;
	if (n>7)  confounds(x,y,z,7) = x0*y0;
	if (n>8)  confounds(x,y,z,8) = x0*x0-y0*y0;
      }
    }
  }
  return confounds;
}


ColumnVector fit_functions(const volume<float>& invol, 
			   const volume4D<float>& confounds)
{			   
  int num = confounds.tsize();
  ColumnVector params(num), xty(num);
  Matrix xtx(num,num);
  volume<float> dotprod;
  // Form the quantities (X' * X) and (X' * Y)
  for (int n=0; n<confounds.tsize(); n++) {
    for (int m=n; m<confounds.tsize(); m++) {
      dotprod = confounds[n] * confounds[m];
      xtx(m+1,n+1) = dotprod.sum();
      xtx(n+1,m+1) = xtx(m+1,n+1);
    }
    dotprod = confounds[n] * invol;
    xty(n+1) = dotprod.sum();
  }
  // Can now get the parameters by solving the GLM
  //   That is:   B = (X' * X)^{-1} * (X' * Y)
  //   Use the pseudo-inverse just in case there is linear dependency
  //    as it should still give a sensible (though not unique) answer
  params = pinv(xtx) * xty;
  return params;
}


volume<float> remove_fitted_functions(const volume<float>& invol,
				      const volume4D<float>& confounds,
				      const ColumnVector& params)
{
  if (confounds.tsize() != params.Nrows() ) {
    cerr << "WARNING:: Do not have the correct number of parameters required" 
	 << endl;
    cerr << "  Removing only a subset of functions" << endl;
  }
  volume<float> outvol = invol;
  int num = Min(confounds.tsize(),params.Nrows());
  for (int n=0; n<num; n++) { 
    outvol -= ((float) params(n+1)) * confounds[n];
  }
  return outvol;
}


/////////////////////////////////////////////////////////////////////////////


int do_work(int argc, char* argv[])
{
  volume<float> invol, outvol, mask;
  volume4D<float> confounds;

  read_volume(invol,inname.value());
  if (maskname.set()) {
    read_volume(mask,maskname.value());
  } else {
    mask = invol*0.0f + 1.0f;
  }
  outvol = invol;


  if (verbose.value()) {cout << "Generating Spherical Harmonic Terms" << endl;}
  confounds = generate_spherical_harmonics(invol,numterms.value());
  for (int n=0; n<confounds.tsize(); n++) { confounds[n] *= mask; }
  if (verbose.value()) {cout << "Fitting Spherical Harmonics" << endl;}
  ColumnVector fits;
  fits = fit_functions(invol*mask,confounds);
  if (verbose.value()) {
    cout << "Spherical Harmonic Amplitudes: ";
    for (int n=1; n<=fits.Nrows(); n++) { cout << fits(n) << "  "; }  
    cout << endl;
  }
  if (unmasked.value()) {
    if (verbose.value()) {cout<<"Regenerating Spherical Harmonic Terms"<<endl;}
    confounds = generate_spherical_harmonics(invol,numterms.value());
  }
  if (verbose.value()) {cout << "Removing Spherical Harmonics" << endl;}
  outvol = remove_fitted_functions(invol,confounds,fits);

  if (!unmasked.value()) { outvol *= mask; }
  save_volume(outvol,outname.value());
  return 0;
}


/////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(outname);
    options.add(maskname);
    options.add(numterms);
    options.add(verbose);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  // Call the local functions

  return do_work(argc,argv);
}

