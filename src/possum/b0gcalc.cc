
// POSSUM

#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#endif //USE_MPI

#include "libprob.h"
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "possumfns.h"
#include "utils/options.h"
#include "newimage/costfns.h"
#include "miscmaths/miscmaths.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace NEWIMAGE;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace Utilities;
//using namespace std;

string title="possum (Version 2.0)\nCopyright(c) 2003, University of Oxford (Ivana Drobnjak)";
string examples="possum -i <input phantom volume> -x <tissue matrix> -r <RFrec> -s <RFtrans> -p <pulse> -f <rf slice profile> -o <output signal matrix> [optional options]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);

//INPUT object and its characteristics (including susc effects on B0 and the RF inhomogeneities)

Option<string> opt_b0(string("-b,--b0p"), string(""),
		  string("3D volume, B0 perturbation due to susceptibility for each voxel (x,y,z) (units in T),in case when rotation Rx or Ry present you need the whole dir"),
		  false, requires_argument);
Option<string> opt_b0timecourse4D(string("--b0time"),string(""),
	          string("1 column b0_timecourse [time(s)]"),
	          false,requires_argument);

Option<int>    opt_level(string("-l,--lev"), 1,
		  string("levels: 1.no motion//basic B0 2.motion//basic B0, 3.motion//full B0, 4.no motion//time changing B0"),
		  true,requires_argument);

int nonoptarg;

/////////////////////////////////////////////////////////////////////////////////////////////////////
int compute_volume(int argc, char *argv[])
{
  int level=opt_level.value();
  ////////////////////////
  //NO MOTION or MOTION in plane
  ///////////////////////
  if (level<=2){
    ////////////////////////
    // B0 PERTURBATION 
    ////////////////////////
    cout<<"Creating 1 B0file together with 3 gradient files..."<<endl;
    volume<double> b0z_dz;
    read_volume(b0z_dz,opt_b0.value()+"z_dz"); 
    print_volume_info(b0z_dz,"b0z_dz");
    volume<double> b0z_dz_gx;
    volume<double> b0z_dz_gy;
    volume<double> b0z_dz_gz;
    calc_gradients(b0z_dz,b0z_dz_gx,b0z_dz_gy,b0z_dz_gz);
    save_volume(b0z_dz_gx,opt_b0.value()+"z_dz_gx");
    save_volume(b0z_dz_gy,opt_b0.value()+"z_dz_gy");
    save_volume(b0z_dz_gz,opt_b0.value()+"z_dz_gz");
    print_volume_info(b0z_dz_gz,"b0z_dz_gz");
 
  ////////////////////////////////////////////////////////////
  //MOTION INVOLVING ROTATION Rx or Ry or both
  ////////////////////////////////////////////////////////////
  if (level==3){
    ////////////////////////////////////////////
    // B0 PERTURBATION 
    ////////////////////////////////////////////
    cout<<"B0 perturbation volumes..."<<endl;
    volume<double> b0x_dx, b0x_dy, b0x_dz, b0y_dx, b0y_dy, b0y_dz, b0z_dx, b0z_dy, b0z_dz;//read in
    volume<double> b0x_dx_gx, b0x_dx_gy, b0x_dx_gz, b0x_dy_gx, b0x_dy_gy, b0x_dy_gz, b0x_dz_gx, b0x_dz_gy, b0x_dz_gz;//calculate from, the calc gradients
    volume<double> b0y_dx_gx, b0y_dx_gy, b0y_dx_gz, b0y_dy_gx, b0y_dy_gy, b0y_dy_gz, b0y_dz_gx, b0y_dz_gy, b0y_dz_gz;
    volume<double> b0z_dx_gx, b0z_dx_gy, b0z_dx_gz, b0z_dy_gx, b0z_dy_gy, b0z_dy_gz, b0z_dz_gx, b0z_dz_gy, b0z_dz_gz;
      cout<<"Reading in the 9 B0 volumes..."<<endl;
      read_volume(b0x_dx,opt_b0.value()+"x_dx");
      read_volume(b0x_dy,opt_b0.value()+"x_dy");
      read_volume(b0x_dz,opt_b0.value()+"x_dz");
      read_volume(b0y_dx,opt_b0.value()+"y_dx");
      read_volume(b0y_dy,opt_b0.value()+"y_dy");
      read_volume(b0y_dz,opt_b0.value()+"y_dz");
      read_volume(b0z_dx,opt_b0.value()+"z_dx");
      read_volume(b0z_dy,opt_b0.value()+"z_dy");
      read_volume(b0z_dz,opt_b0.value()+"z_dz");
      cout<<"Calculating the 27 B0_grad volumes..."<<endl;
      calc_gradients(b0x_dx,b0x_dx_gx,b0x_dx_gy,b0x_dx_gz);
      calc_gradients(b0x_dy,b0x_dy_gx,b0x_dy_gy,b0x_dy_gz);
      calc_gradients(b0x_dz,b0x_dz_gx,b0x_dz_gy,b0x_dz_gz);
      calc_gradients(b0y_dx,b0y_dx_gx,b0y_dx_gy,b0y_dx_gz);
      calc_gradients(b0y_dy,b0y_dy_gx,b0y_dy_gy,b0y_dy_gz);
      calc_gradients(b0y_dz,b0y_dz_gx,b0y_dz_gy,b0y_dz_gz);
      calc_gradients(b0z_dx,b0z_dx_gx,b0z_dx_gy,b0z_dx_gz);
      calc_gradients(b0z_dy,b0z_dy_gx,b0z_dy_gy,b0z_dy_gz);
      calc_gradients(b0z_dz,b0z_dz_gx,b0z_dz_gy,b0z_dz_gz);

      save_volume(b0x_dx_gx,opt_b0.value()+"x_dx_gx");
      save_volume(b0x_dx_gy,opt_b0.value()+"x_dx_gy");
      save_volume(b0x_dx_gz,opt_b0.value()+"x_dx_gz");
      save_volume(b0x_dy_gx,opt_b0.value()+"x_dy_gx");
      save_volume(b0x_dy_gy,opt_b0.value()+"x_dy_gy");
      save_volume(b0x_dy_gz,opt_b0.value()+"x_dy_gz");
      save_volume(b0x_dz_gx,opt_b0.value()+"x_dz_gx");
      save_volume(b0x_dz_gy,opt_b0.value()+"x_dz_gy");
      save_volume(b0x_dz_gz,opt_b0.value()+"x_dz_gz");

      save_volume(b0y_dx_gx,opt_b0.value()+"y_dx_gx");
      save_volume(b0y_dx_gy,opt_b0.value()+"y_dx_gy");
      save_volume(b0y_dx_gz,opt_b0.value()+"y_dx_gz");
      save_volume(b0y_dy_gx,opt_b0.value()+"y_dy_gx");
      save_volume(b0y_dy_gy,opt_b0.value()+"y_dy_gy");
      save_volume(b0y_dy_gz,opt_b0.value()+"y_dy_gz");
      save_volume(b0y_dz_gx,opt_b0.value()+"y_dz_gx");
      save_volume(b0y_dz_gy,opt_b0.value()+"y_dz_gy");
      save_volume(b0y_dz_gz,opt_b0.value()+"y_dz_gz");

      save_volume(b0z_dx_gx,opt_b0.value()+"z_dx_gx");
      save_volume(b0z_dx_gy,opt_b0.value()+"z_dx_gy");
      save_volume(b0z_dx_gz,opt_b0.value()+"z_dx_gz");
      save_volume(b0z_dy_gx,opt_b0.value()+"z_dy_gx");
      save_volume(b0z_dy_gy,opt_b0.value()+"z_dy_gy");
      save_volume(b0z_dy_gz,opt_b0.value()+"z_dy_gz");
      save_volume(b0z_dz_gx,opt_b0.value()+"z_dz_gx");
      save_volume(b0z_dz_gy,opt_b0.value()+"z_dz_gy");
      save_volume(b0z_dz_gz,opt_b0.value()+"z_dz_gz");

   
///////////////////////////////////////////////////////////
  //NO MOTION/
  ///////////////////////////////////////////////////////////
  if (level==4){
    volume4D<double> b0z_dz;
    read_volume4D(b0z_dz,opt_b0.value()+"z_dz"); 
    print_volume_info(b0z_dz,"b0z_dz");
    volume4D<double> b0z_dz_gx;
    volume4D<double> b0z_dz_gy;
    volume4D<double> b0z_dz_gz;
    calc_gradients4D(b0z_dz,b0z_dz_gx,b0z_dz_gy,b0z_dz_gz);
    save_volume4D(b0z_dz_gx,opt_b0.value()+"z_dz_gx");
    save_volume4D(b0z_dz_gy,opt_b0.value()+"z_dz_gy");
    save_volume4D(b0z_dz_gz,opt_b0.value()+"z_dz_gz");
    print_volume_info(b0z_dz_gz,"b0z_dz_gz");
 
    volume4D<double> b0x;
    volume4D<double> b0y;
    volume4D<double> b0z;
    calc_gradients4D(b0,b0x,b0y,b0z);
    //save_volume(b0z,"b0ztest");
    print_volume_info(b0z,"b0z");
   
   /////////////////////////
   //OUTPUT
   /////////////////////////
   for (int i=1;i<=nreadp;i++) {
      signal(1,i)=sreal[i-1]*1e06; //1e06 we need to make signal have larger values so that it more resembles the scanner
      signal(2,i)=simag[i-1]*1e06; //1e06 we need to make signal have larger values so that it more resembles the scanner
    }
    write_binary_matrix(signal,opt_signal.value());
    cout<<"Possum finished generating the signal..."<<endl;
  #endif

  return 0;
}


int main (int argc, char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    options.add(verbose);
    options.add(help);
    options.add(opt_kcoord);
    options.add(opt_object);
    options.add(opt_tissue);
    options.add(opt_motion);
    options.add(opt_pulse);
    options.add(opt_b0);
    options.add(opt_b0timecourse4D);
    options.add(opt_RFrec);
    options.add(opt_RFtrans);
    options.add(opt_activation);
    options.add(opt_timecourse);
    options.add(opt_activation4D);
    options.add(opt_timecourse4D);
    options.add(opt_level);
    options.add(opt_nproc);
    options.add(opt_procid);
    options.add(opt_signal);
    options.add(opt_slcprof);
    options.add(opt_mainmatrix);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if there are less than 2 non-optional args
    //   or the help was requested or a compulsory option was not set
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
  compute_volume(argc, argv);
  return 0;
}
