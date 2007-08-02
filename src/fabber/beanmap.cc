/* BEANMAP main program
 *   Bayes done Efficiently by Approximating Nonlinear Models of ASL Perfusion
 *
 * Implements the Variational Bayes method for approximate inference on
 *  nonlinear models.  Currently, we're applying this to analysis of dual-echo
 *  ASL data, but it would probably be suitable for any type of scan where we 
 *  have a nonlinear forward (predictive) model.  The noise models are also
 *  quite restricted at the moment: AR(1) and a modification of this to two (or
 *  more) correlated timeseries.
 *
 * FMRIB Centre, University of Oxford.  Copyright Adrian Groves, 2007.
 *
 * Last modified:
 * 	$Date: 2007/07/24 09:35:37 $
 * 	$Author: adriang $
 * 	$Revision: 1.31 $
 */

#include <iostream>
#include <exception>
#include <stdexcept>
#include <map>
#include <string>
#include "inference.h"
#include "miscmaths/volumeseries.h"
#include "miscmaths/volume.h"
using namespace std;
using namespace MISCMATHS;

//#include "utils/log.h"
#include "easylog.h"
using namespace Utilities;

/*** Function declarations ***/

void Usage(const string& errorString = "");
ostream& operator<<(ostream& out, const VolumeInfo& info);

/*** Function implementations ***/

int main(int argc, char** argv)
{
  try
    {
      cout << "------------------\n";
      cout << "Welcome to BEANMAP" << endl;

//      map<string,string> args;
//      ParseArguments(argc, argv, args);
      EasyOptions args(argc, argv);

      if (args.ReadBool("help")) 
        { 
            string model = args.ReadWithDefault("model","");
            if (model == "")
                Usage();
            else
                FwdModel::ModelUsageFromName(model, args);
                                 
            return 0; 
        }
      EasyLog::StartLog(args.Read("output", 
        "Must specify an output directory, for example: --output=mytestrun"));
        
      LOG_ERR("Logfile started: " << EasyLog::GetOutputDirectory() 
	          << "/logfile" << endl);

      time_t startTime;
      time(&startTime);
      LOG_ERR("Start time: " << ctime(&startTime));

      // Diagnostic information: software versions
      // This only versions this file... should really use all.
      LOG_ERR("BEANMAP revision: $Id: beanmap.cc,v 1.31 2007/07/24 09:35:37 adriang Exp $\n");
#if 0
      // Slightly more detailed version information
      int ret = system("cvs status | grep Status | grep -v Up-to-date");
      if (ret == 0) LOG_ERR("However, some files appear to have changed.\n");
      else if (ret == 1) LOG_ERR("And all files appear to be up-to-date.\n");
      else LOG_ERR(
        "Return status of 'cvs status | grep Status | grep -v Up-to-date' == " 
        << ret << endl);
#endif
      LOG << "Command line and effective options:\n" << args.Read("") << endl;
      LOG << "--output='" << EasyLog::GetOutputDirectory() << "'" << endl;
      LOG << args << "--------------------" << endl;

      // Start timing/tracing if requested
      bool recordTimings = false;
  
      if (args.ReadBool("debug-timings")) 
        { recordTimings = true; Tracer_Plus::settimingon(); }
      if (args.ReadBool("debug-instant-stack")) 
        { Tracer_Plus::setinstantstackon(); } // instant stack isn't used?
      if (args.ReadBool("debug-running-stack")) 
        { Tracer_Plus::setrunningstackon(); }

      Tracer_Plus tr("BEANMAP main (outer)");
      // can't start it before this or it segfaults if an exception is thown with --debug-timings on.

      // Start a new tracer for timing purposes
      { Tracer_Plus tr2("BEANMAP main()");

      InferenceTechnique* infer = 
        InferenceTechnique::NewFromName(args.Read("method"));

      infer->Setup(args);
      infer->SetOutputFilenames(EasyLog::GetOutputDirectory());
      
      DataSet allData;
      //      Volume mask;
      //      LoadData(args, allData, mask); // we'll need the mask later
      allData.LoadData(args);

      // Arguments should all have been used by now, so complain if there's anything left.
      args.CheckEmpty();   
      
      // Calculations
      infer->DoCalculations(allData);
      infer->SaveResults(allData);
      delete infer;
      
      LOG_ERR("BEANMAP is all done." << endl);

      time_t endTime;
      time(&endTime);
      LOG << "Start time: " << ctime(&startTime);   // Bizarrely, ctime() ends with a \n.
      LOG << "End time: " << ctime(&endTime);
      LOG_ERR("Duration: " << endTime-startTime << " seconds." << endl);

      } // End of timings
     
      if (recordTimings) {
        tr.dump_times(EasyLog::GetOutputDirectory());
        LOG_ERR("Timing profile information recorded to " 
		<< EasyLog::GetOutputDirectory() << "/timings.html" << endl);
      }

      cout << "Logfile was: " << EasyLog::GetOutputDirectory() << "/logfile" << endl;
      EasyLog::StopLog();

      return 0;
    }
  catch (const Invalid_option& e)
    {
      LOG_ERR_SAFE("Invalid_option exception caught in beanmap:\n  " << Exception::what() << endl);
      Usage(Exception::what());
    }
  catch (const exception& e)
    {
      LOG_ERR_SAFE("STL exception caught in beanmap:\n  " << e.what() << endl);
      //Usage(e.what());
    }
  catch (Exception)
    {
      LOG_ERR_SAFE("NEWMAT exception caught in beanmap:\n  " 
	      << Exception::what() << endl);
      //Usage(Exception::what());
    }
  catch (...)
    {
      LOG_ERR_SAFE("Some other exception caught in beanmap!" << endl);
      //Usage("I hate to say this but... an unknown error has occurred.");
    }
  
  if (EasyLog::LogStarted())
    {
      cout << "Logfile was: " << EasyLog::GetOutputDirectory() << "/logfile" << endl;
      EasyLog::StopLog();
    }

  return 1;
}

void Usage(const string& errorString)
{
    cout << "\n\nUsage: beanmap <arguments>\n"
     << "Where most arguments are mandatory except those in [brackets].\n"
     << "Use -@ argfile to read additional arguments from a text file.\n\n";

    cout << "[--help] : drop everything and just print this usage message\n"
     << "--output=/path/to/output : put output here (including logfile\n"
     << "--method={vb|spatialvb} : use VB (or VB with spatial priors)\n"
     << "--max-iterations=NN : maximum (actually, exact) number of iterations\n"
//     << "      [--convergence={maxits|pointzeroone}
     << "--data-order={interleave|concatenate} : data should be interleaved (TE1/TE2) or in-order?\n"
     << "  use --data1=file1, --data2=file2, etc.\n"
     << "--mask=maskfile : inference will only be performed where mask>0\n"
     << "--model={quipss2|simple|grase|eagle} : forward model to use\n"
     << "  Parameters for QUIPSS II model (default values shown):\n"
//     << "[--scan-params=cmdline] : where to obtain scan parameters from\n"
     << "  [--ti1=0.6] and [--ti2=1.5]: TI1 and TI2 of QUIPSS II sequence (s)\n"
     << "  [--t1b=1.66] and [--inv-eff=1] : Assumed T1 of blood and inversion efficiency\n"
     << "  [--tag-pattern=TC] : Repeating pattern of tag & control volumes\n"
     << "  [--te1=9.1] and [--te2=30] : Echo times (ms)\n"
     << "  --bold-basis=vestfile : Basis functions for BOLD effect (omit baseline)\n"
     << "  --cbf-basis=vestfile : Basis functions for CBF (omit baseline regressor)\n"
     << "  --cbv-basis=vestfile : Basis for CBV (omit baseline regressor)\n"
     << "  [--nuisance-basis=vestfile] : Basis for nuisance regressors (defaults to none)\n"
     << "--noise=ar1 : Noise model to use (two cross-linked AR1 models)\n"
     << "--ar1-cross-terms={dual|none|same} : add correlation between TE1 & TE2 to AR noise model"
      //     << "[--noise-prior=hardcoded] : Uses uninformative noise model priors\n"
      //     << "[--fwd-initial-prior] and [--fwd-initial-posterior] : Fwd-model prior and starting point\n"
     << "[--save-model-fit] and [--save-residuals] : Save model fit/residuals files\n"
     << "[--print-free-energy] : Calculate & dump F to the logfile after each update\n"
     << "[--allow-bad-voxels] : Skip to next voxel if a numerical exception occurs (don't stop)\n"
     << endl;


    if (errorString.length() > 0)
        cout << "\nImmediate cause of error: " << errorString << endl;
}
