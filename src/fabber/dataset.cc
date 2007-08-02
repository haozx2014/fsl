/*  dataset.cc - Data-loading class for FABBER

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

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

#include "easylog.h"
#include "dataset.h"

void DumpVolumeInfo(const VolumeInfo& info, string indent = "", ostream& out = LOG)
{
  Tracer_Plus tr("DumpVolumeInfo");
  LOG << indent << "Dimensions: x=" << info.x << ", y=" << info.y 
      << ", z=" << info.z << ", vols=" << info.v << endl;
  LOG << indent << "Voxel size: x=" << info.vx << "mm, y=" << info.vy 
      << "mm, z=" << info.vz << "mm, TR=" << info.tr << " sec\n";
  LOG << indent << "Intents: " << info.intent_code << ", " << info.intent_p1
      << ", " << info.intent_p1 << ", " << info.intent_p2 << ", " 
      << info.intent_p3 << endl;
}


void DataSet::LoadData(ArgsType& args)
{
  Tracer_Plus tr("LoadData");
  string dataOrder = args.ReadWithDefault("data-order","interleave");

  if (dataOrder == "singlefile")
    {
      LOG << "  Loading all data from a single file..." << endl;
      string dataFile = args.Read("data");
      string maskFile = args.Read("mask");

      // Load the files.  Note: these functions don't throw errors if file doesn't exist --
      // they just crash.  Hence the detailed logging before we try anything.

      LOG_ERR("    Loading data from '" << dataFile << "'" << endl);
      VolumeSeries data; data.read(dataFile);
      DumpVolumeInfo(data.getInfo(), "      ");

      LOG_ERR("    Loading mask data from '" + maskFile << "'" << endl);
      mask.read(maskFile);
      DumpVolumeInfo(mask.getInfo(), "      ");

      LOG << "    Applying mask to data..." << endl;
      mask.threshold(1e-16);
      // threshold using mask:
      try {
        data.setPreThresholdPositions(mask.getPreThresholdPositions());
        data.thresholdSeries();
      } catch (Exception) {
        LOG_ERR("\n*** NEWMAT error while thresholding time-series... "
                << "Most likely a dimension mismatch. ***\n");
        DumpVolumeInfo(data.getInfo());
        LOG_ERR("Mask:\n");
        DumpVolumeInfo(mask.getInfo());
        LOG_ERR("\nThis is fatal... rethrowing exception.\n");
        throw;
      }

      voxelData = data;
    }
  else if (dataOrder == "interleave" || dataOrder == "concatenate")
    {
      LOG << "  Loading data from multiple files..." << endl;
      
      vector<VolumeSeries> dataSets;
      int nTimes = -1;
      while (true)
	{
	  int N = dataSets.size() + 1;
	  string datafile = args.ReadWithDefault("data"+stringify(N), "stop!");
	  if (datafile == "stop!") break;

	  // Load the files.  Note: these functions don't throw errors if file doesn't exist --
	  // they just crash.  Hence the detailed logging before we try anything.
	  LOG_ERR("    Loading " << "data"+stringify(N) << " from '" << datafile << "'" << endl);
	  dataSets.push_back(VolumeSeries());
	  dataSets.back().read(datafile);
	  DumpVolumeInfo(dataSets.back().getInfo(), "      ");

	  if (nTimes == -1) 
	    nTimes = dataSets[0].Nrows();
	  else if (nTimes != dataSets.back().Nrows())
	    throw Invalid_option("Data sets must all have the same number of time points");
	}

      int nSets = dataSets.size();
      if (nSets < 1)
	throw Invalid_option("At least one data file is required: --date1=<file1> [--data2=<file2> [...]]");      

      string maskFile = args.Read("mask");
 
      LOG_ERR("    Loading mask data from '" + maskFile << "'" << endl);
      mask.read(maskFile);
      DumpVolumeInfo(mask.getInfo(), "      ");

      LOG << "    Applying mask to all data sets..." << endl;
      mask.threshold(1e-16);
      // threshold using mask:
      for (int i = 0; i < nSets; i++)
	{
	  try {
	    dataSets.at(i).setPreThresholdPositions(mask.getPreThresholdPositions());
	    dataSets.at(i).thresholdSeries();
	  } catch (Exception) {
	    LOG_ERR("\n*** NEWMAT error while thresholding time-series... "
		    << "Most likely a dimension mismatch (more details in logfile) ***\n");
	    LOG << "Data set " << i+1 << ":\n"; 
	    DumpVolumeInfo(dataSets.at(i).getInfo());
	    LOG << "Mask:\n"; 
	    DumpVolumeInfo(mask.getInfo());
	    LOG_ERR("\nThis is fatal... rethrowing exception.\n");
	    throw;
	  }
	}
      
      if (dataOrder == "interleave")
        {
          LOG << "    Combining data into one big matrix by interleaving..." << endl;
          // Interleave:
          voxelData.ReSize(nTimes * nSets, dataSets[0].Ncols());
          for (int i = 0; i < nTimes; i++)
            {
              for (int j = 0; j < nSets; j++)
	        {
	          voxelData.Row(nSets*i+j+1) = dataSets.at(j).Row(i+1);
	        }
	    }
        }
      else
        {
          LOG << "    Combining data into one big matrix by concatenating..." << endl;
          // Concatenate:
          voxelData = dataSets.at(0);
          for (unsigned j = 1; j < dataSets.size(); j++)
            voxelData &= dataSets.at(j);
        }
      
      LOG << "    Done loading data, size = " 
	  << voxelData.Nrows() << " timepoints by "
	  << voxelData.Ncols() << " voxels" << endl;
    }
  else
    throw Invalid_option(("Unrecognized --dataorder: " + dataOrder + " (try interleave or singlefile)").c_str());
}

