#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "EddyHelperClasses.h"
#include "ECScanClasses.h"
#include "DiffusionGP.h"
#include "b0Predictor.h"
#include "EddyUtils.h"
#include "EddyCommandLineOptions.h"

using namespace EDDY;

ReplacementManager Register(const EddyCommandLineOptions&  clo,     // Input
			    ScanType                       st,      // Input
			    ECScanManager&                 sm,      // Input/Output
			    unsigned int                   niter,   // Input
			    NEWMAT::Matrix&                msshist, // Output
			    NEWMAT::Matrix&                phist);  // Output

DiffStatsVector DetectAndReplaceOutliers(// Input
					 const EddyCommandLineOptions& clo,
					 ScanType                      st,
					 // Input/Output
					 ECScanManager&                sm,
					 ReplacementManager&           rm);

boost::shared_ptr<DWIPredictionMaker> LoadPredictionMaker(// Input
							  const EddyCommandLineOptions& clo,
							  ScanType                      st,
							  const ECScanManager&          sm,
							  // Output
							  NEWIMAGE::volume<float>&      mask);

void Diagnostics(const EddyCommandLineOptions&  clo,      // Input
		 unsigned int                   iter,     // Input
		 ScanType                       st,       // Input
		 const ECScanManager&           sm,       // Input
                 const double                   *mss_tmp, // Input
                 const DiffStatsVector&         stats,    // Input
		 NEWMAT::Matrix&                mss,      // Output
		 NEWMAT::Matrix&                phist);   // Output

int main(int argc, char *argv[])
{
    // Parse comand line input
  EddyCommandLineOptions clo(argc,argv); // Command Line Options

  // Read all available info
  if (clo.Verbose()) cout << "Reading images" << endl;
  ECScanManager sm(clo.ImaFname(),clo.MaskFname(),clo.AcqpFname(),clo.TopupFname(),
                   clo.BVecsFname(),clo.BValsFname(),clo.FirstLevelModel(),
		   clo.Indicies(),clo.SessionIndicies(),clo.NoOfSessions()); // Scan Manager
  if (clo.Verbose() && clo.FWHM()) cout << "Smoothing images" << endl;
  sm.SetFWHM(clo.FWHM());
  if (clo.ResamplingMethod() == LSR) {
    if (!sm.CanDoLSRResampling()) throw EddyException("These data do not support least-squares resampling");
  }

  // Set initial parameters. This option is only for testing/debugging/personal use
  if (clo.InitFname() != std::string("")) {
    if (clo.RegisterDWI() && clo.Registerb0()) sm.SetParameters(clo.InitFname(),ANY);
    else if (clo.RegisterDWI()) sm.SetParameters(clo.InitFname(),DWI);
    else sm.SetParameters(clo.InitFname(),B0);
  }

  // Do the registration
  NEWMAT::Matrix dwi_mss, b0_mss, dwi_ph, b0_ph;
  ReplacementManager *dwi_rm;
  if (clo.NIter() && clo.RegisterDWI()) {
    dwi_rm = new ReplacementManager(Register(clo,DWI,sm,clo.NIter(),dwi_mss,dwi_ph));
    // Write outlier information
    if (clo.ReplaceOutliers()) {
      std::vector<unsigned int> i2i = sm.GetDwi2GlobalIndexMapping();
      dwi_rm->WriteReport(i2i,clo.OLReportFname());
    }
  }
  if (clo.NIter() && clo.Registerb0()) ReplacementManager b0_rm = Register(clo,B0,sm,clo.NIter(),b0_mss,b0_ph);

  // Separate field offset from subject movement in PE direction
  if (clo.RegisterDWI()) sm.SeparateFieldOffsetFromMovement();

  // Set reference for location
  if (clo.RegisterDWI()) sm.SetDWIReference();
  if (clo.Registerb0()) sm.Setb0Reference();

  // Write registration parameters
  if (clo.RegisterDWI() && clo.Registerb0()) sm.WriteParameterFile(clo.ParOutFname());
  else if (clo.RegisterDWI()) sm.WriteParameterFile(clo.ParOutFname(),DWI);
  else sm.WriteParameterFile(clo.ParOutFname(),B0);

  // Write registered images
  if (clo.RegisterDWI() && clo.Registerb0()) sm.WriteRegisteredImages(clo.IOutFname(),clo.ResamplingMethod());
  else if (clo.RegisterDWI()) sm.WriteRegisteredImages(clo.IOutFname(),clo.ResamplingMethod(),DWI);
  else sm.WriteRegisteredImages(clo.IOutFname(),clo.ResamplingMethod(),B0);

  // Write EC fields
  if (clo.WriteFields()) {
    if (clo.RegisterDWI() && clo.Registerb0()) sm.WriteECFields(clo.ECFOutFname());
    else if (clo.RegisterDWI()) sm.WriteECFields(clo.ECFOutFname(),DWI);
    else sm.WriteECFields(clo.ECFOutFname(),B0);
  }

  if (clo.NIter() && clo.History()) { 
    if (clo.RegisterDWI()) {
      MISCMATHS::write_ascii_matrix(clo.DwiMssHistoryFname(),dwi_mss); 
      MISCMATHS::write_ascii_matrix(clo.DwiParHistoryFname(),dwi_ph);
    }
    if (clo.Registerb0()) {
      MISCMATHS::write_ascii_matrix(clo.B0MssHistoryFname(),b0_mss); 
      MISCMATHS::write_ascii_matrix(clo.B0ParHistoryFname(),b0_ph);
    }
  }
  exit(EXIT_SUCCESS);
}

ReplacementManager Register(const EddyCommandLineOptions&  clo,     // Input
			    ScanType                       st,      // Input
			    ECScanManager&                 sm,      // Input/Output
			    unsigned int                   niter,   // Input
			    NEWMAT::Matrix&                msshist, // Output
			    NEWMAT::Matrix&                phist)   // Output
{
  msshist.ReSize(niter,sm.NScans(st));
  phist.ReSize(niter,sm.NScans(st)*sm.Scan(0,st).NParam());
  double *mss_tmp = new double[sm.NScans(st)]; 
  ReplacementManager rm(sm.NScans(st),static_cast<unsigned int>(sm.Scan(0,st).GetIma().zsize()),clo.OLNStdev(),clo.OLNVox());
  NEWIMAGE::volume<float> mask = sm.Scan(0,st).GetIma(); EddyUtils::SetTrilinearInterp(mask); mask = 1.0; // Mask in model space

  for (unsigned int iter=0; iter<niter; iter++) {
    // Detect outliers and replace them
    DiffStatsVector stats(sm.NScans(st));
    if (iter && st==DWI && clo.ReplaceOutliers()) stats = DetectAndReplaceOutliers(clo,st,sm,rm);

    // Load prediction maker in model space
    // printf("Starting to load and evaluate predictionmaker\n");
    // clock_t stime = clock();
    boost::shared_ptr<DWIPredictionMaker> pmp = LoadPredictionMaker(clo,st,sm,mask);
    // printf("It took %f sec to load and evaluate predictionmaker\n",double(stime-clock())/double(CLOCKS_PER_SEC));

    // Calculate the parameter updates
    if (clo.Verbose()) cout << "Calculating parameter updates" << endl;
    // int dbl = clo.DebugLevel();  // Can't check DebugLevel inside parallel section :(
    // bool vv = clo.VeryVerbose(); // Compiler crashes if I test clo.VeryVerbose() inside parallel loop :(
# pragma omp parallel for shared(mss_tmp, pmp)
    for (int s=0; s<int(sm.NScans(st)); s++) {
      // Get prediction in model space 
      NEWIMAGE::volume<float> pred = pmp->Predict(s);
      // Update parameters
      if (clo.DebugLevel()) {
	mss_tmp[s] = EddyUtils::param_update_debug(pred,sm.GetSuscHzOffResField(),mask,ALL,true,s,iter,clo.DebugLevel(),sm.Scan(s,st),NULL);
      }
      else mss_tmp[s] = EddyUtils::MovAndECParamUpdate(pred,sm.GetSuscHzOffResField(),mask,true,sm.Scan(s,st));
      if (clo.VeryVerbose()) printf("Iter: %d, scan: %d, mss = %f\n",iter,s,mss_tmp[s]);
    }

    // Print/collect some information that can be used for diagnostics
    Diagnostics(clo,iter,st,sm,mss_tmp,stats,msshist,phist);

    // Maybe use model based EC parameters
    // if () sm.SetPredictedECParam();
  }

  delete [] mss_tmp;
  return(rm);
}

boost::shared_ptr<DWIPredictionMaker> LoadPredictionMaker(// Input
							  const EddyCommandLineOptions& clo,
							  ScanType                      st,
							  const ECScanManager&          sm,
							  // Output
							  NEWIMAGE::volume<float>&      mask)
{
  boost::shared_ptr<DWIPredictionMaker>  pmp;                                 // Prediction Maker Pointer
  if (st==DWI) pmp = boost::shared_ptr<DWIPredictionMaker>(new DiffusionGP);  // Gaussian Process
  else pmp = boost::shared_ptr<DWIPredictionMaker>(new b0Predictor);          // Silly mean predictor
  pmp->SetNoOfScans(sm.NScans(st));
  mask = sm.Scan(0,st).GetIma(); EddyUtils::SetTrilinearInterp(mask); mask = 1.0;

  if (clo.Verbose()) cout << "Loading prediction maker";
  // bool vv = clo.VeryVerbose(); // Compiler crashes if I test clo.VeryVerbose() inside parallel loop :(
  if (clo.VeryVerbose()) cout << endl << "Scan: ";
#pragma omp parallel for shared(pmp,st)
  for (int s=0; s<int(sm.NScans(st)); s++) {
    if (clo.VeryVerbose()) printf(" %d",s);
    NEWIMAGE::volume<float> tmpmask = sm.Scan(s,st).GetIma(); 
    EddyUtils::SetTrilinearInterp(tmpmask); tmpmask = 1.0;
    pmp->SetScan(sm.GetUnwarpedScan(s,tmpmask,st),sm.Scan(s,st).GetDiffPara(),s);
#pragma omp critical
    {
      mask *= tmpmask;
    }
  }
  if (clo.Verbose()) cout << endl << "Evaluating prediction maker model" << endl;
  pmp->EvaluateModel(sm.Mask()*mask);

  return(pmp);
}

void Diagnostics(const EddyCommandLineOptions&  clo,      // Input
		 unsigned int                   iter,     // Input
		 ScanType                       st,       // Input
		 const ECScanManager&           sm,       // Input
                 const double                   *mss_tmp, // Input
                 const DiffStatsVector&         stats,    // Input
		 NEWMAT::Matrix&                mss,      // Output
		 NEWMAT::Matrix&                phist)    // Output
{
  if (clo.Verbose()) {
    double tss=0.0;
    for (unsigned int s=0; s<sm.NScans(st); s++) tss+=mss_tmp[s]; 
    cout << "Iter: " << iter << ", Total mss = " << tss/sm.NScans(st) << endl;
  }
  if (clo.History()) {
    for (unsigned int s=0; s<sm.NScans(st); s++) {
      mss(iter+1,s+1) = mss_tmp[s];
      phist.SubMatrix(iter+1,iter+1,s*sm.Scan(0,st).NParam()+1,(s+1)*sm.Scan(0,st).NParam()) = sm.Scan(s,st).GetParams().t();
    }
  }
  if (clo.WriteSliceStats()) {
    char istring[256];
    sprintf(istring,"EddySliceStatsIteration%02d",iter);
    stats.Write(string(istring));
  }
}

DiffStatsVector DetectAndReplaceOutliers(// Input
					 const EddyCommandLineOptions& clo,
					 ScanType                      st,
					 // Input/Output
					 ECScanManager&                sm,
					 ReplacementManager&           rm)
{
  // printf("Entering DetectAndReplaceOutliers\n");
  // Load up a prediction maker with unwarped and unsmoothed data
  DWIPredictionMaker  *pmp;            // Prediction maker used for the outlier detection
  if (st==DWI) pmp = new DiffusionGP;  // Gaussian Process
  else pmp = new b0Predictor;          // Silly mean predictor
  pmp->SetNoOfScans(sm.NScans(st));
  NEWIMAGE::volume<float> mask(sm.Scan(0,st).GetIma().xsize(),sm.Scan(0,st).GetIma().ysize(),sm.Scan(0,st).GetIma().zsize());
  NEWIMAGE::copybasicproperties(sm.Scan(0,st).GetIma(),mask); mask = 1.0;
  // printf("Starting to load predictionmaker\n");
  // clock_t stime = clock();
#pragma omp parallel for shared(pmp,st)
  for (int s=0; s<int(sm.NScans(st)); s++) {
    NEWIMAGE::volume<float> tmpmask(sm.Scan(s,st).GetIma().xsize(),sm.Scan(s,st).GetIma().ysize(),sm.Scan(s,st).GetIma().zsize());
    NEWIMAGE::copybasicproperties(sm.Scan(s,st).GetIma(),tmpmask); tmpmask = 1.0;
    pmp->SetScan(sm.GetUnwarpedOrigScan(s,tmpmask,st),sm.Scan(s,st).GetDiffPara(),s);
#pragma omp critical
    { mask *= tmpmask; }
  }
  // printf("It took %f sec to load pm\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  // printf("Starting to evaluate the pm model\n");
  // stime = clock();
  pmp->EvaluateModel(sm.Mask()*mask);
  // printf("It took %f sec to evaluate pm\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  // Use these to generate slice-wise stats on difference between observation and prediction
  // printf("Starting to generate slice-wise stats\n");
  // stime = clock();
  DiffStatsVector stats(sm.NScans(st));
#pragma omp parallel for shared(pmp,mask,stats,st)
  for (int s=0; s<int(sm.NScans(st)); s++) {
    NEWIMAGE::volume<float> pred = pmp->Predict(s);
    stats[s] = EddyUtils::GetSliceWiseStats(pred,sm.GetSuscHzOffResField(),mask,sm.Mask(),sm.Scan(s,st));
  }
  // printf("It took %f sec to generate stats\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  // stats.Write("slice_wise_stats");
  // Detect outliers and update replacement manager
  // printf("Updating replacement manager\n");
  // stime = clock();
  rm.Update(stats);
  // printf("It took %f sec to update\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  // Replace outlier slices with their predictions
  // printf("Replacing outliers\n");
  // stime = clock();
#pragma omp parallel for shared(pmp,mask,st)
  for (int s=0; s<int(sm.NScans(st)); s++) {
    // printf("Checking scan %d for outliers\n",s);
    std::vector<unsigned int> ol = rm.OutliersInScan(s);
    if (ol.size()) { // If this scan has outlier slices
      // printf("Scan %d has %d outlier slices\n",s,ol.size());
      NEWIMAGE::volume<float> pred = pmp->Predict(s);
      sm.Scan(s,st).ReplaceSlices(pred,sm.GetSuscHzOffResField(),mask,ol);
    }
  }
  // printf("It took %f sec to replace outliers\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  //  rm.DumpOutlierMap("DebugReplacementManager");
  //  exit(EXIT_SUCCESS);
  delete pmp;
  return(stats);
}
