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

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"


using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace std;



class ZMRISegmentation
{
 public:
   ZMRISegmentation();
  ~ZMRISegmentation()
    {
    }
  bool Segments(bool pve,const NEWIMAGE::volume<float>& pcsf, const NEWIMAGE::volume<float>& pgm, const NEWIMAGE::volume<float>& pwm);
  void TanakaCreate(const NEWIMAGE::volume<float>& image, float fbeta, int nclasses, float nblowpass, bool bbias, int biterationspve, float mixeltypeMRF, int nbiter, int initinitfixed, int winitfixed, int bapused, float hyp, bool verb,string mansegfle,int typeoffile);
  int TanakaMain(NEWIMAGE::volume<float>& pcsf, NEWIMAGE::volume<float>& pgm, NEWIMAGE::volume<float>& pwm);


  volume4D<float> members;
  volume4D<float> m_post;
  volume4D<float> m_pve;
  volume<float> m_Finalbias;
  volume<int> m_Segment;
  volume<int> m_pveSegment;
  volume<int>hardPV;


private:
  NEWIMAGE::volume4D<float> InitclassAlt(int numberofsegs);
  NEWIMAGE::volume4D<float> Initclass(int numberofsegs);
  NEWIMAGE::volume4D<float> pvprobs(int clas);
  NEWIMAGE::volume<float> Convolve(NEWIMAGE::volume<float>& resfieldimage);
  float pvmeans(int clas);
  float pvvar(int clas);
  float ICMCopy(float BB, bool hyperparameter);
  float numberofvoxels();
  float hyper();
  float MeanFieldBetaHyper();
  float PVMRF(int x, int y, int z);
  float HyperparameterEst(float beta, int algorithm);
  float pveConditional(int x, int y, int z);
  float MRFWeightsTotal();
  double MRFWeightsInner(const int x, const int y, const int z,const int c);
  double MRFWeightsAM(const int l, const int m, const int n);
  float pvePosterior(int x, int y, int z);
  float pveConditionalInit(int x, int y, int z);
  float pvePosteriorInit(int x, int y, int z);
  float LoopyMessagePassing(float BB, bool hyperparameter);
  float Loopy2DMessagePassing(float BB, bool hyperparameter);
  float SimplePrior(float BB, bool hyperparameter);
  float logGaussian(float y_i, float mu, float sigma);
  float PVEnergy(int x, int y, int z, float mu, float sigma);
  int qsort();
  int ICM(float BB, bool hyperparameter);
  void TanakaIterations();
  void TanakaHyper();
  void TanakaPriorHyper();
  void Initialise();
  void Initialise(const NEWIMAGE::volume<float>& pcsf, const NEWIMAGE::volume<float>& pgm, const NEWIMAGE::volume<float>& pwm);
  void InitSimple(const NEWIMAGE::volume<float>& pcsf, const NEWIMAGE::volume<float>& pgm, const NEWIMAGE::volume<float>& pwm);
  void Apriorimap(const NEWIMAGE::volume<float>& pcsf, const NEWIMAGE::volume<float>& pgm, const NEWIMAGE::volume<float>& pwm);
  void EMloop();
  void pvrecompmeanvar();
  void onepvefinal();
  void UpdateMembers(NEWIMAGE::volume4D<float>& probability);
  void UpdateWeights();
  void TreeKMeans();
  void WeightedKMeans();
  void pveIterations();
  void pveInitialize();
  void InitWeights();
  void takeexpo();
  void InitKernel(float kernalsize);
  void Classification(int x, int y, int z);
  void Classification();
  void pveClassification(int x, int y, int z);
  void pveClassification();
  void Probabilities(int index);
  void recomputeMeansVariances(int index);
  void maximumVariances(int index);
  void Volumesquant(const NEWIMAGE::volume4D<float>& probs);
  void PVClassificationStep();
  void ICMPV();
  void PVEnergyInit();
  void PVMRFestimation();
  void PVMoreSophisticated();
  void PVestimation();
  void PVMeansVar(const NEWIMAGE::volume4D<float>& probability);
  void MeansVariances(int numberofclasses, NEWIMAGE::volume4D<float>& probability);
  void BiasRemoval();
  void cliques();
  void MeanFieldEstimation();
  void MeanFieldZwiggle();
  void MeanFieldProb();
  void Dimensions();
  void pvbias();
  


  volume4D<float>PVprob;
  volume4D<float> pvprobsinit;
  volume4D<float> talpriors;
  volume4D<float> m_cliquepot;
  volume4D<float> zwiggle;
  volume4D<float> m_postwiggle;
  volume4D<float> m_postwigglepart;
  volume4D<float> m_loopy;
  volume4D<float> m_prob;
  volume4D<float> pvebias;
  volume<double> m_Finalbiastan;
  volume <float>  m_Mricopy;
  volume<float> pveeng;
  volume<float> m_Mri;
  volume<float> pve_eng;
  volume<float> p_bias;
  volume<float> m_brainmask;
  volume<float> p_meaninvcov;
  volume<float> p_resmean;
  volume<int> m_clique;
  volume<int> m_mask;
  volume<int> m_loopyseg;
  ColumnVector kernelx, kernely, kernelz;
  ColumnVector kernelxtan, kernelytan, kernelztan;
  double* globtot;
  double* volumequant;
  float* m_mean;
  float* m_variance;
  float* m_prior;
  float* npve;
  float* weight;
  float* rhs;
  double m_nxdim, m_nydim, m_nzdim;
  float maximum, minimum;
  float m_nbLowpass;
  float TB;
  float maxvariance;
  float beta, pveB;
  float smallnumber;
  float logposterior;
  float nvoxel;
  float pveBmixeltype;
  float Hyper;
  float amx, amy, amz, amxy, amzy, amzx;
  int nomixeltypeiterations;
  int imagetype;
  int neighbour[27];
  int offset3d;
  int m_nSlicesize, m_nWidth, m_nHeight, m_nDepth, m_nbIter;
  int maxvarianceindex;
  int bapusedflag;
  int noclasses;
  int algorithm;
  int inititerations;
  int usingmeanfield;
  int iterationspve;
  int   biassteps,  initfixed;
  bool biasfieldremoval, autoparaflag, weightschanged;
  bool m_b3D;
  bool verboseusage;
  bool loopyflag;
  string mansegfile;
class kmeansexception: public exception
{
 public:
  virtual const char* what() const throw ()
  {
    return "Exception: Not enough classes detected to init KMeans";
  }
} kmeansexc;

};
