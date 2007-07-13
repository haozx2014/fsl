#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"


using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace std;

class ZMRIMULTISegmentation
{
 public:

  ZMRIMULTISegmentation();
  ~ZMRIMULTISegmentation()
    {
    }
  void TanakaCreate(const NEWIMAGE::volume<float>* images, int nclasses, bool b2d, int nbiter,float nblowpass, float fbeta, int bapused, bool pveboolean, int nchan, int altype, bool bbias, int initfixity, bool verb,int biterationspve, int winit, float mixelpB, float Hyp,string mansegfle,int finalT2file);
  int TanakaMain(NEWIMAGE::volume<float>& pcsf, NEWIMAGE::volume<float>& pgm, NEWIMAGE::volume<float>& pwm);

 
  volume4D<float> m_post;
  volume4D<float> members;
  volume4D<float> m_pve;
  volume<float>* m_Finalbias;
  volume<int> m_Segment;
  volume<int> m_pveSegment;
  volume<int> hardPV;
  Matrix coord_minus_mean;
  float maximum, minimum;


private:
  NEWIMAGE::volume<float> Convolve(NEWIMAGE::volume<float>& resfieldimage);
  NEWIMAGE::volume4D<float> InitclassAlt(int);
  NEWIMAGE::volume4D<float>  Initclass(int noclasses);
  Matrix covariancematrix(int classid, volume4D<float> probability);
  float M_2PI(int numberofchan);
  float logpveGaussian(int x, int y, int z, Matrix mu, Matrix sig, float detsig);
  float PVEnergy(int x, int y, int z, Matrix mu, Matrix sigma, float sigdet);
  float logGaussian(int classnumber, int x, int y, int z);
  float pvmeans(int clas);
  float pvvar(int clas);
  void Volumesquant(const NEWIMAGE::volume4D<float>& probs);
  void Initialise();
  void takeexpo();
  float  MRFWeightsTotal();
  double MRFWeightsInner(const int x, const int y, const int z,const int c);
  double MRFWeightsAM(const int l, const int m, const int n);
  void InitWeights();
  void UpdateWeights();
  void InitSimple(const NEWIMAGE::volume<float>& pcsf, const NEWIMAGE::volume<float>& pgm, const NEWIMAGE::volume<float>& pwm);
  void InitAprioriKMeans();
  void BiasRemoval();
  void MeansVariances(int numberofclasses);
  void Dimensions();
  void WeightedKMeans();
  void mresmatrix();
  void UpdateMembers(NEWIMAGE::volume4D<float>& probability);
  void PVClassificationStep();
  void ICMPV();
  void PVMoreSophisticated();
  void PVestimation();
  void PVMeansVar(const NEWIMAGE::volume4D<float>& m_posttemp);
  void Classification(int x, int y, int z);
  void Classification();
  void InitKernel();
  void pveClassification(int x, int y, int z);
  void pveClassification();
  int qsort();
  void TanakaHyper();
  void TanakaPriorHyper();
  void TanakaIterations();
  Matrix*  m_inv_co_variance;
  Matrix* m_co_variance;
  Matrix m_mean;
  volume4D<float> m_prob;
  volume4D<float>PVprob;
  volume4D<float> talpriors;
  volume<float>* m_Mri;
  volume<float>* m_Mricopy;
  volume<float>* p_bias;
  volume<float>* m_meaninvcov;
  volume<float>* m_resmean;
  volume<float> pve_eng;
  volume<float> m_maskc;
  volume<int> m_mask;
  ColumnVector kernelx, kernely, kernelz;
  float amx, amy, amz, amxy, amzx, amzy;
  double* volumequant;
  float* rhs;
  float* globtot;
  float* m_prior;
  float* npve;
  float* weight;
  int imagetype;
  float m_nbLowpass;
  float MSQR2PI;
  float m_nxdim, m_nydim, m_nzdim;
  float beta, pveB, Hyper;
  float pveBmixeltype, mixeltypeMRF, mixelmrf;
  int noclasses;
  int numberofchannels;
  int m_nSlicesize, m_nWidth, m_nHeight, m_nDepth, m_nbIter, initfixed, inititerations, iterationspve;
  bool weightschanged, bapusedflag;
  bool verboseusage;
  bool biasfieldremoval;
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

