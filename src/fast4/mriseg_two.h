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
