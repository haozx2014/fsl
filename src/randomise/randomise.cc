/*  randomise.cc
    Tim Behrens & Steve Smith & Matthew Webster (FMRIB) & Tom Nichols (UMich)
    Copyright (C) 2004-2008 University of Oxford  */
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
#define CLUST_CON 26

// 26 for FSL 18 for SPM

#include "newimage/newimageall.h"
#include "libprob.h"
#include "ranopts.h"
#include <algorithm>

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace RANDOMISE;

vector<Matrix> voxelwise_evs;
vector<int> voxelwise_ev_numbers;

class Permuter
{ 
public:
  bool isFlipping;
  bool isRandom;
  int nGroups;
  int nSubjects;
  double finalPermutation;
  vector<double>       uniquePermutations; //0 is unique for whole design, 1..nGroups is unique per block
  vector<ColumnVector> permutedLabels;
  vector<ColumnVector> originalLocations;
  ColumnVector truePermutation;
  ColumnVector unpermutedVector;
  Permuter();
  ~Permuter();
  void createPermutationGroups(const Matrix& design, const Matrix& groups, const bool oneNonZeroContrast, const long requiredPermutations);
  void initialisePermutationGroups(const ColumnVector& labels, const long requiredPermutations);
  ColumnVector createDesignLabels(const Matrix& design);
  void createTruePermutation(const ColumnVector& labels, ColumnVector copyOldlabels, ColumnVector& permvec);
  ColumnVector nextPermutation(const long perm);
  ColumnVector nextPermutation(const long permutationNumber,vector<ColumnVector>& previousPermutations, const bool printStatus, const bool isStoring);
  bool isPreviousPermutation(const vector<ColumnVector>& previousPermutations,const ColumnVector& newPermutation);
  ColumnVector permutationVector();
  double reportRequiredPermutations(const bool printToScreen);
  ColumnVector returnPreviousTruePermutation(const long permutationNumber,const vector<ColumnVector>& previousPermutations);
private:
  double computeUniquePermutations(const ColumnVector& labels, const bool calculateFlips);
  void nextShuffle(ColumnVector& perm);
  void nextFlip(ColumnVector& mult);
};

class ParametricStatistic
{
public:
  Matrix originalStatistic,uncorrectedStatistic,maximumDistribution,sumStatMat,sumSampMat;
  bool isAveraging,storingUncorrected;
  void   store(const volume<int>& clusterLabels, const ColumnVector& clusterSizes, const volume<float>& mask, const int contrastNo, const unsigned long permNo);
  void   store(const Matrix& parametricMatrix, const unsigned long permNo);
  void   setup(const int nContrasts,const unsigned long nPerms, const int nVoxels, const bool wantAverage, const bool wantUncorrected);
  void   average(const string filename, const float percentileThreshold,const volume<float>& mask);
  ParametricStatistic() { isAveraging=false; }
};

void ParametricStatistic::setup(const int nContrasts,const unsigned long nPerms, const int nVoxels, const bool wantAverage,const bool wantUncorrected=false)
{
  isAveraging=wantAverage;
  storingUncorrected=wantUncorrected;
  maximumDistribution.ReSize(nContrasts,nPerms);
  maximumDistribution=0;
  if ( storingUncorrected ) {
    uncorrectedStatistic.ReSize(1,nVoxels);
    uncorrectedStatistic=0;
  }
  originalStatistic.ReSize(nContrasts,nVoxels);
  originalStatistic=0;
  if ( isAveraging ) {
    sumStatMat=originalStatistic;
    sumSampMat=originalStatistic;
  }
}

void ParametricStatistic::store(const volume<int>& clusterLabels, const ColumnVector& clusterSizes ,const volume<float>& mask, const int contrastNo, const unsigned long permNo)
{
  if ( clusterSizes.Nrows() > 0 ) 
    maximumDistribution(contrastNo,permNo)=clusterSizes.MaximumAbsoluteValue();
  if (permNo==1 || isAveraging) { 
    volume4D<float> parametricImage(mask.xsize(),mask.ysize(),mask.zsize(),1);
    parametricImage=0;
    for(int z=0; z<mask.zsize(); z++)
      for(int y=0; y<mask.ysize(); y++)
	for(int x=0; x<mask.xsize(); x++)
	  if( clusterLabels(x,y,z) ) 
	    parametricImage(x,y,z,0)=clusterSizes(clusterLabels(x,y,z));	  
    if (permNo==1) 
      originalStatistic.Row(contrastNo)=parametricImage.matrix(mask);
    if (isAveraging) {
      sumStatMat.Row(contrastNo)+=parametricImage.matrix(mask);
      sumSampMat.Row(contrastNo)+=SD(parametricImage.matrix(mask),parametricImage.matrix(mask));
    }
  }
}

void ParametricStatistic::store(const Matrix& parametricMatrix, const unsigned long permNo)
{
  maximumDistribution.Column(permNo)=max(parametricMatrix.t()).t();
  if (permNo==1) 
    originalStatistic=parametricMatrix;
  if (storingUncorrected) 
    uncorrectedStatistic += gt(originalStatistic,parametricMatrix);
  if (isAveraging) {
    sumStatMat+=parametricMatrix;
    sumSampMat+=SD(parametricMatrix,parametricMatrix);
  }
  
}

void ParametricStatistic::average(const string filename, const float percentileThreshold,const volume<float>& mask)
{
  if (isAveraging) {
    volume4D<float> temp;
    temp.setmatrix(sumStatMat,mask);
    save_volume4D(temp,filename+"sum");
    temp.setmatrix(sumSampMat,mask);
    save_volume4D(temp,filename+"samp");
    sumStatMat=SD(sumStatMat,sumSampMat);

    if (percentileThreshold>0) {
      temp.setmatrix(sumStatMat,mask);
      float min=temp.percentile(percentileThreshold,mask);
      //cerr << min << " " << percentile((Matrix)tstat_ceav.t(),percentileThreshold*100) << endl;     
      for(int i=1;i<=sumStatMat.Ncols();i++)
	if(sumStatMat(1,i)<min) sumStatMat(1,i)=min;
    }

    temp.setmatrix(sumStatMat,mask);
    save_volume4D(temp,filename+"post");
  }
}

Matrix tfce(const Matrix& tstat, const volume<float>& mask, const float delta, float height_power, float size_power, int connectivity){
  volume4D<float> input_volume;
  input_volume.setmatrix(tstat,mask);
  tfce(input_volume[0],height_power,size_power,connectivity,0,delta);
  return(input_volume.matrix(mask));
}

void checkInput(const short st,const  Matrix& dm,const  Matrix& tc,const  Matrix& fc){
  if (dm.Nrows()!=st) throw Exception("number of rows in design matrix doesn't match number of \"time points\" in input data!"); 
  if (tc.Ncols()!=dm.Ncols()) throw Exception("number of columns in t-contrast matrix doesn't match number of columns in design matrix!");
  if (fc.Ncols() !=0 && fc.Ncols()!=tc.Nrows()) throw Exception("number of columns in f-contrast matrix doesn't match number of rows in t-contrast matrix!");
}

void Initialise(ranopts& opts, volume<float>& mask, Matrix& datam, Matrix& tc, Matrix& dm, Matrix& fc, Matrix& gp)
{
  if (opts.parallelData.value()) opts.verbose.set_value("false");
  if ( opts.randomSeed.set()) srand(opts.randomSeed.value());
  if ( opts.randomSeed.set() && opts.verbose.value() ) cout << "Seeding with " << opts.randomSeed.value() << endl;
  if (opts.verbose.value()) cout << "Loading Data: "; 
  short sx,sy,sz,st;
  {
      FSLIO *IP1;
      IP1 = NewFslOpen(opts.in_fileroot.value(), "r");
      if (IP1==0) throw Exception(("Failed to read volume "+opts.in_fileroot.value()).c_str()); 
      FslGetDim(IP1,&sx,&sy,&sz,&st);
      FslClose(IP1);
  }

  if(opts.one_samp.value())
  {
    dm.ReSize(st,1);
    dm=1;
    tc.ReSize(1,1);
    tc=1;
  }
  else if ( opts.dm_file.value()=="" || opts.tc_file.value()=="" ) throw Exception("Randomise requires a design matrix and contrast as input");
  if (opts.dm_file.value()!="") dm=read_vest(opts.dm_file.value());
  if (opts.tc_file.value()!="") tc=read_vest(opts.tc_file.value());
  if (opts.fc_file.value()!="") fc=read_vest(opts.fc_file.value());
  if (opts.gp_file.value()!="") gp=read_vest(opts.gp_file.value());
  else {
    gp.ReSize(dm.Nrows(),1);
    gp=1;
  }
  checkInput(st,dm,tc,fc);

  if (opts.parallelData.value()) {
    cout << opts.n_perm.value() << " " << tc.Nrows() << " " << opts.out_fileroot.value() << endl;
    exit(0);
  }

  FSLIO *IP1;
  IP1 = NewFslOpen(opts.in_fileroot.value(), "r");
  volume4D<float> data(sx,sy,sz,1);
  float* tbuffer;
  tbuffer = new float[sx*sy*sz];
  for (int t=0;t<st;t++) 
  {
    FslReadBuffer(IP1,tbuffer);
    data[0].reinitialize(sx,sy,sz,tbuffer,false);
    if (t==0)
    {
      if (opts.maskname.value()!="") 
      {
	read_volume(mask,opts.maskname.value());
	if (!samesize(data[0],mask)) throw Exception("Mask dimensions do not match input data dimensions!");
      }
      else mask = data[0];
      set_volume_properties(IP1,mask);
      mask.binarise(0.0001);  
    } 
    if (t!=0) datam&= data.matrix(mask);
    else datam=data.matrix(mask);
    if (opts.verbose.value()) cout << "*" << flush; 
  }
  delete [] tbuffer;
  FslClose(IP1);

  if (opts.demean_data.value()) datam=remmean(datam);
  if (opts.verbose.value()) cout << endl; 

  if (opts.voxelwise_ev_numbers.set() && opts.voxelwise_ev_filenames.set())
  {
    volume4D<float> input;
    voxelwise_ev_numbers=opts.voxelwise_ev_numbers.value();  
    if(opts.voxelwise_ev_filenames.value().size() != voxelwise_ev_numbers.size())
      throw Exception("Number of input voxelwise_ev_filenames must match number of voxelwise_ev_numbers");
    voxelwise_evs.resize(voxelwise_ev_numbers.size());
  
    for(unsigned int i=0; i<voxelwise_ev_numbers.size(); i++)      
      {
	if(voxelwise_ev_numbers[i]>dm.Ncols())
	  throw Exception("voxelwise_ev_numbers option specifies a number greater than number of design EVs)");
	if (opts.verbose.value()) cout << "Loading voxelwise ev: " << opts.voxelwise_ev_filenames.value().at(i) << " for EV " << voxelwise_ev_numbers[i] << endl;
	read_volume4D(input,opts.voxelwise_ev_filenames.value().at(i));
        voxelwise_evs[i]=input.matrix(mask);
      }
  }

  if (opts.verbose.value()) cout << "Data loaded" << endl;
  if (opts.tfce2D.value()) {
    opts.tfce.set_value("true");
    opts.tfce_height.set_value("2");     
    opts.tfce_size.set_value("1");     
    opts.tfce_connectivity.set_value("26");  
  }
}

Matrix PermutedDesign(const Matrix& originalDesign,const ColumnVector& permutation,const bool multiply){
  Matrix output=originalDesign;
  for(int row=1;row<=originalDesign.Nrows();row++)
  {
    if (multiply) output.Row(row)=originalDesign.Row(row)*permutation(row);
    else output.Row(row) << originalDesign.Row(int(permutation(row)));
  }
  return output;
}

void ols(const Matrix& data,const Matrix& des,const Matrix& tc, Matrix& cope,Matrix& varcope,float dof){
  Matrix pdes = pinv(des);
  varcope=diag(tc*pdes*pdes.t()*tc.t());
  pdes*=data;
  cope=tc*pdes;
  Matrix res=data-des*pdes;
  Matrix sigsq=sum(SP(res,res))/dof;
  varcope*=sigsq; 
}

void ols_var_sm(const Matrix& data,const Matrix& des,const Matrix& tc, Matrix& cope,Matrix& varcope,const volume<float>& mask,const volume<float>& mask_sm,float sigma_mm, float dof){
  volume4D<float> sigsqvol;
  
  Matrix pdes = pinv(des);
  Matrix prevar=diag(tc*pdes*pdes.t()*tc.t());
  Matrix pe=pdes*data;
  cope=tc*pe;
  Matrix res=data-des*pe;
  Matrix sigsq=sum(SP(res,res))/dof;
  sigsqvol.setmatrix(sigsq,mask);
  sigsqvol[0]=smooth(sigsqvol[0],sigma_mm);
  sigsqvol[0]/=mask_sm;
  sigsq=sigsqvol.matrix(mask);
  varcope=prevar*sigsq;
}

void OutputStat(const ParametricStatistic input,const volume<float>& mask, const int nPerms,string statLabel,const string fileRoot,const bool outputText, const bool outputRaw=true)
{ 
volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
long nVoxels(input.originalStatistic.Ncols());
Matrix currentStat(1,nVoxels);
 output.setmatrix(input.originalStatistic.Row(1),mask);
 if (outputRaw) save_volume4D(output,fileRoot+statLabel);
 RowVector distribution = input.maximumDistribution.Row(1);    
 if (outputText)
 {
   ofstream output_file((fileRoot+"_corrp"+statLabel+".txt").c_str());
   output_file << distribution.t();
   output_file.close();
 }
 SortAscending(distribution);
 currentStat=0;
 for(int i=1; i<=nVoxels; i++)
   for(int j=nPerms; j>=1; j--)
     if (input.originalStatistic(1,i)>distribution(j))
     {
       currentStat(1,i) = float(j)/nPerms;
       j=0;
     }
 output.setmatrix(currentStat,mask);
 save_volume4D(output,fileRoot+"_corrp"+statLabel);
 if (input.storingUncorrected) {
   output.setmatrix(input.uncorrectedStatistic.Row(1)/float(nPerms),mask);
   save_volume4D(output,fileRoot+"_p"+statLabel);
 }
}

Matrix calculateFStat(Matrix& model,Matrix& data,Matrix& W2,const float dof)
{ 
  Matrix E,Gamma1;
  Gamma1=pinv(model)*data;
  E=data-model*Gamma1;
  E=sum(SP(E,E))/dof;
  Matrix GammaModel=Gamma1.t()*model.t()*model;
  Gamma1=sum(SP(GammaModel.t(),Gamma1));
  return(SD(Gamma1,E)); 
}

void calculatePermutationStatistics(ranopts& opts, const volume<float>& mask, Matrix& datam, Matrix& tc, Matrix& dm,int tstatnum, Matrix& NewW2, float dof, Permuter& permuter)
{
  int nVoxels(datam.Ncols());
  volume4D<float> temp4D,tstat4D(mask.xsize(),mask.ysize(),mask.zsize(),1);
  float tfce_delta(0), clusterThreshold(0), massThreshold(0);
  if (tstatnum>=0) clusterThreshold=opts.cluster_thresh.value();
  else clusterThreshold=opts.f_thresh.value();
  if (tstatnum>=0) massThreshold=opts.clustermass_thresh.value();
  else massThreshold=opts.fmass_thresh.value();
  bool isNormalising( opts.cluster_norm.value() && tstatnum >=0 ), lowram(false);

  string statLabel;
  if (tstatnum<0) statLabel="_fstat"+num2str(-tstatnum);
  else statLabel="_tstat"+num2str(tstatnum);
  // prepare smoothed mask for use (as a convolution renormaliser) in variance smoothing if required
  volume<float> smoothedMask;
  if(opts.var_sm_sig.value()>0) smoothedMask=smooth(mask,opts.var_sm_sig.value());
  // containers for different inference distribution
  ParametricStatistic clusters, clusterMasses, clusterNormals, clusterEnhanced, clusterEnhancedNormals, voxels;
  Matrix dmperm, tstat, tstat_ce, cope, varcope, previousTFCEStat, copesmall, varcopesmall;
  cope.ReSize(1,nVoxels);
  varcope=cope;

  unsigned long nPerms=(unsigned long)permuter.reportRequiredPermutations(opts.verbose.value());
  if ( !((clusterThreshold>0) || (massThreshold>0) || opts.tfce.value() || opts.voxelwiseOutput.value()) )
  {
    cout << "Warning! No output options selected. Outputing raw tstat only" << endl;
    nPerms=1;
  }
  // resize the containers for the relevant inference distributions
  voxels.setup(1,nPerms,nVoxels,false,true);
  if ( clusterThreshold >0 )  
    clusters.setup(1,nPerms,nVoxels,isNormalising);
  if ( massThreshold>0 ) 
    clusterMasses.setup(1,nPerms,nVoxels,false);
  if ( clusters.isAveraging ) 
    clusterNormals.setup(1,nPerms,nVoxels,false);
  if ( opts.tfce.value() )    
    clusterEnhanced.setup(1,nPerms,nVoxels,isNormalising,true);
  if ( clusterEnhanced.isAveraging ) { 
    clusterEnhancedNormals.setup(1,nPerms,nVoxels,false);
    try { previousTFCEStat.ReSize(nPerms,clusterEnhanced.sumStatMat.Ncols());} //between 5e17 - 5e18 values for a 2gb machine
    catch (...) {cerr << "using lowram" << endl; lowram=true;}         
  }
  
  
  ColumnVector clustersizes,permvec;
  vector<ColumnVector> previousPermutations;
  previousPermutations.reserve(nPerms);
  for(unsigned long perm=1; perm<=nPerms; perm++){

    permvec = permuter.nextPermutation(perm,previousPermutations,opts.verbose.value(), isNormalising || opts.outputText.value());    
    dmperm=PermutedDesign(dm,permvec,permuter.isFlipping);

    if (opts.voxelwise_ev_numbers.set() && opts.voxelwise_ev_filenames.set())
    {
      for(int voxel=1;voxel<=datam.Ncols();voxel++)
      {
	Matrix dmtemp=dm;
	for (unsigned int ev=0; ev<voxelwise_evs.size(); ev++)
          dmtemp.Column(voxelwise_ev_numbers[ev])=voxelwise_evs[ev].Column(voxel);
	dof=ols_dof(dmtemp); 
	if (opts.demean_data.value()) dof--;
	dmperm=PermutedDesign(dmtemp,permvec,permuter.isFlipping);
	if(opts.var_sm_sig.value()==0) ols(datam.Column(voxel),dmperm,tc,copesmall,varcopesmall,dof);    
	cope.Column(voxel)=copesmall;
	varcope.Column(voxel)=varcopesmall;
      }
    }
    else
    {
      if(opts.var_sm_sig.value()==0) ols(datam,dmperm,tc,cope,varcope,dof);   
      else ols_var_sm(datam,dmperm,tc,cope,varcope,mask,smoothedMask,opts.var_sm_sig.value(),dof);  
    }

    if ( tstatnum >=0) tstat=SD(cope,sqrt(varcope));
    else tstat=calculateFStat(dmperm,datam,NewW2,dof);
    voxels.store(tstat,perm);
    if (opts.tfce.value())
    {
      if (perm==1) tfce_delta=tstat.Maximum()/100.0;  // i.e. 100 subdivisions of the max input stat height
      tstat_ce=tfce(tstat,mask,tfce_delta,opts.tfce_height.value(),opts.tfce_size.value(),opts.tfce_connectivity.value());
      clusterEnhanced.store(tstat_ce,perm);
      if(!lowram && clusterEnhanced.isAveraging ) previousTFCEStat.Row(perm)=tstat_ce.Row(1);
    }
    if (opts.output_permstat.value()) tstat4D.setmatrix(tstat,mask);
    if (opts.output_permstat.value()) save_volume4D(tstat4D,opts.out_fileroot.value()+"_rawstat" + statLabel + "_" + ((num2str(perm)).insert(0,"00000")).erase(0,num2str(perm).length()));

   if ( clusterThreshold > 0 ) { //cluster thresholding
      tstat4D.setmatrix(tstat,mask);
      tstat4D.binarise( clusterThreshold );
      volume<int> clusterLabels=connected_components(tstat4D[0],clustersizes,CLUST_CON);
      clusters.store(clusterLabels,clustersizes,mask,1,perm);
   } //end of cluster tresholding
     
    
    if ( massThreshold > 0 ) { //cluster mass thresholding
      tstat4D.setmatrix(tstat,mask);
      temp4D=tstat4D;
      tstat4D.binarise(massThreshold);
      volume<int> clusterLabels=connected_components(tstat4D[0],clustersizes,CLUST_CON);
      clustersizes=0;	
      for(int z=0; z<mask.zsize(); z++)
	for(int y=0; y<mask.ysize(); y++)
	  for(int x=0; x<mask.xsize(); x++)
	    if(clusterLabels(x,y,z)>0)
	      clustersizes(clusterLabels(x,y,z))=clustersizes(clusterLabels(x,y,z))+temp4D[0](x,y,z);
      clusterMasses.store(clusterLabels,clustersizes,mask,1,perm);
    }
  }
  //End of Permutations
    
   //Rerun perms for clusternorm
  if ( isNormalising )
  { 
    if ( clusters.isAveraging ) {
      clusters.average(opts.out_fileroot.value()+statLabel+"_clusternorm",0,mask);
      temp4D.setmatrix(clusters.sumStatMat,mask);
    }
    if (clusterEnhanced.isAveraging)
      clusterEnhanced.average(opts.out_fileroot.value()+statLabel+"_tfcenorm",0.02,mask);

    for(unsigned long perm=1; perm<=nPerms; perm++)
    {
      if (opts.verbose.value()) cout << "Starting second-pass " << perm << endl;
      if ( clusters.isAveraging || ( clusterEnhanced.isAveraging && lowram ) ) //Regenerate stats
      { 
	permvec=permuter.returnPreviousTruePermutation(perm,previousPermutations);
	dmperm=PermutedDesign(dm,permvec,permuter.isFlipping);
	if(opts.var_sm_sig.value()==0) ols(datam,dmperm,tc,cope,varcope,dof);   
	else ols_var_sm(datam,dmperm,tc,cope,varcope,mask,smoothedMask,opts.var_sm_sig.value(),dof);   
	tstat=SD(cope,sqrt(varcope));
      }
      if ( clusterEnhanced.isAveraging )
      {
	if (!lowram) tstat_ce=previousTFCEStat.Row(perm);
	else tstat_ce=tfce(tstat,mask,tfce_delta,opts.tfce_height.value(),opts.tfce_size.value(),opts.tfce_connectivity.value());
	tstat_ce=SD(tstat_ce,clusterEnhanced.sumStatMat); 
	clusterEnhancedNormals.store(tstat_ce,perm);
      }
      if ( clusters.isAveraging )
      { 
	tstat4D.setmatrix(tstat,mask);
	tstat4D.binarise(clusterThreshold);
	volume<int> clusterLabels=connected_components(tstat4D[0],clustersizes,CLUST_CON);
	ColumnVector entries,cluster(clustersizes.Nrows());
	cluster=0;
	entries=cluster;
	for(int z=0; z<mask.zsize(); z++)
	  for(int y=0; y<mask.ysize(); y++)
	    for(int x=0; x<mask.xsize(); x++)
	      if (clusterLabels(x,y,z))
	      {
		cluster(clusterLabels(x,y,z))+=temp4D(x,y,z,0);
		entries(clusterLabels(x,y,z))++;
	      }
	clustersizes=SD(clustersizes,SD(cluster,entries));
	clusterNormals.store(clusterLabels,clustersizes,mask,1,perm);
      }
    }
  }
   
  //OUTPUT Routines 
  tstat4D.setmatrix(voxels.originalStatistic.Row(1),mask);
  save_volume4D(tstat4D,opts.out_fileroot.value()+statLabel);
  if ( opts.voxelwiseOutput.value() ) OutputStat(voxels,mask,nPerms,statLabel,opts.out_fileroot.value()+"_vox",opts.outputText.value(),false);
  if ( clusterThreshold > 0 ) OutputStat(clusters,mask,nPerms,statLabel,opts.out_fileroot.value()+"_clustere",opts.outputText.value(),false);
  if ( massThreshold > 0 )    OutputStat(clusterMasses,mask,nPerms,statLabel,opts.out_fileroot.value()+"_clusterm",opts.outputText.value(),false);
  if ( clusters.isAveraging ) OutputStat(clusterNormals,mask,nPerms,statLabel,opts.out_fileroot.value()+"_clustern",opts.outputText.value(),false);  
  if ( opts.tfce.value() )    OutputStat(clusterEnhanced,mask,nPerms,statLabel,opts.out_fileroot.value()+"_tfce",opts.outputText.value());
  if ( clusterEnhanced.isAveraging ) OutputStat(clusterEnhancedNormals,mask,nPerms,statLabel,opts.out_fileroot.value()+"_tfcen",opts.outputText.value());

  if (opts.outputText.value()) {
    ofstream output_file((opts.out_fileroot.value()+"_perm"+statLabel+".txt").c_str());
    for(unsigned long perm=1; perm<=nPerms; perm++) 
      output_file << permuter.returnPreviousTruePermutation(perm,previousPermutations).t();    
    output_file.close();
  }
  
}

void convertContrast(const Matrix& inputDesign,const Matrix& inputContrast,const Matrix& inputData,Matrix& outputModel,Matrix& outputContrast, Matrix& outputData, Matrix& outputConfound )
{
    int r(inputContrast.Nrows()),p(inputContrast.Ncols());
    Matrix tmp=(IdentityMatrix(p)-inputContrast.t()*pinv(inputContrast.t()));
    Matrix U,V;
    DiagonalMatrix D;
    SVD(tmp, D, U, V);
    Matrix c2=U.Columns(1,p-r);
    c2=c2.t();
    Matrix C = inputContrast & c2;
    Matrix W=inputDesign*C.i();
    Matrix W1=W.Columns(1,r);
    Matrix W2=W.Columns(r+1,W.Ncols());

    if (  W2.Nrows() && W2.Ncols() )
    {
      outputModel=W1-W2*pinv(W2)*W1;
      outputData=(IdentityMatrix(W2.Nrows())-W2*pinv(W2))*inputData;
    }
    else
    {
      outputModel=W1;
      outputData=inputData;
    }
    outputContrast=IdentityMatrix(r);
    outputConfound=W2;
}


void analyseContrast(const Matrix& inputContrast,Matrix& dm,Matrix& datam,volume<float>& mask, Matrix& gp,const int& contrastNo,ranopts& opts)
{
  //-ve num for f-stat contrast
  Matrix NewModel,NewCon,NewDataM,NewW2;
  Permuter permuter;

  convertContrast(dm,inputContrast,datam,NewModel,NewCon,NewDataM,NewW2);

  bool oneRegressor( inputContrast.SumAbsoluteValue() == inputContrast.MaximumAbsoluteValue() );

  permuter.createPermutationGroups(remmean(dm)*inputContrast.t(),gp,(contrastNo>0 && oneRegressor),opts.n_perm.value()); 
  if(permuter.isFlipping) cout << "One-sample design detected; sign-flipping instead of permuting." << endl;
  if(opts.verbose.value() || opts.how_many_perms.value()) 
  {
    if(permuter.isFlipping) cout << permuter.uniquePermutations[0] << " sign-flips required for exhaustive test";
    else cout << permuter.uniquePermutations[0] << " permutations required for exhaustive test";
    if (contrastNo>0)  cout << " of t-test " << contrastNo << endl;
    if (contrastNo==0) cout << " of all t-tests " << endl;
    if (contrastNo<0)  cout << " of f-test " << abs(contrastNo) << endl;
    if(opts.how_many_perms.value()) return;
  }  
 
  if (opts.confoundMethod.value()==1) {
    NewCon=inputContrast; 
    NewModel=dm;  
  }
                          
  if (opts.confoundMethod.value()==2 && NewW2.Ncols() > 0)  
  {
    Matrix nuisanceContrast(NewCon.Nrows(),NewW2.Ncols());
    nuisanceContrast=0;
    NewCon = NewCon | nuisanceContrast;
    NewModel = NewModel | NewW2;             
  }

  float dof=ols_dof(dm); 
  if (opts.demean_data.value()) dof--;
  calculatePermutationStatistics(opts,mask,NewDataM,NewCon,NewModel,contrastNo,NewW2,dof,permuter); 
}


void analyseFContrast(Matrix& fc,Matrix& tc,Matrix& model,Matrix& data,volume<float>& mask,Matrix& gp,ranopts& opts)
{   
   for( int fstat=1; fstat<=fc.Nrows() ; fstat++ ) 
   {
      Matrix fullFContrast(0,tc.Ncols());
      for (int tcon=1; tcon<=fc.Ncols() ; tcon++ )
	if (fc(fstat,tcon)==1) fullFContrast &= tc.Row(tcon);
      analyseContrast(fullFContrast,model,data,mask,gp,-fstat,opts);
   }
}

int main(int argc,char *argv[]){
  Log& logger = LogSingleton::getInstance();
  ranopts& opts = ranopts::getInstance();
  opts.parse_command_line(argc,argv,logger);
  Matrix model, Tcontrasts, Fcontrasts, data,grp_lab;
  volume<float> mask;
  try { Initialise(opts,mask,data,Tcontrasts,model,Fcontrasts,grp_lab); 
  bool needsDemean=true;
  for (int i=1;i<=model.Ncols();i++) if ( fabs( (model.Column(i)).Sum() ) > 0.0001 ) needsDemean=false;
  if (needsDemean && !opts.demean_data.value()) cerr << "Warning: All design columns have zero mean - consider using the -D option to demean your data" << endl;
  if (!needsDemean && opts.demean_data.value()) cerr << "Warning: You have demeaned your data, but at least one design column has non-zero mean" << endl;
  if(opts.fc_file.value()!="") analyseFContrast(Fcontrasts,Tcontrasts,model,data,mask,grp_lab,opts); 
  for (int tstat=1; tstat<=Tcontrasts.Nrows() ; tstat++ )  analyseContrast(Tcontrasts.Row(tstat),model,data,mask,grp_lab,tstat,opts); 
  }
  catch(Exception& e) 
  { 
    cerr << "ERROR: Program failed" <<  e.what() << endl << endl << "Exiting" << endl; 
    return 1;
  }
  catch(...) 
  { 
    cerr << "ERROR: Program failed, unknown exception" << endl << endl << "Exiting" << endl; 
    return 1;
  }
  return 0;
}

//Permuter Class
void Permuter::createPermutationGroups(const Matrix& design,const Matrix& groups,const bool oneNonZeroContrast,const long requiredPermutations)
{
  nGroups=int(groups.Maximum())+1;
  nSubjects=design.Nrows();
  ColumnVector labels = createDesignLabels(design);
  isFlipping = ( (labels.Maximum()==1) && oneNonZeroContrast );
  int active=0;  
  originalLocations.resize(nGroups);
  permutedLabels.resize(nGroups);
  for(int row=1;row<=nSubjects;row++) if(groups(row,1)==0 || (design.Row(row).Sum()==0 && !isFlipping)) active++;
  originalLocations[0].ReSize(active);
  permutedLabels[0].ReSize(active);
  for(int row=nSubjects;row>=1;row--) if(groups(row,1)==0 || (design.Row(row).Sum()==0 && !isFlipping)) originalLocations[0](active--)=row;
  for(int group=1;group<=groups.Maximum();group++)
  {
     active=0;
     for(int row=1;row<=nSubjects;row++)
       if(groups(row,1)==group && (design.Row(row).Sum()!=0 || isFlipping)) active++;
     originalLocations[group].ReSize(active);
     permutedLabels[group].ReSize(active);
     for(int row=nSubjects;row>=1;row--) //Now work backwards to fill in the row numbers
       if(groups(row,1)==group && (design.Row(row).Sum()!=0 || isFlipping)) originalLocations[group](active--)=row;
  }
  initialisePermutationGroups(labels,requiredPermutations);
}

ColumnVector Permuter::returnPreviousTruePermutation(const long permutationNumber,const vector<ColumnVector>& previousPermutations)
{
  if (isFlipping) 
    return previousPermutations[permutationNumber-1];
  else {    
    ColumnVector permvec(unpermutedVector);
    for(long perm=1; perm<=permutationNumber; perm++) 
      createTruePermutation(previousPermutations[perm-1],previousPermutations[perm-1-int(perm!=1)],permvec);  
    return permvec;
  }
}

void Permuter::initialisePermutationGroups(const ColumnVector& designLabels,const long requiredPermutations)
{
  truePermutation.ReSize(nSubjects);
  for(int i=1;i<=nSubjects;i++) truePermutation(i)=i;
  if (isFlipping) truePermutation=1;
  unpermutedVector=truePermutation;
  uniquePermutations.resize(nGroups);
  uniquePermutations[0]=1;
  for(int group=0;group<nGroups;group++)
  { 
    for(int row=1;row<=permutedLabels[group].Nrows();row++) 
      permutedLabels[group](row)=designLabels((int)originalLocations[group](row));
    if (group>0) uniquePermutations[group]=computeUniquePermutations(permutedLabels[group],isFlipping);
    uniquePermutations[0]*=uniquePermutations[group];
  }
  isRandom=!(requiredPermutations==0 || requiredPermutations>=uniquePermutations[0]);
  if (isRandom) finalPermutation=requiredPermutations;
  else finalPermutation=uniquePermutations[0];
}

ColumnVector Permuter::permutationVector()
{
ColumnVector newvec(nSubjects); 
   for(int group=0;group<nGroups;group++)
     for(int row=1;row<=permutedLabels[group].Nrows();row++) 
       newvec((int)originalLocations[group](row))=permutedLabels[group](row);
   return newvec;
}

void Permuter::createTruePermutation(const ColumnVector& newLabels,ColumnVector copyOldLabels,ColumnVector& permvec)
{
  if (isFlipping) permvec=permutationVector();
  else
  {
    for(int k=1;k<=newLabels.Nrows();k++)	       
      if(newLabels(k)!=copyOldLabels(k))
        for(int l=1;l<=newLabels.Nrows();l++) 
	  if(newLabels(l)!=copyOldLabels(l) && copyOldLabels(l)==newLabels(k) )
 	   {
	     swap(permvec(l),permvec(k));
	     swap(copyOldLabels(l),copyOldLabels(k));
           }
  }
}

ColumnVector Permuter::nextPermutation(const long permutationNumber)
{
  for(int group=1;group<nGroups;group++)
  {
    long cycle = 1;
    for(int j=group+1;j<nGroups;j++) cycle*=(long)uniquePermutations[j];
    if ( (permutationNumber-1)%cycle==0 ) 
    { 
      if(isFlipping) nextFlip(permutedLabels[group]);
      else nextShuffle(permutedLabels[group]);
    }
  }
  return(permutationVector());
}

ColumnVector Permuter::nextPermutation(const long permutationNumber,vector<ColumnVector>& previousPermutations, const bool printStatus, const bool isStoring)
{
  if (permutationNumber!=1 && printStatus) cout << "Starting permutation " << permutationNumber << endl;
  else if (printStatus) cout << "Starting permutation " << permutationNumber << " (Unpermuted data)" << endl;

  ColumnVector currentLabels=permutationVector();
  ColumnVector newPermutation;
  do
  {
    if (permutationNumber!=1) newPermutation=nextPermutation(permutationNumber);
  } while(isRandom && isPreviousPermutation(previousPermutations,newPermutation));
  if(isStoring || isRandom) previousPermutations.push_back(permutationVector());
  createTruePermutation(permutationVector(),currentLabels,truePermutation);
  return(truePermutation);
}

bool Permuter::isPreviousPermutation(const vector<ColumnVector>& previousPermutations,const ColumnVector& newPermutation){
  for(int i=previousPermutations.size()-1; i>=0; i--)
    if(newPermutation==previousPermutations[i]) return true;
  return false;
  }

void Permuter::nextShuffle(ColumnVector& perm){
   vector<int> temp;
   for (int i=1;i<=perm.Nrows();i++) temp.push_back((int)perm(i));
   if (isRandom) random_shuffle(temp.begin(),temp.end());
   else next_permutation(temp.begin(),temp.end());
   for (int i=1;i<=perm.Nrows();i++) perm(i)=temp[i-1];
}

void Permuter::nextFlip(ColumnVector& flip){
       
  if (isRandom)
  {
    for(int i=1;i<=flip.Nrows();i++)
    {
      float tmp=(float)rand()/RAND_MAX;
      if(tmp > 0.5) flip(i)=1;
      else  flip(i)=-1;
    }     
  }
  else for (int n=flip.Nrows();n>0;n--)
    if(flip(n)==1)
    {
      flip(n)=-1;
      if (n<flip.Nrows()) flip.Rows(n+1,flip.Nrows())=1;
      return;
    } 
}

double Permuter::computeUniquePermutations(const ColumnVector& labels,const bool calculateFlips){
  if (calculateFlips) return std::pow(2.0,labels.Nrows());
  ColumnVector label_counts((int)labels.MaximumAbsoluteValue());
  label_counts=0;
  for(int i=1; i<=labels.Nrows(); i++) label_counts(int(labels(i)))++;
  double yo = lgam(labels.Nrows()+1);
  for(int i=1; i<=labels.MaximumAbsoluteValue(); i++)
    yo -= lgam(label_counts(i)+1);
  return std::floor(exp(yo)+0.5);
}

ColumnVector Permuter::createDesignLabels(const Matrix& design){
  ColumnVector designLabels(design.Nrows());
  vector<RowVector> knownLabels;
  for(int i=1;i<=design.Nrows();i++){
    bool wasExistingLabel=false;
    for(unsigned int l=0;l<knownLabels.size();l++){
      if(design.Row(i)==knownLabels[l]){
	designLabels(i)=l+1;
	wasExistingLabel=true;
      }
    }
    if(!wasExistingLabel){
      knownLabels.push_back(design.Row(i));
      designLabels(i)=knownLabels.size();
    }
  }
  return(designLabels);
}

double Permuter::reportRequiredPermutations(const bool printToScreen)
{
  if (printToScreen)
  {
    if(isRandom) cout<<"Doing " << finalPermutation << " random permutations"<<endl;
    else cout<<"Doing all "<< finalPermutation <<" unique permutations"<<endl;
  }
  return(finalPermutation);
}

Permuter::Permuter()
{
}

Permuter::~Permuter()
{
}
