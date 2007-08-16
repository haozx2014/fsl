/*  randomise.cc
    Tim Behrens & Steve Smith & Matthew Webster (FMRIB) & Tom Nichols (UMich)
    Copyright (C) 2004-2007 University of Oxford  */
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

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "libprob.h"
#include "ranopts.h"
#include <algorithm>

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace RANDOMISE;

class permblock
{ 
public:
bool onesample;
ColumnVector num_perms;
double exhaust_perms;
int num_blocks;
int num_elements;
ColumnVector *permuted_locations;
ColumnVector *original_locations;
permblock();
~permblock();
ColumnVector next_perm(int perm,bool exhaustive);
ColumnVector buildfullperm();
double initpermblocks(bool onesample,ColumnVector labels);
void createpermblocks(const Matrix &regressor,const Matrix &group);
};

Matrix tfce(const Matrix& tstat, const volume<float>& mask, const float delta, float height_power, float size_power, int connectivity){
  volume4D<float> input_volume;
  input_volume.setmatrix(tstat,mask);
  float maxval=input_volume[0].max();
  volume<float> clusterenhance(input_volume[0]);
  clusterenhance=0;
  for (float thresh=delta; thresh<maxval; thresh+=delta)
    {
      volume<float> clusters;
      copyconvert(input_volume[0],clusters);
      clusters.binarise(thresh);

      ColumnVector clustersizes;  
      volume<int>tmpvol=connected_components(clusters,clustersizes,connectivity);
      clustersizes = pow(clustersizes,size_power) * pow(thresh,height_power);
      for(int z=0;z<input_volume[0].zsize();z++)
	for(int y=0;y<input_volume[0].ysize();y++)	    
	  for(int x=0;x<input_volume[0].xsize();x++)
	    if (tmpvol(x,y,z)>0)
	      clusterenhance(x,y,z) += clustersizes(tmpvol(x,y,z));
    }
  copyconvert(clusterenhance,input_volume[0]);
  return(input_volume.matrix(mask));
}


bool check_dims(const short st,const  Matrix& dm,const  Matrix& confounds,const  Matrix& tc,const  Matrix& fc){
  bool  ret=true;
  if (dm.Nrows()!=st){
      cerr << "number of rows in design matrix doesn't match number of \"time points\" in input data!" << endl;
      ret=false;
  }
  if(confounds.Nrows()!=0){
    if(confounds.Nrows()!=dm.Nrows()){
      cerr << "number of rows in confound matrix doesn't match number of rows in design matrix!" << endl;
      ret=false;
    }
  }
  if (tc.Ncols()!=dm.Ncols()){
    cerr << "number of columns in t-contrast matrix doesn't match number of columns in design matrix!" << endl;
    ret=false;
  }
  if(fc.Ncols() !=0){
    if (fc.Ncols()!=tc.Nrows()){
      cerr << "number of columns in f-contrast matrix doesn't match number of rows in t-contrast matrix!" << endl;
      ret=false;
    }
  }
  return ret;
}

void next_mult_vec(ColumnVector& mult, const bool &exhaustive){
  if (exhaustive)
  {
    bool finish=false;
    int n=mult.Nrows();
    while(!finish)
    {
      if(mult(n)==1)
      {
	mult(n)=-1;
	if (n<mult.Nrows()) mult.Rows(n+1,mult.Nrows())=1;
	finish=true;
      }
      else n--;
    }
  }
  else for(int i=1;i<=mult.Nrows();i++)
  {
     float tmp=rand();
     tmp/=RAND_MAX;
     if(tmp > 0.5) mult(i)=1;
     else  mult(i)=-1;
  }
}

void dm_mult(const Matrix& dm, Matrix& dm_new, const ColumnVector& mult){
  for(int i=1;i<=dm.Nrows();i++)
    dm_new.Row(i)=dm.Row(i)*mult(i);
}

void dm_permute(const Matrix& dm_in, Matrix& dm_out, const ColumnVector& r){

  if(r.Nrows() != dm_in.Nrows()){
    cerr<<"permutation vector has wrong number of elements"<<endl;
    exit(-1);
  }
  
  if ( (dm_out.Nrows()!=dm_in.Nrows()) || (dm_out.Ncols()!=dm_in.Ncols()) ){
    dm_out.ReSize(dm_in.Nrows(),dm_in.Ncols());
  }

  for(int row=1;row<=r.Nrows();row++){
    dm_out.Row(row) << dm_in.Row(int(r(row)));
  }

}

void next_perm_vec(ColumnVector& r, const bool &exhaustive){
   int *vec = new int[r.Nrows()+1];
   for (int i=1;i<=r.Nrows();i++) vec[i-1]=(int)r(i);
   if (exhaustive) next_permutation( vec, vec + r.Nrows() );
   else random_shuffle( vec, vec + r.Nrows() );
   for (int i=1;i<=r.Nrows();i++) r(i)=vec[i-1];
   delete [] vec;
}

ColumnVector make_dm_labels(const Matrix& dm){
  ColumnVector labels_orig(dm.Nrows());
  vector<RowVector> label_lu;
  bool found_it=false;
  for(int i=1;i<=dm.Nrows();i++){
    found_it=false;
    for(unsigned int l=0;l<label_lu.size();l++){
      if(dm.Row(i)==label_lu[l]){
	labels_orig(i)=l+1;
	found_it=true;
      }
    }
    if(!found_it){
      label_lu.push_back(dm.Row(i));
      labels_orig(i)=label_lu.size();
    }
  }
  return(labels_orig);
}

bool check_perm(vector<ColumnVector>& oldperms,const ColumnVector& newperm){
  bool isin=false;
  for(int i=oldperms.size()-1; i>=0; i--){
    if(newperm==oldperms[i]){
      isin=true;
      break;
    }
  } 
  return !isin;
}

void ols(const Matrix& data,const Matrix& des,const Matrix& tc, Matrix& cope,Matrix& varcope,float dof){
  // ols
  // data is t x v
  // des is t x ev (design matrix)
  // tc is cons x ev (contrast matrix)
  // cope and varcope will be cons x v
  // but will be resized if they are wrong
  // hence may be passed in uninitialised
  // TB 2004
  if(data.Nrows() != des.Nrows()){
    cerr <<"MISCMATHS::ols - data and design have different number of time points"<<endl;
    exit(-1);
  }
  if(des.Ncols() != tc.Ncols()){
    cerr <<"MISCMATHS::ols - design and contrast matrix have different number of EVs"<<endl;
    exit(-1);
  }  
  Matrix pdes = pinv(des);
  Matrix prevar=diag(tc*pdes*pdes.t()*tc.t());
  Matrix pe=pdes*data;
  cope=tc*pe;
  Matrix res=data-des*pe;
  Matrix sigsq=sum(SP(res,res))/dof;
  varcope=prevar*sigsq; 
}

void ols_var_sm(const Matrix& data,const Matrix& des,const Matrix& tc, Matrix& cope,Matrix& varcope,const volume<float>& mask,const volume<float>& mask_sm,float sigma_mm, float dof){
  // ols_var_sm
  // data is t x v
  // des is t x ev (design matrix)
  // tc is cons x ev (contrast matrix)
  // cope and varcope will be cons x v
  // but will be resized if they are wrong
  // hence may be passed in uninitialised
  // TB 2004

  if(data.Nrows() != des.Nrows()){
    cerr <<"RANDOMISE::ols_var_sm - data and design have different number of time points"<<endl;
    exit(-1);
  }
  if(des.Ncols() != tc.Ncols()){
    cerr <<"RANDOMISE::ols_var_sm - design and contrast matrix have different number of EVs"<<endl;
    exit(-1);
  }  
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

double compute_nperms(const ColumnVector& labels){
    int num_labels=int(labels.MaximumAbsoluteValue());
    ColumnVector label_counts(num_labels);
    label_counts=0;
    for(int i=1; i<=labels.Nrows(); i++){
      label_counts(int(labels(i)))+=1;
    }
    double yo = lgam(labels.Nrows()+1);
    for(int i=1; i<=num_labels; i++)
      yo -= lgam(label_counts(i)+1);
    return std::floor(exp(yo)+0.5);
}

void VoxelSignificance(volume4D <float>&output, const volume <int>&clustermap, RowVector &maxpermsize, const ColumnVector& clustersize, const int n_perms)
{
  output=0;
  SortAscending(maxpermsize);
  for(int currentcluster=1; currentcluster<=clustersize.Nrows(); currentcluster++)
    if(clustersize(currentcluster)>0) 
      for(int passed=n_perms; passed>=1; passed--)
	if(clustersize(currentcluster)>maxpermsize(passed))
	{
	  for(int z=0;z<output.zsize();z++) 
	    for(int y=0;y<output.ysize();y++)
	      for(int x=0;x<output.xsize();x++)
		if(clustermap(x,y,z)==currentcluster)
		  output(x,y,z,0) = float(passed)/n_perms;
          passed=0;
        }
}

void OutputStat(const Matrix &original_stat,const Matrix &input_dist,const volume<float> &mask, const int n_perms,const string name, ranopts& opts,string tstatnumber)
{ 
volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
RowVector tmpvec;
Matrix current_stat(1, original_stat.Ncols());
  for(int row=1; row<=original_stat.Nrows(); row++){ 
    output.setmatrix(original_stat.Row(row),mask);
    if (tstatnumber=="0") tstatnumber=num2str(row);
    save_volume4D(output,opts.out_fileroot.value()+name+tstatnumber);
    tmpvec=input_dist.Row(row);    // max stat
    if (opts.verbose.value())
      {
	ofstream output_file;
	string foo;
	foo = opts.out_fileroot.value()+"_max"+name+tstatnumber+".txt";
	output_file.open(foo.c_str(),ofstream::out);
	output_file << tmpvec.t();
	output_file.close();
      }
    SortAscending(tmpvec);
    current_stat=0;
    for(int i=1; i<=original_stat.Ncols(); i++)
      for(int j=n_perms; j>=1; j--)
	if (original_stat(row,i)>tmpvec(j))
	{
	  current_stat(1,i) = float(j)/n_perms;
          j=0;
        }
    output.setmatrix(current_stat,mask);
    save_volume4D(output,opts.out_fileroot.value()+"_max"+name+tstatnumber);
  }
}

void OutputClusterStat(const Matrix &original_stat,const Matrix &input_dist,const volume<float> &mask, const int n_perms,const string name, ranopts& opts,string tstatnumber,float threshold,volume4D<float> optional_volume,const string mode)
{ 
volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
 output.copyproperties(mask);
volume4D<float> tstat4D(mask.xsize(),mask.ysize(),mask.zsize(),original_stat.Nrows()),tmp_tstat4D;
ColumnVector clustersizes;
  tstat4D.setmatrix(original_stat,mask);
  if (name=="_maxcmass")  tmp_tstat4D=tstat4D;
  tstat4D.binarise(threshold);
  for(int row=1; row<=original_stat.Nrows(); row++){ //for each contrast
    volume<int> clust_label=connected_components(tstat4D[row-1],clustersizes,26);
    if (name=="_maxcn")
    {
       ColumnVector cluster(clustersizes.Nrows()),entries;
       cluster=0;
       entries=cluster;
       for(int z=0; z<clust_label.zsize(); z++)
	 for(int y=0; y<clust_label.ysize(); y++)
	   for(int x=0; x<clust_label.xsize(); x++)
	     if (clust_label(x,y,z))
	     {
	       cluster(clust_label(x,y,z))+=optional_volume(x,y,z,row-1);
                entries(clust_label(x,y,z))++;
	     }
       for(int i=1;i<=clustersizes.Nrows();i++) clustersizes(i)/=(cluster(i)/entries(i));
    }
    if (name=="_maxcmass")
    {
      clustersizes=0;
      for(int z=0; z<clust_label.zsize(); z++)
	for(int y=0; y<clust_label.ysize(); y++)
	  for(int x=0; x<clust_label.xsize(); x++)
	    if(clust_label(x,y,z)>0)
	      clustersizes(clust_label(x,y,z))=clustersizes(clust_label(x,y,z))+tmp_tstat4D[row-1](x,y,z);
    }
    RowVector tmpvec=input_dist.Row(row);
    if (tstatnumber=="0") tstatnumber=num2str(row);
    if (opts.verbose.value())
      {
	ofstream output_file;
	string foo;
	foo = opts.out_fileroot.value()+name+mode+tstatnumber+".txt";
	output_file.open(foo.c_str(),ofstream::out);
	output_file << tmpvec.t();
	output_file.close();
      }
    VoxelSignificance(output,clust_label,tmpvec,clustersizes,n_perms); 
    save_volume4D(output,opts.out_fileroot.value()+name+mode+tstatnumber);
  }
}

void LabelsToPerm(const ColumnVector &labels, ColumnVector & oldlabels,ColumnVector &permvec)
{
   for(int k=1;k<=labels.Nrows();k++)	       
     if(labels(k)!=oldlabels(k))
       for(int l=1;l<=labels.Nrows();l++) 
          if(labels(l)!=oldlabels(l) && oldlabels(l)==labels(k) )
 	  {
 	     double oldvec = permvec(l);
             permvec(l)=permvec(k);
             permvec(k)=oldvec;
             oldvec=oldlabels(l);
             oldlabels(l)=oldlabels(k);
             oldlabels(k)=oldvec;  
          }
}

long int RequiredPerms(double exhaust_perms, ranopts& opts, bool& exhaustive)
{
  exhaustive=1;
  if(opts.n_perm.value()==0){
    if (opts.verbose.value()) cout<<"will do all "<< int(exhaust_perms) <<" unique permutations"<<endl;
    return(int(exhaust_perms));
  }
  else if(opts.n_perm.value()>=exhaust_perms){
    if (opts.verbose.value()) cout<<"will do " << int(exhaust_perms) << " permutations"<<endl;
    return(int(exhaust_perms));
  }
  else{
    if (opts.verbose.value()) cout<<"doing " << opts.n_perm.value() << " permutations as requested"<<endl;
    exhaustive=0;
    return(opts.n_perm.value());
  }  

}

int Initialise(ranopts &opts, volume<float> &mask, Matrix &datam, Matrix &tc, Matrix &dm, Matrix &fc, Matrix &gp)
{
  if (opts.verbose.value()) cout << "Loading Data: "; 
  short sx,sy,sz,st;
  bool lowram=false;
  {
      FSLIO *IP1;
      IP1 = NewFslOpen(opts.in_fileroot.value(), "r");
      if (IP1==0) { imthrow("Failed to read volume "+opts.in_fileroot.value(),22); }
      FslGetDim(IP1,&sx,&sy,&sz,&st);
      FslClose(IP1);
  }

  try
  {
    volume4D<float> *data=new volume4D<float>(sx,sy,sz,(int)(st*1.8)); //1.8 is UPPER bound of mem, 1.7 is too low
     data->destroy();
  }
  catch (...)
  {  
     cerr << "Cannot load full data image, reverting to low-ram usage" << endl;
     lowram=true;
  }
  if (opts.low_ram.value() || lowram)
  {
  for (int t=0;t<st;t++) 
  {
    volume4D<float> data(sx,sy,sz,1);
    read_volume4DROI(data,opts.in_fileroot.value(),0,0,0,t,sx,sy,sz,t);
    if (t==0)
    {
      if (opts.maskname.value()!="") 
      {
         read_volume(mask,opts.maskname.value());
         if (!samesize(data[0],mask))   { cerr << "mask dimensions do not match input data dimensions!" << endl; exit(-1);}
      }
      else mask = data[0];
      mask.binarise(0.0001);  
    }
    if (t!=0) datam=datam & data.matrix(mask);
    else datam=data.matrix(mask);
    if (opts.verbose.value()) cout << "*"; 
  }
  }
  else
  {
   volume4D<float> data;
   read_volume4D(data,opts.in_fileroot.value());
   if (opts.maskname.value()!="") 
   {     
     read_volume(mask,opts.maskname.value());
     if (!samesize(data[0],mask))   { cerr << "mask dimensions do not match input data dimensions!" << endl; exit(-1);}
   }
   else  mask = meanvol(data);
   mask.binarise(0.0001);  
   datam=data.matrix(mask);
  }

  if (opts.verbose.value()) cout << endl;
  if (opts.dm_file.value()!="") dm=read_vest(opts.dm_file.value());
  if (opts.tc_file.value()!="") tc=read_vest(opts.tc_file.value());
  if(opts.one_samp.value())
  {
    dm.ReSize(st,1);
    dm=1;
    tc.ReSize(1,1);
    tc=1;
  }
  else if ( opts.dm_file.value()=="" || opts.tc_file.value()=="" )
  {
    cerr << "Error: Randomise requires a design matrix and contrast as input" << endl;   
    return -1;
  }
  if (opts.fc_file.value()!="") fc=read_vest(opts.fc_file.value());
  if (opts.gp_file.value()!="") gp=read_vest(opts.gp_file.value());
  else {
    gp.ReSize(dm.Nrows(),1);
    gp=1;
  }
  Matrix confound;
  if(opts.confound_file.value()!="") confound=read_vest(opts.confound_file.value());
  if(!check_dims(st,dm,confound,tc,fc)) exit(-1);
  if (opts.demean_data.value()) datam=remmean(datam);
  if(opts.confound_file.value()!="") {
    datam=(Identity(confound.Nrows())-confound*pinv(confound))*datam; //Only here temporarily
    dm=(Identity(confound.Nrows())-confound*pinv(confound))*dm;       //Eventually all confounds will be removed using new method
  }
  if (opts.verbose.value()) cout << "Data loaded" << endl;
  return 0;
}


Matrix generate_fstat(Matrix &model,Matrix &data,Matrix &W2,const float dof)
{ 
  Matrix E;
  Matrix Gamma2=pinv(W2)*data;
  Matrix Y_2=data-W2*Gamma2;
  Matrix Gamma1=pinv(model)*Y_2;
  E=Y_2-model*Gamma1;
  E=sum(SP(E,E));
  E=E/dof;
  Matrix Middlebit=Gamma1.t()*model.t()*model;
  Gamma1=sum(SP(Middlebit.t(),Gamma1));
  return(SD(Gamma1,E)); 
}

void DoPermutationStatistics(ranopts& opts, const volume<float> &mask, Matrix &datam, Matrix &tc, Matrix &dm,int tstatnum, Matrix &NewW2, const float dof, permblock &permbl)
{
  float tfce_delta=0;
  float threshold=0;
  string tstatnumber;
  if (tstatnum<0) tstatnumber=num2str(-tstatnum);
  else tstatnumber=num2str(tstatnum);
  // prepare smoothed mask for use (as a convolution renormaliser) in variance smoothing if required
  volume<float> mask_sm;
  if(opts.var_sm_sig.value()>0) mask_sm=smooth(mask,opts.var_sm_sig.value());
  // containers for different inference distribution
  Matrix tstat, tstat_ce, tstat_ceav, tstat_ce_orig, tstat_orig, maxdist ,maxdistCE, maxdistC,maxdistCN, maxdistCmass, vox_numbigger, cope, varcope, *oldtstatce,tstat_cenorm;
  oldtstatce=&tstat; //just to avoid uninit warning
  volume4D<float> tstat4D;
  bool exhaustive,lowram=opts.low_ram.value();
  long int n_perms=RequiredPerms(permbl.exhaust_perms,opts,exhaustive);
  int n_distrows=1;
  if (tstatnum>=0) n_distrows=tc.Nrows();
  // resize the containers for the relevant inference distributions
  maxdist.ReSize(n_distrows,n_perms);
  maxdist=0;
  vox_numbigger.ReSize(n_distrows,datam.Ncols()); //number of permuted  which are bigger than original
  vox_numbigger=0;
  if ( opts.cluster_thresh.value()>0 || opts.f_thresh.value()>0 ) maxdistC=maxdist;
  if ( opts.clustermass_thresh.value()>0 || opts.fmass_thresh.value()>0 ) maxdistCmass=maxdist;
  if ( opts.cluster_norm.value() )  maxdistCN=maxdist;
  if ( opts.tfce.value() ) maxdistCE=maxdist;
  tstat4D.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),n_distrows);
  volume4D<float> sumclusters,sumcluster_samples(mask.xsize(),mask.ysize(),mask.zsize(),n_distrows);
  sumcluster_samples=0;
  sumclusters=sumcluster_samples;
 
  //Things needed to generate the permutation or multiplication vectors etc.
  Matrix dmperm=dm;
  volume<int> clust_label;
  ColumnVector clustersizes,permvec(dm.Nrows()); 
  for(int i=1;i<=permvec.Nrows();i++) 
    if(!permbl.onesample) permvec(i)=i;//use permvec as true permutation vector
    else permvec(i)=1;                      //use permvec as multiplication vector
  vector<ColumnVector> oldperms,oldlabs;    //store old true perm/mult vecs and labels to avoid redundancy.
  for(int perm=1; perm<=n_perms; perm++){
    if (opts.verbose.value()) cout << "starting permutation " << perm << endl;
    // Get next design matrix (either mult or perm)

    if(permbl.onesample)
    {
       do 
       {
	 if(perm !=1) permvec = permbl.next_perm(perm,exhaustive);
       }while (!exhaustive && !check_perm(oldperms,permvec));
       dm_mult(dm,dmperm,permvec);
    }
    else 
    {
      ColumnVector labels,oldlabels=permbl.buildfullperm();
      do
      {
	  if (perm !=1) labels = permbl.next_perm(perm,exhaustive);
	  else labels=oldlabels;       	          
       } while(!exhaustive && !check_perm(oldlabs,labels));
       LabelsToPerm(labels,oldlabels,permvec);  
       dm_permute(dm,dmperm,permvec);
       if (!exhaustive) oldlabs.push_back(labels);  
    }
    oldperms.push_back(permvec);  
    if(opts.var_sm_sig.value()==0) ols(datam,dmperm,tc,cope,varcope,dof);   
    else ols_var_sm(datam,dmperm,tc,cope,varcope,mask,mask_sm,opts.var_sm_sig.value(),dof);   
    tstat=SD(cope,sqrt(varcope));
    if (opts.tfce.value() && tstatnum>=0 )
    {
      if (perm==1) tfce_delta=tstat.Maximum()/100;  // i.e. 100 subdivisions of the max input stat height
      tstat_ce=tfce(tstat,mask,tfce_delta,opts.tfce_height.value(),opts.tfce_size.value(),opts.tfce_connectivity.value());
      if (perm==1) 
      {  
        tstat_ce_orig=tstat_ce;
	if ( opts.cluster_norm.value() ) 
	{
	  tstat_ceav=tstat_ce;
	  tstat_cenorm=SD(tstat_ce,tstat_ce);
	  if (!lowram)
	  {
	    try { oldtstatce=new Matrix(n_perms,tstat_ceav.Ncols());} //between 5e17 - 5e18 values for a 2gb machine
	    catch (...) {cerr << "using lowram" << endl; lowram=true;}         
	  }
	}
      }
      else if (opts.cluster_norm.value()) 
      {
	tstat_ceav+=tstat_ce;
	tstat_cenorm+=SD(tstat_ce,tstat_ce);
      }
      if(!lowram && opts.cluster_norm.value() ) oldtstatce->Row(perm)=tstat_ce.Row(1);
      maxdistCE.Column(perm)=max(tstat_ce.t()).t();
    }
    if ( tstatnum < 0 ) tstat=generate_fstat(dmperm,datam,NewW2,dof);
    if (perm==1) tstat_orig=tstat;

   maxdist.Column(perm)<<max(tstat.t()).t(); // max stat recording
   vox_numbigger+= geqt(tstat,tstat_orig); // voxelwise stats recording
   if (tstatnum>=0) threshold=opts.cluster_thresh.value();
   else threshold=opts.f_thresh.value();
   if (threshold>0) { //cluster thresholding
      tstat4D.setmatrix(tstat,mask);
      tstat4D.binarise(threshold);
      for(int t=0;t<tstat4D.tsize();t++)
      {
	clust_label=connected_components(tstat4D[t],clustersizes,26);
	if ( clustersizes.Nrows() > 0 ) maxdistC(t+1,perm)=int(clustersizes.MaximumAbsoluteValue());
        if ( opts.cluster_norm.value() && tstatnum>=0)
	{
	  for(int z=0; z<clust_label.zsize(); z++)
	    for(int y=0; y<clust_label.ysize(); y++)
	      for(int x=0; x<clust_label.xsize(); x++)
	      {
		if (clust_label(x,y,z)) 
		{
		  sumclusters(x,y,z,t)+=clustersizes(clust_label(x,y,z));
		  sumcluster_samples(x,y,z,t)++; 
		}
	      }
	} 

      }
   } //end of cluster tresholding
     
    if (tstatnum>=0) threshold=opts.clustermass_thresh.value();
    else threshold=opts.fmass_thresh.value();
    if ( threshold > 0 ) { //cluster mass thresholding
      tstat4D.setmatrix(tstat,mask);
      volume4D<float>tmp_tstat4D(tstat4D);
      tstat4D.binarise(threshold);
      for(int t=0;t<tstat4D.tsize();t++){
	clust_label=connected_components(tstat4D[t],clustersizes,26);
	clustersizes=0;	
	for(int z=0; z<clust_label.zsize(); z++)
	  for(int y=0; y<clust_label.ysize(); y++)
	    for(int x=0; x<clust_label.xsize(); x++)
	      if(clust_label(x,y,z)>0)
		clustersizes(clust_label(x,y,z))=clustersizes(clust_label(x,y,z))+tmp_tstat4D[t](x,y,z);
	if ( clustersizes.Nrows() > 0 ) maxdistCmass(t+1,perm)=int(clustersizes.MaximumAbsoluteValue());
      }
    }
  }

  //Create cluster norm map
 if ( opts.cluster_thresh.value()>0 && opts.cluster_norm.value() && tstatnum>=0) 
 {
   save_volume4D(sumclusters,opts.out_fileroot.value()+"_clusternormpre");
   save_volume4D(sumcluster_samples,opts.out_fileroot.value()+"_clusternormdiv");
   for(int t=0; t<sumclusters.tsize(); t++)
     for(int z=0; z<sumclusters.zsize(); z++)
       for(int y=0; y<sumclusters.ysize(); y++)
         for(int x=0; x<sumclusters.xsize(); x++)
           if(sumcluster_samples(x,y,z,t)>0)
             sumclusters(x,y,z,t)/=sumcluster_samples(x,y,z,t);
 }

 if (opts.tfce.value() && tstatnum>=0 )
 {
   OutputStat(tstat_ce_orig,maxdistCE,mask,n_perms,"_tfce_tstat",opts,tstatnumber);
   if (opts.cluster_norm.value())
   {
     volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
     output.setmatrix(tstat_ceav,mask);
     save_volume4D(output,opts.out_fileroot.value()+"_tfceavnum"+tstatnumber);
     tstat_ceav=SD(tstat_ceav,tstat_cenorm);
     output.setmatrix(tstat_ceav,mask);
     float min=output.percentile(0.02,mask);     
     for(int t=1;t<=tstat_ceav.Ncols();t++)
     {
       if(tstat_ceav(1,t)<min)
       {
	 //cerr << tstat_ceav(1,t) << " " << min << endl;
	 tstat_ceav(1,t)=min;
       }
     }
     output.setmatrix(tstat_ceav,mask);
     save_volume4D(output,opts.out_fileroot.value()+"_tfceavnorm"+tstatnumber);
     output.setmatrix(tstat_cenorm,mask);
     save_volume4D(output,opts.out_fileroot.value()+"_tfceavdenom"+tstatnumber);
   }
 }
    

   //Rerun perms for clusternorm
 if (opts.cluster_norm.value() && tstatnum>=0)
 { 
   for(int perm=1; perm<=n_perms; perm++)
   {
     if (opts.verbose.value()) cout << "starting second-pass permutation " << perm << endl;
     if ( opts.cluster_thresh.value()>0 || ( opts.tfce.value() && lowram ) ) //Regenerate stats
     {
       if(!permbl.onesample) dm_permute(dm,dmperm,oldperms.at(perm-1));
       else dm_mult(dm,dmperm,oldperms.at(perm-1));
       if(opts.var_sm_sig.value()==0) ols(datam,dmperm,tc,cope,varcope,dof);   
       else ols_var_sm(datam,dmperm,tc,cope,varcope,mask,mask_sm,opts.var_sm_sig.value(),dof);   
      tstat=SD(cope,sqrt(varcope));
     }
     if ( opts.tfce.value() )
     {
       if (!lowram) tstat_ce=oldtstatce->Row(perm);
       else tstat_ce=tfce(tstat,mask,tfce_delta,opts.tfce_height.value(),opts.tfce_size.value(),opts.tfce_connectivity.value());
       tstat_ce=SD(tstat_ce,tstat_ceav); 
       maxdistCE.Column(perm)<<max(tstat_ce.t()).t();
     }
     if ( opts.cluster_thresh.value()>0 )
     { 
       tstat4D.setmatrix(tstat,mask);
       tstat4D.binarise(opts.cluster_thresh.value());
       for(int t=0;t<tstat4D.tsize();t++)
       {	
         clust_label=connected_components(tstat4D[t],clustersizes,26);
         ColumnVector entries,cluster(clustersizes.Nrows());
         cluster=0;
         entries=cluster;
         for(int z=0; z<clust_label.zsize(); z++)
	   for(int y=0; y<clust_label.ysize(); y++)
	     for(int x=0; x<clust_label.xsize(); x++)
	       if (clust_label(x,y,z))
	       {
	         cluster(clust_label(x,y,z))+=sumclusters(x,y,z,t);
                 entries(clust_label(x,y,z))++;
	       }
         for(int i=1;i<=clustersizes.Nrows();i++) clustersizes(i)/=(cluster(i)/entries(i));
         if ( clustersizes.Nrows() > 0 ) maxdistCN(t+1,perm)=clustersizes.MaximumAbsoluteValue();
       }
     }
   }
 }
   

  //OUTPUT Routines 
  if (tstatnum>=0) OutputStat(tstat_orig,maxdist,mask,n_perms,"_tstat",opts,tstatnumber);
  if (tstatnum>=0 && opts.tfce.value() && opts.cluster_norm.value() ) tstat_ce_orig=SD(tstat_ce_orig,tstat_ceav); 
  if (tstatnum>=0 && opts.tfce.value() && opts.cluster_norm.value() ) OutputStat(tstat_ce_orig,maxdistCE,mask,n_perms,"_tfceav_tstat",opts,tstatnumber);
  if (tstatnum<0) OutputStat(tstat_orig,maxdist,mask,n_perms,"_fstat",opts,tstatnumber);

    volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
    for(int row=1; row<=n_distrows; row++){
      if (tstatnum==0) tstatnumber=num2str(row);
      Matrix this_t_stat=1-vox_numbigger.Row(row)/float(n_perms);
      output.setmatrix(this_t_stat,mask); 
      if (tstatnum>=0) save_volume4D(output,opts.out_fileroot.value()+"_vox_tstat"+tstatnumber);
      else save_volume4D(output,opts.out_fileroot.value()+"_vox_fstat"+tstatnumber);
    }

    if ( opts.cluster_thresh.value() > 0  && tstatnum >=0 ) OutputClusterStat(tstat_orig,maxdistC,mask,n_perms,"_maxc",opts,tstatnumber,opts.cluster_thresh.value(),sumclusters,"_tstat");
    if ( opts.f_thresh.value() > 0  && tstatnum <0 ) OutputClusterStat(tstat_orig,maxdistC,mask,n_perms,"_maxc",opts,tstatnumber,opts.cluster_thresh.value(),sumclusters,"_fstat");

  if ( opts.cluster_thresh.value() > 0 && opts.cluster_norm.value() && tstatnum >=0 ) 
  {
    OutputClusterStat(tstat_orig,maxdistCN,mask,n_perms,"_maxcn",opts,tstatnumber,opts.cluster_thresh.value(),sumclusters,"_tstat");
    save_volume4D(sumclusters,opts.out_fileroot.value()+"_clusternorm");
  }
  if ( opts.clustermass_thresh.value() > 0 && tstatnum >=0 )  OutputClusterStat(tstat_orig,maxdistCmass,mask,n_perms,"_maxcmass",opts,tstatnumber,opts.clustermass_thresh.value(),sumclusters,"_tstat");
  if ( opts.fmass_thresh.value() > 0 && tstatnum <0 )  OutputClusterStat(tstat_orig,maxdistCmass,mask,n_perms,"_maxcmass",opts,tstatnumber,opts.clustermass_thresh.value(),sumclusters,"_fstat");

  if (opts.cluster_norm.value() && opts.tfce.value() && !lowram) delete oldtstatce;
}

void Convertcontrast(const Matrix &input_x,const Matrix &input_c1,const Matrix &input_datam,Matrix &NewModel,Matrix &NewCon, Matrix &NewDataM, Matrix &NewW2 )
{  
    int r,p;
    Matrix c2;
    r = input_c1.Nrows();
    p = input_c1.Ncols();
    Matrix tmp=(Identity(p)-input_c1.t()*pinv(input_c1.t()));
    Matrix U,V;
    DiagonalMatrix D;
    SVD(tmp, D, U, V);
     c2=U.SubMatrix(1,U.Nrows(),1,p-r);
    c2=c2.t();
    Matrix C = input_c1 & c2;
    Matrix W=input_x*C.i();

    Matrix W1=W.SubMatrix(1,W.Nrows(),1,r);
    Matrix W2=W.SubMatrix(1,W.Nrows(),r+1,W.Ncols());
    NewModel=W1-W2*pinv(W2)*W1;
    NewCon=Identity(r);
    NewDataM=(Identity(W2.Nrows())-W2*pinv(W2))*input_datam;
    NewW2=W2;
}


void doit(Matrix &tc,Matrix &dm,Matrix &datam,volume<float> &mask, Matrix &gp,const int &num,ranopts& opts)
{
  //-ve num for f-stat contrast, 0 for old style "block t-stat"
  Matrix NewModel,NewCon,NewDataM,NewW2,tc2=tc;
  ColumnVector labels;
  permblock permbl;
  if (num<0) Convertcontrast(dm,tc,datam,NewModel,NewCon,NewDataM,NewW2);
  else if (num>0 && tc.Ncols()>1) Convertcontrast(dm,tc.SubMatrix(num,num,1,tc.Ncols()),datam,NewModel,NewCon,NewDataM,NewW2);
  else
  {
        if (num) NewCon=tc.SubMatrix(num,num,1,tc.Ncols());
        else  NewCon=tc;
        NewDataM=datam;
        NewModel=dm;
  }
  if (num>0) tc2=tc.Row(num);
  labels=make_dm_labels(dm*tc2.t());  //dumb regressor
  permbl.createpermblocks(dm*tc2.t(),gp); 
  int  nonzero=0;
  for (int i=1;i<=tc2.Ncols();i++) if (tc2(1,i)) nonzero++;
  permbl.initpermblocks((num>0 && nonzero==1 && labels.Sum() == labels.Nrows()),labels);
  if(permbl.onesample) cout << "One-sample design detected; sign-flipping instead of permuting." << endl;
  if(opts.verbose.value() || opts.how_many_perms.value()) 
  {
    if(!permbl.onesample) cout << permbl.exhaust_perms << " permutations required for exhaustive test";
    else cout << permbl.exhaust_perms << " sign-flipping permutations required for exhaustive test";
    if (num>0)  cout << " of t-test " << num << endl;
    if (num==0) cout << " of all t-tests " << endl;
    if (num<0)  cout << " of f-test " << abs(num) << endl;
    if(opts.how_many_perms.value()) return;
  }  

  float dof=ols_dof(dm); 
  if (opts.demean_data.value()) dof--;
  DoPermutationStatistics(opts,mask,NewDataM,NewCon,NewModel,num,NewW2,dof,permbl); 
}


void fparse(Matrix &fc,Matrix &tc,Matrix &model,Matrix &data,volume<float> &mask,Matrix &gp,ranopts& opts)
{   
    Matrix fstat_newcon(tc.Nrows(),tc.Ncols());
   for(int fstat=1; fstat<=fc.Nrows() ; fstat++ ) 
   {
      int rows=1;
      fstat_newcon.ReSize(tc.Nrows(),tc.Ncols());
      for (int tstat=1; tstat<=fc.Ncols() ; tstat++ )
	if (fc(fstat,tstat)==1) fstat_newcon.Row(rows++)=tc.Row(tstat);
      fstat_newcon=fstat_newcon.SubMatrix(1,rows-1,1,fstat_newcon.Ncols());
      doit(fstat_newcon,model,data,mask,gp,-fstat,opts);
   }
}

int main(int argc,char *argv[]){
  Log& logger = LogSingleton::getInstance();
  ranopts& opts = ranopts::getInstance();
  opts.parse_command_line(argc,argv,logger);
  Matrix model, tstat_con, fstat_con, data,grp_lab;
  volume<float> mask;
  if (Initialise(opts,mask,data,tstat_con,model,fstat_con,grp_lab)==-1) return -1;
  bool needsdemean=true;
  for (int i=1;i<=model.Ncols();i++) if ( fabs( (model.Column(i)).Sum() ) > 0.0001 ) needsdemean=false;
  if (needsdemean && !opts.demean_data.value()) cerr << "Warning: All design columns have zero mean - consider using the -D option to demean your data" << endl;
  if(opts.fc_file.value()!="") fparse(fstat_con,tstat_con,model,data,mask,grp_lab,opts); 
  for (int tstat=1; tstat<=tstat_con.Nrows() ; tstat++ )  doit(tstat_con,model,data,mask,grp_lab,tstat,opts); 
return 0;
}

void permblock::createpermblocks(const Matrix &dmr,const Matrix &groups)
{
  int current=0;
  num_blocks=int(groups.Maximum())+1;
  num_perms.ReSize((int)groups.Maximum());
  num_elements=dmr.Nrows();
  original_locations = new ColumnVector[num_blocks];
  permuted_locations = new ColumnVector[num_blocks];
  //for(int j=1;j<=dmr.Nrows();j++) if((dmr.Row(j)).Sum()==0) current++;
  //original_locations[0].ReSize(current);
  //permuted_locations[0].ReSize(current);
  //current=1;
  //for(int j=1;j<=dmr.Nrows();j++) if(groups(j,1)==0 || (dmr.Row(j)).Sum()==0) original_locations[0](current++)=j;

  for(int i=0;i<=groups.Maximum();i++)
  {
     current=0;
     for(int j=1;j<=dmr.Nrows();j++)
       if(groups(j,1)==i /*&& (dmr.Row(j)).Sum()!=0*/ ) current++;
     original_locations[i].ReSize(current);
     permuted_locations[i].ReSize(current);
     current=1;
     for(int j=1;j<=dmr.Nrows();j++)
       if(groups(j,1)==i /*&& (dmr.Row(j)).Sum()!=0*/ ) original_locations[i](current++)=j;
  }
}

double permblock::initpermblocks(bool onesample_input,ColumnVector labels)
{
   exhaust_perms=1;
   onesample=onesample_input;
   for(int i=0;i<num_blocks;i++)
   { 
     for(int j=1;j<=permuted_locations[i].Nrows();j++) 
     {
       if(!onesample) permuted_locations[i](j)=labels((int)original_locations[i](j));
       else permuted_locations[i](j)=1;
     }
     if (i>0) 
     { 
       if(!onesample) num_perms(i)=compute_nperms(permuted_locations[i]);
       else num_perms(i)=std::pow(2.0,labels.Nrows());
       exhaust_perms*=num_perms(i);
     }
   }
   return(exhaust_perms);
}

ColumnVector permblock::buildfullperm()
{
ColumnVector newvec(num_elements); 
   for(int i=0;i<num_blocks;i++)
     for(int j=1;j<=permuted_locations[i].Nrows();j++) 
       newvec((int)original_locations[i](j))=permuted_locations[i](j);
   return newvec;
}

ColumnVector permblock::next_perm(int perm,bool exhaustive)
{
  for(int i=1;i<num_blocks;i++)
  {
    long int interval = 1;
    for(int j=i+1;j<num_blocks;j++) interval*=(long int)num_perms(j);
    if ( (perm-1)%interval==0) 
    { 
      if(!onesample) next_perm_vec(permuted_locations[i],exhaustive);
      else next_mult_vec(permuted_locations[i],exhaustive);
    }
  }
  return(buildfullperm());
}

permblock::permblock()
{

}

permblock::~permblock()
{
  delete [] original_locations;
  delete [] permuted_locations;
}
