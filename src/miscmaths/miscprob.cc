/*  miscprob.cc

    Christian Beckmann & Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 3.3 (c) 2006, The University of
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

// Miscellaneous maths functions that rely on libprob

#include "miscprob.h"
#include "stdlib.h"
#include "newmatio.h"
#include <iostream>

using namespace NEWMAT;

namespace MISCMATHS {

ReturnMatrix unifrnd(const int dim1, const int dim2, const float start, const float end)
{
  int tdim = dim2;
  double tmpD=1.0;
  if(tdim<0){tdim=dim1;}
  Matrix res(dim1,tdim);
  for (int mc=1; mc<=res.Ncols(); mc++) {
    for (int mr=1; mr<=res.Nrows(); mr++) {
      tmpD = (rand()+1)/double(RAND_MAX+2.0);
      res(mr,mc)=(tmpD)*(end-start)+start;
      //drand(&tmpD);
      //res(mr,mc)=(tmpD-1)*(end-start)+start;
    }
  }

  res.Release();
  return res;
}

ReturnMatrix normrnd(const int dim1, const int dim2, const float mu, const float sigma)
{
  int tdim = dim2;
  double tmpD=1.0;
  if(tdim<0){tdim=dim1;}
  Matrix res(dim1,tdim);
  for (int mc=1; mc<=res.Ncols(); mc++) {
    for (int mr=1; mr<=res.Nrows(); mr++) {
      tmpD = (rand()+1)/double(RAND_MAX+2.0);
      res(mr,mc)=ndtri(tmpD)*sigma+mu ;
      //drand(&tmpD);
      //res(mr,mc)=ndtri(tmpD-1)*sigma+mu ;
    }
  }

  res.Release();

  return res;
}

ReturnMatrix normpdf(const RowVector& vals, const float mu, const float var)
{
  RowVector res(vals);
  for (int mc=1; mc<=res.Ncols(); mc++){
    res(mc) = std::exp(-0.5*(std::pow(vals(mc)-mu,2)/var))*std::pow(2*M_PI*var,-0.5);
  }

  res.Release();
  return res;
}
float normpdf(const ColumnVector& vals, const ColumnVector& mu, const SymmetricMatrix& sigma)
{
  float res;
  LogAndSign ld=(2*M_PI*sigma).LogDeterminant();

  res =  std::exp(-0.5*(
			((vals-mu).t()*sigma.i()*(vals-mu)).AsScalar()
			+ld.LogValue()
			));
  return res;
}
ReturnMatrix normpdf(const Matrix& vals, const ColumnVector& mu, const SymmetricMatrix& sigma)
{
  RowVector res(vals);
  LogAndSign ld = (2*M_PI*sigma).LogDeterminant();
  Matrix isigma = sigma.i();

  for (int mc=1; mc<=res.Ncols(); mc++){
    res(mc) = std::exp( -0.5*(
			      ((vals.Column(mc)-mu).t()*isigma*(vals.Column(mc)-mu)).AsScalar() 
			      +ld.LogValue()
			      ));
  }

  res.Release();
  return res;
}
ReturnMatrix normcdf(const RowVector& vals, const float mu, const float var)
{
  RowVector res(vals);
  RowVector tmp;
  tmp = (vals-mu)/std::sqrt(var);
  for (int mc=1; mc<=res.Ncols(); mc++){
    res(mc) = ndtr(tmp(mc));
  }

  res.Release();
  return res;
}

ReturnMatrix gammacdf(const RowVector& vals, const float mu, const float var)
{
  RowVector res(vals);
  res=0;
  if((mu>0)&&(var>0)){
    float b = std::pow(mu,2)/var;
    float a = mu/var;  
    for (int mc=1; mc<=res.Ncols(); mc++){
      if(vals(mc)>0)
	res(mc) = gdtr(a,b,vals(mc));
    }
  }
  res.Release();
  return res;
}

ReturnMatrix gammapdf(const RowVector& vals, const float mu, const float var)
{
  RowVector res(vals);
  res=0;
  if((mu>0)&&(var>0.00001)){
    float a = std::pow(mu,2)/var;
    float b = mu/var;
    float c = lgam(a);
    if(std::abs(c) < 150){
      for (int mc=1; mc<=res.Ncols(); mc++){
	if(vals(mc)>0.000001){
	  res(mc) = std::exp(a*std::log(b) + 
			     (a-1) * std::log(vals(mc)) 
			     - b*vals(mc) - c);
	}
      }
    }
  }
  res.Release();
  return res;
}


ReturnMatrix normpdf(const RowVector& vals, const RowVector& mu, const RowVector& var)
{
  Matrix res(mu.Ncols(),vals.Ncols());
  for (int mc=1; mc<=res.Ncols(); mc++){
    for (int mr=1; mr<=res.Nrows(); mr++){
      res(mr,mc) = std::exp(-0.5*(std::pow(vals(mc)-mu(mr),2)/var(mr)))*std::pow(2*M_PI*var(mr),-0.5);
    }
  }

  res.Release();
  return res;
}


ReturnMatrix mvnrnd(const RowVector& mu, const SymmetricMatrix& covar, int nsamp) 
{     
//   Matrix eig_vec; 
//   DiagonalMatrix eig_val;
//   EigenValues(covar,eig_val,eig_vec);

//   Matrix ret = ones(nsamp, 1)*mu + dnormrandm(nsamp,mu.Ncols())*sqrt(eig_val)*eig_vec.t();
  Mvnormrandm mvn(mu, covar);

  return mvn.next(nsamp);
}

  /*
// Saad: Wishart and inverseWishart Random Generator
ReturnMatrix wishrnd(const SymmetricMatrix& sigma,const int dof){
  // compute cholesky factor for sigma
  LowerTriangularMatrix L = Cholesky(sigma);

  // for small degrees of freedom, use the definition
  int n = sigma.Nrows();
  Matrix X;
  if(dof <= 81+n ){
    X.ReSize(dof,n);
    X = normrnd(dof,n) * L.t();
  }
  // otherwise, use Smith & Hocking procedure
  else{
    X.ReSize(n,n);
    Matrix A(n,n);
    for(int i=1;i<=n;i++){
      Gamma G((dof-i+1)/2);
      G.Set(rand()/float(RAND_MAX));
      for(int j=1;j<=n;j++){
	if     (i>j) { A(i,j) = 0; }
	else if(i<j) { A(i,j) = normrnd(1,1).AsScalar(); }
	else         { A(i,j) = std::sqrt(2*G.Next()); }
      }
    }
    X = A * L.t();
  }

  SymmetricMatrix res(n);
  res << X.t() * X;

  res.Release();
  return res;
}
ReturnMatrix iwishrnd(const SymmetricMatrix& sigma,const int dof){
  // assumes inv-Wishart(sigma.i(),dof)

  SymmetricMatrix res;
  res = wishrnd(sigma,dof);
  res = res.i();

  res.Release();
  return res;
}

  */

ReturnMatrix perms(const int n){
  if(n<=1){
    Matrix P(1,1);
    P << n;
    P.Release();
    return P;
  }
  Matrix Q = perms(n-1);  // recursive calls
  int m = Q.Nrows();
  Matrix P(n*m,n);
  for(int i=1;i<=m;i++){
    P(i,1)=n;
    for(int j=1;j<=Q.Ncols();j++)
      P(i,j+1)=Q(i,j);
  }
  for(int i=n-1;i>=1;i--){
    int jj=1;
    for(int j=(n-i)*m+1;j<=(n-i+1)*m;j++){
      P(j,1)=i;
      for(int k=1;k<=n-1;k++){
	P(j,k+1)= (Q(jj,k)==i) ? n : Q(jj,k);
      }
      jj++;
    } 
  }
  P.Release();
  return P;
}

}
