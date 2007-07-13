/*
 Brian Patenaude, parts taken from bet2
 FMRIB Image Analysis Group
 */


//#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "math.h"

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "meshclass/meshclass.h"
#include "shapeModel/shapeModel.h"
using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace mesh;
using namespace shapemodel;


string title="overlapMeasures (Version 1.0) University of Oxford (Brian Patenaude)";
string examples="overlapMeasure [options] -i segImage -l goldStandard -k output.txt ";


Option<bool> verbose(string("-v,--verbose"), false, 
					 string("switch on diagnostic messages"), 
					 false, no_argument);
Option<bool> help(string("-h,--help"), false,
				  string("display this message"),
				  false, no_argument);
Option<bool> inputprob(string("--inputprob"), false,
					   string("Input image is probability - do not do any thresholding or smoothing"),
					   false, no_argument);
Option<string> refname(string("-r,--in"), string(""),
					  string("filename of reference image "),
					  false, requires_argument);

Option<string> inname(string("-i,--in"), string(""),
					  string("filename of input image/mesh/bvars"),
					  true, requires_argument);
Option<string> inname2(string("-j,--in"), string(""),
					  string("filename of input image/mesh/bvars"),
					  false, requires_argument);
Option<string> designname(string("-d,--in"), string(""),
					  string("filename of input image/mesh/bvars"),
					  false, requires_argument);
//Option<string> designname2(string("-f,--in"), string(""),
//					  string("filename of input image/mesh/bvars"),
//					  false, requires_argument);
Option<string> modelname(string("-m,--in"), string(""),
					  string("filename of input model"),
					  false, requires_argument);
//Option<string> modelname2(string("-o,--in"), string(""),
//					  string("filename of input model"),
//					  false, requires_argument);

Option<bool> replaceBMV(string("--repbmv"), false,
					   string("replace the .bmv in .bvars file with -m input.bmv"),
					   false, no_argument);

					   
Option<bool> usebvars(string("--usebvars"), false,
					   string("use the barameters is rvm"),
					   false, no_argument);
Option<bool> useRVM(string("--useRVM"), false,
					  string("use relevance vector machien"),
					  false, no_argument);
Option<bool> bvarsnew(string("--bvarsnew"), false,
					  string("rigid reg and new bvars"),
					  false, no_argument);


Option<bool> modeStats(string("--modeStats"), false,
					   string("calculate t values for mode parameters"),
					   false, no_argument);

Option<string> outname(string("-k,--out"), string(""),
					   string("filename of output mesh"),
					   true, requires_argument);
Option<int> numPar(string("-n,--numberparams"),1,
					   string("number of iteratrion"),
					   false, requires_argument);
Option<int> kerntype(string("-a,--kerneltype"),0,
					   string("chan ge kernel"),
					   false, requires_argument);
Option<float> rad(string("-w,--kernel radius"),1.0,
					   string("kernalradius"),
					   false, requires_argument);
Option<float> errval(string("-e,--fractional error"),0.0,
					   string("err in fracrtiona;"),
					   false, requires_argument);
Option<int> polyorder(string("-p,--order of polynomial"),0,
					   string("polynomial order"),
					   false, requires_argument);



int nonoptarg;

////////////////////////////////////////////////////////////////////////////
//global variables

string read_bvars_design(string fname,Matrix* bvars, Matrix* design,vector<string>* vnames, vector<int>* vnvars){
	string stemp;
	string modelNames;
	int N;//number of subjects
		ifstream fin;
		fin.open(fname.c_str());
		//throw away first three lines 
		getline(fin,stemp);//this is bvars file
			getline(fin,modelNames);//modelnames
				fin>>stemp>>N;
				vnames->clear();
							vnvars->clear();
				
				//transform all the bvars
				for (int i=0; i<N;i++){
					fin>>stemp;//read in subject id
					vnames->push_back(stemp);
					
					
					int nvars;//how many vars written for the subject
						fin>>nvars;
						if (i==0){
							bvars->ReSize(nvars,N);
						}
						
							
						vnvars->push_back(nvars);
						//  cout<<"subject "<<i<<" "<<stemp<<" "<<nvars<<endl;
						for (int j=0;j<nvars;j++){
							if (j<nvars){
								float ftemp;
								fin>>ftemp;
								bvars->element(j,i)=ftemp;
							}else{
								bvars->element(j,i)=0;
							}
						}
				}
				//read design line
			int Nx;	
			fin>>stemp>>Nx;
			//cout<<stemp<<endl;
			cout<<"design matrix has "<<Nx<<" regressors"<<endl;
			getline(fin,stemp);
			//cout<<"clear line?"<<stemp<<endl;
			design->ReSize(N,Nx);
			for (int i=0; i<N;i++){
				for (int j=0;j<Nx;j++){
					float ftemp;
					fin>>ftemp;
				//	if (j==0){
				//		ftemp=1.0;
				//	}
					design->element(i,j)=ftemp;
					cout<<ftemp<<" ";
				}
				//cout<<endl;
			}
				
				
				cout<<"design read"<<endl;

				return modelNames;
}
string read_bvars_ModelName(string fname){
  string stemp;
  string modelNames;
   ifstream fin;
  fin.open(fname.c_str());
  //throw away first three lines 
  getline(fin,stemp);//this is bvars file
  getline(fin,modelNames);//modelnames
 
  return modelNames;
}
//************************RELVEVANCE VECTOR MACHINE FUNCTIONS ***************************//
Matrix distSqrd(Matrix* DataTest, Matrix* Data){
	int rows=Data->Nrows();
	int cols=Data->Ncols();
		int colsTest=DataTest->Ncols();

	Matrix MdistSqrd(colsTest,cols);
	//float r=2;
	//float rsq=r*r;
	
			for (int i=0;i<colsTest;i++){
				for (int j=0;j<cols;j++){
					//compute norm dif
					 MdistSqrd.element(i,j)=(DataTest->SubMatrix(1,rows,i+1,i+1)-Data->SubMatrix(1,rows,j+1,j+1)).SumSquare();
					//cout<<exp(-(Data.SubMatrix(1,rows,i+1,i+1)-Data.SubMatrix(1,rows,j+1,j+1)).SumSquare()/r/r)<<" ";
				}	
				//	cout<<endl;
			}
			return MdistSqrd;
}
Matrix compute_kernel(Matrix DataTest, Matrix Data,int kernel_type,  float r){
	int rows=Data.Nrows();
	int cols=Data.Ncols();
	int colsTest=DataTest.Ncols();
	cout<<"Data has "<<rows<<" "<<rows<<" and "<<cols<<" cols"<<endl;
	
	Matrix K(colsTest,cols);
	//float r=2;
	float rsq=r*r;
				ColumnVector Ones(colsTest);
				Matrix L2;
				float order=1;
				switch(kernel_type)
				{
					case 5://polynomia with bias
						order=static_cast<float>(polyorder.value());
						K=DataTest.t()*(Data)/rsq;
						for (int i=0;i<colsTest;i++){
							for (int j=0;j<cols;j++){
								K.element(i,j)+=1;
								K.element(i,j)=powf(K.element(i,j),order);
							}	
						}
							for (int i=0;i<colsTest;i++){
								Ones.element(i)=1;
							}
						K=Ones | K;
						break;		
					case 4:
						//linear with bias
						K.ReSize(Data.Ncols(),1+Data.Nrows());
						Ones.ReSize(Data.Ncols());
						K.SubMatrix(1,Data.Ncols(),2,Data.Nrows()+1)=Data.t();
						//	K=Data.t();
						cout<<"add bias"<<endl;
						for (int i=0;i<Data.Ncols();i++){
							Ones.element(i)=1;
						}
							K.SubMatrix(1,K.Nrows(),1,1)=Ones;
						//K=Ones | K;
						break;
						
					case 3:
						//linear with no bias
						K.ReSize(Data.Nrows(),Data.Ncols());
					//	K=Data.t();
						break; 
					case 2:
						//cubic with bias
						L2=distSqrd(&DataTest,&Data);
						for (int i=0;i<colsTest;i++){
							for (int j=0;j<cols;j++){
								//compute norm dif
								//K.element(i,j)=exp(-(DataTest.SubMatrix(1,rows,i+1,i+1)-Data.SubMatrix(1,rows,j+1,j+1)).SumSquare()/rsq);
								float temp=L2.element(i,j)/rsq;
								K.element(i,j)=temp*sqrt(temp);
								//cout<<exp(-(Data.SubMatrix(1,rows,i+1,i+1)-Data.SubMatrix(1,rows,j+1,j+1)).SumSquare()/r/r)<<" ";
							}	
						}
							for (int i=0;i<colsTest;i++){
								Ones.element(i)=1;
							}
							cout<<"cat oness "<<colsTest<<endl;
						K=Ones | K;
						cout<<"ones cat"<<endl;
						//cout<<K.Nrows();
						break;
					case 1://Gaussian Kernel withoout bias 
						cout<<"L2 to be  calculated"<<endl;
						L2=distSqrd(&DataTest,&Data);
						cout<<"L2 calculated "<<L2.Nrows()<<" "<<L2.Ncols()<<endl;
						for (int i=0;i<colsTest;i++){
							for (int j=0;j<cols;j++){
								//compute norm dif
								//K.element(i,j)=exp(-(DataTest.SubMatrix(1,rows,i+1,i+1)-Data.SubMatrix(1,rows,j+1,j+1)).SumSquare()/rsq);
								K.element(i,j)=exp(-L2.element(i,j)/rsq);
								//cout<<exp(-(Data.SubMatrix(1,rows,i+1,i+1)-Data.SubMatrix(1,rows,j+1,j+1)).SumSquare()/r/r)<<" ";
							}	
						}
							break;
					case 0://Gaussian Kernel with bias 
						L2=distSqrd(&DataTest,&Data);
						for (int i=0;i<colsTest;i++){
							for (int j=0;j<cols;j++){
								//compute norm dif
								//K.element(i,j)=exp(-(DataTest.SubMatrix(1,rows,i+1,i+1)-Data.SubMatrix(1,rows,j+1,j+1)).SumSquare()/rsq);
								K.element(i,j)=exp(-L2.element(i,j)/rsq);
								//cout<<exp(-(Data.SubMatrix(1,rows,i+1,i+1)-Data.SubMatrix(1,rows,j+1,j+1)).SumSquare()/r/r)<<" ";
							}	
						}
							for (int i=0;i<colsTest;i++){
								Ones.element(i)=1;
							}
							cout<<"cat oness "<<colsTest<<endl;
						K=Ones | K;
						cout<<"ones cat"<<endl;
						break;
					default: cerr<<"not valid kernel type"<<endl;
						
				}
				
				bool bias =false;
				if (bias){//add bias to K
					Ones.ReSize(K.Nrows());
					K.SubMatrix(1,Data.Ncols(),2,Data.Nrows()+1)=Data.t();
					for (int i=0;i<Data.Ncols();i++){
						Ones.element(i)=1;
					}
					K.SubMatrix(1,Data.Nrows(),1,1)=Ones;
				}
				
				return K;
}


ColumnVector sigmoid(ColumnVector x){
	//ColumnVector y=x;
	//cout<<"sigmoid "<<endl;
	for (int i=0;i<x.Nrows();i++){
		x.element(i)=1/(1+exp(-x.element(i)));
		//x.element(i)=exp(-x.element(i));
		//cout<<x.element(i)<<endl;
	}
	
	//float y=1/(1+exp(-x));
 return x;
}
ColumnVector calculatePosteriorMode(int maxit, DiagonalMatrix A, Matrix Phi, ColumnVector w, Matrix t, Matrix* Hes){
	float stop_crit_sq=1e-6;
	float lambda_min=2e-8;
	//cout<<"enter post mode finder "<<w.Nrows()<<" "<<Phi.Ncols()<<endl;
	
	//this is required when dealing with classification
	
	ColumnVector Phiw(w.Nrows());
	Phiw=Phi*w;
	///cout<<"Phiz"<<endl;
	int N=Phi.Nrows();
	int M=w.Nrows();
	int d=Phi.Ncols();
	ColumnVector y(M);
	y=sigmoid(Phiw);
	//this is using the negative log posterior
//	cout<<"calc err stuff"<<endl;
	float data_term=0;
	
	
	float reg_term=0;
	
	for (int i=0;i<N;i++){
		//cout<<i<<" "<<N<<" "<<M<<" "<<y.Nrows()<<" "<<t.Nrows()<<" "<<t.Ncols()<<endl;
		if (t.element(i,0)==1){
			data_term+=log(y.element(i));
		}else{
			data_term+=log(1-y.element(i));
		}
		//	cout<<i<<" "<<N<<" "<<M<<" "<<y.Nrows()<<" "<<t.Nrows()<<" "<<t.Ncols()<<endl;
	}	
		data_term*=-1;
	data_term/=N;
	//cout<<"regterm"<<endl;
	for (int i=0;i<M;i++){
//	cout<<i<<" "<<N<<" "<<M<<" "<<y.Nrows()<<" "<<t.Nrows()<<" "<<t.Ncols()<<endl;
//		cout<<"where doe sit break "<<i<<" "<<A.Nrows()<<" "<<A.Ncols()<<" "<<w.Nrows()<<endl;
		reg_term+=A.element(i)*w.element(i)*w.element(i);
	}
//cout<<"finish first loop"<<endl;
	reg_term/=(2.0*N);

	float err_new=data_term+reg_term;
	//cout<<"err_new calced"<<endl;
	//cout<<"err_new "<<err_new<<endl;
	ColumnVector errs(maxit);
	for (int i=0;i<maxit;i++){
		errs.element(i)=0;
	}
	
	ColumnVector vary(y.Nrows());
	vary=y*0;
	Matrix phiv=Phi;
	ColumnVector e;
	ColumnVector g;
	Matrix Hessian;
	ColumnVector delta_w;
	ColumnVector w_new;
//	cout<<"begin ]ations "<<endl;
	for (int iter=0;iter<maxit;iter++){
		
		//update vary
	//	cout<<"vary"<<endl;
		for (int i=0;i<N;i++){
			vary.element(i)=y.element(i)*(1-y.element(i));
			//cout<<vary.element(i)<<endl;
		}
		//cout<<"calc vary"<<endl;
		//update phiv
//		cout<<"phiV"<<endl;
//		for (int i=0;i<N;i++){
//			for (int j=0;j<d;j++){
//				cout<<Phi.element(i,j)<<" ";
//			}
//			cout<<endl;
//		}
		for (int i=0;i<Phi.Nrows();i++){
			for (int j=0;j<Phi.Ncols();j++){
				phiv.element(i,j)=Phi.element(i,j)*vary.element(i);
			//	cout<<phiv.element(i,j)<<" ";
			}
			//cout<<endl;
		}
//		Matrix pvp(5,5);
//		pvp.ReSize(5,5);
//		pvp=(phiv.t())*Phi;
//		cout<<"pvp"<<endl;
//		for (int i=0;i<d;i++){
//			for (int j=0;j<d;j++){
//				//phiv=Phi.element(i,j)*vary.element(i);
//				cout<<pvp.element(i,j)<<" ";
//			}
//			cout<<endl;
//		}
	//	cout<<"calculate phiv"<<endl;
		//update e
		e=t-y;
//		cout<<"e"<<endl;
//		for (int i=0;i<N;i++){
//			
//									cout<<e.element(i)<<endl;
//	}
//			cout<<"A"<<endl;
//		for (int i=0;i<M;i++){
//			
///			cout<<A.element(i)<<endl;
	//	}
//		cout<<"w"<<endl;
//		for (int i=0;i<M;i++){
			
//			cout<<w.element(i)<<endl;
//		}
		g= Phi.t() * e - A*w;
//		cout<<"G"<<endl;
//		for (int i=0;i<M;i++){
			
//			cout<<g.element(i)<<endl;
//		}
		
		Hessian= ( phiv.t() * Phi + A );
		//cout<<"where's the problem"<<endl;
	//	Matrix U;
	//	DiagonalMatrix D;
	//	SVD(Hessian,D,U);
//		for (int i =0; i<Hessian.Nrows();i++){
//		 cout<<D.element(i)<<" ";
//		 if (D.element(i)<1e-10){
//			return w;
//		 }
//		}
//		cout<<endl;
		//have not implement check for well conditioned-ness
		//if ill-conditioned should return and update alpha
		//continue
		//tipping now uses a cholesky decompostion
		//i am not using this 
		errs.element(iter)=err_new;
		//cout<<"ERR "<<err_new<<endl;
	//	cout<<"Stop crit "<<g.SumSquare()/w.Nrows()/w.Nrows()<<" "<<g.SumSquare()<<endl;
		if ((iter>=2)&&(g.SumSquare()/w.Nrows()/w.Nrows())<stop_crit_sq){
			errs=errs.SubMatrix(1,iter,1,1);
			//cout<<"reached stop criterian after "<<iter+1<<" iterations"<<endl;
			break;
		}
		delta_w= (Hessian.t()*Hessian).i()*Hessian.t()*g ;///need to check this ???

	//	cout<<"delta_w"<<endl;
	//	for (int i=0;i<M;i++){
	//		cout<<delta_w.element(i)<<" ";
	//	}
	//	cout<<endl;
		float lambda=1;
		cout<<"reached lambda"<<endl;
//		cout<<"lambda "<<lambda<<" "<<lambda_min<<endl;
		while (lambda>lambda_min){
		//cout<<"lambda "<<lambda<<endl;
			w_new=w+lambda*delta_w;
			Phiw=Phi*w_new;
			y=sigmoid(Phiw);
			data_term=0;
			reg_term=0;
		//	cout<<"update terms"<<endl;
			bool errinf=false;
			for (int i=0;i<N;i++){
				if (((y.element(i)==1)&&(t.element(i,0)==1))|((y.element(i)==0)&&(t.element(i,0)==0))){
					err_new=1e17;
					errinf=true;
				 //cout<<"err_new equal to infintie"<<endl;
					break;
				}else{
					if (t.element(i,0)==1){
						//cout<<"t "<<i<<endl;
						data_term+=log(y.element(i));
					}else{
						data_term+=log(1-y.element(i));
					}
					
					if (i==(N-1)){
						data_term*=-1;
						data_term/=N;
					}
				}
			}
			if (!errinf){
				//cout<<"not errinf"<<N<<" "<<M<<endl;
				for (int i=0;i<M;i++){
					reg_term+=A.element(i)*w_new.element(i)*w_new.element(i);
					if (i==(M-1)){
						reg_term/=(2.0*N);
						err_new=data_term+reg_term;
						//cout<<"update errnew"<<endl;
					}
				}
			}
			//cout<<"lambda "<<lambda<<endl;
			//cout<<"err_new "<<err_new<<" "<<errs.element(iter)<<endl;
			if (err_new>errs.element(iter)){
				lambda/=2;
			}else{
				//cout<<" set w_new"<<endl;
				w=w_new;
				lambda=0;
			}
	
		
		}	
	//g"left lambda"<<endl;
	
	}
	*Hes=Hessian.i();
	cout<<"end posterior mode"<<endl;
	return w;
}


ColumnVector run_rvm(Matrix targ, Matrix des, int niter, Matrix* Tcov, ColumnVector* indsout){
cout<<"rum rvm"<<endl;
	//load target into a matrix
//	Matrix targ(1,vtarg.size());
//	for (unsigned int i=0;i<vtarg.size();i++){
//		targ.element(0,i)=vtarg.at(i);
//	}
//	cout<<"target loaded"<<endl;
	//define parameter
	int N=static_cast<int>(des.Nrows());
	int M=des.Ncols();
	cout<<"Phi size "<<des.Nrows()<<" "<<des.Ncols()<<endl;
	//parameter initialization
	//alpha matrix
	DiagonalMatrix A(M);
	for (int i=0;i<M;i++){
		A.element(i)=1;
	}
	cout<<"alpha initialized"<<endl;
	//alpha matrix
	ColumnVector gam(M);
	for (int i=0;i<M;i++){
		gam.element(i)=1;
	}
	cout<<"gamma initialized"<<endl;
	//noise variance
	float nvar=1;
	
	
	//these keep track of where the indices actually correspond to 
	ColumnVector inds(M);
	for (int i=0;i<M;i++){
		inds.element(i)=i;
	}
	
	
	//initial weigth estimates
	
	//weight covriance
	Matrix cov(M,M);
	Tcov->ReSize(M,M);
	//weight mean
	ColumnVector muw(M);
		for (int j=0;j<(M);j++){
		muw.element(j)=0;
					}


	
		int Nrv=M;//number of relevance vectors
			  //begin iterations
			  	cout<<"Nrv "<<Nrv<<endl;
		for (int j=0;j<Nrv;j++){
		cout<<gam.element(j)<<" ";
		}
		cout<<endl;
		for (int j=0;j<Nrv;j++){
		cout<<A.element(j)<<" ";
		}
		cout<<endl;
		
			//computee mean and covariance
	//		cov= ((des.t()*des)/nvar + A).i();
			//cout<<"cov estimated"<<endl;
	//		muw=cov *des.t() * targ.t() / nvar;
			//cout<<"cov and mu estimated"<<endl;
		
			  cout<<"start iterations"<<endl;
		for (int i=0;i<niter;i++){
		//cout<<"iter "<<i<<endl;
			//prune any zero gammas
			//for (int j=0;j<Nrv;j++){
//				cout<<"alp ";
//			for (int j=0;j<Nrv;j++){
//		cout<<A.element(j)<<" ";
//		}
//		cout<<endl;
//	cout<<"gam ";
//			for (int j=0;j<Nrv;j++){
//		cout<<gam.element(j)<<" ";
//		}
//		cout<<endl;
//			cout<<"gam pre prune "<<endl;
//		for (int j=0;j<Nrv;j++){
//		cout<<gam.element(j)<<" ";
//		}


	
			
			int j=0;
			while (j<Nrv){
				//cout<<"submatric prob"<<endl;
				
				//this loop does the pruning of the irrelevant vectors
				//if (gam.element(j)<0.000000000001){
				if (A.element(j)>1e12){
					//	cout<<"prune "<<j<<endl;
					if (j==Nrv){
				//		cout<<"submatric probA"<<endl;
						Matrix dum=A.SubMatrix(1,Nrv-1,1,Nrv-1);
						A.ReSize(Nrv-1);
						A<<dum;
						
						//			dum=targ.SubMatrix(1,1,1,Nrv-1);
						//			targ.ReSize(1,Nrv-1);
						//			targ<<dum;
						
						dum=gam.SubMatrix(1,Nrv-1,1,1);
						gam.ReSize(Nrv-1);
						gam<<dum;
						
						dum=muw.SubMatrix(1,Nrv-1,1,1);
						muw.ReSize(Nrv-1);
						muw<<dum;
						
						dum=inds.SubMatrix(1,Nrv-1,1,1);
						inds.ReSize(Nrv-1);
						inds<<dum;
						
						dum=des.SubMatrix(1,N,1,Nrv-1);
						des.ReSize(N,Nrv-1);
						des<<dum;
						
						Nrv-=1;
						j--;
					}else if (j>0){
						// j+1 -1
				//		cout<<"submatric probB"<<endl;
						//prun alphas
						Matrix dum=A.SubMatrix(1,j,1,j) | A.SubMatrix(1,j,j+2,Nrv);
						dum=dum & ( A.SubMatrix(j+2,Nrv,1,j) | A.SubMatrix(j+2,Nrv,j+2,Nrv) );
						A.ReSize(Nrv-1);
						A<<dum;
						
						//prune gam
						//			cout<<targ.Nrows()<<" "<<Nrv<<" "<<targ.Ncols()<<endl;
						//			dum=targ.SubMatrix(1,1,1,j);
						//			dum=dum & targ.SubMatrix(1,1,j+2,Nrv);
						//			targ.ReSize(1,Nrv-1);
						//			targ<<dum;
						
						
						//prune gam
						dum=gam.SubMatrix(1,j,1,1);
						dum=dum & gam.SubMatrix(j+2,Nrv,1,1);
						gam.ReSize(Nrv-1);
						gam<<dum;
						
						dum=muw.SubMatrix(1,j,1,1);
						dum=dum & muw.SubMatrix(j+2,Nrv,1,1);
						muw.ReSize(Nrv-1);
						muw<<dum;
						
						
						dum=inds.SubMatrix(1,j,1,1);
						dum=dum & inds.SubMatrix(j+2,Nrv,1,1);
						inds.ReSize(Nrv-1);
						inds<<dum;
						
						
						//prune design matrix
						dum=des.SubMatrix(1,N,1,j) | des.SubMatrix(1,N,j+2,Nrv);
						des.ReSize(N,Nrv-1);
						des<<dum;
						
						
						
						
						Nrv-=1;
						j--;
					}else{
				//		cout<<"submatric probC"<<endl;
						Matrix dum=A.SubMatrix(2,Nrv,2,Nrv);
						A.ReSize(Nrv-1);
						A<<dum;
					//	cout<<targ.Nrows()<<" "<<Nrv<<" "<<targ.Ncols()<<endl;
						
						//			dum=targ.SubMatrix(1,1,2,Nrv);
						//			targ.ReSize(1,Nrv-1);
						//			targ<<dum;
						
						dum=gam.SubMatrix(2,Nrv,1,1);
						gam.ReSize(Nrv-1);
						gam<<dum;
						
						dum=muw.SubMatrix(2,Nrv,1,1);
						muw.ReSize(Nrv-1);
						muw<<dum;
						
						
						dum=inds.SubMatrix(2,Nrv,1,1);
						inds.ReSize(Nrv-1);
						inds<<dum;
						
						
						dum=des.SubMatrix(1,N,2,Nrv);
						des.ReSize(N,Nrv-1);
						des<<dum;
						
						Nrv-=1;
						j--;
					}
				}
				j++;
				}
			//cout<<"submatric prob2"<<endl;
			//						cout<<"gam pruned "<<endl;
			//		for (int j=0;j<Nrv;j++){
			//		cout<<gam.element(j)<<" ";
			//		}
			//					cout<<"A pruned "<<endl;
			//		for (int j=0;j<Nrv;j++){
			//		cout<<A.element(j)<<" ";
			//		}
			//	cout<<"done pruning"<<endl;
			
			//for classification
					bool classif=true;
			if (classif){
				cout<<"calc postmode"<<endl;
	//			for (int j=0;j<Nrv;j++){
	//	cout<<muw.element(j)<<" ";
	//	}
	//	cout<<endl;
		
				muw=calculatePosteriorMode(25, A, des, muw, targ.t(), &cov);
				//	calculatePosteriorMode(50, A, des, muw, targ.t(), &cov);
				//		cout<<"done postmode "<<endl;
			}else{
				//this is only used for regression case
				//computee mean and covariance
				cov= ((des.t()*des)/nvar + A).i();
				//cout<<"cov estimated"<<endl;
				muw=cov *des.t() * targ.t() / nvar;
				//cout<<"cov and mu estimated"<<endl;
				//update gammmas and alphas
			}
		//	cout<<"cov begin"<<endl;
		//	for (int j=0;j<cov.Nrows();j++){
		//		cout<<cov.element(j,j)<<endl;
			//				}
		//	cout<<"cov end"<<endl;
			float sumgam=0;
			for (int j=0;j<Nrv;j++){
					
		//			if (j==(0)){
		//						cout<<"gam "<<gam.element(j)<<" "<<muw.element(j)<<" "<<A.element(j)<<" "<<cov.element(j,j)<<endl;
		//					}
				gam.element(j)=1-A.element(j)*(cov.element(j,j));
				cout<<gam.element(j)<<endl;
		//					if (j==(0)){
		//						cout<<"gam "<<gam.element(j)<<" "<<muw.element(j)<<" "<<A.element(j)<<" "<<cov.element(j,j)<<endl;
		//					}
				sumgam+=gam.element(j);
				//cout<<"i "<<i<<endl;
				if ((i+1)<(0.5*niter)){
				A.element(j)=gam.element(j)/(muw.element(j)*muw.element(j));
				}else{
				
		//		cout<<gam.element(j)/(muw.element(j)*muw.element(j)/gam.element(j)-cov.element(j,j))<<" "<<gam.element(j)/(muw.element(j)*muw.element(j))<<endl;
				A.element(j)=gam.element(j)/(muw.element(j)*muw.element(j)/gam.element(j)-cov.element(j,j));
				//	cout<<"alpha "<<endl;
			
				}
			}
			for (int j=0;j<Nrv;j++){
				if (A.element(j)<1e-13){
					A.element(j)=1e20;
				}
		//	cout<<A.element(j)<<" ";
			}
		//	cout<<endl;
		//	}
			//update noise variance
		//	nvar=(targ.t()-des*muw).SumSquare()/(N-sumgam);
			//cout<<"NRV "<<Nrv<<endl;
			
	//					cout<<"gam "<<endl;
	//		for (int j=0;j<Nrv;j++){
	//			cout<<gam.element(j)<<" ";
	//		}
	//		cout<<endl;
//			
	//		cout<<"muw "<<endl;
	//		for (int j=0;j<Nrv;j++){
	//			cout<<muw.element(j)<<" ";
	//		}
	//		cout<<endl;
	//				
		}
	//	//output gammaa and alpha
		cout<<"relevance "<<endl;
		for (int j=0;j<Nrv;j++){
		cout<<gam.element(j)<<" ";
		}
		cout<<endl;
		cout<<"alpha "<<endl;
		for (int j=0;j<Nrv;j++){
		cout<<A.element(j)<<" ";
		}
		cout<<endl;
		cout<<"mu"<<endl;
		for (int j=0;j<Nrv;j++){
		cout<<muw.element(j)<<" ";
		}
		cout<<endl;
		cout<<"noise variance "<<nvar<<endl;
		cout<<"number of relevance vectors "<<Nrv<<endl;
		*Tcov=cov;
		*indsout=inds;
	return muw;
}
//*******************END OF RVM FUNCTION******************************//
int QDA(Matrix Data, Matrix target, ColumnVector test, int numModes){
//cout<<"enter QDA"<<endl;
//	Matrix DataIn=Data.t();
	//subject per column
	//find number of classes
	//assume they range from 0 to N
	vector<float> desc;
	
	ColumnVector testSub=test.SubMatrix(1,numModes,1,1);
	
	int N=1;//at least one class
	for (int i=0;i< target.Nrows();i++){
		if (target.element(i,0)>=N){ N++;}
	}
//	cout<<"there are "<<N<<" classes"<<endl;
	//now split up the data matrix
	//store all the covariance matrices
	vector< SymmetricMatrix > Vcovk;
	vector< ColumnVector > Vmu;
	
	vector<int> nk;
	for (int cl=0;cl<N;cl++){
		nk.push_back(0);//number of each class
		//Vmu.push_back(0);//mean of each class
	}
//	cout<<Data.Nrows()<<" "<<Data.Ncols()<<endl;
//	cout<<target.Nrows()<<" "<<target.Ncols()<<endl;
	for (int cl=0;cl<N;cl++){
	
		Matrix covk;
		ColumnVector mu(numModes);
		for (int i=1;i<= Data.Ncols();i++){
		//	cout<<"i "<<i<<endl;
			if (target.element(i-1,0)==cl){ 
				//this handles the covariance matrix
				if (nk.at(cl)==0){
					covk=Data.SubMatrix(1,numModes,i,i);
					mu=Data.SubMatrix(1,numModes,i,i);
					nk.at(cl)++;
				}else{
					covk=covk | Data.SubMatrix(1,numModes,i,i);
					mu=mu + Data.SubMatrix(1,numModes,i,i);
					nk.at(cl)++;
				}
			}
		//	cout<<"end i"<<endl;
			
		}
		SymmetricMatrix A(numModes);
		A << (covk*(covk.t())/ Data.Ncols());
		Vcovk.push_back(A);
		Vmu.push_back(mu);
	//	cout<<covk.Nrows()<<" "<<covk.Ncols()<<endl;
	}
	for (int cl=0;cl<N;cl++){
		
		Vmu.at(cl)=Vmu.at(cl)/nk.at(cl);//mean of each class
//		cout<<"class "<<cl<<" mean"<<endl;
//		for (int i=0;i< Vmu.at(cl).Nrows();i++){
//			cout<<Vmu.at(cl).element(i)<<" ";
//		}
//		cout<<endl;
	}
	
	float totalN=0;
	for (int cl=0;cl<N;cl++){
		totalN+=nk.at(cl);
	}
	///CALCULATE DESCRIMINANT, diagonlize for computation
	int Dclass;
	float decisionOld;
	for (int cl=0;cl<N;cl++){
	//	cout<<"diagonalize for class "<<cl<<endl;
		DiagonalMatrix D;
		Matrix V;
		//SymmetricMatrix A=Vcovk.at(cl);
		Jacobi(Vcovk.at(cl),D,V);
		Matrix Ux=V.t()*(testSub-Vmu.at(cl));
		float descMal=(Ux.t()*(D.i())*Ux).AsScalar();
		Matrix ed=(Ux.t()*(D.i())*Ux);
	//	cout<<descMal<<" "<<ed.Nrows()<<" "<<ed.Ncols()<<endl;
		float descDet=0;
		for (int i =0;i<D.Nrows();i++){
			//don't use negtive eigenvalues
			if (D.element(i)>0){
			descDet+=log(D.element(i));
			}
			//cout<<"D "<<i<<" "<<D.element(i)<<endl;
		}
		float descPi=log(static_cast<float>(nk.at(cl))/totalN);
		float decision=-0.5*(descDet + descMal) + descPi;
//		cout<<"decision "<<descDet<<" "<<descPi<<" "<<descMal<<" "<<decision<<endl;
//		cout<<"decision "<<decision<<endl;

		if (cl==0){
			Dclass=0;
			decisionOld=decision;
		}else if (decision>decisionOld){
			Dclass=cl;
			decisionOld=decision;
		}
		
		
	}
	
	
	return Dclass;
}

int QDAsphere(Matrix Data, Matrix target, ColumnVector test, int numModes){
//cout<<"enter QDA"<<endl;
//	Matrix DataIn=Data.t();
	//subject per column
	//find number of classes
	//assume they range from 0 to N
	vector<float> desc;
	
	ColumnVector testSub=test.SubMatrix(1,numModes,1,1);
	
	int N=1;//at least one class
	for (int i=0;i< target.Nrows();i++){
		if (target.element(i,0)>=N){ N++;}
	}
//	cout<<"there are "<<N<<" classes"<<endl;
	//now split up the data matrix
	//store all the covariance matrices
	vector< SymmetricMatrix > Vcovk;
	vector< ColumnVector > Vmu;
	
	vector<int> nk;
	for (int cl=0;cl<N;cl++){
		nk.push_back(0);//number of each class
		//Vmu.push_back(0);//mean of each class
	}
//	cout<<Data.Nrows()<<" "<<Data.Ncols()<<endl;
//	cout<<target.Nrows()<<" "<<target.Ncols()<<endl;
	for (int cl=0;cl<N;cl++){
	
		Matrix covk;
		ColumnVector mu(numModes);
		for (int i=1;i<= Data.Ncols();i++){
		//	cout<<"i "<<i<<endl;
			if (target.element(i-1,0)==cl){ 
				//this handles the covariance matrix
				if (nk.at(cl)==0){
					covk=Data.SubMatrix(1,numModes,i,i);
					mu=Data.SubMatrix(1,numModes,i,i);
					nk.at(cl)++;
				}else{
					covk=covk | Data.SubMatrix(1,numModes,i,i);
					mu=mu + Data.SubMatrix(1,numModes,i,i);
					nk.at(cl)++;
				}
			}
		//	cout<<"end i"<<endl;
			
		}
		SymmetricMatrix A(numModes);
		A << (covk*(covk.t())/ Data.Ncols());
		Vcovk.push_back(A);
		Vmu.push_back(mu);
	//	cout<<covk.Nrows()<<" "<<covk.Ncols()<<endl;
	}
	for (int cl=0;cl<N;cl++){
		
		Vmu.at(cl)=Vmu.at(cl)/nk.at(cl);//mean of each class
//		cout<<"class "<<cl<<" mean"<<endl;
//		for (int i=0;i< Vmu.at(cl).Nrows();i++){
//			cout<<Vmu.at(cl).element(i)<<" ";
//		}
//		cout<<endl;
	}
	
	float totalN=0;
	for (int cl=0;cl<N;cl++){
		totalN+=nk.at(cl);
	}
	///CALCULATE DESCRIMINANT, diagonlize for computation
	int Dclass;
	float decisionOld;
	for (int cl=0;cl<N;cl++){
	//	cout<<"diagonalize for class "<<cl<<endl;
		DiagonalMatrix D;
		Matrix V;
		//SymmetricMatrix A=Vcovk.at(cl);
		Jacobi(Vcovk.at(cl),D,V);
		Matrix Ux=V.t()*(testSub-Vmu.at(cl));
		float descMal=(Ux.t()*(D.i())*Ux).AsScalar();
		Matrix ed=(Ux.t()*(D.i())*Ux);
	//	cout<<descMal<<" "<<ed.Nrows()<<" "<<ed.Ncols()<<endl;
		float descDet=0;
		for (int i =0;i<D.Nrows();i++){
			//don't use negtive eigenvalues
			if (D.element(i)>0){
			descDet+=log(D.element(i));
			}
			//cout<<"D "<<i<<" "<<D.element(i)<<endl;
		}
		float descPi=log(static_cast<float>(nk.at(cl))/totalN);
		float decision=-0.5*(descDet + descMal) + descPi;
//		cout<<"decision "<<descDet<<" "<<descPi<<" "<<descMal<<" "<<decision<<endl;
//		cout<<"decision "<<decision<<endl;

		if (cl==0){
			Dclass=0;
			decisionOld=decision;
		}else if (decision>decisionOld){
			Dclass=cl;
			decisionOld=decision;
		}
		
		
	}
	
	
	return Dclass;
}



Matrix recon_meshesMNI(string refname, string modelname, Matrix bvars, ColumnVector* meanm, Mesh * meshout){
		volume<float> ref;
		read_volume(ref,refname);
	
	
	
		shapeModel* model1= new shapeModel();
		model1->setImageParameters(ref.xsize(), ref.ysize(), ref.zsize(),ref.xdim(), ref.ydim(), ref.zdim());
		cout<<"load model "<<modelname<<endl;
		model1->load_bmv_binaryInfo(modelname,1);
		model1->load_bmv_binary(modelname,1);

			//want to return mean mesh in vector form
		meanm->ReSize(model1->getTotalNumberOfPoints()*3);
		{
		
			int count=0;
			int cumnum=0; 
			for (int sh=0; sh<model1->getNumberOfShapes();sh++){
				Mesh m=model1->getTranslatedMesh(sh);
				*meshout=m;
				for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
					meanm->element(3*cumnum+count)=(*i)->get_coord().X;
					meanm->element(3*cumnum+count+1)=(*i)->get_coord().Y;
					meanm->element(3*cumnum+count+2)=(*i)->get_coord().Z;
					count+=3;
				}
				cumnum+=model1->getNumberOfPoints(sh);	
			}
			
		}
		
		//need number of subjects (to detrmine number of modes)
		int M=model1->getNumberOfSubjects();
		int Tpts=model1->getTotalNumberOfPoints();
			Matrix MeshVerts(3*Tpts,bvars.Ncols());
		//this is different than mnumber of subjects which were used to create model
		//int numSubs=bvars.Ncols();
		//this loads all mesh point into a matrix to do stats on....
		cout<<"generate and load vertices into matrix "<<endl;
		for (int j=0;j<bvars.Ncols();j++){
			//for each subject
			vector<float> vars;
			for (int i=0; i<bvars.Nrows();i++){
				vars.push_back(bvars.element(i,j));
			}
			//keep track of number of points preceding
			int cumnum=0;
			for (int sh=0; sh<model1->getNumberOfShapes();sh++){
				//cout<<"cumnum "<<cumnum<<endl;
				Mesh m=model1->getDeformedMesh(vars,sh,static_cast<int>(vars.size()));	
				
				int count=0;
				for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
					MeshVerts.element(3*cumnum+count,j)=(*i)->get_coord().X;
					MeshVerts.element(3*cumnum+count+1,j)=(*i)->get_coord().Y;
					MeshVerts.element(3*cumnum+count+2,j)=(*i)->get_coord().Z;
					count+=3;
				}	
				//	cout<<"count "<<count<<endl;
				cumnum+=model1->getNumberOfPoints(sh);	
			}	
		}
	return MeshVerts;

}
Matrix recon_meshesNative(string refname, string modelname, Matrix bvars, ColumnVector* meanm, Mesh * meshout){
		volume<float> ref;
		read_volume(ref,refname);
	
	
	
		shapeModel* model1= new shapeModel();
		model1->setImageParameters(ref.xsize(), ref.ysize(), ref.zsize(),ref.xdim(), ref.ydim(), ref.zdim());
		cout<<"load model "<<modelname<<endl;
		model1->load_bmv_binaryInfo(modelname,1);
		model1->load_bmv_binary(modelname,1);

			//want to return mean mesh in vector form
		meanm->ReSize(model1->getTotalNumberOfPoints()*3);
		{
		
			int count=0;
			int cumnum=0; 
			for (int sh=0; sh<model1->getNumberOfShapes();sh++){
				Mesh m=model1->getTranslatedMesh(sh);
				*meshout=m;
				for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
					meanm->element(3*cumnum+count)=(*i)->get_coord().X;
					meanm->element(3*cumnum+count+1)=(*i)->get_coord().Y;
					meanm->element(3*cumnum+count+2)=(*i)->get_coord().Z;
					count+=3;
				}
				cumnum+=model1->getNumberOfPoints(sh);	
			}
			
		}
		
		//need number of subjects (to detrmine number of modes)
		int M=model1->getNumberOfSubjects();
		int Tpts=model1->getTotalNumberOfPoints();
			Matrix MeshVerts(3*Tpts,bvars.Ncols());
		//this is different than mnumber of subjects which were used to create model
		//int numSubs=bvars.Ncols();
		//this loads all mesh point into a matrix to do stats on....
		cout<<"generate and load vertices into matrix "<<endl;
		for (int j=0;j<bvars.Ncols();j++){
			//for each subject
			vector<float> vars;
			for (int i=0; i<bvars.Nrows();i++){
				vars.push_back(bvars.element(i,j));
			}
			//keep track of number of points preceding
			int cumnum=0;
			for (int sh=0; sh<model1->getNumberOfShapes();sh++){
				//cout<<"cumnum "<<cumnum<<endl;
				Mesh m=model1->getDeformedMesh(vars,sh,static_cast<int>(vars.size()));	
				
				int count=0;
				for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
					MeshVerts.element(3*cumnum+count,j)=(*i)->get_coord().X;
					MeshVerts.element(3*cumnum+count+1,j)=(*i)->get_coord().Y;
					MeshVerts.element(3*cumnum+count+2,j)=(*i)->get_coord().Z;
					count+=3;
				}	
				//	cout<<"count "<<count<<endl;
				cumnum+=model1->getNumberOfPoints(sh);	
			}	
		}
	return MeshVerts;

}

Matrix pca_bvars(Matrix Data){
	DiagonalMatrix D;
	Matrix V;
	Matrix U;
	SVD(Data, D,U,V);
	
	return V;
}	
Matrix rigid_pca_bvars(Matrix Data,ColumnVector meanm, Mesh mesh){
	//ColumnVector avgM(sub.Ncols());
	cout<<"Data "<<Data.Nrows()<<" "<<Data.Ncols()<<endl;
	//determine translations
	int Nsub=Data.Ncols();
	int Npoints=Data.Nrows()/3;
	
	//for (int j=0;j<Data.Nrows();j++){
//		cout<<Data.element(j,0)<<endl;
//	}
	
	
	
	//***********CALCULATE CENTROIDS*************//
	//calculate centroid of mean mesh
	float Mxr=0,Myr=0,Mzr=0;
	for (int i=0;i<meanm.Nrows();i=i+3){
		Mxr+=meanm.element(i);
		Myr+=meanm.element(i+1);
		Mzr+=meanm.element(i+2);
	}
	Mxr/=(meanm.Nrows()/3);
	Myr/=(meanm.Nrows()/3);
	Mzr/=(meanm.Nrows()/3);
	
	
	//calculate centroid 
	vector<float> vMx,vMy,vMz;
	for (int i=0;i<Nsub;i++){
	float sx=0,sy=0,sz=0;
	cout<<" i "<<i<<endl;
		for (int j=0;j<Data.Nrows();j=j+3){
			//cout<<"j "<<j<<endl;
			sx+=Data.element(j,i);
			sy+=Data.element(j+1,i);
			sz+=Data.element(j+2,i);
		}
		//cout<<"centroid "<<sx/Npoints<<" "<<sy/Npoints<<" "<<sz/Npoints<<endl;
		vMx.push_back(sx/Npoints);
		vMy.push_back(sy/Npoints);
		vMz.push_back(sz/Npoints);
	}
	Matrix R(3,3);
	for (int subject=0;subject<Data.Ncols();subject++){
			//cout<<"subject "<<subject<<endl;
		
		
		//***********Demena data*************//
		//Deamean the Data and reformat
		//Demean Mean mesh and reformat
		Matrix DataDM(3,Data.Nrows()/3);
		Matrix RefDM(3,Data.Nrows()/3);
		for (int i=0;i<Data.Nrows();i=i+3){
				DataDM.element(0,i/3)=Data.element(i,subject)-vMx.at(subject);
				DataDM.element(1,i/3)=Data.element(i+1,subject)-vMy.at(subject);
				DataDM.element(2,i/3)=Data.element(i+2,subject)-vMz.at(subject);
				RefDM.element(0,i/3)=meanm.element(i)-Mxr;
				RefDM.element(1,i/3)=meanm.element(i+1)-Myr;
				RefDM.element(2,i/3)=meanm.element(i+2)-Mzr;
		}
		
		//***********calculate rotattions*************//
		cout<<"reshaped matrices"<<endl;
		Matrix M=RefDM*(DataDM.t());
		Matrix U;
		DiagonalMatrix D;
		SVD(M.t()*M,D,U);
		//M should always be a 3x3 matrix
		for (int i=0;i<D.Nrows();i++){
		//	if (i==(D.Nrows()-1)){
				cout<<"The singular value is "<<D.element(i)<<endl;
		//	}
			D.element(i)=1/sqrt(D.element(i));
		}
		R=M*(U*D*U.t());
		cout<<"rotation matrix: "<<endl;
	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			cout<<R.element(i,j)<<" ";
		}
		cout<<endl;
		}
		
	}
	
	//*****************APPLY TRANFSORMATIOON TO MESHES**********************//
	//NO SCALE IS CALCULATED
	
	
		Matrix DataNew(Data.Nrows(),Data.Ncols());
	for (int subject=0;subject<Data.Ncols();subject++){
	//	cout<<"subject "<<subject<<endl;
		
		//reshape data
		Matrix DataRS(3,Data.Nrows()/3);
		Matrix RefMean(3,Data.Nrows()/3);
		Matrix DataMean(3,Data.Nrows()/3);
		for (int i=0;i<Data.Nrows();i=i+3){
			//cout<<i/3<<" "<<DataRS.element(i,subject)<<endl;
			DataRS.element(0,i/3)=Data.element(i,subject);
			DataRS.element(1,i/3)=Data.element(i+1,subject);
			DataRS.element(2,i/3)=Data.element(i+2,subject);
	//		if (subject==4){
	//			cout<<i/3<<" "<<DataRS.element(0,i/3)<<endl;
	//		}
			RefMean.element(0,i/3)=Mxr;
			RefMean.element(1,i/3)=Myr;
			RefMean.element(2,i/3)=Mzr;
			DataMean.element(0,i/3)=vMx.at(subject);
			DataMean.element(1,i/3)=vMy.at(subject);
			DataMean.element(2,i/3)=vMz.at(subject);
		}
		///APPLY TRANSFORMATION
		Matrix Reg(3,Data.Nrows()/3);//=R*DataRS+(RefMean-R*DataMean);
		Reg=R*DataRS-(RefMean-R*DataMean);

		for (int i=0;i<Reg.Ncols();i++){
			if (subject==0){
				cout<<i<<" "<<Reg.element(0,i)<<endl;
			}
			DataNew.element(3*i,subject)=Reg.element(0,i);
			DataNew.element(3*i+1,subject)=Reg.element(1,i);
			DataNew.element(3*i+2,subject)=Reg.element(2,i);
		}
	}
	DiagonalMatrix D;
	Matrix V;
	
	SVD(DataNew.t()*DataNew,D,V);
	cout<<V.Nrows()<<" "<<V.Ncols()<<endl;
	Mesh m=mesh;
	int count=0;
	for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
		(*i)->_update_coord.X=DataNew.element(count,0);
		(*i)->_update_coord.Y=DataNew.element(count+1,0);
		(*i)->_update_coord.Z=DataNew.element(count+2,0);
		count+=3;
	}	
				m.update();
				m.save("meshreg.vtk",3);
				
	return V;
}


void do_work_from_bvarsLeave1out(){
	
	vector<int> class_dec;
	vector<float> class_val;
	vector<float> class_fp;
	vector<float> class_fn;

	//read model used to parametrize object
	string mname;
	cout<<"read model name"<<endl;
	mname=read_bvars_ModelName(inname.value() );	
	if (replaceBMV.value()){
		mname=modelname.value();
	}
	
	//cout<<"load the model: "<<mname<<endl;
	//this is used to set the size and dimensions of shape model...not really just effects tranbslation
	volume<short> ref;
	shapeModel* model1 = new shapeModel;
	vector<string> subjectnames;
	Matrix bvars;
	vector<int> vN;
	Matrix target;
	cout<<"read in bvars and target"<<endl;
	read_bvars_design(inname.value(),&bvars,&target,&subjectnames, &vN);
	cout<<bvars.Nrows()<<" "<<bvars.Ncols()<<endl;
	
	if (bvarsnew.value()){
		Matrix verts;
		ColumnVector meanM;
		Mesh m;
		verts=recon_meshesNative(refname.value(),mname,bvars,&meanM,&m);
		
		bvars=rigid_pca_bvars(verts,meanM,m);
	//	Matrix bvars_new;
	//	bvars_new=pca_bvars(verts);
	//	bvars=bvars_new.t();//flir orientation
	}
	
		cout<<bvars.Nrows()<<" "<<bvars.Ncols()<<endl;
	//this is the value we wish to descriminate
	ColumnVector test;
	float sumCorr=0;	
	float sumFP=0;
	float sumFN=0;
	int maxiter=bvars.Ncols();
		for (int iter=1;iter<=maxiter;iter++){//bvars.Nrows();iter++){
				
			Matrix Data;
			
		//	Matrix Kloo;
			Matrix tloo;
			if (iter==1){
			//	cout<<"iter1 "<<iter<<endl; 
				Data=bvars.SubMatrix(1,bvars.Nrows(),2,bvars.Ncols());
				
				tloo=target.SubMatrix(2,target.Nrows(),1,1);
			}else if (iter==bvars.Nrows()){
			//	cout<<"iter2 "<<iter<<endl;
				Data=bvars.SubMatrix(1,bvars.Nrows(),1,iter-1);
				
				tloo=target.SubMatrix(1,target.Nrows()-1,1,1);
			}else{
			//	cout<<"iter3 "<<iter<<endl;
				Data=bvars.SubMatrix(1,bvars.Nrows(),1,iter-1) | bvars.SubMatrix(1,bvars.Nrows(),iter+1,bvars.Ncols());
				tloo=target.SubMatrix(1,iter-1,1,1) & target.SubMatrix(iter+1,target.Nrows(),1,1);
			}
			int dec=0;
			int testTarg=target.element(iter-1,0);
			if (useRVM.value()){
				cout<<"compute kernel"<<endl;
				Matrix K=compute_kernel(Data,Data,kerntype.value(),rad.value());
					Matrix Tcov(K.Ncols(),K.Ncols());
		ColumnVector inds(K.Ncols()); 
				run_rvm(tloo.t(),K,100,&Tcov,&inds);
			}else{
				
				test=bvars.SubMatrix(1,bvars.Nrows(),iter,iter);
				dec=QDA(Data,tloo,test,numPar.value());
			}
			
			
			if(dec==testTarg){
				sumCorr++;
			}else if (dec==0){
				sumFN++;
			}else{
				sumFP++;
			}
			//cout<<"iter "<<iter<<" belongs to class "<<dec<<endl;
		
	}
	cout<<"I got "<<sumCorr<<" out of "<<maxiter<<" "<<sumCorr/maxiter<<endl;
		cout<<"I got "<<sumFP<<" false positive and  "<<sumFN<<" false negatives"<<endl;

}


void do_work_from_bvars(){

	//read model used to parametrize object
	string mname;
	cout<<"read model name"<<endl;
	mname=read_bvars_ModelName(inname.value() );	
	if (replaceBMV.value()){
		mname=modelname.value();
	}

	//cout<<"load the model: "<<mname<<endl;
	//this is used to set the size and dimensions of shape model...not really just effects tranbslation
	volume<short> ref;
		shapeModel* model1 = new shapeModel;
			vector<string> subjectnames;
	Matrix bvars;
	vector<int> vN;
	Matrix target;
	cout<<"read in bvars and target"<<endl;
	read_bvars_design(inname.value(),&bvars,&target,&subjectnames, &vN);
	cout<<bvars.Nrows()<<" "<<bvars.Ncols()<<endl;
	
	
	
	
}

int main(int argc,char *argv[])
{
	
	Tracer tr("main");
	OptionParser options(title, examples);
	
	try {
		// must include all wanted options here (the order determines how
		//  the help message is printed)
		options.add(inname);
		options.add(inname2);
		options.add(rad);
		options.add(bvarsnew);
		options.add(errval);
		options.add(designname);
		options.add(numPar);
		options.add(modelname);
		options.add(refname);
		options.add(replaceBMV);
		options.add(usebvars);
		options.add(useRVM);
	options.add(modeStats);
	options.add(polyorder);
		options.add(outname);
			options.add(kerntype);
		options.add(inputprob);
		options.add(verbose);
		
		options.add(help);
		
		nonoptarg = options.parse_command_line(argc, argv);
		
		// line below stops the program if the help was requested or 
		//  a compulsory option was not set
		if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
		{
			options.usage();
			exit(EXIT_FAILURE);
		}
		
		// Call the local functions
	
		do_work_from_bvarsLeave1out();
		
		//do_work_from_bvars();
		
		
	}  catch(X_OptionError& e) {
		options.usage();
		cerr << endl << e.what() << endl;
		exit(EXIT_FAILURE);
	} catch(std::exception &e) {
		cerr << e.what() << endl;
	} 
	
	return 0;// do_work(argc,argv);
}

