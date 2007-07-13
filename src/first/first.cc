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


string title="first (Version 1.0) University of Oxford (Brian Patenaude)";
string examples="first --baam  -i <input image> -l <flirt matrix> -m <model> -g <number of modes> -y <rob_min> -z <rob_max> ";


Option<bool> verbose(string("-v,--verbose"), false, 
				 string("switch on diagnostic messages"), 
				 false, no_argument);
Option<bool> help(string("-h,--help"), false,
			  string("display this message"),
			  false, no_argument);
Option<bool> inputprob(string("--inputprob"), false,
					   string("Input image is probability - do not do any thresholding or smoothing"),
					   false, no_argument);
Option<bool> loadvars(string("--loadvars"), false,
					   string("load intial parameter estimates from a previous segmentation"),
					   false, no_argument);
Option<bool> shcond(string("--shcond"), false,
					   string("use conditional probability on shape"),
					   false, no_argument);
Option<bool> shcond2(string("--shcond2"), false,
					   string("use conditional probability on 2 shapes"),
					   false, no_argument);
Option<bool> reestimate(string("--modeest"), false,
					   string("reeistmate mode for each iteration"),
					   false, no_argument);
Option<bool> useIntRefModel(string("--useIntRefModel"), false,
					   string("reeistmate mode for each iteration"),
					   false, no_argument);
Option<bool> runAff(string("--runAff"), false,
					   string("run first using affine modes"),
					   false, no_argument);
Option<bool> modeprior(string("--modeprior"), false,
					   string("reeistmate mode for each iteration"),
					   false, no_argument);
//Option<bool> nograd(string("--nograd"), false,
//					   string("allows "),
//					   false, no_argument);
Option<bool> baam(string("--baam"), false,
					   string("use appearaance model"),
					   false, no_argument);
Option<bool> joint(string("--joint"), false,
					   string("use joint model fit"),
					   false, no_argument);
Option<bool> jointN(string("--jointN"), false,
					   string("use jointN model fit"),
					   false, no_argument);
Option<bool> jointN2(string("--jointN2"), false,
					   string("use jointN model fit"),
					   false, no_argument);
Option<bool> histmode(string("--mode"), false,
					   string("use mode of distribution"),
					   false, no_argument);
Option<bool> gmmLesser(string("--lesser"), true,
					   string("use lesser of means in misture model, otherwise use greater"),
					   false, no_argument);
Option<bool> overlap(string("--overlap"), false,
					   string("use overlapcost - can fit mesh to labelled data"),
					   false, no_argument);
Option<bool> useCondMode(string("--useCondMode"), false,
					   string("use conditional mode"),
					   false, no_argument);
Option<bool> useConj3(string("--useConj3"), false,
					   string("use other multifit"),
					   false, no_argument);


Option<bool> singleInit(string("--singleInit"), false,
					   string("initialise structure with single fit"),
					   false, no_argument);
					   
Option<string> inname(string("-i,--in"), string(""),
					  string("filename of input image to be segmented"),
					  true, requires_argument);
Option<string> initname(string("-j,--in"), string(""),
					  string("filename of which structure to initilize"),
					  false, requires_argument);

Option<string> flirtmatname(string("-l,--flirtMatrix"), string(""),
						 string("filename containing flirt matrix (transformattion to MNI space"),
						 true, requires_argument);
Option<string> modelname(string("-m,--modelin"), string(""),
						 string("filename of input model: the structure to be segmented"),
						 true, requires_argument);
Option<string> modelname2(string("-n,--modelin2"), string(""),
						 string("filename of 2nd input model for dual: the structure to be segmented"),
						 false, requires_argument);
Option<int> EMiter(string("-e"), 200,
					   string("maximum number of iterations in the em"),
					   false, requires_argument);
Option<string> bvarsname(string("-o,--name of bvars file your conditioning on"), string(""),
						  string("load bvars"),
						  false, requires_argument);
Option<string> bvarsname2(string("-p,--name of second bvars file your conditioning on"), string(""),
						  string("load bvars2"),
						  false, requires_argument);
Option<string> bmapname(string("-b,--bmapin"), string(""),
						  string("load bmap"),
						  false, requires_argument);
Option<string> bmapname2(string("-d,--bmap2in"), string(""),
						  string("load bmap2"),
						  false, requires_argument);
Option<string> outname(string("-k,--out"), string(""),
					   string("basename of output"),
					   true, requires_argument);
Option<float> stdTrunc(string("-s,--sh"), 5.0,
					   string("number of standard deviation to truncate at"),
					   false, requires_argument);		
Option<float> manMean(string("-a,--mean intensity"), -777.0,
				string("override mean estimation"),
				false, requires_argument);
Option<float> c(string("-c,--cycles"), 1.0,
				string("how many multiples of g"),
				false, requires_argument);
Option<float> f(string("-f,--cycles"), 1.0,
				string("how many cycles"),
				false, requires_argument);
Option<int> g(string("-g,--nmodes"), 10,
			  string("numberOfModes in window"),
			  false, requires_argument);
Option<float> robmin(string("-y,--robmin"), 0.0,
				string("how many multiples of g"),
				false, requires_argument);
Option<float> robmax(string("-z,--robmax"), 255.0,
				string("robust max for normalization"),
				false, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////
//global variables
volume<float> image;
volume<float> Resimage;
int costmodeG=0;
float xdim, ydim, zdim;
float searchDim;
float searchRes=0.2;
volume4D<short> fit;
volume<short> d4tmp;
bool gradzero=false;
const float pi = 3.14159;
//unsigned int M;
vector<float> v_ceigs;
ColumnVector mBx2map, mBx2map2;
Matrix mBx1inv,mBx1inv2;
vector< Matrix > VmBx1inv, VmBx2;
bool GlobShCond=true;
float kpred,kpred2;
float globMean=0;
vector<float> VglobMean;
int gbounds[6]={1000,0,1000,0,1000,0};
float gBeta;
float GmanMean;
float condMode;
bool globfoundmode=false;
vector<bool> Vglobfoundmode;
bool globReest=false;
float globModePriorMean=-777;
	float	globModePriorVar=-777;

int readBmapHalf(string fname, vector<float> b2, Matrix *matBx2){
	int M; //number of subjects
	
	ifstream fin;
	fin.open(fname.c_str());
	string stemp;
	string mapto, ref;
	getline(fin,stemp);
	cout<<stemp<<endl;
	fin>>stemp;
	fin>>M>>mapto>>ref;
	cout<<M<<endl;
	matBx2->ReSize(M,M);
	getline(fin,stemp); 
	getline(fin,stemp);
	//read first matrix 
	//	float ftemp;
	double ftemp;
	//Matrix matBx2(M,M);
	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBx2->element(i,j)=ftemp;
	
		}
	}
	return M;
}
int readBmap(string fname, vector<float> b2, Matrix *matBx2, Matrix * matBcx1){
	int M; //number of subjects
	
	ifstream fin;
	//cout<<"open file "<<fname.c_str()<<endl;
	fin.open(fname.c_str());
//	cout<<"opened"<<endl;
	string stemp;
	
	getline(fin,stemp);
//	cout<<fin<<" "<<stemp<<endl;
	fin>>stemp;
	fin>>M;
//	cout<<M<<endl;
	matBx2->ReSize(M,M);
	matBcx1->ReSize(M,M);
	getline(fin,stemp); 
	getline(fin,stemp);
//	cout<<"bmap: "<<stemp<<endl;
	//read first matrix 
	//	float ftemp;
	double ftemp;
	//Matrix matBx2(M,M);
	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBx2->element(i,j)=ftemp;
		}
	}
	
//	cout<<"readbmap interior "<<endl;
	//cout<<BInv<<endl;
	getline(fin,stemp);
	getline(fin,stemp);
	//Matrix matBcx1(M,M);
	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBcx1->element(i,j)=ftemp;
			
			//if (j==0){
			//	cout<<i<<" "<<ftemp<<" "<<BInv.element(i,j)<<endl;
			     //	     cout<<i<<" "<<ftemp<<" "<<endl;
			     // } 
		}
	}
	getline(fin,stemp);
	getline(fin,stemp);
	float kx2;
	fin>>kx2;
	kpred=kx2;
	
	//try kx2=number of modes
//	kx2=20;
	
	//calculate scale fcator for matBCx1
	float alp=M-1.0/M;
	float gammav=(alp)/(alp-2);
	
	//calculate binner
	float binner=0;
	for (unsigned int i=0;i<b2.size();i++){
		binner+=b2.at(i)*b2.at(i);
	}
	
	float scale=sqrt((alp+binner*gammav)/(alp+kx2));
	//this part causes some practical problems
	*matBcx1=(*matBcx1)*scale;
	//	*matBx2=(*matBx2)*scale;
	getline(fin,stemp);
	getline(fin,stemp);
	//float ftemp;
	v_ceigs.clear();
	
	for (int i=0;i<M;i++){
	  fin>>ftemp;
	  v_ceigs.push_back(ftemp);
	}
	
	return M;
	
	
	
}

int readBmapNoScale(string fname, vector<float> b2, Matrix *matBx2, Matrix * matBcx1){
  int M; //number of subjects
	
  ifstream fin;
  //cout<<"open file "<<fname.c_str()<<endl;
  fin.open(fname.c_str());
  //	cout<<"opened"<<endl;
  string stemp;
	
  getline(fin,stemp);
  //	cout<<fin<<" "<<stemp<<endl;
  fin>>stemp;
  fin>>M;
  cout<<M<<endl;
  matBx2->ReSize(M,M);
  matBcx1->ReSize(M,M);
  getline(fin,stemp); 
  getline(fin,stemp);
  //	cout<<"bmap: "<<stemp<<endl;
  //read first matrix 
  //	float ftemp;
  double ftemp;
  //Matrix matBx2(M,M);
  for (int i=0;i<M;i++){
    for (int j=0;j<M;j++){
      fin>>ftemp;
      matBx2->element(i,j)=ftemp;
    }
  }
	
  //	cout<<"readbmap interior "<<endl;
  //cout<<BInv<<endl;
  getline(fin,stemp);
  getline(fin,stemp);
  //Matrix matBcx1(M,M);
  for (int i=0;i<M;i++){
    for (int j=0;j<M;j++){
      fin>>ftemp;
      matBcx1->element(i,j)=ftemp;
			
      //if (j==0){
      //	cout<<i<<" "<<ftemp<<" "<<BInv.element(i,j)<<endl;
      //	     cout<<i<<" "<<ftemp<<" "<<endl;
      // } 
    }
  }
  getline(fin,stemp);
  getline(fin,stemp);

  float kx2;
  fin>>kx2;
  kpred=kx2;

  //calculate scale fcator for matBCx1
  //float alp=M-1.0/M;
  //	float gammav=(alp)/(alp-2);

  /*	
  //calculate binner
  float binner=0;
  for (unsigned int i=0;i<b2.size();i++){
  binner+=b2.at(i)*b2.at(i);
  }
	
  float scale=sqrt((alp+binner*gammav)/(alp+kx2));
  //this part causes some practical problems
  *matBcx1=(*matBcx1)*scale;
  */

  //	*matBx2=(*matBx2)*scale;
  getline(fin,stemp);
  getline(fin,stemp);
  //float ftemp;
  v_ceigs.clear();
	
  for (int i=0;i<M;i++){
    fin>>ftemp;
    v_ceigs.push_back(ftemp);
  }
	
  return kx2;
}

int readMultiBmapNoScale(string fname, vector<float> b2, vector< Matrix >* vmatBx2, Matrix * matBcx1,  vector<int>* vPredLabel){
  int M; //number of subjects
  int N;//number of predictors
  ifstream fin;
  //cout<<"open file "<<fname.c_str()<<endl;
  fin.open(fname.c_str());
  
  string stemp;
	
  getline(fin,stemp);
  fin>>stemp;
  fin>>M>>N;
  cout<<M<<" "<<N<<endl;
  matBcx1->ReSize(M,M);	
  //need ot read in matruices to calculate gam
  int P=0;
  //theretically N=1 has one only its identity
  if (N==2){
    P=3;
  }

  getline(fin,stemp); 
  getline(fin,stemp);
  cout<<"bmap: "<<stemp<<endl;
  //read in first bmap	
  for (int i=0;i<M;i++){
    for (int j=0;j<M;j++){
      float ftemp;
      fin>>ftemp;
      matBcx1->element(i,j)=ftemp;
    }
  }
  getline(fin,stemp);
  cout<<"getline 2"<<stemp<<endl;
  //dont need to invert because of different implementation in appmodel
			
  //read first matrix 
  //this is the ocnditional mean matrix map
  vmatBx2->clear();
  for (int i=0;i<N;i++){
    getline(fin,stemp);
    cout<<"getline 3"<<stemp<<endl;

    Matrix matBx2(M,M);
    double ftemp;
    //Matrix matBx2(M,M);
    for (int i=0;i<M;i++){
      for (int j=0;j<M;j++){
	fin>>ftemp;
	matBx2.element(i,j)=ftemp;
      }
    }
    vmatBx2->push_back(matBx2);
    getline(fin,stemp);
  }
  //	cout<<"readbmap interior "<<endl;
  //cout<<BInv<<endl;
	
  //Matrix matBcx1(M,M);
  getline(fin,stemp);
  cout<<"getline1 "<<stemp<<endl;

//     vmatG->clear();
//   for (int i=0;i<P;i++){
//     getline(fin,stemp);
//     cout<<"getline 3"<<stemp<<endl;

//     Matrix matG(M,M);
//     double ftemp;
//     //Matrix matBx2(M,M);
//     for (int i=0;i<M;i++){
//       for (int j=0;j<M;j++){
// 	fin>>ftemp;
// 	matG.element(i,j)=ftemp;
//       }
//     }
//     vmatG->push_back(matG);
//     getline(fin,stemp);
//   }

  //Matrix matBcx1(M,M);
  //	getline(fin,stemp);
  //cout<<"getline 3"<<stemp<<endl;

  float kx2;
  fin>>kx2;
  kpred=kx2;
  cout<<kx2<<endl;
  //calculate scale fcator for matBCx1
  //float alp=M-1.0/M;
  //	float gammav=(alp)/(alp-2);

  /*	
  //calculate binner
  float binner=0;
  for (unsigned int i=0;i<b2.size();i++){
  binner+=b2.at(i)*b2.at(i);
  }
	
  float scale=sqrt((alp+binner*gammav)/(alp+kx2));
  //this part causes some practical problems
  *matBcx1=(*matBcx1)*scale;
  */

  //	*matBx2=(*matBx2)*scale;
  getline(fin,stemp);
  getline(fin,stemp);
  cout<<"getline 3 "<<stemp<<endl;
  //getline(fin,stemp);
  //cout<<"getline 4 "<<stemp<<endl;
  //float ftemp;
  //v_ceigs.clear();
  vPredLabel->clear();
  for (int i=0;i<N;i++){
    int itemp;
		
    fin>>itemp;
    cout<<itemp<<" ";
    vPredLabel->push_back(itemp);
  }
  cout<<endl;
  return kx2;
}

void scaleBInvMatrix(vector<float> b2, Matrix * matBinv1, int M, int kx2){


	//calculate scale fcator for matBCx1
	float alp=M-1.0/M;
	float gammav=(alp)/(alp-2);

	//calculate binner
	float binner=0;
	for (unsigned int i=0;i<b2.size();i++){
		binner+=b2.at(i)*b2.at(i);
	}
	
	float scale=sqrt((alp+binner*gammav)/(alp+kx2));
	//this part causes some practical problems
		*matBinv1=(*matBinv1)/scale;	
}
float scaleMultiBInvMatrix(vector< vector<float> > b2, Matrix * matBinv1, int M, int kx2){


	//calculate scale fcator for matBCx1
	float alp=M-1.0/M;
	float gammav=(alp)/(alp-2);

	//calculate binner
	float binner=0;
	
	for (unsigned int i=0;i<b2.size();i++){
		for (unsigned int j=0;j<b2.at(0).size();j++){

		binner+=b2.at(i).at(j)*b2.at(i).at(j);
}
	}
	
	float scale=((alp+binner*gammav)/(alp+kx2));
	//this part causes some practical problems
		*matBinv1=(*matBinv1)/sqrt(scale);
		return scale;	
}


int readBmap2(string fname, vector<float> b2, Matrix *matBx2, Matrix * matBcx1){
//only diffence is the assignment of krped2 versus kpred, anitquate later
	int M; //number of subjects
	
	ifstream fin;
	cout<<"open file "<<fname.c_str()<<endl;
	fin.open(fname.c_str());
	cout<<"opened"<<endl;
	string stemp;
	
	getline(fin,stemp);
	cout<<fin<<" "<<stemp<<endl;
	fin>>stemp;
	fin>>M;
	cout<<M<<endl;
	matBx2->ReSize(M,M);
	matBcx1->ReSize(M,M);
	getline(fin,stemp); 
	getline(fin,stemp);
	cout<<"bmap: "<<stemp<<endl;
	//read first matrix 
	//	float ftemp;
	double ftemp;
	//Matrix matBx2(M,M);
	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBx2->element(i,j)=ftemp;
			
			//cout<<"i "<<i<<" "<<"j"<<" "<<j<<endl;
			//	   if (j==0){
			
			//   cout<<i<<" "<<ftemp<<" "<<endl;
			     //	cout<<i<<" "<<ftemp<<" "<<BInv.element(i,j)<<endl;
			// } 
		}
	}
	
	cout<<"readbmap interior "<<endl;
	//cout<<BInv<<endl;
	getline(fin,stemp);
	getline(fin,stemp);
	//Matrix matBcx1(M,M);
	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBcx1->element(i,j)=ftemp;
			
			//if (j==0){
			//	cout<<i<<" "<<ftemp<<" "<<BInv.element(i,j)<<endl;
			     //	     cout<<i<<" "<<ftemp<<" "<<endl;
			     // } 
		}
	}
	getline(fin,stemp);
	getline(fin,stemp);
	float kx2;
	fin>>kx2;
	kpred2=kx2;
	
	//try kx2=number of modes
//	kx2=20;
	
	//calculate scale fcator for matBCx1
	float alp=M-1.0/M;
	float gammav=(alp)/(alp-2);
	
	//calculate binner
	float binner=0;
	for (unsigned int i=0;i<b2.size();i++){
		binner+=b2.at(i)*b2.at(i);
	}
	
	float scale=sqrt((alp+binner*gammav)/(alp+kx2));
	//this part causes some practical problems
		*matBcx1=(*matBcx1)*scale;
	getline(fin,stemp);
	getline(fin,stemp);
	//float ftemp;
	v_ceigs.clear();
	
	for (int i=0;i<M;i++){
		fin>>ftemp;
		v_ceigs.push_back(ftemp);
	}
	
		return M;
	
	
	
}



void getBounds(Mesh m, int *bounds, float xdim, float ydim, float zdim){
	
	float xmin=1000,xmax=-1000,ymin=1000,ymax=-1000,zmin=1000,zmax=-1000;
	for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
		float tempx=(*i)->get_coord().X;
		float tempy=(*i)->get_coord().Y;
		float tempz=(*i)->get_coord().Z;
		if (tempx<xmin){
			xmin=tempx;
		}
		if (tempx>xmax){
			xmax=tempx;
		}
		if (tempy<ymin){
			ymin=tempy;
		}
		if (tempy>ymax){
			ymax=tempy;
		}
		if (tempz<zmin){
			zmin=tempz;
		}
		if (tempz>zmax){
			zmax=tempz;
		}
	}
	*bounds=static_cast<int>(floor(xmin/xdim)-1);
	*(bounds+1)=static_cast<int>(ceil(xmax/xdim)+1);
	*(bounds+2)=static_cast<int>(floor(ymin/ydim)-1);
	*(bounds+3)=static_cast<int>(ceil(ymax/ydim)+1);
	*(bounds+4)=static_cast<int>(floor(zmin/zdim)-1);
	*(bounds+5)=static_cast<int>(ceil(zmax/zdim)+1);
	
}

void draw_segment(volume<short>& image, const Pt& p1, const Pt& p2, int label)
{
	double xdim = (double) image.xdim();
	double ydim = (double) image.ydim();
	double zdim = (double) image.zdim();

	//in new version of bet2
	double mininc = min(xdim,min(ydim,zdim)) * .5;


	
	Vec n = (p1 - p2);
	double d = n.norm();
	n.normalize();
	//	double l = d*4;
	for (double i=0; i<=d; i+=mininc)
    {
		Pt p = p2 + i* n;
		image((int) floor((p.X)/xdim +.5),(int) floor((p.Y)/ydim +.5),(int) floor((p.Z)/zdim +.5)) = label;
	
    }
}


volume<short> draw_mesh(const volume<short>& image, const Mesh &m, int label)
{
	volume<short> res = image;
	for (list<Triangle*>::const_iterator i = m._triangles.begin(); i!=m._triangles.end(); i++)
    {
		Vec n = (*(*i)->get_vertice(0) - *(*i)->get_vertice(1));
		double d = n.norm();
		n.normalize();
		double l=d*2;
		
		for (int j=0; j<=l /*(floor(1.3*d + 1)) + 1*/; j++)
		{
			Pt p = (*i)->get_vertice(1)->get_coord()  + (double)j*.5 * n;
			draw_segment(res, p, (*i)->get_vertice(2)->get_coord(),label);
		} 
    }
	return res;
}
volume<short> make_mask_from_mesh(const volume<float> & image, const Mesh& m, int label, int* bounds)
{
	//fill outside then inverse	
	float xdim = (float) image.xdim();
	float ydim = (float) image.ydim();
	float zdim = (float) image.zdim();
	
	volume<short> mask;
	copyconvert(image,mask);
	
	mask = 0;
	mask = draw_mesh(mask, m,label);
	// THIS EXCLUDEDS THE ACTUAL MESH
	volume<short> otl=mask;		
		
		getBounds(m,bounds,xdim,ydim,zdim);
		vector<Pt> current;
		current.clear();
		Pt c(bounds[0]-2, bounds[2]-2, bounds[4]-2);
		
		mask.value(static_cast<int>(c.X),static_cast<int>(c.Y),static_cast<int>(c.Z)) = label;
		current.push_back(c);
		int fillCount=0;
		while (!current.empty())
		{
			Pt pc = current.back();
			int x, y, z;
			x=(int) pc.X; y=(int) pc.Y; z=(int) pc.Z;
			current.pop_back();
			fillCount++;
			
			
			if (bounds[0]<=x-1 && mask.value(x-1, y, z)==0) {
				mask.value(x-1, y, z) = label;
				current.push_back(Pt(x-1, y, z));
			}
			if (bounds[2]<=y-1 && mask.value(x, y-1, z)==0) {
				mask.value(x, y-1, z) = label;
				current.push_back(Pt(x, y-1, z));
			}
			if (bounds[4]<=z-1 && mask.value(x, y, z-1)==0) {
				mask.value(x, y, z-1) = label;
				current.push_back(Pt(x, y, z-1));
			}
			if (bounds[1]>=x+1 && mask.value(x+1, y, z)==0){
				mask.value(x+1, y, z) = label;
				current.push_back(Pt(x+1, y, z));
			}
			if (bounds[3]>=y+1 && mask.value(x, y+1, z)==0){
				mask.value(x, y+1, z) = label;
				current.push_back(Pt(x, y+1, z));
			}
			if (bounds[5]>=z+1 && mask.value(x, y, z+1)==0){
				mask.value(x, y, z+1) = label;
				current.push_back(Pt(x, y, z+1)); 
			}
			
		}
for (int i=bounds[0];i<bounds[1];i++){
    for (int j=bounds[2];j<bounds[3];j++){
		for (int k=bounds[4];k<bounds[5];k++){
			if (mask.value(i,j,k)==0){
				otl.value(i,j,k)=label;
			}
		}
    }
}
return otl;
}
volume<short> make_mask_from_meshInOut(const volume<float> & image, const Mesh& m, int label, int* bounds)
{
	
	float xdim = (float) image.xdim();
	float ydim = (float) image.ydim();
	float zdim = (float) image.zdim();
	
	volume<short> mask;
	copyconvert(image,mask);
	
	
	mask = 0;
	mask = draw_mesh(mask, m,label+100);
	
	
	
	// THIS EXCLUDEDS THE ACTUAL MESH
	volume<short> otl=mask;
		getBounds(m,bounds,xdim,ydim,zdim);
		vector<Pt> current;
		current.clear();
		Pt c(bounds[0]-2, bounds[2]-2, bounds[4]-2);
		
		mask.value(static_cast<int>(c.X),static_cast<int>(c.Y),static_cast<int>(c.Z)) = label;
		current.push_back(c);
		int fillCount=0;
		while (!current.empty())
		{
			Pt pc = current.back();
			int x, y, z;
			x=(int) pc.X; y=(int) pc.Y; z=(int) pc.Z;
			
			current.pop_back();
			fillCount++;
			
			
			if (bounds[0]<=x-1 && mask.value(x-1, y, z)==0) {
				mask.value(x-1, y, z) = label;
				current.push_back(Pt(x-1, y, z));
			}
			if (bounds[2]<=y-1 && mask.value(x, y-1, z)==0) {
				mask.value(x, y-1, z) = label;
				current.push_back(Pt(x, y-1, z));
			}
			if (bounds[4]<=z-1 && mask.value(x, y, z-1)==0) {
				mask.value(x, y, z-1) = label;
				current.push_back(Pt(x, y, z-1));
			}
			if (bounds[1]>=x+1 && mask.value(x+1, y, z)==0){
				mask.value(x+1, y, z) = label;
				current.push_back(Pt(x+1, y, z));
			}
			if (bounds[3]>=y+1 && mask.value(x, y+1, z)==0){
				mask.value(x, y+1, z) = label;
				current.push_back(Pt(x, y+1, z));
			}
			if (bounds[5]>=z+1 && mask.value(x, y, z+1)==0){
				mask.value(x, y, z+1) = label;
				current.push_back(Pt(x, y, z+1)); 
			}
			
		}
for (int i=bounds[0];i<bounds[1];i++){
    for (int j=bounds[2];j<bounds[3];j++){
		for (int k=bounds[4];k<bounds[5];k++){
			if (mask.value(i,j,k)==0){
				otl.value(i,j,k)=label;
			}
		}
    }
}
return otl;
}


void addVolumeTo4D(shapeModel* model1, vector<float> vars){
	//Add first translation to 4d analyze file
	int bounds[6]={0,0,0,0,0,0};
//	cout<<"VARSto Add"<<endl;
//	for (int i=0; i<vars.size();i++){
//	  cout<<vars.at(i)<<" ";
//	}
	cout<<endl;
	Mesh mtemp = model1->getDeformedMesh(vars,0, static_cast<int>(vars.size()));
	d4tmp=make_mask_from_meshInOut(image,mtemp, model1->getLabel(0),bounds);
	//	boundaryCleanUp(&d4tmp, mtemp, model1->getLabel(0));
	volume<short> d4tmp2; 
	for (int i=1;i<model1->getNumberOfShapes();i++){
		cout<<"c3 "<<model1->getNumberOfShapes()<<model1->getLabel(i)<<endl;
		
		mtemp = model1->getDeformedMesh(vars,i,static_cast<int>(vars.size()));
		d4tmp2=make_mask_from_mesh(image,mtemp, model1->getLabel(i),bounds);
		//	save_volume(d4tmp2,"prepoutline");
			//	boundaryCleanUp(&d4tmp2, mtemp, model1->getLabel(i));
		  d4tmp=d4tmp+d4tmp2;//* model1->getLabel(i);
	}
	fit.addvolume(d4tmp);
}

vector<float> bTransform(vector<float>* vars,Matrix Bmat,unsigned int M){
  	//load vars 2 (predictor) into a columnvector (newmat)
	ColumnVector Cvars2(M);
	for (unsigned int i =0; i<M;i++){
		if (i < vars->size()){
			Cvars2.element(i)=vars->at(i);
		}else{
			Cvars2.element(i)=0;
		}
	}
	
	//transform vars vars to conditional 
	ColumnVector Cvars1(M);
	Cvars1=Bmat*Cvars2;
	vector<float> v_vars1;
	for (unsigned int i =0; i<M;i++){
		v_vars1.push_back(Cvars1.element(i));
	}
	return v_vars1;
	}
	
vector<float> bTransformFull(vector<float>* vars, ColumnVector v_cmean, Matrix Bmat,unsigned int M){	
	//load vars 2 (predictor) into a columnvector (newmat)
	ColumnVector Cvars2(M);
	for (unsigned int i =0; i<M;i++){
		if (i < vars->size()){
			Cvars2.element(i)=vars->at(i);
		}else{
			Cvars2.element(i)=0;
		}
	}
		//transform vars vars to conditional 
	ColumnVector Cvars1(M);
	Cvars1=v_cmean + Bmat*Cvars2;
	vector<float> v_vars1;
	
	for (unsigned int i =0; i<M;i++){
		v_vars1.push_back(Cvars1.element(i));
	}
	return v_vars1;
	}
ColumnVector bTransformMat(vector<float>* vars,Matrix Bmat,unsigned int M){
	
	//load vars 2 (predictor) into a columnvector (newmat)
	ColumnVector Cvars2(M);
	for (unsigned int i =0; i<M;i++){
		if (i < vars->size()){
			Cvars2.element(i)=vars->at(i);
		}else{
			Cvars2.element(i)=0;
		}
	}
	//transform vars vars to conditional 
	ColumnVector Cvars1(M);
	Cvars1=Bmat*Cvars2;
	
	return Cvars1;
	}
void calculate_meanAndVar(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode, float* mean, float * var){
	volume<short> mask;
	volume<short> maskotl;
	float rdf=0;
	float avgres=0;
	int bounds[6]={0,0,0,0,0,0};
	
	float betas=0;
	
	ColumnVector Best(4);
	vector<float> iprof;
	//new cost
	vector<float> dif;
	for (int i=0;i<model1->getNumberOfShapes()-ex;i++){
		
		float si=0.0,sx=0.0,sy=0.0,sz=0.0;	
		float ssi=0.0, ssx=0.0, ssy=0.0, ssz=0.0;
		float sxy=0.0, sxz=0.0, syz=0.0;		
		float sxi=0.0,syi=0.0,szi=0.0;
		int numpix=0;
		betas=0;
		Mesh m = model1->getDeformedMesh(vars,i,static_cast<int>(vars.size()));
		
		mask=make_mask_from_meshInOut(image,m,model1->getLabel(i),bounds);
		for (int x=bounds[0];x<bounds[1];x++){
			for (int y=bounds[2];y<bounds[3];y++){
				for (int z=bounds[4];z<bounds[5];z++){
				
					//this now operates only on interior to calculate best
					if ((mask.value(x,y,z)==model1->getLabel(i))){
						float tempi;
						tempi=image.value(x,y,z);
						si+=tempi;
						ssi+=tempi*tempi;
						sx+=x;
						sy+=y;
						sz+=z;
						sxi+=x*tempi;
						syi+=y*tempi;
						szi+=z*tempi;
						sxy+=x*y;
						sxz+=x*z;
						syz+=y*z;
						ssx+=x*x;
						ssy+=y*y;
						ssz+=z*z;
						numpix+=1;
						
					}
				}
			}
		}
	//	cout<<"total number of pixels "<<numpix<<endl;
		rdf=numpix-1;
		Matrix A(4,4);											
		Matrix XY(4,1);
		
		A<<numpix<<sx<<sy<<sz
			<<sx<<ssx<<sxy<<sxz
			<<sy<<sxy<<ssy<<syz
			<<sz<<sxz<<syz<<ssz;
		XY<<si<<sxi<<syi<<szi;
		
		Best=A.i()*XY;		
		Matrix val(1,1);
		val=Best.t()*XY;

		avgres=(ssi-val.element(0,0))/numpix;
		
		//begin resdiual profiloe difference portion
		model1->residual(Best, &image, &Resimage,&m, 12);
		float sres=0,ssres=0;
		for (int x=bounds[0];x<bounds[1];x++){
			for (int y=bounds[2];y<bounds[3];y++){
				for (int z=bounds[4];z<bounds[5];z++){
				
					//this now operates only on interior to calculate best
					if ((mask.value(x,y,z)==model1->getLabel(i)+1)){
						//the +1 indicates the coundary
						sres+=Resimage.value(x,y,z);
						ssres+=Resimage.value(x,y,z)*Resimage.value(x,y,z);
							
					}
				}
			}
		}
		///now we change and evaluate probability of voxels belonging to inside
		
		//these define the intensity distribution
		//variance of distyribution inside
		*var=(ssres-(sres*sres/numpix))/(numpix-1);
		*mean=sres/numpix;
}
}



float costfuncInOut(const volume<short> & mask,int* bounds, shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode){
	float rdf=0;
	float avgres=0;
	double cost=0;
	
	float betas=0;
	
	ColumnVector Best(4);
	vector<float> iprof;
	//new cost
	vector<float> dif;
	int bcount=0;
	for (int i=0;i<model1->getNumberOfShapes()-ex;i++){
		
		float si=0.0,sx=0.0,sy=0.0,sz=0.0;	
		float ssi=0.0, ssx=0.0, ssy=0.0, ssz=0.0;
		float sxy=0.0, sxz=0.0, syz=0.0;		
		float sxi=0.0,syi=0.0,szi=0.0;
		int numpix=0;
		betas=0;
			Mesh m = model1->getDeformedMesh(vars,i,static_cast<int>(vars.size()));
		
		for (int x=bounds[0];x<bounds[1];x++){
			for (int y=bounds[2];y<bounds[3];y++){
				for (int z=bounds[4];z<bounds[5];z++){
				
					//this now operates only on interior to calculate best
					if ((mask.value(x,y,z)==model1->getLabel(i))){
						float tempi;
						tempi=image.value(x,y,z);
						si+=tempi;
						ssi+=tempi*tempi;
						sx+=x;
						sy+=y;
						sz+=z;
						sxi+=x*tempi;
						syi+=y*tempi;
						szi+=z*tempi;
						sxy+=x*y;
						sxz+=x*z;
						syz+=y*z;
						ssx+=x*x;
						ssy+=y*y;
						ssz+=z*z;
						numpix+=1;
						
					}
				}
			}
		}
	//	cout<<"total number of pixels "<<numpix<<endl;
		rdf=numpix-1;
		Matrix A(4,4);											
		Matrix XY(4,1);
		
		A<<numpix<<sx<<sy<<sz
			<<sx<<ssx<<sxy<<sxz
			<<sy<<sxy<<ssy<<syz
			<<sz<<sxz<<syz<<ssz;
		XY<<si<<sxi<<syi<<szi;
		
		Best=A.i()*XY;		
		Matrix val(1,1);
		val=Best.t()*XY;

		avgres=(ssi-val.element(0,0))/numpix;
		
		//begin resdiual profiloe difference portion
		model1->residual(Best, &image, &Resimage,&m, 12);
		
		///now we change and evaluate probability of voxels belonging to inside
		
		//these define the intensity distribution
		//variance of distyribution inside
		
		
		float sres=0,ssres=0;
		for (int x=bounds[0];x<bounds[1];x++){
			for (int y=bounds[2];y<bounds[3];y++){
				for (int z=bounds[4];z<bounds[5];z++){
				
					//this now operates only on interior to calculate best
					if ((mask.value(x,y,z)==model1->getLabel(i))){
						//the +1 indicates the coundary
						sres+=Resimage.value(x,y,z);
						ssres+=Resimage.value(x,y,z)*Resimage.value(x,y,z);
							
					}
				}
			}
		}
		float var=(ssres-(sres*sres/numpix))/(numpix-1);
		float mean=sres/numpix;
		
		
		for (int x=bounds[0];x<bounds[1];x++){
			for (int y=bounds[2];y<bounds[3];y++){
				for (int z=bounds[4];z<bounds[5];z++){
				
					//this now operates only on interior to calculate best
					if ((mask.value(x,y,z)==model1->getLabel(i)+1)){
						//the +1 indicates the coundary
						cost+=(Resimage.value(x,y,z)-mean)*(Resimage.value(x,y,z)-mean)/var;
							bcount++;
					}
				}
			}
		}
		
		cost+=bcount*log(var*2*pi);
		} 
		cost*=0.5;
		
		return cost;

}
float logGammaHalf(int Dsh){
	float logGaFunc=0;
	//calculate gamma function for Dsh/2
	if (Dsh%2==0){
		//positive integer into gamma function, calculate factorial of (N/2-1) -1
		for (int d=2;d<=(Dsh/2-1);d++){
			logGaFunc+=log(d);
		}
	}else{
		///gamma(n+0.5)
		int n=(Dsh-1)/2;
		for (int d=3;d<=(2*n-1);d+=2){
			logGaFunc+=log(d);
		}
		logGaFunc-=n*log(2);
		logGaFunc+=log(1.7724538); //mutliply by gam(0.5)= sqrt(pi)
	}
	return logGaFunc;
	
}

//float costfuncProfGrad(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode){
   float costfunc(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode){
	//cout<<"begin cost"<<endl;
	volume<short> mask;
	volume<short> maskotl;
	
	double cost=0;
	
		float gradtot=0.0;
	ColumnVector Best(4);
	vector<float> iprof;
	//new cost
	vector<float> dif;	
	int ipp=5;
	for (int i=0;i<model1->getNumberOfShapes()-ex;i++){
		Mesh m = model1->getDeformedMesh(vars,i,static_cast<int>(vars.size()));
		int count=0;
	
		for (vector<Mpoint*>::iterator k = m._points.begin(); k!=m._points.end(); k++ )
		{
			Pt vertNew;
			Pt vert = (*k)->get_coord();
			Vec normal;
			normal = (*k)->local_normal();
			float avgout=0.0;
			float avgin=0.0;
			for (int j=0;j<ipp;j++){
				int sc=j-(ipp-1)/2;
				vertNew=vert + normal*0.5*sc;
				//calculate probability not residual
								
				if (sc<0){
				  avgin+=image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim);
				}else if(sc>0){
				  avgout+=image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim);
				}
				
			
			}
			gradtot+=abs(avgout-avgin);
			count++;
		}
	}//end of shape interation, all dif in dif
	 //note this is only good if all arjoint
	//cout<<"get Conds"<<endl;
	cost=-gradtot;
return cost;

}



float median(vector<float> vdists){
	int N=vdists.size();
	float median,Q1,Q3;
	if ((N % 2) ==0){
		median=(vdists.at(N/2-1)+vdists.at(N/2))/2;
		if ((N % 4 )==0){
			Q1= (vdists.at(N/4-1)+vdists.at(N/4))/2;
			Q3= (vdists.at(3*N/4-1)+vdists.at(3*N/4))/2;
		}
		else{
			Q1= vdists.at(static_cast<int>(floor(N/4.0)));
			Q3= vdists.at(static_cast<int>(floor(3*N/4.0)));
			
		}
	}else{
		median=vdists.at(static_cast<int>(floor(N/2.0)));
		if (((N-1) % 4 )==0){
			Q1= (vdists.at((N-1)/4-1)+vdists.at((N-1)/4))/2;
			Q3= (vdists.at(3*(N-1)/4)+vdists.at(3*(N-1)/4+1))/2;
		}else{
			Q1= vdists.at(static_cast<int>(floor((N-1)/4.0)));
			Q3= vdists.at(static_cast<int>(ceil(3*(N-1)/4.0)));
			
		}
		
		
	}
	return median;
}	
float mode(vector<float> vdists){
	int N=static_cast<int>(vdists.size());
	int maxcount=0;
	float maxlowint=0;
	float bins=256;
	bins=128;
	float binwidth=(vdists.at(N-1)-0)/bins;
//	cout<<"binwidth"<< binwidth<<" "<<vdists.at(N-1)<<" "<<vdists.at(0)<<endl;
		int count=0;	
		int bincount=1;
			float lowint=0;
	for (int i=0;i<N;i++){
	
	
		if (vdists.at(i)<lowint+binwidth){
			count++;
			if (count>maxcount){
				maxcount=count;
				maxlowint=lowint;
			}
		}else{
			//cout<<count<<" "<<lowint+bincount*binwidth<<" "<<lowint<<endl;
			count=0;
			i--;
			bincount++;
			lowint=bincount*binwidth;
		}
	
	}
	
	//cout<<"bincount "<<bincount<<endl;
	return (maxlowint+binwidth/2.0);
}
float mode(vector<float> vdists, float min, float max){
	int N=static_cast<int>(vdists.size());
	float bins=128;

	float binwidth=(max-min)/bins;
	vector<float> bincounts;
	//innitialize bincounts to zero
	for (int b=0;b<bins;b++){
			bincounts.push_back(0);
		}

	for (int i=0;i<N;i++){
		//search thgrough each bin
		for (int b=0;b<bins;b++){
			//cout<<vdists.at(i)<<endl;
			if (vdists.at(i)<(min+(b+1)*binwidth)){
				//cout<<b<<" "<<min+b*binwidth<<endl;
				bincounts.at(b)++;
				break;
			}
		}
	}
	
	//search for max bin count
	int maxcount=0;
	int maxind=0;
	for (int b=0;b<bins;b++){
		//cout<<bincounts.at(b)<<endl;
			if (bincounts.at(b)>maxcount){
				maxcount=bincounts.at(b);
				maxind=b;
			}
	
	}
	
	return (min+maxind*binwidth+binwidth/2.0);
}

float mode(vector<float> vdists,int *maxcount){
	int N=static_cast<int>(vdists.size());
	*maxcount=0;
	float maxlowint=0;
	float bins=256;
	//bins=128;
	float binwidth=(vdists.at(N-1)-0)/bins;
//	cout<<"binwidth"<< binwidth<<" "<<vdists.at(N-1)<<" "<<vdists.at(0)<<endl;
		int count=0;	
		int bincount=1;
			float lowint=0;
	for (int i=0;i<N;i++){
	
	
		if (vdists.at(i)<lowint+binwidth){
			count++;
			if (count>(*maxcount)){
				*maxcount=count;
				maxlowint=lowint;
			}
		}else{
			//cout<<count<<" "<<lowint+bincount*binwidth<<" "<<lowint<<endl;
			count=0;
			i--;
			bincount++;
			lowint=bincount*binwidth;
		}
	
	}
	
	//cout<<"bincount "<<bincount<<endl;
	return (maxlowint+binwidth/2.0);
}	
float fullwidthhalfmax(vector<float> vdists,float halfmaxval, float *halfmin,float *halfmax){
	int N=static_cast<int>(vdists.size());
	int maxcount=0;
	float maxlowint=0;
	float bins=256;
	bins=128;
	float binwidth=(vdists.at(N-1)-0)/bins;
	int countprev=0;
	
	bool foundmin=false;
	bool foundmax=false;
	
//	cout<<"binwidth"<< binwidth<<" "<<vdists.at(N-1)<<" "<<vdists.at(0)<<endl;
		int count=0;	
		int bincount=1;
			float lowint=0;
	for (int i=0;i<N;i++){
	
	
		if (vdists.at(i)<lowint+binwidth){
			count++;
			if (count>(maxcount)){
				maxcount=count;
				maxlowint=lowint;
			}
		}else{
			//cout<<count<<" "<<lowint+bincount*binwidth<<" "<<lowint<<endl;
			countprev=count;
			if ((count>halfmaxval)&&(!foundmin)){
				*halfmin=lowint+binwidth/2.0;
				foundmin=true;
			}
			if ((count>halfmaxval)){
				*halfmax=lowint+binwidth/2.0;
			}
			count=0;
			i--;
			bincount++;
			lowint=bincount*binwidth;
		}
	
	}
	
	//cout<<"bincount "<<bincount<<endl;
	return (maxlowint+binwidth/2.0);
}	
	//meidan cost func
float costfuncOverlap(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode){
	//cout<<"begin cost"<<endl;
	volume<short> mask;
	volume<short> maskotl;
	
	double cost=0;
	int bounds[6]={0,0,0,0,0,0};
	
	float betas=0;
	
	ColumnVector Best(4);
	vector<float> iprof;
	//new cost
	vector<float> dif;
	
	for (int i=0;i<model1->getNumberOfShapes()-ex;i++){
	  float tp=0,fp=0,fn=0;
		
		betas=0;
		//	cout<<"get deformed mesh "<<vars.at(0)<<" "<<vars.at(1)<<" "<<vars.size()<<endl;
		Mesh m = model1->getDeformedMesh(vars,i,static_cast<int>(vars.size()));
		
		//cout<<"make mask"<<endl;
		mask=make_mask_from_mesh(image,m,model1->getLabel(i),bounds);
		//cout<<gbounds[0]<<" "<<gbounds[1]<<" "<<gbounds[2]<<" "<<gbounds[3]<<" "<<gbounds[4]<<" "<<gbounds[5]<<endl;
		cost+=model1->volumeDistance(&mask,&image,gbounds,&m);
		bool dice=false;
		if (dice){
			//don't wnat to sue dice for now
			for (int x=bounds[0];x<bounds[1];x++){
				for (int y=bounds[2];y<bounds[3];y++){
					for (int z=bounds[4];z<bounds[5];z++){
						if ((mask.value(x,y,z)>0)&&(image.value(x,y,z)>0)){
							tp++;
						}else if ((mask.value(x,y,z)>0)&&(image.value(x,y,z)==0)){
							fp++;
						}else if ((mask.value(x,y,z)==0)&&(image.value(x,y,z)>0)){
							fn++;
						}
						
					}
				}
			}
			cost-=2*tp/(2*tp+fp+fn);
		}
	}
		//cout<<"dice "<<cost<<endl;
			//	cout<<"fd "<<cost<<endl;

	return cost;
}
//float costfuncAppPosterior(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode){

float costfuncApp(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode, float premean){
//cout<<"enter costApp"<<endl;
	//if premean equals -333 then estimate it!
	
	
	//cout<<"begin cost"<<endl;
	//MJ res min vars
	//	float gC1=1/290.0;
	//	float gpreBeta=1/208.018;
	
	
	
	volume<short> mask;
	volume<short> maskotl;
	//float avgres=0;
	//float cost=0;
	double cost=0;
	//double cost2=0;
	int bounds[6]={0,0,0,0,0,0};
	
	float betas=0;
	
	ColumnVector Best(4);
	vector<float> iprof;
	//new cost
	vector<float> dif;
	
	//this is used to store difference from global mean
	vector<float> iprof2;
	//new cost
	vector<float> dif2;
	for (int i=0;i<model1->getNumberOfShapes()-ex;i++){
		float mean=premean;
		Mesh m = model1->getDeformedMesh(vars,i,static_cast<int>(vars.size()));
		if ((premean==-333)&&(GmanMean==-777)&&((!globfoundmode)|(globReest))){
			
	//		float si=0.0;
	//		float ssi=0.0;//this is used for mj cost
						  //added in ssi to look at residual variance for beta estimate
				
	//			int numpix=0;
				betas=0;
				//		cout<<"get deformed mesh "<<vars.size()<<endl;
				
				
				//cout<<"make mask"<<endl;
				mask=make_mask_from_meshInOut(image,m,model1->getLabel(i),bounds);
				//have change away from inout (was used in original build
				//mask=make_mask_from_mesh(image,m,model1->getLabel(i),bounds);
		//		for (int x=bounds[0];x<bounds[1];x++){
		//			for (int y=bounds[2];y<bounds[3];y++){
		//				for (int z=bounds[4];z<bounds[5];z++){
		//					if ((mask.value(x,y,z)>0)){
		//						float tempi;
		//						tempi=image.value(x,y,z);
		//						si+=tempi;
		//						ssi+=tempi*tempi;
		//						numpix+=1;
		//					}
		//				}
		//			}
		//		}
		//		
				//calculate intensity histogram
				float maxint=0,minint=0;
				vector<float> vintens;
				model1->intensityHistMaxMin(&image,&mask,&m,model1->getLabel(i),&vintens, &maxint, &minint);
			
				mean=mode(vintens,minint,maxint);
					if (globMean==mean){
					cout<<"FOund mode"<<endl;
						globfoundmode=true;
					}
				cout<<"mode "<<mean<<endl;
			
				globMean=mean;
			//	cout<<"Fill Mesh and set globMEan "<<globMean<<endl;
			}else if (manMean.value()!=-777){
				globMean=GmanMean;
				//cout<<"set Mean to: "<<manMean.value()<<" "<<GmanMean<<endl;
			}
			
			iprof=model1->getDeformedIprof(vars,i,vars.size());
				vector<float> varstemp;		
			iprof2=model1->getDeformedIprof(varstemp,i,varstemp.size());

			int ipp=model1->getIPP(0);
			
			int count=0;
			for (vector<Mpoint*>::iterator k = m._points.begin(); k!=m._points.end(); k++ )
			{
				Pt vertNew;
				Pt vert = (*k)->get_coord();
				Vec normal;
				normal = (*k)->local_normal();
				for (int j=0;j<ipp;j++){
					int sc=j-(ipp-1)/2;
					vertNew=vert + normal*0.5*sc;
					dif.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-globMean -(iprof.at(count*ipp+j)));
					dif2.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-globMean -(iprof2.at(count*ipp+j)));
				}
				count++;
				//cout<<count<<endl;
			}
			}//end of shape interation, all dif in dif
			 //note this is only good if all arjoint
	
	///%Calculate conditional I | s
	
	vector< vector<float> > prec=model1->getShape(0)->getICondPrec();
	vector<float> eigs=model1->getShape(0)->getICondEigs();
	double probIcond=0;
	for (unsigned int col=0;col<prec.size();col++){
		
		float multemp=0;
		for (unsigned int row=0;row<dif.size();row++){
			multemp+=dif.at(row)*prec.at(col).at(row);
			
		}
		multemp*=multemp;
		probIcond+=multemp/eigs.at(col);
	}
	//the multiplication by n-1 or n is left out becomes constant in log cost
	//probI*=M;//this multiplication is performed later
	//calculates inner product of difference between observed intensity and mean 
	//this is todo 
	float sdif=0;
	for (unsigned int row=0;row<dif.size();row++){
		sdif+=dif.at(row)*dif.at(row);
	}
	//work sonly if all are tied together
	float M=model1->getNumberOfSubjects();
	//sdif*=M/(model1->getShape(0)->getErrs().at(1)*2);
	//the M multiplication is performed later
	sdif*=1/(model1->getShape(0)->getErrs().at(1)*2);
	probIcond+=sdif;
	
	//////%%%%%%%%%%%%%%%%%%%%%%%%Claculate p(I)
	
	
	prec=model1->getShape(0)->getBCondPrec();
	eigs=model1->getShape(0)->getBCondEigs();
	//float probI=0;
	double probI=0;
	//cout<<"caLC1"<<endl;
	for (unsigned int col=0;col<prec.size();col++){
		
		float multemp=0;
		for (unsigned int row=0;row<dif.size();row++){
			multemp+=dif2.at(row)*prec.at(col).at(row);
			
		}
		multemp*=multemp;
		probI+=multemp/eigs.at(col);
	}
	//cout<<"caLC2"<<endl;
	//the multiplication by n-1 or n is left out becomes constant in log cost
	//probI*=M;
	//calculates inner product of difference between observed intensity and mean 
	sdif=0;
	for (unsigned int row=0;row<dif.size();row++){
		sdif+=dif2.at(row)*dif2.at(row);
	}
	//work sonly if all are tied together
	//multiplication by (M-1) is taken care of below		
	sdif*=1/(model1->getShape(0)->getErrs().at(1)*2);
	probI+=sdif;
	
	/////%%%%%%%%%%%%%%%%%%%%%%%%% Now claculate costfunction
	
	
	
	int cumnum=0;
	for (int q=0;q<model1->getNumberOfShapes();q++){
		cumnum+=model1->getNumberOfPoints(q);
	}
	
	double alp=M-1.0/M;
	double binner=0;
	double gammav=alp/(alp-2);
	float k2=3.0*static_cast<float>(cumnum);
	float k1=static_cast<float>(dif.size());
	
	for (unsigned int i=0;i<vars.size();i++){
		binner+=vars.at(i)*vars.at(i);
	}
	
	//much greater than 1 simplification
	//cost+= (alp+k2)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(probI);//-(12/2-1)*log(avgres) +avgres/2;
	cost+= -(k1)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(1+(M-1)*probIcond/(alp+binner*(gammav)));//-(12/2-1)*log(avgres) +avgres/2;
		
		//this is the shaoe prior ...chooeses between no condition, one conditional, or 2 conditionals
		//if (globModePriorMean==-777){
		//cout<<"got to shape prior"<<endl;
		if (true){
			if ((shcond.value())&&(GlobShCond))			{
				
				ColumnVector bx1temp(vars.size());
				for (unsigned int i=0;i<vars.size();i++){
					bx1temp.element(i)=vars.at(i)-mBx2map.element(i);
				}
				//mBx2map=mBx2*Bx2;
				ColumnVector Bcx1=mBx1inv*bx1temp;
				//	mBx1inv=mBcx1.i();
				Matrix Bcinner=Bcx1.t()*Bcx1;
				cost+=(alp+k2+kpred)/2*log(1+Bcinner.element(0,0)/(alp+kpred));
				//and in posteriro bit
				
				//add a secojnd conditional
				if (shcond2.value()){
					for (unsigned int i=0;i<vars.size();i++){
						bx1temp.element(i)=vars.at(i)-mBx2map2.element(i);
					}
					//mBx2map=mBx2*Bx2;
					Bcx1=mBx1inv2*bx1temp;
					//	mBx1inv=mBcx1.i();
					Bcinner=Bcx1.t()*Bcx1;
					cost+=(alp+k2+kpred2)/2*log(1+Bcinner.element(0,0)/(alp+kpred2));
					
				}
			}else{
				
				//p(x) prior--shape prior
				cost+=(alp+k2)/2*log(1+binner*gammav/alp);
			}
		}
		//this bit enables a mode priro
		if (globModePriorMean!=-777){
			cout<<"use prior"<<endl;
			cost-= (M/2)*log(1+(globMean-globModePriorMean)*(globMean-globModePriorMean)/globModePriorVar/(M-1));
			
					}
									//This is the denominator in the posteriro
		//subtract p(i)
	//double negatuve becomes psotivtie
		cost+=(alp+k1)/2*log(1+(M-1)*probI/alp);
		//cout<<"leave cost app"<<endl;
		return cost;
		
	}
	float costfuncAppJ(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode, float premean){
//cout<<"enter costApp"<<endl;
	//if premean equals -333 then estimate it!
	
	
	//cout<<"begin cost"<<endl;
	//MJ res min vars
	//	float gC1=1/290.0;
	//	float gpreBeta=1/208.018;
	
	
	
	volume<short> mask;
	volume<short> maskotl;
	//float avgres=0;
	//float cost=0;
	double cost=0;
	//double cost2=0;
	int bounds[6]={0,0,0,0,0,0};
	
	float betas=0;
	
	ColumnVector Best(4);
	vector<float> iprof;
	//new cost
	vector<float> dif;
	
	//this is used to store difference from global mean
	vector<float> iprof2;
	//new cost
	vector<float> dif2;
	for (int i=0;i<model1->getNumberOfShapes()-ex;i++){
		float mean=premean;
		Mesh m = model1->getDeformedMesh(vars,i,static_cast<int>(vars.size()));
		if ((premean==-333)&&(GmanMean==-777)&&((!globfoundmode)|(globReest))){
			
	//		float si=0.0;
	//		float ssi=0.0;//this is used for mj cost
						  //added in ssi to look at residual variance for beta estimate
				
	//			int numpix=0;
				betas=0;
				//		cout<<"get deformed mesh "<<vars.size()<<endl;
				
				
				//cout<<"make mask"<<endl;
				mask=make_mask_from_meshInOut(image,m,model1->getLabel(i),bounds);
				//have change away from inout (was used in original build
				//mask=make_mask_from_mesh(image,m,model1->getLabel(i),bounds);
		//		for (int x=bounds[0];x<bounds[1];x++){
		//			for (int y=bounds[2];y<bounds[3];y++){
		//				for (int z=bounds[4];z<bounds[5];z++){
		//					if ((mask.value(x,y,z)>0)){
		//						float tempi;
		//						tempi=image.value(x,y,z);
		//						si+=tempi;
		//						ssi+=tempi*tempi;
		//						numpix+=1;
		//					}
		//				}
		//			}
		//		}
		//		
				//calculate intensity histogram
				float maxint=0,minint=0;
				vector<float> vintens;
				model1->intensityHistMaxMin(&image,&mask,&m,model1->getLabel(i),&vintens, &maxint, &minint);
			
				mean=mode(vintens,minint,maxint);
					if (globMean==mean){
					cout<<"FOund mode"<<endl;
						globfoundmode=true;
					}
				cout<<"mode "<<mean<<endl;
			
				globMean=mean;
			//	cout<<"Fill Mesh and set globMEan "<<globMean<<endl;
			}else if (manMean.value()!=-777){
				globMean=GmanMean;
				//cout<<"set Mean to: "<<manMean.value()<<" "<<GmanMean<<endl;
			}
			
			iprof=model1->getDeformedIprof(vars,i,vars.size());
	//			vector<float> varstemp;	
	//			varstemp.push_back(0);	
	//		iprof2=model1->getDeformedIprof(varstemp,i,varstemp.size());

			int ipp=model1->getIPP(0);
			
			int count=0;
			for (vector<Mpoint*>::iterator k = m._points.begin(); k!=m._points.end(); k++ )
			{
				Pt vertNew;
				Pt vert = (*k)->get_coord();
				Vec normal;
				normal = (*k)->local_normal();
				for (int j=0;j<ipp;j++){
					int sc=j-(ipp-1)/2;
					vertNew=vert + normal*0.5*sc;
					dif.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-globMean -(iprof.at(count*ipp+j)));
		//			dif2.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-globMean -(iprof2.at(count*ipp+j)));
				}
				count++;
				//cout<<count<<endl;
			}
			}//end of shape interation, all dif in dif
			 //note this is only good if all arjoint
	
	///%Calculate conditional I | s
	
	vector< vector<float> > prec=model1->getShape(0)->getICondPrec();
	vector<float> eigs=model1->getShape(0)->getICondEigs();
	double probIcond=0;
	for (unsigned int col=0;col<prec.size();col++){
		
		float multemp=0;
		for (unsigned int row=0;row<dif.size();row++){
			multemp+=dif.at(row)*prec.at(col).at(row);
			
		}
		multemp*=multemp;
		probIcond+=multemp/eigs.at(col);
	}
	//the multiplication by n-1 or n is left out becomes constant in log cost
	//probI*=M;//this multiplication is performed later
	//calculates inner product of difference between observed intensity and mean 
	//this is todo 
	float sdif=0;
	for (unsigned int row=0;row<dif.size();row++){
		sdif+=dif.at(row)*dif.at(row);
	}
	//work sonly if all are tied together
	float M=model1->getNumberOfSubjects();
	//sdif*=M/(model1->getShape(0)->getErrs().at(1)*2);
	//the M multiplication is performed later
	sdif*=1/(model1->getShape(0)->getErrs().at(1)*2);
	probIcond+=sdif;
	
	//////%%%%%%%%%%%%%%%%%%%%%%%%Claculate p(I)
	
/*	
	prec=model1->getShape(0)->getBCondPrec();
	eigs=model1->getShape(0)->getBCondEigs();
	//float probI=0;
	double probI=0;
	//cout<<"caLC1"<<endl;
	for (unsigned int col=0;col<prec.size();col++){
		
		float multemp=0;
		for (unsigned int row=0;row<dif.size();row++){
			multemp+=dif2.at(row)*prec.at(col).at(row);
			
		}
		multemp*=multemp;
		probI+=multemp/eigs.at(col);
	}
	//cout<<"caLC2"<<endl;
	//the multiplication by n-1 or n is left out becomes constant in log cost
	//probI*=M;
	//calculates inner product of difference between observed intensity and mean 
	sdif=0;
	for (unsigned int row=0;row<dif.size();row++){
		sdif+=dif2.at(row)*dif2.at(row);
	}
	//work sonly if all are tied together
	//multiplication by (M-1) is taken care of below		
	sdif*=1/(model1->getShape(0)->getErrs().at(1)*2);
	probI+=sdif;
*/	
	/////%%%%%%%%%%%%%%%%%%%%%%%%% Now claculate costfunction
	
	
	
	int cumnum=0;
	for (int q=0;q<model1->getNumberOfShapes();q++){
		cumnum+=model1->getNumberOfPoints(q);
	}
	
	double alp=M-1.0/M;
	double binner=0;
	double gammav=alp/(alp-2);
	float k2=3.0*static_cast<float>(cumnum);
	float k1=static_cast<float>(dif.size());
	
	for (unsigned int i=0;i<vars.size();i++){
		binner+=vars.at(i)*vars.at(i);
	}
	
	//much greater than 1 simplification
	//cost+= (alp+k2)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(probI);//-(12/2-1)*log(avgres) +avgres/2;
	cost+= -(k1)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(1+(M-1)*probIcond/(alp+binner*(gammav)));//-(12/2-1)*log(avgres) +avgres/2;
		
		//this is the shaoe prior ...chooeses between no condition, one conditional, or 2 conditionals
		//if (globModePriorMean==-777){
		//cout<<"got to shape prior"<<endl;
		if (true){
			if ((shcond.value())&&(GlobShCond))			{
				
				ColumnVector bx1temp(vars.size());
				for (unsigned int i=0;i<vars.size();i++){
					bx1temp.element(i)=vars.at(i)-mBx2map.element(i);
				}
				//mBx2map=mBx2*Bx2;
				ColumnVector Bcx1=mBx1inv*bx1temp;
				//	mBx1inv=mBcx1.i();
				Matrix Bcinner=Bcx1.t()*Bcx1;
				cost+=(alp+k2+kpred)/2*log(1+Bcinner.element(0,0)/(alp+kpred));
				//and in posteriro bit
				
				//add a secojnd conditional
				if (shcond2.value()){
					for (unsigned int i=0;i<vars.size();i++){
						bx1temp.element(i)=vars.at(i)-mBx2map2.element(i);
					}
					//mBx2map=mBx2*Bx2;
					Bcx1=mBx1inv2*bx1temp;
					//	mBx1inv=mBcx1.i();
					Bcinner=Bcx1.t()*Bcx1;
					cost+=(alp+k2+kpred2)/2*log(1+Bcinner.element(0,0)/(alp+kpred2));
					
				}
			}else{
				
				//p(x) prior--shape prior
				cost+=(alp+k2)/2*log(1+binner*gammav/alp);
			}
		}
		//this bit enables a mode priro
		if (globModePriorMean!=-777){
			cout<<"use prior"<<endl;
			cost-= (M/2)*log(1+(globMean-globModePriorMean)*(globMean-globModePriorMean)/globModePriorVar/(M-1));
			
					}
									//This is the denominator in the posteriro
		//subtract p(i)
	//	cost-=(alp+k1)/2*log(1+(M-1)*probI/alp);
		//cout<<"leave cost app"<<endl;
		return cost;
		
	}

	
	
	
	float costfuncAppAff(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode, float premean){

	//if premean equals -333 then estimate it!
	
	
//cout<<"begin cost aff"<<endl;
	//MJ res min vars
	//	float gC1=1/290.0;
	//	float gpreBeta=1/208.018;
	
	
	
	volume<short> mask;
	volume<short> maskotl;
	//float avgres=0;
	//float cost=0;
	double cost=0;
	//double cost2=0;
	int bounds[6]={0,0,0,0,0,0};
	
	float betas=0;
	
	ColumnVector Best(4);
	vector<float> iprof;
		vector<float> iprofAff;

		vector<float> iprof2;

	//new cost
	vector<float> dif;
		vector<float> difAff;

	vector<float> dif2;
	for (int i=0;i<model1->getNumberOfShapes()-ex;i++){
		float mean=premean;
		Mesh m = model1->getDeformedMeshAff7(vars,i,static_cast<int>(vars.size()));
		m.save("afftest.vtk",3);
		if ((premean==-333)&&(GmanMean==-777)&&((!globfoundmode)|(globReest))){
			
	//		float si=0.0;
	//		float ssi=0.0;//this is used for mj cost
						  //added in ssi to look at residual variance for beta estimate
				
	//			int numpix=0;
				betas=0;
				//		cout<<"get deformed mesh "<<vars.size()<<endl;
				
				
			//	cout<<"make mask"<<endl;
				mask=make_mask_from_meshInOut(image,m,model1->getLabel(i),bounds);
				//have change away from inout (was used in original build
				//mask=make_mask_from_mesh(image,m,model1->getLabel(i),bounds);
		//		for (int x=bounds[0];x<bounds[1];x++){
		//			for (int y=bounds[2];y<bounds[3];y++){
		//				for (int z=bounds[4];z<bounds[5];z++){
		//					if ((mask.value(x,y,z)>0)){
		//						float tempi;
		//						tempi=image.value(x,y,z);
		//						si+=tempi;
		//						ssi+=tempi*tempi;
		//						numpix+=1;
		//					}
		//				}
		//			}
		//		}
		//		
				//calculate intensity histogram
			//		cout<<"make histograM"<<endl;
				vector<float> vintens;
				model1->intensityHist(&image,&mask,&m,model1->getLabel(i),&vintens);
			//	cout<<"done hist"<<endl;
				//calculate mixture model
				if (histmode.value()){
			//	cout<<"calc mode"<<endl;
				mean=mode(vintens);
			//		cout<<"calced mode"<<endl;
					if (globMean==mean){
					cout<<"FOund mode"<<endl;
						globfoundmode=true;
					}
				cout<<"mode "<<mean<<endl;
				}else{
				mean=model1->EMgmm(&vintens,gmmLesser.value(),EMiter.value());
				}
				globMean=mean;
			//	cout<<"Fill Mesh and set globMEan "<<globMean<<endl;
			}else if (manMean.value()!=-777){
				globMean=GmanMean;
				//cout<<"set Mean to: "<<manMean.value()<<" "<<GmanMean<<endl;
			}
			//cout<<"set up vector"<<endl;
			vector<float> varsSh;
			for (unsigned int k=7;k<vars.size();k++){
				varsSh.push_back(vars.at(k));
			}
		//		cout<<"begin i deform"<<endl;
			iprof=model1->getDeformedIprof(varsSh,i,varsSh.size());
			vector<float> varsAff;
			for (unsigned int k=0;k<7;k++){
				varsAff.push_back(vars.at(k));
			}
			
			iprofAff=model1->getDeformedIprof(varsAff,i,varsAff.size());
			//			iprof=model1->getDeformedIprof(vars,i,vars.size());

		//	cout<<"done i deform"<<endl;
			int ipp=model1->getIPP(0);
			vector<float> varstemp;		
			iprof2=model1->getDeformedIprof(varstemp,i,varstemp.size());
			int count=0;
			
			for (vector<Mpoint*>::iterator k = m._points.begin(); k!=m._points.end(); k++ )
			{
				Pt vertNew;
				Pt vert = (*k)->get_coord();
				Vec normal;
				normal = (*k)->local_normal();
				for (int j=0;j<ipp;j++){
					int sc=j-(ipp-1)/2;
					vertNew=vert + normal*0.5*sc;
					dif.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-globMean -(iprof.at(count*ipp+j)));
					difAff.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-globMean -(iprofAff.at(count*ipp+j)));

					dif2.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-globMean -(iprof2.at(count*ipp+j)));

				}
				count++;
			}
			
			}//end of shape interation, all dif in dif
			 //note this is only good if all arjoint
	
	///%Calculate conditional I | s
	
	vector< vector<float> > prec=model1->getShape(0)->getICondPrec();
	vector<float> eigs=model1->getShape(0)->getICondEigs();
	double probIcond=0;
	for (unsigned int col=0;col<prec.size();col++){
		
		float multemp=0;
		for (unsigned int row=0;row<dif.size();row++){
			multemp+=dif.at(row)*prec.at(col).at(row);
			
		}
		multemp*=multemp;
		probIcond+=multemp/eigs.at(col);
	}
	//the multiplication by n-1 or n is left out becomes constant in log cost
	//probI*=M;//this multiplication is performed later
	//calculates inner product of difference between observed intensity and mean 
	float sdif=0;
	for (unsigned int row=0;row<dif.size();row++){
		sdif+=dif.at(row)*dif.at(row);
	}
//	cout<<"SDIF "<<sdif<<endl;
	//work sonly if all are tied together
	float M=model1->getNumberOfSubjects();
	//sdif*=M/(model1->getShape(0)->getErrs().at(1)*2);
	//the M multiplication is performed later
	sdif*=1/(model1->getShape(0)->getErrs().at(1)*2);
	probIcond+=sdif;
	//	cout<<"COST CPOST icond"<<probIcond<<endl;
	///%Calculate conditional I | s
	//***************Calculate I | aff********************//
	 prec=model1->getShape(0)->getICondPrec();
	eigs=model1->getShape(0)->getICondEigs();
	double probIcondAff=0;
	for (unsigned int col=0;col<prec.size();col++){
		
		float multemp=0;
		for (unsigned int row=0;row<difAff.size();row++){
			multemp+=difAff.at(row)*prec.at(col).at(row);
			
		}
		multemp*=multemp;
		probIcondAff+=multemp/eigs.at(col);
	}
	//the multiplication by n-1 or n is left out becomes constant in log cost
	//probI*=M;//this multiplication is performed later
	//calculates inner product of difference between observed intensity and mean 
	sdif=0;
	for (unsigned int row=0;row<difAff.size();row++){
		sdif+=difAff.at(row)*difAff.at(row);
	}
//	cout<<"SDIF "<<sdif<<endl;
	//work sonly if all are tied together
	//float M=model1->getNumberOfSubjects();
	//sdif*=M/(model1->getShape(0)->getErrs().at(1)*2);
	//the M multiplication is performed later
	sdif*=1/(model1->getShape(0)->getErrs().at(1)*2);
	probIcondAff+=sdif;
	//	cout<<"COST CPOST icond"<<probIcond<<endl
	//////%%%%%%%%%%%%%%%%%%%%%%%%Claculate p(I)
//	
//	
	prec=model1->getShape(0)->getBCondPrec();
	eigs=model1->getShape(0)->getBCondEigs();
	//float probI=0;
	double probI=0;
//	//cout<<"caLC1"<<endl;
	for (unsigned int col=0;col<prec.size();col++){
		
		float multemp=0;
		for (unsigned int row=0;row<dif.size();row++){
			multemp+=dif2.at(row)*prec.at(col).at(row);
			
		}
		multemp*=multemp;
		probI+=multemp/eigs.at(col);
	}
//		cout<<"COST CPOST 1"<<probI<<endl;

	//cout<<"caLC2"<<endl;
	//the multiplication by n-1 or n is left out becomes constant in log cost
	//probI*=M;
	//calculates inner product of difference between observed intensity and mean 
	sdif=0;
	for (unsigned int row=0;row<dif.size();row++){
		sdif+=dif2.at(row)*dif2.at(row);
	}
	//work sonly if all are tied together
	//multiplication by (M-1) is taken care of below		
	sdif*=1/(model1->getShape(0)->getErrs().at(1)*2);
	probI+=sdif;
//	cout<<"COST CPOST 2"<<probI<<endl;
	/////%%%%%%%%%%%%%%%%%%%%%%%%% Now claculate costfunction
	
	
	
	int cumnum=0;
	for (int q=0;q<model1->getNumberOfShapes();q++){
		cumnum+=model1->getNumberOfPoints(q);
	}
	
	double alp=M-1.0/M;
	double binner=0, binnerAff=0;
	double gammav=alp/(alp-2);
	float k2=3.0*static_cast<float>(cumnum);
	float k1=static_cast<float>(dif.size());
	
	//this part changes for Aff....omit first 7 vars
	for (unsigned int i=0;i<7;i++){
		binnerAff+=vars.at(i)*vars.at(i);
	}
	//this part changes for Aff....omit first 7 vars
	for (unsigned int i=7;i<vars.size();i++){
		binner+=vars.at(i)*vars.at(i);
	}
			//	cout<<"COST CPOST "<< -(k1)/2*log((alp+k2)/(alp+binner*(gammav)))<<" "<<(alp+k1+k2)/2*log(1+(M-1)*probIcond/(alp+binner*(gammav)))<<" "<<probIcond<<endl;

	//much greater than 1 simplification
	//cost+= (alp+k2)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(probI);//-(12/2-1)*log(avgres) +avgres/2;
	cost+= -(k1)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(1+(M-1)*probIcond/(alp+binner*(gammav)));//-(12/2-1)*log(avgres) +avgres/2;
	//add aff component
		cost+=(alp+k1+7)/2*log(1+(M-1)*probIcondAff/(alp+binnerAff*(gammav)));
		//this is the shaoe prior ...chooeses between no condition, one conditional, or 2 conditionals
		if (shcond.value()){
			
			ColumnVector bx1temp(vars.size());
			for (unsigned int i=0;i<vars.size();i++){
				bx1temp.element(i)=vars.at(i)-mBx2map.element(i);
			}
			//mBx2map=mBx2*Bx2;
			ColumnVector Bcx1=mBx1inv*bx1temp;
			//	mBx1inv=mBcx1.i();
			Matrix Bcinner=Bcx1.t()*Bcx1;
			cost+=(alp+k2+kpred)/2*log(1+Bcinner.element(0,0)/(alp+kpred));
			//and in posteriro bit
			
			//add a secojnd conditional
			if (shcond2.value()){
				for (unsigned int i=0;i<vars.size();i++){
					bx1temp.element(i)=vars.at(i)-mBx2map2.element(i);
				}
				//mBx2map=mBx2*Bx2;
				Bcx1=mBx1inv2*bx1temp;
				//	mBx1inv=mBcx1.i();
				Bcinner=Bcx1.t()*Bcx1;
				cost+=(alp+k2+kpred2)/2*log(1+Bcinner.element(0,0)/(alp+kpred2));
				
			}
		}else{
		//	cout<<"made it to shape prior"<<endl;
			//p(x) prior--shape prior
			cost+=(alp+k2)/2*log(1+binner*gammav/alp);
		//	cout<<"COST CPOST "<<cost<<endl;
			//need to add in affine components...independent translation, rotation , and sc
			//strat with trans
			float btr=0, brot=0;
			for (unsigned int i=0;i<3;i++){
				btr+=vars.at(i)*vars.at(i);
				brot+=vars.at(i+3)*vars.at(i+3);
			}
			float bsc=vars.at(6)*vars.at(6);
			//gammaTrRot=(M-3)/(M-3-2);
			//normalization becomes alp-2 from gammav/alp
			//cout<<"cost "<<M/2.0*log(1+btr/(M-3-2))<<" "<<M/2.0*log(1+brot/(M-3-2))<<" "<<M/2.0*log(1+bsc/(M-1-2))<<endl;
			cost+=M/2.0*log(1+btr/(M-3-2));
			cost+=M/2.0*log(1+brot/(M-3-2));
			cost+=M/2.0*log(1+bsc/(M-1-2));

			
		}
		
		//This is the denominator in the posteriro
		//subtract p(i)
		cost+=2*(alp+k1)/2*log(1+(M-1)*probI/alp);
		
		//cout<<"exit cost aff"<<endl;
		return cost;
		
	}


float costfuncAppDual(vector< shapeModel* > Vmodel1, vector< vector<float> > Vvars, vector<float> Vpremean){
	//only implemented for the shape model
	///asssumes histmode
	//assumes 2 models sometimes
	//assume 1 structure per model
	
	//if premean equals -333 then estimate it!
	
	
//	cout<<"begin dual cost"<<endl;
	//MJ res min vars
	//	float gC1=1/290.0;
	//	float gpreBeta=1/208.018;
	
	
	
	volume<short> mask;
	volume<short> maskotl;
	//float avgres=0;
	//float cost=0;
	double cost=0;
	//double cost2=0;
	int bounds[6]={0,0,0,0,0,0};
	
	float betas=0;
	
	ColumnVector Best(4);
	vector<float> iprof;
	//new cost
	
	
	int numModels = static_cast<int>(Vmodel1.size());
	//cout<<"costfunc prob1"<<endl;
	//assume a single shape per model
	for (int model=0;model<numModels;model++){
	vector<float> dif;
		//float mean=premean;
//			cout<<"model number "<<model<<endl;
			//cout<<Vpremean.at(model)<<" "<<Vmodel1.at(model)->getLabel(0)<<endl;
		Mesh m = Vmodel1.at(model)->getDeformedMesh(Vvars.at(model),0,static_cast<int>(Vvars.at(model).size()));
//		cout<<"hmmmm costmodel"<<endl;
		//premean
		//	cout<<"costfunc prob4"<<endl;
		if ((Vpremean.at(model)==-333)&&(GmanMean==-777)&&((!Vglobfoundmode.at(model))|(globReest))){
	//		cout<<"costfunc prob3"<<endl;
			betas=0;
//			cout<<"hmmmm make mask"<<endl;
			mask=make_mask_from_meshInOut(image,m,Vmodel1.at(model)->getLabel(0),bounds);
//			cout<<"hmmmm mask made"<<endl;
			vector<float> vintens;
			Vmodel1.at(model)->intensityHist(&image,&mask,&m,Vmodel1.at(model)->getLabel(0),&vintens);
//			cout<<"costfunc prob3.1 "<<vintens.size()<<endl;
			float modeval=mode(vintens);
			Vglobfoundmode.at(model)=true;
//			cout<<"costfunc prob3.2"<<endl;
			VglobMean.at(model)=modeval;
		cout<<"hmmmm find mode "<<modeval<<endl;
		}else if (manMean.value()!=-777){
			VglobMean.at(model)=GmanMean;
			//cout<<"set Mean to: "<<manMean.value()<<" "<<GmanMean<<endl;
		}//else{
			
			//VglobMean=Vpremean;
	//	}
		//	cout<<"costfunc prob2 "<<VglobMean.at(model)<<endl;
		iprof=Vmodel1.at(model)->getDeformedIprof(Vvars.at(model),0,Vvars.at(model).size());
		
		int ipp=Vmodel1.at(model)->getIPP(0);
		
		int count=0;
		
		for (vector<Mpoint*>::iterator k = m._points.begin(); k!=m._points.end(); k++ )
		{
			Pt vertNew;
			Pt vert = (*k)->get_coord();
			Vec normal;
			normal = (*k)->local_normal();
			for (int j=0;j<ipp;j++){
				int sc=j-(ipp-1)/2;
				vertNew=vert + normal*0.5*sc;
				dif.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-VglobMean.at(model) -(iprof.at(count*ipp+j)));
				
			}
			count++;
		}
		//cout<<"hmmmm hmmm"<<endl;
		//end of shape interation, all dif in dif
		//note this is only good if all arjoint
		
		///%Calculate conditional I | s
		
		vector< vector<float> > prec=Vmodel1.at(model)->getShape(0)->getICondPrec();
		vector<float> eigs=Vmodel1.at(model)->getShape(0)->getICondEigs();
		//cout<<"hmmmm hmmm 1.5 "<<dif.size()<<" "<<prec.size()<<" "<<prec.at(0).size()<<endl;
		double probIcond=0;
		for (unsigned int col=0;col<prec.size();col++){
		//	cout<<"col "<<col<<endl;
			float multemp=0;
			for (unsigned int row=0;row<dif.size();row++){
			//cout<<"row "<<row<<endl;
				multemp+=dif.at(row)*prec.at(col).at(row);
				
			}
			multemp*=multemp;
			//cout<<"col2 "<<prec.size()<<" "<<endl;
			probIcond+=multemp/eigs.at(col);
			//cout<<"col3"<<endl;
		}
		//cout<<"hmmmm hmmm 2"<<endl;
		//the multiplication by n-1 or n is left out becomes constant in log cost
		//probI*=M;//this multiplication is performed later
		//calculates inner product of difference between observed intensity and mean 
		float sdif=0;
		for (unsigned int row=0;row<dif.size();row++){
			sdif+=dif.at(row)*dif.at(row);
		}
		//work sonly if all are tied together
		float M=Vmodel1.at(model)->getNumberOfSubjects();
		//sdif*=M/(model1->getShape(0)->getErrs().at(1)*2);
		//the M multiplication is performed later
		sdif*=1/(Vmodel1.at(model)->getShape(0)->getErrs().at(1)*2);
		probIcond+=sdif;
		//	cout<<"hmmmm hmmm P(I)"<<endl;
		//////%%%%%%%%%%%%%%%%%%%%%%%%Claculate p(I)
		
		
		prec=Vmodel1.at(model)->getShape(0)->getBCondPrec();
		eigs=Vmodel1.at(model)->getShape(0)->getBCondEigs();
		//float probI=0;
		double probI=0;
		//cout<<"caLC1"<<endl;
		for (unsigned int col=0;col<prec.size();col++){
			
			float multemp=0;
			for (unsigned int row=0;row<dif.size();row++){
				multemp+=dif.at(row)*prec.at(col).at(row);
				
			}
			multemp*=multemp;
			probI+=multemp/eigs.at(col);
		}
		//cout<<"caLC2"<<endl;
		//the multiplication by n-1 or n is left out becomes constant in log cost
		//probI*=M;
		//calculates inner product of difference between observed intensity and mean 
		sdif=0;
		for (unsigned int row=0;row<dif.size();row++){
			sdif+=dif.at(row)*dif.at(row);
		}
		//work sonly if all are tied together
		//multiplication by (M-1) is taken care of below		
		sdif*=1/(Vmodel1.at(model)->getShape(0)->getErrs().at(1)*2);
		probI+=sdif;
		
		/////%%%%%%%%%%%%%%%%%%%%%%%%% Now claculate costfunction
		
		
		
		int cumnum=0;
		for (int q=0;q<Vmodel1.at(model)->getNumberOfShapes();q++){
			cumnum+=Vmodel1.at(model)->getNumberOfPoints(q);
		}
		
		double alp=M-1.0/M;
		double binner=0;
		double gammav=alp/(alp-2);
		
		float k2=3.0*static_cast<float>(cumnum);
		float k1=static_cast<float>(dif.size());
		
		for (unsigned int i=0;i<Vvars.at(model).size();i++){
			binner+=Vvars.at(model).at(i)*Vvars.at(model).at(i);
		}
		
		//much greater than 1 simplification
		//cost+= (alp+k2)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(probI);//-(12/2-1)*log(avgres) +avgres/2;
		cost+= -(k1)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(1+(M-1)*probIcond/(alp+binner*(gammav)));//-(12/2-1)*log(avgres) +avgres/2;
			//cout<<"cost shape prior"<<endl;
			//this is the shaoe prior ...chooeses between no condition, one conditional, or 2 conditionals
		
			if (model==0){
			
				Matrix mBx2;
				Matrix mBcx1;
		//		cout<<"readbmap "<<bmapname.value()<<endl;
				M=readBmap(bmapname.value(),Vvars.at(model+1),&mBx2,&mBcx1);
		//		cout<<"readbmap done"<<endl;
			//	int nummodes=
				ColumnVector Bx2(Vmodel1.at(model)->getNumberOfModes());
				//load bx2 vars that were read
				
				
				
				//NEED ANOTHER MODEL+1 HAVE TP DOUBLE CHECK
				for (unsigned int i=0;i<Vvars.at(model).size();i++){
					Bx2.element(i)=Vvars.at(model+1).at(i);
				}
				mBx2map=mBx2*Bx2;

				
			
			
				//only implemented for 2 structures.
				ColumnVector bx1temp(Vvars.at(model).size());
				//cout<<bx1temp.Nrows()<<" "<<bx1temp.Ncols()<<endl;
				for (unsigned int i=0;i<Vvars.at(model).size();i++){
					bx1temp.element(i)=Vvars.at(model).at(i)-mBx2map.element(i);
				}
//				cout<<"costfuncAPP hmmm3"<<endl;
				//mBx2map=mBx2*Bx2;
				ColumnVector Bcx1=mBx1inv*bx1temp;
				//	mBx1inv=mBcx1.i();
				Matrix Bcinner=Bcx1.t()*Bcx1;
				cost+=(alp+k2+kpred)/2*log(1+Bcinner.element(0,0)/(alp+kpred));
				//and in posteriro bit
				
			}else{
				
				//p(x) prior--shape prior
				cost+=(alp+k2)/2*log(1+binner*gammav/alp);
			}
				//This is the denominator in the posteriro
		//subtract p(i)
		cost-=(alp+k1)/2*log(1+(M-1)*probI/alp);
		//cout<<"end shape 1"<<endl;
	}
	//cout<<"end dual cost "<<cost<<endl;
	return cost;
	
}
float costfuncAppN(vector< shapeModel* > Vmodel1, vector< vector<float> > Vvars, vector<float> Vpremean, vector<int> vKx2, vector<Matrix> vmBx2,vector<Matrix> vmBx1inv, bool grad, int gradModel){
	//only implemented for the shape model
	///asssumes histmode
	//assume 1 structure per model
	int numModels = static_cast<int>(Vmodel1.size());
	
	cout<<"vvarsCostN in "<<endl;
	for (int model=0;model<numModels;model++){
		
		for (unsigned int i=0;i<Vvars.at(model).size();i++){
			cout<<Vvars.at(model).at(i)<<" ";
			
		}
		cout<<endl;
	}
	
	//if premean equals -333 then estimate it!
	
	volume<short> mask;
	volume<short> maskotl;
	double cost=0;
	int bounds[6]={0,0,0,0,0,0};
	
	float betas=0;
	
	ColumnVector Best(4);
	vector<float> iprof;
	//new cost
	//int M=Vmodel1.at(0)->getNumberOfSubjects();
	
	//assume a single shape per model
	for (int model=0;model<numModels;model++){
		//work sonly if all are tied together
		//these avriables arenned when calculated the ocst function needed for shape priror and for intesndity evaluation
		int M=Vmodel1.at(model)->getNumberOfSubjects();
		double alp=M-1.0/M;
		
		int cumnum=0;
		for (int q=0;q<Vmodel1.at(model)->getNumberOfShapes();q++){
			cumnum+=Vmodel1.at(model)->getNumberOfPoints(q);
		}
		float k2=3.0*static_cast<float>(cumnum);
		
		double binner=0;
		double gammav=alp/(alp-2);
		//if (true){
		if ((!grad)|(gradModel==model)){
			vector<float> dif;
			Mesh m = Vmodel1.at(model)->getDeformedMesh(Vvars.at(model),0,static_cast<int>(Vvars.at(model).size()));
			
			if ((Vpremean.at(model)==-333)&&(GmanMean==-777)&&((!Vglobfoundmode.at(model))|(globReest))){
				betas=0;
				mask=make_mask_from_meshInOut(image,m,Vmodel1.at(model)->getLabel(0),bounds);
				vector<float> vintens;
				Vmodel1.at(model)->intensityHist(&image,&mask,&m,Vmodel1.at(model)->getLabel(0),&vintens);
				float modeval=mode(vintens);
				Vglobfoundmode.at(model)=true;
				VglobMean.at(model)=modeval;
				cout<<"hmmmm find mode "<<modeval<<endl;
			}else if (manMean.value()!=-777){
				VglobMean.at(model)=GmanMean;
				//cout<<"set Mean to: "<<manMean.value()<<" "<<GmanMean<<endl;
			}
			iprof=Vmodel1.at(model)->getDeformedIprof(Vvars.at(model),0,Vvars.at(model).size());
			
			int ipp=Vmodel1.at(model)->getIPP(0);
			
			int count=0;
			
			for (vector<Mpoint*>::iterator k = m._points.begin(); k!=m._points.end(); k++ )
			{
				Pt vertNew;
				Pt vert = (*k)->get_coord();
				Vec normal;
				normal = (*k)->local_normal();
				for (int j=0;j<ipp;j++){
					int sc=j-(ipp-1)/2;
					vertNew=vert + normal*0.5*sc;
					dif.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-VglobMean.at(model) -(iprof.at(count*ipp+j)));
					
				}
				count++;
			}
			//end of shape interation, all dif in dif
			//note this is only good if all arjoint
			
			///%Calculate conditional I | s
			
			vector< vector<float> > prec=Vmodel1.at(model)->getShape(0)->getICondPrec();
			vector<float> eigs=Vmodel1.at(model)->getShape(0)->getICondEigs();
			double probIcond=0;
			for (unsigned int col=0;col<prec.size();col++){
				float multemp=0;
				for (unsigned int row=0;row<dif.size();row++){
					multemp+=dif.at(row)*prec.at(col).at(row);
					
				}
				multemp*=multemp;
				probIcond+=multemp/eigs.at(col);
			}
			//the multiplication by n-1 or n is left out becomes constant in log cost
			//probI*=M;//this multiplication is performed later
			//calculates inner product of difference between observed intensity and mean 
			float sdif=0;
			for (unsigned int row=0;row<dif.size();row++){
				sdif+=dif.at(row)*dif.at(row);
			}
			
			//sdif*=M/(model1->getShape(0)->getErrs().at(1)*2);
			//the M multiplication is performed later
			sdif*=1/(Vmodel1.at(model)->getShape(0)->getErrs().at(1)*2);
			probIcond+=sdif;
			//	cout<<"hmmmm hmmm P(I)"<<endl;
			//////%%%%%%%%%%%%%%%%%%%%%%%%Claculate p(I)
			
			
			prec=Vmodel1.at(model)->getShape(0)->getBCondPrec();
			eigs=Vmodel1.at(model)->getShape(0)->getBCondEigs();
			//float probI=0;
			double probI=0;
			//cout<<"caLC1"<<endl;
			for (unsigned int col=0;col<prec.size();col++){
				
				float multemp=0;
				for (unsigned int row=0;row<dif.size();row++){
					multemp+=dif.at(row)*prec.at(col).at(row);
					
				}
				multemp*=multemp;
				probI+=multemp/eigs.at(col);
			}
			//cout<<"caLC2"<<endl;
			//the multiplication by n-1 or n is left out becomes constant in log cost
			//probI*=M;
			//calculates inner product of difference between observed intensity and mean 
			sdif=0;
			for (unsigned int row=0;row<dif.size();row++){
				sdif+=dif.at(row)*dif.at(row);
			}
			//work sonly if all are tied together
			//multiplication by (M-1) is taken care of below		
			sdif*=1/(Vmodel1.at(model)->getShape(0)->getErrs().at(1)*2);
			probI+=sdif;
			
			/////%%%%%%%%%%%%%%%%%%%%%%%%% Now claculate costfunction
			
			
			
			float k1=static_cast<float>(dif.size());
			
			for (unsigned int i=0;i<Vvars.at(model).size();i++){
				binner+=Vvars.at(model).at(i)*Vvars.at(model).at(i);
			}
			
			//much greater than 1 simplification
			//cost+= (alp+k2)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(probI);//-(12/2-1)*log(avgres) +avgres/2;
			cost+= -(k1)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(1+(M-1)*probIcond/(alp+binner*(gammav)));//-(12/2-1)*log(avgres) +avgres/2;
																														  //subtract p(i)
				cost-=(alp+k1)/2*log(1+(M-1)*probI/alp);
				//cout<<"modelcost "<<model<<" "<<cost<<endl;
		}
		//			cout<<"cost shape prior"<<endl;
		//this is the shaoe prior ...chooeses between no condition, one conditional, or 2 conditionals
		//do fo multiple structures
	//	cout<<"end dualNpreSh cost "<<cost<<endl;

		if (model< static_cast<int>(Vmodel1.size()-1) ){
			cout<<"cond shape"<<endl;
			//assume they are all made form the same number of subjects
			scaleBInvMatrix(Vvars.at(model+1), &vmBx1inv.at(model), M, vKx2.at(model));
			
			ColumnVector Bx2(Vmodel1.at(model)->getNumberOfModes());
			//load bx2 vars that were read
	//		cout<<"vvars"<<endl;
			for (unsigned int i=0;i<Vvars.at(model).size();i++){
				Bx2.element(i)=Vvars.at(model+1).at(i);
	//							cout<<Bx2.element(i)<<" ";

			}
	//		cout<<endl;
			mBx2map=vmBx2.at(model)*Bx2;
			
			
			//only implemented for 2 structures.
			ColumnVector bx1temp(Vvars.at(model).size());
			
			for (unsigned int i=0;i<Vvars.at(model).size();i++){
				bx1temp.element(i)=Vvars.at(model).at(i)-mBx2map.element(i);
	//			cout<<bx1temp.element(i)<<" ";
			}
	//		cout<<endl;
	//		for (unsigned int i=0;i<Vvars.at(model).size();i++){
	//			
	//			cout<<mBx2map.element(i)<<" ";
	//		}
	//		cout<<endl;
			ColumnVector Bcx1=vmBx1inv.at(model)*bx1temp;
			Matrix Bcinner=Bcx1.t()*Bcx1;
		//	cout<<log(1+Bcinner.element(0,0))<<" "<<Bcinner.element(0,0)<<endl;
			cost+=(alp+k2+vKx2.at(model))/2*log(1+Bcinner.element(0,0)/(alp+vKx2.at(model)));
			//and in posteriro bit
			
		}else{
			cout<<"single shape"<<endl;
			//p(x) prior--shape prior
			cost+=(alp+k2)/2*log(1+binner*gammav/alp);
		}
		//This is the denominator in the posteriro
				cout<<"end model dualN cost "<<cost<<endl;

	}
		cout<<"end dualN cost "<<cost<<endl;
	return cost;
	
}
float costfuncAppN2(vector< shapeModel* > Vmodel1, vector< vector<float> > Vvars, vector<float> Vpremean, vector<int> vKx2,vector<  vector<Matrix> > vmBx2,vector<Matrix>  vmBx1inv, vector< vector<int> > vPredLabels, bool grad, int gradModel){
	//only implemented for the shape model
	///asssumes histmode
	//assume 1 structure per model
	//cout<<"enter costfunc appN2"<<endl;
	//if premean equals -333 then estimate it!
//	for (int i =0; i<Vvars.size();i++){
//		for (int j=0; j<Vvars.at(i).size();j++){
//			cout<<Vvars.at(i).at(j)<<" ";
//		}
//		cout<<endl;
//	}
	volume<short> mask;
	volume<short> maskotl;
	double cost=0;
	int bounds[6]={0,0,0,0,0,0};
	
	float betas=0;
	
	ColumnVector Best(4);
	vector<float> iprof;
	//new cost
	//int M=Vmodel1.at(0)->getNumberOfSubjects();
	
	int numModels = static_cast<int>(Vmodel1.size());
	//assume a single shape per model
	for (int model=0;model<numModels;model++){
		//cout<<"MODEL "<<model<<endl;
		//work sonly if all are tied together
		//these avriables arenned when calculated the ocst function needed for shape priror and for intesndity evaluation
		int M=Vmodel1.at(model)->getNumberOfSubjects();
		double alp=M-1.0/M;
		
		int cumnum=0;
		for (int q=0;q<Vmodel1.at(model)->getNumberOfShapes();q++){
			cumnum+=Vmodel1.at(model)->getNumberOfPoints(q);
		}
		float k2=3.0*static_cast<float>(cumnum);
		
		double binner=0;
		double gammav=alp/(alp-2);
		//if (true){
		if ((!grad)|(gradModel==model)){
			vector<float> dif;
			Mesh m = Vmodel1.at(model)->getDeformedMesh(Vvars.at(model),0,static_cast<int>(Vvars.at(model).size()));
			
			if ((Vpremean.at(model)==-333)&&(GmanMean==-777)&&((!Vglobfoundmode.at(model))|(globReest))){
				betas=0;
				mask=make_mask_from_meshInOut(image,m,Vmodel1.at(model)->getLabel(0),bounds);
				vector<float> vintens;
				Vmodel1.at(model)->intensityHist(&image,&mask,&m,Vmodel1.at(model)->getLabel(0),&vintens);
				float modeval=mode(vintens);
				Vglobfoundmode.at(model)=true;
				VglobMean.at(model)=modeval;
				cout<<"hmmmm find mode "<<modeval<<endl;
			}else if (manMean.value()!=-777){
				VglobMean.at(model)=GmanMean;
				////cout<<"set Mean to: "<<manMean.value()<<" "<<GmanMean<<endl;
			}
			iprof=Vmodel1.at(model)->getDeformedIprof(Vvars.at(model),0,Vvars.at(model).size());
			//	cout<<"joitn normalization mean "<<model<<" "<<VglobMean.at(model)<<endl;
			int ipp=Vmodel1.at(model)->getIPP(0);
			
			int count=0;
			
			for (vector<Mpoint*>::iterator k = m._points.begin(); k!=m._points.end(); k++ )
			{
				Pt vertNew;
				Pt vert = (*k)->get_coord();
				Vec normal;
				normal = (*k)->local_normal();
				for (int j=0;j<ipp;j++){
					int sc=j-(ipp-1)/2;
					vertNew=vert + normal*0.5*sc;
					dif.push_back(image.interpolate(vertNew.X/xdim,vertNew.Y/ydim,vertNew.Z/zdim)-VglobMean.at(model) -(iprof.at(count*ipp+j)));
					
				}
				count++;
			}
			//end of shape interation, all dif in dif
			//note this is only good if all arjoint
			
			///%Calculate conditional I | s
			
			vector< vector<float> > prec=Vmodel1.at(model)->getShape(0)->getICondPrec();
			vector<float> eigs=Vmodel1.at(model)->getShape(0)->getICondEigs();
			double probIcond=0;
			for (unsigned int col=0;col<prec.size();col++){
				float multemp=0;
				for (unsigned int row=0;row<dif.size();row++){
					multemp+=dif.at(row)*prec.at(col).at(row);
					
				}
				multemp*=multemp;
				probIcond+=multemp/eigs.at(col);
			}
			//the multiplication by n-1 or n is left out becomes constant in log cost
			//probI*=M;//this multiplication is performed later
			//calculates inner product of difference between observed intensity and mean 
			float sdif=0;
			for (unsigned int row=0;row<dif.size();row++){
				sdif+=dif.at(row)*dif.at(row);
			}
			
			//sdif*=M/(model1->getShape(0)->getErrs().at(1)*2);
			//the M multiplication is performed later
			sdif*=1/(Vmodel1.at(model)->getShape(0)->getErrs().at(1)*2);
			probIcond+=sdif;
			//	//cout<<"hmmmm hmmm P(I)"<<endl;
			//////%%%%%%%%%%%%%%%%%%%%%%%%Claculate p(I)
			
			
			prec=Vmodel1.at(model)->getShape(0)->getBCondPrec();
			eigs=Vmodel1.at(model)->getShape(0)->getBCondEigs();
			//float probI=0;
			double probI=0;
			////cout<<"caLC1"<<endl;
			for (unsigned int col=0;col<prec.size();col++){
				
				float multemp=0;
				for (unsigned int row=0;row<dif.size();row++){
					multemp+=dif.at(row)*prec.at(col).at(row);
					
				}
				multemp*=multemp;
				probI+=multemp/eigs.at(col);
			}
			////cout<<"caLC2"<<endl;
			//the multiplication by n-1 or n is left out becomes constant in log cost
			//probI*=M;
			//calculates inner product of difference between observed intensity and mean 
			sdif=0;
			for (unsigned int row=0;row<dif.size();row++){
				sdif+=dif.at(row)*dif.at(row);
			}
			//work sonly if all are tied together
			//multiplication by (M-1) is taken care of below		
			sdif*=1/(Vmodel1.at(model)->getShape(0)->getErrs().at(1)*2);
			probI+=sdif;
			
			/////%%%%%%%%%%%%%%%%%%%%%%%%% Now claculate costfunction
			
			
			
			float k1=static_cast<float>(dif.size());
			
			for (unsigned int i=0;i<Vvars.at(model).size();i++){
				binner+=Vvars.at(model).at(i)*Vvars.at(model).at(i);
			}
			
			//much greater than 1 simplification
			//cost+= (alp+k2)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(probI);//-(12/2-1)*log(avgres) +avgres/2;
			cost+= -(k1)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(1+(M-1)*probIcond/(alp+binner*(gammav)));//-(12/2-1)*log(avgres) +avgres/2;
																														  //subtract p(i)
				cost-=(alp+k1)/2*log(1+(M-1)*probI/alp);
				////cout<<"modelcost "<<model<<" "<<cost<<endl;
		}
//		}
		//			//cout<<"cost shape prior"<<endl;
		//this is the shaoe prior ...chooeses between no condition, one conditional, or 2 conditionals
		//do fo multiple structures
		//cout<<"Got to the shape prior"<<endl;
		//if (model<Vmodel1.size()-1){
//		if ((!grad)){//this condition ignores joint shape when calculatign gradient
			
			if (model<(numModels-1)){//assumes last model has no conditional
			//	cout<<"Model Number "<<model<<endl;
									 //		for (int model2=0;model2<numModels;model2++){
				if (vPredLabels.at(model).size()>0){
					int bmapcount=0;
					vector< vector<float> > VvarsPred;
			//		cout<<"There are "<<vPredLabels.at(model).size()<<" dependancies "<<Vmodel1.size()<<endl;
					
					while (bmapcount < static_cast<int>(vPredLabels.at(model).size()) ){
						//find indices of desiored labels
			//			cout<<"bmapcoutnt "<<bmapcount<<endl;
						for (unsigned int ibmap=model+1;ibmap<Vmodel1.size();ibmap++){
							
						//	cout<<vPredLabels.at(model).at(bmapcount)<<" "<<Vmodel1.at(ibmap)->getLabel(0)<<endl;
							if (Vmodel1.at(ibmap)->getLabel(0)==vPredLabels.at(model).at(bmapcount)){
								VvarsPred.push_back(Vvars.at(ibmap));
								bmapcount++;
			//					cout<<"FoundN cond struct "<<ibmap<<" "<<Vmodel1.at(ibmap)->getLabel(0)<<endl;
								if (bmapcount==static_cast<int>(vPredLabels.at(model).size())){
									//break if found all dependencies
									break;
								}
							}
							
						}
						
					}
			//		cout<<"Left varspred extracttion"<<endl;
					//assume they are all made form the same number of subjects
					float scale=scaleMultiBInvMatrix(VvarsPred, &vmBx1inv.at(model), M, vKx2.at(model));
					//the implementation fo the bmaps are diffenrewnt form befpore
					//only implemented for 2 structures.
					ColumnVector bx1temp(Vvars.at(model).size());
					//bx1temp=Vvars.at(model);
					//this subtracts the bappede predicted bvars
					for (unsigned int i=0;i<Vvars.at(model).size();i++){
						bx1temp.element(i)=Vvars.at(model).at(i);
					}
					for (unsigned int ibmap=0;ibmap<vPredLabels.at(model).size();ibmap++){
						ColumnVector Bx2(Vmodel1.at(model)->getNumberOfModes());
						//load bx2 vars that were read
						for (unsigned int i=0;i<Vvars.at(model).size();i++){
							Bx2.element(i)=VvarsPred.at(ibmap).at(i);
						}
						mBx2map=vmBx2.at(model).at(ibmap)*Bx2;
						
						for (unsigned int i=0;i<Vvars.at(model).size();i++){
							bx1temp.element(i)-=mBx2map.element(i);
						}
					}	
					ColumnVector Bcx1=vmBx1inv.at(model)*bx1temp;
					Matrix Bcinner=Bcx1.t()*Bcx1;
					//		cout<<"cost pre sh "<<model<<" "<<cost<<" "<<Bcinner.element(0,0)<<" "<<scale<<" "<< k2/2.0*log(scale)<<endl;
					cost+=(alp+k2+vKx2.at(model))/2*log(1+Bcinner.element(0,0)/(alp+vKx2.at(model)));///(alp+vKx2.at(model))
																									 // + k2/2.0*log(scale)
																									 //		cout<<"end sh cost "<<cost<<endl;
																									 //and in posteriro bit
						
			//			cout<<"end"<<endl;
				}else{
					
					//p(x) prior--shape prior
					cost+=(alp+k2)/2*log(1+binner*gammav/alp);
				}
				
			}else{
				//p(x) prior--shape prior
				cost+=(alp+k2)/2*log(1+binner*gammav/alp);
			}
			
//			}else{
//			
//			cost+=(alp+k2)/2*log(1+binner*gammav/alp);
//			}
		//This is the denominator in the posteriro
		
}
	//	//cout<<"end dual cost "<<cost<<endl;
			//cout<<"exit costfunc appN2 "<<cost<<endl;

	return cost;

}
void negGradient(vector<float> *grad, vector<float> p, shapeModel* model1, int ex, int varbeg, int varsize, int& maxind, vector<float> fvals, vector<bool> select,int costtype){
	//costtype=0 gradient ASM
	//costtype=1 AAM
	float searchDim=searchRes;
	vector<float> gradtmp;
	float sumsq=0;
	float costinit=0;
	if (costtype==0){
		costinit=costfunc(model1, p,fvals, ex,0);
	}else if (costtype==1){
		costinit=costfuncApp(model1, p,fvals, ex,0,-333);
	}else if (costtype==2){
		costinit=costfuncOverlap(model1, p,fvals, ex,0);
	}
	float opGrad=0;
	maxind=0;
	float grmaxtmp=0;
	for (int i=0; i<static_cast<int>(p.size()); i++){
		if (select.at(i)){
			
			gradtmp=p;	
			gradtmp.at(i)=gradtmp.at(i)+searchDim;
			//order is reverse because negative gradient
			if (costtype==0){
				grad->at(i)=((costinit-costfunc(model1, gradtmp,fvals, ex,0))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)));
			}else if (costtype==1){
			//if (globReest){
				if (false){
					grad->at(i)=((costinit-costfuncApp(model1, gradtmp,fvals, ex,0,globMean))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)))
					
				}else{
					grad->at(i)=((costinit-costfuncApp(model1, gradtmp,fvals, ex,0,globMean))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)))
				}
					}else if (costtype==2){
						grad->at(i)=((costinit-costfuncOverlap(model1, gradtmp,fvals, ex,0))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)))
	}
			
			if (abs(gradtmp.at(i))>stdTrunc.value()){
				//ignores gradient if goes beyond truncation
				grad->at(i)=grad->at(i)*1e-11;
			}
			//make sure gradient not positive in both directions
			gradtmp.at(i)=gradtmp.at(i)-2*searchDim;
			
			if (costtype==0){
				opGrad=((costinit-costfunc(model1, gradtmp, fvals, ex,0))/(-searchDim*sqrt(model1->getEigenValue(i))));
			}else if (costtype==1){
				//if (globReest){
				if (false){
							opGrad=((costinit-costfuncApp(model1, gradtmp, fvals, ex,0,-333))/(-searchDim*sqrt(model1->getEigenValue(i))));

				}else{
				opGrad=((costinit-costfuncApp(model1, gradtmp, fvals, ex,0,globMean))/(-searchDim*sqrt(model1->getEigenValue(i))));
				}
			}else if (costtype==2){
				opGrad=((costinit-costfuncOverlap(model1, gradtmp, fvals, ex,0))/(-searchDim*sqrt(model1->getEigenValue(i))));
	}
			
			if ((opGrad>0)&&(grad->at(i)<0)){
				grad->at(i)=0;
				
			}else{
				grad->at(i)=(opGrad+grad->at(i))/2;
				
			}
			
			if (abs(grad->at(i))>grmaxtmp){
				grmaxtmp=abs(grad->at(i));
				maxind=i;
			}
			sumsq+=grad->at(i)*grad->at(i);
	}else{
		grad->at(i)=0;
    }
  }
	if (sumsq==0){
		gradzero=true;
		//		cout<<"gradzero"<<endl;
		
	}
	for (unsigned int i=0; i<p.size(); i++){
		grad->at(i)/=sqrt(sumsq);	
		//instead normalize such as max is 1
		//grad->at(i)/=grmaxtmp;
		
		
	}
}
void negGradientAff(vector<float> *grad, vector<float> p, shapeModel* model1, int ex, int varbeg, int varsize, int& maxind, vector<float> fvals, vector<bool> select,int costtype){
	//costtype=0 gradient ASM
	//costtype=1 AAM
	cout<<"enter neg grad"<<endl;
	float searchDim=searchRes;
	vector<float> gradtmp;
	float sumsq=0;
	float costinit=0;
	if (costtype==0){
		costinit=costfunc(model1, p,fvals, ex,0);
	}else if (costtype==1){
		costinit=costfuncAppAff(model1, p,fvals, ex,0,-333);
		//cout<<"cost init "<<costinit<<endl;
	}else if (costtype==2){
		costinit=costfuncOverlap(model1, p,fvals, ex,0);
	}
	float opGrad=0;
	maxind=0;
	float grmaxtmp=0;
	for (int i=0; i<static_cast<int>(p.size()); i++){
		if (select.at(i)){
			
			gradtmp=p;	
			gradtmp.at(i)=gradtmp.at(i)+searchDim;
			//order is reverse because negative gradient
			if (costtype==0){
				grad->at(i)=((costinit-costfunc(model1, gradtmp,fvals, ex,0))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)));
			}else if (costtype==1){
				grad->at(i)=((costinit-costfuncAppAff(model1, gradtmp,fvals, ex,0,globMean))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)))
			//			cout<<"cost initA "<<grad->at(i)<<" "<<gradtmp.at(i)<<endl;

			}else if (costtype==2){
				grad->at(i)=((costinit-costfuncOverlap(model1, gradtmp,fvals, ex,0))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)))
	}
			
			if (abs(gradtmp.at(i))>stdTrunc.value()){
				//ignores gradient if goes beyond truncation
				grad->at(i)=grad->at(i)*1e-11;
			}
			//make sure gradient not positive in both directions
			gradtmp.at(i)=gradtmp.at(i)-2*searchDim;
			
			if (costtype==0){
				opGrad=((costinit-costfunc(model1, gradtmp, fvals, ex,0))/(-searchDim*sqrt(model1->getEigenValue(i))));
			}else if (costtype==1){
				opGrad=((costinit-costfuncAppAff(model1, gradtmp, fvals, ex,0,globMean))/(-searchDim*sqrt(model1->getEigenValue(i))));
				//	cout<<"cost initB "<<opGrad<<" "<<gradtmp.at(i)<<endl;
			}else if (costtype==2){
				opGrad=((costinit-costfuncOverlap(model1, gradtmp, fvals, ex,0))/(-searchDim*sqrt(model1->getEigenValue(i))));
	}
			//cout<<"grad com "<<opGrad<<" "<<grad->at(i)<<endl;
			if ((opGrad>0)&&(grad->at(i)<0)){
				grad->at(i)=0;
				
			}else{
				grad->at(i)=(opGrad+grad->at(i))/2;
				
			}//else if (opGrad<=0){
				//grad->at(i)=opGrad;
			
			//}//
			
			if (abs(grad->at(i))>grmaxtmp){
				grmaxtmp=abs(grad->at(i));
				maxind=i;
			}
			sumsq+=grad->at(i)*grad->at(i);
	}else{
		grad->at(i)=0;
    }
  }
	if (sumsq==0){
		gradzero=true;
		//		cout<<"gradzero"<<endl;
		
	}
	for (unsigned int i=0; i<p.size(); i++){
		grad->at(i)/=sqrt(sumsq);	
	//	cout<<"grad "<<i<<" "<<grad->at(i)<<endl;
		//instead normalize such as max is 1
		//grad->at(i)/=grmaxtmp;
		
		
	}
		cout<<"exit neg grad"<<endl;

}

void negGradientDual(vector <vector<float> >* Vgrad, vector< vector<float> > Vp, vector< shapeModel* > Vmodel1, int varbeg, int varsize, vector<int>* Vmaxind, vector<bool> select){
	//only implemented for appearance model
	//only implemented for same modes of each structure change select if different is desired
	cout<<"Enter NEGGRADIENTDUAL "<<static_cast<int>(Vmodel1.size())<<endl;
	float searchDim=searchRes;

	float sumsq=0;
	
	vector< vector<float> > Vgradtmp;

		vector<float> meantmp;
		for (int model=0; model< static_cast<int>(Vmodel1.size()); model++){
	
			
			meantmp.push_back(-333);
		}
	
		
		
		for (int model=0; model< static_cast<int>(Vmodel1.size()); model++){
			
			float costinit=costfuncAppDual(Vmodel1, Vp ,meantmp);
			float opGrad=0;
			Vmaxind->at(model)=0;
			float grmaxtmp=0;
			cout<<"negGrad "<<Vp.at(model).size()<<" "<<model<<endl;
			for (int i=0; i<static_cast<int>(Vp.at(model).size()); i++){
			//cout<<"i "<<i<<endl;
				if (select.at(i)){
				//	cout<<"i2 "<<i<<endl;
					Vgradtmp=Vp;	
					Vgradtmp.at(model).at(i)=Vgradtmp.at(model).at(i)+searchDim;
					//order is reverse because negative gradient
					
					Vgrad->at(model).at(i)=((costinit-costfuncAppDual(Vmodel1, Vgradtmp,VglobMean))/(searchDim*sqrt(Vmodel1.at(model)->getEigenValue(i))));//(pprev.at(i)))
						
						
						if (abs(Vgradtmp.at(model).at(i))>stdTrunc.value()){
							//ignores gradient if goes beyond truncation
							Vgrad->at(model).at(i)=Vgrad->at(model).at(i)*1e-11;
						}
						//make sure gradient not positive in both directions
						Vgradtmp.at(model).at(i)=Vgradtmp.at(model).at(i)-2*searchDim;
						
						opGrad=((costinit-costfuncAppDual(Vmodel1, Vgradtmp,VglobMean))/(-searchDim*sqrt(Vmodel1.at(model)->getEigenValue(i))));
						
						if ((opGrad>0)&&(Vgrad->at(model).at(i)<0)){
							Vgrad->at(model).at(i)=0;
							
						}else{
							Vgrad->at(model).at(i)=(opGrad+Vgrad->at(model).at(i))/2;
							
						}
						
						if (abs(Vgrad->at(model).at(i))>grmaxtmp){
							grmaxtmp=abs(Vgrad->at(model).at(i));
							Vmaxind->at(model)=i;
						}
						sumsq+=Vgrad->at(model).at(i)*Vgrad->at(model).at(i);
	}else{
		Vgrad->at(model).at(i)=0;
    }
  }
			
  }
		
	if (sumsq==0){
		gradzero=true;
		//		cout<<"gradzero"<<endl;
		
	}
		for (int model=0; model< static_cast<int>(Vmodel1.size()); model++){
			for (unsigned int i=0; i<Vp.at(model).size(); i++){
				//this normalizes by total sum do squares, should I do individual sum of squares?
				Vgrad->at(model).at(i)/=sqrt(sumsq);	
				//instead normalize such as max is 1
				//grad->at(i)/=grmaxtmp;
				
				
			}
		}
		cout<<"END NEGGRADIENTDUAL "<<static_cast<int>(Vmodel1.size())<<endl;
}

void negGradientN(vector <vector<float> >* Vgrad, vector< vector<float> > Vp, vector< shapeModel* > Vmodel1, int varbeg, int varsize, vector<int>* Vmaxind, vector<bool> select, vector<int> vKx2, vector<Matrix> vmBx2,vector<Matrix> vmBx1inv){
	//only implemented for appearance model
	//only implemented for same modes of each structure change select if different is desired
	cout<<"Enter NEGGRADIENTDUAL "<<static_cast<int>(Vmodel1.size())<<endl;
	float searchDim=searchRes;

	float sumsq=0;
	
	vector< vector<float> > Vgradtmp;

		vector<float> meantmp;
		for (int model=0; model< static_cast<int>(Vmodel1.size()); model++){
	
			
			meantmp.push_back(-333);
		}
	
		
		
		for (int model=0; model< static_cast<int>(Vmodel1.size()); model++){
			//float costfuncAppN(vector< shapeModel* > Vmodel1, vector< vector<float> > Vvars, vector<float> Vpremean, vector<int> vKx2, vector<Matrix> vmBx2,vector<Matrix> vmBx1inv, bool grad, int gradModel){

			float costinit=costfuncAppN(Vmodel1, Vp ,meantmp, vKx2, vmBx2, vmBx1inv, true, model);
			float opGrad=0;
			Vmaxind->at(model)=0;
			float grmaxtmp=0;
			cout<<"negGrad "<<Vp.at(model).size()<<" "<<model<<endl;
			for (int i=0; i<static_cast<int>(Vp.at(model).size()); i++){
			//cout<<"i "<<i<<endl;
				if (select.at(i)){
				//	cout<<"i2 "<<i<<endl;
					Vgradtmp=Vp;	
					Vgradtmp.at(model).at(i)=Vgradtmp.at(model).at(i)+searchDim;
					//order is reverse because negative gradient
					
					Vgrad->at(model).at(i)=((costinit-costfuncAppN(Vmodel1, Vgradtmp,VglobMean, vKx2, vmBx2, vmBx1inv, true, model))/(searchDim*sqrt(Vmodel1.at(model)->getEigenValue(i))));//(pprev.at(i)))
						
						
						if (abs(Vgradtmp.at(model).at(i))>stdTrunc.value()){
							//ignores gradient if goes beyond truncation
							Vgrad->at(model).at(i)=Vgrad->at(model).at(i)*1e-11;
						}
						//make sure gradient not positive in both directions
						Vgradtmp.at(model).at(i)=Vgradtmp.at(model).at(i)-2*searchDim;
						
						opGrad=((costinit-costfuncAppN(Vmodel1, Vgradtmp,VglobMean, vKx2, vmBx2, vmBx1inv, true, model))/(-searchDim*sqrt(Vmodel1.at(model)->getEigenValue(i))));
						
						if ((opGrad>0)&&(Vgrad->at(model).at(i)<0)){
							Vgrad->at(model).at(i)=0;
							
						}else{
							Vgrad->at(model).at(i)=(opGrad+Vgrad->at(model).at(i))/2;
							
						}
						
						if (abs(Vgrad->at(model).at(i))>grmaxtmp){
							grmaxtmp=abs(Vgrad->at(model).at(i));
							Vmaxind->at(model)=i;
						}
						sumsq+=Vgrad->at(model).at(i)*Vgrad->at(model).at(i);
	}else{
		Vgrad->at(model).at(i)=0;
    }
  }
			
  }
		
	if (sumsq==0){
		gradzero=true;
		//		cout<<"gradzero"<<endl;
		
	}
		for (int model=0; model< static_cast<int>(Vmodel1.size()); model++){
			for (unsigned int i=0; i<Vp.at(model).size(); i++){
				//this normalizes by total sum do squares, should I do individual sum of squares?
				cout<<"grad N "<<Vgrad->at(model).at(i)<<" "<<sumsq<<endl;
				Vgrad->at(model).at(i)/=sqrt(sumsq);	
				//instead normalize such as max is 1
				//grad->at(i)/=grmaxtmp;
				
				
			}
		}
		cout<<"END NEGGRADIENTDUAL "<<static_cast<int>(Vmodel1.size())<<endl;
}
void negGradientN2(vector <vector<float> >* Vgrad, vector< vector<float> > Vp, vector< shapeModel* > Vmodel1, int varbeg, int varsize, vector<int>* Vmaxind, vector<bool> select, vector<int> vKx2, vector< vector<Matrix> > vmBx2,vector<Matrix> vmBx1inv, vector< vector<int> > vPredLabels){
	//only implemented for appearance model
	//only implemented for same modes of each structure change select if different is desired
	cout<<"Enter NEGGRADIENTDUAL "<<static_cast<int>(Vmodel1.size())<<endl;
	float searchDim=searchRes;

	float sumsq=0;
	
	vector< vector<float> > Vgradtmp;

		vector<float> meantmp;
		for (int model=0; model< static_cast<int>(Vmodel1.size()); model++){
	
			
			meantmp.push_back(-333);
		}
	
	int numModes=0;	
		float grmaxtmp=0;
		for (int model=0; model< static_cast<int>(Vmodel1.size()); model++){
			//float costfuncAppN(vector< shapeModel* > Vmodel1, vector< vector<float> > Vvars, vector<float> Vpremean, vector<int> vKx2, vector<Matrix> vmBx2,vector<Matrix> vmBx1inv, bool grad, int gradModel){

			float costinit=costfuncAppN2(Vmodel1, Vp ,meantmp, vKx2, vmBx2, vmBx1inv,vPredLabels, true, model);
			float opGrad=0;
			Vmaxind->at(model)=0;
			
			cout<<"negGrad "<<Vp.at(model).size()<<" "<<model<<endl;
			for (int i=0; i<static_cast<int>(Vp.at(model).size()); i++){
			//cout<<"i "<<i<<endl;
				if (select.at(i)){
				numModes++;
				//	cout<<"i2 "<<i<<endl;
					Vgradtmp=Vp;	
					Vgradtmp.at(model).at(i)=Vgradtmp.at(model).at(i)+searchDim;
					//order is reverse because negative gradient
					
					Vgrad->at(model).at(i)=((costinit-costfuncAppN2(Vmodel1, Vgradtmp,VglobMean, vKx2, vmBx2, vmBx1inv,vPredLabels, true, model))/(searchDim*sqrt(Vmodel1.at(model)->getEigenValue(i))));//(pprev.at(i)))
						
						
						if (abs(Vgradtmp.at(model).at(i))>stdTrunc.value()){
							//ignores gradient if goes beyond truncation
							Vgrad->at(model).at(i)=Vgrad->at(model).at(i)*1e-11;
						}
						//make sure gradient not positive in both directions
						Vgradtmp.at(model).at(i)=Vgradtmp.at(model).at(i)-2*searchDim;
						
						opGrad=((costinit-costfuncAppN2(Vmodel1, Vgradtmp,VglobMean, vKx2, vmBx2, vmBx1inv,vPredLabels, true, model))/(-searchDim*sqrt(Vmodel1.at(model)->getEigenValue(i))));
						
						
				//		Vgrad->at(model).at(i)=(opGrad+Vgrad->at(model).at(i))/2;
						
						if ((opGrad>0)&&(Vgrad->at(model).at(i)<0)){
							Vgrad->at(model).at(i)=0;
							
						}else{
							Vgrad->at(model).at(i)=(opGrad+Vgrad->at(model).at(i))/2;
							//if (abs(opGrad)>abs(Vgrad->at(model).at(i))){
							//	Vgrad->at(model).at(i)=opGrad;
							//}
							//Vgrad->at(model).at(i)=(opGrad+Vgrad->at(model).at(i))/2;

						}
						
						if (abs(Vgrad->at(model).at(i))>grmaxtmp){
							grmaxtmp=abs(Vgrad->at(model).at(i));
							Vmaxind->at(model)=i;
						}
						sumsq+=Vgrad->at(model).at(i)*Vgrad->at(model).at(i);
	}else{
		Vgrad->at(model).at(i)=0;
    }
  }
			
  }
		
	if (sumsq==0){
		gradzero=true;
		//		cout<<"gradzero"<<endl;
		
	}
		for (int model=0; model< static_cast<int>(Vmodel1.size()); model++){
			for (unsigned int i=0; i<Vp.at(model).size(); i++){
				//this normalizes by total sum do squares, should I do individual sum of squares?
							//	cout<<"grad N "<<Vgrad->at(model).at(i)<<" "<<sumsq<<endl;

				Vgrad->at(model).at(i)/=sqrt(sumsq);///numModes;//static_cast<float>(Vmodel1.size());	
				//instead normalize such as max is 1
			//	Vgrad->at(model).at(i)/=grmaxtmp;
				
				
			}
		}
		cout<<"END NEGGRADIENTDUAL "<<static_cast<int>(Vmodel1.size())<<endl;
}

float conjGradient(shapeModel* model1, vector<float>* vars, int varbeg, int varsize, int ex, vector<float> relStd, vector<float> fvals, vector<bool> select, int costtype, float searchR, float searchRmax){
	cout<<"conjGradient"<<endl;
	//costtype defines whether to use AAM or ASM
	//define variable
	ofstream fcostout;
	if (costtype==2){
		model1->volumeBounds(&image,gbounds);
		//	cout<<gbounds[0]<<" "<<gbounds[1]<<" "<<gbounds[2]<<" "<<gbounds[3]<<" "<<gbounds[4]<<" "<<gbounds[5]<<endl;
		fcostout.open("cost_traj.txt");
	}
	
	vector<float>  svec, res, resPrev;//, pointCost;
		double gamma;
		int n=vars->size();
		double dpResPrev=0; 
		vector<float> gradtmp2=res;
		vector<float> svectmp=*vars;
		vector<float> gradtmp=*vars;
		gradzero=false;
		searchRes=searchR;
		int maxind=0;
		
		//to find cost at truncation
		vector<float> varMax=*vars;
		
		
		//	searchDim=0.001;
		resPrev=*vars;
		res=*vars;
		
		
		negGradient(&res,*vars,model1, ex, varbeg, varsize, maxind, fvals,select,costtype);
		if ((gradzero)){
			cout<<"return nor gradient"<<endl;
			return -777;//if graident == 0
		}
		resPrev=res;
		
		svec=svectmp=res;
		vector<float> tmpVar2;
		for (int i=0;i<n;i++){
			tmpVar2.push_back(0);
		}
		
		vector<bool> vmax;
		vector<float> costVMAX;///used to see if upper bound hit and that is the max
			for (int i=0;i<n;i++){
				vmax.push_back(false);
				costVMAX.push_back(10e6);
			}
			float maxStd=stdTrunc.value();
			for (int iter=0;iter<40;iter++){
				
				//find line mminimizatioon and set svec to vector displavment
				
				//start fitting first mode
			  //	cout<<"iter: "<<iter<<endl;
				
				float tmpVar=0;
				float tmpCost=0;
				vector<float> varstmp=*vars;
				int j=0;
				for (int i=0;i<n;i++){
					tmpVar2.at(i)=0;
				}
				for (int i=0;i<n;i++){
					vmax.at(i)=false;
				}
				int zerocount=0; //to count number of succesive zeros->help speed search
				int searchDist=60;
			
				while ((j<30)){
					
					//enable single select
					for (unsigned int s=0; s<vars->size(); s++){
						if ((abs(varstmp.at(s)+ j*searchRes*svec.at(s))>maxStd)|(abs(gradtmp.at(s)-varstmp.at(s)-j*searchRes*svec.at(s))>relStd.at(s))){
							
							//compare against absoulte allowable maximum std and against relative varaition form start point
							
							//save the value that will render the appropriate max value
							if (!vmax.at(s)){//if not already at its max
								
								if (svec.at(s)==0){//if gradient equals zero, just leave
									tmpVar2.at(s)=0;
								}else{
									if ((abs(varstmp.at(s)+ j*searchRes*svec.at(s))>maxStd)){
										//which max did it hit, handle differently
										if ((varstmp.at(s)+ j*searchRes*svec.at(s))>maxStd){
											tmpVar2.at(s)=(maxStd-varstmp.at(s))/svec.at(s);
										}else{
											tmpVar2.at(s)=(-maxStd-varstmp.at(s))/svec.at(s);
										}
									}else{
										if (gradtmp.at(s)-varstmp.at(s)-j*searchRes*svec.at(s)<relStd.at(s)){
											tmpVar2.at(s)=(gradtmp.at(s)+relStd.at(s)-varstmp.at(s))/svec.at(s);
										}else{
											tmpVar2.at(s)=(gradtmp.at(s)-relStd.at(s)-varstmp.at(s))/svec.at(s);
										}
										
									}
								}
							}
							vmax.at(s)=true;
						}
						
						if (!vmax.at(s)){
							vars->at(s)=varstmp.at(s) + j*searchRes*svec.at(s);
						}else{
							vars->at(s)=varstmp.at(s) + tmpVar2.at(s)*svec.at(s);
						}
					}
					float cost=0;
					if (costtype==0){
						cost=costfunc(model1, *vars, fvals, ex,0);
					}else if (costtype==1){
						cost=costfuncApp(model1, *vars, fvals, ex,0,-333);
					}else if (costtype==2){
						cost=costfuncOverlap(model1,*vars,fvals, ex,0);
						
	}	
					if (j==0){	     
						tmpVar=j*searchRes;
						tmpCost=cost;
					}else if (cost<tmpCost){
						tmpVar=j*searchRes;
						tmpCost=cost;
						zerocount=0;//reset count when tmpvar is updated
					}else{
						zerocount++;
					}
					
					
					if (zerocount>3){//causes break in linear search if nothing found in 5
						j=searchDist;
					}
					
					j++;
				}
				//assign values of position and displacment
				float delta=0;
				
				//this enables single select
				for (unsigned int s=0; s<svectmp.size(); s++){ 
					if ((!vmax.at(s))|((vmax.at(s))&&(tmpVar<=tmpVar2.at(s)))){
						vars->at(s)=varstmp.at(s)+tmpVar*svec.at(s);
						svectmp.at(s)=tmpVar*svectmp.at(s);
						delta+=tmpVar*svec.at(s)*tmpVar*svec.at(s);
					}else{
						vars->at(s)=varstmp.at(s)+tmpVar2.at(s)*svec.at(s);
						svectmp.at(s)=tmpVar2.at(s)*svectmp.at(s);
						delta+=tmpVar2.at(s)*svec.at(s)*tmpVar2.at(s)*svec.at(s);
					}
					
					//delta not actually used in algorithm
				}
				if (delta<0.01){
					tmpVar=0;
				}
				//write new parameters in neggradient call
				
				//incraese search resolution 
					cout<<"vars:"<<endl;
					for (unsigned int i=0;i<vars->size();i++){
					cout<<vars->at(i)<<" ";
				}
				cout<<endl;
				
				
				if ((tmpVar==0)){
					//play with this maaybe no chane in res
					searchRes=searchRes/2;
					//	cout<<"addvolume "<<endl;
					//addVolumeTo4D(model1,*vars);
					if (searchRes<searchRmax){
						gradzero=true;
					}
					
				}
				
				if(!gradzero){
					bool breakloop=false;
					
					while(!breakloop){
						negGradient(&svectmp,*vars,model1, ex, varbeg, varsize, maxind, fvals,select,costtype);
						if (gradzero==true){
							searchRes=searchRes/2;
							if (searchRes<searchRmax){
								breakloop=true;
							}else{
								gradzero=false;
							}
						}else{
							breakloop=true;
						}
					}	
					
				}
				if (gradzero){
					//addVolumeTo4D(model1,*vars);
					if (costtype==0){
					costfunc(model1, *vars, fvals, ex,0);
					}else if (costtype==1){
					costfuncApp(model1, *vars, fvals, ex,0,-333);
					}
					//					cout<<"gbeta1 "<<gBeta<<endl;
					
					gradzero=false;
					//costfunc app will set global means
					cout<<"Leaving conjgraident with mode "<<globMean<<endl;
					return globMean;
				}
				
				//calculate new residual dot product
				double dpRes=0;
				//for (int i=0; i<n; i++){
				for (int i=varbeg; i<varbeg+varsize; i++){
					dpResPrev+=res.at(i)*res.at(i);
					dpRes+=svectmp.at(i)*svectmp.at(i);
					dpRes+=(-svectmp.at(i)+res.at(i))*-svectmp.at(i);
				}
				gamma= dpRes/dpResPrev;
				
				
				//update update vector
				
				resPrev=res;
				res=svectmp;
				
				//so that single selects may be used			
				for (unsigned int i=0; i<svec.size(); i++){
					svectmp.at(i)=svec.at(i)=res.at(i)+gamma*svectmp.at(i);
					//This changes conjugate gradient to gradient
				}
				if (costtype==2){
					fcostout<<iter<<" "<<abs(tmpCost)<<endl;
				}
				
				}
			addVolumeTo4D(model1,*vars);
			}
			
float conjGradientAff(shapeModel* model1, vector<float>* vars, int varbeg, int varsize, int ex, vector<float> relStd, vector<float> fvals, vector<bool> select, int costtype, float searchR, float searchRmax){
	cout<<"conjGradient"<<endl;
	//costtype defines whether to use AAM or ASM
	//define variable
	ofstream fcostout;
	if (costtype==2){
		model1->volumeBounds(&image,gbounds);
		//	cout<<gbounds[0]<<" "<<gbounds[1]<<" "<<gbounds[2]<<" "<<gbounds[3]<<" "<<gbounds[4]<<" "<<gbounds[5]<<endl;
		fcostout.open("cost_traj.txt");
	}
	
	vector<float>  svec, res, resPrev;//, pointCost;
		double gamma;
		int n=vars->size();
		double dpResPrev=0; 
		vector<float> gradtmp2=res;
		vector<float> svectmp=*vars;
		vector<float> gradtmp=*vars;
		gradzero=false;
		searchRes=searchR;
		int maxind=0;
		
		//to find cost at truncation
		vector<float> varMax=*vars;
		
		
		//	searchDim=0.001;
		resPrev=*vars;
		res=*vars;
		
		
		negGradientAff(&res,*vars,model1, ex, varbeg, varsize, maxind, fvals,select,costtype);
		if ((gradzero)){
			cout<<"return nor gradient"<<endl;
			return -777;//if graident == 0
		}
		resPrev=res;
		
		svec=svectmp=res;
		vector<float> tmpVar2;
		for (int i=0;i<n;i++){
			tmpVar2.push_back(0);
		}
		
		vector<bool> vmax;
		vector<float> costVMAX;///used to see if upper bound hit and that is the max
			for (int i=0;i<n;i++){
				vmax.push_back(false);
				costVMAX.push_back(10e6);
			}
			float maxStd=stdTrunc.value();
			for (int iter=0;iter<40;iter++){
				
				//find line mminimizatioon and set svec to vector displavment
					//impse hierarchical between affine and shaep
					float tmpVar=0;
					float tmpCost=0;
					float delta=0;
					for (int groups=0; groups<1;groups++){
					

				//start fitting first mode
			  //	cout<<"iter: "<<iter<<endl;
				
				tmpVar=0;
				tmpCost=0;
				vector<float> varstmp=*vars;
				int j=0;
				for (int i=0;i<n;i++){
					tmpVar2.at(i)=0;
				}
				for (int i=0;i<n;i++){
					vmax.at(i)=false;
				}
				int zerocount=0; //to count number of succesive zeros->help speed search
				int searchDist=60;
			
				while ((j<30)){
					
									//enable single select
					for (unsigned int s=0; s<vars->size(); s++){
						if ((abs(varstmp.at(s)+ j*searchRes*svec.at(s))>maxStd)|(abs(gradtmp.at(s)-varstmp.at(s)-j*searchRes*svec.at(s))>relStd.at(s))){
							
							//compare against absoulte allowable maximum std and against relative varaition form start point
							
							//save the value that will render the appropriate max value
							if (!vmax.at(s)){//if not already at its max
								//the or here is whta allows for hierarchy
								cout<<"groups "<<groups<<" "<<s<<endl;
							//	if ((svec.at(s)==0) | ((groups==0)&&(s<7)) | ((groups==1)&&(s>=7))) {//if gradient equals zero, just leave
								if ((svec.at(s)==0)) {//if gradient equals zero, just leave

											cout<<"set grad 0 "<<s<<endl;
									tmpVar2.at(s)=0;
								}else{
									if ((abs(varstmp.at(s)+ j*searchRes*svec.at(s))>maxStd)){
										//which max did it hit, handle differently
										if ((varstmp.at(s)+ j*searchRes*svec.at(s))>maxStd){
											tmpVar2.at(s)=(maxStd-varstmp.at(s))/svec.at(s);
										}else{
											tmpVar2.at(s)=(-maxStd-varstmp.at(s))/svec.at(s);
										}
									}else{
										if (gradtmp.at(s)-varstmp.at(s)-j*searchRes*svec.at(s)<relStd.at(s)){
											tmpVar2.at(s)=(gradtmp.at(s)+relStd.at(s)-varstmp.at(s))/svec.at(s);
										}else{
											tmpVar2.at(s)=(gradtmp.at(s)-relStd.at(s)-varstmp.at(s))/svec.at(s);
										}
										
									}
								}
							}
							vmax.at(s)=true;
						}
			//			if ( ((groups==0)&&(s<7)) | ((groups==1)&&(s>=7))) {
											//		cout<<"groups "<<groups<<" "<<s<<endl;

							if (!vmax.at(s)){
								//	cout<<"groups "<<groups<<" "<<s<<endl;
								
								vars->at(s)=varstmp.at(s) + j*searchRes*svec.at(s);
								
							}else{
								vars->at(s)=varstmp.at(s) + tmpVar2.at(s)*svec.at(s);
							}
			//			}
					}
				
					float cost=0;
					if (costtype==0){
						cost=costfunc(model1, *vars, fvals, ex,0);
					}else if (costtype==1){
						cost=costfuncAppAff(model1, *vars, fvals, ex,0,-333);
					}else if (costtype==2){
						cost=costfuncOverlap(model1,*vars,fvals, ex,0);
						
	}	
					if (j==0){	     
						tmpVar=j*searchRes;
						tmpCost=cost;
					}else if (cost<tmpCost){
						tmpVar=j*searchRes;
						tmpCost=cost;
						zerocount=0;//reset count when tmpvar is updated
					}else{
						zerocount++;
					}
					
					
					if (zerocount>3){//causes break in linear search if nothing found in 5
						j=searchDist;
					}
					
					j++;
				}
				//assign values of position and displacment
				delta=0;
				
				//this enables single select
				for (unsigned int s=0; s<svectmp.size(); s++){ 
			//	 if ( ((groups==0)&&(s<7)) | ((groups==1)&&(s>=7)) ){
				//	cout<<"groups "<<groups<<" "<<s<<endl;
					if ((!vmax.at(s))|((vmax.at(s))&&(tmpVar<=tmpVar2.at(s)))){
						vars->at(s)=varstmp.at(s)+tmpVar*svec.at(s);
						svectmp.at(s)=tmpVar*svectmp.at(s);
						delta+=tmpVar*svec.at(s)*tmpVar*svec.at(s);
					}else{
						vars->at(s)=varstmp.at(s)+tmpVar2.at(s)*svec.at(s);
						svectmp.at(s)=tmpVar2.at(s)*svectmp.at(s);
						delta+=tmpVar2.at(s)*svec.at(s)*tmpVar2.at(s)*svec.at(s);
					}
			//		}
					//delta not actually used in algorithm
				}
				
					
				//incraese search resolution 
						cout<<"vars Group: "<<groups<<endl;
						for (unsigned int i=0;i<vars->size();i++){
							cout<<vars->at(i)<<" ";
						}
					cout<<endl;
				
				}
				if (delta<0.01){
					tmpVar=0;
				}
				//write new parameters in neggradient call
			
				
				if ((tmpVar==0)){
					//play with this maaybe no chane in res
					searchRes=searchRes/2;
					//	cout<<"addvolume "<<endl;
					//addVolumeTo4D(model1,*vars);
					if (searchRes<searchRmax){
						gradzero=true;
					}
					
				}
				
				if(!gradzero){
					bool breakloop=false;
					
					while(!breakloop){
						cout<<"neg loop"<<endl;
						negGradientAff(&svectmp,*vars,model1, ex, varbeg, varsize, maxind, fvals,select,costtype);
						cout<<"end neg"<<endl;
						if (gradzero==true){
							searchRes=searchRes/2;
							if (searchRes<searchRmax){
								breakloop=true;
							}else{
								gradzero=false;
							}
						}else{
							breakloop=true;
						}
					}	
					
				}
				cout<<"out of loop"<<endl;
				if (gradzero){
					//addVolumeTo4D(model1,*vars);
					costfuncAppAff(model1, *vars, fvals, ex,0,-333);
					//					cout<<"gbeta1 "<<gBeta<<endl;
					
					gradzero=false;
					//costfunc app will set global means
					cout<<"Leaving conjgraident with mode "<<globMean<<endl;
					return globMean;
				}
				
				//calculate new residual dot product
				double dpRes=0;
				//for (int i=0; i<n; i++){
				for (int i=varbeg; i<varbeg+varsize; i++){
					dpResPrev+=res.at(i)*res.at(i);
					dpRes+=svectmp.at(i)*svectmp.at(i);
					dpRes+=(-svectmp.at(i)+res.at(i))*-svectmp.at(i);
				}
				gamma= dpRes/dpResPrev;
				
				
				//update update vector
				
				resPrev=res;
				res=svectmp;
				
				//so that single selects may be used			
				for (unsigned int i=0; i<svec.size(); i++){
					svectmp.at(i)=svec.at(i)=res.at(i)+gamma*svectmp.at(i);
					//This changes conjugate gradient to gradient
				}
				if (costtype==2){
					fcostout<<iter<<" "<<abs(tmpCost)<<endl;
				}
				
				}
			//addVolumeTo4D(model1,*vars);
			}
			
void conjGradientDual(vector< shapeModel* > Vmodel1, vector< vector<float> >* Vvars , int varbeg, int varsize,  vector<float> relStd, vector<bool> select, float searchR, float searchRmax){
	//this will load 2 models and fit simultaneously
	//will need to keep track of as many variables as models
	//though some is general use of conditionals is only for 2 models at this point
	cout<<"conjGradient Dual"<<endl;
	//FOR NOW DEFINE NUMBER OF MODELS AS 2
	int numModels=static_cast<int>(Vmodel1.size());
	//costtype defines whether to use AAM or ASM
	//define variable
	

	//one for each model
	vector< vector<float> >  Vsvec, Vres, VresPrev;//, pointCost;


		vector< double > Vgamma;
		
		//everything is stored in vectors such that it can be generalized to multiple structures
		vector< int > Vn;
		vector< double > VdpResPrev; 
		vector< vector<float> > Vgradtmp;
		vector< vector<float> > Vsvectmp;
		//vector< vector<float> > Vgradtmp2;
		gradzero=false;
		searchRes=searchR;
		
		vector<int> Vmaxind;
		vector< vector<float> > VvarMax;
		
		vector< vector<bool> > Vvmax;
		vector< vector<float> > VcostVMAX;
		vector< vector<float> > VtmpVar2;
		
		//initializations
		float maxStd=stdTrunc.value();
		cout<<"hmmm1 "<<numModels<<endl;
		for (int model=0; model<numModels; model++){
			Vgamma.push_back(0);
			Vn.push_back( Vvars->at(model).size() );
			
			VdpResPrev.push_back(0); 
			
		//	vector<float> gradtmp2=Vres.at(model);
		//	Vgradtmp2.push_back(gradtmp2);
			
			vector<float> svectmp=(*Vvars).at(model);
			Vsvectmp.push_back(svectmp);
			
			vector<float> gradtmp=(*Vvars).at(model);
			Vgradtmp.push_back(gradtmp);
			
			//this stuff is used to handle cost at truncation
			Vmaxind.push_back(0);
		
			//to find cost at truncation
			vector<float> varMax=(*Vvars).at(model);
			VvarMax.push_back(varMax);
		
			//	searchDim=0.001;
			vector<float> resPrev=(*Vvars).at(model);
			VresPrev.push_back(resPrev);
			
			vector<float> res=(*Vvars).at(model);
			Vres.push_back(res);
			
			vector<float> svec=(*Vvars).at(model);
			Vsvec.push_back(svec);
		
		}
		//neggradientdual call should not be in the model loop
			//select is left the same for now....search same modes for each structure
			//cout<<"hmmm2 "<<endl;
			negGradientDual(&Vres,*Vvars,Vmodel1, varbeg, varsize, &Vmaxind, select );
				cout<<"hmmm3 "<<endl;
			if ((gradzero)){
				return;//if graident == 0
			}
			for (int model=0; model<numModels; model++){
				//cout<<"hmmm3.1 "<<endl;
				VresPrev.at(model)=Vres.at(model);
				//	cout<<"hmmm3.1 "<<endl;
				//Vsvec.at(model)=Vsvectmp.at(model)=Vres.at(model);
				Vsvec.at(model)=Vres.at(model);
				Vsvectmp.at(model)=Vres.at(model);
		//		cout<<"hmmm3.1 "<<endl;
				
				vector<float> tmpVar2;
				for (int i=0;i<Vn.at(model);i++){
					tmpVar2.push_back(0);
				}
				VtmpVar2.push_back(tmpVar2);
		//		cout<<"hmmm3.2 "<<endl;
				
				vector<bool> vmax;
				vector<float> costVMAX;///used to see if upper bound hit and that is the max
					
					for (int i=0;i<Vn.at(model);i++){
						vmax.push_back(false);
						costVMAX.push_back(10e6);
					}
					Vvmax.push_back(vmax);
					VcostVMAX.push_back(costVMAX);
			} ///end of initialization for each model
			//STILL NEED TO DEAL WITH GRADZERO
			//	cout<<"hmmm4 "<<endl;
			for (int iter=0;iter<40;iter++){
				
				//find line mminimizatioon and set svec to vector displavment
				
				//start fitting first mode
			  //	cout<<"iter: "<<iter<<endl;
				vector<float> VtmpVar;
				vector<float> VtmpCost;
				vector< vector<float> > Vvarstmp;
				for (int model=0; model<numModels; model++){
					VtmpVar.push_back(0);
					VtmpCost.push_back(0);
					
					vector<float> varstmp=(*Vvars).at(model);
					Vvarstmp.push_back(varstmp);
					for (int i=0;i<Vn.at(model);i++){
						VtmpVar2.at(model).at(i)=0;
					}
					for (int i=0;i<Vn.at(model);i++){
						Vvmax.at(model).at(i)=false;
					}
				
				}
				
			
				int zerocount=0; //to count number of succesive zeros->help speed search
				int searchDist=60;
				int j=0;
				while ((j<30)){
					//used for input to costfuncappdual
					vector<float> vmeantmp;
					
					//enable single select
					for (int model=0; model<numModels; model++){
						for (unsigned int s=0; s<Vvars->at(model).size(); s++){
							if ((abs(Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd)|(abs(Vgradtmp.at(model).at(s)-Vvarstmp.at(model).at(s)-j*searchRes*Vsvec.at(model).at(s))>relStd.at(s))){
								
								//compare against absoulte allowable maximum std and against relative varaition form start point
								
								//save the value that will render the appropriate max value
								if (!Vvmax.at(model).at(s)){//if not already at its max
									
									if (Vsvec.at(model).at(s)==0){//if gradient equals zero, just leave
										VtmpVar2.at(model).at(s)=0;
									}else{
										if ((abs(Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd)){
											//which max did it hit, handle differently
											if ((Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd){
												VtmpVar2.at(model).at(s)=(maxStd-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}else{
												VtmpVar2.at(model).at(s)=(-maxStd-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}
										}else{
											if (Vgradtmp.at(model).at(s)-Vvarstmp.at(model).at(s)-j*searchRes*Vsvec.at(model).at(s)<relStd.at(s)){
												VtmpVar2.at(model).at(s)=(Vgradtmp.at(model).at(s)+relStd.at(s)-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}else{
												VtmpVar2.at(model).at(s)=(Vgradtmp.at(model).at(s)-relStd.at(s)-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}
											
										}
									}
								}
								Vvmax.at(model).at(s)=true;
							}
							
							if (!Vvmax.at(model).at(s)){
								Vvars->at(model).at(s)=Vvarstmp.at(model).at(s) + j*searchRes*Vsvec.at(model).at(s);
							}else{
								Vvars->at(model).at(s)=Vvarstmp.at(model).at(s) + VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
							}
						}
						vmeantmp.push_back(-333);
					}
					float cost=0;
					//only implemented for appearance model
					
					cost=costfuncAppDual(Vmodel1,*Vvars,vmeantmp);
					//cost should only be calculated once for whole model configuration
					
					
					for (int model=0; model<numModels; model++){	
						if (j==0){	     
							VtmpVar.at(model)=j*searchRes;
							VtmpCost.at(model)=cost;
						}else if (cost<VtmpCost.at(model)){
							VtmpVar.at(model)=j*searchRes;
							VtmpCost.at(model)=cost;
							zerocount=0;//reset count when tmpvar is updated
						}else{
							zerocount++;
						}
						
						
						if (zerocount>3){//causes break in linear search if nothing found in 5
							j=searchDist;
						}
					}
					j++;
				}
				//assign values of position and displacment
				vector<float> Vdelta;
				for (int model=0; model<numModels; model++){
					Vdelta.push_back(0);
				}
				for (int model=0; model<numModels; model++){
					//this enables single select
					for (unsigned int s=0; s<Vsvectmp.at(model).size(); s++){ 
						if ((!Vvmax.at(model).at(s))|((Vvmax.at(model).at(s))&&(VtmpVar.at(model)<=VtmpVar2.at(model).at(s)))){
							Vvars->at(model).at(s)=Vvarstmp.at(model).at(s)+VtmpVar.at(model)*Vsvec.at(model).at(s);
							Vsvectmp.at(model).at(s)=VtmpVar.at(model)*Vsvectmp.at(model).at(s);
							Vdelta.at(model)+=VtmpVar.at(model)*Vsvec.at(model).at(s)*VtmpVar.at(model)*Vsvec.at(model).at(s);
						}else{
							Vvars->at(model).at(s)=Vvarstmp.at(model).at(s)+VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
							Vsvectmp.at(model).at(s)=VtmpVar2.at(model).at(s)*Vsvectmp.at(model).at(s);
							Vdelta.at(model)+=VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s)*VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
						}
						
						//delta not actually used in algorithm
					}
					if (Vdelta.at(model)<0.01){
						VtmpVar.at(model)=0;
					}
					
					//write new parameters in neggradient call
					
					//incraese search resolution 
						cout<<model<<" vars:"<<endl;
						for (unsigned int i=0;i<Vvars->at(model).size();i++){
						cout<<Vvars->at(model).at(i)<<" ";
					}
					//cout<<endl;
					
				}//end model iteration
				 //changed to sum of tmpvar this should act on wfole models
				 float sumtmpVar=0;
				 for (int model=0; model<numModels; model++){
					sumtmpVar+=VtmpVar.at(model);
				 }
				if (sumtmpVar==0){
					//play with this maaybe no chane in res
					searchRes=searchRes/2;
					//	cout<<"addvolume "<<endl;
					//	addVolumeTo4D(model1,*vars);
					if (searchRes<searchRmax){
						gradzero=true;
					}
					
				}
				
				if(!gradzero){
					bool breakloop=false;
					
					while(!breakloop){
						cout<<"NEg gradietn here "<<endl;
						negGradientDual(&Vsvectmp,*Vvars,Vmodel1, varbeg, varsize, &Vmaxind,select);
						if (gradzero==true){
							searchRes=searchRes/2;
							if (searchRes<searchRmax){
								breakloop=true;
							}else{
								gradzero=false;
							}
						}else{
							breakloop=true;
						}
					}	
					
				}
				
					cout<<"out of loop "<<endl;
				if (gradzero){
					cout<<"Ive hit the return after "<<iter<<" iterations"<<endl;
				
					//addVolumeTo4D(model1,*vars);
					vector<float> vmeantmp;
					for (int model=0; model<numModels; model++){
						vmeantmp.push_back(-333);
					}
					costfuncAppDual(Vmodel1, *Vvars,vmeantmp);
					//					cout<<"gbeta1 "<<gBeta<<endl;
					gradzero=false;
					
					return ;
				}
				
				//calculate new residual dot product
				//double dpRes=0;
				cout<<"almost have it1"<<endl;
				vector<double> VdpRes;
				 for (int model=0; model<numModels; model++){
					VdpRes.push_back(0);
				 }
				 
				 cout<<"almost have it"<<endl;
				 for (int model=0; model<numModels; model++){
					 //for (int i=0; i<n; i++){
					 for (int i=varbeg; i<varbeg+varsize; i++){
						
						 VdpResPrev.at(model)+=Vres.at(model).at(i)*Vres.at(model).at(i);
						 VdpRes.at(model)+=Vsvectmp.at(model).at(i)*Vsvectmp.at(model).at(i);
						 VdpRes.at(model)+=(-Vsvectmp.at(model).at(i)+Vres.at(model).at(i))*-Vsvectmp.at(model).at(i);
					 }
					  cout<<"almost have it34"<<endl;
					 Vgamma.at(model)= VdpRes.at(model)/VdpResPrev.at(model);
					 
					 cout<<"almost have it3"<<endl;
					 //update update vector
					 
					 VresPrev.at(model)=Vres.at(model);
					 Vres.at(model)=Vsvectmp.at(model);
					 cout<<"HMMMMMMMMMM"<<endl;
					 //so that single selects may be used			
					 for (unsigned int i=0; i<Vsvec.at(model).size(); i++){
						 Vsvectmp.at(model).at(i)=Vsvec.at(model).at(i)=Vres.at(model).at(i)+Vgamma.at(model)*Vsvectmp.at(model).at(i);
						 //This changes conjugate gradient to gradient
					 }
					 }			
				}
			//addVolumeTo4D(model1,*vars);
			}
void conjGradientN(vector< shapeModel* > Vmodel1, vector< vector<float> >* Vvars , vector<Matrix> vmBx2, vector<Matrix> vmBx1inv, vector<int> vKx2,int varbeg, int varsize,  vector<float> relStd, vector<bool> select, float searchR, float searchRmax){
	//this will load 2 models and fit simultaneously
	//will need to keep track of as many variables as models
	//though some is general use of conditionals is only for 2 models at this point
	cout<<"conjGradient Dual"<<endl;
	//FOR NOW DEFINE NUMBER OF MODELS AS 2
	int numModels=static_cast<int>(Vmodel1.size());
	//costtype defines whether to use AAM or ASM
	//define variable
	

	//one for each model
	vector< vector<float> >  Vsvec, Vres, VresPrev;//, pointCost;


		vector< double > Vgamma;
		
		//everything is stored in vectors such that it can be generalized to multiple structures
		vector< int > Vn;
		vector< double > VdpResPrev; 
		vector< vector<float> > Vgradtmp;
		vector< vector<float> > Vsvectmp;
		//vector< vector<float> > Vgradtmp2;
		gradzero=false;
		searchRes=searchR;
		
		vector<int> Vmaxind;
		vector< vector<float> > VvarMax;
		
		vector< vector<bool> > Vvmax;
		vector< vector<float> > VcostVMAX;
		vector< vector<float> > VtmpVar2;
		
		//initializations
		float maxStd=stdTrunc.value();
		cout<<"hmmm1 "<<numModels<<endl;
		for (int model=0; model<numModels; model++){
			Vgamma.push_back(0);
			Vn.push_back( Vvars->at(model).size() );
			
			VdpResPrev.push_back(0); 
			
		//	vector<float> gradtmp2=Vres.at(model);
		//	Vgradtmp2.push_back(gradtmp2);
			
			vector<float> svectmp=(*Vvars).at(model);
			Vsvectmp.push_back(svectmp);
			
			vector<float> gradtmp=(*Vvars).at(model);
			Vgradtmp.push_back(gradtmp);
			
			//this stuff is used to handle cost at truncation
			Vmaxind.push_back(0);
		
			//to find cost at truncation
			vector<float> varMax=(*Vvars).at(model);
			VvarMax.push_back(varMax);
		
			//	searchDim=0.001;
			vector<float> resPrev=(*Vvars).at(model);
			VresPrev.push_back(resPrev);
			
			vector<float> res=(*Vvars).at(model);
			Vres.push_back(res);
			
			vector<float> svec=(*Vvars).at(model);
			Vsvec.push_back(svec);
		
		}
		//neggradientdual call should not be in the model loop
			//select is left the same for now....search same modes for each structure
			//cout<<"hmmm2 "<<endl;

			negGradientN(&Vres,*Vvars,Vmodel1, varbeg, varsize, &Vmaxind, select,vKx2, vmBx2, vmBx1inv );
				cout<<"hmmm3 "<<endl;
			if ((gradzero)){
				return;//if graident == 0
			}
			for (int model=0; model<numModels; model++){
				//cout<<"hmmm3.1 "<<endl;
				VresPrev.at(model)=Vres.at(model);
				//	cout<<"hmmm3.1 "<<endl;
				//Vsvec.at(model)=Vsvectmp.at(model)=Vres.at(model);
				Vsvec.at(model)=Vres.at(model);
				Vsvectmp.at(model)=Vres.at(model);
		//		cout<<"hmmm3.1 "<<endl;
				
				vector<float> tmpVar2;
				for (int i=0;i<Vn.at(model);i++){
					tmpVar2.push_back(0);
				}
				VtmpVar2.push_back(tmpVar2);
		//		cout<<"hmmm3.2 "<<endl;
				
				vector<bool> vmax;
				vector<float> costVMAX;///used to see if upper bound hit and that is the max
					
					for (int i=0;i<Vn.at(model);i++){
						vmax.push_back(false);
						costVMAX.push_back(10e6);
					}
					Vvmax.push_back(vmax);
					VcostVMAX.push_back(costVMAX);
			} ///end of initialization for each model
			//STILL NEED TO DEAL WITH GRADZERO
			//	cout<<"hmmm4 "<<endl;
			for (int iter=0;iter<40;iter++){
				
				//find line mminimizatioon and set svec to vector displavment
				
				//start fitting first mode
			  //	cout<<"iter: "<<iter<<endl;
				vector<float> VtmpVar;
				vector<float> VtmpCost;
				vector< vector<float> > Vvarstmp;
				for (int model=0; model<numModels; model++){
					VtmpVar.push_back(0);
					VtmpCost.push_back(0);
					
					vector<float> varstmp=(*Vvars).at(model);
					Vvarstmp.push_back(varstmp);
					for (int i=0;i<Vn.at(model);i++){
						VtmpVar2.at(model).at(i)=0;
					}
					for (int i=0;i<Vn.at(model);i++){
						Vvmax.at(model).at(i)=false;
					}
				
				}
				
			
				int zerocount=0; //to count number of succesive zeros->help speed search
				int searchDist=60;
				int j=0;
				while ((j<30)){
					//used for input to costfuncappdual
					vector<float> vmeantmp;
					
					//enable single select
					for (int model=0; model<numModels; model++){
						for (unsigned int s=0; s<Vvars->at(model).size(); s++){
							if ((abs(Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd)|(abs(Vgradtmp.at(model).at(s)-Vvarstmp.at(model).at(s)-j*searchRes*Vsvec.at(model).at(s))>relStd.at(s))){
								
								//compare against absoulte allowable maximum std and against relative varaition form start point
								
								//save the value that will render the appropriate max value
								if (!Vvmax.at(model).at(s)){//if not already at its max
									
									if (Vsvec.at(model).at(s)==0){//if gradient equals zero, just leave
										VtmpVar2.at(model).at(s)=0;
									}else{
										if ((abs(Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd)){
											//which max did it hit, handle differently
											if ((Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd){
												VtmpVar2.at(model).at(s)=(maxStd-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}else{
												VtmpVar2.at(model).at(s)=(-maxStd-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}
										}else{
											if (Vgradtmp.at(model).at(s)-Vvarstmp.at(model).at(s)-j*searchRes*Vsvec.at(model).at(s)<relStd.at(s)){
												VtmpVar2.at(model).at(s)=(Vgradtmp.at(model).at(s)+relStd.at(s)-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}else{
												VtmpVar2.at(model).at(s)=(Vgradtmp.at(model).at(s)-relStd.at(s)-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}
											
										}
									}
								}
								Vvmax.at(model).at(s)=true;
							}
							
							if (!Vvmax.at(model).at(s)){
								Vvars->at(model).at(s)=Vvarstmp.at(model).at(s) + j*searchRes*Vsvec.at(model).at(s);
							}else{
								Vvars->at(model).at(s)=Vvarstmp.at(model).at(s) + VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
							}
						}
						vmeantmp.push_back(-333);
					}
					float cost=0;
					//only implemented for appearance model

					cost=costfuncAppN(Vmodel1,*Vvars,vmeantmp,vKx2,vmBx2,vmBx1inv, false, 0);
					//cost should only be calculated once for whole model configuration
					
					
					for (int model=0; model<numModels; model++){	
						if (j==0){	     
							VtmpVar.at(model)=j*searchRes;
							VtmpCost.at(model)=cost;
						}else if (cost<VtmpCost.at(model)){
							VtmpVar.at(model)=j*searchRes;
							VtmpCost.at(model)=cost;
							zerocount=0;//reset count when tmpvar is updated
						}else{
							zerocount++;
						}
						
						
						if (zerocount>3){//causes break in linear search if nothing found in 5
							j=searchDist;
						}
					}
					j++;
				}
				//assign values of position and displacment
				vector<float> Vdelta;
				for (int model=0; model<numModels; model++){
					Vdelta.push_back(0);
				}
				for (int model=0; model<numModels; model++){
					//this enables single select
					for (unsigned int s=0; s<Vsvectmp.at(model).size(); s++){ 
						if ((!Vvmax.at(model).at(s))|((Vvmax.at(model).at(s))&&(VtmpVar.at(model)<=VtmpVar2.at(model).at(s)))){
							Vvars->at(model).at(s)=Vvarstmp.at(model).at(s)+VtmpVar.at(model)*Vsvec.at(model).at(s);
							Vsvectmp.at(model).at(s)=VtmpVar.at(model)*Vsvectmp.at(model).at(s);
							Vdelta.at(model)+=VtmpVar.at(model)*Vsvec.at(model).at(s)*VtmpVar.at(model)*Vsvec.at(model).at(s);
						}else{
							Vvars->at(model).at(s)=Vvarstmp.at(model).at(s)+VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
							Vsvectmp.at(model).at(s)=VtmpVar2.at(model).at(s)*Vsvectmp.at(model).at(s);
							Vdelta.at(model)+=VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s)*VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
						}
						
						//delta not actually used in algorithm
					}
					if (Vdelta.at(model)<0.01){
						VtmpVar.at(model)=0;
					}
					
					//write new parameters in neggradient call
					
					//incraese search resolution 
						cout<<model<<" vars:"<<endl;
						for (unsigned int i=0;i<Vvars->at(model).size();i++){
						cout<<Vvars->at(model).at(i)<<" ";
					}
					//cout<<endl;
					
				}//end model iteration
				 //changed to sum of tmpvar this should act on wfole models
				 float sumtmpVar=0;
				 for (int model=0; model<numModels; model++){
					sumtmpVar+=VtmpVar.at(model);
				 }
				if (sumtmpVar==0){
					//play with this maaybe no chane in res
					searchRes=searchRes/2;
					//	cout<<"addvolume "<<endl;
					//	addVolumeTo4D(model1,*vars);
					if (searchRes<searchRmax){
						gradzero=true;
					}
					
				}
				
				if(!gradzero){
					bool breakloop=false;
					
					while(!breakloop){
						cout<<"NEg gradietn here "<<endl;
						negGradientN(&Vsvectmp,*Vvars,Vmodel1, varbeg, varsize, &Vmaxind,select,vKx2, vmBx2, vmBx1inv );
						if (gradzero==true){
							searchRes=searchRes/2;
							if (searchRes<searchRmax){
								breakloop=true;
							}else{
								gradzero=false;
							}
						}else{
							breakloop=true;
						}
					}	
					
				}
				
					cout<<"out of loop "<<endl;
				if (gradzero){
					cout<<"Ive hit the return after "<<iter<<" iterations"<<endl;
				
					//addVolumeTo4D(model1,*vars);
					vector<float> vmeantmp;
					for (int model=0; model<numModels; model++){
						vmeantmp.push_back(-333);
					}
			

					costfuncAppN(Vmodel1, *Vvars,vmeantmp,vKx2,vmBx2,vmBx1inv, false, 0);
					//					cout<<"gbeta1 "<<gBeta<<endl;
					gradzero=false;
					
					return ;
				}
				
				//calculate new residual dot product
				//double dpRes=0;
				cout<<"almost have it1"<<endl;
				vector<double> VdpRes;
				 for (int model=0; model<numModels; model++){
					VdpRes.push_back(0);
				 }
				 
				 cout<<"almost have it"<<endl;
				 for (int model=0; model<numModels; model++){
					 //for (int i=0; i<n; i++){
					 for (int i=varbeg; i<varbeg+varsize; i++){
						
						 VdpResPrev.at(model)+=Vres.at(model).at(i)*Vres.at(model).at(i);
						 VdpRes.at(model)+=Vsvectmp.at(model).at(i)*Vsvectmp.at(model).at(i);
						 VdpRes.at(model)+=(-Vsvectmp.at(model).at(i)+Vres.at(model).at(i))*-Vsvectmp.at(model).at(i);
					 }
					  cout<<"almost have it34"<<endl;
					 Vgamma.at(model)= VdpRes.at(model)/VdpResPrev.at(model);
					 
					 cout<<"almost have it3"<<endl;
					 //update update vector
					 
					 VresPrev.at(model)=Vres.at(model);
					 Vres.at(model)=Vsvectmp.at(model);
					 cout<<"HMMMMMMMMMM"<<endl;
					 //so that single selects may be used			
				//	 for (unsigned int i=0; i<Vsvec.at(model).size(); i++){
				//		 Vsvectmp.at(model).at(i)=Vsvec.at(model).at(i)=Vres.at(model).at(i)+Vgamma.at(model)*Vsvectmp.at(model).at(i);
						 //This changes conjugate gradient to gradient
				//	 }
					 }			
				}
			//addVolumeTo4D(model1,*vars);
			}
void conjGradientN2(vector< shapeModel* > Vmodel1, vector< vector<float> >* Vvars , vector< vector<Matrix> > vmBx2,  vector<Matrix> vmBx1inv, vector< vector<int> > vPredLabels, vector<int> vKx2,int varbeg, int varsize,  vector<float> relStd, vector<bool> select, float searchR, float searchRmax){
	//this will load 2 models and fit simultaneously
	//will need to keep track of as many variables as models
	//though some is general use of conditionals is only for 2 models at this point
	cout<<"conjGradient Dual"<<endl;
	//FOR NOW DEFINE NUMBER OF MODELS AS 2
	int numModels=static_cast<int>(Vmodel1.size());
	//costtype defines whether to use AAM or ASM
	//define variable
	

	//one for each model
	vector< vector<float> >  Vsvec, Vres, VresPrev;//, pointCost;


		double gamma;
		
		//everything is stored in vectors such that it can be generalized to multiple structures
		vector< int > Vn;
		double dpResPrev; 
		vector< vector<float> > Vgradtmp;
		vector< vector<float> > Vsvectmp;
		//vector< vector<float> > Vgradtmp2;
		gradzero=false;
		searchRes=searchR;
		
		vector<int> Vmaxind;
		vector< vector<float> > VvarMax;
		
		vector< vector<bool> > Vvmax;
		vector< vector<float> > VcostVMAX;
		vector< vector<float> > VtmpVar2;
		
		//initializations
		float maxStd=stdTrunc.value();
		cout<<"hmmm1 "<<numModels<<endl;
		for (int model=0; model<numModels; model++){
			gamma=0;
			Vn.push_back( Vvars->at(model).size() );
			
			dpResPrev=0; 
			
		//	vector<float> gradtmp2=Vres.at(model);
		//	Vgradtmp2.push_back(gradtmp2);
			
			vector<float> svectmp=(*Vvars).at(model);
			Vsvectmp.push_back(svectmp);
			
			vector<float> gradtmp=(*Vvars).at(model);
			Vgradtmp.push_back(gradtmp);
			
			//this stuff is used to handle cost at truncation
			Vmaxind.push_back(0);
		
			//to find cost at truncation
			vector<float> varMax=(*Vvars).at(model);
			VvarMax.push_back(varMax);
		
			//	searchDim=0.001;
			vector<float> resPrev=(*Vvars).at(model);
			VresPrev.push_back(resPrev);
			
			vector<float> res=(*Vvars).at(model);
			Vres.push_back(res);
			
			vector<float> svec=(*Vvars).at(model);
			Vsvec.push_back(svec);
		
		}
		//neggradientdual call should not be in the model loop
			//select is left the same for now....search same modes for each structure
			//cout<<"hmmm2 "<<endl;

			negGradientN2(&Vres,*Vvars,Vmodel1, varbeg, varsize, &Vmaxind, select,vKx2, vmBx2, vmBx1inv, vPredLabels );
				cout<<"hmmm3 "<<endl;
			if ((gradzero)){
				return;//if graident == 0
			}
			for (int model=0; model<numModels; model++){
				//cout<<"hmmm3.1 "<<endl;
				VresPrev.at(model)=Vres.at(model);
				//	cout<<"hmmm3.1 "<<endl;
				//Vsvec.at(model)=Vsvectmp.at(model)=Vres.at(model);
				Vsvec.at(model)=Vres.at(model);
				Vsvectmp.at(model)=Vres.at(model);
		//		cout<<"hmmm3.1 "<<endl;
				
				vector<float> tmpVar2;
				for (int i=0;i<Vn.at(model);i++){
					tmpVar2.push_back(0);
				}
				VtmpVar2.push_back(tmpVar2);
		//		cout<<"hmmm3.2 "<<endl;
				
				vector<bool> vmax;
				vector<float> costVMAX;///used to see if upper bound hit and that is the max
					
					for (int i=0;i<Vn.at(model);i++){
						vmax.push_back(false);
						costVMAX.push_back(10e6);
					}
					Vvmax.push_back(vmax);
					VcostVMAX.push_back(costVMAX);
			} ///end of initialization for each model
			//STILL NEED TO DEAL WITH GRADZERO
			//	cout<<"hmmm4 "<<endl;
			for (int iter=0;iter<40;iter++){
				
				//find line mminimizatioon and set svec to vector displavment
				
				//start fitting first mode
			  //	cout<<"iter: "<<iter<<endl;
				vector<float> VtmpVar;
				vector<float> VtmpCost;
				vector< vector<float> > Vvarstmp;
				for (int model=0; model<numModels; model++){
					VtmpVar.push_back(0);
					VtmpCost.push_back(0);
					
					vector<float> varstmp=(*Vvars).at(model);
					Vvarstmp.push_back(varstmp);
					for (int i=0;i<Vn.at(model);i++){
						VtmpVar2.at(model).at(i)=0;
					}
					for (int i=0;i<Vn.at(model);i++){
						Vvmax.at(model).at(i)=false;
					}
				
				}
				
			
				int zerocount=0; //to count number of succesive zeros->help speed search
				int searchDist=60;
				int j=0;
				while ((j<30)){
					//used for input to costfuncappdual
					vector<float> vmeantmp;
					
					//enable single select
					for (int model=0; model<numModels; model++){
						for (unsigned int s=0; s<Vvars->at(model).size(); s++){
							if ((abs(Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd)|(abs(Vgradtmp.at(model).at(s)-Vvarstmp.at(model).at(s)-j*searchRes*Vsvec.at(model).at(s))>relStd.at(s))){
								
								//compare against absoulte allowable maximum std and against relative varaition form start point
								
								//save the value that will render the appropriate max value
								if (!Vvmax.at(model).at(s)){//if not already at its max
									
									if (Vsvec.at(model).at(s)==0){//if gradient equals zero, just leave
										VtmpVar2.at(model).at(s)=0;
									}else{
										if ((abs(Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd)){
											//which max did it hit, handle differently
											if ((Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd){
												VtmpVar2.at(model).at(s)=(maxStd-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}else{
												VtmpVar2.at(model).at(s)=(-maxStd-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}
										}else{
											if (Vgradtmp.at(model).at(s)-Vvarstmp.at(model).at(s)-j*searchRes*Vsvec.at(model).at(s)<relStd.at(s)){
												VtmpVar2.at(model).at(s)=(Vgradtmp.at(model).at(s)+relStd.at(s)-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}else{
												VtmpVar2.at(model).at(s)=(Vgradtmp.at(model).at(s)-relStd.at(s)-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
											}
											
										}
									}
								}
								Vvmax.at(model).at(s)=true;
							}
							
							if (!Vvmax.at(model).at(s)){
								Vvars->at(model).at(s)=Vvarstmp.at(model).at(s) + j*searchRes*Vsvec.at(model).at(s);
							}else{
								Vvars->at(model).at(s)=Vvarstmp.at(model).at(s) + VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
							}
						}
						vmeantmp.push_back(-333);
					}
					float cost=0;
					//only implemented for appearance model

					cost=costfuncAppN2(Vmodel1,*Vvars,vmeantmp,vKx2,vmBx2,vmBx1inv, vPredLabels, false, 0);
					//cost should only be calculated once for whole model configuration
					for (int model=0; model<numModels; model++){	
						cout<<model<<" varsN2: "<<model<<endl;
						for (unsigned int i=0;i<Vvars->at(model).size();i++){
							cout<<Vvars->at(model).at(i)<<" ";
						}
					}
					for (int model=0; model<numModels; model++){	
						if (j==0){	     
							VtmpVar.at(model)=j*searchRes;
							VtmpCost.at(model)=cost;
						}else if (cost<VtmpCost.at(model)){
							VtmpVar.at(model)=j*searchRes;
							VtmpCost.at(model)=cost;
							zerocount=0;//reset count when tmpvar is updated
						}else{
							zerocount++;
						}
						
						
						if (zerocount>3){//causes break in linear search if nothing found in 5
							j=searchDist;
						}
					}
					j++;
				}
				//assign values of position and displacment
				vector<float> Vdelta;
				for (int model=0; model<numModels; model++){
					Vdelta.push_back(0);
				}
				for (int model=0; model<numModels; model++){
					//this enables single select
					for (unsigned int s=0; s<Vsvectmp.at(model).size(); s++){ 
						if ((!Vvmax.at(model).at(s))|((Vvmax.at(model).at(s))&&(VtmpVar.at(model)<=VtmpVar2.at(model).at(s)))){
							Vvars->at(model).at(s)=Vvarstmp.at(model).at(s)+VtmpVar.at(model)*Vsvec.at(model).at(s);
							Vsvectmp.at(model).at(s)=VtmpVar.at(model)*Vsvectmp.at(model).at(s);
							Vdelta.at(model)+=VtmpVar.at(model)*Vsvec.at(model).at(s)*VtmpVar.at(model)*Vsvec.at(model).at(s);
						}else{
							Vvars->at(model).at(s)=Vvarstmp.at(model).at(s)+VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
							Vsvectmp.at(model).at(s)=VtmpVar2.at(model).at(s)*Vsvectmp.at(model).at(s);
							Vdelta.at(model)+=VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s)*VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
						}
						
						//delta not actually used in algorithm
					}
					if (Vdelta.at(model)<0.01){
						VtmpVar.at(model)=0;
					}
					
					//write new parameters in neggradient call
					
					//incraese search resolution 
						cout<<model<<" vars:"<<endl;
						for (unsigned int i=0;i<Vvars->at(model).size();i++){
						cout<<Vvars->at(model).at(i)<<" ";
					}
					//cout<<endl;
					
				}//end model iteration
				 //changed to sum of tmpvar this should act on wfole models
				 float sumtmpVar=0;
				 for (int model=0; model<numModels; model++){
					sumtmpVar+=VtmpVar.at(model);
				 }
				if (sumtmpVar==0){
					//play with this maaybe no chane in res
					cout<<"NEW RES "<<searchRes<<" to "<<searchRes/2;
					searchRes=searchRes/2;
					//	cout<<"addvolume "<<endl;
					//	addVolumeTo4D(model1,*vars);
					if (searchRes<searchRmax){
						gradzero=true;
					}
					
				}
				
				if(!gradzero){
					bool breakloop=false;
					
					while(!breakloop){
						cout<<"NEg gradietn here "<<endl;
						negGradientN2(&Vsvectmp,*Vvars,Vmodel1, varbeg, varsize, &Vmaxind,select,vKx2, vmBx2, vmBx1inv,vPredLabels );
						if (gradzero==true){
							searchRes=searchRes/2;
							cout<<"NEW RES2 "<<searchRes<<" to "<<searchRes/2;
							if (searchRes<searchRmax){
								breakloop=true;
							}else{
								gradzero=false;
							}
						}else{
							breakloop=true;
						}
					}	
					
				}
				
					cout<<"out of loop "<<endl;
				if (gradzero){
					cout<<"Ive hit the return after "<<iter<<" iterations"<<endl;
				
					//addVolumeTo4D(model1,*vars);
					vector<float> vmeantmp;
					for (int model=0; model<numModels; model++){
						vmeantmp.push_back(-333);
					}
			

					costfuncAppN2(Vmodel1, *Vvars,vmeantmp,vKx2,vmBx2,vmBx1inv, vPredLabels, false, 0);
					//					cout<<"gbeta1 "<<gBeta<<endl;
					gradzero=false;
					
					return ;
				}
				
				//calculate new residual dot product
				//double dpRes=0;
				cout<<"almost have it1"<<endl;
				double dpRes;
				dpRes=0;
		//		  for (int model=0; model<numModels; model++){
		//			VdpRes.push_back(0);
		//		 }
				 
				 cout<<"almost have it"<<endl;
				 for (int model=0; model<numModels; model++){
					 //for (int i=0; i<n; i++){
					 for (int i=varbeg; i<varbeg+varsize; i++){
						
						 dpResPrev+=Vres.at(model).at(i)*Vres.at(model).at(i);
						 dpRes+=Vsvectmp.at(model).at(i)*Vsvectmp.at(model).at(i);
						 dpRes+=(-Vsvectmp.at(model).at(i)+Vres.at(model).at(i))*-Vsvectmp.at(model).at(i);
					 }
					  cout<<"almost have it34"<<endl;
					 gamma= dpRes/dpResPrev;
					 
					 cout<<"almost have it3"<<endl;
					 //update update vector
					 
					 VresPrev.at(model)=Vres.at(model);
					 Vres.at(model)=Vsvectmp.at(model);
					 cout<<"HMMMMMMMMMM"<<endl;
					 //so that single selects may be used			
					 		 for (unsigned int i=0; i<Vsvec.at(model).size(); i++){
					  		 Vsvectmp.at(model).at(i)=Vsvec.at(model).at(i)=Vres.at(model).at(i)+gamma*Vsvectmp.at(model).at(i);
						 //This changes conjugate gradient to gradient
					 	 }
					 }			
				}
			//addVolumeTo4D(model1,*vars);
}
void conjGradientN3(vector< shapeModel* > Vmodel1, vector< vector<float> >* Vvars , vector< vector<Matrix> > vmBx2,  vector<Matrix> vmBx1inv, vector< vector<int> > vPredLabels, vector<int> vKx2,int varbeg, int varsize,  vector<float> relStd, vector<bool> select, float searchR, float searchRmax){
  //this will load 2 models and fit simultaneously
  //will need to keep track of as many variables as models
  //though some is general use of conditionals is only for 2 models at this point
  cout<<"conjGradient Dual"<<endl;
  //FOR NOW DEFINE NUMBER OF MODELS AS 2
  int numModels=static_cast<int>(Vmodel1.size());
  //costtype defines whether to use AAM or ASM
  //define variable
	

  //one for each model
  vector< vector<float> >  Vsvec, Vres, VresPrev;//, pointCost;


  vector<double> Vgamma;
		
  //everything is stored in vectors such that it can be generalized to multiple structures
  vector< int > Vn;
  vector<double> VdpResPrev; 
  vector< vector<float> > Vgradtmp;
  vector< vector<float> > Vsvectmp;
  //vector< vector<float> > Vgradtmp2;
  gradzero=false;
  searchRes=searchR;
		
  vector<int> Vmaxind;
  vector< vector<float> > VvarMax;
		
  vector< vector<bool> > Vvmax;
  vector< vector<float> > VcostVMAX;
  vector< vector<float> > VtmpVar2;
		
  //initializations
  float maxStd=stdTrunc.value();
  cout<<"hmmm1 "<<numModels<<endl;
  for (int model=0; model<numModels; model++){
    Vgamma.push_back(0);
    Vn.push_back( Vvars->at(model).size() );
			
    VdpResPrev.push_back(0); 
			
    //	vector<float> gradtmp2=Vres.at(model);
    //	Vgradtmp2.push_back(gradtmp2);
			
    vector<float> svectmp=(*Vvars).at(model);
    Vsvectmp.push_back(svectmp);
			
    vector<float> gradtmp=(*Vvars).at(model);
    Vgradtmp.push_back(gradtmp);
			
    //this stuff is used to handle cost at truncation
    Vmaxind.push_back(0);
		
    //to find cost at truncation
    vector<float> varMax=(*Vvars).at(model);
    VvarMax.push_back(varMax);
		
    //	searchDim=0.001;
    vector<float> resPrev=(*Vvars).at(model);
    VresPrev.push_back(resPrev);
			
    vector<float> res=(*Vvars).at(model);
    Vres.push_back(res);
			
    vector<float> svec=(*Vvars).at(model);
    Vsvec.push_back(svec);
		
  }
  //neggradientdual call should not be in the model loop
  //select is left the same for now....search same modes for each structure
  //cout<<"hmmm2 "<<endl;

  negGradientN2(&Vres,*Vvars,Vmodel1, varbeg, varsize, &Vmaxind, select,vKx2, vmBx2, vmBx1inv, vPredLabels );
  cout<<"hmmm3 "<<endl;
  if ((gradzero)){
    return;//if graident == 0
  }
  for (int model=0; model<numModels; model++){
    //cout<<"hmmm3.1 "<<endl;
    VresPrev.at(model)=Vres.at(model);
    //	cout<<"hmmm3.1 "<<endl;
    //Vsvec.at(model)=Vsvectmp.at(model)=Vres.at(model);
    Vsvec.at(model)=Vres.at(model);
    Vsvectmp.at(model)=Vres.at(model);
    //		cout<<"hmmm3.1 "<<endl;
				
    vector<float> tmpVar2;
    for (int i=0;i<Vn.at(model);i++){
      tmpVar2.push_back(0);
    }
    VtmpVar2.push_back(tmpVar2);
    //		cout<<"hmmm3.2 "<<endl;
				
    vector<bool> vmax;
    vector<float> costVMAX;///used to see if upper bound hit and that is the max
					
    for (int i=0;i<Vn.at(model);i++){
      vmax.push_back(false);
      costVMAX.push_back(10e6);
    }
    Vvmax.push_back(vmax);
    VcostVMAX.push_back(costVMAX);
  } ///end of initialization for each model
  //STILL NEED TO DEAL WITH GRADZERO
  //	cout<<"hmmm4 "<<endl;
  for (int iter=0;iter<40;iter++){
				
    //find line mminimizatioon and set svec to vector displavment
				
    //start fitting first mode
    //	cout<<"iter: "<<iter<<endl;
    vector<float> VtmpVar;
    vector<float> VtmpCost;
    vector< vector<float> > Vvarstmp;
    vector<float> Vdelta;
    for (int tmodel=0; tmodel<numModels; tmodel++){
      Vdelta.push_back(0);
    }
    for (int model=0; model<numModels; model++){
      VtmpVar.push_back(0);
      VtmpCost.push_back(0);
					
      vector<float> varstmp=(*Vvars).at(model);
      Vvarstmp.push_back(varstmp);
      for (int i=0;i<Vn.at(model);i++){
	VtmpVar2.at(model).at(i)=0;
      }
      for (int i=0;i<Vn.at(model);i++){
	Vvmax.at(model).at(i)=false;
      }
				
    }
				
			
    int zerocount=0; //to count number of succesive zeros->help speed search
    int searchDist=60;
    int j=0;
    while ((j<30)){
      //used for input to costfuncappdual
      vector<float> vmeantmp;
      float cost=0;
      //enable single select
      for (int mtemp=0;mtemp<numModels;mtemp++){
	vmeantmp.push_back(-333);
      }
      for (int model=0; model<numModels; model++){
	for (unsigned int s=0; s<Vvars->at(model).size(); s++){
	  if ((abs(Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd)|(abs(Vgradtmp.at(model).at(s)-Vvarstmp.at(model).at(s)-j*searchRes*Vsvec.at(model).at(s))>relStd.at(s))){
								
	    //compare against absoulte allowable maximum std and against relative varaition form start point
								
	    //save the value that will render the appropriate max value
	    if (!Vvmax.at(model).at(s)){//if not already at its max
									
	      if (Vsvec.at(model).at(s)==0){//if gradient equals zero, just leave
		VtmpVar2.at(model).at(s)=0;
	      }else{
		if ((abs(Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd)){
		  //which max did it hit, handle differently
		  if ((Vvarstmp.at(model).at(s)+ j*searchRes*Vsvec.at(model).at(s))>maxStd){
		    VtmpVar2.at(model).at(s)=(maxStd-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
		  }else{
		    VtmpVar2.at(model).at(s)=(-maxStd-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
		  }
		}else{
		  if (Vgradtmp.at(model).at(s)-Vvarstmp.at(model).at(s)-j*searchRes*Vsvec.at(model).at(s)<relStd.at(s)){
		    VtmpVar2.at(model).at(s)=(Vgradtmp.at(model).at(s)+relStd.at(s)-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
		  }else{
		    VtmpVar2.at(model).at(s)=(Vgradtmp.at(model).at(s)-relStd.at(s)-Vvarstmp.at(model).at(s))/Vsvec.at(model).at(s);
		  }
											
		}
	      }
	    }
	    Vvmax.at(model).at(s)=true;
	  }
							
	  if (!Vvmax.at(model).at(s)){
	    Vvars->at(model).at(s)=Vvarstmp.at(model).at(s) + j*searchRes*Vsvec.at(model).at(s);
	  }else{
	    Vvars->at(model).at(s)=Vvarstmp.at(model).at(s) + VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
	  }
	}
					
					
				
	//only implemented for appearance model
	//	cout<<"costN3 "<<endl;
	cost=costfuncAppN2(Vmodel1,*Vvars,vmeantmp,vKx2,vmBx2,vmBx1inv, vPredLabels, false, 0);
	//	cout<<"costN3 done "<<endl;
					
	//cost should only be calculated once for whole model configuration
	//	for (int tmodel=0; tmodel<numModels; tmodel++){	
	//cout<<tmodel<<" varsN2: "<<tmodel<<endl;
	//  for (unsigned int i=0;i<Vvars->at(tmodel).size();i++){
	//   cout<<Vvars->at(tmodel).at(i)<<" ";
	//  }
	//	}
				
	if (j==0){	     
	  VtmpVar.at(model)=j*searchRes;
	  VtmpCost.at(model)=cost;
	}else if (cost<VtmpCost.at(model)){
	  VtmpVar.at(model)=j*searchRes;
	  VtmpCost.at(model)=cost;
	  zerocount=0;//reset count when tmpvar is updated
	}else{
	  zerocount++;
	}
						
						
	if (zerocount>3){//causes break in linear search if nothing found in 5
	  j=searchDist;
	}
      }
      j++;
    }
    //assign values of position and displacment
				
    for (int tmodel=0; tmodel<numModels; tmodel++){
      Vdelta.at(tmodel)=0;
    }
    for (int model=0; model<numModels; model++){
      //this enables single select
      for (unsigned int s=0; s<Vsvectmp.at(model).size(); s++){ 
	if ((!Vvmax.at(model).at(s))|((Vvmax.at(model).at(s))&&(VtmpVar.at(model)<=VtmpVar2.at(model).at(s)))){
	  Vvars->at(model).at(s)=Vvarstmp.at(model).at(s)+VtmpVar.at(model)*Vsvec.at(model).at(s);
	  Vsvectmp.at(model).at(s)=VtmpVar.at(model)*Vsvectmp.at(model).at(s);
	  Vdelta.at(model)+=VtmpVar.at(model)*Vsvec.at(model).at(s)*VtmpVar.at(model)*Vsvec.at(model).at(s);
	}else{
	  Vvars->at(model).at(s)=Vvarstmp.at(model).at(s)+VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
	  Vsvectmp.at(model).at(s)=VtmpVar2.at(model).at(s)*Vsvectmp.at(model).at(s);
	  Vdelta.at(model)+=VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s)*VtmpVar2.at(model).at(s)*Vsvec.at(model).at(s);
	}
						
	//delta not actually used in algorithm
      }
      if (Vdelta.at(model)<0.01){
	VtmpVar.at(model)=0;
      }
					
      //write new parameters in neggradient call
					
      //incraese search resolution 
      cout<<model<<" vars:"<<endl;
      for (unsigned int i=0;i<Vvars->at(model).size();i++){
	cout<<Vvars->at(model).at(i)<<" ";
      }
      //cout<<endl;
					
    }//end model iteration
    //changed to sum of tmpvar this should act on wfole models
    float sumtmpVar=0;
    for (int model=0; model<numModels; model++){
      sumtmpVar+=VtmpVar.at(model);
    }
    if (sumtmpVar==0){
      //play with this maaybe no chane in res
      cout<<"NEW RES "<<searchRes<<" to "<<searchRes/2;
      searchRes=searchRes/2;
      //	cout<<"addvolume "<<endl;
      //	addVolumeTo4D(model1,*vars);
      if (searchRes<searchRmax){
	gradzero=true;
      }
					
    }
				
    if(!gradzero){
      bool breakloop=false;
					
      while(!breakloop){
	cout<<"NEg gradietn here "<<endl;
	negGradientN2(&Vsvectmp,*Vvars,Vmodel1, varbeg, varsize, &Vmaxind,select,vKx2, vmBx2, vmBx1inv,vPredLabels );
	if (gradzero==true){
	  searchRes=searchRes/2;
	  cout<<"NEW RES2 "<<searchRes<<" to "<<searchRes/2;
	  if (searchRes<searchRmax){
	    breakloop=true;
	  }else{
	    gradzero=false;
	  }
	}else{
	  breakloop=true;
	}
      }	
					
    }
				
    cout<<"out of loop "<<endl;
    if (gradzero){
      cout<<"Ive hit the return after "<<iter<<" iterations"<<endl;
				
      //addVolumeTo4D(model1,*vars);
      vector<float> vmeantmp;
      for (int model=0; model<numModels; model++){
	vmeantmp.push_back(-333);
      }
			

      costfuncAppN2(Vmodel1, *Vvars,vmeantmp,vKx2,vmBx2,vmBx1inv, vPredLabels, false, 0);
      //					cout<<"gbeta1 "<<gBeta<<endl;
      gradzero=false;
					
      return ;
    }
				
    //calculate new residual dot product
    //double dpRes=0;
    cout<<"almost have it1"<<endl;
    vector<double> VdpRes;
    //dpRes=0;
    		  for (int model=0; model<numModels; model++){
    			VdpRes.push_back(0);
    		 }
				 
    cout<<"almost have it"<<endl;
    for (int model=0; model<numModels; model++){
      //for (int i=0; i<n; i++){
      for (int i=varbeg; i<varbeg+varsize; i++){
						
	VdpResPrev.at(model)+=Vres.at(model).at(i)*Vres.at(model).at(i);
	VdpRes.at(model)+=Vsvectmp.at(model).at(i)*Vsvectmp.at(model).at(i);
	VdpRes.at(model)+=(-Vsvectmp.at(model).at(i)+Vres.at(model).at(i))*-Vsvectmp.at(model).at(i);
      }
      cout<<"almost have it34"<<endl;
      Vgamma.at(model)= VdpRes.at(model)/VdpResPrev.at(model);
					 
      cout<<"almost have it3"<<endl;
      //update update vector
					 
      VresPrev.at(model)=Vres.at(model);
      Vres.at(model)=Vsvectmp.at(model);
      cout<<"HMMMMMMMMMM"<<endl;
      //so that single selects may be used			
      for (unsigned int i=0; i<Vsvec.at(model).size(); i++){
	Vsvectmp.at(model).at(i)=Vsvec.at(model).at(i)=Vres.at(model).at(i)+Vgamma.at(model)*Vsvectmp.at(model).at(i);
	//This changes conjugate gradient to gradient
      }
    }			
  }
  //addVolumeTo4D(model1,*vars);
}
float fitModel(shapeModel* model1, string modelName, vector<float>* vars, vector<float>& relStd, int start, int length, int excl, vector<bool> select, int costtype, bool newmodel, float searchR, float searchRmax){
  cout<<"fitting "<<modelName<<endl;
	
	
  if (newmodel){
    model1->clear();
	
    //only use appearnace models now
    model1->load_bmv(modelName,1);
    cout<<"model successfully loaded"<<endl;
  }

  //USe this portion to set a prior on mode
  if (modeprior.value()){
    globModePriorMean=149.2094;
    globModePriorVar=40;//637.0521;
  }
	
  //structure weighting term
  vector<float> fvals;
  for (int i=0; i<model1->getNumberOfShapes();i++){
    fvals.push_back(f.value());
  }
		
  globReest=reestimate.value();
  if (globReest){
    //find mode using less modes of variation
    //this used for wbir	
    int nser=2;
    vector<bool> select2;
    for (int i=0;i<static_cast<int>(vars->size());i++){
      if (i<nser){
	select2.push_back(true);
      }else{
	select2.push_back(false);
      }
    }	
			
    float modeval=conjGradient(model1, vars,0,nser,excl, relStd,fvals,select2, costtype, searchR, searchRmax);
    nser=10;
    //vector<bool> select2;
    for (int i=0;i<static_cast<int>(vars->size());i++){
      if (i<nser){
	select2.at(i)=(true);
      }else{
	select2.at(i)=(false);
      }
    }	
    //GmanMean=-777;
    for (unsigned int i=0; i<vars->size();i++){
      vars->at(i)=0;
      //cout<<vars->at(i)<<" ";
    }

    modeval=conjGradient(model1, vars,0,nser,excl, relStd,fvals,select2, costtype, searchR, searchRmax);
    cout<<"THE MODE IS "<<modeval<<endl;
    //this estimates the mode from the results of the lower mode search
    if (GmanMean!=-777){
      //this is all for single structure
				volume<short> mask;
				int bounds[6]={0,0,0,0,0,0};
				Mesh m = model1->getDeformedMesh(*vars,0,static_cast<int>(vars->size()));
				mask=make_mask_from_meshInOut(image,m,model1->getLabel(0),bounds);
				vector<float> vintens;
				model1->intensityHist(&image,&mask,&m,model1->getLabel(0),&vintens);
			
				
				modeval=mode(vintens);

			cout<<"THE MODE IS "<<modeval<<endl;
		}

		//reset vars to zero
		//set as permanent vlaue
		GmanMean=modeval;
		cout<<"vars "<<endl;
		for (unsigned int i=0; i<vars->size();i++){
				vars->at(i)=0;
		//cout<<vars->at(i)<<" ";
			}
		globReest=false;
			globModePriorMean=-777;
	}
	//globReest=true;
	//GmanMean=-777;
	conjGradient(model1, vars,start,length,excl, relStd,fvals,select, costtype, searchR, searchRmax);

	

 	ofstream fout;
	string name=outname.value()+".bvars";
	fout.open(name.c_str());
	  fout<<"this is a bvars file"<<endl; 
  fout<<modelname.value()<<endl;
  fout<<"NumberOfSubjects "<<1<<endl;
  fout<<inname.value()<<" ";
	fout<<vars->size()<<" ";
	#ifdef PPC64
    int n=0;
#endif
	for (unsigned int i=0;i<vars->size();i++){
		fout<<vars->at(i)<<" ";
		#ifdef PPC64
        if ((n++ % 50) == 0) fout.flush();
		#endif
	}
			fout<<endl;
			fout.close();


}
void fitModelAff(shapeModel* model1, string modelName, vector<float>* vars, vector<float>& relStd, int start, int length, int excl, vector<bool> select, int costtype, bool newmodel, float searchR, float searchRmax){
	cout<<"fitting "<<modelName<<endl;
	
	
	if (newmodel){
		model1->clear();
	
		//only use appearnace models now
		model1->load_bmv(modelName,1);
		cout<<"model successfully loaded"<<endl;
	}
	
	//structure weighting term
	vector<float> fvals;
	for (int i=0; i<model1->getNumberOfShapes();i++){
		fvals.push_back(f.value());
	}
		
	globReest=reestimate.value();
	if (globReest){
		//find mode using less modes of variation
		//this used for wbir	
			int nser=10;
			vector<bool> select2;
			for (int i=0;i<static_cast<int>(vars->size());i++){
				if (i<nser){
					select2.push_back(true);
					}else{
					select2.push_back(false);
					}
			}	
			
		float modeval=conjGradientAff(model1, vars,0,nser,excl, relStd,fvals,select2, costtype, searchR, searchRmax);
		cout<<"THE MODE IS "<<modeval<<endl;
		//this estimates the mode from the results of the lower mode search
		if (GmanMean!=-777){
				//this is all for single structure
				volume<short> mask;
				int bounds[6]={0,0,0,0,0,0};
				Mesh m = model1->getDeformedMesh(*vars,0,static_cast<int>(vars->size()));
				mask=make_mask_from_meshInOut(image,m,model1->getLabel(0),bounds);
				vector<float> vintens;
				model1->intensityHist(&image,&mask,&m,model1->getLabel(0),&vintens);
			
				
				modeval=mode(vintens);

			cout<<"THE MODE IS "<<modeval<<endl;
		}

		//reset vars to zero
		GmanMean=modeval;
		cout<<"vars "<<endl;
		for (unsigned int i=0; i<vars->size();i++){
				//vars->at(i)=0;
		cout<<vars->at(i)<<" ";
			}
		globReest=false;
	}
	conjGradientAff(model1, vars,start,length,excl, relStd,fvals,select, costtype, searchR, searchRmax);

 	ofstream fout;
	string name=outname.value()+".bvars";
	fout.open(name.c_str());
	  fout<<"this is a bvars file"<<endl; 
  fout<<modelname.value()<<endl;
  fout<<"NumberOfSubjects "<<1<<endl;
  fout<<inname.value()<<" ";
	fout<<vars->size()<<" ";
	#ifdef PPC64
    int n=0;
#endif
	for (unsigned int i=0;i<vars->size();i++){
		fout<<vars->at(i)<<" ";
		#ifdef PPC64
        if ((n++ % 50) == 0) fout.flush();
		#endif
	}
			fout<<endl;
			fout.close();


}

	
	
void fitModelDual(vector< shapeModel* > Vmodel1, vector< vector<float> >* Vvars, vector<float>& relStd, int start, int length, vector<bool> select, float searchR, float searchRmax){
	//cout<<"fitting "<<modelName<<endl;
	Vglobfoundmode.clear();
	VglobMean.clear();
	for (unsigned int i=0; i<Vmodel1.size();i++){
		VglobMean.push_back(0);
			Vglobfoundmode.push_back(false);
	}
	
	//dont need newmodel stuff...since we're doing jointly
	//ignore reestimation stuff for now
	conjGradientDual(Vmodel1, Vvars,start,length, relStd, select, searchR, searchRmax);


	//will have to write for multiple vars
	for (unsigned int model=0; model<Vmodel1.size();model++){
		ofstream fout;
		string name=outname.value()+".bvars";
		fout.open(name.c_str());
		fout<<"this is a bvars file"<<endl; 
		fout<<modelname.value()<<endl;
		fout<<"NumberOfSubjects "<<1<<endl;
		fout<<inname.value()<<" ";
		fout<<Vvars->at(model).size()<<" ";
		#ifdef PPC64
		int n=0;
		#endif
		for (unsigned int i=0;i<Vvars->at(model).size();i++){
			fout<<Vvars->at(model).at(i)<<" ";
			#ifdef PPC64
			if ((n++ % 50) == 0) fout.flush();
			#endif
		}
		fout<<endl;
		fout.close();
	}


}
void fitModelN(vector< shapeModel* > Vmodel1, vector< vector<float> >* Vvars, vector<Matrix> vmBx2, vector<Matrix> vmBx1inv, vector<int> vKx2, vector<float>& relStd, int start, int length, vector<bool> select, float searchR, float searchRmax,vector<string> modelNames){
	cout<<"fitting fitmodelN"<<endl;
	Vglobfoundmode.clear();
	VglobMean.clear();
	for (unsigned int i=0; i<Vmodel1.size();i++){
		VglobMean.push_back(0);
			Vglobfoundmode.push_back(false);
	}
	cout<<"enter conjgrad"<<endl;
	
	//dont need newmodel stuff...since we're doing jointly
	//ignore reestimation stuff for now
	conjGradientN(Vmodel1, Vvars, vmBx2,vmBx1inv, vKx2,start,length, relStd, select, searchR, searchRmax);


	//will have to write for multiple vars
	for (unsigned int model=0; model<Vmodel1.size();model++){
		ofstream fout;
		stringstream modelNum;
		modelNum<<model;
		string numstr;
		modelNum>>numstr;
		string name=outname.value()+numstr+".bvars";
		fout.open(name.c_str());
		fout<<"this is a bvars file"<<endl; 
		fout<<modelNames.at(model)<<endl;
		fout<<"NumberOfSubjects "<<1<<endl;
		fout<<inname.value()<<" ";
		fout<<Vvars->at(model).size()<<" ";
		#ifdef PPC64
		int n=0;
		#endif
		for (unsigned int i=0;i<Vvars->at(model).size();i++){
			fout<<Vvars->at(model).at(i)<<" ";
			#ifdef PPC64
			if ((n++ % 50) == 0) fout.flush();
			#endif
		}
		fout<<endl;
		fout.close();
	}


}
void fitModelN2(vector< shapeModel* > Vmodel1, vector< vector<float> >* Vvars, vector< vector<Matrix> > VvmBx2, vector<Matrix> vmBx1inv, vector< vector<int> > vPredLabels, vector<int> vKx2, vector<float>& relStd, int start, int length, vector<bool> select, float searchR, float searchRmax,vector<string> modelNames){
	cout<<"fitting fitmodelN"<<endl;
	Vglobfoundmode.clear();
	VglobMean.clear();
	for (unsigned int i=0; i<Vmodel1.size();i++){
		VglobMean.push_back(0);
			Vglobfoundmode.push_back(false);
	}
	cout<<"enter conjgrad"<<endl;
	
	//dont need newmodel stuff...since we're doing jointly
	//ignore reestimation stuff for now
	if (!useConj3.value()){
	  conjGradientN2(Vmodel1, Vvars, VvmBx2,vmBx1inv, vPredLabels, vKx2,start,length, relStd, select, searchR, searchRmax);}
	else{
	conjGradientN3(Vmodel1, Vvars, VvmBx2,vmBx1inv, vPredLabels, vKx2,start,length, relStd, select, searchR, searchRmax);
	}

	//will have to write for multiple vars
	for (unsigned int model=0; model<Vmodel1.size();model++){
		ofstream fout;
		stringstream modelNum;
		modelNum<<model;
		string numstr;
		modelNum>>numstr;
		string name=outname.value()+numstr+".bvars";
		fout.open(name.c_str());
		fout<<"this is a bvars file"<<endl; 
		fout<<modelNames.at(model)<<endl;
		fout<<"NumberOfSubjects "<<1<<endl;
		fout<<inname.value()<<" ";
		fout<<Vvars->at(model).size()<<" ";
		#ifdef PPC64
		int n=0;
		#endif
		for (unsigned int i=0;i<Vvars->at(model).size();i++){
			fout<<Vvars->at(model).at(i)<<" ";
			#ifdef PPC64
			if ((n++ % 50) == 0) fout.flush();
			#endif
		}
		fout<<endl;
		fout.close();
	}


}

string read_bvars(string fname,vector<float>* bvars,int M){
  string stemp;
  string modelNames;
  int N;//number of subjects
  ifstream fin;
  fin.open(fname.c_str());
  //throw away first three lines 
  getline(fin,stemp);//this is bvars file
  getline(fin,modelNames);//modelnames
  fin>>stemp>>N;
  bvars->clear();
 
  
  //transform all the bvars
  for (int i=0; i<N;i++){
    fin>>stemp;//read in subject id
   
    int nvars;//how many vars written for the subject
    fin>>nvars;
 
       cout<<"subject "<<i<<" "<<stemp<<" "<<nvars<<endl;
    for (int j=0;j<M;j++){
      if (j<nvars){
	float ftemp;
	fin>>ftemp;
	bvars->push_back(ftemp);
      }//else{
      //bvars->element(j,i)=0;
      //}
    }
  }
  return modelNames;
}

int do_work(int argc, char* argv[]) 
{ 
	//load base volume

	read_volume(image,inname.value());
	//normalize image intensities
	cout<<"normalize intensity..."<<endl;

	image=(image-robmin.value())*255/(robmax.value()-robmin.value());
	//cout<<image.getextrapolation()<<endl;
	//image.setpadvalue(0);
	//save_volume(image,"intnorm");
//	if (manMean.value()==-777){
		GmanMean=manMean.value();
//	}else{
//		GmanMean=(manMean.value()-robmin.value())*255/(robmax.value()-robmin.value());
//	}



	
	read_volume(Resimage,inname.value());
	volume<short> brainSeg;
	Resimage=0;
	copyconvert(image,brainSeg);
	brainSeg=0;
	
	const int sizex = image.xsize();
	const int sizey = image.ysize();
	const int sizez = image.zsize();
	xdim=image.xdim();
	ydim=image.ydim();
	zdim=image.zdim();
	
	shapeModel* model1 = new shapeModel;
	//mni152 1mm isotropic
	model1->setImageParameters(182,218, 182, 1, 1, 1);
	//model1->setImageParameters(image.xsize(), image.ysize(), image.zsize(), image.xdim(), image.ydim(), image.zdim());
	cout<<"load model..."<<endl;
		  model1->load_bmv_binaryInfo(modelname.value(),1);
		  model1->load_bmv_binary(modelname.value(),1);	
		  cout<<"register model..."<<endl;
		  //this will set new image parameters
		  model1->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );
		  
		   if (useIntRefModel.value()){
		  shapeModel* modelRef =new shapeModel;
		 	modelRef->setImageParameters(182,218, 182, 1, 1, 1);

			cout<<"use a reference model "<<modelname2.value()<<endl;
			  //this bit of does uses one strcuture as a reference for others.
			  modelRef->load_bmv_binaryInfo(modelname2.value(),1);
			  modelRef->load_bmv_binary(modelname2.value(),1);	
			  cout<<"register model..."<<endl;
			  //this will set new image parameters
			  	  vector<float> varstemp;
			  for (int i=0; i<modelRef->getNumberOfModes();i++){
				  //perturbe the system
				 // cout<<"set tup vars"<<endl;
				  varstemp.push_back(0);
			  }
			  
			  modelRef->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );
			  vector<float> relStd;
			  for (unsigned int i=0;i<varstemp.size();i++){
				  relStd.push_back(stdTrunc.value());
			  }
		
			  //this used for wbir	
			  vector<bool> select;
			  for (unsigned int i=0;i<varstemp.size();i++){
				  if(i<10){
					  select.push_back(true);
				  }else{
					  select.push_back(false);
				  }
			  }		
			  //globshcond enbales one to turn off conditional
			   GlobShCond=false;
		  fitModel(modelRef,modelname2.value(),&varstemp,relStd,0,10,0,select,1, false,0.5, 0.15);
		  GlobShCond=true;
			  Mesh mtemp = modelRef->getDeformedMesh(varstemp,0,static_cast<int>(varstemp.size()));
			  int bounds[6]={0,0,0,0,0,0};
			  volume<short> mask;
			  mask=make_mask_from_meshInOut(image,mtemp,modelRef->getLabel(0),bounds);
			  
			  //calculate intensity histogram
			  float maxint=0,minint=0;
			  vector<float> vintens;
			  modelRef->intensityHistMaxMin(&image,&mask,&mtemp,modelRef->getLabel(0),&vintens, &maxint, &minint);
			  
			  GmanMean=mode(vintens,minint,maxint);
			  cout<<"the mode of the reference distribution is "<<GmanMean<<endl;
			  
		  }
		  

		vector<float> vars;
		for (int i=0; i<model1->getNumberOfModes();i++){
			//perturbe the system
			vars.push_back(0);
		}
	//	addVolumeTo4D(model1,vars);
		searchDim=0.2;
		
		
		
		if (loadvars.value()){
	
		  read_bvars(bvarsname.value(),&vars,model1->getNumberOfModes());
		
		}
		
		if (useCondMode.value()){
		
			shapeModel* model2 = new shapeModel;
			//mni152 1mm isotropic
			model2->setImageParameters(182,218, 182, 1, 1, 1);
			cout<<"load model..."<<endl;
			
			//assumes appearance models
			model2->load_bmv_binaryInfo(modelname2.value(),1);
			model2->load_bmv_binary(modelname2.value(),1);
			cout<<"register model..."<<endl;
			//this will set new image parameters
			model2->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );

			vector<float> varstemp;
			read_bvars(bvarsname.value(),&varstemp,model2->getNumberOfModes());
			
			Mesh mtemp = model1->getDeformedMesh(varstemp,0,static_cast<int>(varstemp.size()));
				int bounds[6]={0,0,0,0,0,0};
			volume<short> mask;
				mask=make_mask_from_meshInOut(image,mtemp,model2->getLabel(0),bounds);
					
				//calculate intensity histogram
				float maxint=0,minint=0;
				vector<float> vintens;
				model2->intensityHistMaxMin(&image,&mask,&mtemp,model2->getLabel(0),&vintens, &maxint, &minint);
			
				float modePri=mode(vintens,minint,maxint);
				
				//evaluate conditional mode (assuem gaussian)
				//for now implement hipp given put
				condMode=149.2094-0.7638*(modePri-184.5554);
			GmanMean=condMode;
			cout<<"conditional mode is "<<GmanMean<<endl;
		}
		
	    
		if (shcond.value()){
		//	cout<<"shcond.value"<<endl;
		int M;
				Matrix mBx2;
				Matrix mBcx1;
		//		cout<<"readbmap "<<bmapname.value()<<endl;
				M=readBmap(bmapname.value(),vars,&mBx2,&mBcx1);
		//		cout<<"readbmap done"<<endl;
			//	int nummodes=
				ColumnVector Bx2(model1->getNumberOfModes());
				//load bx2 vars that were read
				for (unsigned int i=0;i<vars.size();i++){
					Bx2.element(i)=vars.at(i);
				}
				mBx2map=mBx2*Bx2;
				//bmap now directly reflects the transformation t0oo bc
				//mBx1inv=mBcx1;
				mBx1inv=mBcx1.i();
				
				vector<float> v_cmean;
				v_cmean=bTransform(&vars,mBx2,M);
				cout<<"btransfomred"<<endl;
		
				for (int i=0;i<static_cast<int>(vars.size());i++){
				  if (i<g.value()){
					vars.at(i)=v_cmean.at(i);
					cout<<vars.at(i)<<" ";
					}else{
					vars.at(i)=0;
					}
				}
				cout<<endl;
					//vars=v_cmean;
				cout<<"done loading bmap stuff"<<endl;
				if (shcond2.value()){
					vector<float> vars2;
					vars2=vars;
					read_bvars(bvarsname2.value(),&vars,model1->getNumberOfModes());
					cout<<"readbmap2 "<<bmapname2.value()<<endl;
					M=readBmap2(bmapname2.value(),vars,&mBx2,&mBcx1);
					cout<<"readbmap2 done"<<endl;
					//	int nummodes=
					ColumnVector Bx2(model1->getNumberOfModes());
					//load bx2 vars that were read
					for (unsigned int i=0;i<vars.size();i++){
						Bx2.element(i)=vars.at(i);
					}
					mBx2map2=mBx2*Bx2;
					mBx1inv2=mBcx1.i();
					
					vector<float> v_cmean;
					v_cmean=bTransform(&vars,mBx2,M);
					cout<<"btransfomred"<<endl;
					for (int i=0;i<static_cast<int>(vars.size());i++){
						if (i<-1){
							vars.at(i)+=v_cmean.at(i)+vars2.at(i);
						}else{
							vars.at(i)=0;
						}
					}
					//vars=v_cmean;
					cout<<"done loading bmap stuff"<<endl;
					
				}	
				
					
					
		}
	//addVolumeTo4D(model1,vars);
		vector<float> relStd;
		for (unsigned int i=0;i<vars.size();i++){
			relStd.push_back(stdTrunc.value());
		}
		
		string modelName;
		modelName=modelname.value();
			int lb=0,ub=0;
			
			
			//this used for wbir	
			vector<bool> select;
			for (unsigned int i=0;i<vars.size();i++){
				select.push_back(false);
			}		
			
			for (int e=0;e<f.value();e++){
			 
				  for (int shift=0;shift<c.value()*static_cast<int>(g.value());shift+=static_cast<int>(g.value())){
				   				
					lb=0;
					ub=shift+g.value()-1;
					for (int i=0;i<static_cast<int>(vars.size());i++){
						if ((i>=lb)&&(i<=ub)){
							select.at(i)=true;
						}else{
							select.at(i)=false;
						}
					}
				
						costmodeG=0;
						if (baam.value()){
							//use appearance model
					//		fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,1, false,0.15, 0.1);
						fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,1, false,0.5, 0.15);
						//	if (!nograd.value()){
						//	fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,0, false,0.15, 0.1);
							//fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,1, false,0.075, 0.05);
						//	}
						}else if (overlap.value()){
							//uses dice overlap---> work for label data
							fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,2, false,0.15, 0.1);	

						}else{ 
							//use asm 
							fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,0, false,0.15, 0.1);	
 
						}						
				}
			}
		
 			//save_volume4D(fit,"progression");
			addVolumeTo4D(model1,vars);
			save_volume(d4tmp,outname.value());
			Mesh mout = model1->getDeformedMesh(vars,0,static_cast<int>(vars.size()));
			mout.save(outname.value()+".vtk",3);
		//	cout<<"Beta is "<<gBeta<<endl;
					return 0;
		}
int do_work_Aff(int argc, char* argv[]) 
{ 
	//load base volume

	read_volume(image,inname.value());
	//normalize image intensities
	cout<<"normalize intensity... "<<robmax.value()<<" "<<robmin.value()<<endl;

	image=(image-robmin.value())*255/(robmax.value()-robmin.value());
	
	//cout<<image.getextrapolation()<<endl;
	//image.setpadvalue(0);
	//save_volume(image,"intnorm");
//	if (manMean.value()==-777){
		GmanMean=manMean.value();
//	}else{
//		GmanMean=(manMean.value()-robmin.value())*255/(robmax.value()-robmin.value());
//	}



	
	read_volume(Resimage,inname.value());
	volume<short> brainSeg;
	Resimage=0;
	copyconvert(image,brainSeg);
	brainSeg=0;
	
	const int sizex = image.xsize();
	const int sizey = image.ysize();
	const int sizez = image.zsize();
	xdim=image.xdim();
	ydim=image.ydim();
	zdim=image.zdim();
	
	shapeModel* model1 = new shapeModel;
	//mni152 1mm isotropic
	model1->setImageParameters(182,218, 182, 1, 1, 1);
	//model1->setImageParameters(image.xsize(), image.ysize(), image.zsize(), image.xdim(), image.ydim(), image.zdim());
	cout<<"load model..."<<endl;
		  model1->load_bmv_binaryInfo(modelname.value(),6);
		  model1->load_bmv_binary(modelname.value(),6);	
		  cout<<"register model..."<<endl;
		  //this will set new image parameters
		  model1->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );
		  
		   if (useIntRefModel.value()){
		  shapeModel* modelRef =new shapeModel;
		 	modelRef->setImageParameters(182,218, 182, 1, 1, 1);

			cout<<"use a reference model "<<modelname2.value()<<endl;
			  //this bit of does uses one strcuture as a reference for others.
			  modelRef->load_bmv_binaryInfo(modelname2.value(),1);
			  modelRef->load_bmv_binary(modelname2.value(),1);	
			  cout<<"register model..."<<endl;
			  //this will set new image parameters
			  	  vector<float> varstemp;
			  for (int i=0; i<modelRef->getNumberOfModes();i++){
				  //perturbe the system
				 // cout<<"set tup vars"<<endl;
				  varstemp.push_back(0);
			  }
			  
			  modelRef->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );
			  vector<float> relStd;
			  for (unsigned int i=0;i<varstemp.size();i++){
				  relStd.push_back(stdTrunc.value());
			  }
		
			  //this used for wbir	
			  vector<bool> select;
			  for (unsigned int i=0;i<varstemp.size();i++){
				  if(i<10){
					  select.push_back(true);
				  }else{
					  select.push_back(false);
				  }
			  }		
			  //globshcond enbales one to turn off conditional
			   GlobShCond=false;
		  fitModel(modelRef,modelname2.value(),&varstemp,relStd,0,10,0,select,1, false,0.5, 0.15);
		  GlobShCond=true;
			  Mesh mtemp = modelRef->getDeformedMesh(varstemp,0,static_cast<int>(varstemp.size()));
			  int bounds[6]={0,0,0,0,0,0};
			  volume<short> mask;
			  mask=make_mask_from_meshInOut(image,mtemp,modelRef->getLabel(0),bounds);
			  
			  //calculate intensity histogram
			  float maxint=0,minint=0;
			  vector<float> vintens;
			  modelRef->intensityHistMaxMin(&image,&mask,&mtemp,modelRef->getLabel(0),&vintens, &maxint, &minint);
			  
			  GmanMean=mode(vintens,minint,maxint);
			  cout<<"the mode of the reference distribution is "<<GmanMean<<endl;
			  
		  }
		  

		vector<float> vars;
		for (int i=0; i<model1->getNumberOfModes()+7;i++){
			//perturbe the system
			vars.push_back(0);
		}
	//	addVolumeTo4D(model1,vars);
		searchDim=0.2;
		
		
		
		if (loadvars.value()){
	
		  read_bvars(bvarsname.value(),&vars,model1->getNumberOfModes());
		
		}
		
		if (useCondMode.value()){
		
			shapeModel* model2 = new shapeModel;
			//mni152 1mm isotropic
			model2->setImageParameters(182,218, 182, 1, 1, 1);
			cout<<"load model..."<<endl;
			
			//assumes appearance models
			model2->load_bmv_binaryInfo(modelname2.value(),1);
			model2->load_bmv_binary(modelname2.value(),1);
			cout<<"register model..."<<endl;
			//this will set new image parameters
			model2->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );

			vector<float> varstemp;
			read_bvars(bvarsname.value(),&varstemp,model2->getNumberOfModes());
			
			Mesh mtemp = model1->getDeformedMesh(varstemp,0,static_cast<int>(varstemp.size()));
				int bounds[6]={0,0,0,0,0,0};
			volume<short> mask;
				mask=make_mask_from_meshInOut(image,mtemp,model2->getLabel(0),bounds);
					
				//calculate intensity histogram
				float maxint=0,minint=0;
				vector<float> vintens;
				model2->intensityHistMaxMin(&image,&mask,&mtemp,model2->getLabel(0),&vintens, &maxint, &minint);
			
				float modePri=mode(vintens,minint,maxint);
				
				//evaluate conditional mode (assuem gaussian)
				//for now implement hipp given put
				condMode=149.2094-0.7638*(modePri-184.5554);
			GmanMean=condMode;
			cout<<"conditional mode is "<<GmanMean<<endl;
		}
		
	    
		if (shcond.value()){
		//	cout<<"shcond.value"<<endl;
		int M;
				Matrix mBx2;
				Matrix mBcx1;
		//		cout<<"readbmap "<<bmapname.value()<<endl;
				M=readBmap(bmapname.value(),vars,&mBx2,&mBcx1);
		//		cout<<"readbmap done"<<endl;
			//	int nummodes=
				ColumnVector Bx2(model1->getNumberOfModes());
				//load bx2 vars that were read
				for (unsigned int i=0;i<vars.size();i++){
					Bx2.element(i)=vars.at(i);
				}
				mBx2map=mBx2*Bx2;
				//bmap now directly reflects the transformation t0oo bc
				//mBx1inv=mBcx1;
				mBx1inv=mBcx1.i();
				
				vector<float> v_cmean;
				v_cmean=bTransform(&vars,mBx2,M);
				cout<<"btransfomred"<<endl;
		
				for (int i=0;i<static_cast<int>(vars.size());i++){
				  if (i<g.value()){
					vars.at(i)=v_cmean.at(i);
					cout<<vars.at(i)<<" ";
					}else{
					vars.at(i)=0;
					}
				}
				cout<<endl;
					//vars=v_cmean;
				cout<<"done loading bmap stuff"<<endl;
				if (shcond2.value()){
					vector<float> vars2;
					vars2=vars;
					read_bvars(bvarsname2.value(),&vars,model1->getNumberOfModes());
					cout<<"readbmap2 "<<bmapname2.value()<<endl;
					M=readBmap2(bmapname2.value(),vars,&mBx2,&mBcx1);
					cout<<"readbmap2 done"<<endl;
					//	int nummodes=
					ColumnVector Bx2(model1->getNumberOfModes());
					//load bx2 vars that were read
					for (unsigned int i=0;i<vars.size();i++){
						Bx2.element(i)=vars.at(i);
					}
					mBx2map2=mBx2*Bx2;
					mBx1inv2=mBcx1.i();
					
					vector<float> v_cmean;
					v_cmean=bTransform(&vars,mBx2,M);
					cout<<"btransfomred"<<endl;
					for (int i=0;i<static_cast<int>(vars.size());i++){
						if (i<-1){
							vars.at(i)+=v_cmean.at(i)+vars2.at(i);
						}else{
							vars.at(i)=0;
						}
					}
					//vars=v_cmean;
					cout<<"done loading bmap stuff"<<endl;
					
				}	
				
					
					
		}
	//addVolumeTo4D(model1,vars);
		vector<float> relStd;
		for (unsigned int i=0;i<vars.size();i++){
			relStd.push_back(stdTrunc.value());
		}
		
		string modelName;
		modelName=modelname.value();
			int lb=0,ub=0;
			
			
			//this used for wbir	
			vector<bool> select;
			for (unsigned int i=0;i<vars.size();i++){
				select.push_back(false);
			}		
			
			for (int e=0;e<f.value();e++){
			 
				  for (int shift=0;shift<c.value()*static_cast<int>(g.value());shift+=static_cast<int>(g.value())){
				   				
					lb=0;
					ub=shift+g.value()-1;
					for (int i=0;i<static_cast<int>(vars.size());i++){
						if ((i>=lb)&&(i<=ub)){
							select.at(i)=true;
						}else{
							select.at(i)=false;
						}
					}
				
						costmodeG=0;
						if (baam.value()){
							//use appearance model
					//		fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,1, false,0.15, 0.1);
						fitModelAff(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,1, false,0.1, 0.075);
						//	if (!nograd.value()){
						//	fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,0, false,0.15, 0.1);
							//fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,1, false,0.075, 0.05);
						//	}
						}else if (overlap.value()){
							//uses dice overlap---> work for label data
							fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,2, false,0.15, 0.1);	

						}else{ 
							//use asm 
							fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,0, false,0.15, 0.1);	
 
						}						
				}
			}
		
 			//save_volume4D(fit,"progression");
		//	addVolumeTo4D(model1,vars);
			save_volume(d4tmp,outname.value());
			Mesh mout = model1->getDeformedMeshAff7(vars,0,static_cast<int>(vars.size()));
			int bounds[6]={0,0,0,0,0,0};
			d4tmp=make_mask_from_meshInOut(image,mout, model1->getLabel(0),bounds);
			save_volume(d4tmp,outname.value());

			mout.save(outname.value()+".vtk",3);
		//	cout<<"Beta is "<<gBeta<<endl;
					return 0;
		}
		
		

int do_work_joint(int argc, char* argv[]) 
{ 
	//load base volume
	GmanMean=manMean.value();
	read_volume(image,inname.value());
	//normalize image intensities
	cout<<"normalize intensity..."<<endl;

	image=(image-robmin.value())*255/(robmax.value()-robmin.value());
	//save_volume(image,"intnorm");
	
	read_volume(Resimage,inname.value());
	volume<short> brainSeg;
	Resimage=0;
	copyconvert(image,brainSeg);
	brainSeg=0;
	
	const int sizex = image.xsize();
	const int sizey = image.ysize();
	const int sizez = image.zsize();
	xdim=image.xdim();
	ydim=image.ydim();
	zdim=image.zdim();
		
	shapeModel* model1 = new shapeModel;
		//mni152 1mm isotropic
	model1->setImageParameters(182,218, 182, 1, 1, 1);
	//model1->setImageParameters(image.xsize(), image.ysize(), image.zsize(), image.xdim(), image.ydim(), image.zdim());
	cout<<"load model..."<<endl;

	//assumes appearance models
	model1->load_bmv_binaryInfo(modelname.value(),1);
	model1->load_bmv_binary(modelname.value(),1);
		cout<<"register model..."<<endl;
	//this will set new image parameters
	model1->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );


	shapeModel* model2 = new shapeModel;
		//mni152 1mm isotropic
	model2->setImageParameters(182,218, 182, 1, 1, 1);
	//model1->setImageParameters(image.xsize(), image.ysize(), image.zsize(), image.xdim(), image.ydim(), image.zdim());
	cout<<"load model..."<<endl;

	//assumes appearance models
	model2->load_bmv_binaryInfo(modelname2.value(),1);
	model2->load_bmv_binary(modelname2.value(),1);
		cout<<"register model..."<<endl;
	//this will set new image parameters
	model2->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );

	vector< shapeModel* > Vmodels;
	Vmodels.push_back(model1);
	Vmodels.push_back(model2);

		vector<float> vars;
		vector<float> vars2;
		for (int i=0; i<model1->getNumberOfModes();i++){
			//perturbe the system
			vars.push_back(0.00);
			vars2.push_back(0.00);
		}
//if you want to load previous shape
	//	read_bvars(bvarsname.value(),&vars,model1->getNumberOfModes());
//		  read_bvars(bvarsname2.value(),&vars2,model1->getNumberOfModes());
		  
		vector< vector<float> > Vvars;
		Vvars.push_back(vars); 
		Vvars.push_back(vars2);
		 
	//	addVolumeTo4D(model1,vars);
		searchDim=0.2;
		
	///want to load multiple conditionals from a file 
				//load a conditional matrix
//		ifstream fbmaps;
//		fbmaps.open(bmapname.value());
//		while (!fbmaps.eof()){
//				string fnametemp;
//				fbmaps>>fnametemp;
				int M;
				Matrix mBx2;
				Matrix mBcx1;
			//	M=readBmap(fnametemp,vars,&mBx2,&mBcx1);
			M=readBmap(bmapname.value(),vars,&mBx2,&mBcx1);
				//these are the global matrices used
				//bmap now directly reflects the transformation t0oo bc
				//mBx1inv=mBcx1;
				mBx1inv=mBcx1.i();
//				VmBx1inv.push_back(mBx1inv);
//				VmBx2.push_back()
//		}
		//addVolumeTo4D(model1,vars);
		vector<float> relStd;
		for (unsigned int i=0;i<vars.size();i++){
			relStd.push_back(stdTrunc.value());
		}
		
		string modelName;
		modelName=modelname.value();
			int lb=0,ub=0;
			
			
			//this used for wbir	
			vector<bool> select;
			for (unsigned int i=0;i<vars.size();i++){
				select.push_back(false);
			}		
			
			for (int e=0;e<f.value();e++){
			 
				  for (int shift=0;shift<c.value()*static_cast<int>(g.value());shift+=static_cast<int>(g.value())){
				   				
					lb=shift;
					ub=shift+g.value()-1;
					for (int i=0;i<static_cast<int>(vars.size());i++){
						if ((i>=lb)&&(i<=ub)){
							select.at(i)=true;
						}else{
							select.at(i)=false;
						}
					}
				
						costmodeG=0;
					
						fitModelDual(Vmodels,&Vvars,relStd,lb,ub-lb+1,select, 0.5, 0.15);
										
				}
			}
		
 			//save_volume4D(fit,"progression");
			cout<<"final vars"<<endl;
			for (unsigned int i=0;i<Vvars.at(0).size();i++){
			 cout<<Vvars.at(0).at(i)<<" ";
			
			}
			cout<<endl;
			addVolumeTo4D(Vmodels.at(0),Vvars.at(0));
			save_volume(d4tmp,outname.value()+"A");
			addVolumeTo4D(Vmodels.at(1),Vvars.at(1));
				for (unsigned int i=0;i<Vvars.at(0).size();i++){
			 cout<<Vvars.at(1).at(i)<<" ";
			
			}
			save_volume4D(fit,"progression");
			//save_volume(d4tmp,outname.value());
			
			save_volume(d4tmp,outname.value()+"B");
			Mesh mout = model1->getDeformedMesh(vars,0,static_cast<int>(vars.size()));
			mout.save(outname.value()+".vtk",3);
		//	cout<<"Beta is "<<gBeta<<endl;
					return 0;
		}

int do_work_jointN(int argc, char* argv[]) 
{ 
	//load base volume
	GmanMean=manMean.value();
	read_volume(image,inname.value());
	//normalize image intensities
	cout<<"normalize intensity..."<<endl;

	image=(image-robmin.value())*255/(robmax.value()-robmin.value());
	//save_volume(image,"intnorm");
	
	read_volume(Resimage,inname.value());
	volume<short> brainSeg;
	Resimage=0;
	copyconvert(image,brainSeg);
	brainSeg=0;
	
	const int sizex = image.xsize();
	const int sizey = image.ysize();
	const int sizez = image.zsize();
	xdim=image.xdim();
	ydim=image.ydim();
	zdim=image.zdim();
		
		
	vector<shapeModel* > vModels;
	vector< vector<float> > vVars;
	ifstream fmodel;
	fmodel.open(modelname.value().c_str());
	
	
	vector<string> modelNames;
	while (!fmodel.eof()){
		cout<<"reading models "<<endl;
		string mnametemp;
		fmodel>>mnametemp;
		cout<<"manem "<<mnametemp.length()<<" "<<mnametemp<<endl;
		if (mnametemp.length()>0){
			modelNames.push_back(mnametemp);
						shapeModel* model1 = new shapeModel;
			//mni152 1mm isotropic
			model1->setImageParameters(182,218, 182, 1, 1, 1);
			//assumes appearance models
			model1->load_bmv_binaryInfo(mnametemp,1);
			model1->load_bmv_binary(mnametemp,1);
			model1->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );
			vModels.push_back(model1);
			cout<<"add to vector"<<endl;
		}
			
		
			
	}
	
	cout<<"Register the model"<<endl;
	//also setup bvars
	for (unsigned int i=0; i<vModels.size();i++){
		cout<<"i "<<i<<vModels.at(i)<<endl;
		vector<float> vars;
		for (int j=0; j<vModels.at(i)->getNumberOfModes();j++){
			//perturbe the system
		//	if (j<g.value()){
		//	vars.push_back(0.2);
		//	}else{
			vars.push_back(0);
		//	}
		}
		vVars.push_back(vars);
	}
	//this will set new image parameters
	cout<<"model is registered "<<endl;
		
				searchDim=0.2;
		
	///want to load multiple conditionals from a file 
	vector< Matrix > vmBx2;
	vector< Matrix > vmBx1inv;
	vector< int > vKx2;
	ifstream fbmap;
	fbmap.open(bmapname.value().c_str());
	for (unsigned int i=0; i<(vModels.size()-1);i++){
		string btemp;
		fbmap>>btemp;
	//	int M=vModels.at(i)->getNumberOfSubjects();
		Matrix mBx2;
		Matrix mBcx1;
		
		//scale has to be done on the fly
			cout<<"read btemp "<<btemp<<" "<<endl;
		int kx2=readBmapNoScale(btemp,vVars.at(i),&mBx2,&mBcx1);
		
		cout<<"read btemp "<<kx2<<endl;
		vKx2.push_back(kx2);
		vmBx2.push_back(mBx2);
		
		vmBx1inv.push_back(mBcx1.i());
		Matrix bmptest=mBcx1.i();
		for (int i=0;i<317;i++){
		 cout<<bmptest.element(0,i)<<" ";
		}
		
		cout<<endl;
		for (int i=0;i<317;i++){
		 cout<<mBx2.element(0,i)<<" ";
		}
	}
	cout<<"Mapping matrices read"<<endl;
	//addVolumeTo4D(model1,vars);
	vector<float> relStd;
	for (unsigned int i=0;i<vVars.at(0).size();i++){
		relStd.push_back(stdTrunc.value());
	}
	
	int lb=0,ub=0;
	//this used for wbir	
	vector<bool> select;
	for (unsigned int i=0;i<vVars.at(0).size();i++){
		select.push_back(false);
	}		
	
	for (int e=0;e<f.value();e++){
		
		for (int shift=0;shift<c.value()*static_cast<int>(g.value());shift+=static_cast<int>(g.value())){
			
			lb=shift;
			ub=shift+g.value()-1;
			for (int i=0;i<static_cast<int>(vVars.at(0).size());i++){
				if ((i>=lb)&&(i<=ub)){
					select.at(i)=true;
				}else{
					select.at(i)=false;
				}
			}
			
			costmodeG=0;
			cout<<"You are fitting models: "<<endl;
			for (unsigned int i=0; i<vModels.size();i++){
				cout<<vModels.at(i)->getLabel(0)<<" ";
			}
			cout<<endl;
			fitModelN(vModels,&vVars,vmBx2,vmBx1inv, vKx2,relStd,lb,ub-lb+1,select, 0.5, 0.15, modelNames);
			
		}
	}
	
 			//save_volume4D(fit,"progression");
	cout<<"final vars"<<endl;
	for (unsigned int i=0;i<vVars.at(0).size();i++){
		cout<<vVars.at(0).at(i)<<" ";
		
	}
	cout<<endl;
	addVolumeTo4D(vModels.at(0),vVars.at(0));
	save_volume(d4tmp,outname.value()+"A");
	addVolumeTo4D(vModels.at(1),vVars.at(1));
				for (unsigned int i=0;i<vVars.at(0).size();i++){
					cout<<vVars.at(1).at(i)<<" ";
					
				}
				save_volume4D(fit,"progression");
				//save_volume(d4tmp,outname.value());
				
				save_volume(d4tmp,outname.value()+"B");
		//		Mesh mout = vMmodel1->getDeformedMesh(vars,0,static_cast<int>(vars.size()));
		//		mout.save(outname.value()+".vtk",3);
				//	cout<<"Beta is "<<gBeta<<endl;
				return 0;
}
int do_work_jointN2(int argc, char* argv[]) 
{ 
	//load base volume
	GmanMean=manMean.value();
	read_volume(image,inname.value());
	//normalize image intensities
	cout<<"normalize intensity..."<<endl;

	image=(image-robmin.value())*255/(robmax.value()-robmin.value());
	//save_volume(image,"intnorm");
	
	read_volume(Resimage,inname.value());
	volume<short> brainSeg;
	Resimage=0;
	copyconvert(image,brainSeg);
	brainSeg=0;
	
	const int sizex = image.xsize();
	const int sizey = image.ysize();
	const int sizez = image.zsize();
	xdim=image.xdim();
	ydim=image.ydim();
	zdim=image.zdim();
		
		
	vector<shapeModel* > vModels;
	vector< vector<float> > vVars;
	ifstream fmodel;
	fmodel.open(modelname.value().c_str());
	
	
	vector<string> modelNames;
	while (!fmodel.eof()){
		cout<<"reading models "<<endl;
		string mnametemp;
		fmodel>>mnametemp;
		cout<<"manem "<<mnametemp.length()<<" "<<mnametemp<<endl;
		if (mnametemp.length()>0){
			modelNames.push_back(mnametemp);
						shapeModel* model1 = new shapeModel;
			//mni152 1mm isotropic
			model1->setImageParameters(182,218, 182, 1, 1, 1);
			//assumes appearance models
			model1->load_bmv_binaryInfo(mnametemp,1);
			model1->load_bmv_binary(mnametemp,1);
			model1->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );
			vModels.push_back(model1);
			cout<<"add to vector"<<endl;
		}
			
		
			
	}
	
	cout<<"Register the model"<<endl;
	//also setup bvars
	for (unsigned int i=0; i<vModels.size();i++){
		cout<<"i "<<i<<vModels.at(i)<<endl;
		vector<float> vars;
		for (int j=0; j<vModels.at(i)->getNumberOfModes();j++){
			//perturbe the system
		//	if(j<g.value()){
		//	vars.push_back(0.1);
		//	}else{
			vars.push_back(0);
		//	}
		}
		vVars.push_back(vars);
	}
	
	
	
	
	ifstream fbvars;
	fbvars.open(bvarsname.value().c_str());
	if (loadvars.value()){
		for (unsigned int i=0; i<(vModels.size());i++){
			string btemp,stemp;
			fbvars>>btemp;
			getline(fbvars,stemp);
			cout<<"read btemp/stemp bvars"<<btemp<<" "<<stemp<<endl;
			if (strcmp(btemp.c_str(),"none")){
				cout<<"read bvars "<<btemp<<i<<endl;
				read_bvars(btemp,&vVars.at(i),vModels.at(i)->getNumberOfModes());

			}else{
				cout<<"read bvars none "<<btemp<<i<<endl;
			}
			
		}
	}
	
	
	//this will set new image parameters
	cout<<"model is registered "<<endl;
		
				searchDim=0.2;
		
	///want to load multiple conditionals from a file 
	vector< vector< Matrix > > VvmBx2;
	vector<Matrix> vmBcx1inv;
	vector< int > vKx2;
	vector< vector<int> > vPredLabels;
	ifstream fbmap;
	fbmap.open(bmapname.value().c_str());
	for (unsigned int i=0; i<(vModels.size()-1);i++){
		string btemp,stemp;
		fbmap>>btemp;
		getline(fbmap,stemp);
		//	int M=vModels.at(i)->getNumberOfSubjects();
		Matrix mBcx1;
		vector<Matrix> vmBx2;
		
		//scale has to be done on the fly
		cout<<"read btemp/stemp "<<btemp<<" "<<stemp<<endl;
		vector<int> PredLabel;
		int kx2;
		if (strcmp(btemp.c_str(),"none")){
			cout<<"read bmap "<<btemp<<i<<endl;
			kx2=readMultiBmapNoScale(btemp,vVars.at(i),&vmBx2,&mBcx1,&PredLabel);
		}else{
			kx2=0;
			cout<<"read bmap none "<<btemp<<i<<endl;
		}
		
		
		cout<<"pred dim "<<kx2<<endl;
		vKx2.push_back(kx2);
		//dont need to iinvert differen t implementation
		vmBcx1inv.push_back(mBcx1);
//		for (int i=0;i<317;i++){
//		 cout<<mBcx1.element(0,i)<<" ";
//		}
//		cout<<endl;
//		for (int i=0;i<317;i++){
//		 cout<<vmBx2.at(0).element(0,i)<<" ";
//		}
//		cout<<endl;
		VvmBx2.push_back(vmBx2);
		vPredLabels.push_back(PredLabel);
	}
	
	cout<<"check predictive labels"<<endl;
	for (int i=0;i<vPredLabels.size();i++){
		for (int j=0;j<vPredLabels.at(i).size();j++){
		cout<<i<<" "<<j<<" "<<vPredLabels.at(i).at(j)<<endl;
		}
		
		}
	
	
	cout<<"Mapping matrices read"<<endl;
	//addVolumeTo4D(model1,vars);
	vector<float> relStd;
	for (unsigned int i=0;i<vVars.at(0).size();i++){
		relStd.push_back(stdTrunc.value());
	}
	
	int lb=0,ub=0;
	//this used for wbir	
	vector<bool> select;
	for (unsigned int i=0;i<vVars.at(0).size();i++){
		select.push_back(false);
	}		
	
	
	//this tells which models to initialize , should be a numeric file of 0 adn 1
	
// 	ifstream finit;
// 	finit.open(initname.value().c_str());
// 	if (singleInit.value()){
// 	  cout<<"Initialize uisng single structure"<<endl;
		
// 	  string stemp;
			
// 	  //fbvars>>stemp;
// 	  for (unsigned int mod=0; mod<(vModels.size());mod++){
// 	    int Itemp;
		
// 	    finit>>Itemp;
// 	    getline(finit,stemp);
// 	    cout<<"read btemp/stemp bvars "<<Itemp<<" "<<stemp<<endl;
// 	    if (Itemp==1){
// 	      for (int i=0;i<static_cast<int>(vVars.at(mod).size());i++){
// 		if ((i>=0)&&(i<=g.value())){
// 		  select.at(i)=true;
// 		}else{
// 		  select.at(i)=false;
// 		}
// 	      }
// 	      GlobShCond=false;
// 	      fitModel(vModels.at(mod),modelNames.at(mod),&(vVars.at(mod)),relStd,0,g.value(),0,select,1, false,0.5, 0.15);
// 	      //reset the select
// 	      for (int i=0;i<static_cast<int>(vVars.at(mod).size());i++){
// 		select.at(i)=false;
// 	      }
// 	      //reset global found mode
// 	      globfoundmode=false;
// 	    }
// 	    GmanMean=manMean.value();
// 	  }

	
// 	}
	//store which have been initialized

	vector< bool > sinit;
         for (unsigned int mod=0; mod<(vModels.size());mod++){
	   sinit.push_back(false);
	 }

	 //intialize using single structure fit
	ifstream finit;
	finit.open(initname.value().c_str());
	if (singleInit.value()){
	  cout<<"Initialize uisng single structure"<<endl;
		
	  string stemp;
			
	  //fbvars>>stemp;
	  for (unsigned int mod=0; mod<(vModels.size());mod++){
	    int Itemp;
		
	    finit>>Itemp;
	    getline(finit,stemp);
	    cout<<"read btemp/stemp bvars "<<Itemp<<" "<<stemp<<endl;
	    if (Itemp==1){
	      sinit.at(mod)=true;
	      for (int i=0;i<static_cast<int>(vVars.at(mod).size());i++){
		if ((i>=0)&&(i<=g.value())){
		  select.at(i)=true;
		}else{
		  select.at(i)=false;
		}
	      }
	      GlobShCond=false;
	      fitModel(vModels.at(mod),modelNames.at(mod),&(vVars.at(mod)),relStd,0,g.value(),0,select,1, false,0.5, 0.15);
	      //reset the select
	      for (int i=0;i<static_cast<int>(vVars.at(mod).size());i++){
		select.at(i)=false;
	      }
	      //reset global found mode
	      globfoundmode=false;
	    }
	    GmanMean=manMean.value();
	  }

	
	}
	int count=10000;
	while (count<(sinit.size()-1)){
	  cout<<"count "<<count<<endl;
	//run through and intialize neghiibours, this assumes single dependencies
	for (unsigned int mod=0; mod<(vModels.size()-1);mod++){
	  //doesn't initialize last in chain
	  if ((!sinit.at(mod))&&(sinit.at(mod+1))){
	 
	      ColumnVector Bx2(vModels.at(mod+1)->getNumberOfModes());
	      //load bx2 vars that were read
	      for (unsigned int i=0;i<vVars.at(mod).size();i++){
		Bx2.element(i)=vVars.at(mod+1).at(i);
	      }
	      ColumnVector cMean(vVars.at(mod).size());
	      cMean<<VvmBx2.at(mod).at(0)*Bx2;
						
	      cout<<"initialize model "<<mod<<endl;
		for (int i=0;i<g.value()-1;i++){
		  vVars.at(mod).at(i)=cMean.element(i);
		  cout<<cMean.element(i)<<endl;
		}
		sinit.at(mod)=true;
		
	  }
	 
	

	}
	  count=0;
	
	  for (unsigned int i=0;i<(sinit.size()-1);i++){
	    if (sinit.at(i)){
	      count++;
	    } 
	  }
	
	}
	
	
	
	for (int e=0;e<f.value();e++){
		
		for (int shift=0;shift<c.value()*static_cast<int>(g.value());shift+=static_cast<int>(g.value())){
			
			lb=shift;
			ub=shift+g.value()-1;
			for (int i=0;i<static_cast<int>(vVars.at(0).size());i++){
				if ((i>=lb)&&(i<=ub)){
					select.at(i)=true;
				}else{
					select.at(i)=false;
				}
			}
			
			costmodeG=0;
			cout<<"You are fitting models: "<<endl;
			for (unsigned int i=0; i<vModels.size();i++){
				cout<<vModels.at(i)->getLabel(0)<<" ";
			}
			cout<<endl;
			//was nomrally thresh at 0.1
			fitModelN2(vModels,&vVars,VvmBx2,vmBcx1inv,vPredLabels, vKx2,relStd,lb,ub-lb+1,select, 0.5, 0.3, modelNames);
			
		}
	}
	
 			//save_volume4D(fit,"progression");
	cout<<"final vars"<<endl;
	
	for (unsigned int model=0; model<vModels.size();model++){
		ofstream fout;
		stringstream modelNum;
		modelNum<<model;
		string numstr;
		modelNum>>numstr;
		string name=outname.value()+numstr;
		addVolumeTo4D(vModels.at(model),vVars.at(model));
		save_volume(d4tmp,name);
		for (unsigned int i=0;i<vVars.at(model).size();i++){
			cout<<vVars.at(model).at(i)<<" ";
			
		}
	}
	
	

			
		//		save_volume4D(fit,"progression");
				//save_volume(d4tmp,outname.value());
				
		//		save_volume(d4tmp,outname.value()+"B");
		//		Mesh mout = vMmodel1->getDeformedMesh(vars,0,static_cast<int>(vars.size()));
		//		mout.save(outname.value()+".vtk",3);
				//	cout<<"Beta is "<<gBeta<<endl;
				return 0;
}


void testbmap(int argc, char* argv[]) 
{
	Matrix mBx2;
	Matrix mBcx1;
	vector<float> vars;
	vars.push_back(0);
	readBmap(inname.value(),vars,&mBx2,&mBcx1);
}
int main(int argc,char *argv[])
{
	
	Tracer tr("main");
	OptionParser options(title, examples);
	
	try {
		// must include all wanted options here (the order determines how
		//  the help message is printed)
		options.add(inname);
				options.add(initname);

		options.add(singleInit);

		options.add(baam);
		options.add(loadvars);
		options.add(shcond);
		options.add(joint);
		options.add(jointN);
				options.add(jointN2);

		options.add(shcond2);
	//	options.add(nograd);
		
		options.add(bmapname);
			options.add(bmapname2);


		options.add(overlap);
		options.add(outname);
		options.add(stdTrunc);
		options.add(inputprob);
			options.add(verbose);
			options.add(help);
			options.add(EMiter);
		options.add(modelname); 
		options.add(modelname2); 
		options.add(reestimate); 
		options.add(useIntRefModel);
		options.add(runAff);
		options.add(modeprior); 
				options.add(useCondMode); 
	options.add(useConj3); 
options.add(histmode); 
		options.add(bvarsname);
		options.add(bvarsname2);
		options.add(flirtmatname);
		options.add(manMean); 
		options.add(c); 
		options.add(gmmLesser); 
		options.add(f);
		options.add(g);
		options.add(robmin);
		options.add(robmax);
				
		nonoptarg = options.parse_command_line(argc, argv);
		
		// line below stops the program if the help was requested or 
		//  a compulsory option was not set
				if (  (!options.check_compulsory_arguments(true) ))
		{
			options.usage();
			exit(EXIT_FAILURE);
		}
		
		// Call the local functions
		if (runAff.value()){
		do_work_Aff(argc,argv);
		 } else if (joint.value()){
		do_work_joint(argc,argv);
		
		}else if (jointN.value()){
		
		do_work_jointN(argc,argv);
		}else if (jointN2.value()){
		
		do_work_jointN2(argc,argv);
		}else{
			do_work(argc,argv);
		}
		
	}  catch(X_OptionError& e) {
		options.usage();
		cerr << endl << e.what() << endl;
		exit(EXIT_FAILURE);
	} catch(std::exception &e) {
		cerr << e.what() << endl;
	} 
	
	return 0;// do_work(argc,argv);
}

