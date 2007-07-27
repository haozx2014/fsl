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

Option<bool> useIntRefModel(string("--useIntRefModel"), false,
					   string("reeistmate mode for each iteration"),
					   false, no_argument);

Option<bool> baam(string("--baam"), false,
					   string("use appearaance model"),
					   false, no_argument);

Option<bool> overlap(string("--overlap"), false,
					   string("use overlapcost - can fit mesh to labelled data"),
					   false, no_argument);

					   
Option<string> inname(string("-i,--in"), string(""),
					  string("filename of input image to be segmented"),
					  true, requires_argument);

Option<string> flirtmatname(string("-l,--flirtMatrix"), string(""),
						 string("filename containing flirt matrix (transformattion to MNI space"),
						 true, requires_argument);
Option<string> modelname(string("-m,--modelin"), string(""),
						 string("filename of input model: the structure to be segmented"),
						 true, requires_argument);
Option<string> modelname2(string("-n,--modelin2"), string(""),
						 string("filename of 2nd input model for dual: the structure to be segmented"),
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
Option<float> stdTrunc(string("-s,--sh"), 8.0,
					   string("number of standard deviation to truncate at"),
					   false, requires_argument);		
Option<float> manMean(string("-a,--mean intensity"), -777.0,
					  string("override mean estimation"),
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
//volume<float> Resimage;
//int costmodeG=0;
float xdim, ydim, zdim;
//float searchDim;
float searchRes=0.2;
//volume4D<short> fit;
//volume<short> d4tmp;
bool gradzero=false;
//const float pi = 3.14159;
//unsigned int M;
//vector<float> v_ceigs;
ColumnVector mBx2map, mBx2map2;
Matrix mBx1inv,mBx1inv2;
//vector< Matrix > VmBx1inv, VmBx2;
//bool GlobShCond=true;
float kpred,kpred2;
float globMean=0;
//vector<float> VglobMean;
//int gbounds[6]={1000,0,1000,0,1000,0};
//float gBeta;
float GmanMean;
//float condMode;
bool globfoundmode=false;
//vector<bool> Vglobfoundmode;
//bool globReest=false;
//float globModePriorMean=-777;
//	float	globModePriorVar=-777;
//	Matrix bcross;
//	Matrix bmapSh;
//	Matrix bmapAff;


int readBmap(string fname, vector<float> b2, Matrix *matBx2, Matrix * matBcx1){
	int M; //number of subjects
	ifstream fin;
	fin.open(fname.c_str());
	string stemp;
	getline(fin,stemp);
	fin>>stemp;
	fin>>M;
	matBx2->ReSize(M,M);
	matBcx1->ReSize(M,M);
	getline(fin,stemp); 
	getline(fin,stemp);
	double ftemp;

	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBx2->element(i,j)=ftemp;
		}
	}
	
	getline(fin,stemp);
	getline(fin,stemp);

	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBcx1->element(i,j)=ftemp;

		}
	}
	getline(fin,stemp);
	getline(fin,stemp);
	float kx2;
	fin>>kx2;
	kpred=kx2;
	
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
	vector<float> v_ceigs;
	
	for (int i=0;i<M;i++){
	  fin>>ftemp;
	  v_ceigs.push_back(ftemp);
	}
	
	return M;
	
	
	
}

int readBmap2(string fname, vector<float> b2, Matrix *matBx2, Matrix * matBcx1){
//only diffence is the assignment of krped2 versus kpred, anitquate later
	int M; //number of subjects
	
	ifstream fin;

	fin.open(fname.c_str());

	string stemp;
	
	getline(fin,stemp);
	
	fin>>stemp;
	fin>>M;
	
	matBx2->ReSize(M,M);
	matBcx1->ReSize(M,M);
	getline(fin,stemp); 
	getline(fin,stemp);
	
	//read first matrix 

	double ftemp;

	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBx2->element(i,j)=ftemp;

		}
	}
	
	

	getline(fin,stemp);
	getline(fin,stemp);

	for (int i=0;i<M;i++){
		for (int j=0;j<M;j++){
			fin>>ftemp;
			matBcx1->element(i,j)=ftemp;
		}
	}
	getline(fin,stemp);
	getline(fin,stemp);
	float kx2;
	fin>>kx2;
	kpred2=kx2;
	

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
	
	vector<float> v_ceigs;
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
	double mininc = min(xdim,min(ydim,zdim)) * .25;


	
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

  	double xdim = (double) image.xdim();
	double ydim = (double) image.ydim();
	double zdim = (double) image.zdim();

	//in new version of bet2
	double mininc = min(xdim,min(ydim,zdim)) * 0.5;


	volume<short> res = image;
	for (list<Triangle*>::const_iterator i = m._triangles.begin(); i!=m._triangles.end(); i++)
    {
		Vec n = (*(*i)->get_vertice(0) - *(*i)->get_vertice(1));
		double d = n.norm();
		n.normalize();
	       
		
		for (double j=0; j<=d ;  j+=mininc)
		{
			Pt p = (*i)->get_vertice(1)->get_coord()  + (double)j* n;
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



float mode(vector<float> vdists){
	int N=static_cast<int>(vdists.size());
	int maxcount=0;
	float maxlowint=0;
	float bins=256;
	bins=128;
	float binwidth=(vdists.at(N-1)-0)/bins;

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
			
			count=0;
			i--;
			bincount++;
			lowint=bincount*binwidth;
		}
	
	}
	
	
	return (maxlowint+binwidth/2.0);
}
float mode(vector<float> vdists, float min, float max){
	int N=static_cast<int>(vdists.size());
	int bins=128;

	float binwidth=(max-min)/static_cast<float>(bins);
	vector<int> bincounts;
	//innitialize bincounts to zero
	for (int b=0;b<bins;b++){
			bincounts.push_back(0);
		}

	for (int i=0;i<N;i++){
		//search thgrough each bin
		for (int b=0;b<bins;b++){
			
			if (vdists.at(i)<(min+(b+1)*binwidth)){
				
				bincounts.at(b)++;
				break;
			}
		}
	}
	
	//search for max bin count
	int maxcount=0;
	int maxind=0;
	for (int b=0;b<bins;b++){
		
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
			
			count=0;
			i--;
			bincount++;
			lowint=bincount*binwidth;
		}
		
	}
	
	
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
	
	return (maxlowint+binwidth/2.0);
}	
	//meidan cost func
float costfuncOverlap(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode){
	
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
		Mesh m = model1->getDeformedMesh(vars,i,static_cast<int>(vars.size()));
		mask=make_mask_from_mesh(image,m,model1->getLabel(i),bounds);
		
		cost+=model1->volumeDistance(&mask,&image,bounds,&m);
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

	return cost;
}
float costfuncApp(shapeModel* model1, vector<float> vars, vector<float> fvals, int ex, int costmode, float premean){
 	volume<short> mask;
	volume<short> maskotl;
	double cost=0;
	int bounds[6]={0,0,0,0,0,0};
	vector<float> iprof;
	//new cost
	vector<float> dif;
	//new cost
	for (int i=0;i<model1->getNumberOfShapes()-ex;i++){
		float mean=premean;
		Mesh m = model1->getDeformedMesh(vars,i,static_cast<int>(vars.size()));
		if ((premean==-333)&&(GmanMean==-777)&&((!globfoundmode))){
			mask=make_mask_from_meshInOut(image,m,model1->getLabel(i),bounds);
			//calculate intensity histogram
			float maxint=0,minint=0;
			vector<float> vintens;
			model1->intensityHistMaxMin(&image,&mask,&m,model1->getLabel(i),&vintens, &maxint, &minint);
			if (vintens.size()<=1){
			  cout<<"WARNING: NO INTERIOR VOXELS TO ESTIMATE MODE"<<endl;
			}
			mean=mode(vintens,minint,maxint);
			if ((globMean==mean)&&(vintens.size()>=1)){
				if (verbose.value()){
				  cout<<"Found mode "<<mode(vintens,minint,maxint)<<" "<<mean<<endl;
				}
				globfoundmode=true;
			}
			globMean=mean;
		}else if (GmanMean!=-777){
			//can manually set or use in intref
			globMean=GmanMean;
		}
		iprof=model1->getDeformedIprof(vars,i,vars.size());
		int ipp=model1->getIPP(0);
		
		//***************sample the image in this loop*********************//
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
			}
			
			count++;
		}
	}//end of shape interation, all dif in dif
	
	//************Calculate conditional I | s *************************//
	
	
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
	
	//P(I) is exlcuded use proportionailty....doesn't work as well in practice
	
	//the multiplication by n-1 or n is left out becomes constant in log cost
	//probI*=M;
	
	
	/////%%%%%%%%%%%%%%%%%%%%%%%%% Now claculate costfunction
	
	//calculate cumulcative number of points across all meshes
	int cumnum=0;
	for (int q=0;q<model1->getNumberOfShapes();q++){
		cumnum+=model1->getNumberOfPoints(q);
	}
	
	//Define constant
	double alp=M-1.0/M;
	double binner=0;
	double gammav=alp/(alp-2);
	float k2=3.0*static_cast<float>(cumnum);
	float k1=static_cast<float>(dif.size());
	
	//This is the Mahalnobis distance of the shape
	for (unsigned int i=0;i<vars.size();i++){
		binner+=vars.at(i)*vars.at(i);
	}
	
	//this is the intensity portion of the cost function
	cost+= -(k1)/2*log((alp+k2)/(alp+binner*(gammav)))+(alp+k1+k2)/2*log(1+(M-1)*probIcond/(alp+binner*(gammav)));
	
	//this is the shaoe prior ...chooeses between no condition, one conditional, or 2 conditionals
	//add in shape prior term, it can handle 0,1, or 2 conditionals
	if (shcond.value())			{		
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
		
		//p(x) prior--shape prior. No Conditional
		cost+=(alp+k2)/2*log(1+binner*gammav/alp);
	}
	
	return cost;
	
}

void negGradient(vector<float> *grad, vector<float> p, shapeModel* model1, int ex, int varbeg, int varsize, int& maxind, vector<float> fvals, vector<bool> select,int costtype){
	//costtype=0 gradient ASM
	//costtype=1 AAM
	//costtype=2 use Distance weighted Dice
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
	//take gradient with respect to each mode
	for (int i=0; i<static_cast<int>(p.size()); i++){
		//select vector is a bool vector that selecting the mdoes over which to take gradient
		if ( select.at(i)){
			
			//save original mode parameters
			gradtmp=p;	
			//increment the mode parameters
			gradtmp.at(i)=gradtmp.at(i)+searchDim;
			//order is reverse because negative gradient
			//evaluate appropriate cost 
			if (costtype==0){
				grad->at(i)=((costinit-costfunc(model1, gradtmp,fvals, ex,0))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)));
			}else if (costtype==1){
				grad->at(i)=((costinit-costfuncApp(model1, gradtmp,fvals, ex,0,globMean))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)))
			}else if (costtype==2){
				grad->at(i)=((costinit-costfuncOverlap(model1, gradtmp,fvals, ex,0))/(searchDim*sqrt(model1->getEigenValue(i))));//(pprev.at(i)))
			}
			
			//this handles hard max on mode parameters, not needed in practice
			if (abs(gradtmp.at(i))>stdTrunc.value()){
				//ignores gradient if goes beyond truncation
				grad->at(i)=grad->at(i)*1e-11;
			}
			
			//calculate cost difference in opposite direction
			//make sure gradient not positive in both directions
			gradtmp.at(i)=gradtmp.at(i)-2*searchDim;
			
			//evaluate appropriate cost 
			if (costtype==0){
				opGrad=((costinit-costfunc(model1, gradtmp, fvals, ex,0))/(-searchDim*sqrt(model1->getEigenValue(i))));
			}else if (costtype==1){
				opGrad=((costinit-costfuncApp(model1, gradtmp, fvals, ex,0,globMean))/(-searchDim*sqrt(model1->getEigenValue(i))));
			}else if (costtype==2){
				opGrad=((costinit-costfuncOverlap(model1, gradtmp, fvals, ex,0))/(-searchDim*sqrt(model1->getEigenValue(i))));
			}
			
			
			//impose rules 
			if ((opGrad>0)&&(grad->at(i)<0)){
				//if in a valley
				grad->at(i)=0;
			}else{
				//take central derivative
				grad->at(i)=(opGrad+grad->at(i))/2;
			}
			
			
			sumsq+=grad->at(i)*grad->at(i);
		}else{
			//set gradient to zero when not using that mode
			grad->at(i)=0;
		}
  }
	
	//handles cases of zero gradient
	if (sumsq==0){
		gradzero=true;
	}
	//nromalize such that size is one
	for (unsigned int i=0; i<p.size(); i++){
		grad->at(i)/=sqrt(sumsq);	
	}
}
void conjGradient(shapeModel* model1, vector<float>* vars, int varbeg, int varsize, int ex, vector<float> relStd, vector<float> fvals, vector<bool> select, int costtype, float searchR, float searchRmax){
	//costtype defines whether to use AAM or ASM
	//define variable
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
			cout<<"return, no gradient"<<endl;
			return ;//if graident == 0
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
				if (verbose.value()){
					cout<<"vars:"<<endl;
					for (unsigned int i=0;i<vars->size();i++){
						cout<<vars->at(i)<<" ";
					}
					cout<<endl;
				}
				
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
					
					gradzero=false;
					//costfunc app will set global means
				//	cout<<"Leaving conjgraident with mode "<<globMean<<endl;
					return ;
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
								
				}
			//	addVolumeTo4D(model1,*vars);
			}

			
			

void fitModel(shapeModel* model1, string modelName, vector<float>* vars, vector<float>& relStd, int start, int length, int excl, vector<bool> select, int costtype, bool newmodel, float searchR, float searchRmax){
//  cout<<"fitting "<<modelName<<endl;
	
	
  if (newmodel){
    model1->clear();
    //only use appearnace models now
    model1->load_bmv(modelName,1);
 //   cout<<"model successfully loaded"<<endl;
  }

	
  //structure weighting term
  vector<float> fvals;
  for (int i=0; i<model1->getNumberOfShapes();i++){
    fvals.push_back(1);
  }
  
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
	if (verbose.value()){
	cout<<"normalize intensity..."<<endl;
	}
	image=(image-robmin.value())*255/(robmax.value()-robmin.value());
	
		
	if (image.left_right_order()==1){
	  cout<<"Converting to radiological format"<<endl;
		image.makeradiological();
	}
	
	
	
	
	//cout<<image.getextrapolation()<<endl;
	//image.setpadvalue(0);
	//save_volume(image,"intnorm");
		if (manMean.value()==-777){
			GmanMean=manMean.value();
		}else{
			GmanMean=(manMean.value()-robmin.value())*255/(robmax.value()-robmin.value());
		}
	
	
	
	
		
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
	if (verbose.value()){
	cout<<"load model... "<<modelname.value()<<endl;
	}
		  model1->load_bmv_binaryInfo(modelname.value(),1);
		  model1->load_bmv_binary(modelname.value(),1);	
		  if (verbose.value()){
		  cout<<"register model..."<<endl;
		  }
		  //this will set new image parameters
		  model1->modelReg(0, flirtmatname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );
		  
		  if (useIntRefModel.value()){
			  shapeModel* modelRef =new shapeModel;
			  modelRef->setImageParameters(182,218, 182, 1, 1, 1);
			  if (verbose.value()){
			  cout<<"use a reference model "<<modelname2.value()<<endl;
			  }
			  //this bit of does uses one strcuture as a reference for others.
			  modelRef->load_bmv_binaryInfo(modelname2.value(),1);
			  modelRef->load_bmv_binary(modelname2.value(),1);	
			  if (verbose.value()){
			  cout<<"register model..."<<endl;
			  }
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
			  
			  vector<bool> select;
			  for (unsigned int i=0;i<varstemp.size();i++){
				  if(i<10){
					  select.push_back(true);
				  }else{
					  select.push_back(false);
				  }
			  }		
			  //conditional should be off
			  fitModel(modelRef,modelname2.value(),&varstemp,relStd,0,10,0,select,1, false,0.5, 0.15);
			  Mesh mtemp = modelRef->getDeformedMesh(varstemp,0,static_cast<int>(varstemp.size()));
			  int bounds[6]={0,0,0,0,0,0};
			  volume<short> mask;
			  mask=make_mask_from_meshInOut(image,mtemp,modelRef->getLabel(0),bounds);
			  
			  //calculate intensity histogram
			  float maxint=0,minint=0;
			  vector<float> vintens;
			  modelRef->intensityHistMaxMin(&image,&mask,&mtemp,modelRef->getLabel(0),&vintens, &maxint, &minint);
			  GmanMean=mode(vintens,minint,maxint);
			  if (verbose.value()){
			  cout<<"the mode of the reference distribution is "<<GmanMean<<endl;
			  }
		  }
		  
		  
		  vector<float> vars;
		  for (int i=0; i<model1->getNumberOfModes();i++){
			  //perturbe the system
			  vars.push_back(0);
		  }
		  
		  
		  
		  if (loadvars.value()){
			  
			  read_bvars(bvarsname.value(),&vars,model1->getNumberOfModes());
			  
		  }
		  
		  
		  
		  if (shcond.value()){
			  //	cout<<"shcond.value"<<endl;
			  int M;
			  Matrix mBx2;
			  Matrix mBcx1;
			  M=readBmap(bmapname.value(),vars,&mBx2,&mBcx1);
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
			  for (int i=0;i<static_cast<int>(vars.size());i++){
				  if (i<g.value()){
					  vars.at(i)=v_cmean.at(i);
				  }else{
					  vars.at(i)=0;
				  }
			  }
			  
			  //vars=v_cmean;
			  //				cout<<"done loading bmap stuff"<<endl;
			  if (shcond2.value()){
				  vector<float> vars2;
				  vars2=vars;
				  read_bvars(bvarsname2.value(),&vars,model1->getNumberOfModes());
				  M=readBmap2(bmapname2.value(),vars,&mBx2,&mBcx1);
				  ColumnVector Bx2(model1->getNumberOfModes());
				  //load bx2 vars that were read
				  for (unsigned int i=0;i<vars.size();i++){
					  Bx2.element(i)=vars.at(i);
				  }
				  mBx2map2=mBx2*Bx2;
				  mBx1inv2=mBcx1.i();
				  vector<float> v_cmean;
				  v_cmean=bTransform(&vars,mBx2,M);
				  
				  for (int i=0;i<static_cast<int>(vars.size());i++){
					  if (i<-1){
						  vars.at(i)+=v_cmean.at(i)+vars2.at(i);
					  }else{
						  vars.at(i)=0;
					  }
				  }
				  
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
		  
		  
		  lb=0;
		  ub=g.value()-1;
		  for (int i=0;i<static_cast<int>(vars.size());i++){
			  if ((i>=lb)&&(i<=ub)){
				  select.at(i)=true;
			  }else{
				  select.at(i)=false;
			  }
		  }
		  
		  if (baam.value()){
			  fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,1, false,0.5, 0.15);
		  }else if (overlap.value()){
			  fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,2, false,0.15, 0.1);
		  }else{ 
			  //use asm 
			  fitModel(model1,modelName,&vars,relStd,lb,ub-lb+1,0,select,0, false,0.15, 0.1);	
		  }
		  
		  
		Mesh mout = model1->getDeformedMesh(vars,0,static_cast<int>(vars.size()));
		  int bounds[6]={0,0,0,0,0,0};
			volume<short> imout=make_mask_from_meshInOut(image,mout, model1->getLabel(0),bounds);
			
		  save_volume(imout,outname.value());
		  mout.save(outname.value()+".vtk",3);
		  return 0;
}

int main(int argc,char *argv[])
{
	
	Tracer tr("main");
	OptionParser options(title, examples);
	
	try {
		// must include all wanted options here (the order determines how
		//  the help message is printed)
		options.add(inname);
		options.add(baam);
		options.add(loadvars);
		options.add(shcond);
		options.add(shcond2);
		options.add(bmapname);
		options.add(bmapname2);
		options.add(overlap);
		options.add(outname);
		options.add(stdTrunc);
		options.add(inputprob);
		options.add(verbose);
		options.add(help);
		options.add(modelname); 
		options.add(modelname2); 
		options.add(useIntRefModel);
		options.add(bvarsname);
		options.add(bvarsname2);
		options.add(flirtmatname);
		options.add(manMean); 
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
	
			do_work(argc,argv);
		
		
	}  catch(X_OptionError& e) {
		options.usage();
		cerr << endl << e.what() << endl;
		exit(EXIT_FAILURE);
	} catch(std::exception &e) {
		cerr << e.what() << endl;
	} 
	
	return 0;// do_work(argc,argv);
}

