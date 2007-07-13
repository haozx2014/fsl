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
Option<bool> outputzmap(string("--zmap"), false,
					   string("outputs zmap"),
					   false, no_argument);
					   
Option<bool> replaceBMV(string("--overideModel"), false,
					   string("allows you to specify a model other than in .bvars"),
					   false, no_argument);
Option<bool> overlap(string("--overlap"), false,
					   string("calculates overlap"),
					   false, no_argument);
Option<bool> dice(string("--dice"), false,
					   string("do stats on Dice measure"),
					   false, no_argument);
Option<bool> stats(string("--stats"), false,
					   string("do stats"),
					   false, no_argument);
Option<bool> noopen(string("--noopen"), false,
					   string("no morphological open operation"),
					   false, no_argument);
Option<bool> bcd(string("--bcd"), false,
					   string("use boundary corrected Dice"),
					   false, no_argument);
Option<bool> savebcorr(string("--savebcorr"), false,
					   string("save boundary corrected volume"),
					   false, no_argument);
Option<bool> dualbcorr(string("--dualbcorr"), false,
					   string("use to structures to boundary correct"),
					   false, no_argument);
Option<bool> concatbvars(string("--concatbvars"), false,
					   string("use to structures to boundary correct"),
					   false, no_argument);
Option<bool> noname(string("--noname"), false,
					   string("omits file names when writing bvars file"),
					   false, no_argument);
Option<bool> imfrombvars(string("--imfrombvars"), false,
					   string("make image from bavrs"),
					   false, no_argument);
Option<bool> usebvars(string("--usebvars"), false,
					   string("calculate volumes form bvars"),
					   false, no_argument);
Option<bool> doGLM(string("--doGLM"), false,
					   string("doGLM to egenrate tstats"),
					   false, no_argument);
					   
					   
Option<string> inname(string("-i,--in"), string(""),
					  string("filename of input image/mesh/bvars"),
					  true, requires_argument);
Option<string> flirtmatsname(string("-f,--in"), string(""),
					  string("filename of flirt matrices"),
					  false, requires_argument);
Option<string> modelname(string("-m,--in"), string(""),
					  string("filename of input model, if overiding .bvars"),
					  false, requires_argument);
Option<string> cmaname(string("-g,--in"), string(""),
					  string("filename of gold standard image"),
					  false, requires_argument);

Option<string> refname(string("-r,--in"), string(""),
					  string("filename of reference image "),
					  false, requires_argument);

Option<int> meshLabel(string("-l,--meshlabel"), 0,
					   string("If loading a vtk mesh, specify label used."),
					   false, requires_argument);
Option<int> numModes(string("-n,--numModes"), 0,
					   string("number of modes to retaon per structure."),
					   false, requires_argument);
Option<float> thresh(string("-p,--thrsh"), 2.94,
					   string("threshhold for clean up."),
					   false, requires_argument);

Option<string> outname(string("-k,--out"), string(""),
					   string("filename of output mesh"),
					   true, requires_argument);
Option<string> outoverlap(string("-o,--outoverlap"), string(""),
					   string("filename of text file containing overlaps"),
					   false, requires_argument);


int nonoptarg;

////////////////////////////////////////////////////////////////////////////
//global variables

//*************These are the overlap measures function**************************
int findLabel(int label1,vector<int>* vlabels){
	for (unsigned int i=0;i<vlabels->size();i++){
		if (vlabels->at(i)==label1){
			return i;
		}
	}
	return -1;
}

bool findAddLabel(int label1,int label2,int* indseg, vector<int>* vlabels, vector<int>* vTP, vector<int>* vFN, vector<int>* vFP, vector<int>* segImLabels,vector<int>* minInterX,vector<int>* maxInterX,vector<int>* minInterY,vector<int>* maxInterY,vector<int>* minInterZ,vector<int>* maxInterZ){
	int ind2=-1;//ind1=-1, ind2=-1;
	int ind1=-1;
	//return 1 signifies intersection
	
	//cout<<label1<<" "<<label2<<endl;
	*indseg=-1;
	for (unsigned int i=0;i<segImLabels->size();i++){
		if ((segImLabels->at(i)==label1)&&(label1!=0)){
			*indseg=i;
			i=segImLabels->size()+1;
		}
		
		
	} 
	if ((*indseg==-1)&&(label1!=0)){
		
		segImLabels->push_back(label1);
		minInterX->push_back(10000);
		maxInterX->push_back(0);
		minInterY->push_back(10000);
		maxInterY->push_back(0);
		minInterZ->push_back(10000);
		maxInterZ->push_back(0);
		*indseg=segImLabels->size()-1;
	}
	for (unsigned int i=0;i<vlabels->size();i++){	
		// cout<<"labels "<<label1<<" "<<label2<<endl;
		if ((vlabels->at(i)==label1)|(vlabels->at(i)==label2)){ 
			if (label1==label2){
				//	cout<<"found ind1=ind2"<<endl;
				ind1=ind2=i;
				i=vlabels->size()+1;
			}else{
				if (label1==vlabels->at(i)){
					ind1=i;
				}else{
					ind2=i;
				}
			}
		}
		if(i==vlabels->size()-1){
			// cout<<"do i enter?"<<endl;
			if (ind1==-1){
				//  cout<<"end and ind1=-1 "<<label1<<" "<<label2<<endl;
				ind1=vlabels->size();
				vlabels->push_back(label1);
				
				
				
				vTP->push_back(0);
				vFN->push_back(0);
				vFP->push_back(0);
			}
			if (ind2==-1){
				//  cout<<"end and ind2=-1 "<<label1<<" "<<label2<<endl;
				if (label1==label2){
					//	cout<<"end and lab1=lab2"<<endl;
					ind2=ind1;
				}else{
					//		cout<<"end and lab1!=lab2"<<endl;
					ind2=vlabels->size();
					vlabels->push_back(label2);
					vTP->push_back(0);
					vFN->push_back(0);
					vFP->push_back(0);
				}
				
			}
			
			
			
		}
		
	}
	
	//	cout<<"out of loop"<<endl;
	
	if (label1==label2){
		//cout<<"vTP "<<ind1<<endl;
		vTP->at(ind1)+=1;
		//cout<<"vTPout "<<ind1<<endl;
		return true;
	}else{
		//	  cout<<"vFP "<<ind1<<endl;
		vFP->at(ind1)+=1;
		//	cout<<"vFPout "<<ind1<<endl;
		vFN->at(ind2)+=1;
		//	cout<<"vFNout "<<ind2<<endl;
		return false;
	}
	
	//vCount->at(*ind1)+=1;
	
	
}

Matrix overlaps(const volume<short> segIm, const volume<short> gold){
	
	int sizex= segIm.xsize();
	int sizey=segIm.ysize();
	int sizez=segIm.zsize();
	bool inter;
	int indseg;
	//find union and intersection
	vector<int> vlabels,vTP, vFN, vFP,segLabels;
	//these are used to speed up the distance calculation
	vector<int> minInterX,maxInterX,minInterY,maxInterY,minInterZ,maxInterZ;
	vlabels.push_back(0);
	vTP.push_back(0);
	vFN.push_back(0);
	vFP.push_back(0);
	for (int k=0;k<sizez;k++){
		for (int j=0;j<sizey;j++){ 
			for (int i= 0; i<sizex;i++){
				//    cout<<"prefindlabel "<<i<<" "<<j<<" "<<k<<endl;
				inter=findAddLabel(segIm.value(i,j,k), gold.value(i,j,k),&indseg,&vlabels,&vTP, &vFN, &vFP,&segLabels,&minInterX,&maxInterX,&minInterY,&maxInterY,&minInterZ,&maxInterZ);
				//      cout<<"postfindlabel "<<ind1<<endl;
				if ((inter)&&(segIm.value(i,j,k)!=0)){
					//	cout<<minInterX.size()<<" "<<segIm.value(i,j,k)<<" "<<gold.value(i,j,k)<<endl;
					if (i<minInterX.at(indseg)){
						minInterX.at(indseg)=i;
					}
					if (j<minInterY.at(indseg)){
						minInterY.at(indseg)=j;
					}
					if (k<minInterZ.at(indseg)){
						minInterZ.at(indseg)=k;
					}
					if (i>maxInterX.at(indseg)){
						maxInterX.at(indseg)=i;
					}
					if (j>maxInterY.at(indseg)){
						maxInterY.at(indseg)=j;
					}
					if (k>maxInterZ.at(indseg)){
						maxInterZ.at(indseg)=k;
					}
					//	cout<<"leave min max update "<<endl;
				}
			}
		}
	}
	
	
	//this overlap does not give distance weighted stuff
	Matrix simMeasures(static_cast<int>(segLabels.size()),9);
	
	for (unsigned int i=0; i<segLabels.size();i++){
		//cout<<"add measures "<<endl;
		//for each label in segmentation
		//sum intersections
		int ind=findLabel(segLabels.at(i),&vlabels);

		simMeasures.element(i,0)=vlabels.at(ind);
		simMeasures.element(i,1)=vTP.at(ind);
		simMeasures.element(i,2)=vFP.at(ind);
		simMeasures.element(i,3)=vFN.at(ind);
		simMeasures.element(i,4)=static_cast<float>(vTP.at(ind))/(vTP.at(ind)+vFP.at(ind)+vFN.at(ind));
		simMeasures.element(i,5)=2.0*vTP.at(ind)/(2*vTP.at(ind)+vFP.at(ind)+vFN.at(ind));
		simMeasures.element(i,6)=static_cast<float>(vFN.at(ind))/(vFN.at(ind)+vTP.at(ind));
		simMeasures.element(i,7)=static_cast<float>(vFP.at(ind))/(vFN.at(ind)+vTP.at(ind));
		simMeasures.element(i,8)=static_cast<float>(vTP.at(ind))/(vFN.at(ind)+vTP.at(ind));
	}
	return simMeasures;
}

//****************************END OF OVERLAP FUNCTIONS**************************************
//****************************MESH FILLING FUNCTIONS**************************************
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

//****************************END OF MESH FILLING FUNCTIONS**************************************
//****************************BVARS I/O**************************************
string read_bvars(string fname,Matrix* bvars,int M,vector<string>* vnames, vector<int>* vnvars){
	string stemp;
	string modelNames;
	int N;//number of subjects
		ifstream fin;
		fin.open(fname.c_str());
		//throw away first three lines 
		getline(fin,stemp);//this is bvars file
			getline(fin,modelNames);//modelnames
				fin>>stemp>>N;
				bvars->ReSize(M,N);
				vnames->clear();
				vnvars->clear();
				
				//transform all the bvars
				for (int i=0; i<N;i++){
					fin>>stemp;//read in subject id
					vnames->push_back(stemp);
					
					int nvars;//how many vars written for the subject
						fin>>nvars;
						vnvars->push_back(nvars);
						for (int j=0;j<M;j++){
							if (j<nvars){
								float ftemp;
								fin>>ftemp;
								bvars->element(j,i)=ftemp;
							}else{
								bvars->element(j,i)=0;
							}
						}
				}
				return modelNames;
}
void write_bvars(string fname,string modelname,Matrix bvars, int numModes,vector<string> vnames){
	ofstream fout;
	
	fout.open(fname.c_str());
	fout<<"this is a bvars file"<<endl; 
	fout<<modelname<<endl;
	fout<<"NumberOfSubjects "<<bvars.Nrows()<<endl;

	for (int i=0;i<bvars.Nrows();i++){
	if (!noname.value()){
		fout<<vnames.at(i)<<" ";
		fout<<numModes<<" ";
		}
#ifdef PPC64
		int n=0;
#endif
		for (int j=0;j<bvars.Ncols();j++){
			fout<<bvars.element(i,j)<<" ";
#ifdef PPC64
			if ((n++ % 50) == 0) fout.flush();
#endif
		}
		fout<<endl;
	}
	
	fout<<endl;
	fout.close();
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
				
				getline(fin,stemp);
				design->ReSize(N,Nx);
				for (int i=0; i<N;i++){
					for (int j=0;j<Nx;j++){
						float ftemp;
						fin>>ftemp;			
						design->element(i,j)=ftemp;
					}
				}
				
				return modelNames;
}
//****************************END OF BVARS I/O**************************************


//****************************BOUNDARY CORRECTION FUNCTIONS**************************************
float mode(vector<float> vdists, float min, float max){
	int N=static_cast<int>(vdists.size());
	float bins=128;

	float binwidth=(max-min)/bins;
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
	bins=128;
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

float boundaryCorr(volume<short>* mask, volume<float>* ref, int label, float zthresh, int* bounds){
	//returns volume
	//build intensity vector 
	vector<float> vgraylevels;
	vector <float>::iterator Iter;
	float dist=10000;
	
	for (int i=bounds[0];i<bounds[1];i++){
		for (int j=bounds[2];j<bounds[3];j++){
			for (int k=bounds[4];k<bounds[5];k++){		
				
				if (mask->value(i,j,k)==label){
					dist=ref->value(i,j,k);
					
					if (vgraylevels.empty()){
						vgraylevels.push_back(dist);
					}else if (dist>=vgraylevels.back()){
						vgraylevels.push_back(dist);
					}else {
						for ( Iter = vgraylevels.begin( ) ; Iter !=vgraylevels.end( ) ; Iter++ ){
							if (dist<*Iter){
								
								vgraylevels.insert(Iter,dist);
								break;
							}
							
						}
						
						
					}
					
				}
			}
		}
	}
	int maxcount;
	//don't end up using the mode...re-calculate centre from fullwidth half-maximum
	mode(vgraylevels,&maxcount);
	//now find full width half maximum
	//float halfmaxval=maxcount/2.0;
	float halfmin,halfmax;
	fullwidthhalfmax(vgraylevels,maxcount/2.0,&halfmin, &halfmax);
	float mean=(halfmin+halfmax)/2;
	float sdev=abs(halfmax-halfmin)/2.35;
	
	float vol=0;
	
	for (int i=bounds[0];i<bounds[1];i++){
		for (int j=bounds[2];j<bounds[3];j++){
			for (int k=bounds[4];k<bounds[5];k++){
				float z=0.0;
				if (mask->value(i,j,k)==(label+100  )){
					z=(ref->value(i,j,k)-mean)/sdev;
					//maskz.value(i,j,k)=z;
					if (zthresh>=0){
						if (abs(z)>zthresh){
							mask->value(i,j,k)=0;
						}else{
							mask->value(i,j,k)=label;
							vol++;			
						}
					}
				}else if (mask->value(i,j,k)==label){
					vol++;
				}
			}
		}
	}
	return vol;
}

volume<float> boundaryCorrZ(volume<short> mask, volume<float>* ref, int label, int* bounds){
	//returns zstats
	volume<float> zvol;
	copyconvert(mask, zvol);
	zvol=0;
		//build intensity vector 
		vector<float> vgraylevels;
	vector <float>::iterator Iter;
		float dist=10000;
	
		//build up intesity statistics
		
		for (int i=bounds[0];i<bounds[1];i++){
			for (int j=bounds[2];j<bounds[3];j++){
				for (int k=bounds[4];k<bounds[5];k++){		
					
					if (mask.value(i,j,k)==label){
						dist=ref->value(i,j,k);
						
						if (vgraylevels.empty()){
							vgraylevels.push_back(dist);
						}else if (dist>=vgraylevels.back()){
							vgraylevels.push_back(dist);
						}else {
							for ( Iter = vgraylevels.begin( ) ; Iter !=vgraylevels.end( ) ; Iter++ ){
								if (dist<*Iter){
									
									vgraylevels.insert(Iter,dist);
									break;
								}
								
							}
							
							
						}
						
					}
				}
			}
		}
		int maxcount;
		mode(vgraylevels,&maxcount);
		float halfmin,halfmax;
		fullwidthhalfmax(vgraylevels,maxcount/2.0,&halfmin, &halfmax);	
		float mean=(halfmin+halfmax)/2;
		float sdev=abs(halfmax-halfmin)/2.35;

		for (int i=bounds[0];i<bounds[1];i++){
			for (int j=bounds[2];j<bounds[3];j++){
				for (int k=bounds[4];k<bounds[5];k++){
					float z=0.0;
					if (mask.value(i,j,k)==(label+100  )){
						z=(ref->value(i,j,k)-mean)/sdev;
						zvol.value(i,j,k)=z;	
					}
				}
			}
		}
		return zvol;
}
float BoundaryCorrectedDice(volume<short>* mask, volume<short>* ref, int label, int* bounds){
	//returns volume
			float vol=0;
		
		for (int i=bounds[0];i<bounds[1];i++){
			for (int j=bounds[2];j<bounds[3];j++){
				for (int k=bounds[4];k<bounds[5];k++){
					if (mask->value(i,j,k)==(label+100  )){
						if (ref->value(i,j,k)==label){
						mask->value(i,j,k)=ref->value(i,j,k);
						}else{
						mask->value(i,j,k)=0;

						}
					}
				}
			}
		}
		return vol;
}
int findStructLabel(volume<short>* mask, int* bounds){
	//used to find label and set bounds
	int xmin=10000,ymin=10000,zmin=10000;
	int xmax=-1, ymax=-1, zmax=-1;
	int label=999;
	bool found=false;
		for (int i=bounds[0];i<bounds[1];i++){
			for (int j=bounds[2];j<bounds[3];j++){
				for (int k=bounds[4];k<bounds[5];k++){		
						if (mask->value(i,j,k)>0){
							if (xmin>i){ xmin=i; }
							if (ymin>j){ ymin=j; }
							if (zmin>k){ zmin=k; }
							if (xmax<i){ xmax=i; }
							if (ymax<j){ ymax=j; }
							if (zmax<k){ zmax=k; }
						
						}
						if ((mask->value(i,j,k)<100)&&(mask->value(i,j,k)!=0)&&(!found)){
							label=mask->value(i,j,k);
							found=true;
						//	break;
						}
				}
			//	if (found){break;}
			}
			//if (found){break;}
		}
		bounds[0]=xmin-1;
		bounds[1]=xmax+1;
		bounds[2]=ymin-1;
		bounds[3]=ymax+1;
		bounds[4]=zmin-1;
		bounds[5]=zmax+1;
		return label;
}


//****************************END OF BOUNDARY CORRECTION FUNCTIONS**************************************
//*****************************EROSION/DILATION********************************************************
volume<short> myerode(volume<short> mask, int* bounds, int label){
	//returns volume
	//use a 3 cubed neighbpourhood
		volume<short> mero;
		mero=mask;
		mero=mero*0;
		for (int i=bounds[0];i<bounds[1];i++){
			for (int j=bounds[2];j<bounds[3];j++){
				for (int k=bounds[4];k<bounds[5];k++){		
					if (mask.value(i,j,k)!=0){
						bool zero=false;
						for (int l=-1;l<=1;l++){
							for (int m=-1;m<=1;m++){
								for (int n=-1;n<=1;n++){
									if (mask.value(i+l,j+m,k+n)==0){
										zero=true;
										break;
									}
								}
								if (zero){
								break;
								}
							}
							if (zero){break;}
						}
						if (!zero){
							mero.value(i,j,k)=label;
						}
					}
				}
			}
		}
		return mero;
}

volume<short> mydilateM(volume<short> mask,volume<short> maskB, int* bounds,int label, float* vol){
	//returns volume
	//use a 3 cubed neighbpourhood
	//this also adds back the centre
	//mask b has label and 100+label
	*vol=0;
	volume<short> mdil;
		mdil=mask;
		mdil=mdil*0;
		for (int i=bounds[0];i<bounds[1];i++){
			for (int j=bounds[2];j<bounds[3];j++){
				for (int k=bounds[4];k<bounds[5];k++){	
					
					if ((mask.value(i,j,k)==0)&&(maskB.value(i,j,k)!=label)){
					//	short mean=0;
						bool pos=false;
						for (int l=-1;l<=1;l++){
							for (int m=-1;m<=1;m++){
								for (int n=-1;n<=1;n++){
									if (mask.value(i+l,j+m,k+n)!=0){
										pos=true;
										break;
									}
								}
								if (pos){break;}
							}
							if (pos){break;}
						}
						if (pos){
							//mdil.value(i,j,k)=static_cast<short>(round(mean/27));
							mdil.value(i,j,k)=label;
							*vol = *vol +1;
						}
					}else{
						mdil.value(i,j,k)=label;
						*vol = *vol +1;
					}
				}
			}
		}
		return mdil;	
}
//*****************************END OF EROSION/DILATION********************************************************//
//*****************************LINEAR TRANSFORM********************************************************//
Matrix rigid_linear_xfm(Matrix Data,ColumnVector meanm, Mesh mesh){
	//ColumnVector avgM(sub.Ncols());
	cout<<"Data "<<Data.Nrows()<<" "<<Data.Ncols()<<endl;
	//determine translations
	int Nsub=Data.Ncols();
	int Npoints=Data.Nrows()/3;

	
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
	return DataNew;			
	//return V;
}
//********************************End Linear Transforms******************************//

//********************************GLM for stats******************************//
ColumnVector GLM_fit(Matrix G, Matrix D, ColumnVector contrast){

	//start for only well connditioned design matrices
	Matrix A=G.t()*G;
	Matrix Betas(D.Nrows(),D.Ncols());
	Betas=A.i()*G.t()*D;
	
	//calculate residula variance
	ColumnVector avgRes=(D-G*(Betas)).t()*(D-G*(Betas))/(G.Nrows()-G.Ncols());

				//convert to standard error
				avgRes=avgRes*(contrast.t()*A.i()*contrast).AsScalar();		
				Matrix test;
				test=contrast.t()*Betas.SubMatrix(1,Betas.Nrows(),1,1)*Betas.SubMatrix(1,Betas.Nrows(),1,1).t()*contrast;
				
				ColumnVector tstats(avgRes.Nrows());
				for (int i=0; i<avgRes.Nrows();i++){
					tstats.element(i)=(contrast.t()*Betas.SubMatrix(1,Betas.Nrows(),i+1,i+1)/sqrt(avgRes.element(i))).AsScalar();
				}
				return tstats;
}


//******************************EXECUTION FUNCTIONS******************************

//This calculates volumes from a bvars summary file (bounary corrections included).
void do_work_bvars(){
	
	//load the reference volume (mni152_1mm
	volume<float> ref;
	read_volume(ref,refname.value());
	
	//read in bvars
	string mname;
	mname=read_bvars_ModelName(inname.value() );	
	if (replaceBMV.value()){
		mname=modelname.value();
	}
	cout<<"load the model: "<<mname<<endl;
	//load model 
	shapeModel* model1 = new shapeModel;
	model1->setImageParameters(ref.xsize(), ref.ysize(), ref.zsize(), ref.xdim(), ref.ydim(), ref.zdim());
	model1->load_bmv_binaryInfo(mname,1);
	model1->load_bmv_binary(mname,1);
		
	//need number of subjects (to detrmine number of modes)
	int M=model1->getNumberOfSubjects();
	//read bvars into a matrix, pad zeros if bvars not provided
	vector<string> subjectnames;
	Matrix bvars;
	vector<int> vN;
	
	//target is only used for glm
	Matrix target;
	if (doGLM.value()){
		//must include a design matrix with bvars to use glm
		read_bvars_design(inname.value(),&bvars,&target,&subjectnames, &vN);
	}else{
		read_bvars(inname.value(),&bvars,M,&subjectnames, &vN);
	}
	//this is here to replace the appmodel with the sh model
	cout<<"You were using the model: "<<endl;
	cout<<mname<<endl;
	
	
	volume<float> t1im;
	volume<short> segim;
	
	//need flirt matrices
	ifstream flirtmats;
	flirtmats.open(flirtmatsname.value().c_str());
	int bounds[6]={0,0,0,0,0,0};
	
	
	//output volumes
	Matrix Volumes(target.Nrows(),1); //this is only used for GLM
	ofstream fvol_all;
	fvol_all.open(outname.value().c_str());
	
	
	//calculate volume for all subbjects
	for (int subject=0;subject<bvars.Ncols();subject++){
		cout<<"Please make sure images and transformation matrices correspond"<<endl;
		model1->clear();
		model1->setImageParameters(ref.xsize(), ref.ysize(), ref.zsize(), ref.xdim(), ref.ydim(), ref.zdim());
		model1->load_bmv_binaryInfo(mname,1);	
		model1->load_bmv_binary(mname,1);
		
		cout<<"Image: "<<subjectnames.at(subject)<<endl;
		read_volume(t1im,subjectnames.at(subject));
		string flirtname;
		flirtmats>>flirtname;
		
		cout<<"Flirt matrix "<<flirtname<<endl;
		//need to load correspoding flirst matrix
		model1->modelReg(0, flirtname,t1im.xsize(),t1im.ysize(), t1im.zsize(), t1im.xdim(),t1im.ydim() ,t1im.zdim() );
		
		//load bvars into vector, to deform mesh
		vector<float> vars;
		for (int i=0; i<bvars.Nrows();i++){
			vars.push_back(bvars.element(i,subject));
		}
		
		//deform mesh			
		//assume shape 0
		Mesh m = model1->getDeformedMesh(vars,0,static_cast<int>(vars.size()));
		//fill mesh
		segim=make_mask_from_meshInOut(t1im,m,model1->getLabel(0),bounds);
		
		
		stringstream sstemp;
		string outnamess;
		sstemp<<subject;
		sstemp>>outnamess;

		float voltemp;
		volume<short> segimB;
		segimB=segim;
		voltemp=boundaryCorr(&segim, &t1im,model1->getLabel(0), thresh.value(),bounds);
		if (!noopen.value()){
			segim=myerode(segim,bounds, model1->getLabel(0));
			segim=mydilateM(segim,segimB,bounds,model1->getLabel(0),&voltemp);
		}
		fvol_all<<subjectnames.at(subject)<<" "<<voltemp<<" "<<voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim()<<endl;
		//save volume to a matrix if you wish to use in GLM
		if (doGLM.value()){
			Volumes.element(subject,0)=voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim();
		}
		
	}
	
	if (doGLM.value()){
		cout<<"Volumes calculated, now runnning GLM"<<endl;
		ColumnVector tstats;
		ColumnVector contrast(2);
		contrast.element(0)=-1;
		contrast.element(1)=1;
		tstats=GLM_fit(target, Volumes,contrast);
		for (int i =0 ; i <tstats.Nrows();i++){
			cout<<tstats.element(i)<<endl;
		}
	}
	
}


//this function calculates overlap form volumes (bounary corrections included)
void do_work_vols(){
	
	ifstream fims;
	fims.open(inname.value().c_str());
	
	volume<float> t1im;
	volume<short> segim;
	volume<short> cmaim;
	
	//output volumes
	ofstream fvol_all;
	fvol_all.open(outname.value().c_str());
	ofstream foverlap;
	if (overlap.value()){
		foverlap.open(outoverlap.value().c_str());
	}
	int bounds[6]={0,0,0,0,0,0};
	
	
	while (!fims.eof() ){
		string t1name,segname;
		fims>>t1name;
		cout<<"image "<<t1name<<endl;
		read_volume(t1im,t1name);
		fims>>segname;
		read_volume(segim,segname);
		
		//use bounds of whole image
		//need to reset lower bounds as well
		bounds[0]=0;
		bounds[2]=0;
		bounds[4]=0;
		bounds[1]=segim.xsize();
		bounds[3]=segim.ysize();
		bounds[5]=segim.zsize();
		
		float voltemp;
		volume<short> segimB;
		segimB=segim;
		
		//need to find the structure label
		//also sets bounds
		//cout<<"find label "<<" "<<bounds[0]<<" "<<bounds[1]<<" "<<bounds[2]<<" "<<bounds[3]<<" "<<bounds[4]<<" "<<bounds[5]<<" "<<endl;
		
		int label=findStructLabel(&segim, bounds);
		//cout<<"found label do boundary corr "<<label<<" "<<bounds[0]<<" "<<bounds[1]<<" "<<bounds[2]<<" "<<bounds[3]<<" "<<bounds[4]<<" "<<bounds[5]<<" "<<endl;
		if (!bcd.value()){
			voltemp=boundaryCorr(&segim, &t1im,label, thresh.value(),bounds);
			if (!noopen.value()){
				segim=myerode(segim,bounds, label);
				segim=mydilateM(segim,segimB,bounds,label,&voltemp);
			}
			if (savebcorr.value()){
				save_volume(segim,"bcorr");
			}
		}
		Matrix simMeas;
		if (overlap.value()){
			string cmaname;
			fims>>cmaname;
			read_volume(cmaim,cmaname);
			if (bcd.value()){
				voltemp=BoundaryCorrectedDice(&segim, &cmaim,label,bounds);
			}
			simMeas=overlaps(segim,cmaim);
			foverlap<<segname<<" ";
			if (simMeas.Nrows()>1){
				cerr<<"should only have one label!"<<endl;
			}else{
				for (int i=0; i<simMeas.Ncols();i++){
					foverlap<<simMeas.element(0,i)<<" ";
				}
			}
			foverlap<<endl;
		}

		fvol_all<<segname<<" "<<voltemp<<" "<<voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim()<<endl;
	}
	
}
void do_work_dualbcorr(){
	
	//load the reference volume (mni152_1mm
	//	volume<float> ref;
	//	read_volume(ref,refname.value());
	
	ifstream fims;
	fims.open(inname.value().c_str());
	
	volume<float> t1im;
	volume<short> segim;
	volume<short> segim2;
	volume<short> cmaim;
	
	//output volumes
	ofstream fvol_all;
	fvol_all.open(outname.value().c_str());
	ofstream foverlap;
	if (overlap.value()){
		foverlap.open(outoverlap.value().c_str());
	}
	int bounds[6]={0,0,0,0,0,0};
	int bounds2[6]={0,0,0,0,0,0};
	
	while (!fims.eof() ){
		string t1name,segname, segname2;
		fims>>t1name;
		int length=t1name.length();
		cout<<t1name<<" "<<length<<endl;
		if (length==0){
			break;
		}
		
		cout<<"image "<<t1name<<endl;
		read_volume(t1im,t1name);
		fims>>segname;
		cout<<"image "<<segname<<endl;
		read_volume(segim,segname);
		fims>>segname2;
		cout<<"image "<<segname2<<endl;
		read_volume(segim2,segname2);
		
		//use bounds of whole image
		//need to reset lower bounds as well
		bounds[0]=0;
		bounds[2]=0;
		bounds[4]=0;
		bounds[1]=segim.xsize();
		bounds[3]=segim.ysize();
		bounds[5]=segim.zsize();
		
		bounds2[0]=0;
		bounds2[2]=0;
		bounds2[4]=0;
		bounds2[1]=segim2.xsize();
		bounds2[3]=segim2.ysize();
		bounds2[5]=segim2.zsize();
		
		float voltemp,voltemp2;
		//volume<short> segimB;
		//segimB=segim;
		save_volume(segim,"im1a");
		save_volume(segim2,"im2a");
		//need to find the structure label
		//also sets bounds
		//cout<<"find label "<<" "<<bounds[0]<<" "<<bounds[1]<<" "<<bounds[2]<<" "<<bounds[3]<<" "<<bounds[4]<<" "<<bounds[5]<<" "<<endl;
		
		int label=findStructLabel(&segim, bounds);
		int label2=findStructLabel(&segim2, bounds2);
		int bounds3[6]={0,0,0,0,0,0};
		
		
		if (bounds2[0]<bounds[0]) {bounds3[0]=bounds2[0];} else {bounds3[0]=bounds[0];}
		if (bounds2[1]>bounds[1]) {bounds3[1]=bounds2[1];} else {bounds3[1]=bounds[1];}
		if (bounds2[2]<bounds[2]) {bounds3[2]=bounds2[2];} else {bounds3[2]=bounds[2];}
		if (bounds2[3]>bounds[3]) {bounds3[3]=bounds2[3];} else {bounds3[3]=bounds[3];}
		if (bounds2[4]<bounds[4]) {bounds3[4]=bounds2[4];} else {bounds3[4]=bounds[4];}
		if (bounds2[5]>bounds[5]) {bounds3[5]=bounds2[5];} else {bounds3[5]=bounds[5];}
		
		cout<<"find label "<<" "<<bounds3[0]<<" "<<bounds3[1]<<" "<<bounds3[2]<<" "<<bounds3[3]<<" "<<bounds3[4]<<" "<<bounds3[5]<<" "<<endl;
		
		//	voltemp=boundaryCorr(&segim, &t1im,label, thresh.value(),bounds);
		//	voltemp2=boundaryCorr(&segim2, &t1im,label, thresh.value(),bounds);
		//first pass through corrects assigns voxels that are interior to one and exterior to the other	
		for (int i=bounds3[0];i<bounds3[1];i++){
			for (int j=bounds3[2];j<bounds3[3];j++){
				for (int k=bounds3[4];k<bounds3[5];k++){	
					
					//THIS IS THE RULE OF AMYGDALA PRIORITY OF HIPPOCAMPUS
					if ( (segim.value(i,j,k)==label) &&  (segim2.value(i,j,k)==(100+label2)) ){
						segim2.value(i,j,k)=0;
					}else if((segim.value(i,j,k)==(100+label)) &&  (segim2.value(i,j,k)==label2) ){
						segim.value(i,j,k)=0;
					}
				}
			}
		}
		cout<<"ims firstpass done"<<endl;	
		voltemp=boundaryCorr(&segim, &t1im,label, thresh.value(),bounds3);
		voltemp2=boundaryCorr(&segim2, &t1im,label2, thresh.value(),bounds3);
		cout<<"ims firstpass"<<endl;	
		
		//second pass gove applies rules of priority, i.e. choose amygdala over hippocampus	
		for (int i=bounds3[0];i<bounds3[1];i++){
			for (int j=bounds3[2];j<bounds3[3];j++){
				for (int k=bounds3[4];k<bounds3[5];k++){	
					if ( (segim.value(i,j,k)==17) &&  (segim2.value(i,j,k)==(18)) ){
						//segim2.value(i,j,k)=18;
						segim.value(i,j,k)=0;
						
					}else if((segim.value(i,j,k)==18) &&  (segim2.value(i,j,k)==(17)) ){
						segim2.value(i,j,k)=0;
					}
				}
			}
		}
		save_volume(segim,"im1");
		save_volume(segim2,"im2");
		
		Matrix simMeas;
		cout<<"start overlap"<<endl;
		if (overlap.value()){
			string cmaname;
			fims>>cmaname;
			read_volume(cmaim,cmaname);
			if (bcd.value()){
				voltemp=BoundaryCorrectedDice(&segim, &cmaim,label,bounds);
			}
			
			
			simMeas=overlaps(segim,cmaim);
			//out put similarity stats 
			//cout<<simMeas.Ncols()<<" "<<simMeas.Nrows()<<endl;
			foverlap<<segname<<" ";
			if (simMeas.Nrows()>1){
				cerr<<"should only have one label!"<<endl;
			}else{
				for (int i=0; i<simMeas.Ncols();i++){
					foverlap<<simMeas.element(0,i)<<" ";
				}
			}
			foverlap<<endl;
			
			simMeas=overlaps(segim2,cmaim);
			//out put similarity stats 
			//cout<<simMeas.Ncols()<<" "<<simMeas.Nrows()<<endl;
			foverlap<<segname2<<" ";
			if (simMeas.Nrows()>1){
				cerr<<"should only have one label!"<<endl;
			}else{
				for (int i=0; i<simMeas.Ncols();i++){
					foverlap<<simMeas.element(0,i)<<" ";
				}
			}
			foverlap<<endl;
			
			
		}
		//	cout<<"end overlap"<<endl;
		fvol_all<<segname<<" "<<voltemp<<" "<<voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim()<<endl;
	}
	
}
void do_work_dualbcorrInt(){
	
	//load the reference volume (mni152_1mm
//	volume<float> ref;
//	read_volume(ref,refname.value());
	
	ifstream fims;
	fims.open(inname.value().c_str());
	
	volume<float> t1im;
	volume<short> segim;
	volume<short> segim2;
	volume<short> cmaim;
	
	//output volumes
	ofstream fvol_all;
	fvol_all.open(outname.value().c_str());
	ofstream foverlap;
	if (overlap.value()){
		foverlap.open(outoverlap.value().c_str());
	}
	int bounds[6]={0,0,0,0,0,0};
	int bounds2[6]={0,0,0,0,0,0};
	
	while (!fims.eof() ){
		string t1name,segname, segname2;
		fims>>t1name;
		int length=t1name.length();
		cout<<t1name<<" "<<length<<endl;
		if (length==0){
			break;
			}

		cout<<"image "<<t1name<<endl;
		read_volume(t1im,t1name);
		fims>>segname;
		cout<<"image "<<segname<<endl;
		read_volume(segim,segname);
		fims>>segname2;
		cout<<"image "<<segname2<<endl;
		read_volume(segim2,segname2);
		
		//use bounds of whole image
		//need to reset lower bounds as well
		bounds[0]=0;
		bounds[2]=0;
		bounds[4]=0;
		bounds[1]=segim.xsize();
		bounds[3]=segim.ysize();
		bounds[5]=segim.zsize();
		
		bounds2[0]=0;
		bounds2[2]=0;
		bounds2[4]=0;
		bounds2[1]=segim2.xsize();
		bounds2[3]=segim2.ysize();
		bounds2[5]=segim2.zsize();
		
		float voltemp;
		//volume<short> segimB;
		//segimB=segim;
	save_volume(segim,"im1a");
		save_volume(segim2,"im2a");
		//need to find the structure label
		//also sets bounds
		//cout<<"find label "<<" "<<bounds[0]<<" "<<bounds[1]<<" "<<bounds[2]<<" "<<bounds[3]<<" "<<bounds[4]<<" "<<bounds[5]<<" "<<endl;
		
		int label=findStructLabel(&segim, bounds);
	int label2=findStructLabel(&segim2, bounds2);
	int bounds3[6]={0,0,0,0,0,0};
	
			
			if (bounds2[0]<bounds[0]) {bounds3[0]=bounds2[0];} else {bounds3[0]=bounds[0];}
			if (bounds2[1]>bounds[1]) {bounds3[1]=bounds2[1];} else {bounds3[1]=bounds[1];}
			if (bounds2[2]<bounds[2]) {bounds3[2]=bounds2[2];} else {bounds3[2]=bounds[2];}
			if (bounds2[3]>bounds[3]) {bounds3[3]=bounds2[3];} else {bounds3[3]=bounds[3];}
			if (bounds2[4]<bounds[4]) {bounds3[4]=bounds2[4];} else {bounds3[4]=bounds[4];}
			if (bounds2[5]>bounds[5]) {bounds3[5]=bounds2[5];} else {bounds3[5]=bounds[5];}
		
	cout<<"find label "<<" "<<bounds3[0]<<" "<<bounds3[1]<<" "<<bounds3[2]<<" "<<bounds3[3]<<" "<<bounds3[4]<<" "<<bounds3[5]<<" "<<endl;

	//	voltemp=boundaryCorr(&segim, &t1im,label, thresh.value(),bounds);
	//	voltemp2=boundaryCorr(&segim2, &t1im,label, thresh.value(),bounds);
		//first pass through corrects assigns voxels that are interior to one and exterior to the other	
		for (int i=bounds3[0];i<bounds3[1];i++){
			for (int j=bounds3[2];j<bounds3[3];j++){
				for (int k=bounds3[4];k<bounds3[5];k++){	
				
					//THIS IS THE RULE OF AMYGDALA PRIORITY OF HIPPOCAMPUS
					if ( (segim.value(i,j,k)==label) &&  (segim2.value(i,j,k)==(100+label2)) ){
						segim2.value(i,j,k)=0;
					}else if((segim.value(i,j,k)==(100+label)) &&  (segim2.value(i,j,k)==label2) ){
						segim.value(i,j,k)=0;
					}
				}
			}
		}
		cout<<"ims firstpass done"<<endl;	
	volume<float> vz1;
	vz1=boundaryCorrZ(segim, &t1im,label,bounds3);
	volume<float> vz2;
	vz2=boundaryCorrZ(segim2, &t1im,label2,bounds3);
			cout<<"ims firstpass"<<endl;	

	//second pass gove applies rules of priority, i.e. choose amygdala over hippocampus	
		for (int i=bounds3[0];i<bounds3[1];i++){
			for (int j=bounds3[2];j<bounds3[3];j++){
				for (int k=bounds3[4];k<bounds3[5];k++){	
					if ( (segim.value(i,j,k)==17) &&  (segim2.value(i,j,k)==(18)) ){
						//segim2.value(i,j,k)=18;
						//comparing z scores
						if (vz1.value(i,j,k)>=vz2.value(i,j,k)){
							segim.value(i,j,k)=0;
							segim2.value(i,j,k)=18;
						}
					}else if((segim.value(i,j,k)==18) &&  (segim2.value(i,j,k)==(17)) ){
						if (vz1.value(i,j,k)>=vz2.value(i,j,k)){
							segim.value(i,j,k)=0;
							segim2.value(i,j,k)=17;
						}
				
					}else if (segim.value(i,j,k)==18){
							if (vz1.value(i,j,k)>=thresh.value()){
							
							
							}
					}
					
					
					
				}
			}
		}
		save_volume(vz1,"im1z");
		save_volume(vz2,"im2z");
		save_volume(segim,"im1");
		save_volume(segim2,"im2");
					
		Matrix simMeas;
		cout<<"start overlap"<<endl;
		if (overlap.value()){
			string cmaname;
			fims>>cmaname;
			read_volume(cmaim,cmaname);
			if (bcd.value()){
				voltemp=BoundaryCorrectedDice(&segim, &cmaim,label,bounds);
			}
			
			
			simMeas=overlaps(segim,cmaim);
			//out put similarity stats 
			//cout<<simMeas.Ncols()<<" "<<simMeas.Nrows()<<endl;
			foverlap<<segname<<" ";
				if (simMeas.Nrows()>1){
					cerr<<"should only have one label!"<<endl;
				}else{
					for (int i=0; i<simMeas.Ncols();i++){
						foverlap<<simMeas.element(0,i)<<" ";
					}
				}
			foverlap<<endl;
			
			simMeas=overlaps(segim2,cmaim);
			//out put similarity stats 
			//cout<<simMeas.Ncols()<<" "<<simMeas.Nrows()<<endl;
			foverlap<<segname2<<" ";
				if (simMeas.Nrows()>1){
					cerr<<"should only have one label!"<<endl;
				}else{
					for (int i=0; i<simMeas.Ncols();i++){
						foverlap<<simMeas.element(0,i)<<" ";
					}
				}
			foverlap<<endl;
			
			
		}
	
	}
		   
}


void do_work_overlap_stats(){
	
	//flist is a list of overlap files to read
	ifstream flist;
	flist.open(inname.value().c_str());
	

	ofstream fout;
	fout.open(outname.value().c_str());
	ofstream fout2;
	string out2="names"+outname.value();
	fout2.open(out2.c_str());
	
	while (!flist.eof() ){
			float sx=0;
			float sxx=0;
		//name of overlap file
		string overname;
		getline(flist,overname);
		int length=overname.length();
		cout<<overname<<" "<<length<<endl;
		if (length==0){
			break;
			}
		ifstream fover;
		fover.open(overname.c_str());
		int n=0;
		cout<<"overname "<<overname<<endl;
		while (!fover.eof() ){
			//keep track of number of subjects
			n++;
			string imname;
			float overlapmeasure;
			fover>>imname;
			cout<<imname<<endl;
			int measInd=1;
			if (dice.value()){
				measInd=6;
			}
			
			
			//the are 10 outputs in my standard overlap output file
			//this caluclate stats for one
			for (int i=1;i<10;i++){
				fover>>overlapmeasure;
				
				if (i==6){
				//cout<<"Dice "<<overlapmeasure<<endl;
					sx+=overlapmeasure;
					sxx+=overlapmeasure*overlapmeasure;
				}
			}
			getline(fover,imname);
		}
		fover.close();
		fout<<sx/n<<" "<<(sxx-sx*sx/n)/(n-1)<<endl;
		fout2<<overname<<" "<<sx/n<<" "<<(sxx-sx*sx/n)/(n-1)<<endl;
		
	}
	fout.close();	
	fout2.close();	   
}
void do_work_concatbvars(){
	ifstream flist;
	flist.open(inname.value().c_str());
	
	int M=317;
	Matrix bvarsAll;
	int count=0;
	vector<string> subjectnames;
	while (!flist.eof() ){
		//name of overlap file
		string overname;
		getline(flist,overname);
		int length=overname.length();
		cout<<overname<<" "<<length<<endl;
		if (length==0){
			break;
		}
		ifstream fbvars;
		fbvars.open(overname.c_str());
		string stemp;
		
		
		Matrix bvars;
		vector<int> vN;
		
		read_bvars(overname,&bvars,M,&subjectnames, &vN);
		bvars=bvars.t();
		if (count==0){
			bvarsAll=bvars.SubMatrix(1,bvars.Nrows(),1,numModes.value());
		}else{
			bvarsAll=bvarsAll | bvars.SubMatrix(1,bvars.Nrows(),1,numModes.value());
		}
		
		
		
		count++;
	}
	write_bvars(outname.value(),modelname.value(),bvarsAll,numModes.value(),subjectnames);
	
	
}
void do_work_imfrombvars(){

	string mname;
	mname=read_bvars_ModelName(inname.value() );	
	shapeModel* model1 = new shapeModel;
		//mni152 1mm isotropic
	model1->setImageParameters(182,218, 182, 1, 1, 1);
	//model1->setImageParameters(image.xsize(), image.ysize(), image.zsize(), image.xdim(), image.ydim(), image.zdim());
	cout<<"load model..."<<endl;
	
	model1->load_bmv_binaryInfo(mname,1);
	model1->load_bmv_binary(mname,1);

	//read in info from bvars ...this is implemented for single bvar file
	vector<string> subjectnames;
	Matrix bvars;
	vector<int> vN;

	read_bvars(inname.value(),&bvars,model1->getNumberOfSubjects(),&subjectnames, &vN);

	volume<float> image;
	read_volume(image,subjectnames.at(0));
	
	volume<short> brainSeg;
	copyconvert(image,brainSeg);
	brainSeg=0;
	
	const int sizex = image.xsize();
	const int sizey = image.ysize();
	const int sizez = image.zsize();
	
		
	
	
	cout<<"register model..."<<endl;
	//this will set new image parameters
	model1->modelReg(0, flirtmatsname.value(), sizex,sizey, sizez, image.xdim(),image.ydim() ,image.zdim() );
	cout<<bvars.Nrows()<<" "<<bvars.Ncols()<<endl;
	
	vector<float> vars;
	for (int i =0; i<bvars.Nrows();i++){
		vars.push_back(bvars.element(i,0));
	
	}
	
	Mesh m = model1->getDeformedMesh(vars,0, static_cast<int>(vars.size()));
	int bounds[]={0,0,0,0,0,0};
	brainSeg=make_mask_from_meshInOut(image,m,model1->getLabel(0),bounds);

	save_volume(brainSeg,outname.value());
}

void do_work_glm_test(){
	
	
	//use to read in design matrix
	vector<string> subjectnames;
	Matrix bvars;
	vector<int> vN;
	Matrix target;
	cout<<"read in bvars and target"<<endl;
	read_bvars_design(inname.value(),&bvars,&target,&subjectnames, &vN);
	cout<<bvars.Nrows()<<" "<<bvars.Ncols()<<endl;

	//use to read in volumes
	ifstream fin;
	fin.open(refname.value().c_str());
	vector<float> vvols;
	while  (!fin.eof()){
		string stemp;
		float vol;
		fin>>vol;
		//getline(fin,stemp);
		cout<<vol<<" "<<vvols.size()<<endl;
		vvols.push_back(vol);
	}
	cout<<"read volumes"<<endl;
	cout<<static_cast<int>(vvols.size())<<endl;
	Matrix Mvols(static_cast<int>(vvols.size()),1);
	for (int i =0 ; i <static_cast<int>(vvols.size());i++){
		Mvols.element(i,0)=vvols.at(i);
	}
	
	ColumnVector tstats;
	ColumnVector contrast(2);
	contrast.element(0)=-1;
	contrast.element(1)=1;
	cout<<"do glm"<<endl;
	tstats=GLM_fit(target, Mvols,contrast);
	for (int i =0 ; i <tstats.Nrows();i++){
		cout<<tstats.element(i)<<endl;
	}

	
}


int main(int argc,char *argv[])
{
	
	Tracer tr("main");
	OptionParser options(title, examples);
	
	try {
		// must include all wanted options here (the order determines how
		//  the help message is printed)
		options.add(inname);
		options.add(cmaname);
		options.add(refname);
		options.add(flirtmatsname);
		options.add(replaceBMV);
		options.add(overlap);
		options.add(outoverlap);
		options.add(usebvars);
		options.add(doGLM);
		options.add(stats);
		options.add(dice);
		options.add(noopen);
		options.add(dualbcorr);
		options.add(bcd);
		options.add(modelname);
		options.add(outname);
		options.add(outputzmap);
		options.add(thresh);
		options.add(meshLabel);
		options.add(savebcorr);
		options.add(concatbvars);
		options.add(imfrombvars);
		options.add(noname);
		options.add(numModes);
		
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
		if (stats.value()){
			do_work_overlap_stats();
		}else if (dualbcorr.value()){
			do_work_dualbcorrInt();
		}else if(imfrombvars.value()){
			do_work_imfrombvars();
		}else if(concatbvars.value()){
			do_work_concatbvars();
		}else if(usebvars.value()){
			do_work_bvars();
		}else if (doGLM.value()){
			do_work_glm_test();
		}else{
			do_work_vols();
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

