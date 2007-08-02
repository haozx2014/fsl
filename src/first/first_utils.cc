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
#include "miscmaths/miscmaths.h"
using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace mesh;
using namespace shapemodel;
using namespace MISCMATHS;

string title="firt_utils (Version 1.0) University of Oxford (Brian Patenaude)";
string examples="first_utils [options] -i segImage -l goldStandard -k output.txt ";


Option<bool> verbose(string("-v,--verbose"), false, 
					 string("switch on diagnostic messages"), 
					 false, no_argument);
Option<bool> help(string("-h,--help"), false,
				  string("display this message"),
				  false, no_argument);
Option<bool> inputprob(string("--inputprob"), false,
					   string("Input image is probability - do not do any thresholding or smoothing"),
					   false, no_argument);
Option<bool> overlap(string("--overlap"), false,
					 string("calculates overlap"),
					 false, no_argument);

Option<bool> useScale(string("--useScale"), false,
					   string("do stats"),
					   false, no_argument);
Option<bool> vertexAnalysis(string("--vertexAnalysis"), false,
					   string("operate on vertices from bvars"),
					   false, no_argument);
Option<bool> singleBoundaryCorr(string("--singleBoundaryCorr"), false,
					   string("correct boundary voxels of a single structue"),
					   false, no_argument);
Option<bool> bcd(string("--bcd"), false,
			 string("use boundary corrected Dice"),
				 false, no_argument);
Option<bool> usePCAfilter(string("--usePCAfilter"), false,
						   string("pca filter set number of modes to retain"),
						   false, no_argument);
Option<bool> usebvars(string("--usebvars"), false,
						 string("make image from bavrs"),
						 false, no_argument);
Option<bool> useWilks(string("--useWilks"), false,
					   string("generate mesh with vertex wise t-stat vectors"),
					   false, no_argument);

Option<bool> useReconMNI(string("--useReconMNI"), false,
						 string("reconstruct meshes into MNI space"),
						 false, no_argument);

Option<bool> useReconNative(string("--useReconNative"), false,
							string("recomnstruct meshes in native space"),
							false, no_argument);
Option<bool> useRigidAlign(string("--useRigidAlign"), false,
						   string("register meshes using rigid registration"),
						   false, no_argument);

Option<bool> useNorm(string("--useNorm"), false,
					string("use normalization"),
					false, no_argument);

Option<bool> univariateT(string("--univariateT"), false,
						  string("use univariate T in vertex Analysis Mode"),
						  false, no_argument);
Option<bool> useVolumes(string("--useVolumes"), false,
						string("use volume discriminant"),
						false, no_argument);


Option<bool> meshToVol(string("--meshToVol"), false,
					string("convert mesh to a vloume"),
					false, no_argument);
Option<bool> centreOrigin(string("--centreOrigin"), false,
					string("places origin of mesh at the cnere of the image"),
					false, no_argument);

Option<string> inname(string("-i,--in"), string(""),
					  string("filename of input image/mesh/bvars"),
					  true, requires_argument);
Option<string> pathname(string("-a,--in"), string(""),
					  string("path to iimage"),
					  false, requires_argument);
Option<string> flirtmatsname(string("-f,--in"), string(""),
							 string("filename of flirt matrices"),
							 false, requires_argument);
Option<string> modelname(string("-m,--in"), string(""),
						 string("filename of input model, if overiding .bvars"),
						 false, requires_argument);
Option<string> normname(string("-g,--in"), string(""),
						string("filename of gold standard image"),
						false, requires_argument);
Option<string> designname(string("-d,--in"), string(""),
						string("filename of fsl design matrix"),
						false, requires_argument);

Option<string> refname(string("-r,--in"), string(""),
						string("filename of reference image "),
						false, requires_argument);

Option<int> meshLabel(string("-l,--meshlabel"), 1,
					  string("If loading a vtk mesh, specify label used."),
					  false, requires_argument);

Option<float> thresh(string("-p,--thrsh"), 4,
					 string("threshhold for clean up."),
					 false, requires_argument);
Option<int> numModes(string("-n,--numModes"), 0,
					 string("number of modes to retaon per structure."),
					 false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
					   string("filename of output mesh"),
					   true, requires_argument);

int nonoptarg;

////////////////////////////////////////////////////////////////////////////
//global variables

int refXsize=182;
int refYsize=218;
int refZsize=182;
float refXdim=1.0;
float refYdim=1.0;
float refZdim=1.0;







//%%%%%%%%%%%%%%mesh fill
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

  double xdim = (double) image.xdim();
  double ydim = (double) image.ydim();
  double zdim = (double) image.zdim();

  //in new version of bet2
  double mininc = min(xdim,min(ydim,zdim)) * .5;


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
		vTP->at(ind1)+=1;
		return true;
	}else{
		vFP->at(ind1)+=1;
		vFN->at(ind2)+=1;
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
		cout<<"Dice "<<vlabels.at(ind)<<" "<<2.0*vTP.at(ind)/(2*vTP.at(ind)+vFP.at(ind)+vFN.at(ind))<<endl;
	}
	return simMeasures;
}


//****************************BVARS I/O**************************************


string read_bvars(string fname,Matrix* bvars, vector<string>* vnames, vector<int>* vnvars,string impath){
	string stemp;
	string modelNames;
	int N;//number of subjects
		ifstream fin;
		fin.open(fname.c_str());
		//throw away first three lines 
		getline(fin,stemp);//this is bvars file
			getline(fin,modelNames);//modelnames;
				fin>>stemp>>N;
				for (int i=0; i<N;i++){
					fin>>stemp;//read in subject id
					vnames->push_back(impath+stemp);
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
							
				return modelNames;
}

void write_bvars(string fname,string modelname,Matrix bvars, int numModes,vector<string> vnames){
	ofstream fout;
	
	fout.open(fname.c_str());
	fout<<"this is a bvars file"<<endl; 
	fout<<modelname<<endl;
	fout<<"NumberOfSubjects "<<bvars.Nrows()<<endl;
	
	for (int i=0;i<bvars.Nrows();i++){
			fout<<vnames.at(i)<<" ";
			fout<<numModes<<" ";
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
//	cout<<"calc mode"<<endl;
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
		
		//cout<<"imode "<<i<<" "<<N<<endl;
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
	
	//tets for thalamus
//	volume<float> zvol;
//	copyconvert(*mask, zvol);
//	zvol=0;

	
	float vol=0;

	float min=0, max=0;
	//calculates z-value for all lgive a nice place to initialize EM or could act directly oon intesity anduse range to initialize
	for (int i=bounds[0];i<bounds[1];i++){
		for (int j=bounds[2];j<bounds[3];j++){
			for (int k=bounds[4];k<bounds[5];k++){
				float z=0.0;
				if (mask->value(i,j,k)==(label+100  )){
					z=(ref->value(i,j,k)-mean)/sdev;
//					zvol.value(i,j,k)=(z);
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
					z=(ref->value(i,j,k)-mean)/sdev;
//					zvol.value(i,j,k)=(z);
					if (z>max){ max=z; }
					if (z<min){ min=z; }
				}
			}
		}
	}
	
//	save_volume(zvol,outname.value());

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

//********************************GLM for stats******************************//
ColumnVector GLM_fit(Matrix G, Matrix D, ColumnVector contrast){
	
	//start for only well connditioned design matrices
	Matrix A=G.t()*G;
	Matrix Betas(D.Nrows(),D.Ncols());
	Betas=A.i()*G.t()*D;
	Matrix Mres;
	Mres=D-G*(Betas);
	//calculate residula variance
	ColumnVector avgRes(D.Ncols());
	for (int i=0; i<D.Ncols();i++){
		avgRes.element(i)=((Mres.SubMatrix(1,Mres.Nrows(),i+1,i+1)).t()*(Mres.SubMatrix(1,Mres.Nrows(),i+1,i+1))).AsScalar()/(G.Nrows()-G.Ncols());
	} 
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

float MVGLM_fit(Matrix G, Matrix D, Matrix contrast){
	
	
	//Calculate estimnated values
	Matrix Yhat=G*(G.t()*G).i()*G.t()*D;
	//caluclate E covariance matrix
	Matrix E=D-Yhat;
	E=E.t()*E;
	//calculate H, the sum-square /cross square porduct for hypthosis test
	Matrix YhatH= G*contrast.t()*(contrast*G.t()*G*contrast.t()).i()*contrast*G.t()*D;
	Matrix H=D-YhatH;
	//nto effecient but easy to convert to other statistics
	H=H.t()*H-E;
	
	//calculate Wilks Lambda
	int g=D.Ncols();
	float F=0, df2=0,df1=0;
	int p=G.Ncols();//number of dependant
		int N=D.Nrows();//total sampel size
			if (useWilks.value()){
				float wilk=(E.Determinant()/(H+E).Determinant());
				float a=N-g-((p-g+2)/2.0);
				float b=sqrt((p*p*(g-1)*(g-1)-4)/(p*p+(g-1)*(g-1)-5));
				float c=(p*(g-1)-2)/2;
				//F aprox to wilks
				F=( (1-powf(wilk,1/b)) / powf(wilk,1/b)  )*((a*b-c)/(p*(g-1)));
				df2=(a*b-c);//denominator
					df1=(p*(g-1));//numeraotr
						cout<<"Wilk's "<<wilk<<" "<<F<<" "<<df1<<" "<<df2<<endl;
			}else{
				
				float pillai=(H*(H+E).i()).Trace();
				
				int s=1;
				if (p<(g-1)) {s=p;}
				else {s=g-1;}
				float t=(abs(p-g-1)-1)/2.0;
				float u=(N-g-p-1)/2.0;				
				F=((2*u+s+1)/(2*t+s+1))*(pillai/(s-pillai));
				df1=s*(2*t+s+1);
				df2=s*(2*u+s+1);
				cout<<"Pillai F "<<pillai<<" "<<F<<" "<<df1<<" "<<df2<<endl;
			}
				
	return F;
}
//******************************EXECUTION FUNCTIONS******************************


void do_work_SingleClean(){
	//this function is working directly on volumes
	volume<float> t1im;
	volume<short> segim;
	int bounds[6]={0,0,0,0,0,0};
	read_volume(t1im,refname.value());
	read_volume(segim,inname.value());
	//FIND LABEL AND BOUNDS FOR EACH IMAGE
	//need to reset lower bounds as well
	bounds[0]=0;
	bounds[2]=0;
	bounds[4]=0;
	bounds[1]=segim.xsize();
	bounds[3]=segim.ysize();
	bounds[5]=segim.zsize();
	int label=findStructLabel(&segim, bounds);
	volume<float> vz1;

	boundaryCorr(&segim, &t1im,label, thresh.value(), bounds);
		save_volume(segim,outname.value());
	
}

//*****************************LINEA TRANSFORM********************************************************//
Matrix rigid_linear_xfm(Matrix Data,ColumnVector meanm, Mesh mesh, bool writeToFile){
	//ColumnVector avgM(sub.Ncols());
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
		//cout<<" i "<<i<<endl;
		for (int j=0;j<Data.Nrows();j=j+3){
			sx+=Data.element(j,i);
			sy+=Data.element(j+1,i);
			sz+=Data.element(j+2,i);
		}
		vMx.push_back(sx/Npoints);
		vMy.push_back(sy/Npoints);
		vMz.push_back(sz/Npoints);
	}
	vector< Matrix > vR;
	vector< float > vscale;
	
	for (int subject=0;subject<Data.Ncols();subject++){
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
		
		//*************This includes scale calculation ***********
		float scale=1.0;
		if (useScale.value()){
			//***********calculate scale*************//
			//Data is right (in Horn), Ref is left
			//caluclate length vectors and sum
			float sumr=0, suml=0;
			for (int i=0;i<DataDM.Nrows();i++){
				suml+=(DataDM.element(0,i)*DataDM.element(0,i) + DataDM.element(1,i)*DataDM.element(1,i) +DataDM.element(0,i)*DataDM.element(1,i) );
				sumr+=(RefDM.element(0,i)*RefDM.element(0,i) + RefDM.element(1,i)*RefDM.element(1,i) +RefDM.element(0,i)*RefDM.element(1,i) );
			}
			scale=sqrt(sumr/suml);
		}
		vscale.push_back(scale);
		//**********end of scale calculation**************
		//***********calculate rotattions*************//
		//cout<<"reshaped matrices"<<endl;
		Matrix M=RefDM*(DataDM.t());
		Matrix U;
		DiagonalMatrix D;
		SVD(M.t()*M,D,U);
		//M should always be a 3x3 matrix
		for (int i=0;i<D.Nrows();i++){
			D.element(i)=1/sqrt(D.element(i));
		}
		
		Matrix R(3,3);
		R=M*(U*D*U.t());
		vR.push_back(R);
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
									 //difference between 6dof and 7dof
			if (useScale.value()){
				//cout<<"Scale "<<vscale.at(subject)<<endl;
				Reg=vscale.at(subject)*vR.at(subject)*DataRS+(RefMean-vscale.at(subject)*vR.at(subject)*DataMean);
			}else{
				Reg=vR.at(subject)*DataRS+(RefMean-vR.at(subject)*DataMean);
				//Reg=R*DataRS+(RefMean-R*DataMean);
				
			}
			
			for (int i=0;i<Reg.Ncols();i++){
				DataNew.element(3*i,subject)=Reg.element(0,i);
				DataNew.element(3*i+1,subject)=Reg.element(1,i);
				DataNew.element(3*i+2,subject)=Reg.element(2,i);
			}
		
			
			Mesh m=mesh;
			int count=0;
			for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
				(*i)->_update_coord.X=DataNew.element(count,subject);
				(*i)->_update_coord.Y=DataNew.element(count+1,subject);
				(*i)->_update_coord.Z=DataNew.element(count+2,subject);
				count+=3;
			}	
			m.update();
				
			string snum;
			stringstream ssnum;
			ssnum<<subject;
			ssnum>>snum;
			m.save(snum+"reg.vtk",3);
	}
	
	return DataNew;			

}


Matrix recon_meshesMNI( shapeModel* model1, Matrix bvars, ColumnVector* meanm, Mesh * meshout,vector<Mesh>* vMeshes){
	vMeshes->clear();
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
		//int M=model1->getNumberOfSubjects();
		int Tpts=model1->getTotalNumberOfPoints();
			Matrix MeshVerts(3*Tpts,bvars.Ncols());
		//this is different than mnumber of subjects which were used to create model
		//int numSubs=bvars.Ncols();
		//this loads all mesh point into a matrix to do stats on....
		//cout<<"generate and load vertices into matrix "<<endl;
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
				vMeshes->push_back(m);
				int count=0;
				for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
					MeshVerts.element(3*cumnum+count,j)=(*i)->get_coord().X;
					MeshVerts.element(3*cumnum+count+1,j)=(*i)->get_coord().Y;
					MeshVerts.element(3*cumnum+count+2,j)=(*i)->get_coord().Z;
					count+=3;
				}	
				cumnum+=model1->getNumberOfPoints(sh);	
			}	
		}
	return MeshVerts;

}

Matrix recon_meshesNative( string modelname, Matrix bvars, ColumnVector* meanm, Mesh * meshout, vector<string> subjectnames, vector<string> flirtmats,vector<Mesh>* vMeshes){
	vMeshes->clear();
	
	shapeModel* model1= new shapeModel();
	model1->setImageParameters(refXsize,refYsize, refZsize,refXdim, refYdim, refZdim);
	cout<<"load model "<<modelname<<endl;
	model1->load_bmv_binaryInfo(modelname,1);
	model1->load_bmv_binary(modelname,1);
	
	//want to return mean mesh in vector form
	meanm->ReSize(model1->getTotalNumberOfPoints()*3);
	{
		
		int count=0;
		int cumnum=0; 
		for (int sh=0; sh<model1->getNumberOfShapes();sh++){
			//dieferrence between origin specification
			Mesh m=model1->getTranslatedMesh(sh);
//			m.save("targetMesh.vtk",3);
			//Mesh m=model1->getShapeMesh(sh);
			*meshout=m;
			for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
				meanm->element(3*cumnum+count)=(*i)->get_coord().X;
				meanm->element(3*cumnum+count+1)=(*i)->get_coord().Y;
				meanm->element(3*cumnum+count+2)=(*i)->get_coord().Z;
				//				cout<<"coord "<<meanm->element(3*cumnum+count)<<" "<<meanm->element(3*cumnum+count+1)<<" "<<meanm->element(3*cumnum+count+2)<<endl;

				count+=3;
			}
			cumnum+=model1->getNumberOfPoints(sh);	
		}
		
	}
	
	//need number of subjects (to detrmine number of modes)
	int Tpts=model1->getTotalNumberOfPoints();
	Matrix MeshVerts(3*Tpts,bvars.Ncols());
	//this is different than mnumber of subjects which were used to create model
	//this loads all mesh point into a matrix to do stats on....
//	cout<<"generate and load vertices into matrix "<<endl;
	for (int j=0;j<bvars.Ncols();j++){
		volume<float> t1im;
		//for each subject
		vector<float> vars;
		for (int i=0; i<bvars.Nrows();i++){
			vars.push_back(bvars.element(i,j));
		}
		//keep track of number of points preceding
		int cumnum=0;
		for (int sh=0; sh<model1->getNumberOfShapes();sh++){
			Mesh m=model1->getDeformedMesh(vars,sh,static_cast<int>(vars.size()));	
			model1->meshReg(&m, flirtmats.at(j));
			vMeshes->push_back(m);
			int count=0;
			for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
				MeshVerts.element(3*cumnum+count,j)=(*i)->get_coord().X;
				MeshVerts.element(3*cumnum+count+1,j)=(*i)->get_coord().Y;
				MeshVerts.element(3*cumnum+count+2,j)=(*i)->get_coord().Z;
				count+=3;
			}	
			cumnum+=model1->getNumberOfPoints(sh);	
		}	
	}
	return MeshVerts;
	
}

Matrix deMeanMatrix(Matrix M){
	//demean rows
	Matrix Mnew(M.Nrows(),M.Ncols());
	for (int i=0; i<M.Nrows();i++){
		float sum=0;
		for (int j=0; j<M.Ncols();j++){
			sum+=M.element(i,j); 
		}
		//sum becomes mean
		sum/=M.Ncols();
		for (int j=0; j<M.Ncols();j++){
			Mnew.element(i,j)=M.element(i,j)-sum; 
		}
	}
	
	return Mnew;
}


void do_work_bvars(){
	//**********read in bvars and models and lfirt matrices***************//
	string mname;
	mname=read_bvars_ModelName(inname.value() );	
	
	//load model 
	shapeModel* model1 = new shapeModel;
	model1->setImageParameters(refXsize,refYsize, refZsize,refXdim, refYdim, refZdim);
	model1->load_bmv_binaryInfo(mname,1);
	model1->load_bmv_binary(mname,1);
	//need number of subjects (to detrmine number of modes)
	vector<string> subjectnames;
	Matrix bvars;
	vector<int> vN;
	Matrix target;	//target is only used for glm
					//must include a design matrix with bvars to use glm
					//modelname is used to set a path
	
	read_bvars(inname.value(),&bvars,&subjectnames, &vN, pathname.value());
	target=read_vest(designname.value());
				//can filter meshes
	if (usePCAfilter.value()){
		//truncate number of mdoes to recon mehs (smoothing)
		bvars=bvars.SubMatrix(1,numModes.value(),1,bvars.Ncols());
				}
				
				volume<float> t1im;
				volume<short> segim;
				
				vector<Mesh> vMeshes;
				Matrix MeshVerts;//this is used when placing t-stats on a mesh
					
					
					//need flirt matrices
					//load their names into a vector
					ifstream flirtmats;
					vector<string> flirtmatnames;
					
					//**********done reading in bvars and models adn lfirt matrices ***************//
					//****************RECONSTRUCTION AND ALIGNMENT*********************/
					
					//Choose the space in which to reconstruct the meshes
					Mesh modelMeanMesh;
					ColumnVector CVmodelMeanMesh;
					if (useReconNative.value()){
					  	flirtmats.open(flirtmatsname.value().c_str());
						
						for (unsigned int i =0; i<subjectnames.size();i++){
							string stemp;
							flirtmats>>stemp;
							flirtmatnames.push_back(stemp); ///inversion fo flirtmatrix is handled in shape model function when a string is input
						}
						
						//Reconstruct in native spoace of the image (it recons the mni then applies flirt matrix)
						MeshVerts=recon_meshesNative( mname, bvars, &CVmodelMeanMesh, &modelMeanMesh,subjectnames, flirtmatnames, &vMeshes);
					}
					else if(useReconMNI.value()){
						//flirt matrix is not applied
						MeshVerts=recon_meshesMNI( model1, bvars, &CVmodelMeanMesh, &modelMeanMesh,&vMeshes);
					}else{
						//you must choose a method of reconstruction
						cerr<<"choose a mesh reconstruction method"<<endl;
						return;
					}
					//Added in any realignment of meshes here
					if (useRigidAlign.value()){
						//use a least-squares alignment of the meshes (in this case to the mean as defined by the model
						MeshVerts=rigid_linear_xfm(MeshVerts,CVmodelMeanMesh, modelMeanMesh, false);
					}
					//****************END RECONSTRUCTION AND ALIGNMENT*********************/
					//****************DEMEAN DESIGN MATRIX AND ADD ONE COLUMS*********************/
					
					
					
					//if the design has nto been demeaned and you are dooing discrimiant analysis
					//then perform demean of matrix
					
					
						bool isOne=true;
						//checks for mean Column as first column, if it decides there is a column of ones (first column) then 
						//it will not demean the design matrix
						for (int i=0;i<target.Nrows();i++){
							if (target.element(i,0)!=1){
								isOne=false;
							}
						}
						
						
						if(!isOne){
							//create demena deisgn, assume a nomean column at start
							Matrix targTemp(target.Nrows(),target.Ncols()+1);
							for (int i=0;i<targTemp.Ncols();i++){
								if (i==0){
									for (int j=0;j<target.Nrows();j++){
										targTemp.element(j,i)=1;
									}
								}else{
									float mean=0;
									for (int j=0;j<target.Nrows();j++){
										mean+=target.element(j,i-1);
									}	
									mean/=target.Nrows();
									for (int j=0;j<target.Nrows();j++){
										targTemp.element(j,i)=target.element(j,i-1)-mean;
									}
								}
							}
							
							target=targTemp;
						}
						//this displays the demean design matrix
						cout<<"new design matrix"<<endl;
						for (int j=0;j<target.Nrows();j++){
							for (int i=0;i<target.Ncols();i++){	
								cout<<target.element(j,i)<<" ";
							}	
							cout<<endl;
						}
					
					//cout<<"done processing design"<<endl;
					if (!vertexAnalysis.value()){
						
						//calculate volumes
						vector<float> vnorm;
						
						if (useNorm.value()){
							ifstream fnorm;
							fnorm.open(normname.value().c_str());
							for (int subject=0;subject<bvars.Ncols();subject++){
								float temp;
								fnorm>>temp;
								vnorm.push_back(temp);							
							}
						}
						
						Matrix Volumes(target.Nrows(),1); //this is only used for GLM
						ofstream fvol_all;
						fvol_all.open((outname.value()+".vols").c_str());
						for (int subject=0;subject<bvars.Ncols();subject++){
							//this conditions allows you to loads volumes
							if (useVolumes.value()){
								ifstream volIn;
								volIn.open(refname.value().c_str());
								for (int subject=0;subject<bvars.Ncols();subject++){
									float voltemp;
									string stemp;
									volIn>>stemp>>voltemp>>voltemp;
									cout<<"read vol "<<voltemp<<endl;
									
									if (useNorm.value()){
										Volumes.element(subject,0)=voltemp*vnorm.at(subject);
										
									}else{
										Volumes.element(subject,0)=voltemp;
									}
								}
								break;
							}
							//if volumes are load the loop is broken
							read_volume(t1im,pathname.value()+subjectnames.at(subject));
							Mesh m = vMeshes.at(subject);
							
							//fill mesh
					    	int bounds[6]={0,0,0,0,0,0};
							segim=make_mask_from_meshInOut(t1im,m,model1->getLabel(0),bounds);
							stringstream sstemp;
							string outnamess;
							sstemp<<subject;
							sstemp>>outnamess;
							
							float voltemp;
							volume<short> segimB;
							segimB=segim;
							voltemp=boundaryCorr(&segim, &t1im,model1->getLabel(0), thresh.value(),bounds);
							stringstream sstemp2;
							sstemp2<<model1->getLabel(0);
							string lbst;
							sstemp2>>lbst;
							
							save_volume(segim,subjectnames.at(subject)+"FIRSTbcorr_lb"+lbst);
							fvol_all<<subjectnames.at(subject)<<" "<<voltemp<<" "<<voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim()<<endl;
							//save volume to a matrix if you wish to use in GLM
							if (useNorm.value()){
								Volumes.element(subject,0)=voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim()*vnorm.at(subject);
								
							}else{
								Volumes.element(subject,0)=voltemp*t1im.xdim()*t1im.ydim()*t1im.zdim();
							}
							
							
							
							
						}
						//now that volumes are load perform GLM
						cout<<"voluems loaded"<<endl;
						ColumnVector tstats;
						ColumnVector contrast(target.Ncols());
						
						for (int EV=1;EV<target.Ncols();EV++){
							
							stringstream ev2st;
							ev2st<<EV;
							string evnum;
							ev2st>>evnum;
							
							for (int i=0;i<contrast.Nrows();i++){
								if(i==EV){
									contrast.element(i)=1;
								}else{
									contrast.element(i)=0;
								}
								cout<<contrast.element(i)<<" ";
							}
							
							tstats=GLM_fit(target, Volumes,contrast);
							ofstream fTs;
							string ftname=outname.value()+evnum+".tstat";
							fTs.open(ftname.c_str());
							for (int i =0 ; i <tstats.Nrows();i++){
								fTs<<tstats.element(i)<<endl;
							}
						}
						
					}else{
						ColumnVector CVnorm(target.Nrows());
						// if chosen can include a normalization EV (i.e. control for size) ...useful for vertex shape statistics
						if (useNorm.value()){
							ifstream normIn;
							float normtemp;
							normIn.open(normname.value().c_str());
							for (int subject=0;subject<target.Nrows();subject++){
								normIn>>normtemp;
								CVnorm.element(subject)=normtemp;
							}
							target=target | CVnorm ;
						}
						
						
						
						//****************END DEMEAN DESIGN MATRIX AND ADD ONE COLUMS*********************/
						
						Matrix contrast(target.Ncols()-1,target.Ncols());
						int EVmin=1;//ignore mean column first EV to examine
							//int EVmax=1;//EV to include
							int	EVmax=target.Ncols();

								
								for (int EV=EVmin;EV<EVmax;EV++){
									
									
									//append the EV number to the output
									stringstream ev2st;
									ev2st<<EV;
									string evnum;
									ev2st>>evnum;
									
								
									
									
									//update average mesh
									//these vector store tstats to plot on mesh
									vector<float> tstatsx;
									vector<float> tstatsy;
									vector<float> tstatsz;
									
									
									Mesh m=vMeshes.at(0);//mesh.at(0) is used to determien topology/number of verteices
										int count=0;
										for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ )
										{ 
											float meanx1=0,meany1=0,meanz1=0,meanx2=0,meany2=0, meanz2=0;
											// for (int i=0;i<MeshVerts.Nrows();i=i+3){
											int n1=0,n2=0;
											for (int j=0;j<MeshVerts.Ncols();j++){
												if (target.element(j,EVmin) <= 0 ){ //handles demeaned data
													meanx1+=MeshVerts.element(count,j);
													meany1+=MeshVerts.element(count+1,j);
													meanz1+=MeshVerts.element(count+2,j);
													n1++;
												}else{
													meanx2+=MeshVerts.element(count,j);
													meany2+=MeshVerts.element(count+1,j);
													meanz2+=MeshVerts.element(count+2,j);
													n2++;
												}
											}
											// }
											meanx1/=n1;
											meany1/=n1;
											meanz1/=n1;
											
											meanx2/=n2;
											meany2/=n2;
											meanz2/=n2;
											
											//as binary classification was assumed we just finish clauclating the mean for each group
											
											//save vectors point from group1 mean to group2 mean, size 1 (they will be scaled to reflect
											//the test statistic)
											float vecx=meanx2-meanx1;
											float vecy=meany2-meany1;
											float vecz=meanz2-meanz1;
											//	float norm=sqrt(vecx*vecx+vecy*vecy+vecz*vecz);
											
											//scalign occurs later
											tstatsx.push_back(vecx);///norm);
												tstatsy.push_back(vecy);///norm);
													tstatsz.push_back(vecz);///norm);
														
														//for each point calculate mean vertex for each group
														(*i)->_update_coord = Pt(meanx1,meany1,meanz1);
														count+=3;
										}
										m.update();
										model1->setShapeMesh(0,m);
										vector<float> scalarsT;
										
										//this provides the option to do stats on the vertices 
										//create F test contrast matrix (testing each EV separately
											count=0;
											for (int j=0;j<contrast.Ncols();j++){
												if (j!=EV){
													for (int i=0;i<contrast.Ncols();i++){
														if(i==j){
															//if(i==EV){
															contrast.element(count,i)=1;
															}else{
																contrast.element(count,i)=0;
															}
														
														cout<<contrast.element(count,i)<<" ";
														}
													cout<<endl;
													count++;
													}
												}
											
											//use multivariate test on each vertex
											for (int i=0;i<MeshVerts.Nrows();i=i+3){
												//use multivariate multiple regression on each veretx 
												float F=MVGLM_fit(target, MeshVerts.SubMatrix(i+1,i+3,1,MeshVerts.Ncols()).t(),contrast);
												//	tstatsx.at(i/3)*=wilkL;
												//	tstatsy.at(i/3)*=wilkL;
												//	tstatsz.at(i/3)*=wilkL;	
												//this may 
												scalarsT.push_back(F);
												
											}
																												
										for (int sh=0; sh<model1->getNumberOfShapes();sh++){
											model1->setShapeTstatX(sh,tstatsx);
											model1->setShapeTstatY(sh,tstatsy);
											model1->setShapeTstatZ(sh,tstatsz);
											model1->setShapeScalars(sh,scalarsT);
										}
										//save the mesg with vectors
										model1->save(outname.value()+evnum+".vtk",5,0);
									}
					}
					}		


void do_work_meshToVol(){
 
	volume<float> ref;
	volume<short> segim;
	read_volume(ref,inname.value());

	shapeModel* mesh1 = new shapeModel;
	mesh1->setImageParameters(ref.xsize(), ref.ysize(), ref.zsize(), ref.xdim(), ref.ydim(), ref.zdim());
	mesh1->load_vtk(modelname.value(),1);
	Mesh m1;
	if (centreOrigin.value()){
	  m1=mesh1->getTranslatedMesh(0);
	}else{
	  m1=mesh1->getShapeMesh(0);
	}
	int bounds[6]={0,0,0,0,0,0};
	int label=meshLabel.value();	
	
	segim=make_mask_from_meshInOut(ref,m1,label,bounds);
	
	for (int i=bounds[0];i<bounds[1];i++){
		for (int j=bounds[2];j<bounds[3];j++){
			for (int k=bounds[4];k<bounds[5];k++){
				
				if (segim.value(i,j,k)==(label+100  )){
					if (bcd.value()){
						if (ref.value(i,j,k)==label){
							segim.value(i,j,k)=label;
						}else{
							segim.value(i,j,k)=0;
						}
					}else{
						segim.value(i,j,k)=label;
					}
					
				}
			}
		}
	}

	save_volume(segim,outname.value());

}
void do_work_overlap(){
  //this function is working directly on volumes
  

  volume<short> gold;
  volume<short> segim;

  read_volume(gold,refname.value());
  read_volume(segim,inname.value());
  overlaps(segim,gold);
  
}


int main(int argc,char *argv[])
{
	
	Tracer tr("main");
	OptionParser options(title, examples);
	
	try {
		// must include all wanted options here (the order determines how
		//  the help message is printed)
		options.add(inname);
		options.add(normname);
		options.add(refname);
		options.add(pathname);
		options.add(flirtmatsname);
		options.add(useScale);
		options.add(overlap);
		options.add(univariateT);
		options.add(bcd);
		options.add(modelname);
		options.add(outname);
		options.add(thresh);
		options.add(meshLabel);
		options.add(usebvars);
		options.add(useReconMNI);
		options.add(vertexAnalysis);
		options.add(useReconNative);
		options.add(useRigidAlign);
			options.add(designname);

		options.add(useVolumes);
		options.add(meshToVol);
		options.add(centreOrigin);
		options.add(useWilks);
		options.add(inputprob);
		options.add(verbose);
		options.add(usePCAfilter);
		options.add(numModes);
		options.add(help);
		options.add(singleBoundaryCorr);
		nonoptarg = options.parse_command_line(argc, argv);
		
		// line below stops the program if the help was requested or 
		//  a compulsory option was not set
		if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
		{
			options.usage();
			exit(EXIT_FAILURE);
		}
		
		// Call the local functions
		if  (usebvars.value()){
			do_work_bvars();
			
			//	}else if (imfrombvars.value()){
			//	do_work_bvarsToVols();
		}else if (meshToVol.value()){
			do_work_meshToVol();
		}else if(singleBoundaryCorr.value()){
			do_work_SingleClean();
		}else if (overlap.value()){
		  do_work_overlap();
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

