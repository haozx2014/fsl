/*  Copyright (C) 1999-2005 University of Oxford  */

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
#include "shapeModel.h"
#include "mshape.h"

#include "meshclass/meshclass.h"
#include "newimage/newimageall.h"
#include <sstream>
using namespace std;
using namespace NEWIMAGE;
using namespace mesh;
namespace shapemodel {
	
	shapeModel::shapeModel(){
		xsize = 0;
		ysize = 0;
		zsize = 0;
		xdim = 0;
		ydim = 0;
		zdim = 0;
		numShapes=0;
		sumEigVals=0;
		bounds[0]=1000;
		bounds[1]=1000;
		bounds[2]=1000;
		bounds[3]=1000;
		bounds[4]=1000;
		bounds[5]=1000;
		
	}
	shapeModel::~shapeModel(){
		for (unsigned int i=0; i<shapes.size();i++){
			delete shapes.at(i);
		}
	}
	
	void shapeModel::clear(){
		for (unsigned int i=0; i<shapes.size();i++){
			delete shapes.at(i);
		}
		//dont reset image info
		//	xsize = 0;
		//ysize = 0;
		//zsize = 0;
		//xdim = 0;
		//ydim = 0;
		//zdim = 0;
		numShapes=0;
		sumEigVals=0;
		bounds[0]=1000;
		bounds[1]=1000;
		bounds[2]=1000;
		bounds[3]=1000;
		bounds[4]=1000;
		bounds[5]=1000;
		labels.clear();
		npts.clear();
		eigenVals.clear();
		shapes.clear();
		iprofs.clear();
	}
	
	int shapeModel::load_bmv(string s, int type){
	//int type=0;
	//type 0 = normal shape model
	//type 1 = appearance model
	  //cout<<"Loading Shape Model..."<<endl;
	//  cout<<"load bmv2"<<endl;
		string temp;
		int tempInt,tempInt2;
		float tempFloat, vX,vY,vZ;
		ifstream f;
		
		//read in information (not meshes yet!)
		//	cout<<"reading information"<<endl;
		f.open(s.c_str());
		getline(f,temp);
		getline(f,temp);//used to get rid of second ascii line to avoid "modes"
			while (!f.eof()){
				//getline(f,temp);
				f>>temp;
				
				if (!strcmp(temp.c_str(),"numPoints")){
					f>>numShapes>>numShapes>>temp;
					for (int i=0;i<numShapes;i++){
						f>>tempInt;
						npts.push_back(tempInt);
					}
				}else if (!strcmp(temp.c_str(),"eigenValues")){
					f>>tempInt2>>tempInt2>>temp;
					for (int i=0;i<(tempInt2-1);i++){
						f>>tempFloat;
						eigenVals.push_back(tempFloat);
					}
					f>>sumEigVals;
				}else if (!strcmp(temp.c_str(),"labels")){
					
					f>>tempInt2>>tempInt2>>temp;
					for (int i=0;i<tempInt2;i++){
						f>>tempInt;
						labels.push_back(tempInt);
					}
				}else if (!strcmp(temp.c_str(),"numSubjects")){
					
					f>>tempInt2>>tempInt2>>temp;
					//	for (int i=0;i<tempInt2;i++){
					f>>tempInt;
					numSubjects=tempInt;
					//	cout<<numSubjects<<"NUM"<<endl;
						//	}
				}else if (!strcmp(temp.c_str(),"IPP")){
					
					f>>tempInt2>>tempInt2>>temp;
					for (int i=0;i<tempInt2;i++){
						f>>tempInt;
						ipps.push_back(tempInt);
						//cout<<"IPP "<<tempInt<<endl;
						//need to assign it to shape objects stiill
					}
				}else if ((!strcmp(temp.c_str(),"eigenValuesI"))&&(type==1)){
					f>>tempInt2>>tempInt2>>temp;
					for (int i=0;i<(tempInt2-1);i++){
						f>>tempFloat;
						eigenValsI.push_back(tempFloat);
					}
					f>>sumEigValsI;
					}
			}
			f.close();
			//craete shape objects
			for (int i=0;i<numShapes;i++){
				MShape*  Imshape = new MShape();
				shapes.push_back(Imshape);
				if (type==1){
					shapes.at(i)->setIPP(ipps.at(i));
				}
			}
			
			//read in meshes
			//	cout<<"read mesh data"<<endl;
			ifstream fP,f2;
			
			//use f for points and fP for Polygons
			f2.open(s.c_str());
			fP.open(s.c_str());
			//throw away first 5 lines, assumes float
			getline(f2,temp);
			getline(f2,temp);
			getline(f2,temp);
			getline(f2,temp);
			getline(f2,temp);
			do{
				fP>>temp;
			}while (strcmp(temp.c_str(),"POLYGONS"));
			//tempInt is number of Polygons
			fP>>tempInt>>tempInt2;
			int cumPts=0, Pcount=0;    
			int p0=0, p1=0, p2=0;
			//	cout<<"Loading points and Polygons"<<endl;
			for (int sh=0;sh<numShapes;sh++){//create new mesh for each shape
			//	cout<<"new shape"<<endl;							 //Reading the points
				for (int i=0; i<npts.at(sh); i++)
				{//assign points to mesh
					//cout<<"i: "<<i<<endl;
					float x, y, z;
					f2>>x>>y>>z;
					Mpoint * p = new Mpoint(x, y, z, i);
					shapes.at(sh)->pushPoint(p);
				}
			//	cout<<"reading triangles"<<endl;
				
				//reading the triangles
				//tempInt is the number fo polygons   
				bool fpass=true;
				while (true) //break when an poitn index is too large 
				{
					Pcount++;
					if ((sh>0)&&(fpass)){//a point carried over from previous shape
			//			cout<<"npts"<<npts.at(sh-1)<<" "<<p0<<" "<<p1<<" "<<p2<<endl;
						Triangle * t = new Triangle(shapes.at(sh)->getPoint(p0-npts.at(sh-1)), shapes.at(sh)->getPoint(p1-npts.at(sh-1)), shapes.at(sh)->getPoint(p2-npts.at(sh-1)));
			//			cout<<Pcount<<endl;
						shapes.at(sh)->pushTriangle(t);
			//			cout<<Pcount<<"2"<<endl;
					}
					
					fpass=false;
					fP>>p0>>p0>>p1>>p2;
					//cout<<"P: "<<p0<<" "<<p1<<" "<<p2<<" "<<tempInt<<endl;
					p0-=cumPts;
					p1-=cumPts;
					p2-=cumPts;
					
					if ((p0>=(npts.at(sh)))||(p1>=(npts.at(sh)))||(p2>=(npts.at(sh)))||(Pcount>tempInt)){
						//looks for end of shape or end of all shapes
						break;
					}
					Triangle * t = new Triangle(shapes.at(sh)->getPoint(p0), shapes.at(sh)->getPoint(p1), shapes.at(sh)->getPoint(p2));
					shapes.at(sh)->pushTriangle(t);
					
				}
				cumPts+=npts.at(sh);
			}
			fP.close();
			f2.close();
			//completed loading the shapes into meshes (points and polygon)
			//time to read in Vectors
			//first create modeVectors
			//	cout<<"creating mode vectors"<<endl;
			
			ifstream fV;
			fV.open(s.c_str());
			getline(fV,temp);
			getline(fV,temp);//used to get rid of second ascii line to avoid "modes"
				while (!fV.eof()){
					//cout<<"?????"<<endl;
					fV>>temp;
			//	cout<<"problem "<<temp<<endl;
				//cout<<"tempM "<<temp<<endl;
					if (!strcmp(temp.substr(0,4).c_str(),"mode")){
				//		cout<<"MODE FOUND "<<temp<<endl;
						//read in eigenvector (mode)
						string temp2;
						fV>>tempInt>>tempInt>>temp2;
						for (int i=0;i<numShapes;i++){
							vector<Vec> V;
							for (int j=0;j<npts.at(i);j++){
								fV>>vX>>vY>>vZ;
								Vec vtmp(vX,vY,vZ);
								V.push_back(vtmp);
							}
							shapes.at(i)->addModeVector(V);
						}
						
					}else if((!strcmp(temp.substr(0,5).c_str(),"Imean"))&&(type==1)){
						int ipp; 
						string temp2;
						fV>>ipp>>tempInt>>temp2;
						//cout<<"imean"<<endl;
						for (int i=0;i<numShapes;i++){
							vector<float> Vi;
							for (int j=0;j<npts.at(i);j++){
								
								for (int k=0;k<ipp;k++){
									fV>>tempFloat;
									Vi.push_back(tempFloat);
								}
							}
							//fV>>temp2>>temp2;
							shapes.at(i)->setIMean(Vi);
						}
					}else if((!strcmp(temp.substr(0,5).c_str(),"Bmean"))&&(type==1)){
				//	cout<<"Bmean found"<<endl;
						int ipp; 
						string temp2;
						fV>>ipp>>tempInt>>temp2;
						//cout<<"imean"<<endl;
						for (int i=0;i<numShapes;i++){
							vector<float> Vi;
							for (int j=0;j<4;j++){
								fV>>tempFloat;
								Vi.push_back(tempFloat);
							}
							
							//fV>>temp2>>temp2;
							shapes.at(i)->setBMean(Vi);
						}
					}else if ((!strcmp(temp.substr(0,5).c_str(),"Imode"))&&(type==1)){
					//cout<<"FOUND IMODE "<<temp<<endl;
						int ipp; 
						string temp2;
						fV>>ipp>>tempInt>>temp2;
					//	cout<<ipp<<" "<<tempInt<<" "<<temp2<<endl;
						
						for (int i=0;i<numShapes;i++){
							vector<float> Vi;
						//	cout<<"i "<<i<<" "<<npts.at(i)<<endl;
							for (int j=0;j<npts.at(i);j++){
								//cout<<"j "<<j<<endl;
								
								for (int k=0;k<ipp;k++){
									fV>>tempFloat;
									//cout<<"k "<<k<<endl;
									//if (j==1000){
									//	cout<<"ihxdhdh "<<tempFloat<<endl;
									//}
									Vi.push_back(tempFloat);
								}
							//	fV>>temp2;
								//if (j==0){
								//	cout<<"temp "<<temp2<<endl;
							//	}
							}
							//fV>>temp2>>temp2;
						//	cout<<"temp "<<temp2<<endl;
							shapes.at(i)->addIModeVector(Vi);
						//	shapes.at(0)->getIModeVector(0);
						}
					}else if ((!strcmp(temp.substr(0,5).c_str(),"Bmode"))&&(type==1)){
						//	cout<<"FOUND IMODE "<<temp<<endl;
						int ipp; 
						string temp2;
						fV>>ipp>>tempInt>>temp2;
						
						for (int i=0;i<numShapes;i++){
							vector<float> Vi;
							for (int j=0;j<4;j++){
								fV>>tempFloat;
								Vi.push_back(tempFloat);
							}
							
							
							shapes.at(i)->addBModeVector(Vi);
						}
					}else if ((!strcmp(temp.substr(0,5).c_str(),"iCond"))&&(type==1)){
					//cout<<"iconds found"<<endl;
						int sh;
						sh = atoi(temp.substr(9,1).c_str());
						if (!strcmp(temp.substr(0,9).c_str(),"iCondPrec")){
							string temp2;
							int tempInt2;
							//cout<<"enter icondprec"<<endl;
							fV>>tempInt2>>tempInt>>temp2;

							vector< vector<float> > prectemp;							
							//cout<<"number of subjects/cols in u"<<tempInt2<<" "<<tempInt<<endl;
							
							for (int i=0;i<tempInt;i++){
								vector<float> col;
								for (int j=0;j<tempInt2;j++){
									fV>>tempFloat;
									col.push_back(tempFloat);
								//	cout<<"tempFLoat "<<i<<" "<<j<<" "<<tempFloat<<endl;
								}
								prectemp.push_back(col);
							}
							shapes.at(sh)->setICondPrec(prectemp);
							//cout<<"setting icondprec"<<endl;
						}else if (!strcmp(temp.substr(0,9).c_str(),"iCondEigs")){
							//cout<<"enter icondeigs"<<endl;
							string temp2;
							vector<float> eigstemp;
							fV>>tempInt>>tempInt>>temp2;
						//	cout<<"number of subjects/cols in u "<<tempInt<<endl;
							for (int j=0;j<tempInt;j++){
								fV>>tempFloat;
								eigstemp.push_back(tempFloat);
							}
							shapes.at(sh)->setICondEigs(eigstemp);
						}
					//	cout<<"iconds loaded"<<endl;
					}else if ((!strcmp(temp.substr(0,5).c_str(),"bCond"))&&(type==1)){
					//cout<<"found bconds "<<temp<<endl;
						int sh;
						sh = atoi(temp.substr(9,1).c_str());
						if (!strcmp(temp.substr(0,9).c_str(),"bCondPrec")){
					//	cout<<"am i in?"<<endl;
							string temp2;
							int tempInt2;
							//cout<<"enter icondprec"<<endl;
							fV>>tempInt2>>tempInt>>temp2;

							vector< vector<float> > prectemp;							
							//cout<<"number of subjects/cols in u"<<tempInt2<<" "<<tempInt<<endl;
							
							for (int i=0;i<tempInt;i++){
								vector<float> col;
								for (int j=0;j<tempInt2;j++){
									fV>>tempFloat;
									col.push_back(tempFloat);
									//cout<<"tempFLoat "<<i<<" "<<j<<" "<<tempFloat<<endl;
								}
								prectemp.push_back(col);
							}
							shapes.at(sh)->setBCondPrec(prectemp);
							//cout<<"setting icondprec"<<endl;
						}else if (!strcmp(temp.substr(0,9).c_str(),"bCondEigs")){
							//cout<<"enter bcondeigs"<<endl;
							string temp2;
							vector<float> eigstemp;
							fV>>tempInt>>tempInt>>temp2;
						//	cout<<"number of subjects/cols in u "<<tempInt<<endl;
							for (int j=0;j<tempInt;j++){
								fV>>tempFloat;
								eigstemp.push_back(tempFloat);
							}
							shapes.at(sh)->setBCondEigs(eigstemp);
						}
					}else if ((!strcmp(temp.substr(0,9).c_str(),"ErrPriors"))){
					int sh;
						sh = atoi(temp.substr(9,1).c_str());
						string temp2;
						fV>>tempInt>>tempInt>>temp2;
						//cout<<"enter errpiros"<<endl;
							vector<float> verrs;
							for (int j=0;j<2;j++){ //this allows for error in intensity and shape, for sh=0 its is the predictor error
								fV>>tempFloat;
								verrs.push_back(tempFloat);
							//		cout<<"SHAPPPPPPPPEEEEE "<<sh<<" "<<tempFloat<<endl;
							
								}
							shapes.at(sh)->setErrs(verrs);
						}
	//cout<<"problem "<<temp<<endl;
				}

			//	cout<<"Shape Model has finished loading"<<endl;
				return 1;
	}
	int shapeModel::load_vtk(string s, int type){
		clear();
		//type 0 = normal mesh
		//type 1 = intesity profiles joined
		//type 2 = intesnity profiles + planar estimates
		//this supports the loading of a single structure as a shapeModel with 1 mode vector valued at zero
		//read in meshes
		//cout<<"Loading vtk mesh..."<<endl;
		string temp;
		int numPts, tempInt,tempInt2;
		
		
	//		cout<<"read mesh data"<<endl;
			ifstream fP,f2;
			
			//use f for points and fP for Polygons
			f2.open(s.c_str());
			fP.open(s.c_str());
			//throw away first 5 lines, assumes float
			getline(f2,temp);
			getline(f2,temp);
			getline(f2,temp);
			getline(f2,temp);
			f2>>temp>>numPts>>temp;
			getline(f2,temp);
			//cout<<"numPts "<<numPts<<" "<<temp<<endl;
			npts.push_back(numPts);
			do{
				fP>>temp;
			}while (strcmp(temp.c_str(),"POLYGONS"));
			//tempInt is number of Polygons
			fP>>tempInt>>tempInt2;
			int cumPts=0, Pcount=0;    
			int p0=0, p1=0, p2=0;
		//	cout<<"Loading points and Polygons"<<endl;
			
			//create single shape
			MShape*  Imshape = new MShape();
			shapes.push_back(Imshape);
			numShapes=1; 
			for (int sh=0;sh<numShapes;sh++){//create new mesh for each shape
			//	cout<<"new shape"<<endl;							 //Reading the points
				for (int i=0; i<numPts; i++)
				{//assign points to mesh
				//	cout<<"i: "<<i<<endl;
					float x, y, z;
					f2>>x>>y>>z;
					Mpoint * p = new Mpoint(x, y, z, i);
					shapes.at(sh)->pushPoint(p);
				}
		//		cout<<"reading triangles"<<endl;
				
				//reading the triangles
				//tempInt is the number fo polygons   
				bool fpass=true;
				while (true) //break when an poitn index is too large 
				{
					Pcount++;
									
					fpass=false;
					fP>>p0>>p0>>p1>>p2;
					//cout<<"P: "<<p0<<" "<<p1<<" "<<p2<<" "<<tempInt<<endl;
					p0-=cumPts;
					p1-=cumPts;
					p2-=cumPts;
					
					if ((p0>=(npts.at(sh)))||(p1>=(npts.at(sh)))||(p2>=(npts.at(sh)))||(Pcount>tempInt)){
						//looks for end of shape or end of all shapes
						break;
					}
					Triangle * t = new Triangle(shapes.at(sh)->getPoint(p0), shapes.at(sh)->getPoint(p1), shapes.at(sh)->getPoint(p2));
					shapes.at(sh)->pushTriangle(t);
					
				}
				cumPts+=numPts;
			}
			
			
		//	cout<<"Shape Model has finished loading"<<endl;

			
			
			
			fP.close();
			f2.close();
			//completed loading the shapes into meshes (points and polygon)
			//time to read in Vectors
			//first create modeVectors
			//no mode vectors because reading in meshes
			//cout<<"creating mode vectors"<<endl;
			Vec vtmp(0,0,0);
			vector<Vec> V;
			for (int j=0;j<numPts;j++){
				V.push_back(vtmp);
			}
			//cout<<"addshapevector"<<endl;
			shapes.at(0)->addModeVector(V);
			
			eigenVals.push_back(0);
			sumEigVals=0;
			labels.push_back(1);
			if ((type==1)|(type==2)){
				//cout<<"Read intensity profiles"<<endl;
				int ipp;//intensities per point
				ifstream fI;
				fI.open(s.c_str());
				getline(fI,temp);
			//	cout<<temp;
				getline(fI,temp);//used to get rid of second ascii line to avoid "modes"
				//cout<<temp;
				while (!fI.eof()){
					fI>>temp;
					
					//cout<<temp<<endl;
					if (!strcmp(temp.c_str(),"FIELD")){
					//cout<<temp<<endl;
					//	cout<<"found FIELD "<<endl;
						getline(fI,temp);
						fI>>temp;
					//	cout<<temp<<endl;
						fI>>ipp;
						getline(fI,temp);
						//cout<<"There are "<<ipp<<" intensity points"<<endl;
						break;
					}
				}
				float tempI;
				//fI>>tempI;
				//cout<<"IPPPPPP "<<ipp<<endl;
				for (int i=0;i<getTotalNumberOfPoints();i++){
					vector<float> vtemp;
					for (int j=0;j<ipp;j++){
						fI>>tempI;
					//	if (j==0){
					//	cout<<"intenread "<<tempI<<endl;
					//	}
						vtemp.push_back(tempI);
					//	if (j==0){
					//	cout<<"intenread2 "<<vtemp.at(0)<<endl;
					//	}
					}		
					iprofs.push_back(vtemp);		
					//cout<<"tempI "<<tempI<<endl;
				}	
				if (type==2){
					getline(fI,temp);
					getline(fI,temp);
					vector<float> best;
					
					for (int i=0;i<4;i++){
						float tempb;
						fI>>tempb;
						best.push_back(tempb);
					}
					setPlanarParameters(best);
				}
			}
			
			
			return 1;
	
	
	}


	void shapeModel::save(string s, int type, int numModes){
		//numModes currently doesn't do anything
		//type =0 for normal shape model
		//type =1 for appearance mesh
		//type 2 = appearance model
		//type 3 = distance/appearance model
		//type 4= appearance model with planar estimates
		//type 5= mesh and t-stats
		ofstream fshape;  
		fshape.open(s.c_str());
		
		#ifdef PPC64
    int n=0;
#endif
		//calculate total number of points
		int ptTot=0;
		for (int i=0; i<getNumberOfShapes();i++){
			ptTot+=getNumberOfPoints(i);
		}
		if (type==3){
			fshape<<"# vtk DataFile Version 3.0"<<endl<<"IBIM modes of variation file"<<endl<<"ASCII"<<endl<<"DATASET STRUCTURED_POINTS "<<endl<<"DIMENSIONS "<<xsize<<" "<<ysize<<" "<<zsize<<endl;
			fshape<<"ORIGIN "<<origin.at(0)<<" "<<origin.at(1)<<origin.at(2)<<endl<<"SPACING "<<xdim<<" "<<ydim<<" "<<zdim<<endl;
		}else{
			fshape<<"# vtk DataFile Version 3.0"<<endl<<"IBIM modes of variation file"<<endl<<"ASCII"<<endl<<"DATASET POLYDATA"<<endl<<"POINTS "<<ptTot<<"  float"<<endl;
		}
		
		
		if (type!=3){
			//writing points to vtk file
			for (int i=0; i<getNumberOfShapes(); i++) { 
				//	cout<<"wirting point for shape i= "<<i<<endl;
				//	cout<<getNumberOfPoints(i)<<endl;
				for (int j=0; j<getNumberOfPoints(i);j++)
				{    
					fshape<<shapes.at(i)->getPoint(j)->get_coord().X<<" "<<shapes.at(i)->getPoint(j)->get_coord().Y<<" "<<shapes.at(i)->getPoint(j)->get_coord().Z<<endl;
					#ifdef PPC64
						if ((n++ % 20) == 0) fshape.flush();
					#endif
				}
			}
			
			//writing connectuions to vtk file
			int cumNumPoints=0;
			vector<int> Polygons;
			for (int i=0;i<getNumberOfShapes();i++){
				
				Mesh mtemp=shapes.at(i)->getMesh();
				for ( list<Triangle*>::const_iterator t=mtemp._triangles.begin(); t!=mtemp._triangles.end(); t++) {
					
					Polygons.push_back((*t)->get_vertice(0)->get_no()+cumNumPoints);
					Polygons.push_back((*t)->get_vertice(1)->get_no()+cumNumPoints);
					Polygons.push_back((*t)->get_vertice(2)->get_no()+cumNumPoints);	  
				}
				cumNumPoints+=getNumberOfPoints(i);
			}
			
			fshape<<"POLYGONS "<<Polygons.size()/3<<" "<<Polygons.size()/3*4<<endl;
			for (unsigned int i=0;i<Polygons.size();i=i+3){
				fshape<<"3 "<<Polygons.at(i)<<" ";
				fshape<<Polygons.at(i+1)<<" ";
				fshape<<Polygons.at(i+2)<<endl;
					#ifdef PPC64
						if ((n++ % 20) == 0) fshape.flush();
					#endif
			}
		}else if(type==3){
			vector<float> dtemp=shapes.at(0)->getDMean();
			for (unsigned int j=0; j<dtemp.size();j++)
			{    
				fshape<<dtemp.at(j)<<" ";
				#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
			}
			fshape<<endl;
		
		
		}
		//enter modes
		if ((type==0)|(type==2)|(type==3)|(type==4)){	
			if ((type==2)|(type==3)){
				fshape<<"FIELD FieldData "<<2*getNumberOfModes()+5+3*getNumberOfShapes()<<endl; 
			}else if (type==0){
				fshape<<"FIELD FieldData "<<getNumberOfModes()+4<<endl;
			}else if (type==4){
				fshape<<"FIELD FieldData "<<2*getNumberOfModes()+5+3*getNumberOfShapes()+3<<endl;
			}
			if (type!=3){
			//	for (int i=0; i<getNumberOfModes(); i++)
			for (int i=0; i<numModes; i++)  
				{
					fshape<<"mode"<<i<<" "<<3<<" "<<ptTot<<" float"<<endl;
					for (int j=0;j<getNumberOfShapes();j++){
						//cout<<"saving shape "<<j<<endl;
						vector<Vec> tmpmode;
						tmpmode=shapes.at(j)->getModeVector(i);
						for (int k=0; k<getNumberOfPoints(j);k++){
							fshape<<tmpmode.at(k).X<<" ";
							fshape<<tmpmode.at(k).Y<<" ";
							fshape<<tmpmode.at(k).Z<<" "<<endl;
							#ifdef PPC64
						if ((n++ % 20) == 0) fshape.flush();
					#endif
							//	cout<<"mdoe point "<<k<<" mode "<<i<<" "<<getNumberOfPoints(j)<<endl;
						}
					}
				}
			}else{
				for (int i=0; i<getNumberOfModes(); i++)  
				{
					vector<float> tmpmode=shapes.at(0)->getDModeVector(i);
					fshape<<"mode"<<i<<" "<<1<<" "<<tmpmode.size()<<" float"<<endl;
					for (unsigned int k=0; k<tmpmode.size();k++){
						fshape<<tmpmode.at(k)<<" ";
						#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
					}
				
					fshape<<endl;
				}
			}
		}
		if ((type==0)|(type==2)|(type==3)|(type==4)){
			//set intensity mean
			if ((type==2)|(type==4)){
			fshape<<"Imean"<<" "<<shapes.at(0)->getIPP()<<" "<<ptTot<<" float"<<endl;
			for (int j=0;j<getNumberOfShapes();j++){
				//		cout<<"saving imean "<<j<<endl;
				vector<float> tmpmode;
						tmpmode=shapes.at(j)->getIMean();
						for (int k=0; k<getNumberOfPoints(j)*shapes.at(j)->getIPP();k=k+shapes.at(j)->getIPP()){
							for (int l=0;l<shapes.at(j)->getIPP();l++){
								fshape<<tmpmode.at(k+l)<<" ";
								#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
							}
							fshape<<endl;
						}
					}
				
			//	for (int i=0; i<getNumberOfModes(); i++)  
			for (int i=0; i<numModes; i++)
					{				
					//there is	an assumption here that all shapes need the same number of ipp
					fshape<<"Imode"<<i<<" "<<shapes.at(0)->getIPP()<<" "<<ptTot<<" float"<<endl;
					for (int j=0;j<getNumberOfShapes();j++){
					//cout<<"saving imode "<<i<<" "<<j<<endl;
						vector<float> tmpmode;
						tmpmode=shapes.at(j)->getIModeVector(i);
					//	cout<<"found imode "<<getNumberOfPoints(j)<<" "<<shapes.at(j)->getIPP()<<" "<<tmpmode.size()<<endl;
						for (int k=0; k<getNumberOfPoints(j)*shapes.at(j)->getIPP();k=k+shapes.at(j)->getIPP()){
							for (int l=0;l<shapes.at(j)->getIPP();l++){
							//	cout<<"hmm"<<endl;
								fshape<<tmpmode.at(k+l)<<" ";
								#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
							//	cout<<"hmmm "<<k<<endl;
							}
							fshape<<endl;
						}
					}
				}
			}else if (type==3){
				vector<float> tmpmode;
				tmpmode=shapes.at(0)->getIMean();
				fshape<<"Imean"<<" "<<shapes.at(0)->getIPP()<<" "<<tmpmode.size()<<" float"<<endl;
				//		cout<<"saving imean "<<j<<endl;
				
				for (unsigned int k=0; k<tmpmode.size();k++){
					fshape<<tmpmode.at(k)<<" ";
					#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
				}
				fshape<<endl;
				
				for (int i=0; i<getNumberOfModes(); i++)  
				{				
					//there is	an assumption here that all shapes need the same number of ipp
					fshape<<"Imode"<<i<<" "<<shapes.at(0)->getIPP()<<" "<<shapes.at(0)->getIMean().size()<<" float"<<endl;
					for (int j=0;j<getNumberOfShapes();j++){
						//cout<<"saving imode "<<j<<endl;
						vector<float> tmpmode;
						tmpmode=shapes.at(j)->getIModeVector(i);
						for (unsigned int k=0; k<tmpmode.size();k++){
							fshape<<tmpmode.at(k)<<" ";
							#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
						}
						fshape<<endl;
					}
				}
				
			}
			////////////////write Best mode
			if (type==4){
				//cout<<"Bmean"<<endl;
				//this is actually used for intensity covariance
		//		fshape<<"Bmean"<<" "<<4<<" "<<getNumberOfShapes()<<" float"<<endl;
		//		for (int j=0;j<getNumberOfShapes();j++){
		//			//		cout<<"saving imean "<<j<<endl;
		//			vector<float> tmpmode;
		//			tmpmode=shapes.at(j)->getBMean();
		//			for (int k=0; k<4;k++){
		//				fshape<<tmpmode.at(k)<<" ";
		//				#ifdef PPC64
		//				if ((n++ % 50) == 0) fshape.flush();
		//			#endif
		//			}
		//			fshape<<endl;
		//		}
				
		//		for (int i=0; i<getNumberOfModes(); i++)  
		//		{				
		//			//there is	an assumption here that all shapes need the same number of ipp
		//			fshape<<"Bmode"<<i<<" "<<shapes.at(0)->getIPP()<<" "<<shapes.at(0)->getIMean().size()<<" float"<<endl;
		//			for (int j=0;j<getNumberOfShapes();j++){
		//				cout<<"saving bmode "<<j<<endl;
		//				vector<float> tmpmode;
		//				tmpmode=shapes.at(j)->getBModeVector(i);
		//				cout<<"found bmode "<<i<<endl;
		//				for (int k=0; k<4;k++){
		//					fshape<<tmpmode.at(k)<<" ";
		//					cout<<"k "<<k<<endl;
		//				}
		//				fshape<<endl;
		//			}
		//		}
				
			
			
			}
			
			//write covariance matrix of predictive structure...
			for (int i=0; i < 1;i++){//getNumberOfShapes();i++){
				
				
				if (type!=0){
					//cout<<"Shape "<<i<<" end"<<endl;
					
					
					//////////////////////////////enter Intensity conditionals
					vector< vector<float> > precimat = shapes.at(i)->getICondPrec();
					fshape<<"iCondPrec"<<i<<" "<<precimat.at(0).size()<<" "<<precimat.size()<<" float"<<endl;
					for (unsigned int j=0;j<precimat.size();j++){
						for (unsigned int k=0; k<precimat.at(j).size();k++){
							fshape<<precimat.at(j).at(k)<<" ";
							#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
						}
						fshape<<endl;
					}
					
					
					fshape<<"iCondEigs"<<i<<" "<<1<<" "<<shapes.at(i)->getICondEigs().size()<<" float"<<endl;
					for (unsigned int j=0;j<shapes.at(i)->getICondEigs().size();j++){
						
						fshape<<shapes.at(i)->getICondEigs().at(j)<<" ";
						#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
					}
					fshape<<endl;
					//cout<<"Shape "<<i<<" end"<<endl;
					
				//	if (type==6){
					//////////////////////////////enter Intensity covariance
					precimat = shapes.at(i)->getBCondPrec();
					fshape<<"bCondPrec"<<i<<" "<<precimat.at(0).size()<<" "<<precimat.size()<<" float"<<endl;
					for (unsigned int j=0;j<precimat.size();j++){
					//cout<<"bcond "<<j<<" "<<precimat.at(j).size()<<endl;
						for (unsigned int k=0; k<precimat.at(j).size();k++){
				//		if ((k==0)|(k==(precimat.at(j).size()-1))){
				//		cout<<precimat.at(j).at(k)<<" "<<endl;
				//		}
							fshape<<precimat.at(j).at(k)<<" ";
							#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
						}
						fshape<<endl;
					}
					
					
					fshape<<"bCondEigs"<<i<<" "<<1<<" "<<shapes.at(i)->getBCondEigs().size()<<" float"<<endl;
					for (unsigned int j=0;j<shapes.at(i)->getICondEigs().size();j++){
						
						fshape<<shapes.at(i)->getBCondEigs().at(j)<<" ";
						#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
					}
					fshape<<endl;
					
					
					
				}
					
					
				
				fshape<<"ErrPriors"<<i<<" "<<1<<" "<<shapes.at(i)->getErrs().size()<<" float"<<endl;
				for (unsigned int j=0;j<shapes.at(i)->getErrs().size();j++){
					fshape<<shapes.at(i)->getErrs().at(j)<<" ";
					#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
				}
				fshape<<endl;
				if (type==4){
						//////////////////////////////enter Best conditionals
			/*			vector< vector<float> > condmat = shapes.at(i)->getBCondMat();
						fshape<<"bCondMat"<<i<<" "<<condmat.at(0).size()<<" "<<condmat.size()<<" float"<<endl;
						for (unsigned int j=0;j<condmat.size();j++){
							vector<float> vtemp;
							for (unsigned int i=0;i<condmat.size();i++){
								fshape<<condmat.at(i).at(j)<<" ";
							}
							fshape<<endl;
						}
			*/
			
							//			precimat = shapes.at(i)->getBCondPrec();
			//			fshape<<"bCondPrec"<<i<<" "<<precimat.at(0).size()<<" "<<precimat.size()<<" float"<<endl;
			//			for (unsigned int j=0;j<precimat.size();j++){
			//				for (unsigned int k=0; k<precimat.at(j).size();k++){
			//					fshape<<precimat.at(j).at(k)<<" ";
			//				}
			//				fshape<<endl;
			//			}
						
						
			//			fshape<<"bCondEigs"<<i<<" "<<1<<" "<<shapes.at(i)->getBCondEigs().size()<<" float"<<endl;
			//			for (unsigned int j=0;j<shapes.at(i)->getBCondEigs().size();j++){
			//				
			//				fshape<<shapes.at(i)->getBCondEigs().at(j)<<" ";
			//			}
			//			fshape<<endl;
						
				}
				//cout<<"Shape "<<i<<" end"<<endl;
				
			}
			if (type!=0){
				//encoding eigenvalues
				fshape<<"eigenValuesI 1 "<<getNumberOfModes()+1<<" float"<<endl;
				for (int i=0;i<getNumberOfModes();i++)
				{
					fshape<<eigenValsI.at(i)<<" ";
					#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
				}
				fshape<<sumEigValsI<<endl;
			}
			
			//encoding number of points per structure
			fshape<<"numPoints 1 "<<getNumberOfShapes()<<" int"<<endl;
			for (int i=0;i<getNumberOfShapes();i++)
			{
				//cout<<getNumberOfPoints(0)<<endl;
				fshape<<getNumberOfPoints(i)<<" ";//getNumberOfPoints.at(i)<<" ";
				#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
			}
			fshape<<endl;
			//encoding eigenvalues
			fshape<<"eigenValues 1 "<<getNumberOfModes()+1<<" float"<<endl;
			for (int i=0;i<getNumberOfModes();i++)
			{
				fshape<<eigenVals.at(i)<<" ";
				#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
			}
			fshape<<sumEigVals<<endl;
			
			//encoding scalar labels
			fshape<<"labels 1 "<<getNumberOfShapes()<<" int"<<endl;
			for (int i=0;i<getNumberOfShapes();i++)
			{
			//	cout<<"labels "<<i<<" "<<labels.size()<<" "<<labels.at(i)<<endl;
				fshape<<labels.at(i)<<" ";
				#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
			}
			fshape<<endl;	
			
			fshape<<"numSubjects 1 1 int"<<endl;
			fshape<<numSubjects<<endl;
			if (type!=0){
				//IPP for each shape
				fshape<<"IPP 1 "<<getNumberOfShapes()<<" int"<<endl;
				for (int i=0;i<getNumberOfShapes();i++)
				{
					fshape<<getShape(i)->getIPP()<<" ";
					#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
				}
				fshape<<endl;	
			}
		}else if(type==1){
			//write intensity information
			fshape<<"FIELD FieldData 2"<<endl;
			//cout<<"IPROFS "<<iprofs.at(0).size()<<endl;
			fshape<<"intensity_profiles "<<iprofs.at(0).size()<<" "<<getNumberOfPoints(0)<<" float"<<endl;
			for (int i =0;i<getNumberOfPoints(0);i++){
				for (unsigned int j=0; j<iprofs.at(i).size();j++){
					fshape<<iprofs.at(i).at(j)<<" ";
					#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
				}
				fshape<<endl;
			}
			fshape<<"planarEstimates "<<vbest.size()<<" "<<1<<" float"<<endl;
			for (unsigned int i =0;i<vbest.size();i++){
					fshape<<vbest.at(i)<<" ";
					#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
			}
				fshape<<endl;
		}else if(type==5){
			//write t-stats
			fshape<<"POINT_DATA "<<getTotalNumberOfPoints()<<endl;
			fshape<<"VECTORS tstast float"<<endl;
			for (int sh=0;sh<getNumberOfShapes();sh++){
				vector<float> tstatx=shapes.at(sh)->getTstatX();
				vector<float> tstaty=shapes.at(sh)->getTstatY();
				vector<float> tstatz=shapes.at(sh)->getTstatZ();

				for (int i=0;i<getNumberOfPoints(sh);i++){
					fshape<<tstatx.at(i)<<" "<<tstaty.at(i)<<" "<<tstatz.at(i)<<endl;
					#ifdef PPC64
						if ((n++ % 20) == 0) fshape.flush();
					#endif
				}
			}
		
			
		}
	//	cout<<"closing shape/appearance model"<<endl;
		fshape.close();
	//	cout<<"closed"<<endl;
	}
	
	void shapeModel::save_binary(string s, int type, int numModes){
	//this is implemented for  types 0 and 1
		//numModes currently doesn't do anything
		//type =0 for normal shape model
		//type =1 for appearance mesh
		//type 2 = appearance model
		//type 3 = distance/appearance model
		//type 4= appearance model with planar estimates
		//type 5= mesh and t-stats
		//type 6= Affmodes
		ofstream fshape;  
		fshape.open(s.c_str(),ios::out|ios::binary);
		//fshape.open(s.c_str());

		#ifdef PPC64
    int n=0;
#endif
		//calculate total number of points
		int ptTot=0;
		for (int i=0; i<getNumberOfShapes();i++){
			ptTot+=getNumberOfPoints(i);
		}
	
		fshape<<"# vtk DataFile Version 3.0"<<endl;
		//add in a binary number to test for endianess
		unsigned int bintest=42;
	//	cout<<"size of unsigned int "<<sizeof(unsigned int)<<endl;
		fshape.write(reinterpret_cast<char *>(&bintest),sizeof(bintest));
		//cout<<"cast 42 "<<reinterpret_cast<char *>(&bintest)<<endl;
		fshape<<"IBIM modes of variation file "<<endl<<"BINARY"<<endl<<"DATASET POLYDATA"<<endl<<"POINTS "<<ptTot<<"  float"<<endl;
	
		
		
		if (type!=3){
			//writing points to vtk file
			for (int i=0; i<getNumberOfShapes(); i++) { 
				//	cout<<"wirting point for shape i= "<<i<<endl;
				//	cout<<getNumberOfPoints(i)<<endl;
				for (int j=0; j<getNumberOfPoints(i);j++)
				{    
				//	cout<<"j "<<j<<endl;
					float coord;
					coord=shapes.at(i)->getPoint(j)->get_coord().X;
					fshape.write(reinterpret_cast<char *>(&coord),sizeof(coord));
					coord=shapes.at(i)->getPoint(j)->get_coord().Y;
					fshape.write(reinterpret_cast<char *>(&coord),sizeof(coord));
					coord=shapes.at(i)->getPoint(j)->get_coord().Z;
					fshape.write(reinterpret_cast<char *>(&coord),sizeof(coord));
					#ifdef PPC64
						if ((n++ % 20) == 0) fshape.flush();
					#endif
				}
			}
			
			//writing connectuions to vtk file
			int cumNumPoints=0;
			vector<int> Polygons;
			for (int i=0;i<getNumberOfShapes();i++){
				
				Mesh mtemp=shapes.at(i)->getMesh();
				for ( list<Triangle*>::const_iterator t=mtemp._triangles.begin(); t!=mtemp._triangles.end(); t++) {
					
					Polygons.push_back((*t)->get_vertice(0)->get_no()+cumNumPoints);
					Polygons.push_back((*t)->get_vertice(1)->get_no()+cumNumPoints);
					Polygons.push_back((*t)->get_vertice(2)->get_no()+cumNumPoints);	  
				}
				cumNumPoints+=getNumberOfPoints(i);
			}
			
			fshape<<"POLYGONS "<<Polygons.size()/3<<" "<<Polygons.size()/3*4<<endl;
			for (unsigned int i=0;i<Polygons.size();i=i+3){
				int thr=3;
				int poly;
				fshape.write(reinterpret_cast<char *>(&thr),sizeof(thr));
				poly=Polygons.at(i);
				fshape.write(reinterpret_cast<char *>(&poly),sizeof(poly));
				poly=Polygons.at(i+1);
				fshape.write(reinterpret_cast<char *>(&poly),sizeof(poly));
				poly=Polygons.at(i+2);
				fshape.write(reinterpret_cast<char *>(&poly),sizeof(poly));
				
					#ifdef PPC64
						if ((n++ % 20) == 0) fshape.flush();
					#endif
			}
		}	//	cout<<"closing shape/appearance model"<<endl;
		//enter modes
	

				//is chnaged to 6 because intesiyt covaricne
				
				
		//	fshape<<"FIELD FieldData "<<2*getNumberOfModes()+6+3*getNumberOfShapes()<<endl; 
			//start with shape model only
	//		fshape<<"FIELD FieldData "<<getNumberOfModes()<<endl;
	if ((type!=6)&&(type!=7)){
			fshape<<"FIELD FieldData "<<2*getNumberOfModes()+6+7<<endl;
		}else{
					fshape<<"FIELD FieldData "<<2*getNumberOfModes()+6+7+7+8<<endl;

		}
	//	cout<<"nummodes "<<getNumberOfModes()<<endl;
			//	for (int i=0; i<getNumberOfModes(); i++)
			//for (int i=0; i<numModes; i++)  
			for (int i=0; i<getNumberOfModes(); i++)
				{
					fshape<<"mode"<<i<<" "<<3<<" "<<ptTot<<" float"<<endl;
					for (int j=0;j<getNumberOfShapes();j++){
						//cout<<"saving shape "<<j<<endl;
						vector<Vec> tmpmode;
						tmpmode=shapes.at(j)->getModeVector(i);
						for (int k=0; k<getNumberOfPoints(j);k++){
							float coord;
							coord=tmpmode.at(k).X;
							fshape.write(reinterpret_cast<char *>(&coord),sizeof(coord));
							coord=tmpmode.at(k).Y;
							fshape.write(reinterpret_cast<char *>(&coord),sizeof(coord));
							coord=tmpmode.at(k).Z;
							fshape.write(reinterpret_cast<char *>(&coord),sizeof(coord));
							
							#ifdef PPC64
							if ((n++ % 20) == 0) fshape.flush();
							#endif
							//	cout<<"mdoe point "<<k<<" mode "<<i<<" "<<getNumberOfPoints(j)<<endl;
						}
					}
				}
			
		
			//set intensity mean
		if (type!=0){
			fshape<<"Imean"<<" "<<shapes.at(0)->getIPP()<<" "<<ptTot<<" float"<<endl;
			for (int j=0;j<getNumberOfShapes();j++){
				//		cout<<"saving imean "<<j<<endl;
				vector<float> tmpmode;
						tmpmode=shapes.at(j)->getIMean();
						for (int k=0; k<getNumberOfPoints(j)*shapes.at(j)->getIPP();k=k+shapes.at(j)->getIPP()){
							for (int l=0;l<shapes.at(j)->getIPP();l++){
								float val;
								val=tmpmode.at(k+l);
								//cout<<"write imean "<<val<<endl;
								fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
								
								#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
							}
						//	fshape<<endl;
						}
					}
				
				for (int i=0; i<getNumberOfModes(); i++)  
		//	for (int i=0; i<numModes; i++)
					{				
					//there is	an assumption here that all shapes need the same number of ipp
					fshape<<"Imode"<<i<<" "<<shapes.at(0)->getIPP()<<" "<<ptTot<<" float"<<endl;
					for (int j=0;j<getNumberOfShapes();j++){
						vector<float> tmpmode;
						tmpmode=shapes.at(j)->getIModeVector(i);
						for (int k=0; k<getNumberOfPoints(j)*shapes.at(j)->getIPP();k=k+shapes.at(j)->getIPP()){
							for (int l=0;l<shapes.at(j)->getIPP();l++){
							float val;
								val=tmpmode.at(k+l);
								fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
			
								#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
							//	cout<<"hmmm "<<k<<endl;
							}
						//	fshape<<endl;
						}
					}
				}
				
			}
				
				
				
			//write covariance matrix of predictive structure...
			for (int i=0; i < 1;i++){//getNumberOfShapes();i++){
				
				
				if (type!=0){
					//cout<<"Shape "<<i<<" end"<<endl;
					
					
					//////////////////////////////enter Intensity conditionals
					vector< vector<float> > precimat = shapes.at(i)->getICondPrec();
					fshape<<"iCondPrec"<<i<<" "<<precimat.at(0).size()<<" "<<precimat.size()<<" float"<<endl;
					for (unsigned int j=0;j<precimat.size();j++){
						for (unsigned int k=0; k<precimat.at(j).size();k++){
							float val;
							val=precimat.at(j).at(k);
							fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
#ifdef PPC64
							if ((n++ % 50) == 0) fshape.flush();
#endif
						}
						//fshape<<endl;
					}
					
					
					fshape<<"iCondEigs"<<i<<" "<<1<<" "<<shapes.at(i)->getICondEigs().size()<<" float"<<endl;
					for (unsigned int j=0;j<shapes.at(i)->getICondEigs().size();j++){
						float val;
						val=shapes.at(i)->getICondEigs().at(j);
						fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
						
#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
#endif
					}
					//fshape<<endl;
					
					//////////////////////////////enter Intensity covariance
					precimat = shapes.at(i)->getBCondPrec();
					fshape<<"bCondPrec"<<i<<" "<<precimat.at(0).size()<<" "<<precimat.size()<<" float"<<endl;
					for (unsigned int j=0;j<precimat.size();j++){
						//cout<<"bcond "<<j<<" "<<precimat.at(j).size()<<endl;
						for (unsigned int k=0; k<precimat.at(j).size();k++){
							float val;
							val=precimat.at(j).at(k);
							fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
							
#ifdef PPC64
							if ((n++ % 50) == 0) fshape.flush();
#endif
						}
						//fshape<<endl;
					}
					
					
					fshape<<"bCondEigs"<<i<<" "<<1<<" "<<shapes.at(i)->getBCondEigs().size()<<" float"<<endl;
					for (unsigned int j=0;j<shapes.at(i)->getICondEigs().size();j++){
						float val;
						val=shapes.at(i)->getBCondEigs().at(j);
						fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
						
#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
#endif
					}
					//fshape<<endl;
					
					
					
				}
					
					
				
				fshape<<"ErrPriors"<<i<<" "<<1<<" "<<shapes.at(i)->getErrs().size()<<" float"<<endl;
				for (unsigned int j=0;j<shapes.at(i)->getErrs().size();j++){
					float err;
					err=shapes.at(i)->getErrs().at(j);
					fshape.write(reinterpret_cast<char *>(&err),sizeof(err));
					#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
				}
				//fshape<<endl;
		
				//cout<<"Shape "<<i<<" end"<<endl;
				
			}


			//encoding number of points per structure
			fshape<<"numPoints 1 "<<getNumberOfShapes()<<" int"<<endl;
			for (int i=0;i<getNumberOfShapes();i++)
			{
				int num;
				num=getNumberOfPoints(i);
				fshape.write(reinterpret_cast<char *>(&num),sizeof(num));
			
				#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
			}
			//fshape<<endl;
			
			//encoding eigenvalues Intesnity
			{
				fshape<<"eigenValuesI 1 "<<getNumberOfModes()+1<<" float"<<endl;
				for (int i=0;i<getNumberOfModes();i++)
				{
					float eig;
					eig=eigenValsI.at(i);
					fshape.write(reinterpret_cast<char *>(&eig),sizeof(eig));
					#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
				}
				//fshape<<sumEigValsI;//<<endl;
				float eig;
			eig=sumEigValsI;
			fshape.write(reinterpret_cast<char *>(&eig),sizeof(eig));
			}
			{
			////encoding eigenvalues shape
			fshape<<"eigenValues 1 "<<getNumberOfModes()+1<<" float"<<endl;
			for (int i=0;i<getNumberOfModes();i++)
			{
				float eig;
				eig=eigenVals.at(i);
			//cout<<"eig "<<eig<<endl;
				fshape.write(reinterpret_cast<char *>(&eig),sizeof(eig));
				#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
			}
			float eig;
			eig=sumEigVals;
			fshape.write(reinterpret_cast<char *>(&eig),sizeof(eig));
			}
			
			
			//encoding scalar labels
			fshape<<"labels 1 "<<getNumberOfShapes()<<" int"<<endl;
			for (int i=0;i<getNumberOfShapes();i++)
			{
			//	cout<<"labels "<<i<<" "<<labels.size()<<" "<<labels.at(i)<<endl;
				int lab;
				lab=labels.at(i);
				fshape.write(reinterpret_cast<char *>(&lab),sizeof(lab));
				
				#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
					#endif
			}
			//fshape<<endl;	
			
			fshape<<"numSubjects 1 1 int"<<endl;
			
			fshape.write(reinterpret_cast<char *>(&numSubjects),sizeof(numSubjects));
			
			//IPP for each shape
			fshape<<"IPP 1 "<<getNumberOfShapes()<<" int"<<endl;
			for (int i=0;i<getNumberOfShapes();i++)
			{
				int ipp;
				ipp=getShape(i)->getIPP();
				fshape.write(reinterpret_cast<char *>(&ipp),sizeof(ipp));
				
#ifdef PPC64
				if ((n++ % 50) == 0) fshape.flush();
#endif
			}
			//fshape<<endl;	
			if ((type==6)|(type==7)){
		//		cout<<"load aff model stuff"<<endl;
			//write out affine modes
			//write translations
			//get all relevant data
			
			
			fshape<<"transVec 3 3 float"<<endl;
		//	cout<<"first trvec"<<shapes.at(0)->getAfftxVecs().at(0)<<endl;
				fshape.write(reinterpret_cast<char *>(&(shapes.at(0)->getAfftxVecs().at(0))),sizeof(float));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAfftxVecs().at(1)),sizeof(shapes.at(0)->getAfftxVecs().at(1)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAfftxVecs().at(2)),sizeof(shapes.at(0)->getAfftxVecs().at(2)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAfftyVecs().at(0)),sizeof(shapes.at(0)->getAfftyVecs().at(0)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAfftyVecs().at(1)),sizeof(shapes.at(0)->getAfftyVecs().at(1)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAfftyVecs().at(2)),sizeof(shapes.at(0)->getAfftyVecs().at(2)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAfftzVecs().at(0)),sizeof(shapes.at(0)->getAfftzVecs().at(0)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAfftzVecs().at(1)),sizeof(shapes.at(0)->getAfftzVecs().at(1)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAfftzVecs().at(2)),sizeof(shapes.at(0)->getAfftzVecs().at(2)));
				fshape<<"transEigs 1 3 float"<<endl;
				float temp=shapes.at(0)->getAfftxEig();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAfftxEig()));
				temp=shapes.at(0)->getAfftyEig();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAfftyEig()));
				temp=shapes.at(0)->getAfftzEig();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAfftzEig()));
				fshape<<"transMean 1 3 float"<<endl;
				temp=shapes.at(0)->getAfftxMean();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAfftxMean()));
				temp=shapes.at(0)->getAfftyMean();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAfftyMean()));
				temp=shapes.at(0)->getAfftzMean();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAfftzMean()));
				//write rotations
			fshape<<"RotationVecs 3 3 float"<<endl;
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAffrxVecs().at(0)),sizeof(shapes.at(0)->getAffrxVecs().at(0)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAffrxVecs().at(1)),sizeof(shapes.at(0)->getAffrxVecs().at(1)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAffrxVecs().at(2)),sizeof(shapes.at(0)->getAffrxVecs().at(2)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAffryVecs().at(0)),sizeof(shapes.at(0)->getAffryVecs().at(0)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAffryVecs().at(1)),sizeof(shapes.at(0)->getAffryVecs().at(1)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAffryVecs().at(2)),sizeof(shapes.at(0)->getAffryVecs().at(2)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAffrzVecs().at(0)),sizeof(shapes.at(0)->getAffrzVecs().at(0)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAffrzVecs().at(1)),sizeof(shapes.at(0)->getAffrzVecs().at(1)));
				fshape.write(reinterpret_cast<char *>(&shapes.at(0)->getAffrzVecs().at(2)),sizeof(shapes.at(0)->getAffrzVecs().at(2)));
				fshape<<"RotationEigs 1 3 float"<<endl;
				temp=shapes.at(0)->getAffrxEig();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAffrxEig()));
				temp=shapes.at(0)->getAffryEig();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAffryEig()));
				temp=shapes.at(0)->getAffrzEig();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAffrzEig()));
				fshape<<"RotationMean 1 3 float"<<endl;
				temp=shapes.at(0)->getAffrxMean();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAffrxMean()));
				temp=shapes.at(0)->getAffryMean();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAffryMean()));
				temp=shapes.at(0)->getAffrzMean();
				fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAffrzMean()));
					//write scale
					if (type==7){
					fshape<<"scMeanVar 1 2 float"<<endl;
					temp=shapes.at(0)->getAffscMean();
					fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAffscMean()));
					temp=shapes.at(0)->getAffscEig();
					fshape.write(reinterpret_cast<char *>(&temp),sizeof(shapes.at(0)->getAffscEig()));
					}
					
					
					///deal with intensity predicttions at the moment
					int numAff=7;
					if (type==6){
						numAff=6;
					}
					for (int i=0; i<numAff; i++)  
						//	for (int i=0; i<numModes; i++)
					{		//	cout<<"write aff icond "<<i<<endl;	
						//there is	an assumption here that all shapes need the same number of ipp
						fshape<<"AffImode"<<i<<" "<<shapes.at(0)->getIPP()<<" "<<ptTot<<" float"<<endl;
						for (int j=0;j<getNumberOfShapes();j++){
							vector<float> tmpmode;
							tmpmode=shapes.at(j)->getAffIModeVector(i);
							for (int k=0; k<getNumberOfPoints(j)*shapes.at(j)->getIPP();k=k+shapes.at(j)->getIPP()){
								for (int l=0;l<shapes.at(j)->getIPP();l++){
									float val;
									val=tmpmode.at(k+l);
									fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
									
#ifdef PPC64
									if ((n++ % 50) == 0) fshape.flush();
#endif
									//	cout<<"hmmm "<<k<<endl;
								}
								//	fshape<<endl;
							}
						}
					}
					
					//////////////////////////////enter Intensity conditionals
					for (int i=0; i<1;i++){//shapes
					vector< vector<float> > precimat = shapes.at(i)->getICondPrecAff();
					//cout<<"write aff icondMat "<<i<<endl;	
					fshape<<"iAffCondPrec"<<i<<" "<<precimat.at(0).size()<<" "<<precimat.size()<<" float"<<endl;
					for (unsigned int j=0;j<precimat.size();j++){
						for (unsigned int k=0; k<precimat.at(j).size();k++){
							float val;
							val=precimat.at(j).at(k);
							fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
#ifdef PPC64
							if ((n++ % 50) == 0) fshape.flush();
#endif
						}
						//fshape<<endl;
					}
					
					//cout<<"write aff icondEig"<<i<<endl;	
					fshape<<"iAffCondEigs"<<i<<" "<<1<<" "<<shapes.at(i)->getICondEigsAff().size()<<" float"<<endl;
					for (unsigned int j=0;j<shapes.at(i)->getICondEigsAff().size();j++){
						float val;
						val=shapes.at(i)->getICondEigsAff().at(j);
						fshape.write(reinterpret_cast<char *>(&val),sizeof(val));
						
#ifdef PPC64
						if ((n++ % 50) == 0) fshape.flush();
#endif
					}
					}
					//fshape<<endl;
					
					
					
					
			}
			
			
			
		
		fshape.close();
	//	cout<<"closed"<<endl;
	}
	
#define BINFLAG 42

void shapeModel::load_bmv_binaryInfo(string s, int type){
	//int type=0;
	//type 0 = normal shape model
	//type 1 = appearance model
	//type 2= affine component(doesn't chnager anything)
	//cout<<"Loading BinaryInfo Shape Model..."<<endl;
	
			string temp;
		int tempInt,tempInt2;
		float tempFloat, vX,vY,vZ;
		ifstream f;
		int numModes=0;
		
		//read in information (not meshes yet!)
		//	cout<<"reading information"<<endl;
		f.open(s.c_str(), ios::in | ios::binary);
		getline(f,temp);
		//getline(f,temp);//used to get rid of second ascii line to avoid "modes"
		//read in line upto binary test
		//f>>temp>>temp>>temp>>temp>>temp;
		//cout<<temp<<endl;				
			
		unsigned int testval;
		// test for byte swapping
		f.read(reinterpret_cast<char *>(&testval),sizeof(testval));
		
		bool swapbytes=false;

		if (testval!=BINFLAG) {
			swapbytes = true;
			 Swap_Nbytes(1,sizeof(testval),&testval);
		//	cout<<testval<<" "<<BINFLAG<<endl;
			//cout<<"failed endian test"<<endl;
		}//else{
			//cout<<"endianess is ok"<<endl;
		//	}
			
		if (testval!=BINFLAG) { 
			cerr << "Unrecognised binary matrix file format" << endl;
		}
		

		//eliminate space
	//	char * ctemp;
	//	ctemp=new char[50];
	//	cout<<"ctemp size "<<sizeof(ctemp)<<endl;
	//	f.read(ctemp,sizeof(ctemp));
	//	for (int i=0;i<50;i++){
	//		cout<<*ctemp;
	//		ctemp++;
	//	}
	//	cout<<endl;
		//cout<<"Thhis should be a space"<<*ctemp<<"is there a space"<<endl;
		
		//endiantest
		//bool swapbytes = false;
		
		
			
		
		getline(f,temp);//get rid of endl;
		//	cout<<"this should be blank "<<temp<<endl;
		getline(f,temp);
		getline(f,temp);//this is polydata line
		//	cout<<temp<<endl;
			int totPoints;
			f>>temp>>totPoints;//read in number of points
				getline(f,temp); //assumes float, so throws away rest of line
				//	cout<<temp<<endl;
				
				//read points
				//	cout<<"number of points: "<<totPoints<<endl;
				for (int i=0;i<3*totPoints;i++){
				//don't need to byyte swap here because not using values
					f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				}
				f>>temp>>tempInt;//stores number of polygons
				//	cout<<temp<<endl;
					getline(f,temp); 
					//	cout<<temp<<endl;
					int numPolys=tempInt;
					for (int i=0;i<numPolys;i++){
						//don't need to byyte swap here because not using values
						f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
						//throw "3" away, we knwo we're dealing with triangles
						f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
						f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
						f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
					
					}
					
					//read FIELD Data
					int numField; 
					f>>temp>>temp>>numField;
					getline(f,temp);
					
					//	cout<<"should be blank: "<<temp<<endl;//this should be 
						vector<float> verrs;
						//			cout<<"numFiedl "<<numField<<endl;
						for (int i=0;i<numField;i++){
							f>>temp;
							//cout<<"temp "<<temp<<endl;
							if (!strcmp(temp.substr(0,4).c_str(),"mode")){
							  //cout<<"found "<<temp<<endl;
								string temp2;
								getline(f,temp);//need getlines before reads
									
									for (int j=0;j<totPoints;j++){
										f.read(reinterpret_cast<char*>(&vX),sizeof(float));
										f.read(reinterpret_cast<char*>(&vY),sizeof(float));
										f.read(reinterpret_cast<char*>(&vZ),sizeof(float));
										
									}
									
							}else if(!strcmp(temp.substr(0,5).c_str(),"Imean")){
							  //cout<<"found "<<temp<<endl;

								int ipp; 
								string temp2;
								f>>ipp>>tempInt>>temp2;
								getline(f,temp); 
								for (int j=0;j<totPoints*ipp;j++){
									f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
								//	cout<<"imean read "<<tempFloat<<endl;
								}
							
							}else if (!strcmp(temp.substr(0,5).c_str(),"Imode")){
							  //		cout<<"found "<<temp<<endl;
								string temp2;
								int ipp;
								f>>ipp>>tempInt>>temp2;
								getline(f,temp);
						//		for (int i=0;i<numShapes;i++){
						//			for (int j=0;j<npts.at(i);j++){
						//				for (int k=0;k<ipp;k++){
						for (int j=0;j<totPoints*ipp;j++){

										f.read(reinterpret_cast<char*>(&vX),sizeof(float));
										}
						//				}
						//			}
						//		}
								
							}else if (!strcmp(temp.substr(0,5).c_str(),"iCond")){
							  //	cout<<"found "<<temp<<endl;
								//cout<<"iconds found"<<endl;
								int sh;
								sh = atoi(temp.substr(9,1).c_str());
								if (!strcmp(temp.substr(0,9).c_str(),"iCondPrec")){
								  //		cout<<"found "<<temp<<endl;
									string temp2;
									int tempInt2;
									f>>tempInt2>>tempInt>>temp2;
									getline(f,temp); 
									for (int i=0;i<tempInt;i++){
										for (int j=0;j<tempInt2;j++){
											f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
										}
									}
									
								}else if (!strcmp(temp.substr(0,9).c_str(),"iCondEigs")){
								  //		cout<<"found "<<temp<<endl;
									//cout<<"enter icondeigs"<<endl;
									string temp2;
									vector<float> eigstemp;
									f>>tempInt>>tempInt>>temp2;
									getline(f,temp); 
									for (int j=0;j<tempInt;j++){
										f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
									}
								}
								//		cout<<"iconds loaded"<<endl;
							}else if ((!strcmp(temp.substr(0,5).c_str(),"bCond"))){
							  //	cout<<"found bconds "<<temp<<endl;
								int sh;
								sh = atoi(temp.substr(9,1).c_str());
								if (!strcmp(temp.substr(0,9).c_str(),"bCondPrec")){
								  //		cout<<"found "<<temp<<endl;
									string temp2;
									int tempInt2;
									f>>tempInt2>>tempInt>>temp2;	
									getline(f,temp); 				
									vector< vector<float> > prectemp;							
									for (int i=0;i<tempInt;i++){
										vector<float> col;
										for (int j=0;j<tempInt2;j++){
											f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
										}
									}
									
									
								}else if (!strcmp(temp.substr(0,9).c_str(),"bCondEigs")){
								  //		cout<<"found "<<temp<<endl;
									//cout<<"enter bcondeigs"<<endl;
									string temp2;
									vector<float> eigstemp;
									f>>tempInt>>tempInt>>temp2;
									getline(f,temp); 
									
									for (int j=0;j<tempInt;j++){
										f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
									}
									
								}
								
							}else if (!strcmp(temp.c_str(),"numPoints")){
															//cout<<"found "<<temp<<endl;

								f>>numShapes>>numShapes>>temp;
								getline(f,temp);
								for (int i=0;i<numShapes;i++){
									f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
									if (swapbytes){
									Swap_Nbytes(1,sizeof(tempInt),&tempInt);
									}
									npts.push_back(tempInt);
								}
							}else if (!strcmp(temp.c_str(),"eigenValues")){
															//cout<<"found "<<temp<<endl;

								f>>tempInt2>>tempInt2>>temp;
								getline(f,temp);
								//set the number of modes by the number of eigenvalues
								numModes=tempInt2-1;
								for (int i=0;i<(tempInt2-1);i++){
									f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
									if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
									eigenVals.push_back(tempFloat);
								}
								f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
								if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
								sumEigVals=tempFloat;
							}else if (!strcmp(temp.c_str(),"labels")){
							  //			cout<<"found "<<temp<<endl;
								f>>tempInt2>>tempInt2>>temp;
								getline(f,temp);
								for (int i=0;i<tempInt2;i++){
									f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
									if (swapbytes){
									Swap_Nbytes(1,sizeof(tempInt),&tempInt);
									}
									labels.push_back(tempInt);
								}
							}else if (!strcmp(temp.c_str(),"numSubjects")){
							  //			cout<<"found "<<temp<<endl;
								f>>tempInt2>>tempInt2>>temp;
								getline(f,temp);
								//	for (int i=0;i<tempInt2;i++){
								f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
								if (swapbytes){
									Swap_Nbytes(1,sizeof(tempInt),&tempInt);
									}
								numSubjects=tempInt;
								//	cout<<numSubjects<<"NUM"<<endl;
								//	}
							}else if (!strcmp(temp.c_str(),"IPP")){
							  //		cout<<"found "<<temp<<endl;
								f>>tempInt2>>tempInt2>>temp;
								//		cout<<"IPP "<<tempInt2<<endl;
								getline(f,temp);
								for (int i=0;i<tempInt2;i++){
									f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
									if (swapbytes){
									Swap_Nbytes(1,sizeof(tempInt),&tempInt);
									}
									ipps.push_back(tempInt);
									//cout<<"IPP "<<tempInt<<endl;
									//need to assign it to shape objects stiill
								}
							}else if (!strcmp(temp.c_str(),"eigenValuesI")){
							  //	cout<<"found "<<temp<<endl;
								f>>tempInt2>>tempInt2>>temp;
								getline(f,temp);
								for (int i=0;i<(tempInt2-1);i++){
									f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
									if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
									eigenValsI.push_back(tempFloat);
								}
								f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
								if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
								sumEigValsI=tempFloat;
							}else if (!strcmp(temp.substr(0,9).c_str(),"ErrPriors")){
							  //		cout<<"found "<<temp<<endl;
								int sh;
								sh = atoi(temp.substr(9,1).c_str());
								string temp2;
								f>>tempInt>>tempInt>>temp2;
								getline(f,temp);
								//		cout<<"errPri "<<temp<<endl;
								//cout<<"enter errpiros"<<endl;
								
								for (int j=0;j<2;j++){ //this allows for error in intensity and shape, for sh=0 its is the predictor error
									f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
								
								}
								
							}
							
							
							
							
		}
					
						
	f.close();
	}
int shapeModel::load_bmv_binary(string s, int type){
	//int type=0;
	//type 0 = normal shape model
	//type 1 = appearance model
	//typw 2 = affine component
	//cout<<"Loading Binary Shape Model..."<<endl;
	
	
	//		cout<<"constructing shape Model"<<endl;
	//Data is read now set up shape model
	//craete shape objects
	for (int i=0;i<numShapes;i++){
	  //		cout<<"adding shape "<<i;
		MShape*  Imshape = new MShape();
		shapes.push_back(Imshape);
		if (type>=1){
			shapes.at(i)->setIPP(ipps.at(i));
			//			cout<<"with IPP "<<ipps.at(i);
		}
		//cout<<endl;
	}

	
	
	
	
	
			string temp;
		int tempInt,tempInt2;
		float tempFloat, vX,vY,vZ;
		ifstream f;
		int numModes=0;
		
		//read in information (not meshes yet!)
		//	cout<<"reading information"<<endl;
		f.open(s.c_str(), ios::in | ios::binary);
		getline(f,temp);
		//getline(f,temp);//used to get rid of second ascii line to avoid "modes"
		//read in line upto binary test
		//f>>temp>>temp>>temp>>temp>>temp;
		//cout<<temp<<endl;				
			
		unsigned int testval;
		// test for byte swapping
		f.read(reinterpret_cast<char *>(&testval),sizeof(testval));
		
		bool swapbytes=false;

		if (testval!=BINFLAG) {
			swapbytes = true;
			 Swap_Nbytes(1,sizeof(testval),&testval);
			//cout<<testval<<" "<<BINFLAG<<endl;
			//cout<<"failed endian test"<<endl;
		}//else{
			//cout<<"endianess is ok"<<endl;
			//}
			
		if (testval!=BINFLAG) { 
			cerr << "Unrecognised binary matrix file format" << endl;
		}
		

		//eliminate space
	//	char * ctemp;
	//	ctemp=new char[50];
	//	cout<<"ctemp size "<<sizeof(ctemp)<<endl;
	//	f.read(ctemp,sizeof(ctemp));
	//	for (int i=0;i<50;i++){
	//		cout<<*ctemp;
	//		ctemp++;
	//	}
	//	cout<<endl;
		//cout<<"Thhis should be a space"<<*ctemp<<"is there a space"<<endl;
		
		//endiantest
		//bool swapbytes = false;
		
		
			
		
		getline(f,temp);//get rid of endl;
		//		cout<<"this should be blank "<<temp<<endl;
		getline(f,temp);
		getline(f,temp);//this is polydata line
		//	cout<<temp<<endl;

		int totPoints;
			f>>temp>>totPoints;//read in number of points
				getline(f,temp); //assumes float, so throws away rest of line
				//	cout<<temp<<endl;
				
				//read points

	
	
			vector<float> points;
		vector<int> polys;
		
		//read points
		//		cout<<"number of points: "<<totPoints<<endl;
		for (int i=0;i<3*totPoints;i++){
			f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
			if (swapbytes){
							Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);		
									}
			points.push_back(tempFloat);
		}
		f>>temp>>tempInt;//stores number of polygons
		//			cout<<temp<<endl;
		getline(f,temp); 
		//	cout<<temp<<endl;
		int numPolys=tempInt;
		for (int i=0;i<numPolys;i++){
			f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
			//throw "3" away, we knwo we're dealing with triangles
			f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
			if (swapbytes){
							Swap_Nbytes(1,sizeof(tempInt),&tempInt);		
									}
			polys.push_back(tempInt);
			f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
			if (swapbytes){
							Swap_Nbytes(1,sizeof(tempInt),&tempInt);		
									}
			polys.push_back(tempInt);
			f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
			if (swapbytes){
							Swap_Nbytes(1,sizeof(tempInt),&tempInt);		
									}
			polys.push_back(tempInt);
		}
		
		//read FIELD Data
		int numField; 
		f>>temp>>temp>>numField;
		getline(f,temp);
		
		//		cout<<"should be blank: "<<temp<<endl;//this should be 
	//int ipp=13;
		//this read all fields store modes seperately to divide by shape at the end
		vector< vector<Vec> > vmodes; 
		//vector< vector<float> > vimodes; 
		vector<float> imean_all; 
		vector<float> verrs;//only use a single set
			//	cout<<"numFiedl "<<numField<<endl;
		for (int i=0;i<numField;i++){
			f>>temp;
				//	cout<<temp<<endl;
			if (!strcmp(temp.substr(0,4).c_str(),"mode")){
			  //			cout<<"found "<<temp<<endl;
				string temp2;
			//	f>>tempInt>>tempInt>>temp2;
				getline(f,temp);//need getlines before reads
				vector<Vec> V;
					for (int j=0;j<totPoints;j++){
						f.read(reinterpret_cast<char*>(&vX),sizeof(float));
						if (swapbytes){
								Swap_Nbytes(1,sizeof(vX),&vX);	
									}
						f.read(reinterpret_cast<char*>(&vY),sizeof(float));
						if (swapbytes){
								Swap_Nbytes(1,sizeof(vY),&vY);	
									}
						f.read(reinterpret_cast<char*>(&vZ),sizeof(float));
						if (swapbytes){
								Swap_Nbytes(1,sizeof(vZ),&vZ);	
									}
						Vec vtmp(vX,vY,vZ);
						V.push_back(vtmp);
					}
					vmodes.push_back(V);
				
				
			}else if(!strcmp(temp.substr(0,5).c_str(),"Imean")){
			  //			cout<<"found "<<temp<<endl;
						int ipp; 
						string temp2;
						f>>ipp>>tempInt>>temp2;
							getline(f,temp); 
							//						cout<<"This should be empty "<<temp<<endl;
							//	cout<<"IMEAN "<<ipp<<" "<<tempInt<<" "<<temp2<<endl;
							//	cout<<"totpoint "<<totPoints<<endl;
						//	for (int j=0;j<totPoints;j++){
						for (int i=0;i<getNumberOfShapes();i++){
						vector<float> Vi;
							for (int j=0;j<npts.at(i);j++){

								for (int k=0;k<ipp;k++){
									//fV>>tempFloat;
									//cout<<"imean read"<<endl;
									f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
									if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
									Vi.push_back(tempFloat);
								}
								}
								shapes.at(i)->setIMean(Vi);
							}
						//				cout<<"done Imean"<<endl;
				}else if (!strcmp(temp.substr(0,5).c_str(),"Imode")){
			  //		cout<<"found "<<temp<<endl;
				string temp2;
				int ipp;
				f>>ipp>>tempInt>>temp2;
				getline(f,temp);
				for (int i=0;i<numShapes;i++){
					vector<float> Vi;
					for (int j=0;j<npts.at(i);j++){
						for (int k=0;k<ipp;k++){
							f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
							if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
							Vi.push_back(tempFloat);
						}
					}
					shapes.at(i)->addIModeVector(Vi);
					//vimodes.push_back(Vi);
					
				}
			}else if (!strcmp(temp.substr(0,5).c_str(),"iCond")){
				//cout<<"iconds found"<<endl;
				int sh;
				sh = atoi(temp.substr(9,1).c_str());
				if (!strcmp(temp.substr(0,9).c_str(),"iCondPrec")){
					string temp2;
					int tempInt2;
					f>>tempInt2>>tempInt>>temp2;
					getline(f,temp);
					vector< vector<float> > prectemp;							
					for (int i=0;i<tempInt;i++){
						vector<float> col;
						for (int j=0;j<tempInt2;j++){
							//fV>>tempFloat;
							f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
							if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
							col.push_back(tempFloat);
						}
						prectemp.push_back(col);
					}
					shapes.at(sh)->setICondPrec(prectemp);
				}else if (!strcmp(temp.substr(0,9).c_str(),"iCondEigs")){
					//cout<<"enter icondeigs"<<endl;
					string temp2;
					vector<float> eigstemp;
					f>>tempInt>>tempInt>>temp2;
					getline(f,temp);
					for (int j=0;j<tempInt;j++){
						f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
						if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
						eigstemp.push_back(tempFloat);
					}
					shapes.at(sh)->setICondEigs(eigstemp);
				}
				//			cout<<"iconds loaded"<<endl;
			}else if ((!strcmp(temp.substr(0,5).c_str(),"bCond"))&&(type>=1)){
			  //			cout<<"found bconds "<<temp<<endl;
				int sh;
				sh = atoi(temp.substr(9,1).c_str());
				if (!strcmp(temp.substr(0,9).c_str(),"bCondPrec")){
					string temp2;
					int tempInt2;
					f>>tempInt2>>tempInt>>temp2;	
					getline(f,temp);				
					vector< vector<float> > prectemp;							
					for (int i=0;i<tempInt;i++){
						vector<float> col;
						for (int j=0;j<tempInt2;j++){
							f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
							if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
							col.push_back(tempFloat);
						}
						prectemp.push_back(col);
					}
					shapes.at(sh)->setBCondPrec(prectemp);
				}else if (!strcmp(temp.substr(0,9).c_str(),"bCondEigs")){
				  //				cout<<"enter bcondeigs"<<endl;
					string temp2;
					vector<float> eigstemp;
					f>>tempInt>>tempInt>>temp2;
					getline(f,temp);
					for (int j=0;j<tempInt;j++){
						f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
						if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
						eigstemp.push_back(tempFloat);
					}
					shapes.at(sh)->setBCondEigs(eigstemp);
				}
			}
			//for info field which are already read in remove the storage, allows for rearranging of fields
			else if (!strcmp(temp.c_str(),"numPoints")){
			  //			cout<<"found numPoints"<<endl;
				f>>numShapes>>numShapes>>temp;
				getline(f,temp);
				for (int i=0;i<numShapes;i++){
					f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
				}
			}else if (!strcmp(temp.c_str(),"eigenValues")){
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				//set the number of modes by the number of eigenvalues
				numModes=tempInt2-1;
				for (int i=0;i<(tempInt2-1);i++){
					f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
					eigenVals.push_back(tempFloat);
				}
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
			}else if (!strcmp(temp.c_str(),"labels")){
				
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				for (int i=0;i<tempInt2;i++){
					f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
				}
			}else if (!strcmp(temp.c_str(),"numSubjects")){
				
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);

				f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
		
			}else if (!strcmp(temp.c_str(),"IPP")){
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				for (int i=0;i<tempInt2;i++){
					f.read(reinterpret_cast<char*>(&tempInt),sizeof(int));
				}
			}else if ((!strcmp(temp.c_str(),"eigenValuesI"))&&(type>=1)){
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				for (int i=0;i<(tempInt2-1);i++){
					f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				}
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
			
			}else if (!strcmp(temp.substr(0,9).c_str(),"ErrPriors")){
				int sh;
				sh = atoi(temp.substr(9,1).c_str());
				string temp2;
				f>>tempInt>>tempInt>>temp2;
				getline(f,temp);
				for (int j=0;j<2;j++){ //this allows for error in intensity and shape, for sh=0 its is the predictor error
					f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
						if (swapbytes){
							Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
						}
						verrs.push_back(tempFloat);
				}
				shapes.at(sh)->setErrs(verrs);
			}else if (!strcmp(temp.c_str(),"transVec")){
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
			//	cout<<tempInt2<<" "<<temp<<endl;
				vector<float> tx;
				vector<float> ty;
				vector<float> tz;
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				tx.push_back(tempFloat);
				//cout<<"transVec "<<tempFloat<<endl;
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				tx.push_back(tempFloat);
								//cout<<"transVec "<<tempFloat<<endl;

				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				tx.push_back(tempFloat);
								//cout<<"transVec "<<tempFloat<<endl;

				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				ty.push_back(tempFloat);
								//cout<<"transVec "<<tempFloat<<endl;

				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				ty.push_back(tempFloat);
								//cout<<"transVec "<<tempFloat<<endl;

				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				ty.push_back(tempFloat);
								
								//cout<<"transVec "<<tempFloat<<endl;

				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				tz.push_back(tempFloat);
								//cout<<"transVec "<<tempFloat<<endl;

				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				tz.push_back(tempFloat);
								//cout<<"transVec "<<tempFloat<<endl;

				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				tz.push_back(tempFloat);
								//cout<<"transVec "<<tempFloat<<endl;

				shapes.at(0)->setAffModeVectorsTr(tx,ty,tz);
			}else if (!strcmp(temp.c_str(),"transEigs")){
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				float tx;
				f.read(reinterpret_cast<char*>(&tx),sizeof(float));
				float ty;
				f.read(reinterpret_cast<char*>(&ty),sizeof(float));
				float tz;
				f.read(reinterpret_cast<char*>(&tz),sizeof(float));
				shapes.at(0)->setAffEigsTr(tx,ty,tz);
								//cout<<"tE "<<tx<<" "<<ty<<" "<<tz<<endl;

			}else if (!strcmp(temp.c_str(),"transMean")){
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				float tx;
				f.read(reinterpret_cast<char*>(&tx),sizeof(float));
				float ty;
				f.read(reinterpret_cast<char*>(&ty),sizeof(float));
				float tz;
				f.read(reinterpret_cast<char*>(&tz),sizeof(float));
				shapes.at(0)->setAffMeanTr(tx,ty,tz);
								//cout<<"tM "<<tx<<" "<<ty<<" "<<tz<<endl;

			}else if (!strcmp(temp.c_str(),"RotationVecs")){
				//cout<<"tmep "<<temp<<endl;
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				vector<float> rx;
				vector<float> ry;
				vector<float> rz;
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				rx.push_back(tempFloat);
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				rx.push_back(tempFloat);
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				rx.push_back(tempFloat);
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				ry.push_back(tempFloat);
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				ry.push_back(tempFloat);
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				ry.push_back(tempFloat);
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				rz.push_back(tempFloat);
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				rz.push_back(tempFloat);
				f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
				rz.push_back(tempFloat);
								shapes.at(0)->setAffModeVectorsRot(rx,ry,rz);

			}else if (!strcmp(temp.c_str(),"RotationEigs")){
			//cout<<"tmep "<<temp<<endl;
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				float tx;
				f.read(reinterpret_cast<char*>(&tx),sizeof(float));
				float ty;
				f.read(reinterpret_cast<char*>(&ty),sizeof(float));
				float tz;
				f.read(reinterpret_cast<char*>(&tz),sizeof(float));
								shapes.at(0)->setAffEigsRot(tx,ty,tz);
				//cout<<"rotE "<<tx<<" "<<ty<<" "<<tz<<endl;

			}else if (!strcmp(temp.c_str(),"RotationMean")){
			//cout<<"tmep "<<temp<<endl;
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				float tx;
				f.read(reinterpret_cast<char*>(&tx),sizeof(float));
				float ty;
				f.read(reinterpret_cast<char*>(&ty),sizeof(float));
				float tz;
				f.read(reinterpret_cast<char*>(&tz),sizeof(float));
				shapes.at(0)->setAffMeanRot(tx,ty,tz);
				//cout<<"rotM "<<tx<<" "<<ty<<" "<<tz<<endl;

			}else if (!strcmp(temp.c_str(),"scMeanVar")){
			//cout<<"tmep "<<temp<<endl;
				f>>tempInt2>>tempInt2>>temp;
				getline(f,temp);
				float tmean;
				f.read(reinterpret_cast<char*>(&tmean),sizeof(float));
				float tvar;
				f.read(reinterpret_cast<char*>(&tvar),sizeof(float));
				shapes.at(0)->setAffEigsSc(tvar);
				shapes.at(0)->setAffMeanSc(tmean);
				//cout<<"sc "<<tvar<<" "<<tmean<<endl;
			}else if (!strcmp(temp.substr(0,8).c_str(),"AffImode")){
			  	//	cout<<"found "<<temp<<endl;
				string temp2;
				int ipp;
				f>>ipp>>tempInt>>temp2;
				getline(f,temp);
				for (int i=0;i<numShapes;i++){
					vector<float> Vi;
					for (int j=0;j<npts.at(i);j++){
						for (int k=0;k<ipp;k++){
							f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
							if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
							Vi.push_back(tempFloat);
						}
					}
					shapes.at(i)->addAffIModeVector(Vi);
					//vimodes.push_back(Vi);
					
				}
			}else if (!strcmp(temp.substr(0,8).c_str(),"iAffCond")){
				//cout<<"iAffconds found"<<endl;
				int sh;
				sh = atoi(temp.substr(9,1).c_str());
				if (!strcmp(temp.substr(0,12).c_str(),"iAffCondPrec")){
					//cout<<"found AffCOndI"<<endl;
					string temp2;
					int tempInt2;
					f>>tempInt2>>tempInt>>temp2;
					getline(f,temp);
					vector< vector<float> > prectemp;							
					for (int i=0;i<tempInt;i++){
						vector<float> col;
						for (int j=0;j<tempInt2;j++){
							//fV>>tempFloat;
							f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
							if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
							col.push_back(tempFloat);
						}
						prectemp.push_back(col);
					}
					shapes.at(sh)->setICondPrecAff(prectemp);
				}else if (!strcmp(temp.substr(0,12).c_str(),"iAffCondEigs")){
				//cout<<"found AffCOndIEiug"<<endl;
					//cout<<"enter icondeigs"<<endl;
					string temp2;
					vector<float> eigstemp;
					f>>tempInt>>tempInt>>temp2;
					getline(f,temp);
					for (int j=0;j<tempInt;j++){
						f.read(reinterpret_cast<char*>(&tempFloat),sizeof(float));
						if (swapbytes){
									Swap_Nbytes(1,sizeof(tempFloat),&tempFloat);
									}
						eigstemp.push_back(tempFloat);
					}
					shapes.at(sh)->setICondEigsAff(eigstemp);
				}
				//			cout<<"iconds loaded"<<endl;
			}

			
			
		}
					
						
	f.close();
//	cout<<"constructing shape Model"<<endl;
//	//Data is read now set up shape model
//	//craete shape objects
///	for (int i=0;i<numShapes;i++){
	//	MShape*  Imshape = new MShape();
//		shapes.push_back(Imshape);
//		if (type==1){
//			shapes.at(i)->setIPP(ipps.at(i));
//		}
//	}
	//load in Eignenvalues
	
	
	
	//	cout<<"shape added"<<endl;
	vector<int> cumnum;
	cumnum.push_back(0);
	for (unsigned int i=0;i<npts.size();i++){
		
	}
	
	
	for (int sh=0;sh<numShapes;sh++){//create new mesh for each shape
	  //			cout<<"Shape "<<sh<<endl;	
								
		for (int i=0; i<npts.at(sh); i++)
		{//assign points to mesh
			int offset=0;
			if (sh>0){
				offset=cumnum.at(sh);
			}
				Mpoint * p = new Mpoint(points.at(3*i+offset), points.at(3*i+1+offset), points.at(3*i+2+offset), i);
			
			shapes.at(sh)->pushPoint(p);
		}
		//		cout<<"points added"<<endl;
		//reading the triangles
		
		bool fpass=true;
		int cumPts=0, Pcount=0;    
		int p0=0, p1=0, p2=0;

		while (true) //break when an poitn index is too large 
		{
			Pcount++;
		     
			if(Pcount>numPolys){
				break;
			}
			if ((sh>0)&&(fpass)){//a point carried over from previous shape
								 //			cout<<"npts"<<npts.at(sh-1)<<" "<<p0<<" "<<p1<<" "<<p2<<endl;
				Triangle * t = new Triangle(shapes.at(sh)->getPoint(p0-npts.at(sh-1)), shapes.at(sh)->getPoint(p1-npts.at(sh-1)), shapes.at(sh)->getPoint(p2-npts.at(sh-1)));
				//			cout<<Pcount<<endl;
				shapes.at(sh)->pushTriangle(t);
				//			cout<<Pcount<<"2"<<endl;
			}
			
			fpass=false;
			//cout<<Pcount<<" "<<polys.size()<<endl;
			p0=polys.at(3*(Pcount-1));
			p1=polys.at(3*(Pcount-1)+1);
			p2=polys.at(3*(Pcount-1)+2);
			//cout<<"P: "<<p0<<" "<<p1<<" "<<p2<<" "<<tempInt<<endl;
			p0-=cumPts;
			p1-=cumPts;
			p2-=cumPts;
			
			if ((p0>=(npts.at(sh)))||(p1>=(npts.at(sh)))||(p2>=(npts.at(sh)))||(Pcount>numPolys)){
				//looks for end of shape or end of all shapes
				break;
			}
			Triangle * t = new Triangle(shapes.at(sh)->getPoint(p0), shapes.at(sh)->getPoint(p1), shapes.at(sh)->getPoint(p2));
			shapes.at(sh)->pushTriangle(t);
			
		}
		
		//		cout<<"polygons added"<<endl;
		cumPts+=npts.at(sh);
	}
	//int nummodes=static_cast<int>(vmodes.size());
	//	cout<<"numModes "<<numModes<<endl;
	for (int i =0;i<numModes;i++){
		for (int sh=0;sh<numShapes;sh++){
			vector< Vec > V;
			for (int j=cumnum.at(sh);j<cumnum.at(sh)+npts.at(sh);j++){
				V.push_back(vmodes.at(i).at(j));
			}
			shapes.at(sh)->addModeVector(V);
			
		}
	}

				return 1;
	}

	
	
	int shapeModel::getNumberOfModes() {
		return shapes.at(0)->getNumberOfModes();
	}
	
	int shapeModel::getNumberOfShapes(){
		return shapes.size();
	}
	int shapeModel::getNumberOfSubjects(){
		return numSubjects;
	}
  void shapeModel::setNumberOfSubjects(int N) {
		numSubjects=N;
	}
	int shapeModel::getNumberOfPoints(int shape){
		return shapes.at(shape)->getNumberOfPoints();
	}	
 
	int shapeModel::getTotalNumberOfPoints(){
		int N=0;
		for (int i=0;i<getNumberOfShapes();i++){
			N+=shapes.at(i)->getNumberOfPoints();
			//cout<<"N: "<<N<<endl;
			}
			return N;
	}
	int shapeModel::getLabel(int shape){
		return labels.at(shape);
	}
	int shapeModel::getIPP(int shape){
		return ipps.at(shape);
	}
	vector<int> shapeModel::getLabels(){
		return labels;
	}
	
	void shapeModel::setLabels(vector<int> labs){
		labels=labs;
	}
	float shapeModel::getEigenValue(int mode){ 
		return eigenVals.at(mode);
	}
	float shapeModel::getEigenValueI(int mode){ 
		return eigenValsI.at(mode);
	}
	
	void shapeModel::setEigenValues(vector<float> vals){ 
		eigenVals=vals;
		//set the sum of the eigenvalues
		sumEigVals=0;
		for (int i=0;i<static_cast<int>(vals.size());i++){
			sumEigVals+=eigenVals.at(i);
		}
	}
	void shapeModel::setEigenValuesI(vector<float> vals){ 
		eigenValsI=vals;
		//set the sum of the eigenvalues
		sumEigValsI=0;
		for (int i=0;i<static_cast<int>(vals.size());i++){
			sumEigValsI+=eigenValsI.at(i);
		}
	}
	void shapeModel::setIntensityProfiles(vector< vector<float> > profiles){
		iprofs=profiles;
	}
	vector< vector<float> > shapeModel::getIntensityProfiles(){
		return iprofs;
	}
	float shapeModel::getSumEigenValues(){ 
		
		return sumEigVals;
	}
	void shapeModel::setICondPrecMatrix(vector< vector<float> > precmat, vector<float> Eigs, int shape){
		shapes.at(shape)->setICondPrecEigs(precmat, Eigs );
	}

	
	
	Mesh shapeModel::getDeformedMesh( vector<float> var, int shape, int numModes){
		//get shape mesh shapes.at(shape);
	//cout<<"get shape"
		Mesh m = getShapeMesh(shape);//shapes.at(shape)->getMesh();
	
		//cout<<"varsize "<<var.size()<<endl;
		for (unsigned int mode=0; mode<var.size();mode++){
			//for each mode
		//	cout<<"get shape mode "<<mode<<endl;
			const vector<Vec> mvec = shapes.at(shape)->getModeVector(mode);
	//	cout<<"mode got"<<endl;
			int count=0;
			
			for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
				(*i)->_update_coord = (*i)->get_coord() + (var.at(mode)*sqrt(eigenVals.at(mode))*mvec.at(count));
				count++;
			}
		//	cout<<"reached end p"<<endl;
			m.update();
		}
		//translates mesh to centre of image
		Vec trans((xsize-1)/2.0*abs(xdim),(ysize-1)/2.0*abs(ydim),(zsize-1)/2.0*abs(zdim));
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			(*i)->_update_coord = (*i)->get_coord() + trans;
		}
		
		m.update();
		return m;
	}
	
	
	Mesh shapeModel::getDeformedMeshAff6( vector<float> var, int shape, int numModes){
		//this function assumes that the first 7 modes are the affine parameters
		///they should be used in the specific order
		//must have at least 7 modes
		

		//get shape mesh shapes.at(shape);
		Mesh m = getShapeMesh(shape);//shapes.at(shape)->getMesh();
//cout<<"vars size "<<var.size()<<endl;
		//get the translation vectors
		vector< float > tr1=shapes.at(0)->getAfftxVecs();
		vector< float > tr2=shapes.at(0)->getAfftyVecs();
		vector< float > tr3=shapes.at(0)->getAfftzVecs();
		float Etr1=shapes.at(0)->getAfftxEig();
		float Etr2=shapes.at(0)->getAfftyEig();
		float Etr3=shapes.at(0)->getAfftzEig();
		//cout<<"got translation paarametersw "<<endl;
		//cout<<var.at(0)<<" "<<var.at(1)<<" "<<var.at(2)<<" "<<var.at(3)<<" "<<var.at(4)<<" "<<var.at(5)<<" "<<var.at(6)<<" "<<" "<<endl;
	//	cout<<Etr1<<" "<<Etr2<<" "<<Etr3<<" "<< shapes.at(0)->getAffrxEig()<<" "<< shapes.at(0)->getAffryEig()<<" "<< shapes.at(0)->getAffrzEig()<<endl;
	//	cout<<shapes.at(0)->getAfftxMean()<<" "<<shapes.at(0)->getAfftyMean()<<" "<< shapes.at(0)->getAfftzMean()<<" "<<shapes.at(0)->getAffrxMean()<<" "<<shapes.at(0)->getAffryMean()<<" "<< shapes.at(0)->getAffrzMean()<<endl;
		
		Vec trOffSet(var.at(0)*sqrt(Etr1)*tr1.at(0)+var.at(1)*sqrt(Etr2)*tr2.at(0)+var.at(2)*sqrt(Etr3)*tr3.at(0) + shapes.at(0)->getAfftxMean() , \
					 var.at(0)*sqrt(Etr1)*tr1.at(1)+var.at(1)*sqrt(Etr2)*tr2.at(1)+var.at(2)*sqrt(Etr3)*tr3.at(1) + shapes.at(0)->getAfftyMean(), \
					 var.at(0)*sqrt(Etr1)*tr1.at(2)+var.at(1)*sqrt(Etr2)*tr2.at(2)+var.at(2)*sqrt(Etr3)*tr3.at(2) + shapes.at(0)->getAfftzMean()); 
		//cout<<"trOffSet "<<trOffSet.X<<" "<<trOffSet.Y<<" "<<trOffSet.Z<<endl;
		
		Pt MeanPt;
		//calculate centroid
		int count=0;
		
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			MeanPt+=(*i)->get_coord() ;
			count++;
		}
		MeanPt/=static_cast<float>(count);
	//	cout<<"centroid "<<MeanPt.X<<" "<<MeanPt.Y<<" "<<MeanPt.Z<<endl;
		
		
		//calculate scale
	//	float scale =var.at(6)*sqrt( shapes.at(0)->getAffscEig())+ shapes.at(0)->getAffscMean();
	float scale=1;
	//	cout<<"scale "<<scale<<endl;
		//calculate rotation
		vector< float > rr1= shapes.at(0)->getAffrxVecs();
		vector< float > rr2= shapes.at(0)->getAffryVecs();
		vector< float > rr3= shapes.at(0)->getAffrzVecs();
		float rx=var.at(3)*sqrt( shapes.at(0)->getAffrxEig())*rr1.at(0)+var.at(4)*sqrt( shapes.at(0)->getAffryEig())*rr2.at(0)+var.at(5)*sqrt( shapes.at(0)->getAffrzEig())*rr3.at(0) + shapes.at(0)->getAffrxMean();
		float ry=var.at(3)*sqrt( shapes.at(0)->getAffrxEig())*rr1.at(1)+var.at(4)*sqrt( shapes.at(0)->getAffryEig())*rr2.at(1)+var.at(5)*sqrt( shapes.at(0)->getAffrzEig())*rr3.at(1) + shapes.at(0)->getAffryMean();
		float rz=var.at(3)*sqrt( shapes.at(0)->getAffrxEig())*rr1.at(2)+var.at(4)*sqrt( shapes.at(0)->getAffryEig())*rr2.at(2)+var.at(5)*sqrt( shapes.at(0)->getAffrzEig())*rr3.at(2) + shapes.at(0)->getAffrzMean();

	//	cout<<"Angles "<<rx<<" "<<ry<<" "<<rz<<endl;
		//calculate rotation matrix
		float cX=cos(rx);
		float sX=sin(rx);
		float cY=cos(ry);
		float sY=sin(ry);
		float cZ=cos(rz);
		float sZ=sin(rz);


		//calculate rotation matrix
		float r11=cY*cZ;
		float r12= cY*sZ;
		float r13=-sY;
		float r21=sX*sY*cZ-cX*sZ;
		float r22=sX*sY*sZ+cX*cZ;
		float r23=sX*cY;
		float r31=cX*sY*cZ+sX*sZ;
		float r32=cX*sY*sZ-sX*cZ;
		float r33= cX*cY;
		
	//	m.rotation(cY*cZ, cY*sZ, -sY, sX*sY*cZ-cX*sZ, sX*sY*sZ+cX*cZ , sX*cY , cX*sY*cZ+sX*sZ, cX*sY*sZ-sX*cZ, cX*cY, MeanPt.X,MeanPt.Y,MeanPt.Z);
	//	m.update();
		//m.rescale(scale, MeanPt.X,MeanPt.Y,MeanPt.Z);
		//m.update();
		//m.translation(trOffSet);
	//	m.update();
	//	cout<<"varsize "<<var.size()-7<<endl;
	
	
	
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
				Vec cen= (*i)->get_coord() - Pt(MeanPt.X,MeanPt.Y,MeanPt.Z);
				(*i)->_update_coord = Pt(MeanPt.X,MeanPt.Y,MeanPt.Z) + scale*Vec( (cen.X*r11+cen.Y*r12+cen.Z*r13) , (cen.X*r21+cen.Y*r22+cen.Z*r23) , (cen.X*r31+cen.Y*r32+cen.Z*r33)) + trOffSet;
			//	(*i)->_update_coord = Pt(MeanPt.X,MeanPt.Y,MeanPt.Z) +scale*cen + trOffSet;

				//	cout<<(*i)->_update_coord.X<<endl;
				count++;
			}
	
		m.update();
		for (unsigned int mode=0; mode<(var.size()-6);mode++){
			//for each mode
		//cout<<"get shape mode "<<mode<<endl;
			const vector<Vec> mvec = shapes.at(shape)->getModeVector(mode);
		//	cout<<"mode got"<<endl;
			int count=0;
			
			for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
				(*i)->_update_coord = (*i)->get_coord() + (var.at(mode+6)*sqrt(eigenVals.at(mode))*mvec.at(count));
				//	cout<<(*i)->_update_coord.X<<endl;
				count++;
			}
		//	cout<<"reached end p"<<endl;
			m.update();
		}
		//cout<<"done deform1"<<endl;
		//translates mesh to centre of image
		Vec trans((xsize-1)/2.0*abs(xdim),(ysize-1)/2.0*abs(ydim),(zsize-1)/2.0*abs(zdim));
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			(*i)->_update_coord = (*i)->get_coord() + trans;
		//	cout<<(*i)->_update_coord.X<<endl;
		}
		
		m.update();
			//	cout<<"done deform2"<<endl;

		return m;
	}
		
	Mesh shapeModel::getDeformedMeshAff7( vector<float> var, int shape, int numModes){
		//this function assumes that the first 7 modes are the affine parameters
		///they should be used in the specific order
		//must have at least 7 modes
		

		//get shape mesh shapes.at(shape);
		Mesh m = getShapeMesh(shape);//shapes.at(shape)->getMesh();

		//get the translation vectors
		vector< float > tr1=shapes.at(0)->getAfftxVecs();
		vector< float > tr2=shapes.at(0)->getAfftyVecs();
		vector< float > tr3=shapes.at(0)->getAfftzVecs();
		float Etr1=shapes.at(0)->getAfftxEig();
		float Etr2=shapes.at(0)->getAfftyEig();
		float Etr3=shapes.at(0)->getAfftzEig();
		
	//	cout<<var.at(0)<<" "<<var.at(1)<<" "<<var.at(1)<<" "<<var.at(2)<<" "<<var.at(3)<<" "<<var.at(4)<<" "<<var.at(5)<<" "<<var.at(6)<<" "<<endl;
	//	cout<<Etr1<<" "<<Etr2<<" "<<Etr3<<" "<< shapes.at(0)->getAffrxEig()<<" "<< shapes.at(0)->getAffryEig()<<" "<< shapes.at(0)->getAffrzEig()<<endl;
	//	cout<<shapes.at(0)->getAfftxMean()<<" "<<shapes.at(0)->getAfftyMean()<<" "<< shapes.at(0)->getAfftzMean()<<" "<<shapes.at(0)->getAffrxMean()<<" "<<shapes.at(0)->getAffryMean()<<" "<< shapes.at(0)->getAffrzMean()<<endl;
		
		Vec trOffSet(var.at(0)*sqrt(Etr1)*tr1.at(0)+var.at(1)*sqrt(Etr2)*tr2.at(0)+var.at(2)*sqrt(Etr3)*tr3.at(0) + shapes.at(0)->getAfftxMean() , \
					 var.at(0)*sqrt(Etr1)*tr1.at(1)+var.at(1)*sqrt(Etr2)*tr2.at(1)+var.at(2)*sqrt(Etr3)*tr3.at(1) + shapes.at(0)->getAfftyMean(), \
					 var.at(0)*sqrt(Etr1)*tr1.at(2)+var.at(1)*sqrt(Etr2)*tr2.at(2)+var.at(2)*sqrt(Etr3)*tr3.at(2) + shapes.at(0)->getAfftzMean()); 
		//cout<<"trOffSet "<<trOffSet.X<<" "<<trOffSet.Y<<" "<<trOffSet.Z<<endl;
		
		Pt MeanPt;
		//calculate centroid
		int count=0;
		
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			MeanPt+=(*i)->get_coord() ;
			count++;
		}
		MeanPt/=static_cast<float>(count);
	//	cout<<"centroid "<<MeanPt.X<<" "<<MeanPt.Y<<" "<<MeanPt.Z<<endl;
		
		
		//calculate scale
		float scale =var.at(6)*sqrt( shapes.at(0)->getAffscEig())+ shapes.at(0)->getAffscMean();
	//	cout<<"scale "<<scale<<endl;
		//calculate rotation
		vector< float > rr1= shapes.at(0)->getAffrxVecs();
		vector< float > rr2= shapes.at(0)->getAffryVecs();
		vector< float > rr3= shapes.at(0)->getAffrzVecs();
		float rx=var.at(3)*sqrt( shapes.at(0)->getAffrxEig())*rr1.at(0)+var.at(4)*sqrt( shapes.at(0)->getAffryEig())*rr2.at(0)+var.at(5)*sqrt( shapes.at(0)->getAffrzEig())*rr3.at(0) + shapes.at(0)->getAffrxMean();
		float ry=var.at(3)*sqrt( shapes.at(0)->getAffrxEig())*rr1.at(1)+var.at(4)*sqrt( shapes.at(0)->getAffryEig())*rr2.at(1)+var.at(5)*sqrt( shapes.at(0)->getAffrzEig())*rr3.at(1) + shapes.at(0)->getAffryMean();
		float rz=var.at(3)*sqrt( shapes.at(0)->getAffrxEig())*rr1.at(2)+var.at(4)*sqrt( shapes.at(0)->getAffryEig())*rr2.at(2)+var.at(5)*sqrt( shapes.at(0)->getAffrzEig())*rr3.at(2) + shapes.at(0)->getAffrzMean();

	//	cout<<"Angles "<<rx<<" "<<ry<<" "<<rz<<endl;
		//calculate rotation matrix
		float cX=cos(rx);
		float sX=sin(rx);
		float cY=cos(ry);
		float sY=sin(ry);
		float cZ=cos(rz);
		float sZ=sin(rz);


		//calculate rotation matrix
		float r11=cY*cZ;
		float r12= cY*sZ;
		float r13=-sY;
		float r21=sX*sY*cZ-cX*sZ;
		float r22=sX*sY*sZ+cX*cZ;
		float r23=sX*cY;
		float r31=cX*sY*cZ+sX*sZ;
		float r32=cX*sY*sZ-sX*cZ;
		float r33= cX*cY;
		
	//	m.rotation(cY*cZ, cY*sZ, -sY, sX*sY*cZ-cX*sZ, sX*sY*sZ+cX*cZ , sX*cY , cX*sY*cZ+sX*sZ, cX*sY*sZ-sX*cZ, cX*cY, MeanPt.X,MeanPt.Y,MeanPt.Z);
	//	m.update();
		//m.rescale(scale, MeanPt.X,MeanPt.Y,MeanPt.Z);
		//m.update();
		//m.translation(trOffSet);
	//	m.update();
	//	cout<<"varsize "<<var.size()-7<<endl;
	
	
	
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
				Vec cen= (*i)->get_coord() - Pt(MeanPt.X,MeanPt.Y,MeanPt.Z);
				(*i)->_update_coord = Pt(MeanPt.X,MeanPt.Y,MeanPt.Z) + scale*Vec( (cen.X*r11+cen.Y*r12+cen.Z*r13) , (cen.X*r21+cen.Y*r22+cen.Z*r23) , (cen.X*r31+cen.Y*r32+cen.Z*r33)) + trOffSet;
			//	(*i)->_update_coord = Pt(MeanPt.X,MeanPt.Y,MeanPt.Z) +scale*cen + trOffSet;

				//	cout<<(*i)->_update_coord.X<<endl;
				count++;
			}
	
		m.update();
		for (unsigned int mode=0; mode<(var.size()-7);mode++){
			//for each mode
	//		cout<<"get shape mode "<<mode<<endl;
			const vector<Vec> mvec = shapes.at(shape)->getModeVector(mode);
	//	cout<<"mode got"<<endl;
			int count=0;
			
			for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
				(*i)->_update_coord = (*i)->get_coord() + (var.at(mode+7)*sqrt(eigenVals.at(mode))*mvec.at(count));
				//	cout<<(*i)->_update_coord.X<<endl;
				count++;
			}
		//	cout<<"reached end p"<<endl;
			m.update();
		}
		//cout<<"done deform1"<<endl;
		//translates mesh to centre of image
		Vec trans((xsize-1)/2.0*abs(xdim),(ysize-1)/2.0*abs(ydim),(zsize-1)/2.0*abs(zdim));
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			(*i)->_update_coord = (*i)->get_coord() + trans;
		//	cout<<(*i)->_update_coord.X<<endl;
		}
		
		m.update();
			//	cout<<"done deform2"<<endl;

		return m;
	}

	
	Mesh shapeModel::getDeformedMeshP3dof( vector<float> var, int shape, int numModes,float tx, float ty, float tz){
		//get shape mesh shapes.at(shape);
	//cout<<"get shape"
		Mesh m = getShapeMesh(shape);//shapes.at(shape)->getMesh();
		Vec toffset(tx,ty,tz);
		//cout<<"varsize "<<var.size()<<endl;
		for (unsigned int mode=0; mode<var.size();mode++){
			//for each mode
		//	cout<<"get shape mode "<<mode<<endl;
			const vector<Vec> mvec = shapes.at(shape)->getModeVector(mode);
	//	cout<<"mode got"<<endl;
			int count=0;
			
			for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
				(*i)->_update_coord = (*i)->get_coord() + (var.at(mode)*sqrt(eigenVals.at(mode))*mvec.at(count))+toffset;
				count++;
			}
		//	cout<<"reached end p"<<endl;
			m.update();
		}
		//translates mesh to centre of image
		Vec trans((xsize-1)/2.0*abs(xdim),(ysize-1)/2.0*abs(ydim),(zsize-1)/2.0*abs(zdim));
		
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			(*i)->_update_coord = (*i)->get_coord() + trans;
		}
		
		m.update();
		return m;
	}

	
	
	vector<float> shapeModel::getDeformedIprof( vector<float> var, int shape, int numModes){
	
		vector<float> imean=shapes.at(shape)->getIMean();
		for (unsigned int mode=0; mode<var.size();mode++){
			//for each mode
			//cout<<"mode "<<mode<<" "<<var.at(mode)<<endl;
			//cout<<"shape "<<shape<<" "<<shapes.size()<<" "<<endl;
			vector<float> mvec = shapes.at(shape)->getIModeVector(mode);
		//	cout<<"helpiprof"<<endl;
			for (unsigned int i=0;i<imean.size();i++){
					//cout<<"helpiprof "<<i<<endl;
					imean.at(i)=imean.at(i)+var.at(mode)*sqrt(eigenVals.at(mode))*mvec.at(i);
			}
		//	cout<<"helpiprof2"<<endl;
		}
		return imean;
	
	}
		vector<float> shapeModel::getDeformedIprofAff7( vector<float> var, int shape, int numModes){
	
		vector<float> imean=shapes.at(shape)->getIMean();
		for (unsigned int mode=0; mode<var.size();mode++){
			//for each mode
			vector<float> mvec;
			//cout<<"mode "<<mode<<" "<<var.at(mode)<<endl;
			//cout<<"shape "<<shape<<" "<<shapes.size()<<" "<<endl;
			//should not have more than 7 modes anyways
			if (mode<7){
				mvec= shapes.at(shape)->getAffIModeVector(mode);
				//	cout<<"helpiprof"<<endl;
				float eig=0;
				switch (mode){
					case 0:
						eig=shapes.at(shape)->getAfftxEig();
						break;
						case 1:
						eig=shapes.at(shape)->getAfftyEig();
						break;
							case 2:
						eig=shapes.at(shape)->getAfftzEig();
						break;
							case 3:
						eig=shapes.at(shape)->getAffrxEig();
						break;
							case 4:
						eig=shapes.at(shape)->getAffryEig();
						break;
							case 5:
						eig=shapes.at(shape)->getAffrzEig();
						break;
							case 6:
						eig=shapes.at(shape)->getAffscEig();
						break;
					
				}
				for (unsigned int i=0;i<imean.size();i++){
					//cout<<"helpiprof "<<i<<endl;
					imean.at(i)=imean.at(i)+var.at(mode)*sqrt(eig)*mvec.at(i);
				}
		}
		
		//	cout<<"helpiprof2"<<endl;
		}
		return imean;
	
	}
vector<float> shapeModel::getDeformedIprofAff6( vector<float> var, int shape, int numModes){
	
		vector<float> imean=shapes.at(shape)->getIMean();
		for (unsigned int mode=0; mode<var.size();mode++){
			//for each mode
			vector<float> mvec;
			//cout<<"mode "<<mode<<" "<<var.at(mode)<<endl;
			//cout<<"shape "<<shape<<" "<<shapes.size()<<" "<<endl;
			//should not have more than 7 modes anyways
			if (mode<6){
				mvec= shapes.at(shape)->getAffIModeVector(mode);
				//	cout<<"helpiprof"<<endl;
				float eig=0;
				switch (mode){
					case 0:
						eig=shapes.at(shape)->getAfftxEig();
						break;
						case 1:
						eig=shapes.at(shape)->getAfftyEig();
						break;
							case 2:
						eig=shapes.at(shape)->getAfftzEig();
						break;
							case 3:
						eig=shapes.at(shape)->getAffrxEig();
						break;
							case 4:
						eig=shapes.at(shape)->getAffryEig();
						break;
							case 5:
						eig=shapes.at(shape)->getAffrzEig();
						break;
					
				}
				for (unsigned int i=0;i<imean.size();i++){
					//cout<<"helpiprof "<<i<<endl;
					imean.at(i)=imean.at(i) +var.at(mode)*sqrt(eig)*mvec.at(i);
				}
		}
		
		//	cout<<"helpiprof2"<<endl;
		}
		return imean;
	
	}

//	Mesh shapeModel::getTranslatedMesh( vector<float> var, int shape){
//		
//		//get shape mesh
//		Mesh m = shapes.at(shape)->getMesh();
//		
//		for (unsigned int mode=0; mode<var.size();mode++){
//			//for each mode
//			const vector<Vec> mvec = shapes.at(shape)->getModeVector(mode);
//			int count=0;
///			for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
	//			(*i)->_update_coord = (*i)->get_coord() + (var.at(mode)*sqrt(eigenVals.at(mode))*mvec.at(count));
//				count++;
//			}
//			m.update();
//		}
//		//translates mesh to centre of image
//		Vec trans((xsize-1)/2.0*abs(xdim),(ysize-1)/2.0*abs(ydim),(zsize-1)/2.0*abs(zdim));
//		cout<<"Xsize "<<xsize<<endl;
//		//trans.Y=trans.Y+var;
		//trans.Y=trans.Y+3.4;
//		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
//			(*i)->_update_coord = (*i)->get_coord() + trans;
//		}
		
//		m.update();
//		return m;
//	}
	Mesh shapeModel::getTranslatedMesh( int shape){
		
		//get shape mesh
		Mesh m = shapes.at(shape)->getMesh();
		
		//translates mesh to centre of image
		Vec trans((xsize-1)/2.0*abs(xdim),(ysize-1)/2.0*abs(ydim),(zsize-1)/2.0*abs(zdim));
		//cout<<"Xsize "<<xsize<<endl;
		//trans.Y=trans.Y+var;
		//trans.Y=trans.Y+3.4;
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			(*i)->_update_coord = (*i)->get_coord() + trans;
		}
		
		m.update();
		return m;
	}
	vector<float> shapeModel::getProjectModeParameters( Mesh m , int shape, int trunc, float* res ){
		vector<float> vars;
		vars.push_back(0);
		//get MeanMesh
		Mesh m1=getDeformedMesh(vars,shape,static_cast<int>(vars.size()));

		//translates mesh to centre of image
//		Vec trans(-(xsize-1)/2.0*abs(xdim),-(ysize-1)/2.0*abs(ydim),-(zsize-1)/2.0*abs(zdim));
		vector< float > vtrans;
		//load mean into vector
		for (vector<Mpoint*>::iterator i = m1._points.begin(); i!=m1._points.end(); i++ ){
			vtrans.push_back( (*i)->get_coord().X );
			vtrans.push_back( (*i)->get_coord().Y );
			vtrans.push_back( (*i)->get_coord().Z );
		}

		//trans.Y=trans.Y+var;
		//trans.Y=trans.Y+3.4;
		vector< float > vdmean;
		vector< float > vmesh;

		int count=0;
		//demean and store in vector
	//	cout<<"sp 1"<<endl;
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			vdmean.push_back((*i)->get_coord().X-vtrans.at(count));
			vdmean.push_back((*i)->get_coord().Y-vtrans.at(count+1));
			vdmean.push_back((*i)->get_coord().Z-vtrans.at(count+2));
			
			vmesh.push_back((*i)->get_coord().X);
			vmesh.push_back((*i)->get_coord().Y);
			vmesh.push_back((*i)->get_coord().Z);
			
			count+=3;
		}
	//			cout<<"sp 2"<<endl;

		
		//now project caculate D^(-1)U^tvdean
		vector<float> bvars;
	//	cout<<eigenVals.size()<<endl;
		for (int i=0;i<getNumberOfModes();i++){
	//						cout<<"sp 31 "<<i<<" "<<eigenVals.at(i)<<endl;

			vector< Vec > vmode = shapes.at(shape)->getModeVector(i);
			float vari=0;
			count=0;
	//				cout<<"sp 3 "<<i<<endl;
			if(i<trunc){
			for (unsigned int j=0;j<vmode.size();j++){
				vari+=vdmean.at(count)*vmode.at(j).X;
				vari+=vdmean.at(count+1)*vmode.at(j).Y;
				vari+=vdmean.at(count+2)*vmode.at(j).Z;
				count+=3;
			}
			}
			bvars.push_back(vari/sqrt(eigenVals.at(i)));
		}
		
		
		m1=getDeformedMesh(bvars,shape,static_cast<int>(bvars.size()));
	//	cout<<"lengthsa "<<vdmean.size()<<" "<<vtrans.size()<<" "<<bvars.size()<<endl;
		vector< float > vdif;
		 count=0;
		//demean and store in vector
	//	cout<<"sp 1"<<endl;
		float ssx=0;
		for (vector<Mpoint*>::iterator i = m1._points.begin(); i!=m1._points.end(); i++ ){
			float difx=vmesh.at(count)- ((*i)->get_coord().X);
			float dify=vmesh.at(count+1)- ((*i)->get_coord().Y);
			float difz=vmesh.at(count+2)- ((*i)->get_coord().Z);
			ssx+=difx*difx+dify*dify+difz*difz;
			count+=3;
		}

		*res=ssx;
		return bvars;
	}

	
	Mesh shapeModel::getInverseTranslatedMesh( Mesh m ){
		
		//translates mesh to centre of image
		Vec trans(-(xsize-1)/2.0*abs(xdim),-(ysize-1)/2.0*abs(ydim),-(zsize-1)/2.0*abs(zdim));
		//trans.Y=trans.Y+var;
		//trans.Y=trans.Y+3.4;
		for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
			(*i)->_update_coord = (*i)->get_coord() + trans;
		}
		
		m.update();
		return m;
	}
	
	
	void shapeModel::addShape(MShape* shape){
		shapes.push_back(shape);
	}
	
	
		void shapeModel::setOrigin(int orgx,int orgy,int orgz){
		origin.push_back(orgx);
		origin.push_back(orgy);
		origin.push_back(orgz);
	}
	
	void shapeModel::setImageParameters(int sx, int sy, int sz, float dx, float dy, float dz){
		xsize = sx;
		ysize = sy;
		zsize = sz;
		xdim = dx;
		ydim = dy;
		zdim = dz;
	}
	
	MShape* shapeModel::getShape(int ind){
		return shapes.at(ind);
	}
	int shapeModel::getShapeIndex(int label){
		int ind=-1;
		for (unsigned int i=0; i<labels.size();i++){
			if (label==labels.at(i)){
				ind=i;
			}
		}
		return ind;
	}
	
	Mesh shapeModel::getShapeMesh(int ind){
		return shapes.at(ind)->getMesh();
	}
	void shapeModel::setShapeMesh(int ind, Mesh m){
		return shapes.at(ind)->setMesh(m);
	}
	void shapeModel::centreAndSetShapeMesh(Mesh m,int ind){
		shapes.at(ind)->setMesh(getInverseTranslatedMesh(m));
	}
	
	vector<Vec> shapeModel::getShapeMode(int shape, int mode){
		return shapes.at(shape)->getModeVector(mode);
	}
	void shapeModel::setShapeTstatX(int shape, vector<float> v){
		shapes.at(shape)->setTstatX(v);
	}
	void shapeModel::setShapeTstatY(int shape, vector<float> v){
		shapes.at(shape)->setTstatY(v);
	}

	void shapeModel::setShapeTstatZ(int shape, vector<float> v){
		shapes.at(shape)->setTstatZ(v);
	}

	void shapeModel::setShapeMode(int shape, int mode, vector<Vec> v){
		shapes.at(shape)->setModeVector(v,mode);
	}
	
	bool shapeModel::getIntersection(){
		bool inter=false;
		for (int i=0;i<getNumberOfShapes();i++){
			if (shapes.at(i)->getMesh().real_self_intersection()){
				inter=true;
				break;
			}
			//cout<<"inteertest"<<i<<endl;
		}
		return inter;
	}
	Vec shapeModel::getModelGlobalTrans(int mode){
		float xmin=1,ymin=1,zmin=1;
		Vec tr(0,0,0);
		Vec tmp;
		for (int i=0;i<numShapes;i++){
			tmp=shapes.at(i)->getShapeGlobalTrans(mode);
			if (i==0){
				xmin=tmp.X;
				ymin=tmp.Y;
				zmin=tmp.Z;
			}else{ 
				if (abs(tmp.X)<abs(xmin)){
					xmin=tmp.X;
				}
				if (abs(tmp.Y)<abs(ymin)){
					ymin=tmp.Y;
				}
				if (abs(tmp.Z)<abs(zmin)){
					zmin=tmp.Z;
				}
				
			}
			tr.X=xmin;
			tr.Y=ymin;
			tr.Z=zmin;
		}
		return tr;
		
	}
	vector<float> shapeModel::projectShape (const int label, vector<float> vars,  shapeModel* newModel, int Nproj) {
		int shapeIndex=newModel->getShapeIndex(label);
		//cout<<"SHAPE INDEX: "<<shapeIndex<<endl;
		//find total translation
		vector<Vec> oldVec =  getShapeMode(getShapeIndex(label),0);//shapes.at(getShapeIndex(label))->getModeVector(0);
		for (unsigned int i=0; i<oldVec.size(); i++ ){
			oldVec.at(i)=oldVec.at(i)*(vars.at(0)*sqrt(eigenVals.at(0)));
		}
		//intialize frist mode
		//add other modes
		float normsq=0;
		for (unsigned int mode=1;mode<vars.size();mode++){//getNumberOfModes();mode++){
					vector<Vec> tmpvec = getShapeMode(getShapeIndex(label),mode);//shapes.at(shapeIndex)->getModeVector(mode);
		
			for (unsigned int i=0; i<oldVec.size(); i++ ){
				oldVec.at(i)+=(vars.at(mode)*sqrt(eigenVals.at(mode))*tmpvec.at(i));
				//				oldVec.at(i)+=tmpvec.at(i);
				normsq+=oldVec.at(i)|oldVec.at(i);
			}	
			
		}
		
		//calculate norm of shape deformat6ion vector, is unity if model has only one shape!
		
		vector<float> newvars;
		//will truncate modes to number in original model, not if using decopled translations
		for (int mode=0;mode<Nproj;mode++){//mode<newModel->getNumberOfModes();mode++){
			//cout<<"Mode "<<mode<<endl;
			vector<Vec> newModeVec = newModel->getShapeMode(shapeIndex, mode);
			float tmp=0;
			for (unsigned int i=0; i<oldVec.size(); i++ ){
				tmp+=oldVec.at(i)|newModeVec.at(i);
			}
			//assumes unit vector, i.e. single shape model
			//sqrt(newModel->getEigenValue(mode));
			//cout<<normsq<<endl;
			if (newModel->getEigenValue(mode)>0.01){
				tmp/=sqrt(newModel->getEigenValue(mode));//sqrt(normsq);
				if (abs(tmp)<1e-5){
					tmp=0;
					}
				}else{
				tmp=0;//prevent rounding errors
			}
			for (unsigned int i=0; i<oldVec.size(); i++ ){
			if (i==0){
				//cout<<"old: "<<oldVec.at(i).X<<" "<<oldVec.at(i).Y<<" "<<oldVec.at(i).Z<<" "<<sqrt(newModel->getEigenValue(mode))<<endl;
				//cout<<"new: "<<newModeVec.at(i).X<<" "<<newModeVec.at(i).Y<<" "<<newModeVec.at(i).Z<<endl;
				}
				
					oldVec.at(i)=oldVec.at(i)-tmp*sqrt(newModel->getEigenValue(mode))*newModeVec.at(i);
			
					
				
				//if (i==0){
				//	cout<<oldVec.at(i).X<<" "<<oldVec.at(i).Y<<" "<<oldVec.at(i).Z<<endl;
				//}
			}
			
			//cout<<"NEW VARS "<<tmp<<endl;
			newvars.push_back(tmp);
		}
		
		return newvars;
	}
	vector<float> shapeModel::projectVectors(const int label, vector<Vec> gvec, shapeModel* newModel, int beg, int Nproj) {
		vector<float> modematch;
		int shapeIndex=newModel->getShapeIndex(label);
		//cout<<"SHAPE INDEX: "<<shapeIndex<<endl;
		//find total translation
				
		//calculate norm of shape deformat6ion vector, is unity if model has only one shape!
				//will truncate modes to number in original model, not if using decopled translations
		for (int i=0;i<beg;i++){
			modematch.push_back(0);
		}
		for (int mode=beg;mode<beg+Nproj-1;mode++){//mode<newModel->getNumberOfModes();mode++){
			//cout<<"Mode "<<mode<<endl;
			vector<Vec> newModeVec = newModel->getShapeMode(shapeIndex, mode);
			float tmp=0;
		//	float normsq=0;
			for (unsigned int i=0; i<gvec.size(); i++ ){
				tmp+=(gvec.at(i)|newModeVec.at(i));
				
			}
			//assumes unit vector, i.e. single shape model
			//sqrt(newModel->getEigenValue(mode));
			//cout<<normsq<<endl;
			modematch.push_back(tmp);
			for (unsigned int i=0; i<gvec.size(); i++ ){
				gvec.at(i)=gvec.at(i)-tmp*newModeVec.at(i);
				
			}
			
			//cout<<"NEW VARS "<<tmp<<endl;
			//	newvars.push_back(tmp);
		}
		
		return modematch;
	}
void shapeModel::residual( volume<float>* image, const volume<short>* mask, Mesh* m, int label, int extent){
	//this doesn't explictly use the shape model, however is placed in here because a maks may be created using the shape model
	//this has been created for training purposes

	int bounds[6]={0,0,0,0,0,0};
	
	getBounds(m,bounds,0);
	
	ColumnVector B(4);
	
	float si=0.0,sx=0.0,sy=0.0,sz=0.0;	
	float ssi=0.0, ssx=0.0, ssy=0.0, ssz=0.0;
	float sxy=0.0, sxz=0.0, syz=0.0;
	float sxi=0.0,syi=0.0,szi=0.0;
	int numpix=0;
	
	//cout<<bounds[0]<<" "<<bounds[1]<<bounds[2]<<" "<<bounds[3]<<" "<<bounds[4]<<" "<<bounds[5]<<endl;
	for (int x=bounds[0];x<bounds[1];x++){
		for (int y=bounds[2];y<bounds[3];y++){
			for (int z=bounds[4];z<bounds[5];z++){
				
				if ((mask->value(x,y,z)==label)){
					float tempi;
					tempi=image->value(x,y,z);
					
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
	//	  cout<<" set up matrices "<<endl;
	Matrix A(4,4);											
	Matrix XY(4,1);
	
	A<<numpix<<sx<<sy<<sz
		<<sx<<ssx<<sxy<<sxz
		<<sy<<sxy<<ssy<<syz
		<<sz<<sxz<<syz<<ssz;
	XY<<si<<sxi<<syi<<szi;
	//look at slopes
		
	getBounds(m,bounds,6);
	
	B=A.i()*XY;
	//correct B(0) to be mean intensity
	B.element(0)=B.element(0)+B.element(1)*(bounds[0]+bounds[1])/2+B.element(2)*(bounds[2]+bounds[3])/2+B.element(3)*(bounds[4]+bounds[5])/2;	
	vbest.clear();//clear any existing
	for (int i=0;i<4;i++){
		vbest.push_back(B.element(i));
	}

	
	
	
	ColumnVector Y((bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]));
	ColumnVector Res((bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]));
	Matrix X((bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]),4);
	//cout<<"calculating residual at points"<<endl;
	int ptcount=0;
	for (int x=bounds[0];x<bounds[1];x++){
		for (int y=bounds[2];y<bounds[3];y++){
			for (int z=bounds[4];z<bounds[5];z++){
			//	cout<<"ptcount "<<x<<" "<<y<<" "<<z<<" "<<image.value(x,y,z)<<endl;
				X.element(ptcount,0)=1;
				X.element(ptcount,1)=x;
				X.element(ptcount,2)=y;
				X.element(ptcount,3)=z;
				Y.element(ptcount)=image->value(x,y,z);
				
				ptcount++;
			}
		}
	}
	
	
	Matrix val(1,1);
	//val=(A.i()*XY).t()*XY;
	//cout<<"rescal"<<endl;
	Res=Y-X*A.i()*XY;
		

	ptcount=0;
	for (int x=bounds[0];x<bounds[1];x++){
		for (int y=bounds[2];y<bounds[3];y++){
			for (int z=bounds[4];z<bounds[5];z++){
				image->value(x,y,z)=Res.element(ptcount);
	//			cout<<Res.element(ptcount)<<endl;
				ptcount++;
			}
		}
	}
	
	
	//cout<<"save"<<endl;
	//save_volume(*image,"restest");
	
	
	
	//cumVar += (ssi - val.element(0,0);	
	
	
}

void shapeModel::residualMeanOnly( volume<float>* image, const volume<short>* mask, Mesh* m, int label, int extent){
	//this doesn't explictly use the shape model, however is placed in here because a maks may be created using the shape model
	//this has been created for training purposes

	int bounds[6]={0,0,0,0,0,0};
	
	getBounds(m,bounds,4);
	
	ColumnVector B(4);
	
	float si=0.0;
	int numpix=0;
	
	//cout<<bounds[0]<<" "<<bounds[1]<<bounds[2]<<" "<<bounds[3]<<" "<<bounds[4]<<" "<<bounds[5]<<endl;
	for (int x=bounds[0];x<bounds[1];x++){
		for (int y=bounds[2];y<bounds[3];y++){
			for (int z=bounds[4];z<bounds[5];z++){
				
				if ((mask->value(x,y,z)==label)){
					float tempi;
					tempi=image->value(x,y,z);
					
					si+=tempi;
					
					numpix+=1;
					
				}
			}
		}
	}
	//	  cout<<" set up matrices "<<endl;

	getBounds(m,bounds,extent);
	

	//correct B(0) to be mean intensity
				
float mean=si/numpix;
cout<<"mean "<<mean<<endl;
	for (int x=bounds[0];x<bounds[1];x++){
		for (int y=bounds[2];y<bounds[3];y++){
			for (int z=bounds[4];z<bounds[5];z++){
				image->value(x,y,z)=image->value(x,y,z)-mean;
					//			cout<<Res.element(ptcount)<<endl;
				
			}
		}
	}
	
	
	//cout<<"save"<<endl;
	//save_volume(*image,"restest");
	
	
	
	//cumVar += (ssi - val.element(0,0);	
	
	
}
void shapeModel::residualMeanOnly( float mean, volume<float>* image, volume<float>* resimage, Mesh* m, int extent){
	//this doesn't explictly use the shape model, however is placed in here because a maks may be created using the shape model
	//this has been created for training purposes

	int bounds[6]={0,0,0,0,0,0};
	
	getBounds(m,bounds,extent);
	
	for (int x=bounds[0];x<bounds[1];x++){
		for (int y=bounds[2];y<bounds[3];y++){
			for (int z=bounds[4];z<bounds[5];z++){
				resimage->value(x,y,z)=image->value(x,y,z)-mean;
					//			cout<<Res.element(ptcount)<<endl;
				
			}
		}
	}
	
	
//	cout<<"save"<<endl;
//	save_volume(*resimage,"restest");
	
	
	
	//cumVar += (ssi - val.element(0,0);	
	
	
}
int shapeModel::cumResidual( volume<float>* image, const volume<short>* mask,  int label, int extent, float* totres){
	//this doesn't explictly use the shape model, however is placed in here because a maks may be created using the shape model
	//this has been created for training purposes
	
	
	float si=0.0,sx=0.0,sy=0.0,sz=0.0;	
	float ssi=0.0, ssx=0.0, ssy=0.0, ssz=0.0;
	float sxy=0.0, sxz=0.0, syz=0.0;
	float sxi=0.0,syi=0.0,szi=0.0;
	int numpix=0;
	
	for (int x=0;x<image->xsize();x++){
		for (int y=0;y<image->ysize();y++){
			for (int z=0;z<image->zsize();z++){
				if ((mask->value(x,y,z)==label)){
					float tempi;
					tempi=image->value(x,y,z);
					
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
	Matrix A(4,4);											
	Matrix XY(4,1);
	
	A<<numpix<<sx<<sy<<sz
		<<sx<<ssx<<sxy<<sxz
		<<sy<<sxy<<ssy<<syz
		<<sz<<sxz<<syz<<ssz;
	XY<<si<<sxi<<syi<<szi;
	//look at slopes
	
	
	
	ColumnVector Y(numpix);
	ColumnVector Res(numpix);
	Matrix X(numpix,4);
	//cout<<"calculating residual at points"<<endl;
	int ptcount=0;
//	for (int x=bounds[0];x<bounds[1];x++){
//		for (int y=bounds[2];y<bounds[3];y++){
//			for (int z=bounds[4];z<bounds[5];z++){
			for (int x=0;x<image->xsize();x++){
		for (int y=0;y<image->ysize();y++){
			for (int z=0;z<image->zsize();z++){
				if ((mask->value(x,y,z)==label)){
					//	cout<<"ptcount "<<x<<" "<<y<<" "<<z<<" "<<image.value(x,y,z)<<endl;
					X.element(ptcount,0)=1;
					X.element(ptcount,1)=x;
					X.element(ptcount,2)=y;
					X.element(ptcount,3)=z;
					Y.element(ptcount)=image->value(x,y,z);
					
					ptcount++;
				}
			}
		}
	}
	if (ptcount!=numpix){
		cerr<<"mismatch between in number of voxels"<<endl;
	}
	//cout<<"NUMPIX "<<numpix<<" ptcount "<<ptcount<<endl;
	Matrix val(1,1);
	val=Y.t()*Y-(A.i()*XY).t()*XY;
	
	//cout<<"Residual Calculation "<<val<<endl;
	//cout<<"rescal"<<endl;
	Res=Y-X*A.i()*XY;
	Matrix restot(1,1);
	restot=(Y-X*A.i()*XY).t()*(Y-X*A.i()*XY);
	float test=0.0;
		double sumres=0.0;
	for (int i=0; i<ptcount;i++){
		sumres+=Res.element(i)*Res.element(i);
		test=test+Res.element(i)*Res.element(i);
	//	cout<<Res.element(i)*Res.element(i)<<" "<<sumres<<" "<<restot.element(0,0)<<" "<<test<<endl;
		
	}
	
	*totres=restot.element(0,0);////val.element(0,0);
		//cout<<"res2w "<<sumres<<endl;
		return numpix;
}
void shapeModel::residual( ColumnVector Best, const volume<float>* image, volume<float>* Resimage, Mesh* m, int extent){
	//this doesn't explictly use the shape model, however is placed in here because a maks may be created using the shape model
	//this has been created for training purposes
	//cout<<"cost begin"<<endl;
	
	
	int bounds[6]={0,0,0,0,0,0};
	
	getBounds(m,bounds,extent);
	
	ColumnVector Y((bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]));
	ColumnVector Res((bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]));
	Matrix X((bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]),4);
//	cout<<"calculating residual at points"<<endl;
	int ptcount=0;
	for (int x=bounds[0];x<bounds[1];x++){
		for (int y=bounds[2];y<bounds[3];y++){
			for (int z=bounds[4];z<bounds[5];z++){
			//	cout<<"ptcount "<<x<<" "<<y<<" "<<z<<" "<<image.value(x,y,z)<<endl;
				X.element(ptcount,0)=1;
				X.element(ptcount,1)=x;
				X.element(ptcount,2)=y;
				X.element(ptcount,3)=z;
				Y.element(ptcount)=image->interpolate(x,y,z);
				
				ptcount++;
			}
		}
	}
	
	
	Matrix val(1,1);
	//val=(A.i()*XY).t()*XY;
	//cout<<"rescal"<<endl;
	Res=Y-X*Best;
	
	ptcount=0;
	for (int x=bounds[0];x<bounds[1];x++){
		for (int y=bounds[2];y<bounds[3];y++){
			for (int z=bounds[4];z<bounds[5];z++){
				Resimage->value(x,y,z)=Res.element(ptcount);
	//			cout<<Res.element(ptcount)<<endl;
				ptcount++;
			}
		}
	}
	
	
	//cout<<"save"<<endl;
	//save_volume(*image,"restest");
	
	
	
	//cumVar += (ssi - val.element(0,0));	
	
	
}

void shapeModel::residualMedianOnly(float median,volume<float>* image,volume<float>* resimage, Mesh* m, int extent){
	//this doesn't explictly use the shape model, however is placed in here because a maks may be created using the shape model
	//this has been created for training purposes
	
	int bounds[6]={0,0,0,0,0,0};
	
	getBounds(m,bounds,extent);
//cout<<bounds[0]<<" "<<bounds[1]<<" "<<bounds[2]<<" "<<bounds[3]<<" "<<bounds[4]<<" "<<bounds[5]<<endl;
//cout<<"median "<<median<<endl;
	for (int x=bounds[0];x<=bounds[1];x++){
		for (int y=bounds[2];y<=bounds[3];y++){
			for (int z=bounds[4];z<bounds[5];z++){
				resimage->value(x,y,z)=image->value(x,y,z)-median;
			}
		}
	}
	
	
	//cout<<"save"<<endl;
	//save_volume(*resimage,"medtest");
	
	
	
	//cumVar += (ssi - val.element(0,0);	
	
	
}
void shapeModel::getBounds(Mesh* m, int* bounds, int extra){
	
	float xmin=1000,xmax=-1000,ymin=1000,ymax=-1000,zmin=1000,zmax=-1000;
	for (vector<Mpoint*>::iterator i = m->_points.begin(); i!=m->_points.end(); i++ ){
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
	
	
	*bounds=static_cast<int>(floor(xmin/xdim)-1)-extra;
	if (*bounds<0){
		*bounds=0;
		}
	*(bounds+1)=static_cast<int>(ceil(xmax/xdim)+1+extra);
	*(bounds+2)=static_cast<int>(floor(ymin/ydim)-1-extra);
	if (*(bounds+2)<0){
		*(bounds+2)=0;
		}
	*(bounds+3)=static_cast<int>(ceil(ymax/ydim)+1+extra);
	*(bounds+4)=static_cast<int>(floor(zmin/zdim)-1-extra);
	if (*(bounds+4)<0){
		*(bounds+4)=0;
		}
	*(bounds+5)=static_cast<int>(ceil(zmax/zdim)+1+extra);

}


void shapeModel::draw_segment(volume<short>& image, const Pt& p1, const Pt& p2, int label)
{
	double xdim = (double) image.xdim();
	double ydim = (double) image.ydim();
	double zdim = (double) image.zdim();
	
	double d = (p1 - p2).norm();
	Vec n = (p1 - p2);
	n.normalize();
	double l = d*2;
	for (int i=0; i<=l; i++)
    {
		Pt p = p2 + i*.5 * n;
		image((int) floor((p.X)/xdim +.5),(int) floor((p.Y)/ydim +.5),(int) floor((p.Z)/zdim +.5)) = label;
		
	}
}


volume<short> shapeModel::draw_mesh(const volume<short>* image, const Mesh* m, int label)
{
	volume<short> res = *image;
	for (list<Triangle*>::const_iterator i = m->_triangles.begin(); i!=m->_triangles.end(); i++)
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

volume<short> shapeModel::make_mask_from_mesh(const volume<short> * image, Mesh* m, int label, bool binout)
{
//	cout<<"INVmake_mask_from_mesh begins"<<endl;
	//	cout<<"INVmake_mask_from_mesh begins"<<endl;
	
//	float xdim = (float) image->xdim();
//	float ydim = (float) image->ydim();
//	float zdim = (float) image->zdim();
	
	volume<short> mask=*image;
//	copyconvert(image,mask);
//	cout<<"copy"<<endl;
//	int xsize = mask.xsize();
//	int ysize = mask.ysize();
//	int zsize = mask.zsize();
	
	mask = 0;
//	cout<<"draw mesh"<<endl;
	if (binout){
	mask = draw_mesh(&mask, m,label+1);
	}else{
		mask = draw_mesh(&mask, m,label);

	}
//	cout<<"mesh drawn"<<endl;


	//save_volume(mask,"MASK");	

// THIS EXCLUDEDS THE ACTUAL MESH
	volume<short> otl=mask;//(mask-label)*-1;
	//otl=mask*2;
//	otl=0;	


//	save_volume(mask,"MASKIm");
//	float Rx=0,Ry=0,Rz=0;
//cout<<"set bounds "<<image.xsize()<<endl;
	getBounds(m,bounds,3);
	//cout<<bounds[0]<<" "<<bounds[1]<<" "<<bounds[2]<<" "<<bounds[3]<<" "<<bounds[4]<<" "<<bounds[5]<<endl;
	vector<Pt> current;
	current.clear();
	Pt c(bounds[0], bounds[2], bounds[4]);
	//cout<<"fill init "<<c.X<<" "<<c.Y<<" "<<c.Z<<" "<<mask.value(static_cast<int>(c.X),static_cast<int>(c.Y),static_cast<int>(c.Z))<<endl;
	
	mask.value(static_cast<int>(c.X),static_cast<int>(c.Y),static_cast<int>(c.Z)) = label;
//	cout<<"fill"<<endl;
			current.push_back(c);
		int fillCount=0;
	
		while (!current.empty())
		{
			//cout<<"looping"<<endl;
			Pt pc = current.back();
			int x, y, z;
			x=(int) pc.X; y=(int) pc.Y; z=(int) pc.Z;
			//cout<<x<<" "<<y<<" "<<z<<endl;
			current.pop_back();
			//if (mask.value(x, y, z)!=label){
			fillCount++;
			//}
		
	
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
		//save_volume(mask,"outmask");
		for (int i=bounds[0];i<bounds[1];i++){
			for (int j=bounds[2];j<bounds[3];j++){
				for (int k=bounds[4];k<bounds[5];k++){
					if (mask.value(i,j,k)==0){
						otl.value(i,j,k)=label;
					}
				}
			}
		}
	//cout<<"endmake mask"<<endl;
return otl;
}

int shapeModel::frac_overlap(const volume<short> segIm, const volume<short> gold,string fname)
{
	//calculate fractional overlap from distance maps
	ofstream foverlap;
	foverlap.open(fname.c_str());
	
	int numShapes=getNumberOfShapes();
	int label=0;
	
	//float unin, FN,TN,TP,FP;	
	for (int i=0;i<numShapes; i++){
		label=labels.at(i);
		
		int xsize= segIm.xsize();
		int ysize=segIm.ysize();
		int zsize=segIm.zsize();
		
		//find union and intersection
		float unin=0.0;
		float TP=0.0;
		float FN=0.0;
		float TN=0.0;
		float FP=0.0;
		for (int k=0;k<zsize;k++){
			for (int j=0;j<ysize;j++){ 
				for (int i= 0; i<xsize;i++){
					
				  if ((segIm(i,j,k)!=label)&&(gold(i,j,k)!=label)){
						TN++;
						}
					
					if ((segIm(i,j,k)==label)|(gold(i,j,k)==label)){
						unin++;
						if ((segIm(i,j,k)==label)&&(gold(i,j,k)==label)){
							TP++;
						}else if ((segIm(i,j,k)==label)&&(gold(i,j,k)!=label)){
							FP++;
						} else if ((segIm(i,j,k)!=label)&&(gold(i,j,k)==label)){
							FN++;
						}
					}
				}
			}
		}
		foverlap<<"label "<<label<<endl;
		foverlap<<"Fractional_Overlap "<<TP/unin<<endl;
		foverlap<<"Dice "<<2*TP/(2*TP+FP+FN)<<endl;
		foverlap<<"FNegative Fractional_Overlap "<<1-(FN)/unin<<endl;
		foverlap<<"FPositive Fractional_Overlap "<<1-(FP)/unin<<endl;
		foverlap<<"Sensitivity "<<TP/(TP+FN)<<endl;
		foverlap<<"Specificity "<<TN/(TN+FP)<<endl;
	}
	
	return 0;
}
int shapeModel::frac_overlap(Mesh m, const volume<short> gold, int label, string fname)
{
	//calculate fractional overlap from distance maps
	ofstream foverlap;
	foverlap.open(fname.c_str());
	volume<float> goldf;
	copyconvert(gold,goldf);
//	cout<<"copied "<<endl;
	volume<short> segIm=make_mask_from_mesh(&gold, &m, label,false);
//	cout<<"filled"<<endl;
	//int numShapes=getNumberOfShapes();
	//int label=0;
	
	//float unin, FN,TN,TP,FP;	
	//for (int i=0;i<numShapes; i++){
	//	label=labels.at(i);
		
		int xsize= segIm.xsize();
		int ysize=segIm.ysize();
		int zsize=segIm.zsize();
		
		//find union and intersection
		float unin=0.0;
		float TP=0.0;
		float FN=0.0;
		float TN=0.0;
		float FP=0.0;
		for (int k=0;k<zsize;k++){
			for (int j=0;j<ysize;j++){ 
				for (int i= 0; i<xsize;i++){
					
				  if ((segIm(i,j,k)!=label)&&(gold(i,j,k)!=label)){
						TN++;
						}
					
					if ((segIm(i,j,k)==label)|(gold(i,j,k)==label)){
						unin++;
						if ((segIm(i,j,k)==label)&&(gold(i,j,k)==label)){
							TP++;
						}else if ((segIm(i,j,k)==label)&&(gold(i,j,k)!=label)){
							FP++;
						} else if ((segIm(i,j,k)!=label)&&(gold(i,j,k)==label)){
							FN++;
						}
					}
				}
			}
		}
		foverlap<<"label "<<label<<endl;
		foverlap<<"Fractional_Overlap "<<TP/unin<<endl;
		foverlap<<"Dice "<<2*TP/(2*TP+FP+FN)<<endl;
		foverlap<<"FNegative Fractional_Overlap "<<1-(FN)/unin<<endl;
		foverlap<<"FPositive Fractional_Overlap "<<1-(FP)/unin<<endl;
		foverlap<<"Sensitivity "<<TP/(TP+FN)<<endl;
		foverlap<<"Specificity "<<TN/(TN+FP)<<endl;
	
	
	return 0;
}
int shapeModel::volumeDistance(const volume<short> segIm, const volume<short> gold,  string fname)
{
	//calculate fractional overlap from distance maps
	
		ofstream fs;
		fs.open(fname.c_str());
	
		//segIm and gold should have same dimensions
	int numShapes=getNumberOfShapes();
	int label=0;
//	cout<<"Xsize "<<xsize<<" "<<ysize<<" "<<zsize<<endl;
	int minxI=xsize,minyI=ysize,minzI=zsize;
	int maxxI=0,maxyI=0,maxzI=0;
	
	int minxG=xsize,minyG=ysize,minzG=zsize;
	int maxxG=0,maxyG=0,maxzG=0;
	
	int minxS=xsize,minyS=ysize,minzS=zsize;
	int maxxS=0,maxyS=0,maxzS=0;
	
	
	//float unin, FN,TN,TP,FP;	
//	cout<<"numSahpes"<<numShapes<<endl;
	
	
	
	for (int i=0;i<numShapes; i++){
	
		
		label=labels.at(i);
		//cout<<"label "<<label<<endl;
		//find union and intersection
		int TP=0;
		int FN=0;
		int FP=0;
		
		float TSvox=0;
		float TGvox=0;
		for (int k=0;k<zsize;k++){
			for (int j=0;j<ysize;j++){ 
				for (int i= 0; i<xsize;i++){
						if ((segIm(i,j,k)==label)&&(gold(i,j,k)==label)){
							TSvox++;
							TGvox++;
							//find bounds of intersection
							if (minxI>i){
								minxI=i;
							}
							if (maxxI<i){
								maxxI=i;
							}
							if (minyI>j){
								minyI=j;
							}
							if (maxyI<j){
								maxyI=j;
							}
							if (minzI>k){
								minzI=k;
							}
							if (maxzI<k){
								maxzI=k;
							}
					
							TP++;
						}else if ((segIm(i,j,k)==label)&&(gold(i,j,k)!=label)){
							TSvox++;
							//find bounds of seg
							if (minxS>i){
								minxS=i;
							}
							if (maxxS<i){
								maxxS=i;
							}
							if (minyS>j){
								minyS=j;
							}
							if (maxyS<j){
								maxyS=j;
							}
							if (minzS>k){
								minzS=k;
							}
							if (maxzS<k){
								maxzS=k;
							}
						
							FP++;
						} else if ((segIm(i,j,k)!=label)&&(gold(i,j,k)==label)){
						TGvox++;
							if (minxG>i){
								minxG=i;
							}
							if (maxxG<i){
								maxxG=i;
							}
							if (minyG>j){
								minyG=j;
							}
							if (maxyG<j){
								maxyG=j;
							}
							if (minzG>k){
								minzG=k;
							}
							if (maxzG<k){
								maxzG=k;
							}
							
							FN++;
						}
					
				}
			}
		}
		//cout<<"new bounds founds...calculate distance"<<endl;
		//calculate FNdistance
		//volume<float> FNdmap;
		//copyconvert(segIm, FNdmap);
		//FNdmap=0;
		float FNd=0;
		float FNdSq=0;
		int FNvox=0;
		for (int k=minzG;k<=maxzG;k++){
			for (int j=minyG;j<=maxyG;j++){ 
				for (int i= minxG; i<=maxxG;i++){
					float dist=10000;
				//	cout<<"dist"<<endl;
					if ((segIm(i,j,k)!=label)&&(gold(i,j,k)==label)){
						//calculate distance to nearest intersection point
						for (int n=minzI;n<=maxzI;n++){
							for (int m=minyI;m<=maxyI;m++){ 
								for (int l= minxI; l<=maxxI;l++){
							//		cout<<segIm(l,m,n)<<" "<<gold(l,m,n)<<endl;
									if ((segIm(l,m,n)==label)&&(gold(l,m,n)==label)){
											float disttemp=(((l-i)*xdim*(l-i)*xdim) + ((m-j)*ydim*(m-j)*ydim) +((n-k)*zdim*(n-k)*zdim));
											if (disttemp<dist){
												dist=disttemp;
				//								FNdmap(i,j,k)=sqrt(dist);
											}
									}	
								
								}
							}
						}
					//	cout<<"dist "<<dist<<endl;
						
							FNdSq+=(dist);
							FNd+=sqrt(dist);
							FNvox++;
						
					}
					
				}
			}
		}
		//cout<<"saving FPdmap"<<endl;
		//save_volume(FNdmap,"FNdmap");	
		//volume<float> FPdmap;
		//copyconvert(segIm, FPdmap);
		//FPdmap=0;
		float FPd=0;
		float FPdSq=0;
		int FPvox=0;
		for (int k=minzS;k<=maxzS;k++){
			for (int j=minyS;j<=maxyS;j++){ 
				for (int i= minxS; i<=maxxS;i++){
					float dist=10000;
				//	cout<<"dist"<<endl;
					if ((segIm(i,j,k)==label)&&(gold(i,j,k)!=label)){
						//calculate distance to nearest intersection point
						for (int n=minzI;n<=maxzI;n++){
							for (int m=minyI;m<=maxyI;m++){ 
								for (int l= minxI; l<=maxxI;l++){
							//		cout<<segIm(l,m,n)<<" "<<gold(l,m,n)<<endl;
									if ((segIm(l,m,n)==label)&&(gold(l,m,n)==label)){
											float disttemp=(((l-i)*xdim*(l-i)*xdim) + ((m-j)*ydim*(m-j)*ydim) +((n-k)*zdim*(n-k)*zdim));
											if (disttemp<dist){
												dist=disttemp;
										//		FPdmap(i,j,k)=sqrt(dist);
												
											}
									}	
								
								}
							}
						}
						//	cout<<"dist "<<dist<<endl;
						
						FPdSq+=(dist);
						FPd+=sqrt(dist);
						FPvox++;
						
						
						
					}
					
				}
			}
		}
		//cout<<"saving FPdmap"<<endl;
		//save_volume(FPdmap,"FPdmap");	
			
		
		fs<<label<<" "<<TGvox<<" "<<TSvox<<endl;
		fs<<FPd<<" "<<FPdSq<<" "<<FPvox<<endl;
		fs<<FNd<<" "<<FNdSq<<" "<<FNvox<<endl;
		
		

	}
	fs.close();
	//foverlap.close();
	return 0;
}
float shapeModel::volumeDistance(const volume<short>* segIm, const volume<float>* gold)
{
	//calculate fractional overlap from distance maps
		
		//segIm and gold should have same dimensions
	//int numShapes=getNumberOfShapes();
	//assumes 1 shape
	int minxI=xsize,minyI=ysize,minzI=zsize;
	int maxxI=0,maxyI=0,maxzI=0;
	
	int minxG=xsize,minyG=ysize,minzG=zsize;
	int maxxG=0,maxyG=0,maxzG=0;
	
	int minxS=xsize,minyS=ysize,minzS=zsize;
	int maxxS=0,maxyS=0,maxzS=0;
	

	
	//for (int i=0;i<numShapes; i++){
		//find union and intersection
		int TP=0;
		int FN=0;
		int FP=0;
		
		float TSvox=0;
		float TGvox=0;
		for (int k=0;k<zsize;k++){
			for (int j=0;j<ysize;j++){ 
				for (int i= 0; i<xsize;i++){
						if ((segIm->value(i,j,k)>0)&&(gold->value(i,j,k)>0)){
							TSvox++;
							TGvox++;
							//find bounds of intersection
							if (minxI>i){
								minxI=i;
							}
							if (maxxI<i){
								maxxI=i;
							}
							if (minyI>j){
								minyI=j;
							}
							if (maxyI<j){
								maxyI=j;
							}
							if (minzI>k){
								minzI=k;
							}
							if (maxzI<k){
								maxzI=k;
							}
					
							TP++;
						}else if ((segIm->value(i,j,k)>0)&&(gold->value(i,j,k)==0)){
							TSvox++;
							//find bounds of seg
							if (minxS>i){
								minxS=i;
							}
							if (maxxS<i){
								maxxS=i;
							}
							if (minyS>j){
								minyS=j;
							}
							if (maxyS<j){
								maxyS=j;
							}
							if (minzS>k){
								minzS=k;
							}
							if (maxzS<k){
								maxzS=k;
							}
						
							FP++;
						} else if ((segIm->value(i,j,k)==0)&&(gold->value(i,j,k)>0)){
						TGvox++;
							if (minxG>i){
								minxG=i;
							}
							if (maxxG<i){
								maxxG=i;
							}
							if (minyG>j){
								minyG=j;
							}
							if (maxyG<j){
								maxyG=j;
							}
							if (minzG>k){
								minzG=k;
							}
							if (maxzG<k){
								maxzG=k;
							}
							
							FN++;
						}
					
				}
			}
		}
		//cout<<"new bounds founds...calculate distance"<<endl;
		//calculate FNdistance
		
		float FNd=0;
				for (int k=minzG;k<=maxzG;k++){
			for (int j=minyG;j<=maxyG;j++){ 
				for (int i= minxG; i<=maxxG;i++){
					float dist=10000;
					if ((segIm->value(i,j,k)==0)&&(gold->value(i,j,k)>0)){
						//calculate distance to nearest intersection point
						for (int n=minzI;n<=maxzI;n++){
							for (int m=minyI;m<=maxyI;m++){ 
								for (int l= minxI; l<=maxxI;l++){
									if ((segIm->value(l,m,n)>0)&&(gold->value(l,m,n)>0)){
											float disttemp=(((l-i)*xdim*(l-i)*xdim) + ((m-j)*ydim*(m-j)*ydim) +((n-k)*zdim*(n-k)*zdim));
											if (disttemp<dist){
												dist=disttemp;
											}
									}	
								
								}
							}
						}
						
							FNd+=sqrt(dist);
						
					}
					
				}
			}
		}
		
		
		float FPd=0;
	
		for (int k=minzS;k<=maxzS;k++){
			for (int j=minyS;j<=maxyS;j++){ 
				for (int i= minxS; i<=maxxS;i++){
					float dist=10000;
					if ((segIm->value(i,j,k)>0)&&(gold->value(i,j,k)==0)){
						//calculate distance to nearest intersection point
						for (int n=minzI;n<=maxzI;n++){
							for (int m=minyI;m<=maxyI;m++){ 
								for (int l= minxI; l<=maxxI;l++){
									if ((segIm->value(l,m,n)>0)&&(gold->value(l,m,n)>0)){
											float disttemp=(((l-i)*xdim*(l-i)*xdim) + ((m-j)*ydim*(m-j)*ydim) +((n-k)*zdim*(n-k)*zdim));
											if (disttemp<dist){
												dist=disttemp;
												
											}
									}	
								
								}
							}
						}
						
						FPd+=sqrt(dist);
	
					}
					
				}
			}
		}
			return FPd+FNd;
}

float shapeModel::volumeDistance(const volume<short>* segIm, const volume<float>* gold, int *gbounds, Mesh* m)
{
	int sbounds[6]={0,0,0,0,0,0};
	getBounds(m,sbounds,2);	
	
	int minxI=xsize,minyI=ysize,minzI=zsize;
	int maxxI=0,maxyI=0,maxzI=0;
	
	
	int minxS=sbounds[0],minyS=sbounds[2],minzS=sbounds[4];
	int maxxS=sbounds[1],maxyS=sbounds[3],maxzS=sbounds[5];
	
	int minxG=gbounds[0],minyG=gbounds[2],minzG=gbounds[4];
	int maxxG=gbounds[1],maxyG=gbounds[3],maxzG=gbounds[5];
	 int minxsg,minysg,minzsg; 
		int maxxsg,maxysg,maxzsg; 
		
		//set mins and maxes to search for intersection minmax
		if (minxG>minxS){
			minxsg=minxS;
		}else{
			minxsg=minxG;
		}
		if (minyG>minyS){
			minysg=minyS;
		}else{
			minysg=minyG;
		}
		if (minzG>minzS){
			minzsg=minzS;
		}else{
			minzsg=minzG;
		}

		if (maxxG<maxxS){
			maxxsg=maxxS;
		}else{
			maxxsg=maxxG;
		}
		if (maxyG<maxyS){
			maxysg=maxyS;
		}else{
			maxysg=maxyG;
		}
		if (maxzG<maxzS){
			maxzsg=maxzS;
		}else{
			maxzsg=maxzG;
		}

				
		for (int k=minzsg;k<maxzsg;k++){
			for (int j=minysg;j<maxysg;j++){ 
				for (int i= minxsg; i<maxxsg;i++){
						if ((segIm->value(i,j,k)>0)&&(gold->value(i,j,k)>0)){
							
							//find bounds of intersection
							if (minxI>i){
								minxI=i;
							}
							if (maxxI<i){
								maxxI=i;
							}
							if (minyI>j){
								minyI=j;
							}
							if (maxyI<j){
								maxyI=j;
							}
							if (minzI>k){
								minzI=k;
							}
							if (maxzI<k){
								maxzI=k;
							}
					
							
						}
						}
						}
						}		
		//cout<<"new bounds founds...calculate distance"<<endl;
		//calculate FNdistance
		
		float FNd=0;
		for (int k=minzG;k<=maxzG;k++){
			for (int j=minyG;j<=maxyG;j++){ 
				for (int i= minxG; i<=maxxG;i++){
					float dist=10000;
					if ((segIm->value(i,j,k)==0)&&(gold->value(i,j,k)>0)){
						//calculate distance to nearest intersection point
						for (int n=minzI;n<=maxzI;n++){
							for (int m=minyI;m<=maxyI;m++){ 
								for (int l= minxI; l<=maxxI;l++){
									if ((segIm->value(l,m,n)>0)&&(gold->value(l,m,n)>0)){
											float disttemp=(((l-i)*xdim*(l-i)*xdim) + ((m-j)*ydim*(m-j)*ydim) +((n-k)*zdim*(n-k)*zdim));
											if (disttemp<dist){
												dist=disttemp;
											}
									}	
								
								}
							}
						}
						
							FNd+=sqrt(dist);
						
					}
					
				}
			}
		}
		
		
		float FPd=0;
	
		for (int k=minzS;k<=maxzS;k++){
			for (int j=minyS;j<=maxyS;j++){ 
				for (int i= minxS; i<=maxxS;i++){
					float dist=10000;
					if ((segIm->value(i,j,k)>0)&&(gold->value(i,j,k)==0)){
						//calculate distance to nearest intersection point
						for (int n=minzI;n<=maxzI;n++){
							for (int m=minyI;m<=maxyI;m++){ 
								for (int l= minxI; l<=maxxI;l++){
									if ((segIm->value(l,m,n)>0)&&(gold->value(l,m,n)>0)){
											float disttemp=(((l-i)*xdim*(l-i)*xdim) + ((m-j)*ydim*(m-j)*ydim) +((n-k)*zdim*(n-k)*zdim));
											if (disttemp<dist){
												dist=disttemp;
												
											}
									}	
								
								}
							}
						}
						
						FPd+=sqrt(dist);
	
					}
					
				}
			}
		}
			return FPd+FNd;
}


void shapeModel::volumeBounds( const volume<float>* gold, int *gbounds)
{
	
	for (int k=0;k<zsize;k++){
		for (int j=0;j<ysize;j++){ 
			for (int i=0; i<xsize;i++){
				if (gold->value(i,j,k)>0){
					
					//find bounds of intersection
					if (gbounds[0]>i){
						gbounds[0]=i;
					}
					if (gbounds[1]<i){
						gbounds[1]=i;
					}
					if (gbounds[2]>j){
						gbounds[2]=j;
					}
					if (gbounds[3]<j){
						gbounds[3]=j;
					}
					if (gbounds[4]>k){
						gbounds[4]=k;
					}
					if (gbounds[5]<k){
						gbounds[5]=k;
					}
					
					
				}
			}
		}
	}		
	//cout<<"new bounds founds...calculate distance"<<endl;
}



int shapeModel::meshDistance(const volume<short>* gold, int shape, int tol, vector<float>* vdists)
{
	//returns a sorted min distance for each vertex 
	vdists->clear();
	float dist=10000;
	int bounds[6]={0,0,0,0,0,0};
	Mesh m=getTranslatedMesh(shape);
	//Mesh m=shapes.at(shape)->getMesh();
	int label=getLabel(shape);
	getBounds(&m,bounds,tol);
	vector <float>::iterator Iter;
	int lastminz=bounds[4],lastminy=bounds[2],lastminx=bounds[0];
	//cout<<bounds[4]<<" "<<bounds[2]<<" "<<bounds[0]<<endl;
	float ldist=0;
	for (vector<Mpoint*>::iterator i = m._points.begin(); i!=m._points.end(); i++ ){
		float tempx=(*i)->get_coord().X;
		float tempy=(*i)->get_coord().Y;
		float tempz=(*i)->get_coord().Z;
		
		
		
		int tminz=lastminz,tminy=lastminy,tminx=lastminx;
			int tmaxz=bounds[5],tmaxy=bounds[3],tmaxx=bounds[1];
		
		//do in form of while loops and update end limit
		if (dist!=10000){
			//calculates distance to lastpoints
			ldist=sqrt((tempx-lastminx*xdim)*(tempx-lastminx*xdim) + (tempy-lastminy*ydim)*(tempy-lastminy*ydim) + (tempz-lastminz*zdim)*(tempz-lastminz*zdim));
			tminz=static_cast<int>(floor(tempz-ldist)/zdim);
			tmaxz=static_cast<int>(ceil(tempz+ldist)/zdim);
			tmaxy=static_cast<int>(ceil(tempy+ldist)/ydim);
			tminy=static_cast<int>(floor(tempy-ldist)/ydim);
			tmaxx=static_cast<int>(ceil(tempx+ldist)/xdim);
			tminx=static_cast<int>(floor(tempx-ldist)/xdim);
		}else{
		 tminz=bounds[4];
		 tminy=bounds[2];
		 tminx=bounds[0];
			tmaxz=bounds[5];
			tmaxy=bounds[3];
			tmaxx=bounds[1];
		
		}
		dist=10000;
		//cout<<tmaxz<<" "<<tmaxy<<" "<<tmaxx<<endl;
		int n=tminz;
		while (n<=tmaxz){
			//for (int n=bounds[4];n<=bounds[5];n++){
			int m=tminy; 
			while (m<=tmaxy){
			//		for (int m=bounds[2];m<=bounds[3];m++){ 
			int l=tminx;
				while (l<=tmaxx){
				//	for (int l=bounds[0]; l<=bounds[1];l++){
				//cout<<n<<" "<<m<<" "<<l<<endl;
					if (gold->value(l,m,n)==label){
						float disttemp= sqrt((tempx-l*xdim)*(tempx-l*xdim) + (tempy-m*ydim)*(tempy-m*ydim) + (tempz-n*zdim)*(tempz-n*zdim));
						if (disttemp<dist){
							dist=(disttemp);
							lastminz=n;
							lastminy=m;
							lastminx=l;
						}
					}	
					l++;
					}
				m++;
				}
			n++;
				}
		//		cout<<"add to vector"<<endl;
		//dist=sqrt(dist);
		if (vdists->empty()){
		//	cout<<"empty"<<endl;
			vdists->push_back(dist);
		}else if (dist>=vdists->back()){
	//	cout<<"back"<<endl;
			vdists->push_back(dist);
		}else {
		//	cout<<"insert"<<endl;
			for ( Iter = vdists->begin( ) ; Iter != vdists->end( ) ; Iter++ ){
				if (dist<*Iter){
					
					vdists->insert(Iter,dist);
					break;
				}
			
			}
			
			
		}
	//	cout<<"vec updated "<<endl;
		
			}
	//			for ( Iter = vdists->begin( ) ; Iter != vdists->end( ) ; Iter++ ){
	//		cout<<*Iter<<" ";
	//			Iter++;
	//		}
	//		cout<<endl;

	
	return 0;
			}

int shapeModel::intensityHist(const volume<float>* image, const volume<short>* mask, Mesh* m,int label , vector<float>* vgraylevels)
{
	//returns a sorted min distance for each vertex 
	vgraylevels->clear();
	float dist=10000;
//	int bounds[6]={0,0,0,0,0,0};
	//Mesh m=getTranslatedMesh(shape);
	//Mesh m=shapes.at(shape)->getMesh();
	//int label=getLabel(shape);
//getBounds(m,bounds,tol);
	int bounds[6]={0,0,0,0,0,0};
	
	getBounds(m,bounds,5);
	
	vector <float>::iterator Iter;
	//cout<<bounds[4]<<" "<<bounds[2]<<" "<<bounds[0]<<endl;
	for (int n=bounds[4];n<=bounds[5];n++){
		for (int m=bounds[2];m<=bounds[3];m++){ 
			for (int l=bounds[0]; l<=bounds[1];l++){
				if (mask->value(l,m,n)==label){
					dist=image->value(l,m,n);
					if (vgraylevels->empty()){
					//	cout<<"empty"<<endl;
						vgraylevels->push_back(dist);
					}else if (dist>=vgraylevels->back()){
					//	cout<<"back"<<endl;
						vgraylevels->push_back(dist);
					}else {
					//	cout<<"insert"<<endl;
						for ( Iter = vgraylevels->begin( ) ; Iter !=vgraylevels->end( ) ; Iter++ ){
							if (dist<*Iter){
								
								vgraylevels->insert(Iter,dist);
								break;
							}
							
						}
						
						
					}
					
					
				
				
				}
			}
		}
	}
		
	
	return 0;
			}
int shapeModel::intensityHistMaxMin(const volume<float>* image, const volume<short>* mask, Mesh* m,int label , vector<float>* vgraylevels, float* maxint, float* minint)
{
	//this does not sort the intesnities but save max min intensity
	//returns a sorted min distance for each vertex 
	*maxint=-1000000;
	*minint=1000000;
	vgraylevels->clear();
	float dist=10000;
//	int bounds[6]={0,0,0,0,0,0};
	//Mesh m=getTranslatedMesh(shape);
	//Mesh m=shapes.at(shape)->getMesh();
	//int label=getLabel(shape);
//getBounds(m,bounds,tol);
	int bounds[6]={0,0,0,0,0,0};
	
	getBounds(m,bounds,5);
	
	vector <float>::iterator Iter;
	//cout<<bounds[4]<<" "<<bounds[2]<<" "<<bounds[0]<<endl;
	for (int n=bounds[4];n<=bounds[5];n++){
		for (int m=bounds[2];m<=bounds[3];m++){ 
			for (int l=bounds[0]; l<=bounds[1];l++){
				//cout<<"in bounds?"<<endl;
				if (mask->in_bounds(l,m,n)){
					if (mask->value(l,m,n)==label){
						dist=image->value(l,m,n);
						vgraylevels->push_back(dist);
						if ((*maxint)<dist){ *maxint=dist; }
						if ((*minint)>dist){ *minint=dist; }
						
						
					}
				}
			}
		}
	}
		
	//cout<<"done intesn hist"<<endl;
	return 0;
			}
int shapeModel::intensityHistMult(const volume<float>* image, const volume<short>* mask, Mesh* m,vector<int> vlabel , vector<float>* vgraylevels)
{
	//thsi is used when training from multiple labels
	//returns a sorted min distance for each vertex 
	vgraylevels->clear();
	float dist=10000;
//	int bounds[6]={0,0,0,0,0,0};
	//Mesh m=getTranslatedMesh(shape);
	//Mesh m=shapes.at(shape)->getMesh();
	//int label=getLabel(shape);
//getBounds(m,bounds,tol);
	int bounds[6]={0,0,0,0,0,0};
	
	getBounds(m,bounds,5);
	
	vector <float>::iterator Iter;
	//cout<<bounds[4]<<" "<<bounds[2]<<" "<<bounds[0]<<endl;
	for (int n=bounds[4];n<=bounds[5];n++){
		for (int m=bounds[2];m<=bounds[3];m++){ 
			for (int l=bounds[0]; l<=bounds[1];l++){
				int label=mask->value(l,m,n);
				bool foundLabel=false;
				for (unsigned int i=0; i<vlabel.size();i++){
					if (label==vlabel.at(i)){
						foundLabel=true;
					}
				}
				if (foundLabel){
					dist=image->value(l,m,n);
					if (vgraylevels->empty()){
					//	cout<<"empty"<<endl;
						vgraylevels->push_back(dist);
					}else if (dist>=vgraylevels->back()){
					//	cout<<"back"<<endl;
						vgraylevels->push_back(dist);
					}else {
					//	cout<<"insert"<<endl;
						for ( Iter = vgraylevels->begin( ) ; Iter !=vgraylevels->end( ) ; Iter++ ){
							if (dist<*Iter){
								
								vgraylevels->insert(Iter,dist);
								break;
							}
							
						}
						
						
					}
					
					
				
				
				}
			}
		}
	}
		
	
	return 0;
			}
			
float shapeModel::EMgmm3(vector<float>* vgl,bool lesser,int niter, vector<float>* mu, vector<float>* vvars, float init0, float init1, float init2){
//implement a mixture of 3 gaussians .... csf,grey matter, white matter
//returns lower mean
	//implement a mixture of 2 gaussians
	vector <float>::iterator Iter;
	float si=0.0, ssi=0.0;
	int N=0;
//		float muest1=10;//vgl->at(static_cast<int>(floor(static_cast<float>(vgl->size())/3.0)));
//	float muest2=vgl->at(static_cast<int>(floor(vgl->size()/2.0)));
//	float muest3=vgl->at(static_cast<int>(floor(vgl->size()*2.0/3.0)));;
	
		float muest1=init0;//vgl->at(static_cast<int>(floor(static_cast<float>(vgl->size())/3.0)));
	float muest2=init1;//vgl->at(static_cast<int>(floor(vgl->size()/2.0)));
	float muest3=init2;//vgl->at(static_cast<int>(floor(vgl->size()*2.0/3.0)));;
	for ( Iter = vgl->begin( ) ; Iter !=vgl->end( ) ; Iter++ ){
	//	cout<<"Iter "<<*Iter<<endl;
		si+=*Iter;
		ssi+=(*Iter)*(*Iter);
		N++;
	}
	//cout<<"SSI "<<si<<" "<<ssi<<endl;
	float varest1=1;//(ssi-si*si/N)/(N-1);
	float varest2=1;//varest1;
	float varest3=1;//varest1;
	//cout<<"init estimates "<<muest1<<" "<<muest2<<" "<<varest1<<" "<<varest2<<endl;
	//int niter=10;
	
	//initial misture weigthing
	float f1=0.5;//1/3.0;
	float f2=0.25;//1/3.0;
	float f3=1-f1-f2;
	//cout<<"Niter "<<niter<<endl;
	//float mulast=0;
	 
	for (int i=0; i<niter; i++){
		float C1=1/sqrt(2.0*3.14*varest1);
		float C2=1/sqrt(2.0*3.14*varest2);
		float C3=1/sqrt(2.0*3.14*varest3);
		//vector<float> resp1,resp2;
	//	vector<float>resp;
		float sumfresp1=0.0,sumfresp2=0.0,sumfresp3=0.0;
		//update means and mixture
		float muesttemp1=0.0, muesttemp2=0.0, muesttemp3=0.0;
			float varesttemp1=0.0, varesttemp2=0.0, varesttemp3=0.0;
			float fresp1=0,fresp2=0,fresp3=0;
		for ( Iter = vgl->begin( ) ; Iter !=vgl->end( ) ; Iter++ ){
			float resp1=C1*exp(-0.5*(*Iter-muest1)*(*Iter-muest1)/varest1);
			float resp2=C2*exp(-0.5*(*Iter-muest2)*(*Iter-muest2)/varest2);
			float resp3=C3*exp(-0.5*(*Iter-muest3)*(*Iter-muest3)/varest3);
			float Sresp=f1*resp1+f2*resp2+f3*resp3;
			//cout<<"resp "<<resp1<<" "<<resp2<<" "<<resp3<<" "<<Sresp<<endl;

			fresp1=f1*resp1/Sresp;
			fresp2=f2*resp2/Sresp;
			fresp3=f3*resp3/Sresp;
			
			sumfresp1+=fresp1;
			sumfresp2+=fresp2;
			sumfresp3+=fresp3;
	
			muesttemp1+= fresp1*(*Iter);
			muesttemp2+= fresp2*(*Iter);
			muesttemp3+= fresp3*(*Iter);
		//				cout<<"fresp "<<fresp1<<" "<<fresp2<<" "<<fresp3<<" "<<muesttemp1<<" "<<muesttemp2<<" "<<muesttemp3<<endl;

			
			//this should really be in a seperate loop
				varesttemp1+=fresp1*(*Iter-muest1)*(*Iter-muest1);
			varesttemp2+=fresp2*(*Iter-muest2)*(*Iter-muest2);
			varesttemp3+=fresp3*(*Iter-muest3)*(*Iter-muest3);
			
		}
		muest1=muesttemp1/sumfresp1;
		muest2=muesttemp2/sumfresp2;
		muest3=muesttemp3/sumfresp3;
//cout<<"estimeate "<<muest1<<" "<<sumfresp1<<" "<<muest2<<" "<<sumfresp2<<" "<<endl;
	//	for ( Iter = vgl->begin( ) ; Iter !=vgl->end( ) ; Iter++ ){
		
	
	//		}

			
					varest1=varesttemp1/sumfresp1;
		varest2=varesttemp2/sumfresp2;
		varest3=varesttemp3/sumfresp3;
		f1=sumfresp1/N;
		f2=sumfresp2/N;
		f3=sumfresp3/N;
		//cout<<"estimeate "<<muest1<<" "<<muest2<<" "<<muest3<<" "<<varest1<<" "<<varest2<<" "<<varest3<<" "<<f1<<" "<<f2<<" "<<f3<<endl;
		/*
		//implement a convergence criteria
		if (lesser){
			float dmu=abs(mulast-muest1);
			if (dmu<0.05){
			cout<<"estimeate "<<muest1<<" "<<muest2<<" "<<varest1<<" "<<varest2<<" "<<f1<<" "<<i<<endl;
				break;
			}else{
				mulast=muest1;
			}
		}else{
			float dmu=abs(mulast-muest2);
			if (dmu<0.05){
				break;
			}else{
				mulast=muest2;
			}
			
		}
	}
	
	if (muest1<=muest2){
		if (lesser){
			//cout<<"setmuest 1 to 40"<<endl;
			//muest1=40;
			*var=varest1;
			return muest1;
		}else{
			return muest2;
		}
	}else{
		if (lesser){
			return muest2;
		}else{
			return muest2;
		}
	}*/
}

mu->clear();
vvars->clear();
if (muest1<muest2){	
	if (muest1<muest3){
		mu->push_back(muest1);
		vvars->push_back(varest1);
//		vf->push_back(f1);
		if(muest3<muest2){
			mu->push_back(muest3);
//			vf->push_back(f1);
			vvars->push_back(varest3);
			
			mu->push_back(muest2);
//			vf->push_back(f1);
			vvars->push_back(varest2);
			
		}else{
			mu->push_back(muest2);
//			vf->push_back(f1);
			vvars->push_back(varest2);
			
			mu->push_back(muest3);
//			vf->push_back(f1);
			vvars->push_back(varest3);
			
		}
	}
	
}else if (muest2<muest3){
	mu->push_back(muest2);
//	vf->push_back(f1);
	vvars->push_back(varest2);
	
	if(muest3<muest1){
		mu->push_back(muest3);
//		vf->push_back(f1);
		vvars->push_back(varest3);
		
		mu->push_back(muest1);
//		vf->push_back(f1);
		vvars->push_back(varest1);
		
	}else{
		mu->push_back(muest2);
//		vf->push_back(f1);
		vvars->push_back(varest2);
		
		mu->push_back(muest1);
//		vf->push_back(f1);
		vvars->push_back(varest1);
		
	}
}else{
	mu->push_back(muest3);
//	vf->push_back(f1);
	vvars->push_back(varest3);
	
	mu->push_back(muest2);
//	vf->push_back(f1);
	vvars->push_back(varest2);
	
	mu->push_back(muest1);
//	vf->push_back(f1);
	vvars->push_back(varest1);
	
}
		//cout<<"estimeate "<<mu->at(0)<<" "<<mu->at(1)<<" "<<mu->at(2)<<" "<<vvars->at(0)<<" "<<vvars->at(1)<<" "<<vvars->at(2)<<" "<<f1<<" "<<f2<<" "<<f3<<endl;

return -1; 	
		}


float shapeModel::EMgmm(vector<float>* vgl,bool lesser,int niter){
//returns lower mean
	//implement a mixture of 2 gaussians
	vector <float>::iterator Iter;
	float si=0.0, ssi=0.0;
	int N=0;
		float muest1=vgl->at(static_cast<int>(floor(static_cast<float>(vgl->size())/4.0)));
	float muest2=vgl->at(static_cast<int>(floor(vgl->size()*3.0/4.0)));
	
	for ( Iter = vgl->begin( ) ; Iter !=vgl->end( ) ; Iter++ ){
	//	cout<<"Iter "<<*Iter<<endl;
		si+=*Iter;
		ssi+=(*Iter)*(*Iter);
		N++;
	}
	
	float varest1=(ssi-si*si/N)/(N-1);
	float varest2=varest1;
	//cout<<"init estimates "<<muest1<<" "<<muest2<<" "<<varest1<<" "<<varest2<<endl;
	//int niter=10;
	
	//initial misture weigthing
	float f=0.5;
	//cout<<"Niter "<<niter<<endl;
	float mulast=0;
	 
	for (int i=0; i<niter; i++){
		float C1=1/sqrt(2.0*3.14*varest1);
		float C2=1/sqrt(2.0*3.14*varest2);
		//vector<float> resp1,resp2;
	//	vector<float>resp;
		float sumresp=0.0,invsumresp=0.0;
		//update means and mixture
		float muesttemp1=0.0, muesttemp2=0.0;
			float varesttemp1=0.0, varesttemp2=0.0;
		for ( Iter = vgl->begin( ) ; Iter !=vgl->end( ) ; Iter++ ){
			float resp1=C1*exp(-0.5*(*Iter-muest1)*(*Iter-muest1)/varest1);
			float resp2=C2*exp(-0.5*(*Iter-muest2)*(*Iter-muest2)/varest2);
			float resp=f*resp2/( (1-f)*resp1 + f*resp2 );
			sumresp+=resp;
			invsumresp+=1-resp;
			muesttemp1+= (1-resp)*(*Iter);
			muesttemp2+= resp*(*Iter);
			varesttemp1+=(1-resp)*(*Iter-muest1)*(*Iter-muest1);
			varesttemp2+=resp*(*Iter-muest2)*(*Iter-muest2);
			
		}
		muest1=muesttemp1/invsumresp;
		muest2=muesttemp2/sumresp;
		varest1=varesttemp1/invsumresp;
		varest2=varesttemp2/sumresp;
		f=sumresp/N;
	//cout<<"estimeate "<<muest1<<" "<<muest2<<" "<<varest1<<" "<<varest2<<" "<<f<<endl;
		//implement a convergence criteria
		if (lesser){
			float dmu=abs(mulast-muest1);
			if (dmu<0.05){
			//cout<<"estimeate "<<muest1<<" "<<muest2<<" "<<varest1<<" "<<varest2<<" "<<f<<" "<<i<<endl;
				break;
			}else{
				mulast=muest1;
			}
		}else{
			float dmu=abs(mulast-muest2);
			if (dmu<0.05){
				break;
			}else{
				mulast=muest2;
			}
			
		}
	}
	
	if (muest1<=muest2){
		if (lesser){
		//	cout<<"setmuest 1 to 40"<<endl;
		//muest1=60;
		//	*var=varest1;

			//cout<<"returen lesser "<<muest1<<endl;
			return muest1;
		}else{
		  //*var=varest2;
		//cout<<"returen greater "<<muest2<<endl;
			return muest2;
		}
	}else{
		if (lesser){
			//cout<<"returen lesser "<<muest2<<endl;
			//	*var=varest2;
			return muest2;
		}else{
		  //	*var=varest2;
						//cout<<"returen greater "<<muest1<<endl;

			return muest1;
		}
	}
	return -1; 	
}

			
//int shapeModel::intensityHistInter(const volume<float>* image, const volume<short>* mask, Mesh* m,int label , vector<float>* vgraylevels)
//{
//	//returns a sorted min distance for each vertex 
//	vgraylevels->clear();
//	float dist=10000;
//	int bounds[6]={0,0,0,0,0,0};
//	//Mesh m=getTranslatedMesh(shape);
//	//Mesh m=shapes.at(shape)->getMesh();
//	//int label=getLabel(shape);
//getBounds(m,bounds,tol);
//	int bounds[6]={0,0,0,0,0,0};
	
//	getBounds(m,bounds,5);
	
//	vector <float>::iterator Iter;
//	//cout<<bounds[4]<<" "<<bounds[2]<<" "<<bounds[0]<<endl;
//	for (int n=bounds[4];n<=bounds[5];n++){
//		for (int m=bounds[2];m<=bounds[3];m++){ 
//			for (int l=bounds[0]; l<=bounds[1];l++){
//				if (mask->value(l,m,n)==label){
//					dist=image->value(l,m,n);
//					if (vgraylevels->empty()){
//					//	cout<<"empty"<<endl;
//						vgraylevels->push_back(dist);
//					}else if (dist>=vgraylevels->back()){
//					//	cout<<"back"<<endl;
//						vgraylevels->push_back(dist);
//					}else {
//					//	cout<<"insert"<<endl;
//						for ( Iter = vgraylevels->begin( ) ; Iter !=vgraylevels->end( ) ; Iter++ ){
//							if (dist<*Iter){
//								
//								vgraylevels->insert(Iter,dist);
//								break;
//							}
//							
//						}
//						
//						
//					}
					
					
				
				
//				}
//			}
//		}
//	}
//		
//	
//	return 0;
//			}
						
void shapeModel::setPlanarParameters(vector<float> best){
	vbest=best;
}
vector<float> shapeModel::getPlanarParameters(){
	return vbest;
}

void shapeModel::modelReg(int appmode, string flirtmatname, int refxsize, int refysize, int refzsize, float refxdim, float refydim, float refzdim){
	//refsize is actually target image
	//appmode==0 is for shape mdoel
	//appmode==1 is for meshes
	//appmode==2 is for conditional shape
	//make sure you've set the image dimensions/resolution prior to this cal
	vector<float> vars;	
	
	vars.push_back(0);
	
	//vector< vector<float> > iprofs;
	int numPoints=getTotalNumberOfPoints();
	Matrix MeshPts(4,numPoints);
	Matrix NewMeshPts(4,numPoints);
	Matrix flirtmat(4,4);
	
	//read in flirt matrix, uses ascii
	
	ifstream fmat;
	fmat.open(flirtmatname.c_str());
	float tmpfloat=0;
	for (int i=0; i<4;i++){
		for (int j=0; j<4;j++){
			fmat>>tmpfloat;
			flirtmat.element(i,j)=tmpfloat;
		}
	}
				flirtmat=flirtmat.i();
	//test reading of flirt matrices
//		for (int i=0; i<4;i++){
//			for (int j=0; j<4;j++){
//				cout<<flirtmat.element(i,j);
//			}
//			cout<<endl;
//		}
	
	//need to add for multiple shapes
	for (int s=0; s<getNumberOfShapes();s++){
		Mesh m1;
//		cout<<"reg shape "<<s<<endl;
	//	cout<<"get deformed mesh"<<endl;
		//deformed matruxbrigs it into flirt space (world coordinates)
		m1=getDeformedMesh(vars,s,static_cast<int>(vars.size()));
		int count=0;
	//	cout<<"transform mesh points... "<<xdim<<" "<<ydim<<" "<<zdim<<endl;
		
		
		for (vector<Mpoint*>::iterator i = m1._points.begin(); i!=m1._points.end(); i++ ){
		//	cout<<"count "<<count<<endl;	
			MeshPts.element(0,count)=(*i)->get_coord().X;
			MeshPts.element(1,count)=(*i)->get_coord().Y;
			MeshPts.element(2,count)=(*i)->get_coord().Z;
	//		MeshPts.element(0,count)=(*i)->get_coord().X/xdim;
	//		MeshPts.element(1,count)=(*i)->get_coord().Y/ydim;
	//		MeshPts.element(2,count)=(*i)->get_coord().Z/zdim;
			MeshPts.element(3,count)=1;
			
			//Pt newPt(jointAvgShapes.element(count),jointAvgShapes.element(count+1),jointAvgShapes.element(count+2));
			//(*i)->_update_coord = newPt;
			//count=count+3;
			count++;
		}
	//cout<<"mesh points loaded into matrix... "<<refxdim<<" "<<refydim<<" "<<refzdim<<endl;
		NewMeshPts=flirtmat*MeshPts;
		count=0;
		for (vector<Mpoint*>::iterator i = m1._points.begin(); i!=m1._points.end(); i++ ){
			//Pt newPt(NewMeshPts.element(0,count)*refxdim,NewMeshPts.element(1,count)*refydim,NewMeshPts.element(2,count)*refzdim);
						
						
			Pt newPt(NewMeshPts.element(0,count),NewMeshPts.element(1,count),NewMeshPts.element(2,count));

			(*i)->_update_coord = newPt;
			count++;
		}
		m1.update();
		//ste new image parameters
		setImageParameters(refxsize, refysize, refzsize, refxdim, refydim, refzdim);
		centreAndSetShapeMesh(m1,s);//actually inver translation
	//	cout<<"...transformation complete"<<endl;
		if ((appmode==0)|(appmode==2)){
	//	cout<<"register modes sahpe"<<s<<endl;
			for (int m=0; m<getNumberOfModes();m++){
				//cout<<"mode "<<m<<endl;
				//this registers modes (i.e. covariance matrix)
				vector<Vec> vmode=getShapeMode(s, m);
				for (int p=0; p<static_cast<int>(vmode.size());p++){
				//	if (m==1){
				//		cout<<"p "<<p<<endl;
				//	}
					ColumnVector mode(4);
					ColumnVector newmode(4);
					mode.element(0)=vmode.at(p).X;
					mode.element(1)=vmode.at(p).Y;
					mode.element(2)=vmode.at(p).Z;
					mode.element(3)=0;
				//	cout<<" "<<vmode.at(p).X;
					newmode=flirtmat*mode;
					vmode.at(p).X=newmode.element(0);
					vmode.at(p).Y=newmode.element(1);
					vmode.at(p).Z=newmode.element(2);
				//	cout<<" "<<vmode.at(p).X<<endl;
				}
		//		cout<<s<<" "<<m<<endl;
				vector<Vec> vmode2=getShapeMode(s, m);
				setShapeMode(s,m,vmode);
				vector<Vec> vmode3=getShapeMode(s, m);
		//		cout<<vmode2.at(0).X<<" "<<vmode3.at(0).X<<endl;

				
			}
					
		}
	
		
	}
	//only implemented for 2 structures
	if (appmode==2){//for conditonakl shape ---anitquated
		//registers the conditional matrix for I given, assumes only 1 condmat
			vector<vector<float> > condMat;
		condMat=getShape(0)->getICondPrec();
		//int rows=static_cast<int>(condMat.at(0).size());
		int cols=static_cast<int>(condMat.size());
	//	cout<<"reg rows and cols "<<rows<<" "<<cols<<endl;
		int snum=getNumberOfPoints(1); 
			//go through each row and register the covariance matrix
		
		for (int j=0;j<cols;j++){
			for (int i =0;i<snum;i=i+3){
				ColumnVector mode(4);
				ColumnVector newmode(4);
				mode.element(0)=condMat.at(j).at(i);
				mode.element(1)=condMat.at(j).at(i+1);
				mode.element(2)=condMat.at(j).at(i+2);
				mode.element(3)=0;
				newmode=flirtmat*mode;
				condMat.at(j).at(i)=newmode.element(0);
				condMat.at(j).at(i+1)=newmode.element(1);
				condMat.at(j).at(i+2)=newmode.element(2);
				
				
			}
		}
		getShape(0)->setICondPrec(condMat);
	}
//if ((appmode.value()==0)|(appmode.value()==2)){
//	model1->save(outname.value().c_str(),2,0);
//	model1->save(outname.value().c_str(),4,nummodes.value());
//}else{
//	model1->save(outname.value().c_str(),0,nummodes.value());
//}

  //save("modelRegtest",2,142);
	
}
//void shapeModel::meshReg(Mesh* m, string flirtmatname){
void shapeModel::meshReg(Mesh* m, Matrix flirtmat){

	//refsize is actually target image
	int numPoints=m->nvertices();
	Matrix MeshPts(4,numPoints);
	Matrix NewMeshPts(4,numPoints);
//	Matrix flirtmat(4,4);
	
	//read in flirt matrix, uses ascii
	
//	ifstream fmat;
//	fmat.open(flirtmatname.c_str());
//	float tmpfloat=0;
//	for (int i=0; i<4;i++){
//		for (int j=0; j<4;j++){
//			fmat>>tmpfloat;
//			flirtmat.element(i,j)=tmpfloat;
//		}
//	}
//				flirtmat=flirtmat.i();
	
	
    //	cout<<"transform mesh points..."<<endl;
	int count=0;
	for (vector<Mpoint*>::iterator i = m->_points.begin(); i!=m->_points.end(); i++ ){
		//	cout<<"count "<<count<<endl;	
		
		MeshPts.element(0,count)=(*i)->get_coord().X;
		MeshPts.element(1,count)=(*i)->get_coord().Y;
		MeshPts.element(2,count)=(*i)->get_coord().Z;
		MeshPts.element(3,count)=1;
		
		count++;
	}
	//		cout<<"mesh points loaded into matrix..."<<endl;
	NewMeshPts=flirtmat*MeshPts;
	count=0;
	for (vector<Mpoint*>::iterator i = m->_points.begin(); i!=m->_points.end(); i++ ){
		Pt newPt(NewMeshPts.element(0,count),NewMeshPts.element(1,count),NewMeshPts.element(2,count));
		(*i)->_update_coord = newPt;
		count++;
	}
	m->update();
	
	
}

void  shapeModel::worldToVoxelCoords(Mesh* m){
	for (vector<Mpoint*>::iterator i = m->_points.begin(); i!=m->_points.end(); i++ ){
		Pt newPt((*i)->get_coord().X / xdim,(*i)->get_coord().Y / ydim,(*i)->get_coord().Z /zdim);
		(*i)->_update_coord = newPt;
		
	}
	m->update();
}

void  shapeModel::voxelToWorldCoords(Mesh* m){
	for (vector<Mpoint*>::iterator i = m->_points.begin(); i!=m->_points.end(); i++ ){
		Pt newPt((*i)->get_coord().X * xdim,(*i)->get_coord().Y * ydim,(*i)->get_coord().Z *zdim);
		(*i)->_update_coord = newPt;
		
	}
	m->update();
}

void shapeModel::meshReg(Mesh* m, string flirtmatname){
//void shapeModel::meshReg(Mesh* m, Matrix flirtmat){

	//refsize is actually target image
	int numPoints=m->nvertices();
	Matrix MeshPts(4,numPoints);
	Matrix NewMeshPts(4,numPoints);
	Matrix flirtmat(4,4);
	
	//read in flirt matrix, uses ascii
	
	ifstream fmat;
	fmat.open(flirtmatname.c_str());
	float tmpfloat=0;
	for (int i=0; i<4;i++){
		for (int j=0; j<4;j++){
			fmat>>tmpfloat;
			flirtmat.element(i,j)=tmpfloat;
		//	cout<<flirtmat.element(i,j)<<" ";
		}
		//cout<<endl;
	}
				flirtmat=flirtmat.i();
	
	
    //	cout<<"transform mesh points..."<<endl;
	int count=0;
	for (vector<Mpoint*>::iterator i = m->_points.begin(); i!=m->_points.end(); i++ ){
		//	cout<<"count "<<count<<endl;	
		
		MeshPts.element(0,count)=(*i)->get_coord().X;
		MeshPts.element(1,count)=(*i)->get_coord().Y;
		MeshPts.element(2,count)=(*i)->get_coord().Z;
		MeshPts.element(3,count)=1;
		
		count++;
	}
	//		cout<<"mesh points loaded into matrix..."<<endl;
	NewMeshPts=flirtmat*MeshPts;
	count=0;
	for (vector<Mpoint*>::iterator i = m->_points.begin(); i!=m->_points.end(); i++ ){
		Pt newPt(NewMeshPts.element(0,count),NewMeshPts.element(1,count),NewMeshPts.element(2,count));
		(*i)->_update_coord = newPt;
		count++;
	}
	m->update();
	
	
}

// Uninteresting byte swapping functions

typedef struct { unsigned char a,b ; } TWObytes ;

void shapeModel::Swap_2bytes( int n , void *ar )    /* 2 bytes at a time */
{
  register int ii ;
   register TWObytes *tb = (TWObytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].b ; tb[ii].b = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

typedef struct { unsigned char a,b,c,d ; } FOURbytes ;

void shapeModel::Swap_4bytes( int n , void *ar )    /* 4 bytes at a time */
{
   register int ii ;
   register FOURbytes *tb = (FOURbytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].d ; tb[ii].d = tt ;
     tt = tb[ii].b ; tb[ii].b = tb[ii].c ; tb[ii].c = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

typedef struct { unsigned char a,b,c,d , D,C,B,A ; } EIGHTbytes ;

void shapeModel::Swap_8bytes( int n , void *ar )    /* 8 bytes at a time */
{
   register int ii ;
   register EIGHTbytes *tb = (EIGHTbytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].A ; tb[ii].A = tt ;
     tt = tb[ii].b ; tb[ii].b = tb[ii].B ; tb[ii].B = tt ;
     tt = tb[ii].c ; tb[ii].c = tb[ii].C ; tb[ii].C = tt ;
     tt = tb[ii].d ; tb[ii].d = tb[ii].D ; tb[ii].D = tt ;
   }
   return ;
}

/*---------------------------------------------------------------------------*/

typedef struct { unsigned char a,b,c,d,e,f,g,h ,
                               H,G,F,E,D,C,B,A  ; } SIXTEENbytes ;

void shapeModel::Swap_16bytes( int n , void *ar )    /* 16 bytes at a time */
{
   register int ii ;
   register SIXTEENbytes *tb = (SIXTEENbytes *)ar ;
   register unsigned char tt ;

   for( ii=0 ; ii < n ; ii++ ){
     tt = tb[ii].a ; tb[ii].a = tb[ii].A ; tb[ii].A = tt ;
     tt = tb[ii].b ; tb[ii].b = tb[ii].B ; tb[ii].B = tt ;
     tt = tb[ii].c ; tb[ii].c = tb[ii].C ; tb[ii].C = tt ;
     tt = tb[ii].d ; tb[ii].d = tb[ii].D ; tb[ii].D = tt ;

     tt = tb[ii].e ; tb[ii].e = tb[ii].E ; tb[ii].E = tt ;
     tt = tb[ii].f ; tb[ii].f = tb[ii].F ; tb[ii].F = tt ;
     tt = tb[ii].g ; tb[ii].g = tb[ii].G ; tb[ii].G = tt ;
     tt = tb[ii].h ; tb[ii].h = tb[ii].H ; tb[ii].H = tt ;
   }
   return ;
}


//This was taken directly from miscmaths
void shapeModel::Swap_Nbytes( int n , int siz , void *ar )  /* subsuming case */
{
   switch( siz ){
     case 2:  Swap_2bytes ( n , ar ) ; break ;
     case 4:  Swap_4bytes ( n , ar ) ; break ;
     case 8:  Swap_8bytes ( n , ar ) ; break ;
     case 16: Swap_16bytes( n , ar ) ; break ;
   }
   return ;
}


}
	
	
