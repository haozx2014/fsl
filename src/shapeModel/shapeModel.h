#ifndef _shapeModel
#define _shapeModel

#include <iostream>

#include <string>
#include <fstream>
#include <stdio.h>
#include <algorithm>

#include "meshclass/meshclass.h"
#include "mshape.h"
using namespace std;
using namespace mesh;


namespace shapemodel{
		class shapeModel {
		
public:
		shapeModel();
		~shapeModel();
		void clear();
			
		int load_bmv(string s="manual_input", int type=0);
		int load_bmv_binary(string s="manual_input", int type=0);
		void load_bmv_binaryInfo(string s="manual_input", int type=0);
		int load_vtk(string s="manual_input",int type=0);
		
		int getNumberOfModes() ;
		int getNumberOfShapes() ;
		int getNumberOfSubjects() ;
		void setNumberOfSubjects(int N);
		int getNumberOfPoints(int shape) ;
		int getTotalNumberOfPoints() ;
		int getLabel(int shape) ;
		int getIPP(int shape) ;
		vector<int> getLabels() ;
		void setLabels(vector<int> labs);		
		float getEigenValue(int mode) ;
		float getEigenValueI(int mode) ;
		void setEigenValues(vector<float> vals);
		void setEigenValuesI(vector<float> vals);
		float getSumEigenValues() ;
		void setIntensityProfiles( vector< vector<float> > profiles);
		vector< vector<float> > getIntensityProfiles();
		Mesh getDeformedMesh(vector<float> var, int shape, int numModes);
		Mesh getDeformedMeshAff7( vector<float> var, int shape, int numModes);
				Mesh getDeformedMeshAff6( vector<float> var, int shape, int numModes);

		//Mesh getDeformedMeshP3dof( vector<float> var, int shape, int numModes,float tx, float ty, float tz);
		Mesh getDeformedMeshP3dof( vector<float> var, int shape, int numModes,float tx, float ty, float tz);
		vector<float> getDeformedIprof(vector<float> var, int shape,int numModes);
		vector<float> getDeformedIprofAff7(vector<float> var, int shape,int numModes);
				vector<float> getDeformedIprofAff6(vector<float> var, int shape,int numModes);

	//	Mesh getTranslatedMesh( vector<float> var, int shape);
		Mesh getTranslatedMesh( int shape);
		Mesh getInverseTranslatedMesh( Mesh m );
				vector<float> getProjectModeParameters( Mesh m, int shape, int trunc, float* res );

		//void testRead();
		void setImageParameters(int sx, int sy, int sz, float dx, float dy, float dz);
		int getShapeIndex(int label) ;
		bool getIntersection();
		void addShape(MShape* shape);
		MShape* getShape(int ind) ;
		Mesh getShapeMesh(int ind);
		void setShapeMesh(int ind, Mesh m);
		void centreAndSetShapeMesh(Mesh m, int ind);
		vector<Vec> getShapeMode(int shape, int mode);
		void setShapeMode(int shape, int mode, vector<Vec> v);
		void setShapeTstatX(int shape, vector<float> v);
		void setShapeTstatY(int shape, vector<float> v);
		void setShapeTstatZ(int shape, vector<float> v);
		void save(string s, int type, int numModes);
		void save_binary(string s, int type, int numModes);

		
		Vec getModelGlobalTrans(int mode);
		vector<float> projectShape(int label, vector<float> vars, shapeModel* newModel, int Nproj) ;
		//void getBounds(int shape, int *bounds);//, float *xmax, float *ymin, float *ymax, float *zmin, float *zmax);
		vector<float> projectVectors(const int label, vector<Vec> gvec, shapeModel* newModel, int beg, int Nproj);
		void residual(volume<float>* image, const volume<short>* mask, Mesh* m, int label, int extent);
			void residualMeanOnly( float mean, volume<float>* image, volume<float>* resimage, Mesh* m, int extent);

					void residualMeanOnly(volume<float>* image, const volume<short>* mask, Mesh* m, int label, int extent);
		int cumResidual(volume<float>* image, const volume<short>* mask, int label, int exten, float* totres);

		void residual( ColumnVector Best, const volume<float>* image, volume<float>* Resimage, Mesh* m, int extent);
		void getBounds(Mesh* m, int *bounds, int extra);
		void draw_segment(volume<short>& image, const Pt& p1, const Pt& p2, int label);
		volume<short> draw_mesh(const volume<short>* image, const Mesh* m, int label);
		volume<short> make_mask_from_mesh(const volume<short>* image, Mesh* m, int label, bool binout);
		int frac_overlap(const volume<short> segIm, const volume<short> gold, string fname);
		int frac_overlap(Mesh m, const volume<short> gold, int label, string fname);
		int volumeDistance(const volume<short> segIm, const volume<short> gold,  string fname);
		float volumeDistance(const volume<short>* segIm, const volume<float>* gold);
		float volumeDistance(const volume<short>* segIm, const volume<float>* gold, int *gbounds, Mesh* M);
		void volumeBounds(const volume<float>* gold, int *gbounds);
		float EMgmm(vector<float>* vgl,bool lesser,int niter);

		float EMgmm3(vector<float>* vgl,bool lesser,int niter,vector<float>* vmu, vector<float>* vvars, float init0, float init1, float init2);
		int meshDistance(const volume<short> gold,  int shape);
		void setICondPrecMatrix(vector< vector<float> > iconds, vector<float> Eigs, int shape);
		void setPlanarParameters(vector<float> best);
		vector<float> getPlanarParameters();
		int meshDistance(const volume<short>* gold, int shape, int tol, vector<float>* vdists);
		void setOrigin(int orgx, int orgy, int orgz);
						void meshReg(Mesh* m, string flirtmatname);
						void worldToVoxelCoords(Mesh* m);
						void voxelToWorldCoords(Mesh* m);

		void modelReg(int appmode, string flirtmatname,int refxsize, int refysize, int refzsize, float refxdim, float refydim, float refzdim);
		void meshReg(Mesh* m, Matrix flirtmat);

		int intensityHist(const volume<float>* image, const volume<short>* mask, Mesh* m,int label , vector<float>* vgraylevels);
				int intensityHistMaxMin(const volume<float>* image, const volume<short>* mask, Mesh* m,int label , vector<float>* vgraylevels, float* maxint, float*minint);

				int intensityHistMult(const volume<float>* image, const volume<short>* mask, Mesh* m,vector<int> vlabel , vector<float>* vgraylevels);

		void residualMedianOnly(float median,volume<float>* image,volume<float>* resimage, Mesh* m, int extent);
		//this is taken from miscmaths
						  // Uninteresting byte swapping functions
		  //taken from mismaths
  void Swap_2bytes ( int n , void *ar ) ;
  void Swap_4bytes ( int n , void *ar ) ;
  void Swap_8bytes ( int n , void *ar ) ;
  void Swap_16bytes( int n , void *ar ) ;
  void Swap_Nbytes ( int n , int siz , void *ar ) ;
		
		
private:
		int xsize, ysize, zsize;
		float xdim, ydim, zdim;
		int bounds[6];
		//float xmin,xmax,ymin,ymax,zmin,zmax;
		int numShapes;
		int numSubjects;
		int currentMode;
		float sumEigVals;
		float sumEigValsI;
		vector<int> labels;
		vector<int> npts;
		vector<int> ipps;
		vector< vector<float> > iprofs;
		vector<float> eigenVals;
		vector<float> eigenValsI;
		vector<MShape*> shapes; 
		vector<float> vbest;
		
		vector<int> origin;
	
	};
}
#endif
