#ifndef _modeVectors
#define _modeVectors

#include "meshclass/meshclass.h"

using namespace std;
using namespace mesh;


namespace shapemodel{

class MShape {
 
public:
  MShape();
  ~MShape();
  //mesh functions
  Mesh getMesh();
  void setMesh(Mesh mset);
  void pushPoint(Mpoint* p);
  void pushTriangle(Triangle* t);
  Mpoint* getPoint(int x);
  Vec getShapeGlobalTrans(int mode);
  void saveMesh(string s);
  int getNumberOfPoints();
  void setShapeIntensityProfiles(vector< vector<float> > profiles);
  void setIPP(int num);
  int getIPP();
  void setErrs(vector<float> num);
  vector<float> getErrs();
  
  const vector<float> getIMean();
  const vector<float> getDMean();
  const vector<float> getBMean();
  
  void setIMean(vector<float> v);
  void setDMean(vector<float> v);
  void setBMean(vector<float> v);
  
  //vectors functions
  void insertModeVector(vector<Vec> v, int ind);
  void setModeVector(vector<Vec> v, int ind);  
  void addModeVector(vector<Vec> v);
  void addIModeVector(vector<float> v);
  void addDModeVector(vector<float> v);
  void addBModeVector(vector<float> v);
    void addAffIModeVector(vector<float> v);

    void setAffModeVectors(vector<float> tx, vector<float> ty,vector<float> tz,vector<float> rx,vector<float> ry,vector<float> rz);  
      void setAffModeVectorsTr(vector<float> tx, vector<float> ty,vector<float> tz);  
    void setAffModeVectorsRot(vector<float> rx,vector<float> ry,vector<float> rz);  

  void setAffMean( float tx, float ty,float tz,float sc,float rx,float ry,float rz );
  void setAffMeanTr( float tx, float ty,float tz);
  void setAffMeanRot( float rx,float ry,float rz );
  void setAffMeanSc( float sc );
  
  void setAffEigs(float tx, float ty,float tz,float sc,float rx,float ry,float rz);
 void setAffEigsTr( float tx, float ty,float tz);
  void setAffEigsRot( float rx,float ry,float rz );
  void setAffEigsSc( float sc );
  int getNumberOfModes();	
  
  const vector<Vec> getModeVector(int mode);
  const vector<float> getIModeVector(int mode);
    const vector<float> getAffIModeVector(int mode);

  const vector<float> getDModeVector(int mode);
  const vector<float> getBModeVector(int mode);
  
  void setICondPrecEigs(vector< vector<float> > precmat, vector<float> Eigs);
  void setICondPrecEigsAff(vector< vector<float> > precmat, vector<float> Eigs);
  void setICondPrec(vector< vector<float> > precmat);
  void setICondEigs( vector<float> Eigs);
   void setICondPrecAff(vector< vector<float> > precmat);
  void setICondEigsAff( vector<float> Eigs);
  
  void setBCondPrecEigs(vector< vector<float> > precmat, vector<float> Eigs);
  void setBCondPrec(vector< vector<float> > precmat);
  void setBCondMat(vector< vector<float> > condmat);
  void setBCondEigs( vector<float> Eigs);
  
  vector< vector<float> > getICondPrec();
  vector<float> getICondEigs();
  vector< vector<float> > getICondPrecAff();
  vector<float> getICondEigsAff();
  
  vector< vector<float> > getBCondPrec();
  vector< vector<float> > getBCondMat();
  vector<float> getBCondEigs();
  
   vector<float> getTstatX();
	vector<float> getTstatY();
	 vector<float> getTstatZ();
	 void setTstatX(vector<float> v);
	 void setTstatY(vector<float> v);
	 void setTstatZ(vector<float> v);
	 
	 vector<float> getAfftxVecs();
	  vector<float> getAfftyVecs();
	   vector<float> getAfftzVecs(); 
	   vector<float> getAffrxVecs();
	    vector<float> getAffryVecs(); 
		vector<float> getAffrzVecs();
		
		float getAfftxEig();
	  float getAfftyEig();
	   float getAfftzEig(); 
	   float getAffrxEig();
	    float getAffryEig(); 
		float getAffrzEig();
				float getAffscEig();

		
				float getAfftxMean();
	  float getAfftyMean();
	   float getAfftzMean(); 
	   float getAffrxMean();
	    float getAffryMean(); 
		float getAffrzMean();
  	float getAffscMean();
	
private:
  Mesh* m;
  int numModes;
  int ipp;
  vector<float> errs;
  vector<float> tstatx;
  vector<float> tstaty;
  vector<float> tstatz;

  vector< vector<float> > iprofs;
  vector< vector<Vec> > modeVecs;
  vector< vector<float> > dmodeVecs;
  vector< vector<float> > imodeVecs;
  vector< vector<float> > bmodeVecs;
    vector< vector<float> > affimodeVecs;

  //include affine vecs
 
  
  
  vector<float> imean;
  vector<float> dmean;
    vector<float> bmean;

	vector< vector<float> > icondprec;
	vector< float > icondEigs;
	vector< vector<float> > icondprecAff;
	vector< float > icondEigsAff;
	vector< vector<float> > bcondprec;
	vector< vector<float> > bcondmat;
	vector< float > bcondEigs;
	 vector< float > txVecs;
  vector< float > tyVecs;
  vector< float > tzVecs;
 
  vector< float > rxVecs;
  vector< float > ryVecs;
  vector< float > rzVecs;

	float txmean;
  float tymean;
  float tzmean;
   float scmean;
  float rxmean;
  float rymean;
  float rzmean;
  	float txEig;
  float tyEig;
  float tzEig;
   float scEig;
  float rxEig;
  float ryEig;
  float rzEig;
};
}
#endif
