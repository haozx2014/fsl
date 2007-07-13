
#include "mshape.h"
#include "meshclass/meshclass.h"

using namespace std;
using namespace mesh;
namespace shapemodel {
  

MShape::MShape(){
	numModes=0;
	m=new Mesh();
}

MShape::~MShape(){
	delete m;
}


Mesh MShape::getMesh() {
	return *m;
}

void MShape::setMesh(Mesh mset){
	*m=mset;
}

void MShape::addModeVector( vector<Vec> v ){
  modeVecs.push_back(v);
  numModes++;
}

void MShape::setModeVector(  vector<Vec> v , int ind){
  modeVecs.at(ind)=v;
}

void MShape::setAffModeVectors( vector<float> tx, vector<float> ty,vector<float> tz,vector<float> rx,vector<float> ry,vector<float> rz ){
  txVecs=tx;
   tyVecs=ty;
   tzVecs=tz;
   rxVecs=rx;
   ryVecs=ry;
   rzVecs=rz;
}
void MShape::setAffModeVectorsTr( vector<float> tx, vector<float> ty,vector<float> tz ){
  txVecs=tx;
   tyVecs=ty;
   tzVecs=tz;
  
}
void MShape::setAffModeVectorsRot(vector<float> rx,vector<float> ry,vector<float> rz ){

   rxVecs=rx;
   ryVecs=ry;
   rzVecs=rz;
}

void MShape::setAffMean( float tx, float ty,float tz,float sc,float rx,float ry,float rz ){
  txmean=tx;
   tymean=ty;
   tzmean=tz;
   scmean=sc;
   rxmean=rx;
   rymean=ry;
   rzmean=rz;
}
void MShape::setAffMeanTr( float tx, float ty,float tz ){
  txmean=tx;
   tymean=ty;
   tzmean=tz;
 
}
void MShape::setAffMeanRot( float rx,float ry,float rz ){

   rxmean=rx;
   rymean=ry;
   rzmean=rz;
}
void MShape::setAffMeanSc( float sc){

   scmean=sc;
  
}
void MShape::setAffEigs( float tx, float ty,float tz,float sc,float rx,float ry,float rz ){
  txEig=tx;
   tyEig=ty;
   tzEig=tz;
   scEig=sc;
   rxEig=rx;
   ryEig=ry;
   rzEig=rz;
}
void MShape::setAffEigsTr( float tx, float ty,float tz){
  txEig=tx;
   tyEig=ty;
   tzEig=tz;
  
}
void MShape::setAffEigsRot( float rx,float ry,float rz ){

   rxEig=rx;
   ryEig=ry;
   rzEig=rz;
}
void MShape::setAffEigsSc( float sc ){
 
   scEig=sc;
  
}
vector<float>  MShape::getAfftxVecs(){
	return txVecs;
}
vector<float>  MShape::getAfftyVecs(){
	return tyVecs;
}
vector<float>  MShape::getAfftzVecs(){
	return tzVecs;
}
vector<float>  MShape::getAffrxVecs(){
	return rxVecs;
}
vector<float>  MShape::getAffryVecs(){
	return ryVecs;
}
vector<float>  MShape::getAffrzVecs(){
	return rzVecs;
}
float  MShape::getAfftxEig(){
	return txEig;
}
float  MShape::getAfftyEig(){
	return tyEig;
}
float  MShape::getAfftzEig(){
	return tzEig;
}
float  MShape::getAffrxEig(){
	return rxEig;
}
float  MShape::getAffryEig(){
	return ryEig;
}
float  MShape::getAffrzEig(){
	return rzEig;
}
float  MShape::getAffscEig(){
	return scEig;
}

float  MShape::getAfftyMean(){
	return txmean;
}
float  MShape::getAfftzMean(){
	return tymean;
} 
float  MShape::getAffrxMean(){
	return tzmean;
}
float  MShape::getAfftxMean(){
	return rxmean;
}

float  MShape::getAffryMean(){
	return rymean;
} 
float  MShape::getAffrzMean(){
	return rzmean;
}
float  MShape::getAffscMean(){
	return scmean;
}


 
 
 
 
 
 
 
void MShape::insertModeVector( vector<Vec> v, int ind ){
  vector< vector<Vec> >::iterator iter = modeVecs.begin();
  modeVecs.insert(iter+ind,v);
  numModes++;
}

const vector<Vec> MShape::getModeVector( int mode ){
  return modeVecs.at(mode);
}

const vector<float> MShape::getAffIModeVector( int mode ){
	//cout<<"imode vectors (mshape) "<<imodeVecs.size()<<endl;
  return affimodeVecs.at(mode);
}

const vector<float> MShape::getIModeVector( int mode ){
	//cout<<"imode vectors (mshape) "<<imodeVecs.size()<<endl;
  return imodeVecs.at(mode);
}
const vector<float> MShape::getDModeVector( int mode ){
	//cout<<"imode vectors (mshape) "<<imodeVecs.size()<<endl;
  return dmodeVecs.at(mode);
}
const vector<float> MShape::getBModeVector( int mode ){
	//cout<<"imode vectors (mshape) "<<imodeVecs.size()<<endl;
  return bmodeVecs.at(mode);
}

int MShape::getNumberOfModes(){
	return numModes;
}

int MShape::getNumberOfPoints(){
	return m->nvertices();
}
void MShape::pushPoint(Mpoint* p){
	m->_points.push_back(p);
}

Mpoint* MShape::getPoint(int x){
	return m->get_point(x);
}

Vec MShape::getShapeGlobalTrans(int mode){
	float xmin=1,ymin=1,zmin=1;
	Vec tr(0,0,0);
	vector<Vec> vtmp=modeVecs.at(mode);
	for (int i=0; i< static_cast<int> (vtmp.size()); i++){//search for mins
		if (abs(vtmp.at(i).X)<abs(xmin)){
			xmin=vtmp.at(i).X;
		}
		if (abs(vtmp.at(i).Y)<abs(ymin)){
			ymin=vtmp.at(i).Y;
		}
		if (abs(vtmp.at(i).Z)<abs(zmin)){
			zmin=vtmp.at(i).Z;
		}
	}
	tr.X=xmin;
	tr.Y=ymin;
	tr.Z=zmin;
	return tr;
}

void MShape::pushTriangle(Triangle* t){
	m->_triangles.push_back(t);
}

void MShape::saveMesh(string s){
	m->save(s.c_str(),3);
}
void MShape::setShapeIntensityProfiles( vector< vector<float> > profiles){
	iprofs=profiles;
}

void MShape::setIPP(int num){
	ipp=num;
}
int MShape::getIPP(){
	return ipp;
}

void MShape::setErrs(vector<float> num){
	errs=num;
}
vector<float> MShape::getErrs(){
	return errs;
}

void MShape::addIModeVector(vector<float> v){
	imodeVecs.push_back(v);
}
void MShape::addAffIModeVector(vector<float> v){
	affimodeVecs.push_back(v);
}
void MShape::addBModeVector(vector<float> v){
	bmodeVecs.push_back(v);
}
void MShape::addDModeVector(vector<float> v){
	dmodeVecs.push_back(v);
	numModes++;//beware can only use mode vector or DMode or else will double count modes
}
const vector<float> MShape::getIMean(){
	return imean;
}

const vector<float> MShape::getDMean(){
	return dmean;
}
const vector<float> MShape::getBMean(){
	return bmean;
}
void MShape::setIMean(vector<float> v){
	imean=v;
	}
void MShape::setDMean(vector<float> v){
	dmean=v;
	}
	
void MShape::setBMean(vector<float> v){
	bmean=v;
	}
void MShape::setICondPrecEigs(vector< vector<float> > precmat, vector<float> Eigs){
	icondprec=precmat;
	icondEigs=Eigs;
}
void MShape::setICondPrecEigsAff(vector< vector<float> > precmat, vector<float> Eigs){
	icondprecAff=precmat;
	icondEigsAff=Eigs;
}
void MShape::setICondPrec(vector< vector<float> > precmat){
	icondprec=precmat;
}
void MShape::setICondEigs( vector<float> Eigs){
	icondEigs=Eigs;
}
void MShape::setICondPrecAff(vector< vector<float> > precmat){
	icondprecAff=precmat;
}
void MShape::setICondEigsAff( vector<float> Eigs){
	icondEigsAff=Eigs;
}

void MShape::setBCondPrecEigs(vector< vector<float> > precmat, vector<float> Eigs){
	bcondprec=precmat;
	bcondEigs=Eigs;
}
void MShape::setBCondMat(vector< vector<float> > condmat){
	bcondmat=condmat;
}
void MShape::setBCondPrec(vector< vector<float> > precmat){
	bcondprec=precmat;
}
void MShape::setBCondEigs( vector<float> Eigs){
	bcondEigs=Eigs;
}

vector< vector<float> > MShape::getICondPrec(){
	return icondprec;

}

vector<float> MShape::getICondEigs(){
	return icondEigs;
}
vector< vector<float> > MShape::getICondPrecAff(){
	return icondprecAff;

}

vector<float> MShape::getICondEigsAff(){
	return icondEigsAff;
}

vector< vector<float> > MShape::getBCondPrec(){
	return bcondprec;

}
vector< vector<float> > MShape::getBCondMat(){
	return bcondmat;

}
vector<float> MShape::getBCondEigs(){
	return bcondEigs;
}
vector<float> MShape::getTstatX(){
	return tstatx;
}
vector<float> MShape::getTstatY(){
	return tstaty;
}
vector<float> MShape::getTstatZ(){
	return tstatz;
}
void MShape::setTstatX(vector<float> v){
	tstatx=v;
}
void MShape::setTstatY(vector<float> v){
	tstaty=v;
}
void MShape::setTstatZ(vector<float> v){
	tstatz=v;
}

}
