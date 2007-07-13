#include "maxlen.h"

struct fdf
{
  int nx;
  int ny;
  int nz;
  int nv;
  int bytes;
  float vx;
  float vy;
  float vz;
  float px;
  float py;
  float pz;
  float phi;
  float theta;
  float psi;
  float dmin;
  float dmax;
  float tr;
  float te;
  int array_index;
  int array_dim;
  char seqfil[MAXLEN];
  char studyid[MAXLEN];
  char storage[MAXLEN];
  float orient[9];
};

void fdf_set_vox(struct fdf *header,float vx,float vy,float vz);
void fdf_set_dim(struct fdf *header,int nx, int ny, int nz, int nv);
void fdf_set_type(struct fdf *header,const char *storage,int bytes);
void fdf_set_maxmin(struct fdf *header,float dmax,float dmin);
int write_fdf_volume(const void *data,char *fidpath,char *outpath,struct fdf *header);
void fdf_set_params(struct fdf *header,char *fidpath);
