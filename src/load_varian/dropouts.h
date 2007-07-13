void calc_dropouts(float *drop,float *errcount,int ppl,int lpi,int ipv,int nvol,int cvol);
float dropout_xtrt(float *data,int x,int y,int z,int t,int xsize,int ysize,int zsize,int tsize);
int count_dropouts(float *dropouts,int ipv,int nt);
