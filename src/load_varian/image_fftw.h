void image_fftw_r(float*data,int ppl,int lpi,int ipv,int dir);
void image_fftw_c(float*data,int ppl,int lpi,int ipv,int dir);
void image_fftw_s(float *data,int ppl,int lpi,int ipv,float p,int dir);
void fftw_baseline(float *data,int points,int step);

#define image_oned_fftw_r(data,ppl,lpi,ipv) \
     image_fftw_r(data,ppl,lpi,ipv,1)
#define image_oned_ifftw_r(data,ppl,lpi,ipv) \
     image_fftw_r(data,ppl,lpi,ipv,-1)

#define image_oned_fftw_c(data,ppl,lpi,ipv) \
     image_fftw_c(data,ppl,lpi,ipv,1)
#define image_oned_ifftw_c(data,ppl,lpi,ipv) \
     image_fftw_c(data,ppl,lpi,ipv,-1)

#define image_oned_fftw_s(data,ppl,lpi,ipv) \
     image_fftw_s(data,ppl,lpi,ipv,0.0,1)
#define image_oned_fftw_sp(data,ppl,lpi,ipv,p) \
     image_fftw_s(data,ppl,lpi,ipv,p,1)
#define image_oned_ifftw_s(data,ppl,lpi,ipv) \
     image_fftw_s(data,ppl,lpi,ipv,0.0,-1)

