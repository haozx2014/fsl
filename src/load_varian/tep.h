#include "newmat.h"
using namespace NEWMAT;
float tep_cost_func(const ColumnVector &input);
void tep_interpolate(float *out,float *in,long points,long lines,float shift);
float optimum_tep(float *data,int ppl,int lpi,int ipv);
float groa_calc(float *data,int ppl,int lpi,int ipv);
float grof_calc(float *data,int ppl,int lpi,int ipv);
void tep_interpolate_complex(float *data,int ppl,int lines,float shift);
void auto_centre_kspace(float *data,int ppl,int lpi,int ipv);
int tep_first_shift(float *data,int ppl,int lpi,int ipv);
