#include "newmat.h"
using namespace NEWMAT;
float entropy_cost_function(const ColumnVector &input);
float iepi_cost_function(const ColumnVector& input);
void iepi_rephase(float *in,float *out,int ppl,int lpi,int nseg, const ColumnVector &phase);
double calc_entropy(float *data,long imsize);
double calc_entropy1(float *data,long imsize);
double calc_entropy2(float *data,long imsize);
double calc_entropy3(float *data,int ppl,int lpi);
void calc_gradient(float *data,int ppl,int lpi);
float calc_max(float *data,long imsize);
float calc_min(float *data,long imsize);
void calc_hist(float *data,long imsize,float min,float max,int bins,int *hist);
float calc_sumsq(float *data,long imsize);
void phase_correct(float *in, float *out, int ppl, int lpi,int nseg,float ph1);
void entropy_ghost_reduction(float *data,int ppl,int lpi,int ipv,int nseg);
void output_entropy_cost_function(float *data,int ppl,int lpi,int ipv,int nseg);
void fix_phase_file(float *cal, int ipv);
void entropy_ghost_reduction_array(float *data,int ppl,int lpi,int ipv,
				  int nseg,float *array);
void entropy_ghost_reduction_apply(float *data,int ppl,int lpi,int ipv,
				  int nseg,float *array);

void iepi_entropy_ghost_reduction(float *data,int ppl,int lpi,int ipv,int nseg);
