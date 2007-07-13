int pss_to_array(char* filename,float** pss_array);
void pss_reorder(float* data,long ppi,int ipv,float* array);
void copyslice(float *in,float *out,long pts);
int output_pss(char* infile,char* outfile);
int mean_pss(char* filename,float* pss);
