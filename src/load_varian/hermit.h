int hks_herm(float *data,int ppl,int lpi,int ipv);
int hks_herm_ref(float *data,float *ref,int ppl,int lpi,int ipv);
int hks_centre_kspace(float *data,int ppl,int lpi,int ipv);
float hks_phase_images(float *data,int ppl,int lpi,int *partial_lpi,int ipv);
float hks_phase_images_ref(float *data,float *ref,int ppl,int lpi,int *partial_lpi,int ipv);
void hks_hermitian_conjugation(float *data,int ppl,int lpi,int partial_lpi,int ipv);

