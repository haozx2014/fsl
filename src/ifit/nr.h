void covsrt(float **covar, int ma, int ia[], int mfit);
int gaussj(float **a, int n, float **b, int m);
void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	    int ma, float **alpha, float beta[], float *chisq,
	    void (*funcs)(float, float [], float *, float [], int));
int mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	    int ma, float **covar, float **alpha, float *chisq,
	    void (*funcs)(float, float [], float *, float [], int),
	    float *alamda);
float brent(float ax, float bx, float cx, float (*f)(float), float tol,
            float *xmin);
float f1dim(float x);
void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []));
void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
            float (*func)(float));
void powell_lim(float p[], float **xi, int n, float ftol, int *iter,
                float *fret,float (*func)(float []),int itmax);

