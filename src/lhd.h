// likelihood functionals

#ifndef __LHDMOD_H__
#define __LHDMOD_H__

double intercept(int n, double *e, double *v, double *z, double vsum);
double sse(int n, double a, double *e, double *y, double *v);
double bin_nllhd(int n, double a, double *e, double *y, double *v);
double po_nllhd(int n, double a, double *e, double *y, double *v);

double bin_reweight(int n, double a, double *e, 
					double *y, double *v, double *z, int *vzf);
double po_reweight(int n, double a, double *e, 
					double *y, double *v, double *z, int *vzf);

#endif
