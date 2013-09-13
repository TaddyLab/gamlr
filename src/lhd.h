// likelihood functionals

#ifndef __LHDMOD_H__
#define __LHDMOD_H__

double grad(int n, double *x, int *o, 
            double a, double *e, double *w, double *z,
            int N, double xm);
double curve(int n, double *x, int *o, double xm,  
          double *w, double wsum, double *wxm);
double intercept(int n, double *e, double *w, double *z);

double lin_nllhd(int n, double a, double *e, double *y);
double bin_nllhd(int n, double a, double *e, double *y);
double po_nllhd(int n, double a, double *e, double *y);

double bin_reweight(int n, double a, double *e, 
					double *y, double *w, double *z);
double po_reweight(int n, double a, double *e, 
					double *y, double *w, double *z);

#endif
