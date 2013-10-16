// likelihood functionals

#ifndef __LHDMOD_H__
#define __LHDMOD_H__

double grad(int n, double *x, int *o, 
			double vxsum, double vxz,
            double a, double *e, double *v);
double curve(int n, double *x, int *o, double xm,  
          double *v, double vsum, double *vxm);
double intercept(int n, double *e, double *v, double *z, double vsum);

double lin_nllhd(int n, double a, double *e, double *y);
double bin_nllhd(int n, double a, double *e, double *y);
double po_nllhd(int n, double a, double *e, double *y);

double bin_reweight(int n, double a, double *e, 
					double *y, double *v, double *z);
double po_reweight(int n, double a, double *e, 
					double *y, double *v, double *z);

#endif
