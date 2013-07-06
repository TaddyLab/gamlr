#ifndef __LHDMOD_H__
#define __LHDMOD_H__

double lin_nllhd(int n, double *e, double *y);
double lin_grad(int n, double *x, int *o, double *e, double *xy);
double lin_intercept(int n, double *e, double *ysum);

double bin_nllhd(int n, double *e, double *y);
double bin_grad(int n, double *x, int *o, double *e, double *xy);
double bin_curve(int n, double *x, int *o, double *e, double *d);
double bin_intercept(int n, double *e, double *ysum);

double po_nllhd(int n, double *e, double *y);
double po_grad(int n, double *x, int *o, double *e, double *xy);
double po_curve(int n, double *x, int *o, double *e, double *d);
double po_intercept(int n, double *e, double *ysum);

#endif
