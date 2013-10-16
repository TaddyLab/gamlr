// likelihood functionals

#include <stdlib.h>
#include <stdio.h>
#include <Rmath.h>
#include "lhd.h"
#include "vec.h"

// Weighted least squares functions
double grad(int n, double *x, int *o, 
            double vxsum, double vxz,
            double a, double *e, double *v){
  double g = 0.0;
  double vi = 1.0;

  g = -vxz + a*vxsum;
  for(int i=0; i<n; i++){
    if(v) vi = v[o[i]];
    g += vi*x[i]*e[o[i]];
  }

  return g;
}


double curve(int n, double *x, int *o, double xm,  
          double *v, double vsum, double *vxm){
  double h = 0.0;
  double vi = 1.0;
  double vxs = 0.0;

  for(int i=0; i<n; i++){
    if(v) vi = v[o[i]];
    vxs += x[i]*vi;
    h += x[i]*vi*x[i];
  }
  
  // center
  h += xm*xm*vsum - 2.0*vxs*xm;

  *vxm = vxs/vsum;

  return h;
}

double intercept(int n, double *e, double *v, double *z, double vsum){
  int i;

  double rsum = 0.0;
  if(v){
    for(i=0; i<n; i++)
      rsum += v[i]*(z[i]-e[i]); 
  }
  else{ // shouldn't be called 
    rsum = sum_dvec(z,n) - sum_dvec(e,n); }

  return rsum/vsum;
}

// Negative Log LHD

double lin_nllhd(int n, double a, double *e, double *y){
  double l = 0.0;
  double r;
  for(int i=0; i<n; i++){
    r = (y[i] - a - e[i]);
    l += 0.5*r*r; 
  }
  return l;
}

double bin_nllhd(int n, double a, double *e, double *y){
  double l = 0.0;
  double f;
  for(int i=0; i<n; i++){
    f = a + e[i];
    l += -y[i]*f + log(1 + exp(f));
  }
  return l;
}

double po_nllhd(int n, double a, double *e, double *y){
  double l = 0.0;
  double f;
  for(int i=0; i<n; i++){
    f = a + e[i];
    l += exp(f) - y[i]*f;
  }
  return l;
}

// Re-weightings  
double bin_reweight(int n, double a, double *e, 
            double *y, double *v, double *z){
  double q, ee, vs, f;
  vs = 0.0;
  if(!isfinite(a)) return 0.0;
  for(int i=0; i<n; i++){
    f = a + e[i];
    ee = exp(f);
    q = ee/(1.0+ee);
    v[i] = 1.0/(2.0 + 1.0/ee + ee);
    if(v[i]==0.0) return 0.0;
    z[i] = f + (y[i]-q)/v[i];
    vs += v[i];
  }
  return vs;
}

double po_reweight(int n, double a, double *e, 
            double *y, double *v, double *z){
  double vs, f;
  vs = 0.0;
  for(int i=0; i<n; i++){
    f = a + e[i];
    v[i] = exp(f);
    z[i] = f + y[i]/v[i] - 1.0;
    if(v[i]==0.0) return 0.0;
    vs += v[i];
  }
  return vs;
}
