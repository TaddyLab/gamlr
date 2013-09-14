// likelihood functionals

#include <stdlib.h>
#include <stdio.h>
#include <Rmath.h>
#include "lhd.h"
#include "vec.h"

// Weighted least squares functions
double grad(int n, double *x, int *o, 
            double a, double *e, double *w, double *z,
            int N, double xm){
  double g = 0.0;
  double wi = 0.5;

  for(int i=0; i<n; i++){
    if(w) wi = w[o[i]];
    g += -x[i]*wi*(z[o[i]] - a - e[o[i]]);
  }

  g *= 2.0;
  return g;
}

double curve(int n, double *x, int *o, double xm,  
          double *w, double wsum, double *wxm){
  double h = 0.0;
  double wi = 0.5;
  double wxs = 0.0;

  for(int i=0; i<n; i++){
    if(w) wi = w[o[i]];
    wxs += x[i]*wi;
    h += x[i]*wi*x[i];
  }
  
  // center
  h += xm*xm*wsum - 2.0*wxs*xm;

  *wxm = wxs/wsum;

  h *= 2.0;
  return h;
}

double intercept(int n, double *e, double *w, double *z){
  double wsum, rsum, alpha;
  int i;

  wsum = rsum = 0.0;
  if(w){
    for(i=0; i<n; i++){
      wsum += w[i];
      rsum += w[i]*(z[i]-e[i]); 
     }
  }
  else{ // shouldn't be called 
    wsum = (double) n;
    rsum = sum_dvec(z,n) - sum_dvec(e,n);
  }

  return rsum/wsum;
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
    l += exp(f) - y[i]*(f);
  }
  return l;
}

// Re-weightings  

double bin_reweight(int n, double a, double *e, 
            double *y, double *w, double *z){
  double q, ee, ws, f;
  ws = 0.0;
  for(int i=0; i<n; i++){
    f = a + e[i];
    ee = exp(f);
    q = ee/(1.0+ee);
    w[i] = q*(1.0-q);
    if(w[i]==0.0){
      ws = 0.0; break; }
    z[i] = f + (y[i]-q)/w[i];
    ws += w[i];
  }
  return ws;
}

double po_reweight(int n, double a, double *e, 
            double *y, double *w, double *z){
  double ws, f;
  ws = 0.0;
  for(int i=0; i<n; i++){
    f = a + e[i];
    w[i] = exp(f);
    z[i] = f + y[i]/w[i] - 1.0;
    ws += w[i];
  }
  return ws;
}
