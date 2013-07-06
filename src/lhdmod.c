/* Model-specific probability and likelihood functions */

#include <stdlib.h>
#include <Rmath.h>
#include "rhelp.h"
#include "lhdmod.h"


// Linear
double lin_nllhd(int n, double *e, double *y){
  double l = 0.0;
  for(int i=0; i<n; i++) l += 0.5*(y[i] - e[i])*(y[i] - e[i]);
  return l;
}

double lin_grad(int n, double *x, int *o, double *e, double *xy){
  double g = -xy[0];
  for(int i=0; i<n; i++)
    g += x[i]*e[o[i]];
  return g;
}

double lin_intercept(int n, double *e, double *ysum){  
  double iG = ysum[0];  
  double dbet;
  int i;
  for(i=0; i<n; i++) iG += -e[i];
  dbet = iG/((double) n);
  for(i=0; i<n; i++) e[i] += dbet;
  return dbet;
}

// binomial
double bin_nllhd(int n, double *e, double *y){
  double l = 0.0;
  for(int i=0; i<n; i++)
    l += -y[i]*log(e[i]) + log(1 + e[i]);
  return l;
}

double bin_grad(int n, double *x, int *o, double *e, double *xy){
  double g = -xy[0];
  for(int i=0; i<n; i++) 
    g += x[i]*e[o[i]]/(1.0+e[o[i]]);
  return g;
}

double bin_intercept(int n, double *e, double *ysum){  
  int i; 
  double iG, iH, dbet, q;
  iG = -ysum[0];
  iH = 0.0;

  for(i=0; i<n; i++){ 
    q = e[i]/(1.0 + e[i]);
    iG += q;
    iH += q*(1-q); }
  dbet = -iG/iH;
  for(i=0; i<n; i++) e[i] *= exp(dbet);
  return dbet;
}

double bin_curve(int n, double *x, int *o, double *e, double *d){
  double h = 0.0;
  double ed, edx;
  int tr = isfinite(*d);
  for(int i=0; i<n; i++){
    if(tr){
      edx = exp(d[0]*fabs(x[i]));
      if(e[o[i]] > 1.0) ed = fmax(1.0, e[o[i]]/edx);
      else ed = fmax(1.0, edx/e[o[i]]); 
    } else ed = e[o[i]];
    h += x[i]*x[i]/(2.0 + ed + 1.0/ed);
  }
  return h;
 }

// Poisson 
double po_nllhd(int n, double *e, double *y){
  double l = 0.0;
  for(int i=0; i<n; i++) l += e[i] - y[i]*log(e[i]);
  return l;
}

double po_grad(int n, double *x, int *o, double *e, double *xy){
  double g = -xy[0];
  for(int i=0; i<n; i++) g += x[i]*e[o[i]];
  return g;
}

double po_curve(int n, double *x, int *o, double *e, double *d){
  double h = 0.0;
  if(isfinite(*d)) for(int i=0; i<n; i++) 
          h += x[i]*x[i]*e[o[i]]*exp(d[0]*fabs(x[i]));
  else for(int i=0; i<n; i++) 
          h += x[i]*x[i]*e[o[i]];
  return h;
}

double po_intercept(int n, double *e, double *ysum){  
  double es, dbet;
  es = 0.0;  
  int i; 
  for(i=0; i<n; i++) es += e[i];
  dbet = log(ysum[0]) - log(es);
  for(i=0; i<n; i++) e[i] *= exp(dbet);
  return dbet;
}
