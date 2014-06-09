// likelihood functionals

#include <stdlib.h>
#include <stdio.h>
#include <Rmath.h>
#include "lhd.h"
#include "vec.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

// Weighted least squares functions
double grad(int n, double *x, int *o, 
            double vxsum, double vxz,
            double a, double *e, double *v){
  double vi = 1.0;
  int irls = (v[0]!=0);
  double g = -vxz + a*vxsum;
  double *vxe = NULL;

#ifdef _OPENMP
  int nt = omp_get_max_threads();
  int step = (n + nt + 1)/nt;
  //printf("nt %d, step %d, n %d\n",nt,step,n);

  if(step>1e3){
    vxe = new_dzero(nt);
    #pragma omp parallel
    {
      int moi=omp_get_thread_num();
      int debut=moi*step;
      int fin=debut+step;
      if(fin>n) fin = n;
      //printf("moi %d, debut %d, fin %d of %d\n",moi,debut,fin,n);
      for(int i=debut; i<fin; i++){
        if(irls) vi = v[o[i]];
        vxe[moi] += vi*x[i]*e[o[i]]; }
    }
  }
#endif

  if(vxe){
    //printf("omp\n");
    for(int k=0; k<nt; k++)
      g+=vxe[k];
    free(vxe);
  } 
  else{
    //printf("serial\n");
    for(int i=0; i<n; i++){
      if(irls) vi = v[o[i]];
      g += vi*x[i]*e[o[i]];
    }
  }

  return g;
}


double curve(int n, double *x, int *o, double xm,  
          double *v, double vsum, double *vxs){
  double h = 0.0;
  *vxs = 0.0;

  for(int i=0; i<n; i++){
    *vxs += x[i]*v[o[i]];
    h += x[i]*v[o[i]]*x[i];
  }
  
  // center
  h += xm*xm*vsum - 2.0*vxs[0]*xm;

  return h;
}

double intercept(int n, double *e, double *v, double *z, double vsum){
  int i;

  double rsum = 0.0;
  for(i=0; i<n; i++)
      rsum += v[i]*(z[i]-e[i]); 
  return rsum/vsum;
}

// Negative Log LHD
double sse(int n, double a, double *e, double *y, double *v){
  double l = 0.0;
  double r;
  double vi = 1.0;

  for(int i=0; i<n; i++){
    if(v[0]!=0) vi = v[i];
    r = (y[i] - a - e[i]);
    l += 0.5*r*r*vi; 
  }
  return l;
}

double bin_nllhd(int n, double a, double *e, double *y, double *v){
  double l = 0.0;
  double f;
  for(int i=0; i<n; i++){
    f = a + e[i];
    l += -y[i]*f + log(1 + exp(f));
  }
  return l;
}

double po_nllhd(int n, double a, double *e, double *y, double *v){
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
            double *y, double *v, double *z, int *vzf){
  double q, ee, vs, f;
  vs = 0.0;
  if(!isfinite(a)) return 0.0;
  for(int i=0; i<n; i++){
    f = a + e[i];
    ee = exp(f);
    q = ee/(1.0+ee);
    v[i] = 1.0/(2.0 + 1.0/ee + ee);
    if(v[i]<1e-12){ // perfect separation
      v[i] = 0.0;
      z[i] = y[i]; 
      *vzf = 1;
    }else{
      z[i] = f + (y[i]-q)/v[i];
      vs += v[i];
    }
  }
  return vs;
}

double po_reweight(int n, double a, double *e, 
            double *y, double *v, double *z, int *vzf){
  double vs, f;
  vs = 0.0;
  for(int i=0; i<n; i++){
    f = a + e[i];
    v[i] = exp(f);
   if(v[i]<1e-12){ // perfect fit
      v[i] = 0.0;
      z[i] = y[i]; 
      *vzf = 1;
    }else{
      vs += v[i];
      z[i] = f + y[i]/v[i] - 1.0;
    }
  }
  return vs;
}
