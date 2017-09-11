// array tools
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "vec.h"

double *new_dvec(int n)
{
  double *v;
  if(n == 0) return NULL;
  v = (double*)  malloc(sizeof(double) * (unsigned) n);
  assert(v);
  return v;
}

double* new_dseq(double from, double to, int n)
{
  int i;

  assert( from <= to);
  double *v = new_dvec(n);
  v[0] = from;
  double by = (to-from)/((double) n-1);
  for(i=1; i<n; i++) v[i] = v[i-1] + by;
  return v;
}

double* new_dzero(int n)
{
  int i;
  double *v = new_dvec(n);
  for(i=0; i<n; i++) v[i] = 0;
  return v;
}

double* new_dup_dvec(double *v, int n)
{
  double* dv_new = new_dvec(n);
  copy_dvec(dv_new, v, n);
  return dv_new;
}

double sum_dvec(double *v, int n)
{
  int i;
  if(n==0) return 0.0;
  assert(v);
  double s = 0.0;
  for(i=0; i<n; i++) s += v[i];
  return(s);
}

void zero_dvec(double *v, int n)
{
  int i;
  for(i=0; i<n; i++) v[i] = 0;
}

void copy_dvec(double *v, double *orig, int n)
{
  int i;
  if(n > 0) assert(v && orig);
  for(i=0; i<n; i++) v[i] = orig[i];
}


double* drep(double val, int n)
{
  int i;
  double *v = new_dvec(n);
  for(i=0; i<n; i++) v[i] = val;
  return v;
}

double dmin(double *v, int n){
  double m = FP_INFINITE;
  int i;
  for(i=0; i<n; i++) if(v[i]<m) m=v[i];
  return m;
}

double dmax(double *v, int n){
  double m = -FP_INFINITE;
  int i;
  for(i=0; i<n; i++) if(v[i]>m) m=v[i];
  return m;
}

double dabsmin(double *v, int n){
  double m = FP_INFINITE;
  int i;
  for(i=0; i<n; i++) if(fabs(v[i])<m) m = fabs(v[i]);
  return m;
}

double dabsmax(double *v, int n){
  double m = 0.0;
  int i;
  for(i=0; i<n; i++) if(fabs(v[i])>m) m = fabs(v[i]);
  return m;
}

int *new_ivec(int n)
{
  int *iv;
  if(n == 0) return NULL;
  iv = (int*)  malloc(sizeof(int) * (unsigned) n);
  assert(iv);
  return iv;
}

int* new_iseq(int from, int to)
{
  int n,i;

  assert( from <= to);
  n =  (to - from) + 1;
  int *v = new_ivec(n);
  v[0] = from;
  for(i=1; i<n; i++) v[i] = v[i-1] + 1;
  return v;
}

int* new_izero(int n)
{
  int i;
  int *v = new_ivec(n);
  for(i=0; i<n; i++) v[i] = 0;
  return v;
}

int* new_dup_ivec(int *v, int n)
{
  int* iv_new = new_ivec(n);
  copy_ivec(iv_new, v, n);
  return iv_new;
}

int sum_ivec(int *v, int n)
{
  int i;
  if(n==0) return 0;
  assert(v);
  int s = 0;
  for(i=0; i<n; i++) s += v[i];
  return(s);
}

void zero_ivec(int *v, int n)
{
  int i;
  for(i=0; i<n; i++) v[i] = 0;
}

void copy_ivec(int *v, int *orig, int n)
{
  int i;
  if(n > 0) assert(v && orig);
  for(i=0; i<n; i++) v[i] = orig[i];
}

int* irep(int val, int n)
{
  int i;
  int *v = new_ivec(n);
  for(i=0; i<n; i++) v[i] = val;
  return v;
}
