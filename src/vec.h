// array tools

#ifndef __VEC_H__
#define __VEC_H__ 

#include <stdio.h>
#include <stdlib.h>

// 1D double array tools.

double* new_dvec(int n);
double* new_dseq(double from, double to, int n);
double* new_dzero(int n);
double* new_dup_dvec(double *v, int n);
double sum_dvec(double *v, int n);
void zero_dvec(double *v, int n);
void copy_dvec(double *v, double *orig, int n);
double* drep(double val, int n);
double dmin(double *v, int n);
double dmax(double *v, int n);
double dabsmin(double *v, int n);
double dabsmax(double *v, int n);

// 1D int array tools.

int* new_ivec(int n);
int* new_iseq(int from, int to);
int* new_izero(int n);
int* new_dup_ivec(int *v, int n);
int sum_ivec(int *v, int n);
void zero_ivec(int *v,  int n);
void copy_ivec(int *v, int *orig, int n);
int* irep(int val, int n);

#endif

