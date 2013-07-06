/*  Linear algebra tools (wrappers for lapack and blas) 
 *  and some basic array tools.
 *  Everything is COLUMN MAJOR for communication with fortran.
 *  
 */

#ifndef __LATOOLS_H__
#define __LATOOLS_H__ 

#include <stdio.h>
#include <stdlib.h>

/* These will all be defined if you're compiling with R */

// fns from lapack
extern void dpotrf_(char*, int*, double*, int*, int*);
extern void dposv_(char*, int*, int*, double*, int*, double*, int*, int*);
extern void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
extern void dsysv_(char*, int*,int *,double *,int*,int*,double*,int*,double*,int*,int*);

// fns from blas
extern void dgemm_(char*, char*,  int*, int*, int*, double*,
		  double*, int*, double*, int*, double*, double*, int*);
extern void dsymm_(char*, char*, int*, int*, double*,
		  double*, int*, double*, int*, double*, double*, int*);

// friendlier versions.  

void la_dgemm(int tA, int tB, int Arow, int Acol, int Brow, int Bcol, int Crow, int Ccol,
	      double *A, double *B, double *C, double alpha, double beta);
void la_dsymm(int Alhs, int Arow, int Acol, int Brow, int Bcol, int Crow, int Ccol,
	      double *A, double *B, double *C, double alpha, double beta);

int la_dposv(int Arow, int Bcol, double *A, double *B);
int la_dgesv(int Arow, int Bcol, double *A, double *B);
int la_dsysv(int Arow, int Bcol, double *A, double *B);
int la_dpotrf(int dim, double *A);

// 2D double array tools.

double** new_mat(int nr, int nc);
double** new_mat_fromv(int nr, int nc, double *v);
int** new_imat(int nr, int nc);
int** new_imat_fromv(int nr, int nc, int *v);
double** new_zero_mat(int nr, int nc);
double** new_val_mat(int nr, int nc, double val);
double** new_id_mat(int dim);
double** new_dup_mat(int nr, int nc, double **orig);
void delete_mat(double **mat);
void delete_imat(int **mat);
void copy_mat(int nr, int nc, double** mat, double** orig);
void print_mat(int nr, int nc, double **mat, FILE *outfile);
void print_imat(int nr, int nc, int **mat, FILE *outfile);
void zero_mat(double **mat, int nr, int nc);
void id_mat(double **mat, int dim);
void diag_mat(double **mat, double *v, int dim);

// 1D double array tools.

double* new_dvec(int n);
double* new_dseq(double from, double to, int n);
double* new_dzero(int n);
double* new_dup_dvec(double *v, int n);
double sum_dvec(double *v, int n);
void zero_dvec(double *v, int n);
void copy_dvec(double *v, double *orig, int n);
void print_dvec(double *v, int n, FILE *outfile);
double* drep(double val, int n);
double normalize(double *v, int n);
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
void print_ivec(int *v, int n, FILE *outfile);
int* irep(int val, int n);

#endif

