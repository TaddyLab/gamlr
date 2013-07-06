#include "latools.h"
#include "rhelp.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


/*
 * la_dgemm:  blas wrapper
 *
 * C := alpha*op(A)*op(B) +  beta*C.
 */

void la_dgemm(int tA, int tB, int Arow, int Acol, int Brow, int Bcol, int Crow, int Ccol,
		  double *A, double *B, double *C, double alpha, double beta)
{
	char opA, opB;
	int m, n, k, lda, ldb, ldc;

	if(tA){ opA = 'T'; m = Acol; k = Arow; lda = k; }
	else{ opA = 'N'; m = Arow; k = Acol; lda = m; }

	if(tB){ opB = 'T'; n = Brow; assert(Bcol==k); ldb = n; }
	else{ opB = 'N'; n = Bcol; assert(Brow==k); ldb = k; }

	ldc = m; assert(Crow==m); assert(Ccol==n);

	dgemm_(&opA,&opB,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}


/*
 * la_dsymm:  blas wrapper
 *
 * C := alpha*A*B + beta*C or alpha*B*A + beta*C, where A is symmetric
 */


void la_dsymm(int Alhs, int Arow, int Acol, int Brow, int Bcol, int Crow, int Ccol,
		  double *A, double *B, double *C, double alpha, double beta)
{
  int m, n, lda, ldb, ldc;
  char Atri = 'L';  // reference the lower triangle of column major A.
  char Aside;

  assert(Arow=Acol);

  m = Crow; n = Ccol; lda = Arow; ldb = Brow; ldc = Crow;
  assert(ldb == m);

  if(Alhs){  Aside='L';  assert(lda == m); }
  else{  Aside='R';  assert(lda == n); }

  dsymm_(&Aside,&Atri,&m,&n,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}


/*
 * la_dposv:  lapack wrapper
 *
 * Solves the system A*X = B, where A is symmetric AND positive definite.
 * Upon return, B is the solution and A contains its cholesky decomp.
 */

int la_dposv(int Arow, int Bcol, double *A, double *B)
{
  char cholTri = 'L'; // return cholesky decomp of col-major A in lower triangle
  int info;
  dposv_(&cholTri, &Arow, &Bcol, A, &Arow, B, &Arow, &info);
  return info; // if info = -i, i'th arg is wrong.  if info > A is not pos-def.
}


/*
 * la_dgesv:  lapack wrapper
 *
 * Solves the system A*X = B, where A is a general square matrix.
 * Upon return, B is the solution and A contains L and U from
 * the decomposition A = P * L * U (unit diag for L are not returned)
 */

int la_dgesv(int Arow, int Bcol, double *A, double *B)
{
  int info;
  int *ip = new_ivec(Arow);  /* pivot indices which define P;
			      row i was interchanged with row ip[i]*/
  dgesv_(&Arow, &Bcol, A, &Arow, ip, B, &Arow, &info);
  return info; // if info = -i, i'th arg is wrong.  if info > 0, A is not pos-def.
}


/*
 * la_dsysv:  lapack wrapper
 *
 * Solves the system A*X = B, where A is a symmetric matrix.
 * Upon return, B is block diagonal matrix D and the
 * multipliers used to obtain the factor U or L from the
 * factorization A = U*D*U**T or A = L*D*L**T as computed by DSYTRF.
 */

int la_dsysv(int Arow, int Bcol, double *A, double *B)
{
  int info;
  double wkopt;
  int lwork = -1;
  int *ipiv = new_ivec(Arow);
  char uplo = 'L';

  dsysv_(&uplo, &Arow, &Bcol, A, &Arow, ipiv, B, &Arow, &wkopt, &lwork, &info);

  lwork = (int) wkopt;
  double *work = new_dvec(lwork);

  dsysv_(&uplo, &Arow, &Bcol, A, &Arow, ipiv, B, &Arow, work, &lwork, &info);

  free(work);
  free(ipiv);
  return info;
}


/*
 * la_dpotrf:  lapack wrapper
 *
 * Replaces col-major A with the lower triangle cholesky L s.t. A = LL'
 * Upper triangle is untouched.
 *
 */

int la_dpotrf(int dim, double *A)
{
  int info;
  char cholTri = 'L'; // return cholesky decomp of col-major A in lower triangle
  dpotrf_(&cholTri,&dim,A,&dim,&info);
  return info; // if info = -i, i'th arg is wrong.  if info > 0, A is not pos-def.
}


/*
 * I'll take one brand-new ncxnr array of doubles.
 *
 */

double** new_mat(int nr, int nc){

  if(nr == 0 || nc == 0) return NULL;
  int j;
  double **mat = (double**) malloc(sizeof(double*) * (unsigned) nc);
  assert(mat);
  mat[0] = (double*) malloc(sizeof(double) *  (unsigned) (nc*nr));
  assert(mat[0]);

  for(j=1; j<nc; j++) mat[j] = mat[j-1] +  (unsigned) nr;

  return mat;
}


double** new_mat_fromv(int nr, int nc, double *v){
  int i, j;
  double **mat = new_mat(nr, nc);
  for(j=0; j<nc; j++) for(i=0; i<nr; i++) mat[j][i] = v[j*nr + i];

  return mat;
}


/*
 * and for ints
 *
 */

int** new_imat(int nr, int nc){

  if(nr == 0 || nc == 0) return NULL;
  int j;
  int **mat = (int**) malloc(sizeof(int*) * (unsigned) nc);
  assert(mat);
  mat[0] = (int*) malloc(sizeof(int) *  (unsigned) (nc*nr));
  assert(mat[0]);

  for(j=1; j<nc; j++) mat[j] = mat[j-1] +  (unsigned) nr;

  return mat;
}


int** new_imat_fromv(int nr, int nc, int *v){
  int i, j;
  int **mat = new_imat(nr, nc);
  for(j=0; j<nc; j++) for(i=0; i<nr; i++) mat[j][i] = v[j*nr + i];

  return mat;
}

/*
 * Gimme a new ncxnr array full of zeros.
 *
 */

double** new_zero_mat(int nr, int nc){
  int i, j;
  double **mat = new_mat(nr, nc);
  for(j=0; j<nc; j++) for(i=0; i<nr; i++) mat[j][i] = 0.0;
  return mat;
}


/*
 * Gimme a new ncxnr array full of val
 *
 */

double** new_val_mat(int nr, int nc, double val){
  int i, j;
  double **mat = new_mat(nr, nc);
  for(j=0; j<nc; j++) for(i=0; i<nr; i++) mat[j][i] = val;
  return mat;
}

/*
 * How about a new dimxdim ID array?
 *
 */

double** new_id_mat(int dim){
  int j;
  double **mat = new_zero_mat(dim, dim);
  for(j=0; j<dim; j++) mat[j][j] = 1.0;
  return mat;
}

/*
 * I'll have the same.
 *
 */

double** new_dup_mat(int nr, int nc, double **orig){
  int i, j;
  double **mat = new_mat(nr, nc);
  for(j=0; j<nc; j++) for(i=0; i<nr; i++) mat[j][i] = orig[j][i];
  return mat;
}


/*
 * That 2d array I created above?  I don't need it any more.
 *
 */


void delete_mat(double **mat){
 if(mat == NULL) return;
  assert(*mat);
  free(*mat);
  assert(mat);
  free(mat);
}


void delete_imat(int **mat){
 if(mat == NULL) return;
  assert(*mat);
  free(*mat);
  assert(mat);
  free(mat);
}


void zero_mat(double **mat, int nr, int nc){
  int i, j;
  for(j=0; j<nc; j++) for(i=0; i<nr; i++) mat[j][i] = 0.0;
}


void id_mat(double **mat, int dim){
  int j;
  zero_mat(mat, dim, dim);
  for(j=0; j<dim; j++) mat[j][j] = 1.0;
}

void diag_mat(double **mat, double *v, int dim){
  int j;
  zero_mat(mat, dim, dim);
  for(j=0; j<dim; j++) mat[j][j] = v[j];
}

/*
 * copy orig to mat
 */

void copy_mat(int nr, int nc, double** mat, double** orig)
{
  int i, j;
  for(j=0; j<nc; j++) for(i=0; i<nr; i++) mat[j][i] = orig[j][i];
}


/*
 * print an n x col matrix allocated as above out an opened outfile.
 * actually, this routine can print any double**
 */

void print_mat(int nr, int nc, double **mat, FILE *outfile)
{
  int i,j;
  for(i=0; i<nr; i++) {
    for(j=0; j<nc; j++) {
      if(j==nc-1) myprintf(outfile, "%g\n", mat[j][i]);
      else myprintf(outfile, "%g ", mat[j][i]);
    }
  }
}

void print_imat(int nr, int nc, int **mat, FILE *outfile)
{
  int i,j;
  for(i=0; i<nr; i++) {
    for(j=0; j<nc; j++) {
      if(j==nc-1) myprintf(outfile, "%d\n", mat[j][i]);
      else myprintf(outfile, "%d ", mat[j][i]);
    }
  }
}


/************* double vec tools ********/

/*
 * new vec of integers of length n
 */

double *new_dvec(int n)
{
  double *v;
  if(n == 0) return NULL;
  v = (double*)  malloc(sizeof(double) * (unsigned) n);
  assert(v);
  return v;
}


/*
 * allocate and return an array containing
 * the double n-seqence [from...to]
 */

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


/*
 * allocates a new int array of size n
 * and fills it with zeros
 */

double* new_dzero(int n)
{
  int i;
  double *v = new_dvec(n);
  for(i=0; i<n; i++) v[i] = 0;
  return v;
}

/*
 * allocate a new double vec of length n and copy the integer
 * contents of dv into it
 */

double* new_dup_dvec(double *v, int n)
{
  double* dv_new = new_dvec(n);
  copy_dvec(dv_new, v, n);
  return dv_new;
}


/*
 * sum:
 *
 * return the sum of the contents of the dvec
 */

double sum_dvec(double *v, int n)
{
  int i;
  if(n==0) return 0.0;
  assert(v);
  double s = 0.0;
  for(i=0; i<n; i++) s += v[i];
  return(s);
}


/*
 * zeros out v
 *
 */

void zero_dvec(double *v, int n)
{
  int i;
  for(i=0; i<n; i++) v[i] = 0;
}


/*
 * copy orig into v
 */

void copy_dvec(double *v, double *orig, int n)
{
  int i;
  if(n > 0) assert(v && orig);
  for(i=0; i<n; i++) v[i] = orig[i];
}


/*
 * printing an double vector out to outfile
 */

void print_dvec(double *v, int n, FILE *outfile)
{
  int i;
  for(i=0; i<n; i++) myprintf(outfile, "%g ", v[i]);
  myprintf(outfile, "\n");
}


double* drep(double val, int n)
{
  int i;
  double *v = new_dvec(n);
  for(i=0; i<n; i++) v[i] = val;
  return v;
}

double normalize(double *v, int n){
  int i;
  double vsum = 0.0;
  for(i=0; i<n; i++) vsum += v[i];
  if(vsum == 0.0)
    { for(i=0; i<n; i++) v[i] = 1.0/((double) n);
      return 0.0; }
  for(i=0; i<n; i++)
    v[i] = v[i]/vsum;
  return vsum;
}

double dmin(double *v, int n){
  double m = FP_INFINITE;
  for(int i=0; i<n; i++) if(v[i]<m) m=v[i];
  return m;
}

double dmax(double *v, int n){
  double m = -FP_INFINITE;
  for(int i=0; i<n; i++) if(v[i]>m) m=v[i];
  return m;
}

double dabsmin(double *v, int n){
  double m = FP_INFINITE;
  for(int i=0; i<n; i++) if(fabs(v[i])<m) m = fabs(v[i]);
  return m;
}

double dabsmax(double *v, int n){
  double m = 0.0;
  for(int i=0; i<n; i++) if(fabs(v[i])>m) m = fabs(v[i]);
  return m;
}


/*****integer vector tools******/


/*
 * new vec of integers of length n
 */

int *new_ivec(int n)
{
  int *iv;
  if(n == 0) return NULL;
  iv = (int*)  malloc(sizeof(int) * (unsigned) n);
  assert(iv);
  return iv;
}


/*
 * allocate and return an array containing
 * the integer seqence [from...to]
 */

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




/*
 * allocates a new int array of size n
 * and fills it with zeros
 */

int* new_izero(int n)
{
  int i;
  int *v = new_ivec(n);
  for(i=0; i<n; i++) v[i] = 0;
  return v;
}

/*
 * allocate a new integer vec of length n and copy the integer
 * contents of iv into it
 */

int* new_dup_ivec(int *v, int n)
{
  int* iv_new = new_ivec(n);
  copy_ivec(iv_new, v, n);
  return iv_new;
}


/*
 * sumiv:
 *
 * return the sum of the contents of the ivec
 */

int sum_ivec(int *v, int n)
{
  int i;
  if(n==0) return 0;
  assert(v);
  int s = 0;
  for(i=0; i<n; i++) s += v[i];
  return(s);
}


/*
 * zeros out v
 *
 */

void zero_ivec(int *v, int n)
{
  int i;
  for(i=0; i<n; i++) v[i] = 0;
}


/*
 * copy orig into v
 */

void copy_ivec(int *v, int *orig, int n)
{
  int i;
  if(n > 0) assert(v && orig);
  for(i=0; i<n; i++) v[i] = orig[i];
}


/*
 * printing an integer vector out to outfile
 */

void print_ivec(int *v, int n, FILE *outfile)
{
  int i;
  for(i=0; i<n; i++) myprintf(outfile, "%d ", v[i]);
  myprintf(outfile, "\n");
}


int* irep(int val, int n)
{
  int i;
  int *v = new_ivec(n);
  for(i=0; i<n; i++) v[i] = val;
  return v;
}
