/* log penalized regression -- Matt Taddy 2013 */

#include <stdlib.h>
#include <string.h>
#include <Rmath.h>
#include <time.h>

#include "latools.h"
#include "rhelp.h"
#include "lhdmod.h"

/**** global variables ****/

// argument variables
unsigned int dirty = 0;
int n, p, l;
double nd, pd;

// pointers to arguments
double *y = NULL;
double *xv = NULL;
int *xi = NULL;
int *xp = NULL;
double *V = NULL;

// variables to be created
int fam;
double *par = NULL;
double A, NLLHD;
double *B = NULL;
double *E = NULL;
double trbnd;
int itertotal, npass;

double gam;
unsigned int fixlam;
unsigned int subsel;

double ysum,ybar;
double *xm = NULL;
double *xs = NULL;
double *xy = NULL;

double *H = NULL;
double *D = NULL;
double *G = NULL;
double *ag0 = NULL;

double *QN0 = NULL;
double *QN1 = NULL;
double *qnE = NULL;

// function pointers
double (*myexp)(double) = (*exp);
double (*calcL)(int, double*, double*) = NULL;
double (*calcG)(int, double*, int*, double*, double*) = NULL;
double (*calcH)(int, double*, int*, double*, double*) = NULL;
double (*Imove)(int, double*, double*) = NULL;

/* global cleanup function */
void gamlr_cleanup(){
  if(!dirty) return;

  if(xm){ free(xm); xm = NULL; }
  if(xs){ free(xs); xs = NULL; }
  if(xy){ free(xy); xy = NULL; }

  if(par){ free(par); par = NULL; }
  if(B){ free(B); B = NULL; }
  if(E){ free(E); E = NULL; }
  if(D){ free(D); D = NULL; }
  if(G){ free(G); G = NULL; }
  if(H){ free(H); H = NULL; }
  if(ag0){ free(ag0); ag0 = NULL; }

  if(QN0){ free(QN0); QN0 = NULL; }
  if(QN1){ free(QN1); QN1 = NULL; }
  if(qnE){ free(qnE); qnE = NULL; }

  dirty = 0;
}

void checkdata(int standardize){
  int i,j;
  ysum = sum_dvec(y,n); 
  ybar = ysum/nd;
  xm = new_dzero(p);
  H = new_dzero(p);
  xs = new_dzero(p);
  xy = new_dzero(p);

  // centers
  for(j=0; j<p; j++){
    for(i=xp[j]; i<xp[j+1]; i++) xm[j] += xv[i];
    xm[j] = xm[j]/nd; }

  // momments
  for(j=0; j<p; j++)
    for(i=xp[j]; i<xp[j+1]; i++)
    {  H[j] += xv[i]*xv[i]; xy[j] += xv[i]*y[xi[i]]; }
      
  // check for no var
  for(j=0; j<p; j++){
    xs[j] = sqrt( (H[j]/nd - xm[j]*xm[j]) );
    if(xs[j]==0.0){
      V[j] = INFINITY; 
      xs[j] = 1.0; }
  }   

  // to scale or not to scale
  if(!standardize) for(j=0; j<p; j++) xs[j] = 1.0;
}

// calculates degrees of freedom, as well as
// other gradient dependent statistics.
double dof(double *mu){
  double df =  1.0;
  int j;
  
  // subset selection
  if(subsel){
    int jnext = 0;
    double gnext = 0.0;
    for(j=0; j<p; j++){
      if(!isfinite(V[j])){
        G[j] = (*calcG)(xp[j+1]-xp[j], &xv[xp[j]], 
                        &xi[xp[j]], E, &xy[j]); 
        ag0[j] = fabs(G[j])/xs[j];
        if(ag0[j]>gnext){
          V[jnext] = INFINITY;
          jnext = j;
          V[j] = 0.0;
          gnext = G[j]; }
      }
      else df++; }
    mu[0] = exp(-df);
    return df; }

  // lasso or log pen initialization  
  if(!isfinite(mu[0])){
    for(j=0; j<p; j++){
      ag0[j] = fabs(G[j])/xs[j];  
      if(V[j] == 0.0) df++; }
    mu[0] = dmax(ag0,p);///nd; 
    return df; }

  // lasso 
  if(fixlam){
    for(j=0; j<p; j++)
      if( (B[j]!=0.0) | (V[j]==0.0) ) df ++;
    return df;
  }
  
  // log penalty
  double s,r;
  for(j=0; j<p; j++){
    if(V[j]>0.0){
      if(!isfinite(V[j])) continue;
      if(B[j]==0.0) ag0[j] = fabs(G[j])/xs[j];
      s = par[0]*V[j];
      r = par[1]*V[j]; 
      df += pgamma(ag0[j], s, 1.0/r, 1, 0); 
    } else df++;
  }

  return df;
}

/* penalty cost function */
double calcC(double *b){
  double cost =  0.0;
  int j;
  double s,r;
  for(j=0; j<p; j++)
      if( (V[j] > 0.0) & isfinite(V[j]) ){
        if(fixlam) cost += par[0]*fabs(b[j])*xs[j]*V[j];
        else{
          s = par[0]*V[j];
          r = par[1]*V[j]; 
          cost += s*(1.0 - log(s/(r+fabs(b[j])*xs[j]))); }
      }
  return cost;
}

/* The gradient descent move for given direction */
double Bmove(int j)
{
  double dbet, ghb, l1pen;
  myassert(H[j] != 0.0); // happens for all zero xj

  // unpenalized
  if(V[j]==0.0) dbet = -G[j]/H[j]; 
  else{
    if(fixlam) l1pen = xs[j]*par[0]*V[j];
    else l1pen = xs[j]*par[0]*V[j]/(par[1]*V[j]+fabs(B[j])*xs[j]); 
    ghb = (G[j] - H[j]*B[j]);
    if(fabs(ghb) < l1pen) dbet = -B[j];
    else dbet = -(G[j]-sign(ghb)*l1pen)/H[j];
  }
  if(fabs(dbet) > trbnd) dbet = sign(dbet)*trbnd;
 
  return dbet;
}

/* Quasi-Newton update (meta: inlcudes previous qn steps)
   z0 is altered upon return as your qn solution  */
void quasi_newton(int dim, double *z0, double *z1, double *z2){
  double u,v,w;
  for(int j=0; j<dim; j++){
    u = z1[j]-z0[j];
    v = z2[j]-z1[j];
    if((u!=0) & (u!=v))
      { w = u/(u-v);
        z0[j] = (1.0-w)*z1[j] + w*z2[j]; }
        else z0[j] = z2[j];
  }
}

void QNmove(double *P){  
  int i,j;
  double qnA, Lqn, Pqn;

  // move beta
  quasi_newton(p, QN0, QN1, B);
  if(fam==1){
    for(i=0; i<n; i++) qnE[i] = A;
    for(j=0; j<p; j++)
        for(i=xp[j]; i<xp[j+1]; i++)
          qnE[xi[i]] += xv[i]*QN0[j]; }
  else{
    for(i=0; i<n; i++) qnE[i] = myexp(A);
    for(j=0; j<p; j++)
        for(i=xp[j]; i<xp[j+1]; i++)
          qnE[xi[i]] *= myexp(xv[i]*QN0[j]); }


  // update intercept and add to qnE
  qnA = A + Imove(n, qnE, &ysum);

  // check posterior
  Lqn = calcL(n, qnE, y);
  Pqn = Lqn + calcC(QN0);
      
  // printf("QN log lhd jump: %g\n", *P-Pqn);
  // printf("current: %g | ", A); print_dvec(B,p,mystdout);
  // printf("qn move: %g | ", qnA); print_dvec(QN0,p,mystdout);

  if(Pqn < *P)
    { 
      *P = Pqn;
      NLLHD = Lqn;
      copy_dvec(B, QN0, p);
      copy_dvec(E, qnE, n); 
      A = qnA;
    }
}

/* coordinate descent for log penalized regression */
int cdsolve(double tol, int M, int qn)
{
  int t,i,j,dozero,exitstat,dopen; 
  double Pnew,Pold,dbet,Pdiff,Bdiff;

  // initialize
  dopen = isfinite(par[0]);
  exitstat=0;
  dozero=1;
  t=0;
  Pdiff = tol*2.0;
  if(isfinite(trbnd)) trbnd = 1.0;
  Pnew = NLLHD = calcL(n, E, y);
  if(dopen) Pnew += calcC(B);

  // CD loop
  while( ( (Pdiff > tol) | dozero ) & (t < M) ){

    Pold = Pnew;
    Bdiff = 0.0;

    // possibly store quasi-newton stuff
    if(qn & dopen)
     {  void *tmp = QN0; 
        QN0 = QN1; 
        QN1 = tmp; 
        copy_dvec(QN1, B, p); }

    // loop through coefficients
    for(j=0; j<p; j++){

      // always skip the zero sd var
      if(!isfinite(V[j])) continue;

      // skip the in-active set unless 'dozero'
      if(!dozero & (B[j]==0.0) & (V[j]>0.0)) continue;

      // update gradient 
      G[j] = (*calcG)(xp[j+1]-xp[j], 
                    &xv[xp[j]], &xi[xp[j]], 
                    E, &xy[j]);

      // for null model skip penalized variables
      if(!dopen & (V[j]>0.0)) continue;

      // update curvature
      if((fam!=1) & dozero)
        H[j] = (*calcH)(xp[j+1]-xp[j], 
                      &xv[xp[j]], &xi[xp[j]], 
                      E, &trbnd); 

      // calculate the move and update
      dbet = Bmove(j);
      if(dbet!=0.0){ 
        Bdiff += fabs(dbet);
        B[j] += dbet;
        if(fam==1)
          for(i=xp[j]; i<xp[j+1]; i++) E[xi[i]] += xv[i]*dbet;
        else
          for(i=xp[j]; i<xp[j+1]; i++) E[xi[i]] *= myexp(xv[i]*dbet);
      }
    }

    // break for intercept only linear model
    if( (fam==1) & (Bdiff==0.0) & dozero ) break;

    // draw the intercept
    if(fam!=0)
      A += Imove(n, E, &ysum);

    // iterate and update objective 
    t++;
    Pnew = NLLHD = calcL(n, E, y);
    if(dopen) Pnew += calcC(B);

    // trust region
    if(isfinite(trbnd)) trbnd = fmax(trbnd/2.0, Bdiff/pd);

    // possibly accelerate
    if(qn & (t>2) & (t%3 == 0) & dopen) QNmove(&Pnew);

    // posterior diff
    Pdiff = Pold - Pnew;

    //printf("t = %d: log posterior drop = %g\n", t, Pdiff);
    //printf("param: %g | ", A);
    //print_dvec(B,p, mystdout);
    
    // check for irregular exits
    if((Pnew!=Pnew) | !isfinite(Pnew)){
      warning("Stopped due to infinite posterior. \n");
      exitstat = 1;
      break;
    }
    if((Pdiff < 0.0) &  (fabs(Pdiff) > tol)){
      if(!isfinite(trbnd) & (fam!=1)){
        trbnd = Bdiff/(pd+1.0); 
        Pdiff = fabs(Pdiff); }
      else{
        warning("Stopped due to divergent optimization. \n");
        exitstat = 1;
        break; } 
    }

    // check for active set update
    if(dozero == 1) dozero = 0;
    else if(Pdiff < tol) dozero = 1; 

  } 
  npass = t;
  return exitstat;
}

/*
 * Main Function: Rgamlr
 *
 * path estimation of penalized coefficients
 *
 */

 void R_gamlr(int *famid, 
              int *n_in, int *p_in, int *l_in, 
              int *xi_in, int *xp_in, double *xv_in, 
              double *y_in, double *weight, int *standardize,
              int *npen, double *pminratio, double *varpen,  
              double *thresh, int *maxit, int *qn,  
              double *mu, double *deviance, double *df, 
              double *alpha,  double *beta, 
              int *exits, int *verb)
 {
  dirty = 1; // flag to say the function has been called
  // time stamp for periodic R interaction
  time_t itime = time(NULL);  

  /** Build global variables **/
  fam = *famid;
  n = *n_in;
  p = *p_in;
  nd = (double) n;
  pd = (double) p;
  l = *l_in;
  V = weight;
  y = y_in;
  xi = xi_in;
  xp = xp_in;
  xv = xv_in;

  checkdata(*standardize);

  trbnd = INFINITY;
  npass = itertotal = 0;

  switch( fam )
  {
    case 1:
      calcL = &lin_nllhd;
      calcG = &lin_grad;
      Imove = &lin_intercept;
      A = ybar;
      E = drep(A, n);
      break;
    case 2:
      calcL = &bin_nllhd;
      calcG = &bin_grad;
      Imove = &bin_intercept;
      calcH = &bin_curve;
      A = log(ybar/(1-ybar));
      E = drep(exp(A),n);
      break;
    case 3:
      calcL = &po_nllhd;
      calcG = &po_grad;
      Imove = &po_intercept;
      calcH = &po_curve;
      A = log(ybar);
      E = drep(exp(A), n);
      break;
    default: error( "unrecognized family type." );
  }

  B = new_dzero(p);
  G = new_dzero(p);
  par = new_dvec(2); 
  ag0 = new_dvec(p);

  gam = *varpen;
  fixlam = (gam == 0.0); // lasso
  subsel = !isfinite(gam); // subset selection
  if(subsel)
    for(int j=0; j<p; j++) 
      if(V[j]>0) V[j] = INFINITY;

  if(*qn) 
    { QN0 = new_dvec(p); 
      QN1 = new_dvec(p); 
      qnE = new_dvec(n); }

  if(*verb)
    myprintf(mystdout,
      "*** regression for %d observations and %d covariates ***\n", 
      n, p);

  double delta = exp( log(*pminratio)/((double) *npen-1) );
  double Lold = INFINITY;
  NLLHD = 0.0;
  int s;

  for(s=0; s<*npen; s++){

    if(s>0)
      mu[s] = mu[s-1]*delta;

    par[0] = mu[s];//*nd;
    if(!fixlam & !subsel){ 
      par[1] = mu[s]/gam;
      par[0] *= par[1]; }

    exits[s] = cdsolve(*thresh,*maxit,*qn);
    itertotal += npass;


    if(exits[s] | (Lold < NLLHD) | (npass>=*maxit)){ 
      myprintf(mystderr, 
        "Terminating path: Did you choose the wrong response family?\n");
      *npen = s; break; }

    deviance[s] = 2.0*NLLHD;
    alpha[s] = A;
    copy_dvec(&beta[s*p],B,p);
    if(s==0)
      *thresh *= fabs(deviance[0]); 

    df[s] = dof(&mu[s]);

    if(*verb) 
      myprintf(mystdout, 
          "segment %d: mu = %.4g, dev = %.4g, npasses = %d\n", 
          s+1, mu[s], deviance[s], npass);

    itime = my_r_process_events(itime); 
  }

  *maxit = itertotal;
  gamlr_cleanup();
}










