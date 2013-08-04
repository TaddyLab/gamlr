/* The Gamma Lasso -- Matt Taddy 2013 */

#include <stdlib.h>
#include <string.h>
#include <Rmath.h>
#include <time.h>

#include "vec.h"
#include "lhd.h"
#include "gui.h"

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
double A;
double *B = NULL;
double *E = NULL;
double trust;
unsigned int itertotal,npass;

double gam;
double L1pen;

double ysum,ybar;
double *xm = NULL;
double *xs = NULL;
double *xy = NULL;

double *H = NULL;
double *G = NULL;
double *ag0 = NULL;

// function pointers
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

  if(B){ free(B); B = NULL; }
  if(G){ free(G); G = NULL; }
  if(H){ free(H); H = NULL; }
  if(ag0){ free(ag0); ag0 = NULL; }

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

/* penalty cost function */
double calcC(void){
  if(!isfinite(L1pen)) return 0.0;
  double cost = 0.0;
  for(int j=0; j<p; j++)
    if( (V[j] > 0.0) & isfinite(V[j]))
      cost += L1pen*V[j]*fabs(B[j])*xs[j]; 
  return cost;
}

// calculates degrees of freedom, as well as
// other gradient dependent statistics.
double dof(double *lam, double L){
  double df =  1.0;
  int j;
  
  // initialization  
  if(!isfinite(lam[0])){
    for(j=0; j<p; j++){
      ag0[j] = fabs(G[j])/xs[j];  
      if(V[j] == 0.0) df++; }
    lam[0] = dmax(ag0,p)/nd; 
    return df; }

  // lasso 
  if(gam==0.0){
    for(j=0; j<p; j++)
      if( (B[j]!=0.0) | (V[j]==0.0) ) df ++;
    return df;
  }
  // gamma lasso
  double s,phi;
  if(fam==1) phi = L*2/nd; 
  else phi = 1.0;
  for(j=0; j<p; j++){
    if(V[j]==0.0) df++;
    else if(isfinite(V[j])){
      if(!isfinite(V[j])) continue;
      if(B[j]==0.0) ag0[j] = fabs(G[j])/xs[j];
      s = L1pen/gam;
      df += pgamma(ag0[j], 
                    s/phi,//+1.0, 
                    phi*gam,//*V[j], 
                    1, 0); 
    }
  }

  return df;
}

/* The gradient descent move for given direction */
double Bmove(int j)
{
  double dbet;
  if(H[j]==0) return -B[j];

  // unpenalized
  if(V[j]==0.0) dbet = -G[j]/H[j]; 
  else{
    double pen,ghb;
    pen = xs[j]*L1pen*V[j];
    ghb = (G[j] - H[j]*B[j]);
    if(fabs(ghb) < pen) dbet = -B[j];
    else dbet = -(G[j]-sign(ghb)*pen)/H[j];
  }
  if(fabs(dbet) > trust) dbet = sign(dbet)*trust;
 
  return dbet;
}

/* coordinate descent for log penalized regression */
int cdsolve(double tol, int M)
{
  int t,i,j,dozero,exitstat,dopen; 
  double dbet,Bdiff;
  double POST,Pold,Pdiff;

  // initialize
  dopen = isfinite(L1pen);
  POST=INFINITY;
  Bdiff = Pdiff = INFINITY;
  if(isfinite(trust)) trust = 1.0;
  exitstat=0;
  dozero=1;
  t=0;

  // CD loop
  while( ( (Pdiff > tol) | dozero ) & (t < M) ){
    Bdiff = 0.0;

    /****** loop through coefficients ******/
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
                      E, &trust); 

      // calculate the move and update
      dbet = Bmove(j);
      if(dbet!=0.0){ 
        B[j] += dbet;
        if(fam==1)
          for(i=xp[j]; i<xp[j+1]; i++) E[xi[i]] += xv[i]*dbet;
        else
          for(i=xp[j]; i<xp[j+1]; i++) E[xi[i]] *= exp(xv[i]*dbet);
        Bdiff += fabs(dbet);
      }
    }

    // draw the intercept
    A += Imove(n, E, &ysum);
    
    // break for intercept only linear model
    if( (fam==1) & (Bdiff==0.0) & dozero ) break;

    /****  iterate forward *****/

    t++;
    // trust region
    if(isfinite(trust)) trust = fmax(trust/2.0, Bdiff/pd);

    //  check objective 
    Pold = POST;
    POST = calcL(n, E, y);
    if(dopen) POST += calcC();
    Pdiff = Pold - POST;
    // check for irregular exits
    if((POST!=POST) | !isfinite(POST)){
      shout("Stopped descent due to infinite posterior. \n");
      exitstat = 1;
      break;
    }
    if((Pdiff < 0.0) &  (fabs(Pdiff) > 0.01)){
      shout("Stopped descent due to divergent optimization. \n");
      exitstat = 1;
      break; 
    }
    if(t == M-1){
      shout("Hit max CD iterations.  "); 
      exitstat = 1;
      break;
    }
    // printf("t = %d: diff = %g\n", t, Pdiff);
    // printf("param: %g | ", A);
    // print_dvec(B,5, mystdout);

    // check for active set update
    if(dozero == 1) dozero = 0;
    else if(Pdiff < tol) dozero = 1; 

  } 
  npass = t;
  return exitstat;
}

/*
 * Main Function: gamlr
 *
 * path estimation of penalized coefficients
 *
 */

 void gamlr(int *famid, // 1 gaus, 2 bin, 3 pois
            int *n_in, // nobs 
            int *p_in, // nvar
            int *l_in, // length of nonzero x entries
            int *xi_in, // length-l row ids for nonzero x
            int *xp_in, // length-p+1 pointers to each column start
            double *xv_in, // nonzero x entry values
            double *y_in, // length-n y
            double *eta, // length-n fixed shifts
            double *weight, // length-p weights
            int *standardize, // whether to scale penalty by sd(x_j)
            int *nlam, // length of the path
            double *minratio, // lam_nlam/lam_1
            double *penscale,  // gamma in the GL paper
            double *thresh,  // cd convergence
            int *maxit, // cd max iterations 
            double *lam, // output lambda
            double *deviance, // output deviance
            double *df, // output df
            double *alpha,  // output intercepts
            double *beta, // output coefficients
            int *exits, // exit status.  0 is normal
            int *verb) // talk? 
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
  E = eta;
  y = y_in;
  xi = xi_in;
  xp = xp_in;
  xv = xv_in;

  checkdata(*standardize);

  trust = INFINITY;
  npass = itertotal = 0;
  
  A=0.0;
  B = new_dzero(p);
  G = new_dzero(p);
  ag0 = new_dvec(p);
  gam = *penscale;

  switch( fam )
  {
    case 2:
      calcL = &bin_nllhd;
      calcG = &bin_grad;
      Imove = &bin_intercept;
      calcH = &bin_curve;
      A = log(ybar/(1-ybar));
      for(int i=0; i<n; i++) E[i] *= exp(A);
      break;
    case 3:
      calcL = &po_nllhd;
      calcG = &po_grad;
      Imove = &po_intercept;
      calcH = &po_curve;
      A = log(ybar);
      for(int i=0; i<n; i++) E[i] *= exp(A);
      break;
    default: 
      fam = 1; // if it wasn't already
      calcL = &lin_nllhd;
      calcG = &lin_grad;
      Imove = &lin_intercept;
      A = ybar;
      for(int i=0; i<n; i++) E[i] += A;
  }

  if(*verb)
    speak("*** n=%d observations and p=%d covariates ***\n", n,p);

  double delta = exp( log(*minratio)/((double) *nlam-1) );
  double Lold = INFINITY;
  double NLLHD =  calcL(n, E, y);
  double D0;
  int s;

  for(s=0; s<*nlam; s++){

    if(s>0)
      lam[s] = lam[s-1]*delta;

    L1pen = lam[s]*nd;

    exits[s] = cdsolve(*thresh,*maxit);
    itertotal += npass;

    NLLHD =  calcL(n, E, y);
    deviance[s] = 2.0*NLLHD;
    if(s==0){
      D0 = deviance[0];
      *thresh *= D0; }
    df[s] = dof(&lam[s], NLLHD);

    // exit checks
    if((fam==2) & (deviance[s]/D0 < 0.05)){
      exits[s] = 1;
      shout("Near perfect fit.  "); }
    if((Lold - NLLHD < -0.05)){
      exits[s] = 1;
      shout("Path divergence.  "); }
    if(df[s] >= nd){
      exits[s] = 1;
      shout("Saturated model.  "); }
    if(exits[s]){
      shout("Terminating path.\n");
      *nlam = s; break; }
    
    Lold = NLLHD;
    alpha[s] = A;
    copy_dvec(&beta[s*p],B,p);

    for(int j=0; j<p; j++) 
      if(isfinite(gam)){
        if( (V[j]>0.0) & isfinite(V[j]) )
          V[j] = 1.0/(1.0+gam*fabs(B[j]));
      } else if(B[j]!=0.0){
        V[j] = 0.0;
      }


    if(*verb) 
      speak("segment %d: lam = %.4g, dev = %.4g, npasses = %d\n", 
          s+1, lam[s], deviance[s], npass);

    itime = interact(itime); 
  }

  *maxit = itertotal;
  gamlr_cleanup();
}










