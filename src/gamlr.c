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
double *Y = NULL;
double *xv = NULL;
int *xi = NULL;
int *xp = NULL;
double *W = NULL;

// variables to be created
unsigned int fam;
double A;
double *B = NULL;
double *E = NULL;
double *Z = NULL;
double *V = NULL;
double *vxbar = NULL;

unsigned int itertotal,npass;

double gam;
double L1pen;

double ysum,ybar;
double *xbar = NULL;
double *xsd = NULL;

double *H = NULL;
double *G = NULL;
double *ag0 = NULL;
double df0;

// function pointers
double (*nllhd)(int, double, double*, double*) = NULL;
double (*reweight)(int, double, double*, 
                double *, double*, double*) = NULL;

/* global cleanup function */
void gamlr_cleanup(){
  if(!dirty) return;

  if(xbar){ free(xbar); xbar = NULL; }
  if(xsd){ free(xsd); xsd = NULL; }

  if(B){ free(B); B = NULL; }
  if(G){ free(G); G = NULL; }
  if(H){ free(H); H = NULL; }
  if(ag0){ free(ag0); ag0 = NULL; }

  if(V){
    free(V); V = NULL;
    free(Z); Z = NULL;
    free(vxbar); vxbar = NULL;
  }

  dirty = 0;
}

void checkdata(int standardize){
  int i,j;

  // means
  ysum = sum_dvec(Y,n); 
  ybar = ysum/nd;
  xbar = new_dzero(p);
  for(j=0; j<p; j++){
    for(i=xp[j]; i<xp[j+1]; i++) 
      xbar[j] += xv[i];
    xbar[j] *= 1.0/nd;
  }

  // dispersion
  H = new_dvec(p);
  xsd = new_dvec(p);
  for(j=0; j<p; j++){
    H[j] = -nd*xbar[j]*xbar[j];
    for(i=xp[j]; i<xp[j+1]; i++) 
      H[j] += xv[i]*xv[i]; 
    if(H[j]==0.0){
      W[j] = INFINITY; 
      xsd[j] = 1.0; 
    }
    else xsd[j] = sqrt(H[j]/nd);
  }

  // to scale or not to scale
  if(!standardize) for(j=0; j<p; j++) xsd[j] = 1.0;

}

// calculates degrees of freedom, as well as
// other gradient dependent statistics.
double dof(int s, double *lam, double L){
  int j;
  
  // initialization  
  if(s==0){
    df0 = 1.0;
    for(j=0; j<p; j++){
      ag0[j] = fabs(G[j])/xsd[j];  
      if(W[j] == 0.0) df0++; 
    }
    if(!isfinite(lam[0])){  
      lam[0] = dmax(ag0,p)/nd;
      return df0; 
    }
  }

  double df = df0;

  // lasso 
  if(gam==0.0){
    for(j=0; j<p; j++)
      if( (B[j]!=0.0) | (W[j]==0.0) ) df ++;
    return df;
  }
  // gamma lasso
  double shape,phi;
  if(fam==1) phi = L*2/nd; 
  else phi = 1.0;
  for(j=0; j<p; j++){
    if(W[j]==0.0) df++;
    else if(isfinite(W[j])){
      if(!isfinite(W[j])) continue;
      if(B[j]==0.0) ag0[j] = fabs(G[j])/xsd[j];
      shape = L1pen/gam;
      df += pgamma(ag0[j], 
                    shape/phi, 
                    phi*gam, 
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
  if(W[j]==0.0) dbet = -G[j]/H[j]; 
  else{
    // penalty is L1pen*W[j]*fabs(B[j])*xsd[j].
    double pen,ghb;
    pen = xsd[j]*L1pen*W[j];
    ghb = (G[j] - H[j]*B[j]);
    if(fabs(ghb) < pen) dbet = -B[j];
    else dbet = -(G[j]-sign(ghb)*pen)/H[j];
  }

  return dbet;
}

/* coordinate descent for log penalized regression */
int cdsolve(double tol, int M)
{
  int t,i,j,dozero,dopen; 
  double dbet,imove,Bdiff,vsum,exitstat;

  // initialize
  vsum = nd;
  dopen = isfinite(L1pen);
  Bdiff = INFINITY;
  exitstat = 0;
  dozero = 1;
  t = 0;

  // CD loop
  while( ( (Bdiff > tol) | dozero ) & (t < M) ){

    Bdiff = 0.0;
    imove = 0.0;
    if(dozero){
      if(V){
        if((t>0) | V[0]==0.0){
          vsum = reweight(n, A, E, Y, V, Z);
          if(vsum==0.0){ // perfect separation
            shout("Infinite Likelihood.   ");
            exitstat = 1;
            break; }
          for(j=0; j<p; j++)
            H[j] = curve(xp[j+1]-xp[j], 
                  &xv[xp[j]], &xi[xp[j]], xbar[j],
                  V, vsum, &vxbar[j]);
          dbet = intercept(n, E, V, Z)-A;
          A += dbet;
          Bdiff += fabs(dbet);
        }
      }
      //speak("A[%d]=%g,Bdiff=%g\n",t,A,Bdiff);
    }


    /****** cycle through coefficients ******/
    for(j=0; j<p; j++){

      // always skip the zero sd var
      if(!isfinite(W[j])) continue;

      // skip the in-active set unless 'dozero'
      if(!dozero & (B[j]==0.0) & (W[j]>0.0)) continue;

      // update gradient 
      G[j] = grad(xp[j+1]-xp[j], 
              &xv[xp[j]], &xi[xp[j]], 
              A, E, V, Z, n, xbar[j]);

      // for null model skip penalized variables
      if(!dopen & (W[j]>0.0)){ dbet = 0.0; continue; }

      // calculate the move and update
      dbet = Bmove(j);
      if(dbet!=0.0){ 
        B[j] += dbet;
        Bdiff += fabs(dbet);
        for(i=xp[j]; i<xp[j+1]; i++)
          E[xi[i]] += xv[i]*dbet; 
        A += -vxbar[j]*dbet;
      }

    }

    // break for intercept only linear model
    if( (fam==1) & (Bdiff==0.0) & dozero ) break;

    // iterate
    t++;

    // check for irregular exits
    if(t == M){
      shout("Hit max CD iterations.  "); 
      exitstat = 1;
      break;
    }

    // check for active set update
    if(dozero == 1) dozero = 0;
    else if(Bdiff < tol) dozero = 1; 

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
            double *eta, // on input, length-n fixed shifts
            double *weight, // length-p weights
            int *standardize, // whether to scale penalty by sd(x_j)
            int *nlam, // length of the path
            double *delta, // path stepsize
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
  W = weight;
  E = eta;
  Y = y_in;
  xi = xi_in;
  xp = xp_in;
  xv = xv_in;

  checkdata(*standardize);
  npass = itertotal = 0;

  A=0.0;
  B = new_dzero(p);
  G = new_dzero(p);
  ag0 = new_dvec(p);
  gam = *penscale;

  // some local variables
  double Lold, NLLHD; 
  int s;

  // family dependent settings
  switch( fam )
  {
    case 2:
      nllhd = &bin_nllhd;
      reweight = &bin_reweight;
      A = log(ybar/(1-ybar));
      break;
    case 3:
      nllhd = &po_nllhd;
      reweight = &po_reweight;
      A = log(ybar);
      break;
    default: 
      fam = 1; // if it wasn't already
      nllhd = &lin_nllhd;
      A = (ysum - sum_dvec(eta,n))/nd;
  }
  if(fam!=1){
    Z = new_dvec(n);
    V = new_dzero(n);
    vxbar = new_dvec(n);
  }
  else{ 
    Z = Y;
    vxbar = xbar; }

  Lold = INFINITY;
  NLLHD =  nllhd(n, A, E, Y);

  if(*verb)
    speak("*** n=%d observations and p=%d covariates ***\n", n,p);

  // move along the path
  for(s=0; s<*nlam; s++){

    // deflate the penalty
    if(s>0)
      lam[s] = lam[s-1]*(*delta);
    L1pen = lam[s]*nd;

    // run descent
    exits[s] = cdsolve(*thresh,*maxit);

    // update parameters and objective
    itertotal += npass;
    Lold = NLLHD;
    NLLHD =  nllhd(n, A, E, Y);
    deviance[s] = 2.0*NLLHD;
    df[s] = dof(s, lam, NLLHD);
    alpha[s] = A;
    copy_dvec(&beta[s*p],B,p);

    // exit checks
    if(df[s] >= nd){
      exits[s] = 1;
      shout("Saturated model.  "); 
    }
    if(exits[s]){
      shout("Terminating path.\n");
      *nlam = s; break; }
    
    // gamma lasso updateing
    for(int j=0; j<p; j++) 
      if(isfinite(gam)){
        if( (W[j]>0.0) & isfinite(W[j]) )
          W[j] = 1.0/(1.0+gam*fabs(B[j]));
      } else if(B[j]!=0.0){
        W[j] = 0.0;
      }

    if(*verb) 
      speak("segment %d: lam = %.4g, dev = %.4g, npass = %d\n", 
          s+1, lam[s], deviance[s], npass);

    itime = interact(itime); 
  }

  *maxit = itertotal;
  gamlr_cleanup();
}










