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
int n, p, N;
double nd, pd;

// pointers to arguments
double *Y = NULL;
double *xv = NULL;
int *xi = NULL;
int *xp = NULL;
double *W = NULL;

double *xxv = NULL;
int doxx;

// variables to be created
unsigned int fam;
double A;
double *B = NULL;
double *E = NULL;
double *Z = NULL;
double *V = NULL;
double *vxbar = NULL;
double *vxz = NULL;
double vsum;

unsigned int itertotal,npass;

double gam,l1pen;

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
  if(vxz){ free(vxz); vxz = NULL; }

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
  xsd = new_dvec(p);
  for(j=0; j<p; j++){
    H[j] = -nd*xbar[j]*xbar[j];
    if(doxx)
      H[j] += xxv[j*(j+1)/2 + j];
    else
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
  
  //calculate absolute grads
  for(j=0; j<p; j++)
    if(isfinite(W[j]) & (B[j]==0.0))
      ag0[j] = fabs(G[j])/xsd[j];

  // initialization  
  if(s==0){
    df0 = 1.0;
    for(j=0; j<p; j++) 
      if(W[j] == 0.0) df0++; 
    if(!isfinite(lam[0]))  
      lam[0] = dmax(ag0,p)/nd;
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
      shape = lam[s]*nd/gam;
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
    // penalty is lam[s]*nd*W[j]*fabs(B[j])*xsd[j].
    double pen,ghb;
    pen = xsd[j]*l1pen*W[j];
    ghb = (G[j] - H[j]*B[j]);
    if(fabs(ghb) < pen) dbet = -B[j];
    else dbet = -(G[j]-sign(ghb)*pen)/H[j];
  }

  return dbet;
}

/* coordinate descent for log penalized regression */
int cdsolve(double tol, int M)
{
  int t,i,j,k,dozero,dopen; 
  double dbet,imove,Bdiff,exitstat;

  // initialize
  dopen = isfinite(l1pen);
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
        if((t>0) | (V[0]==0.0)){
          vsum = reweight(n, A, E, Y, V, Z);
          if(vsum==0.0){ // perfect separation
            shout("Warning: infinite likelihood.   ");
            exitstat = 1;
            break; }
          for(j=0; j<p; j++){
            H[j] = curve(xp[j+1]-xp[j], 
                  &xv[xp[j]], &xi[xp[j]], xbar[j],
                  V, vsum, &vxbar[j]);
            vxz[j] = 0.0;
            for(i=xp[j]; i<xp[j+1]; i++)
              vxz[j] += V[xi[i]]*xv[i]*Z[xi[i]];
          }
          dbet = intercept(n, E, V, Z, vsum)-A;
          A += dbet;
          Bdiff = fabs(vsum*dbet*dbet);
        }
      }
    }


    /****** cycle through coefficients ******/
    for(j=0; j<p; j++){

      // always skip the zero sd var
      if(!isfinite(W[j])) continue;

      // skip the in-active set unless 'dozero'
      if(!dozero & (B[j]==0.0) & (W[j]>0.0)) continue;

      // update gradient
      if(doxx){
        if(V){ 
          shout("Error: only doxx for gaussians.\n"); 
          return 1; 
        }
        G[j] = -vxz[j] + A*vxbar[j]*vsum;
        int jind = j*(j+1)/2;
        for(k=0; k<j; k++)
          G[j] += xxv[jind+k]*B[k];
        for(k=j; k<p; k++)
          G[j] += xxv[k*(k+1)/2 + j]*B[k];
      } 
      else{  
        G[j] = grad(xp[j+1]-xp[j], 
              &xv[xp[j]], &xi[xp[j]], 
              vxbar[j]*vsum, vxz[j],
              A, E, V); }


      // for null model skip penalized variables
      if(!dopen & (W[j]>0.0)){ dbet = 0.0; continue; }

      // calculate the move and update
      dbet = Bmove(j);
      if(dbet!=0.0){ 
        B[j] += dbet;
        if(!doxx)
          for(i=xp[j]; i<xp[j+1]; i++)
            E[xi[i]] += xv[i]*dbet; 
        A += -vxbar[j]*dbet;
        Bdiff = fmax(Bdiff,H[j]*dbet*dbet);
        //speak("%d %d dbet=%g, Bdiff=%g\n",t,j,dbet,Bdiff);
      }
    }

    // break for intercept only linear model
    if( (fam==1) & (Bdiff==0.0) & dozero ) break;

    // iterate
    t++;

    // check for irregular exits
    if(t == M){
      shout("Warning: hit max CD iterations.  "); 
      exitstat = 1;
      break;
    }

    // check for active set update
    if(dozero == 1) dozero = 0;
    else if(Bdiff < tol) dozero = 1; 

  }

  // calc preds if they were skipped
  if(doxx){
    zero_dvec(E,n);
    for(j=0; j<p; j++)
      if(B[j]!=0)
        for(i=xp[j];i<xp[j+1];i++)
          E[xi[i]] += xv[i]*B[j];
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
            int *N_in, // length of nonzero x entries
            int *xi_in, // length-l row ids for nonzero x
            int *xp_in, // length-p+1 pointers to each column start
            double *xv_in, // nonzero x entry values
            double *y_in, // length-n y
            int *prexx, // indicator for pre-calc xx
            double *xxv_in, // dense columns of upper tri for xx
            double *eta, // length-n fixed shifts (assumed zero for gaussian)
            double *weight, // length-p weights
            int *standardize, // whether to scale penalty by sd(x_j)
            int *nlam, // length of the path
            double *delta, // path stepsize
            double *penscale,  // gamma in the GL paper
            double *thresh,  // cd convergence
            int *maxit, // cd max iterations 
            double *lambda, // output lambda
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
  N = *N_in;
  W = weight;
  E = eta;
  Y = y_in;
  xi = xi_in;
  xp = xp_in;
  xv = xv_in;

  doxx = *prexx;
  xxv = xxv_in;
  H = new_dvec(p);

  checkdata(*standardize);

  A=0.0;
  B = new_dzero(p);
  G = new_dzero(p);
  ag0 = new_dzero(p);
  gam = *penscale;

  npass = itertotal = 0;


  // some local variables
  double Lold, NLLHD, Lsat;
  int s;

  // family dependent settings
  switch( fam )
  {
    case 2:
      nllhd = &bin_nllhd;
      reweight = &bin_reweight;
      A = log(ybar/(1-ybar));
      Lsat = 0.0;
      break;
    case 3:
      nllhd = &po_nllhd;
      reweight = &po_reweight;
      A = log(ybar);
      // nonzero saturated deviance
      Lsat = ysum;
      for(int i=0; i<n; i++)
        if(Y[i]!=0) Lsat += -Y[i]*log(Y[i]);
      break;
    default: 
      fam = 1; // if it wasn't already
      nllhd = &lin_nllhd;
      A = (ysum - sum_dvec(eta,n))/nd;
      Lsat=0.0;
  }
  if(fam!=1){
    Z = new_dvec(n);
    V = new_dzero(n);
    vxbar = new_dvec(p);
    vxz = new_dvec(p);
  }
  else{ 
    Z = Y;
    vxbar = xbar; 
    vxz = new_dzero(p);
    vsum = nd;
    for(int j=0; j<p; j++)
      for(int i=xp[j]; i<xp[j+1]; i++)
          vxz[j] += xv[i]*Y[xi[i]]; 
  }

  l1pen = INFINITY;
  Lold = INFINITY;
  NLLHD =  nllhd(n, A, E, Y);

  if(*verb)
    speak("*** n=%d observations and p=%d covariates ***\n", n,p);

  // move along the path
  for(s=0; s<*nlam; s++){

    // deflate the penalty
    if(s>0)
      lambda[s] = lambda[s-1]*(*delta);
    l1pen = lambda[s]*nd;

    // run descent
    exits[s] = cdsolve(*thresh,*maxit);

    // update parameters and objective
    itertotal += npass;
    Lold = NLLHD;
    NLLHD =  nllhd(n, A, E, Y);
    deviance[s] = 2.0*(NLLHD - Lsat);
    df[s] = dof(s, lambda, NLLHD);
    alpha[s] = A;
    copy_dvec(&beta[s*p],B,p);

    if(s==0) *thresh *= deviance[0]; // relativism
    
    // gamma lasso updating
    for(int j=0; j<p; j++) 
      if(isfinite(gam)){
        if( (W[j]>0.0) & isfinite(W[j]) )
          W[j] = 1.0/(1.0+gam*fabs(B[j]));
      } else if(B[j]!=0.0){
        W[j] = 0.0;
      }

    // verbalize
    if(*verb) 
      speak("segment %d: lambda = %.4g, dev = %.4g, npass = %d\n", 
          s+1, lambda[s], deviance[s], npass);

    // exit checks
    if(deviance[s]<0.0){
      exits[s] = 1;
      shout("Warning: negative deviance.  ");
    }
    if(df[s] >= nd){
      exits[s] = 1;
      shout("Warning: saturated model.  "); 
    }
    if(exits[s]){
      shout("Finishing path early.\n");
      *nlam = s; break; }

    itime = interact(itime); 
  }

  *maxit = itertotal;
  gamlr_cleanup();
}










