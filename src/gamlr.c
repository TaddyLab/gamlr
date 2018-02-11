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
double *V = NULL;
double *gam = NULL;
double *dof = NULL;

int prexx;
double *xbar = NULL;
double *vxsum = NULL;
double *vxz = NULL;
double *vxx = NULL;

// variables to be created
double *omega = NULL;
unsigned int fam;
double A;
double *B = NULL;
double *E = NULL;
double *Z = NULL;
double vsum;

double *l1fixedcost = NULL;

unsigned int npass,nrw;
double ntimeslam;

double ysum,ybar;
double *H = NULL;
double *G = NULL;
double *ag0 = NULL;

// function pointers
double (*nllhd)(int, double, double*, double*, double*) = NULL;
double (*reweight)(int, double, double*, 
                double *, double*, double*, int*) = NULL;

// old print functions, unused
// void printArray(int *ptr, size_t length)          
// {         
//   //for statment to print values using array             
//   size_t i = 0;
//   for( ; i < length; ++i )      
//     printf("%d", ptr[i]);        
// }   

// void printString(const char *ptr)          
// {         
//   //for statment to print values using array             
//   for( ; *ptr!=NULL; ++ptr)        
//     printf("%c", *ptr);        
// }         


// void printvec(double *ptr, int num){
//   double *ptrval = &ptr[0];
//   int idx;
//   for(idx = 0; idx < num; idx++){
//     printf("%f, %d \n", *(ptrval+idx), (ptrval+idx));
//   }
// }



/* global cleanup function */
void R_gamlr_cleanup(){
  if(!dirty) return;
  if(Z){ free(Z); Z = NULL; }
  if(B){ free(B); B = NULL; }
  if(G){ free(G); G = NULL; }
  if(H){ free(H); H = NULL; }
  if(ag0){ free(ag0); ag0 = NULL; }
  if(omega){ free(omega); omega = NULL; }

  dirty = 0;
}

/* The gradient descent move for given direction */
double Bmove(int j)
{
  double dbet, pen, ghb;
  if(H[j]==0) return -B[j];

  pen = l1fixedcost[j]*nd;
  if(W[j] > 0.0) pen += ntimeslam*W[j]*omega[j];
  ghb = (G[j] - H[j]*B[j]);
  if(fabs(ghb) < pen) dbet = -B[j];
  else dbet = -(G[j]-sign(ghb)*pen)/H[j];

  return dbet;
}

void donullgrad(void){
  int j = 0;
  for(j=0; j<p; j++)
    if( (W[j]>0.0) & isfinite(W[j]) & (B[j]==0.0) ){
      ag0[j] = fabs(G[j])/W[j] - l1fixedcost[j]*nd;
      if(ag0[j]<0.0) ag0[j] = 0.0; 
    }
}

double getdf(double L){
  int j;
  double dfs = 1.0;

  // penalized bit
  double shape,phi;
  if(fam==1) phi = L*2/nd; 
  else phi = 1.0;
  for(j=0; j<p; j++)
    if(isfinite(W[j])){
      if( (gam[j]==0.0) | (W[j]==0.0) ){  
        if( (B[j]!=0.0) ) dfs++;
      } else{ 
        shape = ntimeslam/gam[j];
        dfs += pgamma(ag0[j], shape/phi, phi*gam[j], 1, 0); 
      }
  }
  return(dfs);
}

void doxbar(void){
  int j = 0;
  for(j=0; j<p; j++){
      xbar[j] = 0.0;
      int i;
      for(i=xp[j]; i<xp[j+1]; i++) 
        xbar[j] += xv[i];
      xbar[j] *= 1.0/nd; }
}

void docurve(void){
  double vx;  
  int j;
  for(j=0; j<p; j++){
    H[j] = vxsum[j] = vxz[j] = 0.0;
    int i;
    for(i=xp[j]; i<xp[j+1]; i++){
      vx = V[xi[i]]*xv[i];
      vxsum[j] += vx;
      vxz[j] += vx*Z[xi[i]];
      H[j] += vx*xv[i];
    }
    H[j] += xbar[j]*(xbar[j]*vsum - 2.0*vxsum[j]); 
  }
}

void dograd(int j){
  int k;
  G[j] = -vxz[j] + A*vxsum[j]; 
  if(prexx){
    int jind = j*(j+1)/2;
    for(k=0; k<j; k++)
      G[j] += vxx[jind+k]*B[k];
    for(k=j; k<p; k++)
      G[j] += vxx[k*(k+1)/2 + j]*B[k];
  } 
  else{
    int i;
    for(i=xp[j]; i<xp[j+1]; i++) 
        G[j] += V[xi[i]]*xv[i]*E[xi[i]];
  }
}

/* coordinate descent for log penalized regression */
int cdsolve(double tol, int M, int RW)
{
  int rw,t,i,j,dozero,dopen,exitstat; 
  double dbet,Bdiff;

  // initialize
  dopen = isfinite(ntimeslam);
  Bdiff = INFINITY;
  exitstat = 0;
  dozero = 1;
  t = 0;
  rw = 0;

  // CD loop
  while( (Bdiff > tol) | dozero ){

    Bdiff = 0.0;
    if(dozero)
      if( (fam!=1) & (RW>rw) ){
          rw +=1;
          vsum = reweight(n, A, E, Y, V, Z, &exitstat);
          docurve();
          dbet = intercept(n, E, V, Z, vsum)-A;
          A += dbet;
          Bdiff = fabs(vsum*dbet*dbet);
      }

    /****** cycle through coefficients ******/
    for(j=0; j<p; j++){

      // always skip the zero sd var
      if(!isfinite(W[j])) continue;

      // skip the in-active set unless 'dozero'
      if(!dozero & (B[j]==0.0) & (W[j]>0.0)) continue;

      // update gradient
      dograd(j);

      // for null model skip penalized variables
      if(!dopen & (W[j]>0.0)){ dbet = 0.0; continue; }

      // calculate the move and update
      dbet = Bmove(j);
      if(dbet!=0.0){ 
        B[j] += dbet;
        if(!prexx)
          for(i=xp[j]; i<xp[j+1]; i++)
            E[xi[i]] += xv[i]*dbet; 
        A += -vxsum[j]*dbet/vsum;
        Bdiff = fmax(Bdiff,H[j]*dbet*dbet);
      }
    }

    // break for intercept only linear model
    if( (fam==1) & (Bdiff==0.0) & dozero ) break;

    // iterate
    t++;

    // check for max iterations
    if(t == M){
      shout("Warning: hit max CD iterations.  "); 
      exitstat = 2;
      break;
    }

    // check for active set update
    if(dozero == 1) dozero = 0;
    else if(Bdiff < tol) dozero = 1; 

  }

  // calc preds if they were skipped
  if(prexx & (N>0)){
    // got to figure out what to do with this if N=0
    zero_dvec(E,n);
    for(j=0; j<p; j++)
      if(B[j]!=0)
        for(i=xp[j];i<xp[j+1];i++)
          E[xi[i]] += xv[i]*B[j];
  }

  npass = t;
  nrw = rw; 
  return exitstat;
}

/*
 * Main Function: _gamlr
 * path estimation of adaptively penalized coefficients
 *
 */

 void R_gamlr(int *famid, // 1 gaus, 2 bin, 3 pois
            int *n_in, // nobs 
            int *p_in, // nvar
            int *N_in, // length of nonzero x entries
            int *xi_in, // length-l row ids for nonzero x
            int *xp_in, // length-p+1 pointers to each column start
            double *xv_in, // nonzero x entry values
            double *y_in, // length-n y
            int *prexx_in, // indicator for pre-calculated covariances
            double *xbar_in,  // un-weighted covariate means
            double *vxsum_in, // weighted sums of x values
            double *vxx_in, // dense columns of upper tri for XVX
            double *vxy_in, // weighted correlation between x and y
            double *eta, // length-n fixed shifts (assumed zero for gaussian)
            double *varweight, // length-p weights
            double *obsweight, // length-n weights
            int *standardize, // whether to scale penalty by sd(x_j)
            int *nlam, // length of the path
            double *delta, // path stepsize
            double *gamvec,  // gamma in the GL paper
            double *fixedcost, // additional fixed cost
            double *thresh,  // cd convergence
            int *maxit, // cd max iterations 
            int *maxrw, // max irls reweights
            double *lambda, // output lambda
            double *deviance, // output deviance
            double *dofvec, // output dof
            double *alpha,  // output intercepts
            double *beta, // output coefficients
            int *exits, // exit status.  0 is normal, 1 warn, 2 break path
            int *verb) // talk? 
 {
  dirty = 1; // flag to say the function has been called
  // time stamp for periodic R interaction
  time_t itime = time(NULL);  
//  printf("The color: %s\n", "blue");
  l1fixedcost = fixedcost;

  /** Build global variables **/
  fam = *famid;
  n = *n_in;
//  printf("%d\n", n);
  p = *p_in;
//  printf("%d\n", p);
  nd = (double) n;
  pd = (double) p;
  N = *N_in;
//  printf("%d\n", N);
  
  
  
  E = eta;
  Y = y_in;
  Z = new_dup_dvec(Y,n); // duplicates the array in Y at another address via Z
  ysum = sum_dvec(Y,n);  // take sum of all the elements of the array with pointer Y
// printf("%f\n", ysum);
  ybar = ysum/nd; // take the mean of all the elements of the array with pointer Y

  xi = xi_in;
  xp = xp_in;
  xv = xv_in;
  prexx = *prexx_in;
  xbar = xbar_in;
  vxsum = vxsum_in;
  vxx = vxx_in;
  vxz = vxy_in;

  H = new_dvec(p); // create a vector H of all 0s of length p 
  W = varweight; // loading the variable weight vector (this is just the pointer)
  omega = drep(1.0,p);  // gamma lasso adaptations, a vector of all 1s of length p
 // printvec(omega, p);
  V = obsweight; // observation weight vector- defaults to all 1 usually - length n
  //printvec(V, 10);
  vsum = sum_dvec(V,n);  // take sum of the observation weights across all samples 
 // printf("%f\n", vsum);

  if(prexx){
    int j;
    for(j=0; j<p; j++)
      H[j] = vxx[j*(j+1)/2+j] 
        + xbar[j]*(xbar[j]*vsum - 2.0*vxsum[j]); 
  } 
  else{
    doxbar();
    if(*standardize | (fam==1)) docurve(); 
  }
  //printvec(H, p);

  if(*standardize){
    int j;
    for(j=0; j<p; j++){
      if(fabs(H[j])<1e-10){ H[j]=0.0; W[j] = INFINITY; }
      else W[j] *= sqrt(H[j]/vsum);
    }
  }
  
  //printvec(W, p);

  A=0.0;
  B = new_dzero(p);
  G = new_dzero(p);
  ag0 = new_dzero(p);
  gam = gamvec;
  dof = dofvec;

  // some local variables
  double NLLHD, NLsat;
  int s;

  // family dependent settings
  switch( fam )
  {
    case 2:
      nllhd = &bin_nllhd;
      reweight = &bin_reweight;
      A = log(ybar/(1-ybar));
      NLsat = 0.0;
      break;
    case 3:
      nllhd = &po_nllhd;
      reweight = &po_reweight;
      A = log(ybar);
      // nonzero saturated negative log likelihood
      NLsat = ysum;
      int i;
      for(i=0; i<n; i++)
        if(Y[i]!=0) NLsat += -Y[i]*log(Y[i]);
      break;
    default: 
      fam = 1; // if it wasn't already
      nllhd = &sse;
      NLsat=0.0;
      A = intercept(n, E, V, Z, vsum);
      int j;
      for(j=0; j<p; j++) dograd(j);
  }

  NLLHD =  nllhd(n, A, E, Y, V);
  if(*verb)
    speak("*** n=%d observations and p=%d covariates ***\n", n,p);

  // move along the path
  for(s=0; s<*nlam; s++){

    // deflate the penalty
    if(s>0)
      lambda[s] = lambda[s-1]*(*delta);
    ntimeslam = lambda[s]*nd;

    // run descent
    exits[s] = cdsolve(*thresh,maxit[s],maxrw[s]);

    // update parameters and objective
    maxit[s] = npass;
    maxrw[s] = nrw;
    if( (s==0) | (N>0) ) 
      NLLHD =  nllhd(n, A, E, Y, V);
    deviance[s] = 2.0*(NLLHD - NLsat);

    donullgrad();
    if( (s==0) & !isfinite(lambda[0])){ 
      lambda[0] = dmax(ag0,p)/nd;
      ntimeslam = lambda[0]*nd; }

    dof[s] = getdf(NLLHD); 
    alpha[s] = A;
    copy_dvec(&beta[s*p],B,p);

    if(s==0) *thresh *= deviance[0]; // relativism
    
    // gamma lasso updating
    int j;
    for(j=0; j<p; j++) 
      if(gam[j]>0.0){
        if(isfinite(gam[j])){
          if( (W[j]>0.0) & isfinite(W[j]) ){
            omega[j] = 1.0/(1.0+gam[j]*fabs(B[j])); 
          }
        } 
        else if(B[j]!=0.0) omega[j] = 0.0; 
      }
      
    // verbalize
    if(*verb) 
      speak("segment %d: lambda = %.4g, dev = %.4g, npass = %d\n", 
          s+1, lambda[s], deviance[s], npass);

    // exit checks
    if(deviance[s]!=deviance[s]){
      exits[s] = 2;
      shout("Warning: NaN deviance.  ");
    }
    if(deviance[s]<0.0){
      exits[s] = 2;
      shout("Warning: negative deviance.  ");
    }
    if(exits[s]==2){
      shout("Finishing path early.\n");
      *nlam = s; break; }

    itime = interact(itime); 
  }

  // deviance calcs are wrong for null X
  // so we just make the last model look best
  if( (N==0) & (s>0) ) deviance[*nlam-1] = 0.0;
  R_gamlr_cleanup();
}










