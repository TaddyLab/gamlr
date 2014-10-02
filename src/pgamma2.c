/*Nelson Auner: taken from http://compbio.mit.edu/spimap/pub/spimap/src/gamma.cpp
and edited Oct. 2014 */ 

/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gamma distribution

=============================================================================*/

// c++ headers
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// gsl
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>


// spidir headers
#include "common.h"
#include "gamma.h"

namespace spidir {


extern "C" {

/* natural log of the gamma function
   gammln as implemented in the
 * first edition of Numerical Recipes in C */
double gammln(double xx)
{
    double x, tmp, ser;
    const static double cof[6]={76.18009172947146,    -86.50532032941677,
                                24.01409824083091,    -1.231739572450155,
                                0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.5066282746310005*ser);
}

//
// Lanczos approximation to the gamma function. 
// 
// found on http://www.rskey.org/gamma.htm   
//
double gamm(double x) 
{
    double ret = (1.000000000190015 + 
                 76.18009172947146 / (x + 1) +  
                 -86.50532032941677 / (x + 2) + 
                 24.01409824083091 / (x + 3) +  
                 -1.231739572450155 / (x + 4) + 
                 1.208650973866179e-3 / (x + 5) + 
                 -5.395239384953e-6 / (x + 6));
    
    return ret * sqrt(2*M_PI)/x * pow(x + 5.5, x+.5) * exp(-x-5.5);
}


// natural log of the gamma distribution PDF
double gammalog(double x, double a, double b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return -x * b + (a - 1.0) * log(x) + a * log(b) - gammln(a);
}

// the gamma distribution PDF
double gammaPdf(double x, double a, double b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return exp(-x * b) * pow(x, a - 1.0) * pow(b, a) / gamm(a);
}


// Inverse gamma distribution PDF
double invgammaPdf(double x, double a, double b)
{
    return exp(loginvgammaPdf(x, a, b));
    //return pow(b, a) / gamm(a) * pow(1/x, a+1) * exp(-b/x);
}


// Log Inverse gamma distribution PDF
double loginvgammaPdf(double x, double a, double b)
{
    return a*log(b) - gammln(a) + (a+1)*log(1/x) - b/x;
}


double invgammaDerivA(double x, double a, double b)
{
    return gammaDerivA(1/x, a, b) / x / x;
}


double invgammaDerivB(double x, double a, double b)
{
    return gammaDerivB(1/x, a, b) / x / x;
}


double invgammaDerivG(double x, double d)
{
    double logd = log(d);
    double logix = log(1/x);
    double p = - gammln(1+d) + d*logd - d/x + (3+d)*logix;
    return exp(p) * (d * (-1 + x) + x + d*x*(logd + logix) - 
                     d*x*gsl_sf_psi_n(0, 1+d));
    /*
    return (1/gamm(1+d) * pow(d,d) * exp(-d/x) * pow(1/x, 3+d)) * 
            (d * (-1 + x) + x + d*x*(log(d) + log(1/x)) - 
             d*x*gsl_sf_psi_n(0, 1+d));
    */
}

double invgammaDerivG2(double x, double d)
{

    double logd = log(d);
    double logix = log(1/x);
    double dxlogdix = d*x*(logd + logix);
    double psi1 = gsl_sf_psi_n(0, 1 + d);
    double psi2 = gsl_sf_psi_n(1, 1 + d);

    double p = - gammln(1 + d) + d*log(d) -d/x + (4+d)*log(1/x);

    return exp(p) * (d - 2*(1 + d)*x + (3 + d)*x*x + 
         x*(logd + logix)*(2*(d*(-1 + x) + x) + dxlogdix) + 
         x*(-2*(d*(-1 + x) + x + dxlogdix)*
            psi1 + d*x*psi1*psi1 - d*x*psi2));
    return exp(p);
    
    /*
    return (1/gamm(1 + d))*pow(d,d)*exp(-d/x)*pow(1/x, 4 + d) *
        (d - 2*(1 + d)*x + (3 + d)*x*x + 
         x*(log(d) + log(1/x))*(2*(d*(-1 + x) + x) + 
                                d*x*(log(d) + log(1/x))) + 
         x*(-2*(d*(-1 + x) + x + d*x*(log(d) + log(1/x)))*
            gsl_sf_psi_n(0, 1 + d) + 
            d*x*pow(gsl_sf_psi_n(0, 1 + d), 2) - 
            d*x*gsl_sf_psi_n(1, 1 + d)));
    */
}


// Derivative of Gamma distribution with respect to x
double gammaDerivX(double x, double a, double b)
{ 
    return pow(b, a) / gamm(a) * exp(-b*x) * pow(x, a-2) * (a - b*x - 1);
}
    
// Derivative of Gamma distribution with respect to a
double gammaDerivA(double x, double a, double b)
{
    return exp(-b*x) * pow(x,a-1) * pow(b, a) / gamm(a) * 
        (log(b) + log(x) - gsl_sf_psi_n(0, a));
}

// Derivative of Gamma distribution with respect to b
double gammaDerivB(double x, double a, double b)
{
    return pow(x, a-1) / gamm(a) * exp(-b*x) * pow(b, a-1) * (a - x*b);
}


// Derivative of Gamma distribution with respect to nu (its variance)
double gammaDerivV(double x, double v)
{
    float q = 1.0 / v;
    
    double lnq = log(q);
    double lnx = log(x);
    double w = lnq*(q+2) + lnx*(q-1) - x * q - gammln(q);
    double z = 1 - x + lnq + lnx - gsl_sf_psi_n(0, q);

    if (z > 0)
	return -exp(w + log(z));
    else
	return exp(w + log(-z));
    

    //return - pow(q, q+2) * pow(x, q-1) * exp(-x * q) / gamm(q) * 
    //	(1 - x + log(q) + log(x) - gsl_sf_psi_n(0, q));
}

// Second Derivative of Gamma distribution with respect to nu (its variance)
double gammaDerivV2(double x, double v)
{
    double q = 1.0 / v;
    double lnx = log(x);
    double lnq = log(q);
    double A = 1.0 + v - x + lnx;
    double psi0 = gsl_sf_psi_n(0, q);
    double psi1 = gsl_sf_psi_n(1, q);

    double w = -x*q + lnq*(q+4) + lnx*(q-1) - gammln(q);
    double z = 1 + 3*v - 2*x - 2*v*x + x*x + lnq*lnq - lnx*lnx +
	2*A*lnx + 2*A*lnq -
	2*(A+lnq)*psi0 + psi0*psi0 - psi1;
    double y;

    if (z > 0.0)
	y = exp(w + log(z));
    else
	y = -exp(w + log(-z));
    
    //printf("x=%f, v=%f, y=%f %f; %f %f\n", x, v, y, gammln(q), w, z);

    //assert(y < INFINITY);
    //assert(y > -INFINITY);
    //assert(!isnan(y));

    return y;

    /*
    return exp(-x*q) * pow(q,q+4)*pow(x,q-1)/gamm(q) *
	(1 + 3*v - 2*x - 2*v*x + x*x + lnq*lnq - lnx*lnx +
	 2*A*lnx + 2*A*lnq -
	 2*(A+lnq)*psi0 + psi0*psi0 - psi1);
    */
}


// PDF of a sum of n gamma variables
double gammaSumPdf(double y, int n, float *alpha, float *beta, 
		   float tol)
{
    const int nterms = 100; // maximum number of terms to compute

    // y must be positive
    if (y <= 0.0)
        return 0.0;

    // handle case where a single gamma is present
    if (n == 1)
	return gammaPdf(y, alpha[0], beta[0]);

    double b1 = 1.0 / beta[0];
    double beta2[n];
    for (int i=0; i<n; i++) {

	// convert betas into their reciprocals
	beta2[i] = 1.0 / beta[i];

	// find min beta (b1)
	if (beta2[i] < b1)
	    b1 = beta2[i];
    }

    // compute terms C and p
    double C = 1.0;
    double p = 0.0;
    for (int i=0; i<n; i++) {
	C *= pow(b1 / beta2[i], alpha[i]);
	p += alpha[i];	
    }


    double prob = 0.0, prob2 = 0.0;
    double gammas[nterms]; // originally 1-indexed
    double deltas[nterms]; // zero-indexed
    deltas[0] = 1.0;
    for (int i=0; i<nterms-1; i++) {
        int k = i + 1;
        gammas[i] = 0.0;
	for (int j=0; j<n; j++)
            gammas[i] += alpha[j] * pow(1.0 - b1 / beta2[j], k) / k;
	
        k = i;
        deltas[k+1] = 0.0;
	for (int j=1; j<k+2; j++)
	    deltas[k+1] += j * gammas[j-1] * deltas[k+1-j];
	deltas[k+1] /= (k+1);

        prob2 = deltas[k] * pow(y, p+k-1) * exp(-y/b1) /
	    (gamm(p+k) * pow(b1, p+k));

        prob += prob2;

	// test for convergence
	if (prob2 / prob < tol)
	    break;
    }

    return C * prob;
}




/* =========================================================================
// UNNEEDED

double negbinomPdf(int k, double r, double p)
{
    return gsl_ran_negative_binomial_pdf(k, p, r);
}


// Log Derivative of Negative Binomial distribution with respect to r
double negbinomDerivR(int k, double r, double p)
{
    // (1-p)^k p^r G(k+r) / (k! G(r)) ( psi^(0)(k+r) - psi^(0)(r) + log(p) )
    double lnp = log(p);
    double A = k * log(1-p) + r*lnp + gammln(k+r) - gammln(k+1) - gammln(r);
    double B = gsl_sf_psi_n(0, k+r) - gsl_sf_psi_n(0, r) + lnp;
    return exp(A) * B;
}


// Derivative of Negative Binomial distribution with respect to p
double negbinomDerivP(int k, double r, double p)
{
    // (1-p)^k p^r G(k+r) / (k! G(r)) (k/(p-1) + r/p)
    double lnp = log(p);
    double A = k * log(1-p) + r*lnp + gammln(k+r) - gammln(k+1) - gammln(r);
    double B = k/(p-1) + r/p;
    return exp(A) * B;
}

*/

} // extern "C"



} // namespace spidir