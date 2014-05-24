<<<<<<< HEAD
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// Rmath.h: Rcpp R/C++ interface class library -- Wrappers for R's Rmath API
//
// Copyright (C) 2012 Dirk Eddelbuettel and Romain Francois
//
// This file is part of Rcpp.
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or           
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#ifndef Rcpp_Rmath_h
#define Rcpp_Rmath_h

namespace R {

    // see R's Rmath.h as well as Writing R Extension

    /* Random Number Generators */
    inline double norm_rand(void) 	{ return ::norm_rand(); }
    inline double unif_rand(void)	{ return ::unif_rand(); }
    inline double exp_rand(void)	{ return ::exp_rand(); } 

    /* Normal Distribution */
    inline double dnorm(double x, double mu, double sigma, int lg)              { return ::Rf_dnorm4(x, mu, sigma, lg); }
    inline double pnorm(double x, double mu, double sigma, int lt, int lg)      { return ::Rf_pnorm5(x, mu, sigma, lt, lg); }
    inline double qnorm(double p, double mu, double sigma, int lt, int lg)      { return ::Rf_qnorm5(p, mu, sigma, lt, lg); }
    inline double rnorm(double mu, double sigma)                                { return ::Rf_rnorm(mu, sigma); }
    inline void	pnorm_both(double x, double *cum, double *ccum, int lt, int lg) { return ::Rf_pnorm_both(x, cum, ccum, lt, lg); }

    /* Uniform Distribution */
    inline double dunif(double x, double a, double b, int lg)		{ return ::Rf_dunif(x, a, b, lg); }
    inline double punif(double x, double a, double b, int lt, int lg)   { return ::Rf_punif(x, a, b, lt, lg); }
    inline double qunif(double p, double a, double b, int lt, int lg)   { return ::Rf_qunif(p, a, b, lt, lg); }
    inline double runif(double a, double b)                             { return ::Rf_runif(a, b); }

    /* Gamma Distribution */
    inline double dgamma(double x, double shp, double scl, int lg)	   { return ::Rf_dgamma(x, shp, scl, lg); }
    inline double pgamma(double x, double alp, double scl, int lt, int lg) { return ::Rf_pgamma(x, alp, scl, lt, lg); }
    inline double qgamma(double p, double alp, double scl, int lt, int lg) { return ::Rf_qgamma(p, alp, scl, lt, lg); }
    inline double rgamma(double a, double scl)                             { return ::Rf_rgamma(a, scl); }

    inline double log1pmx(double x)                  { return ::Rf_log1pmx(x); }
    inline double log1pexp(double x)                 { return ::log1pexp(x); }  // <-- ../nmath/plogis.c
    inline double lgamma1p(double a)                 { return ::Rf_lgamma1p(a); }
    inline double logspace_add(double lx, double ly) { return ::Rf_logspace_add(lx, ly); }
    inline double logspace_sub(double lx, double ly) { return ::Rf_logspace_sub(lx, ly); }

    /* Beta Distribution */
    inline double dbeta(double x, double a, double b, int lg)         { return ::Rf_dbeta(x, a, b, lg); }
    inline double pbeta(double x, double p, double q, int lt, int lg) { return ::Rf_pbeta(x, p, q, lt, lg); }
    inline double qbeta(double a, double p, double q, int lt, int lg) { return ::Rf_qbeta(a, p, q, lt, lg); }
    inline double rbeta(double a, double b)                           { return ::Rf_rbeta(a, b); }

    /* Lognormal Distribution */
    inline double dlnorm(double x, double ml, double sl, int lg)	 { return ::Rf_dlnorm(x, ml, sl, lg); }
    inline double plnorm(double x, double ml, double sl, int lt, int lg) { return ::Rf_plnorm(x, ml, sl, lt, lg); }
    inline double qlnorm(double p, double ml, double sl, int lt, int lg) { return ::Rf_qlnorm(p, ml, sl, lt, lg); }
    inline double rlnorm(double ml, double sl)                           { return ::Rf_rlnorm(ml, sl); }

    /* Chi-squared Distribution */
    inline double dchisq(double x, double df, int lg)          { return ::Rf_dchisq(x, df, lg); }
    inline double pchisq(double x, double df, int lt, int lg)  { return ::Rf_pchisq(x, df, lt, lg); }
    inline double qchisq(double p, double df, int lt, int lg)  { return ::Rf_qchisq(p, df, lt, lg); }	
    inline double rchisq(double df)                            { return ::Rf_rchisq(df); }

    /* Non-central Chi-squared Distribution */
    inline double dnchisq(double x, double df, double ncp, int lg)          { return ::Rf_dnchisq(x, df, ncp, lg); }
    inline double pnchisq(double x, double df, double ncp, int lt, int lg)  { return ::Rf_pnchisq(x, df, ncp, lt, lg); }
    inline double qnchisq(double p, double df, double ncp, int lt, int lg)  { return ::Rf_qnchisq(p, df, ncp, lt, lg); }	
    inline double rnchisq(double df, double lb)                             { return ::Rf_rnchisq(df, lb); }

    /* F Distibution */
    inline double df(double x, double df1, double df2, int lg)		{ return ::Rf_df(x, df1, df2, lg); }
    inline double pf(double x, double df1, double df2, int lt, int lg)	{ return ::Rf_pf(x, df1, df2, lt, lg); }
    inline double qf(double p, double df1, double df2, int lt, int lg)	{ return ::Rf_qf(p, df1, df2, lt, lg); }
    inline double rf(double df1, double df2)				{ return ::Rf_rf(df1, df2); }

    /* Student t Distibution */
    inline double dt(double x, double n, int lg)			{ return ::Rf_dt(x, n, lg); }
    inline double pt(double x, double n, int lt, int lg)		{ return ::Rf_pt(x, n, lt, lg); }
    inline double qt(double p, double n, int lt, int lg)		{ return ::Rf_qt(p, n, lt, lg); }
    inline double rt(double n)						{ return ::Rf_rt(n); }

    /* Binomial Distribution */
    inline double dbinom(double x, double n, double p, int lg)	  	{ return ::Rf_dbinom(x, n, p, lg); }
    inline double pbinom(double x, double n, double p, int lt, int lg)  { return ::Rf_pbinom(x, n, p, lt, lg); }
    inline double qbinom(double p, double n, double m, int lt, int lg)  { return ::Rf_qbinom(p, n, m, lt, lg); }
    inline double rbinom(double n, double p)				{ return ::Rf_rbinom(n, p); }

    /* Multnomial Distribution */
    inline void rmultinom(int n, double* prob, int k, int* rn)		{ return ::rmultinom(n, prob, k, rn); }

    /* Cauchy Distribution */
    inline double dcauchy(double x, double lc, double sl, int lg)		{ return ::Rf_dcauchy(x, lc, sl, lg); }
    inline double pcauchy(double x, double lc, double sl, int lt, int lg)	{ return ::Rf_pcauchy(x, lc, sl, lt, lg); }
    inline double qcauchy(double p, double lc, double sl, int lt, int lg)	{ return ::Rf_qcauchy(p, lc, sl, lt, lg); }
    inline double rcauchy(double lc, double sl)					{ return ::Rf_rcauchy(lc, sl); }

    /* Exponential Distribution */
    inline double dexp(double x, double sl, int lg)		{ return ::Rf_dexp(x, sl, lg); }
    inline double pexp(double x, double sl, int lt, int lg)	{ return ::Rf_pexp(x, sl, lt, lg); }
    inline double qexp(double p, double sl, int lt, int lg)	{ return ::Rf_qexp(p, sl, lt, lg); }
    inline double rexp(double sl)				{ return ::Rf_rexp(sl); }

    /* Geometric Distribution */
    inline double dgeom(double x, double p, int lg)		{ return ::Rf_dgeom(x, p, lg); }
    inline double pgeom(double x, double p, int lt, int lg)	{ return ::Rf_pgeom(x, p, lt, lg); }
    inline double qgeom(double p, double pb, int lt, int lg)	{ return ::Rf_qgeom(p, pb, lt, lg); }
    inline double rgeom(double p)				{ return ::Rf_rgeom(p); }

    /* Hypergeometric Distibution */
    inline double dhyper(double x, double r, double b, double n, int lg)		{ return ::Rf_dhyper(x, r, b, n, lg); }
    inline double phyper(double x, double r, double b, double n, int lt, int lg)	{ return ::Rf_phyper(x, r, b, n, lt, lg); }
    inline double qhyper(double p, double r, double b, double n, int lt, int lg)	{ return ::Rf_qhyper(p, r, b, n, lt, lg); }
    inline double rhyper(double r, double b, double n)					{ return ::Rf_rhyper(r, b, n); }

    /* Negative Binomial Distribution */
    inline double dnbinom(double x, double sz, double pb, int lg)		{ return ::Rf_dnbinom(x, sz, pb, lg); }
    inline double pnbinom(double x, double sz, double pb, int lt, int lg)	{ return ::Rf_pnbinom(x, sz, pb, lt, lg); }
    inline double qnbinom(double p, double sz, double pb, int lt, int lg)	{ return ::Rf_qnbinom(p, sz, pb, lt, lg); }
    inline double rnbinom(double sz, double pb)					{ return ::Rf_rnbinom(sz, pb); }

    inline double dnbinom_mu(double x, double sz, double mu, int lg)		{ return ::Rf_dnbinom(x, sz, mu, lg); }
    inline double pnbinom_mu(double x, double sz, double mu, int lt, int lg)	{ return ::Rf_pnbinom(x, sz, mu, lt, lg); }
    inline double qnbinom_mu(double x, double sz, double mu, int lt, int lg)	{ return ::Rf_qnbinom(x, sz, mu, lt, lg); }
    inline double rnbinom_mu(double sz, double mu)				{ return ::Rf_rnbinom(sz, mu); }

    /* Poisson Distribution */
    inline double dpois(double x, double lb, int lg)		{ return ::Rf_dpois(x, lb, lg); }
    inline double ppois(double x, double lb, int lt, int lg)	{ return ::Rf_ppois(x, lb, lt, lg); }
    inline double qpois(double p, double lb, int lt, int lg)	{ return ::Rf_qpois(p, lb, lt, lg); }
    inline double rpois(double mu)				{ return ::Rf_rpois(mu); }

    /* Weibull Distribution */
    inline double dweibull(double x, double sh, double sl, int lg)		{ return ::Rf_dweibull(x, sh, sl, lg); }
    inline double pweibull(double x, double sh, double sl, int lt, int lg)	{ return ::Rf_pweibull(x, sh, sl, lt, lg); }
    inline double qweibull(double p, double sh, double sl, int lt, int lg)	{ return ::Rf_qweibull(p, sh, sl, lt, lg); }
    inline double rweibull(double sh, double sl)				{ return ::Rf_rweibull(sh, sl); }

    /* Logistic Distribution */
    inline double dlogis(double x, double lc, double sl, int lg)		{ return ::Rf_dlogis(x, lc, sl, lg); }
    inline double plogis(double x, double lc, double sl, int lt, int lg)	{ return ::Rf_plogis(x, lc, sl, lt, lg); }
    inline double qlogis(double p, double lc, double sl, int lt, int lg)	{ return ::Rf_qlogis(p, lc, sl, lt, lg); }
    inline double rlogis(double lc, double sl)					{ return ::Rf_rlogis(lc, sl); }

    /* Non-central Beta Distribution */
    inline double dnbeta(double x, double a, double b, double ncp, int lg)		{ return ::Rf_dnbeta(x, a, b, ncp, lg); }
    inline double pnbeta(double x, double a, double b, double ncp, int lt, int lg)	{ return ::Rf_pnbeta(x, a, b, ncp, lt, lg); }
    inline double qnbeta(double p, double a, double b, double ncp, int lt, int lg)	{ return ::Rf_qnbeta(p, a, b, ncp, lt, lg); }
    inline double rnbeta(double a, double b, double np)					{ return ::Rf_rnbeta(a, b, np); }

    /* Non-central F Distribution */
    inline double dnf(double x, double df1, double df2, double ncp, int lg)		{ return ::Rf_dnf(x, df1, df2, ncp, lg); }
    inline double pnf(double x, double df1, double df2, double ncp, int lt, int lg)	{ return ::Rf_pnf(x, df1, df2, ncp, lt, lg); }
    inline double qnf(double p, double df1, double df2, double ncp, int lt, int lg)	{ return ::Rf_qnf(p, df1, df2, ncp, lt, lg); }

    /* Non-central Student t Distribution */
    inline double dnt(double x, double df, double ncp, int lg)		{ return ::Rf_dnt(x, df, ncp, lg); }
    inline double pnt(double x, double df, double ncp, int lt, int lg)	{ return ::Rf_pnt(x, df, ncp, lt, lg); }
    inline double qnt(double p, double df, double ncp, int lt, int lg)	{ return ::Rf_qnt(p, df, ncp, lt, lg); }

    /* Studentized Range Distribution */
    inline double ptukey(double q, double rr, double cc, double df, int lt, int lg)	{ return ::Rf_ptukey(q, rr, cc, df, lt, lg); }
    inline double qtukey(double p, double rr, double cc, double df, int lt, int lg)	{ return ::Rf_qtukey(p, rr, cc, df, lt, lg); }

    /* Wilcoxon Rank Sum Distribution */
    inline double dwilcox(double x, double m, double n, int lg)		{ return ::Rf_dwilcox(x, m, n, lg); }
    inline double pwilcox(double q, double m, double n, int lt, int lg)	{ return ::Rf_pwilcox(q, m, n, lt, lg); }
    inline double qwilcox(double x, double m, double n, int lt, int lg)	{ return ::Rf_qwilcox(x, m, n, lt, lg); }
    inline double rwilcox(double m, double n)				{ return ::Rf_rwilcox(m, n); }

    /* Wilcoxon Signed Rank Distribution */
    inline double dsignrank(double x, double n, int lg)			{ return ::Rf_dsignrank(x, n, lg); }
    inline double psignrank(double x, double n, int lt, int lg)		{ return ::Rf_psignrank(x, n, lt, lg); }
    inline double qsignrank(double x, double n, int lt, int lg)		{ return ::Rf_qsignrank(x, n, lt, lg); }
    inline double rsignrank(double n)					{ return ::Rf_rsignrank(n); }

    /* Gamma and Related Functions */
    inline double gammafn(double x)			{ return ::Rf_gammafn(x); }
    inline double lgammafn(double x)			{ return ::Rf_lgammafn(x); }
    inline double lgammafn_sign(double x, int *sgn)	{ return ::Rf_lgammafn_sign(x, sgn); }
    inline void   dpsifn(double x, int n, int kode, int m, double *ans, int *nz, int *ierr)	{ return ::Rf_dpsifn(x, n, kode, m, ans, nz, ierr); }
    inline double psigamma(double x, double deriv)	{ return ::Rf_psigamma(x, deriv); }
    inline double digamma(double x)	{ return ::Rf_digamma(x); }
    inline double trigamma(double x)	{ return ::Rf_trigamma(x); }
    inline double tetragamma(double x)	{ return ::Rf_tetragamma(x); }
    inline double pentagamma(double x)	{ return ::Rf_pentagamma(x); }
    
    inline double beta(double a, double b)	{ return ::Rf_beta(a, b); }
    inline double lbeta(double a, double b)	{ return ::Rf_lbeta(a, b); }

    inline double choose(double n, double k)	{ return ::Rf_choose(n, k); }
    inline double lchoose(double n, double k)	{ return ::Rf_lchoose(n, k); }

    /* Bessel Functions */
    inline double bessel_i(double x, double al, double ex)	{ return ::Rf_bessel_i(x, al, ex); }
    inline double bessel_j(double x, double al)			{ return ::Rf_bessel_j(x, al); }
    inline double bessel_k(double x, double al, double ex)	{ return ::Rf_bessel_k(x, al, ex); }
    inline double bessel_y(double x, double al)			{ return ::Rf_bessel_y(x, al); }
    inline double bessel_i_ex(double x, double al, double ex, double *bi)	{ return ::Rf_bessel_i_ex(x, al, ex, bi); }
    inline double bessel_j_ex(double x, double al, double *bj)			{ return ::Rf_bessel_j_ex(x, al, bj); }
    inline double bessel_k_ex(double x, double al, double ex, double *bk)	{ return ::Rf_bessel_k_ex(x, al, ex, bk); }
    inline double bessel_y_ex(double x, double al, double *by)			{ return ::Rf_bessel_y_ex(x, al, by); }

    /* General Support Functions */
#ifndef HAVE_HYPOT
    inline double hypot(double a, double b)	{ return ::Rf_hypot(a, b); }
#endif
    inline double pythag(double a, double b)	{ return ::Rf_pythag(a, b); }
#ifndef HAVE_EXPM1
    inline double expm1(double x); /* = exp(x)-1 {care for small x} */	{ return ::Rf_expm1(x); }
#endif
#ifndef HAVE_LOG1P
    inline double log1p(double x); /* = log(1+x) {care for small x} */ { return ::Rf_log1p(x); }
#endif
    inline int imax2(int x, int y)		{ return ::Rf_imax2(x, y); }
    inline int imin2(int x, int y)		{ return ::Rf_imin2(x, y); }
    inline double fmax2(double x, double y)	{ return ::Rf_fmax2(x, y); }
    inline double fmin2(double x, double y)	{ return ::Rf_fmin2(x, y); }
    inline double sign(double x)		{ return ::Rf_sign(x); }
    inline double fprec(double x, double dg)	{ return ::Rf_fprec(x, dg); }
    inline double fround(double x, double dg)	{ return ::Rf_fround(x, dg); }
    inline double fsign(double x, double y)	{ return ::Rf_fsign(x, y); }
    inline double ftrunc(double x)		{ return ::Rf_ftrunc(x); }
    
}

#endif
=======
/* -*- C -*-
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2012  The R Core Team
 *  Copyright (C) 2004       The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *

 * Rmath.h  should contain ALL headers from R's C code in `src/nmath'
   -------  such that ``the Math library'' can be used by simply

   ``#include <Rmath.h> ''

   and nothing else.

   It is part of the API and supports 'standalone Rmath'.

*/
#ifndef RMATH_H
#define RMATH_H

#define MATHLIB_STANDALONE

/* Note that on some systems we need to include math.h before the
   defines below, to avoid redefining ftrunc */
#ifndef NO_C_HEADERS
# include <math.h>
#endif

/*-- Mathlib as part of R --  define this for standalone : */
/* #undef MATHLIB_STANDALONE */

#define R_VERSION_STRING "3.0.2"

#ifndef HAVE_EXPM1
# define HAVE_EXPM1 1
#endif

#ifndef HAVE_HYPOT
# define HAVE_HYPOT 1
#endif

#ifndef HAVE_LOG1P
# define HAVE_LOG1P 1
#endif

#ifndef HAVE_WORKING_LOG1P
# define HAVE_WORKING_LOG1P 1
#endif

#if defined(HAVE_LOG1P) && !defined(HAVE_WORKING_LOG1P)
/* remap to avoid problems with getting the right entry point */
double  Rlog1p(double);
#define log1p Rlog1p
#endif


	/* Undo SGI Madness */

#ifdef ftrunc
# undef ftrunc
#endif
#ifdef qexp
# undef qexp
#endif
#ifdef qgamma
# undef qgamma
#endif


/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2)) */
#endif


#ifdef MATHLIB_STANDALONE
# ifndef R_EXT_BOOLEAN_H_
/* "copy-paste" R_ext/Boolean.h if not already included: */
 #define R_EXT_BOOLEAN_H_
 #undef FALSE
 #undef TRUE
 typedef enum { FALSE = 0, TRUE } Rboolean;
# endif
#else
# include <R_ext/Boolean.h>
#endif


#ifndef MATHLIB_STANDALONE
#define bessel_i	Rf_bessel_i
#define bessel_j	Rf_bessel_j
#define bessel_k	Rf_bessel_k
#define bessel_y	Rf_bessel_y
#define bessel_i_ex	Rf_bessel_i_ex
#define bessel_j_ex	Rf_bessel_j_ex
#define bessel_k_ex	Rf_bessel_k_ex
#define bessel_y_ex	Rf_bessel_y_ex
#define beta		Rf_beta
#define choose		Rf_choose
#define dbeta		Rf_dbeta
#define dbinom		Rf_dbinom
#define dcauchy		Rf_dcauchy
#define dchisq		Rf_dchisq
#define dexp		Rf_dexp
#define df		Rf_df
#define dgamma		Rf_dgamma
#define dgeom		Rf_dgeom
#define dhyper		Rf_dhyper
#define digamma		Rf_digamma
#define dlnorm		Rf_dlnorm
#define dlogis		Rf_dlogis
#define dnbeta		Rf_dnbeta
#define dnbinom		Rf_dnbinom
#define dnbinom_mu	Rf_dnbinom_mu
#define dnchisq		Rf_dnchisq
#define dnf		Rf_dnf
#define dnorm4		Rf_dnorm4
#define dnt		Rf_dnt
#define dpois		Rf_dpois
#define dpsifn		Rf_dpsifn
#define dsignrank	Rf_dsignrank
#define dt		Rf_dt
#define dtukey		Rf_dtukey
#define dunif		Rf_dunif
#define dweibull	Rf_dweibull
#define dwilcox		Rf_dwilcox
#define fmax2		Rf_fmax2
#define fmin2		Rf_fmin2
#define fprec		Rf_fprec
#define fround		Rf_fround
#define ftrunc		Rf_ftrunc
#define fsign		Rf_fsign
#define gammafn		Rf_gammafn
#define imax2		Rf_imax2
#define imin2		Rf_imin2
#define lbeta		Rf_lbeta
#define lchoose		Rf_lchoose
#define lgammafn	Rf_lgammafn
#define lgammafn_sign	Rf_lgammafn_sign
#define lgamma1p	Rf_lgamma1p
#define log1pmx		Rf_log1pmx
#define logspace_add	Rf_logspace_add
#define logspace_sub	Rf_logspace_sub
#define pbeta		Rf_pbeta
#define pbeta_raw	Rf_pbeta_raw
#define pbinom		Rf_pbinom
#define pcauchy		Rf_pcauchy
#define pchisq		Rf_pchisq
#define pentagamma	Rf_pentagamma
#define pexp		Rf_pexp
#define pf		Rf_pf
#define pgamma		Rf_pgamma
#define pgeom		Rf_pgeom
#define phyper		Rf_phyper
#define plnorm		Rf_plnorm
#define plogis		Rf_plogis
#define pnbeta		Rf_pnbeta
#define pnbinom		Rf_pnbinom
#define pnbinom_mu     	Rf_pnbinom_mu
#define pnchisq		Rf_pnchisq
#define pnf		Rf_pnf
#define pnorm5		Rf_pnorm5
#define pnorm_both	Rf_pnorm_both
#define pnt		Rf_pnt
#define ppois		Rf_ppois
#define psignrank	Rf_psignrank
#define psigamma	Rf_psigamma
#define pt		Rf_pt
#define ptukey		Rf_ptukey
#define punif		Rf_punif
#define pythag		Rf_pythag
#define pweibull	Rf_pweibull
#define pwilcox		Rf_pwilcox
#define qbeta		Rf_qbeta
#define qbinom		Rf_qbinom
#define qcauchy		Rf_qcauchy
#define qchisq		Rf_qchisq
#define qchisq_appr	Rf_qchisq_appr
#define qexp		Rf_qexp
#define qf		Rf_qf
#define qgamma		Rf_qgamma
#define qgeom		Rf_qgeom
#define qhyper		Rf_qhyper
#define qlnorm		Rf_qlnorm
#define qlogis		Rf_qlogis
#define qnbeta		Rf_qnbeta
#define qnbinom		Rf_qnbinom
#define qnbinom_mu     	Rf_qnbinom_mu
#define qnchisq		Rf_qnchisq
#define qnf		Rf_qnf
#define qnorm5		Rf_qnorm5
#define qnt		Rf_qnt
#define qpois		Rf_qpois
#define qsignrank	Rf_qsignrank
#define qt		Rf_qt
#define qtukey		Rf_qtukey
#define qunif		Rf_qunif
#define qweibull	Rf_qweibull
#define qwilcox		Rf_qwilcox
#define rbeta		Rf_rbeta
#define rbinom		Rf_rbinom
#define rcauchy		Rf_rcauchy
#define rchisq		Rf_rchisq
#define rexp		Rf_rexp
#define rf		Rf_rf
#define rgamma		Rf_rgamma
#define rgeom		Rf_rgeom
#define rhyper		Rf_rhyper
#define rlnorm		Rf_rlnorm
#define rlogis		Rf_rlogis
#define rnbeta		Rf_rnbeta
#define rnbinom		Rf_rnbinom
#define rnchisq		Rf_rnchisq
#define rnf		Rf_rnf
#define rnorm		Rf_rnorm
#define rnt		Rf_rnt
#define rpois		Rf_rpois
#define rsignrank	Rf_rsignrank
#define rt		Rf_rt
#define rtukey		Rf_rtukey
#define runif		Rf_runif
#define rweibull	Rf_rweibull
#define rwilcox		Rf_rwilcox
#define sign		Rf_sign
#define tetragamma	Rf_tetragamma
#define trigamma	Rf_trigamma
#endif

#define	rround	fround
#define	prec	fprec
#undef trunc
#define	trunc	ftrunc

#ifdef  __cplusplus
extern "C" {
#endif
	/* R's versions with !R_FINITE checks */

double R_pow(double x, double y);
double R_pow_di(double, int);

	/* Random Number Generators */

double	norm_rand(void);
double	unif_rand(void);
double	exp_rand(void);
#ifdef MATHLIB_STANDALONE
void	set_seed(unsigned int, unsigned int);
void	get_seed(unsigned int *, unsigned int *);
#endif

	/* Normal Distribution */

#define pnorm pnorm5
#define qnorm qnorm5
#define dnorm dnorm4

double	dnorm(double, double, double, int);
double	pnorm(double, double, double, int, int);
double	qnorm(double, double, double, int, int);
double	rnorm(double, double);
void	pnorm_both(double, double *, double *, int, int);/* both tails */

	/* Uniform Distribution */

double	dunif(double, double, double, int);
double	punif(double, double, double, int, int);
double	qunif(double, double, double, int, int);
double	runif(double, double);

	/* Gamma Distribution */

double	dgamma(double, double, double, int);
double	pgamma(double, double, double, int, int);
double	qgamma(double, double, double, int, int);
double	rgamma(double, double);

double  log1pmx(double);
double  log1pexp(double); // <-- ../nmath/plogis.c
double  lgamma1p(double);
double  logspace_add(double, double);
double  logspace_sub(double, double);

	/* Beta Distribution */

double	dbeta(double, double, double, int);
double	pbeta(double, double, double, int, int);
double	qbeta(double, double, double, int, int);
double	rbeta(double, double);

	/* Lognormal Distribution */

double	dlnorm(double, double, double, int);
double	plnorm(double, double, double, int, int);
double	qlnorm(double, double, double, int, int);
double	rlnorm(double, double);

	/* Chi-squared Distribution */

double	dchisq(double, double, int);
double	pchisq(double, double, int, int);
double	qchisq(double, double, int, int);
double	rchisq(double);

	/* Non-central Chi-squared Distribution */

double	dnchisq(double, double, double, int);
double	pnchisq(double, double, double, int, int);
double	qnchisq(double, double, double, int, int);
double	rnchisq(double, double);

	/* F Distibution */

double	df(double, double, double, int);
double	pf(double, double, double, int, int);
double	qf(double, double, double, int, int);
double	rf(double, double);

	/* Student t Distibution */

double	dt(double, double, int);
double	pt(double, double, int, int);
double	qt(double, double, int, int);
double	rt(double);

	/* Binomial Distribution */

double	dbinom(double, double, double, int);
double	pbinom(double, double, double, int, int);
double	qbinom(double, double, double, int, int);
double	rbinom(double, double);

	/* Multnomial Distribution */

void	rmultinom(int, double*, int, int*);

	/* Cauchy Distribution */

double	dcauchy(double, double, double, int);
double	pcauchy(double, double, double, int, int);
double	qcauchy(double, double, double, int, int);
double	rcauchy(double, double);

	/* Exponential Distribution */

double	dexp(double, double, int);
double	pexp(double, double, int, int);
double	qexp(double, double, int, int);
double	rexp(double);

	/* Geometric Distribution */

double	dgeom(double, double, int);
double	pgeom(double, double, int, int);
double	qgeom(double, double, int, int);
double	rgeom(double);

	/* Hypergeometric Distibution */

double	dhyper(double, double, double, double, int);
double	phyper(double, double, double, double, int, int);
double	qhyper(double, double, double, double, int, int);
double	rhyper(double, double, double);

	/* Negative Binomial Distribution */

double	dnbinom(double, double, double, int);
double	pnbinom(double, double, double, int, int);
double	qnbinom(double, double, double, int, int);
double	rnbinom(double, double);

double	dnbinom_mu(double, double, double, int);
double	pnbinom_mu(double, double, double, int, int);
double	qnbinom_mu(double, double, double, int, int);
double	rnbinom_mu(double, double);

	/* Poisson Distribution */

double	dpois(double, double, int);
double	ppois(double, double, int, int);
double	qpois(double, double, int, int);
double	rpois(double);

	/* Weibull Distribution */

double	dweibull(double, double, double, int);
double	pweibull(double, double, double, int, int);
double	qweibull(double, double, double, int, int);
double	rweibull(double, double);

	/* Logistic Distribution */

double	dlogis(double, double, double, int);
double	plogis(double, double, double, int, int);
double	qlogis(double, double, double, int, int);
double	rlogis(double, double);

	/* Non-central Beta Distribution */

double	dnbeta(double, double, double, double, int);
double	pnbeta(double, double, double, double, int, int);
double	qnbeta(double, double, double, double, int, int);
double	rnbeta(double, double, double);

	/* Non-central F Distribution */

double  dnf(double, double, double, double, int);
double	pnf(double, double, double, double, int, int);
double	qnf(double, double, double, double, int, int);

	/* Non-central Student t Distribution */

double	dnt(double, double, double, int);
double	pnt(double, double, double, int, int);
double	qnt(double, double, double, int, int);

	/* Studentized Range Distribution */

double	ptukey(double, double, double, double, int, int);
double	qtukey(double, double, double, double, int, int);

	/* Wilcoxon Rank Sum Distribution */

double dwilcox(double, double, double, int);
double pwilcox(double, double, double, int, int);
double qwilcox(double, double, double, int, int);
double rwilcox(double, double);

	/* Wilcoxon Signed Rank Distribution */

double dsignrank(double, double, int);
double psignrank(double, double, int, int);
double qsignrank(double, double, int, int);
double rsignrank(double);

	/* Gamma and Related Functions */
double	gammafn(double);
double	lgammafn(double);
double	lgammafn_sign(double, int*);
void    dpsifn(double, int, int, int, double*, int*, int*);
double	psigamma(double, double);
double	digamma(double);
double	trigamma(double);
double	tetragamma(double);
double	pentagamma(double);

double	beta(double, double);
double	lbeta(double, double);

double	choose(double, double);
double	lchoose(double, double);

	/* Bessel Functions */

double	bessel_i(double, double, double);
double	bessel_j(double, double);
double	bessel_k(double, double, double);
double	bessel_y(double, double);
double	bessel_i_ex(double, double, double, double *);
double	bessel_j_ex(double, double, double *);
double	bessel_k_ex(double, double, double, double *);
double	bessel_y_ex(double, double, double *);


	/* General Support Functions */

#ifndef HAVE_HYPOT
double 	hypot(double, double);
#endif
double 	pythag(double, double);
#ifndef HAVE_EXPM1
double  expm1(double); /* = exp(x)-1 {care for small x} */
#endif
#ifndef HAVE_LOG1P
double  log1p(double); /* = log(1+x) {care for small x} */
#endif
int	imax2(int, int);
int	imin2(int, int);
double	fmax2(double, double);
double	fmin2(double, double);
double	sign(double);
double	fprec(double, double);
double	fround(double, double);
double	fsign(double, double);
double	ftrunc(double);

double  log1pmx(double); /* Accurate log(1+x) - x, {care for small x} */
double  lgamma1p(double);/* accurate log(gamma(x+1)), small x (0 < x < 0.5) */

/* Compute the log of a sum or difference from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 * or  log (exp (logx) - exp (logy))
 *
 * without causing overflows or throwing away too much accuracy:
 */
double  logspace_add(double logx, double logy);
double  logspace_sub(double logx, double logy);




/* ----------------- Private part of the header file ------------------- */

	/* old-R Compatibility */

#ifdef OLD_RMATH_COMPAT
# define snorm	norm_rand
# define sunif	unif_rand
# define sexp	exp_rand
#endif

#if defined(MATHLIB_STANDALONE) && !defined(MATHLIB_PRIVATE_H)
/* second is defined by nmath.h */

/* If isnan is a macro, as C99 specifies, the C++
   math header will undefine it. This happens on OS X */
# ifdef __cplusplus
  int R_isnancpp(double); /* in mlutils.c */
#  define ISNAN(x)     R_isnancpp(x)
# else
#  define ISNAN(x)     (isnan(x)!=0)
# endif

# define R_FINITE(x)    R_finite(x)
int R_finite(double);

# ifdef WIN32  /* not Win32 as no config information */
#  ifdef RMATH_DLL
#   define R_EXTERN extern __declspec(dllimport)
#  else
#   define R_EXTERN extern
#  endif
R_EXTERN double NA_REAL;
R_EXTERN double R_PosInf;
R_EXTERN double R_NegInf;
R_EXTERN int N01_kind;
#  undef R_EXTERN
#else
extern int N01_kind;
# endif

#endif /* MATHLIB_STANDALONE */

#ifdef  __cplusplus
}
#endif

#endif /* RMATH_H */
>>>>>>> b4b15d3... pulled in R source code
