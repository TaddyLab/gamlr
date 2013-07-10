\name{gamlr}
\alias{gamlr}
\alias{predict.gamlr}
\alias{plot.gamlr}
\alias{coef.gamlr}
\alias{logLik.gamlr}
\title{Gamma-Lasso regression}
\description{ Log penalized regression estimation. }
\usage{
gamlr(x, y, 
   family=c("gaussian","binomial","poisson"),
   varpen=0,
   npen=100, pen.start=Inf,  
   pen.min.ratio=0.01, weight=NULL, 
   standardize=TRUE, verb=FALSE,
   thresh=1e-6, maxit=1e5, qn=FALSE)

\method{plot}{gamlr}(x, against=c("pen","dev"), col="navy", ...)
\method{coef}{gamlr}(object, select=which.min(BIC(object)), ...)
\method{predict}{gamlr}(object, newdata,
            select=which.min(BIC(object)),
            type = c("link", "response"), ...)
\method{logLik}{gamlr}(object, ...)
}
\arguments{
  \item{x}{ A dense \code{matrix} 
      or sparse \code{Matrix} of covariates,
      with \code{ncol(x)} variables and 
      \code{nrow(x)==length(y)} observations.
      This should not include the intercept.}
  \item{y}{A vector of response values. 
      There is almost no argument checking, 
      so be careful to match \code{y} with the appropriate \code{family}. 
      }
  \item{family}{ Response model type; 
  either "gaussian", "poisson", or "binomial".  
  Note that for "binomial", \code{y} is in \eqn{[0,1]}. }
  \item{varpen}{ Prior variance of the 
      L1 cost parameter \eqn{\lambda}.  
      The default of zero corresponds to the lasso,
      and higher values correspond to a more concave penalty.  
      See details. }
  \item{npen}{ Number of regularization path segments. }
  \item{pen.start}{ Initial penalty value.  Default of \code{Inf}
  implies the suppremum penalty that returns all zero
  coefficients.  This is the largest absolute coefficient gradient at the null model. }
  \item{pen.min.ratio}{ The smallest penalty weight 
    (expected L1 cost) as a ratio of the path start value.  
    Our default is always 0.01; note that this differs from \code{glmnet}
    whose default depends upon the dimension of \code{x}. }
  \item{weight}{ Parameter weights; for \code{varpen>0} 
      this is an up-scaling of 
      \code{varpen} that is specific to each covariate (see details). 
      For the lasso (\code{varpen=0}), 
      these are multipliers on the L1 penalty.  The default is 
      all 1.0, and a \code{weight} entry of 0 always
      means that the corresponding coefficient is unpenalized. }
  \item{standardize}{ Whether to standardize 
    the coefficients to have standard deviation of one.  
    This is equivalent to multiplying the L1 penalty 
    by each coefficient standard deviation. }
  \item{verb}{ Whether to print some output for each path segment. }
  \item{thresh}{ Optimization convergence tolerance for each 
    inner coordinate-descent loop.  This is measured against the 
    drop in negative log posterior after a full parameter set update, 
    relative to the null model negative log likelihood. }
  \item{maxit}{ Max iterations for a single segment
         coordinate descent routine. }
  \item{qn}{ Whether to us quasi-newton acceleration. 
  This is useful for high dimensions with high multicollinearity, 
  but the overhead is too costly otherwise. }
  \item{object}{ A gamlr object.}
  \item{against}{ Whether to plot paths 
  against log penalty or deviance.}
\item{col}{ The single, or vector of \code{ncol(x)} colors for each
    regularization path. }
\item{newdata}{ New \code{x} data for prediction.}
\item{select}{ Index of the path segment 
  for which you want coefficients or prediction. }
\item{type}{ Either "link" for the linear equation, 
or "response" for predictions transformed 
to the same domain as \code{y}.}
\item{...}{ Extra arguments to each method. }
}
\details{ Finds the posterior modes along a regularization path
		   of \emph{prior expected L1 penalties} using coordinate descent.

  Given L1 penalty \eqn{\lambda}, each path segment minimizes
 -logLHD\eqn{(\beta_1 ... \beta_p) 
  + (\lambda/\phi) \sum |\beta_j|}, with \eqn{\phi} the exponential
  family dispersion parameter (\eqn{\sigma^2} for \code{family="gaussian"},
  one otherwise).
  Via the \code{varpen} argument, 
  \eqn{\lambda} is either fixed (\code{varpen=0}, i.e. the lasso) 
  or it is assigned a \eqn{Gamma(s, r)} prior, 
  with expectation \code{pen = s/r}
  and variance \code{varpen = s/r^2}. 
  The latter \code{varpen>0} case is equivalent to a log penalty \eqn{s*log[r + |\beta|]}, and the path for \code{pen} is then in terms of \eqn{E[\lambda]}.

  Note that the \code{penalty} here, \eqn{E[\lambda]}, is not scaled by 
  the number of observations as is done in \code{glmnet}.  Thus to get 
  values comparable to \code{glmnet}'s \code{lambda} you divide 
  \code{fit$penalty/fit$nobs} from a \code{gamlr} object. 

  Finally, for completeness, setting \code{varpen=Inf} yields a stepwise
      subset selection algorithm: each step adds (unpenalized)
      the coefficient with highest absolute gradient.  The reported
      \code{penalty} is just \code{exp(-df)}, where \code{df} is 
      the number of nonzero coefficients.
}
\value{
  \item{penalty}{The path of fitted \emph{prior expected} L1 penalties.}
  \item{nobs}{ The number of observations.}
  \item{alpha}{Intercepts.}
  \item{beta}{Regression coefficients.}
  \item{df}{Approximate degrees of freedom.}
  \item{deviance}{Fitted deviance: 
  -2( logLHD.fitted - logLHD.saturated). }
  \item{iterations}{Total number of coordinate descent cycles. }
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{

### a low-D test (highly multi-collinear)

n <- 1000
p <- 3
xvar <- matrix(0.9, nrow=p,ncol=p)
diag(xvar) <- 1
x <- matrix(rnorm(p*n), nrow=n)\%*\%chol(xvar)
y <- 4 + 3*x[,1] + -1*x[,2] + rnorm(n)

fitlasso <- gamlr(x, y, varpen=0) # lasso
fitgl <- gamlr(x, y, varpen=n/5) # log penalty
fitglbv <- gamlr(x, y, varpen=n) # big variance log penalty

par(mfrow=c(1,3))
ylim = range(c(fitglbv$beta@x))
plot(fitlasso, ylim=ylim)
plot(fitgl, ylim=ylim, col="maroon")
plot(fitglbv, ylim=ylim, col="darkorange")

 }
\seealso{cv.gamlr, hockey}
