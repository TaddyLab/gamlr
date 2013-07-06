\name{gamlr}
\alias{gamlr}
\alias{print.gamlr}
\alias{summary.gamlr}
\alias{predict.gamlr}
\alias{plot.gamlr}
\alias{coef.gamlr}
\alias{logLik.gamlr}
\alias{family.gamlr}
\title{Gamma-Lasso regression}
\description{ Log penalized regression estimation. }
\usage{
gamlr(x, y, 
   family=c("gaussian","binomial","poisson"),
   varpen=0,
   npen=100, pen.min.ratio=0.01, 
   pen.start=Inf,  weight=NULL, 
   standardize=TRUE, verb=FALSE,
   thresh=1e-6, maxit=1e5, qn=FALSE, ...)
plot.gamlr(fit, against=c("pen","dev"), col="navy", ...)
coef.gamlr(object, select=which.min(BIC(object)),
            type = c("link", "response"), ...)
predict.gamlr(object, newdata,
            select=which.min(BIC(object)),
            type = c("link", "response"))
summary.gamlr(object, ...)
print.gamlr(object, ...)
logLik.gamlr(object, ...)
family.gamlr(object, ...)
}
\arguments{
  \item{x}{ A dense \code{matrix} 
      or sparse \code{Matrix} of covariates,
      with \code{ncol(x)} variables and 
      \code{nrow(x)==length(y)} observations.
      This should not include the intercept.}
  \item{y}{A vector of response values. 
      Note that there is no argument checking; 
      if you choose the wrong 
      \code{family} is could cause problems. }
  \item{family}{ Response model type; 
  either "gaussian", "poisson", or "binomial".  
  Note that for "binomial", \code{y} is in \eqn{[0,1]}. }
  \item{varpen}{ Prior variance of the 
      L1 cost parameter \eqn{\lambda}.  
      The default of zero corresponds to the lasso,
      and higher values correspond to a MORE concave penalty.  
      See details. }
  \item{npen}{ Number of regularization path segments. }
  \item{pen.start}{ Initial penalty value.  Default of \code{Inf}
  implies the suppremum penalty that returns all zero
  coefficients.  This is the largest absolute coefficient gradient at the null model. }
  \item{pen.min.ratio}{ The smallest penalty weight 
    (expected L1 cost) as a ratio of the path start value. }
  \item{weight}{ Parameter weights; for \code{varpen>0} 
      this is an up-scaling of 
      \code{varpen} that is specific to each covariate (see details). 
      For the lasso (\code{varpen=0}), 
      these are multipliers on the L1 penalty.  The default is 
      all \code{1.0}, and a \code{weight} entry of \code{0} always
      means that the corresponding coefficient is unpenalized. }
  \item{standardize}{ Whether to standardize 
    the coefficients to have standard deviation of one.  
    This is equivalent to multiplying the L1 penalty 
    by each coefficient standard deviation. }
  \item{verb}{ Whether to print some output for each path segment. }
  \item{thresh}{ Optimization convergence tolerance for each 
    inner coordinate-descent loop.  This is measured as the 
    drop in negative log posterior after a full parameter set update, 
    relative to the null model negative log likelihood. }
  \item{maxit}{ Max iterations for a segment
         coordinate descent routine. }
  \item{qn}{ Whether to us quasi-newton acceleration. 
  This is useful for high dimensions with high multicollinearity, 
  but the overhead is too costly otherwise. }
}
\details{ Finds the posterior modes along a regularization path
		   for regression parameters using coordinate descent.

  \emph{Coefficient penalization} is based upon the precision parameters
  \eqn{\lambda} of independent L1 penalties. Via the
       \code{varpen} argument, this precision is either fixed (\code{varpen=0}, i.e. the lasso) or it is assigned a
       \eqn{Gamma(s, r)} prior, with expectation \code{pen = s/r} 
       and variance \code{varpen = s/r^2}. This is equivalent to a log penalty \eqn{s*log[r + |b|]}.
}
\value{
  The returned \code{gamlr} object list including some of the
  input variables (possibly cleaned or updated), as well as
  \item{penalty}{The path of fitted expected L1 penalties.}
  \item{nobs}{ The number of observations.}
  \item{alpha}{Intercepts.}
  \item{beta}{Regression coefficients.}
  \item{df}{Approximate degrees of freedom.}
  \item{deviance}{Fitted deviance: 
  -2( logLHD.fitted - logLHD.saturated). }
  \item{iterations}{Total number of passes through updates for the
    full parameter set. }
}
\references{
  Taddy 2013.
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

fitlasso <- gamlr(x, y, penvar=0) # lasso
fitgl1 <- gamlr(x, y, penvar=1) # log penalty var 1
fitgl10 <- gamlr(x, y, penvar=10) # log penalty var 10

par(mfrow=c(1,3))
ylim = range(c(fitgl10$beta$x))
plot(fitlasso, ylim=ylim)
plot(fitgl1, ylim=ylim, col="maroon")
plot(fitgl10, ylim=ylim, col="darkorange")

 }
\seealso{cv.gamlr,hockey}
