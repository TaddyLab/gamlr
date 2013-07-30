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
   gamma=0,
   nlambda=100, lambda.start=Inf,  
   lambda.min.ratio=0.01, weight=NULL, 
   standardize=TRUE, 
   thresh=1e-6, maxit=1e5,
   verb=FALSE)

\method{plot}{gamlr}(x, against=c("pen","dev"), 
    col=rgb(0,0,.5,.75), select=TRUE, df=TRUE, ...)
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
  \item{gamma}{ Prior scale on the 
      L1 cost parameter \eqn{\lambda}.  
      The default of zero corresponds to the lasso,
      and higher values correspond to a more concave penalty.  
      See details. }
  \item{nlambda}{ Number of regularization path segments. }
  \item{lambda.start}{ Initial penalty value.  Default of \code{Inf}
  implies the suppremum lambda that returns all zero
  coefficients.  This is the largest absolute coefficient gradient at the null model. }
  \item{lambda.min.ratio}{ The smallest penalty weight 
    (expected L1 cost) as a ratio of the path start value.  
    Our default is always 0.01; note that this differs from \code{glmnet}
    whose default depends upon the dimension of \code{x}. }
  \item{weight}{ Parameter weights.  
      The default is 
      all 1.0, and a \code{weight} entry of 0 always
      means that the corresponding coefficient is unpenalized.
      Other values have penalty-dependent behavior: for \code{gamma>0} 
      \code{weight} is an up-scaling of 
      \code{gamma} that is specific to each covariate (see details), and 
      for the lasso (\code{gamma=0}), 
      these are multipliers on the L1 penalty 
      (e.g., for the adaptive lasso).   }
  \item{standardize}{ Whether to standardize 
    the coefficients to have standard deviation of one.  
    This is equivalent to multiplying the L1 penalty 
    by each coefficient standard deviation. }
  \item{verb}{ Whether to print some output for each path segment. }
  \item{thresh}{ Optimization convergence tolerance for each 
    inner coordinate-descent loop.  This is measured against the 
    largest coordinate change after a full parameter-set update, 
    relative to the null model negative log likelihood. }
  \item{maxit}{ Max iterations for a single segment
         coordinate descent routine. }
   \item{object}{ A gamlr object.}
  \item{against}{ Whether to plot paths 
  against log penalty or deviance.}
\item{col}{ A single plot color (default is 3/4 opaque navy), 
  or vector of length \code{ncol(x)} colors for each coefficient
    regularization path. }
\item{select}{ Index of the path segment 
  for which you want coefficients or prediction, or in \code{plot} 
  whether to add lines marking AIC and BIC selected models.}
\item{df}{ Whether to add to the plot degrees of freedom along the top axis.}
\item{newdata}{ New \code{x} data for prediction.}
\item{type}{ Either "link" for the linear equation, 
or "response" for predictions transformed 
to the same domain as \code{y}.}
\item{...}{ Extra arguments to each method. }
}
\details{ Finds the posterior modes along a regularization path
		   of \emph{prior expected L1 penalties} using coordinate descent.

  Each path segment minimizes the log-cost penalized objective 
  -\eqn{\phi}logLHD\eqn{(\beta_1 ... \beta_p) + \sum s log[r + |\beta_j|]},
  where \eqn{\phi} is the exponential
  family dispersion parameter (\eqn{\sigma^2} for \code{family="gaussian"},
  one otherwise).  

  This log penalty arises through joint estimation for coefficients and their
  unique L1 penalties \eqn{\omega_1 ... \omega_p} under a \eqn{Gamma(s/\phi, r/\phi)} prior with expectation \eqn{E[\omega] = s/r}
  and variance \eqn{var[\omega] = \phi s/r^2}.
  Reparameterizing \eqn{E[\omega] = n\lambda} and
  \eqn{var[\omega] = n \phi\gamma}, our penalty path is a grid
  in \code{lambda} = \eqn{\lambda} under scale \code{gamma} = \eqn{\gamma}.
  For \code{gamma=0} the L1 penalties are all fixed at 
  \eqn{\omega = n\lambda} (i.e. equivalent to \code{glmnet}'s lasso).

    For completeness, setting \code{gamma=Inf} yields a stepwise
      subset selection algorithm: each step adds (unpenalized)
      the coefficient with highest absolute gradient.  The reported
      \code{lambda} is just \code{exp(-df)}, where \code{df} is 
      the number of nonzero coefficients 
      (i.e. not the true degrees of freedom, which is essentially \code{ncol(x)}).

      \code{plot.gamlr} can be used to graph the results: it 
      shows the regularization paths for \emph{beta} and marks
      the minimum AIC and BIC models.
}
\value{
  \item{lambda}{The path of fitted \emph{prior expected} L1 penalties.}
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
\references{Taddy (2013), The Gamma Lasso.}

\examples{

### a low-D test (highly multi-collinear)

n <- 1000
p <- 3
xvar <- matrix(0.9, nrow=p,ncol=p)
diag(xvar) <- 1
x <- matrix(rnorm(p*n), nrow=n)\%*\%chol(xvar)
y <- 4 + 3*x[,1] + -1*x[,2] + rnorm(n)

fitlasso <- gamlr(x, y, gamma=0) # lasso
fitgl <- gamlr(x, y, gamma=1/6) # log penalty
fitglbv <- gamlr(x, y, gamma=6) # big variance log penalty

par(mfrow=c(1,3))
ylim = range(c(fitglbv$beta@x))
plot(fitlasso, ylim=ylim)
plot(fitgl, ylim=ylim, col="maroon")
plot(fitglbv, ylim=ylim, col="darkorange")

 }
\seealso{cv.gamlr, hockey}
