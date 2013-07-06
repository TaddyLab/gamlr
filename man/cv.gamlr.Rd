\name{cv.gamlr}
\alias{cv.gamlr}
\alias{plot.cv.gamlr}
\alias{coef.cv.gamlr}
\alias{predict.cv.gamlr}
\alias{summary.cv.gamlr}
\alias{print.cv.gamlr}
\alias{family.cv.gamlr}
\title{ Cross Validation for gamlr}
\description{Cross validation for log penalty selection. }
\usage{
cv.gamlr(x, y, nfold=5, verb=FALSE, ...)
plot.cv.gamlr(x, ...)
coef.cv.gamlr(object, seg=c("1se","min"), ...)
predict.cv.gamlr(object, newdata, seg=c("1se","min"), ...)
summary.cv.gamlr(object, ...)
print.cv.gamlr(x, ...)
}
\arguments{
  \item{x}{ A dense \code{matrix} 
      or sparse \code{Matrix} of covariates,
      with \code{ncol(x)} variables and 
      \code{length(y)} observations.
      This should not include the intercept.}
  \item{y}{A vector of response values. 
      Note that there is no argument checking; 
      if you choose the wrong 
      \code{family} is could cause problems. }
  \item{nfold}{ The number of cross validation folds. }
  \item{verb}{ Whether to print progress through folds. }
  \item{seg}{ Whether to select the minimum or 1se OOS deviance segment;
        see details.}
  \item{...}{ Arguments to \code{gamlr}. }
}
\details{ Fits a \code{gamlr} log penalized regression to the full dataset, and then performs \code{nfold} 
cross validation to evaluate out-of-sample (OOS)
performance  for different penalty weights.

\code{plot.cv.gamlr} can be used to plot the results: it 
shows mean OOS deviance with 1se error bars as a function
of the penalty weight.

The other S3 methods return information relevant to the `best'
penalty: \code{seg="min"} is that with minimum average OOS deviance,
and  \code{seg="1se"} is that whose average OOS deviance is
no more than 1 standard error away from the minimum.
} 
\value{
  \item{gamlr}{ The full-data fitted \code{gamlr} object.}
  \item{family}{ The linear model exponential family. }
  \item{nfold}{ The number of CV folds. }
  \item{cvm}{ Mean OOS deviance by \code{gamlr\$penalty} }
  \item{cvs}{ The standard errors on \code{cvm}.}
  \item{seg.min}{ The index of minimum \code{cvm}. }
  \item{seg.1se}{ The index of \code{1se} \code{cvm} (see details). }
  \item{pen.min}{ Penalty at minimum \code{cvm}. }
  \item{pen.1se}{ Penalty at \code{1se} \code{cvm}. }
 }
\references{
  Taddy 2013.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{

n <- 100
p <- 100

xvar <- matrix(ncol=p,nrow=p)
for(i in 1:p) for(j in i:p) xvar[i,j] <- 0.5^{abs(i-j)}
x <- matrix(rnorm(p*n), nrow=n)\%*\%chol(xvar)
beta <- matrix( (-1)^(1:p)*exp(-(1:p)/10) )
mu = x\%*\%beta
y <- mu + rnorm(n,sd=sd(as.matrix(mu))/2)

cvfit <- cv.gamlr(x, y, verb=TRUE)

coef(cvfit)[1:3,] # 1se default
coef(cvfit, seg="min")[1:3,] # min OOS deviance

par(mfrow=c(1,2))
plot(cvfit)
plot(cvfit$gamlr)

}

\seealso{gamlr}
