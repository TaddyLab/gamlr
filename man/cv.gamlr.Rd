\name{cv.gamlr}
\alias{cv.gamlr}
\alias{plot.cv.gamlr}
\alias{coef.cv.gamlr}
\alias{predict.cv.gamlr}
\title{ Cross Validation for gamlr}
\description{Cross validation for gamma lasso penalty selection. }
\usage{
cv.gamlr(x, y, nfold=5, foldid=NULL, verb=FALSE, cl=NULL, ...)
\method{plot}{cv.gamlr}(x, select=TRUE, ...)
\method{coef}{cv.gamlr}(object, select=c("1se","min"), ...)
\method{predict}{cv.gamlr}(object, newdata, select=c("1se","min"), ...)
}
\arguments{
  \item{x}{ Covariates; see \code{gamlr}. }
  \item{y}{ Response; see \code{gamlr}. }
  \item{nfold}{ The number of cross validation folds. }
  \item{foldid}{ An optional length-n vector of fold memberships for each observation.  If specified, this dictates \code{nfold}.}
  \item{verb}{ Whether to print progress through folds. }
  \item{cl}{possible \code{parallel} library cluster.  If this is not-\code{NULL}, the CV folds are executed in parallel.  This copies the data \code{nfold} times, so make sure you have the memory space.}
  \item{...}{ Arguments to \code{gamlr}. }
  \item{object}{ A gamlr object.}
  \item{newdata}{ New \code{x} data for prediction.}
  \item{select}{ In prediction and coefficient extraction, 
  select which "best" model to return: 
  \code{select="min"} is that with minimum average OOS deviance,
  and  \code{select="1se"} is that whose average OOS deviance is
  no more than 1 standard error away from the minimum. In \code{plot}, whether to draw these selections.}
}
\details{ Fits a \code{gamlr} regression to the full dataset, and then performs \code{nfold} 
cross validation to evaluate out-of-sample (OOS)
performance  for different penalty weights.

\code{plot.cv.gamlr} can be used to plot the results: it 
shows mean OOS deviance with 1se error bars.
} 
\value{
  \item{gamlr}{ The full-data fitted \code{gamlr} object.}
  \item{nfold}{ The number of CV folds. }
  \item{foldid}{ The length-n vector of fold memberships. }
  \item{cvm}{ Mean OOS deviance by \code{gamlr\$lambda} }
  \item{cvs}{ The standard errors on \code{cvm}.}
  \item{seg.min}{ The index of minimum \code{cvm}. }
  \item{seg.1se}{ The index of \code{1se} \code{cvm} (see details). }
  \item{lambda.min}{ Penalty at minimum \code{cvm}. }
  \item{lambda.1se}{ Penalty at \code{1se} \code{cvm}. }
 }
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\references{Taddy (2013), The Gamma Lasso, http://arxiv.org/abs/1308.5623}

\examples{

n <- 100
p <- 100

xvar <- matrix(ncol=p,nrow=p)
for(i in 1:p) for(j in i:p) xvar[i,j] <- 0.5^{abs(i-j)}
x <- matrix(rnorm(p*n), nrow=n)\%*\%chol(xvar)
beta <- matrix( (-1)^(1:p)*exp(-(1:p)/10) )
mu = x\%*\%beta
y <- mu + rnorm(n,sd=sd(as.vector(mu))/2)

## fit with gamma=1 concavity
cvfit <- cv.gamlr(x, y, gamma=1, verb=TRUE)

coef(cvfit)[1:3,] # 1se default
coef(cvfit, select="min")[1:3,] # min OOS deviance

predict(cvfit, x[1:2,], select="min")
predict(cvfit$gamlr, x[1:2,], select=cvfit$seg.min)

par(mfrow=c(1,2))
plot(cvfit)
plot(cvfit$gamlr)

}

\seealso{gamlr, hockey}
