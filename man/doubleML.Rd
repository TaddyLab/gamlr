\name{doubleML}
\alias{doubleML}
\title{double ML}
\description{double (i.e., double) Machine Learning for treatment effect estimation }
\usage{
doubleML(x, d, y, nfold=2, foldid=NULL, family="gaussian", cl=NULL, ...)}
\arguments{
  \item{x}{ Covariates; see \code{gamlr}. }
  \item{d}{ The matrix of treatment variables.  Each column is used as a response by \code{gamlr} during the residualization procedure.}
  \item{y}{ Response; see \code{gamlr}. }
  \item{nfold}{ The number of cross validation folds. }
  \item{foldid}{ An optional length-n vector of fold memberships for each observation.  If specified, this dictates \code{nfold}.}
  \item{family}{ Response model type for the treatment prediction; 
  either "gaussian", "poisson", or "binomial".  This can be either be a single family shared by all columns of \code{d} or a vector of families of length \code{ncol(d)}}
  \item{cl}{possible \code{parallel} library cluster.  If this is not-\code{NULL}, the CV folds are executed in parallel.  This copies the data \code{nfold} times, so make sure you have the memory space.}
  \item{...}{ Arguments to all the \code{gamlr} regressions. }
}
\details{ Performs the double ML procedure of Chernozhukov et al. (2017) to produce an unbiased estimate of the average linear treatment effects of \code{d} on \code{y}.  This procedure uses \code{gamlr} to regress \code{y} and each column of \code{d} onto \code{x}.  In the cross-fitting routine described in Taddy (2019), these regressions are trained on a portion of the data and the out-of-sample residuals are calculated on the left-out fold.  Model selection for these residualization steps is based on the AICc selection rule.  The response residuals are then regressed onto the treatment residuals using \code{lm} and the resulting estimates and standard errors are unbiased for the treatment effects under the assumptions of Chernozhukov et al.
} 
\value{
  A fitted \code{lm} object estimating the treatment effect of \code{d} on \code{y}.  The \code{lm} function has been called with \code{x=TRUE, y=TRUE} such that this object contains the residualized \code{d} as \code{x} and residualized \code{y} as \code{y}.
 }
\author{
  Matt Taddy \email{mataddy@gmail.com}
}
\references{Chernozhukov, Victor and Chetverikov, Denis and Demirer, Mert and Duflo, Esther and Hansen, Christian and Newey, Whitney and Robins, James (The Econometrics Journal, 2017), Double/debiased machine learning for treatment and structural parameters

Matt Taddy, 2019.  Business Data Science, McGraw-Hill}

\examples{

data(hockey)
who <- which(colnames(player)=="SIDNEY_CROSBY")
s <- sample.int(nrow(player),10000) # subsample for a fast example
doubleML(x=player[s,-who], d=player[s,who], y=goal$homegoal[s], standardize=FALSE)

}

\seealso{gamlr, hockey, AICc}
