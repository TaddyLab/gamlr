\name{AICc}
\alias{AICc}
\title{Corrected AIC}
\description{ Corrected AIC calculation. }
\usage{
AICc(object, k=2)
}
\arguments{
  \item{object}{ Some model object that you can call \code{logLik} on (such as a \code{gamlr} or \code{glm} fit). }
\item{k}{ The usual \code{AIC} complexity penalty.  \code{k} defaults to 2. }
}
\details{ This works just like usual AIC, but instead calculates the small sample (or high dimensional) corrected version  from Hurvich and Tsai
\deqn{AICc = -2\log LHD + k*df*\frac{n}{n-df-1}.}
}
\value{
  A numeric value for every model evaluated.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\references{Hurvich, C. M. and C-L Tsai, 1989.  "Regression and Time Series Model Selection in Small Samples", Biometrika 76.}

\seealso{gamlr, hockey}
