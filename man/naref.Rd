\name{naref}
\alias{naref}
\title{NA reference level}
\description{ Set NA as the reference level for factor variables. }
\usage{
naref(x)
}
\arguments{
\item{x}{ A data frame or vector of factor variables. }
}
\details{ \code{naref} sets \code{NA} as the reference level for a factor variable or for all factor variables in a data frame.  This is useful to build model matrices for regularized regression, and for dealing with missing values, as in Taddy 2019.
}
\value{
  The factor vector if \code{x} is a vector, or the data frame containing factor columns if \code{x} is a dataframe, are returned with \code{NA} as the reference level(s).  
}
\author{
  Matt Taddy \email{mataddy@gmail.com}
}
\references{Matt Taddy, 2019.  "Business Data Science", McGraw-Hill}
\examples{
x <- factor(c(1,2,3))
levels(x)
naref(x)
}
