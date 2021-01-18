\name{naref}
\alias{naref}
\title{NA reference level}
\description{ Set NA as the reference level for factor variables and do imputation on missing values for numeric variables.  This is useful to build model matrices for regularized regression, and for dealing with missing values, as in Taddy 2019.}
\usage{
naref(x, impute=FALSE, pzero=0.5)
}
\arguments{
\item{x}{ A data frame. }
\item{impute}{ Logical, whether to impute missing values in numeric columns.}
\item{pzero}{ If \code{impute==TRUE}, then if more than \code{pzero} of the values in a column are zero do zero imputation, else do mean imputation.}
}
\details{ For every \code{factor} or \code{character} column in \code{x}, \code{naref} sets \code{NA} as the reference level for a \code{factor} variable. Columns coded as \code{character} class are first converted to factors via \R{factor(x)}.  If \code{impute=TRUE} then the numeric columns are converted to two columns, one appended \code{.x} that contains imputed values and another appended \code{.miss} which is a binary variable indicating whether the original value was missing.  Numeric columns are returned without change if \code{impute=FALSE} or if they do not contain any missing values.  
}
\value{
  A data frame where the factor and character columns have been converted to factors with reference level \code{NA}, and if \code{impute=TRUE} the missing values in numeric columns have been imputed and a flag for missingness has been added.  See details.   
}
\author{
  Matt Taddy \email{mataddy@gmail.com}
}
\references{Matt Taddy, 2019.  "Business Data Science", McGraw-Hill}
\examples{
( x <- data.frame(a=factor(c(1,2,3)),b=c(1,NA,3)) )
naref(x, impute=TRUE)
}
