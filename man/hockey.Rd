\name{hockey}
\alias{hockey}

\title{ NHL hockey data }

\description{Every NHL goal from fall 2002 through the 2013 cup finals.}

\details{ The data comprise of information about 
  play configuration and the
  players on ice  (including goalies) for every
  goal from 2002-03 to 2012-13 NHL seasons. 
  Collected using A. C. Thomas's \code{nlhscrapr} package. 
  See the Chicago hockey analytics project at \code{github.com/mataddy/hockey}.}
\value{
   \item{goal}{ Info about each goal scored.  }
   \item{onice}{ Matrix with entries for who was on the ice for each goal: +1 for a home team player, -1 for an away
	 team player, zero otherwise. }
   \item{config}{ Special teams info. For example,
   \code{S5v4} is a 5 on 4 powerplay,
   +1 if it is for the home-team and -1 for the away team.}
   \item{player}{ Player info: names, most common position, and plus-minus.}
 }

\references{Gramacy, Jensen, and Taddy (2013): "Estimating Player
  Contribution in Hockey with Regularized Logistic Regression."
  http://arxiv.org/abs/1209.5026.
}

\author{
Matt Taddy, \email{taddy@chicagobooth.edu}
}

\examples{
data(hockey)
x <- cBind(config,onice)
y <- as.numeric(goal$who=="HOME")
fit <- gamlr(x=x, y=y, family="binomial", nlam=10,
            standardize=FALSE, varpen=1e-2, verb=TRUE)

## BIC selected coefficients
sum( (B <- as.matrix(coef(fit))) !=0)
round(B[order(-abs(B))[1:20],],2)

## fit plot
plot(fit, col=c(rep("slategrey",ncol(config)),
                rep("darkorange",ncol(onice))))

## predict the last goal of the 2013 playoffs
goal[fit$n,]
predict(fit,x[fit$n,,drop=FALSE],type="response")
}
\seealso{gamlr}
