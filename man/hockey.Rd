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
y <- as.numeric(goal$whoscored=="HOME")

## fit a lambda (gamma=0) model
plot(fit <- gamlr(x=x, y=y, family="binomial",  
      free=1:ncol(config), standardize=FALSE, 
      thresh=1e-8))

## BIC selected player coefficients
B <- coef(fit, s=which.min(BIC(fit)))[colnames(onice),] 
length(B <- B[B!=0])
## print nonzero player effects on who-scored log odds
round(B[order(-B)],3)

## predict the last goal of the 2013 playoffs (this time use AIC)
goal[fit$n,]
predict(fit,x[fit$n,,drop=FALSE],type="response",k=2)
}
\seealso{gamlr}
