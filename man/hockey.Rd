\name{hockey}
\alias{hockey}

\title{ NHL hockey data }

\description{Every even-strength goal from 2007-2011.  Collected from
  \code{nhl.com} by Andrew Thomas.  }

\details{ The data comprise of information about teams
  and the players on ice (including goalies) for every
  even strength goal during the four regular seasons of 2007-2008
  through 2010-2011.  There were  1467 players involved in  18154
  goals. All matrices are stored as slam simple triplet matrices.
}

\value{
   \item{Ygoals}{ Response: +1 for a home goal, -1 for an away goal.}
   \item{Xteam}{ Team info: +1 for the home team, -1 for the away team,
	 zero otherwise. }
   \item{Xplayer}{ Player info: +1 for a home team player, -1 for an away
	 team player, zero otherwise. }
   \item{plusminus}{ A list containing team and player plus-minus over
	 the four seasons.}
 }

\references{Gramacy, Jensen, and Taddy (2013): "Estimating Player
  Contribution in Hockey with Regularized Logistic Regression."
  http://arxiv.org/abs/1209.5026.

  Thomas, Ventura, Jensen, and Ma (2012):
  "Competing Process Hazard Function Models for Player Ratings in Ice Hockey." http://arxiv.org/abs/1209.5026.
}

\author{
Matt Taddy, \email{taddy@chicagobooth.edu}
}

\examples{
data(hockey)

## fit a regularization path for 20 significant players
## team effects are unpenalized
fit <- gamlr(covars=cbind(Xteam,Xplayer), response=Ygoals,
               free=1:30, verb=1, family="binomial", scale=FALSE,
               stoparg=list(rule="maxvar", val=20))	
plot(fit)
summary(fit)

}