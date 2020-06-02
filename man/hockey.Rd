\name{hockey}
\alias{hockey}
\alias{config}
\alias{goal}
\alias{team}
\alias{player}

\title{ NHL hockey data }

\description{Every NHL goal from fall 2002 through the 2014 cup finals.}

\details{ The data comprise of information about 
  play configuration and the
  players on ice  (including goalies) for every
  goal from 2002-03 to 2012-14 NHL seasons. 
  Collected using A. C. Thomas's \code{nlhscrapr} package. 
  See the Chicago hockey analytics project at \code{github.com/mataddy/hockey}.}
\value{
   \item{goal}{ Info about each goal scored, including \code{homegoal} -- an indicator for the home team scoring.}
   \item{player}{ Sparse Matrix with entries for who was on the ice for each goal: +1 for a home team player, -1 for an away team player, zero otherwise. }
   \item{team}{ Sparse Matrix with indicators for each team*season interaction: +1 for home team, -1 for away team. }
   \item{config}{ Special teams info. For example,
   \code{S5v4} is a 5 on 4 powerplay,
   +1 if it is for the home-team and -1 for the away team.}

 }

\references{Gramacy, Jensen, and Taddy (2013): "Estimating Player
  Contribution in Hockey with Regularized Logistic Regression", the Journal of Quantitative Analysis in Sport.

  Gramacy, Taddy, and Tian (2015): "Hockey Player Performance via Regularized Logistic Regression", the Handbook of statistical methods for design and analysis in sports.
}

\author{
Matt Taddy, \email{mataddy@gmail.com}
}

\examples{
## design 
data(hockey)
x <- cbind(config,team,player)
y <- goal$homegoal

## fit the plus-minus regression model
## (non-player effects are unpenalized)

fit <- gamlr(x, y, 
  lambda.min.ratio=0.05, nlambda=40, ## just so it runs in under 5 sec
  free=1:(ncol(config)+ncol(team)),
  standardize=FALSE, family="binomial")
plot(fit)

## look at estimated player [career] effects
B <- coef(fit)[colnames(player),]
sum(B!=0) # number of measurable effects (AICc selection)
B[order(-B)[1:10]] # 10 biggest

## convert to 2013-2014 season partial plus-minus
now <- goal$season=="20132014"
pm <- colSums(player[now,names(B)]*c(-1,1)[y[now]+1]) # traditional plus minus
ng <- colSums(abs(player[now,names(B)])) # total number of goals
# The individual effect on probability that a
# given goal is for vs against that player's team
p <- 1/(1+exp(-B)) 
# multiply ng*p - ng*(1-p) to get expected plus-minus
ppm <- ng*(2*p-1)

# organize the data together and print top 20
effect <- data.frame(b=round(B,3),ppm=round(ppm,3),pm=pm)
effect <- effect[order(-effect$ppm),]
print(effect[1:20,])
}
\seealso{gamlr}
