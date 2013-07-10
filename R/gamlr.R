##########################################################
##### path estimation for log penalized regression  ######
##########################################################

## Wrapper function; most happens in c
gamlr <- function(x, y, 
            family=c("gaussian","binomial","poisson"),
            varpen=0, npen=100, 
            pen.start=Inf,  
            pen.min.ratio=0.01, 
            weight=NULL, standardize=TRUE, verb=FALSE,
            thresh=1e-6, maxit=1e5, qn=FALSE)
{
  on.exit(.C("gamlr_cleanup", PACKAGE = "gamlr"))

  ## integer family codes
  family=match.arg(family)
  famid = switch(family, 
    "gaussian"=1, "binomial"=2, "poisson"=3)

  ## data formatting
  y <- drop(y)
  stopifnot(is.null(dim(y)))
  if(is.factor(y)&family=="binomial") y <- as.numeric(y)-1
  y <- as.double(y)
  n <- length(y)

  if(is.data.frame(x)) x <- as.matrix(x)
  x=as(x,"dgCMatrix") 
  p <- ncol(x)
  if(is.null(colnames(x))) colnames(x) <- 1:p
  stopifnot(nrow(x)==n) 


  ## precision weights
  if(is.null(weight)) weight <- rep(1,p)
  stopifnot(all(weight>=0) & length(weight)==p)
  weight <- as.double(weight)

  ## check and clean all arguments
  stopifnot(pen.min.ratio<=1)
  stopifnot(all(c(npen,pen.min.ratio)>0))
  stopifnot(all(c(pen.start,varpen)>=0))
  stopifnot(all(c(thresh,maxit)>0))
  if(is.infinite(varpen)){
    npen=min(npen,sum(weight!=0)+1)
    pen.start=Inf }
  mu <- double(npen)
  mu[1] <- pen.start

  ## drop it like it's hot
  fit <- .C("R_gamlr",
            famid=as.integer(famid), 
            n=n,
            p=p,
            l=length(x@i),
            xi=x@i,
            xp=x@p,
            xv=as.double(x@x),
            y=y,
            weight=weight,
            standardize=as.integer(standardize>0),
            npen=as.integer(npen),
            pminratio=as.double(pen.min.ratio),
            varpen=as.double(varpen),
            thresh=as.double(thresh),
            maxit=as.integer(maxit),
            qn=as.integer(qn>0),
            mu=mu,
            deviance=double(npen),
            df=double(npen),
            alpha=as.double(rep(0,npen)),
            beta=as.double(rep(0,npen*p)),
            exits=integer(npen), 
            verb=as.integer(verb>0),
            PACKAGE="gamlr",
            NAOK=TRUE,
            dup=FALSE)

  ## coefficients
  npen <- fit$npen
  if(npen == 0) stop("could not converge for any penalty.")
  alpha <- head(fit$alpha,npen)
  names(alpha) <- paste0('seg',(1:npen))
  beta <- Matrix(head(fit$beta,npen*p),
                    nrow=p, ncol=npen, 
                    dimnames=list(colnames(x),names(alpha)),
                    sparse=TRUE)
  ## path stats
  pen <- head(fit$mu,npen)
  dev <- head(fit$deviance,npen)
  df <- head(fit$df,npen)
  exits <- head(fit$exits,npen)
  names(df) <- names(dev) <- names(pen) <- names(alpha)

  ## nonzero saturated poisson deviance
  if(family=="poisson")
    dev <- dev + 2*sum(ifelse(y>0,y*log(y),0) - y) 

  if(is.infinite(varpen)) 
    fit$weight <- weight+fit$weight

  ## build return object and exit
  out <- list(penalty=pen, 
             varpen=fit$varpen,
             weight=fit$weight, 
             nobs=fit$n,
             family=family,
             alpha=alpha,
             beta=beta, 
             df=df,
             deviance=dev,
             iterations=fit$maxit,
             call=match.call()) 

  class(out) <- "gamlr"
  invisible(out)
}
 

#### S3 method functions

plot.gamlr <- function(x, against=c("pen","dev"), 
                      col="navy",...)
{
  npen <- ncol(x$beta)
  p <- nrow(x$beta)
  nzr <- unique(x$beta@i)+1
  nzr <- nzr[x$weight[nzr]!=0 & is.finite(x$weight[nzr])]
  beta <- as.matrix(x$beta[nzr,,drop=FALSE])

  if(length(col)==1) col <- rep(col,p)
  col <- col[nzr]

  against=match.arg(against)
  if(against=="pen"){
      xv <- log(x$penalty)
      xvn <- "log penalty"
  } else if(against=="dev"){
      xv <- x$dev
      xvn <- "deviance"
  } else
    stop("unrecognized 'against' argument.")

  argl = list(...)
  if(is.null(argl$ylim)) argl$ylim=range(beta)
  if(is.null(argl$ylab)) argl$ylab="coefficient"
  if(is.null(argl$xlab)) argl$xlab=xvn
  if(is.null(argl$lty)) argl$lty=1
  do.call(plot, c(list(x=xv, y=rep(0,npen), col="grey70", type="l"), argl))

  matplot(xv, t(beta), col=col, add=TRUE, type="l", lty=argl$lty)

  dfi <- unique(round(
    seq(1,npen,length=ceiling(length(axTicks(1))))))
  axis(3,at=log(x$penalty[dfi]), 
    labels=round(x$df[dfi],1),tick=FALSE, line=-.5)

  abline(v=log(x$penalty[which.min(BIC(x))]), 
    lty=3, col="grey50")
}

coef.gamlr <- function(object, 
  select=which.min(BIC(object)), ...){
  if(is.null(select)) select <- 1:ncol(object$beta)
  return(rBind(intercept=object$alpha[select], 
              object$beta[,select,drop=FALSE]))
}

predict.gamlr <- function(object, newdata,
                    select=which.min(BIC(object)),
                    type = c("link", "response"), ...)
{
  stopifnot(inherits(newdata, c("matrix","Matrix")))
  if(is.null(select)) select <- 1:ncol(object$beta)
  
  eta <- matrix(object$alpha[select],
            nrow=nrow(newdata),
            ncol=length(object$alpha[select]),
            byrow=TRUE)
  eta <- eta + newdata%*%object$beta[,select]
                              
  type=match.arg(type)
  if(object$family=="binomial" & type=="response")
    return(apply(eta,2,function(c) 1/(1+exp(-c))))
  if(object$family=="poisson" & type=="response")
    return(apply(eta,2,exp))

  return(eta)
}

summary.gamlr <- function(object, ...){
  print(object)

  return(data.frame(
    pen=object$penalty,
    par=diff(object$b@p)+1,
    df=object$df,
    r2=1-object$dev/object$dev[1],
    bic=BIC(object)))
}

print.gamlr <- function(x, ...){
  cat("\n")
  cat(sprintf(
    "%s gamlr with %d inputs and %d segments.", 
    x$family, nrow(x$beta), ncol(x$beta)))
  cat("\n\n")
}

logLik.gamlr <- function(object, ...){
  ll <- -0.5*object$dev
  attr(ll,"nobs") = object$nobs
  attr(ll,"df") = object$df
  class(ll) <- "logLik"
  ll
}

family.gamlr <- function(object, ...) object$family



