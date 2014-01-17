##########################################################
##### path estimation for log penalized regression  ######
##########################################################

## Wrapper function; most happens in c
gamlr <- function(x, y, 
            family=c("gaussian","binomial","poisson"),
            gamma=0, 
            nlambda=100, 
            lambda.start=Inf,  
            lambda.min.ratio=0.01, 
            free=NULL, 
            standardize=TRUE, 
            doxx=n>0.5*(p^2),  
            tol=1e-7, 
            maxit=1e4,
            verb=FALSE, ...)
{
  on.exit(.C("gamlr_cleanup", PACKAGE = "gamlr"))

  ## integer family codes
  family=match.arg(family)
  famid = switch(family, 
    "gaussian"=1, "binomial"=2, "poisson"=3)

  ## data formatting (more follows after doxx)
  y <- drop(y)
  stopifnot(is.null(dim(y)))
  if(is.factor(y)&family=="binomial") y <- as.numeric(y)-1
  y <- as.double(y)
  n <- length(y)

  if(inherits(x,"numeric")) x <- matrix(x)
  if(inherits(x,"data.frame")) x <- as.matrix(x)
  if(inherits(x,"simple_triplet_matrix"))
    x <- sparseMatrix(i=x$i,j=x$j,x=x$v,
              dims=dim(x),dimnames=dimnames(x))
  p <- ncol(x)

  ## extras
  xtr = list(...)

  ## alias from glmnet terminology
  if(!is.null(xtr$thresh)) tol = xtr$thresh

  ## fixed shifts for poisson
  eta <- rep(0.0,n)
  if(!is.null(xtr$fix)){
    if(family=="gaussian") y = y-xtr$fix
    else eta <- xtr$fix   } 
  stopifnot(length(eta)==n)
  eta <- as.double(eta)

  ## observation weights
  if(!is.null(xtr$obsweight)){
    obsweight <- xtr$obsweight
    stopifnot(all(obsweight>0))
    stopifnot(length(obsweight)==n)
  } else obsweight <- as.double(rep(1,n))

  ## precalc of x'x
  doxx = doxx & (family=="gaussian")
  if(doxx){
    if(is.null(xtr$xx))
      xtr$xx <- as(
        tcrossprod(t(x*sqrt(obsweight))),
        "matrix")
    xx <- as(xtr$xx,"dspMatrix")
    if(xx@uplo=="L") xx <- t(xx)
    xx <- as.double(xx@x)
  } else xx <- double(0) 

  ## final x formatting
  x=as(x,"dgCMatrix") 
  if(is.null(colnames(x))) 
    colnames(x) <- 1:p
  stopifnot(nrow(x)==n) 

  ## variable weights
  if(!is.null(xtr$varweight)){
    varweight <- xtr$varweight
    stopifnot(all(varweight>=0))
    stopifnot(length(varweight)==p)
  } else{ varweight <- rep(1,p) }
  varweight[free] <- 0
  varweight <- as.double(varweight)

  ## check and clean all arguments
  stopifnot(lambda.min.ratio<=1)
  stopifnot(all(c(nlambda,lambda.min.ratio)>0))
  stopifnot(all(c(lambda.start)>=0))
  stopifnot(all(c(tol,maxit)>0))
  lambda <- double(nlambda)
  lambda[1] <- lambda.start

  ## stepsize
  delta <- exp( log(lambda.min.ratio)/(nlambda-1) )

  ## drop it like it's hot
  fit <- .C("gamlr",
            famid=as.integer(famid), 
            n=n,
            p=p,
            l=length(x@i),
            xi=x@i,
            xp=x@p,
            xv=as.double(x@x),
            y=y,
            doxx=as.integer(doxx),
            xx=xx,
            eta=eta,
            varweight=varweight,
            obsweight=obsweight,
            standardize=as.integer(standardize>0),
            nlambda=as.integer(nlambda),
            delta=as.double(delta),
            gamma=as.double(gamma),
            tol=as.double(tol),
            maxit=as.integer(maxit),
            lambda=as.double(lambda),
            deviance=double(nlambda),
            df=double(nlambda),
            alpha=as.double(rep(0,nlambda)),
            beta=as.double(rep(0,nlambda*p)),
            exits=integer(nlambda), 
            verb=as.integer(verb>0),
            PACKAGE="gamlr",
            NAOK=TRUE,
            dup=FALSE)

  ## coefficients
  nlambda <- fit$nlambda
  if(nlambda == 0){
    warning("could not converge for any lambda.")
    nlambda <- 1
    fit$deviance <- NA
  }

  alpha <- head(fit$alpha,nlambda)
  names(alpha) <- paste0('seg',(1:nlambda))
  beta <- Matrix(head(fit$beta,nlambda*p),
                    nrow=p, ncol=nlambda, 
                    dimnames=list(colnames(x),names(alpha)),
                    sparse=TRUE)

  ## path stats
  lambda <- head(fit$lambda,nlambda)
  dev <- head(fit$deviance,nlambda)
  df <- head(fit$df,nlambda)
  names(df) <- names(dev) <- names(lambda) <- names(alpha)

  ## build return object and exit
  out <- list(lambda=lambda, 
             gamma=fit$gamma,
             nobs=fit$n,
             family=family,
             alpha=alpha,
             beta=beta, 
             df=df,
             deviance=dev,
             totalpass=fit$maxit,
             free=free,
             call=sys.call(1)) 

  class(out) <- "gamlr"
  invisible(out)
}
 

#### S3 method functions

plot.gamlr <- function(x, against=c("pen","dev"), 
                      col="navy", 
                      select=TRUE, df=TRUE, ...)
{
  nlambda <- ncol(x$beta)
  if(!is.null(x$free)) 
    x$beta <- x$beta[-x$free,,drop=FALSE]
  p <- nrow(x$beta) 
  nzr <- unique(x$beta@i)+1
  if(length(nzr)==0) return("nothing to plot")
  beta <- as.matrix(x$beta[nzr,,drop=FALSE])

  if(!is.null(col)){
    if(length(col)==1) col <- rep(col,p)
    col <- col[nzr] }
  else col <- 1:6 # matplot default

  against=match.arg(against)
  if(against=="pen"){
      xv <- log(x$lambda)
      xvn <- "log lambda"
  } else if(against=="dev"){
      xv <- x$dev
      xvn <- "deviance"
  } else
    stop("unrecognized 'against' argument.")

  if(!is.finite(xv[1])) stop("refusing to plot an unconverged fit")

  argl = list(...)
  if(is.null(argl$ylim)) argl$ylim=range(beta)
  if(is.null(argl$ylab)) argl$ylab="coefficient"
  if(is.null(argl$xlab)) argl$xlab=xvn
  if(is.null(argl$lty)) argl$lty=1
  if(is.null(argl$bty)) argl$bty="n"
  do.call(plot, c(list(x=xv, y=rep(0,nlambda), col="grey70", type="l"), argl))

  matplot(xv, t(beta), col=col, add=TRUE, type="l", lty=argl$lty)

  if(df){
    dfi <- unique(round(
      seq(1,nlambda,length=ceiling(length(axTicks(1))))))
    axis(3,at=xv[dfi], labels=round(x$df[dfi],1),tick=FALSE, line=-.5) }

  if(select){
    abline(v=xv[which.min(BIC(x))], lty=3, col="grey20")
    abline(v=xv[which.min(AICc(x))], lty=3, col="grey20")
  }
}

coef.gamlr <- function(object, select=NULL, k=2, ...)
{
  if(length(select)==0)
    select <- which.min(AICc(object,k=k))
  else if(select==0)
   select <- 1:ncol(object$beta)

  select[select>ncol(object$beta)] <- ncol(object$beta)
  return(rBind(intercept=object$alpha[select], 
              object$beta[,select,drop=FALSE]))
}

predict.gamlr <- function(object, newdata,
                    type = c("link", "response"), ...)
{
  if(inherits(newdata,"data.frame")) 
    newdata <- as.matrix(newdata)
  if(inherits(newdata,"simple_triplet_matrix"))
    newdata <- sparseMatrix(
                    i=newdata$i,
                    j=newdata$j,
                    x=newdata$v,
                    dims=dim(newdata),
                    dimnames=dimnames(newdata))

  B <- coef(object, ...)
  eta <- matrix(B[1,],
            nrow=nrow(newdata),
            ncol=ncol(B),
            byrow=TRUE)
  eta <- eta + newdata%*%B[-1,,drop=FALSE]
                              
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
    lambda=object$lambda,
    par=diff(object$b@p)+1,
    df=object$df,
    r2=1-object$dev/object$dev[1],
    bic=BIC(object)))
}

print.gamlr <- function(x, ...){
  cat("\n")
  cat(sprintf(
    "gamma = %g %s gamlr with %d inputs and %d segments.", 
    x$gamma, x$family, nrow(x$beta), ncol(x$beta)))
  cat("\n\n")
}

logLik.gamlr <- function(object, ...){
  ll <- -0.5*object$dev
  attr(ll,"nobs") = object$nobs
  attr(ll,"df") = object$df
  class(ll) <- "logLik"
  ll
}




