##########################################################
##### path estimation for log penalized regression  ######
##########################################################

## Wrapper function; most happens in c
gamlr <- function(x, y, 
            family=c("gaussian","binomial","poisson"),
            gamma=0, nlambda=100, 
            lambda.start=Inf,  
            lambda.min.ratio=0.01, 
            free=NULL, standardize=TRUE, 
            thresh=1e-7, maxit=1e3,
            verb=FALSE, ...)
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

  if(inherits(x,"data.frame")) x <- as.matrix(x)
  if(inherits(x,"simple_triplet_matrix"))
    x <- sparseMatrix(i=x$i,j=x$j,x=x$v,
              dims=dim(x),dimnames=dimnames(x))
  x=as(x,"dgCMatrix") 
  if(is.null(colnames(x))) 
    colnames(x) <- 1:ncol(x)
  stopifnot(nrow(x)==n) 
  p <- ncol(x)

  ## weight
  weight <- rep(1,p)
  weight[free] <- 0
  weight <- as.double(weight)

  ## check and clean all arguments
  stopifnot(lambda.min.ratio<=1)
  stopifnot(all(c(nlambda,lambda.min.ratio)>0))
  stopifnot(all(c(lambda.start)>=0))
  stopifnot(all(c(thresh,maxit)>0))
  lambda <- double(nlambda)
  lambda[1] <- lambda.start

  ## stepsize
  delta <- exp( log(lambda.min.ratio)/(nlambda-1) )

  ## extras
  xtr = list(...)

  ## fixed shifts
  if(!is.null(xtr$fix))
    eta <- xtr$fix
  else
    eta <- rep(0.0,n)
  stopifnot(length(eta)==n)
  eta <- as.double(eta)

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
            eta=eta,
            weight=weight,
            standardize=as.integer(standardize>0),
            nlambda=as.integer(nlambda),
            delta=as.double(delta),
            gamma=as.double(gamma),
            thresh=as.double(thresh),
            maxit=as.integer(maxit),
            lambda=lambda,
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
  if(nlambda == 0) stop("could not converge for any lambda.")
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

  ## nonzero saturated poisson deviance
  if(family=="poisson")
    dev <- dev + 2*sum(ifelse(y>0,y*log(y),0) - y) 

  ## build return object and exit
  out <- list(lambda=lambda, 
             gamma=fit$gamma,
             weight=fit$weight, 
             nobs=fit$n,
             family=family,
             alpha=alpha,
             beta=beta, 
             df=df,
             deviance=dev,
             totalpass=fit$maxit,
             call=match.call()) 

  class(out) <- "gamlr"
  invisible(out)
}
 

#### S3 method functions

plot.gamlr <- function(x, against=c("pen","dev"), 
                      col=rgb(0,0,.5,.75), 
                      select=TRUE, df=TRUE, ...)
{
  nlambda <- ncol(x$beta)
  p <- nrow(x$beta)
  nzr <- unique(x$beta@i)+1
  beta <- as.matrix(x$beta[nzr,,drop=FALSE])

  if(length(col)==1) col <- rep(col,p)
  col <- col[nzr]

  against=match.arg(against)
  if(against=="pen"){
      xv <- log(x$lambda)
      xvn <- "log lambda"
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
  if(is.null(argl$bty)) argl$bty="n"
  do.call(plot, c(list(x=xv, y=rep(0,nlambda), col="grey70", type="l"), argl))

  matplot(xv, t(beta), col=col, add=TRUE, type="l", lty=argl$lty)

  if(df){
    dfi <- unique(round(
      seq(1,nlambda,length=ceiling(length(axTicks(1))))))
    axis(3,at=xv[dfi], labels=round(x$df[dfi],1),tick=FALSE, line=-.5) }

  if(select){
    abline(v=xv[which.min(BIC(x))], lty=3, col="grey20")
    abline(v=xv[which.min(AIC(x))], lty=3, col="grey20") }
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
  if(inherits(newdata,"data.frame")) 
    newdata <- as.matrix(newdata)
  if(inherits(newdata,"simple_triplet_matrix"))
    newdata <- sparseMatrix(
                    i=newdata$i,
                    j=newdata$j,
                    x=newdata$v,
                    dims=dim(newdata),
                    dimnames=dimnames(newdata))

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

family.gamlr <- function(object, ...) object$family



