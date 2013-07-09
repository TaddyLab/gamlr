###########################################################
##### cross-validation for log penalized regression  ######
###########################################################

## just an R loop that calls gamlr
cv.gamlr <- function(x, y, nfold=5, verb=FALSE, ...){
  
  full <- gamlr(x,y,...)
  fam <- family(full)

  nfold <- min(nfold,full$nobs)
  rando <- sample.int(full$n)
  foldsize <- floor(full$nobs/nfold)

  argl <- list(...)
  argl$pen.start <- full$pen[1]
  argl$npen <- length(full$pen)
  argl$pen.min.ratio <- full$pen[argl$npen]/argl$pen.start

  oos <- matrix(Inf, nrow=nfold, ncol=length(full$pen),
                dimnames=list(NULL,names(full$pen)))

  if(verb) cat("fold ")
  for(i in 1:nfold){
    train <- rando[-(foldsize*(i-1) + 1:foldsize)]
    fit <- do.call(gamlr, 
      c(list(x=x[train,],y=y[train]), argl))
    eta <- predict(fit, x[-train,], select=NULL)

    dev <- apply(eta,2, 
      function(e) 
        mean(switch(fam, 
          "gaussian" = (e-y[-train])^2, 
          "binomial" = -2*(y[-train]*e - log(1+exp(e))),
          "poisson" = -2*(y[-train]*e - exp(e)))))

    if(fam=="poisson"){
      satnllhd <- mean(ifelse(y[-train]>0,
                      y[-train]*log(y[-train]),
                      0.0) - y[-train]) 
      dev <- dev + 2*satnllhd }
    oos[i,1:length(fit$pen)] <- dev 
    if(verb) cat(sprintf("%d,",i))
  }
  
  cvm <- apply(oos,2,mean)
  cvs <- apply(oos,2,sd)/sqrt(nfold-1)

  seg.min <- which.min(cvm)
  pen.min = full$pen[seg.min]

  cv1se <- (cvm[seg.min]+cvs[seg.min])-cvm
  seg.1se <- min((1:length(cvm))[cv1se>=0])
  pen.1se = full$pen[seg.1se]

  if(verb) cat("done.\n")
  out <- list(gamlr=full,
          family=fam,
          nfold=nfold,
          cvm=cvm,
          cvs=cvs,
          seg.min=seg.min,
          seg.1se=seg.1se,
          pen.min=pen.min,
          pen.1se=pen.1se)
  class(out) <- "cv.gamlr"
  invisible(out)
}

## S3 method functions

plot.cv.gamlr <- function(x, ...){
  argl = list(...)
  if(is.null(argl$xlab)) argl$xlab="log penalty"
  if(is.null(argl$ylab)) argl$ylab=sprintf("%s deviance",x$family)
  if(is.null(argl$pch)) argl$pch=20
  if(is.null(argl$col)) argl$col=4

  cvlo <- x$cvm-x$cvs
  cvhi <- x$cvm+x$cvs

  if(is.null(argl$ylim)) argl$ylim=range(c(cvlo,cvhi))

  argl$x <- log(x$gamlr$pen)
  argl$y <- x$cvm
  argl$type <- "n"

  suppressWarnings(do.call(plot, argl))
  segments(x0=argl$x, y0=cvlo, y1=cvhi, col="grey70")
  argl$type <- NULL
  suppressWarnings(do.call(points, argl))

  abline(v=log(x$pen.min), lty=3, col="grey20")
  abline(v=log(x$pen.1se), lty=3, col="grey20")

  dfi <- unique(round(
    seq(1,length(argl$x),length=ceiling(length(axTicks(1))))))
  axis(3,at=argl$x[dfi], 
    labels=round(x$gamlr$df[dfi],1),tick=FALSE, line=-.5)

}

coef.cv.gamlr <- function(object, 
                          select=c("1se","min"), ...){
  seg = paste("seg",match.arg(select),sep=".")
  coef(object$gamlr, select=object[[seg]])
}

predict.cv.gamlr <- function(object, newdata,
                          select=c("1se","min"), ...){
  seg = paste("seg",match.arg(select),sep=".")
  predict.gamlr(object$gamlr, newdata, select=object[[seg]], ...)
}

summary.cv.gamlr <- function(object, ...){
  print(object)

  return(data.frame(
    pen=object$gamlr$penalty,
    par=diff(object$gamlr$b@p)+1,
    oos.r2=1-object$cvm/object$cvm[1]))
}

print.cv.gamlr <- function(x, ...){
  cat("\n")
  cat(sprintf(
    "%d-fold %s cv.gamlr object", 
    x$nfold, x$gamlr$family))
  cat("\n\n")
}


family.cv.gamlr <- function(object, ...) object$family













