# double ML with cross-fitting

doubleML <- function(x, d, y, nfold=2, foldid=NULL, 
					family="gaussian", cl=NULL, ...)
{
	# check arguments
	nobs <- length(y)
	if(is.null(dim(x))) x <- as.matrix(x)
	if(is.null(dim(d)) | is.data.frame(d)) d <- as.matrix(d)
	if(nrow(x)!=nobs | nrow(d)!=nobs) 
		stop("mismatch in number of observations in x, y, d")
	if(length(family)!=ncol(d)){
		if(length(family)==1) family <- rep(family, ncol(d))
		else stop("argument family length doesn't match number of d columns")
	}
	argl <- list(...)

	# randomly split data into folds
  	if(is.null(foldid)){
    	foldid <- rep.int(1:nfold, 
    		times = ceiling(nobs/nfold))[sample.int(nobs)]
  	} else stopifnot(length(foldid)==nobs)
    I <- split(1:nobs, foldid)
  	nfold <- length(I)

 	# function to run each OOS orthogonalization
	getresids <- function(k){
		# response
		yfit <- do.call(gamlr, 
			c(list(x=x[-I[[k]],], y=y[-I[[k]]]), argl))
		yhat <- predict(yfit, x[I[[k]],], type="response")
		ytil <- drop(y[I[[k]]] - yhat)
		# treatments
		dtil <- c()
		for(j in 1:ncol(d)){
			dfit <- do.call(gamlr,
				c(list(x=x[-I[[k]],], y=d[-I[[k]],j], family=family[j]), argl))
			dhat <- predict(dfit, x[I[[k]],], type="response")
			dtil <- cbind(dtil, drop(d[I[[k]],j] - dhat))
		}
		return(cbind(ytil, dtil))
	}

	# apply the getresids function
  	if(!is.null(cl)){
   	 if (requireNamespace("parallel", quietly = TRUE)) {
    	  parallel::clusterExport(cl,
        	c("x","d","y","I","argl","family"), 
        	envir=environment())
      	resids <- parallel::parLapply(cl,1:nfold,getresids)
    	} else {
     		warning("cl is not NULL, but parallel package unavailable.")
      		cl <- NULL
   	 	}
  	}
  	if(is.null(cl)) resids <- lapply(1:nfold,getresids)

	# rename and return
	for(k in 1:length(resids)){
		y[I[[k]]] <- resids[[k]][,1]
		d[I[[k]],] <- resids[[k]][,-1]
	}
	lm( y ~ d - 1, x=TRUE, y=TRUE ) 
}

