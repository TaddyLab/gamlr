## xnaref: make missing (NA) the reference level of a factor
xnaref <- function(x){
	if(is.character(x)) x <- factor(x)
	if(is.factor(x))
		if(!is.na(levels(x)[1]))
			x <- factor(x,levels=c(NA,levels(x)),exclude=NULL)
	return(x) 
}

## ximpute: mean or zero imputation for missing values
ximpute <- function(x, pzero){
	if(is.numeric(x)){
		nas <- is.na(x)
		if(any(nas)){
			if(mean(x==0,na.rm=TRUE) > pzero) val <- 0
			else val <- mean(x, na.rm=TRUE)
			miss <- is.na(x)
			x[is.na(x)] <- val
			x <- cbind(x,miss)
		}
	} 
	return(x)
}

## wrapper to apply both functions to a dataframe
naref <- function(x, impute=FALSE, pzero=0.5){
	if(!is.data.frame(x)) 
		stop("You need to give me a data.frame.")
	x <- lapply(x, xnaref)
	if(impute) x <- lapply(x, ximpute, pzero=pzero)
	return(as.data.frame(x))
}