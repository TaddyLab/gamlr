## naref: make missing (NA) the reference level of a factor
xnaref <- function(x){
	if(is.factor(x))
		if(!is.na(levels(x)[1]))
			x <- factor(x,levels=c(NA,levels(x)),exclude=NULL)
	return(x) }

naref <- function(x){
	if(is.null(dim(x))) return(xnaref(x))
	if(!is.data.frame(x)) 
		stop("You need to give me a data.frame or a factor.")
	x <- lapply(x, xnaref)
	return(as.data.frame(x))
}