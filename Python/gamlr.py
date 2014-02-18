#---------- BEGIN example --------------------------------#
##Following the basic steps from gamlr.R wrapper:
#make x a matrix (numpy) array? 
##observation weights (error check them!)

import numpy 
import SparseMatrix as sm
import numpy.ctypeslib as npct
from ctypes import *

##Here, I walk through the basic demo of gamlr.R to put to gamlr.C

n = 1000
p = 3
#xvar = matrix(0.9, nrow=p,ncol=p)
xvar=numpy.ones(shape=(p,p))*.9
#diag(xvar) <- 1
numpy.fill_diagonal(xvar,1)

#x <- matrix(rnorm(p*n), nrow=n)%*%chol(xvar)
x = numpy.random.standard_normal(p*n)
x.shape = (n,p)
numpy.dot(x,numpy.linalg.cholesky(xvar))
#y <- 4 + 3*x[,1] + -1*x[,2] + rnorm(n)
y = 4+3*x[:,0]+-1*x[:,1]+numpy.random.standard_normal(n)
#---------END example------------------------------------#


#------- BEGIN gamlr.R translation ------------------------
##########################################################
##### path estimation for log penalized regression  ######
##########################################################

## Wrapper function; most happens in c

# gamlr <- function(x, y, 
#             family=c("gaussian","binomial","poisson"),
#             gamma=0, 
#             nlambda=100, 
#             lambda.start=Inf,  
#             lambda.min.ratio=0.01, 
#             free=NULL, 
#             standardize=TRUE, 
#             doxx=n>0.5*(p^2),  
#             tol=1e-7, 
#             maxit=1e4,
#             verb=FALSE, ...)
#{


def gamlr(x, y, family="gaussian",
gamma=0,nlambda=100,
lambda_start = float("inf"),
lambda_min_ratio=.01,
free=None,
standardize=True,
doxx = None, #no idea how he does this in R
tol = 1e-7,
maxit = 1e4,
verb=False):
  #  on.exit(.C("gamlr_cleanup", PACKAGE = "gamlr")) NOT IMPLEMENTED
  ## integer family codes - ADD THESE FEATURES IN LATER N.A.
  #family=match.arg(family)
  #famid = switch(family, 
  #  "gaussian"=1, "binomial"=2, "poisson"=3)
  famid = 1 #We'll have to work this out later
  # # data formatting (more follows after doxx)
  # y <- drop(y)
  # stopifnot(is.null(dim(y)))
  # if(is.factor(y)&family=="binomial") y <- as.numeric(y)-1
  # y <- as.double(y)
  # n <- length(y)
  n = y.shape[0]
  # if(inherits(x,"numeric")) x <- matrix(x)
  # if(inherits(x,"data.frame")) x <- as.matrix(x)
  # if(inherits(x,"simple_triplet_matrix"))
  #   x <- sparseMatrix(i=x$i,j=x$j,x=x$v,
  #             dims=dim(x),dimnames=dimnames(x))
  #p <- ncol(x)
  p = x.shape[1] #second dimension, ie. # columns

  ## extras
  #xtr = list(...)

  ## alias from glmnet terminology
  #if(!is.null(xtr$thresh)) tol = xtr$thresh

  ## fixed shifts for poisson
  #eta <- rep(0.0,n)
  eta = [0.0]*n
  # if(!is.null(xtr$fix)){
  #   if(family=="gaussian") y = y-xtr$fix
  #   else eta <- xtr$fix   } 
  #stopifnot(length(eta)==n)
  assert(len(eta)==n)
  #eta <- as.double(eta)
  eta = numpy.array(eta)

  ## observation weights
  # if(!is.null(xtr$obsweight)){
  #   obsweight <- xtr$obsweight
  #   stopifnot(all(obsweight>0))
  #   stopifnot(length(obsweight)==n)
  # } else obsweight <- as.double(rep(1,n))

  ## precalc of x'x
  # doxx = doxx & (family=="gaussian")
  # if(doxx){
  #   if(is.null(xtr$xx))
  #     xtr$xx <- as(
  #       tcrossprod(t(x*sqrt(obsweight))),
  #       "matrix")
  #   xx <- as(xtr$xx,"dspMatrix")
  #   if(xx@uplo=="L") xx <- t(xx)
  #   xx <- as.double(xx@x)
  # } else xx <- double(0) 

  ## final x formatting
  ###########
  #x=as(x,"dgCMatrix") 
  x = sm.makeSM(x)
  #This turns x into a class SM, with x.ex (matrix) x.i (=x@i) and x.p (=x@p from R)
  #if(is.null(colnames(x))) 
  #  colnames(x) <- 1:p
  #stopifnot(nrow(x)==n) 
  assert(x.ex.shape[0]==n) #since x is now of class makeSM
  ## variable weights
  # if(!is.null(xtr$varweight)){
  #   varweight <- xtr$varweight
  #   stopifnot(all(varweight>=0))
  #   stopifnot(length(varweight)==p)
  # } else{ varweight <- rep(1,p) }
  varweight = numpy.array([1]*p)
  # varweight[free] <- 0
  # varweight <- as.double(varweight)

  ## check and clean all arguments
  #stopifnot(lambda.min.ratio<=1)
  assert(lambda_min_ratio<=1)
  #stopifnot(all(c(nlambda,lambda.min.ratio)>0))
  assert(nlambda>0 and lambda_min_ratio > 0 )
  #stopifnot(all(c(lambda.start)>=0))
  assert(lambda_start >=0)
  #stopifnot(all(c(tol,maxit)>0))
  assert(tol >0 and maxit >0)
  #lambda <- double(nlambda)
  Lambda = numpy.zeros(nlambda) # "lambda" is a function name!!
  #lambda[1] <- lambda_start
  Lambda[0] = lambda_start

  ## stepsize
  #delta <- exp( log(lambda.min.ratio)/(nlambda-1) )
  delta = numpy.exp( numpy.log(lambda_min_ratio)/(nlambda-1) )
  return(delta)
  
  ################
#  _libgamlr = ctypes.CDLL('./gamlr.so')
  
  #alternative method?
_libgamlr = numpy.ctypeslib.load_library('./gamlr.so','.')
_libgamlr.gamlr.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int),
                            POINTER(c_double),POINTER(c_double),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),
                            POINTER(c_int),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_int),
                            POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),
                            POINTER(c_int),POINTER(c_int)]
print("hello! I finished")    
  #
  #fit = _libgamlr.gamlr(
  #  byref(c_int(famid)),
  #  byref(c_int(n)),
  #  byref(c_int(p)),
  #  byref(len(x.i))
  #  byref(c_int(x.i))
  #  byref(c_int(x.v))

#We're going to have to do a lot of work to make sure that all of the moving parts come through ok. 
  ## drop it like it's hot
  # fit <- .C("gamlr",
  #           famid=as.integer(famid), 
  #           n=n,
  #           p=p,
  #           l=length(x@i),
  #           xi=x@i,
  #           xp=x@p,
  #           xv=as.double(x@x),
  #           y=y,
  #           doxx=as.integer(doxx),
  #           xx=xx,12
  #           eta=eta,
  #           varweight=varweight,
  #           obsweight=obsweight,
  #           standardize=as.integer(standardize>0),
  #           nlambda=as.integer(nlambda),
  #           delta=as.double(delta),
  #           gamma=as.double(gamma),
  #           tol=as.double(tol),
  #           maxit=as.integer(maxit),
  #           lambda=as.double(lam2bda),
  #           deviance=double(nlambda),
  #           df=double(nlambda),
  #           alpha=as.double(rep(0,nlambda)),
  #           beta=as.double(rep(0,nlambda*p)),
  #           exits=integer(nlambda), 
  #           verb=as.integer(verb>0),
  #           PACKAGE="gamlr",
  #           NAOK=TRUE,
  #           dup=FALSE)
