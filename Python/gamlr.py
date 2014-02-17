##Following the basic steps from gamlr.R wrapper:

##data formatting:
#make x a matrix (numpy) array? 

##observation weights (error check them!)



##Here, I walk through the basic demo of gamlr.R to put to gamlr.C
import numpy 

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

#OK, looks good. 
#now, to run through the rest of the 