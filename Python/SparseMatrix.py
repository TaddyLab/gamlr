import numpy
class makeSM(object):
	def __init__(self,x):
			#Check that we were given a numpy array
			assert(isinstance(x,numpy.ndarray))
			#self._as_parameter_ = x sigh..this doesn't work
			self.ex = x
			self.x = []
			self.i = []
			self.p = [0]
			#we'll fill these in shortly
			pcount=0
			for j in xrange(x.shape[1]):
				for i in xrange(x.shape[0]):
					if x[i,j] != 0.0:
						pcount+=1
						self.i.append(i)
						self.x.append(x[i,j])
				self.p.append(pcount)
