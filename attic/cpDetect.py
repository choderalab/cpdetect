
# function to quickly calculate the means and means sums of squares of
# all the partitions of a set of points

import cPickle
import numpy
from math import pi
from scipy import array, zeros
from mpmath import log, gamma, mpf # arbitrary float precision!
import sys

inv_log10=1.0/log(10)

def load_testdata( filename = "trajectory.dat" ):
	FILE = open( filename )
	u = cPickle.Unpickler( FILE )
	data = u.load()
	FILE.close()
	return data

def findPoissonChangePoint( data, factorial ):
	# data is a list of counts in each time period, uniformly spaced

	# the denominator (including both P(D|H1) and constant parts of P(D|H2) )
	C = data.sum()
	N = mpf(len(data))
	denominator = factorial[C-1] * pi / ( 2 * N**C )

	# the numerator (trickier)
	# this needs to be averaged over the possible change points 
	weights = zeros(N,dtype=object)
	CA = 0
	CB = C
	for i in range(1,N) :
		# points up through i are in data set A; the rest are in B
		datapoint = data[i-1]	
		NA = mpf(i)   ; CA += datapoint
		NB = mpf(N-i) ; CB -= datapoint
		
		fraction_num = factorial[CA] * factorial[CB] 
		fraction_den = NA**(CA+1) * NB**(CB+1) * ( (CA/NA)**2 + (CB/NB)**2 )
		#weights.append( fraction_num/fraction_den )
		weights[i-1] = mpf(fraction_num)/fraction_den

	numerator = weights.mean()
	lognum= inv_log10 * log( numerator )
	logden= inv_log10 * log( denominator )
	logodds = lognum - logden
	print "num:",numerator, "log num:", lognum, "| denom:", denominator, "log denom:", logden, "|| log odds:", logodds 

	# If there is a change point, then logodds will be greater than 0
	if logodds < 0 : return None
	return ( weights.argmax(), logodds ) 

def findGaussianChangePoint( data, gammatable ):
	N = len( data )
	if N<6 : return None # can't find a cp in data this small

	# the denominator. This is the easy part.
	denom = (pi**1.5) * mpf(( N*data.var() ))**( -N/2.0 + 0.5 ) * gammatable[N]

	# BEGIN weight calculation
	# the numerator. A little trickier.
	weights=[0,0,0] # the change cannot have occurred in the last 3 points
	data2=data**2

	#initialize
	dataA=data[0:3] ; dataA2=data2[0:3] ; NA = len(dataA)
	dataB=data[3:] ; dataB2=data2[3:] ;  NB = len(dataB)
	sumA=dataA.sum() ; sumsqA=dataA2.sum()
	sumB=dataB.sum()  ; sumsqB=dataB2.sum()

	# first data point--this could be done in the loop but it's okay here
	meanA=sumA/NA ; meansumsqA = sumsqA/NA ; meanA2 = meanA**2 ; sA2=meansumsqA-meanA2
	meanB=sumB/NB ; meansumsqB = sumsqB/NB ; meanB2 = meanB**2 ; sB2=meansumsqB-meanB2

	wnumf1 = mpf(NA)**(-0.5*NA + 0.5 ) * mpf(sA2)**(-0.5*NA + 1) * gammatable[NA]
	wnumf2 = mpf(NB)**(-0.5*NB + 0.5 ) * mpf(sB2)**(-0.5*NB + 1) * gammatable[NB]
	wdenom = (sA2 + sB2) * (meanA2*meanB2)
	weights.append( (wnumf1*wnumf2)/wdenom ) 

	for i in range( 3, N-3 ):
		NA += 1	; NB -= 1
		next = data[i]
		sumA += next	; sumB -= next
		nextsq = data2[i]
		sumsqA += nextsq; sumsqB -= nextsq
		meanA=sumA/NA ; meansumsqA = sumsqA/NA ; meanA2 = meanA**2 ; sA2=meansumsqA-meanA2
		meanB=sumB/NB ; meansumsqB = sumsqB/NB ; meanB2 = meanB**2 ; sB2=meansumsqB-meanB2
		wnumf1 = mpf(NA)**(-0.5*NA + 0.5 ) * mpf(sA2)**(-0.5*NA + 1) * gammatable[NA]
		wnumf2 = mpf(NB)**(-0.5*NB + 0.5 ) * mpf(sB2)**(-0.5*NB + 1) * gammatable[NB]
		wdenom = (sA2 + sB2) * (meanA2*meanB2)
		weights.append( (wnumf1*wnumf2)/wdenom) 
	weights.extend( [0,0] ) # the change cannot have occurred at the last 2 points
	weights=array(weights)
	# END weight calculation

	num = 2.0**2.5 * abs(data.mean()) * weights.mean()
	logodds = log( num ) - log( denom ) 	
	print "num:", num, "log num:", log(num), "| denom:", denom, "log denom:", log(denom), "|| log odds:", logodds 
	
	# If there is a change point, then logodds will be greater than 0
	if logodds < 0 : return None
	return ( weights.argmax(), logodds ) 

class ChangePointDetector:
	def __init__( self, data, function, table=None ):
		self.data = data
		self.datalen = len( self.data )
		self.function = function
		self.table=table
		self.changepoints = []
		self.logodds = {}
		self.niter = 0
		self.maxiter = 1000000 # just in case

	def nchangepoints( self ):
		return len( self.changepoints )

	def split_init( self, verbose=False ):
		self.split( 0, self.datalen, verbose )

	def split( self, start, end, verbose=False ):
		if self.niter > self.maxiter :
			print "Change point detection error: number of iterations exceeded"
			print "If this is the right result, you may need to increase"
			print "ChangePointDetector.maxiter (currently %d)" % self.maxiter
			return
		self.niter += 1
		if verbose:
			print "\nIteration %d" % self.niter
			print "Trying to split the segment:", self.data[start:end], "(data from %d to %d)" % ( start, end)
			print self.data[start:end]

		# try to find a change point in the data segment 
		try:
			result = self.function( self.data[ start: end ], self.table )	
		except TypeError: 
			print "trying to test data from %d to %d failed" % ( start,end )
			#print self.data
			raise

		# otherwise, store the cp and call self.split on the two ends
		if result is not None :
			try: # fails if only one value is returned
				logodds = result[1]
				self.logodds[ start+result[0] ] = logodds
				result = start+result[0]
			except TypeError: # must mean it's one number?
				result += start

			if verbose: print "!! change point detected at %d !!" % result
			self.changepoints.append( result )
			self.split( start, result, verbose )
			self.split( result+1, end, verbose )

	def sort( self ): self.changepoints.sort()

	# display the change points
	def show( self ): print self.changepoints

	# show the change points along with the log odds
	def showall( self ):
		for i in range( len( self.changepoints ) ) :
			changepoint = self.changepoints[i]
			try: 
				logodds = self.logodds[ changepoint ]
			except KeyError:
				logodds = None
			print "%d (%f)" % ( changepoint, logodds )

	def largest_logodds( self ):
		return array( self.logodds.values() ).max()
