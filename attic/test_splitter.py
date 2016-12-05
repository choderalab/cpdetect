
#!/usr/bin/python

from cpDetect import ChangePointDetector
from numpy import array, ones, column_stack
from random import randint

NONE = None

datalen=20
data = ones( datalen )
for i in range( randint(1,10 ) ) :
	data[ randint(0, datalen-1) ] = randint(2, 100)
data[17]=101.0
print "data length", len(data)
timepts = array( range(datalen) )

print "arg max:", data.argmax()

def function( data ):
	npts = len(data)
	if npts > 5 :
		max=3+data[3:-2].argmax()
		min=3+data[3:-2].argmin()
		if data[max] > data[min] : 
			return max
		else: return NONE
	else: return NONE

def function2( data ):
	# if the 'data' array is not empty, return 
	if len(data)>0:
		amax = data.argmax()
	else: return NONE 
	if data[amax] > 1: return amax
	else: return NONE
		
cpd= ChangePointDetector( data, function )
cpd.split( 0, len( data ), verbose=True )

print column_stack( (data, timepts) )
#cpd.changepoints.sort()
print cpd.changepoints
