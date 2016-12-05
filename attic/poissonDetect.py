
#!/usr/bin/python

import matplotlib.pyplot as plt
from time import time
from numpy.random import poisson
from numpy import concatenate, ones
from cpDetect import *
from mpmath import gamma

l1=5	; l2=8
n1=150	; n2=100
d1=poisson( l1, n1 )
d2=poisson( l2, n2 )
trajectory = concatenate( (d1,d2,d1,d1,d2) )

# build factorial table
factorial=[]
maxc = trajectory.sum()+2
for i in range(1,maxc): 
	factorial.append( gamma(i) )

cpd = ChangePointDetector( trajectory,findPoissonChangePoint,factorial)

# the part we want to time
t0 = time()
print "splitting"
cpd.split_init()
t1 = time()
print "took %3.3f seconds" % ( t1-t0 )
cpd.sort()
cpd.showall()
print "found %d change points" % cpd.nchangepoints()

times = array( range( len(trajectory) ) )

ymax= trajectory.max() * 1.01 
ymin= trajectory.min() * 0.99 

# get the shade scale. Any cp will have log B > 0, so scale to the biggest one
for i in range( cpd.nchangepoints() ):
	changepoint = cpd.changepoints[i] 
	y = ( ymin, ymax )
	x = changepoint * ones( len(y) ) 
	logodds = cpd.logodds[changepoint]

	if logodds > 1: color=(0,0,0) # black
	elif logodds > 0.477: color=(.5,.5,.5) # dark gray
	else : color=(.7,.7,.7) # light gray

	print "point ", changepoint, "log odds", logodds, "color", color
	plt.plot( x, y, "-", color=color, linewidth=1.5 )

plt.plot( times, trajectory, "b-" )
plt.savefig( "figure.eps" )
plt.show()
