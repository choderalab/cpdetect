
#!/usr/bin/python

from time import time
from old2cpDetect import *
import matplotlib.pyplot as plt
import scipy
from numpy import ones
from mpmath import mpf, gamma
import cPickle

FILE=open( "trajectory-1.dat", "r" )
u = cPickle.Unpickler( FILE )
trajectory = scipy.array( u.load() )

cpd = ChangePointDetector( trajectory,findGaussianChangePoint )

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
