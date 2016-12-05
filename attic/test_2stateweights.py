
#!/usr/bin/python

import scipy
from cpDetect import *
import matplotlib.pyplot as plt
from mpmath import log, mpf

FILE=open( "trajectory-1.dat", "r" )
u = cPickle.Unpickler( FILE )
trajectory = scipy.array( u.load() )

wts = calc_twostate_weights(trajectory)
for i in range( len(wts) ): print i, wts[i]

# normalize the weights
tot = mpf( 0.0 )
for elem in wts: tot+=elem
print tot

logwts = []
for elem in wts : logwts.append( log(elem,10 ))

for i in range( len(logwts) ): print i, logwts[i]

plt.plot( logwts )
plt.show()

