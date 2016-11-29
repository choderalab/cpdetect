
#!/usr/bin/python

import sys
import cPickle
import matplotlib.pyplot as plt

filename = sys.argv[1]
file = open( filename )
u = cPickle.Unpickler( file )
trajectory=u.load()

plt.plot( trajectory, "o", color=(0,0,0) )
plt.show()
