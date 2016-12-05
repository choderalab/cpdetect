
#!/usr/bin/python

from cpDetect import *

data = load_testdata()
#print "data:", data
print "number of points:", len(data)

mean_var_array=calc_mean_mss( data )

print numpy.array( mean_var_array )
print "number of points:", len(mean_var_array)
