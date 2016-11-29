
#!/usr/bin/python
from cPickle import Unpickler
file=open( "trajectory-1.dat" )
out = open( "data.txt", "w" )
u = Unpickler( file )
data = u.load()
for elem in data : out.write( "%f\n" % elem )
file.close()
out.close()
