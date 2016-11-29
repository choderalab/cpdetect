
import cPickle
import random

maxnpts=5000;
nstates = 3 # int( random.uniform( 5, 15) )
davg_max = random.uniform( 5, 15 )
dsig_max = 5 # random.uniform( .1, .3 )

# initial values
avg=50
sigma=1
ptrans = 0.1 # random.uniform( 0.001, .05 )

# build the trajectory
npts = 0
state = 0
trajectory=[]
statemod=200
while True :
	
	sample = random.gauss( avg, sigma )
	trajectory.append( sample )
	npts += 1

	# quit if the maximum has been reached
	if npts >= maxnpts : break
	"""
	# test for a new state
	if random.random() < ptrans :
		print "!! change point at %d !!" % npts 
		ptrans = random.uniform( 0.001, .05 )
		avg = avg+random.uniform( -davg_max, davg_max )
		sigma = abs( sigma+random.uniform( -dsig_max, dsig_max ) )
		state += 1
	
	# test if we need a new state
	if state > nstates: break
	"""
	try:
		if (npts%statemod) == 0:
			print "change of state!"
			statemod=statemod+random.randint(-100,100)
			avg += 10
			sigma += 1
	except ZeroDivisionError: pass

	
FILE=open("trajectory-1.dat","w")
p=cPickle.Pickler( FILE )
p.dump( trajectory )
