
Optimization:

1. '../fast-means/maybeMeans.py' seemed faster than '../fast-means/slowerMeans.py.' The latter just called mean() and var()
methods repeatedly, while the former keeps track of stacks of numbers (which data is in which pile) and calculates only
N averages, average squares, and variances. However, when implemented in cpDetect it did not seem that much faster.

result: move cpDetect.py to old-cpDetect.py

2. tabulated gamma.
	N	with table	without table
	1000	2.116 s		3.827 s
	10000	31.329 s	60.942 s	
