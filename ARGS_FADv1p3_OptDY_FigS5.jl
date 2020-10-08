mrw  = (pOp  = [:mU,:mW,:eP],	# Parameters to optimize
		pMin = [-3,-3,-3],		# Minimum parameter value to explore (log10)
		pMax = [3,3,3],			# Maximum parameter value to explore (log10)
		runs = 1,			# Number of optimization runs
		iter = 1000,			# Number of iterations per optimization run
		cov  = [0.1,0.1,0.1],	# Covariance to calculate parameter random walk
		M    = 10,				# "Mutation step size" for multiplicative random walk
		temp = 0,				# Flag for simulated annealing (if 0, MRW)
		prtW = 1,				# Flag for printing each walk step
		prtD = 1);				# Flag for printing full DY curve