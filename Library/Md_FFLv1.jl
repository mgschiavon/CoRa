# Feedback with feedforward loop (v01)

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * (U + W))       - ((g + gY) * Y)
		dU  = (mU * (kD/(Y + kD))) - ((g + gU) * U)
		dW  = (mW * U)             - ((g + gW) * W)
	end g mY gY mU kD gU mW gW mYs gYs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * (U + W))        - ((g + gY) * Y)
		dU  = (mU * (kD/(Ys + kD))) - ((g + gU) * U)
		dW  = (mW * U)              - ((g + gW) * W)
		dYs =    mYs                - ((g + gYs) * Ys)
	end g mY gY mU kD gU mW gW mYs gYs;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mYs] = p[:mY] * (ss[2] + ss[3]);
		p[:gYs] = p[:gY];
	end;
end