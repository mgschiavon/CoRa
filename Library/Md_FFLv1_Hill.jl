# Feedback with feedforward loop (v01)
# !! Using Hill function to describe negative regulation

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * (U + W))                      - ((g + gY) * Y)
		dU  = (mU * ((kD^nH)/((Y^nH) + (kD^nH)))) - ((g + gU) * U)
		dW  = (mW * U)                            - ((g + gW) * W)
	end g mY gY mU kD nH gU mW gW mYs gYs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * (U + W))                       - ((g + gY) * Y)
		dU  = (mU * ((kD^nH)/((Ys^nH) + (kD^nH)))) - ((g + gU) * U)
		dW  = (mW * U)                             - ((g + gW) * W)
		dYs =    mYs                               - ((g + gYs) * Ys)
	end g mY gY mU kD nH gU mW gW mYs gYs;

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