# Feedback with buffering (v01)
#   with U inducing Y synthesis
# !! Using Hill function to describe negative regulation

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * U)                            - ((g + gY) * Y)
		dU  = (mU * ((kD^nH)/((Y^nH) + (kD^nH)))) - ((g + gU) * U)  - (b * U) + (bs * Us)
		dUs =                                     - ((g + gU) * Us) + (b * U) - (bs * Us)
	end g mY gY mU kD nH gU b bs mYs gYs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * U)                             - ((g + gY) * Y)
		dU  = (mU * ((kD^nH)/((Ys^nH) + (kD^nH)))) - ((g + gU) * U)  - (b * U) + (bs * Us)
		dUs =                                      - ((g + gU) * Us) + (b * U) - (bs * Us)
		dYs =    mYs    - ((g + gYs) * Ys)
	end g mY gY mU kD nH gU b bs mYs gYs;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mYs] = p[:mY] * ss[2];
		p[:gYs] = p[:gY];
	end;
end