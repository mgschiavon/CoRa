# Feedback with buffering (v01)
#   with U inducing Y synthesis
# !! Using linear negative regulation as in Hancock et al (2017)

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * U)      - ((g + gY)  * Y)
		dU  =  mU - (k * Y) - ((g + gU)  * U)  - (b * U) + (bs * Us)
		dUs =               - ((g + gUs) * Us) + (b * U) - (bs * Us)
	end g mY gY mU k gU gUs b bs mYs gYs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * U)      - ((g + gY)  * Y)
		dU  = mU - (k * Ys) - ((g + gU)  * U)  - (b * U) + (bs * Us)
		dUs =               - ((g + gUs) * Us) + (b * U) - (bs * Us)
		dYs =    mYs        - ((g + gYs) * Ys)
	end g mY gY mU k gU gUs b bs mYs gYs;

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