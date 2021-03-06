# Feedback by active degradation + positive feedback (v02)
#   where only W is degraded in the complex,
#   W induces its own synthesis
#   and active W in complex form

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY = (mY * (W + C))                  - ((g + gY) * Y)
		dU = (mU * Y)                        - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)
		dW = (mW * ((W + C)/((W + C) + kD))) - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		dC =                                 - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
	end g mY gY mU gU mW kD gW e0 eP eM mYs gYs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * (W + C))                  - ((g + gY) * Y)
		dU  = (mU * Ys)                       - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)
		dW  = (mW * ((W + C)/((W + C) + kD))) - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		dC  =                                 - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
		dYs =    mYs                          - ((g + gYs) * Ys)
	end g mY gY mU gU mW kD gW e0 eP eM mYs gYs;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mYs] = p[:mY] * (ss[3] + ss[4]);
		p[:gYs] = p[:gY];
	end;
end