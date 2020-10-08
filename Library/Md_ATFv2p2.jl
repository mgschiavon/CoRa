# Antithetic feedback (v02 p02)
#   with active W in complex form
#   considering a more complex process

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * (kD/(Y1 + kD))) - ((g + gY) * Y)
		dU  = (mU * Y)              - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)
		dW  =    mW                 - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		dC  =                       - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
		dY0 = (m0 * (W + C))        - ((g + gY) * Y0)
		dY1 = (m1 * (k1/(Y0 + k1))) - ((g + gY) * Y1)
	end g mY kD gY mU gU mW gW e0 eP eM m0 m1 k1 mYs gYs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * (kD/(Y1 + kD))) - ((g + gY) * Y)
		dU  = (mU * Ys)             - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)
		dW  =    mW                 - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		dC  =                       - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
		dY0 = (m0 * (W + C))        - ((g + gY) * Y0)
		dY1 = (m1 * (k1/(Y0 + k1))) - ((g + gY) * Y1)
		dYs =    mYs    - ((g + gYs) * Ys)
	end g mY kD gY mU gU mW gW e0 eP eM m0 m1 k1 mYs gYs;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mYs] = p[:mY] * (p[:kD]/(ss[6] + p[:kD]));
		p[:gYs] = p[:gY];
	end;
end