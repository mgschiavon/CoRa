# Brink motif feedback (v02)
#   with Y represses A synthesis,

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * U)              - ((g + gY) * Y)
		dA  = (mA * (kD/(kD + Y)))  - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		dB  =  mB                   - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		dC  =                       - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		dU  =  mU                   - (g * U) + (bA * A * Us) - (bI * B * U)
		dUs =                      - (g * Us) - (bA * A * Us) + (bI * B * U)
	end g mY gY mA kD mB mU eP e0 bA bI mYs gYs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * U)              - ((g + gY) * Y)
		dA  = (mA * (kD/(kD + Ys))) - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		dB  =  mB                   - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		dC  =                       - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		dU  =  mU                   - (g * U) + (bA * A * Us) - (bI * B * U)
		dUs =                      - (g * Us) - (bA * A * Us) + (bI * B * U)
		dYs = mYs - ((g + gYs) * Ys)
	end g mY gY mA kD mB mU eP e0 bA bI mYs gYs;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mYs] = p[:mY] * ss[5];
		p[:gYs] = p[:gY];
	end;
end