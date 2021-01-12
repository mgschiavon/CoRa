# Lineages Renewal Control (v01)
#   inspired by Buzi et al (2015)

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY = (2 * X * (kY * ((Y/oY)^n))/(kZ + kX + (kY * ((Y/oY)^n)))) - (g * Y)
		dX = (2 * X *                kX/(kZ + kX + (kY * ((Y/oY)^n)))) - X
		dZ = (2 * X *                kZ/(kZ + kX + (kY * ((Y/oY)^n)))) - (g * Z)
	end g kX kY kZ oY n mYs gYs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY = (2 * X * (kY * ((Ys/oY)^n))/(kZ + kX + (kY * ((Ys/oY)^n)))) - (g * Y)
		dX = (2 * X *                 kX/(kZ + kX + (kY * ((Ys/oY)^n)))) - X
		dZ = (2 * X *                 kZ/(kZ + kX + (kY * ((Ys/oY)^n)))) - (g * Z)
		dYs = mYs - (gYs * Ys)
	end g kX kY kZ oY n mYs gYs;

	# Define system's output (e.g. total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mYs] = (2 * ss[2] * (p[:kY] * ((ss[1]/p[:oY])^p[:n]))/(p[:kZ] + p[:kX] + (p[:kY] * ((ss[1]/p[:oY])^p[:n]))));
		p[:gYs] = p[:g];
	end;
end
