# Lineages Renewal Control (v01)
#   inspired by Buzi et al (2015)

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY = (2 * X * (kY * ((Y/oY)^n))/(kA + kX + (kY * ((Y/oY)^n)))) - (g * Y)
		dX = (2 * X *                kX/(kA + kX + (kY * ((Y/oY)^n)))) - X
		dZ = (2 * X *                kA/(kA + kX + (kY * ((Y/oY)^n)))) - (g * Z)
	end g kX kY kA oY n n mYs gYs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY = (2 * X * (kY * ((Y/oY)^n))/(kA + kX + (kY * ((YS/oY)^n)))) - (g * Y)
		dX = (2 * X *                kX/(kA + kX + (kY * ((YS/oY)^n)))) - X
		dZ = (2 * X *                kA/(kA + kX + (kY * ((YS/oY)^n)))) - (g * Z)
		dYs = mYs - (gYs * Y)
	end g kX kY kA oY n mYs gYs;

	# Define system's output (e.g. total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mYs] = (2 * ss[2] * (kY * ((ss[1]/oY)^n))/(kA + kX + (kY * ((ss[1]/oY)^n))));
		p[:gYs] = p[:g];
	end;
end
