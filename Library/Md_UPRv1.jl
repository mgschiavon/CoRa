# Unfolded protein response (v01)
#   intracellular signaling pathway
#   as modeled in Pincus et al. (2010)

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dU   = sU - (cB * U * B) + (gUB * UB) - (cD * U) + (cE * E * Ud) + (gB * UB)
		dUB  = (cB * U * B) - (gUB * UB) - (cD * UB) + (cE * E * UdB) - (gB * UB) - (gF * UB)
		dUd  = - (cB * Ud * B) + (gUB * UdB) + (cD * U) - (cE * E * Ud) + (gB * UdB)
		dUdB = (cB * Ud * B) - (gUB * UdB) + (cD * UB) - (cE * E * UdB) - (gB * UdB)
		dI   = - (cB * B * I) - (cA * (U + Ud + UdB) * I) + (gIB * IB) + (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dIB  = (cB * B * I) - (gIB * IB)
		dIa  = (cA * (U + Ud + UdB) * I) - (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dHu  = bHu - (gHs * Hu) - (bHs * min(Ia,Hu))
		dHs  = - (gHs * Hs) + (bHs * min(Ia,Hu))
		dB   = (bB * (1 + (nB * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gB * B) - (cB * U * B) + (gUB * UB) - (cB * Ud * B) + (gUB * UdB) + (gF * UB)
		dE   = (bE * (1 + (nE * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gE * E)
	end sU cB gUB cD cE gB gF cA gIB gIA kI nI bHu gHs bHs bB nB a0 a1 gB bE nE gE uT;
	# ODE system without feedback
	odeNF = @ode_def begin
		dU   = sU - (cB * U * B) + (gUB * UB) - (cD * U) + (cE * E * Ud) + (gB * UB)
		dUB  = (cB * U * B) - (gUB * UB) - (cD * UB) + (cE * E * UdB) - (gB * UB) - (gF * UB)
		dUd  = - (cB * Ud * B) + (gUB * UdB) + (cD * U) - (cE * E * Ud) + (gB * UdB)
		dUdB = (cB * Ud * B) - (gUB * UdB) + (cD * UB) - (cE * E * UdB) - (gB * UdB)
		dI   = - (cB * B * I) - (cA * uT * I) + (gIB * IB) + (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dIB  = (cB * B * I) - (gIB * IB)
		dIa  = (cA * uT * I) - (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dHu  = bHu - (gHs * Hu) - (bHs * min(Ia,Hu))
		dHs  = - (gHs * Hs) + (bHs * min(Ia,Hu))
		dB   = (bB * (1 + (nB * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gB * B) - (cB * U * B) + (gUB * UB) - (cB * Ud * B) + (gUB * UdB) + (gF * UB)
		dE   = (bE * (1 + (nE * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gE * E)
	end sU cB gUB cD cE gB gF cA gIB gIA kI nI bHu gHs bHs bB nB a0 a1 gB bE nE gE uT;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1] + ss[3] + ss[4];
	end;
	function outNF(ss)
		return ss[1] + ss[3] + ss[4];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:uT] = ss[1] + ss[3] + ss[4];
	end;

	# dUTs = dU + dUd + dUdB = sU - (cB * U * B) + (gUB * UB) + (gB * UB) + (cD * UB) - (cE * E * UdB)
end
