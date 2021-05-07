# Unfolded protein response (v05)
#   with explicit B & E mRNA dynamics + reporter
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
		dI   = mI - (gI * I) - (cBI * B * I) - (cA * (U + Ud + UdB) * I) + (gIB * IB) + (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dIB  = - (gI * IB) + (cBI * B * I) - (gIB * IB)
		dIa  = - (gI * Ia) + (cA * (U + Ud + UdB) * I) - (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dHu  = bHu - (gHs * Hu) - ((bHo * min(I,Hu)) + (bHs * min(Ia,Hu)))
		dHs  = - (gHs * Hs) + ((bHo * min(I,Hu)) + (bHs * min(Ia,Hu)))
		dBm  = (bBm * (1 + (nB * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gBm * Bm)
		dB   = (bB * Bm) - (gB * B) - (cB * U * B) + (gUB * UB) - (cB * Ud * B) + (gUB * UdB) + (gF * UB) + (gI * IB)
		dEm  = (bEm * (1 + (nE * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gEm * Em)
		dE   = (bE * Em) - (gE * E)
		dRm  = - (gHs * Rm) + ((bHo * min(I,Hu)) + (bHs * min(Ia,Hu)))
		dRs  = (0.008333333333333 * Rm) - (0.00003472222222222222 * Rs)
	end sU cB cBI gUB cD cE gB gF mI gI cA gIB gIA kI nI bHu gHs bHs bHo bBm nB a0 a1 gBm bB gB bEm nE gEm bE gE uT;
	# ODE system without feedback
	odeNF = @ode_def begin
		dU   = sU - (cB * U * B) + (gUB * UB) - (cD * U) + (cE * E * Ud) + (gB * UB)
		dUB  = (cB * U * B) - (gUB * UB) - (cD * UB) + (cE * E * UdB) - (gB * UB) - (gF * UB)
		dUd  = - (cB * Ud * B) + (gUB * UdB) + (cD * U) - (cE * E * Ud) + (gB * UdB)
		dUdB = (cB * Ud * B) - (gUB * UdB) + (cD * UB) - (cE * E * UdB) - (gB * UdB)
		dI   = mI - (gI * I) - (cBI * B * I) - (cA * uT * I) + (gIB * IB) + (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dIB  = - (gI * IB) + (cBI * B * I) - (gIB * IB)
		dIa  = - (gI * Ia) + (cA * uT * I) - (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dHu  = bHu - (gHs * Hu) - ((bHo * min(I,Hu)) + (bHs * min(Ia,Hu)))
		dHs  = - (gHs * Hs) + ((bHo * min(I,Hu)) + (bHs * min(Ia,Hu)))
		dBm  = (bBm * (1 + (nB * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gBm * Bm)
		dB   = (bB * Bm) - (gB * B) - (cB * U * B) + (gUB * UB) - (cB * Ud * B) + (gUB * UdB) + (gF * UB) + (gI * IB)
		dEm  = (bEm * (1 + (nE * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gEm * Em)
		dE   = (bE * Em) - (gE * E)
		dRm  = - (gHs * Rm) + ((bHo * min(I,Hu)) + (bHs * min(Ia,Hu)))
		dRs  = (0.008333333333333 * Rm) - (0.00003472222222222222 * Rs)
	end sU cB cBI gUB cD cE gB gF mI gI cA gIB gIA kI nI bHu gHs bHs bHo bBm nB a0 a1 gBm bB gB bEm nE gEm bE gE uT;

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
