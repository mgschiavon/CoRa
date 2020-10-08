# Steady state & DY calculation functions

# Julia v.1.1.1

module fn
	# Required libraries
	using DifferentialEquations
	using Statistics

	# Steady state function for a given system
	# INPUT: syst - Handle for the ODE system (@ode_def)
	#        p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	#        rtol - Tolerance value for ODE solver
	#        uns  - 1 to use a slower, more stable ODE solver
	# OUPUT: ss   - Vector of steady state of the ODE system
	function SS(syst, p, x0, rtol, uns)
		pV = [p[i] for i in syst.params];
		if (uns == 1)
			s0 = x0;
			t = 0.
			while (any([abs(x) for x in syst(s0, pV, 0.)] .> (rtol * s0)) && t<1e6)
				ss = solve(ODEProblem(syst,s0,10000.,pV),AutoTsit5(Rosenbrock23()),reltol=rtol);
				t += 100.;
				s0 = last(ss.u);
				#println("Time: ",t,", ",s0)
				end;
			return s0
		else
			ss = solve(SteadyStateProblem(syst, x0, pV), SSRootfind());
			# Verify this is the steady state:
			if any(ss.u .< 0)
				ss = solve(SteadyStateProblem(syst, x0, pV), DynamicSS(Rodas5(); reltol=rtol));
				return ss.u
			end
			if any([abs(x) for x in syst(ss.u, pV, 0.)] .> (rtol * ss.u))
				ss = solve(SteadyStateProblem(syst, ss.u, pV), DynamicSS(Rodas5(); reltol=rtol));
				return ss.u
			end
			return ss.u
			end
	end;

	# ODE dynamics for a given system
	# INPUT: syst - Handle for the ODE system (@ode_def)
	#        p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	#        tspan- Time to simulate
	# OUPUT: xD   - Vector of steady state of the ODE system
	function Dyn(syst, p, x0, tspan)
		pV = [p[i] for i in syst.params];
		xD = solve(ODEProblem(syst,x0,tspan,pV),AutoTsit5(Rosenbrock23()),reltol=1e-6);
		return xD
	end;

	# DY function
	# INPUT: Y    - Output of full system before perturbation
	#        YD   - Output of full system after perturbation
	#        Ynf  - Output of non-feedback system before perturbation (Ynf:=Y)
	#        YnfD - Output of non-feedback system after perturbation
	# OUPUT:      - DY value
	function DY(Y,YD,Ynf,YnfD)
		if abs(log10(YnfD/Ynf)) < 1e-4
			return NaN
		end
		return log10(YD/Y)/log10(YnfD/Ynf);
	end;

	# DY curve
	# INPUT: p     - Dictionary function with the ODE parameters & values
	#        pert  - Handle for the perturbation details
	#        motif - Handle for the considered motif
	#        uns  - 1 to use a slower, more stable ODE solver
	# OUPUT: DYs   - Vector of DY values for the range of parameters
	function DYc(p, pert, motif, uns)
		r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
		DYs = Array{Float64}(undef,length(r));
		DYs /= 0;
		p[pert.p] = pert.c;
		for i in 1:length(r)
			p[pert.p] *= r[i];
			rtol = 1e-6;
			flg1 = 1;
			ssR = ones(length(motif.odeFB.syms));
			soR = ones(length(motif.odeNF.syms));
			while(rtol >= 1e-24)
				# Reference steady state:
				ssR = SS(motif.odeFB, p, ssR, rtol, uns);
				# Locally analogous system reference steady state:
				motif.localNF(p,ssR);
				soR = SS(motif.odeNF, p, soR, rtol, uns);
				if(abs(motif.outFB(ssR) - motif.outNF(soR)) > 1e-4)
					rtol *= 1e-3;
					if(rtol < 1e-24)
						println("ERROR: Check NF system (reltol=",rtol*1e3,").")
						println(vcat(pert.p,i,[p[i] for i in motif.odeFB.params],motif.outFB(ssR),motif.outNF(soR)))
						#throw(DomainError("x-("))
						if(abs(motif.outFB(ssR) - motif.outNF(soR))/motif.outFB(ssR) > 0.01)
							flg1 = 0;
							println("SS results excluded!")
						end
					end
				else
					break
				end
			end
			# Perturbation:
			p[pert.p] *= pert.d;
			if(flg1==1)
				ssD = SS(motif.odeFB, p, ssR, rtol, uns);
				soD = SS(motif.odeNF, p, soR, rtol, uns);
				DYs[i] = DY(motif.outFB(ssR), motif.outFB(ssD), motif.outNF(soR), motif.outNF(soD));
			end
			p[pert.p] /= pert.d;
			p[pert.p] /= r[i];
		end
		return DYs
	end;

	# DY "metrics"
	# INPUT: DYs   - Vector of DY values for the range of parameters
	#        pert  - Handle for the perturbation details
	# OUPUT:       - [Range of DY<=eps,start,end,min(DY),optimal rho]
	function DYm(DYs, pert)
		r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
		i = DYs .<= pert.eps;
		j = findall(i);
		x = copy(DYs);
		x[x .=== NaN] .= Inf;
		if isempty(j)
			return [sum(DYs[i])./length(DYs[i]), NaN, NaN, minimum(x), r[argmin(x)]]
		end
		return [sum(DYs[i])./length(DYs[i]), pert.c * r[j[1]], pert.c * r[j[end]], minimum(x), r[argmin(x)]]
	end;

	# SSs for "optimal" control
	# INPUT: p     - Dictionary function with the ODE parameters & values
	#        pert  - Handle for the perturbation details
	#        motif - Handle for the considered motif
	#        DYs   - Vector of DY values for the range of parameters
	#        uns  - 1 to use a slower, more stable ODE solver
	# OUPUT:       - [Optimal rho, steady state for the full system before & after perturbation, and for the non-feedback system before & after perturbation]
	function SSopt(p, pert, motif, DYs, uns)
		r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
		x = copy(DYs);
		x[x .=== NaN] .= Inf;
		p[pert.p] = pert.c;
		p[pert.p] *= r[argmin(x)];
			rtol = 1e-6;
			ssR = zeros(length(motif.odeFB.syms));
			soR = zeros(length(motif.odeNF.syms));
			while(rtol >= 1e-24)
				# Reference steady state:
				ssR = SS(motif.odeFB, p, ssR, rtol, uns);
				# Locally analogous system reference steady state:
				motif.localNF(p,ssR);
				soR = SS(motif.odeNF, p, soR, rtol, uns);
				if(abs(motif.outFB(ssR) - motif.outNF(soR)) > 1e-4)
					rtol *= 1e-3;
					if(rtol < 1e-24)
						println("ERROR: Check NF system (reltol=",rtol*1e3,").")
						println(vcat(pert.p,i,[p[i] for i in motif.odeFB.params],motif.outFB(ssR),motif.outNF(soR)))
						#throw(DomainError("x-("))
					end
				else
					break
				end
			end
			# Perturbation:
			p[pert.p] *= pert.d;
			ssD = SS(motif.odeFB, p, ssR, rtol, uns);
			soD = SS(motif.odeNF, p, soR, rtol, uns);
		p[pert.p] /= pert.d;
		p[pert.p] /= r[argmin(x)];
		return [vcat(p[pert.p] * r[argmin(x)],ssR,ssD,soR,soD)]
	end;
end