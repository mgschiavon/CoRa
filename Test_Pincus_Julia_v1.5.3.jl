using DifferentialEquations
using Statistics
using DelimitedFiles
using Distributions
using Plots

cd("C:\\Users\\mgsch\\Dropbox (MGS-UCSF)\\MODEL - Quantifying feedback\\Paper\\CoRa\\")
iARG = (mm = "UPRv5",  # Label for motif file
     ex = "Ex01",      # Label for parameters file
     pp = :cD,         # Label for perturbation type
     ax = :cD,         # Label for condition/environment
     an = "ExSSs");    # Chose analysis type (Options: ExSSs, ExDyn, DYms, OptDY)


# Load functions & parameters:
include(string("InputFiles\\ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))	# Core parameters
mm = include(string("Library\\Md_",iARG.mm,".jl"));
fn = include(string("Library\\FN_DYs.jl"));
pO = copy(p);

# Find initial conditions as 0nM DTT steady state (i.e. pre-disturbance):
p[:cD] = 0.0;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
x0 = zeros(length(mm.odeFB.syms));
x0[5] = 256;                # :I, I total = 256 mol
x0[8] = p[:bHu]/p[:gHs];    # :Hu, basal Hac1 = 200 mol
x0[11] = 430000;            # :B, basal BiP = 430,000 mol
ss = try
        solve(ODEProblem(mm.odeFB,x0,1e6,pV),alg_hints=[:stiff],isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
    catch
        try
            solve(ODEProblem(mm.odeFB,x0,1e6,pV),isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
        catch
            println("WARNING: Error in steady state calculation.")
        end
    end;
x0 = ss.u[end];
# Choose variable to plot:
iS = 15;

# Simulate and plot different perturbation values:
### :cD - Enzymatic rate of disulfide bond breaking
### [0.0015 (1/mol)(1/s)] x DTT concentration [0 - 6500 mol = 0-5 mM]:
### [0 - 9710.475 (1/s)]
p[:cD] = 0.0;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = try
        solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff],isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
    catch
        try
            solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
        catch
            println("WARNING: Error in steady state calculation.")
        end
    end;
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," mM"),lw=3,legend=:topleft)

p[:cD] = 9710.475*0.66/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = try
        solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff],isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
    catch
        try
            solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
        catch
            println("WARNING: Error in steady state calculation.")
        end
    end;
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," mM"),lw=3)

p[:cD] = 9710.475*1/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = try
        solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff],isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
    catch
        try
            solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
        catch
            println("WARNING: Error in steady state calculation.")
        end
    end;
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," mM"),lw=3)

p[:cD] = 9710.475*1.5/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = try
        solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff],isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
    catch
        try
            solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
        catch
            println("WARNING: Error in steady state calculation.")
        end
    end;
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," mM"),lw=3)

p[:cD] = 9710.475*2.2/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = try
        solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff],isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
    catch
        try
            solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
        catch
            println("WARNING: Error in steady state calculation.")
        end
    end;
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," mM"),lw=3)

p[:cD] = 9710.475*3.3/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = try
        solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff],isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
    catch
        try
            solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
        catch
            println("WARNING: Error in steady state calculation.")
        end
    end;
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," mM"),lw=3)

p[:cD] = 9710.475;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = try
        solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff],isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
    catch
        try
            solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),isoutofdomain=(u,p,t) -> any(x -> (x < 0), u));
        catch
            println("WARNING: Error in steady state calculation.")
        end
    end;
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," mM"),lw=3)

    xlabel!("Minutes")
    ylabel!(string(mm.odeFB.syms[iS]," molecules"))
