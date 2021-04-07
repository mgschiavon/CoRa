using DifferentialEquations
using Statistics
using DelimitedFiles
using Distributions
using Plots

cd("C:\\Users\\mgsch\\Dropbox (MGS-UCSF)\\MODEL - Quantifying feedback\\Paper\\CoRa\\")
iARG = (mm = "UPRv1",  # Label for motif file
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
x0 = zeros(length(mm.odeFB.syms));
x0[5] = 256;    # :I, I total = 256 mol
x0[8] = 200;    # :Hu, basal Hac1 = 200 mol
x0[10] = 430000;# :B, basal BiP = 430,000 mol
ss = solve(ODEProblem(mm.odeFB,x0,1e6,pV),alg_hint=[:stiff],reltol=1e-4,callback=TerminateSteadyState());
x0 = ss.u[end];

# Choose variable to plot:
iS = 9;

# Simulate and plot different perturbation values:
### :cD - Enzymatic rate of disulfide bond breaking
### [0.0015 (1/mol)(1/s)] x DTT concentration [0 - 6500 mol = 0-5 mM]:
### [0 - 9.75 (1/s)]
p[:cD] = 0.0;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,12000.0,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot(ss.t,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

p[:cD] = 9.75/2;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,12000.0,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

p[:cD] = 9.75;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,12000.0,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

    xlabel!("time")
    ylabel!(string(mm.odeFB.syms[iS]," molecules"))
