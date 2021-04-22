using DifferentialEquations
using Statistics
using DelimitedFiles
using Distributions
using Plots

cd("C:\\Users\\mgsch\\Dropbox (MGS-UCSF)\\MODEL - Quantifying feedback\\Paper\\CoRa\\")
iARG = (mm = "UPRv5",  # Label for motif file
     ex = "mutH",      # Label for parameters file
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
pV = [p[i] for i in mm.odeFB.params];
x0 = zeros(length(mm.odeFB.syms));
x0[5] = 256;    # :I, I total = 256 mol
x0[8] = p[:bHu]/p[:gHs];    # :Hu, basal Hac1 = 200 mol
x0[11] = 430000;# :B, basal BiP = 430,000 mol
#ss = solve(ODEProblem(mm.odeFB,x0,1e7,pV),alg_hint=[:stiff],reltol=1e-5,callback=TerminateSteadyState());
#x0 = ss.u[end];
########## For UPRv5 | Ex01 ##########
#ss = solve(ODEProblem(mm.odeFB,x0,1e6,pV),alg_hint=[:stiff]);
#x0 = [28430.812810605, 372148.85954383,      0.0,      0.0,     18.133988230549,    237.36823020938,      0.6128894545444,    199.37654960166,      1.1036424751692,    247.5388675555,  73302.29924954,   1664.6746405859,      2994019.1377436,      1.1036424751692,    264.9378053966];
########## For UPRv5 | mutH ##########
#ss = solve(ODEProblem(mm.odeFB,x0,1e6,pV),alg_hint=[:stiff]);
x0 = [31449.290697777767
 372148.8595413453
      0.0
      0.0
     19.89902676170431
    235.47213383742258
      0.743947411157102
    199.14055088050873
      1.3396411963220092
    243.62818590704651
  66266.80291661603
   1619.1904047976013
      2.9122129582690676e6
      1.3396411963220092
    321.59109966476416];

# Choose variable to plot:
iS = 15;

# Simulate and plot different perturbation values:
### :cD - Enzymatic rate of disulfide bond breaking
### [0.0015 (1/mol)(1/s)] x DTT concentration [0 - 6500 mol = 0-5 mM]:
### [0 - 9710.475 (1/s)]
p[:cD] = 0.0;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9710.475," nM"),lw=3)

p[:cD] = 9710.475*0.66/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9710.475," nM"),lw=3)

p[:cD] = 9710.475*1/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9710.475," nM"),lw=3)

p[:cD] = 9710.475*1.5/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9710.475," nM"),lw=3)

p[:cD] = 9710.475*2.2/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9710.475," nM"),lw=3)

p[:cD] = 9710.475*3.3/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9710.475," nM"),lw=3)

p[:cD] = 9710.475;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9710.475," nM"),lw=3)

    xlabel!("Minutes")
    ylabel!(string(mm.odeFB.syms[iS]," molecules"))
