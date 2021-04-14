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
#x0 = [28430.69339981689,372148.7146545400,0.0,0.0,18.1339171204673,237.3682039471226,0.6128844770696,199.3765691647196,1.103622912111158,247.5385343531709,73302.5785928124,1664.670765186364,0.00000299];  #For UPRv3 & UPRv4
ss = solve(ODEProblem(mm.odeFB,x0,1e6,pV),alg_hint=[:stiff]);
x0 = [28430.81281
 372148.85954
      0.0
      0.0
     18.13398
    237.36823
      0.61288
    199.37654
      1.10364
    247.53886
  73302.29924
   1664.67464
      0.0000029
      1.103642475169266
    264.937805396655];

# Choose variable to plot:
iS = 15;

# Simulate and plot different perturbation values:
### :cD - Enzymatic rate of disulfide bond breaking
### [0.0015 (1/mol)(1/s)] x DTT concentration [0 - 6500 mol = 0-5 mM]:
### [0 - 9.75 (1/s)]
p[:cD] = 0.0;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

p[:cD] = 9.75*0.66/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

p[:cD] = 9.75*1/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

p[:cD] = 9.75*1.5/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

p[:cD] = 9.75*2.2/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

p[:cD] = 9.75*3.3/5;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

p[:cD] = 9.75;
pV = [p[i] for i in mm.odeFB.params];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hint=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",p[:cD]*5.0/9.75," nM"),lw=3)

    xlabel!("Minutes")
    ylabel!(string(mm.odeFB.syms[iS]," molecules"))
