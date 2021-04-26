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
x0[5] = 256;    # :I, I total = 256 mol
x0[8] = p[:bHu]/p[:gHs];    # :Hu, basal Hac1 = 200 mol
x0[11] = 430000;# :B, basal BiP = 430,000 mol
ss = solve(ODEProblem(mm.odeFB,x0,1e6,pV),alg_hints=[:stiff]);
x0 = ss.u[end];
########## For UPRv5 | Ex01 ##########
#ss = solve(ODEProblem(mm.odeFB,x0,1e6,pV),alg_hints=[:stiff]);
#x0 = [28430.812810605, 372148.85954383,      0.0,      0.0,     18.133988230549,    237.36823020938,      0.6128894545444,    199.37654960166,      1.1036424751692,    247.5388675555,  73302.29924954,   1664.6746405859,      2994019.1377436,      1.1036424751692,    264.9378053966];
########## For UPRv5 | mutH ##########
#ss = solve(ODEProblem(mm.odeFB,x0,1e6,pV),alg_hints=[:stiff]);
#x0 = [31449.290697777767, 372148.8595413453,      0.0,      0.0,     19.89902676170431,    235.47213383742258,      0.743947411157102,    199.14055088050873,      1.3396411963220092,    243.62818590704651,  66266.80291661603,   1619.1904047976013,      2.9122129582690676e6,      1.3396411963220092,    321.59109966476416];
########## For UPRv5 | mutI ##########
#ss = solve(ODEProblem(mm.odeFB,x0,1e6,pV),alg_hints=[:stiff]);
#x0 = [9938.6135044051, 372148.85954332,      0.0,      0.0,    253.124484807358,      0.0,      2.99062310503571,    195.09491637774,      5.3852756990812,    323.50330569012, 209691.618214591,   2548.19949027465,      4.5830926081760e6,      5.3852756990812,   1292.77655884681];
# Choose variable to plot:
iS = 15;

# Simulate and plot different perturbation values:
### :cD - Enzymatic rate of disulfide bond breaking
### [0.0015 (1/mol)(1/s)] x DTT concentration [0 - 6500 mol = 0-5 mM]:
### [0 - 9710.475 (1/s)]
p[:cD] = 0.0;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," nM"),lw=3,legend=:topleft)

p[:cD] = 9710.475*0.66/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," nM"),lw=3)

p[:cD] = 9710.475*1/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," nM"),lw=3)

p[:cD] = 9710.475*1.5/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff]);
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," nM"),lw=3)

p[:cD] = 9710.475*2.2/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
#ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff]);
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV));
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," nM"),lw=3)

p[:cD] = 9710.475*3.3/5;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
#ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff]);
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV));
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," nM"),lw=3)

p[:cD] = 9710.475;
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.odeFB.sys.ps];
#ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV),alg_hints=[:stiff]);
ss = solve(ODEProblem(mm.odeFB,x0,240.0*60,pV));
x = zeros(length(ss.u));
for i in 1:length(ss.u)
    x[i] = ss.u[i][iS];
end
plot!(ss.t/60,x,label=string("DTT = ",round(100*p[:cD]*5.0/9710.475)/100," nM"),lw=3)

    xlabel!("Minutes")
    ylabel!(string(mm.odeFB.syms[iS]," molecules"))