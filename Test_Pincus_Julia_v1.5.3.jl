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
x0[1] = 0;                  # :U, unfolded proteins
x0[2] = 1000;               # :UB, unfolded proteins:BiP complex
x0[5] = 256;                # :I, Ire1
x0[8] = 190;                # :Hu, Hac1 unspliced mRNA
x0[10] = 12.9287;           # :Bm, BiP mRNA
x0[11] = 75000;             # :B, BiP
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
DTT = [0,0.6600,0.9900,1.4800,2.2200,3.5500,5.0000];    # DTT concentration
### :cD - Enzymatic rate of disulfide bond breaking x DTT number of molecules in the ER
for d in DTT
    p[:cD] = p[:cD0] * (p[:mMc] * (d * p[:ERv]));
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
    if d==0
        plot(ss.t/60,x,label=string("DTT = ",d," mM"),lw=3,legend=:topleft)
    else
        plot!(ss.t/60,x,label=string("DTT = ",d," mM"),lw=3,legend=:topleft)
    end
end

    xlabel!("Minutes")
    ylabel!(string(mm.odeFB.syms[iS]," molecules"))
