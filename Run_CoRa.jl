## Running in julia terminal
	cd("C:\\Users\\mgsch\\Dropbox (MGS-UCSF)\\MODEL - Quantifying feedback\\Paper\\CoRa\\")
	using Pkg; Pkg.activate(".");
	iARG = (mm = "ATFv1",  # Label for motif file
       ex = "Ex01",      # Label for parameters file
       pp = :mY,         # Label for perturbation type
       ax = :mY,         # Label for condition/environment
       an = "OptDY");    # Chose analysis type (Options: ExSSs, ExDyn, DYms, OptDY)
	include("CoRa_Main_Julia_v_1_5_3.jl")
