# Kinetic parameters
p = Dict([
    :sU  => 310.0,     # Source rate for protein unfolding [mol/s]
    :cB  => 0.0327,    # Attachment rate of single BiP molecule to unfolded proteins [(1/mol)(1/s)]
    :cBI => 0.0327,    # Attachment rate of single BiP molecule to Ire1 [(1/mol)(1/s)]
    :gUB => 195.9896,  # Dissociation rate of BiP from folding complex [1/s]
    #:cD0 => 0.001516231919278, # Enzymatic rate of disulfide bond breaking [(1/mol)(1/s)]
    #:mMc => 6.022e+05,         # To convert DTT concetration (mM) to expected number of molecules
    #:ERv => 2.15434934,        # ER volume
    :cD  => 6000,      # Enzymatic rate of disulfide bond breaking x DTT number of molecules in the ER: (DTT [mM])*(cD0*mMc*ERv); 3000 ~ 1.53mM DTT
    :cE  => 0.0015,    # Enzymatic rate of disulfide bond formation [(1/mol)(1/s)]
    :gF  => 0.00083333,# Folding rate of protein in the folding complex [1/s]
    :mI  => 0.0356,    # Ire1 synthesis rate to keep the population of Ire1 approximately 256 mol [mol/s]
    :gI  => 0.000139,  # Decay rate of Ire1 (assumed the same as BiP and Ero1) [1/s]
    :cA  => 0.00021777,# Attachment rate of Ire1 molecule to an unfolded protein [(1/mol)(1/s)]
    :gIB => 195.9896,  # Dissociation rate of BiP from the inactive complex [1/s]
    :gIA => 195.9896,  # Part of non-linear (active Ire1 cooperative) decay rate [1/s]
    :kI  => 45.0,      # Part of non-linear (active Ire1 cooperative) decay rate [mol]
    :nI  => 4.5,       # Part of non-linear (active Ire1 cooperative) decay rate
    :bHu => 0.1667,    # Transcription rate of Hac1 unspliced mRNA [mol/s]
    :gHs => 0.00083333,# Decay rate of Hac1 mRNA [1/s]
    :bHs => 0.0015,    # Splicing rate of Hac1 mRNA [1/s]
    :bHo => 1.5152e-9, # Basal splicing rate of Hac1 mRNA [1/s]
    :bBm => 0.1625,    # BiP mRNA basal synthesis rate [mol/s]
    :nB  => 3.0,       # Part of BiP synthesis rate function
    :a0  => 296.5118,  # Part of BiP synthesis rate function
    :a1  => 5.2558,    # Part of BiP synthesis rate function
    :gBm => 0.00066667,# Decay rate of BiP mRNA (1/s)
    :bB  => 0.25,      # BiP synthesis rate [1/s]
    :gB  => 0.00013889,# Decay rate of BiP [1/s]
    :bEm => 1.0772,    # Ero1 mRNA basal synthesis rate [mol/s]
    :nE  => 6.0,       # Part of Ero1 synthesis rate function
    :gEm => 0.00066667,# Decay rate of Ero1 mRNA [1/s]
    :bE  => 0.25,      # Ero1 basal synthesis rate [1/s]
    :gE  => 0.00013889,# Decay rate of Ero1 [1/s]
    :uT  => NaN,       # LOCAL: Constitutive I activation "boost" (to get the locally anologous system)
]);

# Initial conditions
x0 = zeros(length(mm.odeFB.syms));
x0[1] = 0;                  # :U, unfolded proteins
x0[2] = 1000;               # :UB, unfolded proteins:BiP complex
x0[5] = 256;                # :I, Ire1
x0[8] = 190;                # :Hu, Hac1 unspliced mRNA
x0[10] = 12.9287;           # :Bm, BiP mRNA
x0[11] = 75000;             # :B, BiP
