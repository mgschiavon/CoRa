# Kinetic parameters
p = Dict([
    :sU  => 310.0,     # Source rate for protein unfolding [mol/s]
    :cB  => 0.0350,    # Attachment rate of single BiP molecule to unfolded proteins [(1/mol)(1/s)]
    :cBI => 0.0350,    # Attachment rate of single BiP molecule to Ire1 [(1/mol)(1/s)]
    :gUB => 196.0,     # Dissociation rate of BiP from folding complex [1/s]
    :cD  => 9710.475,  # Enzymatic rate of disulfide bond breaking [0.0015 (1/mol)(1/s)] x DTT concentration [0 - 6473650 mol = 0-5 mM]: [0 - 9710.475 (1/s)]
    :cE  => 0.0015,    # Enzymatic rate of disulfide bond formation [(1/mol)(1/s)]
    :gB  => 0.000139,  # Decay rate of BiP [1/s]
    :gF  => 0.000833,  # Folding rate of protein in the folding complex [1/s]
    :mI  => 0.0356,    # Ire1 synthesis rate to keep the population of Ire1 approximately 256 mol [mol/s]
    :gI  => 0.000139,  # Decay rate of Ire1 (assumed the same as BiP and Ero1) [1/s]
    :cA  => 0.000233,  # Attachment rate of Ire1 molecule to an unfolded protein [(1/mol)(1/s)]
    :gIB => 196.0,     # Dissociation rate of BiP from the inactive complex [1/s]
    :gIA => 196.0,     # Part of non-linear (active Ire1 cooperative) decay rate [1/s]
    :kI  => 45.0,      # Part of non-linear (active Ire1 cooperative) decay rate [mol]
    :nI  => 4.5,       # Part of non-linear (active Ire1 cooperative) decay rate
    :bHu => 0.167,     # Transcription rate of Hac1 unspliced mRNA [mol/s]
    :gHs => 0.000833,  # Decay rate of Hac1 mRNA [1/s]
    :bHs => 0.0015,    # Splicing rate of Hac1 mRNA [1/s]
    :bBm => 0.1625,    # BiP mRNA basal synthesis rate [mol/s]
    :nB  => 4.0,       # Part of BiP synthesis rate function
    :a0  => 296.5,     # Part of BiP synthesis rate function
    :a1  => 5.26,      # Part of BiP synthesis rate function
    :gBm => 0.000667,  # Decay rate of BiP mRNA (1/s)
    :bB  => 0.25,      # BiP synthesis rate [1/s]
    :gB  => 0.000139,  # Decay rate of BiP (1/s)
    :bEm => 1.08,      # Ero1 mRNA basal synthesis rate [mol/s]
    :nE  => 7.0,       # Part of Ero1 synthesis rate function
    :gEm => 0.000667,  # Decay rate of Ero1 mRNA [1/s]
    :bE  => 0.25,      # Ero1 basal synthesis rate [1/s]
    :gE  => 0.000139,  # Decay rate of Ero1 [1/s]
    :uT  => NaN,       # LOCAL: Constitutive I activation "boost" (to get the locally anologous system)
]);
