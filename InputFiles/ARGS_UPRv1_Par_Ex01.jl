# Kinetic parameters
p = Dict([
    :sU  => 310.0,     # Source rate for protein unfolding [mol/s]
    :cB  => 0.0350,    # Attachment rate of single BiP molecule [(1/(mol^2)(1/s)]
    :gUB => 196.0,     # Dissociation rate of BiP from folding complex [(1/mol)(1/s)]
    :cD  => 0.0015,    # Enzymatic rate of disulfide bond breaking [(1/(mol^2)(1/s)] x DTT concentration [mol]
    :cE  => 0.0015,    # Enzymatic rate of disulfide bond formation [(1/(mol^2)(1/s)]
    :gB  => 0.000139,  # Decay rate of BiP [(1/mol)(1/s)]
    :gF  => 0.000833,  # Folding rate of protein in the folding complex [(1/mol)(1/s)]
    :cA  => 0.01,      #
    :gIB => 0.01,      #
    :gIA => 0.01,      #
    :kI  => 0.01,      #
    :nI  => 0.01,      #
    :bHu => 0.01,      #
    :gHs => 0.01,      #
    :bHs => 0.01,      #
    :bB  => 0.01,      #
    :nB  => 0.01,      #
    :a0  => 0.01,      #
    :a1  => 0.01,      #
    :gB  => 0.01,      #
    :bE  => 0.01,      #
    :nE  => 0.01,      #
    :gE  => 0.01,      #
    :uT  => NaN,       # LOCAL: Constitutive I activation "boost" (to get the locally anologous system)
]);
