# Kinetic parameters
p = Dict([
    :g   => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 0.125,     # Y synthesis rate dependent of W (nM/min)
    :gY  => 0.1,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mA  => 0.372,     # A maximum synthesis rate repressed by Y (nM/min)
    :kD  => 1,         # Y repression KD (nM)
    :nH  => 1,         # Hill coefficient for Y repression
    :mB  => 0.125,     # B constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :mU  => 0.1,       # U synthesis rate dependent of Y (nM/min)
    :e0  => 0.0001,    # U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
    :eP  => 0.05,      # U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    :bA  => 0.5,       # U to Us transition rate (1/(nM min))
    :bI  => 0.5,       # U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    :mYs => NaN,       # LOCAL: Ys constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :gYs => NaN,       # LOCAL: Ys constitutive degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
]);