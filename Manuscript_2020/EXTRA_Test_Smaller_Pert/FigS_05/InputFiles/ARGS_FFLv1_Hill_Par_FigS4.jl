# Kinetic parameters
p = Dict([
    :g   => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 0.125,     # Y synthesis rate dependent of W (nM/min)
    :gY  => 0.1,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 0.0334,    # U maximum synthesis rate repressed by Y (nM/min)
    :kD  => 1,         # Y repression KD (nM)
    :nH  => 1,         # Hill coefficient for Y repression
    :gU  => 0.01,      # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mW  => 0.125,     # W synthesis rate dependent of U (1/min)
    :gW  => 0.01,      # W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mYs => NaN,       # LOCAL: Ys constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :gYs => NaN,       # LOCAL: Ys constitutive degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
]);