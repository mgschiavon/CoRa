# Kinetic parameters
p = Dict([
    :g   => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 0.125,     # Y synthesis rate dependent of W (nM/min)
    :gY  => 0.1,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 2,         # U maximum synthesis rate repressed by Y (nM/min)
    :kD  => 1,         # Y repression KD (nM)
    :nH  => 1,         # Hill coefficient for Y repression
    :gU  => 0.0001,    # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :b   => 0.1565,    # U to Us transition rate (1/min)
    :bs  => 0.0001,    # Us to U transition rate (1/min)
    :mYs => NaN,       # LOCAL: Ys constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :gYs => NaN,       # LOCAL: Ys constitutive degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
]);