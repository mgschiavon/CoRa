# Kinetic parameters
p = Dict([
    :g   => 0.01,      # Cell-death rate (1/reproduction rate)
    :kY  => 0.5,       # Propensity of X to differentiate to Y
    :kX  => 1,         # Propensity of Y to self-renew
    :kZ  => 0.5,       # Propensity of X to differentiate to Z when Y = oY
    :n   => 1,         # Feedback strength (i.e. nonlinearity)
    :oY  => 100,       # "Optimal" Y number (tissue size)
    :mYs => NaN,       # LOCAL: Ys constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :gYs => NaN,       # LOCAL: Ys constitutive degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
]);
