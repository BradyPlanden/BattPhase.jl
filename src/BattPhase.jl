module BattPhase

    using LinearAlgebra
    export y₀1!, y₀4!
    export Eqs11!, Jac!
    export Upwind!, WENO!, ENO!
    export Ef!, Ef2!, Ef3!
    export Φ₊!, Φ̄₊, YStore, PhiSwitch!, MidPred!
    export rk3solve, rk3asolve, pdrk3solve, pdrk3asolve, ps

    include("custom_functions.jl")
    include("PardisoSolve.jl")
    include("OpenBLASSolve.jl")

end #module
