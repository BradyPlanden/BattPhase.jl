module BattPhase

    using LinearAlgebra
    export y₀!, Eqs11!, Jac!, Upwind!, Φ₊!, Φ̄₊, Ef!
    export Ef2!, Ef3!, YStore, PhiSwitch!, MidPred!

    include("custom_functions.jl")

end #module
