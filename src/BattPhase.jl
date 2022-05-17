module BattPhase

    using LinearAlgebra, SparseArrays
    export y₀1!, y₀4!
    export Eqs11!, Jac!
    export Upwind!, WENO!, ENO!
    export Ef!, Ef2!, Ef3!
    export Φ₊!, Φ̄₊, YStore, PhiSwitch!, MidPred!
    export rk3solve, rk3asolve, pdrk3solve, pdrk3asolve, ps

    include("custom_functions.jl")
    include("PardisoSolve.jl")
    include("OpenBLASSolve.jl")


"""
    flatten_(a::Tuple, b...) 
    Flattens input Tuple "a" and inserts "b" 
"""
function flatten! end
flatten!() = ()
flatten!(a::Tuple) = Tuple(a)
flatten!(a) = (a,)
flatten!(a::Tuple, b...) = tuple(a..., flatten!(b...)...)
flatten!(a, b...) = tuple(a, flatten!(b...)...)
flatten_tuple(x::Tuple) = flatten!(x...)


end #module
