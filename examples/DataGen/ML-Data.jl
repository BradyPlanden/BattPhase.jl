using BattPhase, HDF5

NN = MM = Int(320)
δ = 0.1
tt = 0.
tf = 2.
ν = ki₀ = 1.
Ydata = tuple()

for i ∈ 1:100
    κ = randn()
    Ydata₀ = Semicircle(NN,MM,κ,δ,ν,ki0,tt,tf)
    Ydata = flatten!(Ydata,Ydata₀)
end


A = collect(reshape(1:120, 15, 8))
h5write("/tmp/test2.h5", "mygroup2/A", A)