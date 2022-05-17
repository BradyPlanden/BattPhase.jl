using BattPhase, HDF5, Infiltrator

NN = MM = Int(20)
δ = 0.1
tt = 0.
tf = 2.
StepNum = 60
Δₜ = (tf-tt)/StepNum
ν = ki₀ = 1.
ψ = 10
Ydata₁ = Array{Float64}(undef,ψ,StepNum+2,NN+2,MM+2)
for i ∈ 1:ψ
    κ = rand(0.1:0.00001:0.4)
    Ydata₁[i,:,:,:] = Semicircle(NN,MM,κ,δ,ν,ki₀,tt,tf,Δₜ)
end


# A = collect(reshape(1:120, 15, 8))
# h5write("/tmp/test2.h5", "mygroup2/A", A)