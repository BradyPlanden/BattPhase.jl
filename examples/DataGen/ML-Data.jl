using BattPhase, HDF5, Infiltrator

NN = MM = Int(100)
δ = 0.1
tt = 0.
tf = 2.
StepNum = 60
Δₜ = (tf-tt)/StepNum
ν = ki₀ = 1.
ψ = 500
Ydata₁ = Array{Float64}(undef,ψ,StepNum+2,NN+2,MM+2)
κ₁ = Vector{Float64}(undef,ψ)
for i ∈ 1:ψ
    κ = rand(0.1:1e-5:0.4)
    Ydata₁[i,:,:,:], κ₁[i] = Semicircle(NN,MM,κ,δ,ν,ki₀,tt,tf,Δₜ)
end


save("pde_60_100.h5", OrderedDict("Ydata₁"=>Ydata₁,"κ₁"=>κ₁))