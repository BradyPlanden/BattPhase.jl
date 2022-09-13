using BattPhase, HDF5, OrderedCollections, FileIO, Infiltrator

function Data(ρ, StepNum, ti, tf, ψ, ζ, ξ, ν, ki₀)
    PhaseOut = κout = tuple()
    Δₜ = (tf-ti)/StepNum
    Z = ζ.^2
    δ = rand(0.1:1e-3:2,ψ)

    for η ∈ ζ
        MM = η
        NN = η
        γ = (NN)*(MM)
        LoopStepNum = StepNum
        Ydata₁ = Array{Float64}(undef,ψ,LoopStepNum,NN*MM)
        κ₁ = Vector{Float64}(undef,ψ)

        for i ∈ 1:ψ
            κ = rand(0.1:1e-5:1)
            Ydata₁[i,:,:], κ₁[i] = Seed1D(NN,MM,κ,δ[i],ν,ki₀,ti,tf,Δₜ)
        end

        @show size(Ydata₁)

        PhaseOut = flatten!(PhaseOut,permutedims(Ydata₁,[3,2,1])) # Align the output .h5 for python 
        κout = flatten!(κout,κ₁)
        
    end


    # Training Dataset 
    f = h5open("CE_$(ρ)_E1.h5","w")  #"r+"? / CE_train_$(Int(tf*60)).h5
    create_group(f, "$ρ")
    g = f["$ρ"] 
    g["beta"] = zeros(ψ)
    g["gamma"] = zeros(ψ)
    g["alpha"] = δ 
    for i ∈ 1:length(ζ)
        ϵ = Vector(range(0,ξ,step=(ξ/(Z[i]-1)))) 
        g["pde_$(StepNum)-$(Z[i])"] = PhaseOut[i]
        attributes(g["pde_$(StepNum)-$(Z[i])"])["dt"] = Δₜ
        attributes(g["pde_$(StepNum)-$(Z[i])"])["dx"] = ξ/(Z[i]-1)
        attributes(g["pde_$(StepNum)-$(Z[i])"])["nt"] = StepNum
        attributes(g["pde_$(StepNum)-$(Z[i])"])["nx"] = ζ[i]
        attributes(g["pde_$(StepNum)-$(Z[i])"])["tmax"] = tf
        attributes(g["pde_$(StepNum)-$(Z[i])"])["tmin"] = ti
        attributes(g["pde_$(StepNum)-$(Z[i])"])["x"] = ϵ
    end
    close(f)

return PhaseOut, δ

end

ξ = 40 # Physical Spacial Range
ti = 0.
tf = 1.
StepNum = 100
ν = ki₀ = 1.
ψ = [16,8,8]#[1024,256,256]
ζ = [86, 43]

ρ = ["train" "valid" "test"]

for i ∈ 1:length(ρ)
    Phase, δ = Data(ρ[i], StepNum, ti, tf, ψ[i], ζ, ξ, ν, ki₀)
end

#Dict("train"=>Dict{String,Any}("Phase1"=>[2,3,5],"Phase3"=>rand(3)))
