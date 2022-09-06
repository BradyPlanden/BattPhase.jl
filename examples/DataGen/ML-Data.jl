using BattPhase, HDF5, OrderedCollections, FileIO, Infiltrator

function Data()
    PhaseOut = κout = tuple()
    ξ = 40 # Physical Spacial Range
    ti = 0.
    tf = 1.
    StepNum = 250
    Δₜ = (tf-ti)/StepNum
    ν = ki₀ = 1.
    ψ = 16
    ζ = [50, 25]
    Z = ζ.^2
    δ = rand(0.01:1e-5:2,ψ)

    for ρ ∈ ["train" "valid" "test"]
        #PhaseFieldData = tuple()
        for η ∈ ζ
            MM = η
            NN = η #5 #Only 5 horizontal points currently
            γ = (NN)*(MM)
            LoopStepNum = StepNum+3
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
        g["alpha"] = δ #zeros(ψ)
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

    end

    # ω = OrderedDict(
    #     "PhaseFieldData_$(ζ[1])"=>PhaseFieldData[1],
    #     "PhaseFieldData_$(ζ[2])"=>PhaseFieldData[2],
    #     "PhaseFieldData_$(ζ[3])"=>PhaseFieldData[3],
    #     "PhaseFieldData_$(ζ[4])"=>PhaseFieldData[4],
    #     )
    # save("CE_train_$(Int(tf*60)).h5", ω)

    return PhaseOut, δ
end

Phase, δ = Data()


#Dict("train"=>Dict{String,Any}("Phase1"=>[2,3,5],"Phase3"=>rand(3)))
