using BattPhase, HDF5, OrderedCollections, FileIO, Infiltrator

function Data()

    PhaseFieldData = κ = tuple()
    δ = 0.1
    ti = 0.
    tf = 1.
    StepNum = 30
    Δₜ = (tf-ti)/StepNum
    ν = ki₀ = 1.
    ψ = 300
    ζ = [10, 22]

    for α ∈ ["train" "valid" "test"]

        for η ∈ ζ
            NN = MM = η
            Ydata₁ = Array{Float64}(undef,ψ,StepNum+3,NN+2,MM+2)
            κ₁ = Vector{Float64}(undef,ψ)

            for i ∈ 1:ψ
                κ = rand(0.1:1e-5:0.4)
                Ydata₁[i,:,:,:], κ₁[i] = Semicircle(NN,MM,κ,δ,ν,ki₀,ti,tf,Δₜ)
            end

            PhaseFieldData = flatten!(PhaseFieldData, permutedims(Ydata₁,[3,4,2,1])) # Align the output .h5 for python 
            κ = flatten!(κ,κ₁)
        end

        # Training Dataset 
        f = h5open("CE_$(α)_E1.h5","w")  #"r+"? / CE_train_$(Int(tf*60)).h5
        create_group(f, "$α")
        g = f["$α"]
        for i ∈ 1:length(ζ)
        g["pde_$(StepNum)-$(ζ[i])"] = PhaseFieldData[i]
        attributes(g["pde_$(StepNum)-$(ζ[i])"])["dt"] = Δₜ
        attributes(g["pde_$(StepNum)-$(ζ[i])"])["dx"] = 1/ζ[i]
        attributes(g["pde_$(StepNum)-$(ζ[i])"])["nt"] = StepNum
        attributes(g["pde_$(StepNum)-$(ζ[i])"])["nx"] = ζ[i]
        attributes(g["pde_$(StepNum)-$(ζ[i])"])["tmax"] = tf
        attributes(g["pde_$(StepNum)-$(ζ[i])"])["tmin"] = ti
        attributes(g["pde_$(StepNum)-$(ζ[i])"])["x"] = Vector(range(0,ζ[i])) 
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
end

Phase = Data()


#Dict("train"=>Dict{String,Any}("Phase1"=>[2,3,5],"Phase3"=>rand(3)))
