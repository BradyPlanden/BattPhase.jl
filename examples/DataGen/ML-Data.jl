using BattPhase, HDF5, OrderedCollections, FileIO, Infiltrator

function Data()

    PhaseFieldData = κ = tuple()
    δ = 0.1
    tt = 0.
    tf = 1.
    StepNum = 30
    Δₜ = (tf-tt)/StepNum
    ν = ki₀ = 1.
    ψ = 300
    ζ = [10, 20, 40, 80]

    for η ∈ ζ
        NN = MM = η
        Ydata₁ = Array{Float64}(undef,ψ,StepNum+3,NN+2,MM+2)
        κ₁ = Vector{Float64}(undef,ψ)

        for i ∈ 1:ψ
            κ = rand(0.1:1e-5:0.4)
            Ydata₁[i,:,:,:], κ₁[i] = Semicircle(NN,MM,κ,δ,ν,ki₀,tt,tf,Δₜ)
        end

        PhaseFieldData = flatten!(PhaseFieldData, Ydata₁)
        κ = flatten!(κ,κ₁)
    end

    ω = OrderedDict(
        "PhaseFieldData_$(ζ[1])"=>PhaseFieldData[1],
        "PhaseFieldData_$(ζ[2])"=>PhaseFieldData[2],
        "PhaseFieldData_$(ζ[3])"=>PhaseFieldData[3],
        "PhaseFieldData_$(ζ[4])"=>PhaseFieldData[4],
        )
    save("CE_train_$(Int(tf*60)).h5", ω)
end

Phase = Data()


#Dict("train"=>Dict{String,Any}("Phase1"=>[2,3,5],"Phase3"=>rand(3)))
