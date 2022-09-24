using BattPhase, HDF5, OrderedCollections, FileIO, Infiltrator, LatinHypercubeSampling

function Data(ρ, StepNum, ti, tf, ψ, δ, κ, ζ, ξ, ν, ki₀)
    PhaseOut = κout = tuple()
    Δₜ = (tf-ti)/StepNum
    Z = ζ.^2
    #δ = rand(0.1:1e-3:2,ψ)

    for η ∈ ζ
        MM = η
        NN = η
        γ = (NN)*(MM)
        LoopStepNum = StepNum
        Ydata₁ = Array{Float64}(undef,ψ,LoopStepNum,NN*MM)
        κ₁ = Vector{Float64}(undef,ψ)

        for i ∈ 1:ψ
            #κ = rand(0.1:1e-5:1)
            #@infiltrate cond=true
            Ydata₁[i,:,:], κ₁[i] = Seed1D(NN,MM,κ[i],δ[i],ν[i],ki₀[i],ti,tf,Δₜ) #Seed1D(NN,MM,κ[i],δ[i],ν[i],ki₀[i],ti,tf,Δₜ)
        end

        @show size(Ydata₁)

        PhaseOut = flatten!(PhaseOut,permutedims(Ydata₁,[3,2,1])) # Align the output .h5 for python 
        κout = flatten!(κout,κ₁)
        
    end


    # Training Dataset 
    f = h5open("CE_$(ρ)_E1.h5","w")  #"r+"? / CE_train_$(Int(tf*60)).h5
    create_group(f, "$ρ")
    g = f["$ρ"] 
    g["beta"] = ki₀
    g["gamma"] = ν
    g["alpha"] = δ
    g["kappa"] = κ 
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

return PhaseOut

end

function LatinHyper()
    Output = δ_out = κ_out = ν_out = ki₀_out = tuple()
    ξ = 40 # Physical Spacial Range
    ti = 0.
    tf = 1/3
    StepNum = 100
    #ν = ki₀ = 1.
    ψ = [256,128,128]#[1024,128,128]#[512,64,64]
    ζ = [80, 40]
    Mw = 6.941
    F = 96485
    R = 8.314
    T = 298.15
    ρₛ = 5.34e5
    n = 1
    σₑ = 0.05
    Iₐ = 10

    ρ = ["train" "valid" "test"]

    for i ∈ 1:length(ρ)

        plan, _ = LHCoptim((ψ[i])÷2,4,1000)
        δ = plan[:,1]./maximum(plan[:,1]).*-1.0
        κ = plan[:,2]./maximum(plan[:,2]).*0.05 .+0.94
        ki₀ = plan[:,4]./maximum(plan[:,4]).*0.2 .+ 1.46


        plan, _ = LHCoptim(ψ[i]÷2,4,1000)
        δ = [δ; plan[:,1]./maximum(plan[:,1]).*1.0]
        κ = [κ; plan[:,2]./maximum(plan[:,2]).*0.05]
        ki₀ = [ki₀; plan[:,4]./maximum(plan[:,4]).*0.2 .+ 1.46]

        ν = @. (Mw*3600*Iₐ)/(ρₛ*F*n*ξ*1e-6*δ)

        

        #@infiltrate cond=true

        Out = Data(ρ[i], StepNum, ti, tf, ψ[i], δ, κ, ζ, ξ, ν, ki₀)
        Output = flatten!(Output,Out)
        κ_out = flatten!(κ_out,κ)
        δ_out = flatten!(δ_out,δ)
        ν_out = flatten!(ν_out,ν)
        ki₀_out = flatten!(ki₀_out,ki₀)
    end
    
    return Output, δ_out, κ_out, ν_out, ki₀_out
end

Output, δ, κ, ν, ki₀ = LatinHyper()