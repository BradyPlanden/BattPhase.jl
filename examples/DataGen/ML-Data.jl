using BattPhase, HDF5, OrderedCollections, FileIO, Infiltrator, LatinHypercubeSampling

function Data(ρ, StepNum, ti, tf, ψ, δ, κ, ζ, ξ, ν, ki₀, dims, γ)
    PhaseOut = κout = tuple()
    Δₜ = (tf-ti)/StepNum

    if dims == 2
        Z = ζ.^2
    else
        Z = ζ
    end

    for η ∈ ζ
        if dims == 2
            MM = η
            NN = η
            Ydata₁ = Array{Float64}(undef,ψ,StepNum,NN*MM)
            
        elseif dims == 1
            MM = η
            NN = η
            Ydata₁ = Array{Float64}(undef,ψ,StepNum,NN,MM)
        end
        
        κ₁ = Vector{Float64}(undef,ψ)

        for i ∈ 1:ψ
            # if iseven(ψ) == true 
            #     @show i
            # end
            if dims == 2
                Ydata₁[i,:,:], κ₁[i] = Seed1D(NN,MM,κ[i],δ[i],ν[i],ki₀[i],ti,tf,Δₜ,dims,γ[i]) #Seed1D(NN,MM,κ[i],δ[i],ν[i],ki₀[i],ti,tf,Δₜ)
            elseif dims == 1
                Ydata₁[i,:,:,:], κ₁[i] = Seed1D(NN,MM,κ[i],δ[i],ν[i],ki₀[i],ti,tf,Δₜ,dims,γ[i])
            end

        end

        @show size(Ydata₁)
        if dims == 2 
            PhaseOut = flatten!(PhaseOut,permutedims(Ydata₁,[3,2,1])) # Align the output .h5 for python 
        else dims == 1
            PhaseOut = flatten!(PhaseOut,permutedims(Ydata₁,[3,4,2,1])) # Align the output .h5 for python 
        end
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
        if dims == 2
            g["pde_$(StepNum)-$(Z[i])"] = PhaseOut[i]
        elseif dims == 1
            g["pde_$(StepNum)-$(Z[i])"] = PhaseOut[i][size(PhaseOut[i],1),:,:,:]
        end
        attributes(g["pde_$(StepNum)-$(Z[i])"])["dt"] = Δₜ
        attributes(g["pde_$(StepNum)-$(Z[i])"])["dx"] = ξ/(Z[i]-1)
        attributes(g["pde_$(StepNum)-$(Z[i])"])["nt"] = StepNum
        attributes(g["pde_$(StepNum)-$(Z[i])"])["nx"] = ζ[i]
        attributes(g["pde_$(StepNum)-$(Z[i])"])["tmax"] = tf
        attributes(g["pde_$(StepNum)-$(Z[i])"])["tmin"] = ti
        attributes(g["pde_$(StepNum)-$(Z[i])"])["x"] = Vector(range(0,ξ,step=(ξ/(Z[i]-1)))) 
    end
    close(f)

return PhaseOut

end

function LatinHyper()
    Output = δ_out = κ_out = ν_out = ki₀_out = γ_out = tuple()
    dims = 2
    ξ = 40 # Physical Spacial Range
    ti = 0.
    tf = 1.
    StepNum = 50
    #ν = ki₀ = 1.
    ψ = [512,64,64]#[1024,128,128]#[512,64,64]
    ζ = [40]
    Mw = 6.941
    F = 96485
    R = 8.314
    T = 298.15
    ρₛ = 5.34e5
    n = 1
    σₑ = 0.05
    Iₐ = 10

    ρ = ["train" "valid" "test"]

    if dims == 2
    for i ∈ 1:length(ρ)

        plan, _ = LHCoptim((ψ[i])÷2,4,1000)
        δ = plan[:,1]./maximum(plan[:,1]).*3.0#-0.6
        κ = plan[:,2]./maximum(plan[:,2]).*0.75#.*0.3 .+ 0.65 
        γ = plan[:,3]./maximum(plan[:,3]).*0.75 .+ 0.125 
        #ki₀ = plan[:,4]./maximum(plan[:,4]).*0.2 .+ 1.455


        plan, _ = LHCoptim(ψ[i]÷2,4,1000)
        δ = [δ; plan[:,1]./maximum(plan[:,1]).*-3.0]#0.6]
        κ = [κ; plan[:,2]./maximum(plan[:,2]).*0.6 .+ 0.25] #For grid size = 40 (max(κ) = 13.5 for full plated grid), grid = 80, κ = 27
        γ = [γ; plan[:,3]./maximum(plan[:,3]).*0.75 .+ 0.125]
        #ki₀ = [ki₀; plan[:,4]./maximum(plan[:,4]).*0.2 .+ 1.455]

        ν = @. (Mw*3600*Iₐ)/(ρₛ*F*n*ξ*1e-6*abs(δ))
        ki₀ = ones(ψ[i]).*1.557
        

        #@infiltrate cond=true

        Out = Data(ρ[i], StepNum, ti, tf, ψ[i], δ, κ, ζ, ξ, ν, ki₀, dims,γ)
        Output = flatten!(Output,Out)
        κ_out = flatten!(κ_out,κ)
        δ_out = flatten!(δ_out,δ)
        ν_out = flatten!(ν_out,ν)
        ki₀_out = flatten!(ki₀_out,ki₀)
        γ_out = flatten!(γ_out,γ)
    end

elseif dims == 1
    for i ∈ 1:length(ρ)

        plan, _ = LHCoptim((ψ[i])÷2,4,1000)
        δ = plan[:,1]./maximum(plan[:,1]).*-0.3#-3.0
        κ = plan[:,2]./maximum(plan[:,2]).*0.6 .+ 0.35
        γ = plan[:,3]./maximum(plan[:,3]).*0.75 .+ 0.125 
        #ki₀ = plan[:,4]./maximum(plan[:,4]).*0.2 .+ 1.455


        plan, _ = LHCoptim(ψ[i]÷2,4,1000)
        δ = [δ; plan[:,1]./maximum(plan[:,1]).*0.3]
        κ = [κ; plan[:,2]./maximum(plan[:,2]).*0.75] #For grid size = 40 (max(κ) = 13.5 for full plated grid), grid = 80, κ = 27
        γ = [γ; plan[:,3]./maximum(plan[:,3]).*0.75 .+ 0.125]
        #ki₀ = [ki₀; plan[:,4]./maximum(plan[:,4]).*0.2 .+ 1.455]

        ν = ones(ψ[i]) #@. (Mw*3600*Iₐ)/(ρₛ*F*n*ξ*1e-6*δ)
        ki₀ = ones(ψ[i]) #ones(ψ[i]).*1.557
        

        #@infiltrate cond=true

        Out = Data(ρ[i], StepNum, ti, tf, ψ[i], δ, κ, ζ, ξ, ν, ki₀, dims,γ)
        Output = flatten!(Output,Out)
        κ_out = flatten!(κ_out,κ)
        δ_out = flatten!(δ_out,δ)
        ν_out = flatten!(ν_out,ν)
        ki₀_out = flatten!(ki₀_out,ki₀)
        γ_out = flatten!(γ_out,γ)
    end
    
end
    
    return Output, δ_out, κ_out, ν_out, ki₀_out, γ_out
end

Output, δ, κ, ν, ki₀, γ = LatinHyper()