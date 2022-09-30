using BattPhase, LinearAlgebra, BenchmarkTools, SparseArrays, Plots, Infiltrator#, ProfileView, Infiltrator
  
function Seed1D(NN,MM,κ,δ,ν,ki₀,tt,tf,Δₜ,dims)

    ## Discretisation Parameters
    N = Int(NN+2)
    h = 1/NN
    M = Int(MM+2)
    Ntot = N*M

    # Pre-Allocations
    Y₀ = Array{Float64}(undef,N,M) .=0
    F₀ = Array{Float64}(undef,N,M) .=0
    Φ₀ = Vector{Float64}(undef,Ntot) .= 0
    ff₀ = Vector{Float64}(undef,Ntot) .= 0
    j₀ = spzeros(Ntot,Ntot)
    ymid = Array{Float64}(undef,N,M) .=0
    Φₜ₀ = Vector{Float64}(undef,Ntot) .= 0
    Φₘ₀ = Vector{Float64}(undef,Ntot) .= 0

    # Initialisations
    y₀4!(NN,MM,Y₀,κ)
    Eqs11!(Y₀,Φ₀,δ,ki₀,ff₀,N,M,h)
    Jac!(Y₀,ki₀,j₀,N,M,h)
    dᵦ₀ = j₀\-ff₀
    Φ₀ += dᵦ₀


    vv₀ = max(abs(Φ₀[Ntot-2*N+1]),abs(Φ₀[Ntot-2*N+N÷2]),abs(Φ₀[Ntot-2*N+N÷2+1]),abs(Φ₀[Ntot-N]))
    dt₀ = min(h/vv₀/ν/ki₀,tf-tt)
    Nₜ = ceil(Int64,(tf/Δₜ)+2)
    V2 = Vector{Float64}(undef,Nₜ) .= 0
    V1 = Vector{Float64}(undef,Nₜ) .= 0
    V = Vector{Float64}(undef,Nₜ) .= 0
    TT = Vector{Float64}(undef,Nₜ) .= 0
    Φₐ₀ = Vector{Float64}(undef,Nₜ) .= 0
    V2[1] = (Φ₀[Ntot-2*N+N÷2]+Φ₀[Ntot-2*N+N÷2+1])/2+h/2*δ
    V1[1] = Φ₀[2*N]
    V[1] = Φ₀[Ntot-N]+0.5*h*δ

    Φₐ₀[1] = Φ̄₊(Y₀,N,M)
    ymid₀ = copy(Y₀)
    Ydata = Array{Float64}(undef,Nₜ+1,N,M) .= 0
    YStore(Y₀,Ydata,N,M,1)

    ## Single Run ##
    # @time Ydata_rk3, V_rk3, V1_rk3, V2_rk3, TT_rk3, Φₐ_rk3 = rk3solve(Y₀,Φ₀,F₀,dt₀,N,M,δ,ki₀,ymid₀,j₀,Ntot,tt,tf,TT,V,V1,V2,ff₀,dᵦ₀,ν,vv₀,h,Φₐ₀,Ydata)
    Ydata_rk3a, V_rk3a, V1_rk3a, V2_rk3a, TT_rk3a, Φₐ_rk3a = rk3asolve(Y₀,Φ₀,Φₜ₀,Φₘ₀,F₀,dt₀,N,M,δ,ki₀,ymid₀,j₀,Ntot,tt,tf,TT,V,V1,V2,ff₀,dᵦ₀,ν,vv₀,h,Φₐ₀,Ydata,Δₜ)  

    #@infiltrate cond=true

    #2D to 1D
    if dims == 2
        Y2D = Array{Float64}(undef,size(Ydata_rk3a,1)-3,NN*MM) .= 0
        for k ∈ 1:size(Ydata_rk3a,1)-3
            Ψ = Array{Float64}(undef,0)
            for j ∈ 2:size(Ydata_rk3a,3)-1
                if isodd(j) == true
                    Ψ = [Ψ; reverse(Ydata_rk3a[k,j,2:end-1])]
                else
                    Ψ = [Ψ; Ydata_rk3a[k,j,2:end-1]]
                end
            end
        Y2D[k,:] = Ψ
        end
        return Y2D, κ
    elseif dims == 1
        Y2D = Array{Float64}(undef,size(Ydata_rk3a,1)-3,NN,MM) .= 0
        for k ∈ 1:size(Ydata_rk3a,1)-3
            for j ∈ 2:size(Ydata_rk3a,3)-1
                Y2D[k,:,j-1] = Ydata_rk3a[k,2:end-1,j]
            end
        end
        return Y2D, κ
    end
                

end