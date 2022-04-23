using BattPhase, Plots, BenchmarkTools, SparseArrays, ProfileView, Infiltrator, DifferentialEquations, MKL

    ## Discretisation Parameters
    Nx = NN = MM = Ny = Int64(160)
    δ = 0.1
    tf = 2.
    ν = ki₀ = 1.

    N = Int64(NN+2)
    h = 1/NN
    M = Int64(MM+2)
    Ntot = Int64(N*M)


    Y₀ = Array{Float64}(undef,N,M)
    eq = Array{Float64}(undef,N,M)
    F₀ = Array{Float64}(undef,N,M)
    Φ₀ = Vector{Float64}(undef,Ntot) .= 0
    #Φ₀ = sparsevec(Φ₀)
    ff = Vector{Float64}(undef,Ntot) .= 0
    #j₀ = Array{Float64}(undef,Ntot,Ntot) #Modify to be Sparse
    j₀ = spzeros(Ntot,Ntot)
    ymid = Array{Float64}(undef,N,M)

    #if δ<0 then read("Y₀data.m"):end: ??
    #if δ<0 then read("tdata.m"):end: ??

    #plotly()
    #surface(Y₀,camera = (0, 90))


    #Calls
    #kk = @benchmark Eqs11!($Y₀,$Φ₀,$δ,$ki₀,$ff,$N,$M)
    #ttq = @benchmark Jac!($Y₀,$Φ₀,$δ,$ki₀,$j₀,$N,$M)
    y₀!(NN,MM,Y₀)
    Eqs11!(Y₀,Φ₀,δ,ki₀,ff,N,M,h)
    Jac!(Y₀,Φ₀,δ,ki₀,j₀,N,M,h)

    dᵦ = j₀\-ff
    Φ₀ += dᵦ


    vv₀ = max(abs(Φ₀[Ntot-2*N+1]),abs(Φ₀[Ntot-2*N+N÷2]),abs(Φ₀[Ntot-2*N+N÷2+1]),abs(Φ₀[Ntot-N]))

    #@infiltrate cond=true

    tt = 0.
    dt = min(h/vv₀/ν/ki₀,tf-tt)
    Nₜ = round(Int64,tf/dt)+50
    V = Vector{Float64}(undef,Nₜ)
    TT = Vector{Float64}(undef,Nₜ) .= 0
    Φₐ = Vector{Float64}(undef,Nₜ)
    V[1] = (Φ₀[Int(Ntot-2*N+N/2)]/2+Φ₀[Int(Ntot-2*N+N/2+1)]/2)+h/2*δ
    #TT[1] = 0.


    Φₐ[1] = Φᵥ₊(Y₀,N,M,Ny)
    ymid = copy(Y₀)
    Ydata = Array{Float64}(undef,Nₜ+1,N,M)
    YStore(Y₀,Ydata,N,M,1)

    #qqw = @benchmark Ef!($Y₀,$ymid,$Φ₀,$F₀,$dt,$ν,$ki₀,$N,$M)
    #rrt = @benchmark Upwind!($Y₀,$Φ₀,$F₀,$dt,$N,$M,$δ)
    #fgv = @benchmark $j₀\-$ff
    ii = 1

    @inline function solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀::Float64,ymid,ii,j₀,Ntot,tt::Float64,tf::Float64,Φₐ,TT,V,Nₜ,ff,dᵦ,ν::Float64,vv₀::Float64,h::Float64)
        while tt < tf

            #Ef!
            Upwind!(Y₀,Φ₀,F₀,dt,N,M,δ,h)
            Ef!(Y₀,ymid,Φ₀,F₀,dt,ν,ki₀,N,M)
            Eqs11!(ymid,Φ₀,δ,ki₀,ff,N,M,h)
            Jac!(ymid,Φ₀,δ,ki₀,j₀,N,M,h)
            dᵦ = j₀\-ff
            Φ₊!(Ntot,Φ₀,dᵦ)

            #Ef2!
            Upwind!(ymid,Φ₀,F₀,dt,N,M,δ,h)
            Ef2!(Y₀,ymid,Φ₀,F₀,dt,ν,ki₀,N,M)
            Eqs11!(ymid,Φ₀,δ,ki₀,ff,N,M,h)
            Jac!(ymid,Φ₀,δ,ki₀,j₀,N,M,h)
            dᵦ = j₀\-ff
            Φ₊!(Ntot,Φ₀,dᵦ)

            #Ef3!
            Upwind!(ymid,Φ₀,F₀,dt,N,M,δ,h)
            Ef3!(Y₀,ymid,Φ₀,F₀,dt,ν,ki₀,N,M)
            Eqs11!(Y₀,Φ₀,δ,ki₀,ff,N,M,h)
            Jac!(Y₀,Φ₀,δ,ki₀,j₀,N,M,h)
            dᵦ = j₀\-ff
            Φ₊!(Ntot,Φ₀,dᵦ)


            @show ii += 1
            V[ii] = (Φ₀[Ntot-2*N+N÷2]/2+Φ₀[Ntot-2*N+N÷2+1]/2)+h/2*δ
            TT[ii] = TT[ii-1]+dt
            tt += dt
            YStore(Y₀,Ydata,N,M,ii+1)
            Φₐ[ii] = Φᵥ₊(Y₀,N,M,Ny)
            vv₀ = max(abs(Φ₀[Ntot-2*N+1]),abs(Φ₀[Ntot-2*N+N÷2]),abs(Φ₀[Ntot-2*N+N÷2+1]),abs(Φ₀[Ntot-N]))
            dt = min(h/vv₀/ν/ki₀,tf-tt)
            #YdatStore(Y₀,Ydat,N,M,i+1);
            #Φₐ[i]:=Φᵥ₊(Y₀,N,M,Ny);

        end
        return Ydata, V, TT
    end

    #@time solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀,ymid,ii,j₀,Ntot,tt,tf,Φₐ,TT,V,Nₜ)

    #return kk, ttq, qqw, rrt, fgv


@time solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀,ymid,ii,j₀,Ntot,tt,tf,Φₐ,TT,V,Nₜ,ff,dᵦ,ν,vv₀,h)
#kk,ttq, qqw, rrt, fgv = Circle()
#@profview solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀,ymid,ii,j₀,Ntot,tt,tf,Φₐ,TT,V,Nₜ,ff,dᵦ,ν,vv₀,h)
#@code_warntype solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀,ymid,ii,j₀,Ntot,tt,tf,Φₐ,TT,V,Nₜ,ff,dᵦ,ν,vv₀,h)
#Ydata, V, TT = solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀,ymid,ii,j₀,Ntot,tt,tf,Φₐ,TT,V,Nₜ,ff,dᵦ,ν,vv₀,h)