using BattPhase, LinearAlgebra, Plots, BenchmarkTools, SparseArrays#, ProfileView, Infiltrator, DifferentialEquations, MKL

    ## Discretisation Parameters
    Nx = NN = MM = Ny = Int(160)
    δ = 0.1
    tt = 0.
    tf = 2.
    ν = ki₀ = 1.
    
    N = Int(NN+2)
    h = 1/NN
    M = Int64(MM+2)
    Ntot = N*M

    # Pre-Allocations
    Y₀ = Array{Float64}(undef,N,M)
    F₀ = Array{Float64}(undef,N,M)
    Φ₀ = Vector{Float64}(undef,Ntot) .= 0
    ff₀ = Vector{Float64}(undef,Ntot) .= 0
    j₀ = spzeros(Ntot,Ntot)
    ymid = Array{Float64}(undef,N,M)
    Φₜ₀ = Vector{Float64}(undef,Ntot) .= 0
    Φₘ₀ = Vector{Float64}(undef,Ntot) .= 0

    y₀!(NN,MM,Y₀)
    Eqs11!(Y₀,Φ₀,δ,ki₀,ff₀,N,M,h)
    Jac!(Y₀,Φ₀,δ,ki₀,j₀,N,M,h)

    dᵦ₀ = j₀\-ff₀
    Φ₀ += dᵦ₀

    vv₀ = max(abs(Φ₀[Ntot-2*N+1]),abs(Φ₀[Ntot-2*N+N÷2]),abs(Φ₀[Ntot-2*N+N÷2+1]),abs(Φ₀[Ntot-N]))
    dt₀ = min(h/vv₀/ν/ki₀,tf-tt)
    Nₜ = round(Int64,tf/dt₀)
    V2 = Vector{Float64}(undef,Nₜ) .= 0
    V1 = Vector{Float64}(undef,Nₜ) .= 0
    V = Vector{Float64}(undef,Nₜ) .= 0
    TT = Vector{Float64}(undef,Nₜ) .= 0
    Φₐ₀ = Vector{Float64}(undef,Nₜ) .= 0
    V[1] = (Φ₀[Ntot-2*N+N÷2]/2+Φ₀[Ntot-2*N+N÷2+1]/2)+h/2*δ

    Φₐ₀[1] = Φ̄₊(Y₀,N,M,Ny)
    ymid₀ = copy(Y₀)
    Ydata = Array{Float64}(undef,Nₜ+1,N,M)
    YStore(Y₀,Ydata,N,M,1)


   
    # SSP-RK3 Numerical Loop
    @inline function rk3solve!(Y,Φ,F,dt,N,M,δ,ki,ymid,j,Ntot,tt,tf,TT,V,ff,dᵦ,ν,vv,h)
        i = 1
       
        while tt < tf

            #Ef!
            Upwind!(Y,Φ,F,dt,N,M,δ,h)
            Ef!(Y,ymid,Φ,F,dt,ν,ki,N,M)
            Eqs11!(ymid,Φ,δ,ki,ff,N,M,h)
            Jac!(ymid,Φ,δ,ki,j,N,M,h)
            dᵦ = j\-ff
            Φ₊!(Ntot,Φ,dᵦ)

            #Ef2!
            Upwind!(ymid,Φ,F,dt,N,M,δ,h)
            Ef2!(Y,ymid,Φ,F,dt,ν,ki,N,M)
            Eqs11!(ymid,Φ,δ,ki,ff,N,M,h)
            Jac!(ymid,Φ,δ,ki,j,N,M,h)
            dᵦ = j\-ff
            Φ₊!(Ntot,Φ,dᵦ)

            #Ef3!
            Upwind!(ymid,Φ,F,dt,N,M,δ,h)
            Ef3!(Y,ymid,Φ,F,dt,ν,ki,N,M)
            Eqs11!(Y,Φ,δ,ki,ff,N,M,h)
            Jac!(Y,Φ,δ,ki,j,N,M,h)
            dᵦ = j\-ff
            Φ₊!(Ntot,Φ,dᵦ)

            i += 1
            V[i] = (Φ[Ntot-2*N+N÷2]/2+Φ[Ntot-2*N+N÷2+1]/2)+h/2*δ
            TT[i] = TT[i-1]+dt
            tt += dt
            #YStore(Y,Ydata,N,M,i+1) 
            #Φₐ[i] = Φ̄₊(Y,N,M,Ny)
            vv = max(abs(Φ[Ntot-2*N+1]),abs(Φ[Ntot-2*N+N÷2]),abs(Φ[Ntot-2*N+N÷2+1]),abs(Φ[Ntot-N]))
            dt = min(h/vv/ν/ki,tf-tt)

        end
        return Ydata, V, TT
    end


    # SSP-RK3a Numerical Loop
    @inline function rk3asolve!(Y,Φ,Φₜ,Φₘ,F,dt,N,M,δ,ki,ymid,j,Ntot,tt,tf,TT,V,V1,V2,ff,dᵦ,ν,vv,h)
        i = 1

        while tt < tf

            #Ef!
            PhiSwitch!(N,Φ,Φₜ)
            Upwind!(Y,Φ,F,dt,N,M,δ,h)
            Ef!(Y,ymid,Φ,F,dt,ν,ki,N,M)
            Eqs11!(ymid,Φ,δ,ki,ff,N,M,h)
            Jac!(ymid,Φ,δ,ki,j,N,M,h)
            dᵦ = j\-ff
            Φ₊!(Ntot,Φ,dᵦ)

            #Ef2!
            MidPred!(N,Φ,Φₜ,Φₘ)
            Upwind!(ymid,Φ,F,dt,N,M,δ,h)
            Ef2!(Y,ymid,Φ,F,dt,ν,ki,N,M)

            #Ef3!
            Upwind!(ymid,Φₘ,F,dt,N,M,δ,h)
            Ef3!(Y,ymid,Φₘ,F,dt,ν,ki,N,M)

            i += 1
            V2[i] = (Φ[Ntot-2*N+N÷2]/2+Φ[Ntot-2*N+N÷2+1]/2)+h/2*δ
            V1[i] = Φ[2*N]
            V[i] = Φ[Ntot-N]+0.5*h*0.1
            TT[i] = TT[i-1]+dt
            tt += dt
            #YStore(Y,Ydata,N,M,i+1) 
            #Φₐ[i] = Φ̄₊(Y,N,M,Ny)
            vv = max(abs(Φ[Ntot-2*N+1]),abs(Φ[Ntot-2*N+N÷2]),abs(Φ[Ntot-2*N+N÷2+1]),abs(Φ[Ntot-N]))
            dt = min(h/vv/ν/ki,tf-tt)

        end
        return Ydata, V, V1, V2, TT
    end


    # Benchmark Rk3
    t =  @benchmarkable rk3solve!(a,b,c,d,e,f,g,hh,ii,j,k,l,m,o,p,r,s,t,u,v) setup = begin; a=copy($Y₀);b=copy($Φ₀);c=copy($F₀);d=copy($dt₀);e=copy($N);f=copy($M);g=copy($δ);hh=copy($ki₀);ii=copy($ymid₀);j=copy($j₀);k=copy($Ntot);l=copy($tt);m=copy($tf);o=copy($TT);p=copy($V);r=copy($ff₀);s=copy($dᵦ₀);t=copy($ν);u=copy($vv₀);v=copy($h) end
    run(t, evals=1, seconds=500.0, samples = 7)

    # Benchmark Rk3a
    # t =  @benchmarkable rk3asolve!(a,b,ba,bb,c,d,e,f,g,hh,ii,j,k,l,m,o,p,pa,pb,r,s,t,u,v) setup = begin; a=copy($Y₀);b=copy($Φ₀);ba=copy($Φₜ);bb=copy($Φₘ);c=copy($F₀);d=copy($dt₀);e=copy($N);f=copy($M);g=copy($δ);hh=copy($ki₀);ii=copy($ymid₀);j=copy($j₀);k=copy($Ntot);l=copy($tt);m=copy($tf);o=copy($TT);p=copy($V);pa=copy($V1);pb=copy($V2);r=copy($ff₀);s=copy($dᵦ₀);t=copy($ν);u=copy($vv₀);v=copy($h) end
    # run(t, evals=1, seconds=500.0, samples = 7)

    # Run
    #Ydata_rk3, V_rk3, TT_rk3 = rk3solve!(Y₀,Φ₀,F₀,dt₀,N,M,δ,ki₀,ymid₀,j₀,Ntot,tt,tf,TT,V,ff₀,dᵦ₀,ν,vv₀,h)
    #Ydata_rk3a, V_rk3a, V1_rk3a, V2_rk3a, TT_rk3a = rk3asolve!(Y₀,Φ₀,Φₜ₀,Φₘ₀,F₀,dt₀,N,M,δ,ki₀,ymid₀,j₀,Ntot,tt,tf,TT,V,V1,V2,ff₀,dᵦ₀,ν,vv₀,h)
    
    ## Performance Debug Tools ##
    #@profview solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀,ymid,ii,j₀,Ntot,tt,tf,Φₐ,TT,V,Nₜ,ff,dᵦ,ν,vv₀,h)
    #@code_warntype solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀,ymid,ii,j₀,Ntot,tt,tf,Φₐ,TT,V,Nₜ,ff,dᵦ,ν,vv₀,h)


    # Plotting
    # plotly()
    # surface(Y₀,camera = (0, 90))