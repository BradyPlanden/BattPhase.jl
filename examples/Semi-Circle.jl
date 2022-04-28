using BattPhase, LinearAlgebra, BenchmarkTools, SparseArrays, Plots#, ProfileView, Infiltrator

    ## Discretisation Parameters
    Nx = NN = MM = Ny = Int(320)
    δ = 0.1
    tt = 0.
    tf = 2.
    ν = ki₀ = 1.
    
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
    y₀!(NN,MM,Y₀)
    Eqs11!(Y₀,Φ₀,δ,ki₀,ff₀,N,M,h)
    Jac!(Y₀,ki₀,j₀,N,M,h)

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
    V2[1] = (Φ₀[Ntot-2*N+N÷2]+Φ₀[Ntot-2*N+N÷2+1])/2+h/2*δ
    V1[1] = Φ₀[2*N]
    V[1] = Φ₀[Ntot-N]+0.5*h*δ

    Φₐ₀[1] = Φ̄₊(Y₀,N,M,Ny)
    ymid₀ = copy(Y₀)
    Ydata = Array{Float64}(undef,Nₜ+1,N,M) .=0
    YStore(Y₀,Ydata,N,M,1)

   
    # SSP-RK3 Numerical Loop
    @inline function rk3solve(Y,Φ,F,dt,N,M,δ,ki,ymid,j,Ntot,tt,tf,TT,V,V1,V2,ff,dᵦ,ν,vv,h,Φₐ)
        i = 1
       
        while tt < tf

            #Ef!
            Upwind!(Y,F,N,M,δ,h)
            Ef!(Y,ymid,Φ,F,dt,ν,ki,N,M)
            Eqs11!(ymid,Φ,δ,ki,ff,N,M,h)
            Jac!(ymid,ki,j,N,M,h)
            dᵦ = j\-ff
            Φ₊!(Ntot,Φ,dᵦ)

            #Ef2!
            Upwind!(ymid,F,N,M,δ,h)
            Ef2!(Y,ymid,Φ,F,dt,ν,ki,N,M)
            Eqs11!(ymid,Φ,δ,ki,ff,N,M,h)
            Jac!(ymid,ki,j,N,M,h)
            dᵦ = j\-ff
            Φ₊!(Ntot,Φ,dᵦ)

            #Ef3!
            Upwind!(ymid,F,N,M,δ,h)
            Ef3!(Y,ymid,Φ,F,dt,ν,ki,N,M)
            Eqs11!(Y,Φ,δ,ki,ff,N,M,h)
            Jac!(Y,ki,j,N,M,h)
            dᵦ = j\-ff
            Φ₊!(Ntot,Φ,dᵦ)

            i += 1
            V[i] = (Φ[Ntot-2*N+N÷2]/2+Φ[Ntot-2*N+N÷2+1]/2)+h/2*δ
            V1[i] = Φ[2*N]
            V[i] = Φ[Ntot-N]+0.5*h*δ
            TT[i] = TT[i-1]+dt
            tt += dt
            YStore(Y,Ydata,N,M,i) 
            Φₐ[i] = Φ̄₊(Y,N,M,Ny)
            vv = max(abs(Φ[Ntot-2*N+1]),abs(Φ[Ntot-2*N+N÷2]),abs(Φ[Ntot-2*N+N÷2+1]),abs(Φ[Ntot-N]))
            dt = min(h/vv/ν/ki,tf-tt)

            if TT[i] >= tf/2
                δ = -0.1
            end

        end

        ## Shrink Vectors ##
        V2 = V2[1:i]
        V1 = V1[1:i]
        V = V[1:i]
        TT = TT[1:i]

        return Ydata, V2, V1, V2, TT, Φₐ
    end


    # SSP-RK3a Numerical Loop
    @inline function rk3asolve(Y,Φ,Φₜ,Φₘ,F,dt,N,M,δ,ki,ymid,j,Ntot,tt,tf,TT,V,V1,V2,ff,dᵦ,ν,vv,h,Φₐ)
        i = 1
    
        while tt < tf
            
            #Ef!
            PhiSwitch!(Ntot,Φ,Φₜ)
            Upwind!(Y,F,N,M,δ,h)
            Ef!(Y,ymid,Φ,F,dt,ν,ki,N,M)
            Eqs11!(ymid,Φ,δ,ki,ff,N,M,h)
            Jac!(ymid,ki,j,N,M,h)
            dᵦ = j\-ff
            Φ₊!(Ntot,Φ,dᵦ)

            #Ef2!
            MidPred!(Ntot,Φ,Φₜ,Φₘ)
            Upwind!(ymid,F,N,M,δ,h)
            Ef2!(Y,ymid,Φ,F,dt,ν,ki,N,M)

            #Ef3!
            Upwind!(ymid,F,N,M,δ,h)
            Ef3!(Y,ymid,Φₘ,F,dt,ν,ki,N,M)

            i += 1
            V2[i] = (Φ[Ntot-2*N+N÷2]+Φ[Ntot-2*N+N÷2+1])/2+h/2*δ
            V1[i] = Φ[2*N]
            V[i] = Φ[Ntot-N]+0.5*h*δ
            TT[i] = TT[i-1]+dt
            tt += dt
            YStore(Y,Ydata,N,M,i) 
            Φₐ[i] = Φ̄₊(Y,N,M,Ny)
            vv = max(abs(Φ[Ntot-2*N+1]),abs(Φ[Ntot-2*N+N÷2]),abs(Φ[Ntot-2*N+N÷2+1]),abs(Φ[Ntot-N]))
            dt = min(h/vv/ν/ki,tf-tt)

            if TT[i] >= tf/2
                δ = -0.1
            end

        end

        ## Shrink Vectors ##
        V2 = V2[1:i]
        V1 = V1[1:i]
        V = V[1:i]
        TT = TT[1:i]

        return Ydata, V, V1, V2, TT, Φₐ
    end


    ## Benchmark Rk3 ##
    # t =  @benchmarkable rk3solve(a,b,c,d,e,f,g,hh,ii,j,k,l,m,o,p,pa,pb,r,s,t,u,v,w) setup = begin; a=copy($Y₀);b=copy($Φ₀);c=copy($F₀);d=copy($dt₀);e=copy($N);f=copy($M);g=copy($δ);hh=copy($ki₀);ii=copy($ymid₀);j=copy($j₀);k=copy($Ntot);l=copy($tt);m=copy($tf);o=copy($TT);p=copy($V);pa=copy($V1);pb=copy($V2);r=copy($ff₀);s=copy($dᵦ₀);t=copy($ν);u=copy($vv₀);v=copy($h);w=copy($Φₐ₀) end
    # run(t, evals=1, seconds=500.0, samples = 7)

    ## Benchmark Rk3a ##
    #t =  @benchmarkable rk3asolve(a,b,ba,bb,c,d,e,f,g,hh,ii,j,k,l,m,o,p,pa,pb,r,s,t,u,v,w) setup = begin; a=copy($Y₀);b=copy($Φ₀);ba=copy($Φₜ₀);bb=copy($Φₘ₀);c=copy($F₀);d=copy($dt₀);e=copy($N);f=copy($M);g=copy($δ);hh=copy($ki₀);ii=copy($ymid₀);j=copy($j₀);k=copy($Ntot);l=copy($tt);m=copy($tf);o=copy($TT);p=copy($V);pa=copy($V1);pb=copy($V2);r=copy($ff₀);s=copy($dᵦ₀);t=copy($ν);u=copy($vv₀);v=copy($h);w=copy($Φₐ₀) end
    #run(t, evals=1, seconds=200.0, samples = 7)

    ## Run ##
    # @time Ydata_rk3, V_rk3, V1_rk3, V2_rk3, TT_rk3, Φₐ_rk3 = rk3solve(Y₀,Φ₀,F₀,dt₀,N,M,δ,ki₀,ymid₀,j₀,Ntot,tt,tf,TT,V,V1,V2,ff₀,dᵦ₀,ν,vv₀,h,Φₐ₀)
     @time Ydata_rk3a, V_rk3a, V1_rk3a, V2_rk3a, TT_rk3a, Φₐ_rk3a = rk3asolve(Y₀,Φ₀,Φₜ₀,Φₘ₀,F₀,dt₀,N,M,δ,ki₀,ymid₀,j₀,Ntot,tt,tf,TT,V,V1,V2,ff₀,dᵦ₀,ν,vv₀,h,Φₐ₀)
    
    ## Performance Debug Tools ##
    # @profview solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀,ymid,ii,j₀,Ntot,tt,tf,Φₐ,TT,V,Nₜ,ff,dᵦ,ν,vv₀,h)
    # @code_warntype solve!(Y₀,Φ₀,F₀,dt,N,M,δ,ki₀,ymid,ii,j₀,Ntot,tt,tf,Φₐ,TT,V,Nₜ,ff,dᵦ,ν,vv₀,h)


    ## Plotting ##
    #  plotly()
    #  surface((Ydata[43,:,:]-Ydata[1,:,:]),camera = (0, 90))
    #  surface((Ydata[1,:,:]),camera = (0, 90))
    #  surface((Ydata[60,:,:]),camera = (0, 90))

    TT_trunc = trunc.(TT_rk3a, digits=1)
    kk = Vector{String}(undef,length(V2_rk3a))
    for i ∈ 1:length(V2_rk3a)
        if TT[i]<=tf/2
            kk[i] = "Plating"
        else
            kk[i] = "Stripping"
        end
    end

    anim = @animate for i ∈ 1:length(V2_rk3a)
        heatmap(Ydata_rk3a[i,2:end-1,2:end-1], annotations = (300, 300, Plots.text(kk[i], :center)), box=:on, c = :davos,bottom_margin=5Plots.mm, left_margin = 7.5Plots.mm, right_margin = 0Plots.mm, top_margin = 5Plots.mm, ylabel = "Position (μm)", xlabel = "Position (μm)",title="Lithium Metal Anode Evolution\nfor Semi-Circle Seed", size=(1280,720))#GnBu_3
        annotate!(300, 20, Plots.text("$(TT_trunc[i]) Hr", :center))
    end
    gif(anim, "anim_fps15.gif", fps = 15)
    #plot(heatmap(Ydata_rk3[10,2:end-1,2:end-1], annotations = (300, 300, Plots.text(kk[2], :left)), box=:on, c = :davos,bottom_margin=5Plots.mm, left_margin = 7.5Plots.mm, right_margin = 0Plots.mm, top_margin = 5Plots.mm, ylabel = "Position (μm)", xlabel = "Position (μm)",title="Lithium Metal Anode Evolution\nfor Semi-Circle Seed", size=(1920,1080)))#GnBu_3)
    