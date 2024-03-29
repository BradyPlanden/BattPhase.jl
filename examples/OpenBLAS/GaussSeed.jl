using BattPhase, LinearAlgebra, BenchmarkTools, SparseArrays, Plots
   
    ## Discretisation Parameters
    Nx = NN = MM = Ny = Int(320)
    δ = -0.1
    tt = 0.
    tf = 1.

    Mw = 6.941
    F = 96485
    R = 8.314
    T = 298.15
    L = 40e-6
    ρ = 5.34e5
    n = 1
    Iₐ = 10
    i₀ = 50
    κₑ = 0.05

    ki₀ = (i₀*L*F)/(κₑ*R*T)
    ν = (3600*Mw*Iₐ)/(ρ*F*L*n*δ)
    
    
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
    y₀4!(NN,MM,Y₀)
    Eqs11!(Y₀,Φ₀,δ,ki₀,ff₀,N,M,h)
    Jac!(Y₀,ki₀,j₀,N,M,h)
    dᵦ₀ = j₀\-ff₀
    Φ₀ += dᵦ₀


    vv₀ = max(abs(Φ₀[Ntot-2*N+1]),abs(Φ₀[Ntot-2*N+N÷2]),abs(Φ₀[Ntot-2*N+N÷2+1]),abs(Φ₀[Ntot-N]))
    dt₀ = min(h/vv₀/abs(ν)/ki₀,tf-tt)
    Nₜ = ceil(Int64,tf/dt₀+1)
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


    ## Benchmark Rk3 ##
    # t =  @benchmarkable rk3solve(a,b,c,d,e,f,g,hh,ii,j,k,l,m,o,p,pa,pb,r,s,t,u,v,w,x)
     # setup = begin;
     #      a=copy($Y₀);
     #      b=copy($Φ₀);
     #      c=copy($F₀);
     #      d=copy($dt₀);
     #      e=copy($N);
     #      f=copy($M);
     #      g=copy($δ);
     #      hh=copy($ki₀);
     #      ii=copy($ymid₀);
     #      j=copy($j₀);
     #      k=copy($Ntot);
     #      l=copy($tt);
     #      m=copy($tf);
     #      o=copy($TT);
     #      p=copy($V);
     #      pa=copy($V1);
     #      pb=copy($V2);
     #      r=copy($ff₀);
     #      s=copy($dᵦ₀);
     #      t=copy($ν);
     #      u=copy($vv₀);
     #      v=copy($h);
     #      w=copy($Φₐ₀)
     #      x=copy($Ydata)
     # end
    # run(t, evals=1, seconds=500.0, samples = 7)

    ## Benchmark Rk3a ##
    #t =  @benchmarkable rk3asolve(a,b,ba,bb,c,d,e,f,g,hh,ii,j,k,l,m,o,p,pa,pb,r,s,t,u,v,w,x) 
     # setup = begin;
     #      a=copy($Y₀);
     #      b=copy($Φ₀);
     #      ba=copy($Φₜ₀);
     #      bb=copy($Φₘ₀);
     #      c=copy($F₀);
     #      d=copy($dt₀);
     #      e=copy($N);
     #      f=copy($M);
     #      g=copy($δ);
     #      hh=copy($ki₀);
     #      ii=copy($ymid₀);
     #      j=copy($j₀);
     #      k=copy($Ntot);
     #      l=copy($tt);
     #      m=copy($tf);
     #      o=copy($TT);
     #      p=copy($V);
     #      pa=copy($V1);
     #      pb=copy($V2);
     #      r=copy($ff₀);
     #      s=copy($dᵦ₀);
     #      t=copy($ν);
     #      u=copy($vv₀);
     #      v=copy($h);
     #      w=copy($Φₐ₀) 
     #      x=copy($Ydata)
     # end
     # run(t, evals=1, seconds=200.0, samples = 7)

    ## Single Run ##
     #@time Ydata_rk3, V_rk3, V1_rk3, V2_rk3, TT_rk3, Φₐ_rk3 = rk3solve(Y₀,Φ₀,F₀,dt₀,N,M,δ,ki₀,ymid₀,j₀,Ntot,tt,tf,TT,V,V1,V2,ff₀,dᵦ₀,ν,vv₀,h,Φₐ₀,Ydata)
     @time Ydata_rk3a, V_rk3a, V1_rk3a, V2_rk3a, TT_rk3a, Φₐ_rk3a = rk3asolve(Y₀,Φ₀,Φₜ₀,Φₘ₀,F₀,dt₀,N,M,δ,ki₀,ymid₀,j₀,Ntot,tt,tf,TT,V,V1,V2,ff₀,dᵦ₀,ν,vv₀,h,Φₐ₀,Ydata) 


    ## Plotting ##
    TT_trunc = trunc.(TT_rk3a, digits=1)
    kk = Vector{String}(undef,length(V2_rk3a))
    for i ∈ 1:length(V2_rk3a)
        if TT[i]<=tf/2
            kk[i] = "Plating"
        else
            kk[i] = "Stripping"
        end
    end
    Lₛ = 40
    x = y = Lₛ/NN:Lₛ/NN:Lₛ

    gr()
    anim = @animate for i ∈ 1:length(V2_rk3a)
        heatmap(x,y,Ydata_rk3a[i,2:end-1,2:end-1]', annotations = (38, 38, Plots.text(kk[i], :center)), box=:on, c = :davos,bottom_margin=5Plots.mm, left_margin = 7.5Plots.mm, right_margin = 0Plots.mm, top_margin = 5Plots.mm, ylabel = "Position (μm)", xlabel = "Position (μm)",title="Lithium Metal Anode Evolution\nfor Gaussian Seed at $(TT_trunc[i]) Hr", xtickfontsize=16,ytickfontsize=16,xguidefontsize=16,yguidefontsize=16, size=(1280,720))#GnBu_3
         annotate!(20, 400, Plots.text("$(TT_trunc[i]) Hr", :center))
    end
    gif(anim, "anim_fps15.gif", fps = 15) 

    plotly()
    l = @layout [a b c{0.33w}]
    p1 = plot(heatmap(x,y,Ydata_rk3a[1,2:end-1,2:end-1]',clims=(0, 1), color = :davos), colorbar=false)
    p2 = plot(heatmap(x,y,Ydata_rk3a[50,2:end-1,2:end-1]',clims=(0, 1), color = :davos), colorbar=false)
    p3 = plot(heatmap(x,y,Ydata_rk3a[100,2:end-1,2:end-1]',clims=(0, 1), color = :davos))
    plot(p1,p2,p3, color = :davos, layout=l, bottom_margin=5Plots.mm, left_margin = 7.5Plots.mm, right_margin = 0Plots.mm, top_margin = 5Plots.mm, ylabel = "Position (μm)", xlabel = "Position (μm)", xtickfontsize=16, ytickfontsize=16, xguidefontsize=16, yguidefontsize=16,size=(1920,1080))