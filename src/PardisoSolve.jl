using Pardiso
ps = MKLPardisoSolver()

# Pardiso SSP-RK3 Numerical Loop
@inline function pdrk3solve(Y,Φ,F,dt,N,M,δ,ki,ymid,j,Ntot,tt,tf,TT,V,V1,V2,ff,dᵦ,ν,vv,h,Φₐ,Ydata)
    i = 1
   
    while tt < tf

        #Ef!
        Upwind!(Y,F,N,M,δ,h)
        Ef!(Y,ymid,Φ,F,dt,ν,ki,N,M)
        Eqs11!(ymid,Φ,δ,ki,ff,N,M,h)
        Jac!(ymid,ki,j,N,M,h)
        solve!(ps,dᵦ,j,-ff)
        Φ₊!(Ntot,Φ,dᵦ)

        #Ef2!
        Upwind!(ymid,F,N,M,δ,h)
        Ef2!(Y,ymid,Φ,F,dt,ν,ki,N,M)
        Eqs11!(ymid,Φ,δ,ki,ff,N,M,h)
        Jac!(ymid,ki,j,N,M,h)
        solve!(ps,dᵦ,j,-ff)
        Φ₊!(Ntot,Φ,dᵦ)

        #Ef3!
        Upwind!(ymid,F,N,M,δ,h)
        Ef3!(Y,ymid,Φ,F,dt,ν,ki,N,M)
        Eqs11!(Y,Φ,δ,ki,ff,N,M,h)
        Jac!(Y,ki,j,N,M,h)
        solve!(ps,dᵦ,j,-ff)
        Φ₊!(Ntot,Φ,dᵦ)

        i += 1
        V2[i] = (Φ[Ntot-2*N+N÷2]/2+Φ[Ntot-2*N+N÷2+1]/2)+h/2*δ
        V1[i] = Φ[2*N]
        V[i] = Φ[Ntot-N]+0.5*h*δ
        TT[i] = TT[i-1]+dt
        tt += dt
        YStore(Y,Ydata,N,M,i) 
        Φₐ[i] = Φ̄₊(Y,N,M)
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
    Φₐ = Φₐ[1:i]

    return Ydata, V2, V1, V2, TT, Φₐ
end


# Pardiso - SSP-RK3a Numerical Loop
@inline function pdrk3asolve(Y,Φ,Φₜ,Φₘ,F,dt,N,M,δ,ki,ymid,j,Ntot,tt,tf,TT,V,V1,V2,ff,dᵦ,ν,vv,h,Φₐ,Ydata)
    i = 1
    

    while tt < tf
        
        #Ef!
        PhiSwitch!(Ntot,Φ,Φₜ)
        Upwind!(Y,F,N,M,δ,h)
        Ef!(Y,ymid,Φ,F,dt,ν,ki,N,M)
        Eqs11!(ymid,Φ,δ,ki,ff,N,M,h)
        Jac!(ymid,ki,j,N,M,h)
        solve!(ps,dᵦ,j,-ff)
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
        Φₐ[i] = Φ̄₊(Y,N,M)
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
    Φₐ = Φₐ[1:i]

    return Ydata, V, V1, V2, TT, Φₐ
end