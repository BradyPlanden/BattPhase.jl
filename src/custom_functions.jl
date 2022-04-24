## Model 1 Geometry
function y₀!(NN,MM,Y₀)
    N = NN+2
    h = 1.0/NN
    w = h/2
    M = Ny = MM+2
    for j ∈ 2:M-1 @simd for i ∈ 2:N-1
        rr = ((i-3/2)*h)^2+((j-3/2)*h)^2 #combining two variables
        Y₀[i,j] = max(eps(),0.5*tanh(sqrt(2)*(sqrt(rr)-0.3)/w)+0.5) 

        for i ∈ 1:N 
            Y₀[i,1] = Y₀[i,2]
            Y₀[i,M] = Y₀[i,M-1]
            for j ∈ 1:M 
                Y₀[1,j] = Y₀[2,j]
                Y₀[N,j] = Y₀[N-1,j]
            end
        end
    end end
end  

## System Eqs
function Eqs11!(Y₀,y,δ,ki₀,ff,N,M,h)

    @inbounds for j ∈ 2:M-1 @simd for i ∈ 2:N-1 
        i1 = i+(j-1)*N
        @views ff[i1] = (Y₀[i,j]+Y₀[i,j+1])*(y[i1+N]-y[i1]) -
        (Y₀[i,j]+Y₀[i,j-1])*(y[i1]-y[i1-N]) +
        (Y₀[i+1,j]+Y₀[i,j])*(y[i1+1]-y[i1]) -
        (Y₀[i,j]+Y₀[i-1,j])*(y[i1]-y[i1-1]) -
        ki₀*y[i1]*h*(1e-24+(Y₀[i+1,j]-Y₀[i-1,j])^2+(Y₀[i,j+1]-Y₀[i,j-1])^2)^(1/2)
    end end

    @inbounds for j ∈ 1:M @simd for i ∈ 1:1 
        i1 = i+(j-1)*N
        @views ff[i1] = -y[i1]+y[i1+1]
    end end

    @inbounds for j ∈ 1:M @simd for i ∈ N:N 
        i1 = i+(j-1)*N
        @views ff[i1] = -y[i1]+y[i1-1]
    end end

    @inbounds for j ∈ 1:1 @simd for i ∈ 1:N 
        i1 = i+(j-1)*N
        @views ff[i1] = -y[i1]+y[i1+N]
    end end

    @inbounds for j ∈ M:M @simd for i ∈ 1:N 
        i1 = i+(j-1)*N
        @views ff[i1] = -y[i1]+y[i1-N]+δ*h
    end end
end

## Jacobian
function Jac!(Y₀,y,δ,ki₀,j₀,N,M,h)

    @inbounds for j ∈ 2:M-1 @simd for i ∈ 2:N-1 
        i1 = i+(j-1)*N
        @views j₀[i1,i1] = -4*Y₀[i,j]-Y₀[i,j+1]-Y₀[i,j-1]-Y₀[i+1,j]-Y₀[i-1,j]-
        ki₀*h*(eps()+(Y₀[i+1,j]-Y₀[i-1,j])^2+(Y₀[i,j+1]-Y₀[i,j-1])^2)^(1/2)
        @views j₀[i1,i1+1] = Y₀[i+1,j]+Y₀[i,j]
        @views j₀[i1,i1-1] = Y₀[i-1,j]+Y₀[i,j]
        @views j₀[i1,i1+N] = Y₀[i,j+1]+Y₀[i,j]
        @views j₀[i1,i1-N] = Y₀[i,j-1]+Y₀[i,j]
    end end

    @inbounds for j ∈ 1:M @simd for i ∈ 1:1 
        i1 = i+(j-1)*N
        @views j₀[i1,i1] = -1.
        @views j₀[i1,i1+1] =1.
        #ff[i1] = -y[i1]+y[i1+1]
    end end

    @inbounds for j ∈ 1:M @simd for i ∈ N:N 
        i1 = i+(j-1)*N
        @views j₀[i1,i1] = -1.
        @views j₀[i1,i1-1] = 1.
        #ff[i1]:=-y[i1]+y[i1-1]:
    end end

    @inbounds for j ∈ 1:1 @simd for i ∈ 1:N 
        i1 = i+(j-1)*N
        #ff[i1] = -y[i1]+y[i1+N]
        @views j₀[i1,i1] = -1.
        @views j₀[i1,i1+N] =1.
    end end

    @inbounds for j ∈ M:M @simd for i ∈ 1:N 
        i1 = i+(j-1)*N
        @views j₀[i1,i1] = -1.
        @views j₀[i1,i1-N] = 1.
        #ff[i1] = -y[i1]+y[i1-N]+δ*h
    end end 
end


## Upwind
function Upwind!(Y₀,Φ₀,F₀,dt,N,M,v₀,h)

    vv₀ = v₀
    @views for i ∈ 1:N
        Y₀[i,1] = Y₀[i,2]
        Y₀[i,M] = Y₀[i,M-1]
    end

    @views for j ∈ 1:M
        Y₀[1,j] = Y₀[2,j]
        Y₀[N,j] = Y₀[N-1,j]
    end

    @inbounds for j ∈ 2:M-1 @simd for i ∈ 2:N-1 
        vx = vy =0.
        @views vx1 = (Y₀[i,j]-Y₀[i-1,j])/h
        @views vx2 = (Y₀[i+1,j]-Y₀[i,j])/h
        @views vy1 = (Y₀[i,j]-Y₀[i,j-1])/h
        @views vy2 = (Y₀[i,j+1]-Y₀[i,j])/h

    if v₀ >= 0.
        vx1 = max(vx1,0) 
        vx2 = -min(vx2,0)
    else
        vx1 = -min(vx1,0)
        vx2 = max(vx2,0)
    end

    if v₀ >= 0.
        vy1 = max(vy1,0)
        vy2 = -min(vy2,0)
    else
        vy1 = -min(vy1,0)
        vy2 = max(vy2,0)
    end
        nx = sqrt(max(vx1,vx2)^2+max(vy1,vy2)^2)
        F₀[i,j] = nx
    end end

end


# Φ+
function Φ₊!(N,Φ,dᵦ) #Phi Add
    @views for i ∈ 1:N 
        Φ[i] += dᵦ[i]
    end 
end


function Φ̄₊(Y₀,N,M,Ny) #Phi Average Add
    Φᵥ= 0.
    @inbounds for j ∈ 2:M-1 @simd for i ∈ 2:N-1
        @views Φᵥ += Y₀[i,j]
    end end
    return Φᵥ/(M-2)/(N-2)
end


## EF Add 
function Ef!(Y₀,ymid,Φ₀,F₀,dt,ν,ki₀,N,M)
    #for i from 2 to N-1 do for j from 2 to M-1 do ymid[i,j]:=Y₀[i,j]-dt*F₀[i,j]*Φ₀[i+(j-1)*N]:od:od:
    @inbounds for j ∈ 2:M-1 @simd for i ∈ 2:N-1  
        @views ymid[i,j] = max(eps(),Y₀[i,j]-dt*ν*ki₀*F₀[i,j]*Φ₀[i+(j-1)*N])
    end end

    @inbounds for i ∈ 1:N
        @views ymid[i,1] = ymid[i,2]
        @views ymid[i,M] = ymid[i,M-1]
    end

    @inbounds for j ∈ 1:M
        @views ymid[1,j] = ymid[2,j]
        @views ymid[N,j] = ymid[N-1,j]
    end
end


function Ef2!(Y₀,ymid,Φ₀,F₀,dt,ν,ki₀,N,M)
    @inbounds for j ∈ 2:M-1 @simd for i ∈ 2:N-1 
        @views ymid[i,j] = max(eps(),Y₀[i,j]*3/4+ymid[i,j]/4-dt/4*ν*ki₀*F₀[i,j]*Φ₀[i+(j-1)*N])
    end end

    @inbounds @views for i ∈ 1:N
         ymid[i,1] = ymid[i,2]
         ymid[i,M] = ymid[i,M-1]
    end

    @inbounds @views for j ∈ 1:M
         ymid[1,j] = ymid[2,j]
         ymid[N,j] = ymid[N-1,j]
    end
end


function Ef3!(Y₀,ymid,Φ₀,F₀,dt,ν,ki₀,N,M)
    @inbounds for j ∈ 2:M-1 @simd for i ∈ 2:N-1 
        @views Y₀[i,j] = max(eps(),Y₀[i,j]*1/3+ymid[i,j]*2/3-dt*2/3*ν*ki₀*F₀[i,j]*Φ₀[i+(j-1)*N])
    end end

        @inbounds @views for i ∈ 1:N
             Y₀[i,1] = Y₀[i,2]
             Y₀[i,M] = Y₀[i,M-1]
        end

        @inbounds @views for j ∈ 1:M
            Y₀[1,j] = Y₀[2,j]
            Y₀[N,j] = Y₀[N-1,j]
        end
end

function YStore(Y₀,Ydata,N,M,jj)
    @inbounds for j ∈ 2:M-1 @simd for i ∈ 2:N-1 
        @views Ydata[jj,i,j] = Y₀[i,j]
    end end
end


function PhiSwitch!(N,Φ₀,Φₜ)
    @inbounds @views for i ∈ 1:N
        Φₜ[i] = Φ₀[i]
    end
end


function MidPred!(N,Φ₀,Φₜ,Φₘ)
    @inbounds @views for i ∈ 1:N
        Φₘ[i] = (Φₜ[i]+Φ₀[i])/2
    end
end