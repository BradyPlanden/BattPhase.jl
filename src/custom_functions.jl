## Model 1 Geometry (Semi-Circle)
function y₀1!(NN,MM,Y₀,κ)
    N = NN+2
    h = 1.0/NN
    w = h/2
    M = MM+2
    @views for j ∈ 2:M-1 for i ∈ 2:N-1
        rr = ((i-3/2)*h)^2+((j-3/2)*h)^2 #combining two variables
        Y₀[i,j] = max(eps(),0.5*tanh(√2*(√rr-κ)/w)+0.5) # κ = 0.3

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

# Model 4 Geometry (Gauss Seed)
function y₀4!(NN,MM,Y₀)
    N = NN+2
    h = 1.0/NN
    w = h/2
    M = MM+2
    @views for j ∈ 2:M-1 for i ∈ 2:N-1
        xₓ = (i-3/2)*h
        y = (j-3/2)*h
        rᵣ₁ = 0.4+0.3*exp(-100*(xₓ-0.5)^2)
        rᵣ₂ = 0.4+0.2*exp(-300*(xₓ-0.2)^2)
        rᵣ₃ = 0.4+0.1375*exp(-300*(xₓ-0.8)^2)
        Y₀[i,j] = max(eps(),0.5+min(0.5*tanh((y-rᵣ₁)/w/√2),0.5*tanh((y-rᵣ₂)/w/√2),0.5*tanh((y-rᵣ₃)/w/√2)))#,0.5+0.5*tanh((y-rᵣ₂)/w/√2))

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

    @inbounds @views for j ∈ 2:M-1 for i ∈ 2:N-1 
        i1 = i+(j-1)*N
        ff[i1] = (Y₀[i,j]+Y₀[i,j+1])*(y[i1+N]-y[i1]) -
        (Y₀[i,j]+Y₀[i,j-1])*(y[i1]-y[i1-N]) +
        (Y₀[i+1,j]+Y₀[i,j])*(y[i1+1]-y[i1]) -
        (Y₀[i,j]+Y₀[i-1,j])*(y[i1]-y[i1-1]) -
        ki₀*y[i1]*h*(eps()+(Y₀[i+1,j]-Y₀[i-1,j])^2 +
        (Y₀[i,j+1]-Y₀[i,j-1])^2)^(1/2)
    end end


    ## Boundary Conditions
    @inbounds @views for j ∈ 1:M
        i1 = 1+(j-1)*N
        ff[i1] = -y[i1]+y[i1+1]
    end

    @inbounds @views for j ∈ 1:M
        i1 = N+(j-1)*N
        ff[i1] = -y[i1]+y[i1-1]
    end

    @inbounds @views for i ∈ 1:N 
        ff[i] = -y[i]+y[i+N]
    end

    @inbounds @views for i ∈ 1:N 
        i1 = i+(M-1)*N
        ff[i1] = -y[i1]+y[i1-N]+abs(δ)*h
    end
end

## Jacobian
function Jac!(Y₀,ki₀,j₀,N,M,h)

    @inbounds @views for j ∈ 2:M-1 for i ∈ 2:N-1
        i1 = i+(j-1)*N
        j₀[i1,i1] = -4*Y₀[i,j]-Y₀[i,j+1]-Y₀[i,j-1]-Y₀[i+1,j]-Y₀[i-1,j]-
        ki₀*h*(eps()+(Y₀[i+1,j]-Y₀[i-1,j])^2+(Y₀[i,j+1]-Y₀[i,j-1])^2)^(1/2)
        j₀[i1,i1+1] = Y₀[i+1,j]+Y₀[i,j]
        j₀[i1,i1-1] = Y₀[i-1,j]+Y₀[i,j]
        j₀[i1,i1+N] = Y₀[i,j+1]+Y₀[i,j]
        j₀[i1,i1-N] = Y₀[i,j-1]+Y₀[i,j]
    end end

    @inbounds @views for j ∈ 1:M
        i1 = 1+(j-1)*N
        j₀[i1,i1] = -1.
        j₀[i1,i1+1] = 1.
    end

    @inbounds @views for j ∈ 1:M 
        i1 = N+(j-1)*N
        j₀[i1,i1] = -1.
        j₀[i1,i1-1] = 1.

    end

    @inbounds @views for i ∈ 1:N 
        j₀[i,i] = -1.
        j₀[i,i+N] = 1.
    end

    @inbounds @views for i ∈ 1:N 
        i1 = i+(M-1)*N
        j₀[i1,i1] = -1.
        j₀[i1,i1-N] = 1.
    end
end

## Upwind
function Upwind!(Y₀,F₀,N,M,δ,h)

    @inbounds @views for i ∈ 1:N
        Y₀[i,1] = Y₀[i,2]
        Y₀[i,M] = Y₀[i,M-1]
    end

    @inbounds @views for j ∈ 1:M
        Y₀[1,j] = Y₀[2,j]
        Y₀[N,j] = Y₀[N-1,j]
    end

    @inbounds @views for j ∈ 2:M-1 for i ∈ 2:N-1 
        vx1 = (Y₀[i,j]-Y₀[i-1,j])/h
        vx2 = (Y₀[i+1,j]-Y₀[i,j])/h
        vy1 = (Y₀[i,j]-Y₀[i,j-1])/h
        vy2 = (Y₀[i,j+1]-Y₀[i,j])/h

    if δ >= zero(δ)
        vx1 = max(vx1,zero(vx1))
        vx2 = -min(vx2,zero(vx2))
        vy1 = max(vy1,zero(vy1))
        vy2 = -min(vy2,zero(vy2))
    else
        vx1 = -min(vx1,zero(vx1))
        vx2 = max(vx2,zero(vx2))
        vy1 = -min(vy1,zero(vy1))
        vy2 = max(vy2,zero(vy2))
    end
        nx = sqrt(max(vx1,vx2)^2+max(vy1,vy2)^2)
        F₀[i,j] = nx
    end end

end

## ENO2
function ENO!(Y₀,F₀,N,M,δ,h)

    @inbounds @views for i ∈ 1:N
        Y₀[i,1] = Y₀[i,2]
        Y₀[i,M] = Y₀[i,M-1]
    end

    @inbounds @views for j ∈ 1:M
        Y₀[1,j] = Y₀[2,j]
        Y₀[N,j] = Y₀[N-1,j]
    end

    @inbounds @views for j ∈ 2:M-1 for i ∈ 2:N-1 
        sdx = (Y₀[i+1,j]-2*Y₀[i,j]+Y₀[i-1,j])/h
        sdy = (Y₀[i,j+1]-2*Y₀[i,j]+Y₀[i,j-1])/h

        if i == 2 
            sdxb = (Y₀[i,j]-2*Y₀[i-1,j]+Y₀[i-1,j])/h
        else
            sdxb = (Y₀[i,j]-2*Y₀[i-1,j]+Y₀[i-2,j])/h
        end

        if sdx*sdxb >= zero(sdx)
            s1x = 1.0 
        else 
            s1x = 0.0
        end

        vx1 = (Y₀[i,j]-Y₀[i-1,j])/h+0.5*sign(sdx)*s1x*min(abs(sdx),abs(sdxb))

        if i == N-1
            sdxf = (Y₀[i+1,j]-2*Y₀[i+1,j]+Y₀[i,j])/h
        else
            sdxf = (Y₀[i+2,j]-2*Y₀[i+1,j]+Y₀[i,j])/h
        end

        if sdx*sdxf >= zero(sdx)
            s1x = 1.0
        else 
            s1x = 0.0
        end

        vx2 = (Y₀[i+1,j]-Y₀[i,j])/h-0.5*sign(sdx)*s1x*min(abs(sdx),abs(sdxf))
        
        if j == 2
            sdyb = (Y₀[i,j]-2*Y₀[i,j-1]+Y₀[i,j-1])/h
        else
            sdyb = (Y₀[i,j]-2*Y₀[i,j-1]+Y₀[i,j-2])/h 
        end
        
        if sdy*sdyb >= zero(sdy)
            s1y = 1.0 
        else 
            s1y = 0.0
        end

        vy1 = (Y₀[i,j]-Y₀[i,j-1])/h+0.5*sign(sdy)*s1y*min(abs(sdy),abs(sdyb))

        if j == M-1
            sdyf = (Y₀[i,j+1]-2*Y₀[i,j+1]+Y₀[i,j])/h
        else
            sdyf = (Y₀[i,j+2]-2*Y₀[i,j+1]+Y₀[i,j])/h
        end

        if sdy*sdyf >= zero(sdy)
            s1y = 1.0
        else 
            s1y = 0.0
        end

        vy2 = (Y₀[i,j+1]-Y₀[i,j])/h+0.5*sign(sdy)*s1y*min(abs(sdy),abs(sdyf))

        if δ >= zero(δ)
            vx1 = max(vx1,zero(vx1))
            vx2 = -min(vx2,zero(vx2))
            vy1 = max(vy1,zero(vy1))
            vy2 = -min(vy2,zero(vy2))
        else
            vx1 = -min(vx1,zero(vx1))
            vx2 = max(vx2,zero(vx2))
            vy1 = -min(vy1,zero(vy1))
            vy2 = max(vy2,zero(vy2))
        end

        nx = sqrt(max(vx1,vx2)^2+max(vy1,vy2)^2)
        F₀[i,j] = nx
    end end
end


## WENO3
function WENO!(Y₀,F₀,N,M,δ,h)
    e1 = eps()
    @inbounds @views for i ∈ 1:N
        Y₀[i,1] = Y₀[i,2]
        Y₀[i,M] = Y₀[i,M-1]
    end

    @inbounds @views for j ∈ 1:M
        Y₀[1,j] = Y₀[2,j]
        Y₀[N,j] = Y₀[N-1,j]
    end

    @inbounds @views for i ∈ 2:N-1 for j ∈ 2:M-1
        phix = (Y₀[i+1,j]-Y₀[i-1,j])/2/h
        phiy = (Y₀[i,j+1]-Y₀[i,j-1])/2/h 

        if i == 2
            sdb =Y₀[i,j]-2*Y₀[i-1,j]+Y₀[i-1,j] 
        else 
            sdb = Y₀[i,j]-2*Y₀[i-1,j]+Y₀[i-2,j]
        end
        
        sd = Y₀[i+1,j]-2*Y₀[i,j]+Y₀[i-1,j]

        if i == N-1
            sdf = Y₀[i,j]-2*Y₀[i+1,j]+Y₀[i+1,j]
        else 
            sdf = Y₀[i+2,j]-2*Y₀[i+1,j]+Y₀[i,j]
        end

        r1 = (e1+sdb^2)/(e1+sd^2)
        w1 = 1/(1+2*r1^2)
        r2 = (e1+sdf^2)/(e1+sd^2)
        w2 = 1/(1+2*r2^2)
        vx1 = phix-0.5*w1/h*(sd-sdb)
        vx2 = phix-0.5*w2/h*(sdf-sd)

        if j == 2
            sdb = Y₀[i,j]-2*Y₀[i,j-1]+Y₀[i,j-1]
        else 
            sdb = Y₀[i,j]-2*Y₀[i,j-1]+Y₀[i,j-2]
        end

        sd = Y₀[i,j+1]-2*Y₀[i,j]+Y₀[i,j-1]

        if j == M-1
            sdf = Y₀[i,j]-2*Y₀[i,j+1]+Y₀[i,j+1]
        else
            sdf = Y₀[i,j+2]-2*Y₀[i,j+1]+Y₀[i,j]
        end

        r1 = (e1+sdb^2)/(e1+sd^2)
        w1 = 1/(1+2*r1^2)
        r2 = (e1+sdf^2)/(e1+sd^2)
        w2 = 1/(1+2*r2^2)
        vy1 = phiy-0.5*w1/h*(sd-sdb)
        vy2 = phiy-0.5*w2/h*(sdf-sd)

        if δ >= zero(δ)
            vx1 = max(vx1,zero(vx1))
            vx2 = -min(vx2,zero(vx2))
            vy1 = max(vy1,zero(vy1))
            vy2 = -min(vy2,zero(vy2))
        else
            vx1 = -min(vx1,zero(vx1))
            vx2 = max(vx2,zero(vx2))
            vy1 = -min(vy1,zero(vy1))
            vy2 = max(vy2,zero(vy2))
        end
        nx = sqrt(max(vx1,vx2)^2+max(vy1,vy2)^2)
        F₀[i,j] = nx
    end end
end


# Φ+
function Φ₊!(N,Φ,dᵦ) #Phi Add
    @inbounds @views for i ∈ 1:N
        Φ[i] += dᵦ[i]
    end 
end


function Φ̄₊(Y₀,N,M) #Phi Average Add
    Φᵥ= zero(Y₀[1])
    @inbounds @views for j ∈ 2:M-1 for i ∈ 2:N-1
        Φᵥ += Y₀[i,j]
    end end
    return Φᵥ/(M-2)/(N-2)
end


## EF Add 
function Ef!(Y₀,ymid,Φ₀,F₀,dt,ν,ki₀,N,M)
    @inbounds @views for j ∈ 2:M-1 for i ∈ 2:N-1  
        ymid[i,j] = max(eps(),Y₀[i,j]-dt*ν*ki₀*F₀[i,j]*Φ₀[i+(j-1)*N])
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


function Ef2!(Y₀,ymid,Φ₀,F₀,dt,ν,ki₀,N,M)
    @inbounds @views for j ∈ 2:M-1 for i ∈ 2:N-1 
        ymid[i,j] = max(eps(),Y₀[i,j]*3/4+ymid[i,j]/4-dt/4*ν*ki₀*F₀[i,j]*Φ₀[i+(j-1)*N])
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
    @inbounds @views for j ∈ 2:M-1  for i ∈ 2:N-1 
        Y₀[i,j] = max(eps(),Y₀[i,j]*1/3+ymid[i,j]*2/3-dt*2/3*ν*ki₀*F₀[i,j]*Φ₀[i+(j-1)*N])
    end end

        @inbounds @views  for i ∈ 1:N
            Y₀[i,1] = Y₀[i,2]
            Y₀[i,M] = Y₀[i,M-1]
        end

        @inbounds @views  for j ∈ 1:M
            Y₀[1,j] = Y₀[2,j]
            Y₀[N,j] = Y₀[N-1,j]
        end
end

function YStore(Y₀,Ydata,N,M,jj)
    @inbounds @views for j ∈ 2:M-1 for i ∈ 2:N-1 
        Ydata[jj,i,j] = Y₀[i,j]
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