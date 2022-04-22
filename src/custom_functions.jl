## Discretisation Parameters
Nx = NN = MM = Ny =  100
δ = 0.1
tf =2.
ν = 1.
ki₀ = 1.

N = NN+2
h = 1/NN
M = MM+2

Nₜ = N*M


## Model 1 Geometry
function y₀!(NN,MM,Y₀)
    N = NN+2
    h = 1.0/NN
    w = h/2
    M = Ny = MM+2
    for i ∈ 2:N-1 for j ∈ 2:M-1
        rr = ((i-3/2)*h)^2+((j-3/2)*h)^2 #combining two variables
        Y₀[i,j] = max(eps(),0.5*tanh(sqrt(2)*(sqrt(rr)-0.3)/w)+0.5) #eps() replaces 1e-9

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


Y₀ = Array{Float64}(undef,N,M)
eq = Array{Float64}(undef,N,M)
F₀ = Array{Float64}(undef,N,M)
y₀(NN,MM,Y₀)
Φ₀ = Vector{Float64}(undef,Nₜ)
ff = copy(Φ₀)
j₀ = Array{Float64}(undef,Nₜ,Nₜ) #Modify to be Sparse

#if delta<0 then read("Y₀data.m"):end: ??
#if delta<0 then read("tdata.m"):end: ??

plotly()
surface(Y₀,camera = (0, 90))

## System Eqs
function Eqs11!(Y₀,y,δ,ki₀,ff,N,M)
    h=1.0/(N-2)

    for i ∈ 2:N-1 for j ∈ 2:M-1
        i1 = i+(j-1)*N
        ff[i1] = (Y₀[i,j]+Y₀[i,j+1])*(y[i1+N]-y[i1]) -
        (Y₀[i,j]+Y₀[i,j-1])*(y[i1]-y[i1-N]) +
        (Y₀[i+1,j]+Y₀[i,j])*(y[i1+1]-y[i1]) -
        (Y₀[i,j]+Y₀[i-1,j])*(y[i1]-y[i1-1]) -
        ki₀*y[i1]*h*(1e-24+(Y₀[i+1,j]-Y₀[i-1,j])^2+(Y₀[i,j+1]-Y₀[i,j-1])^2)^(1/2)
    end end

    for i ∈ 1:1 for j ∈ 1:M
        i1 = i+(j-1)*N
        ff[i1] = -y[i1]+y[i1+1]
    end end

    for i ∈ N:N for j ∈ 1:M
        i1 = i+(j-1)*N
        ff[i1] = -y[i1]+y[i1-1]
    end end

    for i ∈ 1:N for j ∈ 1:1
        i1 = i+(j-1)*N
        ff[i1] = -y[i1]+y[i1+N]
    end end

    for i ∈ 1:N for j ∈ M:M
        i1 = i+(j-1)*N
        ff[i1] = -y[i1]+y[i1-N]+δ*h
    end end
end

## Jacobian
function Jac!(Y₀,y,δ,ki₀,j₀)

    for i ∈ 2:N-1 for j ∈ 2:M-1
        i1 = i+(j-1)*N
        j₀[i1,i1] = -4*Y₀[i,j]-Y₀[i,j+1]-Y₀[i,j-1]-Y₀[i+1,j]-Y₀[i-1,j]-
        ki₀*h*(1e-24+(Y₀[i+1,j]-Y₀[i-1,j])^2+(Y₀[i,j+1]-Y₀[i,j-1])^2)^(1/2)
        j₀[i1,i1+1] = Y₀[i+1,j]+Y₀[i,j]
        j₀[i1,i1-1] = Y₀[i-1,j]+Y₀[i,j]
        j₀[i1,i1+N] = Y₀[i,j+1]+Y₀[i,j]
        j₀[i1,i1-N] = Y₀[i,j-1]+Y₀[i,j]
    end end

    for i ∈ 1:1 for j ∈ 1:M
        i1 = i+(j-1)*N
        j₀[i1,i1] = -1.
        j₀[i1,i1+1] =1.
        #ff[i1] = -y[i1]+y[i1+1]
    end end

    for i ∈ N:N for j ∈ 1:M
        i1 = i+(j-1)*N
        j₀[i1,i1] = -1.
        j₀[i1,i1-1] = 1.
        #ff[i1]:=-y[i1]+y[i1-1]:
    end end

    for i ∈ 1:N for j ∈ 1:1
        i1 = i+(j-1)*N
        #ff[i1] = -y[i1]+y[i1+N]
        j₀[i1,i1] = -1.
        j₀[i1,i1+N] =1.
    end end

    for i ∈ 1:N for j ∈ M:M
        i1 = i+(j-1)*N
        j₀[i1,i1] = -1.
        j₀[i1,i1-N] = 1.
        #ff[i1] = -y[i1]+y[i1-N]+delta*h
    end end 
end


## Upwind
function Upwind!(Y₀,Φ₀,F₀,dₜ,M,v₀)
    h=1/(N-2)

    vv₀ = v₀
    for i ∈ 1:N
        Y₀[i,1] = Y₀[i,2]
        Y₀[i,M] = Y₀[i,M-1]
    end

    for j ∈ 1:M
        Y₀[1,j] = Y₀[2,j]
        Y₀[N,j] = Y₀[N-1,j]
    end

    for i ∈ 2:N-1 for j ∈ 2:M-1   
        vx = vy =0.
        vx1 = (Y₀[i,j]-Y₀[i-1,j])/h
        vx2 = (Y₀[i+1,j]-Y₀[i,j])/h
        vy1 = (Y₀[i,j]-Y₀[i,j-1])/h
        vy2 = (Y₀[i,j+1]-Y₀[i,j])/h

    if v0 >= 0.
        vx1 = max(vx1,0) 
        vx2 = -min(vx2,0)
    else
        vx1 = -min(vx1,0)
        vx2 = max(vx2,0)
    end

    if v0 >= 0.
        vy1 = max(vy1,0)
        vy2 = -min(vy2,0)
    else
        vy1 = -min(vy1,0)
        vy2 = max(vy2,0)
    end
        nx = sqrt(max(vx1,vx2)^2+max(vy1,vy2)^2)
        F0[i,j] = nx
    end end

end



# Φ+
function Φ₊!(N,Φ₀,dᵦ)
    for i ∈ 1:N 
        Φ₀[i] += dᵦ[i]
    end 
end

#Calls
Eqs11!(Y₀,Φ₀,δ,ki₀,ff,N,M)
Jac!(Y₀,Φ₀,δ,ki₀,j₀)

dᵦ = j₀\-ff



