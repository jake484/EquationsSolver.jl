using LinearAlgebra

A=rand(10,10)+5*I
r=rand(10)
n=10
m=5

function ArnoldiProcess(A::Matrix{Float64},r::Vector{Float64},n::Int64,m::Int64)
    V=Array{Float64,2}(undef,n,m+1)
    V[:,1]=r/norm(r)
    H=UpperHessenberg(Array{Float64}(undef,m+1,m))
    v_next=Vector{Float64}(undef,n)
    for k=1:m-1
        v_next = A*V[:,k]
        for i=1:k
            H[i,k]=dot(v_next, V[:,i])
        end
        for i=1:k
            v_next -= H[i,k] * V[:,i]
        end
        H[k+1,k] = norm(v_next)
        V[:,k+1] = v_next / H[k+1,k]
    end
    for i=1:m
        H[i,m]=dot(A * V[:,m], V[:,i])
    end
    v_next = A*V[:,m]
    for i=1:m
        v_next -= H[i,m] * V[:,i]
    end
    H[m+1,m] = norm(v_next)
    V[:,m+1] = (H[m+1,m]==0) ? zeros(Float64, n) : v_next/H[m+1,m]
    return V, H
end

V,H=ArnoldiProcess(A,r,n,m)