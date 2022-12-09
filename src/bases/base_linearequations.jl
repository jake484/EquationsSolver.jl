
function LU_solve(F::LU, b::Vector)
    b = b[F.p]
    n = size(F.factors, 1)
    #y=inv(L)*b
    y::Vector{Float64} = zeros(n)
    y[1] = b[1]
    x = zeros(n)
    for i = 2:n
        y[i] = b[i]
        for j = 1:i-1
            y[i] -= F.factors[i, j] * y[j]
        end
    end
    x[n] = y[n] / F.factors[n, n]
    for i = n-1:-1:1
        x[i] = y[i]
        for j = i+1:n
            x[i] -= F.factors[i, j] * x[j]
        end
        x[i] /= F.factors[i, i]
    end
    return x
end
### end LU

### 共轭梯度法 Conjugate_Gradient
"""
CG_solve(A::Symmetric{Float64}, b::Vector{Float64}, guessValue::Vector{Float64}, abstol=1e-8)

当A为n×n对称正定矩阵时,可以采用共轭梯度法求解Ax=b的解

终止条件为残余矢量的模小于ep或迭代次数不小于4*n

对于良态矩阵或忽略摄入误差的情况下，一般迭代n次即可找到解

对稀疏矩阵，迭代次数常常小于n

返回x
"""
function CG_solve(A::Symmetric{Float64}, b::Vector{Float64}, guessValue::Vector{Float64}, abstol=1e-8)
    n = size(A, 1)
    x = guessValue
    r = b
    d = r
    flag = 0
    temp = sum(abs2.(r))
    while (temp > abstol) & (flag < n * 4)
        alpha = temp / (d' * A * d)
        x = x + alpha * d
        r = b - A * x
        beta = sum(abs2.(r)) / temp
        temp = sum(abs2.(r))
        d = r + beta * d
        flag += 1
    end
    return x
end
### end 共轭梯度法

### GMRES
"""
generalized minimum residual (GMRES) algorithm
"""
function GMRES_restarted(A::Matrix{Float64}, b::Vector{Float64}, guessValue::Vector{Float64},m, abstol=1e-8)
    n = size(A, 1)
    r=b-A*guessValue
    beta=norm(r)
    V,H=ArnoldiProcess(A,r,n,m)
    
    return x
end

function HessenbergQR()
    
end

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

### end GMRES