
function LU_solve(F::LU, b::Vector{Float64})
    b = b[F.p]
    y = F.L \ b
    x = F.U \ y
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
function CG_solve(A::Symmetric{Float64}, b::Vector{Float64}, guessValue::Vector{Float64}, abstol::Float64=1e-8)
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
    return x, temp
end
### end 共轭梯度法

### GMRES
"""
generalized minimum residual (GMRES) algorithm
"""
function GMRES_restarted(A::AbstractMatrix{Float64}, b::Vector{Float64}, guessValue::Vector{Float64}, m::Int64, maxiter::Int64, abstol::Float64=1e-8)
    n = size(A, 1)
    error::Float64 = Inf
    iter::Int64 = 0
    r::Vector{Float64} = Vector{Float64}(undef, n)
    c::Vector{Float64} = Vector{Float64}(undef, m + 1)
    while error > abstol && iter < maxiter
        r = b - A * guessValue
        V, H = ArnoldiProcess(A, r, n, m)
        #=
        这段的优化空间：
        1. 对于Hessenberg矩阵直接用m次Givens变换，而不用自带的Householder分解
        2. 矩阵H计算后不再使用，可以直接在H的内存计算R
        3. 直接采用Julia自带的最小二乘问题求解替换QR分解求最小二乘问题是否更快
        4. 即便做QR分解，矩阵Q也不用全部求出，只需求出第一行
        5. 将H记录为m*m，最后一个元素单独记录是否能提高效率
        =#
        Q, R = qr(H)
        c = Q[1, :] * norm(r)
        y = R \ c[1:m]

        guessValue += V[:, 1:m] * y
        error = abs(c[m+1])
        iter += 1
    end
    return guessValue, error
end

function HessenbergQR()

end

function ArnoldiProcess(A::AbstractMatrix{Float64}, r::Vector{Float64}, n::Int64, m::Int64)
    V = Array{Float64,2}(undef, n, m + 1)
    V[:, 1] = r / norm(r)
    H = UpperHessenberg(Array{Float64}(undef, m + 1, m))
    v_next = Vector{Float64}(undef, n)
    for k = 1:m-1
        v_next = A * V[:, k]
        for i = 1:k
            H[i, k] = dot(v_next, V[:, i])
        end
        for i = 1:k
            v_next -= H[i, k] * V[:, i]
        end
        H[k+1, k] = norm(v_next)
        V[:, k+1] = v_next / H[k+1, k]
    end
    for i = 1:m
        H[i, m] = dot(A * V[:, m], V[:, i])
    end
    v_next = A * V[:, m]
    for i = 1:m
        v_next -= H[i, m] * V[:, i]
    end
    H[m+1, m] = norm(v_next)
    V[:, m+1] = (H[m+1, m] == 0) ? zeros(Float64, n) : v_next / H[m+1, m]
    return V, H
end

### end GMRES