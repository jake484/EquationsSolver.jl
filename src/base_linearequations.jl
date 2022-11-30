###LU分解求线性方程组的根
"""
LU_factorization(A; sparsed=false)

输入矩阵A，返回矩阵A的LUP分解L,U,p

使得A[p,:]=L*U

返回 L,U,p
"""
function LU_factorization(A; sparsed=false)
    #下三角L，Doolittle分解
    L, U, p = lu(A)
    return sparsed ? (sparsed(L), sparsed(U), p) : (L, U, p)
end

"""
LU_solve(A, b)

使用LU分解法求解方程Ax=b

返回x
"""
function LU_solve(A, b)
    n = size(A, 1)
    _, U0, _ = LU_factorization(hcat(A, b))
    U = U0[:, 1:end-1]
    y = U0[:, end]
    x = zeros(n)
    x[n] = y[n] / U[n, n]
    for i = n-1:-1:1
        x[i] = (y[i] - U[i, i+1:n]' * x[i+1:n]) / U[i, i]
    end
    return x
end

"""
LU_solve(L, U, p, b)

使用LU分解法求解方程Ax=b

如果已经将矩阵A完成LUP分解，即PA=LU或A[p,:]=L*U

那么该方法可直接使用分解结果计算Ax=b的解

返回x
"""
function LU_solve(L, U, p, b)
    b = b[p]
    n = size(L, 1)
    #y=inv(L)*b
    y = zeros(n)
    y[1] = b[1]
    x = zeros(n)
    for i = 2:n
        y[i] = b[i] - L[i, 1:i-1]' * y[1:i-1]
    end
    x[n] = y[n] / U[n, n]
    for i = n-1:-1:1
        x[i] = (y[i] - U[i, i+1:n]' * x[i+1:n]) / U[i, i]
    end
    return x
end
### end LU

### 共轭梯度法 Conjugate_Gradient
"""
CG_solve(A, b; ep=1e-5)

当A为n×n对称正定矩阵时,可以采用共轭梯度法求解Ax=b的解

终止条件为残余矢量的模小于ep或迭代次数不小于4*n

对于良态矩阵或忽略摄入误差的情况下，一般迭代n次即可找到解

返回x
"""
function CG_solve(A, b; ep=1e-8)
    n = size(A, 1)
    x = zeros(n)
    r = b
    d = r
    flag = 0
    temp=sum(abs2.(r))
    while (temp > ep) & (flag<n*4)
        alpha = temp / (d' * A * d)
        x = x + alpha * d
        r=b-A*x
        beta=sum(abs2.(r))/temp
        temp=sum(abs2.(r))
        d=r+beta*d
        flag+=1
    end
    return x, flag
end
### end 共轭梯度法


function linearSolve(A, b, method::Symbol=:LU)

end