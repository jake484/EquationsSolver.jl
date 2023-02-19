
# 牛顿法
function NLNewton(f::Function, ja::Function, guessValue::Vector{Float64}, maxiter::Int64, abstol::Float64; lpMethod::AbstractLinearMethod=Direct(), lpSolverSetting::Dict{Symbol,T}=Dict()) where {T<:Real}
    this_ja(jaValue::Matrix{Float64}, x::Vector{Float64}) = Base.invokelatest(ja, jaValue, x)
    this_f(x::Vector{Float64}) = Base.invokelatest(f, x)
    n = length(guessValue)
    # 存结果
    newValue::Vector{Float64} = guessValue[:]
    lastValue::Vector{Float64} = Inf * ones(n)
    # 存函数值
    eqsValue::Vector{Float64} = zeros(n)
    # 存Jacobian矩阵
    jaValue::Matrix{Float64} = zeros(n, n)

    count::Int64 = 0

    ep::Float64 = Inf

    while count < maxiter && ep > abstol
        lastValue = newValue
        this_ja(jaValue, newValue)
        eqsValue = this_f(newValue)
        lp = LinearProblem(jaValue, -eqsValue)
        dx = solve(lp, lpMethod; lpSolverSetting...)
        newValue += dx
        ep = maximum(abs.(dx))
        count += 1
    end

    return newValue, ep
end

function NLJacobianIter(f::Function, guessValue::Vector{Float64}, maxiter::Int64, abstol::Float64)
    this_f(x::Vector{Float64}) = Base.invokelatest(f, x) + x
    iter::Int64 = 0
    error::Float64 = Inf
    while error > abstol && iter < maxiter
        temp = guessValue[:]
        guessValue = this_f(guessValue)
        error = maximum(abs.(temp - guessValue))
        iter += 1
    end
    return guessValue, error
end

# 对方程组简单迭代法加速后和简化牛顿法相同
function NLSimplifiedNewton(f::Function, guessValue::Vector{Float64},jaValue::Matrix{Float64}, maxiter::Int64, abstol::Float64)
    this_f(x::Vector{Float64}) = Base.invokelatest(f, x)
    n = length(guessValue)
    # 存结果
    newValue::Vector{Float64} = guessValue[:]
    lastValue::Vector{Float64} = Inf * ones(n)
    # 存函数值
    eqsValue::Vector{Float64} = zeros(n)
    # 存Jacobian矩阵
    F = lu(jaValue)

    count::Int64 = 0

    ep::Float64 = Inf

    while count < maxiter && ep > abstol
        lastValue = newValue
        eqsValue = this_f(newValue)
        dx = LU_solve(F, -eqsValue)
        newValue += dx
        ep = maximum(abs.(dx))
        count += 1
    end

    if count >= maxiter
        error("Error: iterate times exceed maxiter")
    end

    return newValue, ep
end

# 弦割法
function NLSecant(f::Function, guessValue::Vector{Float64}, h::Float64, maxiter::Int64, abstol::Float64; lpMethod::AbstractLinearMethod=Direct(), lpSolverSetting::Dict{Symbol,T}=Dict()) where {T<:Real}
    this_f(x::Vector{Float64}) = Base.invokelatest(f, x)
    n::Int64 = length(guessValue)
    count::Int64 = 0
    A::Matrix{Float64} = zeros(n, n)
    newValue::Vector{Float64} = guessValue[:]
    lastValue::Vector{Float64} = Inf * ones(n)

    while count < maxiter && h > abstol
        lastValue = newValue

        newValue[1] += h
        A[:, 1] = this_f(newValue)
        for i = 2:n
            newValue[i] += h
            newValue[i-1] -= h
            A[:, i] = this_f(newValue)
        end
        newValue[end] -= h

        lp = LinearProblem(A, this_f(lastValue))
        z = solve(lp, lpMethod; lpSolverSetting...)
        dx = h * z / (sum(z) - 1)
        newValue += dx
        h = maximum(abs.(dx))

        count += 1
    end

    return newValue, h
end

function NLBroydenR1(f::Function, guessValue::Vector{Float64}, A::Matrix{Float64}, maxiter::Int64, abstol::Float64)
    this_f(x::Vector{Float64}) = Base.invokelatest(f, x)
    n::Int64 = length(guessValue)
    count::Int64 = 0

    lastValue::Vector{Float64} = guessValue[:]

    temp1::Vector{Float64} = this_f(lastValue)
    dx::Vector{Float64} = -A * temp1
    newValue::Vector{Float64} = lastValue + dx
    temp2::Vector{Float64} = this_f(newValue)
    y::Vector{Float64} = temp2 - temp1
    err::Float64 = Inf

    while count < maxiter && err > abstol
        A += (dx - A * y) * dx' * A / (dx' * A * y)
        dx = -A * temp2
        err = maximum(abs.(dx))
        lastValue = newValue
        newValue += dx
        temp1 = temp2
        temp2 = this_f(newValue)
        y = temp2 - temp1
    end

    return newValue, err
end
