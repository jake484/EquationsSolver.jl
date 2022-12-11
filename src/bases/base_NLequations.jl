
# 牛顿法
function NLNewton(f::Function, ja::Function, guessValue::Vector{Float64}, maxiter, abstol, lpMethod::AbstractLinearMethod=Direct(),lpiter::Int64=10000)
    this_ja(jaValue::Matrix{Float64}, x::Vector{Float64})=Base.invokelatest(ja,jaValue,x)
    this_f(x::Vector{Float64})=Base.invokelatest(f,x)
    n=length(guessValue)
    # 存结果
    newValue::Vector{Float64} = guessValue[:]
    lastValue::Vector{Float64} = Inf*ones(n)
    # 存函数值
    eqsValue::Vector{Float64} = zeros(n)
    # 存Jacobian矩阵
    jaValue::Matrix{Float64} = zeros(n,n)

    count::Int64 = 0

    while count < maxiter && maximum(abs.(newValue .- lastValue)) > abstol
        lastValue = newValue
        this_ja(jaValue, newValue)
        eqsValue=this_f(newValue)
        lp=LinearProblem(jaValue,-eqsValue)
        newValue += solve(lp,lpMethod,lpiter)
        count += 1
    end

    if count >= maxiter
        error("Error: iterate times exceed maxiter")
    end

    return newValue
end

function NLJacobianIter(f::Function, guessValue::Vector{Float64}, maxiter, abstol)
    this_f(x::Vector{Float64})=Base.invokelatest(f,x)
    temp=guessValue[:]
    guessValue=this_f(guessValue)
    iter=0
    while maximum(abs.(temp-guessValue))>abstol && iter<maxiter
        guessValue=f(guessValue)
        iter+=1
    end
    return newValue
end