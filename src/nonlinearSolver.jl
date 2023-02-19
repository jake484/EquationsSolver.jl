struct NonlinearProblem <: AbstractEquationsProblem
    f::Function
    guessValue::Vector{Float64}
    NonlinearProblem(f::Function,guessValue::Vector{Float64})=new(f,guessValue)
end

struct NonlinearProblem_functiontype <: AbstractEquationsProblem
    NP::NonlinearProblem
    vars::Vector{Num}
    NonlinearProblem_functiontype(NP::NonlinearProblem,vars::Vector{Num})=new(NP,vars)
end

function NonlinearProblem(eqs::Any, vars::Dict)
    eqs = check_eqs(eqs)
    res = check_vars(eqs, vars)
    vars = promote_vars(vars)
    guessValue=collect(values(vars))
    vars=collect(keys(vars))
    f_expr = build_function(eqs, vars)
    Base.remove_linenums!.(f_expr)
    f = eval(f_expr[1])

    NP=NonlinearProblem(f,guessValue)

    return NonlinearProblem_functiontype(NP, vars)
end

function Base.getproperty(NP::NonlinearProblem,d::Symbol)
    if d==:jacobian
        N=length(NP.guessValue)
        @variables x[1:N]
        x=collect(x)
        ja_symbolic=Symbolics.jacobian(NP.f(x), x)
        ja_expr = build_function(ja_symbolic, x)
        Base.remove_linenums!.(ja_expr)
        ja=eval(ja_expr[2])
        return ja
    else
        return getfield(NP,d)
    end
end


struct Newton <: AbstractNLMethod end
struct NLJacobian <: AbstractNLMethod end #简单迭代，simple iteration method
struct SimplifiedNewton <: AbstractNLMethod end
struct Secant <: AbstractNLMethod end
struct BroydenRank1 <: AbstractNLMethod end

# 牛顿法
function solve(problem::NonlinearProblem,::Newton; maxiter=100, abstol=1e-6, lpMethod::AbstractLinearMethod=Direct(),lpSolverSetting::Dict{Symbol,T}=Dict{Symbol,Real}()) where {T<:Real}
    x,err=NLNewton(problem.f, problem.jacobian,problem.guessValue, maxiter, abstol; lpMethod=lpMethod, lpSolverSetting=lpSolverSetting)
    if err>abstol
        @warn string("Unsuccessful solve, absolute error = ",err,".")
    end
    return x
end
solve(problem::NonlinearProblem)=solve(problem,Newton())

# 对符号方程类型求解
function solve(problem::NonlinearProblem_functiontype,method::AbstractNLMethod;lpMethod::AbstractLinearMethod=Direct(),lpSolverSetting::Dict{Symbol,T}=Dict{Symbol,Real}()) where {T<:Real}
    res=solve(problem.NP,method;lpMethod=lpMethod,lpSolverSetting=lpSolverSetting)
    return Dict(zip(problem.vars,res))
end
solve(problem::NonlinearProblem_functiontype)=solve(problem,Newton())

function solve(problem::NonlinearProblem, ::NLJacobian; maxiter=1000, abstol=1e-8)
    x,err = NLJacobianIter(problem.f, problem.guessValue, maxiter, abstol)
    if err>abstol
        @warn string("Unsuccessfully solved, absolute error = ",err,". Try larger iteration times, or the iteration form may not be convergence.")
    end
    return x
end

function solve(problem::NonlinearProblem,::SimplifiedNewton; maxiter=100, abstol=1e-6)
    jaValue = ForwardDiff.jacobian(problem.f,problem.guessValue)
    x,err=NLSimplifiedNewton(problem.f, problem.guessValue, jaValue, maxiter, abstol)
    if err>abstol
        @warn string("Unsuccessful solve, absolute error = ",err,".")
    end
    return x
end

function solve(problem::NonlinearProblem,::Secant;h::Float64=1e-3, maxiter=100, abstol=1e-6, lpMethod::AbstractLinearMethod=Direct(),lpSolverSetting::Dict{Symbol,T}=Dict{Symbol,Real}()) where {T<:Real}
    x,err=NLSecant(problem.f, problem.guessValue,h, maxiter, abstol;lpMethod=lpMethod, lpSolverSetting=lpSolverSetting)
    if err>abstol
        @warn string("Unsuccessful solve, absolute error = ",err,".")
    end
    return x
end

function solve(problem::NonlinearProblem,::BroydenRank1; maxiter=100, abstol=1e-6)
    A = inv(ForwardDiff.jacobian(problem.f,problem.guessValue))
    x,err=NLBroydenR1(problem.f, problem.guessValue,A, maxiter, abstol)

    if err>abstol
        @warn string("Unsuccessful solve, absolute error = ",err,".")
    end
    return x
end


