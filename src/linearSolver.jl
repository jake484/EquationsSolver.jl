"""
LinearProblem type has 
"""
struct LinearProblem <: AbstractEquationsProblem
    A::AbstractMatrix{Float64}
    b::AbstractVector{Float64}
    guessValue::Vector{Float64}
    LinearProblem(A::AbstractMatrix{Float64}, b::AbstractVector{Float64},guessValue::Vector{Float64}) = new(A, b, guessValue)
end

LinearProblem(A::Matrix{Float64}, b::Vector{Float64})=LinearProblem(A, b,zeros(length(b)))

struct LinearProblem_functiontype <: AbstractEquationsProblem
    lp::LinearProblem
    vars::Vector{Num}
    LinearProblem_functiontype(lp::LinearProblem,vars::Vector{Num})=new(lp,vars)
end

"""
Form a LinearProblem_functiontype with equations
"""
function LinearProblem(eqs::Any, vars::Dict)
    checked_eqs = check_eqs(eqs)
    res = check_vars(checked_eqs, vars)
    dict = Dict(key => 0.0 for key in keys(vars))
    A::Matrix{Float64} = Symbolics.value.(Symbolics.jacobian(checked_eqs, collect(keys(vars))))
    b::Vector{Float64} = -Symbolics.value.(substitute.(checked_eqs, (dict,)))
    lp=LinearProblem(A,b)
    return LinearProblem_functiontype(lp, collect(keys(vars)))
end

struct Direct <: AbstractLinearMethod end               #Julia自带方法
struct LUFactorized <: AbstractLinearMethod end         #已经进行LU分解的问题采用LU回带
struct ConjugateGradient <: AbstractLinearMethod end    #共轭梯度法，CG为缩写
struct CG <: AbstractLinearMethod end
struct GMRESM <: AbstractLinearMethod end


"""
Solve a LinearProblem with method '\\' originally in Julia
"""
function solve(problem::LinearProblem)
    problem.A \ problem.b
end
solve(problem::LinearProblem, ::Direct) = solve(problem::LinearProblem)
solve(problem::LinearProblem, ::Direct, lpiter::Int64) = solve(problem::LinearProblem)
function solve(lpf::LinearProblem_functiontype)
    Dict(zip(lpf.vars,solve(lpf.lp)))
end

"""
Solve a LinearProblem whose coefficient matrix has been LUP factorized
"""
function solve(problem::LinearProblem, ::LUFactorized, LUFactor::LU)
    LU_solve(LUFactor, problem.b)
end

"""
Use Conjugate Gradient method to solve LinearProblem Ax=b where A is symmetric positive definit matrix.
This method is very efficient, thus is recommended if A suit the case.
"""
function solve(problem::LinearProblem, ::CG; abstol=1e-8)
    CG_solve(Symmetric(problem.A), problem.b, problem.guessValue, abstol)
end
@deprecate solve(problem::LinearProblem, ::ConjugateGradient; abstol=1e-8) solve(problem, CG(); abstol=abstol)

function  solve(problem::LinearProblem, ::GMRESM, m=10; maxiter=1000, abstol=1e-8)
    x, error=GMRES_restarted(problem.A, problem.b, problem.guessValue,m , maxiter, abstol)
    if error>abstol
        @warn string("Unsuccessful solve, error = ",error,". Try larger m or larger iteration times.")
    end
    return x
end
