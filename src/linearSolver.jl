struct LinearProblem <: AbstractEquationsProblem
    A::Matrix{Float64}
    b::Vector{Float64}
    maxiters::Int
    LinearProblem(A, b, maxiters=10000) = new(A, b, maxiters)
end

struct LinearProblem_functiontype <: AbstractEquationsProblem
    lp::LinearProblem
    vars::Vector{Num}
    LinearProblem_functiontype(lp::LinearProblem,vars::Vector{Num})=new(lp,vars)
end

function LinearProblem(eqs::Any, vars::Dict, maxiters=10000)
    checked_eqs = check_eqs(eqs)
    res = check_vars(checked_eqs, vars)
    dict = Dict(key => 0.0 for key in keys(vars))
    A::Matrix{Float64} = Symbolics.value.(Symbolics.jacobian(checked_eqs, collect(keys(vars))))
    b::Vector{Float64} = -Symbolics.value.(substitute.(checked_eqs, (dict,)))
    lp=LinearProblem(A,b,maxiters)
    return LinearProblem_functiontype(lp, collect(keys(vars)))
end

struct Direct <: AbstractLinearMethod end               #Julia自带方法
struct LUFactorized <: AbstractLinearMethod end         #已经进行LU分解的问题采用LU回带
struct ConjugateGradient <: AbstractLinearMethod end    #共轭梯度法，CG为缩写
struct CG <: AbstractLinearMethod end

#solve_liner = Symbolics.solve_for

#采用Julia自带的求解方法
function solve(problem::LinearProblem)
    problem.A \ problem.b
end
solve(problem::LinearProblem, ::Direct) = solve(problem::LinearProblem)
solve(problem::LinearProblem, ::Direct, lpiter::Int64) = solve(problem::LinearProblem)
function solve(lpf::LinearProblem_functiontype)
    Dict(zip(lpf.vars,solve(lpf.lp)))
end

#对于已经进行LU分解的问题采用LU分解的回带过程
function solve(problem::LinearProblem, ::LUFactorized, LUFactor::LU)
    LU_solve(LUFactor, problem.b)
end

#采用共轭梯度法，这里直接将矩阵A改为对称矩阵。如果A本来就对称，则没有影响
function solve(problem::LinearProblem, ::CG; abstol=1e-8)
    CG_solve(Symmetric(problem.A), problem.b; abstol=abstol)
end
@deprecate solve(problem::LinearProblem, ::ConjugateGradient; abstol=1e-8) solve(problem, CG(); abstol=abstol)