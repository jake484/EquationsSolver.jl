struct NonlinearProblem <: AbstractEquationsProblem
    f::Function
    guessValue::Vector{Float64}
    maxiters::Int
    abstol::Float64
    NonlinearProblem(f::Function,guessValue::Vector{Float64},maxiters::Int,abstol::Float64)=new(f,guessValue,maxiters,abstol)
end

struct NonlinearProblem_functiontype <: AbstractEquationsProblem
    NP::NonlinearProblem
    vars::Vector{Num}
    NonlinearProblem_functiontype(NP::NonlinearProblem,vars::Vector{Num})=new(NP,vars)
end

function NLProblem(eqs::Any, vars::Dict; maxiters=10000, abstol=1.0E-6)
    eqs = check_eqs(eqs)
    res = check_vars(eqs, vars)
    vars = promote_vars(vars)
    guessValue=collect(values(vars))
    vars=collect(keys(vars))

    f_expr = build_function(eqs, vars)
    Base.remove_linenums!.(f_expr)
    f = eval(f_expr[1])

    NP=NonlinearProblem(f,guessValue,maxiters,abstol)

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


function solve(problem::NonlinearProblem,::Newton;lpMethod::AbstractLinearMethod=Direct(),lpiter::Int64=10000)
    NLNewton(problem.f, problem.jacobian, problem.guessValue, problem.maxiters,problem.abstol, lpMethod, lpiter)
end
solve(problem::NonlinearProblem)=solve(problem,Newton)

function solve(problem::NonlinearProblem_functiontype,method::AbstractNLMethod;lpMethod::AbstractLinearMethod=Direct(),lpiter::Int64=10000)
    res=solve(problem.NP,method;lpMethod=lpMethod,lpiter=lpiter)
    return Dict(zip(problem.vars,res))
end
solve(problem::NonlinearProblem_functiontype)=solve(problem,Newton())