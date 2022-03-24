module EquationsSolver

using Symbolics

include("base.jl")

abstract type EquationsProblem end
struct LinearProblem <: EquationsProblem
    eqs::Vector{Equation}
    vars::Vector{Num}
    maxiters::Int
end

struct NonlinearProblem <: EquationsProblem
    eqs::Vector{Num}
    vars::Vector{Num}
    maxiters::Int
end

function LinearProblem(eqs::Any, vars::Dict, maxiters=10000)
    checked_eqs = check_eqs(eqs)
    res = check_vars(checked_eqs, vars)
    return LinearProblem(eqs, collect(keys(vars)), maxiters)
end
function NonlinearProblem(eqs::Any, vars::Dict, maxiters=10000)
    eqs = check_eqs(eqs)
    res = check_vars(eqs, vars)
    return LinearProblem(eqs, vars, maxiters)
end

export LinearProblem
export NonlinearProblem
export solve
include("solver.jl")
end
