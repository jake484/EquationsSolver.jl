module EquationsSolver

using Symbolics
import Symbolics: @variables, Equation

include("base.jl")

abstract type EquationsProblem end

struct LinearProblem <: EquationsProblem
    eqs::Vector{Equation}
    vars::Vector{Num}
    maxiters::Int
end

struct NonlinearProblem <: EquationsProblem
    eqs::Vector{Num}
    vars::Dict
    maxiters::Int
    abstol::Float64
end

function LinearProblem(eqs::Any, vars::Dict, maxiters=10000)
    checked_eqs = check_eqs(eqs)
    res = check_vars(checked_eqs, vars)
    return LinearProblem(eqs, collect(keys(vars)), maxiters)
end

function NLProblem(eqs::Any, vars::Dict; maxiters=10000, abstol=1.0E-6)
    eqs = check_eqs(eqs)
    res = check_vars(eqs, vars)
    vars = promote_vars(vars)
    return NonlinearProblem(eqs, vars, maxiters, abstol)
end

export LinearProblem
export NLProblem
export solve
export @variables, Equation
include("solver.jl")
end
