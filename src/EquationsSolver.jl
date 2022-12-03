module EquationsSolver

using Symbolics
using LinearAlgebra, SparseArrays
import Symbolics: @variables, Equation

include("base_preprocess.jl")
include("base_linearequations.jl")

abstract type AbstractSolveMethod end
abstract type AbstractEquationsProblem end

abstract type AbstractLinearMethod end
struct Direct <: AbstractLinearMethod end
struct CG <: AbstractLinearMethod end

abstract type AbstractNLMethod end
struct Newton <: AbstractNLMethod end

struct LinearProblem <: AbstractEquationsProblem
    A::Matrix{Float64}
    b::Vector{Float64}
    maxiters::Int
    LinearProblem(A,b,maxiters=10000)=new(A,b,maxiters)
end

struct NonlinearProblem <: AbstractEquationsProblem
    eqs::Vector{Num}
    vars::Dict
    maxiters::Int
    abstol::Float64
end

function LinearProblem(eqs::Any, vars::Dict, maxiters=10000)
    checked_eqs = check_eqs(eqs)
    res = check_vars(checked_eqs, vars)
    dict=Dict(key=>0.0 for key in keys(vars))
    A::Matrix{Float64} = Symbolics.value.(Symbolics.jacobian(checked_eqs, collect(keys(vars))))
    b::Vector{Float64} = -Symbolics.value.(substitute.(checked_eqs, (dict,)))
    return LinearProblem(A, b, maxiters), keys(vars)
end

function NLProblem(eqs::Any, vars::Dict; maxiters=10000, abstol=1.0E-6)
    eqs = check_eqs(eqs)
    res = check_vars(eqs, vars)
    vars = promote_vars(vars)
    return NonlinearProblem(eqs, vars, maxiters, abstol)
end

include("linearSolver.jl")
include("nonlinearSolver.jl")

export LinearProblem
export NLProblem
export solve
export @variables, Equation

end
