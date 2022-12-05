module EquationsSolver

using Symbolics
using LinearAlgebra, SparseArrays
import Symbolics: @variables, Equation

abstract type AbstractSolveMethod end
abstract type AbstractEquationsProblem end

abstract type AbstractLinearMethod <: AbstractSolveMethod end
abstract type AbstractNLMethod <: AbstractSolveMethod end

include("./bases/base_preprocess.jl")
include("./bases/base_linearequations.jl")
include("./bases/base_NLequations.jl")

include("linearSolver.jl")
include("nonlinearSolver.jl")

export
    # Problem types
    LinearProblem, NonlinearProblem,
    # Type function
    NLProblem,
    # Solve function
    solve,
    # Linear methods
    Direct,
    LUFactorized,
    CG, ConjugateGradient,
    # Nonlinear methods
    Newton

export @variables, Equation
    
end
