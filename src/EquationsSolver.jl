module EquationsSolver

using Symbolics, ForwardDiff
using LinearAlgebra, SparseArrays
import Symbolics: @variables, Equation

abstract type AbstractSolveMethod end
abstract type AbstractEquationsProblem end

abstract type AbstractLinearMethod <: AbstractSolveMethod end
abstract type AbstractNLMethod <: AbstractSolveMethod end
abstract type AbstractPolyMethod <: AbstractSolveMethod end

include("./bases/base_preprocess.jl")
include("./bases/base_linearequations.jl")
include("./bases/base_NLequations.jl")
include("./bases/base_polynomial.jl")

include("linearSolver.jl")
include("nonlinearSolver.jl")
include("polynomialSolver.jl")

export
    # Problem types
    LinearProblem, NonlinearProblem, PolyProblem,
    # Solve function
    solve,
    # Linear methods
    Direct,
    LUFactorized,
    CG, ConjugateGradient,
    GMRESM,
    # Nonlinear methods
    Newton, NLJacobian, SimplifiedNewton, Secant, BroydenRank1,
    # Polynomial methods
    Bairstov

export @variables, Equation

end
