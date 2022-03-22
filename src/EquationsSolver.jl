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
    eqs::Vector{Equation}
    vars::Vector{Num}
    maxiters::Int
end
solve_liner = Symbolics.solve_for
function EquationsProblem(eqs::Any, vars, maxiters=10000)
    eqs = check_eqs(eqs)
    res = check_vars(eqs, vars)
    if islinear(eqs, vars)
        return LinearProblem(eqs, vars, maxiters)
    else
        return NonlinearProblem(eqs, vars, maxiters)
    end
end

export EquationsProblem

end
