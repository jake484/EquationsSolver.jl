using EquationsSolver
using Test
using Symbolics
#=
function LinearProblem(eqs::Any, vars::Dict, maxiters=10000)
    checked_eqs = check_eqs(eqs)
    res = check_vars(checked_eqs, vars)
    dict=Dict(key=>0.0 for key in keys(vars))
    A=Symbolics.jacobian(checked_eqs, vars)
    b=-Symbolics.value.(substitute.(checked_eqs, (dict,)))
    return LinearProblem(A, b, maxiters)
end
=#

A = rand(4, 4)
b = rand(4)
lp = LinearProblem(A, b)

@variables x, y
eqs = [
    x - y ~ 0,
    x + y ~ 1
]
vars = Dict(x => 1.0, y => 1.0)

checked_eqs = EquationsSolver.check_eqs(eqs)
res = EquationsSolver.check_vars(checked_eqs, vars)
dict = Dict(key => 0.0 for key in keys(vars))
A = Symbolics.jacobian(checked_eqs, collect(keys(vars)))
b = -Symbolics.value.(substitute.(checked_eqs, (dict,)))
LinearProblem(A, b, maxiters)