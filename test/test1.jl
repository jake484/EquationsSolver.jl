using EquationsSolver
using Test
using Symbolics
@testset "LinearProblem" begin
    @variables x, y
    eqs = [
        x - y ~ 0,
        x + y ~ 1
    ]
    vars = Dict(x => 1.0, y => 1.0)
    res, varVector = LinearProblem(eqs, vars)
    @test typeof(res) == LinearProblem
end
