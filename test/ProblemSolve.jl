using EquationsSolver
using Test
using Symbolics

@testset "test: x + 5 - exp(x) ~ 0" begin
    @variables x, y, z
    eqs = [
        x + 5 - exp(x)
    ]
    vars = Dict(x => 2.0)
    jac = Symbolics.jacobian(eqs, [x])
    res = EquationsSolver._solve(vars, eqs, jac, 1000, 1.0E-6)
    @test length(res) == 1
end
