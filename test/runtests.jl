using EquationsSolver
using Test
using Symbolics


@testset "check_eqs" begin
    @variables x
    @test typeof(EquationsSolver.check_eqs(x + 1 ~ 1)) == typeof([x + 1 ~ 1])
    @test typeof(EquationsSolver.check_eqs([x + 1 ~ 1])) == typeof([x + 1 ~ 1])
end
@variables x, y
@testset "get_vars" begin
    res = EquationsSolver.get_vars(x + y)
    @test EquationsSolver.isin(x, res) == true
    @test EquationsSolver.isin(x, res) == true
end

@testset "islinear_sym" begin
    @test Symbolics.islinear(x - y, [x, y])
    @test Symbolics.islinear(x, [x, y])
    @test Symbolics.islinear(y, [x, y])
end
dic = Dict(x => 1.0, y => 1.0)

@test EquationsSolver.isin(x, EquationsSolver.get_Num(dic))

@testset "get_vars" begin
    res = EquationsSolver.get_vars(x + y)
    @test length(res) == length(EquationsSolver.get_Num(dic))
end

@testset "islinear" begin
    @variables x, y
    @test EquationsSolver.islinear([x - y ~ 0, x + y ~ 1], Dict(x => 1.0, y => 1.0))
end
expr = EquationsSolver.toexpr([x - y ~ 0, x + y ~ 1])
EquationsSolver.get_all_vars(expr)
prob = EquationsProblem([x - y ~ 0, x + y ~ 0], Dict(x => 1.0, y => 1.0))