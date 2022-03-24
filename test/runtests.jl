using EquationsSolver
using Test
using Symbolics


@testset "check_eqs" begin
    @variables x, y
    res = EquationsSolver.check_eqs(x + y ~ 1)
    @test isequal(res, [x + y - 1])
    res = EquationsSolver.check_eqs([x + y ~ 1])
    @test isequal(res, [x + y - 1])
end

@testset "get_vars" begin
    @variables x, y
    res = EquationsSolver.get_vars(x + y)
    @test EquationsSolver.isin(x, res) == true
    @test EquationsSolver.isin(x, res) == true
end

@testset "islinear_sym" begin
    @variables x, y
    @test Symbolics.islinear(x - y, [x, y])
    @test Symbolics.islinear(x, [x, y])
    @test Symbolics.islinear(y, [x, y])
end

@testset "isin" begin
    dic = Dict(x => 1.0, y => 1.0)
    @test EquationsSolver.isin(x, EquationsSolver.get_Num(dic))
end

@testset "get_vars" begin
    res = EquationsSolver.get_vars(x + y)
    @test length(res) == length(EquationsSolver.get_Num(dic))
end

@testset "islinear" begin
    @variables x, y
    @test EquationsSolver.islinear([x - y ~ 0, x + y ~ 1], Dict(x => 1.0, y => 1.0))
end

@testset "get_all_vars" begin
    @variables x, y, z
    expr1 = EquationsSolver.toexpr([x - y ~ 0, x + y ~ 1])
    expr2 = EquationsSolver.toexpr([x - y ~ 0, x + y - z ~ 1])
    res1 = EquationsSolver.get_all_vars(expr1)
    res2 = EquationsSolver.get_all_vars(expr2)
    @test isequal(res1, [x, y])
    @test isequal(res2, [x, y, z])
end

@testset "check_vars" begin
    @variables x, y
    eqs = [
        x - y ~ 0,
        x + y ~ 1
    ]
    vars = Dict(x => 1.0, y => 1.0)
    @test EquationsSolver.check_vars(eqs, vars)
end


@testset "LinearProblem" begin
    @variables x, y
    eqs = [
        x - y ~ 0,
        x + y ~ 1
    ]
    vars = Dict(x => 1.0, y => 1.0)
    res = LinearProblem(eqs, vars)
    @test typeof(res) == LinearProblem
end

@testset "solve" begin
    @variables x, y
    eqs = [
        x + y ~ 3.0,
        x + 2y ~ 5.0
    ]
    vars = Dict(x => 1.0, y => 1.0)
    res = EquationsSolver.solve(LinearProblem(eqs, vars))
    @test res == Dict(x => 1.0, y => 2.0)
end
