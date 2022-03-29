using EquationsSolver
using Test
# using Symbolics


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



@testset "one_cal" begin
    @variables x
    eqs = [
        x + 5 - exp(x)
    ]
    vars = Dict(x => 2.0)
    jac = Symbolics.jacobian(eqs, [x])
    res = EquationsSolver.one_cal([1.0], vars, eqs, jac)
    @test length(res) == 1
end



@testset "get_dict" begin
    @variables x, y, z
    value = [1.0, 2.0, 3.0]
    vars = Dict(x => 2.0, y => 1.0, z => 3.0)
    res = EquationsSolver.get_dict(value, vars, [x, y, z])
    @test res == Dict(x => 1.0, y => 2.0, z => 3.0)
end



@testset "_solve" begin
    @variables x, y, z
    eqs = [
        x + 5 - exp(x)
    ]
    vars = Dict(x => 2.0)
    jac = Symbolics.jacobian(eqs, [x])
    res = EquationsSolver._solve(vars, eqs, jac, 1000, 1.0E-6, [x])
    @test length(res) == 1

    eqs = [
        x + y + z - 3,
        x + 4y + 9z - 14,
        x + 2y + 3z - 6
    ]
    vars = Dict(x => 1.0, y => 0.0, z => 0.0)
    jac = Symbolics.jacobian(eqs, [x, y, z])
    res = EquationsSolver._solve(vars, eqs, jac, 10000, 1.0E-6, [x, y, z])
    @test round.(collect(values(res)), digits=1) == [1.0,1.0,1.0]
end


@testset "solve-NLProblem" begin
    @variables x
    eqs = [
        x + 5 ~ exp(x)
    ]
    vars = Dict(x => 2.0)
    pro = NLProblem(eqs,vars)
    res = solve(pro)
    @test length(res) == 1
end


@testset "solve-LinearProblem" begin
    @variables x,y,z
    eqs = [
        x + y + z ~ 3,
        x + 4y + 9z ~ 14,
        x + 2y + 3z ~ 6
    ]
    vars = Dict(x => 1.0, y => 0.0, z => 0.0)
    pro = LinearProblem(eqs,vars)
    res = solve(pro)
    @test round.(collect(values(res)), digits=1) == [1.0,1.0,1.0]
end

