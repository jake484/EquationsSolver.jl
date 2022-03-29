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
        x - y - 0,
        x + y - 1
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


@testset "NLProblem_1: x , y" begin
    @variables x, y
    eqs = [
        cos(x^2 + 0.4 * y) + x^2 + y^2 ~ 1.6
        1.5 * x^2 - 1 / 0.36 * y^2 ~ 1.0
    ]
    vars = Dict(x => 2.0, y => 1.0)
    pro = NLProblem(eqs, vars)
    res = solve(pro)
    @test round(res[x], digits=2) == 1.04 &&
          round(res[y], digits=2) == 0.47
end

@testset "NLProblem_2: x , y" begin
    @variables x, y
    eqs = [
        4 * x + 0.1 * exp(x) ~ y + 1
        4 * y + x^2 / 8 ~ x
    ]
    vars = Dict(x => 2.0, y => 1.0)
    pro = NLProblem(eqs, vars)
    res = solve(pro)
    @test round(res[x], digits=2) == 0.23 &&
          round(res[y], digits=2) == 0.06
end

@testset "NLProblem: x , y , z" begin
    @variables x, y, z
    eqs = [
        x^2 + y^2 + z^2 ~ 1.0
        2 * x^2 + y^2 - 4 * z ~ 0
        3 * x^2 - 4 * y^2 + z^2 ~ 0
    ]
    vars = Dict(x => 2.0, y => 1.0, z => 1.0)
    pro = NLProblem(eqs, vars)
    res = solve(pro)
    @test round(res[x], digits=2) == 0.70 &&
          round(res[y], digits=2) == 0.63 &&
          round(res[z], digits=2) == 0.34
end

@testset "NLProblem: x , y , z" begin
    @variables x, y, z
    eqs = [
        x^2 + y^2 + z^2 ~ 1.0
        2 * x^2 + y^2 - 4 * z ~ 0
        3 * x^2 - 4 * y^2 + z^2 ~ 0
    ]
    vars = Dict(x => 2.0, y => 1.0, z => 1.0)
    pro = NLProblem(eqs, vars)
    res = solve(pro)
    @test round(res[x], digits=2) == 0.70 &&
          round(res[y], digits=2) == 0.63 &&
          round(res[z], digits=2) == 0.34
end

@testset "solve-LinearProblem" begin
    N = 20
    @variables x[N]
    eqs = Vector{Equation}([])
    for i in 1:N
        rand_num = rand(1:20, Vector{Int}, 3)
        s = sum(rand_num)
        if i == 1
            push!(eqs, sum([rand_num[j-i+1] * x[j] for j in i:i+2]) ~ s)
        elseif i == N
            push!(eqs, sum([rand_num[j-i+3] * x[j] for j in i-2:i]) ~ s)
        else
            push!(eqs, sum([rand_num[j-i+2] * x[j] for j in i-1:i+1]) ~ s)
        end
    end
    vars=Dict([x[i]=>0.0 for i in 1:N])
    pro = LinearProblem(eqs,vars)
    res = solve(pro)
    @test round.(collect(values(res)), digits=1) == [1.0 for i in 1:N]
end



@testset "solve-LinearProblem by NLProblem" begin
    N = 100
    @variables x[N]
    eqs = Vector{Equation}([])
    for i in 1:N
        rand_num = rand(1:20, Vector{Int}, 3)
        s = sum(rand_num)
        if i == 1
            push!(eqs, sum([rand_num[j-i+1] * x[j] for j in i:i+2]) ~ s)
        elseif i == N
            push!(eqs, sum([rand_num[j-i+3] * x[j] for j in i-2:i]) ~ s)
        else
            push!(eqs, sum([rand_num[j-i+2] * x[j] for j in i-1:i+1]) ~ s)
        end
    end
    vars=Dict([x[i]=>0.0 for i in 1:N])
    pro = NLProblem(eqs,vars)
    res = solve(pro)
    @test round.(collect(values(res)), digits=1) == [1.0 for i in 1:N]
end

