using EquationsSolver
using Test

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

