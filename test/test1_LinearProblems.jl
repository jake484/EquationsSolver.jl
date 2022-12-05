using EquationsSolver
using Test
using Symbolics
using LinearAlgebra

@testset "LinearProblem" begin
    @variables x, y
    eqs = [
        x - y ~ 0,
        x + y ~ 1
    ]
    vars = Dict(x => 1.0, y => 1.0)
    res = LinearProblem(eqs, vars)
    @test typeof(res) == EquationsSolver.LinearProblem_functiontype
end


@testset "linear solver" begin
    A = rand(4, 4)
    b = rand(4)
    lp = LinearProblem(A, b)
    x0 = A \ b
    @test solve(lp) ≈ x0
    @test solve(lp, Direct()) ≈ x0
    F = lu(A)
    @test solve(lp, LUFactorized(), F) ≈ x0
    A = A + A' + 2 * I
    b = 2 * b .+ 1
    lp = LinearProblem(A, b)
    x0 = A \ b
    @test solve(lp, CG()) ≈ x0
    @test solve(lp, ConjugateGradient(); abstol=1e-8) ≈ x0
end