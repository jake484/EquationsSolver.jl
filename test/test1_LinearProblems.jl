using EquationsSolver
using Test
using Symbolics
using LinearAlgebra, SparseArrays

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

@testset "GMRES method" begin
    n=100000
    k=100.0
    A=sprand(n,n,k/n)
    dig=A*ones(n)
    ia,ja,va=findnz(A)
    for i=1:n
        push!(ia,i)
        push!(ja,i)
        push!(va,dig[i])
    end
    A=sparse(ia,ja,va)
    b=rand(n)
    guessValue=zeros(n)
    m=10
    maxiters=100
    abstol=1e-6

    lp=LinearProblem(A,b,guessValue)

    @time x,err,iter=solve(lp,GMRESM(),m;maxiter=maxiters)
    @test abs(err) < abstol
end
