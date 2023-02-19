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
    @test solve(lp, LUFactorized();LUFactor = F) ≈ x0
    A = A + A' + 2 * I
    b = 2 * b .+ 1
    lp = LinearProblem(A, b)
    x0 = A \ b
    @test solve(lp, CG()) ≈ x0
    @test solve(lp, ConjugateGradient(); abstol=1e-8) ≈ x0
end

@testset "GMRES method" begin
    n=1000
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

    println("\n GMRES time:")
    @time x=solve(lp,GMRESM();m=m,maxiter=maxiters)
    println("\n")
    println("\n Direct time:")
    @time x0=A\b
    println("\n")
    @test x ≈ x0
    print(maximum(abs.(x-x0)))
end

@testset "Set before solve" begin
    A = rand(4, 4)
    b = rand(4)
    lp = LinearProblem(A, b)
    x0 = A \ b
    
    lpM1=Direct()
    setting1::Dict{Symbol,Real}=Dict()

    lpM2=GMRESM()
    setting2=Dict(:m=>4)

    F=lu(A)
    lpM3=LUFactorized()
    setting3=Dict(:LUFactor=>F)

    x1=solve(lp,lpM1;setting1...)
    @test isapprox(x0,x1)

    x2=solve(lp,lpM2;setting2...)
    @test isapprox(x0,x2)

    x3=solve(lp,lpM3;setting3...)
    @test isapprox(x0,x3)

    setting4=Dict(
        :method=>GMRESM(),
        :m=>4
    )
    x4=solve(lp;lpSolverSetting=setting4)
    @test isapprox(x0,x4)

end
