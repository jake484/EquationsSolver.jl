using EquationsSolver
using Test
using Symbolics, ForwardDiff

@testset "NonlinearProblem generate" begin
    @variables x, y, z

    eqs = [
        x^2 + y^2 + z^2 ~ 29,
        x + 2 * y + 3 * z ~ 20,
        x * y - y * z + 2 * z * x ~ 10
    ]
    #x=2,y=3,z=4
    val = Dict(
        x => 5,
        y => 2,
        z => 3
    )

    nlp = NonlinearProblem(eqs, val)
    @test typeof(nlp) == EquationsSolver.NonlinearProblem_functiontype

    nlp.vars

    thisNP = nlp.NP

    @test typeof(thisNP) == NonlinearProblem

    ja = thisNP.jacobian

    @test typeof(ja) <: Function

    res = solve(thisNP, Newton())

    x0 = Dict(x => 2, y => 3, z => 4)

    x1 = Dict(zip(nlp.vars, res))

    @test x1[x] ≈ x0[x] &&
          x1[y] ≈ x0[y] &&
          x1[z] ≈ x0[z]

    @test solve(nlp) == x1
end

@testset "function input" begin
    function myfun(x)
        return [x[1] - cos(x[2]), x[2] - 2 * cos(x[1])]
    end

    guessValue = [0.0, 0.0]

    NLP = NonlinearProblem(myfun, guessValue)

    x = solve(NLP, Newton(); maxiter=100, abstol=1e-8)
    @test maximum(abs.(myfun(x))) < 1e-8
end

@testset "NLJacobian" begin
    function myfun(x)
        return [cos(x[2]) - x[1], 2 * cos(x[1]) - x[2]]
    end

    guessValue = [0.0, 0.0]

    NLP = NonlinearProblem(myfun, guessValue)

    x = solve(NLP, NLJacobian(); maxiter=1000, abstol=1e-8)
    @test maximum(abs.(myfun(x))) < 1e-8
end

@testset "SimplifiedNewton" begin
    function myfun(x)
        return [cos(x[2]) - x[1], 2 * cos(x[1]) - x[2]]
    end

    guessValue = [0.0, 0.0]

    NLP = NonlinearProblem(myfun, guessValue)

    x = solve(NLP, SimplifiedNewton(); maxiter=1000, abstol=1e-8)
    @test maximum(abs.(myfun(x))) < 1e-8


    function myfun2(x)
        [x[1] - cos(x[1])]
    end
    guessValue = [0.0]
    NLP2 = NonlinearProblem(myfun2, guessValue)
    x = solve(NLP2, SimplifiedNewton(); maxiter=1000, abstol=1e-8)
    @test maximum(abs.(myfun2(x))) < 1e-8
end

@testset "Secant" begin
    function myfun1(x)
        [
            x[1]^2 + x[2]^2 + x[3]^2 - 29,
            x[1] + 2 * x[2] + 3 * x[3] - 20,
            x[1] * x[2] - x[2] * x[3] + 2 * x[3] * x[1] - 10
        ]
    end

    guessValue = [1.0, 0.2, 3.4]

    NLP = NonlinearProblem(myfun1, guessValue)

    lpM = GMRESM()
    lpSetting = Dict(
        :m => 3
    )

    nlSetting = Dict(
        :h => 1.0,
        :lpMethod => lpM,
        :lpSolverSetting => lpSetting
    )

    x = solve(NLP, Secant(); nlSetting...)
    @test isapprox(myfun1(x), zeros(3); atol=1e-5, rtol=0.1)
end

@testset "BroydenRank1" begin
    function myfun1(x)
        [
            x[1]^2 + x[2]^2 + x[3]^2 - 29,
            x[1] + 2 * x[2] + 3 * x[3] - 20,
            x[1] * x[2] - x[2] * x[3] + 2 * x[3] * x[1] - 10
        ]
    end

    #guessValue = [1.0, 2.0, 3.0] 此处雅克比矩阵奇异
    guessValue = [1.5, 2.0, 3.0]

    NLP = NonlinearProblem(myfun1, guessValue)

    x = solve(NLP, BroydenRank1())
    @test isapprox(myfun1(x), zeros(3); atol=1e-5, rtol=0.1)
end
