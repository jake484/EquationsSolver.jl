using EquationsSolver
using Test
using Symbolics

@testset "NonlinearProblem generate" begin
    @variables x,y,z

    eqs=[
        x^2+y^2+z^2~29,
        x+2*y+3*z~20,
        x*y-y*z+2*z*x~10
    ]
    #x=2,y=3,z=4
    val=Dict(
        x=>5,
        y=>2,
        z=>3
    )

    nlp=NonlinearProblem(eqs,val)
    @test typeof(nlp)==EquationsSolver.NonlinearProblem_functiontype

    nlp.vars

    thisNP=nlp.NP

    @test typeof(thisNP)==NonlinearProblem

    ja=thisNP.jacobian

    @test typeof(ja) <: Function 

    res=solve(thisNP,Newton())

    x0=Dict(x=>2,y=>3,z=>4)

    x1=Dict(zip(nlp.vars,res))

    @test x1[x]≈x0[x] &&
          x1[y]≈x0[y] &&
          x1[z]≈x0[z]

    @test solve(nlp)==x1
end

@testset "function input" begin
    function myfun(x)
        return [x[1]-cos(x[2]),x[2] - 2*cos(x[1])]
    end
    
    guessValue=[0.0,0.0]
    
    NLP=NonlinearProblem(myfun,guessValue)
    
    x=solve(NLP,Newton(); maxiter=100, abstol=1e-8)
    @test maximum(abs.(myfun(x))) < 1e-8
end

@testset "NLJacobian" begin
    function myfun(x)
        return [cos(x[2])-x[1],2*cos(x[1])-x[2]]
    end
    
    guessValue=[0.0,0.0]
    
    NLP=NonlinearProblem(myfun,guessValue)
    
    x=solve(NLP,NLJacobian();maxiter=1000,abstol=1e-8)
    @test maximum(abs.(myfun(x))) < 1e-8
end
