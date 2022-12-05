using EquationsSolver
using Test
using Symbolics

@testset "NLProblem generate" begin
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

    nlp=NLProblem(eqs,val)
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