using EquationsSolver
using Test

function poly(x,a)
    s=a[1]
    n=length(a)
    for i=2:n
        s=s*x+a[i]
    end
    return s
end

@testset "real root functions" begin
    a=[1,-6,13,-12,4]#(x-1)(x-2)^2
    prob=PolyProblem(a)
    realRoots,complexRoots=solve(prob)
    @test maximum(abs.(map(x->poly(x,a),realRoots)))<1e-3
    @test isempty(complexRoots)
end

@testset "complex root functions" begin
    a=[1,0,3,0,2]#(x^+1)(x^2+2)
    prob=PolyProblem(a)
    realRoots,complexRoots=solve(prob)
    @test maximum(abs.(map(x->poly(x,a),complexRoots)))<1e-3
    @test isempty(realRoots)
end

@testset "hybrid root functions" begin
    a=[1,0,-2,-1,-3,2,0,3]#(x^2+1)(x^2-3)(x^3-1)
    prob=PolyProblem(a)
    realRoots,complexRoots=solve(prob,Bairstov();abstol=1e-20)
    @test maximum(abs.(map(x->poly(x,a),complexRoots)))<1e-3
    @test maximum(abs.(map(x->poly(x,a),realRoots)))<1e-3
end
