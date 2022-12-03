solve_liner = Symbolics.solve_for

#采用Julia自带的求解方法
function solve(problem::LinearProblem)
    problem.A \ problem.b
end
solve(problem::LinearProblem, ::Direct) = solve(problem::LinearProblem)

#对于已经进行LU分解的问题采用LU分解的回带过程
function solve(problem::LinearProblem, ::LUFactorized, LUFactor::LU)
    LU_solve(LUFactor, problem.b)
end

#采用共轭梯度法，这里直接将矩阵A改为对称矩阵。如果A本来就对称，则没有影响
function solve(problem::LinearProblem, ::CG; abserr=1e-8)
    CG_solve(Symmetric(problem.A), problem.b; abserr=abserr)
end
@deprecate solve(problem::LinearProblem, ::ConjugateGradient; abserr=1e-8) solve(problem, CG(); abserr=abserr)