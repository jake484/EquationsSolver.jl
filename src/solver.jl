solve_liner = Symbolics.solve_for

function solve(problem::LinearProblem)
    eqs = problem.eqs
    vars = problem.vars
    sol = solve_liner(eqs,vars)
    return Dict([(vars[i],sol[i]) for i in 1:length(vars) ])
end

function solve(problem::NonlinearProblem)
    eqs = problem.eqs
    vars = problem.vars
    u0 = collect(values(vars))
    sol = solve_liner(eqs,vars)
    return Dict([(vars[i],sol[i]) for i in 1:length(vars) ])
end
