solve_liner = Symbolics.solve_for


function solve(problem::LinearProblem)
    eqs = problem.eqs
    vars = problem.vars
    sol = solve_liner(eqs, vars)
    return Dict([(vars[i], sol[i]) for i in 1:length(vars)])
end


function one_cal(array, dict, f, ja)
    return array .- inv(Symbolics.value.(substitute.(ja, (dict,)))) * Symbolics.value.(substitute.(f, (dict,)))
end


function get_dict(array, dict, vars)
    for i in 1:length(vars)
        dict[vars[i]] = array[i]
    end
    return dict
end


function _solve(dict, f, ja, maxiter, abstol, vars)
    last_v = Inf
    new_v = collect(values(dict))
    dict = dict
    count = 0

    while count < maxiter && maximum(abs.(new_v .- last_v)) > abstol
        last_v = new_v
        new_v = one_cal(new_v, dict, f, ja)
        dict = get_dict(new_v, dict, vars)
        count += 1
    end

    if count >= maxiter
        error("Error: iterate times exceed maxiter")
    end

    return dict
end



function solve(problem::NonlinearProblem)
    eqs = problem.eqs
    dict_u = problem.vars
    vars = collect(keys(problem.vars))
    maxiters = problem.maxiters
    abstol = problem.abstol
    jac = Symbolics.jacobian(eqs, vars)
    sol = _solve(dict_u, eqs, jac, maxiters, abstol,vars)
    return sol
end
