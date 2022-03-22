
get_vars = Symbolics.get_variables

function isin(var::Num, vars)
    for varible in vars
        if isequal(varible, var)
            return true
        end
    end
    return false
end

function get_Num(vars::Dict)
    return collect(keys(vars))
end

function check_eqs(eqs)
    if typeof(eqs) == Equation
        return Vector{Equation}([eqs])
    elseif typeof(eqs) == Vector{Equation}
        return eqs
    else
        error("Error: type of eqs must be Equation or Vector{Equation}")
    end
end

function toexpr(eqs)
    expr = Vector{Num}([])
    for eq in eqs
        push!(expr,eq.lhs - eq.rhs)
    end
    return expr
end

function get_all_vars(eqs_expr)
    return get_vars(sum(eqs_expr))
end
function check_vars(eqs, vars)
    eqs_expr = toexpr(eqs)
    eqs_var = get_all_vars(eqs_expr)
    vars = collect(keys(vars))
    print(length(eqs_var))
    if length(eqs_var) != length(vars)
        error("Error: number of varialbes in vars don't equal to 
        number of variables in eqs")
    end
    for var in eqs_var
        try
            pop!(vars, var)
        catch
            error("Error: $(var) not in vars")
        end
    end
    rest_var = collect(keys(vars))
    if length(rest_var) > 0
        error("$(rest_var) not in eqs")
    end
    return true
end

function islinear(eqs, vars)
    for eq in eqs
        if !Symbolics.islinear(eq.lhs, get_Num(vars))
            return false
        end
    end
    return true
end