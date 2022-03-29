# Solve LinearProblem

There some LinearProblem example.

## Simple linear problem

solve function:

$$
\left\{\begin{matrix} x+y+z = 3\\ 
x + 4y + 9z = 14 \\ 
x + 2y + 3z = 6
\end{matrix}\right.    
$$

julia code:

```julia
@variables x,y,z
eqs = [
    x + y + z ~ 3,
    x + 4y + 9z ~ 14,
    x + 2y + 3z ~ 6
]
vars = Dict(x => 1.0, y => 0.0, z => 0.0)
pro = LinearProblem(eqs,vars)
res = solve(pro)

res = Dict{Num, Float64}(
    y => 0.9999999999999996, 
    z => 1.0000000000000002, 
    x => 1.0000000000000002
) 
```

## Using NLProblem to solve LP 

Actually, the way EquationsSolver using to solve NLP is Newton iterative method. It also can solve LP, because LP can be treated as a special NLP.

We can make a problem like this. The value of all symbols is 1.0.

```julia
N = 100
@variables x[N]
eqs = Vector{Equation}([])
for i in 1:N
    rand_num = rand(1:20, Vector{Int}, 3)
    s = sum(rand_num)
    if i == 1
        push!(eqs, sum([rand_num[j-i+1] * x[j] for j in i:i+2]) ~ s)
    elseif i == N
        push!(eqs, sum([rand_num[j-i+3] * x[j] for j in i-2:i]) ~ s)
    else
        push!(eqs, sum([rand_num[j-i+2] * x[j] for j in i-1:i+1]) ~ s)
    end
end
vars=Dict([x[i]=>0.0 for i in 1:N])
pro = NLProblem(eqs,vars)
res = solve(pro)
```