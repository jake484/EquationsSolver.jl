# Solve NonlinearProblem

There some NonlinearProblem example.

## NLP examples

solve function:

$$
\left\{\begin{matrix} x^2 + y^2 + z^2 = 1.0\\ 
2  x^2 + y^2 - 4  z = 0 \\ 
 3  x^2 - 4  y^2 + z^2 = 0
\end{matrix}\right.
$$

julia code:

```julia
@variables x, y, z
eqs = [
    x^2 + y^2 + z^2 ~ 1.0
    2 * x^2 + y^2 - 4 * z ~ 0
    3 * x^2 - 4 * y^2 + z^2 ~ 0
]
vars = Dict(x => 2.0, y => 1.0, z => 1.0)
pro = NLProblem(eqs, vars)
res = solve(pro)
@show res

res = Dict{Num, Float64}(
    y => 0.6285242979602138, 
    z => 0.3425641896895694, 
    x => 0.6982886099715139
    )
```

solve function:

$$
\left\{\begin{matrix}
    cos(x^2 + 0.4  y) + x^2 + y^2 = 1.6\\
    1.5  x^2 - 1 / 0.36 * y^2 = 1.0
\end{matrix}\right.
$$

julia code:

```julia
@variables x, y
eqs = [
    cos(x^2 + 0.4 * y) + x^2 + y^2 ~ 1.6
    1.5 * x^2 - 1 / 0.36 * y^2 ~ 1.0
]
vars = Dict(x => 2.0, y => 1.0)
pro = NLProblem(eqs, vars)
res = solve(pro)

res = Dict{Num, Float64}(
    y => 0.47172595266055173, 
    x => 1.0386292376770578
    )
```

solve function:

$$
\left\{\begin{matrix}
    4  x + 0.1 * e^x = y + 1\\
    4 y + x^2 / 8 = x
\end{matrix}\right.
$$

julia code:

```julia
@variables x, y
eqs = [
    4 * x + 0.1 * exp(x) ~ y + 1
    4 * y + x^2 / 8 ~ x
]
vars = Dict(x => 2.0, y => 1.0)
pro = NLProblem(eqs, vars)
res = solve(pro)

res = Dict{Num, Float64}(
    y => 0.056451519652141575, 
    x => 0.23256700509067185
    ) 
```
