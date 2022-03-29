# EquationsSolver

Documentation for [EquationsSolver](https://github.com/jake484/EquationsSolver.jl).

## Installation

To install EquationsSolver.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("EquationsSolver")
```

## Solve Problem

There two steps to solve a problem-Define and Solve.

### Define problem

Defining like this.

```julia
using EquationsSolver
@variables x
eqs = [
    x + 5 ~ exp(x)
]
vars = Dict(x => 2.0)
pro = NLProblem(eqs,vars)
```

Using @variables to define a symbol variable. It is from Symbolics.jl

```julia
@variables x
```

And then write equations with x. Only one equation is ok. It's better to write equation in a vector.

```julia
eqs = [
    x + 5 ~ exp(x)
]
```

Next, give symbol x a initial value by Dict.

```julia
vars = Dict(x => 2.0)
```

Finally, define problem-**LinearProblem** or **NonlinearProblem**.

```julia
pro = NLProblem(eqs,vars)
```

### Solve problem

Sovling problem is very easy.Just use solve function and get the result.

```julia
res = solve(pro)
```

If print the res, we will get

```julia
Dict{Num, Float64}(x => 1.9368474072202186)
```
