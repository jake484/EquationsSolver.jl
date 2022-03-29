# EquationsSolver

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jake484.github.io/EquationsSolver.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jake484.github.io/EquationsSolver.jl/dev)
[![Build Status](https://travis-ci.com/jake484/EquationsSolver.jl.svg?branch=main)](https://travis-ci.com/jake484/EquationsSolver.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jake484/EquationsSolver.jl?svg=true)](https://ci.appveyor.com/project/jake484/EquationsSolver-jl)
[![Build Status](https://api.cirrus-ci.com/github/jake484/EquationsSolver.jl.svg)](https://cirrus-ci.com/github/jake484/EquationsSolver.jl)
[![Coverage](https://codecov.io/gh/jake484/EquationsSolver.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jake484/EquationsSolver.jl)

EquationsSolver is a little user-friendly tool to solve linear equations and nonlinear equations.

It is based on Symbolics.jl

It can test your little problems very fast and easily.

For example,

```julia
@variables x
eqs = [
    x + 5 ~ exp(x)
]
vars = Dict(x => 2.0)
pro = NLProblem(eqs,vars)
res = solve(pro)
```

