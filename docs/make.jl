using EquationsSolver
using Documenter

DocMeta.setdocmeta!(EquationsSolver, :DocTestSetup, :(using EquationsSolver); recursive=true)

makedocs(;
    modules=[EquationsSolver],
    authors="yjy <522432938@qq.com> and contributors",
    repo="https://github.com/jake484/EquationsSolver.jl/blob/{commit}{path}#{line}",
    sitename="EquationsSolver.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jake484.github.io/EquationsSolver.jl",
        assets=["assets/logo.png"],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorials/LinearProblem.md",
            "tutorials/NLProblem.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/jake484/EquationsSolver.jl",
    devbranch="main",
)
