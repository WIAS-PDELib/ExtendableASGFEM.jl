using Documenter
using ExtendableASGFEM

makedocs(
    modules = [ExtendableASGFEM],
    sitename = "ExtendableASGFEM.jl",
    authors = "Christian Merdon, Martin Eigel",
    format = Documenter.HTML(; repolink = "https://github.com/WIAS-PDELib/ExtendableASGFEM.jl", mathengine = MathJax3()),
    clean = true,
    checkdocs = :all,
    warnonly = false,
    doctest = true,
    pages = [
        "Home" => "index.md",
        "SGFEM Overview" => "sgfem.md",
        "Model problems" => [
            "modelproblems.md",
            "coefficients.md",
        ],
        "Stochastic discretization" => [
            "orthogonal_polynomials.md",
            "onbasis.md",
            "tonbasis.md",
        ],
        "Solving" => [
            "solvers.md",
            "estimators.md",
        ],
        "Plotting" => [
            "plots.md",
        ],
        "Scripts" => [
            "poisson_script.md",
        ]
    ],
)


deploydocs(
    repo = "github.com/WIAS-PDELib/ExtendableASGFEM.jl",
)
