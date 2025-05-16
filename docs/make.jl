using Documenter
using ExtendableASGFEM

makedocs(
		modules = [ExtendableASGFEM],
		sitename = "ExtendableASGFEM.jl",
		authors = "Christian Merdon, Martin Eigel",
		format = Documenter.HTML(; repolink = "https://github.com/WIAS-PDELib/ExtendableASGFEM.jl", mathengine = MathJax3()),
        clean = false,
		checkdocs = :all,
		warnonly = false,
		doctest = true,
		pages = [
			"Home" => "index.md"
			"KLE expansions" => [
				"coefficients.md",
			]
			"Stochastic discretization" => [
				"orthogonal_polynomials.md",
				"onbasis.md",
				"tonbasis.md",
			]
			"Solvers" => [
				"sgfevector.md",
				"modelproblems.md",
			]
			"Adaptivity" => [
				"estimators.md"
			]
			"Plotting" => [
				"plots.md"
			]
		],
	)


deploydocs(
    repo = "github.com/WIAS-PDELib/ExtendableASGFEM.jl",
)