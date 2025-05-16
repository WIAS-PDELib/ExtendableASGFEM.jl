using Documenter
using ExtendableASGFEM

push!(LOAD_PATH,"../src/")

makedocs(
		modules = [ExtendableASGFEM],
		sitename = "ExtendableASGFEM.jl",
		authors = "Christian Merdon, Martin Eigel",
		repo = "https://lab.wias-berlin.de/merdon/ASGFEMJulia",
		# repolink = "https://lab.wias-berlin.de/merdon/ASGFEMJulia", 
		format = Documenter.HTML(; mathengine = MathJax3()),
		clean = false,
		source  = "src",
		build   = "build",
		checkdocs = :all,
		warnonly = true,
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