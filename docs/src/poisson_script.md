# Script on Poisson Problems

The script `scripts/poisson.jl` provides an implementation of adaptive stochastic Galerkin FEM for various formulations of the parametric Poisson problem.

## Available Problem Types

- `PoissonProblemPrimal`: Standard Poisson problem with linear coefficient `a`
- `LogTransformedPoissonProblemPrimal`: Log-transformed Poisson problem with exponential coefficient `exp(a)`
- `LogTransformedPoissonProblemDual`: Dual formulation of the log-transformed problem (WIP)

## Usage

The main function to run experiments is:
```julia
run(; problem = problem_type, kwargs...)
```

Results can be loaded and visualized using:
```julia
show_results(; kwargs...)    # show results and statistics
produce_plots(; kwargs...)   # generate convergence and other plots
```

## Parameters

The following parameters can be specified (with their defaults):

### Problem Configuration
- `problem = LogTransformedPoissonProblemPrimal`: Problem type to solve
- `domain = "square"`: Domain shape ("square" or "lshape")
- `order = 1`: Polynomial order of finite elements
- `initial_refs = 1`: Initial uniform mesh refinements

### Stochastic Parameters
- `C = StochasticCoefficientCosinus`: Type of stochastic coefficient
- `decay = 2`: Decay rate of the coefficient expansion
- `mean = 0`: Mean value of the coefficient
- `maxm = 150`: Maximum number of terms in coefficient expansion
- `initial_modes = [[0]]`: Initial active stochastic modes

### Adaptive Parameters
- `θ_stochastic = 0.5`: Marking parameter for stochastic refinement
- `θ_spatial = 0.5`: Marking parameter for spatial refinement
- `factor_tail = 1`: Factor for tail estimator comparison
- `tail_extension = [10, 2]`: Extension parameters for [0] mode and others
- `maxdofs = 1.0e4`: Maximum degrees of freedom
- `nsamples = 150`: Number of samples for error computation
- `use_equilibration_estimator = false`: Use equilibration error estimator instead of standard residual estimator

### Solver Options
- `use_iterative_solver = true`: Use iterative solver (instead of a direct solver with full matrix)
- `bonus_quadorder_a = 2`: Additional quadrature order for coefficient
- `bonus_quadorder_f = 0`: Additional quadrature order for right-hand side

### Right-Hand Side
- `f = (result, qpinfo) -> (result[1] = 1)`: Right-hand side function

### Visualization
- `Plotter = nothing`: Plotting backend (e.g., CairoMakie)

## Results Storage

Results are automatically stored using [DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl) in the `data` folder. The filename includes most problem parameters for reproducibility:

```
data/[problem]/[domain]/order=[order]_maxdofs=[maxdofs]_decay=[decay]_mean=[mean]_θ=([θ_spatial],[θ_stochastic])_tail=[tail_extension]
```

If the equilibration estimator is used, `_eq` is appended to the filename.