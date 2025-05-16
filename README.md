# ExtendableASGFEM

This package implements the stochastic Galerkin finite element method for certain model problems,
including iterative solvers and a posteriori error control in two dimensions.

## Installation

Instantiate this project after cloning or updating via
```
$ julia --project=.
julia> using Pkg
julia> Pkg.instantiate()
```

## Scripts/Examples

- `scripts/poisson.jl`: Adaptive stochastic Galerkin FEM for the Poisson problem with stochastic diffusion coefficient
