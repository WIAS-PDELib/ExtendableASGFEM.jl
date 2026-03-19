[![Build status](https://github.com/WIAS-PDELib/ExtendableASGFEM.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/WIAS-PDELib/ExtendableASGFEM.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://wias-pdelib.github.io/ExtendableASGFEM.jl/stable/index.html)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wias-pdelib.github.io/ExtendableASGFEM.jl/dev/index.html)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

# ExtendableASGFEM

This package implements the stochastic Galerkin finite element method for certain two dimensional
model problems involving KLE of stochastic coefficients. The rather huge systems have a tensorized
structure and are solved by iterative solvers. A posteriori error estimators steer the spatial
and stochastic refinement.

The spatial discretization is based on the finite element package [ExtendableFEM.jl](https://github.com/WIAS-PDELib/ExtendableFEM.jl)/[ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl)

## Scripts/Examples

- `scripts/poisson.jl`: Adaptive stochastic Galerkin FEM for (different formulations of) the Poisson problem with stochastic diffusion coefficient
