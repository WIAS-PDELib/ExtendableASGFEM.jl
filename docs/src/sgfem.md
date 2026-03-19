# Overview

This page gives a concise overview of the main building blocks 
of the stochastic Galerkin (SG) methods implemented in this repository.

## Key blocks
- Represent uncertain coefficients/inputs by a finite parametric expansion (e.g. Karhunen–Loève / KL expansion).
- Expand the solution unknowns in a tensor product of a spatial finite element basis and a stochastic polynomial basis (Galerkin ansatz).
- Galerkin project onto the combined basis to obtain a coupled deterministic system for the expansion coefficients.
- Solve the system with solvers that exploit the tensor/block structure.
- Adaptivity in space and in the stochastic index set to reduce costs.

## Example: Parametric Poisson Problem

To illustrate the workflow, consider a Poisson problem with parametric diffusion coefficient:

```math
-\mathrm{div}(a(y,x) \nabla u(y,x)) = f(x) \quad \text{for } (y,x) \in \Gamma \times D
```

where the coefficient has the Karhunen-Loève expansion:

```math
a(y,x) = a_0(x) + \sum_{m=1}^M \sqrt{\lambda_m} \phi_m(x) y_m
```

with mean field $a_0$, eigenpairs $(\lambda_m,\phi_m)$ and parameters $y_m \in [-1,1]$.

### Stochastic Galerkin Formulation

1. Expand the solution in tensorized basis functions:
   ```math
   u(y,x) = \sum_{\mu} u_\mu(x) H_\mu(y)
   ```
   where $H_\mu(y) = \prod_{m=1}^M H_{\mu_m}(y_m)$ are multivariate Legendre polynomials.

2. Galerkin projection yields the weak form:
   ```math
   \sum_{\mu} \int_\Gamma a(y,x) \nabla u_\mu(x) \cdot \nabla v(x) H_\mu(y) H_\nu(y)\,dy = \int_\Gamma f(x)v(x) H_\nu(y)\,dy
   ```
   for all test functions $v(x)$ and indices $\nu$.

3. This results in a coupled block system $\mathbf{A}\mathbf{u} = \mathbf{b}$ where:
   - Each block $\mathbf{A}_{\mu,\nu}$ involves the mean and KL terms of $a(y,x)$
   - The tensor structure allows efficient matrix-free operations


## Where to find documentation/implementations of the key blocks
1. Parametric model / KL representation of the random coefficient.
   - See:
     - [Stochastic coefficients](coefficients.md) — available random coefficients and polynomial chaos expansions.
     - [Model problems](modelproblems.md) — available model problems that involve random coefficients. See also source: `src/modelproblems/`.
2. Choose stochastic basis, orthogonal polynomials, suitable for the parameter space
   - See:
     - [Orthogonal polynomials and recurrence relations](orthogonal_polynomials.md) — Legendre / Hermite polynomials and recurrence coefficients.
     - [ONBasis (one-dimensional orthonormal basis)](onbasis.md) — construction and utilities for evaluating 1D orthonormal polynomial bases (norms, quadrature, evaluations).
     - [Tensorized / multivariate basis (TensorizedBasis)](tonbasis.md) — assembling multivariate bases from ONBasis instances and precomputing triple products.
1. Build spatial FE spaces and blocks of system matrix.
   - Use FESpace types (H1, HDIV, ...) from the ExtendableFEM ecosystem.
   - See: [ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl) — basis finite-element spaces, basis evaluations and low-level FE utilities like standard interpolations.
   - See: [ExtendableFEM.jl](https://github.com/WIAS-PDELib/ExtendableFEM.jl) — high-level deterministic operator assembly and helpers used throughout the codebase.
2. Assemble and solve the full SG system
   - Galerkin projection yields a coupled deterministic block system. Briefly:
     - Expand u(y,x)=∑_μ u_μ(x) H_μ(y), test with H_ν ⇒ block matrix with entries a_{μ,ν} A(·).
   - Solver choices:
     - direct assembly + dense solver (for debugging / small problems)
     - matrix‑free iterative solvers exploiting tensor/block structure + preconditioners
   - See problem solver implementations: `src/modelproblems/solvers_*.jl`.
3. Postprocessing, error estimation & adaptivity
   - Error estimators and marking criteria implemented in the script driver `scripts/poisson.jl` for the available Poisson model problems.
   - Error estimators are problem-dependent and can be currently found in `src/estimate.jl`
   - Spatial (mesh refinement) uses refinement routines from [ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl)
   - Stochastic refinement (enrich multi‑index set) uses functions from `src/mopcontrol.jl`
   - see documentation page on [Estimators](estimators.md) for some more details
   - Results, parameters and reproducible outputs are stored with DrWatson (see scripts/poisson.jl for naming pattern).
