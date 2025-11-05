# Overview

This page gives a concise overview of the main building blocks 
of the stochastic Galerkin (SG) methods implemented in this repository.

## Key blocks
- Represent uncertain coefficients/inputs by a finite parametric expansion (e.g. Karhunen–Loève / KL expansion).
- Expand the solution unknowns in a tensor product of a spatial finite element basis and a stochastic polynomial basis (Galerkin ansatz).
- Galerkin project onto the combined basis to obtain a coupled deterministic system for the expansion coefficients.
- Solve the system with solvers that exploit the tensor/block structure.
- Adaptivity in space and in the stochastic index set to reduce costs.

## Where to find documentation/implementations of the key blocks
1. Parametric model / KL representation of the random coefficient.
   - See also:
     - [Stochastic coefficients](coefficients.md) — available random coefficients
     - and polynomial chaos expansions.
     - [Model problems](modelproblems.md) — available model problems that involve random coefficients. See also source: `src/modelproblems/`.

2. Choose stochastic basis, orthogonal polynomials, suitable for the parameter space
   - See also:
     - [Orthogonal polynomials and recurrence relations](orthogonal_polynomials.md) — Legendre / Hermite polynomials and recurrence coefficients.
     - [ONBasis (one-dimensional orthonormal basis)](onbasis.md) — construction and utilities for evaluating 1D orthonormal polynomial bases (norms, quadrature, evaluations).
     - [Tensorized / multivariate basis (TensorizedBasis)](tonbasis.md) — assembling multivariate bases from ONBasis instances and precomputing triple products.

3. Build spatial FE spaces and blocks of system matrix.
   - Use FESpace types (H1, HDIV, ...) from the ExtendableFEM ecosystem.
   - See: [ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl) — basis finite-element spaces, basis evaluations and low-level FE utilities like standard interpolations.
   - See: [ExtendableFEM.jl](https://github.com/WIAS-PDELib/ExtendableFEM.jl) — high-level deterministic operator assembly and helpers used throughout the codebase.
  
4. Assemble and solve the full SG system
   - Galerkin projection yields a coupled deterministic block system. Briefly:
     - Expand u(y,x)=∑_μ u_μ(x) H_μ(y), test with H_ν ⇒ block matrix with entries a_{μ,ν} A(·).
   - Solver choices:
     - direct assembly + dense solver (for debugging / small problems)
     - matrix‑free iterative solvers exploiting tensor/block structure + preconditioners
   - See problem solver implementations: `src/modelproblems/solvers_*.jl`.

5. Postprocessing, error estimation & adaptivity
   - Error estimators and marking criteria implemented in the script driver `scripts/poisson.jl` for the available Poisson model problems.
   - Error estimators are problem-dependent and can be currently found in `src/estimate.jl`
   - Spatial (mesh refinement) uses refinement routines from [ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl)
   - Stochastic refinement (enrich multi‑index set) uses functions from 
   `src/mopcontrol.jl`
   - see dcomumentation page on [Estimators](estimators.md) for some more details
   - Results, parameters and reproducible outputs are stored with DrWatson (see scripts/poisson.jl for naming pattern).
