# Solver

Solving a problem requires both spatial and stochastic discretization. These are combined into a specialized vector structure, which is then passed to a solver function that executes an iterative algorithm tailored to each model problem.

## SGFEVector

The spatial discretization is defined by a single finite element space from [ExtendableFEM.jl](https://github.com/WIAS-PDELib/ExtendableFEM.jl), while the stochastic discretization uses a tensorized basis for the parameter space of the stochastic coefficient. Both components must be set up in advance.

!!! note
    Currently, it is not possible to use different finite element spaces for different multi-indices. This feature may be added in the future.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["sgfevector.jl"]
Order   = [:type, :function]
```

## Solve Dispatchers

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["modelproblems/modelproblems.jl"]
Order   = [:type, :function]
```

## Poisson Primal Solvers

The primal stochastic Poisson problem can be solved both iteratively (matrix-free) and directly (full assembly). See [Iterative Solution of the Primal Poisson Problem](solvers_poisson_primal.md) for a detailed explanation of the iterative algorithm.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["modelproblems/solvers_poisson_primal.jl"]
Order   = [:type, :function]
```

## Log-Transformed Poisson Solvers

The matrix-free solvers of the log-transformed problems follow the same pattern as the primal solver:
the system operator only implements `LinearAlgebra.mul!`, with the nonzero coupling coefficients of `G`
precomputed into a list so that vanishing couplings are never iterated over.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["modelproblems/solvers_logpoisson_primal.jl", "modelproblems/solvers_logpoisson_dual.jl"]
Order   = [:type, :function]
```

## Initial Guess for Adaptive Loops

In an adaptive refinement loop, the iterative solvers above can start from the solution
of the previous level instead of zero. `interpolate_initial_guess!` transfers the spatial
coefficient blocks of the previously active stochastic modes: they are copied directly if
the grid is unchanged, and pointwise interpolated (`lazy_interpolate!`) onto the dof points
of the new grid after spatial refinement. The transferred vector is passed to the solver as
part of the solution vector; the solver eliminates the guess and lets GMRES solve for the
correction (see [Iterative Solution of the Primal Poisson Problem](solvers_poisson_primal.md)).
The stopping criterion stays the cold-start one, $\|b - \mathbf{A}u\| \le \mathrm{atol} + \mathrm{rtol}\,\|b\|$,
so that iteration counts benefit from a good guess while accuracy does not depend on it.
This is used in the AFEM loop of `scripts/poisson.jl` (see [Poisson Script](poisson_script.md)).
