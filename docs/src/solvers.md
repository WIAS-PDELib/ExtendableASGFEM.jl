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

## Solve Function

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["modelproblems/modelproblems.jl"]
Order   = [:type, :function]
```
