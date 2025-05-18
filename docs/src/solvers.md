# Solver

Solving requires a spatial and stochastic discretization.
Both are connected in a special vector structure
that is passed to a solve function that runs a special
iterative solver for each model problem.

## SGFEVector

The spatial discretization is represented by
s single finite element space from [ExtendableFEM.jl](https://github.com/WIAS-PDELib/ExtendableFEM.jl),
while the stochastic discretization is represented by a tensorized basis
for the parameter space of the stochastic coefficient. Both have to be prepared in
advance.

!!! note

    Currently it is not possible to use different finite element spaces for different multi-indices,
    but this feature might be added in the future.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["sgfevector.jl"]
Order   = [:type, :function]
```

## Solve function

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["modelproblems/modelproblems.jl"]
Order   = [:type, :function]
```
