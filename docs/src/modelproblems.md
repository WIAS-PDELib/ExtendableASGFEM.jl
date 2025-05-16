# Model Problems

The following model problems are available and are
dispatched via the following types.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["modelproblems/poisson_primal.jl",
         "modelproblems/logpoisson_primal.jl",
         "modelproblems/logpoisson_dual.jl"]
Order   = [:type, :function]
```

## Solver

Each of the above model problems has an iterative solver
implemented that is used as an default in the following
solve method.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["modelproblems/modelproblems.jl"]
Order   = [:type, :function]
```

## A posteriori error estimator
