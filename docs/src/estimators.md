# Adaptivity and error control

## Residual-based a posteriori error estimation
Spatial error estimation refers to classical residual-based error estimation
for the zero-th multi index that refers to the mean value.
Stochastic error control has to estimate which stochastic mode needs to be refined
in the sense that either the polynomial degree is increased or neighbouring
stochastic modes are activated. Both is represented by the multi-indices.
Then unified error control allows to perform residual-based error estimation for
the subresiduals that are associated to each multi-index and depend on the model
problem (see references on the main page for details).


```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["estimate.jl"]
Order   = [:type, :function]
```

## Multi index management

Depending on the model problem and stochastic coefficient the
amount of multi indices that should be added to the error estimation
varies.
Here are some methods that help with enriching the set of multi-indices.


```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["mopcontrol.jl"]
Order   = [:type, :function]
```

## Monte carlo sampling estimator

There is also a hierarchical Monte carlo error estimator available that
compares the solution with a higher order discrete solution for sampled
deterministic problems. This is merely intended as a way to compute the
reference error to assess the efficiency of the residual-based error estimator.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["sampling_error.jl"]
Order   = [:type, :function]
```
