# Adaptivity and error control

## Residual-based a posteriori error estimation
Spatial error estimation refers to classical residual-based error estimation for the zero-th multi-index, which corresponds to the mean value.
Stochastic error control determines which stochastic mode should be refined, either by increasing the polynomial degree or by activating neighboring stochastic modes. Both are represented by multi-indices.
Unified error control enables residual-based error estimation for the subresiduals associated with each multi-index, depending on the model problem (see references on the main page for details).


```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["estimate.jl"]
Order   = [:type, :function]
```

## Multi-index management

Depending on the model problem and stochastic coefficient, the number of multi-indices to be added for error estimation varies.
The following methods assist in enriching the set of multi-indices.


```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["mopcontrol.jl"]
Order   = [:type, :function]
```

## Monte Carlo sampling estimator

A hierarchical Monte Carlo error estimator is also available. It compares the solution with a higher-order discrete solution for sampled deterministic problems. This is mainly intended to compute a reference error to assess the efficiency of the residual-based error estimator.


```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["sampling_error.jl"]
Order   = [:type, :function]
```
