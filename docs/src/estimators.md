# A posteriori error estimation

## Spatial refinement
Spatial error estimation refers to classical residual-based error estimation
for the zero-th multi index that refers to the mean value.
Stochastic error control has to estimate which stochastic mode needs to be refined
in the sense that either the polynomial degree is increased or neighbouring
stochastic modes are activated. Both is represented by the multi-indices.
Then unified error control allows to perform residual-based error estimation for
the subresiduals that are associated to each multi-index
(see references on the main page for details).


```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["estimators.jl"]
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