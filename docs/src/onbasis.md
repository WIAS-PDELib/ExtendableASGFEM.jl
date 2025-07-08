# ONBasis

An ONBasis (orthonormal basis) encapsulates information about the orthogonal polynomials associated with a given probability distribution. This includes norms, quadrature rules, and cached evaluations at quadrature points. The ONBasis serves as a fundamental building block for constructing the tensorized basis linked to the multi-indices in the stochastic discretization.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["onbasis.jl"]
Order   = [:type, :function]
```
