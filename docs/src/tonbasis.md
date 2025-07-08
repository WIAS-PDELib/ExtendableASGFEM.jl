# TensorizedBasis

Each multi-index $\mu = [\mu_1, \mu_2, \ldots, \mu_M]$ defines a tensorized basis function for the parameter space of the form $H_\mu = \prod_{k=1}^M H_k$, where each $H_k$ is an orthogonal polynomial.

The `TensorizedBasis` object collects all information required to evaluate these basis functions, including the set of multi-indices and the triple products of the form $(y_m H_\mu, H_\lambda)$ for each $m$ and $\mu, \lambda$ in the set of multi-indices, stored as a sparse matrix. While analytic formulas exist to compute these triple products using recurrence coefficients, storing them can significantly speed up evaluations.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["tensorizedbasis.jl"]
Order   = [:type, :function]
```
