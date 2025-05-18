# TensorizedBasis

Each multi-index ``\mu = [\mu_1,\mu_2,\ldots,\mu_M]``
encodes a tensorized basis function for the parameter space
of the form ``H_\mu = \prod_{k=1}^M H_k`` where the
``H_k`` are the orthogonal polynomials.
The TensorizedBasis collects all information necessary
to evaluate those basis functions, i.e. the set of multi-indices
and the triple products of the form ``(y_mH_\mu, H_\lambda)``
for each ``m`` and ``\mu, \lambda`` in the set of multi-indices
as a sparse matrix. There are analytic formulas to evaluate
these triple products in terms of recurrence coefficients, but it makes
sense to store them for faster evaluation times.




```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["tensorizedbasis.jl"]
Order   = [:type, :function]
```
