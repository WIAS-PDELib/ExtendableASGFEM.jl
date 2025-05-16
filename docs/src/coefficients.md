# Stochastic Coefficients

Coefficients are assumed to be represented by a Karhunen-Loeve expansion (KLE)
that have the following general structure and API.

## API
```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["coefficients/coefficients.jl"]
Order   = [:type, :function]
```

## StochasticCoefficientCosinus

An example of a particular KLE with cosinus type basis function
is given by the following subtype.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["coefficients/cosinus.jl"]
Order   = [:type, :function]
```
