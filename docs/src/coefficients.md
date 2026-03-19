# Stochastic Coefficients

Stochastic coefficients play a central role in uncertainty quantification and stochastic finite element methods. In this package, coefficients are represented using a Karhunen-Loève expansion (KLE), which allows for the efficient representation of random fields with prescribed covariance structure.

## Overview

The Karhunen-Loève expansion expresses a stochastic process as a series of orthogonal functions weighted by uncorrelated random variables:

```math
a(x, \omega) = a_0(x) + \sum_{n=1}^N \sqrt{\lambda_n}\, \phi_n(x)\, \xi_n(\omega)
```

where $a_0(x)$ is the mean, $\lambda_n$ and $\phi_n(x)$ are the eigenvalues and eigenfunctions of the covariance operator, and $\xi_n$ are independent standard normal random variables.

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
