# Orthogonal Polynomials

In the stochastic discretization of random variables, global polynomials that are orthogonal with respect to the probability distribution of each random variable $y_m$ are used. These orthogonal polynomials can be generated via recurrence relations, with coefficients determined by the underlying distribution.

## Recurrence Relations
Orthogonal polynomials $H_n$ with respect to a weight function $\omega$ satisfy

```math
\int_\Gamma \omega(y) H_{n}(y) H_m(y) \,dy = N^2_{nm}\delta_{nm}
```
where
```math
N_{nn} := \| H_n \|_{\omega}^2 := \int_\Gamma \omega(y) H_{n}(y) H_n(y) \,dy
```
The polynomials satisfy the three-term recurrence relation:
```math
\begin{aligned}
  H_{n+2}(y) & = (a_n y - b_n) H_{n+1}(y) - c_n H_{n}(y)
\end{aligned}
```
with initial values $H_0 = 0$ and $H_1 = 1$.

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["orthogonal_polynomials/orthogonal_polynomials.jl"]
Order   = [:type, :function]
```

## Legendre Polynomials (Uniform Distribution)

For the weight function $\omega(y) = 1/2$ on the interval $[-1,1]$ (uniform distribution), the recurrence coefficients are $a_n = (2n+1)/(n+1)$, $b_n = 0$, and $c_n = n/(n+1)$. The norms of the resulting Legendre polynomials are
```math
    \| H_n \|^2_\omega = \frac{2}{2n+1}
```

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["orthogonal_polynomials/Legendre_uniform.jl"]
Order   = [:type, :function]
```

## Hermite Polynomials (Normal Distribution)

For the weight function $\omega(y) = \exp(-y^2/2)/(2\pi)$ (normal distribution), the recurrence coefficients are $a_n = 1$, $b_n = 0$, and $c_n = n$. The first six Hermite polynomials are:
```math
\begin{aligned}
H_0 & = 0\\
H_1 & = 1\\
H_2 & = y\\
H_3 & = y^2 - 1\\
H_4 & = y^3 - 3y\\
H_5 & = y^4 - 6y^2 + 3\\
\end{aligned}
```
Their norms are given by
```math
    \| H_n \|^2_\omega = n!
```

```@autodocs
Modules = [ExtendableASGFEM]
Pages = ["orthogonal_polynomials/Hermite_normal.jl"]
Order   = [:type, :function]
```
