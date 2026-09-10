"""
$(TYPEDEF)

Type for dispatching Legendre polynomials, which are orthogonal with respect to the uniform distribution on [-1, 1], i.e., the weight function `w = 0.5`.
"""
abstract type LegendrePolynomials <: OrthogonalPolynomialType end

"""
$(TYPEDSIGNATURES)

Returns the recurrence coefficients `(a, b, c)` for the k-th Legendre polynomial, corresponding to the three-term recurrence relation:

    P_{k+1}(x) = (a_k + b_k * x) P_k(x) - c_k P_{k-1}(x)

For Legendre polynomials:
- `a = 0`
- `b = (2k + 1) / (k + 1)`
- `c = k / (k + 1)`
- `P_{-1} = 0`
- `P_0 = 1`
"""
recurrence_coefficients(::Type{LegendrePolynomials}, k::Integer) =
    0, (2 * k + 1) / (k + 1), k / (k + 1)

issymmetric(::Type{LegendrePolynomials}) = true

"""
$(TYPEDSIGNATURES)

Returns the norm of the k-th Legendre polynomial, i.e.,

    ||P_k|| = sqrt(1 / (2k + 1))
"""
norm_basis(::Type{LegendrePolynomials}, k) = sqrt(1 / (2 * k + 1))

"""
$(TYPEDSIGNATURES)

Returns the probability distribution associated with the Legendre polynomials, which is the uniform distribution on [-1, 1].
"""
distribution(::Type{LegendrePolynomials}) = Uniform(-1, 1)
