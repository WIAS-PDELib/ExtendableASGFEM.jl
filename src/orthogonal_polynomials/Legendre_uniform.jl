"""
$(TYPEDEF)

Type for dispatching Legendre Polynomials
"""
abstract type LegendrePolynomials <: OrthogonalPolynomialType end


"""
$(TYPEDSIGNATURES)

Returns the recurrence coefficients for the k-th Legendre polynomial
"""
recurrence_coefficients(::Type{LegendrePolynomials}, k::Integer) =
    0, (2 * k + 1) // (k + 1), k // (k + 1)

issymmetric(::Type{LegendrePolynomials}) = true

"""
$(TYPEDSIGNATURES)

Returns the norm of the k-th Legendre polynomial
"""
norms(::Type{LegendrePolynomials}, k) = sqrt(2 / (2 * k + 1))

"""
$(TYPEDSIGNATURES)

Returns the distribution associated to the Legendre polynomial
"""
distribution(::Type{LegendrePolynomials}) = Uniform(-1, 1)
