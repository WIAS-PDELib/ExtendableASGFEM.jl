"""
$(TYPEDEF)

Type for dispatching Hermite Polynomials
"""
abstract type HermitePolynomials <: OrthogonalPolynomialType end

"""
$(TYPEDSIGNATURES)

Returns the recurrence coefficients for the k-th Legendre polynomial
"""
recurrence_coefficients(::Type{HermitePolynomials}, k::Integer) =
    0, 1, k

issymmetric(::Type{HermitePolynomials}) = true

"""
$(TYPEDSIGNATURES)

Returns the norm of the k-th Hermite polynomial
"""
norms(::Type{HermitePolynomials}, k) = sqrt.(factorial.(big(k)))

"""
$(TYPEDSIGNATURES)

Returns the distribution associated to the Hermite polynomial
"""
distribution(::Type{HermitePolynomials}) = Normal(0, 1)
