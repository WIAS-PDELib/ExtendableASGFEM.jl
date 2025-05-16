abstract type LaguerrePolynomials{α <: Real} <: OrthogonalPolynomialType end

recurrence_coefficients(::Type{LaguerrePolynomials{α}}, k::Integer) where {α <: Real} =
    (2k + 1 + α) / (k + 1),
    -1 / (k + 1),
    (k + α) / (k + 1)

issymmetric(::Type{<:LaguerrePolynomials}) = false
