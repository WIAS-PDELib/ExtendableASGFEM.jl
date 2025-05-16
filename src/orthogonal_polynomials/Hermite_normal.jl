abstract type HermitePolynomials <: OrthogonalPolynomialType end


recurrence_coefficients(::Type{HermitePolynomials}, k::Integer) =
    0, 1, k

issymmetric(::Type{HermitePolynomials}) = true

norms(::Type{HermitePolynomials}, k) = sqrt.(factorial.(big(k)))

distribution(::Type{HermitePolynomials}) = Normal(0, 1)


# abstract type HermitePolynomialsNormalized <: OrthogonalPolynomialType end

# recurrence_coefficients(::Type{HermitePolynomialsNormalized}, k::Integer) =
#   0, 1/sqrt(k+2), sqrt((k+1)/(k+2))

# issymmetric(::Type{HermitePolynomialsNormalized}) = true

# norms(::Type{HermitePolynomialsNormalized}, k) = 1

# distribution(::Type{HermitePolynomialsNormalized}) = Normal(0,1)
