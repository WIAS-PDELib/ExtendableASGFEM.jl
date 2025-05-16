abstract type LegendrePolynomials <: OrthogonalPolynomialType end

recurrence_coefficients(::Type{LegendrePolynomials}, k::Integer) =
  0, (2*k+1) // (k+1), k // (k+1)

issymmetric(::Type{LegendrePolynomials}) = true

norms(::Type{LegendrePolynomials}, k) = sqrt(2/(2*k+1))

distribution(::Type{LegendrePolynomials}) = Uniform(-1,1)