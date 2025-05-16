abstract type ChebyshevTPolynomials <: OrthogonalPolynomialType end

recurrence_coefficients(::Type{ChebyshevTPolynomials}, k::Integer) =
  0, (k==0 ? 1 : 2), 1

issymmetric(::Type{ChebyshevTPolynomials}) = true


abstract type ChebyshevUPolynomials <: OrthogonalPolynomialType end

recurrence_coefficients(::Type{ChebyshevUPolynomials}, k::Integer) =
  0, 2, 1

issymmetric(::Type{ChebyshevUPolynomials}) = true