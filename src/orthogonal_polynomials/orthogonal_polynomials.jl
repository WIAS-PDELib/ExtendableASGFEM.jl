### taken/modified from : https://github.com/aleadev/orthopoly.jl

abstract type OrthogonalPolynomialType end

include("Hermite_normal.jl")
include("Legendre_uniform.jl")

eltypes(::Type{Tuple{T1}}) where {T1} = (T1,)
eltypes(::Type{Tuple{T1, T2}}) where {T1, T2} = (T1, T2)
eltypes(::Type{Tuple{T1, T2, T3}}) where {T1, T2, T3} = (T1, T2, T3)
eltypes(::Type{Tuple{T1, T2, T3, T4}}) where {T1, T2, T3, T4} = (T1, T2, T3, T4)

function flip_tuple_array(a::Array{T}) where {T}
    types = eltypes(T)
    N = length(a)
    b = [zeros(t, N) for t in types]
    for i in 1:length(types)
        for j in 1:N
            b[i][j] = a[j][i]
        end
    end
    return (b)
end


function rc_array(basis::Type{<:OrthogonalPolynomialType}, n::Integer)
    arr = [recurrence_coefficients(basis, k) for k in 0:n]
    return flip_tuple_array(arr)
end

function rc_array_monic(basis::Type{<:OrthogonalPolynomialType}, n::Integer)
    r = rc_array(basis, n)
    # extract columns
    a = -r[1]
    b = r[3][2:end]
    c = r[2][:]

    # convert to monic polynomials
    α = a ./ c
    β = b ./ (c[1:(end - 1)] .* c[2:end])
    return α, β
end


"""
$(TYPEDSIGNATURES)

Evaluates the first n+1 orthogonal polynomials at a single x and write result into y
"""
function evaluate!(y, basis::Type{<:OrthogonalPolynomialType}, n::Integer, x::Real)
    y[1] = 1
    for k in 0:(n - 1)
        a, b, c = recurrence_coefficients(basis, k)
        y[k + 2] = (a + b * x) * y[k + 1]
        if k > 0
            y[k + 2] = y[k + 2] - c * y[k]
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Evaluates the first n+1 orthogonal polynomials at a single x
"""
function evaluate(basis::Type{<:OrthogonalPolynomialType}, n::Integer, x::Real)
    y = zeros(typeof(x), n + 1)
    y[1] = 1
    for k in 0:(n - 1)
        a, b, c = recurrence_coefficients(basis, k)
        y[k + 2] = (a + b * x) * y[k + 1]
        if k > 0
            y[k + 2] = y[k + 2] - c * y[k]
        end
    end
    return y
end

"""
$(TYPEDSIGNATURES)

Evaluates the first n+1 orthogonal polynomials at a vector of x
"""
function evaluate(basis::Type{<:OrthogonalPolynomialType}, n::Integer, x::AbstractVector{T}) where {T}
    y = ones(T, length(x), n + 1)
    for k in 0:(n - 1)
        a, b, c = recurrence_coefficients(basis, k)
        for j in 1:length(x)
            y[j, k + 2] = (a + b * x[j]) * y[j, k + 1]
            if k > 0
                y[j, k + 2] = y[j, k + 2] - c * y[j, k]
            end
        end
    end
    return y
end


"""
$(TYPEDSIGNATURES)

computes the norms of the first n+1 polynomials by Gauss quadrature
"""
function LinearAlgebra.norm(basis::Type{<:OrthogonalPolynomialType}, n; gr = gauss_rule(basis, 2 * n))

    function norm_internal!(norms, qw, y)
        for m in 1:length(norms)
            norms[m] = dot(view(y, :, m), qw)
        end
        return
    end

    T = eltype(gr[1])
    y = evaluate(basis, n, gr[1]) .^ 2
    norms = zeros(T, n + 1)
    norm_internal!(norms, gr[2], y)
    return sqrt.(norms)
end

function gauss_rule(basis::Type{<:OrthogonalPolynomialType}, n::Integer; T = Float64)
    α, β = rc_array_monic(basis, n - 1)
    α = Vector{T}(α)
    β = Vector{T}(β)

    return gauss_rule((α, β))
end

function gauss_rule(monic_rc::Tuple{Vector{S}, Vector{T}}) where {S, T}
    # W. GAUTSCHI, ORTHOGONAL POLYNOMIALS AND QUADRATURE, Electronic
    # Transactions on Numerical Analysis, Volume 9, 1999, pp. 65-76.

    α = monic_rc[1]
    β = monic_rc[2]

    # set up Jacobi matrix and compute eigenvalues
    J = SymTridiagonal(α, sqrt.(β))
    #J = float(diagm(α) + diagm(√β,1) + diagm(√β,-1))

    x = eigvals(J)
    V = eigvecs(J)
    w = vec(V[1, :]' .^ 2)
    w = w / sum(w)

    # symmetrise
    if all(α == 0)
        x = 0.5 * (x - reverse(x))
        w = 0.5 * (w + reverse(w))
    end
    return (x, w)
end

function sqnorms_from_rc(OBT::Type{<:OrthogonalPolynomialType}, k)
    (a0, b0, c0) = recurrence_coefficients(OBT, 0)
    (an, bn, cn) = recurrence_coefficients(OBT, k)
    h = b0 / bn
    for i in 1:(k + 1)
        (ai, bi, ci) = recurrence_coefficients(OBT, i)
        h *= ci
    end
    return h
end


"""
$(TYPEDSIGNATURES)

Changes the recurrence coefficients, such that basis functions are normalized.
"""
function normalise_recurrence_coefficients(OBT::Type{<:OrthogonalPolynomialType}, k)
    (a, b, c) = recurrence_coefficients(OBT, k)
    h2 = norms(OBT, k + 1)
    h1 = norms(OBT, k)
    if k > 0
        h0 = norms(OBT, k - 1)
    else
        h0 = 0
    end
    return (a * h1 / h2, b * h1 / h2, c * h0 / h2)
end
