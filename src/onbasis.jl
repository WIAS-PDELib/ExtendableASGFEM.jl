"""
$(TYPEDEF)

Structure that stores information of an orthogonal basis, i.e.,
their norms, quadrature rules and values at reference points,
and the result of the last evaluate call.
"""
struct ONBasis{T <: Real, OBT <: OrthogonalPolynomialType, npoly, nquad}
    norms::Vector{T}                                     # norms of all basis functions
    val::Vector{T}                                              # heap for saving last evaluations
    gauss_rule::Tuple{SVector{nquad,T},SVector{nquad,T}}        # quadrature rule points and weights
    vals4xref::SMatrix{nquad,npoly,T}                           # evaluations of all ONB functions at quadrature points
end

vals4qp(ONB::ONBasis, k) = view(ONB.vals4xref,k,:)
vals4poly(ONB::ONBasis, p) = view(ONB.vals4xref,:,p+1)
qw(ONB::ONBasis) = ONB.gauss_rule[2]
qp(ONB::ONBasis) = ONB.gauss_rule[1]
vals4xref(ONB::ONBasis) = ONB.vals4xref
norm4poly(ONB::ONBasis, p) = getindex(ONB.norms, p+1)
distribution(ONB::ONBasis{T,OBT}) where {T,OBT} = distribution(OBT)
OrthogonalPolynomialType(ONB::ONBasis{T,OBT}) where {T,OBT} = OBT


"""
$(TYPEDSIGNATURES)

Constructs an ONBasis for the given OrthogonalPolynomialType (and the associated weight function)
and the maximal polynomial order.
"""
function ONBasis(OBT::Type{<:OrthogonalPolynomialType}, maxorder, maxquadorder = 2*maxorder; T = Float64)
    ## find quadrature rule and evaluate basis at all quadrature points
    gr = gauss_rule(OBT, maxquadorder; T = T)
    val = zeros(T, maxorder+1)
    vals4xref = evaluate(OBT, maxorder, gr[1])
    
    ## normalize
    qw = gr[2]
    nqp = length(qw)
    norms = zeros(T, maxorder+1) #norm(OBT, maxorder; gr = gr).^(-1/2)
    for m = 1 : maxorder+1, k = 1 : nqp
        norms[m] += vals4xref[k,m]^2 * qw[k]
    end
    norms = sqrt.(norms)

    return ONBasis{T, OBT, maxorder+1, length(gr[1])}(norms, val, gr, vals4xref)
end

function _triple_product(vals4xref,w,j,k,l)
    T = eltype(vals4xref)
    val::T = 0.0
    for q = 1 : length(w)
        val += vals4xref[q,j+1] * vals4xref[q,k+1] * vals4xref[q,l+1] * w[q]
    end
    return val
end

"""
$(TYPEDSIGNATURES)

Evaluates the triple product of three basis functions.
"""
function triple_product(ONB::ONBasis, j, k, l; normalize = true)
    if normalize
        return _triple_product(vals4xref(ONB), qw(ONB), j, k, l) / norm4poly(ONB, j) / norm4poly(ONB, k) / norm4poly(ONB, l)
    else
        return _triple_product(vals4xref(ONB), qw(ONB), j, k, l)
    end
end

function _triple_product_y(vals4xref,w,qp,j,k)
    T = eltype(vals4xref)
    val::T = 0.0
    for q = 1 : length(w)
        val += vals4xref[q,j+1] * vals4xref[q,k+1] * qp[q] * w[q]
    end
    return val
end

"""
$(TYPEDSIGNATURES)

Evaluates the triple product of two basis functions and y.
"""
function triple_product_y(ONB::ONBasis, j, k; normalize = true)
    if normalize
        return _triple_product_y(vals4xref(ONB), qw(ONB), qp(ONB), j, k) / norm4poly(ONB, j) / norm4poly(ONB, k)
    else
        return _triple_product_y(vals4xref(ONB), qw(ONB), qp(ONB), j, k)
    end
end

function _scalar_product(vals4xref,w,j,k)
    T = eltype(vals4xref)
    val::T = 0.0
    for q = 1 : length(w)
        val += vals4xref[q,j+1] * vals4xref[q,k+1] * w[q]
    end
    return val
end

"""
$(TYPEDSIGNATURES)

Evaluates the scalar product of two basis functions.
"""
function scalar_product(ONB::ONBasis, j, k; normalize = true)
    if normalize
        return _scalar_product(vals4xref(ONB), qw(ONB), j, k) / norm4poly(ONB, j) / norm4poly(ONB, k) 
    else
        return _scalar_product(vals4xref(ONB), qw(ONB), j, k)
    end
end


function _integral(vals4xref,w,j)
    T = eltype(vals4xref)
    val::T = 0.0
    for q = 1 : length(w)
        val += vals4xref[q,j+1] * w[q]
    end
    return val
end

"""
$(TYPEDSIGNATURES)

Evaluates the integral of a single basis function.
"""
function integral(ONB::ONBasis, j; normalize = true)
    if normalize
        return _integral(vals4xref(ONB), qw(ONB), j) / norm4poly(ONB, j)
    else
        return _integral(vals4xref(ONB), qw(ONB), j)
    end
end

"""
$(TYPEDSIGNATURES)

Evaluates all basis polynomials at point x
"""
function evaluate(ONB::ONBasis{T, OBT, maxorder}, x; normalize = true) where {T,OBT,maxorder}
    evaluate!(ONB.val, OBT, maxorder-1, x)
    if normalize
        ONB.val ./= ONB.norms
    end
    return ONB.val
end