"""
$(TYPEDEF)

A structure representing a tensorized orthogonal polynomial basis for stochastic Galerkin methods, associated with a given set of multi-indices.

# Fields
- `ONB::ONBType`: The underlying univariate orthogonal basis (e.g., Hermite, Legendre).
- `vals::Vector{Vector{T}}`: Cached evaluations of the univariate basis functions at the most recent sample(s).
- `nmodes::Int`: Number of multi-indices (i.e., the number of stochastic modes).
- `G::ExtendableSparseMatrix{T, Int64}`: Coupling matrix encoding triple products or recurrence relations between basis functions.
- `multi_indices::MIType`: Collection of multi-indices (typically an array of integer arrays), each representing a multi-dimensional polynomial degree.

# Description
The `TensorizedBasis` type encapsulates all information required to evaluate and manipulate a tensor-product polynomial basis in multiple stochastic dimensions. It supports efficient evaluation, storage of basis values, and access to coupling coefficients for use in stochastic Galerkin finite element methods (SGFEM).

# Example
```julia
# Construct a tensorized Hermite basis with 3 variables and total degree 2
TB = TensorizedBasis(HermitePolynomials, 3, 2, 2)

# Evaluate all basis functions at a sample point
set_sample!(TB, [0.1, -0.2, 0.3])

# Get the value of the 5th basis function at the current sample
val = evaluate(TB, 5)
```
"""
struct TensorizedBasis{T <: Real, ONBType <: ONBasis, MIType}
    ONB::ONBType
    vals::Vector{Vector{T}}                 # heap for saving last evaluations of M samples
    nmodes::Int
    G::ExtendableSparseMatrix{T, Int64}      ## saves coupling indices
    multi_indices::MIType                   # e.g. Array{Array{Int,1},1}
end

maxlength(A::Array{Array{T, 1}, 1}) where {T} = maximum([length(A[j]) for j in 1:length(A)])

"""
$(TYPEDSIGNATURES)

returns the number of multi-indices
"""
num_multiindices(TB::TensorizedBasis) = TB.nmodes

"""
$(TYPEDSIGNATURES)

returns the maximal length of the stored multi-indices
"""
maxlength_multiindices(TB::TensorizedBasis{T, ONBT}) where {T, ONBT} = maxlength(TB.multi_indices)

"""
$(TYPEDSIGNATURES)

returns distribution associated to the orthogonal basis
"""
distribution(TB::TensorizedBasis{T, ONBT}) where {T, ONBT} = distribution(TB.ONB)


"""
$(TYPEDSIGNATURES)

returns the j-th multi-index.
"""
get_multiindex(TB::TensorizedBasis, j) = TB.multi_indices[j]


"""
$(TYPEDSIGNATURES)

returns the triple products between ``y_m`` and ``H_j`` and ``H_k``
for two multi-indices `j` and `k`.
"""
get_coupling_coefficient(TB::TensorizedBasis{T}, m, j, k) where {T} = TB.G[(m - 1) * TB.nmodes + j, k]::T # ⟨ ξ_m ψ_mi(j) ψ_mi(k) ⟩


"""
$(TYPEDSIGNATURES)

Generate samples and weights
for the distribution of the ON basis that the tensorized basis is based upon
(that can be used for a Monte-Carlo estimator).
"""
function sample_distribution(TB::TensorizedBasis, nsamples; M = maxlength_multiindices(TB), Mweights = M)
    dist = distribution(TB)
    Samples::Array{Float64, 2} = zeros(Float64, M, nsamples)
    Random.seed!(123)
    rand!(dist, Samples)
    weights = [prod([pdf(dist, Samples[j, s]) for j in 1:Mweights]) for s in 1:nsamples]
    return Samples, weights
end


"""
$(TYPEDSIGNATURES)

shows information on the tensorized basis.
"""
function Base.show(io::IO, TB::TensorizedBasis)
    # Summary section
    nvars = maxlength_multiindices(TB)
    nmodes = num_multiindices(TB)
    dist = distribution(TB)
    Gnnz = length(TB.G.cscmatrix.nzval)
    println(io, "TENSORIZED BASIS")
    println(io, "----------------")
    println(io, "Number of modes: $nmodes")
    println(io, "Max multi-index length: $nvars")
    println(io, "Underlying basis: $(typeof(TB.ONB))")
    println(io, "Distribution: $dist")
    println(io, "Coupling matrix nonzeros: $Gnnz\n")

    # Multi-indices display
    println(io, "Multi-indices:")
    nshow = 8
    nend = 2
    total = length(TB.multi_indices)
    function print_mi_row(startidx, endidx)
        for j in startidx:endidx
            mi = TB.multi_indices[j]
            print(io, lpad(string(mi), 12))
            if (j - startidx + 1) % 4 == 0
                println(io, "")
            end
        end
        return
    end
    if total <= nshow + nend
        print_mi_row(1, total)
        println(io, "")
    else
        print_mi_row(1, nshow)
        println(io, "  ...")
        print_mi_row(total - nend + 1, total)
        println(io, "")
    end

    # ON-BASIS section
    ONB = TB.ONB
    norms = getfield(ONB, :norms)
    println(io, "ON-BASIS")
    println(io, "--------")
    println(io, "Type: $(typeof(ONB))")
    print(io, "Norms: ")
    maxn = min(length(norms), 8)
    for j in 1:maxn
        @printf(io, "%g", norms[j])
        if j < maxn
            print(io, ", ")
        end
    end
    if length(norms) > maxn
        print(io, ", ...")
    end
    println(io, "")
    return
end

"""
$(TYPEDSIGNATURES)

constructor for a tensorized basis for the given OrthogonalPolynomialType. If no multi-indices are provided it automatically generates
all multi-indices up to support length M and polynomial order maxorder.
"""
function TensorizedBasis(OBT::Type{<:OrthogonalPolynomialType}, M, order, maxorder, maxquadorder = 2 * maxorder; T = Float64, multi_indices = "full")
    @assert order <= maxorder
    ONB = ONBasis(OBT, maxorder, maxquadorder; T = T)
    vals = Vector{Vector{T}}(undef, M)
    for j in 1:M
        vals[j] = zeros(T, length(ONB.val))
    end
    if multi_indices == "full"
        multi_indices = generate_multiindices(M, order)
    end

    nmodes = length(multi_indices)

    G = get_tensor_multiplication_with_ym(T, multi_indices, ONB)

    return TensorizedBasis{T, typeof(ONB), typeof(multi_indices)}(ONB, vals, nmodes, G, multi_indices)
end

function triple_product(TB::TensorizedBasis, j, k, l; normalize = true)
    val = 1.0
    ONB = TB.ONB
    multi_indices = TB.multi_indices
    for d in 1:maxlength_multiindices(TB)
        val *= triple_product(ONB, multi_indices[j][d], multi_indices[k][d], multi_indices[l][d]; normalize = normalize)
    end
    return val
end

function get_tensor_multiplication_with_ym(T, multi_indices, ONB)
    #a = zeros(Int,M)
    M::Int = maxlength(multi_indices)
    nmodes::Int = length(multi_indices)
    OBT = OrthogonalPolynomialType(ONB)
    G = ExtendableSparseMatrix{T, Int}(M * nmodes, nmodes) # quick and dirty solution for a 3D sparse array, G[a,b,c] = G[(a-1)*M+b,c]

    for j in 1:nmodes
        for m in 1:M
            (a, b, c) = normalise_recurrence_coefficients(OBT, multi_indices[j][m])
            mu1 = deepcopy(multi_indices[j])
            mu2 = deepcopy(multi_indices[j])
            mu1[m] += 1
            mu2[m] -= 1
            for k in 1:nmodes
                if all(multi_indices[k] .== mu1)
                    G[(m - 1) * nmodes + j, k] = 1 / b # = beta + 1
                elseif all(multi_indices[k] .== mu2)
                    G[(m - 1) * nmodes + j, k] = c / b # = beta - 1
                end
            end
        end
    end
    flush!(G)
    return G
end

"""
$(TYPEDSIGNATURES)

evaluates all basis functions at sample vector x; call this
before using evaluate!
"""
function set_sample!(TB::TensorizedBasis, x; normalize = true)
    vals = TB.vals
    for j in 1:length(x)
        vals[j] .= evaluate(TB.ONB, x[j]; normalize = normalize)
    end
    for j in (length(x) + 1):length(vals)
        fill!(vals[j], 0)
        vals[j][1] = 1
    end
    return TB.vals
end


"""
$(TYPEDSIGNATURES)

evaluates the basis function for the j-th multi-index at the sample that was set with set_sample.
"""
function evaluate(TB::TensorizedBasis{T}, j) where {T}
    uni_vals = TB.vals
    multi_index = TB.multi_indices[j]
    prod::T = 1
    for m in 1:length(multi_index)
        prod *= uni_vals[m][multi_index[m] + 1]
    end
    return prod
end
