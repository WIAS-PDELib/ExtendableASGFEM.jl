
"""
$(TYPEDEF)

Structure that stores information of a tensorized
orthogonal for a certain set of multi-indices.
"""
struct TensorizedBasis{T <: Real, ONBType <: ONBasis, MIType}
    ONB::ONBType
    vals::Vector{Vector{T}}                 # heap for saving last evaluations of M samples
    nmodes::Int
    G::ExtendableSparseMatrix{T,Int64}      ## saves coupling indices
    multi_indices::MIType                   # e.g. Array{Array{Int,1},1}
end

maxlength(A::Array{Array{T,1},1}) where {T} = maximum([length(A[j]) for j in 1 : length(A)])
num_multiindices(TB::TensorizedBasis) = TB.nmodes
maxlength_multiindices(TB::TensorizedBasis{T, ONBT}) where {T, ONBT} = maxlength(TB.multi_indices)
distribution(TB::TensorizedBasis{T, ONBT}) where {T, ONBT} = distribution(TB.ONB)
get_multiindex(TB::TensorizedBasis, j) = TB.multi_indices[j]
get_coupling_coefficient(TB::TensorizedBasis{T}, m, j, k) where {T} = TB.G[(m-1)*TB.nmodes+j,k]::T # ⟨ ξ_m ψ_mi(j) ψ_mi(k) ⟩


"""
$(TYPEDSIGNATURES)

Generate samples and weights
for the distribution of the ON basis that the tensorized basis is based upon
(that can be used for a Monte-Carlo estimator).
"""
function sample_distribution(TB::TensorizedBasis, nsamples; M = maxlength_multiindices(TB), Mweights = M)
    dist = distribution(TB)
    Samples::Array{Float64,2} = zeros(Float64,M,nsamples)
    Random.seed!(123)
    rand!(dist, Samples)
    weights = [prod([pdf(dist, Samples[j,s]) for j = 1 : Mweights]) for s = 1 : nsamples]
    return Samples, weights
end


"""
$(TYPEDSIGNATURES)

shows information on the tensorized basis.
"""
function Base.show(io::IO, TB::TensorizedBasis)
    println(io, "TENSORIZED BASIS")
    println(io, "nmodes = $(size(TB.multi_indices))")
    println(io, "")
    for j = 1 : min(40, length(TB.multi_indices))
        print(io, "$(TB.multi_indices[j]) ")
        if j % 4 == 0
            println(io, "");
        end
    end
    if length(TB.multi_indices) > 50
        println(io, "...")
    end

    ONB = TB.ONB
    qp = ExtendableASGFEM.qp(ONB)
    qw = ExtendableASGFEM.qw(ONB)
    norms = ONB.norms
    vals = ONB.vals4xref
    println(io, "ON-BASIS")
    println(io, "type = $(typeof(ONB))")
    println(io, "")
    #for j = 1 : length(qp)
    #    println(io, "$(qp[j]) $(vals[j,:])")
    #end
    for j = 1 : length(norms)
        println(io, "||H_$(j-1)||_ω =$(norms[j])")
    end
end

"""
$(TYPEDSIGNATURES)

constructor for a tensorized basis for the given OrthogonalPolynomialType. If no multi-indices are provided it automatically generates
all multi-indices up to support length M and polynomial order maxorder.
"""
function TensorizedBasis(OBT::Type{<:OrthogonalPolynomialType}, M, order, maxorder, maxquadorder = 2*maxorder; T = Float64, multi_indices = "full")
    @assert order <= maxorder
    ONB = ONBasis(OBT, maxorder, maxquadorder; T = T)
    vals = Vector{Vector{T}}(undef, M)
    for j = 1 : M
        vals[j] = zeros(T, length(ONB.val))
    end
    if multi_indices == "full"
        multi_indices = generate_multiindices(M, order)
    end

    nmodes = length(multi_indices)

    G = get_tensor_multiplication_with_ym(T, multi_indices, ONB)

    return TensorizedBasis{T, typeof(ONB), typeof(multi_indices)}(ONB, vals, nmodes, G, multi_indices)
end


function multiplication_with_ym(TB::TensorizedBasis, j, k, m; normalize = true)
    val = 1.0
    ONB = TB.ONB
    multi_indices = TB.multi_indices
    #for d = 1 : maxlength_multiindices(TB)
    #    if d == m
    #        val *= triple_product_y(ONB, multi_indices[j][d], multi_indices[k][d]; normalize = normalize)
    #    else
    #        val *= scalar_product(ONB, multi_indices[j][d], multi_indices[k][d]; normalize = normalize)
    #    end
    #end
    #val = triple_product_y(ONB, multi_indices[j][m], multi_indices[k][m]; normalize = normalize)
    if multi_indices[k][m] == multi_indices[j][m] + 1
        val = norm4poly(ONB, multi_indices[j][m] + 1) / norm4poly(ONB, multi_indices[j][m])
    elseif multi_indices[k][m] == multi_indices[j][m] - 1
        if multi_indices[j][m] == 0
            val = 0
        else
            val = (multi_indices[j][m] - 1) * norm4poly(ONB, multi_indices[j][m] - 1) / norm4poly(ONB, multi_indices[j][m])
        end
    end
    
    return val
end

function triple_product(TB::TensorizedBasis, j, k, l; normalize = true)
    val = 1.0
    ONB = TB.ONB
    multi_indices = TB.multi_indices
    for d = 1 : maxlength_multiindices(TB)
        val *= triple_product(ONB, multi_indices[j][d], multi_indices[k][d], multi_indices[l][d]; normalize = normalize)
    end
    return val
end

function get_tensor_multiplication_with_ym(T, multi_indices, ONB)
    #a = zeros(Int,M)
    M::Int = maxlength(multi_indices)
    nmodes::Int = length(multi_indices)
    OBT = OrthogonalPolynomialType(ONB)
    G = ExtendableSparseMatrix{T,Int}(M*nmodes,nmodes) # quick and dirty solution for a 3D sparse array, G[a,b,c] = G[(a-1)*M+b,c]

    prod::T = 0.0
    for j = 1 : nmodes
        for m = 1 : M
            (a,b,c) = normalise_recurrence_coefficients(OBT, multi_indices[j][m])
            mu1 = deepcopy(multi_indices[j])
            mu2 = deepcopy(multi_indices[j])
            mu1[m] += 1
            mu2[m] -= 1
            for k = 1 : nmodes
                if all(multi_indices[k] .== mu1)
                    #@info "MI $k=$(multi_indices[k]) is +1 neigbour of MI $j=$(multi_indices[j]) with factor $(1/b)"
                    G[(m-1)*nmodes+j,k] = 1/b # = beta + 1
                elseif all(multi_indices[k] .== mu2)
                    #@info "MI $k=$(multi_indices[k]) is -1 neigbour of MI $j=$(multi_indices[j]) with factor $(c/b)"
                    G[(m-1)*nmodes+j,k] = c/b # = beta - 1
                end
            end
        end
    end
    #@assert false
    flush!(G)
    return G
end

## todo : rewrite
function get_tensor(TB::TensorizedBasis{T, ONBT}) where {T, ONBT}
    @assert false # never use this function
    #a = zeros(Int,M)
    multi_indices = TB.multi_indices
    M::Int = maxlength_multiindices(TB)
    nmodes::Int = length(multi_indices)
    G = ExtendableSparseMatrix{Float64,Int}(M*nmodes,nmodes) # quick and dirty solution for a 3D sparse array, G[a,b,c] = G[(a-1)*M+b,c]

    prod::T = 0.0
    for j = 1 : nmodes, k = 1 : nmodes, e = 1 : M
        prod = triple_product(TB, e, j, k; normalize = true)
        if abs(prod) > 1e-12
            G[(e-1)*nmodes+j,k] = prod
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
    for j = 1 : length(x)
        vals[j] .= evaluate(TB.ONB, x[j]; normalize = normalize)
    end
    for j = length(x)+1:length(vals)
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
    for m = 1 : length(multi_index)
        prod *= uni_vals[m][multi_index[m]+1]
    end
    return prod
end

