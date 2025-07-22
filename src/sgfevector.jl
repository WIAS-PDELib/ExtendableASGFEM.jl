"""
$(TYPEDEF)

A structure that extends `ExtendableFEMBase.FEVector` to include information about the stochastic discretization, specifically the associated `TensorizedBasis`.

# Fields
- `FES_space`: Array of finite element spaces for all stochastic modes.
- `TB`: The tensorized basis used for the stochastic discretization.
- `active_modes`: Indices of the active stochastic modes (with respect to `TB.multi_indices`).
- `length4modes`: Offsets for each mode in the global vector.
- `FEVectorBlocks`: Vector block mask for each multi-index (stochastic mode).
- `entries`: The full coefficient vector containing all degrees of freedom.
- `last_sample`: The last sample used for evaluation (for efficient repeated evaluation).
- `FEV`: An `FEVector` used for evaluating the SGFEVector at a given sample.

This structure enables efficient storage, evaluation, and manipulation of solutions in stochastic Galerkin finite element methods.
"""
struct SGFEVector{T, Tv, Ti, ONBType <: ONBasis, MIType}
    FES_space::Array{FESpace{Tv, Ti}, 1}        ## finite element space for all modes
    TB::TensorizedBasis{Tv, ONBType, MIType}    ## tensorized basis with
    active_modes::Array{Ti, 1}                   ## indices of active modes (w.r.t. TB.mulit_indices)
    length4modes::Array{Ti, 1}                   ## offsets for each mode
    FEVectorBlocks::Array{FEVectorBlock{T, Tv, Ti}, 1} ## Vector block mask for each multi-index
    entries::Array{T, 1}                         ## complete coefficient vector
    last_sample::Array{Tv, 1}                    ## last sample used for evaluation
    FEV::FEVector{T, Tv, Ti}                      ## Vector used for evaluating
end

"""
$(TYPEDSIGNATURES)

Evaluates the `SGFEVector` at the given sample `S` and stores the result in `SGFEV.FEV`.

# Arguments
- `SGFEV`: The stochastic Galerkin finite element vector to evaluate.
- `S`: A vector representing the sample (values for the stochastic variables).

# Details
- The sample `S` is stored in `SGFEV.last_sample` (truncated or padded as needed).
- The tensorized basis is evaluated at `S`, and the resulting coefficients are used to assemble the spatial solution in `SGFEV.FEV`.
- This enables efficient evaluation of the SGFEM solution at arbitrary points in the stochastic parameter space.
"""
function set_sample!(SGFEV::SGFEVector, S::AbstractVector)
    last_sample = SGFEV.last_sample
    fill!(last_sample, 0)
    if length(S) <= length(last_sample)
        last_sample[1:length(S)] .= S
    else
        last_sample .= view(S, 1:length(last_sample))
        @warn "ignoring end of sample, since SGFEVector only has multi-indices up to m = $(length(last_sample))"
    end

    ## evaluate SGFEVector at sample
    FEV = SGFEV.FEV
    TB = SGFEV.TB
    nmodes = num_multiindices(SGFEV)
    set_sample!(TB, S; normalize = true)
    nunknowns = length(FEV)
    for u in 1:nunknowns
        fill!(FEV[u], 0)
        for k in 1:nmodes
            r = ExtendableASGFEM.evaluate(TB, k) # = P_k(s)
            if r != 0
                view(FEV[u]) .+= r * view(SGFEV[(u - 1) * nmodes + k])
            end
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Constructs an `SGFEVector` for the given spatial finite element spaces `FES` and the tensorized basis `TB` representing the stochastic discretization.

# Arguments
- `FES`: Array of finite element spaces (one for each unknown).
- `TB`: The tensorized basis for the stochastic discretization.
- `active_modes`: Indices of the active stochastic modes to include (default: all modes in `TB`).
- `T`: The floating-point type for the coefficient vector (default: `Tv`).
- `unames`: Names or identifiers for the unknowns (default: `1:length(FES)`).

# Returns
An `SGFEVector` object that stores the spatial and stochastic discretization, the coefficient vector, and all necessary metadata for efficient evaluation and manipulation in stochastic Galerkin FEM.
"""
function SGFEVector(FES::Array{<:FESpace{Tv, Ti}, 1}, TB::TensorizedBasis{Tv, ONBType, MIType}; active_modes = :all, T = Tv, unames = 1:length(FES)) where {Tv, Ti, ONBType, MIType}
    if active_modes == :all
        active_modes = 1:nmodes(TB)
    end
    length4modes = [0]
    feblocks = Array{FEVectorBlock{T, Tv, Ti}, 1}(undef, 0)
    entries = zeros(T, 0)
    nunknowns = length(FES)
    for j in 1:length(FES)
        fes = FES[j]
        ndofs = fes.ndofs
        append!(entries, zeros(Float64, ndofs * length(active_modes)))
        for m in active_modes
            name = nunknowns == 1 ? "$(TB.multi_indices[m])" : "$(unames[j]) $(TB.multi_indices[m])"
            push!(feblocks, FEVectorBlock{T, Tv, Ti, typeof(entries), eltype(fes), ExtendableFEMBase.assemblytype(fes)}(name, fes, length4modes[end], length4modes[end] + ndofs, entries))
            push!(length4modes, length4modes[end] + ndofs)
        end
    end

    return SGFEVector{T, Tv, Ti, ONBType, MIType}(FES, TB, active_modes, length4modes, feblocks, entries, zeros(Tv, maxlength_multiindices(TB)), FEVector(FES))
end

function SGFEVector(FES::FESpace{Tv, Ti}, TB::TensorizedBasis{Tv}; kwargs...) where {Tv, Ti}
    return SGFEVector(Array{FESpace{Tv, Ti}, 1}([FES]), TB; kwargs...)
end

"""
$(TYPEDSIGNATURES)

returns the `i`-th stochastic mode

"""
Base.getindex(SGFEV::SGFEVector, i::Int) = SGFEV.FEVectorBlocks[i]
"""
$(TYPEDSIGNATURES)

returns the FEVectorBlock for the `i`-th stochastic mode of the `u`-th unknown

"""
Base.getindex(SGFEV::SGFEVector, u::Int, i::Int) = SGFEV.FEVectorBlocks[(u - 1) * length(SGFEV.active_modes) + i]
"""
$(TYPEDSIGNATURES)

returns a tuple with the number of active stochastic modes and the number of spatial degrees of freedom

"""
Base.size(SGFEV::SGFEVector) = (length(SGFEV.active_modes), SGFEV.FES_space.ndofs)

"""
$(TYPEDSIGNATURES)

returns the length of the full vector, i.e., the total number of degrees of freedom

"""
Base.length(SGFEV::SGFEVector) = length(SGFEV.entries)


"""
$(TYPEDSIGNATURES)

returns the number of active modes (that are used from the stored tensorized basis)

"""
num_multiindices(SGFEV::SGFEVector) = length(SGFEV.active_modes)

"""
$(TYPEDSIGNATURES)

returns the finite element types of the (spatial) FE spaces

"""
fetypes(SGFEV::SGFEVector) = [eltype(SGFEV.FES_space[j]) for j in 1:length(SGFEV.FES_space)]

"""
$(TYPEDSIGNATURES)

Custom `show` function for `FEVector` that prints some information on its blocks.
"""
function Base.show(io::IO, SGFEV::SGFEVector)
    println(io, "\nSGFEVector information")
    println(io, "======================")
    println(io, "FETypes = ", fetypes(SGFEV))
    println(io, "")
    ndofs = SGFEV.FES_space[1].ndofs
    nblocks = num_multiindices(SGFEV)
    # Compute blockwise L2 norms
    block_norms = [norm(SGFEV[j]) for j in 1:nblocks]
    total_norm = norm(block_norms)
    # Table header
    print(io, "   block  |                 MI                 |   ndofs  |     min      |     max      |    norm   ")
    print(io, "\n-----------------------------------------------------------------------------------------------------")
    for j in 1:nblocks
        mi = SGFEV.TB.multi_indices[SGFEV.active_modes[j]]
        ext = extrema(view(SGFEV[j]))
        nrm = block_norms[j]
        @printf(
            io, "\n [%5d]  | %s | %8d | %12.4e | %12.4e | %10.4e ",
            j, lpad(string(mi), 34), ndofs, ext[1], ext[2], nrm
        )
    end
    println(io, "\n-----------------------------------------------------------------------------------------------------")
    @printf(io, "Total L2 norm = %.4e\n", total_norm)
    return @printf(io, "Total ndofs = %d\n", length(SGFEV))
end


function norms(SGFEV::SGFEVector{T}, p::Real = 2) where {T}
    nmodes = num_multiindices(SGFEV)
    n = zeros(T, nmodes)
    for j in 1:nmodes
        n[j] = norm(SGFEV[j], p)
    end
    return n
end
