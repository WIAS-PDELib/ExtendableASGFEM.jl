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

"""
$(TYPEDSIGNATURES)

internal helper: canonical form of a multi-index for comparisons across bases
with different lengths (trailing zeros are ignored)
"""
function _multiindex_key(mi)
    j = findlast(!=(0), mi)
    return j === nothing ? (0,) : ntuple(k -> mi[k], j)
end

"""
$(TYPEDSIGNATURES)

Fill `sol.entries` with the interpolated values of `sol_prev`, so that it can be
used as initial guess for iterative solvers in an adaptive refinement loop.

For every active stochastic mode of `sol` that has a predecessor in `sol_prev`, the
corresponding spatial coefficient blocks are transferred:
- if the grid is unchanged (stochastic refinement), the blocks are copied directly
- if the grid was refined (spatial refinement), the block of `sol_prev` is pointwise
  evaluated at the dof points of `sol` via `lazy_interpolate!`

Modes of `sol` without predecessor in `sol_prev` are left untouched (usually zero).
Requires matching finite element types (per unknown) in both `SGFEVector`s.

# Keywords
- `match_modes`: how to associate modes of `sol` with modes of `sol_prev`
  - `:append` (default): assume that both mode lists coincide up to the length of the
    shorter one, i.e. new modes were only appended at the end (as done by the AFEM
    loops in `scripts/`). A consistency check of the multi-indices (ignoring trailing
    zeros) is performed; if it fails, a warning is issued and `:multiindex` is used.
  - `:multiindex`: match modes by their multi-index (ignoring trailing zeros).

Returns the number of transferred stochastic modes.
"""
function interpolate_initial_guess!(sol::SGFEVector, sol_prev::SGFEVector; match_modes = :append)
    nunknowns = length(sol.FES_space)
    if nunknowns != length(sol_prev.FES_space)
        @warn "number of unknowns differs, skipping initial guess interpolation"
        return 0
    end
    for u in 1:nunknowns
        if eltype(sol.FES_space[u]) != eltype(sol_prev.FES_space[u])
            @warn "finite element types differ, skipping initial guess interpolation"
            return 0
        end
    end

    nmodes = num_multiindices(sol)
    nmodes_prev = num_multiindices(sol_prev)
    keys_new = [_multiindex_key(sol.TB.multi_indices[sol.active_modes[k]]) for k in 1:nmodes]
    keys_prev = [_multiindex_key(sol_prev.TB.multi_indices[sol_prev.active_modes[j]]) for j in 1:nmodes_prev]

    ## determine pairs (j,k) of matching modes: sol_prev mode j -> sol mode k
    pairs = Tuple{Int, Int}[]
    if match_modes == :append
        ncommon = min(nmodes, nmodes_prev)
        if all(keys_new[j] == keys_prev[j] for j in 1:ncommon)
            append!(pairs, [(j, j) for j in 1:ncommon])
        else
            @warn "mode lists do not coincide up to appended modes, falling back to multi-index matching"
            empty!(pairs)
        end
    end
    if isempty(pairs)
        if match_modes in (:append, :multiindex)
            prev_modes = Dict{Tuple, Int}()
            for j in 1:nmodes_prev
                prev_modes[keys_prev[j]] = j
            end
            for k in 1:nmodes
                j = get(prev_modes, keys_new[k], 0)
                j == 0 && continue
                push!(pairs, (j, k))
            end
        else
            error("unknown match_modes option $(match_modes)")
        end
    end

    nmoved = 0
    for (j, k) in pairs
        for u in 1:nunknowns
            tgt = sol[u, k]
            src = sol_prev[u, j]
            if tgt.FES.ndofs == src.FES.ndofs && tgt.FES.xgrid[Coordinates] == src.FES.xgrid[Coordinates]
                copyto!(view(tgt), view(src))
            else
                lazy_interpolate!(tgt, [src])
            end
        end
        nmoved += 1
    end
    return nmoved
end
