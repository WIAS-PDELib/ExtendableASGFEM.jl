function get_neighbours(OBT, multi_indices)
    M = length(multi_indices[1])
    nmodes = length(multi_indices)
    idx4mode = Dict{NTuple{M, Int}, Int}()
    for k in 1:nmodes
        idx4mode[(multi_indices[k]...,)] = k
    end
    mneighboursPLUS = zeros(Int, M, nmodes)
    mneighboursMINUS = zeros(Int, M, nmodes)
    for j in 1:nmodes
        mu = (multi_indices[j]...,)
        for m in 1:M
            k = get(idx4mode, ntuple(mm -> mu[mm] + (mm == m), M), 0)
            k > 0 && (mneighboursPLUS[m, j] = k)
            k = get(idx4mode, ntuple(mm -> mu[mm] - (mm == m), M), 0)
            k > 0 && (mneighboursMINUS[m, j] = k)
        end
    end
    return mneighboursPLUS, mneighboursMINUS
end

function get_next_tail(multi_indices, mode, startm, maxm)
    for m in startm:maxm
        mu1 = deepcopy(mode)
        if length(mu1) < m
            append!(mu1, zeros(Int, m - length(mu1)))
        end
        mu1[m] += 1
        if !(mu1 in multi_indices)
            return mu1, m
        end
    end
    return false
end

## Cache for prepare_extended_modes (see its docstring). The entry type is Any because the
## concrete types of TB_extended (ONB parameterized by polynomial family and degree) depend on
## the basis and cannot be declared for a heterogeneous Dict; type stability for the callers is
## recovered by type-annotating the fields when reading the entry back (see below).
const _extended_modes_cache_maxlen = 4
const _extended_modes_cache = Dict{Any, Any}()

## Cache for the interpolations of <e^-a, H_nu> of the log-transformed estimator: they depend
## only on (extended multi-indices, coefficient, interpolation order, grid), not on the current
## solution vector. The grid enters via objectid and a hash of its coordinates so that in-place
## mesh refinement (which the FEEvaluator's dofmap would not see) invalidates the entry.
## Same Any-storage + ::-annotated readback scheme as _extended_modes_cache.
const _lambda_interp_cache_maxlen = 2
const _lambda_interp_cache = Dict{Any, Any}()

"""
    prepare_extended_modes(sol::SGFEVector; tail_extension = [10, 2])

Shared estimator preparation: extend the multi-index set of `sol` by boundary/tail modes,
rebuild the associated tensorized basis and precompute the mode neighbour tables.

Since `TensorizedBasis` is immutable, the multi-index set plus `tail_extension` uniquely
determine the result; it is therefore cached (few slots) and recomputed only once per
basis extension, e.g. within an adaptive loop that calls the estimators repeatedly.

# Returns
`(multi_indices_extended, TB_extended, G, mneighboursPLUS, mneighboursMINUS, M_extended, nmodes_extended)`
"""
function prepare_extended_modes(sol::SGFEVector; tail_extension = [10, 2])
    TB = sol.TB
    ## cache key: (polynomial family, multi-index content, tail_extension).
    ## The polynomial family must be part of the key since TB_extended (and the triple
    ## products in its G matrix) depend on it; the multi-index content is safe to use as
    ## part of the key because TensorizedBasis is immutable, so it cannot go stale.
    cachekey = (OrthogonalPolynomialType(TB.ONB), Tuple.(TB.multi_indices)..., Tuple(tail_extension))
    if !haskey(_extended_modes_cache, cachekey)
        M = maxlength_multiindices(TB)
        OBT = OrthogonalPolynomialType(TB.ONB)

        ## extend multi_indices (note: add_boundary_modes mutates its argument)
        multi_indices_extended = add_boundary_modes(deepcopy(TB.multi_indices); tail_extension = tail_extension)
        M_extended = length(multi_indices_extended[1])
        maxorder = maximum(maximum(multi_indices_extended[j]) for j in 1:length(multi_indices_extended))
        TB_extended = TensorizedBasis(OBT, M + 1, maxorder, 2 * maxorder, 2 * maxorder, multi_indices = multi_indices_extended)
        nmodes_extended = length(multi_indices_extended)
        G = TB_extended.G

        ## prepare neighbours of modes
        mneighboursPLUS, mneighboursMINUS = get_neighbours(OBT, multi_indices_extended)

        length(_extended_modes_cache) >= _extended_modes_cache_maxlen && empty!(_extended_modes_cache)
        _extended_modes_cache[cachekey] = (multi_indices_extended, TB_extended, G, mneighboursPLUS, mneighboursMINUS, M_extended, nmodes_extended)
    end
    ## The cache stores Any, so read the entry back with ::-annotations: these are cheap
    ## runtime type checks that simultaneously give the compiler concrete types for the
    ## destructured return values (a bare return of the cache entry would type the callers'
    ## G, mneighbours..., etc. as Any and break the type stability of the estimator loops).
    cacheval = _extended_modes_cache[cachekey]
    return (
        cacheval[1]::Vector{Vector{Int}}, cacheval[2]::TensorizedBasis{Float64}, cacheval[3]::ExtendableSparseMatrix{Float64, Int64},
        cacheval[4]::Matrix{Int}, cacheval[5]::Matrix{Int}, cacheval[6]::Int, cacheval[7]::Int,
    )
end

"""
$(TYPEDSIGNATURES)

Compute the residual-based a posteriori error estimator for a stochastic Galerkin solution vector `sol`.

# Arguments
- `::Type{AbstractModelProblem}`: The model problem type for which the estimator is called. Used for dispatching to the appropriate estimator implementation.
- `sol::SGFEVector`: The current stochastic Galerkin solution vector.
- `C::AbstractStochasticCoefficient`: The stochastic coefficient (random field or parameterization).
- `rhs`: (Optional) Right-hand side function for the PDE (default: `nothing`).
- `bonus_quadorder`: (Optional) Additional quadrature order for integration (default: 1).
- `tail_extension`: (Optional) Two-element vector controlling the boundary mode extension (default: `[10, 2]`, see `add_boundary_modes`).
- `kwargs...`: Additional keyword arguments passed to the estimator.

# Returns
A tuple containing:
- `eta4modes::Vector{Float64}`: Total error estimator for each multi-index (stochastic mode), corresponding to the enriched set of multi-indices (with current active modes first).
- `eta4cell::Matrix{Float64}`: Error estimator for each cell in the spatial grid (for spatial refinement), for each multi-index.
- `multi_indices_extended`: The enriched set of multi-indices used in the computation (including boundary extensions).
- `ζ_data::Float64`: Estimated data truncation error (0.0 if the estimator does not provide one).

# Description
This function computes a residual-based a posteriori error estimator for the current SGFEM solution. It supports both spatial and stochastic adaptivity by providing error indicators for each cell and each stochastic mode. The estimator is tailored to the model problem and the stochastic coefficient, and can be extended to include additional boundary modes for improved reliability.

If no specialized estimator is available for the given model problem type, an error is raised and empty arrays are returned.
```
"""
function estimate(::Type{AbstractModelProblem}, sol::SGFEVector, C::AbstractStochasticCoefficient; kwargs...)
    @error "no error estimator for the model problem type available"
    return zeros(Float64, 0), zeros(Float64, 0, 0), [[0]], 0.0
end


function estimate(::Type{LogTransformedPoissonProblemPrimal}, sol::SGFEVector, C::AbstractStochasticCoefficient; rhs = nothing, bonus_quadorder = 1, tail_extension = [10, 2])
    if rhs === nothing
        error("estimate: a right-hand side `rhs` is needed for the residual computation")
    end

    FES = sol.FES_space[1]
    FEType = eltype(FES)
    xgrid = FES.xgrid
    ncells = num_cells(xgrid)
    EG = xgrid[UniqueCellGeometries][1]
    TB = sol.TB
    nmodes = TB.nmodes
    order = get_polynomialorder(FEType, EG)

    ## extend multi_indices, rebuild tensor basis and prepare mode neighbours
    multi_indices_extended, TB_extended, G, mneighboursPLUS, mneighboursMINUS, M_extended, nmodes_extended = prepare_extended_modes(sol; tail_extension = tail_extension)

    Mcoeff = maxm(C)
    ## coefficient modes beyond the coefficient's own expansion are identically zero
    Mcoup = min(Mcoeff, M_extended)

    ## prepare quadrature rule
    quadorder = 2 * (order - 1) + bonus_quadorder
    qf = QuadratureRule{Float64, EG}(quadorder)
    weights::Vector{Float64} = qf.w
    xref::Vector{Vector{Float64}} = qf.xref
    nweights::Int = length(weights)
    cellvolumes = xgrid[CellVolumes]

    ## prepare FE basis evaluator and dofmap
    FEBasis_∇ = FEEvaluator(FES, Gradient, qf)
    FEBasis_Δ = FEEvaluator(FES, Laplacian, qf)
    ∇vals = FEBasis_∇.cvals
    Δvals = FEBasis_Δ.cvals
    L2G = L2GTransformer(EG, xgrid, ON_CELLS)
    celldofs = FES[CellDofs]
    ndofs4cell::Int = get_ndofs(ON_CELLS, FEType, EG)
    ndofs = FES.ndofs
    coeffs = sol.entries

    ## interpolate <e^-a, H_nu> for all extended modes
    ## This depends only on the extended tensor basis, the coefficient, the quadrature order
    ## and the grid (not on the solution values), so the interpolations and the identity
    ## FE-basis evaluator are reused across estimator calls until one of these changes.
    ## As in _extended_modes_cache the entry is stored as Any and unpacked with ::-checks
    ## below to keep the hot loop in barrier() type stable.
    interpkey = (quadorder, multi_indices_extended, C, objectid(xgrid), hash(xgrid[Coordinates]))
    if !haskey(_lambda_interp_cache, interpkey)
        @info "interpolating <e^-a, H_nu> for all modes with quadorder = $quadorder"
        ## expansion of the coefficient, lambda_μ!(result, x, μ) evaluates <e^-a, H_μ>(x)
        expa_PCE!, lambda_μ! = expa_PCE_mop(TB_extended, C; factor = -1.0)
        FES_interp = FESpace{H1Pk{1, 2, quadorder}}(xgrid)
        FEBasis_id = FEEvaluator(FES_interp, ExtendableFEMBase.Identity, qf)
        expaf_interpolations = FEVector([FES_interp for j in 1:nmodes_extended])
        for j in 1:nmodes_extended
            interpolate!(expaf_interpolations[j], (result, qpinfo) -> lambda_μ!(result, qpinfo.x, j); quadorder = 2 * quadorder)
        end
        length(_lambda_interp_cache) >= _lambda_interp_cache_maxlen && empty!(_lambda_interp_cache)
        _lambda_interp_cache[interpkey] = (FEBasis_id, expaf_interpolations, FES_interp[CellDofs], FES_interp.ndofs, get_ndofs(ON_CELLS, H1Pk{1, 2, quadorder}, EG))
    end
    cacheval = _lambda_interp_cache[interpkey]
    FEBasis_id = cacheval[1]::ExtendableFEMBase.SingleFEEvaluator{Float64}
    idvals = FEBasis_id.cvals
    expaf_interpolations = cacheval[2]::FEVector{Float64}
    celldofs_interp = cacheval[3]::Union{VariableTargetAdjacency{Int32}, SerialVariableTargetAdjacency{Int32}, Array{Int32, 2}}
    offset_interp = cacheval[4]::Int
    ndofs4cell_interp = cacheval[5]::Int
    coeffs_interp = expaf_interpolations.entries

    ## precompute the mode coupling structure of sigma_h: for each mode j, the list of
    ## (coefficient index m, solution mode k, triple-product weight g) with g != 0 such
    ## that the sigma term of the mode-j residual reads sum_(m,k,g) g * grad(a_m) * grad(u_k)
    neighbours4mode = [Tuple{Int, Int, Float64}[] for j in 1:nmodes_extended]
    for j in 1:nmodes_extended
        for m in 1:Mcoup
            for k in (mneighboursPLUS[m, j], mneighboursMINUS[m, j])
                if 0 < k <= nmodes
                    g = G[(m - 1) * nmodes_extended + j, k]
                    if g != 0
                        push!(neighbours4mode[j], (m, k, g))
                    end
                end
            end
        end
    end
    ## modes whose residual gains the Δu_h term (requires second derivatives of the FE basis)
    add_lap4j = [order > 1 && j <= nmodes for j in 1:nmodes_extended]

    ## compute volume terms
    f4modes = zeros(Float64, nmodes_extended)
    eta4cell = zeros(Float64, ncells, nmodes_extended)
    eta4modes = zeros(Float64, nmodes_extended)
    ζ_data1 = 0.0
    ζ_data2 = 0.0
    function barrier(EG, L2G::L2GTransformer)
        gradam = zeros(Float64, 2)
        gradamvals = zeros(Float64, 2, Mcoup)
        am = zeros(Float64, 1)
        kmL2 = 0.0
        ftemp = zeros(Float64, 1)
        x = zeros(Float64, 2)

        for cell in 1:ncells
            update_basis!(FEBasis_∇, cell)
            update_basis!(FEBasis_id, cell)
            if order > 1
                update_basis!(FEBasis_Δ, cell)
            end
            update_trafo!(L2G, cell)

            for qp in 1:nweights
                eval_trafo!(x, L2G, xref[qp])

                ## evaluate all coefficient gradients once per quadrature point
                for m in 1:Mcoup
                    get_gradam!(gradam, x, m, C)
                    gradamvals[1, m] = gradam[1]
                    gradamvals[2, m] = gradam[2]
                end

                # compute ζ_data
                kmL2 = 0.0
                for m in 1:Mcoeff
                    am[1] = 0
                    get_am!(am, x, m, C)
                    kmL2 += am[1]^2
                end
                rhs(ftemp, x)
                ζ_data1 += ftemp[1]^2 * exp(2 * kmL2) * weights[qp] * cellvolumes[cell]
                fval = ftemp[1]

                ## h_T|| f_nu + \sigma_\nu ||
                for j in 1:nmodes_extended
                    ## evaluate <e^-a, H_nu>(x_qp) from its interpolation
                    # (instead of the much more expensive lambda_μ!(result, x, j))
                    lambda_val = 0.0
                    for d in 1:ndofs4cell_interp
                        dof = (j - 1) * offset_interp + celldofs_interp[d, cell]
                        lambda_val += idvals[1, d, qp] * coeffs_interp[dof]
                    end
                    ζ_data2 += lambda_val^2 * fval^2 * weights[qp] * cellvolumes[cell]

                    ## residual f_nu = <e^-a f, H_nu> + Δu_h (the latter only for P_k, k > 1)
                    res = lambda_val * fval
                    if add_lap4j[j]
                        for d in 1:ndofs4cell
                            dof = (j - 1) * ndofs + celldofs[d, cell]
                            res += coeffs[dof] * Δvals[1, d, qp]
                        end
                    end

                    ## sigma term via the precomputed coupling lists
                    sigma_val = 0.0
                    for (m, k, g) in neighbours4mode[j]
                        gax = gradamvals[1, m]
                        gay = gradamvals[2, m]
                        koff = (k - 1) * ndofs
                        dudg = 0.0
                        for d in 1:ndofs4cell
                            dudg += coeffs[koff + celldofs[d, cell]] * (gax * ∇vals[1, d, qp] + gay * ∇vals[2, d, qp])
                        end
                        sigma_val += g * dudg
                    end
                    eta4cell[cell, j] += (res + sigma_val)^2 * weights[qp]
                    f4modes[j] += res^2 * weights[qp] * cellvolumes[cell]
                end
            end

            ## boundary modes enjoy no Galerkin orthogonality and therefore have no additional h power
            for j in 1:nmodes_extended
                if j <= nmodes
                    eta4cell[cell, j] *= cellvolumes[cell]^2
                else
                    eta4cell[cell, j] *= cellvolumes[cell]
                end
            end
        end
        for j in 1:nmodes_extended
            eta4modes[j] = sqrt(sum(view(eta4cell, :, j)))
        end
        return
    end

    barrier(EG, L2G)

    for j in 1:nmodes_extended
        @info "mode = $(multi_indices_extended[j]) \t ||f_nu||^2 = $(f4modes[j]))"
    end

    ## compute normal jumps
    sol_j = FEVector(FES)
    cellfaces = xgrid[CellFaces]
    JumpIntegrator = ItemIntegratorDG((result, input, qp) -> (result[1] = dot(input, input)), [jump(grad(1))]; resultdim = 1, entities = ON_IFACES)
    for j in 1:nmodes
        sol_j.entries .= view(sol[j])
        jumps4face = ExtendableFEM.evaluate(JumpIntegrator, sol_j)
        jumps4face[1, :] .*= xgrid[FaceVolumes]
        for cell in 1:ncells, f in 1:size(cellfaces, 1)
            eta4cell[cell, j] += jumps4face[cellfaces[f, cell]]
        end
        ## merge jump contribution in quadratic norm
        eta4modes[j] = sqrt(eta4modes[j]^2 + sum(view(jumps4face, :)))
    end

    ζ_data = ζ_data1 - ζ_data2
    if ζ_data < 0
        @warn "negative data truncation error estimate = $(ζ_data) (ζ_data1 = $(ζ_data1), ζ_data2 = $(ζ_data2))
            This is likely due to insufficient quadrature order for the interpolation of <e^-a, H_nu>."
    else
        @info "estimated data truncation error = $(ζ_data) (ζ_data1 = $(ζ_data1), ζ_data2 = $(ζ_data2))"
    end

    return eta4modes, eta4cell, multi_indices_extended, ζ_data1 - ζ_data2 #, f4modes
end


function estimate(::Type{PoissonProblemPrimal}, sol::SGFEVector, C::AbstractStochasticCoefficient; rhs = nothing, bonus_quadorder = 1, tail_extension = [10, 2])
    if rhs === nothing
        error("estimate: a right-hand side `rhs` is needed for the residual computation")
    end

    FES = sol.FES_space[1]
    FEType = eltype(FES)
    xgrid = FES.xgrid
    ncells = num_cells(xgrid)
    EG = xgrid[UniqueCellGeometries][1]
    TB = sol.TB
    nmodes = TB.nmodes
    order = get_polynomialorder(FEType, EG)

    ## extend multi_indices, rebuild tensor basis and prepare mode neighbours
    multi_indices_extended, TB_extended, G, mneighboursPLUS, mneighboursMINUS, M_extended, nmodes_extended = prepare_extended_modes(sol; tail_extension = tail_extension)

    ## coefficient modes beyond the coefficient's own expansion are identically zero
    Mcoup = min(maxm(C), M_extended)

    ## prepare quadrature rule
    quadorder = 2 * (order - 1) + bonus_quadorder
    qf = QuadratureRule{Float64, EG}(quadorder)
    weights::Vector{Float64} = qf.w
    xref::Vector{Vector{Float64}} = qf.xref
    nweights::Int = length(weights)
    cellvolumes = xgrid[CellVolumes]

    ## prepare FE basis evaluator and dofmap
    FEBasis_Δ = FEEvaluator(FES, Laplacian, qf)
    Δvals = FEBasis_Δ.cvals
    L2G = L2GTransformer(EG, xgrid, ON_CELLS)
    celldofs = FES[CellDofs]
    ndofs4cell::Int = get_ndofs(ON_CELLS, FEType, EG)
    ndofs = FES.ndofs
    coeffs = sol.entries

    ## precompute the mode-wise coupling structure of div(a grad u_h): for each mode j,
    ## the list of (coefficient index m, solution mode k, triple-product weight g) such
    ## that the mode-j residual reads f_nu + sum_m a_m * sum_(m,k,g) g * Laplace(u_k)
    neighbours4mode = [Tuple{Int, Int, Float64}[] for j in 1:nmodes_extended]
    for j in 1:nmodes_extended
        if j <= nmodes
            push!(neighbours4mode[j], (0, j, 1.0))
        end
        for m in 1:Mcoup
            for k in (mneighboursPLUS[m, j], mneighboursMINUS[m, j])
                if 0 < k <= nmodes
                    g = G[(m - 1) * nmodes_extended + j, k]
                    if g != 0
                        push!(neighbours4mode[j], (m, k, g))
                    end
                end
            end
        end
    end

    ## compute volume terms
    eta4cell = zeros(Float64, ncells, nmodes_extended)
    eta4modes = zeros(Float64, nmodes_extended)
    ## local wrapper function: keeps the hot loop in a type-stable local scope
    function barrier(EG, L2G::L2GTransformer)
        amtmp = zeros(Float64, 1)
        amvals = zeros(Float64, Mcoup + 1)
        ftemp = zeros(Float64, 1)
        x = zeros(Float64, 2)
        lap4mode = zeros(Float64, nmodes, nweights)

        for cell in 1:ncells
            if order > 1
                update_basis!(FEBasis_Δ, cell)
            end
            update_trafo!(L2G, cell)
            if order > 1
                ## Laplace(u_k)(x_qp) depends only on (cell, k, qp): precompute once per cell
                for k in 1:nmodes
                    koff = (k - 1) * ndofs
                    for qp in 1:nweights
                        lap = 0.0
                        for d in 1:ndofs4cell
                            lap += coeffs[koff + celldofs[d, cell]] * Δvals[1, d, qp]
                        end
                        lap4mode[k, qp] = lap
                    end
                end
                for qp in 1:nweights
                    eval_trafo!(x, L2G, xref[qp])

                    ## evaluate all coefficient functions a_m(x) once per quadrature point
                    for m in 0:Mcoup
                        get_am!(amtmp, x, m, C)
                        amvals[m + 1] = amtmp[1]
                    end
                    rhs(ftemp, x)

                    ## h_T|| f_nu + div \sigma_\nu ||
                    for j in 1:nmodes_extended
                        res = j == 1 ? ftemp[1] : 0.0
                        for (m, k, g) in neighbours4mode[j]
                            res += amvals[m + 1] * g * lap4mode[k, qp]
                        end
                        eta4cell[cell, j] += res^2 * weights[qp]
                    end
                end
            else
                ## P1: the element-wise Laplacian vanishes, only the deterministic rhs contributes
                for qp in 1:nweights
                    eval_trafo!(x, L2G, xref[qp])
                    rhs(ftemp, x)
                    eta4cell[cell, 1] += ftemp[1]^2 * weights[qp]
                end
            end

            ## mesh-size scaling (h_T ~ sqrt(|T|) in 2D): the reference quadrature weights only
            ## provide the cell mean of the squared residual, and interior modes gain two
            ## additional powers of h_T from Galerkin orthogonality;
            ## boundary modes enjoy no Galerkin orthogonality and therefore have no additional h power
            for j in 1:nmodes_extended
                if j <= nmodes
                    eta4cell[cell, j] *= cellvolumes[cell]^3
                else
                    eta4cell[cell, j] *= cellvolumes[cell]
                end
            end
        end
        for j in 1:nmodes_extended
            eta4modes[j] = sqrt(sum(view(eta4cell, :, j)))
        end
        return
    end

    barrier(EG, L2G)

    ## compute normal jumps of a grad u_h.
    ## Linearity trick: instead of assembling the flux sum_m a_m * sum_k g*u_k of each mode j
    ## (one FaceInterpolator sweep per mode and coefficient), evaluate a_m * jump(grad u_k)
    ## at the face quadrature points only for the distinct (m, k) pairs that actually occur
    ## in the coupling lists neighbours4mode = (m, k, g), and form the mode-j flux values as
    ## linear combination of these stored pair values (the (0, j, 1.0) entries cover a_0*grad u_j)
    sol_j = FEVector(FES)
    cellfaces = xgrid[CellFaces]
    jumps4face = zeros(Float64, 1, size(xgrid[FaceNodes], 2))

    m_pointer = [0]
    JumpEvaluator = FaceInterpolator((result, input, qpinfo) -> (get_am!(result, qpinfo.x, m_pointer[1], C); result .= result[1] * input), [jump(grad(1))], quadorder = quadorder)
    ExtendableFEM.build_assembler!(JumpEvaluator, [sol_j[1]])
    sol_jumps = deepcopy(JumpEvaluator.value)
    nquadvals = length(JumpEvaluator.value.entries)

    ## enumerate the needed (coefficient, solution mode) pairs
    idx4pair = zeros(Int, Mcoup + 1, nmodes)
    pairs = Tuple{Int, Int}[]
    for j in 1:nmodes_extended
        for (m, k, g) in neighbours4mode[j]
            if idx4pair[m + 1, k] == 0
                push!(pairs, (m, k))
                idx4pair[m + 1, k] = length(pairs)
            end
        end
    end
    pairvalues = zeros(Float64, length(pairs), nquadvals)
    for (i, (m, k)) in enumerate(pairs)
        sol_j.entries .= view(sol[k])
        m_pointer[1] = m
        ExtendableFEM.evaluate!(JumpEvaluator, sol_j)
        pairvalues[i, :] .= JumpEvaluator.value.entries
    end

    JumpIntegrator = L2NormIntegrator([id(1)]; entities = ON_FACES)
    bfaces = xgrid[BFaceFaces]
    sigmavalues = zeros(Float64, nquadvals)

    for j in 1:nmodes_extended
        ## pointwise flux values sum_(m,k,g) g * a_m * jump(grad u_k) at the face quadrature points
        fill!(sigmavalues, 0)
        for (m, k, g) in neighbours4mode[j]
            prow = view(pairvalues, idx4pair[m + 1, k], :)
            for ii in 1:nquadvals
                sigmavalues[ii] += g * prow[ii]
            end
        end
        sol_jumps.entries .= sigmavalues
        # @error "Mode:", multi_indices_extended[j]
        # @error "Normal:", extrema(sol_jumps.entries)
        jumps4face .= sum(ExtendableFEM.evaluate(JumpIntegrator, sol_jumps), dims = 1)
        jumps4face[1, bfaces] .= 0

        ## jump scaling with edge size h_F ~ |F|: interior modes get one power of h_F,
        ## tail modes (no Galerkin orthogonality) an inverse power
        if j <= nmodes
            jumps4face[1, :] .*= xgrid[FaceVolumes]
        else
            jumps4face[1, :] ./= xgrid[FaceVolumes]
        end

        for cell in 1:ncells
            for f in 1:size(cellfaces, 1)
                eta4cell[cell, j] += jumps4face[cellfaces[f, cell]]
            end
        end
        eta4modes[j] = sqrt(eta4modes[j]^2 + sum(view(jumps4face, :)))
    end

    ## no data truncation estimate for this problem (4th slot keeps the return signature uniform)
    return eta4modes, eta4cell, multi_indices_extended, 0.0
end

function estimate(::Type{StokesProblemPrimal}, sol::SGFEVector, C::AbstractStochasticCoefficient; rhs = nothing, rhs_curl = nothing, bonus_quadorder = 1, tail_extension = 5)
    if isnothing(rhs)
        error("estimate: a right-hand side `rhs` is needed for the residual computation")
    end
    if isnothing(rhs_curl)
        error("estimate: the curl `rhs_curl` of the `rhs` is needed for the residual computation")
    end

    FES = sol.FES_space[1]
    FEType = eltype(FES)
    xgrid = FES.xgrid
    ncells = num_cells(xgrid)
    EG = xgrid[UniqueCellGeometries][1]
    TB = sol.TB
    nmodes = TB.nmodes

    order = 1 # ToDo: Insert get_polynomialorder(FEType, EG) again

    ## extend multi_indices, rebuild tensor basis and prepare mode neighbours
    multi_indices_extended, TB_extended, G, mneighboursPLUS, mneighboursMINUS, M_extended, nmodes_extended = prepare_extended_modes(sol; tail_extension = tail_extension)

    ## coefficient modes beyond the coefficient's own expansion are identically zero
    Mcoup = min(maxm(C), M_extended)

    ## prepare quadrature rule
    quadorder = 2 * (order - 1) + bonus_quadorder
    qf = QuadratureRule{Float64, EG}(quadorder)
    weights::Vector{Float64} = qf.w
    xref::Vector{Vector{Float64}} = qf.xref
    nweights::Int = length(weights)
    cellvolumes = xgrid[CellVolumes]

    ## prepare FE basis evaluator and dofmap
    FEBasis_Δ = FEEvaluator(FES, Laplacian, qf)
    Δvals = FEBasis_Δ.cvals
    L2G = L2GTransformer(EG, xgrid, ON_CELLS)
    celldofs = FES[CellDofs]
    ndofs4cell::Int = get_ndofs(ON_CELLS, FEType, EG)
    ndofs = FES.ndofs
    coeffs = sol.entries

    ## precompute the mode-wise coupling structure of div(a grad u_h): for each mode j,
    ## the list of (coefficient index m, solution mode k, triple-product weight g) such
    ## that the mode-j residual reads f_nu + sum_m a_m * sum_(m,k,g) g * Laplace(u_k)
    neighbours4mode = [Tuple{Int, Int, Float64}[] for j in 1:nmodes_extended]
    for j in 1:nmodes_extended
        if j <= nmodes
            push!(neighbours4mode[j], (0, j, 1.0))
        end
        for m in 1:Mcoup
            for k in (mneighboursPLUS[m, j], mneighboursMINUS[m, j])
                if 0 < k <= nmodes
                    g = G[(m - 1) * nmodes_extended + j, k]
                    if g != 0
                        push!(neighbours4mode[j], (m, k, g))
                    end
                end
            end
        end
    end

    # ToDo:
    #     η_curl(σ_μ)   = ||h_T^2 curl(f_μ + div σ_μ)||_{L^2(D)}
    #     η_jump(σ_μ)   = ||sqrt(h_E) [σ_μn_E]||_{L^2(E^∘)}
    #     η_jump,2(σ_μ) = ||h_E^(3/2) [(f_μ + div σ_μ) ⋅ τ_E]||_{L^2(E^∘)}
    #     η_cons,1(σ_μ) = ||div σ_μ ∘ (1 - Π)||_{V'_{h, 0}}

    ## compute volume terms
    eta4cell = zeros(Float64, ncells, nmodes_extended)
    eta4modes = zeros(Float64, nmodes_extended)
    ## local wrapper function: keeps the hot loop in a type-stable local scope
    function barrier(EG, L2G::L2GTransformer)
        amtmp = zeros(Float64, 1)
        amvals = zeros(Float64, Mcoup + 1)
        ftemp = zeros(Float64, 2)
        x = zeros(Float64, 2)
        lap4mode = zeros(Float64, nmodes, nweights)

        # update_basis!(FEBasis_NormalFlux, cell)
        # update_basis!(FEBasis_TangentFlux, cell)

        for cell in 1:ncells
            if order > 1
                update_basis!(FEBasis_Δ, cell)
            end
            update_trafo!(L2G, cell)
            if order > 1
                error("The a posteriori error estimator for the StokesProblemPrimal is only implemented for FE pairs of order = 1")
            else
                #     η_cons,1(σ_μ) = ||div σ_μ ∘ (1 - Π)||_{V'_{h, 0}} <- result .= qpinfo.volume .* sum(Δu .^ 2)
                #                   = ||Δu||
                #     Volume term -> Implement in this loop

                for qp in 1:nweights
                    #     η_curl(σ_μ)   = ||h_T^2 curl(f_μ + div σ_μ)||_{L^2(D)}
                    eval_trafo!(x, L2G, xref[qp])
                    ## If velocity space is order 1: ||h_T^2 curl(f_μ + div σ_μ)|| = ||h_T^2 curl f_μ||
                    rhs_curl(ftemp, x)
                    eta4cell[cell, 1] += ftemp[1]^2 * weights[qp]
                    ## Approximation for η_jump,2
                end
            end

            ## mesh-size scaling (h_T ~ sqrt(|T|) in 2D): the reference quadrature weights only
            ## provide the cell mean of the squared residual, and interior modes gain two
            ## additional powers of h_T from Galerkin orthogonality;
            ## boundary modes enjoy no Galerkin orthogonality and therefore have no additional h power
            for j in 1:nmodes_extended
                if j <= nmodes
                    eta4cell[cell, j] *= cellvolumes[cell]^3
                else
                    eta4cell[cell, j] *= cellvolumes[cell]
                end
            end
        end
        for j in 1:nmodes_extended
            eta4modes[j] = sqrt(sum(view(eta4cell, :, j)))
        end
        return
    end

    barrier(EG, L2G)

    ## compute normal jumps of a grad u_h.
    ## Linearity trick: instead of assembling the flux sum_m a_m * sum_k g*u_k of each mode j
    ## (one FaceInterpolator sweep per mode and coefficient), evaluate a_m * jump(grad u_k)
    ## at the face quadrature points only for the distinct (m, k) pairs that actually occur
    ## in the coupling lists neighbours4mode = (m, k, g), and form the mode-j flux values as
    ## linear combination of these stored pair values (the (0, j, 1.0) entries cover a_0*grad u_j)

    ## New compute tangential jumps of grad u_h?

    sol_j = FEVector(FES)
    cellfaces = xgrid[CellFaces]
    jumps4face = zeros(Float64, 1, size(xgrid[FaceNodes], 2))
    jumps4face_laplace = zeros(Float64, 1, size(xgrid[FaceNodes], 2))

    m_pointer = [0]

    #     η_jump(σ_μ)   = ||sqrt(h_E) [σ_μn_E]||_{L^2(E^∘)}
    #     η_jump,2(σ_μ) = ||h_E^(3/2) [(f_μ + div σ_μ) ⋅ τ_E]||_{L^2(E^∘)}
    JumpEvaluator = FaceInterpolator(
        (result, input, qpinfo) -> (get_am!(result, qpinfo.x, m_pointer[1], C); result .= result[1] .* [input[1] input[2]; input[3] input[4]] * qpinfo.normal), [jump(grad(1))];
        resultdim = 2, # quadorder = 0
    )
    ExtendableFEM.build_assembler!(JumpEvaluator, [sol_j[1]])
    sol_jumps = deepcopy(JumpEvaluator.value)

    LaplaceJumpEvaluator = FaceInterpolator(
        (result, input, qpinfo) -> (
            get_am!(result, qpinfo.x, m_pointer[1], C);
            result .= result[1] * input
        ), [jump(Δ(1))] #quadorder = 0
    )
    ExtendableFEM.build_assembler!(LaplaceJumpEvaluator, [sol_j[1]])
    sol_laplace_jumps = deepcopy(LaplaceJumpEvaluator.value)

    JumpIntegrator = L2NormIntegrator([id(1)]; entities = ON_FACES)
    bfaces = xgrid[BFaceFaces]

    TangentJumpIntegrator = ItemIntegrator(
        (result, input, qpinfo) -> (result[1] = (input[1] * qpinfo.normal[2] - input[2] * qpinfo.normal[1])^2), [id(1)];
        entities = ON_FACES, resultdim = 1, quadorder = 0
    )

    @time for j in 1:nmodes_extended
        fill!(sol_jumps[1], 0)
        fill!(sol_laplace_jumps[1], 0)
        for m in 0:M_extended
            fill!(sol_j[1], 0)
            if m == 0 && j <= nmodes
                sol_j.entries .+= view(sol[j])
            elseif m > 0
                if mneighboursPLUS[m, j] > 0 && mneighboursPLUS[m, j] <= nmodes
                    sol_j.entries .+= view(sol[mneighboursPLUS[m, j]]) * G[(m - 1) * nmodes_extended + j, mneighboursPLUS[m, j]]
                end
                if mneighboursMINUS[m, j] > 0 && mneighboursMINUS[m, j] <= nmodes
                    sol_j.entries .+= view(sol[mneighboursMINUS[m, j]]) * G[(m - 1) * nmodes_extended + j, mneighboursMINUS[m, j]]
                end
            end
            m_pointer[1] = m
            ExtendableFEM.evaluate!(JumpEvaluator, sol_j)
            sol_jumps.entries .+= JumpEvaluator.value.entries
            ExtendableFEM.evaluate!(LaplaceJumpEvaluator, sol_j)
            sol_laplace_jumps.entries .+= LaplaceJumpEvaluator.value.entries
        end
        @error "Mode:", multi_indices_extended[j]
        @error "Normal:", extrema(sol_jumps.entries)
        @error "Laplacian:", extrema(sol_laplace_jumps.entries)

        jumps4face .= sum(ExtendableFEM.evaluate(JumpIntegrator, sol_jumps), dims = 1)
        jumps4face[1, bfaces] .= 0

        jumps4face_laplace .= sum(ExtendableFEM.evaluate(TangentJumpIntegrator, sol_laplace_jumps), dims = 1)
        jumps4face_laplace[1, bfaces] .= 0

        if j <= nmodes
            jumps4face[1, :] .*= xgrid[FaceVolumes]
            jumps4face_laplace[1, :] .*= xgrid[FaceVolumes] .^ 3
        else
            jumps4face[1, :] ./= xgrid[FaceVolumes]
            jumps4face_laplace[1, :] .*= xgrid[FaceVolumes]
        end

        for cell in 1:ncells
            for f in 1:3
                eta4cell[cell, j] += jumps4face[cellfaces[f, cell]]
                eta4cell[cell, j] += jumps4face_laplace[cellfaces[f, cell]]
            end
        end
        eta4modes[j] = sqrt(eta4modes[j]^2 + sum(view(jumps4face, :)) + sum(view(jumps4face_laplace, :)))
    end

    ## no data truncation estimate for this problem (4th slot keeps the return signature uniform)
    return eta4modes, eta4cell, multi_indices_extended, 0.0
end


# function estimate(::Type{LogTransformedPoissonProblemDual}, sol::SGFEVector, C::AbstractStochasticCoefficient; problem = LogTransformedPoissonProblemPrimal, rhs = nothing, bonus_quadorder = 0)

#     ## read FE spaces
#     FES_p = sol.FES_space[1]
#     FES_u = sol.FES_space[2]
#     FEType_p = eltype(FES_p)
#     FEType_u = eltype(FES_u)

#     ## read grid
#     xgrid = FES_u.xgrid
#     ncells = num_cells(xgrid)
#     EG = xgrid[UniqueCellGeometries][1]
#     order = get_polynomialorder(FEType_u, EG)

#     ## read tensor basis
#     TB = sol.TB
#     M = maxlength_multiindices(TB)
#     nmodes = TB.nmodes
#     multi_indices = TB.multi_indices
#     OBT = OrthogonalPolynomialType(TB.ONB)

#     ## interpolate u into H1 space
#     FEType_uint = H1Pk{1,2,order+1}
#     FES_uint = FESpace{FEType_uint}(xgrid)
#     sol_uj = FEVector(FES_u)
#     sol_uint = FEVector([FES_uint for j = 1 : nmodes])
#     for j = 1 : nmodes
#         sol_uj.entries .= view(sol[nmodes+j])
#         lazy_interpolate!(sol_uint[j], sol_uj, [id(1)])
#     end

#     ## extend multi_indices
#     multi_indices_extended = add_boundary_modes(deepcopy(multi_indices))
#     M_extended = length(multi_indices_extended[1])
#     maxorder = maximum(maximum(multi_indices_extended[j]) for j = 1 : length(multi_indices_extended))
#     TB_extended = TensorizedBasis(OBT, M+1, maxorder, 2*maxorder, 2*maxorder, multi_indices = multi_indices_extended)
#     nmodes_extended = length(multi_indices_extended)
#     G = TB_extended.G

#     ## prepare neighbours of modes
#     mneighboursPLUS, mneighboursMINUS = get_neighbours(OBT, multi_indices_extended)

#     ## prepare quadrature rule
#     qf = QuadratureRule{Float64, EG}(2 * (order - 1) + bonus_quadorder)
#     weights::Vector{Float64} = qf.w
#     xref::Vector{Vector{Float64}} = qf.xref
#     nweights::Int = length(weights)
#     cellvolumes = xgrid[CellVolumes]

#     ## prepare FE basis evaluator and dofmap
#     FEBasis_∇uint = FEEvaluator(FES_uint, Gradient, qf)
#     FEBasis_uint = FEEvaluator(FES_uint, ExtendableFEMBase.Identity, qf)
#     FEBasis_p = FEEvaluator(FES_p, ExtendableFEMBase.Identity, qf)
#     FEBasis_divp = FEEvaluator(FES_p, Divergence, qf)
#     ∇uint_vals = FEBasis_∇uint.cvals
#     uint_vals = FEBasis_uint.cvals
#     p_vals = FEBasis_p.cvals
#     divp_vals = FEBasis_divp.cvals
#     L2G = L2GTransformer(EG, xgrid, ON_CELLS)
#     celldofs_p = FES_p[CellDofs]
#     celldofs_u = FES_u[CellDofs]
#     celldofs_uint = FES_uint[CellDofs]
#     ndofs4cell_p::Int = get_ndofs(ON_CELLS, FEType_p, EG)
#     ndofs4cell_uint::Int = get_ndofs(ON_CELLS, FEType_uint, EG)
#     ndofs_p = FES_p.ndofs
#     ndofs_uint = FES_uint.ndofs
#     coeffs = sol.entries
#     coeffs_uint = sol_uint.entries

#     ## prepare expansion of coefficient
#     expa_PCE!, lambda_μ! = expa_PCE_mop(TB_extended, C)


#     @info "...starting error estimation"

#     ## compute deterministic best-approximation matrix
#     E = FEMatrix(FES_uint)
#     L = FEMatrix(FES_p,FES_uint)
#     bη = FEVector(FES_uint)
#     vμ = FEVector(FES_uint)
#     vη = FEVector([FES_uint for j = 1 : nmodes])
#     assemble!(E, BilinearOperator([grad(1)]))

#     ## get boundary dofs and penalize them
#     bfacedofs = FES_uint[BFaceDofs]
#     nbfaces = num_sources(bfacedofs)
#     bdofs = []
#     for bface = 1 : nbfaces
#         append!(bdofs, view(bfacedofs,:,bface))
#     end
#     unique!(bdofs)
#     for dof in bdofs
#         E.entries[dof,dof] = 1e60
#     end

#     ## get LU decomposition
#     flush!(E.entries)
#     LUfacE = lu(E.entries.cscmatrix)

#     ## get coefficients for (e^-a,H_μ), evaluate at x with lambda_μ!(result,x,μ)
#     expa_PCE!, lambda_μ! = expa_PCE_mop(TB_extended, C)


#     ## best--approximate each component of q_h
#     function λ_kernel!(result,input,qpinfo)
#         lambda_μ!(result,qpinfo.x,qpinfo.params[1],-1.0)
#         result[2] = result[1] * input[2]
#         result[1] = result[1] * input[1]
#         return nothing
#     end
#     for μ = 1 : nmodes
#         # update bη (bnd data and right-hand side)
#         fill!(bη.entries,0)
#         for λ = 1 : nmodes_extended
#             fill!(L.entries.cscmatrix.nzval,0)
#             assemble!(L, BilinearOperator(λ_kernel!, [id(1)], [grad(1)]; params = [λ]))
#             for γ = 1 : nmodes
#                 ## need triple product here !!!!
#                 #g = G[(λ-1)*nmodes+γ,μ]
#                 g = triple_product(TB_extended, λ, γ, μ; normalize = true)
#                 if abs(g) > 1e-13
#                     addblock_matmul!(bη[1], L[1,1], sol[γ]; factor = -g, transposed = true)
#                 end
#             end
#         end

#         # solve
#         vμ.entries .= LUfacE \ bη.entries

#         # copy into long vector
#         addblock!(vη[μ], vμ.entries)

#         println(stdout, unicode_scalarplot(vη[μ]; title = "v for $(multi_indices_extended[μ])"))

#     end

#     ## TODO evaluate estimator η_μ := || <e^{a} p_h, H_μ> - ∇ v_μ || for each μ

#     ## compute volume terms
#     eta4cell = zeros(Float64, ncells, nmodes_extended)
#     eta4modes = zeros(Float64, nmodes_extended)
#     function barrier2(EG, L2G::L2GTransformer)
#         lambda_temp = zeros(Float64, 1)
#         gradam = zeros(Float64, 2)
#         ftemp = zeros(Float64, 1)
#         sigmatemp = zeros(Float64, 2)
#         x = zeros(Float64, 2)

#         for cell = 1 : ncells
#             update_basis!(FEBasis_∇uint, cell)
#             update_basis!(FEBasis_p, cell)
#             update_trafo!(L2G, cell)
#             for j = 1 : nmodes_extended

#                     eta4cell[cell,j] += (sigmatemp[1]^2 + sigmatemp[2]^2) * weights[qp]
#                 end

#                 eta4cell[cell,j] *= cellvolumes[cell]
#             end
#         end
#         for j = 1 : nmodes_extended
#             eta4modes[j] = sqrt(sum(view(eta4cell,:,j)))
#         end
#     end

#     barrier2(EG, L2G)


#     return eta4modes, eta4cell, multi_indices_extended
# end


################################
### EQUILIBRATION ESTIMATORS ###
################################


## kernel for equilibration error estimator
function eqestimator_kernel!(result, input, qpinfo)
    σ_h, divσ_h, ∇u_h = view(input, 1:2), input[3], view(input, 4:5)
    result[1] = norm(σ_h .- ∇u_h)^2 + divσ_h^2
    return nothing
end

## this function computes the local equilibrated fluxes
## by solving local problems on (disjunct groups of) node patches
function estimate_equilibration(::Type{PoissonProblemPrimal}, sol::SGFEVector, C::AbstractStochasticCoefficient; FETypeDual = nothing, rhs = nothing, tail_extension = [10, 2], bonus_quadorder = 0)
    if rhs === nothing
        error("estimate_equilibration: a right-hand side `rhs` is needed for the residual computation")
    end

    ## needed grid stuff
    FES = sol.FES_space[1]
    ndofs = FES.ndofs
    xgrid = FES.xgrid
    xCellNodes::Array{Int32, 2} = xgrid[CellNodes]
    xCellVolumes::Array{Float64, 1} = xgrid[CellVolumes]
    xNodeCells::Adjacency{Int32} = atranspose(xCellNodes)
    nnodes::Int = num_sources(xNodeCells)

    order = get_polynomialorder(eltype(FES), Triangle2D)
    if FETypeDual === nothing
        FETypeDual = HDIVRTk{2, order}
    end

    ## get node patch groups that can be solved in parallel
    group4node = xgrid[NodePatchGroups]

    ## init equilibration space (and Lagrange multiplier space)
    FESDual = FESpace{FETypeDual}(xgrid)
    celldofs_dual::Union{VariableTargetAdjacency{Int32}, SerialVariableTargetAdjacency{Int32}, Array{Int32, 2}} = FESDual[CellDofs]
    celldofs::Union{VariableTargetAdjacency{Int32}, SerialVariableTargetAdjacency{Int32}, Array{Int32, 2}} = FES[CellDofs]

    ## extract tensor basis
    TB = sol.TB
    nmodes = TB.nmodes

    ## extend multi_indices, rebuild tensor basis and prepare mode neighbours
    multi_indices_extended, TB_extended, G, mneighboursPLUS, mneighboursMINUS, M_extended, nmodes_extended = prepare_extended_modes(sol; tail_extension = tail_extension)

    ## coefficient modes beyond the coefficient's own expansion are identically zero
    Mcoup = min(maxm(C), M_extended)

    ## append block in solution vector for equilibrated fluxes
    sol_eq = SGFEVector(FESDual, TB_extended; active_modes = 1:length(multi_indices_extended))

    ## partition of unity and their gradients = P1 basis functions
    POUFES = FESpace{H1P1{1}}(xgrid)
    POUqf = QuadratureRule{Float64, Triangle2D}(0)

    ## quadrature formulas
    qf = QuadratureRule{Float64, Triangle2D}(2 * get_polynomialorder(FETypeDual, Triangle2D))
    xref::Vector{Vector{Float64}} = qf.xref
    weights::Array{Float64, 1} = qf.w

    ## some constants
    offset::Int = 0
    ncells::Int = num_cells(xgrid)
    div_penalty::Float64 = 1.0e5      # divergence constraint is realized by penalisation
    bnd_penalty::Float64 = 1.0e60     # penalty for non-involved dofs of a group
    maxdofs::Int = max_num_targets_per_source(celldofs_dual)
    maxdofs_uh::Int = max_num_targets_per_source(celldofs)
    coeffs = sol.entries

    ## redistribute groups for more equilibrated thread load (first groups are larger)
    maxgroups = maximum(group4node)
    groups = Array{Int, 1}(1:maxgroups)
    for j::Int in 1:floor(maxgroups / 2)
        a = groups[j]
        groups[j] = groups[2 * j]
        groups[2 * j] = a
    end
    X = Array{Array{Float64, 1}, 2}(undef, maxgroups, nmodes_extended)

    f4modes = zeros(Float64, nmodes_extended)
    function solve_patchgroup!(group, mode)
        ## temporary variables
        graduh::Array{Float64, 1} = zeros(Float64, 2)
        coeffs_uh::Array{Float64, 1} = zeros(Float64, maxdofs_uh)
        f_temp::Array{Float64, 1} = zeros(Float64, 1)
        am::Array{Float64, 1} = zeros(Float64, 1)
        x::Array{Float64, 1} = zeros(Float64, 2)
        Alocal = zeros(Float64, maxdofs, maxdofs)
        blocal = zeros(Float64, maxdofs)

        ## init system
        A = ExtendableSparseMatrix{Float64, Int64}(FESDual.ndofs, FESDual.ndofs)
        b = zeros(Float64, FESDual.ndofs)

        ## init FEBasiEvaluators
        FEE_∇φ = FEEvaluator(POUFES, Gradient, POUqf)
        FEE_xref = FEEvaluator(POUFES, ExtendableFEMBase.Identity, qf)
        FEE_∇u = FEEvaluator(FES, Gradient, qf)
        FEE_div = FEEvaluator(FESDual, Divergence, qf)
        FEE_id = FEEvaluator(FESDual, ExtendableFEMBase.Identity, qf)
        L2G = L2GTransformer(Triangle2D, xgrid, ON_CELLS)
        idvals = FEE_id.cvals
        divvals = FEE_div.cvals
        xref_vals = FEE_xref.cvals
        ∇φvals = FEE_∇φ.cvals

        ## find dofs at boundary of current node patches
        ## and in interior of cells outside of current node patch group
        is_noninvolveddof = zeros(Bool, FESDual.ndofs)
        outside_cell::Bool = false
        for cell in 1:ncells
            outside_cell = true
            for k in 1:3
                if group4node[xCellNodes[k, cell]] == group
                    outside_cell = false
                    break
                end
            end
            if (outside_cell) # mark interior dofs of outside cell
                for j in 1:maxdofs
                    is_noninvolveddof[celldofs_dual[j, cell]] = true
                end
            end
        end

        dofplus::Int = 0
        dofminus::Int = 0
        coeffMINUS::Float64 = 0
        coeffPLUS::Float64 = 0
        for node in 1:nnodes
            if group4node[node] == group
                for c in 1:num_targets(xNodeCells, node)
                    cell::Int = xNodeCells[c, node]

                    ## find local node number of global node z
                    ## and evaluate (constant) gradient of nodal basis function phi_z
                    localnode = 1
                    while xCellNodes[localnode, cell] != node
                        localnode += 1
                    end
                    FEE_∇φ.citem[] = cell
                    update_basis!(FEE_∇φ)

                    ## update other FE evaluators
                    FEE_∇u.citem[] = cell
                    FEE_div.citem[] = cell
                    FEE_id.citem[] = cell
                    update_basis!(FEE_∇u)
                    update_basis!(FEE_div)
                    update_basis!(FEE_id)
                    update_trafo!(L2G, cell)


                    ## assembly on this cell
                    for i in eachindex(weights)
                        eval_trafo!(x, L2G, xref[i])

                        weight = weights[i] * xCellVolumes[cell]

                        ## read coefficients for discrete flux
                        fill!(coeffs_uh, 0)
                        if mode <= nmodes
                            get_am!(am, x, 0, C)
                            for d in 1:maxdofs_uh
                                dof = (mode - 1) * ndofs + celldofs[d, cell]
                                coeffs_uh[d] += coeffs[dof] * am[1]
                            end
                        end
                        for m in 1:Mcoup
                            get_am!(am, x, m, C)
                            for d in 1:maxdofs_uh
                                dofplus = mneighboursPLUS[m, mode] <= nmodes ? (mneighboursPLUS[m, mode] - 1) * ndofs + celldofs[d, cell] : 0
                                dofminus = mneighboursMINUS[m, mode] <= nmodes ? (mneighboursMINUS[m, mode] - 1) * ndofs + celldofs[d, cell] : 0
                                coeffPLUS = dofplus > 0 ? coeffs[dofplus] * G[(m - 1) * nmodes_extended + mode, mneighboursPLUS[m, mode]] : 0.0
                                coeffMINUS = dofminus > 0 ? coeffs[dofminus] * G[(m - 1) * nmodes_extended + mode, mneighboursMINUS[m, mode]] : 0.0
                                coeffs_uh[d] += (coeffPLUS + coeffMINUS) * am[1]
                            end
                        end

                        ## evaluate grad(u_h) and nodal basis function at quadrature point
                        fill!(graduh, 0)
                        eval_febe!(graduh, FEE_∇u, coeffs_uh, i)

                        ## evaluate rhs
                        f_temp[1] = 0
                        if mode == 1
                            rhs(f_temp, x)
                        end

                        ## compute residual -f*phi_z + grad(u_h) * grad(phi_z) at quadrature point i
                        temp2 = div_penalty * sqrt(xCellVolumes[cell]) * weight
                        temp = temp2 * (-f_temp[1] * xref_vals[1, localnode, i] + dot(graduh, view(∇φvals, :, localnode, 1)))
                        for dof_i in 1:maxdofs
                            ## right-hand side for best-approximation (grad(u_h)*phi)
                            blocal[dof_i] += dot(graduh, view(idvals, :, dof_i, i)) * xref_vals[1, localnode, i] * weight
                            ## mass matrix Hdiv
                            for dof_j in dof_i:maxdofs
                                Alocal[dof_i, dof_j] += dot(view(idvals, :, dof_i, i), view(idvals, :, dof_j, i)) * weight
                            end
                            ## div-div matrix Hdiv * penalty (quick and dirty to avoid Lagrange multiplier)
                            blocal[dof_i] += temp * divvals[1, dof_i, i]
                            temp3 = temp2 * divvals[1, dof_i, i]
                            for dof_j in dof_i:maxdofs
                                Alocal[dof_i, dof_j] += temp3 * divvals[1, dof_j, i]
                            end
                        end
                    end

                    ## write into global A and b
                    for dof_i in 1:maxdofs
                        dofi = celldofs_dual[dof_i, cell]
                        b[dofi] += blocal[dof_i]
                        for dof_j in 1:maxdofs
                            dofj = celldofs_dual[dof_j, cell]
                            if dof_j < dof_i # use that Alocal is symmetric
                                _addnz(A, dofi, dofj, Alocal[dof_j, dof_i], 1)
                            else
                                _addnz(A, dofi, dofj, Alocal[dof_i, dof_j], 1)
                            end
                        end
                    end

                    ## reset local A and b
                    fill!(Alocal, 0)
                    fill!(blocal, 0)
                end
            end
        end

        ## penalize dofs that are not involved
        for j in 1:FESDual.ndofs
            if is_noninvolveddof[j]
                A[j, j] = bnd_penalty
                b[j] = 0
            end
        end

        ## solve local problem
        return A \ b

    end

    #Threads.@threads
    for group in groups
        grouptime = @elapsed begin
            @info "  Starting equilibrating patch group $group on thread $(Threads.threadid())... "

            for mode in 1:nmodes_extended
                X[group, mode] = solve_patchgroup!(group, mode)
            end

        end

        @info "Finished equilibration patch group $group on thread $(Threads.threadid()) in $(grouptime)s "
    end

    ## write local solutions to global vector
    for group in 1:maxgroups, mode in 1:nmodes_extended
        view(sol_eq[mode]) .+= X[group, mode]
    end
    coeffs_eq = sol_eq.entries
    ndofs_eq = FESDual.ndofs

    ## compute volume terms
    function compute_volume_terms()

        ## define error estimator : || σ_h - a∇u_h ||^2 + || f + div σ_h ||^2
        #EQIntegrator = ItemIntegrator(eqestimator_kernel!, [grad(1), id(2), div(2)]; resultdim = 1, quadorder = 2 * order)

        eta4cell = zeros(Float64, ncells, nmodes_extended)
        lambda_temp = zeros(Float64, 1)
        f_temp = zeros(Float64, 1)
        am = zeros(Float64, 1)
        sigmatemp = zeros(Float64, 2)
        x = zeros(Float64, 2)

        FEE_div2 = FEEvaluator(FESDual, Divergence, qf)
        FEE_id2 = FEEvaluator(FESDual, ExtendableFEMBase.Identity, qf)
        FEE_∇u2 = FEEvaluator(FES, Gradient, qf)
        ∇uvals = FEE_∇u2.cvals
        idvals = FEE_id2.cvals
        divvals = FEE_div2.cvals
        L2G = L2GTransformer(Triangle2D, xgrid, ON_CELLS)

        for cell in 1:ncells
            FEE_∇u2.citem[] = cell
            FEE_div2.citem[] = cell
            FEE_id2.citem[] = cell
            update_basis!(FEE_∇u2)
            update_basis!(FEE_div2)
            update_basis!(FEE_id2)
            update_trafo!(L2G, cell)

            for j in 1:nmodes_extended

                for qp in eachindex(weights)
                    fill!(sigmatemp, 0)
                    eval_trafo!(x, L2G, xref[qp])

                    ## evaluate flux of discrete solution
                    if j <= nmodes
                        get_am!(am, x, 0, C)
                        for d in 1:maxdofs_uh
                            dof = (j - 1) * ndofs + celldofs[d, cell]
                            for k in 1:2
                                sigmatemp[k] += coeffs[dof] * am[1] * ∇uvals[k, d, qp]
                            end
                        end
                    end
                    for m in 1:Mcoup
                        get_am!(am, x, m, C)
                        for d in 1:maxdofs_uh
                            dofplus = mneighboursPLUS[m, j] <= nmodes ? (mneighboursPLUS[m, j] - 1) * ndofs + celldofs[d, cell] : 0
                            dofminus = mneighboursMINUS[m, j] <= nmodes ? (mneighboursMINUS[m, j] - 1) * ndofs + celldofs[d, cell] : 0
                            coeffPLUS = dofplus > 0 ? coeffs[dofplus] * G[(m - 1) * nmodes_extended + j, mneighboursPLUS[m, j]] : 0.0
                            coeffMINUS = dofminus > 0 ? coeffs[dofminus] * G[(m - 1) * nmodes_extended + j, mneighboursMINUS[m, j]] : 0.0
                            for k in 1:2
                                sigmatemp[k] += (coeffPLUS + coeffMINUS) * am[1] * ∇uvals[k, d, qp]
                            end
                        end
                    end


                    ## subtract equilibrated flux
                    div_sigma = 0.0
                    if j == 1
                        f_temp[1] = 0
                        rhs(f_temp, x)
                        div_sigma += f_temp[1]
                    end
                    for d in 1:maxdofs
                        dof = (j - 1) * ndofs_eq + celldofs_dual[d, cell]
                        for k in 1:2
                            sigmatemp[k] -= coeffs_eq[dof] * idvals[k, d, qp]
                        end
                        div_sigma += coeffs_eq[dof] * divvals[1, d, qp]
                    end

                    eta4cell[cell, j] += (dot(sigmatemp, sigmatemp) + div_sigma^2 / pi^2) * weights[qp] * xCellVolumes[cell]
                end
            end
        end
        return eta4cell
    end

    eta4cell = compute_volume_terms()

    eta4modes = zeros(Float64, nmodes_extended)
    for j in 1:nmodes_extended
        eta4modes[j] = sqrt(sum(view(eta4cell, :, j)))
    end

    ## no data truncation estimate (4th slot keeps the return signature uniform)
    return eta4modes, eta4cell, multi_indices_extended, 0.0
end
