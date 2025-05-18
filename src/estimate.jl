function get_neighbours(OBT, multi_indices)
    M = length(multi_indices[1])
    nmodes = length(multi_indices)
    mneighboursPLUS = zeros(Int, M, nmodes)
    mneighboursMINUS = zeros(Int, M, nmodes)
    for j in 1:nmodes
        for m in 1:M
            mu1 = deepcopy(multi_indices[j])
            mu2 = deepcopy(multi_indices[j])
            mu1[m] += 1
            mu2[m] -= 1
            for k in 1:nmodes
                if all(multi_indices[k] .== mu1)
                    mneighboursPLUS[m, j] = k
                elseif all(multi_indices[k] .== mu2)
                    mneighboursMINUS[m, j] = k
                end
            end
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

"""
$(TYPEDSIGNATURES)

computes the residual-based a posteriori error estimator for the current solution vector `sol`
for the given model problem and stochastic coefficient `C` and returns:
- `eta4modes` = array with total error estimators for each multi-index (corresponding to the enriched set of multi-indices with the current active modes coming first)
- `eta4cell` = array with total error estimators for each cell in grid (for spatial refinement)
- `multi_indices_extended` = enriched set of multi-indices used for the computations

"""
function estimate(::Type{AbstractModelProblem}, sol::SGFEVector, C::AbstractStochasticCoefficient; kwargs...)
    @error "no error estimator for the model problem type available"
    return Vector{Float64}, Matrix{Float64}, [[0]]
end


function estimate(::Type{LogTransformedPoissonProblemPrimal}, sol::SGFEVector, C::AbstractStochasticCoefficient; rhs = nothing, bonus_quadorder = 1, tail_extension = 5)

    FES = sol.FES_space[1]
    FEType = eltype(FES)
    xgrid = FES.xgrid
    ncells = num_cells(xgrid)
    EG = xgrid[UniqueCellGeometries][1]
    TB = sol.TB
    M = maxlength_multiindices(TB)
    nmodes = TB.nmodes
    multi_indices = TB.multi_indices
    OBT = OrthogonalPolynomialType(TB.ONB)
    order = get_polynomialorder(FEType, EG)

    ## extend multi_indices
    multi_indices_extended = add_boundary_modes(deepcopy(multi_indices); tail_extension = tail_extension)
    M_extended = length(multi_indices_extended[1])
    maxorder = maximum(maximum(multi_indices_extended[j]) for j in 1:length(multi_indices_extended))
    TB_extended = TensorizedBasis(OBT, M + 1, maxorder, 2 * maxorder, 2 * maxorder, multi_indices = multi_indices_extended)
    nmodes_extended = length(multi_indices_extended)
    G = TB_extended.G

    ## prepare neighbours of modes
    mneighboursPLUS, mneighboursMINUS = get_neighbours(OBT, multi_indices_extended)

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

    ## prepare expansion of coefficient
    expa_PCE!, lambda_μ! = expa_PCE_mop(TB_extended, C; factor = -1.0)

    ## interpolate <e^-a f, H_nu>
    FES_interp = FESpace{H1Pk{1, 2, order}}(xgrid)
    FEBasis_id = FEEvaluator(FES_interp, ExtendableFEMBase.Identity, qf)
    idvals = FEBasis_id.cvals
    expaf_interpolations = FEVector([FES_interp for j in 1:nmodes_extended])
    for j in 1:nmodes_extended
        interpolate!(expaf_interpolations[j], (result, qpinfo) -> lambda_μ!(result, qpinfo.x, j); quadorder = quadorder)
    end
    offset_interp = FES_interp.ndofs
    celldofs_interp = FES_interp[CellDofs]
    ndofs4cell_interp = get_ndofs(ON_CELLS, H1Pk{1, 2, order}, EG)
    coeffs_interp = expaf_interpolations.entries

    ## compute volume terms
    f4modes = zeros(Float64, nmodes_extended)
    eta4cell = zeros(Float64, ncells, nmodes_extended)
    eta4modes = zeros(Float64, nmodes_extended)
    function barrier(EG, L2G::L2GTransformer)
        lambda_temp = zeros(Float64, 1)
        gradam = zeros(Float64, 2)
        ftemp = zeros(Float64, 1)
        sigmatemp = zeros(Float64, 1)
        x = zeros(Float64, 2)

        for cell in 1:ncells
            update_basis!(FEBasis_∇, cell)
            update_basis!(FEBasis_id, cell)
            if order > 1
                update_basis!(FEBasis_Δ, cell)
            end
            update_trafo!(L2G, cell)
            for j in 1:nmodes_extended
                ## h_T|| f_nu + \sigma_\nu ||

                ## evaluate <e^-a f, H_nu> = lambda_nu * f
                for qp in 1:nweights
                    eval_trafo!(x, L2G, xref[qp])
                    lambda_temp = 0.0
                    for d in 1:ndofs4cell_interp
                        dof = (j - 1) * offset_interp + celldofs_interp[d, cell]
                        lambda_temp += idvals[1, d, qp] * coeffs_interp[dof]
                    end
                    #lambda_μ!(lambda_temp, x, j) # most expensive line
                    rhs(ftemp, x)
                    ftemp[1] = lambda_temp[1] * ftemp[1]
                    if order > 1 && j <= nmodes
                        for d in 1:ndofs4cell
                            dof = (j - 1) * ndofs + celldofs[d, cell]
                            ftemp[1] += coeffs[dof] * Δvals[1, d, qp]
                        end
                    end

                    sigmatemp[1] = 0
                    for m in 1:M_extended
                        get_gradam!(gradam, x, m, C)
                        for d in 1:ndofs4cell
                            dofplus = mneighboursPLUS[m, j] <= nmodes ? (mneighboursPLUS[m, j] - 1) * ndofs + celldofs[d, cell] : 0
                            dofminus = mneighboursMINUS[m, j] <= nmodes ? (mneighboursMINUS[m, j] - 1) * ndofs + celldofs[d, cell] : 0
                            coeffPLUS = dofplus > 0 ? coeffs[dofplus] * G[(m - 1) * nmodes_extended + j, mneighboursPLUS[m, j]] : 0.0
                            coeffMINUS = dofminus > 0 ? coeffs[dofminus] * G[(m - 1) * nmodes_extended + j, mneighboursMINUS[m, j]] : 0.0
                            sigmatemp[1] += (coeffPLUS + coeffMINUS) * dot(gradam, view(∇vals, :, d, qp))
                        end
                    end
                    eta4cell[cell, j] += (ftemp[1] + sigmatemp[1])^2 * weights[qp]
                    f4modes[j] += ftemp[1]^2 * weights[qp] * cellvolumes[cell]
                end

                ## boundary modes enjoy no Galerkin orthogonality and therefore have no additional h power
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
        for cell in 1:ncells, f in 1:3
            eta4cell[cell, j] += jumps4face[cellfaces[f, cell]]
        end
        eta4modes[j] += sqrt(eta4modes[j]^2 + sum(view(jumps4face, :)))
    end

    return eta4modes, eta4cell, multi_indices_extended #, f4modes
end


function estimate(::Type{PoissonProblemPrimal}, sol::SGFEVector, C::AbstractStochasticCoefficient; rhs = nothing, bonus_quadorder = 1, tail_extension = 5)

    FES = sol.FES_space[1]
    FEType = eltype(FES)
    xgrid = FES.xgrid
    ncells = num_cells(xgrid)
    EG = xgrid[UniqueCellGeometries][1]
    TB = sol.TB
    M = maxlength_multiindices(TB)
    nmodes = TB.nmodes
    multi_indices = TB.multi_indices
    OBT = OrthogonalPolynomialType(TB.ONB)
    order = get_polynomialorder(FEType, EG)

    ## extend multi_indices
    multi_indices_extended = add_boundary_modes(deepcopy(multi_indices); tail_extension = tail_extension)
    M_extended = length(multi_indices_extended[1])
    maxorder = maximum(maximum(multi_indices_extended[j]) for j in 1:length(multi_indices_extended))
    TB_extended = TensorizedBasis(OBT, M + 1, maxorder, 2 * maxorder, 2 * maxorder, multi_indices = multi_indices_extended)
    nmodes_extended = length(multi_indices_extended)
    G = TB_extended.G

    ## prepare neighbours of modes
    mneighboursPLUS, mneighboursMINUS = get_neighbours(OBT, multi_indices_extended)

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

    ## compute volume terms
    eta4cell = zeros(Float64, ncells, nmodes_extended)
    eta4modes = zeros(Float64, nmodes_extended)
    f4modes = zeros(Float64, nmodes_extended)
    function barrier(EG, L2G::L2GTransformer)
        lambda_temp = zeros(Float64, 1)
        am = zeros(Float64, 1)
        ftemp = zeros(Float64, 1)
        sigmatemp = zeros(Float64, 1)
        x = zeros(Float64, 2)

        for cell in 1:ncells
            update_basis!(FEBasis_∇, cell)
            if order > 1
                update_basis!(FEBasis_Δ, cell)
            end
            update_trafo!(L2G, cell)
            for j in 1:nmodes_extended
                ## h_T|| f_nu + div \sigma_\nu ||

                for qp in 1:nweights
                    eval_trafo!(x, L2G, xref[qp])

                    ftemp[1] = 0
                    if j == 1
                        rhs(ftemp, x)
                    end
                    if order > 1
                        if j <= nmodes
                            get_am!(am, x, 0, C)
                            for d in 1:ndofs4cell
                                dof = (j - 1) * ndofs + celldofs[d, cell]
                                ftemp[1] += coeffs[dof] * am[1] * Δvals[1, d, qp]
                            end
                        end
                        for m in 1:M_extended
                            get_am!(am, x, m, C)
                            for d in 1:ndofs4cell
                                dofplus = mneighboursPLUS[m, j] <= nmodes ? (mneighboursPLUS[m, j] - 1) * ndofs + celldofs[d, cell] : 0
                                dofminus = mneighboursMINUS[m, j] <= nmodes ? (mneighboursMINUS[m, j] - 1) * ndofs + celldofs[d, cell] : 0
                                coeffPLUS = dofplus > 0 ? coeffs[dofplus] * G[(m - 1) * nmodes_extended + j, mneighboursPLUS[m, j]] : 0.0
                                coeffMINUS = dofminus > 0 ? coeffs[dofminus] * G[(m - 1) * nmodes_extended + j, mneighboursMINUS[m, j]] : 0.0
                                ftemp[1] += (coeffPLUS + coeffMINUS) * am[1] * Δvals[1, d, qp]
                            end
                        end
                    end
                    eta4cell[cell, j] += ftemp[1]^2 * weights[qp]
                    f4modes[j] += ftemp[1]^2 * weights[qp] * cellvolumes[cell]
                end

                ## boundary modes enjoy no Galerkin orthogonality and therefore have no additional h power
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

    ## compute normal jumps, todo: am is missing
    sol_j = FEVector(FES)
    cellfaces = xgrid[CellFaces]
    jumps4face = zeros(Float64, 1, size(xgrid[FaceNodes], 2))

    m_pointer = [0]
    JumpEvaluator = FaceInterpolator((result, input, qpinfo) -> (get_am!(result, qpinfo.x, m_pointer[1], C); result .= result[1] * input), [jump(grad(1))], quadorder = quadorder)
    ExtendableFEM.build_assembler!(JumpEvaluator, [sol_j[1]])
    sol_jumps = deepcopy(JumpEvaluator.value)
    JumpIntegrator = L2NormIntegrator([id(1)]; entities = ON_FACES)
    bfaces = xgrid[BFaceFaces]

    @time for j in 1:nmodes_extended
        fill!(sol_jumps[1], 0)
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
        end
        jumps4face .= sum(ExtendableFEM.evaluate(JumpIntegrator, sol_jumps), dims = 1)
        jumps4face[1, bfaces] .= 0

        if j <= nmodes
            jumps4face[1, :] .*= xgrid[FaceVolumes]
        else
            jumps4face[1, :] ./= xgrid[FaceVolumes]
        end

        for cell in 1:ncells
            for f in 1:3
                eta4cell[cell, j] += jumps4face[cellfaces[f, cell]]
            end
        end
        eta4modes[j] = sqrt(eta4modes[j]^2 + sum(view(jumps4face, :)))
    end

    return eta4modes, eta4cell, multi_indices_extended #, f4modes
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
function estimate_equilibration(::Type{PoissonProblemPrimal}, sol::SGFEVector, C::AbstractStochasticCoefficient; FETypeDual = nothing, rhs = nothing, tail_extension = 5, bonus_quadorder = 0)
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
    M = maxlength_multiindices(TB)
    nmodes = TB.nmodes
    multi_indices = TB.multi_indices
    OBT = OrthogonalPolynomialType(TB.ONB)

    ## extend multi_indices
    multi_indices_extended = add_boundary_modes(deepcopy(multi_indices); tail_extension = tail_extension)
    M_extended = length(multi_indices_extended[1])
    maxorder = maximum(maximum(multi_indices_extended[j]) for j in 1:length(multi_indices_extended))
    TB_extended = TensorizedBasis(OBT, M + 1, maxorder, 2 * maxorder, 2 * maxorder, multi_indices = multi_indices_extended)
    nmodes_extended = length(multi_indices_extended)
    G = TB_extended.G

    ## prepare neighbours of modes
    mneighboursPLUS, mneighboursMINUS = get_neighbours(OBT, multi_indices_extended)

    ## append block in solution vector for equilibrated fluxes
    sol_eq = FEVector(FESDual)
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
                        for m in 1:M_extended
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
                    for m in 1:M_extended
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

    return eta4modes, eta4cell, multi_indices_extended #, f4modes
end
