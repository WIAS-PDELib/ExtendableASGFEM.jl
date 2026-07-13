module StokesSimple

using ExtendableASGFEM
using ExtendableFEM
using ExtendableFEMBase
using ExtendableGrids
using GridVisualize
using Symbolics

function prepare_data()
    @variables x1, x2

    u = [0, (1 / 2) * (1 - x1^2)]
    eval_u! = build_function(u, x1, x2, expression = Val{false})[2]

    return (result, _) -> (result .= 0),
        (result, qpinfo) -> (eval_u!(result, qpinfo.x...))
end

function main(;
        problem = StokesProblemPrimal,
        nrefs = 3,                      ## number of uniform refinements of the initial grid
        order = 2,                      ## polynomial order of the FEspaces
        domain = "square",              ## domain, e.g., "square" or "lshape"
        modes = [[0], [1]],             ## initial multi-indices for stochastic basis
        use_iterative_solver = true,    ## use iterative solver? (otherwise direct)
        nsamples = 50,
        Plotter = nothing,
    )
    ## prepare stochastic coefficient
    C = StochasticCoefficientConstants()

    ## prepare grid
    (xgrid, boundary_regions) = if domain == "square"
        (uniform_refine(grid_unitsquare(Triangle2D), nrefs), 1:4)
    elseif domain == "lshape"
        (uniform_refine(grid_lshape(Triangle2D), nrefs), 1:8)
    else
        error("unknown domain: $domain")
    end

    ## prepare stochastic basis
    multi_indices = Array{Array{Int, 1}, 1}(modes)
    prepare_multi_indices!(multi_indices)
    M = maximum(length.(multi_indices))
    OBType = LegendrePolynomials
    ansatz_deg = maximum([maximum(multi_indices[k]) for k in 1:length(multi_indices)]) + 4
    TensorBasis = TensorizedBasis(OBType, M, ansatz_deg, 2 * ansatz_deg, 2 * ansatz_deg, multi_indices = multi_indices)

    ## prepare synthetic problem data
    f!, exact_u! = prepare_data()
    bonus_quadorder_f = 2

    ## prepare FE spaces
    FETypes = (H1BR{2}, L2P0{1})
    FES = [FESpace{FETypes[1]}(xgrid), FESpace{FETypes[2]}(xgrid)]
    unames = ["u", "p"]

    ## create solution vector
    sol = SGFEVector(FES, TensorBasis; active_modes = 1:length(multi_indices), unames)

    ## solve problem
    @info "Solving..."
    solve!(problem, sol, f!, C; (exact_boundary!) = (exact_u!), bonus_quadorder_f, use_iterative_solver)

    ## plot solution
    if !isnothing(Plotter)
        p = plot_modes(sol; Plotter = Plotter, ncols = 4)
        display(p)
    end

    ## compute exact error (by MC sampling)
    weightederrorH1, weightederrorL2u, weightederrorL2p, uniformerrorH1, uniformerrorL2u, uniformerrorL2p =
        calculate_sampling_error_2(
        sol,
        f!,
        C;
        problem,
        metrics_configurations = stokes_metrics_configuration,
        (exact_boundary!) = (exact_u!),
        boundary_regions,
        order = order + 1,
        nsamples,
    )

    @info "RESULTS
        || ∇(u - u_h) || (w,u) = $(sqrt(weightederrorH1[end])), $(sqrt(uniformerrorH1[end]))
        || u - u_h || (w,u) = $(sqrt(weightederrorL2u[end])), $(sqrt(uniformerrorL2u[end]))
        || p - p_h || (w,u) = $(sqrt(weightederrorL2p[end])), $(sqrt(uniformerrorL2p[end]))"

    return sol
end

end
