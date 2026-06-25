module StokesSimple

using ExtendableASGFEM
using ExtendableFEM
using ExtendableFEMBase
using ExtendableGrids
using GridVisualize

function f!(result, qpinfo)
    x = qpinfo.x[1]
    y = qpinfo.x[2]
    result[1] = 5(x^4) + 12(x^2) * ((-1 + x)^2) * (-1 + y) + 12(x^2) * ((-1 + x)^2) * y + 4(x^2) * y * ((-1 + y)^2) + 4(x^2) * (-1 + y) * (y^2) + 4((-1 + x)^2) * (-1 + y) * (y^2) + 4((-1 + x)^2) * y * ((-1 + y)^2) + 16x * (-1 + x) * (-1 + y) * (y^2) + 16x * (-1 + x) * y * ((-1 + y)^2)
    result[2] = 5(y^4) - 4x * ((-1 + x)^2) * ((-1 + y)^2) - 4(-1 + x) * (x^2) * ((-1 + y)^2) - 4x * ((-1 + x)^2) * (y^2) - 4(-1 + x) * (x^2) * (y^2) - 16(-1 + x) * (x^2) * y * (-1 + y) - 16x * ((-1 + x)^2) * y * (-1 + y) - 12x * (y^2) * ((-1 + y)^2) - 12(-1 + x) * (y^2) * ((-1 + y)^2)
    return nothing
    #result[1] = 5(x^4) - 0.001(-12(x^2) * ((-1 + x)^2) * (-1 + y) - 12(x^2) * ((-1 + x)^2) * y - 4(x^2) * y * ((-1 + y)^2) - 4(x^2) * (-1 + y) * (y^2) - 4((-1 + x)^2) * (-1 + y) * (y^2) - 4((-1 + x)^2) * y * ((-1 + y)^2) - 16x * (-1 + x) * (-1 + y) * (y^2) - 16x * (-1 + x) * y * ((-1 + y)^2))
    #result[2] = 5(y^4) - 0.001(4x * ((-1 + x)^2) * ((-1 + y)^2) + 4(-1 + x) * (x^2) * ((-1 + y)^2) + 4x * ((-1 + x)^2) * (y^2) + 4(-1 + x) * (x^2) * (y^2) + 16(-1 + x) * (x^2) * y * (-1 + y) + 16x * ((-1 + x)^2) * y * (-1 + y) + 12x * (y^2) * ((-1 + y)^2) + 12(-1 + x) * (y^2) * ((-1 + y)^2))
    #result[1] = 1000 * (1 - x) * x * y * (1 - y)
    #result[2] = 1000 * (1 - x) * x * y * (1 - y)
end

function main(;
        problem = StokesProblemPrimal,
        nrefs = 3,      # number of uniform refinements of the initial grid
        order = 2,      # polynomial order of the FEspaces
        decay = 2.0,    # decay factor for the random coefficient
        mean = 1, # mean value of coefficient
        domain = "square",  # domain, e.g., "square" or "lshape"
        initial_modes = [[0], [1, 0], [0, 1], [2, 0], [0, 0, 1]],   # initial multi-indices for stochastic basis
        (f!) = (f!),       # right-hand side function
        use_iterative_solver = true,    # use iterative solver ? (otherwise direct)
        nsamples = 20,
        Plotter = nothing,
        debug = false,
    )
    ## prepare stochastic coefficient
    C = StochasticCoefficientCosinus(; τ = 0.9, decay = decay, mean = mean)

    ## prepare grid
    xgrid = if domain == "square"
        uniform_refine(grid_unitsquare(Triangle2D), nrefs)
    elseif domain == "lshape"
        uniform_refine(grid_lshape(Triangle2D), nrefs)
    else
        error("unknown domain: $domain")
    end

    ## prepare stochastic basis
    multi_indices = Array{Array{Int, 1}, 1}(initial_modes)
    prepare_multi_indices!(multi_indices)
    M = maximum(length.(multi_indices))
    OBType = LegendrePolynomials
    ansatz_deg = maximum([maximum(multi_indices[k]) for k in 1:length(multi_indices)]) + 4
    TensorBasis = TensorizedBasis(OBType, M, ansatz_deg, 2 * ansatz_deg, 2 * ansatz_deg, multi_indices = multi_indices)

    ## prepare FE spaces
    FETypes = (H1BR{2}, L2P0{1})
    FES = [FESpace{FETypes[1]}(xgrid), FESpace{FETypes[2]}(xgrid)]
    unames = ["u", "p"]

    ## create solution vector
    sol = SGFEVector(FES, TensorBasis; active_modes = 1:length(multi_indices), unames = unames)

    ## solve problem
    @info "Solving..."
    solve!(problem, sol, C; rhs = (f!), use_iterative_solver = use_iterative_solver)

    ## plot solution
    if !isnothing(Plotter)
        p = plot_modes(sol; Plotter = Plotter, ncols = 4)
        display(p)
    end

    ## compute exact error (by MC sampling)
    weightederrorH1, weightederrorL2u, weightederrorL2p, uniformerrorH1, uniformerrorL2u, uniformerrorL2p = calculate_sampling_error_2(
        sol, C; problem = problem, metrics_configurations = stokes_metrics_configuration, rhs = (f!), order = order + 1, nsamples, debug
    )

    @info "RESULTS
        || ∇(u-u_h) || (w,u) = $(sqrt(weightederrorH1[end])), $(sqrt(uniformerrorH1[end]))
        || u - u_h || (w,u) = $(sqrt(weightederrorL2u[end])), $(sqrt(uniformerrorL2u[end]))
        || p - p_h || (w,u) = $(sqrt(weightederrorL2p[end])), $(sqrt(uniformerrorL2p[end]))"

    return sol
end

end
