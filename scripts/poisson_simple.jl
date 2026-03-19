#= 
([source code](SOURCE_URL))

minimalistic script, that solves a stochastic Poisson problem on a uniform mesh
 
usage:
- run main: main(; problem = problem, kwargs...)
=#

module PoissonSimple

using ExtendableASGFEM
using ExtendableFEM
using ExtendableFEMBase
using ExtendableGrids
using GridVisualize

function main(;
        problem = PoissonProblemPrimal,
        nrefs = 3,      # number of uniform refinements of the initial grid
        order = 2,      # polynomial order of the FEspaces
        decay = 2.0,    # decay factor for the random coefficient
        mean = problem == PoissonProblemPrimal ? 1.0 : 0.0, # mean value of coefficient
        domain = "square",  # domain, e.g., "square" or "lshape"
        initial_modes = [[0], [1,0], [0,1], [2,0], [0,0,1]],   # initial multi-indices for stochastic basis
        f! = (result, qpinfo) -> (result[1] = 1),       # right-hand side function
        use_iterative_solver = true,    # use iterative solver ? (otherwise direct)
        Plotter = nothing)

    ## prepare stochastic coefficient
    τ = (problem <: PoissonProblemPrimal) ? 0.9 : 1.0
    if problem <: PoissonProblemPrimal
        @assert mean >= 1 "coefficient needs to be at least 1 to ensure ellipticity"
    end
    C = StochasticCoefficientCosinus(; τ = τ, decay = decay, mean = mean)
    
    ## prepare grid
    xgrid = if domain == "square"
        uniform_refine(grid_unitsquare(Triangle2D), nrefs)
    elseif domain == "lshape"
        uniform_refine(grid_lshape(Triangle2D), nrefs)
    else
        error("unknown domain: $domain")
    end

    ## prepare stochastic basis
    multi_indices = Array{Array{Int,1},1}(initial_modes)
    prepare_multi_indices!(multi_indices)
    M = maximum(length.(multi_indices))
    OBType = problem <: PoissonProblemPrimal ? LegendrePolynomials : HermitePolynomials
    ansatz_deg = maximum([maximum(multi_indices[k]) for k in 1:length(multi_indices)]) + 4
    TensorBasis = TensorizedBasis(OBType, M, ansatz_deg, 2*ansatz_deg, 2*ansatz_deg, multi_indices = multi_indices)

    ## prepare FE spaces
    if problem <: LogTransformedPoissonProblemDual
        FEType = [HDIVRTk{2,order}, order == 0 ? L2P0{1} : H1Pk{1,2,order}]
        FES = [FESpace{FEType[1]}(xgrid), FESpace{FEType[2]}(xgrid; broken = true)]
        unames = ["p", "u"]
    else
        FEType = H1Pk{1,2,order}
        FES = FESpace{FEType}(xgrid)
        unames = ["u"]
    end
    
    ## create solution vector
    sol = SGFEVector(FES, TensorBasis; active_modes = 1:length(multi_indices), unames = unames)

    ## solve problem
    @info "Solving..."
    solve!(problem, sol, C; rhs = f!, use_iterative_solver = use_iterative_solver)

    ## compute exact error (by MC sampling)
    weightederrorH1, weightederrorL2, uniformerrorH1, uniformerrorL2 = calculate_sampling_error(sol, C; problem = problem, rhs = f!, order = order + 1, nsamples = 50)

    ## plot solution
    if !isnothing(Plotter)
        p = plot_modes(sol; Plotter = Plotter, ncols = 4)
        display(p)
    end

    @info "RESULTS
        || ∇(u-u_h) || (w,u) = $(sqrt(weightederrorH1[end])), $(sqrt(uniformerrorH1[end]))
        || u - u_h || (w,u) = $(sqrt(weightederrorL2[end])), $(sqrt(uniformerrorL2[end]))"

    return sol
end

end # module