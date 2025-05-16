
"""
$(TYPEDEF)

Poisson problem with linear stochastic coefficient.
"""
abstract type PoissonProblemPrimal <: AbstractModelProblem end

include("solvers_poisson_primal.jl")

## deterministic problem description
function deterministic_problem(::Type{PoissonProblemPrimal}, C::AbstractStochasticCoefficient, sample_pointer; get_a! = get_a!, rhs = nothing, bonus_quadorder_a = 2, bonus_quadorder_f = 0)
    ## kernel for diffusion operator (= exp(a))

    a! = get_a!(C; factor = 1)

    function kernel_diffusion_a!(result, input, qpinfo)
        a!(result, qpinfo.x, sample_pointer)
        result .= result[1] * input
        return nothing
    end

    PD = ProblemDescription("Poisson problem (primal)")
    u = Unknown("u"; name = "potential")
    assign_unknown!(PD, u)
    assign_operator!(PD, BilinearOperator(kernel_diffusion_a!, [grad(u)]; bonus_quadorder = bonus_quadorder_a))
    if rhs !== nothing
        assign_operator!(PD, LinearOperator(rhs, [id(u)]; store = true, bonus_quadorder = bonus_quadorder_f))
    end
    assign_operator!(PD, HomogeneousBoundaryData(u; regions = 1:8))

    return PD, u
end

## solver for stochastic Galerkin problem
function solve!(
    ::Type{PoissonProblemPrimal},
    sol::SGFEVector,                        ## target SGFEM vector
    C::AbstractStochasticCoefficient;       ## stochastic coefficient a
    bonus_quadorder_a = 2,                  ## additional quadrature order for grad(am)
    bonus_quadorder_f = 0,                  ## additional quadrature order for rhs f
    rhs = nothing,
    debug = false,
	use_iterative_solver = true)

	## create FESpace for space discretisation
	FES = sol.FES_space[1]
    TB = sol.TB
    xgrid = FES.xgrid
    dim = size(xgrid[Coordinates], 1)

	## Laplacian: A = (a∇u,∇v) with building blocks Am = (a_m∇u,∇v)
    A0 = FEMatrix(FES, FES)
    assemble!(A0, BilinearOperator(get_am_x(0, C), [grad(1)], [grad(1)]; bonus_quadorder = bonus_quadorder_a))
	A = []
	for m = 1 : maxlength_multiindices(TB)
		Am = FEMatrix(FES, FES)
		assemble!(Am, BilinearOperator(get_am_x(m, C), [grad(1)], [grad(1)]; bonus_quadorder = bonus_quadorder_a))
		push!(A, Am)
	end
    
    ## assemble right-hand side: b = (f, v)
    if rhs !== nothing
        b = FEVector(FES)
        assemble!(b, LinearOperator(rhs, [id(1)]; bonus_quadorder = bonus_quadorder_f))
    else
        @error "need right-hand side"
    end

    ## solve
	if use_iterative_solver
		@time bdofs = solve_primal!(sol,A0,A,b,TB.G,TB.nmodes,1)
	else
		@time bdofs = solve_full_primal!(sol,A0,A,b,TB.G,TB.nmodes,1)
	end

    return bdofs
end