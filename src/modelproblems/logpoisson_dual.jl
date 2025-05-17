"""
$(TYPEDEF)

Dual log-transformed formulation of the Poisson problem with exponential stochastic coefficient (WIP).
"""
abstract type LogTransformedPoissonProblemDual <: AbstractModelProblem end

#########################
### ASSEMBLY & SOLVER ###
#########################

include("solvers_logpoisson_dual.jl")

## deterministic problem description
function deterministic_problem(::Type{LogTransformedPoissonProblemDual}, C::AbstractStochasticCoefficient, sample_pointer; rhs = nothing, bonus_quadorder_a = 2, bonus_quadorder_f = 0)
    grada! = get_grada!(C; factor = -1)

    function u_times_grada!(result, input, qpinfo)
        grada!(result, qpinfo.x, sample_pointer)
        result .*= input[1]
        return nothing
    end

    PD = ProblemDescription("log-transformed Poisson problem (dual)")
    p = Unknown("p"; name = "stress p = -∇ũ - ∇aũ = - epx(-a) ∇u")
    u = Unknown("u"; name = "potential ũ = exp(-a) u")
    assign_unknown!(PD, p)
    assign_unknown!(PD, u)
    assign_operator!(PD, BilinearOperator([id(p)]; store = true))
    assign_operator!(PD, BilinearOperator([div(p)], [id(u)]; store = true, factor = -1, transposed_copy = 1))
    assign_operator!(PD, BilinearOperator(u_times_grada!, [id(p)], [id(u)]; factor = 1, bonus_quadorder = bonus_quadorder_a))
    if rhs !== nothing
        assign_operator!(PD, LinearOperator(rhs, [id(u)]; factor = -1, store = true, bonus_quadorder = bonus_quadorder_f))
    end
    #assign_operator!(PD, HomogeneousBoundaryData(u; regions = 1:4))

    return PD, u
end

function solve!(
        ::Type{LogTransformedPoissonProblemDual},
        sol::SGFEVector,                        ## target SGFEM vector
        C::AbstractStochasticCoefficient;       ## stochastic coefficient a
        bonus_quadorder_a = 2,                  ## additional quadrature order for grad(am)
        bonus_quadorder_f = 0,                  ## additional quadrature order for rhs f
        rhs = nothing,
        debug = false,
        use_iterative_solver = true
    )

    ## create FESpace for space discretisation
    FES = sol.FES_space
    TensorBasis = sol.TB
    xgrid = FES[1].xgrid
    dim = size(xgrid[Coordinates], 1)

    ## assemble mass matrix: A = (p,q)
    A = FEMatrix(FES[1])
    B = FEMatrix(FES[1], FES[2])
    assemble!(A, BilinearOperator([id(1)]; factor = 1))
    assemble!(B, BilinearOperator([div(1)], [id(1)]; factor = -1))

    ## assemble convection term: N[m] = - (u ∇a[m], ∇v)
    N0 = FEMatrix(FES[1], FES[2])
    N = []
    for m in 1:maxlength_multiindices(TensorBasis)
        Nm = FEMatrix(FES[1], FES[2])
        assemble!(Nm, BilinearOperator(get_gradam_x_u(m, C), [id(1)], [id(1)]; factor = -1, bonus_quadorder = bonus_quadorder_a))
        push!(N, Nm)
    end

    ## assemble right-hand side: b[μ] = (λ[μ], v)
    if rhs !== nothing
        b0 = FEVector(FES[2])
        assemble!(b0, LinearOperator(rhs, [id(1)]; factor = -1, bonus_quadorder = bonus_quadorder_f))

    else
        @error "need right-hand side"
    end

    ## get G matrix
    G = TensorBasis.G
    nmodes = num_multiindices(TensorBasis)

    ## solve
    if use_iterative_solver
        @time bdofs = solve_logpoisson_dual!(sol, A, B, N0, N, b0, G, nmodes, 1)
    else
        @time bdofs = solve_logpoisson_dual_full!(sol, A, B, N0, N, b0, G, nmodes, 1)
    end

    return bdofs
end
