"""
$(TYPEDEF)

Log-transformed formulation of the Poisson problem with exponential stochastic coefficient.
"""
abstract type LogTransformedPoissonProblemPrimal <: AbstractModelProblem end

include("solvers_logpoisson_primal.jl")

## deterministic problem description
function deterministic_problem(::Type{LogTransformedPoissonProblemPrimal}, C::AbstractStochasticCoefficient, sample_pointer; rhs = nothing, bonus_quadorder_a = 2, bonus_quadorder_f = 0)
    ## kernel for diffusion operator (= exp(a))

    expa! = get_expa!(C; factor = 1)

    function kernel_diffusion_expa!(result, input, qpinfo)
        expa!(result, qpinfo.x, sample_pointer)
        result .= result[1] * input
        return nothing
    end

    PD = ProblemDescription("log-transformed Poisson problem (primal)")
    u = Unknown("u"; name = "potential")
    assign_unknown!(PD, u)
    assign_operator!(PD, BilinearOperator(kernel_diffusion_expa!, [grad(u)]; bonus_quadorder = bonus_quadorder_a))
    if rhs !== nothing
        assign_operator!(PD, LinearOperator(rhs, [id(u)]; store = true, bonus_quadorder = bonus_quadorder_f))
    end
    assign_operator!(PD, HomogeneousBoundaryData(u; regions = 1:8))

    return PD, u
end

function deterministic_problem2(::Type{LogTransformedPoissonProblemPrimal}, C::AbstractStochasticCoefficient, sample_pointer; rhs = nothing, bonus_quadorder_a = 2, bonus_quadorder_f = 0)
    ## kernel for diffusion operator (= exp(a))

    expa! = get_expa!(C; factor = -1)
    grad_a! = get_grada!(C; factor = 1)

    grad_a = zeros(Float64, 2)
    function kernel_convection!(result, input, qpinfo)
        fill!(grad_a, 0)
        grad_a!(grad_a, qpinfo.x, sample_pointer)
        result .= dot(grad_a, input)
        return nothing
    end

    PD = ProblemDescription("log-transformed Poisson problem (primal)")
    u = Unknown("u"; name = "potential")
    assign_unknown!(PD, u)
    assign_operator!(PD, BilinearOperator([grad(u)]; bonus_quadorder = bonus_quadorder_a))
    assign_operator!(PD, BilinearOperator(kernel_convection!, [id(1)], [grad(1)]; factor = -1, bonus_quadorder = bonus_quadorder_a))

    if rhs !== nothing
        f_eval = zeros(Float64, 1)
        function kernel_fexp(result, qpinfo)
            result[1] = 0
            expa!(result, qpinfo.x, sample_pointer)
            rhs(f_eval, qpinfo.x)
            return result[1] = result[1] * f_eval[1]
        end
        assign_operator!(PD, LinearOperator(kernel_fexp, [id(1)]; bonus_quadorder = max(bonus_quadorder_a, bonus_quadorder_f)))
    end
    assign_operator!(PD, HomogeneousBoundaryData(u; regions = 1:8))

    return PD, u
end

## solver for stochastic Galerkin problem
function solve!(
        ::Type{LogTransformedPoissonProblemPrimal},
        sol::SGFEVector,                        ## target SGFEM vector
        C::AbstractStochasticCoefficient;       ## stochastic coefficient a
        bonus_quadorder_a = 2,                  ## additional quadrature order for grad(am)
        bonus_quadorder_f = 0,                  ## additional quadrature order for rhs f
        rhs = nothing,
        debug = false,
        use_iterative_solver = true
    )

    ## create FESpace for space discretisation
    FES = sol.FES_space[1]
    TensorBasis = sol.TB
    xgrid = FES.xgrid
    dim = size(xgrid[Coordinates], 1)

    ## assemble Laplacian: A = (∇u,∇v)
    A = FEMatrix(FES, FES)
    assemble!(A, BilinearOperator([grad(1)]; factor = 1))

    ## assemble convection term: N[m] = - (∇a[m] v, ∇u)
    N = []
    for m in 1:maxlength_multiindices(TensorBasis)
        Nm = FEMatrix(FES, FES)
        assemble!(Nm, BilinearOperator(get_gradam_x_sigma(dim, m, C), [id(1)], [grad(1)]; factor = -1, bonus_quadorder = bonus_quadorder_a))
        push!(N, Nm)
    end

    ## assemble right-hand side: b[μ] = (λ[μ], v)
    if rhs !== nothing
        expa_PCE!, lambda_μ! = expa_PCE_mop(TensorBasis, C; factor = -1.0)
        f_eval = zeros(Float64, 1)
        function kernel_fexp(result, qpinfo)
            result[1] = 0
            lambda_μ!(result, qpinfo.x, qpinfo.params[1])
            rhs(f_eval, qpinfo.x)
            return result[1] = result[1] * f_eval[1] # --> (λ_μ f, •)
        end
        b = []
        nmodes = num_multiindices(TensorBasis)
        for m in 1:nmodes
            bm = FEVector(FES)
            assemble!(bm, LinearOperator(kernel_fexp, [id(1)]; params = [m], bonus_quadorder = bonus_quadorder_f))
            #@info "b[$m] = $(bm.entries)"
            push!(b, bm)
            if debug
                Ifλ = FEVector(FES)
                interpolate!(Ifλ[1], kernel_fexp; params = [m])
                println(stdout, unicode_scalarplot(Ifλ[1]; title = "f*λ[μ] for m=$m"))
            end
        end
    else
        @error "need right-hand side"
    end

    ## get G matrix
    G = TensorBasis.G

    ## solve
    if use_iterative_solver
        @time bdofs = solve_logpoisson_primal!(sol, A, FEMatrix(FES), N, b, G, nmodes, 1)
    else
        @time bdofs = solve_logpoisson_primal_full!(sol, A, FEMatrix(FES), N, b, G, nmodes, 1)
    end

    return bdofs
end
