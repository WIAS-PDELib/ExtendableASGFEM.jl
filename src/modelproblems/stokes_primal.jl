"""

Stokes problem with linear stochastic coefficient ``ν`` that seeks ``(u, p)`` such that

``-\\nu\\Delta u + \\nabla p = f \\quad \\text{for } (y,x) \\in \\Gamma \\times D``
``-\\nu\\mathrm{div} u = 0 \\quad \\text{for } (y,x) \\in \\Gamma \\times D``

"""
abstract type StokesProblemPrimal <: AbstractModelProblem end

include("solvers_stokes_primal.jl")

## deterministic problem description
function deterministic_problem(
        ::Type{StokesProblemPrimal},
        rhs!::Function,
        C::AbstractStochasticCoefficient,
        sample_pointer;
        (get_a!) = (get_a!),
        (exact_boundary!) = nothing,
        reconstruct = true,
        boundary_regions = 1:4,
        bonus_quadorder_a = 2,
        bonus_quadorder_f = 0,
    )
    get_ν! = get_a!(C)
    ν = zeros(1)
    function stokes_kernel!(result, input, qpinfo)
        get_ν!(ν, qpinfo.x, sample_pointer)
        result[1] = ν[1] * input[1] - input[5]
        result[2] = ν[1] * input[2]
        result[3] = ν[1] * input[3]
        result[4] = ν[1] * input[4] - input[5]
        result[5] = -(input[1] + input[4])
        return nothing
    end

    PD = ProblemDescription("Stokes problem (primal)")
    u = Unknown("u", name = "velocity", dim = 2)
    p = Unknown("p", name = "pressure", dim = 1)
    assign_unknown!(PD, u)
    assign_unknown!(PD, p)
    assign_operator!(PD, BilinearOperator(stokes_kernel!, [grad(u), id(p)]; bonus_quadorder = bonus_quadorder_a))
    assign_operator!(
        PD,
        LinearOperator(
            rhs!, [reconstruct ? apply(u, Reconstruct{HDIVRT1{2}, Identity}) : id(u)];
            bonus_quadorder = bonus_quadorder_f
        )
    )
    # Is a non-zero boundary condition or some non-standard definition required?
    if !isnothing(exact_boundary!)
        assign_operator!(
            PD,
            InterpolateBoundaryData(
                u, exact_boundary!;
                regions = boundary_regions, bonus_quadorder = bonus_quadorder_f
            )
        )
    else # Define zero boundary if no additional constraints are defined
        assign_operator!(PD, HomogeneousBoundaryData(u; regions = boundary_regions))
    end
    assign_restriction!(PD, ZeroMeanValueRestriction(p))

    return PD, u
end

## solver for stochastic Galerkin problem
function solve!(
        ::Type{StokesProblemPrimal},
        sol::SGFEVector,                     ## target SGFEM vector
        rhs!::Function,                      ## kernel for right-hand side calculation
        C::AbstractStochasticCoefficient;    ## stochastic coefficient ν
        (exact_boundary!) = nothing,         ## kernel function for boundary values
        bonus_quadorder_a = 2,               ## additional quadrature order for grad(am)
        bonus_quadorder_f = 0,               ## additional quadrature order for rhs f
        use_iterative_solver = true,         ## use iterative solver? (otherwise direct)
        reconstruct = false,                 ## apply HDIVRT0{2} reconstruction operator to rhs?
    )
    FES = sol.FES_space
    TB = sol.TB

    ## Laplacian: A = (a∇u,∇v) with building blocks Am = (a_m∇u,∇v)
    A0 = FEMatrix(FES[1])
    assemble!(A0, BilinearOperator(get_am_x(0, C), [grad(1)], [grad(1)]; bonus_quadorder = bonus_quadorder_a))
    A = []

    for m in 1:maxlength_multiindices(TB)
        Am = FEMatrix(FES[1])
        assemble!(Am, BilinearOperator(get_am_x(m, C), [grad(1)], [grad(1)]; bonus_quadorder = bonus_quadorder_a))
        push!(A, Am)
    end

    ## Deterministic B = (∇ ⋅ v, q) term
    B = FEMatrix(FES[1], FES[2])
    assemble!(B, BilinearOperator([div(1)], [id(1)]; factor = -1))

    ## assemble right-hand side (f, v)
    b0 = FEVector(FES[1])
    assemble!(
        b0,
        LinearOperator(
            rhs!, [reconstruct ? apply(1, Reconstruct{HDIVRT0{2}, Identity}) : id(1)];
            bonus_quadorder = bonus_quadorder_f
        )
    )

    ## calculate bdofs for boundary conditions and future use
    bfacedofs = FES[1][BFaceDofs]
    nbfaces = num_sources(bfacedofs)
    bdofs = []
    for bface in 1:nbfaces
        append!(bdofs, view(bfacedofs, :, bface))
    end
    unique!(bdofs)

    ## set appropriate boundary on the right-hand side
    if !isnothing(exact_boundary!)
        interpolate!(b0[1], ON_BFACES, exact_boundary!; bonus_quadorder = bonus_quadorder_f)
    end

    ## get G matrix
    G = TB.G
    nmodes = num_multiindices(TB)

    ## solve
    @time (
        use_iterative_solver
            ? solve_stokes_primal!
            : solve_stokes_primal_full!
    )(sol, A0, A, B, b0, G, nmodes, bdofs)

    return bdofs
end
