abstract type AbstractModelProblem end

"""
$(TYPEDSIGNATURES)

Solves the specified model problem using the provided stochastic coefficient `C` and right-hand side `rhs`, writing the solution into `sol`.

The `sol` vector communicates both the spatial and stochastic discretization, as well as any initial data required for the iterative solver.

- If `use_iterative_solver` (default: `true`) is set, an iterative solver is used. Otherwise, the full system matrix is assembled and solved directly (note: this is very slow for large systems).
- The parameters `bonus_quadorder_f` (default: `0`) and `bonus_quadorder_a` (default: `2`) allow you to increase the quadrature order for terms involving the right-hand side or the stochastic coefficient, respectively.
- Additional keyword arguments can be passed via `kwargs`.

If no solver is implemented for the given model problem, an error is thrown.
"""
function solve!(
        ::Type{AbstractModelProblem},
        sol::SGFEVector,                        ## target SGFEM vector
        C::AbstractStochasticCoefficient;
        rhs = nothing,
        use_iterative_solver = true,
        bonus_quadorder_f = 0,
        bonus_quadorder_a = 0,
        kwargs...
    )
    return @error "This Model problem seems to have no solver implemented!"
end

include("logpoisson_primal.jl")
include("logpoisson_dual.jl")
include("poisson_primal.jl")
