abstract type AbstractModelProblem end

"""
$(TYPEDSIGNATURES)

solves the specified model problem with the given stochastic coefficient `C` and right-hand side `rhs`
and writes the solution into `sol`. Via this `sol` vector the spatial and stochastic discretization
is communicated as well as initial data for the iterative solver.
The boolean `use_iterative_solver` (default is true) determines if the iterative solver is used
or if the full matrix is assembled and solved by a direct solver (very slow for larger systems).
The parameters `bonus_quadorder_f` (default is 0) and `bonus_quadorder_a` (default is 2) can be used
to increase the quadrature order in terms that involve the rhs or the stochastic coefficient, respectively.


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
