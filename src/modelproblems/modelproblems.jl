abstract type AbstractModelProblem end

"""
$(TYPEDSIGNATURES)

solves the specified model problem with the given stochastic coefficient `C` and right-hand side `rhs`
and writes the solution into `sol`. Via this `sol` vector the spatial and stochastic discretization
is communicated as well as initial data for the iterative solver.

"""
function solve!(
    ::Type{AbstractModelProblem},
    sol::SGFEVector,                        ## target SGFEM vector
    C::AbstractStochasticCoefficient;
    kwargs...)
    @error "This Model problem seems to have no solver implemented!"
end

include("logpoisson_primal.jl")
include("logpoisson_dual.jl")
include("poisson_primal.jl")