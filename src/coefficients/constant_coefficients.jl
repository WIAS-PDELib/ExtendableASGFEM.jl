"""
$(TYPEDEF)

Expansion of the form

    a(x,y) = a_0(x) + ∑_m y_m a_m(x)

where the `a_m` are of cosinus type.

"""
struct StochasticCoefficientConstants{T} <: AbstractStochasticCoefficient{T}
    constants::Vector{T}
end

maxm(SC::StochasticCoefficientConstants) = length(SC.constants)
meanvalue(SC::StochasticCoefficientConstants) = SC.constants[1]


"""
$(TYPEDSIGNATURES)

constructor for StochasticCoefficientConstants of type `T` (default = Float64),
where `decay` (default = 2) steers the decay of the coefficient basis functions (the larger the faster),
`mean` (default = 0) is the mean value of the coefficient, `maxm` (default = 100) is the maximal number of stochastic random variables,
and `τ` (default = 1) is a uniform scaling factors

"""
function StochasticCoefficientConstants(; constants = [1.0, 0.2])
    return StochasticCoefficientConstants{eltype(constants)}(constants)
end


function get_am!(result, x, m, SC::StochasticCoefficientConstants)
    result[1] = m < length(SC.constants) ? SC.constants[m + 1] : 0
    return nothing
end

function get_gradam!(result, x, m, SC::StochasticCoefficientConstants)
    fill!(result, 0)
    return nothing
end

function Base.show(io::IO, SC::StochasticCoefficientConstants)
    println(io, "COEFFICIENT DATA")
    return println(io, "constants = $(SC.constants)")
end
