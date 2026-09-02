"""
$(TYPEDEF)

"""
struct StochasticCoefficientConstants{T} <: AbstractStochasticCoefficient{T}
    constants::Vector{T}
end

maxm(SC::StochasticCoefficientConstants) = length(SC.constants)
meanvalue(SC::StochasticCoefficientConstants) = SC.constants[1]


"""
$(TYPEDSIGNATURES)

constructor for StochasticCoefficientConstants of type `T` (default = Float64)

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
