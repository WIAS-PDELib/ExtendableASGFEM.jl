
"""
$(TYPEDEF)

Expansion of the form

    a(x,y) = a_0(x) + ∑_m y_m a_m(x)

where the `a_m` are of cosinus type.

"""
struct StochasticCoefficientCosinus{T} <: AbstractStochasticCoefficient{T}
    decay::T
    mean_value::T
    decay_factors::Array{T,1}
    b1::Array{Int,1}
    b2::Array{Int,1}
end

maxm(SC::StochasticCoefficientCosinus) = length(SC.decay_factors)
meanvalue(SC::StochasticCoefficientCosinus) = SC.mean_value


"""
$(TYPEDSIGNATURES)

constructor for StochasticCoefficientCosinus of type `T` (default = Float64),
where `decay` (default = 2) steers the decay of the coefficient basis functions (the larger the faster),
`mean` (default = 0) is the mean value of the coefficient, `maxm` (default = 100) is the maximal number of stochastic random variables,
and `τ` (default = 1) is a uniform scaling factors

"""
function StochasticCoefficientCosinus(; T = Float64, τ = 1, start = 2, decay = 2, mean = 0, maxm = 100)
    decay_factors = zeros(T, maxm)
    b1 = zeros(Int, maxm)
    b2 = zeros(Int, maxm)
    k = 0
    j = 0
    for m = 1 : maxm+1
        if m > 1
            decay_factors[m-1] = Float64(m-2+start)^-decay
            b1[m-1] = j #floor((m+2)/2)-1
            b2[m-1] = k #ceil((m+2)/2)-1
        end
        if k > 0
            j += 1
            k -= 1
        else
            k = (j+k)+1
            j = 0
        end 
    end
    amp = τ/zeta(decay, start) # Hurwitz zeta function (limit value of geometric sequence)
    decay_factors .*= amp
    return StochasticCoefficientCosinus{T}(decay, mean, decay_factors, b1, b2)
end


function get_am!(result, x, m, SC::StochasticCoefficientCosinus)
    if m == 0
        result[1] = SC.mean_value
    else
        result[1] = SC.decay_factors[m]*cos(π*SC.b1[m]*x[1])*cos(π*SC.b2[m]*x[2])
    end
    return nothing
end

function get_gradam!(result, x, m, SC::StochasticCoefficientCosinus)
    if m == 0
        fill!(result,0)
    else
        result[1] = -SC.b1[m]*π*sin(SC.b1[m]*π*x[1])*cos(SC.b2[m]*π*x[2])
        result[2] = -SC.b2[m]*π*cos(SC.b1[m]*π*x[1])*sin(SC.b2[m]*π*x[2])
        result .*= SC.decay_factors[m]
    end
    return nothing
end

function Base.show(io::IO, SC::StochasticCoefficientCosinus)
    println(io, "COEFFICIENT DATA")
    println(io, "decay = $(SC.decay)")
    println(io, "mean_value = $(SC.mean_value)")
    println(io, "")
    for j = 1 : length(SC.decay_factors)
        println(io, "$(SC.b1[j]) | $(SC.b2[j]) | $(SC.decay_factors[j])")
    end
    println(io, "sum(coeffs) = $(sum(SC.decay_factors))")
end
