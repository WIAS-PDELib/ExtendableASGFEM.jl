"""
$(TYPEDEF)

A stochastic coefficient is assumed to have the Karhunen-Loeve expansion form

    a(x,y) = a_0(x) + ∑_m y_m a_m(x)

with (centered independent) random variables `y_m` and
basis functions `a_m(x)` that need to be specified (together with their gradients)
and expectation value `a_0`
(in general they steam from a spectral analysis of the covariance operator
of a).


"""
abstract type AbstractStochasticCoefficient{T} end
include("cosinus.jl")


"""
$(TYPEDSIGNATURES)

prepares a function of interface

    a!(result, x, y)

that evaluates the coefficient (times a factor)
at space coordinates x and random variables y
into result (a vector of length 1).

"""
function get_a!(SC::AbstractStochasticCoefficient{T}; factor = 1) where {T}
    eval_am = zeros(T, 1)
    return function closure(result, x, y)
        result[1] = meanvalue(SC) # a[0]
        for m in 1:length(y)
            get_am!(eval_am, x, m, SC)
            result[1] += y[m] * eval_am[1]
        end
        result .*= factor
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

prepares a function of interface

    get_grada!(result, x, y)

that evaluates the spatial
gradient of the coefficient (times a factor)
at space coordinates x and random variables y
into result (a vector of same length as x).

"""
function get_grada!(SC::AbstractStochasticCoefficient{T}; factor = 1) where {T}
    eval_gradam = zeros(T, 2)
    return function closure(result, x, y)
        for m in 1:length(y)
            get_gradam!(eval_gradam, x, m, SC)
            result .+= y[m] * eval_gradam
        end
        result .*= factor
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

prepares a function of interface

    get_expa!(result, x, y)

that evaluates the exponential of
the coefficient (times a factor)
at space coordinates x and random variables y
into result (a vector of length 1).

"""
function get_expa!(SC::AbstractStochasticCoefficient; factor = 1)
    eval_a! = get_a!(SC)
    return function closure(result, x, y)
        eval_a!(result, x, y)
        result[1] = exp(factor * result[1])
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

returns a scalarplot of the m-th coefficient functiona `a_m`
interpolated on the given grid.
The Plotter backend can be changed with the `Plotter` argument
(e.g. GLMakie, CairoMakie, PyPlot, Plots).

"""
function plot_am(xgrid::ExtendableGrid, m, SC::AbstractStochasticCoefficient; Plotter = nothing, kwargs...)
    FES = FESpace{H1P1{1}}(xgrid)
    I = FEVector(FES)
    ExtendableFEMBase.interpolate!(I[1], (result, qpinfo) -> get_am!(result, qpinfo.x, m, SC))
    return scalarplot(xgrid, I.entries; Plotter, kwargs...)
end


"""
$(TYPEDSIGNATURES)

prepares a function of interface

    get_gradam_x_u!(result, input, qpinfo)

that evaluates `∇a_m(qpinfo.x) u(qpinfo.x)`
(used to define ExtendableFEM operator kernels
and input is expected to be some scalar quantity).

"""
function get_gradam_x_u(m, SC::AbstractStochasticCoefficient{T}) where {T}
    return function closure(result, input, qpinfo)
        get_gradam!(result, qpinfo.x, m, SC)
        result .*= input[1]
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

prepares a function of interface

    get_am_x!(result, input, qpinfo)

that evaluates `a_m(qpinfo.x) input`
(used to define ExtendableFEM operator kernels).

"""
function get_am_x(m, SC::AbstractStochasticCoefficient{T}) where {T}
    return function closure(result, input, qpinfo)
        get_am!(result, qpinfo.x, m, SC)
        result .= result[1] * input
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

prepares a function of interface

    get_gradam_x_sigma!(result, input, qpinfo)

that evaluates `∇a_m(qpinfo.x) ⋅ σ(qpinfo.x)`
(used to define ExtendableFEM operator kernels
and `σ` is expected to be some vector-valued
quantity of same length).

"""
## helper function that dot-products grad(a_m) with a vector valued input
function get_gradam_x_sigma(dim, m, SC::AbstractStochasticCoefficient{T}) where {T}
    ∇a::Array{T, 1} = zeros(T, dim)
    return function closure(result, input, qpinfo)
        get_gradam!(∇a, qpinfo.x, m, SC)
        result[1] = dot(∇a, input)
        return nothing
    end
end


"""
$(TYPEDSIGNATURES)

prepares a function of interface

    get_grada_x_sigma!(result, input, qpinfo)

that evaluates `∇a(qpinfo.x) ⋅ σ(qpinfo.x)`
(used to define ExtendableFEM operator kernels
and `σ` is expected to be some vector-valued
quantity of same length).

"""
function get_grada_x_sigma(dim, SC::AbstractStochasticCoefficient{T}, y) where {T}
    ∇a::Array{T, 1} = zeros(T, dim)
    eval_gradam::Array{T, 1} = zeros(T, dim)
    return function closure(result, input, qpinfo)
        for m in 1:length(y)
            get_gradam!(eval_gradam, qpinfo.x, m, SC)
            ∇a .+= y[m] * eval_gradam
        end
        result[1] = dot(∇a, input)
        return nothing
    end
end


"""
$(TYPEDSIGNATURES)

prepares two functions

    lambda_mu!(result, input, qpinfo)
    expa!(result, input, qpinfo)

that calculate `λ_μ`, which are orthogonal decomposition coefficient functions of `exp(a) ≈ ∑ λ_μ H_μ` w.r.t.
to the multi-indices `μ` and their associated orthogonal basis functions `H_μ` in TB, as well as
(an approximation based on this decomposition of) exp(a).

"""
function expa_PCE_mop(TB::TensorizedBasis{T, ONBType}, SC::AbstractStochasticCoefficient{T}; N_truncate::Int = maxm(SC), factor = 1) where {T, ONBType}
    eval_am::Array{T, 1} = zeros(T, 1)
    eval_amu_storage::T = 0.0
    eval_Hmu_storage::T = 0.0
    multi_indices = TB.multi_indices
    M = length(multi_indices[1])
    N = maximum([maximum(multi_indices[k]) for k in 1:length(multi_indices)])
    eval_aHmu_storage::T = 0.0
    mu_fac_storage::T = 0
    nmodes::Int = length(multi_indices)
    @debug "init PCE with M=$M and $nmodes modes"
    last_x = zeros(Float64, 2)

    function lambda_mu(result, x, μ) # for e^{a*factor}

        ## calculate exp(1/2 ∑a_m^2)
        sum_am = 0
        for m in 1:N_truncate
            get_am!(eval_am, x, m, SC) # gives a_m
            sum_am += eval_am[1]^2
        end
        sum_am = exp(sum_am / 2) * exp(meanvalue(SC) * factor)

        ## evaluate a_mu
        eval_amu::Float64 = eval_amu_storage
        mu_fac::Float64 = mu_fac_storage
        eval_amu = 1.0
        mu_fac = 1
        for d in 1:length(multi_indices[μ])
            get_am!(eval_am, x, d, SC)
            eval_amu *= eval_am[1]^multi_indices[μ][d]
            mu_fac *= factorial(multi_indices[μ][d])
        end
        eval_amu /= sqrt(mu_fac) * factor^sum(multi_indices[μ])

        result[1] = eval_amu * sum_am
        return nothing
    end

    function closure(result, x, y)
        eval_aHmu::Float64 = eval_aHmu_storage
        eval_Hmu::Float64 = eval_Hmu_storage
        eval_aHmu = 0.0

        ## evaluate all op basis polynomials once
        set_sample!(TB, y)

        for μ in 1:nmodes
            ## evaluate H_mu
            eval_Hmu = evaluate(TB, μ)

            ## evaluate lambda_mu := (e^a, H_μ) (coefficient for H_μ)
            lambda_mu(result, x, μ)
            #@show multi_indices[μ], result[1]

            eval_aHmu += result[1] * eval_Hmu
        end

        result[1] = eval_aHmu
        return nothing
    end

    return closure, lambda_mu
end


###### testing area
struct SingleStochasticCoefficient{T, k} <: AbstractStochasticCoefficient{T}
    coefficient::T
end

meanvalue(SC::SingleStochasticCoefficient) = 0
maxm(SC::SingleStochasticCoefficient{T, k}) where {T, k} = k

function get_am!(result, x, m, SC::SingleStochasticCoefficient{T, k}) where {T, k}
    return if m == k
        result[1] = SC.coefficient * x[1]
    else
        result[1] = 0
    end
end

function get_gradam!(result, x, m, SC::SingleStochasticCoefficient{T, k}) where {T, k}
    return if m == k
        result[1] = SC.coefficient
        result[2] = 0
    else
        fill!(result, 0)
    end
end
