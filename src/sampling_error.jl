struct Metric
    name::String
    ids::Vector{Tuple{Int, Int}}
    integrator::ExtendableFEM.ItemIntegrator
    solution::Vector{FEVectorBlock}
end

struct MetricConfiguration
    name::String
    closure::Function
    resultdim::Union{Int, Nothing}
    kwargs::Union{Dict{Symbol, Any}, Nothing}
end

function create_squared_l2_dist_closure!(; kwargs...)
    inputdim = get(kwargs, :inputdim, nothing)
    if !isnothing(inputdim) && isodd(inputdim)
        throw(DomainError("inputdim must be even for l2 error calculation"))
    end

    return function (result, input, _)
        N = (!isnothing(inputdim) ? inputdim : length(input)) >> 1
        result[1] = sum((input[1:N] - input[(N + 1):(N << 1)]) .^ 2)
        return nothing
    end
end

## kernel for L2 error of stress, i.e || grad(u) - grad(u_h) ||
function data_error_stress(
        ::Union{Type{<:PoissonProblemPrimal}, Type{<:LogTransformedPoissonProblemPrimal}},
        C::AbstractStochasticCoefficient, sample_pointer; kwargs...
    )
    return create_squared_l2_dist_closure!(; kwargs...), [grad(1), grad(2)], [(1, 1), (2, 1)]
end

function data_error_stress(
        ::Type{<:LogTransformedPoissonProblemDual},
        C::AbstractStochasticCoefficient, sample_pointer; kwargs...
    )
    return create_squared_l2_dist_closure!(; kwargs...), [id(1), id(2)], [(1, 1), (2, 1)]
end

## kernel for L2 error || u - u_h ||
function data_error_u(
        ::Union{Type{<:PoissonProblemPrimal}, Type{<:LogTransformedPoissonProblemPrimal}},
        C::AbstractStochasticCoefficient, sample_pointer; kwargs...
    )
    return create_squared_l2_dist_closure!(; kwargs...), [id(1), id(2)], [(1, 1), (2, 1)]
end
function data_error_u(
        ::Type{<:LogTransformedPoissonProblemDual},
        C::AbstractStochasticCoefficient, sample_pointer; kwargs...
    )
    return create_squared_l2_dist_closure!(; kwargs...), [id(1), id(2)], [(1, 2), (2, 2)]
end

u_with_stress_metric_configuration::Vector{MetricConfiguration} = [
    MetricConfiguration("L2stress", data_error_stress, 2, nothing),
    MetricConfiguration("L2u", data_error_u, 1, nothing),
]

FES4sampling(::Type{PoissonProblemPrimal}, dim, xgrid, order) = [FESpace{H1Pk{1, dim, order}}(xgrid)]
FES4sampling(::Type{LogTransformedPoissonProblemPrimal}, dim, xgrid, order) = [FESpace{H1Pk{1, dim, order}}(xgrid)]
FES4sampling(::Type{LogTransformedPoissonProblemDual}, dim, xgrid, order) = [
    FESpace{HDIVRTk{dim, order}}(xgrid), FESpace{H1Pk{1, dim, order}}(xgrid; broken = true),
]

function setup_metrics(
        problem,
        metrics_configuration::Vector{MetricConfiguration},
        C::AbstractStochasticCoefficient,
        sample_pointer;
        order::Int = 2
    )
    ## collect metrics from metrics_configuration in a single list
    metrics = []
    for configuration in metrics_configuration
        kwargs = isnothing(configuration.kwargs) ?
            Dict() :
            configuration.kwargs
        resultdim = isnothing(configuration.resultdim) ? 1 : configuration.resultdim

        kernel_closure!, input, ids = configuration.closure(
            problem, C, sample_pointer; kwargs...
        )

        push!(
            metrics,
            Metric(
                configuration.name,
                ids,
                ItemIntegrator(
                    kernel_closure!, input; resultdim, quadorder = 2 * order + 1
                ),
                Vector{FEVectorBlock}(undef, 2),
            ),
        )
    end

    return metrics
end

function sample_reference_solution(
        problem,
        FESSampling,
        Samples,
        C::AbstractStochasticCoefficient;
        rhs = nothing,
        bonus_quadorder_a = 2,
        bonus_quadorder_f = 0,
        nsamples = 100,
        parallel_sampling = true,
        debug = false,
    )
    ## prepare array with deterministic solutions
    sol_det = Vector{FEVector}(undef, nsamples)

    nthreads = parallel_sampling ? Threads.nthreads() : 1
    @info "Estimating exact error by MC sampling (with nthreads = $nthreads)"

    ## compute deterministic solutions (in parallel)
    Threads.@threads for s in 1:nsamples
        ## deterministic problem description
        PD, _ = deterministic_problem(
            problem,
            C,
            Samples[:, s];
            rhs = rhs,
            bonus_quadorder_a = bonus_quadorder_a,
            bonus_quadorder_f = bonus_quadorder_f,
        )

        sol_det[s] = ExtendableFEM.solve(PD, FESSampling; verbosity = debug ? 0 : -1, timeroutputs = :none)
        print(".")
    end
    println(" MC samples solved")

    return sol_det
end

function apply_error_integrators_from_metrics!(
        SolutionSGFEM::SGFEVector,
        sol_det,
        metrics,
        Samples,
        weights,
        M,
        Msamples;
        dimensionwise_error = true,
    )
    sol_sgfem = SolutionSGFEM.FEV
    M0 = dimensionwise_error ? 0 : M
    nsamples = size(Samples, 2)
    metric_names = [metric.name for metric in metrics]

    error_per_run = Dict()
    for metric_name in metric_names
        error_per_run[metric_name] = zeros(Float64, M + 1, nsamples)
    end

    nsamples = size(Samples, 2)
    for s in 1:nsamples
        for metric in metrics
            for (i, entry) in enumerate(metric.ids)
                if entry[1] == 1
                    metric.solution[i] = sol_det[s][entry[2]]
                else
                    metric.solution[i] = sol_sgfem[entry[2]]
                end
            end
        end

        csample = Samples[:, s]
        for m in M0:M
            y = view(csample, 1:m)
            set_sample!(SolutionSGFEM, y)
            for metric in metrics
                sol = ExtendableFEM.evaluate(metric.integrator, metric.solution)
                error_per_run[metric.name][m + 1, s] = sum(view(sol, 1, :))
            end
        end

        @info "SAMPLE $s of $nsamples
            sample = $(view(csample, 1:min(Msamples, 2 * M)))
            weight = $(weights[s])
            $(join(map(name -> "error($name) = $(error_per_run[name][end, s])", metric_names), "\n            "))"
    end

    error = Dict()
    weightsum = sum(weights)
    for name in metric_names
        name_weighted = "$(name)_weighted"
        name_uniform = "$(name)_uniform"
        error[name_weighted] = zeros(Float64, M + 1)
        error[name_uniform] = zeros(Float64, M + 1)
        for s in 1:nsamples
            slice = view(error_per_run[name], :, s)
            error[name_uniform] .+= slice
            error[name_weighted] .+= slice * weights[s]
        end
        for m in 1:(M + 1)
            error[name_uniform][m] /= nsamples
            error[name_weighted][m] /= weightsum
        end
    end

    return [error["$(name)_weighted"] for name in metric_names]...,
        [error["$(name)_uniform"] for name in metric_names]...
end


"""
$(TYPEDSIGNATURES)

Estimates the error for the given model problem by Monte Carlo sampling. For each sample, a deterministic finite element solution is computed for a fixed (sampled) coefficient with polynomial order `order` and compared to the provided stochastic Galerkin solution `SolutionSGFEM`.

# Returns
A tuple of 2n many arrays (each of length `M+1`, where `M` is the maximum order of the multi-indices) containing the
results of the n many configured metrics. The first n arrays are weighted by the distribution, while the latter n
arrays are weighted uniformly. The order of the calculated errors coincides with the order of the argument `metrics_configuration`.

The m-th component of these arrays gives the error when only multi-indices up to order m are included.

# Arguments
- `SolutionSGFEM`: The stochastic Galerkin solution (SGFEVector).
- `C`: The stochastic coefficient.
- `problem`: The model problem type (default: `LogTransformedPoissonProblemPrimal`).
- `metrics_configuration`: The metrics to be calculated (default: `u_with_stress_metric_configuration`).
- `bonus_quadorder_a`: Additional quadrature order for the coefficient (default: 10).
- `bonus_quadorder_f`: Additional quadrature order for the right-hand side (default: 0).
- `order`: Polynomial order for the deterministic reference solution (default: 2).
- `rhs`: Right-hand side function (optional).
- `dim`: Spatial dimension (default: inferred from the solution).
- `Msamples`: Number of random variables to sample (default: `maxm(C)`).
- `parallel_sampling`: Whether to use parallel sampling (default: true).
- `dimensionwise_error`: If true, computes errors dimensionwise (default: false).
- `debug`: If true, enables debug output and plotting (default: false).
- `nsamples`: Number of Monte Carlo samples (default: 100).

# Details
This function computes deterministic reference solutions for each sample, then compares them to the SGFEM solution to estimate the error. Both weighted and unweighted (uniform) averages are returned for stress and solution errors.
"""
function calculate_sampling_error(
        SolutionSGFEM::SGFEVector,
        C::AbstractStochasticCoefficient;
        problem = LogTransformedPoissonProblemPrimal,
        metrics_configuration = u_with_stress_metric_configuration,
        bonus_quadorder_a = 10,
        bonus_quadorder_f = 0,
        order = 2,
        rhs = nothing,
        dim = size(SolutionSGFEM.FES_space[1].xgrid[Coordinates], 1),
        Msamples = maxm(C),
        parallel_sampling = true,
        dimensionwise_error = true,
        debug = false,
        nsamples = 100
    )
    FES = SolutionSGFEM.FES_space
    xgrid = FES[1].xgrid
    TensorBasis = SolutionSGFEM.TB

    ## generate samples
    Samples, weights = sample_distribution(TensorBasis, nsamples; M = Msamples, Mweights = Msamples)

    ## prepare array with deterministic solutions
    ## calculate MC samples
    FES4Sampling = FES4sampling(problem, dim, xgrid, order)
    sol_det = sample_reference_solution(
        problem, FES4Sampling, Samples, C;
        rhs, bonus_quadorder_a, bonus_quadorder_f, nsamples, debug, parallel_sampling
    )

    ## create closure functions for each specified metric
    metrics = setup_metrics(problem, metrics_configuration, C, Samples; order)

    ## Calculate error for each metric
    M = maxlength_multiindices(TensorBasis)
    return apply_error_integrators_from_metrics!(
        SolutionSGFEM, sol_det, metrics, Samples, weights, M, Msamples; dimensionwise_error
    )
end
