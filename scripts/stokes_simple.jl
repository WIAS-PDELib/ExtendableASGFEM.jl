module StokesSimple

using CairoMakie
using DataFrames
using DelimitedFiles
using ExtendableASGFEM
using ExtendableFEM
using ExtendableFEMBase
using ExtendableGrids
using GridVisualize
using LaTeXStrings
using PythonPlot
using Symbolics

## Structs for the description of numerical examples
struct DomainConfiguration
    type::String
    shift::Union{Vector{Float64}, Nothing}
    scale::Union{Vector{Float64}, Nothing}
end

struct StochasticConfiguration{T1 <: OrthogonalPolynomialType, T2 <: AbstractStochasticCoefficient}
    OBType::Type{T1}
    SCType::Type{T2}
end

struct RunConfiguration{
        TOBType <: OrthogonalPolynomialType, TSCType <: AbstractStochasticCoefficient, Tkwargs,
    }
    generate_data::Function
    domain_configuration::Union{DomainConfiguration, Nothing}
    stochastic_configuration::Union{StochasticConfiguration{TOBType, TSCType}, Nothing}
    initial_modes::Vector{Vector{Int}}
    runner::Function
    kwargs_list::Vector{Tkwargs}
end

## Helper methods for creating functions from symbolic expressions
function map_expression_to_eval_func(expr, variables)
    compiled = build_function(expr, variables..., expression = Val{false})
    return isa(compiled, Tuple) ? compiled[2] : (result, values...) -> (result[1] = compiled(values...))
end

function map_expression_to_deterministic_kernel(expr, variables)
    eval_func! = map_expression_to_eval_func(expr, variables)
    return (result, qpinfo) -> eval_func!(result, qpinfo.x...)
end

function map_expression_to_stochastic_kernel(
        expr,
        variables,
        sample_pointer_map = nothing
    )
    eval_func! = map_expression_to_eval_func(expr, variables)
    if !isnothing(sample_pointer_map)
        return function (sample_pointer)
            return function (result, qpinfo)
                eval_func!(result, sample_pointer_map(sample_pointer)..., qpinfo.x...)
                return nothing
            end
        end
    else # assume that expr encodes a deterministic function -> first closure ignores sample_pointer
        return (_) -> (result, qpinfo) -> eval_func!(result, qpinfo.x...)
    end
end

## Helper methods for the generation of the synthetic problem data and for the construction
## of data structures for the spatial and stochastic discretization
function prepare_data_for_example_with_analytical_solution(rc::RunConfiguration; kwargs...)
    exact_f!, exact_u!, exact_∇u!, exact_p! = rc.generate_data(; kwargs...)

    function error_configuration_closure(exact_func!::Function, oa_args, ids)
        return function (_...) # Ignore data_error_* args in sampling_error.jl
            function exact_func_closure(sample_pointer)
                exact_func_for_fixed_stochastic_input! = exact_func!(sample_pointer)
                return function (result, input, qpinfo)
                    exact_func_for_fixed_stochastic_input!(result, qpinfo)
                    result[1] = sum((result - input) .^ 2)
                    return nothing
                end
            end

            return exact_func_closure, oa_args, ids
        end
    end

    error_closures = map(
        error_configuration_closure,
        (exact_u!, exact_∇u!, exact_p!),
        ([id(1)], [grad(1)], [id(1)]),
        ([(2, 1)], [(2, 1)], [(2, 2)])
    )

    return exact_f!, exact_u!, error_closures...
end

function prepare_domain(rc::RunConfiguration, nrefs::Int)
    dc = rc.domain_configuration
    shift = !isnothing(dc.shift) ? dc.shift : [0, 0]
    scale = !isnothing(dc.scale) ? dc.scale : [1, 1]

    return if dc.type == "square"
        (
            uniform_refine(
                grid_unitsquare(Triangle2D; shift, scale), nrefs
            ), 1:4,
        )
    elseif dc.type == "lshape"
        (
            uniform_refine(
                grid_lshape(Triangle2D; shift, scale), nrefs
            ), 1:8,
        )
    else
        error("unknown domain type: $(dc.type)")
    end
end

function prepare_stochastic_data_structures(rc::RunConfiguration)
    sc = rc.stochastic_configuration
    multi_indices = Vector{Vector{Int}}(rc.initial_modes)
    prepare_multi_indices!(multi_indices)
    M = maximum(length.(multi_indices))
    ansatz_deg = maximum([maximum(multi_indices[k]) for k in 1:length(multi_indices)]) + 4
    TensorBasis = TensorizedBasis(
        sc.OBType, M, ansatz_deg, 2 * ansatz_deg, 2 * ansatz_deg, multi_indices = multi_indices
    )

    return multi_indices, TensorBasis
end

function create_plot(
        xdata, ydata; xscale = identity, yscale = identity, xlabel::LaTeXStrings.LaTeXString,
        ylabel::LaTeXStrings.LaTeXString, title::LaTeXStrings.LaTeXString, position = :rt,
        plot_path::String, labels::Vector{LaTeXStrings.LaTeXString}, show_optimal_rate = false,
    )
    colors = Makie.wong_colors()
    p = Makie.Figure(fontsize = 18, size = (900, 600))
    ax = Axis(
        p[1, 1], xscale = xscale, yscale = yscale, xlabel = xlabel,
        ylabel = ylabel, title = title, xminorticksvisible = true,
        xminorgridvisible = true, yminorgridvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(10)
    )
    for i in 1:length(ydata[1])
        scatterlines!(
            xdata, getindex.(ydata, i), color = colors[i], label = labels[i],
            markersize = 10, linewidth = 2, linestyle = :solid
        )
    end
    if show_optimal_rate
        linestyles = [:dot, :dashdot, :dashdotdot]
        for i in 1:3
            scatterlines!(
                xdata, xdata .^ (-i / 2), label = L"\mathcal{O}(h^{%$(-i)})", color = :gray,
                markersize = 0, linewidth = 2, linestyle = linestyles[i],
            )
        end
    end
    axislegend(
        ax, merge = true, labelsize = 16, orientation = :horizontal, position = position
    )
    mkpath(dirname(plot_path))
    CairoMakie.save(plot_path, p)

    return nothing
end

## Methods for generating synthetic problem data
function generate_data_for_simple_example(; kwargs...)
    ν = get(kwargs, :ν, [1.0, 0.2])

    @variables x1, x2, y1

    u = [0, (1 / 2) * (1 - x1^2)]
    ∇u = Symbolics.jacobian(u, [x1, x2])
    ∇u_reshaped = [∇u[1, 1], ∇u[1, 2], ∇u[2, 1], ∇u[2, 2]]
    p = - ν[1] * x2 - ν[2] * x2 * y1
    f = 0 * [x1, x2]

    return map_expression_to_deterministic_kernel(f, [x1, x2]),
        map(
            map_expression_to_stochastic_kernel,
            (u, ∇u_reshaped, p),
            ([x1, x2], [x1, x2], [y1, x1, x2]),
            (nothing, nothing, y -> [y[1]])
        )...
end

function generate_data_for_simple_example_with_modified_rhs(; kwargs...)
    ν = get(kwargs, :ν, [1.0, 0.2])

    @variables x1, x2, y1

    _, exact_u!, exact_∇u!, _ = generate_data_for_simple_example(; kwargs...)
    f = [x1, x2]
    # D = [-1, 1] × [-1, 1] and ∫_D p = 4/3 => subtract 1/3
    p = - ν[1] * x2 - ν[2] * x2 * y1 + 0.5 * (x1^2 + x2^2) - 1 / 3

    exact_f! = map_expression_to_deterministic_kernel(f, [x1, x2])
    exact_p! = map_expression_to_stochastic_kernel(p, [y1, x1, x2], y -> [y[1]])

    return exact_f!, exact_u!, exact_∇u!, exact_p!
end

function generate_data_for_complex_example(; kwargs...)
    ν = get(kwargs, :ν, [1.0, 0.2])

    @variables x1, x2

    ψ = x1^2 * (x1 - 1)^2 * x2^2 * (x2 - 1)^2
    ∇ψ = Symbolics.gradient(ψ, [x1, x2])
    curl_ψ = [-∇ψ[2], ∇ψ[1]]
    ∇u = Symbolics.jacobian(curl_ψ, [x1, x2])
    Δu = [
        Symbolics.derivative(∇u[1, 1], x1) + Symbolics.derivative(∇u[1, 2], x2),
        Symbolics.derivative(∇u[2, 1], x1) + Symbolics.derivative(∇u[2, 2], x2),
    ]
    p = x1^3
    ∇p = Symbolics.gradient(p, [x1, x2])
    f = - ν[1] * Δu + ∇p
    boundary_data = 0 * [x1, x2]

    return map_expression_to_deterministic_kernel(f, [x1, x2]),
        map_expression_to_stochastic_kernel(boundary_data, [x1, x2])
end

## Code containing the actual numerical examples
function runner_plot_modes_for_solution(rc::RunConfiguration; kwargs...)
    α = get(kwargs, :α, 1)
    ν = α * get(kwargs, :ν, [1.0, 0.2])
    reconstruct = get(kwargs, :reconstruct, false)
    FETypes = get(kwargs, :FETypes, (H1P2{2, 2}, H1P1{1}))
    nrefs = get(kwargs, :nrefs, 5)

    SCType = rc.stochastic_configuration.SCType
    C = SCType(ν)

    exact_f!, exact_boundary! = rc.generate_data(; ν = ν)

    xgrid, _ = prepare_domain(rc, nrefs)
    multi_indices, TensorBasis = prepare_stochastic_data_structures(rc)

    unames = ["u", "p"]
    FES = [FESpace{FETypes[1]}(xgrid), FESpace{FETypes[2]}(xgrid)]

    SolutionSGFEM = SGFEVector(FES, TensorBasis; active_modes = 1:length(multi_indices), unames)
    solve!(
        # Stochastic boundary data are not supported by the codebase. Thus, exact_u! is required to be
        # deterministic for this function to work. If exact_u! is a deterministic function, then the
        # boundary data can be obtained by evaluating exact_u! at any point in the stochastic domain.
        StokesProblemPrimal, SolutionSGFEM, exact_f!, C; (exact_boundary!) = (exact_boundary!([0, 0])),
        bonus_quadorder_f = 10, use_iterative_solver = true, reconstruct = reconstruct,
    )

    display(plot_modes(SolutionSGFEM; Plotter = PythonPlot, ncols = 4, normalizePattern = ["p"]))

    return SolutionSGFEM
end

function runner_error_for_taylor_hood_element(rc::RunConfiguration; kwargs...)
    exact_f!, exact_u!, error_u_closure, error_∇u_closure, error_p_closure =
        prepare_data_for_example_with_analytical_solution(rc)
    C = rc.stochastic_configuration.SCType()
    multi_indices, TensorBasis = prepare_stochastic_data_structures(rc)
    M::Int = maxlength_multiindices(TensorBasis)

    nsamples, Msamples = 100, maxm(C)
    Samples, weights = sample_distribution(
        TensorBasis, nsamples; M = Msamples, Mweights = Msamples
    )
    metrics_configuration = [
        Dict("name" => "L2u", "closure" => error_u_closure, "dim" => 2),
        Dict("name" => "L2stress", "closure" => error_∇u_closure, "dim" => 4),
        Dict("name" => "L2p", "closure" => error_p_closure, "dim" => 1),
    ]
    metrics = setup_metrics(StokesProblemPrimal, metrics_configuration)

    df = DataFrame(
        ndofs = zeros(Int, 0),
        ∇u_error = zeros(Float64, 0),
        u_error = zeros(Float64, 0),
        p_error = zeros(Float64, 0),
    )

    FETypes, unames = (H1P2{2, 2}, H1P1{1}), ["u", "p"]

    SolutionSGFEM = nothing
    for nrefs in 1:6
        xgrid, _ = prepare_domain(rc, nrefs)
        FES = [FESpace{FETypes[1]}(xgrid), FESpace{FETypes[2]}(xgrid)]

        SolutionSGFEM = SGFEVector(FES, TensorBasis; active_modes = 1:length(multi_indices), unames)
        solve!(
            # Stochastic boundary data are not supported by the codebase. Thus, exact_u! is required to be
            # deterministic for this function to work. If exact_u! is a deterministic function, then the
            # boundary data can be obtained by evaluating exact_u! at any point in the stochastic domain.
            StokesProblemPrimal, SolutionSGFEM, exact_f!, C; (exact_boundary!) = (exact_u!([0, 0])),
            bonus_quadorder_f = 2, use_iterative_solver = true, reconstruct = false,
        )

        u_error, ∇u_error, p_error = apply_error_integrators_from_metrics!(
            SolutionSGFEM, nothing, metrics, Samples, weights, M, Msamples; dimensionwise_error = false
        )

        push!(
            df, [(FES[1].ndofs + FES[2].ndofs) * length(multi_indices) ∇u_error[end] u_error[end] p_error[end]]
        )
    end

    if haskey(kwargs, :table_path)
        writedlm(kwargs[:table_path], Dict("result" => df))
    end

    return SolutionSGFEM
end

function runner_non_pressure_robustness_error(rc::RunConfiguration; kwargs...)
    nrefs = get(kwargs, :nrefs, 4)
    reconstruct = get(kwargs, :reconstruct, false)

    xgrid, _ = prepare_domain(rc, nrefs)
    SCType = rc.stochastic_configuration.SCType
    multi_indices, TensorBasis = prepare_stochastic_data_structures(rc)
    ν = [1.0, 0.2]

    FETypes, unames = (H1BR{2}, L2P0{1}), ["u", "p"]
    FES = [FESpace{FETypes[1]}(xgrid), FESpace{FETypes[2]}(xgrid)]

    SolutionSGFEM, metrics = nothing, nothing
    exponent_range, error = 1:-1:-9, Vector{Vector{Float64}}(undef, 0)
    for i in exponent_range
        α = 10.0^i
        C = SCType(α * ν)
        exact_f!, exact_u!, error_u_closure, error_∇u_closure, error_p_closure =
            prepare_data_for_example_with_analytical_solution(rc; ν = α * ν)

        SolutionSGFEM = SGFEVector(FES, TensorBasis; active_modes = 1:length(multi_indices), unames)
        solve!(
            # Stochastic boundary data are not supported by the codebase. Thus, exact_u! is required to be
            # deterministic for this function to work. If exact_u! is a deterministic function, then the
            # boundary data can be obtained by evaluating exact_u! at any point in the stochastic domain.
            StokesProblemPrimal, SolutionSGFEM, exact_f!, C; (exact_boundary!) = (exact_u!([0, 0])),
            bonus_quadorder_f = 2, use_iterative_solver = true, reconstruct = reconstruct,
        )

        nsamples, Msamples = 100, maxm(C)
        Samples, weights = sample_distribution(
            TensorBasis, nsamples; M = Msamples, Mweights = Msamples
        )
        metrics_configuration = [
            Dict("name" => "L2u", "closure" => error_u_closure, "dim" => 2),
            Dict("name" => "L2stress", "closure" => error_∇u_closure, "dim" => 4),
            Dict("name" => "L2p", "closure" => error_p_closure, "dim" => 1),
        ]
        metrics = setup_metrics(StokesProblemPrimal, metrics_configuration)

        u_error, ∇u_error, p_error = apply_error_integrators_from_metrics!(
            SolutionSGFEM, nothing, metrics, Samples, weights, 1, Msamples; dimensionwise_error = false
        )

        push!(error, [∇u_error[end]; u_error[end]; p_error[end]])
    end

    if haskey(kwargs, :plot_path)
        create_plot(
            10.0 .^ (exponent_range), error; xscale = log10, yscale = log10, xlabel = L"\alpha",
            ylabel = L"error measured in $L^2$-norm", plot_path = kwargs[:plot_path],
            title = L"""Error for $-\nabla\cdot((\alpha \cdot \nu)\nabla\mathbf{u}) + \nabla p = \mathbf{x}$,
            with Taylor-Hood pair $\mathcal{P}_2/\mathcal{P}_1$""",
            labels = [L"\|∇(u - u_h)\|_{L^2}", L"\|u - u_h\|_{L^2}", L"\|p - p_h\|_{L^2}"],
        )
    end

    return SolutionSGFEM
end

function runner_convergence_plot(rc::RunConfiguration; kwargs...)
    α = get(kwargs, :α, 1)
    ν = α * get(kwargs, :ν, [1.0, 0.2])
    reconstruct = get(kwargs, :reconstruct, false)
    FETypes = get(kwargs, :FETypes, (H1P2{2, 2}, H1P1{1}))
    nrefs_max = get(kwargs, :nrefs_max, 5)

    SCType = rc.stochastic_configuration.SCType
    C = SCType(ν)
    multi_indices, TensorBasis = prepare_stochastic_data_structures(rc)

    exact_f!, exact_boundary! = rc.generate_data(; ν = ν)

    unames = ["u", "p"]

    SolutionSGFEM = nothing
    ndofs, error = Vector{Float64}(undef, 0), Vector{Vector{Float64}}(undef, 0)
    for nrefs in 1:nrefs_max
        xgrid, _ = prepare_domain(rc, nrefs)
        FES = [FESpace{FETypes[1]}(xgrid), FESpace{FETypes[2]}(xgrid)]

        SolutionSGFEM = SGFEVector(FES, TensorBasis; active_modes = 1:length(multi_indices), unames)
        solve!(
            # Stochastic boundary data are not supported by the codebase. Thus, exact_u! is required to be
            # deterministic for this function to work. If exact_u! is a deterministic function, then the
            # boundary data can be obtained by evaluating exact_u! at any point in the stochastic domain.
            StokesProblemPrimal, SolutionSGFEM, exact_f!, C; (exact_boundary!) = (exact_boundary!([0, 0])),
            bonus_quadorder_f = 10, use_iterative_solver = true, reconstruct = reconstruct,
        )

        l2stress_error, l2u_error, l2p_error = calculate_sampling_error_2(
            SolutionSGFEM, exact_f!, C; problem = StokesProblemPrimal, nsamples = 5, order = 2,
            metrics_configuration = stokes_metrics_configuration, dimensionwise_error = true,
        )

        push!(ndofs, (FES[1].ndofs + FES[2].ndofs) * length(multi_indices))
        push!(error, [l2stress_error[end]; l2u_error[end]; l2p_error[end]])
    end

    if haskey(kwargs, :plot_path)
        create_plot(
            ndofs, error; xscale = log10, yscale = log10, xlabel = L"\text{number of degrees of freedom}",
            ylabel = L"error measured in $L^2$-norm", plot_path = kwargs[:plot_path], position = :rt,
            title = reconstruct ? L"\text{Convergence plot for the complex example with reconstruction operator}"
                : L"\text{Convergence plot for the complex example without reconstruction operator}",
            labels = [L"\|∇(u - u_h)\|_{L^2}", L"\|u - u_h\|_{L^2}", L"\|p - p_h\|_{L^2}"], show_optimal_rate = true,
        )
    end

    return SolutionSGFEM
end

## Declarative definition of numerical examples
domain_configuration_standard_square = DomainConfiguration("square", nothing, nothing)
domain_configuration_large_square = DomainConfiguration("square", [-0.5, -0.5], [2, 2])

stochastic_configuration_constant_coefficients = StochasticConfiguration(
    LegendrePolynomials, StochasticCoefficientConstants,
)

run_configuration_show_modes = RunConfiguration(
    generate_data_for_simple_example,
    domain_configuration_large_square,
    stochastic_configuration_constant_coefficients,
    [[0], [1]],
    runner_plot_modes_for_solution,
    [],
)

run_configuration_error_for_taylor_hood_element = RunConfiguration(
    generate_data_for_simple_example,
    domain_configuration_large_square,
    stochastic_configuration_constant_coefficients,
    [[0], [1]],
    runner_error_for_taylor_hood_element,
    [Dict(:table_path => "data/taylor_hood_error_for_simple_example.txt")],
)

run_configuration_non_pressure_robust_example = RunConfiguration(
    generate_data_for_simple_example_with_modified_rhs,
    domain_configuration_large_square,
    stochastic_configuration_constant_coefficients,
    [[0], [1]],
    runner_non_pressure_robustness_error,
    [Dict(:plot_path => "plots/non_pressure_robust_error_for_simple_example.svg")],
)

run_configuration_show_modes_complex_example = RunConfiguration(
    generate_data_for_complex_example,
    domain_configuration_standard_square,
    stochastic_configuration_constant_coefficients,
    [[0], [1], [2], [3]],
    runner_plot_modes_for_solution,
    [
        Dict(:α => 10.0^(-6), :FETypes => (H1BR{2}, L2P0{1}), :nrefs => 3, :reconstruct => false),
        Dict(:α => 10.0^(-6), :FETypes => (H1BR{2}, L2P0{1}), :nrefs => 3, :reconstruct => true),
    ]
)

run_configuration_convergence_plot_complex_example = RunConfiguration(
    generate_data_for_complex_example,
    domain_configuration_standard_square,
    stochastic_configuration_constant_coefficients,
    [[0], [1], [2], [3]],
    runner_convergence_plot,
    [
        Dict(
            :α => 10.0^(-4), :FETypes => (H1BR{2}, L2P0{1}), :reconstruct => false,
            :plot_path => "plots/convergence_plot_for_complex_example_non_pressure_robust.svg"
        ),
        Dict(
            :α => 10.0^(-4), :FETypes => (H1BR{2}, L2P0{1}), :reconstruct => true,
            :plot_path => "plots/convergence_plot_for_complex_example_pressure_robust.svg"
        ),
    ]
)

run_configurations = [
    run_configuration_show_modes,
    run_configuration_error_for_taylor_hood_element,
    run_configuration_non_pressure_robust_example,
    run_configuration_show_modes_complex_example,
    run_configuration_convergence_plot_complex_example,
]

## Helper methods for executing numerical examples, which are represented as RunConfiguration structs
function execute_runner(rc::RunConfiguration)
    return isempty(rc.kwargs_list) ?
        rc.runner(rc) :
        map(kwargs -> rc.runner(rc; kwargs...), rc.kwargs_list)
end

function execute_all()
    return map(execute_runner, run_configurations)
end

end
