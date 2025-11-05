#= 

([source code](SOURCE_URL))

runs AFEM loop for stochastic Poisson problem
 
usage:
- run experiment: run(; problem = problem, kwargs...)
- load results  : show_results(; kwargs...)
- produce plots : produce_plots(; kwargs...)

possible values for problem are
- PoissonProblemPrimal = Poisson problem with linear coefficient a
- LogTransformedPoissonProblemPrimal = log-transformed Poisson problem with exponential coefficient exp(a)
- LogTransformedPoissonProblemDual = dual formulation of the log-transformed Poisson problem

=#

module Poisson

using ExtendableASGFEM
using ExtendableFEM
using ExtendableFEMBase
using ExtendableGrids
using DataFrames
using DrWatson
using UnicodePlots
using LaTeXStrings
using GridVisualize
using CairoMakie
using DelimitedFiles
using Pkg

function filename(data; folder = "data", add = "", makepath = false)
    problem = data["problem"]
    decay = data["decay"]
    mean = data["mean"]
    order = data["order"]
    domain = data["domain"]
    use_equilibration_estimator = data["use_equilibration_estimator"]
    θ_stochastic = data["θ_stochastic"]
    θ_spatial = data["θ_spatial"]
    tail_extension = data["tail_extension"]
    maxdofs = data["maxdofs"]
    filename = "$(folder)/$(problem)/$(domain)/"
    if makepath
        mkpath(filename)
    end
    filename *= add * "order=$(order)_maxdofs=$(maxdofs)_decay=$(decay)_mean=$(mean)_θ=($(θ_spatial),$(θ_stochastic))_tail=$(tail_extension)"
    if use_equilibration_estimator
        filename *= "_eq"
    end
    return filename
end

default_args = Dict(
    "initial_refs" => 1,
    "domain" => "square", ## or lshape
    "order" => 1,
    "problem" => LogTransformedPoissonProblemPrimal,
    "C" => StochasticCoefficientCosinus,
    "decay" => 2,
    "mean" => 0,
    "maxm" => 150,
    "bonus_quadorder_a" => 2,
    "f" => (result, qpinfo) -> (result[1] = 1),
    "bonus_quadorder_f" => 0,
    "θ_stochastic" => 0.5,
    "θ_spatial" => 0.5,
    "factor_tail" => 1,
    "tail_extension" => [10, 2], # for [0] mode, and all others
    "maxdofs" => 1.0e4,
    "initial_modes" => [[0]],
    "nsamples" => 150,
    "use_equilibration_estimator" => false,
    "use_iterative_solver" => true,
    "Plotter" => nothing,
)


## that is the main function to run and save results
function run(; force = false, kwargs...)
    data = deepcopy(default_args)

    for (k, v) in kwargs
        data[String(k)] = v
    end
    Plotter = data["Plotter"]

    ## load/produce data
    data, ~ = produce_or_load(_main, data, filename = filename, force = force)

    ## plot final solution
    results = data["results"]
    sol = data["solution"]
    multi_indices = data["multi_indices"]
    repair_grid!(sol.FES_space[1].xgrid)
    if Plotter === nothing
        for j in 1:num_multiindices(sol)
            @show extrema(view(sol[j]))
            println(stdout, unicode_scalarplot(sol[j]; title = "MI $(multi_indices[j])"))
        end
    else
        plot_modes(sol; Plotter = Plotter, ncols = 4)
    end

    @info "final multi_indices = "
    for m in multi_indices
        println("$m")
    end
    @show results
    return data
end


# loads (or produces) a result file and prints the data
function show_results(; force = false, mode_history_up_to_level = 15, kwargs...)
    data = deepcopy(default_args)

    for (k, v) in kwargs
        data[String(k)] = v
    end

    data, ~ = produce_or_load(main, data, filename = filename, force = force)

    results = data["results"]
    nmodes = results[!, "nmodes"]
    nlevels = length(nmodes)

    ## multi_index statisticts
    MI = data["multi_indices"]
    istart = 1
    for j in 1:min(mode_history_up_to_level, nlevels)
        istart = j == 1 ? 1 : nmodes[j - 1] + 1
        iend = nmodes[j]
        if iend >= istart
            @info "NEW MODES on level $j : "
            maxlength = maximum([findlast(!=(0), MI[k]) for k in istart:iend])
            if maxlength === nothing
                maxlength = 1
            end
            for k in istart:iend
                @info MI[k][1:maxlength]
            end
        end
    end


    return data
end

# the main AFEM loop for the given problem and data
# (this function does not check for existing results, better use run !)
function _main(
        data = nothing;
        debug = false,
        plot_solution = false,
        Plotter = nothing,
        kwargs...
    )

    if isnothing(data)
        data = deepcopy(default_args)
    end

    for (k, v) in kwargs
        data[String(k)] = v
    end

    ## read arguments
    problem = data["problem"]
    nrefs = data["initial_refs"]
    decay = data["decay"]
    mean = data["mean"]
    τ = 1
    if problem <: PoissonProblemPrimal
        @assert mean >= 1 "coefficient for $problem needs to be at least 1 to ensure ellipticity"
        τ = 0.9
    end
    maxm = data["maxm"]
    C = data["C"](; τ = τ, decay = decay, mean = mean, maxm = maxm)
    θ_stochastic = data["θ_stochastic"]
    θ_spatial = data["θ_spatial"]
    multi_indices = Array{Array{Int, 1}, 1}(data["initial_modes"])
    nsamples = data["nsamples"]
    f! = data["f"]
    bonus_quadorder_a = data["bonus_quadorder_a"]
    bonus_quadorder_f = data["bonus_quadorder_f"]
    factor_tail = data["factor_tail"]
    maxdofs = data["maxdofs"]
    order = data["order"]
    use_iterative_solver = data["use_iterative_solver"]
    use_equilibration_estimator = data["use_equilibration_estimator"]

    prepare_multi_indices!(multi_indices)
    @info "initial multi_indices = $(multi_indices)"

    ############################
    ## spatial discretisation ##
    ############################
    dim = 2
    domain = data["domain"]
    if dim == 1
        xgrid = uniform_refine(reference_domain(Edge1D), nrefs)
    elseif dim == 2
        if domain == "square"
            xgrid = uniform_refine(grid_unitsquare(Triangle2D), nrefs)
        elseif domain == "lshape"
            xgrid = uniform_refine(grid_lshape(Triangle2D), nrefs)
        else
            @error "domain $domain not known"
        end
    else
        @error "dimension $dim not available yet"
    end


    df = DataFrame(
        ndofs_space = zeros(Int, 0),
        nmodes = zeros(Int, 0),
        exact_error_stress = zeros(Float64, 0),
        exact_error_u = zeros(Float64, 0),
        exact_error_stress_uni = zeros(Float64, 0),
        exact_error_u_uni = zeros(Float64, 0),
        estimate = zeros(Float64, 0),
        estimate_active = zeros(Float64, 0),
        estimate_tail = zeros(Float64, 0),
        time_solve = zeros(Float64, 0),
        time_estimate = zeros(Float64, 0)
    )

    sol = nothing
    lvl = 0
    while (true)

        push!(df, [0 0 0 0 0 0 0 0 0 0 0])

        ###############################
        ## stochastic discretisation ##
        ###############################
        lvl += 1
        @info "LEVEL $lvl"
        prepare_multi_indices!(multi_indices)
        M = maximum(length.(multi_indices))
        OBType = problem <: PoissonProblemPrimal ? LegendrePolynomials : HermitePolynomials
        ansatz_deg = maximum([maximum(multi_indices[k]) for k in 1:length(multi_indices)]) + 4
        TensorBasis = TensorizedBasis(OBType, M, ansatz_deg, 2 * ansatz_deg, 2 * ansatz_deg, multi_indices = multi_indices)
        if debug
            @info C, TensorBasis
            #plot_basis(TensorBasis; Plotter = Plotter)
        end

        if num_nodes(xgrid) < 200
            println(stdout, unicode_gridplot(xgrid))
        end

        ############
        ## solver ##
        ############
        ## generate solution object
        if problem <: LogTransformedPoissonProblemDual
            FEType = [HDIVRTk{2, order}, order == 0 ? L2P0{1} : H1Pk{1, dim, order}]
            FES = [FESpace{FEType[1]}(xgrid), FESpace{FEType[2]}(xgrid; broken = true)]
            unames = ["p", "u"]
        elseif problem <: LogTransformedPoissonProblemPrimal || problem <: PoissonProblemPrimal
            FEType = H1Pk{1, dim, order}
            FES = FESpace{FEType}(xgrid)
            unames = ["u"]
        else
            @error "problem type not available"
        end
        sol = SGFEVector(FES, TensorBasis; active_modes = 1:length(multi_indices), unames = unames)

        ## call solve
        time_solve = @elapsed bdofs = solve!(problem, sol, C; rhs = f!, bonus_quadorder_a = bonus_quadorder_a, bonus_quadorder_f = bonus_quadorder_f, use_iterative_solver = use_iterative_solver, debug = debug)
        df[lvl, :time_solve] = time_solve

        ############
        ## report ##
        ############
        ## plot modes
        if Plotter == nothing && plot_solution
            for j in 1:num_multiindices(sol)
                if debug
                    @show extrema(view(sol[j]))
                end
                println(stdout, unicode_scalarplot(sol[j]; title = "MI $(multi_indices[j])"))
            end
        elseif plot_solution
            plot_modes(sol; Plotter = Plotter, ncols = 4)
        end

        ############
        ## errors ##
        ############

        ## exact error by hierarchic sampling
        weightederrorH1, weightederrorL2, uniformerrorH1, uniformerrorL2 = calculate_sampling_error(sol, C; problem = problem, rhs = f!, order = order + 1, debug = debug, bonus_quadorder_a = bonus_quadorder_a, bonus_quadorder_f = bonus_quadorder_f, nsamples = nsamples)

        ## error estimator
        tail_extension = data["tail_extension"]
        if use_equilibration_estimator
            time_estimate = @elapsed η4modes, η4cells, multi_indices_extended = estimate_equilibration(problem, sol, C; rhs = f!, bonus_quadorder = max(bonus_quadorder_f, bonus_quadorder_a), tail_extension = tail_extension)
        else
            time_estimate = @elapsed η4modes, η4cells, multi_indices_extended = estimate(problem, sol, C; rhs = f!, bonus_quadorder = max(bonus_quadorder_f, bonus_quadorder_a), tail_extension = tail_extension)
        end
        for j in 1:length(multi_indices_extended)
            @info "mode = $(multi_indices_extended[j]) | error = $(η4modes[j])" # \t| f4mode = $(j <= length(multi_indices) ? f4modes[j] : 0)"
        end
        df[lvl, :time_estimate] = time_estimate

        ## compute norms of solution modes and grad(am)
        # nmodes = length(multi_indices)
        # normsuh4mode = zeros(Float64, nmodes)
        # M = FEMatrix(FES, FES)
        # assemble!(M, BilinearOperator([id(1)]))
        # for j = 1 : nmodes
        #     normsuh4mode[j] = lrmatmul(view(sol[j]), M.entries, view(sol[j]))
        # end
        # @show normsuh4mode

        # norms4gradam = zeros(Float64, maxm)
        # for m = 1 : maxm
        #     norms4gradam[m] = norm(integrate(xgrid, ON_CELLS, (result, qpinfo)->(get_gradam!(result, qpinfo.x, m, C); result.= result.^2;), 2))
        # end
        # @show norms4gradam

        # ## alternative tail estimator
        # for j = 1 : nmodes
        #     normsuh4mode[j] = lrmatmul(view(sol[j]), M.entries, view(sol[j]))
        #     tail_mode, tail_m = get_next_tail(multi_indices, multi_indices[j], 2, maxm)
        #     tail_estimate = f4modes[j] + normsuh4mode[j] * norms4gradam[tail_m]
        #     @show j, tail_mode, tail_m, tail_estimate
        # end

        #########################
        ## ADAPTIVE REFINEMENT ##
        #########################
        inactive_else, inactive_bnd, inactive_bnd2, active_bnd, active_int = classify_modes(multi_indices_extended, multi_indices_extended[1:TensorBasis.nmodes])

        sorted_id = sortperm(η4modes; rev = true)
        EstimatorInterior = 0
        EstimatorBoundary = 0
        actives = union(active_int, active_bnd)
        for j in 1:length(multi_indices_extended)
            print("$(multi_indices_extended[sorted_id[j]]) | ")
            if sorted_id[j] in inactive_else
                print(" [inactive far]")
            elseif sorted_id[j] in inactive_bnd
                print("[inactive bnd1]")
            elseif sorted_id[j] in inactive_bnd2
                print("[inactive bnd2]")
            elseif sorted_id[j] in active_int
                print("  [active int] ")
            elseif sorted_id[j] in active_bnd
                print("  [active bnd] ")
            end
            println(" | $(η4modes[sorted_id[j]])")
            if sorted_id[j] in actives
                EstimatorInterior += η4modes[sorted_id[j]]^2
            else
                EstimatorBoundary += η4modes[sorted_id[j]]^2
            end
        end
        @show actives

        ## save RESULTS
        @info "FINAL RESULTS
        || ∇(u-u_h) || (w,u) = $(sqrt.(weightederrorH1)), $(sqrt.(uniformerrorH1))
        estimator_total = $(sqrt(sum(η4modes .^ 2)))
        estimator_space = $(sqrt(EstimatorInterior))
        estimator_stoch = $(sqrt(EstimatorBoundary))
        || u - u_h || (w,u) = $(sqrt.(weightederrorL2)), $(sqrt.(uniformerrorL2))"

        df[lvl, :ndofs_space] = sum([sol.FES_space[j].ndofs for j in 1:length(sol.FES_space)])
        df[lvl, :nmodes] = sol.TB.nmodes
        df[lvl, :exact_error_stress] = sqrt(weightederrorH1[end])
        df[lvl, :exact_error_u] = sqrt(weightederrorL2[end])
        df[lvl, :exact_error_stress_uni] = sqrt(uniformerrorH1[end])
        df[lvl, :exact_error_u_uni] = sqrt(uniformerrorL2[end])
        df[lvl, :estimate] = sqrt(sum(η4modes .^ 2))
        df[lvl, :estimate_active] = sqrt(EstimatorInterior)
        df[lvl, :estimate_tail] = sqrt(EstimatorBoundary)
        data["results"] = df
        data["solution"] = sol
        data["multi_indices"] = multi_indices

        if length(sol.entries) >= maxdofs
            break
        end

        if EstimatorInterior > factor_tail * EstimatorBoundary
            println("...spatial refinement")
            if θ_spatial >= 1
                ## uniform mesh refinement
                xgrid = uniform_refine(xgrid)
            else
                ## refine by red-green-blue refinement (incl. closuring)
                facemarker = bulk_mark(xgrid, view(sum(view(η4cells, :, actives), dims = 2), :), θ_spatial; indicator_AT = ON_CELLS)
                # facemarker = bulk_mark(xgrid, view(η4cells,:,1), θ_spatial; indicator_AT = ON_CELLS)
                xgrid = RGB_refine(xgrid, facemarker)
            end
        else # stochastic refinement
            accumulated_stochastic_error = 0
            println("...stochastic refinement")
            for j in 1:length(multi_indices_extended)
                if (sorted_id[j] in inactive_bnd || sorted_id[j] in inactive_bnd2)
                    percentage = η4modes[sorted_id[j]]^2 / EstimatorBoundary
                    println("... adding mode $(multi_indices_extended[sorted_id[j]]) with error $(η4modes[sorted_id[j]]) (=$(Float16(percentage * 100))% of tail error)")
                    accumulated_stochastic_error += η4modes[sorted_id[j]]^2
                    push!(multi_indices, multi_indices_extended[sorted_id[j]])
                    if accumulated_stochastic_error >= θ_stochastic * EstimatorBoundary
                        break
                    end
                end
            end
        end

        @info "MODE-CLASSIFICATION
            active_interior = $(active_int)
            active_boundary = $(active_bnd)
            inactive_bnd = $(inactive_bnd)
            inactive_bnd2 = $(inactive_bnd2)
            inactive_else = $(inactive_else)"
    end

    filename_params = filename(data; add = "", folder = "data", makepath = true) * ".txt"

    data["version"] = pkgversion(ExtendableASGFEM)
    writedlm(filename_params, data, "=")

    return data
end


function produce_plots(;
        order = default_args["order"],
        decay = default_args["decay"],
        legend_position = :lb,
        xscale = log10,
        yscale = log10,
        force = false,
        scaling_factor = true,
        markersize = 10,
        plotM = 80,
        maxdegree = 12,
        template = :convergence, # :custom = specify your own cols
        cols = [:estimate, :exact_error_stress],
        xlabel = "ndofs",
        show_optimal_rate = false,
        kwargs...
    )

    colors = Makie.wong_colors()
    linestyles = [:solid, :dot, :dashdot, :dashdotdot]

    basedata = deepcopy(default_args)
    for (k, v) in kwargs
        basedata[String(k)] = v
    end

    if template == :convergence
        cols = [:estimate, :exact_error_stress]
    elseif template == :convergence_with_L2
        cols = [:estimate, :exact_error_stress, :exact_error]
    elseif template == :active_vs_tail
        cols = [:estimate_active, :estimate_tail]
    elseif template == :dofs
        cols = [:ndofs_space, :nmodes, :ndofs_all]
        xscale = nothing
        legend_position = :lt
        xlabel = "level"
    end


    ## load/produce data
    if typeof(order) <: Real
        order = [order]
    end
    if typeof(decay) <: Real
        decay = [decay]
    end
    data = Array{Dict{String, Any}, 2}(undef, length(order), length(decay))
    for o in 1:length(order), d in 1:length(decay)
        basedata["order"] = order[o]
        basedata["decay"] = decay[d]
        data[o, d], ~ = produce_or_load(main, basedata, filename = filename, force = force)
        xgrid = data[o, d]["solution"].FES_space[1].xgrid
        repair_grid!(xgrid)
        xgrid[BFaceRegions] .= 1

        ## gridplot
        gplt = GridVisualize.gridplot(xgrid; Plotter = CairoMakie, linewidth = 1, colorbar = :none, legend = :none)
        filename_plot = filename(basedata; add = "grid_", folder = "plots", makepath = true) * ".png"
        CairoMakie.save(filename_plot, gplt)
        @info "grid plot for $cols saved under $filename_plot"

    end

    ## load results
    for d in 1:length(decay)
        ## convergence plot
        ylabel = ""
        if length(cols) == 1
            if cols[1] == :exact_error_u
                ylabel = L"|| u - u_h ||_{L^2}"
            elseif cols[1] == :exact_error_stress
                ylabel = L"|| u - u_h ||_A"
            elseif cols[1] == :estimate
                ylabel = L"\eta"
            elseif cols[1] == :estimate_tail
                ylabel = L"\eta_\text{tail}"
            elseif cols[1] == :estimate_active
                ylabel = L"\eta_\text{active}"
            end
        end
        title = ""
        f_conv = Figure(fontsize = 18, size = (900, 600))
        if xscale !== nothing
            ax_conv = Axis(
                f_conv[1, 1], title = title, xlabel = xlabel, ylabel = ylabel, yscale = yscale, xscale = xscale,
                xminorticksvisible = true, xminorgridvisible = true, yminorgridvisible = true, yminorticksvisible = true,
                xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(10)
            )
        else
            ax_conv = Axis(
                f_conv[1, 1], title = title, xlabel = xlabel, ylabel = ylabel, yscale = yscale,
                xminorticksvisible = true, xminorgridvisible = true, yminorgridvisible = true, yminorticksvisible = true,
                xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(10)
            )
        end

        @info "producing convergence history..."
        for o in 1:length(order)
            df = data[o, d]["results"]
            nlevels = size(df, 1)
            ndofs = nothing
            ndofs_space = nothing
            ndofs_stochastic = nothing
            for c in 1:length(cols)
                color = colors[o]
                ndofs_space = df[!, :ndofs_space]
                ndofs_stochastic = df[!, :nmodes]
                if template == :dofs
                    ndofs = 1:length(ndofs_space)
                else
                    ndofs = df[!, :ndofs_space] .* df[!, :nmodes]
                end
                if cols[c] == :ndofs_all
                    plotdata = df[!, :ndofs_space] .* df[!, :nmodes]
                else
                    plotdata = df[!, cols[c]]
                end

                label = ""
                if length(cols) > 1
                    if cols[c] == :exact_error_u
                        label = L"|| u - u_h ||_{L^2}"
                    elseif cols[c] == :exact_error_stress
                        label = L"|| u - u_h ||_A"
                    elseif cols[c] == :estimate
                        if scaling_factor !== nothing
                            plotdata .*= scaling_factor
                        end
                        label = L"\eta"
                    elseif cols[c] == :estimate_tail
                        if scaling_factor !== nothing
                            plotdata .*= scaling_factor
                        end
                        label = L"\eta(\partial_h \Lambda)"
                    elseif cols[c] == :estimate_active
                        if scaling_factor !== nothing
                            plotdata .*= scaling_factor
                        end
                        label = L"\eta(\Lambda)"
                    elseif cols[c] == :ndofs_space
                        label = L"| V_h |"
                    elseif cols[c] == :nmodes
                        label = L"| \Lambda |"
                    elseif cols[c] == :ndofs_all
                        label = L"| \mathcal{V}_h |"
                    end
                end
                label = L"%$label order = %$(order[o])"
                scatterlines!(ndofs, plotdata, label = label, color = color, markersize = markersize, linewidth = 2, linestyle = linestyles[c])
            end
            if show_optimal_rate
                label = L"\mathcal{O}(N_\text{space}^{%$(-order[o]/2)})"
                scatterlines!(ndofs, ndofs .^ (-order[o] / 2), label = label, color = :gray, markersize = 0, linewidth = 2, linestyle = linestyles[length(cols) + 1])
            end
        end
        axislegend(ax_conv, position = legend_position, merge = true, labelsize = 16, orientation = :horizontal, nbanks = length(order) * length(cols) + 1)

        basedata["order"] = order
        basedata["decay"] = decay[d]
        filename_plot = filename(basedata; add = "$(template)_", folder = "plots", makepath = true) * ".png"
        CairoMakie.save(filename_plot, f_conv)
        @info "convergence plot for $cols saved under $filename_plot"
    end

    ## bar plot
    for o in 1:length(order)
        max4mode = zeros(length(decay), plotM)
        cat = zeros(Int, plotM * length(decay))
        height = zeros(Float64, plotM * length(decay))
        groups = zeros(Int, plotM * length(decay))
        for d in 1:length(decay)
            multi_indices = data[o, d]["multi_indices"]
            nmodes = length(multi_indices)
            M = min(maximum(length.(multi_indices)), plotM)
            max4mode[d, 1:M] = [maximum([multi_indices[j][k] for j in 1:nmodes]) for k in 1:M]
            cat[d:length(decay):end] .= 1:plotM
            groups[d:length(decay):end] .= d
            height[d:length(decay):end] .= max4mode[d, :]
        end


        # Plot
        replace!(height, 0 => 0.1)
        if maxdegree == 0
            #  maxdegree = maximum(max4mode[:])
        end
        yticks = (1:maxdegree, ["$(Int(j))" for j in 1:maxdegree])
        xticks = (1:plotM, [j % 2 == 1 ? "$(Int(j))" : "" for j in 1:plotM])
        f_bar = Figure(fontsize = 12, size = (900, 250))
        Axis(f_bar[1, 1], limits = (0, plotM, 0, maxdegree), yticks = yticks, xticks = xticks, xlabel = "stochastic dimension", ylabel = "max degree", xticklabelrotation = pi / 2)
        CairoMakie.barplot!(
            cat, height,
            dodge = AbstractVector{Integer}(groups),
            color = colors[groups]
        )

        # Legend
        labels = ["decay $d" for d in decay]
        elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
        title = ""

        Legend(f_bar[1, 2], elements, labels, title)

        basedata["decay"] = decay
        basedata["order"] = order[o]
        filename_plot = filename(basedata; add = "barplot_", folder = "plots", makepath = true) * ".png"
        CairoMakie.save(filename_plot, f_bar)
        @info "modes bar plot for $cols saved under $filename_plot"
    end


    return
end

function repair_grid!(xgrid::ExtendableGrid)
    xgrid.components[CellGeometries] = VectorOfConstants{ElementGeometries, Int}(xgrid.components[CellGeometries][1], num_cells(xgrid))
    xgrid.components[FaceGeometries] = VectorOfConstants{ElementGeometries, Int}(xgrid.components[FaceGeometries][1], length(xgrid.components[FaceGeometries]))
    xgrid.components[BFaceGeometries] = VectorOfConstants{ElementGeometries, Int}(xgrid.components[BFaceGeometries][1], length(xgrid.components[BFaceGeometries]))

    xgrid.components[UniqueCellGeometries] = Vector{ElementGeometries}([xgrid.components[CellGeometries][1]])
    xgrid.components[UniqueFaceGeometries] = Vector{ElementGeometries}([xgrid.components[FaceGeometries][1]])
    return xgrid.components[UniqueBFaceGeometries] = Vector{ElementGeometries}([xgrid.components[BFaceGeometries][1]])
end

end # module
