## kernel for L2 error of stress, i.e || grad(u) - grad(u_h) ||
function data_error_stress(::Union{Type{<:PoissonProblemPrimal}, Type{<:LogTransformedPoissonProblemPrimal}}, dim, C::AbstractStochasticCoefficient, sample_pointer)
    function closure(result, input, qpinfo)
        result[1] = 1
        if dim == 1
            result[2] = (input[1] - input[2])^2
        elseif dim == 2
            result[2] = (input[1] - input[3])^2 + (input[2] - input[4])^2
        end
        result[1] *= result[2]
        return nothing
    end
    return closure, [grad(1), grad(2)], [(1, 1), (2, 1)]
end

function data_error_stress(::Type{<:LogTransformedPoissonProblemDual}, dim, C::AbstractStochasticCoefficient, sample_pointer)
    function closure(result, input, qpinfo)
        result[1] = 1
        if dim == 1
            result[2] = (input[1] - input[2])^2
        elseif dim == 2
            result[2] = (input[1] - input[3])^2 + (input[2] - input[4])^2
        end
        result[1] *= result[2]
        return nothing
    end
    return closure, [id(1), id(2)], [(1, 1), (2, 1)]
end

## kernel for L2 error || u - u_h ||
function data_error_u(::Union{Type{<:PoissonProblemPrimal}, Type{<:LogTransformedPoissonProblemPrimal}}, dim, C::AbstractStochasticCoefficient, sample_pointer)
    function closure(result, input, qpinfo)
        result[1] = (input[1] - input[2])^2
        return nothing
    end
    return closure, [id(1), id(2)], [(1, 1), (2, 1)]
end
function data_error_u(::Type{<:LogTransformedPoissonProblemDual}, dim, C::AbstractStochasticCoefficient, sample_pointer)
    function closure(result, input, qpinfo)
        result[1] = (input[1] - input[2])^2
        return nothing
    end
    return closure, [id(1), id(2)], [(1, 2), (2, 2)]
end

FES4sampling(::Type{PoissonProblemPrimal}, dim, xgrid, order) = [FESpace{H1Pk{1, dim, order}}(xgrid)]
FES4sampling(::Type{LogTransformedPoissonProblemPrimal}, dim, xgrid, order) = [FESpace{H1Pk{1, dim, order}}(xgrid)]
FES4sampling(::Type{LogTransformedPoissonProblemDual}, dim, xgrid, order) = [FESpace{HDIVRTk{dim, order}}(xgrid), FESpace{H1Pk{1, dim, order}}(xgrid; broken = true)]


"""
$(TYPEDSIGNATURES)

Estimates the error for the given model problem by Monte Carlo sampling. For each sample, a deterministic finite element solution is computed for a fixed (sampled) coefficient with polynomial order `order` and compared to the provided stochastic Galerkin solution `SolutionSGFEM`.

# Returns
A tuple of arrays (each of length `M+1`, where `M` is the maximum order of the multi-indices):
- `totalerrorL2stress_weighted`: Mean L2 error of the stress, weighted by the distribution.
- `totalerrorL2u_weighted`: Mean L2 error of the solution, weighted by the distribution.
- `totalerrorL2stress_uniform`: Mean L2 error of the stress, weighted uniformly.
- `totalerrorL2u_uniform`: Mean L2 error of the solution, weighted uniformly.

The m-th component of these arrays gives the error when only multi-indices up to order m are included.

# Arguments
- `SolutionSGFEM`: The stochastic Galerkin solution (SGFEVector).
- `C`: The stochastic coefficient.
- `problem`: The model problem type (default: `LogTransformedPoissonProblemPrimal`).
- `bonus_quadorder_a`: Additional quadrature order for the coefficient (default: 10).
- `bonus_quadorder_f`: Additional quadrature order for the right-hand side (default: 0).
- `order`: Polynomial order for the deterministic reference solution (default: 2).
- `rhs`: Right-hand side function (optional).
- `dim`: Spatial dimension (default: inferred from the solution).
- `Msamples`: Number of random variables to sample (default: `maxm(C)`).
- `parallel_sampling`: Whether to use parallel sampling (default: true).
- `dimensionwise_error`: If true, computes errors dimensionwise (default: false).
- `energy_norm`: If true, uses the energy norm for stress error (default: true).
- `debug`: If true, enables debug output and plotting (default: false).
- `nsamples`: Number of Monte Carlo samples (default: 100).

# Details
This function computes deterministic reference solutions for each sample, then compares them to the SGFEM solution to estimate the error. Both weighted and unweighted (uniform) averages are returned for stress and solution errors.
"""
function calculate_sampling_error(
        SolutionSGFEM::SGFEVector,
        C::AbstractStochasticCoefficient;
        problem = LogTransformedPoissonProblemPrimal,
        bonus_quadorder_a = 10,
        bonus_quadorder_f = 0,
        order = 2,
        rhs = nothing,
        dim = size(SolutionSGFEM.FES_space[1].xgrid[Coordinates], 1),
        Msamples = maxm(C),
        parallel_sampling = true,
        dimensionwise_error = true,
        energy_norm = true,
        debug = false,
        nsamples = 100
    )

    nthreads = Threads.nthreads()
    @info "Estimating exact error by MC sampling (with nthreads = $nthreads)"
    FES = SolutionSGFEM.FES_space
    xgrid = FES[1].xgrid
    sol_sgfem = SolutionSGFEM.FEV
    TensorBasis = SolutionSGFEM.TB
    multi_indices::Array{Array{Int, 1}, 1} = TensorBasis.multi_indices
    M::Int = maxlength_multiindices(TensorBasis)
    nmodes = num_multiindices(TensorBasis)

    ## generate samples
    Samples, weights = sample_distribution(TensorBasis, nsamples; M = Msamples, Mweights = Msamples)

    ## prepare array with deterministic solutions
    sol_det = Array{FEVector, 1}(undef, nsamples)

    ## compute deterministic solutions (in parallel)
    Threads.@threads for s in 1:nsamples

        ## deterministic problem description
        PD, u = deterministic_problem(problem, C, Samples[:, s]; rhs = rhs, bonus_quadorder_a = bonus_quadorder_a, bonus_quadorder_f = bonus_quadorder_f)

        ## solve problem for the current sample
        FESSampling = FES4sampling(problem, dim, xgrid, order)
        sol_det[s] = ExtendableFEM.solve(PD, FESSampling; verbosity = debug ? 0 : -1, timeroutputs = :none)

        print(".")
    end
    println(" MC samples solved")

    ## compute errors for each sample (sequentially)

    ## construct error estimation kernels
    csample = Samples[:, 1]
    kernel_stress!, input_stress, ids_stress = data_error_stress(problem, dim, C, csample)
    kernel_u!, input_u, ids_u = data_error_u(problem, dim, C, csample)
    ErrorIntegratorL2stress = ItemIntegrator(kernel_stress!, input_stress; resultdim = dim, quadorder = 2 * order)
    ErrorIntegratorL2u = ItemIntegrator(kernel_u!, input_u; quadorder = 2 * order)


    errorL2stress = zeros(Float64, M + 1, nsamples)
    errorL2stress2 = zeros(Float64, M + 1, nsamples)
    errorL2u = zeros(Float64, M + 1, nsamples)
    r::Float64 = 1.0
    input_stress_sol = Array{FEVectorBlock, 1}(undef, 2)
    input_u_sol = Array{FEVectorBlock, 1}(undef, 2)
    M0 = dimensionwise_error ? 0 : M

    for s in 1:nsamples
        csample = Samples[:, s]
        for entry in ids_stress
            if entry[1] == 1
                input_stress_sol[1] = sol_det[s][entry[2]]
            else
                input_stress_sol[2] = sol_sgfem[entry[2]]
            end
        end
        for entry in ids_u
            if entry[1] == 1
                input_u_sol[1] = sol_det[s][entry[2]]
            else
                input_u_sol[2] = sol_sgfem[entry[2]]
            end
        end

        ## compare with the SGFEM solution
        for m in M0:M
            ## evaluate SGFEM solution = ∑_k SolutionSGFEM[k] * P_k(s)
            set_sample!(SolutionSGFEM, view(csample, 1:m))

            ## compute errors
            result = ExtendableFEM.evaluate(ErrorIntegratorL2stress, input_stress_sol)
            errorL2stress[m + 1, s] = sum(view(result, energy_norm ? 1 : 2, :))
            errorL2stress2[m + 1, s] = sum(view(result, energy_norm ? 2 : 1, :))
            errorL2u[m + 1, s] = sum(view(ExtendableFEM.evaluate(ErrorIntegratorL2u, input_u_sol), 1, :))
        end

        @info "SAMPLE $s of $nsamples
            sample = $(view(csample, 1:min(Msamples, 2 * M)))
            weight = $(weights[s])
            error(L2stress) = $(errorL2stress[end, s])
            error(L2stress2) = $(errorL2stress2[end, s])
            error(L2u) = $(errorL2u[end, s])"

        if debug ## plot deterministic solution and evaluation of SGFEM solution for sample
            @show extrema(view(sol_det.entries, 1:num_nodes(xgrid)))
            println(stdout, unicode_scalarplot(input_u_sol[1]; title = "det. reference solution for sample $s"))
            @show extrema(view(sol_sgfem.entries, 1:num_nodes(xgrid)))
            println(stdout, unicode_scalarplot(input_u_sol[2]; title = "SGFEM evaluation for sample $s"))
        end
    end

    ## average
    totalerrorL2stress_weighted = zeros(Float64, M + 1)
    totalerrorL2stress2_weighted = zeros(Float64, M + 1)
    totalerrorL2u_weighted = zeros(Float64, M + 1)
    totalerrorL2stress_uniform = zeros(Float64, M + 1)
    totalerrorL2u_uniform = zeros(Float64, M + 1)
    for s in 1:nsamples
        totalerrorL2stress_uniform .+= view(errorL2stress, :, s)
        totalerrorL2u_uniform .+= view(errorL2u, :, s)
        totalerrorL2stress_weighted .+= view(errorL2stress, :, s) * weights[s]
        totalerrorL2stress2_weighted .+= view(errorL2stress2, :, s) * weights[s]
        totalerrorL2u_weighted .+= view(errorL2u, :, s) * weights[s]
    end

    weightsum = sum(weights)
    for m in 1:length(totalerrorL2stress_weighted)
        totalerrorL2stress_uniform[m] = totalerrorL2stress_uniform[m] / nsamples
        totalerrorL2u_uniform[m] = totalerrorL2u_uniform[m] / nsamples
        totalerrorL2stress_weighted[m] = totalerrorL2stress_weighted[m] / weightsum
        totalerrorL2stress2_weighted[m] = totalerrorL2stress2_weighted[m] / weightsum
        totalerrorL2u_weighted[m] = totalerrorL2u_weighted[m] / weightsum
    end

    @info totalerrorL2stress_weighted, totalerrorL2stress2_weighted

    return totalerrorL2stress_weighted, totalerrorL2u_weighted, totalerrorL2stress_uniform, totalerrorL2u_uniform
end
