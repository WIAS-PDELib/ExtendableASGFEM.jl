
"""
$(TYPEDEF)

structure that builds upon ExtendableFEMBase.FEVector and adds information
on the stochastic discretization, i.e. the used TensoriedBasis

"""
struct SGFEVector{T,Tv, Ti, ONBType <: ONBasis, MIType}
    FES_space::Array{FESpace{Tv, Ti},1}        ## finite element space for all modes
    TB::TensorizedBasis{Tv, ONBType, MIType}    ## tensorized basis with
    active_modes::Array{Ti,1}                   ## indices of active modes (w.r.t. TB.mulit_indices)
    length4modes::Array{Ti,1}                   ## offsets for each mode
	FEVectorBlocks::Array{FEVectorBlock{T, Tv, Ti}, 1} ## Vector block mask for each multi-index
    entries::Array{T,1}                         ## complete coefficient vector
    last_sample::Array{Tv,1}                    ## last sample used for evaluation
    FEV::FEVector{T,Tv,Ti}                      ## Vector used for evaluating 
end

"""
$(TYPEDSIGNATURES)

evaluates the SGVector at that sample and stores the results in SGFEV.FEV.

"""
function set_sample!(SGFEV::SGFEVector, S::AbstractVector)
    last_sample = SGFEV.last_sample
    fill!(last_sample, 0)
    if length(S) <= length(last_sample)
        last_sample[1:length(S)] .= S
    else
        last_sample .= view(S, 1:length(last_sample))
        @warn "ignoring end of sample, since SGFEVector only has multi-indices up to m = $(length(last_sample))"
    end

    ## evalute SGFEVector at sample
    FEV = SGFEV.FEV
    TB = SGFEV.TB
    nmodes = num_multiindices(SGFEV)
    set_sample!(TB, S; normalize = true)
    nunknowns = length(FEV)
    for u = 1 : nunknowns
        fill!(FEV[u],0)
        for k = 1 : nmodes
            r = ExtendableASGFEM.evaluate(TB, k) # = P_k(s)
            if r != 0
                view(FEV[u]) .+= r * view(SGFEV[(u-1)*nmodes + k])
            end
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

contructs an SGFEVector for the given (spatial) finite element space `FES` and the
(tensorized basis of the stochastic discretization) `TB`. 

"""
function SGFEVector(FES::Array{<:FESpace{Tv,Ti},1}, TB::TensorizedBasis{Tv, ONBType, MIType}; active_modes = :all, T = Tv, unames = 1:length(FES)) where {Tv, Ti, ONBType, MIType}
    if active_modes == :all
        active_modes = 1:nmodes(TB)
    end
    length4modes = [0]
    feblocks = Array{FEVectorBlock{T, Tv, Ti}, 1}(undef, 0)
    entries = zeros(T, 0)
    nunknowns = length(FES)
    for j = 1 : length(FES)
        fes = FES[j]
        ndofs = fes.ndofs
        append!(entries, zeros(Float64, ndofs*length(active_modes)))
        for m in active_modes
            name = nunknowns == 1 ? "$(TB.multi_indices[m])" : "$(unames[j]) $(TB.multi_indices[m])"
            push!(feblocks, FEVectorBlock{T, Tv, Ti, eltype(fes), ExtendableFEMBase.assemblytype(fes)}(name, fes, length4modes[end], length4modes[end] + ndofs, entries))
            push!(length4modes, length4modes[end]+ndofs)
        end
    end

    return SGFEVector{T, Tv, Ti, ONBType, MIType}(FES, TB, active_modes, length4modes, feblocks, entries, zeros(Tv, maxlength_multiindices(TB)), FEVector(FES))
end

function SGFEVector(FES::FESpace{Tv, Ti}, TB::TensorizedBasis{Tv}; kwargs...) where {Tv,Ti}
    return SGFEVector(Array{FESpace{Tv,Ti},1}([FES]), TB; kwargs...)
end

"""
$(TYPEDSIGNATURES)

returns the `i`-th stochastic mode

"""
Base.getindex(SGFEV::SGFEVector, i::Int) = SGFEV.FEVectorBlocks[i]
"""
$(TYPEDSIGNATURES)

returns the FEVectorBlock for the `i`-th stochastic mode of the `u`-th unknown

"""
Base.getindex(SGFEV::SGFEVector, u::Int, i::Int) = SGFEV.FEVectorBlocks[(u-1)*length(SGFEV.active_modes) + i]
"""
$(TYPEDSIGNATURES)

returns a tuple with the number of active stochastic modes and the number of spatial degrees of freedom

"""
Base.size(SGFEV::SGFEVector) = (length(SGFEV.active_modes), SGFEV.FES_space.ndofs)

"""
$(TYPEDSIGNATURES)

returns the length of the full vector, i.e., the total number of degrees of freedom

"""
Base.length(SGFEV::SGFEVector) = length(SGFEV.entries)


"""
$(TYPEDSIGNATURES)

returns the number of active modes (that are used from the stored tensorized basis)

"""
num_multiindices(SGFEV::SGFEVector) = length(SGFEV.active_modes)

"""
$(TYPEDSIGNATURES)

returns the finite element types of the (spatial) FE spaces

"""
fetypes(SGFEV::SGFEVector) = [eltype(SGFEV.FES_space[j]) for j = 1 : length(SGFEV.FES_space)]

"""
$(TYPEDSIGNATURES)

Custom `show` function for `FEVector` that prints some information on its blocks.
"""
function Base.show(io::IO, SGFEV::SGFEVector)
	println(io, "\nSGFEVector information")
	println(io, "======================")
    println(io, "FETypes = ", fetypes(SGFEV))
    println(io, "")
    ndofs = [SGFEV.FES_space[j].ndofs for j = 1 : length(SGFEV.FES_space)]
	print(io, "   block  |     MI\t|    ndofs \t|    min  /  max")
	for j = 1 : num_multiindices(SGFEV)
		@printf(io, "\n [%5d]  | ", j)
		@printf(io, "  %s\t| ", "$(SGFEV.TB.multi_indices[SGFEV.active_modes[j]])")
		@printf(io, " %8d\t|", sum(ndofs))
		ext = extrema(SGFEV[j])
		@printf(io, " %.2e/%.2e", ext[1], ext[2])
	end
    @printf(io, "\n\t\t\t %10d (total)", length(SGFEV))

end


function norms(SGFEV::SGFEVector{T}, p::Real = 2) where {T}
    nmodes = num_multiindices(SGFEV)
	n = zeros(T, nmodes)
	for j ∈ 1:nmodes
		n[j] = norm(SGFEV[j], p)
	end
	return n
end
