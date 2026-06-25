####################
### SOLVER STUFF ###
####################

struct StokesPrimal{Tv, MT, VT, GT}
    A0::MT
    A::Vector{MT}
    B::MT
    G::GT
    bdofs::Vector{Int}
    nmodes::Int
    vsize::Array{Int, 1}
end

struct StokesPrimalPreconditioner{Tv, FAC}
    LUS::FAC
    DA::Array{Tv, 1}
    temp::Array{Tv, 1}
    bdofs::Vector{Int}
    nmodes::Int
    vsize::Array{Int, 1}
end

function stokesPrimalPreconditioner(A0::ExtendableSparseMatrix{Tv, Ti}, B::ExtendableSparseMatrix{Tv, Ti}, bdofs, nmodes, vsize) where {Tv, Ti}
    #for dof in bdofs
    #    A0[dof, dof] = 1.0e60
    #end
    #flush!(A0)

    DA::Array{Tv, 1} = zeros(Tv, size(A0, 1))
    for j in 1:length(DA)
        DA[j] = A0[j, j]
    end

    # compute S = B inv(A_diag) B'
    S = ExtendableSparseMatrix{Tv, Ti}(size(B, 2), size(B, 2))
    cscmat::SparseMatrixCSC{Tv, Ti} = B.cscmatrix
    rows::Array{Ti, 1} = rowvals(cscmat)
    valsB::Array{Tv, 1} = cscmat.nzval
    value::Tv = 0
    row::Ti = 0
    for i in 1:size(B, 2), j in 1:size(B, 2)
        for r in nzrange(cscmat, i), r2 in nzrange(cscmat, j)
            if rows[r] == rows[r2]
                row = rows[r]
                value = valsB[r] * valsB[r2] / DA[row]
                _addnz(S, i, j, value, 1)
            end
        end
    end

    # compute LU factorisation of S
    flush!(S)
    LUS = lu(S.cscmatrix)
    #LUA = lu(A0.cscmatrix)

    # temporary storage array for solver
    temp = zeros(Tv, size(B, 2))

    return StokesPrimalPreconditioner{Tv, typeof(LUS)}(LUS, DA, temp, bdofs, nmodes, vsize)
end

@inline LinearAlgebra.ldiv!(C::StokesPrimalPreconditioner, b) = ldiv!(b, C, b)
@inline function LinearAlgebra.ldiv!(y, C::StokesPrimalPreconditioner{Tv, FAC}, b) where {Tv, FAC}
    a::Int = 0
    c::Int = 0
    DA::Array{Tv, 1} = C.DA
    temp::Array{Tv, 1} = C.temp
    nmodes::Int = C.nmodes
    vsize::Array{Int, 1} = C.vsize
    for mu in 1:nmodes
        # upper left block of preconditioner (I ⊗ A_diag)
        a = (mu - 1) * vsize[1] + 1
        c = mu * vsize[1]
        for i in a:c
            y[i] = b[i] / DA[i - a + 1]
        end
        a = nmodes * vsize[1] + (mu - 1) * vsize[2] + 1
        c = nmodes * vsize[1] + mu * vsize[2]
        if y !== b
            ldiv!(view(y, a:c), C.LUS, view(b, a:c))
        else
            ldiv!(temp, C.LUS, view(b, a:c))
            y[a:c] .= temp
        end
    end

    return y
end

@inline function LinearAlgebra.:\(C::StokesPrimalPreconditioner, b)
    y = zero(b)
    ldiv!(y, C, b)
    return y
end

function LinearAlgebra.mul!(Ax::Vector{Tv}, S::StokesPrimal{Tv, MT, VT, GT}, x) where {Tv, MT, VT, GT}
    fill!(Ax, 0)
    g::Tv = 0
    vsize::Array{Int, 1} = S.vsize
    nmodes::Int = S.nmodes
    bdofs::Vector{Int} = S.bdofs
    G::GT = S.G
    A::Vector{MT} = S.A
    M::Int = length(A) # size(G,1) / nmodes
    a::Int = 0
    b::Int = 0
    a2::Int = 0
    b2::Int = 0
    A0::MT = S.A0
    B::MT = S.B
    for mu in 1:nmodes
        # deterministic part
        a = (mu - 1) * vsize[1] + 1
        b = mu * vsize[1]
        a2 = a
        b2 = b
        addblock_matmul!(view(Ax, a:b), A0[1, 1], view(x, a2:b2))
        a2 = nmodes * vsize[1] + (mu - 1) * vsize[2] + 1
        b2 = nmodes * vsize[1] + mu * vsize[2]
        addblock_matmul!(view(Ax, a:b), B[1, 1], view(x, a2:b2))
        a2 = a
        b2 = b
        a = nmodes * vsize[1] + (mu - 1) * vsize[2] + 1
        b = nmodes * vsize[1] + mu * vsize[2]
        addblock_matmul!(view(Ax, a:b), B[1, 1], view(x, a2:b2); transposed = true)

        # stochastic part
        a = (mu - 1) * vsize[1] + 1
        b = mu * vsize[1]
        for nu in 1:nmodes, e in 1:M
            g = G[(e - 1) * nmodes + mu, nu]
            if abs(g) > 1.0e-12
                a2 = (nu - 1) * vsize[1] + 1
                b2 = nu * vsize[1]
                addblock_matmul!(view(Ax, a:b), A[e][1, 1], view(x, a2:b2); factor = g)
            end
        end

        for dof in bdofs
            Ax[a + dof - 1] = 0
        end
    end

    return nothing
end

Base.eltype(S::StokesPrimal) = typeof(S).parameters[1]
Base.size(S::StokesPrimal) = S.nmodes .* (size(S.A0.entries) .+ size(S.B.entries)[2])

function solve_stokes_primal!(SolutionSGFEM::SGFEVector, A0, A, B, b0, G, nmodes, bfac; atol = 1.0e-14, rtol = 1.0e-14)
    ## create fullmatrix-free matrix evaluator
    @info "Solving StochasticFEM iteratively and matrix-free (ndofs = $(length(SolutionSGFEM.entries)))..."
    vsize = [SolutionSGFEM[1].FES.ndofs, SolutionSGFEM[nmodes + 1].FES.ndofs]
    ## boundary data
    bfacedofs = SolutionSGFEM.FES_space[1][BFaceDofs]
    nbfaces = num_sources(bfacedofs)
    bdofs = []
    for bface in 1:nbfaces
        append!(bdofs, view(bfacedofs, :, bface))
    end
    bdofs = unique(bdofs)

    S = StokesPrimal{eltype(G), typeof(A0), typeof(b0), typeof(G)}(A0, A, B, G, bdofs, nmodes, vsize)
    @info "...initializing Preconditioner"
    @time P = stokesPrimalPreconditioner(A0.entries, B.entries, bdofs, nmodes, vsize)

    ## right-hand side
    b = deepcopy(SolutionSGFEM)
    fill!(b.entries, 0)
    addblock!(b[1], b0[1]; factor = bfac)
    for m in 1:nmodes
        for dof in bdofs
            b[m][dof] = 0
        end
    end

    ## solve
    @info "...starting right-conditioned GMRES"
    x, history = Krylov.gmres(S, b.entries, SolutionSGFEM.entries; ldiv = true, atol = atol, rtol = rtol, M = P)
    SolutionSGFEM.entries .= x
    @show history

    ## Pressure mean von pressure modes abziehen, damit Vorfaktor wie bei use_iterative_solver=false ist.
    ## https://wias-pdelib.github.io/ExtendableFEM.jl/stable/module_examples/Example252_NSEPlanarLatticeFlow/
    xgrid = SolutionSGFEM[1].FES.xgrid
    for i in 1:nmodes
        pintegrate = ItemIntegrator([id(1)])
        pmean = sum(ExtendableFEM.evaluate(pintegrate, [SolutionSGFEM.FEVectorBlocks[nmodes + i]])) / sum(xgrid[CellVolumes])
        view(SolutionSGFEM.FEVectorBlocks[nmodes + i]) .-= pmean
    end

    ## check residual
    Ax = zero(SolutionSGFEM.entries)
    mul!(Ax, S, SolutionSGFEM.entries)
    @info "solver residual = $(sqrt(sum((Ax - b.entries) .^ 2)))"

    return bdofs
end

function solve_stokes_primal_full!(SolutionSGFEM::SGFEVector, A0, A, B, b0, G, nmodes, rhsfac)
    M::Int = length(A) # size(G,1) / nmodes
    FES = SolutionSGFEM.FES_space
    multi_indices = SolutionSGFEM.TB.multi_indices
    nmodes = num_multiindices(SolutionSGFEM)
    bigFES = Array{FESpace{Float64, Int32}, 1}([FES[1] for j in 1:nmodes])
    append!(bigFES, [FES[2] for j in 1:nmodes])

    bigS = FEMatrix(bigFES)
    bigb = FEVector(bigFES)

    for j in 1:nmodes
        addblock!(bigS[j, j], A0[1, 1])
        addblock!(bigS[j, nmodes + j], B[1, 1])
        addblock!(bigS[nmodes + j, j], B[1, 1]; transpose = true)
    end

    ## right-hand side
    addblock!(bigb[1], b0[1])

    g::Float64 = 0
    for j in 1:nmodes, k in 1:nmodes
        for e in 1:M
            g = G[(e - 1) * nmodes + j, k] # ⟨ ξ_m ψ_mi(j) ψ_mi(k) ⟩
            if abs(g) > 1.0e-12
                addblock!(bigS[j, k], A[e][1, 1]; factor = g)
            end
        end
    end
    flush!(bigS.entries)

    ## boundary data
    bfacedofs = FES[1][BFaceDofs]
    nbfaces = num_sources(bfacedofs)
    bdofs = []
    for bface in 1:nbfaces
        append!(bdofs, view(bfacedofs, :, bface))
    end
    unique!(bdofs)

    for m in 1:nmodes
        for dof in bdofs
            bigS[m, m][dof, dof] = 1.0e60
            bigb[m][dof] = 0
        end
    end
    flush!(bigS.entries)

    @info "Solving StochasticFEM with full matrix..."
    SolutionSGFEM.entries .= bigS.entries \ bigb.entries

    xgrid = SolutionSGFEM[1].FES.xgrid
    for i in 1:nmodes
        pintegrate = ItemIntegrator([id(1)])
        pmean = sum(ExtendableFEM.evaluate(pintegrate, [SolutionSGFEM.FEVectorBlocks[nmodes + i]])) / sum(xgrid[CellVolumes])
        view(SolutionSGFEM.FEVectorBlocks[nmodes + i]) .-= pmean
    end

    residual = bigS.entries * SolutionSGFEM.entries .- bigb.entries
    for m in 1:nmodes
        residual[FES[1].ndofs * (m - 1) .+ bdofs] .= 0
    end
    println("linear residual = $(sqrt(sum(residual .^ 2)))")

    return bdofs
end
