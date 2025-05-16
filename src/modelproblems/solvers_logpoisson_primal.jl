####################
### SOLVER STUFF ###
####################


### Preconditioner needs to implement these functions
#struct Identity end
#
#\(::Identity, x) = copy(x)
#ldiv!(::Identity, x) = x
#ldiv!(y, ::Identity, x) = copyto!(y, x)


struct MySystemLogPrimal{Tv, MT, VT, GT}
    A::MT
    N0::MT
    Nm::Vector{MT}
    G::GT
    bdofs::Vector{Int}
    nmodes::Int
end
Base.size(S::MySystemLogPrimal) = S.nmodes .* size(S.A.entries)

struct MyPreconditionerLogPrimal{Tv, FAC}
    LUA::FAC
    DA::Array{Tv, 1}
    bdofs::Vector{Int}
    nmodes::Int
end

function MyPreconditionerLogPrimal(A::ExtendableSparseMatrix{Tv, Ti}, bdofs, nmodes) where {Tv, Ti}
    DA::Array{Tv, 1} = zeros(Tv, size(A, 1))
    for j in 1:length(DA)
        DA[j] = A[j, j]
    end

    # compute LU factorisation of S
    for dof in bdofs
        A[dof, dof] = 1.0e60
    end
    flush!(A)
    LUA = LUFactorization(A)

    return MyPreconditionerLogPrimal{Tv, typeof(LUA)}(LUA, DA, bdofs, nmodes)
end


#@inline LinearAlgebra.ldiv!(x::AbstractArray, C::ExtendableSparseMatrix, b::AbstractArray) = x = C\b
@inline LinearAlgebra.ldiv!(C::MyPreconditionerLogPrimal, b) = ldiv!(b, C, b)
@inline function LinearAlgebra.ldiv!(y, C::MyPreconditionerLogPrimal{Tv, FAC}, b) where {Tv, FAC}
    #copyto!(y, b)
    #return
    a::Int = 0
    c::Int = 0
    DA::Array{Tv, 1} = C.DA
    #bdofs = C.bdofs
    vsize::Int = length(DA)
    nmodes::Int = C.nmodes
    for mu::Int in 1:nmodes
        # upper left block of preconditioner (I ⊗ A_diag)
        a = (mu - 1) * vsize + 1
        c = mu * vsize
        if (false)
            for i in a:c
                y[i] = b[i] / DA[i - a + 1]
            end
        else
            if y !== b
                ldiv!(view(y, a:c), C.LUA, view(b, a:c))  # y = A\b
            else
                temp = zeros(Tv, vsize)
                ldiv!(temp, C.LUA, view(b, a:c))
                y[a:c] .= temp
            end
        end
        #if dof in bdofs
        # y[a+dof-1] = 0
        #end
    end
    return y
end
@inline function LinearAlgebra.:\(C::MyPreconditionerLogPrimal, b)
    y = zero(b)
    ldiv!(y, C, b)
    return y
end


function LinearAlgebra.mul!(Ax, S::MySystemLogPrimal{Tv, MT, VT, GT}, x) where {Tv, MT, VT, GT}
    fill!(Ax, 0)
    g::Tv = 0
    A::MT = S.A
    vsize::Int = size(A.entries, 1)
    nmodes::Int = S.nmodes
    bdofs::Vector{Int} = S.bdofs
    G::GT = S.G
    Nm::Vector{MT} = S.Nm
    M::Int = length(Nm) # size(G,1) / nmodes
    a::Int = 0
    b::Int = 0
    a2::Int = 0
    b2::Int = 0
    N0::MT = S.N0

    for mu::Int in 1:nmodes
        # deterministic part
        a = (mu - 1) * vsize + 1
        b = mu * vsize
        a2 = a
        b2 = b
        addblock_matmul!(view(Ax, a:b), A[1, 1], view(x, a2:b2))
        addblock_matmul!(view(Ax, a:b), N0[1, 1], view(x, a2:b2))

        # stochastic part
        for nu::Int in 1:nmodes, e::Int in 1:M
            g = G[(e - 1) * nmodes + mu, nu]
            if g != 0
                a2 = (nu - 1) * vsize + 1
                b2 = nu * vsize
                addblock_matmul!(view(Ax, a:b), Nm[e][1, 1], view(x, a2:b2); factor = g)
            end
        end

        for dof in bdofs
            Ax[a + dof - 1] = 0
        end
    end
    return
end

Base.eltype(S::MySystemLogPrimal) = typeof(S).parameters[1]
Base.size(S::MySystemLogPrimal, d::Int) = S.nmodes * size(S.A.entries, 1)


function solve_logpoisson_primal!(SolutionSGFEM::SGFEVector, A, N0, Nm, b0, G, nmodes, bfac; atol = 1.0e-14, rtol = 1.0e-14)

    ## create fullmatrix-free matrix evaluator
    @info "Solving StochasticFEM iteratively and matrix-free (ndofs = $(length(SolutionSGFEM)))..."

    ## boundary data
    bfacedofs = SolutionSGFEM.FES_space[1][BFaceDofs]
    nbfaces = num_sources(bfacedofs)
    bdofs = []
    for bface in 1:nbfaces
        append!(bdofs, view(bfacedofs, :, bface))
    end
    bdofs = unique(bdofs)

    S = MySystemLogPrimal{eltype(G), typeof(A), typeof(b0), typeof(G)}(A, N0, Nm, G, bdofs, nmodes)
    @info "...initializing Preconditioner"
    @time P = MyPreconditionerLogPrimal(A.entries, bdofs, nmodes)

    ## right-hand side
    b = deepcopy(SolutionSGFEM)
    for m in 1:nmodes
        addblock!(b[m], b0[m][1])
        for dof in bdofs
            b[m][dof] = 0
        end
    end

    ## solve
    @info "...starting preconditioned GMRES"
    #x, stats = IterativeSolvers.gmres!(SolutionSGFEM.entries, S, b.entries; log = true, Pr = P, Pl = P)
    x, stats = Krylov.gmres(S, b.entries, SolutionSGFEM.entries; ldiv = true, atol = atol, rtol = rtol, M = P)
    SolutionSGFEM.entries .= x
    @show stats

    ## check residual
    Ax = zero(SolutionSGFEM.entries)
    mul!(Ax, S, SolutionSGFEM.entries)
    @info "solver residual = $(sqrt(sum((Ax - b.entries) .^ 2)))"
    return bdofs
end


function solve_logpoisson_primal_full!(SolutionSGFEM::SGFEVector, A, N0, N, b, G, nmodes, rhsfac)

    M::Int = length(N) # size(G,1) / nmodes

    FES = SolutionSGFEM.FES_space[1]
    nmodes = num_multiindices(SolutionSGFEM)
    bigFES = [FES for j in 1:nmodes]
    x::Vector{Float64} = zeros(Float64, 2)

    bigS = FEMatrix(bigFES)
    bigb = FEVector(bigFES)

    ## add AA and BB and N0
    for j in 1:nmodes
        addblock!(bigS[j, j], A[1, 1])
        if N0 !== nothing
            addblock!(bigS[j, j], N0[1, 1])
        end
    end

    ## add Nm blocks of NN
    g::Float64 = 0
    for j in 1:nmodes, k in 1:nmodes
        for e in 1:M
            g = G[(e - 1) * nmodes + j, k] # ⟨ ξ_m ψ_mi(j) ψ_mi(k) ⟩
            if abs(g) > 1.0e-12
                #@show g, [e,j,k]
                addblock!(bigS[j, k], N[e][1, 1]; factor = g)
            end
        end
    end

    ## right-hand side
    for m in 1:nmodes
        addblock!(bigb[m], b[m][1])
    end

    ## boundary data
    bfacedofs = FES[BFaceDofs]
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
    #gmres!(SolutionSGFEM.entries, bigS.entries, bigb.entries)

    residual = bigS.entries * SolutionSGFEM.entries .- bigb.entries
    println("linear residual = $(sqrt(sum(residual .^ 2)))")
    return bdofs
end
