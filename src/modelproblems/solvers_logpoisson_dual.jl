####################
### SOLVER STUFF ###
####################


### Preconditioner needs to implement these functions
#struct Identity end
#
#\(::Identity, x) = copy(x)
#ldiv!(::Identity, x) = x
#ldiv!(y, ::Identity, x) = copyto!(y, x)


struct MySystemLogDual{Tv, MT, VT, GT}
    A::MT
    B::MT
    N0::MT
    Nm::Vector{MT}
    G::GT
    nmodes::Int
    vsize::Array{Int, 1}
    couplings::Vector{Tuple{Int, Int, Int, Tv}}
end
Base.size(S::MySystemLogDual) = S.nmodes .* (size(S.A.entries) .+ size(S.B.entries)[2])

"""
    MySystemLogDual(A, B, N0, Nm, G, nmodes, vsize)

Convenience constructor that precomputes the nonzero coupling list `(mu, nu, e, g)` from `G`.
"""
function MySystemLogDual(A, B, N0, Nm, G, nmodes, vsize)
    Tv = eltype(G)
    M = length(Nm)
    couplings = Tuple{Int, Int, Int, Tv}[]
    for mu in 1:nmodes, nu in 1:nmodes, e in 1:M
        g = G[(e - 1) * nmodes + mu, nu]
        if g != 0
            push!(couplings, (mu, nu, e, g))
        end
    end
    return MySystemLogDual{Tv, typeof(A), typeof(couplings), typeof(G)}(A, B, N0, Nm, G, nmodes, vsize, couplings)
end

struct MyPreconditionerLogDual{Tv, FAC}
    LUS::FAC
    DA::Array{Tv, 1}
    temp::Array{Tv, 1}
    nmodes::Int
    vsize::Array{Int, 1}
end

function MyPreconditionerLogDual(A::ExtendableSparseMatrix{Tv, Ti}, B::ExtendableSparseMatrix{Tv, Ti}, nmodes, vsize) where {Tv, Ti}
    DA::Array{Tv, 1} = zeros(Tv, size(A, 1))
    for j in 1:length(DA)
        DA[j] = A[j, j]
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

    # temporary storage array for solver
    temp = zeros(Tv, size(B, 2))

    return MyPreconditionerLogDual{Tv, typeof(LUS)}(LUS, DA, temp, nmodes, vsize)
end


@inline LinearAlgebra.ldiv!(C::MyPreconditionerLogDual, b) = ldiv!(b, C, b)
@inline function LinearAlgebra.ldiv!(y, C::MyPreconditionerLogDual{Tv, FAC}, b) where {Tv, FAC}
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
@inline function LinearAlgebra.:\(C::MyPreconditionerLogDual, b)
    y = zero(b)
    ldiv!(y, C, b)
    return y
end


function LinearAlgebra.mul!(Ax::Vector{Tv}, S::MySystemLogDual{Tv, MT, VT, GT}, x) where {Tv, MT, VT, GT}
    fill!(Ax, 0)
    A::MT = S.A
    N0::MT = S.N0
    B::MT = S.B
    Nm::Vector{MT} = S.Nm
    vsize::Array{Int, 1} = S.vsize
    nmodes::Int = S.nmodes
    couplings = S.couplings
    a::Int = 0
    b::Int = 0
    a2::Int = 0
    b2::Int = 0
    for mu in 1:nmodes
        # deterministic part
        a = (mu - 1) * vsize[1] + 1
        b = mu * vsize[1]
        a2 = a
        b2 = b
        addblock_matmul!(view(Ax, a:b), A[1, 1], view(x, a2:b2))
        addblock_matmul!(view(Ax, a:b), N0[1, 1], view(x, a2:b2))
        a2 = nmodes * vsize[1] + (mu - 1) * vsize[2] + 1
        b2 = nmodes * vsize[1] + mu * vsize[2]
        addblock_matmul!(view(Ax, a:b), B[1, 1], view(x, a2:b2))
        a2 = a
        b2 = b
        a = nmodes * vsize[1] + (mu - 1) * vsize[2] + 1
        b = nmodes * vsize[1] + mu * vsize[2]
        addblock_matmul!(view(Ax, a:b), B[1, 1], view(x, a2:b2); transposed = true)
    end

    # stochastic part
    for (mu, nu, e, g) in couplings
        a = (mu - 1) * vsize[1] + 1
        b = mu * vsize[1]
        a2 = nmodes * vsize[1] + (nu - 1) * vsize[2] + 1
        b2 = nmodes * vsize[1] + nu * vsize[2]
        addblock_matmul!(view(Ax, a:b), Nm[e][1, 1], view(x, a2:b2); factor = g)
    end
    return Ax
end

Base.eltype(::MySystemLogDual{Tv}) where {Tv} = Tv
Base.size(S::MySystemLogDual, d::Int) = S.nmodes * (S.vsize[1] + S.vsize[2])

function solve_logpoisson_dual!(SolutionSGFEM::SGFEVector, A, B, N0, Nm, b0, G, nmodes, bfac; atol = 1.0e-14, rtol = 1.0e-14)

    ## create fullmatrix-free matrix evaluator
    @info "Solving StochasticFEM iteratively and matrix-free (ndofs = $(length(SolutionSGFEM.entries)))..."
    vsize = [SolutionSGFEM[1].FES.ndofs, SolutionSGFEM[nmodes + 1].FES.ndofs]
    S = MySystemLogDual(A, B, N0, Nm, G, nmodes, vsize)
    @info "...initializing Preconditioner"
    @time P = MyPreconditionerLogDual(A.entries, B.entries, nmodes, vsize)

    ## right-hand side (only has one deterministic block)
    b = deepcopy(SolutionSGFEM)
    fill!(b.entries, 0)
    addblock!(b[nmodes + 1], b0[1]; factor = bfac)

    ## solve via LinearSolve's higher-level API while keeping the matrix-free operator
    @info "...starting right-preconditioned GMRES via LinearSolve"
    Aop = FunctionOperator(
        (y, x, u, p, t) -> mul!(y, S, x),
        zero(b.entries),
        zero(b.entries);
        T = eltype(b.entries),
        islinear = true,
        isconstant = true,
        ifcache = false,
    )
    prob = LinearProblem(Aop, b.entries)
    sol = solve(
        prob;
        alg = KrylovJL_GMRES(),
        abstol = atol,
        reltol = rtol,
        Pl = P,
        Pr = P,
        verbose = false,
    )
    SolutionSGFEM.entries .= sol.u
    @show sol.stats

    ## check residual
    Ax = zero(SolutionSGFEM.entries)
    mul!(Ax, S, SolutionSGFEM.entries)
    return @info "solver residual = $(sqrt(sum((Ax - b.entries) .^ 2)))"
end


function solve_logpoisson_dual_full!(SolutionSGFEM::SGFEVector, A, B, N0, N, b0, G, nmodes)

    M::Int = length(N) # size(G,1) / nmodes
    FES = SolutionSGFEM.FES_space
    bigFES = Array{FESpace{Float64, Int32}, 1}([FES[1] for j in 1:nmodes])
    append!(bigFES, [FES[2] for j in 1:nmodes])

    bigS = FEMatrix(bigFES)
    bigb = FEVector(bigFES)

    ## add AA and BB and N0
    for j in 1:nmodes
        addblock!(bigS[j, j], A[1, 1])
        addblock!(bigS[j, nmodes + j], N0[1, 1])
        addblock!(bigS[j, nmodes + j], B[1, 1])
        addblock!(bigS[nmodes + j, j], B[1, 1]; transpose = true)
    end

    ## rhs
    addblock!(bigb[nmodes + 1], b0[1])

    ## add Nm blocks of NN
    g::Float64 = 0
    for j in 1:nmodes, k in 1:nmodes
        for e in 1:M
            g = G[(e - 1) * nmodes + j, k]
            if abs(g) > 1.0e-12
                addblock!(bigS[j, nmodes + k], N[e][1, 1]; factor = g)
            end
        end
    end

    @info "Solving StochasticFEM with full matrix..."
    SolutionSGFEM.entries .= bigS.entries \ bigb.entries
    #gmres!(SolutionSGFEM.entries, bigS.entries, bigb.entries)

    residual = bigS.entries * SolutionSGFEM.entries .- bigb.entries
    println("linear residual = $(sqrt(sum(residual .^ 2)))")
    return []
end
