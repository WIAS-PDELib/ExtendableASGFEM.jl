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
    couplings::Vector{Tuple{Int, Int, Int, Tv}}
end
Base.size(S::MySystemLogPrimal) = S.nmodes .* size(S.A.entries)

"""
    MySystemLogPrimal(A, N0, Nm, G, bdofs, nmodes)

Convenience constructor that precomputes the nonzero coupling list `(mu, nu, e, g)` from `G`.
"""
function MySystemLogPrimal(A, N0, Nm, G, bdofs, nmodes)
    Tv = eltype(G)
    M = length(Nm)
    couplings = Tuple{Int, Int, Int, Tv}[]
    for mu in 1:nmodes, nu in 1:nmodes, e in 1:M
        g = G[(e - 1) * nmodes + mu, nu]
        if g != 0
            push!(couplings, (mu, nu, e, g))
        end
    end
    return MySystemLogPrimal{Tv, typeof(A), typeof(couplings), typeof(G)}(A, N0, Nm, G, bdofs, nmodes, couplings)
end

struct MyPreconditionerLogPrimal{Tv, FAC}
    LUA::FAC
    DA::Array{Tv, 1}
    temp::Array{Tv, 1}
    bdofs::Vector{Int}
    nmodes::Int
end

function MyPreconditionerLogPrimal(A::ExtendableSparseMatrix{Tv, Ti}, bdofs, nmodes) where {Tv, Ti}
    DA::Array{Tv, 1} = zeros(Tv, size(A, 1))
    for j in 1:length(DA)
        DA[j] = A[j, j]
    end

    # compute LU factorisation of A
    for dof in bdofs
        A[dof, dof] = 1.0e60
    end
    flush!(A)
    LUA = lu(A.cscmatrix)

    # temporary storage array for solver
    temp = zeros(Tv, size(A, 1))

    return MyPreconditionerLogPrimal{Tv, typeof(LUA)}(LUA, DA, temp, bdofs, nmodes)
end

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
                ldiv!(C.temp, C.LUA, view(b, a:c))
                y[a:c] .= C.temp
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
    A::MT = S.A
    N0::MT = S.N0
    Nm::Vector{MT} = S.Nm
    vsize::Int = size(A.entries, 1)
    nmodes::Int = S.nmodes
    bdofs::Vector{Int} = S.bdofs
    couplings = S.couplings
    a::Int = 0
    b::Int = 0
    a2::Int = 0
    b2::Int = 0

    for mu::Int in 1:nmodes
        # deterministic part
        a = (mu - 1) * vsize + 1
        b = mu * vsize
        addblock_matmul!(view(Ax, a:b), A[1, 1], view(x, a:b))
        addblock_matmul!(view(Ax, a:b), N0[1, 1], view(x, a:b))
    end

    # stochastic part
    for (mu, nu, e, g) in couplings
        a = (mu - 1) * vsize + 1
        b = mu * vsize
        a2 = (nu - 1) * vsize + 1
        b2 = nu * vsize
        addblock_matmul!(view(Ax, a:b), Nm[e][1, 1], view(x, a2:b2); factor = g)
    end

    for mu::Int in 1:nmodes
        a = (mu - 1) * vsize + 1
        for dof in bdofs
            Ax[a + dof - 1] = 0
        end
    end
    return Ax
end

Base.eltype(::MySystemLogPrimal{Tv}) where {Tv} = Tv
Base.size(S::MySystemLogPrimal, d::Int) = S.nmodes * size(S.A.entries, 1)


function solve_logpoisson_primal!(SolutionSGFEM::SGFEVector, A, N0, Nm, b0, G, nmodes; atol = 1.0e-14, rtol = 1.0e-14)

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

    S = MySystemLogPrimal(A, N0, Nm, G, bdofs, nmodes)
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

    ## solve via LinearSolve's higher-level API while keeping the matrix-free operator
    @info "...starting preconditioned GMRES via LinearSolve"
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
    @info "solver residual = $(sqrt(sum((Ax - b.entries) .^ 2)))"
    return bdofs
end


function solve_logpoisson_primal_full!(SolutionSGFEM::SGFEVector, A, N0, N, b, G, nmodes)

    M::Int = length(N) # size(G,1) / nmodes

    FES = SolutionSGFEM.FES_space[1]
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
