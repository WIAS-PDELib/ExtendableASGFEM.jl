####################
### SOLVER STUFF ###
####################


### Preconditioner needs to implement these functions
#struct Identity end
#
#\(::Identity, x) = copy(x)
#ldiv!(::Identity, x) = x
#ldiv!(y, ::Identity, x) = copyto!(y, x)


"""
    MySystemPrimal{Tv, MT, VT, GT}

Matrix-free block system evaluator for the primal stochastic Poisson problem.

Represents the block-structured linear operator arising from the stochastic Galerkin
discretisation with a Karhunen-Loeve expanded diffusion coefficient. The operator acts
on vectors partitioned by stochastic mode and implements `mul!` for matrix-vector
products without assembling the full block matrix.

# Fields
- `A0::MT`: Mean diffusion block.
- `Am::Vector{MT}`: KL perturbation blocks (one per random dimension).
- `G::GT`: Coupling tensor with entries `G[e, mu, nu]` mapping KL indices `e` to mode pairs.
- `bdofs::Vector{Int}`: Boundary dofs for homogeneous Dirichlet conditions.
- `nmodes::Int`: Number of stochastic modes.
- `couplings::Vector{Tuple{Int, Int, Int, Tv}}`: Nonzero entries of `G` precomputed as
  `(mu, nu, e, g)` tuples, so that `mul!` never iterates over vanishing couplings.
"""
struct MySystemPrimal{Tv, MT, VT, GT}
    A0::MT
    Am::Vector{MT}
    G::GT
    bdofs::Vector{Int}
    nmodes::Int
    couplings::Vector{Tuple{Int, Int, Int, Tv}}
end
Base.size(S::MySystemPrimal) = S.nmodes .* size(S.A0.entries)

"""
    MySystemPrimal(A0, Am, G, bdofs, nmodes)

Convenience constructor that precomputes the nonzero coupling list `(mu, nu, e, g)` from `G`.
"""
function MySystemPrimal(A0, Am, G, bdofs, nmodes)
    Tv = eltype(G)
    M = length(Am)
    couplings = Tuple{Int, Int, Int, Tv}[]
    for mu in 1:nmodes, nu in 1:nmodes, e in 1:M
        g = G[(e - 1) * nmodes + mu, nu]
        if g != 0
            push!(couplings, (mu, nu, e, g))
        end
    end
    return MySystemPrimal{Tv, typeof(A0), typeof(couplings), typeof(G)}(A0, Am, G, bdofs, nmodes, couplings)
end

"""
    MyPreconditionerPrimal{Tv, FAC}

Block-diagonal LU preconditioner for the primal SG Poisson system.

Each diagonal block is the inverse of the mean stiffness matrix ``A_0``, obtained via
LU factorisation. Boundary dofs are handled by stiffening those diagonal entries in the
pre-factored ``A_0`` so that forward/backward substitution implicitly enforces zero
Dirichlet conditions.

# Fields
- `LUA::FAC`: LU factorisation of ``A_0`` with stiffened boundary dofs.
- `DA::Vector{Tv}`: Diagonal of ``A_0`` (kept for compatibility, not used in current path).
- `temp::Vector{Tv}`: Scratch vector reused by the in-place path of `ldiv!`.
- `bdofs::Vector{Int}`: Boundary dofs.
- `nmodes::Int`: Number of stochastic modes.
"""
struct MyPreconditionerPrimal{Tv, FAC}
    LUA::FAC
    DA::Array{Tv, 1}
    temp::Array{Tv, 1}
    bdofs::Vector{Int}
    nmodes::Int
end

"""
    MyPreconditionerPrimal(A::ExtendableSparseMatrix, bdofs, nmodes)

Construct a block-diagonal LU preconditioner from the mean matrix ``A_0``.

Stiffens boundary diagonal entries to ``1e60`` and precomputes the LU factorisation of
the modified matrix. The factorisation is reused for every preconditioner-vector
product in the GMRES iteration.
"""
function MyPreconditionerPrimal(A::ExtendableSparseMatrix{Tv, Ti}, bdofs, nmodes) where {Tv, Ti}
    DA::Array{Tv, 1} = zeros(Tv, size(A, 1))
    for j in 1:length(DA)
        DA[j] = A[j, j]
    end

    # compute LU factorisation of S
    for dof in bdofs
        A[dof, dof] = 1.0e60
    end
    flush!(A)
    LUA = lu(A.cscmatrix)

    # temporary storage array for solver
    temp = zeros(Tv, size(A, 1))

    return MyPreconditionerPrimal{Tv, typeof(LUA)}(LUA, DA, temp, bdofs, nmodes)
end

@inline LinearAlgebra.ldiv!(C::MyPreconditionerPrimal, b) = ldiv!(b, C, b)
@inline function LinearAlgebra.ldiv!(y, C::MyPreconditionerPrimal{Tv, FAC}, b) where {Tv, FAC}
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
@inline function LinearAlgebra.:\(C::MyPreconditionerPrimal, b)
    y = zero(b)
    ldiv!(y, C, b)
    return y
end


"""
    LinearAlgebra.mul!(Ax, S::MySystemPrimal, x)

Matrix-free matmul for the block-structured SG Poisson operator.

Computes the product ``Ax`` in-place without assembling the full block matrix. For each
stochastic mode ``\\mu``, the deterministic diffusion ``A_0`` is applied to mode ``\\mu``
of ``x``, and every KL perturbation ``A_e`` is applied to mode ``\\nu`` of ``x`` weighted
by the coupling coefficient ``G_{e,\\mu,\\nu}``. Only the nonzero couplings, precomputed
in `S.couplings` at construction time, are traversed. Boundary rows are zeroed out after
accumulation.
"""
function LinearAlgebra.mul!(Ax, S::MySystemPrimal{Tv, MT, VT, GT}, x) where {Tv, MT, VT, GT}
    fill!(Ax, 0)
    A0::MT = S.A0
    Am::Vector{MT} = S.Am
    vsize::Int = size(A0.entries, 1)
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
        addblock_matmul!(view(Ax, a:b), A0[1, 1], view(x, a:b))
    end

    # stochastic part
    for (mu, nu, e, g) in couplings
        a = (mu - 1) * vsize + 1
        b = mu * vsize
        a2 = (nu - 1) * vsize + 1
        b2 = nu * vsize
        addblock_matmul!(view(Ax, a:b), Am[e][1, 1], view(x, a2:b2); factor = g)
    end

    for mu::Int in 1:nmodes
        a = (mu - 1) * vsize + 1
        for dof in bdofs
            Ax[a + dof - 1] = 0
        end
    end
    return Ax
end

Base.eltype(::MySystemPrimal{Tv}) where {Tv} = Tv
Base.size(S::MySystemPrimal, d::Int) = S.nmodes * size(S.A0.entries, 1)


"""
    solve_primal!(SolutionSGFEM::SGFEVector, A0, Am, b0, G, nmodes; atol, rtol)

Solve the primal SG Poisson system iteratively using matrix-free preconditioned GMRES.

Builds a `MySystemPrimal` operator and a `MyPreconditionerPrimal` from
the mean matrix ``A_0``, KL perturbation blocks ``A_m``, coupling tensor ``G``, and
boundary info extracted from the ``SGFEVector``. The right-hand side is assembled by
copying the solution vector, adding the deterministic force block ``b_0`` to mode 1,
and zeroing boundary entries on all modes.

Krylov.jl's ``Krylov.gmres`` is used as the solver. Default tolerances are 1e-14.

# Arguments
- `SolutionSGFEM`: Output ``SGFEVector`` (modified in-place).
- `A0`: Mean stiffness block.
- `Am`: Vector of KL-perturbation stiffness blocks.
- `b0`: Deterministic right-hand side block (applied to mode 1 only).
- `G`: Coupling tensor ``G[e, mu, nu]``.
- `nmodes`: Number of stochastic modes.

# Keywords
- `atol`: Absolute GMRES tolerance (default 1.0e-14).
- `rtol`: Relative GMRES tolerance (default 1.0e-14).
"""
function solve_primal!(SolutionSGFEM::SGFEVector, A0, Am, b0, G, nmodes; atol = 1.0e-14, rtol = 1.0e-14)

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

    S = MySystemPrimal(A0, Am, G, bdofs, nmodes)
    @info "...initializing Preconditioner"
    @time P = MyPreconditionerPrimal(A0.entries, bdofs, nmodes)

    ## right-hand side
    b = deepcopy(SolutionSGFEM)
    addblock!(b[1], b0[1])
    for m in 1:nmodes
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


"""
    solve_full_primal!(SolutionSGFEM::SGFEVector, A0, A, b, G, nmodes)

Build the full block SG Poisson matrix and solve it with a direct backslash.

Assembles the complete block system into a ``FEMatrix`` with ``nmodes * nmodes`` blocks:
the diagonal receives contributions from the mean operator ``A_0``, and off-diagonal
and diagonal blocks are filled by the KL perturbations ``A_e`` weighted by the coupling
tensor ``G``. Boundary dofs are stiffened and RHS entries set to zero. The system is
solved in a single line via Julia's ``\\`` operator.

This function is primarily useful for verification against the matrix-free solver
``solve_primal!`` since it scales as ``O((nmodes * ndofs)^3)``.

# Arguments
- `SolutionSGFEM`: Output ``SGFEVector`` (modified in-place).
- `A0`: Mean stiffness block.
- `A`: Vector of KL-perturbation blocks.
- `b`: Right-hand side blocks.
- `G`: Coupling tensor.
- `nmodes`: Number of stochastic modes.
"""
function solve_full_primal!(SolutionSGFEM::SGFEVector, A0, A, b, G, nmodes)

    M::Int = length(A) # size(G,1) / nmodes

    FES = SolutionSGFEM.FES_space[1]
    @show FES
    multi_indices = SolutionSGFEM.TB.multi_indices
    bigFES = [FES for j in 1:nmodes]
    x::Vector{Float64} = zeros(Float64, 2)

    bigS = FEMatrix(bigFES)
    bigb = FEVector(bigFES)

    ## add AA and BB and N0
    for j in 1:nmodes
        addblock!(bigS[j, j], A0[1, 1])
    end

    ## add Nm blocks of NN
    g::Float64 = 0
    for j in 1:nmodes, k in 1:nmodes
        for e in 1:M
            g = G[(e - 1) * nmodes + j, k] # ⟨ ξ_m ψ_mi(j) ψ_mi(k) ⟩
            if abs(g) > 1.0e-14
                addblock!(bigS[j, k], A[e][1, 1]; factor = g)
            end
        end
    end

    ## right-hand side
    addblock!(bigb[1], b[1])

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
