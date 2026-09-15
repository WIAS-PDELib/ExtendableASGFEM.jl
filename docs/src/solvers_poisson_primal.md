# Iterative Solution of the Primal Poisson Problem

This page describes how the primal stochastic Poisson system is solved using a matrix-free preconditioned GMRES method. A complementary direct solver is also available for verification purposes.

## The Block system

The stochastic Galerkin discretisation of the Poisson equation with a Karhunen-Loève expanded diffusion coefficient yields a block-structured linear system

```math
\mathbf{A}\,\mathbf{u} = \mathbf{b} \,,
```

where the unknown and right-hand side vectors are partitioned by stochastic mode:

```math
\mathbf{u} = \begin{bmatrix} u_1 \\ u_2 \\ \vdots \\ u_{n_{\text{modes}}} \end{bmatrix}, \quad
\mathbf{b} = \begin{bmatrix} b_1 \\ 0 \\ \vdots \\ 0 \end{bmatrix}.
```

The matrix $\mathbf{A}$ is composed of $n_{\text{modes}} \times n_{\text{modes}}$ blocks $\mathbf{A}_{\mu,\nu}$, each of size $n_{\text{dofs}} \times n_{\text{dofs}}$. Block $(\mu,\nu)$ is given by

```math
\mathbf{A}_{\mu,\nu} = A_0 + \sum_{e=1}^{M} G_{e,\mu,\nu}\, A_e \,,
```

where $A_0$ is the mean stiffness matrix, $A_e$ are the KL-perturbation matrices, and $G_{e,\mu,\nu} = \langle \phi_e H_\mu H_\nu \rangle$ are the coupling coefficients obtained by integrating the stochastic basis functions against the $e$-th KL eigenfunction.

## Matrix-free approach

Assembling $\mathbf{A}$ explicitly is prohibitively expensive for all but the smallest problems. Instead, the package provides the `MySystemPrimal` type, which implements `LinearAlgebra.mul!` for in-place matrix-vector products without storing the full system matrix.

```julia
mul!(Ax, S::MySystemPrimal, x)
```

For each stochastic mode $\mu$, the deterministic diffusion $A_0$ is applied to mode $\mu$ of $x$, and every KL perturbation $A_e$ is applied to mode $\nu$ of $x$ weighted by the coupling coefficient $G_{e,\mu,\nu}$. Boundary rows are zeroed after accumulation. The overall cost per matmul is $O(n_{\text{modes}}^2 \cdot M \cdot n_{\text{dofs}}^2)$ FLOPs, but avoids storing the full $n_{\text{modes}} n_{\text{dofs}} \times n_{\text{modes}} n_{\text{dofs}}$ matrix.

## Preconditioner

The preconditioner `MyPreconditionerPrimal` is block-diagonal: each diagonal block is the inverse of the mean stiffness matrix $A_0$. The construction:

1. Stiffens diagonal entries at boundary dofs to `1e60`, enforcing homogeneous Dirichlet conditions implicitly.
2. Computes an LU factorisation of the modified $A_0$.
3. During a preconditioner–vector product, the LU solve is applied independently to each stochastic mode block.

Because the factorisation of $A_0$ is computed only once, each preconditioner application costs $O(n_{\text{modes}} \cdot n_{\text{dofs}}^2)$ (forward/backward substitution per block).

## Iterative solve: `solve_primal!`

The entry point `solve_primal!` handles the full solve:

```julia
solve_primal!(SolutionSGFEM::SGFEVector, A0, Am, b0, G, nmodes, bfac; atol=1e-14, rtol=1e-14)
```

1. Extract boundary information from the `SGFEVector`'s finite element space.
2. Build the matrix-free operator and preconditioner.
3. Assemble the right-hand side $b$ by adding the deterministic force block $b_0$ to mode 1 and zeroing boundary entries on all modes.
4. Solve with Krylov.jl's GMRES: `Krylov.gmres(S, b, x; ldiv=true, M=P)`, where $P$ is the block-diagonal preconditioner.
5. Report the residual norm $\|\mathbf{Ax} - \mathbf{b}\|_2$.

The default tolerances are `1e-14`.

## Direct solver: `solve_full_primal!`

The `solve_full_primal!` function assembles the full block matrix into a `FEMatrix` and solves it with Julia's backslash operator. This is primarily a verification tool:

```julia
solve_full_primal!(SolutionSGFEM::SGFEVector, A0, A, b, G, nmodes, rhsfac)
```

The full assembly has storage and computational cost of $O((nmodes \cdot ndofs)^2)$ and $O((nmodes \cdot ndofs)^3)$, respectively, so it is only practical for very small problems.

## Summary

| Type / Function | Purpose |
|---|---|
| `MySystemPrimal` | Matrix-free block operator for matvec products |
| `MyPreconditionerPrimal` | Block-diagonal preconditioner (LU on $A_0$) |
| `solve_primal!` | Matrix-free preconditioned GMRES |
| `solve_full_primal!` | Direct solver on assembled full block matrix |
