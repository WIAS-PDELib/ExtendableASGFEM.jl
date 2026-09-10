using ExtendableASGFEM
using LinearAlgebra
using DoubleFloats
using ExplicitImports
using Aqua
using Test

include("../scripts/poisson_simple.jl")

function main()

    @testset "Aqua.jl" begin
        Aqua.test_all(
            ExtendableASGFEM;
            ambiguities = false,
        )
        #Aqua.test_ambiguities(ExtendableASGFEM)
    end

    @testset "ExplicitImports" begin
        @test ExplicitImports.check_no_implicit_imports(ExtendableASGFEM) === nothing
        @test ExplicitImports.check_no_stale_explicit_imports(ExtendableASGFEM) === nothing
    end

    if isdefined(Docs, :undocumented_names) # >=1.11
        #@testset "UndocumentedNames" begin
        #    @test isempty(Docs.undocumented_names(ExtendableASGFEM))
        #end
    end

    T = Double64
    tol = 1e-15

    @testset "norm_basis" begin
        # Hermite: ||H_k|| = sqrt(k!)
        @test norm_basis(HermitePolynomials, 0) ≈ 1.0
        @test norm_basis(HermitePolynomials, 1) ≈ 1.0
        @test norm_basis(HermitePolynomials, 2) ≈ sqrt(2)
        @test norm_basis(HermitePolynomials, 3) ≈ sqrt(6)
        @test norm_basis(HermitePolynomials, 4) ≈ sqrt(24)
        @test norm_basis(HermitePolynomials, 5) ≈ sqrt(120)

        # Legendre: ||P_k|| = sqrt(1/(2k+1))
        @test norm_basis(LegendrePolynomials, 0) ≈ 1.0
        @test norm_basis(LegendrePolynomials, 1) ≈ sqrt(1/3)
        @test norm_basis(LegendrePolynomials, 2) ≈ sqrt(1/5)
        @test norm_basis(LegendrePolynomials, 3) ≈ sqrt(1/7)
    end

    @testset "evaluate recurrence -- Hermite" begin
        # Known polynomial values at specific points
        # H_0=1, H_1=y, H_2=y^2-1, H_3=y^3-3y, H_4=y^4-6y^2+3, H_5=y^5-10y^3+15y
        x = T(0.5)
        n = 5
        vals = evaluate(HermitePolynomials, n, x)

        exact = [T(1),
                 x,
                 x^2 - T(1),
                 x^3 - 3x,
                 x^4 - 6x^2 + T(3),
                 x^5 - 10x^3 + 15x]
        @test all(abs.(vals .- exact) .< tol)

        # H_n(0): odd → 0, even → (-1)^{k/2}*(k-1)!!
        x0 = zero(T)
        vals0 = evaluate(HermitePolynomials, 5, x0)
        @test vals0[1] ≈ 1    # H_0
        @test vals0[2] ≈ 0    # H_1(0)
        @test vals0[3] ≈ -1   # H_2(0)
        @test vals0[4] ≈ 0    # H_3(0)
        @test vals0[5] ≈ 3    # H_4(0)
        @test vals0[6] ≈ 0    # H_5(0)
    end

    @testset "evaluate recurrence -- Legendre" begin
        # Known polynomial values at y=0.5
        # P_0=1, P_1=y, P_2=(3y^2-1)/2, P_3=(5y^3-3y)/2, P_4=(35y^4-30y^2+3)/8, P_5=(63y^5-70y^3+15y)/8
        x = T(0.5)
        n = 5
        vals = evaluate(LegendrePolynomials, n, x)

        exact = [T(1),
                 x,
                 (3x^2 - T(1)) / T(2),
                 (5x^3 - 3x) / T(2),
                 (35x^4 - 30x^2 + T(3)) / T(8),
                 (63x^5 - 70x^3 + 15x) / T(8)]
        @test all(abs.(vals .- exact) .< tol)

        # Orthogonality: P_n(0) for even/odd
        x0 = zero(T)
        vals0 = evaluate(LegendrePolynomials, 5, x0)
        @test vals0[1] ≈ 1       # P_0(0)=1
        @test vals0[2] ≈ 0       # P_1(0)=0
        @test vals0[3] ≈ -0.5    # P_2(0)=-1/2
        @test vals0[4] ≈ 0       # P_3(0)=0
    end

    @testset "evaluate array" begin
        xs = [T(-0.8), T(-0.4), T(0.0), T(0.4), T(0.8)]
        vals = evaluate(HermitePolynomials, 3, xs)
        @test size(vals) == (length(xs), 4)

        # Check P_0 is all ones
        @test all(abs.(vals[:, 1] .- 1) .< tol)

        # Check P_1 = x
        @test all(abs.(vals[:, 2] .- xs) .< tol)

        # Check P_2 = x^2 - 1
        @test all(abs.(vals[:, 3] .- (xs .^ 2 .- 1)) .< tol)
    end

    @testset "gauss_rule: nodes and weights" begin
        # Gauss-Legendre(n=3) with 6 nodes (n=2*3 in gauss_rule call)
        gr = gauss_rule(LegendrePolynomials, 6; T = T)
        nodes, weights = gr

        @test length(nodes) == 6
        @test length(weights) == 6
        @test all(>(zero(T)), weights)
        @test abs(sum(weights) - 1.0) < 1e-10  # integral of w=1/2 on [-1,1] is 2
        @test sum(nodes) ≈ zero(T) atol = 1e-10  # symmetric about 0

        # Gauss-Hermre with 6 nodes, weights sum to sqrt(pi) for standard Gauss-Hermite.
        # For normal weight it should sum to 1.
        gr_h = gauss_rule(HermitePolynomials, 6; T = T)
        wh, weightsh = gr_h
        @test length(wh) == 6
        @test abs(sum(weightsh) - 1.0) < 1e-10  # probability measure sums to 1
        @test sum(wh) ≈ zero(T) atol = 1e-10  # centred about 0
    end

    @testset "evaluate against gauss quadrature" begin
        # For both bases, integrate p_n^2 using Gauss quadrature and compare to norm_basis(n)^2
        for (basis, maxn) in [(HermitePolynomials, 6), (LegendrePolynomials, 8)]
            gr = gauss_rule(basis, 2 * maxn; T = T)
            nodes, weights = gr

            for n in 0:maxn
                pn = evaluate(basis, maxn, nodes)[:, n+1]  # column n+1 (Julia 1-based)
                quad_norm2 = sum(pn .^ 2 .* weights)
                expected2 = norm_basis(basis, n)^2
                @test abs(quad_norm2 - expected2) < tol * maxn
            end

            # Orthogonality check: integral of p_n * p_m for n ≠ m should be 0
            for n in 0:(maxn-1), m in (n+1):maxn
                pn = evaluate(basis, maxn, nodes)[:, n+1]
                pm = evaluate(basis, maxn, nodes)[:, m+1]
                integrand = sum(pn .* pm .* weights)
                @test abs(integrand) < tol
            end
        end
    end

    @testset "poisson_simple integration" begin
        sol = PoissonSimple.main(
            nrefs = 1,
            order = 1,
            decay = 2.0,
            mean = 1.0,
            domain = "square",
            initial_modes = [[0], [1, 0], [0, 1]],
            use_iterative_solver = false,
            Plotter = nothing,
        )
        @test !isnothing(sol)
    end

    return
end

main()
