using ExtendableASGFEM
using LinearAlgebra
using DoubleFloats
using ExplicitImports
using Aqua
using Test

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
    
    ## configure test precision
    T = Double64 # higher precision instead of Float64
    normalize = true

    for T in (Float64, Double64)
        if normalize
            tolerance = T == Float64 ? 1e-12 : 1e-12
        else
            tolerance = T == Float64 ? 1e-8 : 1e-19
        end
        val::T = 0.0
        val_compare::T = 0.0

        #########################
        ## Hermite polynomials ##
        #########################

        ## generate ONB
        OBType = HermitePolynomials
        order = T == Float64 ? 6 : 8
        ONB = ONBasis(OBType, order, 3*order; T = T)

        @info "Checking integrals..."
            for j = 0 : order
                val = integral(ONB, j; normalize = normalize)
                val_compare = j == 0 ? 1 : 0
                @info "$([j]) = $val | $(val_compare) | error = $(val-val_compare)"
                @test abs(val - val_compare) < tolerance
            end

        @info "Checking scalar products (orthonormality)..."

            ## test if (quadratured) norm of k-th Hermite polynomials equals hard-coded norm
            for j = 0 : order
                val = norm4poly(ONB, j)
                val_compare = norms(OBType, j)
                @info "$([j]) = $(val_compare) | error = $(val-val_compare)"
                @test abs(val - val_compare) < tolerance
            end

            ## check ortho-normality
            for j = 0 : order, k = j : order
                val = scalar_product(ONB, j, k; normalize = normalize)
                val_compare = j == k ? (normalize ? 1 : norms(OBType, j) * norms(OBType, k)) : 0
                @test abs(val - val_compare) < tolerance
                if j == k
                    @info "$([j,k]) = $(val_compare) | error = $(val-val_compare)"
                end
            end

        ## check triple products
        @info "Checking triple products..."
            for j = 0 : order, k = j : order, l = k : order
                val = triple_product(ONB, j, k, l; normalize = normalize)
                if iseven(j + k + l) && minimum([j+k-l,k+l-j,l+j-k]) >= 0
                    val_compare = factorial(l) * factorial(j) * factorial(k) / (factorial(Int((j+k-l)/2)) * factorial(Int((k+l-j)/2)) * factorial(Int((j+l-k)/2)))
                    if normalize
                        val_compare /= norms(OBType, j) * norms(OBType, k) * norms(OBType, l)
                    end
                    @info "$([j,k,l]) = $(val_compare) | error = $(val-val_compare)"
                    @test abs(val - val_compare) < tolerance
                else
                    @test abs(val) < tolerance
                end
            end

        ## check triple products
        @info "Checking allocations..."
            ## test that evaluation of polynomials is allocation-free
            evaluate(ONB, 0.5) # run once for compilation
            a = @allocated evaluate(ONB, 0.5) #now measure allocation
            @test a == 0 
    end


    tolerance = 1e-12

    OBType = HermitePolynomials
    TB = TensorizedBasis(OBType, 3, 3, 6)
    nmi = num_multiindices(TB)
    multi_indices = TB.multi_indices

    @info "Checking triple products..."
    for j = 1: nmi, k = 1 : nmi, l = 1 : nmi
        val = triple_product(TB, j, k, l; normalize = true)
        val_compare = 1
        for d = 1 : maxlength_multiindices(TB)
            a = multi_indices[j][d]
            b = multi_indices[k][d]
            c = multi_indices[l][d]
            if iseven(a + b + c) && minimum([a+b-c,b+c-a,a+c-b]) >= 0
                val_compare *= factorial(c) * factorial(a) * factorial(b) / (factorial(Int((a+b-c)/2)) * factorial(Int((b+c-a)/2)) * factorial(Int((a+c-b)/2))) / (norms(OBType, a) * norms(OBType, b) * norms(OBType, c))
            else
                val_compare = 0
            end
        end
        @info "$([multi_indices[j],multi_indices[k],multi_indices[l]]) = $(val_compare) | $(val) | error = $(val-val_compare)"
        @test abs(val - val_compare) < tolerance
    end


end

main()
