module ExtendableASGFEM

using Distributed: Distributed
using Distributions: Distributions, Normal, Uniform, dim, pdf
using DocStringExtensions: DocStringExtensions, TYPEDEF, TYPEDSIGNATURES
using ExtendableFEM: ExtendableFEM, BilinearOperator, FaceInterpolator,
    HomogeneousBoundaryData, ItemIntegrator, ItemIntegratorDG,
    L2NormIntegrator, LinearOperator, ProblemDescription,
    Unknown, assemble!, assign_operator!, assign_unknown!,
    div, grad, id, jump, plot
using ExtendableFEMBase: ExtendableFEMBase, BFaceDofs, CellDofs, Divergence,
    FEEvaluator, FEMatrix, FESpace, FEVector,
    FEVectorBlock, Gradient, H1P1, H1Pk, HDIVRTk,
    Laplacian, QuadratureRule, _addnz, addblock!,
    addblock_matmul!, eval_febe!, fill!, get_ndofs,
    get_polynomialorder, unicode_scalarplot,
    update_basis!
using ExtendableGrids: ExtendableGrids, Adjacency, BFaceFaces, CellFaces,
    CellNodes, CellVolumes, Coordinates, ExtendableGrid,
    FaceNodes, FaceVolumes, L2GTransformer, NodePatchGroups,
    ON_CELLS, ON_FACES, ON_IFACES,
    SerialVariableTargetAdjacency, Triangle2D,
    UniqueCellGeometries, VariableTargetAdjacency, append!,
    atranspose, eval_trafo!, interpolate!,
    max_num_targets_per_source, num_cells, num_nodes,
    num_sources, num_targets, unique, update_trafo!
using ExtendableSparse: ExtendableSparse, ExtendableSparseMatrix,
    LUFactorization, flush!
using GridVisualize: GridVisualize, GridVisualizer, scalarplot, scalarplot!
using IterativeSolvers: IterativeSolvers
using Krylov: Krylov
using LinearAlgebra: LinearAlgebra, SymTridiagonal, dot, eigvals, eigvecs,
    ldiv!, mul!, norm
using Printf: Printf, @printf
using Random: Random, rand!
using SparseArrays: SparseArrays, SparseMatrixCSC, nzrange, rowvals
using SpecialFunctions: SpecialFunctions, zeta
using StaticArrays: StaticArrays, SMatrix, SVector

include("mopcontrol.jl")
export generate_multiindices, add_boundary_modes, classify_modes
export get_neighbourhood_relation_matrix

include("orthogonal_polynomials/orthogonal_polynomials.jl")
export OrthogonalPolynomialType
export HermitePolynomials, LegendrePolynomials #, LaguerrePolynomials, ChebyshevTPolynomials, ChebyshevUPolynomials
export evaluate, evaluate!, gauss_rule
export norm, norms
export distribution
export normalise_recurrence_coefficients

include("onbasis.jl")
export ONBasis
export evaluate
export vals4qp, vals4poly, qw, qp, norm4poly
export triple_product, scalar_product, integral
export set_sample!
export get_multiindex, prepare_multi_indices!
export get_coupling_coefficient

include("tensorizedbasis.jl")
export TensorizedBasis
export num_multiindices, maxlength_multiindices
export get_tensor_multiplication_with_ym
export sample_distribution

include("sgfevector.jl")
export SGFEVector
export set_sample!

include("coefficients/coefficients.jl")
export AbstractStochasticCoefficient
export StochasticCoefficientCosinus, SingleStochasticCoefficient
export get_am!, get_gradam!
export get_a!, get_expa!
export get_am_x, get_gradam_x_sigma, get_gradam_x_u
export expa_PCE_mop, maxm

include("plots.jl")
export plot_modes
export plot_basis

include("modelproblems/modelproblems.jl")
export solve!
export deterministic_problem, deterministic_problem2
export LogTransformedPoissonProblemPrimal
export LogTransformedPoissonProblemDual
export PoissonProblemPrimal

include("sampling_error.jl")
export calculate_sampling_error

include("estimate.jl")
export estimate
export estimate_equilibration
export get_next_tail

greet() = print("Hello World!")

end # module
