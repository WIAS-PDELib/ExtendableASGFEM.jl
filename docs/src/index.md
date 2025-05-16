# ExtendableASGFEM

This package implements the stochastic Galerkin finite element method for certain two dimensional
model problems involving KLE of stochastic coefficients. The rather huge systems have a tensorized
structure and are solved by iterative solvers. A posteriori error estimators steer the spatial
and stochastic refinement.

The spatial discretization is based on the finite element package [ExtendableFEM.jl](https://github.com/WIAS-PDELib/ExtendableFEM.jl) and [ExtendableFEMBase.jl](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl)


## References

- [1]   "Adaptive stochastic Galerkin FEM"\
        M. Eigel, C.J. Gittelson, C. Schwab, E. Zander\
        CMAME 270, 1 (2014), 247--269\
        [>Journal-Link<](https://www.doi.org/10.1016/J.CMA.2013.11.015)
        [>Preprint-Link<](https://www.research-collection.ethz.ch/handle/20.500.11850/154941)
- [2]   "A posteriori error control for stochastic Galerkin FEM with high-dimensional random parametric PDEs",\
        M. Eigel, C. Merdon\
        to be published in: Error Control, Adaptive Discretizations, and Applications, Part 3, Academic Press\
        [>Preprint-Link<](https://www.wias-berlin.de/publications/wias-publ/run.jsp?template=abstract&type=Preprint&year=&number=3174)
- [3]   "Local equilibration error estimators for guaranteed error control in adaptive higher-order stochastic  Galerkin finite element methods"\
        M. Eigel and C. Merdon\
        SIAM/ASA J. Uncertainty Quantification 4(1) (2016), 1372--1397"\
        [>Journal-Link<](https://epubs.siam.org/doi/10.1137/15M102188X)
        [>Preprint-Link<](http://www.wias-berlin.de/publications/wias-publ/run.jsp?template=abstract&type=Preprint&year=2014&number=1997)