"""
$(TYPEDSIGNATURES)

Plot scalar plots of the stochastic modes of an `SGFEVector` solution using `ExtendableFEM.plots.GridVisualize`.

# Arguments
- `sol::SGFEVector`: The stochastic Galerkin finite element solution vector whose modes are to be visualized.
- `unknown`: (Optional, default: 1) Index of the unknown to plot (for multi-unknown systems).
- `Plotter`: (Optional) Plotting backend to use (e.g., `GLMakie`, `CairoMakie`, `PyPlot`, `Plots`).
- `nmodes`: (Optional) Number of stochastic modes to plot (default: all modes in the tensorized basis).
- `ncols`: (Optional, default: 3) Number of columns in the plot grid.
- `width`: (Optional) Total width of the plot grid in pixels (default: `400 * ncols`).
- `sort`: (Optional, default: `false`) If `true`, modes are sorted by their L2 norm (largest first).

# Example
```julia
plot_modes(sol; Plotter=GLMakie, ncols=4, sort=true)
```
"""
function plot_modes(sol::SGFEVector; unknown = 1, Plotter = nothing, nmodes = num_multiindices(sol.TB), ncols = 3, width = 400 * ncols, sort = false)
    nrows = Int(ceil(nmodes / ncols))
    FES = sol.FES_space
    nunknowns = length(FES)
    TensorBasis = sol.TB
    norm4modes = norms(sol)
    nmodes = TensorBasis.nmodes
    if sort && nunknowns == 1
        modes = sortperm(norm4modes; rev = true)
        modes = modes[1:nmodes]
    else
        modes = Array{Int, 1}(1:nmodes)
    end
    for u in 2:nunknowns
        append!(modes, (u - 1) * nmodes .+ modes)
    end
    p = plot([id(j) for j in modes], [sol[j] for j in modes]; Plotter = Plotter)
    return p
end


function plot_basis(TB::TensorizedBasis; kwargs...)
    return plot_basis(TB.ONB; kwargs...)
end

"""
$(TYPEDSIGNATURES)

plots the basis functions of the ONBasis `ONB` via GridVisualize.
The `Plotter` argument determines the backend (e.g. GLMakie, CairoMakie, PyPlot, Plots).
"""
function plot_basis(ONB::ONBasis{T, OBT, npoly, nquad}; Plotter = nothing, resolution = (500, 500), kwargs...) where {T, OBT, npoly, nquad}
    vals4xref = ONB.vals4xref
    gr = ONB.gauss_rule
    xgrid = ExtendableGrid{Float64, Int}()
    xgrid[Coordinates] = reshape(gr[1], (1, length(gr[1])))
    xgrid[CellNodes] = [1:(nquad - 1) 2:nquad]
    p = GridVisualizer(; Plotter = Plotter, layout = (1, 1), clear = true, resolution = resolution, kwargs...)
    for j in 1:npoly
        scalarplot!(p[1, 1], xgrid, vals4xref; clear = false, Plotter = Plotter, label = "$j", xlimits = (-1, 1), limits = (-4, 4))
    end
    return
end
