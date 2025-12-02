
# # Example using Quadratic Flux Minimising (QFM) coordinates

using ChaoticShearAlfvenWaves
using CSAWViz

# QFM surfaces are generated in the companion package CSAWCantori.
# We have some precomputed used for testing.
surf_dir = abspath(joinpath(pathof(ChaoticShearAlfvenWaves), "../../test/data/"))
surfs = surfaces_from_file(joinpath(surf_dir, "benchmark_surfaces.jld2"));


# This case uses a very large non-resonant perturbation
isl = init_island(m0=3, n0=1, A=0.1)
geo = init_geometry()
fields = init_fields(:ψ, q=cantori_q, isl=isl)
prob = init_problem(geometry=geo, fields=fields)

#we can view the poincare section to understand this perturbation
poincare_plot(prob, 500, LinRange(0.15, 0.9, 50), zeros(50), color=:black)
# and overlay the qfm surfaces
plot_surfaces!(surfs, color=:red, linewidth=1.0)

# We use finite elements for mapping below.
# The unrealistic perturbation causes issues at the boundaries so we shrink the size.
sgrid = init_grid(:s, 30, start=0.15, stop=0.9)
ϑgrid = init_grid(:ϑ, 5, pf=1)
ζgrid = init_grid(:ζ, 1, pf=-1)

grids = init_grids(sgrid, ϑgrid, ζgrid)

#solve over the full spectrum
solver = init_solver(prob=prob, full_spectrum=true)
#we now also need to pass the QFM surfaces in.
evals, ϕ, ϕft = compute_spectrum(prob, grids, solver, surfs, deriv=true);

#despite the perturbation the continuum appears nowrmal
#as we have 'straightened' the flux surfaces.
continuum_plot(evals)

#we look at specific solutions
#a TAE and a typical continuum solution
tae_ind = find_ind(evals, 0.25) 
cont_ind = find_ind(evals, 0.20)
harmonic_plot(ϕft, grids, tae_ind)
harmonic_plot(ϕft, grids, cont_ind)

# To understand the QFM coordinates, we can map the solutions back to toroidal coordinates
# Note that this mapping can be slow for large grids.
# The radial grid is reduced again to prevent issues with interpolation at the edges.
ψgrid = init_grid(:ψ, 80, start=0.25, stop=0.8)
θgrid = init_grid(:θ, 30)
φgrid = init_grid(:φ, 10)
tor_grids = init_grids(ψgrid, θgrid, φgrid)

# Map the solutions
tor_evals, ϕ_tor, ϕft_tor = qfm_spectrum_to_tor(evals, ϕ, ϕft, grids, tor_grids, surfs);

#harmonics show the more complicated structure in toroidal coordinates.
harmonic_plot(ϕft_tor, tor_grids, tae_ind, label_max=0.5)
harmonic_plot(ϕft_tor, tor_grids, cont_ind, label_max=0.5)

#Looking at the continuum solution, we can see the solution conforms to the modified flux surfaces
contour_plot(ϕ_tor, tor_grids, cont_ind)
