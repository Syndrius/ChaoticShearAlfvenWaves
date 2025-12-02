# # Example with a chaotic magnetic field using QFM coordinates
# See Thomas et al. 2026

using ChaoticShearAlfvenWaves
using CSAWViz

# QFM surfaces are generated in the companion package CSAWCantori.
# We have some precomputed used for testing.
surf_dir = abspath(joinpath(pathof(ChaoticShearAlfvenWaves), "../../test/data/"))
surfs = surfaces_from_file(joinpath(surf_dir, "chaos_surfaces.jld2"));

# Chaotic trajectories are generated with overlapping islands.
# Surfaces are generated for specific problem setup below.
k = 0.0012
geo = init_geometry(:cyl, R0=1.0)
isl1 = init_island(:ψ, m0=3, n0=-2, A=k/3)
isl2 = init_island(:ψ, m0=4, n0=-3, A=k/4)

#Multiple islands are input as a list
isls = [isl1, isl2]
fields = init_fields(:ψ, q=cantori_q, isls=isls)
prob = init_problem(geometry=geo, fields=fields)

#we can view the poincare section to understand this perturbation
poincare_plot(prob, 500, LinRange(0.4, 0.8, 100), zeros(100), color=:black, ylimits=(0.5, 0.667))
# and overlay the qfm surfaces
plot_surfaces!(surfs, color=:red, linewidth=1.0)

sgrid = init_grid(:s, 100, start=0.15, stop=0.9, sep1=0.5, sep2=0.667, frac=0.5)
ϑgrid = init_grid(:sm, 3, start=1)
ζgrid = init_grid(:sm, 2, start=-2)

grids = init_grids(sgrid, ϑgrid, ζgrid)

# we solve with multiple shift and invert solves, slicing up the spectrum.
solver = init_solver(prob=prob, targets=[0.15, 0.25])
#we now also need to pass the QFM surfaces in.
evals, ϕ, ϕft = compute_spectrum(prob, grids, solver, surfs, deriv=true);

# Continuum shows typical structure, excpet at the islands at s=0.5 and s=0.667
continuum_plot(evals)

#looking at a specific solution
ind = find_ind(evals, 0.147)
# we can see the standard continuum structure
harmonic_plot(ϕft, grids, ind)

# Further chaotic cases require significantly higher resolution (see CSAWParallel)
# These cases are discussed in Thomas et al. 2026
