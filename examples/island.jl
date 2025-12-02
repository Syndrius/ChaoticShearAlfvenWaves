# Example using Straight field line island coordinates.
# See Qu and Hole 2023, Konies et al 2024 and Thomas et al 2026.

using ChaoticShearAlfvenWaves
using CSAWViz

# geometry is set to island with :κ
geo = init_geometry(:κ, R0=1.0)
# we form the island, specifying the radial coordinate, :κ
# this requires that we also set ψ0, the radial location and qp the derivative of the qprofile at ψ0
#these values allow the q-profile required for analytical coordinates to be computed.
isl = init_island(:κ, m0=1, n0=-1, w=0.1, ψ0=0.5, qp=1.0)

#specify the radial variable :κ, with the specific island q-profile.
fields = init_fields(:κ, q=island_q, isl=isl)
prob = init_problem(geometry=geo, fields=fields)

#radial grid is set to island, :κ, to specify default values
#specifically sets the boundary at 0.999 to avoid the sepratrix.
κgrid = init_grid(:κ, 20)
# we use finite elements in every dimension for mapping below.
ᾱgrid = init_grid(:ᾱ, 5, pf=1)
τgrid = init_grid(:τ, 1, pf=0)
grids = init_grids(κgrid, ᾱgrid, τgrid)

#solve for the full spectrum
solver = init_solver(prob=prob, full_spectrum=true)

#we set deriv=true so that we also return the derivative of ϕ
#which is needed for mapping.
evals, ϕ, ϕft = compute_spectrum(prob, grids, solver, deriv=true);

# we can see the island structure
continuum_plot(evals, ylimits=(0.0, 0.05))

# and a MiAE at 0.015
ind = find_ind(evals, 0.015)

# We can see the global structure.
harmonic_plot(ϕft, grids, ind)

# We can also map the solutions from the island coordinates back to toroidal coordinates.

# This requires that we create a new toroidal grid
#we restrict it to only the island region.
ψgrid = init_grid(:ψ, 80, start=0.45, stop=0.55)
θgrid = init_grid(:θ, 20)
φgrid = init_grid(:φ, 10)
tor_grids = init_grids(ψgrid, θgrid, φgrid)

# Mapping the spectrum
evals_tor, ϕ_tor, ϕ_torft = isl_spectrum_to_tor(evals, ϕ, ϕft, grids, tor_grids, prob);

#now we can plot the same MiAE solution
#but in toroidal coordinates
#we shift x3=φ to centre the island mode.
contour_plot(ϕ_tor, tor_grids, ind, x3=6)
