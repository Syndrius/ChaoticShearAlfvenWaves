# Example showing cylindrical geometry.

using ChaoticShearAlfvenWaves
using CSAWViz

#the problem is initiated with cylindrical geoemtry
geo = init_geometry(:cyl)
fields = init_fields()
prob = init_problem(geometry=geo, fields=fields)

# Creating grids for computing the continuum
ψgrid = init_grid(:cont, 100)
θgrid = init_grid(:sm, 2, start=1)
φgrid = init_grid(:sm, 1, start=-1)
grids = init_grids(ψgrid, θgrid, φgrid)

# Cylindrical geometry does not have mode coupling 
# So we do not see a gap.
evals = compute_spectrum(prob, grids);

continuum_plot(evals)

# With cylindrical geoemtry, we can compute the continuum analytically
evals_an = analytical_spectrum(prob, grids);

#overlaying the analytical continuum we can see the agreement.
continuum_plot!(evals_an, markersize=2.0)

# We can also find a GAE in cylindrical geometry.
# This follows the problem outined in Van Rij et al. 1985.

# Fields are set with the geometric radius as the radial coordinate
# with specific density profile to create a minimum in the continuum.
fields = init_fields(:r, q=gae_q, dens=gae_dens)
prob = init_problem(geometry=geo, fields=fields)
solver = init_solver(prob=prob, full_spectrum=true)

#the density and q-profiles can be viewed.
dens_plot(gae_dens, grids)
q_plot(gae_q, grids)

#create the grids, using the spectral method in poloidal and toroidal dimensions.
rgrid = init_grid(:r, 100)
θgrid = init_grid(:sm, 1, start=1)
ζgrid = init_grid(:sm, 1, start=0)
grids = init_grids(rgrid, θgrid, ζgrid)

# solve the problem
evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);

#looking at the continuum
continuum_plot(evals)

#we can see a GAE below the continuum.
ind = find_ind(evals, 0.6)

#looking closer, we verify the global structure.
harmonic_plot(ϕft, grids, ind)

