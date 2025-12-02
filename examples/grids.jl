# # Example showing the different kind of Grids used.
# The choice of grids dictates the method of solution.

using ChaoticShearAlfvenWaves
using CSAWViz

# define the problem 
geo = init_geometry()
fields = init_fields()
prob = init_problem(geometry=geo, fields=fields)

# We can specify the radial grid to consider only the second derivative on flux surfaces
# This efficiently gives the continuum
ψgrid = init_grid(:cont, 100)
θgrid = init_grid(:sm, 2, start=1)
φgrid = init_grid(:sm, 1, start=-1)
grids = init_grids(ψgrid, θgrid, φgrid)

# For the continuum we do not need to specify the solver
# And we only get the eigenvalues
evals = compute_spectrum(prob, grids);

continuum_plot(evals)

# To compute the full spectrum we require a solver object
solver = init_solver(prob=prob, full_spectrum=true)

# We can use the Fourier spectral method in both poloidal and toroidal
# With finite elements in the radial dimension
ψgrid = init_grid(:ψ, 100)
θgrid = init_grid(:sm, 2, start=1)
φgrid = init_grid(:sm, 1, start=-1)
grids = init_grids(ψgrid, θgrid, φgrid)

evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);

continuum_plot(evals)

# Now we can look at individual solutions.
# For example, the TAE at ω≈0.3
ind = find_ind(evals, 0.3)

#then look at the harmonics
harmonic_plot(ϕft, grids, ind)

# Alternatively, we can use finite elements in ψ and θ, with the spectral method in φ.
# Noting that the matrices are slower to construct with finite elements
# So here we use a smaller grid
ψgrid = init_grid(:ψ, 40)
#specifying :θ here defaults the finite element grid to be 2π periodic.
#pf (phase factor) centres the Fourier transform on a specific mode
#giving higher accuracy
θgrid = init_grid(:θ, 10, pf=1)
φgrid = init_grid(:sm, 1, start=-1)
grids = init_grids(ψgrid, θgrid, φgrid)

evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);

continuum_plot(evals)

#look for the same TAE.
ind = find_ind(evals, 0.3)

#then look at the harmonics
harmonic_plot(ϕft, grids, ind)


# Finally, we can use finite elements in all dimensions.
# This is a bit impractical in serial (see CSAWParallel)
# But is most appropriate for cases with broken symmetry
ψgrid = init_grid(:ψ, 30)
θgrid = init_grid(:θ, 5, pf=1)
φgrid = init_grid(:φ, 1, pf=-1)
grids = init_grids(ψgrid, θgrid, φgrid)

evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);

continuum_plot(evals)

#look for the same TAE.
ind = find_ind(evals, 0.3)

#then look at the harmonics
harmonic_plot(ϕft, grids, ind)
