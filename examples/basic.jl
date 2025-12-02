# # Example of Basic usage.

using ChaoticShearAlfvenWaves

# This main function, compute_spectrum, requires 3 inputs.

# The first defines the problem to be solved 
# This requires information about the geometry and the fields.

#by default, the geometry is initialised as toroidal with R0=4.0
geo = init_geometry()

#by default, the fields are intialised with a quadratic q-profile
#and a uniform density using the toroidal flux, ψ, as the radial coordinate
fields = init_fields()

#we can then create the problem from these.
prob = init_problem(geometry=geo, fields=fields)

# We next define the grids used.

# These are defined in each dimension separately, then combined into a 3d grid.
# The type of grids constructed dictate the method of solution.

#the radial grid, :ψ specifies the type, and 100 is the number of points.
ψgrid = init_grid(:ψ, 100)
#:sm specifies that the Fourier Spectral Method is used in the poloidal direction
#with 2 modes, m=1 and m=2.
θgrid = init_grid(:sm, 2, start=1)
#n=-1
φgrid = init_grid(:sm, 1, start=-1)

#we then combine the grids
grids = init_grids(ψgrid, θgrid, φgrid)

# Finally, we specify the solver to use.

#the solver is set to compute the entire spectrum
#which is only practical for small scale problems.
#the prob input is needed here for normalisation.
solver = init_solver(prob=prob, full_spectrum=true)

# We now compute the spectrum.

#this returns an object storing the eigenvalue information
#the perturbation to the electric potential
#and the Fourier transform in θ and φ of this potential.
evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);


# We have created a companion package, CSAWViz, specifcally for viewing these outputs.
using CSAWViz

#first the continuum
continuum_plot(evals)

#we can then find a specific solution to look at.
#for example, the TAE located at ω≈0.3
ind = find_ind(evals, 0.3)

#then look at the harmonics
harmonic_plot(ϕft, grids, ind)

#or a contour
contour_plot(ϕ, grids, ind)

#or a surface plot
surface_plot(ϕ, grids, ind)

#Now the TAE frequncy is known, we can change the solver to target this specific frequency.
#This allows us to increase the resolution and solve much larger matrices.

#sets the solver to perform shift and inert, getting the 100 solutions closets to 0.3.
solver = init_solver(prob=prob, target=0.3, nev=100)

#grid can be increase in resolution.
ψgrid = init_grid(:ψ, 1000)
θgrid = init_grid(:sm, 2, start=1)
φgrid = init_grid(:sm, 1, start=-1)

grids = init_grids(ψgrid, θgrid, φgrid)

evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);

# Looking at the solution again

#first the continuum
continuum_plot(evals)

#looking at the same TAE
ind = find_ind(evals, 0.3)

#then look at the harmonics
harmonic_plot(ϕft, grids, ind)

contour_plot(ϕ, grids, ind)

surface_plot(ϕ, grids, ind)

