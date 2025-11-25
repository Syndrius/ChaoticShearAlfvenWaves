using MID
#unsure if this is better
#it may even be worse.

"""
This main function, compute_spectrum, requires 3 inputs.

The first defines the problem to be solved and is requires information about the geometry and the fields.
"""

#by default, the geometry is initialised as toroidal with R0=4.0
geo = init_geometry()

#by default, the fields are intialised with a quadratic q-profile
#and a uniform density using the toroidal flux, ψ, as the radial coordinate
fields = init_fields()

#we can then create the problem from these.
prob = init_problem(geometry=geo, fields=fields)

"""
We next define the grids used.

These are defined in each dimension separately, then combined into a 3d grid.
The type of grids constructed dictate the method of solution.
"""

#the radial grid, :ψ specifies the type, and 100 is the number of points.
ψgrid = init_grid(:ψ, 100)
#:sm specifies that the Fourier Spectral Method is used in the poloidal direction
#with 2 modes, m=1 and m=2.
θgrid = init_grid(:sm, 2, start=1)
#n=-1
φgrid = init_grid(:sm, 1, start=-1)

#we then combine the grids
grids = init_grids(ψgrid, θgrid, φgrid)

"""
Finally, we specify the solver to use.
"""

#the solver is set to compute the entire spectrum
#which is only practical for small scale problems.
#the prob input is needed here for normalisation.
solver = init_solver(prob=prob, full_spectrum=true)

"""
We now compute the spectrum.
"""

#this returns an object storing the eigenvalue information
#the perturbation to the electric potential
#and the Fourier transform in θ and φ of this potential.
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);

"""
The companion Package CSAWViz allows for simple plotting.
"""
using MIDViz

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

"""
Now the TAE frequncy is known, we can change the solver to target this specific frequency.
This allows us to increase the resolution and solve much larger matrices.
"""

#sets the solver to perform shift and inert, getting the 100 solutions closets to 0.3.
solver = init_solver(prob=prob, target=0.3, nev=100)

#grid can be increase in resolution.
ψgrid = init_grid(:ψ, 1000)
θgrid = init_grid(:sm, 2, start=1)
φgrid = init_grid(:sm, 1, start=-1)

grids = init_grids(ψgrid, θgrid, φgrid)

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);


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
#%%
#Below we should move to the grids.jl
#TODO
#example of interpolation in here.
"""
By changing the radial grid, we can compute only the continuum.
"""

ψgrid_c = init_grid(:cont, 100)
grids_c = init_grids(ψgrid_c, θgrid, φgrid)

#ok so this can't have no kwargs, that is silly af.
evals_c = compute_spectrum(prob, grids_c);

continuum_plot(evals_c)

"""
Or we can change the poloidal and/or the toroidal grid to use Finite Elements in every dimension.
Note that this is more demanding.
"""

ψgrid = init_grid(:ψ, 30)
#setting phase factor (pf) focusses on a specific branch, giving more accuracy.
θgrid = init_grid(:θ, 5, pf=1)
φgrid = init_grid(:φ, 2, pf=-1)

grids = init_grids(ψgrid, θgrid, φgrid)

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);

continuum_plot(evals)

ind = find_ind(evals, 0.31)

harmonic_plot(ϕft, grids, ind)
#%%
#should demonstrate basic functionality,
#perhaps this case should consider different geometries and fields?
#who knows.
#not sure if this shoudl cover all the different scenarious we can have?
#perhaps instead of basic and ffs etc
#we can have a basic one that just gets this working

#then a general one that covers heaps of different cases.



#this example may be able to cover all the grid cases.

#met = MID.Geometry.flux_toroidal_metric!
geo = init_geometry()#, met)

fields = init_fields(:r)

prob = init_problem(geometry=geo, fields=fields)
#%%

#the way this is done is cooked beyond beleif.
#think we should have an init_fem_grid, init_spectral_grid and init_grid which defaults to fem.
#and make the options more like the newer version.
x1grid = init_grid(:ψ, 15)
x2grid = init_grid(:θ, 3, pf=1)
#x2grid = init_grid(:sm, 2, start=1)
x3grid = init_grid(:ζ, 1, pf=-1)
#x3grid = init_grid(:spectral, 1, start=-1) #perhaps remove start from kwargs?

grids = init_grids(x1grid, x2grid, x3grid)
#%%

#prob here is suboptimal!
#ideally, the normalisation is handled somewhere else.
#perhaps solver can have targets and normalised targets? or something?
#we defs have to do something about this!
solver = init_solver(prob=prob, full_spectrum=true)
#%%

@time evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);
scatter(evals.x1, real.(evals.ω), ylimits=(0.0, 1.0))

find_ind(evals, 0.301)
