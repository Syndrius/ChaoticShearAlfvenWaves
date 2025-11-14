using MID
using Plots; plotlyjs()
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
