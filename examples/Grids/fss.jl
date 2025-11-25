using MID
using MIDViz
#%%
#this example may be able to cover all the grid cases.

#met = MID.Geometry.flux_toroidal_metric!
geo = init_geometry()#, met)
fields = init_fields()

prob = init_problem(geometry=geo, fields=fields)
#%%
x1grid = init_grid(:ψ, 100)
x2grid = init_grid(:sm, 2, start=1)
x3grid = init_grid(:spectral, 1, start=-1) #perhaps remove start from kwargs?

grids = init_grids(x1grid, x2grid, x3grid)
#%%

#prob here is suboptimal!
#ideally, the normalisation is handled somewhere else.
#perhaps solver can have targets and normalised targets? or something?
#we defs have to do something about this!
solver = init_solver(prob=prob, full_spectrum=true)
#%%

@time evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);

continuum_plot(evals)
