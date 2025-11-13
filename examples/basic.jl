using MID
using Plots
#%%
#this example may be able to cover all the grid cases.

#met = MID.Geometry.flux_toroidal_metric!
geo = MID.Geometry.init_geometry(4.0)#, met)

q = MID.Fields.fu_dam_q
dens = MID.Fields.uniform_dens
isls = MID.Fields.IslandT[]
#fields = init_fields(q=q, dens=dens, isls=isls)
fields = init_fields()

prob = init_problem(geo=geo, fields=fields)
#%%

#the way this is done is cooked beyond beleif.
#think we should have an init_fem_grid, init_spectral_grid and init_grid which defaults to fem.
#and make the options more like the newer version.
x1grid = init_grid(:ψ, 50)
#x2grid = init_grid(:θ, 4)
x2grid = init_grid(:sm, 3, start=1)
#x3grid = init_grid(:ζ, 2)
x3grid = init_grid(:spectral, 2, start=-2) #perhaps remove start from kwargs?

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
