using MID
using MIDViz
#%%
#don't think this should be its own file.
#move to basic!
#ideally we could specify cyl with :ψ or whatever we actually use!
geo = init_geometry(:cyl)

#need to change this to allow a single isl to be input.
fields = init_fields(:r, q=gae_q, dens=gae_dens)
#%%

#should make this grid auto change the start/stop.
rgrid = init_grid(:r, 60)
θgrid = init_grid(:sm, 3, start=0)
ζgrid = init_grid(:sm, 2, start=0)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
dens_plot(gae_dens, grids)
q_plot(gae_q, grids)
#%%
prob = init_problem(geometry=geo, fields=fields)
#island has not been instantiated properly!
#so metric is cooked af.
#%%
solver = init_solver(prob=prob, full_spectrum=true)
#%%
evals, ϕ, ϕft = compute_spectrum(prob=prob, solver=solver, grids=grids);
#%%
#so we can see basic global mode structure, I guess that is all that matters tbh!
continuum_plot(evals)


ind = find_ind(evals, 0.6)

harmonic_plot(ϕft, grids, ind)

