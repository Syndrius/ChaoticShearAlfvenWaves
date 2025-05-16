
using MID
using MIDViz
#%%
#Integration has been massivly sped up, but this is still slow af, probbaly requires multiproc.
#~90% of the time is spent numerically integrating. Wonder if there is anything we can do???
Nr = 30
Nθ = 6
Nζ = 2
rgrid = init_grid(type=:rf, N=Nr);
θgrid = init_grid(type=:af, N=Nθ, pf=2);
ζgrid = init_grid(type=:af, N=Nζ, pf=-2);
grids = init_grids(rgrid, θgrid, ζgrid);
#%%

#first define the problem
geo = init_geo(R0=4.0)

prob = init_problem(q=fu_dam_q, geo=geo); 
#%%
solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4, 0.5])
#%%

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%
#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))
continuum_plot(evals)


ind = find_ind(evals, 0.33)
#ind = 348
potential_plot(ϕft, grids, ind)
