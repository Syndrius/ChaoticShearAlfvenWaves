
using MID
using MIDViz
#%%
#first define the problem
geo = init_geo(R0=4.0)

#rgrid = collect(LinRange(0, 1, N));
prob = init_problem(q=MID.Equilibrium.flux_fu_dam_q, geo=geo, type=:flux); 
#%%
#then create the grids
Nr = 30;
rgrid = init_grid(type=:rf, N=Nr, stop=0.5)
θgrid = init_grid(type=:as, N = 2, start = 1)
ζgrid = init_grid(type=:as, N = 1, start = -1)
grids = init_grids(rgrid, θgrid, ζgrid);

#%%
#then define the solver
solver = init_solver(full_spectrum=true, prob=prob)
#solver = init_solver(target=0.33, prob=prob)
#solver = init_solver(prob=prob, nshifts=4)
#solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4, 0.5])
#%%

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=false);
#%%
#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))
continuum_plot(evals)


#freq is a bit different, but otherwise this is the same
ind = find_ind(evals, 0.284)
#ind = 348
potential_plot(ϕft, grids, ind)
