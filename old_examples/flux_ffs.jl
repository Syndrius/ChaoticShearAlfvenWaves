
using MID
using MIDViz
using Plots; plotlyjs()
#%%

#first define the problem
geo = init_geo(R0=4.0)

#rgrid = collect(LinRange(0, 1, N));
prob = init_problem(q=MID.Equilibrium.flux_fu_dam_q, geo=geo, type=:flux); 
#%%
#then create the grids
Nr = 20;
Nθ = 4;
rgrid = init_grid(type=:rf, N=Nr, stop=0.5)
θgrid = init_grid(type=:af, N=Nθ, pf=1)
ζgrid = init_grid(type=:as, N = 1, start = -1)
grids = init_grids(rgrid, θgrid, ζgrid);

#%%

solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4, 0.5])
solver = init_solver(prob=prob, full_spectrum=true)
#%%
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);

#%%
#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))
continuum_plot(evals)

#still sort of works, but not as well
ind = find_ind(evals, 0.294)
evals.ω[ind]
potential_plot(ϕft, grids, ind)
