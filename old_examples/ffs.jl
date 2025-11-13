
using MID
using MIDViz
#%%

#first define the problem
geo = init_geo(R0=4.0)

#rgrid = collect(LinRange(0, 1, N));
prob = init_problem(q=fu_dam_q, geo=geo); 
#%%
#then create the grids
Nr = 30;
Nθ = 5;
rgrid = init_grid(type=:rf, N=Nr)
θgrid = init_grid(type=:af, N=Nθ)
ζgrid = init_grid(type=:as, N = 1, start = -1)
grids = init_grids(rgrid, θgrid, ζgrid);

#%%

solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4, 0.5])
#%%
@time evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);

#%%
#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))
continuum_plot(evals)


ind = find_ind(evals, 0.302)
potential_plot(ϕft, grids, ind)
