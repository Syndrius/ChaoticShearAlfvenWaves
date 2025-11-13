
using MID
using MIDViz
using Plots; plotlyjs()
#%%
#Integration has been massivly sped up, but this is still slow af, probbaly requires multiproc.
#~90% of the time is spent numerically integrating. Wonder if there is anything we can do???
Nr = 20
Nθ = 4
Nζ = 2
rgrid = init_grid(type=:rf, N=Nr);
θgrid = init_grid(type=:af, N=Nθ, pf=1);
ζgrid = init_grid(type=:af, N=Nζ, pf=-1);
grids = init_grids(rgrid, θgrid, ζgrid);
#%%

#first define the problem
geo = init_geo(R0=4.0)

prob = init_problem(q=fu_dam_q, geo=geo); 
#%%
#solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4, 0.5])
solver = init_solver(prob=prob, full_spectrum=true)
#%%

@time evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%
#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))

continuum_plot(evals, n=-1)


ind = find_ind(evals, 0.300)
#ind = 348
potential_plot(ϕft, grids, ind)

zeros(Matrix{3}{undef})

Array{2}(undef, 3, 3, 3, 3)
