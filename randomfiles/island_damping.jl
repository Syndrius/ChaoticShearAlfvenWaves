
#return to damping just in case in randomly works again
using MID
using MIDViz
using Plots; plotlyjs()
#%%

function damp_q(r::Float64)
    a = 1.1
    b = 0.625 #this looks ok, will have a 3/2 island at 0.8 for damping comparison.
    return a+b*r^2, 2*b*r
end
#%%

rgrid = init_grid(type=:rf, N=80)
θgrid = init_grid(type=:as, N=5, start=1)
ζgrid = init_grid(type=:as, N=2, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
geo = init_geo(R0=4.0)
prob = init_problem(q = damp_q, geo=geo)
solver = init_solver(prob=prob, full_spectrum=true)
#%%

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);


#%%

continuum_plot(evals)
#%%
ind = find_ind(evals, 0.2794)

potential_plot(ϕft, grids, ind)


