
#basis case of computing the continuum

using MID
using MIDViz
#%%

Nr = 100;
geo = init_geo(R0=4.0)

prob = init_problem(q=fu_dam_q, geo=geo)#, met=no_delta_metric!); 
#%%
#then create the grids
Nr = 100;
rgrid = init_grid(type=:rc, N=Nr)
θgrid = init_grid(type=:as, N = 2, start = 1)
ζgrid = init_grid(type=:as, N = 1, start = -1)
grids = init_grids(rgrid, θgrid, ζgrid);
#%%
ω_cont = compute_continuum(prob, grids);
#%%
continuum_plot(ω_cont, grids)
