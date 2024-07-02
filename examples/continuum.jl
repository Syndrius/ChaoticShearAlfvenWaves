
#basis case of computing the continuum

using MID
using Plots; plotlyjs()


Nr = 100;
geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo)#, met=no_delta_metric!); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=2, count=2)
ζgrid = init_sm_grid(start=-2, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)


ω_cont = continuum(prob = prob, grids=grids);

plot_continuum(ω=ω_cont, grids=grids)

