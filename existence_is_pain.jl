
#can we atleast get the same results as a simple case of Axel's
using MID

using Plots; plotlyjs()


#start very small, matrix scales much more extremly
Nr = 100;

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo)#, met=no_delta_metric!); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=2, count=2)
ζgrid = init_sm_grid(start=-2, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)


evals, ϕ, ϕ_ft = compute_spectrum(grids=grids, prob=prob, full_spectrum=true);

plot_continuum(evals)

ind = find_ind(evals, 0.11075)

plot_potential(ϕ_ft, grids, ind)



#this took forever lol.
Nr = 80;
Nθ = 15

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo)#, met=no_delta_metric!); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ)
ζgrid = init_sm_grid(start=-2, count=1)
grids = init_grids(rgrid, θgrid, ζgrid)


evals, ϕ, ϕ_ft = compute_spectrum(grids=grids, prob=prob, full_spectrum=true);

plot_continuum(evals)

ind = find_ind(evals, 0.11075)

plot_potential(ϕ_ft, grids, ind)