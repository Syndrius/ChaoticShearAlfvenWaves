

using MID
using Plots; plotlyjs()


#start very small, matrix scales much more extremly
Nr = 500;
Nθ = 20

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo, dens=axel_dens, δ=-4.0e-7); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=2)
ζgrid = init_sm_grid(start=-2, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)



ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=0.387);


ϕms = mode_structure(ϕ, grids);
reconstruct_continuum(ω, ϕms, grids)#, filename="data/fu_dam_spectrum.png")

#unsure why it is flipped... hopefully a resolution problemo
#this is an extremly different tae frequency... real good lol
tae_ind = find_ind(ω, 0.387)
tae_freq = ω[tae_ind]
plot_potential(ϕms, grids, tae_ind, 1)