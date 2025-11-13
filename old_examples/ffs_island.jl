
using MID

using Plots; plotlyjs()


#1e-4 seems to be to big to make sense, this predicts no gap
#this could be because we don't have enough modes, hard to tell
#5e-5 still predicts a gap, but the upshift is significantly bigger than island_continuum predicts
#again this could be because of lack of modes, hard to tell.

Nr = 600;
Nθ = 15

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=2)
ζgrid = init_sm_grid(start=-2, count = 1)
grids = init_grids(rgrid, θgrid, ζgrid);

isl = IslandT(A=0.0e-4, m0=5, n0=4);
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=island_damping_q, isl=isl, geo=geo);

#ω_cont = continuum(prob=prob, grids=grids);
#plot_continuum(ω = ω_cont, grids=grids, n=-2)

#gapmin = maximum(ω_cont[:, 1])
#gapmax = minimum(ω_cont[:, 2])

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=0.389, nev=100);


ϕms = mode_structure(ϕ, grids);
reconstruct_continuum(ω, ϕms, grids)

tae_ind = find_ind(ω, 0.389)
tae_freq = ω[tae_ind]
plot_potential(ϕms, grids, tae_ind, 1)