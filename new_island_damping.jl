
using MID
using Plots; plotlyjs()


Nr = 100;
Nθ = 15

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=2)
#θgrid = init_sm_grid(start=2, count = 2)
ζgrid = init_sm_grid(start=-3, count = 3)
grids = init_grids(rgrid, θgrid, ζgrid);

#no tangible interaction between tae and island... yay :(   
isl = IslandT(A=3.0e-4, m0=3, n0=2);
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=test_q, isl=isl, geo=geo);

#ω_cont = continuum(prob=prob, grids=grids);
#plot_continuum(ω = ω_cont, grids=grids, n=-2)

#gapmin = maximum(ω_cont[:, 1])
#gapmax = minimum(ω_cont[:, 2])

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=0.387, nev=50);


ϕms = mode_structure(ϕ, grids);
reconstruct_continuum(ω, ϕms, grids)

tae_ind = find_ind(ω, 0.387)

tae_freq = ω[tae_ind]
plot_potential(ϕms, grids, 2, 3)