
#basic test of island damping
using MID

N = 100;
grids = init_grids(N=N, mstart=2, mcount=2, nstart=-2, ncount=1);

isl = IslandT(A=0e-4, m0=5, n0=4);
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=island_damping_q, isl=isl, geo=geo);

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);


reconstruct_continuum(ω = ω.^2, ϕ = ϕ, grids = grids)

tae_ind = find_ind(ω.^2, geo.R0^2*0.001465)

#tae_freq=0.0014653584 may need to confirm this more accuratly!

display(ω[tae_ind])

plot_potential(grids=grids, ϕ=ϕ, ind=tae_ind, n=1)