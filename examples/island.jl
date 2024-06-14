
#basic test of island damping
using MID

N = 300;
grids = init_grids(N=N, mstart=2, mcount=2, nstart=-2, ncount=1);

isl = IslandT(A=0.0e-4, m0=5, n0=4);
geo = GeoParamsT(R0=3.0);
prob = init_problem(q=island_damping_q, isl=isl, geo=geo);

ω_cont = continuum(prob=prob, grids=grids);
plot_continuum(ω = ω_cont, grids=grids, n=-2)

gapmin = maximum(ω_cont[:, 1])
gapmax = minimum(ω_cont[:, 2])

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);


reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

tae_ind = find_ind(ω, 0.449)

#tae_freq=0.0014653584 may need to confirm this more accuratly!

display(ω[tae_ind])

plot_potential(grids=grids, ϕ=ϕ, ind=tae_ind, n=1)




N = 1000;
grids = init_grids(N=N, mstart=-3, mcount=12, nstart=-6, ncount=3, nincr=4);

isl = IslandT(A=1.0e-5, m0=5, n0=4);
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=island_damping_q, isl=isl, geo=geo, δ=-4e-7);

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=0.001465);


#reconstruct_continuum(ω = ω.^2, ϕ = ϕ, grids = grids)

#tae_ind = find_ind(ω.^2, geo.R0^2*0.001465)

#tae_freq=0.0014653584 may need to confirm this more accuratly!

display(ω[2])
display(imag(ω[1])/real(ω[1]))

plot_potential(grids=grids, ϕ=ϕ, ind=1, n=1)