
#test for newer q, designed to have tae with n=3, and island with m0=3, so that hopefully n_I=1 modes will overlap with tae.

using MID
using Plots; plotlyjs()

function nq(r)
    a = 1.4
    b = 6-4*a

    q = a + b*r^2
    dq = 2 * b * r

    return q, dq
end

#with individual parts.
#different ways of doing it are basically identical.
#10.992975 seconds (15.38 M allocations: 1.195 GiB, 4.40% gc time)

#start very small, matrix scales much more extremly
Nr = 50;
Nθ = 8

geo = GeoParamsT(R0=5.0)

isl = IslandT(A=0.0e-4, m0=3, n0=2);

prob = init_problem(q=nq, geo=geo, isl=isl); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=3)
#θgrid = init_sm_grid(start=1, count=5)
ζgrid = init_sm_grid(start=-3, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)

#cont_grids = init_grids(Nr, θgrid, ζgrid)

#ω_cont = continuum(prob, cont_grids)

#plot_continuum(ω_cont, cont_grids)

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);

ϕms = mode_structure(ϕ, grids)
reconstruct_continuum(ω, ϕms, grids)

tae_ind = find_ind(ω, 0.301)
tae_freq = ω[tae_ind]
plot_potential(ϕ, grids, tae_ind)


#perhaps with more res?

Nr = 80
Nθ = 30

geo = GeoParamsT(R0=5.0)

isl = IslandT(A=8.0e-4, m0=3, n0=2);

prob = init_problem(q=nq, geo=geo, isl=isl, δ=-4e-7); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=3)
#θgrid = init_sm_grid(start=1, count=5)
ζgrid = init_sm_grid(start=-3, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ = 0.301, nev=100);

ϕms = mode_structure(ϕ, grids)
reconstruct_continuum(ω, ϕms, grids)

#perf, in theory, island should overlap here!
#does seem to be no overlap, may a resolution issue, but stil quite surpriseing!
#maybe we need more toroidal mode numbers? island n should be 1, but we are only considering n=3??
tae_ind = find_ind(abs.(ω), 0.319)
tae_freq = ω[tae_ind]
plot_potential(ϕms, grids, tae_ind)

display(ω[1:5])

display(imag(ω[tae_ind])/real(ω[tae_ind]))

plot_phi_surface(ϕ, grids, tae_ind)