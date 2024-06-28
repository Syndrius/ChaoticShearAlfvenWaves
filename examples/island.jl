
#basic test of island damping
using MID
using Plots; plotlyjs()


#1e-4 seems to be to big to make sense, this predicts no gap
#this could be because we don't have enough modes, hard to tell
#5e-5 still predicts a gap, but the upshift is significantly bigger than island_continuum predicts
#again this could be because of lack of modes, hard to tell.

N = 1000;
grids = init_grids(N=N, mstart=1, mcount=4, nstart=-2, ncount=1, nincr=4);

isl = IslandT(A=0.0e-4, m0=5, n0=4);
geo = GeoParamsT(R0=3.0);
prob = init_problem(q=island_damping_q, isl=isl, geo=geo, met=no_delta_metric!);

#ω_cont = continuum(prob=prob, grids=grids);
#plot_continuum(ω = ω_cont, grids=grids, n=-2)

#gapmin = maximum(ω_cont[:, 1])
#gapmax = minimum(ω_cont[:, 2])

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.35267/geo.R0)^2, nev=100);

#0.3846 is what island continuum predicts as upshift for A=1e-4
reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

tae_ind = find_ind(ω, 0.352)

#tae_freq=0.0014653584 may need to confirm this more accuratly!

display(ω[tae_ind]^2/geo.R0^2)

plot_potential(grids=grids, ϕ=ϕ, ind=tae_ind, n=1)


scatter(ones(length(ω)), real.(ω))


z = construct_surface(ϕ, length(ω), grids, 3);

plot_surface(z, grids, tae_ind)


island_width(isl, 5/4, 0.8)



N = 1000;
grids = init_grids(N=N, mstart=-3, mcount=12, nstart=-6, ncount=3, nincr=4);

isl = IslandT(A=2.5e-4, m0=5, n0=4);
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=island_damping_q, isl=isl, geo=geo, δ=-4e-7);

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=0.001465);


#reconstruct_continuum(ω = ω.^2, ϕ = ϕ, grids = grids)

#tae_ind = find_ind(ω.^2, geo.R0^2*0.001465)

#tae_freq=0.0014653584 may need to confirm this more accuratly!

display(ω[2])
display(imag(ω[1])/real(ω[1]))

plot_potential(grids=grids, ϕ=ϕ, ind=7, n=2)

#not really sure what this is telling is...
plot_sum_potential(grids=grids, ϕ=ϕ, ind=1)


z = construct_surface(ϕ, length(ω), grids, π/8);

#no notable island structure...
#may need to try this with higher res results.
#maybe we need to try and plot θ vs ζ at r=0.5?? unsure..
#or we need to think about what ζ values we should be using...
plot_surface(z, grids, 7)