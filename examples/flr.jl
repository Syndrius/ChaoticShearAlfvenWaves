
#test for flr effects

using MID
using Plots; plotlyjs()

geo = GeoParamsT(R0=12.5)

rgrid = init_fem_grid(N=200)
θgrid = init_sm_grid(start=2, count=2)
#θgrid = init_fem_grid(N=8, pf=2)
ζgrid = init_sm_grid(start=-2, count=1)

grids = init_grids(rgrid, θgrid, ζgrid)

#not sure why kwargs not working here but they do for geo...
flr = FLRT(δ=0.0, ρ_i=0.002, δ_e=0.1)

prob = init_problem(q=flr_q, geo=geo, flr=flr)

evals, ϕ, ϕft = compute_spectrum(grids=grids, prob=prob, full_spectrum=true);

#display(size(ϕ))
#reconstruct_continuum(ω, ϕ, grids)
plot_continuum(evals)

#display(grids.θ.count)


tae_ind = find_ind(evals, 0.4260)

evals.ω[tae_ind]

plot_potential(ϕft, grids, tae_ind)

(evals.ω[tae_ind])^2

