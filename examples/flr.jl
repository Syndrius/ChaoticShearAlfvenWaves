
#test for flr effects

using MID
using Plots; plotlyjs()

geo = GeoParamsT(R0=12.5)

rgrid = init_fem_grid(N=40)
#θgrid = init_sm_grid(start=2, count=2)
θgrid = init_fem_grid(N=8, pf=2)
ζgrid = init_sm_grid(start=-2, count=1)

grids = init_grids(rgrid, θgrid, ζgrid)

#not sure why kwargs not working here but they do for geo...
flr = FLRT(δ=0.0, ρ_i=1e-12, δ_e=0.1)

prob = init_problem(q=flr_q, geo=geo, flr=flr)

ω, ϕ = construct_and_solve(grids=grids, prob=prob, full_spectrum=true);

#display(size(ϕ))
reconstruct_continuum(ω, ϕ, grids)

#display(grids.θ.count)

display(ω[1:4])

tae_ind = find_ind(real.(ω), 0.389)

ω[tae_ind]

ϕms = mode_structure(ϕ, grids);

plot_potential(ϕms, grids, tae_ind)