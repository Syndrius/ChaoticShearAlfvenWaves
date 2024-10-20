
using MID
using Plots; plotlyjs()



geo = GeoParamsT(R0 = 10)
rgrid = rfem_grid(N=60)
θgrid = asm_grid(start=0, N=3)
ζgrid = asm_grid(start=-1, N=2)

grids = init_grids(rgrid, θgrid, ζgrid)

prob = init_problem(q = island_mode_21, geo=geo)#, met=cylindrical_metric!)

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

continuum_plot(evals)#, n=-2)

#println(real.(evals.ω[91:100]))

#so to summarise, gams do whatever the fk they want.
potential_plot(ϕft, grids, ind)

ind = find_ind(evals, 0.14737)

potential_plot(ϕft, grids, ind)