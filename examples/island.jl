
#basic test of island damping
using MID
using Plots; plotlyjs()



geo = GeoParamsT(R0 = 1000)
rgrid = rfem_grid(N=100, start=0.0, stop=1.0, sep1=0.4, sep2=0.6, frac=0.5)#, left_bc=false)
θgrid = asm_grid(start=-2, N=5)
ζgrid = asm_grid(start=-1, N=3)

grids = init_grids(rgrid, θgrid, ζgrid)

isl = IslandT(m0=2, n0=-1, r0=0.5, w=0.03, qp=2.0)
#prob = init_problem(q = inside_island_q, geo=geo, met=MID.Geometry.Axel_island_metric!)
#prob = init_problem(q = inside_island_q, geo=geo, met=island_metric!)

prob = init_problem(q=island_mode_21, geo=geo, isl=isl)



evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);#, target_freq=10);

continuum_plot(evals)#, ymax=0.08)#, ymax



island_om = evals.ω[0.42 .< evals.r .< 0.58]


isl_ind = 31

#display(island_om[6]) #why the fk is there a huge imaginary one????
tae_ind = find_ind(evals, island_om[isl_ind])

potential_plot(ϕft, grids, tae_ind)


contour_plot(ϕ, grids, ind=tae_ind)
