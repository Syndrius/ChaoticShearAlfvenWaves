
using MID
using Plots; plotlyjs()





Nr = 100;
#a = 0.523 -> rescale R0 from 5.71 to equivalent with a=0
geo = GeoParamsT(R0=1/0.523 * 5.71)

prob = init_problem(q=Axel_island_q, geo=geo)#, met=no_delta_metric!); 

#rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=0, count=10)
ζgrid = init_sm_grid(start=-8, count=8)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(Nr, θgrid, ζgrid)


ω_cont = continuum(prob, grids);

plot_continuum(ω_cont, grids, ymax=1)


using MIDViz


isl = IslandT(m0=5, n0=-5, A=0.001)


Ntraj = 40
flux_list = LinRange(0.01, 0.5, Ntraj)
rlist = @. sqrt(2*flux_list)
rp, θp = poincare_plot(prob.q, slab_to_plot, 500, Ntraj, 0, 0.0, 0.0, 0.0, prob.geo.R0, isl, rlist);

#Axel's island is ~4.3cm full width.
0.043/0.523
#scales to 0.08 when a-> 1.
#A=0.001 seems to be a pretty darn good estimate lol.
#may want to be more precise.
#our width is 