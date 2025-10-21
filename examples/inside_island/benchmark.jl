
#replicating the benchmark from Zhisongs paper, figure 2.
using MID
using MIDIslands
using Plots
using MIDViz

#%%
isl = init_island(m0=5, n0=-2, A=1e-5, qp=4.0, ψ0=0.125, coords=true)
geo = init_geo(R0=1000.0)
prob = init_problem(type=:island, geo=geo, isl=isl, q=island_q)
#%%

#χlist = - isl.A + 0.01 * isl.A * np.exp.(-1 * LinRange(0, 16, 100))
κgrid = init_grid(type=:rc, start=0.995, stop=0.999999, N=100)
ᾱgrid = init_grid(type=:as, start=-7, N=15)
τgrid = init_grid(type=:as, start=0, N=1)

grids = init_grids(κgrid, ᾱgrid, τgrid)
#%%

evals = compute_continuum(prob, grids)
#%%

continuum_plot(evals, ylimits=(-0.01, 0.03))

