
using MID
using MIDIslands
using Plots
using MIDViz

#%%

isl = init_island(m0=4, n0=-3, w=0.1, qp=2.0, ψ0=0.5, coords=true)
isl = MID.inst_island(isl)
met = MID.MetT()
#%%

#ok, we are computing, probs garbage.
#but this is enough to check if the continuum is computed accuratly
#jks need Jacobian, and maybe its deriv?? -> if we don't need deriv we can get away with taking the determinant
#unsure what the continuum calculation actually requires.
MID.Geometry.island_metric!(met, 0.001, π/2, 0.0, 1.0, isl)

met.gu
#%%
#compare continuums.

#so N=4 for the poloidal grid just cooks this...
θgrid = init_grid(type=:as, N=4, start=-1)
ζgrid = init_grid(type=:as, N=1, start=0)

κlist = LinRange(0.001, 0.999, 100)
χlist = @. isl.A - κlist * 2 * isl.A 

ψland = PsiIslandT(isl.m0, isl.n0, isl.A, isl.q0, isl.qp, isl.ψ0)

geo = init_geo(R0=1.0)

quom = island_continuum(χlist, θgrid, ζgrid, geo, ψland, 0)

κp = repeat(κlist, 1, 3)
scatter(κp, sqrt.(abs.(quom .* geo.R0^2)), legend=false, markersize=1.0)#, ylimits=(0, 0.08))

#%%
κc = init_grid(type=:rc, N=10, start=0.001, stop=0.999)
grids = init_grids(κc, θgrid, ζgrid)

prob = init_problem(geo=geo, q=island_q, type=:island, isl=isl)

evals = compute_continuum(prob, grids)

continuum_plot(evals)

scatter!(κp, sqrt.(abs.(quom .* geo.R0^2)), legend=false, markersize=1.0)#, ylimits=(0, 0.08))


