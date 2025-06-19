
#just a simple test case for inside the island, making sure it works at least a wee bit
#this is maybe just a test that we can compute the same spectrum based on all of our different islands!
#cool, confident this is working now!
#ideally we should be able to make it work with rad case interchangably
#we know we have rad working, and flux working so it should be ok
using MID
using MIDViz
using MIDIslands
#%%

geo = init_geo(R0=4.0)

#is flux island at 0.5 going to cause problemos?
#perhaps we should have done the flux case at 0.125? and gone to 0.5?
#island has no issue with it, which is probably all that matters.
isl21a = init_island(m0=2, n0=-1, w=0.1, ψ0=0.5, qp=2.0, flux=true)
#start with no islands
#prob = init_problem(geo=geo, q=island_q, met=:cylinder, isl=isl21a, type=:flux)
prob = init_problem(geo=geo, q=island_q, met=:island, isl=isl21a, type=:island)
#%%

κgrid = init_grid(type=:rc, N=100, start=0.0, stop=0.999)
ᾱgrid = init_grid(type=:as, N=3, start=-1)
τgrid = init_grid(type=:as, N=1, start=0)

grids = init_grids(κgrid, ᾱgrid, τgrid)
#%%

evals_cont = compute_continuum(prob, grids);
#%%
continuum_plot(evals_cont, grids, ymax=0.1)

#%%
#doesn't look good lol
chiland = MID.inst_island(isl21a)
chiland = PsiIslandT(2, -1, chiland.A, 2.0, 2.0, 0.5)
χlist = collect(chiland.A .- MID.inst_grid(κgrid) .* (2*chiland.A))

#lots of negatives lol
ω2list = island_continuum(χlist, ᾱgrid, τgrid, geo, chiland, 0)

#%%


ψ̄list = compute_ψ̄(chiland, χlist, 0);
r = repeat(4*sqrt.(2 .* ψ̄list), 1,  ᾱgrid.N * τgrid.N);
κplot = repeat(MID.inst_grid(κgrid), 1,  ᾱgrid.N * τgrid.N);

#%%

#other than x axis (again!) this actually looks like it matches. which is wild af.
#scatter!(r, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.5, ylimits=(0, 0.08))
scatter!(κplot, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.5, ylimits=(0, 0.08))
#%%

κgrid = init_grid(type=:rf, N=100, start=0.0, stop=0.999, left_bc=false)
ᾱgrid = init_grid(type=:as, N=3, start=-1)
τgrid = init_grid(type=:as, N=1, start=0)

grids = init_grids(κgrid, ᾱgrid, τgrid)
#%%
solver = init_solver(prob=prob, full_spectrum=true)
#%%
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%
#pretty fkn bang on, clearly some cooked af behaviour at κ=0, for m=0, perhaps the gradient needs to be zero? who knows.
#perhaps this is why the (0, 0) branch has been ignored!
continuum_plot(evals, ymax=0.1)
