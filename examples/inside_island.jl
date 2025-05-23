

#this works now (well enough! will need some serious cleaning, also unsure how this will be put into MID. Should be possible,
#but we may need to consider a separate package.
#unsure if this functionailty should be moved to MIDIslands
#think either that or the creation of the problem etc needs to be a bit more streamlined
#perhaps we just pick the island metric then the q-profile doesn't matter?
#and that creates an island problemo?


using MID
using MIDViz
using Plots; plotlyjs()
#%%
#fss seems to be incapable of finding any gap modes???
#bit odd.


geo = init_geo(R0 = 10.0)
rgrid = init_grid(type=:rf, N=80, start=0.0, stop=0.999, left_bc=false)
#θgrid = asm_grid(start=-4, N=9)
θgrid = init_grid(type=:as, N=6, start=-2)
ζgrid = init_grid(type=:as, start=-1, N=2)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
k = 0.002
#think this matches the island_mode_21, we will need further verification.
isl = init_island(m0=2, n0=-1, A=k/2, r0=0.5, qp=2.0)
#prob = init_problem(q = inside_island_q, geo=geo, met=MID.Geometry.Axel_island_metric!)
#prob = init_problem(q = inside_island_q, geo=geo, met=island_metric!)

prob = MID.Structures.init_isl_problem(geo=geo, isl=isl)
display(prob.isl)

MID.Structures.inst_island(isl)
#%%

solver = init_solver(prob=prob, full_spectrum=true)
#%%

#so ϕ for fss is not working... gee wiz.
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);#, target_freq=10);
#%%
continuum_plot(evals, ymax=0.3, ymin=-0.01)#, ymax=10)#, n=-2)
#%%
#obvs not working, check equivalence between ours and Axel, by considering coord transformation from κ to κ^2
#also may need to check elliptic functions are doing what we want...

#(-1, 1) modes seem to be perfectly coupled, bit of a problemo..
#perhaps this is the expected behaviour????
#maybe this will change with fff???
#based on gae testing, this may just be a problem with fss.

ind = find_ind(evals, 0.00665294)

ind = 234

potential_plot(ϕft, grids, ind)

contour_plot(ϕ, grids, ind=ind)


#current hard codes island for this case.
isl = IslandT(2, -1, 0.00015625000000000003, 2.0, 2.0, 0.5, 0.05)

#pmd = asm_grid(start=-12, N=26, incr=1)
#tmd = asm_grid(start=-2, N=5, incr=1)
#tmd = MID.ModeDataT(start=-8, count=10, incr=2)


#this makes far more sense lol.
κlist = LinRange(0.000001, 0.999, 100)

χlist = @. -(2*isl.A*κlist - isl.A)



ω2list = island_continuum(χlist, θgrid, ζgrid, geo, isl, 0);


scatter!(κlist, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.1)#, ylimits=(0.3, 0.5))#.2, 0.6))
#3.997713 = 4
#0.0049998 = 0.005
