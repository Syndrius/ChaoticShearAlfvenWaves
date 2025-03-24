

#this works now (well enough! will need some serious cleaning, also unsure how this will be put into MID. Should be possible,
#but we may need to consider a separate package.


using MID
using Plots; plotlyjs()

#fss seems to be incapable of finding any gap modes???
#bit odd.


geo = GeoParamsT(R0 = 1000)
rgrid = rfem_grid(N=80, start=0.0, stop=0.999, left_bc=false)
#θgrid = asm_grid(start=-4, N=9)
θgrid = afem_grid(N=6, pf=1)
ζgrid = asm_grid(start=-1, N=2)

grids = init_grids(rgrid, θgrid, ζgrid)

isl = IslandT(m0=2, n0=-1, r0=0.5, w=0.03, qp=2.0)
#prob = init_problem(q = inside_island_q, geo=geo, met=MID.Geometry.Axel_island_metric!)
#prob = init_problem(q = inside_island_q, geo=geo, met=island_metric!)

prob = MID.Structures.init_isl_problem(geo=geo, isl=isl)


#so ϕ for fss is not working... gee wiz.
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);#, target_freq=10);

continuum_plot(evals, ymax=0.08, ymin=-0.01)#, ymax=10)#, n=-2)

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
