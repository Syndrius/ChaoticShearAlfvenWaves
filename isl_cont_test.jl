

#see if we can actually undertsand the island continuum a bit.
#ideally inside and out...


using MID
using Plots#; plotlyjs()
using LaTeXStrings



#copied directly from earlier problemo.
#isl = IslandT(2, -1, 0.00015625000000000003, 2.0, 2.0, 0.5, 0.05)
isl = IslandT(2, -1, 0.00015625000000000003, 2.0, 2.0, 0.5, 0.05)
#I guess we should now consider our island...
#will need to understand wot the hek is going on with χ lol
#may also want R0 to be smaller so gap is larger

#w = 0.03 case, seems to be best resolved example.
isl = IslandT(2, -1, 5.625e-5, 2.0, 2.0, 0.5, 0.03)

#ie 1e-5 is the part we pass into our other code.
#so with 1e-6 we actually expect the tae to not experience damping
#but 1e-5 perhaps not? there is like 1 `strand` that it could interact with. But majority is below.
#similar with 5e-5, perhaps we just don't have enough modes?
#A = 4e-4 / 4 #4 to reflect r(1-r) factor

#so with small islands the structures are almost independent?

#educated guess lol.
#isl = ContIslandT(5, 4, A, 5/4, 0.8, 0.125)

geo = GeoParamsT(R0=1000.0)

#start with this???
pmd = asm_grid(start=-15, N=31, incr=1)
tmd = asm_grid(start=-5, N=11, incr=1)
#tmd = MID.ModeDataT(start=-8, count=10, incr=2)


#this makes far more sense lol.
κlist = LinRange(0.000001, 0.999, 300)

χlist = @. -(2*isl.A*κlist - isl.A)


ω2list = island_continuum(χlist, pmd, tmd, geo, isl, 0);


scatter(sqrt.(κlist), sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.9, ylimits=(-0.005, 0.08),xlabel=L"\sqrt{\kappa}", ylabel=L"\omega/\omega_A", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600)




scatter!(ones(size(reduce(vcat, ω2list))) .+ 0.05, sqrt.(abs.(reduce(vcat, ω2list) .* geo.R0^2)), label="Isl")

savefig("aapps_pics/isl_cont.png")


prob = init_problem(q=island_mode_21, geo=geo)

rgrid = MID.ContGridDataT(N=100)

grids = init_grids(rgrid, pmd, tmd)

evals_cont = cyl_cont(prob, grids);

ω_cont = continuum(prob, grids);

continuum_plot(evals_cont)


#v confused by this!
#island gaps do not seem to match...
#i guess we are considering the wrong mode numbers or whatevre.

continuum_plot(ω_cont, grids)





using Elliptic

xvals = LinRange(0, 0.8, 1000)

plot(xvals, Elliptic.K.(xvals))


#perhaps we can consider the zero -ellipticicty limit for the island to get a grasp of the continuum, this should be essentially the cylindrical limit of the island case, so same modes but without coupling...

function inside_island_q(r)
    om0 = sqrt(isl.A*isl.qp/(isl.q0^2*isl.r0))
    return -2*Elliptic.K(r) / (isl.m0*π * om0) #assume we don't need dq for cont??
end


#so now these match, we should be able to understand the island continuum...
#looks like this is basically just the normal cont in Qu, but coupling is very high, so gap width is huge.
#this seems more like a fact of the island q-profile rather than the island itself. But I guess we can always taylor expand q-profile, so as long as island is quite small, should be accurate
#so maybe the more meaningful question is why does the q-profile need to look like an elliptic K function??? pretty wild.

#defs m and m+2 modes are coupling, almost seems like an easier way to model this would be to do an elliptic thing with the above q-profile?? although I guess that is what island coords are???
#I think one difference is that the ellipticity would be a function of κ so this would probably be v difficult to do. 
#also, the contours are not actually ellipses, probbaly pretty close for κ ≈ 0, but would be a better approx than circles. Probably not much easier than just using the island coords tbh.

Nstart = 0
Ncount = 1
rgrid = MID.ContGridDataT(N=100)
θgrid = asm_grid(start=0, N=10)
ζgrid = asm_grid(start=Nstart, N=Ncount)
#sweeet!
φgrid = asm_grid(start=Nstart*isl.m0, N=Ncount, incr=isl.m0)#*isl.m0)


grids = init_grids(rgrid, θgrid, ζgrid)

prob = init_problem(q = inside_island_q, geo = geo)


evals_cont = cyl_cont(prob, grids);

continuum_plot(evals_cont, ymax=0.08, ymin=-0.01)


ω2list = island_continuum(χlist, θgrid, φgrid, geo, isl, 0);


scatter!(κlist, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.1)#, ylimits=(-0.01, 0.08))



κplot = LinRange(0, 1, 300)



plot(sqrt.(κplot), inside_island_q.(κplot) ./ inside_island_q(0), ylimits=(0.99, 2.51),xlabel=L"\sqrt{\kappa}", ylabel=L"q(\kappa)", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=24, xguidefontsize=25, xtickfontsize=10, ytickfontsize=10, dpi=600, legend=false)

savefig("aapps_pics/isl_q.png")