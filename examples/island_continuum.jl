

using MID
using Plots; plotlyjs()

display(prob.isl)

#copied directly from earlier problemo.
isl = IslandT(2, -1, 0.00015625000000000003, 2.0, 2.0, 0.5, 0.05)
#I guess we should now consider our island...
#will need to understand wot the hek is going on with χ lol
#may also want R0 to be smaller so gap is larger

#ie 1e-5 is the part we pass into our other code.
#so with 1e-6 we actually expect the tae to not experience damping
#but 1e-5 perhaps not? there is like 1 `strand` that it could interact with. But majority is below.
#similar with 5e-5, perhaps we just don't have enough modes?
#A = 4e-4 / 4 #4 to reflect r(1-r) factor

#so with small islands the structures are almost independent?

#educated guess lol.
#isl = ContIslandT(5, 4, A, 5/4, 0.8, 0.125)

geo = GeoParamsT(R0=10.0)

#start with this???
pmd = asm_grid(start=-12, N=26, incr=1)
tmd = asm_grid(start=-2, N=5, incr=1)
#tmd = MID.ModeDataT(start=-8, count=10, incr=2)


#this makes far more sense lol.
κlist = LinRange(0.000001, 0.999, 100)

χlist = @. -(2*isl.A*κlist - isl.A)

#χlist = LinRange(-A+A*0.05, A-A*0.05, 50) #this is fked, think we need to cluster near the spratrix.

#not sure how we will get this to conform to our other work....
#χlist1 =  1 .- 0.01 * exp.( -1 .* collect(LinRange(0, 12, 11)))
#χlist2 = LinRange(0, sqrt(χlist1[1]), 191)[1:end-1] .^2 .+ 1e-5
#ok this matches Zhisong, again no fkn idea what the hell this list is.
#χlist = isl.A .- vcat(χlist2, χlist1) .* 2 .* isl.A

ω2list = island_continuum(χlist, pmd, tmd, geo, isl, 0);

#width = 4 * sqrt(isl.A * isl.q0^2/isl.qp)
#ψ_isl = 2 * width / (π * isl.m0)
#ψ̄m = ψ_isl
#ψ̄list = MID.Continuum.compute_ψ̄(isl, χlist, 0)

#κlist = @. (-χlist + isl.A) / (2 * isl.A)

#this is such a fkn awful x axis. by gollly this will be annoying to change.
#ideally our x-axis will be κ, ideally χ would die lol.
#r = repeat(sqrt.(2 .* ψ̄list), 1,  pmd.count * tmd.count)

#we should normalise the x-axis, it is going from centre of island to edge I think.

#ok so this kind of matches our case a bit.
#we probably need to match our island cases to island coords to see what is happening.
#probably requires that we fix up the other stuff.
#we will probably need to do it for inside and outside separatly.
scatter(κlist, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.1)#, ylimits=(0.3, 0.5))#.2, 0.6))

hline!([0.425, 0.425]) #top, 
hline!([0.3764, 0.3764]) #bottom
hline!([0.380, 0.380]) #tae

savefig("data/9e-5_continuum.png")
#our gap has max at 0.424
#min at 0.3764
#tae at ~0.380




A = 1e-4

isl = ContIslandT(5, 2, A, 5/2, 4.0, 0.125)

geo = GeoParamsT(R0=10.0)


pmd = MID.ModeDataT(start=-10, count=2, incr=1)
tmd = MID.ModeDataT(start=-6, count=2, incr=2)

#@. doesn't work with exp of linrange for some reason!
χlist1 =  1 .- 0.01 * exp.( -1 .* collect(LinRange(0, 12, 2)))[1:end-1]
χlist2 = LinRange(0, sqrt(χlist1[1]), 2)[1:end-1] .^2 .+ 1e-5
#ok this matches Zhisong, again no fkn idea what the hell this list is.
χlist = A .- vcat(χlist2, χlist1) .* 2 .* A

#may actually be different!
#not ideal, will need to test with less modes. takes like ~5mins.
ω2list = island_continuum(χlist, pmd, tmd, geo, isl, 0);

ω2list_old = passing_continuum(χlist, pmd, tmd, geo, isl);

#seeing tiny differences, I think form actually computing gu vs taking inverse.
display(ω2list)# .- 
display(ω2list_old)



χlistplus = @. -A - LinRange(0, 1, 401)[2:end]^2 * 0.005
#χlistplus = @. -noA - LinRange(0, 1, 2)[2:end]^2 * 0.044
χlistminus = @. -A - LinRange(0, 1, 401)[2:end]^2 * 0.003




#I guess we need to export modedataT..
pmd = MID.ModeDataT(start=0, count=20, incr=1)
tmd = MID.ModeDataT(start=-5, count=10, incr=2)

ω2list_plus_old = passing_continuum(χlistplus, pmd, tmd, geo, isl, 1);
#ω2list_minus_old = passing_continuum(χlistminus, pmd, tmd, geo, isl, -1);


ω2list_plus = island_continuum(χlistplus, pmd, tmd, geo, isl, 1);
#ω2list_minus = island_continuum(χlistminus, pmd, tmd, geo, isl, -1);


#thinkt they are the same now, hard to be certain though!
display(ω2list_plus .- ω2list_plus_old)